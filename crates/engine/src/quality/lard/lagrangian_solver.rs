// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/quality/lard/LagrangianSolver.hpp and reactSpeciesBlock
// from src/engine/transport/components/ReactionModule/ReactionLegacyBinding.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin,
// licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The LAGRANGIAN quality engine (QUALITY_SOLVER LAGRANGIAN): per-link segment transport. Each
// substep ages the parcels, drains every conduit's back into its downstream node ledger, mixes
// the nodes in topological order (with [TREATMENT] and the zero-volume passthrough of the
// non-conduit links), releases new front segments, applies the RWPT dispersion, decays and
// reacts; the routing step then applies the heat fluxes and the bed and publishes the conduit
// means.
//
// Differences in form, never in result:
//   * step() returns Result: the reaction stage, the booster sources and the heat flux modules
//     can meet code not translated yet (StepError::NotTranslated).
//   * reactSpeciesBlock's process-wide static BindingScratch (the integrator workspace and the
//     gather block) is a per-call local: every block entry is written before it is read and
//     the integrator writes every workspace entry it reads.
//   * thermalRegions, speciesRegions, applyThermalIncrements, applySpeciesIncrements and the
//     checkpoint serializers serve the C API's state accessors and engine checkpoints, which
//     are not on the 1D run path; they are not translated.
//
// NOT YET TRANSLATED (StepError::NotTranslated where met; SWMMEngine::stepRouting refuses such
// a deck first): applyFluxesAndBed past its guard (the surface, radiative and bed heat fluxes
// and the bed solute exchange), which needs the heat flux modules.

use crate::core::cmath::{cmax, cvt_i32};
use crate::core::context::SimulationContext;
use crate::core::engine_step_routing::StepError;
use crate::data::heat::HeatSource;
use crate::data::link::{LinkType, XsectShape};
use crate::data::node::NodeType;
use crate::data::reaction_tokens::RxHydVar;
use crate::data::water_age::WaterAgeSource;
use crate::hydraulics::node::{store_volume, store_volume_at};
use crate::quality::lard::rwpt_dispersion::{rwpt_hyd_radius, RwptDispersion, K_RWPT_PARTICLES_PER_LINK};
use crate::quality::lard::segment_store::SegmentStore;
use crate::quality::quality_routing::{apply_node_treatment, book_negative_source_clamp, node_is_reactor, update_untreated_storage_hrt, ZERO_VOLUME};
use crate::transport::heat_legacy::{apply_initial_temp_overrides, bed_exchange_enabled};
use crate::transport::reaction_legacy_binding::{ensure_msx_state, legacy_reactions_active};
use crate::transport::reaction_legacy_step::{clamp_rate_species, integrator_step, note_negative_clamp, reaction_area_volume, reaction_booster_mass, Ws};
use crate::transport::transport_policy::network_1d_enables;
use crate::transport::water_age_legacy::{apply_initial_age_overrides, book_negative_age_clamp};

/// cfs; below this a link moves nothing
pub const K_TINY_FLOW: f64 = 1.0e-8;

/// lard::SpeciesRowLayout: pollutants, then age, temperature and the reactions species.
#[derive(Debug, Clone, Copy)]
pub struct SpeciesRowLayout {
    pub np: i32,
    pub age_row: i32,
    pub temp_row: i32,
    pub msx_first: i32,
    pub ns: i32,
}

/// lard::rowLayout
pub fn row_layout(ctx: &SimulationContext) -> SpeciesRowLayout {
    let en = network_1d_enables(ctx);
    let mut l = SpeciesRowLayout { np: en.n_pollut, age_row: -1, temp_row: -1, msx_first: -1, ns: 0 };
    l.ns = l.np;
    if en.age {
        l.age_row = l.ns;
        l.ns += 1;
    }
    if en.temperature {
        l.temp_row = l.ns;
        l.ns += 1;
    }
    if legacy_reactions_active(ctx) {
        l.msx_first = l.ns;
        l.ns += en.n_msx;
    }
    l
}

/// transport::reactSpeciesBlock: one element's reactions species through the integrator.
#[allow(clippy::too_many_arguments)]
fn react_species_block(ctx: &mut SimulationContext, tank: bool, dt: f64, species_block: &mut [f64], pollut: Option<&[f64]>, hrt_seconds: f64, temp_c: f64, element: i32) -> Result<(), StepError> {
    let ns = ctx.reactions.n_species();
    if ns == 0 {
        return Ok(());
    }
    let uns = ns as usize;
    let mut ws = Ws::new(&ctx.reactions);
    let mut block = vec![0.0f64; uns];
    let mut coef_scratch: Vec<f64> = Vec::new();
    let mut hydvar = [0.0f64; RxHydVar::Count as usize];
    hydvar[RxHydVar::Hrt as usize] = hrt_seconds;
    hydvar[RxHydVar::Temp as usize] = if temp_c.is_finite() { temp_c } else { ctx.reactions.default_temp_c };
    if ctx.reactions.wall_geometry && element >= 0 {
        let v = if tank { ctx.nodes.volume[element as usize] } else { ctx.links.volume[element as usize] };
        hydvar[RxHydVar::Av as usize] = reaction_area_volume(ctx, !tank, element, v);
    }
    block[..uns].copy_from_slice(&species_block[..uns]);
    let rep = {
        let rx = &ctx.reactions;
        let coefs: Vec<f64> = rx.coefficients(!tank, element, &mut coef_scratch).to_vec();
        let (forms, spans) = if tank { (&rx.tank_form, &rx.tank_expr) } else { (&rx.pipe_form, &rx.pipe_expr) };
        integrator_step(rx, forms, spans, dt, &mut block, &mut hydvar, &mut ws, pollut, None, &coefs)?
    };
    if rep.ok {
        if clamp_rate_species(&ctx.reactions, tank, &mut block, None) > 0.0 {
            note_negative_clamp(ctx, b"a LARD element");
        }
        species_block[..uns].copy_from_slice(&block[..uns]);
    } else if !ctx.reactions.warned_react_failure {
        ctx.reactions.warned_react_failure = true;
        let mut w = b"Reaction step failed on a LARD element (element state left unchanged): ".to_vec();
        w.extend_from_slice(&rep.error);
        ctx.warnings.push(w);
    }
    Ok(())
}

/// lard::LagrangianSolver
#[derive(Debug, Clone, Default)]
pub struct LagrangianSolver {
    pub store_: SegmentStore,
    pub rwpt_: RwptDispersion,
    pub release_vol_: Vec<f64>,
    pub substep_counter_: u64,
    pub topo_: Vec<i32>,
    pub node_out_links_: Vec<Vec<i32>>,
    pub flow_sign_: Vec<i8>,
    pub node_mass_in_: Vec<f64>,
    pub node_vol_in_: Vec<f64>,
    pub scratch_: Vec<f64>,
    pub treat_cin_: Vec<f64>,
    pub msx_block_: Vec<f64>,
    pub msx_poll_: Vec<f64>,
    pub initialized_: bool,
}

fn assign<T: Clone>(v: &mut Vec<T>, n: usize, x: T) {
    v.clear();
    v.resize(n, x);
}

impl LagrangianSolver {
    pub fn initialized(&self) -> bool {
        self.initialized_
    }

    /// step: one routing step (QUALITY_STEP substeps, then the heat fluxes and the bed, then
    /// publish).
    pub fn step(&mut self, ctx: &mut SimulationContext, dt_routing: f64) -> Result<(), StepError> {
        let l = row_layout(ctx);
        let ns = l.ns;
        if ns <= 0 {
            return Ok(());
        }
        if !self.initialized_ {
            self.init(ctx);
        }
        let nl = ctx.n_links();
        // 0. Flow reversal, once per routing step
        let mut topo_dirty = false;
        for li in 0..nl {
            let ul = li as usize;
            let q = ctx.links.flow[ul];
            if q.abs() <= K_TINY_FLOW {
                continue;
            }
            let sign: i8 = if q >= 0.0 { 1 } else { -1 };
            if sign != self.flow_sign_[ul] {
                if ctx.links.type_[ul] == LinkType::Conduit {
                    self.store_.reverse(li);
                }
                self.flow_sign_[ul] = sign;
                topo_dirty = true;
            }
        }
        if topo_dirty || self.topo_.is_empty() {
            self.compute_topo_order(ctx);
        }
        // X3a: QUALITY_STEP substepping
        let dtq = ctx.options.quality_step;
        let nsub = if dtq > 0.0 && dtq < dt_routing { cvt_i32((dt_routing / dtq).ceil()) } else { 1 };
        let dt = dt_routing / nsub as f64;
        let frac = 1.0 / nsub as f64;
        for _ in 0..nsub {
            self.substep(ctx, dt, frac)?;
        }
        self.apply_fluxes_and_bed(ctx, dt_routing)?;
        self.publish(ctx, dt_routing);
        Ok(())
    }

    /// substep: age, drain, mix, release, disperse, decay, react.
    pub fn substep(&mut self, ctx: &mut SimulationContext, dt: f64, frac: f64) -> Result<(), StepError> {
        let lay = row_layout(ctx);
        let np = lay.np;
        let ns = lay.ns;
        let age = lay.age_row >= 0;
        let heat = lay.temp_row >= 0;
        let nn = ctx.n_nodes();
        let nl = ctx.n_links();

        // AGE (before transport)
        if age {
            self.store_.add_species(lay.age_row, dt);
            for n in 0..nn as usize {
                ctx.water_age_state.node_age[n] += dt;
            }
        }

        // 1. DRAIN (all conduits, before any node mixes)
        assign(&mut self.scratch_, ns as usize, 0.0);
        assign(&mut self.treat_cin_, (if np > 0 { np } else { 0 }) as usize, 0.0);
        for li in 0..nl {
            let ul = li as usize;
            if ctx.links.type_[ul] != LinkType::Conduit {
                continue;
            }
            let q = ctx.links.flow[ul].abs();
            let dn = self.downstream_node(ctx, li);
            let up = self.upstream_node(ctx, li);
            if q > K_TINY_FLOW && dn >= 0 {
                self.scratch_.iter_mut().for_each(|x| *x = 0.0);
                let drained = self.store_.drain_back(li, q * dt, &mut self.scratch_);
                self.add_to_ledger(dn, drained, ns);
            }
            let v_new = ctx.links.volume[ul];
            let v_rem = self.store_.total_volume(li);
            if v_rem > v_new && up >= 0 {
                self.scratch_.iter_mut().for_each(|x| *x = 0.0);
                let shed = self.store_.drain_front(li, v_rem - v_new, &mut self.scratch_);
                self.add_to_ledger(up, shed, ns);
            }
        }

        // 2. MIX in topo order, passthrough zero-volume links
        for t in 0..self.topo_.len() {
            let n = self.topo_[t];
            let un = n as usize;
            let v_old = store_volume_at(&ctx.nodes, n, ctx.nodes.old_depth[un], ctx.nodes.old_volume[un]);
            let v_in = self.node_vol_in_[un];
            let zero_bf = ctx.options.outfall_backflow_zero && ctx.nodes.type_[un] == NodeType::Outfall && v_in + ctx.nodes.qual_vol_in[un] * frac <= 0.0;
            for s in 0..ns {
                let is_age = s == lay.age_row;
                let is_temp = s == lay.temp_row;
                let is_msx = lay.msx_first >= 0 && s >= lay.msx_first;
                let xi = if is_msx { un * (ns - lay.msx_first) as usize + (s - lay.msx_first) as usize } else { 0 };
                let li = un * ns as usize + s as usize;
                let ci = un * np as usize + s as usize;
                let st_old = if is_age {
                    ctx.water_age_state.node_age[un]
                } else if is_temp {
                    ctx.heat_state.node_temp[un]
                } else if is_msx {
                    ctx.reactions.msx_node_conc[xi]
                } else {
                    ctx.nodes.conc[ci]
                };
                let m_ext = if is_age {
                    ctx.water_age_state.node_age_vol_in[un] * dt
                } else if is_temp {
                    ctx.heat_state.node_temp_vol_in[un] * dt
                } else if is_msx {
                    if xi < ctx.reactions.msx_ext_mass_in.len() {
                        ctx.reactions.msx_ext_mass_in[xi] * dt
                    } else {
                        0.0
                    }
                } else {
                    ctx.nodes.qual_mass_in[ci] * dt
                };
                let mut m = st_old * v_old + self.node_mass_in_[li] + m_ext;
                if !is_age && !is_temp && !is_msx && (s as usize) < self.treat_cin_.len() {
                    self.treat_cin_[s as usize] = self.node_mass_in_[li] + m_ext;
                }
                if m < 0.0 && !is_temp {
                    if is_age {
                        book_negative_age_clamp(ctx, n);
                    } else if !is_msx {
                        book_negative_source_clamp(ctx, n, s, -m);
                    }
                    m = 0.0;
                }
                let denom = v_old + v_in + ctx.nodes.qual_vol_in[un] * frac;
                if is_msx && !zero_bf {
                    m += reaction_booster_mass(ctx, n, s - lay.msx_first, denom, v_in + ctx.nodes.qual_vol_in[un] * frac)?;
                }
                let st_new = if zero_bf {
                    0.0
                } else if denom > 1.0e-12 {
                    m / denom
                } else {
                    st_old
                };
                if is_age {
                    ctx.water_age_state.node_age[un] = st_new;
                } else if is_temp {
                    ctx.heat_state.node_temp[un] = st_new;
                } else if is_msx {
                    ctx.reactions.msx_node_conc[xi] = st_new;
                } else {
                    ctx.nodes.conc[ci] = st_new;
                }
                self.node_mass_in_[li] = 0.0;
            }
            self.node_vol_in_[un] = 0.0;
            // P2.3: [TREATMENT] on the freshly mixed node state
            if np > 0 && dt > 0.0 && un < ctx.treatment.has_treatment.len() && ctx.treatment.has_treatment[un] {
                let v_in_total = v_in + ctx.nodes.qual_vol_in[un] * frac;
                for p in 0..np as usize {
                    self.treat_cin_[p] = if v_in_total > 1.0e-12 { cmax(0.0, self.treat_cin_[p]) / v_in_total } else { 0.0 };
                }
                apply_node_treatment(ctx, n, dt, v_in_total / dt, Some(&self.treat_cin_));
            }
            // Zero-volume passthrough of the outgoing non-conduit links
            for k in 0..self.node_out_links_[un].len() {
                let lk = self.node_out_links_[un][k];
                let ul = lk as usize;
                if ctx.links.type_[ul] == LinkType::Conduit {
                    continue;
                }
                let q = ctx.links.flow[ul].abs();
                if q <= K_TINY_FLOW {
                    continue;
                }
                let dn = self.downstream_node(ctx, lk);
                if dn < 0 {
                    continue;
                }
                for p in 0..np as usize {
                    let ni = un * np as usize + p;
                    let lp = ul * np as usize + p;
                    self.scratch_[p] = ctx.nodes.conc[ni];
                    ctx.links.conc[lp] = ctx.nodes.conc[ni];
                }
                if age {
                    self.scratch_[lay.age_row as usize] = ctx.water_age_state.node_age[un];
                    ctx.water_age_state.link_age[ul] = ctx.water_age_state.node_age[un];
                }
                if heat {
                    self.scratch_[lay.temp_row as usize] = ctx.heat_state.node_temp[un];
                    ctx.heat_state.link_temp[ul] = ctx.heat_state.node_temp[un];
                }
                if lay.msx_first >= 0 {
                    let nsp = (ns - lay.msx_first) as usize;
                    for sp in 0..nsp {
                        let c = ctx.reactions.msx_node_conc[un * nsp + sp];
                        self.scratch_[lay.msx_first as usize + sp] = c;
                        ctx.reactions.msx_link_conc[ul * nsp + sp] = c;
                    }
                }
                self.add_to_ledger_rate(dn, q * dt, ns);
            }
        }
        // 3. RELEASE new front segments
        for li in 0..nl {
            let ul = li as usize;
            if ctx.links.type_[ul] != LinkType::Conduit {
                continue;
            }
            let up = self.upstream_node(ctx, li);
            let need = ctx.links.volume[ul] - self.store_.total_volume(li);
            self.release_vol_[ul] = if need > 0.0 { need } else { 0.0 };
            if need <= 0.0 || up < 0 {
                continue;
            }
            let uu = up as usize;
            for p in 0..np as usize {
                self.scratch_[p] = ctx.nodes.conc[uu * np as usize + p];
            }
            if age {
                self.scratch_[lay.age_row as usize] = ctx.water_age_state.node_age[uu];
            }
            if heat {
                self.scratch_[lay.temp_row as usize] = ctx.heat_state.node_temp[uu];
            }
            if lay.msx_first >= 0 {
                let nsp = (ns - lay.msx_first) as usize;
                for sp in 0..nsp {
                    self.scratch_[lay.msx_first as usize + sp] = ctx.reactions.msx_node_conc[uu * nsp + sp];
                }
            }
            self.store_.push_front(li, need, &self.scratch_);
        }

        // 3b. RWPT dispersion
        if ctx.options.lard_rwpt {
            self.substep_counter_ = self.substep_counter_.wrapping_add(1);
            for li in 0..nl {
                let ul = li as usize;
                if ctx.links.type_[ul] != LinkType::Conduit {
                    continue;
                }
                let q = ctx.links.flow[ul].abs();
                if q <= K_TINY_FLOW {
                    continue;
                }
                let row = ctx.link_subtypes.conduit_row(li);
                if row < 0 {
                    continue;
                }
                let ur = row as usize;
                let len = ctx.link_subtypes.conduits.length[ur];
                let vol = ctx.links.volume[ul];
                if len <= 0.0 || vol <= 0.0 {
                    continue;
                }
                let a_flow = vol / len;
                let ubar = q / a_flow;
                let h = ctx.links.depth[ul];
                let circ = ctx.links.xsect_shape[ul] == XsectShape::Circular;
                let rh = rwpt_hyd_radius(a_flow, h, ctx.links.xsect_geom1[ul], circ);
                let rough = ctx.link_subtypes.conduits.roughness[ur];
                self.rwpt_.disperse(&mut self.store_, li, ubar, h, rh, rough, self.release_vol_[ul], q * dt, dt, self.substep_counter_, ctx.options.rwpt_seed as u64, &mut self.scratch_);
            }
        }

        // 4. DECAY (exact exponential)
        for p in 0..np {
            let k = ctx.pollutants.k_decay[p as usize];
            if k == 0.0 {
                continue;
            }
            let f = (-k * dt).exp(); // LIBM
            let mut removed = self.store_.decay_species(p, f);
            for n in 0..nn {
                let idx = n as usize * np as usize + p as usize;
                let v = store_volume(&ctx.nodes, n);
                removed += ctx.nodes.conc[idx] * (1.0 - f) * v;
                ctx.nodes.conc[idx] *= f;
            }
            if (p as usize) < ctx.mass_balance.qual_routing_reacted.len() {
                ctx.mass_balance.qual_routing_reacted[p as usize] += removed;
            }
        }

        // 4b. REACT the reactions species per segment and per node store
        if lay.msx_first >= 0 {
            update_untreated_storage_hrt(ctx, dt);
            let nsp = ns - lay.msx_first;
            if (self.msx_block_.len() as i32) < nsp {
                assign(&mut self.msx_block_, nsp as usize, 0.0);
            }
            if (self.msx_poll_.len() as i32) < np && np > 0 {
                assign(&mut self.msx_poll_, np as usize, 0.0);
            }
            let mut block = std::mem::take(&mut self.msx_block_);
            let mut poll = std::mem::take(&mut self.msx_poll_);
            for li in 0..nl {
                let ul = li as usize;
                if ctx.links.type_[ul] != LinkType::Conduit {
                    continue;
                }
                let cnt = self.store_.count(li);
                for i in 0..cnt {
                    for sp in 0..nsp {
                        block[sp as usize] = self.store_.seg_conc(li, i, lay.msx_first + sp);
                    }
                    for p in 0..np {
                        poll[p as usize] = self.store_.seg_conc(li, i, p);
                    }
                    let seg_temp_c = if lay.temp_row >= 0 { self.store_.seg_conc(li, i, lay.temp_row) } else { f64::NAN };
                    let r = react_species_block(ctx, false, dt, &mut block, if np > 0 { Some(&poll) } else { None }, 0.0, seg_temp_c, li);
                    if let Err(e) = r {
                        self.msx_block_ = block;
                        self.msx_poll_ = poll;
                        return Err(e);
                    }
                    for sp in 0..nsp {
                        self.store_.set_seg_conc(li, i, lay.msx_first + sp, block[sp as usize]);
                    }
                }
            }
            for n in 0..nn {
                let un = n as usize;
                if !node_is_reactor(ctx, n) || ctx.nodes.volume[un] <= ZERO_VOLUME {
                    continue;
                }
                let unsp = nsp as usize;
                for sp in 0..unsp {
                    block[sp] = ctx.reactions.msx_node_conc[un * unsp + sp];
                }
                for p in 0..np as usize {
                    poll[p] = ctx.nodes.conc[un * np as usize + p];
                }
                let hrt = if un < ctx.nodes.hrt.len() { ctx.nodes.hrt[un] } else { 0.0 };
                let node_temp_c = if lay.temp_row >= 0 && un < ctx.heat_state.node_temp.len() { ctx.heat_state.node_temp[un] } else { f64::NAN };
                let r = react_species_block(ctx, true, dt, &mut block, if np > 0 { Some(&poll) } else { None }, hrt, node_temp_c, n);
                if let Err(e) = r {
                    self.msx_block_ = block;
                    self.msx_poll_ = poll;
                    return Err(e);
                }
                for sp in 0..unsp {
                    ctx.reactions.msx_node_conc[un * unsp + sp] = block[sp];
                }
            }
            self.msx_block_ = block;
            self.msx_poll_ = poll;
        }
        Ok(())
    }

    /// publish: the conduit means into links.conc (and age, temperature, species).
    pub fn publish(&mut self, ctx: &mut SimulationContext, dt_routing: f64) {
        let lay = row_layout(ctx);
        let np = lay.np;
        let age = lay.age_row >= 0;
        let heat = lay.temp_row >= 0;
        let nl = ctx.n_links();
        for li in 0..nl {
            let ul = li as usize;
            if ctx.links.type_[ul] != LinkType::Conduit {
                continue;
            }
            self.store_.mean_conc(li, &mut self.scratch_);
            for p in 0..np as usize {
                ctx.links.conc[ul * np as usize + p] = self.scratch_[p];
            }
            if age {
                if self.store_.count(li) > 0 {
                    ctx.water_age_state.link_age[ul] = self.scratch_[lay.age_row as usize];
                } else {
                    ctx.water_age_state.link_age[ul] += dt_routing;
                }
            }
            if heat && self.store_.count(li) > 0 {
                ctx.heat_state.link_temp[ul] = self.scratch_[lay.temp_row as usize];
            }
            if lay.msx_first >= 0 && self.store_.count(li) > 0 {
                let nsp = (lay.ns - lay.msx_first) as usize;
                for sp in 0..nsp {
                    ctx.reactions.msx_link_conc[ul * nsp + sp] = self.scratch_[lay.msx_first as usize + sp];
                }
            }
        }
        ctx.links.conc_old.clone_from(&ctx.links.conc);
        ctx.nodes.conc_old.clone_from(&ctx.nodes.conc);
    }

    /// applyFluxesAndBed: its guard (the heat flux modules are not translated yet).
    fn apply_fluxes_and_bed(&mut self, ctx: &mut SimulationContext, dt: f64) -> Result<(), StepError> {
        let lay = row_layout(ctx);
        let heat_on = lay.temp_row >= 0;
        let bed_on = bed_exchange_enabled(ctx);
        let flux_on = heat_on && (ctx.heat_config.surface_exchange || ctx.heat_config.radiative_exchange);
        if (!flux_on && !bed_on) || !(dt > 0.0) {
            return Ok(());
        }
        Err(StepError::NotTranslated("LARD heat fluxes and bed"))
    }

    /// init: size the store and ledgers, seed age and temperature, one segment per wet conduit.
    pub fn init(&mut self, ctx: &mut SimulationContext) {
        let lay = row_layout(ctx);
        let np = lay.np;
        let ns = lay.ns;
        let age = lay.age_row >= 0;
        let heat = lay.temp_row >= 0;
        let nn = ctx.n_nodes();
        let nl = ctx.n_links();
        self.store_.resize(nl, ns, std::cmp::max(2, ctx.options.max_segments_per_link));
        assign(&mut self.flow_sign_, nl as usize, 1);
        assign(&mut self.node_mass_in_, nn as usize * ns as usize, 0.0);
        assign(&mut self.node_vol_in_, nn as usize, 0.0);
        assign(&mut self.scratch_, ns as usize, 0.0);
        assign(&mut self.treat_cin_, (if np > 0 { np } else { 0 }) as usize, 0.0);
        assign(&mut self.node_out_links_, nn as usize, Vec::new());
        assign(&mut self.release_vol_, nl as usize, 0.0);
        if ctx.options.lard_rwpt {
            self.rwpt_.resize(nl, K_RWPT_PARTICLES_PER_LINK);
        }
        if lay.msx_first >= 0 {
            ensure_msx_state(ctx);
        }
        if age {
            let nsc = ctx.n_subcatches();
            let ws = &mut ctx.water_age_state;
            if ws.node_age.len() != nn as usize {
                ws.resize(nn, nl, nsc);
            }
            if !ws.hotstart_loaded {
                let a0 = ctx.water_age_config.global_age[WaterAgeSource::InitialState as usize];
                if a0 > 0.0 {
                    ws.node_age.iter_mut().for_each(|x| *x = a0);
                    ws.link_age.iter_mut().for_each(|x| *x = a0);
                }
            }
            apply_initial_age_overrides(ctx);
        }
        if heat {
            let t0 = ctx.heat_config.global_temp[HeatSource::InitialState as usize];
            let hstate = &mut ctx.heat_state;
            if hstate.node_temp.len() != nn as usize {
                hstate.resize(nn, nl, t0);
            }
            if !hstate.legacy_seeded {
                hstate.node_temp.iter_mut().for_each(|x| *x = t0);
                hstate.link_temp.iter_mut().for_each(|x| *x = t0);
                apply_initial_temp_overrides(ctx);
                ctx.heat_state.legacy_seeded = true;
            }
        }
        for li in 0..nl {
            let ul = li as usize;
            self.store_.clear_link(li);
            if ctx.links.type_[ul] != LinkType::Conduit {
                continue;
            }
            let v = ctx.links.volume[ul];
            if v <= 0.0 {
                continue;
            }
            for p in 0..np as usize {
                self.scratch_[p] = ctx.links.conc[ul * np as usize + p];
            }
            if age {
                self.scratch_[lay.age_row as usize] = ctx.water_age_state.link_age[ul];
            }
            if heat {
                self.scratch_[lay.temp_row as usize] = ctx.heat_state.link_temp[ul];
            }
            if lay.msx_first >= 0 {
                let nsp = (ns - lay.msx_first) as usize;
                for sp in 0..nsp {
                    self.scratch_[lay.msx_first as usize + sp] = ctx.reactions.msx_link_conc[ul * nsp + sp];
                }
            }
            self.store_.push_front(li, v, &self.scratch_);
            self.flow_sign_[ul] = if ctx.links.flow[ul] >= 0.0 { 1 } else { -1 };
        }
        self.compute_topo_order(ctx);
        self.initialized_ = true;
    }

    fn upstream_node(&self, ctx: &SimulationContext, l: i32) -> i32 {
        let ul = l as usize;
        if self.flow_sign_[ul] >= 0 {
            ctx.links.node1[ul]
        } else {
            ctx.links.node2[ul]
        }
    }
    fn downstream_node(&self, ctx: &SimulationContext, l: i32) -> i32 {
        let ul = l as usize;
        if self.flow_sign_[ul] >= 0 {
            ctx.links.node2[ul]
        } else {
            ctx.links.node1[ul]
        }
    }

    /// addToLedger(n, vol, scratch_, stride): drained mass into node n's ledger. Upstream's
    /// fourth parameter, named np, is the row stride; every call passes ns.
    fn add_to_ledger(&mut self, n: i32, vol: f64, np: i32) {
        let un = n as usize;
        self.node_vol_in_[un] += vol;
        for p in 0..np as usize {
            self.node_mass_in_[un * np as usize + p] += self.scratch_[p];
        }
    }
    /// addToLedgerRate(n, vol, scratch_ concentrations, stride)
    fn add_to_ledger_rate(&mut self, n: i32, vol: f64, np: i32) {
        let un = n as usize;
        self.node_vol_in_[un] += vol;
        for p in 0..np as usize {
            self.node_mass_in_[un * np as usize + p] += vol * self.scratch_[p];
        }
    }

    /// computeTopoOrder: Kahn's order over the flowing links; cycle leftovers appended.
    fn compute_topo_order(&mut self, ctx: &SimulationContext) {
        let nn = ctx.n_nodes();
        let nl = ctx.n_links();
        let mut indeg = vec![0i32; nn as usize];
        for v in self.node_out_links_.iter_mut() {
            v.clear();
        }
        for li in 0..nl {
            let ul = li as usize;
            let up = self.upstream_node(ctx, li);
            let dn = self.downstream_node(ctx, li);
            if up < 0 || dn < 0 {
                continue;
            }
            self.node_out_links_[up as usize].push(li);
            if ctx.links.flow[ul].abs() > K_TINY_FLOW {
                indeg[dn as usize] += 1;
            }
        }
        self.topo_.clear();
        self.topo_.reserve(nn as usize);
        let mut q: Vec<i32> = Vec::new();
        for n in 0..nn {
            if indeg[n as usize] == 0 {
                q.push(n);
            }
        }
        let mut seen = vec![0u8; nn as usize];
        let mut qi = 0usize;
        while qi < q.len() {
            let n = q[qi];
            qi += 1;
            if seen[n as usize] != 0 {
                continue;
            }
            seen[n as usize] = 1;
            self.topo_.push(n);
            for k in 0..self.node_out_links_[n as usize].len() {
                let li = self.node_out_links_[n as usize][k];
                if ctx.links.flow[li as usize].abs() <= K_TINY_FLOW {
                    continue;
                }
                let dn = self.downstream_node(ctx, li);
                if dn >= 0 {
                    indeg[dn as usize] -= 1;
                    if indeg[dn as usize] == 0 {
                        q.push(dn);
                    }
                }
            }
        }
        for n in 0..nn {
            if seen[n as usize] == 0 {
                self.topo_.push(n);
            }
        }
    }
}
