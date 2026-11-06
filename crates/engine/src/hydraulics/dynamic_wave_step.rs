// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/DynamicWave.cpp and DynamicWave.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The dynamic-wave solver's routing step (the run-time half of DWSolver; dynamic_wave.rs holds
// the state and init): execute, initNodeStates, findBypassedLinks, computeLinkGeometry,
// recomputeConduitLossOne, momentumKernels, processDryLink, applyFlowLimits,
// processManningLink, processForceMainLink, updateNodeFlows, gatherConduitNodeFlows,
// limitConduitLosses, computeAASkipFlags, updateNodeDepthsTeam, setNodeDepth,
// commitNodeDepthState, getRoutingStep, getLinkStep, accumulateLinkStepStats,
// getCrownCutoff, getSlotWidth and the file-local conveyArea and computeYnorm.
//
// Differences in form, never in result:
//   * The persistent OpenMP team runs as a team of one: every parallel loop writes disjoint
//     single-producer elements and the only cross-thread combine is an integer count, so
//     upstream's result is the same at any thread count (its own BIT-EXACTNESS note), and the
//     serial order here is that result. `omp single` / `omp for` blocks run in their textual
//     order.
//   * `groups_` (upstream a pointer into the owning Router) is passed in as `groups`.
//   * The non-conduit callback (a std::function capturing the DWSolver upstream) receives the
//     solver's per-node arrays it uses through NcView (nodeNewSurfAreaDataMut, nodeSumDqdh,
//     isBypassed, isInitialized).
//   * Environment switches read with getenv upstream take their unset defaults:
//     SWMM_DISABLE_BYPASS_MASK (mask on), SWMM_VJ_M3 / SWMM_VJ_M4 (on); the trace hooks
//     (SWMM_TRACE_LINK, SWMM_TRACE_NODE, OPENSWMM_DT_TRACE) print only and are left out.
//   * SWMM_FAST_MANNING_POW is not defined in the reference build: std::pow is used.
//
// NOT YET TRANSLATED (execute returns StepError::NotTranslated before changing anything):
// DYNAMIC_SLOT (applyDPSGeometry, updateDPSState, spatialSmoothP), TPA (updateTpaLatch,
// applyTpaGeometry), physical junction areas (NodeData::storage_area > 0: chamberEndArea,
// chamberPipeVolumeChange and the chamber branches), and node LID storage (lidnode ports).

use crate::core::cmath::{cmax, cmin};
use crate::core::constants::{self, EXTRAN_CROWN_CUTOFF, FUDGE, GRAVITY, MAX_VELOCITY, MIN_TIMESTEP, SLOT_CROWN_CUTOFF};
use crate::core::context::SimulationContext;
use crate::core::engine_step_routing::StepError;
use crate::core::options::NodeContinuity;
use crate::core::unit_conversion as ucf;
use crate::data::forcing::ForcingMode;
use crate::data::link::{FlowClass, LinkType, XsectShape};
use crate::data::node::NodeType;
use crate::hydraulics::culvert;
use crate::hydraulics::dynamic_wave::{build_xsp, DWSolver, MomentumCategory, NodeTile, SurchargeMethod};
use crate::hydraulics::force_main;
use crate::hydraulics::hyd_structures::StructureSolver;
use crate::hydraulics::node;
use crate::hydraulics::outfall;
use crate::hydraulics::xsect_batch::{XSectGroups, XSectParams};
use crate::hydraulics::xsection as xsect;
use crate::hydrology::lid_node_step as lidnode;

/// The solver state the non-conduit callback reaches through its captured DWSolver upstream.
pub struct NcView<'a> {
    /// nodeNewSurfAreaDataMut()
    pub new_surf_area: &'a mut Vec<f64>,
    /// nodeSumDqdh(n)
    pub sumdqdh: &'a mut Vec<f64>,
    /// isBypassed(j)
    pub bypassed: &'a [u8],
    /// isInitialized()
    pub initialized: bool,
}

impl NcView<'_> {
    pub fn is_bypassed(&self, j: i32) -> bool {
        let uj = j as usize;
        uj < self.bypassed.len() && self.bypassed[uj] != 0
    }
}

/// DWSolver::NonConduitFlowFunc: (ctx, dt, Picard iteration, the solver's node arrays, the
/// engine's StructureSolver).
pub type NonConduitFn<'a> = dyn FnMut(&mut SimulationContext, f64, i32, &mut NcView, &mut StructureSolver) -> Result<(), StepError> + 'a;

/// SurfaceExchange.hpp boundSurfaceOutfallFlow.
pub fn bound_surface_outfall_flow(ctx: &SimulationContext, link: i32, q: f64) -> f64 {
    if ctx.surface_outfall_link_limit.is_empty() || q == 0.0 {
        return q;
    }
    let donor = if q > 0.0 { ctx.links.node1[link as usize] } else { ctx.links.node2[link as usize] };
    if donor < 0 || donor as usize >= ctx.surface_outfall_link_limit.len() {
        return q;
    }
    let cap = ctx.surface_outfall_link_limit[donor as usize];
    if q > 0.0 {
        cmin(q, cap)
    } else {
        cmax(q, -cap)
    }
}

/// AquiferLinkExchange.hpp aquiferSeepFactor.
pub fn aquifer_seep_factor(depth: f64, gw_head_rel: f64, d_c: f64) -> f64 {
    if !(d_c > 0.0) {
        return 1.0;
    }
    let h_eff = if gw_head_rel > 0.0 { gw_head_rel } else { 0.0 };
    cmin((depth - h_eff) / d_c, 1.0)
}

/// AquiferLinkExchange.hpp capAquiferGain.
pub fn cap_aquifer_gain(seep_rate: f64, gain_max: f64) -> f64 {
    if seep_rate < -gain_max {
        -gain_max
    } else {
        seep_rate
    }
}

/// AquiferLinkExchange.hpp applySeepageForcing.
pub fn apply_seepage_forcing(ctx: &SimulationContext, link: usize, seep_rate: &mut f64) {
    let f = &ctx.forcing;
    if link >= f.link_seepage_mode.len() {
        return;
    }
    match f.link_seepage_mode[link] {
        ForcingMode::Override => *seep_rate = f.link_seepage_value[link],
        ForcingMode::Add => *seep_rate += f.link_seepage_value[link],
        _ => {}
    }
}

/// isOpenShape (legacy xsect_isOpen, storage shape codes).
fn is_open_shape(shape: XsectShape) -> bool {
    matches!(
        shape,
        XsectShape::Dummy
            | XsectShape::RectOpen
            | XsectShape::Trapezoidal
            | XsectShape::Triangular
            | XsectShape::Parabolic
            | XsectShape::Power
            | XsectShape::Irregular
            | XsectShape::StreetXsect
    )
}

/// conveyArea: the momentum area without the Preissmann slot (min(a, a_full)).
#[inline]
fn convey_area(a: f64, a_full: f64) -> f64 {
    if a > a_full {
        a_full
    } else {
        a
    }
}

/// computeYnorm: Manning section factor inverted to a normal depth.
fn compute_ynorm(xs: &XSectParams, beta: f64, q_max: f64, mut q: f64) -> f64 {
    if beta <= 0.0 || q <= 0.0 {
        return 0.0;
    }
    if q_max > 0.0 && q > q_max {
        q = q_max;
    }
    let s = q / beta;
    let a = xsect::get_a_of_s(xs, s);
    xsect::get_y_of_a(xs, a)
}

impl DWSolver {
    pub fn get_crown_cutoff(&self) -> f64 {
        if self.surcharge_method == SurchargeMethod::Slot || self.surcharge_method == SurchargeMethod::DynamicSlot || self.surcharge_method == SurchargeMethod::Tpa {
            return SLOT_CROWN_CUTOFF;
        }
        EXTRAN_CROWN_CUTOFF
    }

    /// Op-for-op legacy dwflow.c getSlotWidth (0 unless SLOT).
    pub fn get_slot_width(&self, y: f64, y_full: f64, w_max: f64, shape: XsectShape) -> f64 {
        if self.surcharge_method != SurchargeMethod::Slot {
            return 0.0;
        }
        if y_full <= 0.0 || is_open_shape(shape) {
            return 0.0;
        }
        let y_norm = y / y_full;
        if y_norm < SLOT_CROWN_CUTOFF {
            return 0.0;
        }
        if y_norm > 1.78 {
            return 0.01 * w_max;
        }
        w_max * 0.5423 * (-y_norm.powf(2.4)).exp() // LIBM pow, exp
    }

    /// The paths execute does not translate yet (see the file notes), checked before any state
    /// changes.
    fn not_translated(&self, ctx: &SimulationContext) -> Option<&'static str> {
        if self.surcharge_method == SurchargeMethod::DynamicSlot {
            return Some("DYNAMIC_SLOT surcharge method");
        }
        if self.surcharge_method == SurchargeMethod::Tpa {
            return Some("TPA surcharge method");
        }
        if ctx.nodes.storage_area.iter().any(|&a| a > 0.0) {
            return Some("physical junction areas");
        }
        if !ctx.coupled_pond_area.is_empty() {
            return Some("2D-coupled ponds");
        }
        None
    }

    /// DWSolver::execute: one routing step's Picard iterations. Returns the iteration count.
    pub fn execute(&mut self, groups: &XSectGroups, ctx: &mut SimulationContext, dt: f64, structures: &mut StructureSolver, mut non_conduit_fn: Option<&mut NonConduitFn>) -> Result<i32, StepError> {
        if let Some(w) = self.not_translated(ctx) {
            return Err(StepError::NotTranslated(w));
        }

        // Per-timestep constant
        self.dt_gravity_ = dt * GRAVITY;

        // Zero-loss latch (flips permanently the first time a nonzero loss becomes possible).
        if self.evap_rate > 0.0 || self.any_conduit_seep_ {
            self.losses_all_zero_ = false;
        }

        // area_mid_ from the PREVIOUS timestep for the unsteady momentum term.
        let n_am = self.area_mid_.len();
        self.area_old_[..n_am].copy_from_slice(&self.area_mid_);

        // Clear bypass flags at the start of each timestep
        self.bypassed_.iter_mut().for_each(|b| *b = 0);

        // Refresh the node-invariant tile for this step.
        {
            let nd = &ctx.nodes;
            if self.node_tile_.len() != self.n_nodes_ as usize {
                self.node_tile_.resize(self.n_nodes_ as usize, NodeTile::default());
            }
            self.unit_sys_ = ucf::get_unit_system(ctx.options.flow_units as i32);
            for i in 0..self.n_nodes_ {
                let ui = i as usize;
                let t = &mut self.node_tile_[ui];
                t.full_depth = nd.full_depth[ui];
                t.y_crown = nd.crown_elev[ui] - nd.invert_elev[ui];
                t.invert_elev = nd.invert_elev[ui];
                t.ponded_area = nd.ponded_area[ui];
                t.sur_depth = nd.sur_depth[ui];
                t.full_volume = nd.full_volume[ui];
                t.rpt_full_volume = nd.rpt_full_volume[ui];
                t.degree = nd.degree[ui];
                t.is_storage = if nd.type_[ui] == NodeType::Storage || nd.storage_area[ui] > 0.0 { 1 } else { 0 };
                t.is_outfall = if nd.type_[ui] == NodeType::Outfall { 1 } else { 0 };
                t.is_virtual = if ui < nd.is_virtual.len() && nd.is_virtual[ui] != 0 { 1 } else { 0 };
            }
        }

        let has_nc = non_conduit_fn.is_some();
        let mut t_steps = 0;
        let mut t_converged = false;
        let mut unconv_shared: i32;

        while t_steps < self.max_trials {
            // Step 0: init node states
            self.init_node_states(ctx);

            // Step 1: cross-section geometry for every conduit
            self.compute_link_geometry(groups, ctx);

            // Steps 1.5 + 2: losses, STEP E overrides, momentum, flow commit
            self.momentum_kernels(ctx, dt, t_steps);

            // (omp single)
            {
                if !ctx.surface_outfall_link_limit.is_empty() {
                    for ci in 0..self.n_conduits_ as usize {
                        let j = self.tile_uj_[ci];
                        let uj = j as usize;
                        let accepted = bound_surface_outfall_flow(ctx, j, self.new_flow_[uj]);
                        if accepted != self.new_flow_[uj] {
                            self.new_flow_[uj] = accepted;
                            ctx.links.flow[uj] = accepted;
                            self.q1_[uj] = accepted / self.tile_barrels_d_[ci];
                            self.dqdh_[uj] = 0.0;
                            ctx.links.dqdh[uj] = 0.0;
                        }
                    }
                }
                lidnode::reset_ports(ctx)?;
                for ci in 0..self.n_conduits_ as usize {
                    let j = self.tile_uj_[ci];
                    let uj = j as usize;
                    if !lidnode::active(ctx, self.tile_n1_[ci]) && !lidnode::active(ctx, self.tile_n2_[ci]) {
                        continue;
                    }
                    let before = self.new_flow_[uj];
                    let after = lidnode::exchange_ports(ctx, j, before, dt)?;
                    self.new_flow_[uj] = after;
                    ctx.links.flow[uj] = after;
                    self.q1_[uj] = after / self.tile_barrels_d_[ci];
                    if before != 0.0 && before != after {
                        self.dqdh_[uj] *= after / before;
                        ctx.links.dqdh[uj] = self.dqdh_[uj];
                    }
                }
            }

            // Step 3a: scatter conduit flows to nodes
            if has_nc {
                self.gather_conduit_node_flows(ctx);
            }

            // Serial scatter/structures phase
            {
                if !has_nc {
                    self.update_node_flows(ctx);
                }
                if let Some(f) = non_conduit_fn.as_mut() {
                    let initialized = self.is_initialized();
                    let mut view = NcView { new_surf_area: &mut self.xnode_.new_surf_area, sumdqdh: &mut self.xnode_.sumdqdh, bypassed: &self.bypassed_, initialized };
                    f(ctx, dt, t_steps, &mut view, structures)?;
                }
                if !self.losses_all_zero_ {
                    self.limit_conduit_losses(ctx, dt);
                }
                outfall::set_all_outfall_depths(ctx, dt);
                if self.anderson_accel {
                    self.compute_aa_skip_flags(ctx);
                }
                unconv_shared = 0;
            }

            // Step 6: node depths and convergence
            self.update_node_depths_team(ctx, dt, t_steps, &mut unconv_shared);
            t_converged = unconv_shared == 0;

            t_steps += 1;

            if t_steps > 1 {
                if t_converged {
                    break;
                }
                self.find_bypassed_links(ctx);
            }
        }

        let steps = t_steps;
        let converged = t_converged;

        // Post-Picard: per-node non-convergence counts (outfalls excluded)
        if !converged {
            for i in 0..self.n_nodes_ {
                let ui = i as usize;
                if self.node_tile_[ui].is_outfall == 0 && self.xnode_.converged[ui] == 0 {
                    ctx.nodes.stat_non_converged_count[ui] += 1;
                }
            }
        }

        // (DYNAMIC_SLOT updateDPSState: refused at entry)
        self.vj_accumulate_residuals(ctx);

        self.last_converged_ = converged;
        Ok(steps)
    }

    /// initNodeStates: reset per-node accumulators and the surface-area baseline.
    fn init_node_states(&mut self, ctx: &mut SimulationContext) {
        let unit_sys = ucf::get_unit_system(ctx.options.flow_units as i32);
        for i in 0..self.n_nodes_ {
            let ui = i as usize;
            self.xnode_.converged[ui] = 0;
            self.xnode_.sumdqdh[ui] = 0.0;

            let node_can_pond = ctx.options.allow_ponding || (ui < ctx.coupled_node.len() && ctx.coupled_node[ui] != 0);
            let t = self.node_tile_[ui];
            if t.is_storage == 0 {
                self.xnode_.new_surf_area[ui] = if node_can_pond && ctx.nodes.depth[ui] > t.full_depth && t.ponded_area > 0.0 { t.ponded_area } else { 0.0 };
                // (coupled_pond_area is empty in a 1D build; execute refuses it otherwise)
            } else {
                self.xnode_.new_surf_area[ui] = if node_can_pond {
                    lidnode_free_ponded_area(ctx, i, ctx.nodes.depth[ui], unit_sys)
                } else {
                    node::get_surf_area(&ctx.nodes, i, ctx.nodes.depth[ui], Some(&ctx.tables), unit_sys, Some(&ctx.node_subtypes))
                };
            }

            // (storage_area > 0: refused in execute)

            let nodes = &mut ctx.nodes;
            nodes.inflow[ui] = 0.0;
            nodes.outflow[ui] = nodes.losses[ui];
            let lat = nodes.lat_flow[ui];
            if lat >= 0.0 {
                nodes.inflow[ui] += lat;
            } else {
                nodes.outflow[ui] -= lat;
            }
        }
    }

    /// findBypassedLinks: links whose both end nodes converged skip the next iteration.
    fn find_bypassed_links(&mut self, ctx: &SimulationContext) {
        let links = &ctx.links;
        for j in 0..self.n_links_ {
            let uj = j as usize;
            let n1 = links.node1[uj];
            let n2 = links.node2[uj];
            if n1 < 0 || n2 < 0 {
                self.bypassed_[uj] = 0;
                continue;
            }
            let b = if self.xnode_.converged[n1 as usize] != 0 && self.xnode_.converged[n2 as usize] != 0 { 1 } else { 0 };
            self.bypassed_[uj] = b;
        }
        if !self.vjunc_.is_empty() {
            for p in self.vjunc_.iter() {
                if p.link_a < 0 || p.link_b < 0 {
                    continue;
                }
                let ua = p.link_a as usize;
                let ub = p.link_b as usize;
                if self.bypassed_[ua] != self.bypassed_[ub] {
                    self.bypassed_[ua] = 0;
                    self.bypassed_[ub] = 0;
                }
            }
        }
    }

    /// computeLinkGeometry: depths, widths, flow classes, surface areas, areas and hydraulic
    /// radii of every non-bypassed conduit.
    fn compute_link_geometry(&mut self, groups: &XSectGroups, ctx: &mut SimulationContext) {
        groups.set_bypass_mask(Some(&self.bypassed_));

        let slot_mode = self.surcharge_method == SurchargeMethod::Slot || self.surcharge_method == SurchargeMethod::DynamicSlot;

        // STEP A + STEP B prep
        for ci in 0..self.n_conduits_ as usize {
            let uci = ci;
            let j = self.tile_uj_[uci];
            let uj = j as usize;
            if self.bypassed_[uj] != 0 {
                continue;
            }
            let un1 = self.tile_n1_[uci];
            let un2 = self.tile_n2_[uci];

            let inv1 = self.tile_inv1_elev_[uci];
            let inv2 = self.tile_inv2_elev_[uci];
            let z1 = inv1 + self.tile_z1_off_[uci];
            let z2 = inv2 + self.tile_z2_off_[uci];
            let h1 = cmax(lidnode::port_depth(ctx, un1, self.tile_z1_off_[uci]) + inv1, z1);
            let h2 = cmax(lidnode::port_depth(ctx, un2, self.tile_z2_off_[uci]) + inv2, z2);

            let mut y1 = cmax(h1 - z1, FUDGE);
            let mut y2 = cmax(h2 - z2, FUDGE);

            let yf = self.tile_y_full_[uci];
            if !slot_mode {
                y1 = cmin(y1, yf);
                y2 = cmin(y2, yf);
            }
            let y_mid = 0.5 * (y1 + y2);

            self.depth1_[uj] = y1;
            self.depth2_[uj] = y2;
            self.depth_mid_[uj] = y_mid;
            self.h1_[uj] = h1;
            self.h2_[uj] = h2;
            self.fasnh_[uj] = 1.0;

            if !slot_mode {
                let bs = self.tile_xsect_batch_shape_[uci];
                let is_open = xsect::is_open(bs);
                let y_cap = if !is_open && yf > 0.0 { EXTRAN_CROWN_CUTOFF * yf } else { 1e30 };
                self.wcap_d1_[uj] = cmin(y1, y_cap);
                self.wcap_d2_[uj] = cmin(y2, y_cap);
                self.wcap_dm_[uj] = cmin(y_mid, y_cap);
            }
        }

        // STEP B: widths
        if !slot_mode {
            groups.compute_widths_triple_team(&self.wcap_d1_, &self.wcap_d2_, &self.wcap_dm_, &mut self.width1_, &mut self.width2_, &mut self.width_mid_, self.n_links_, 0, 1);
        } else {
            groups.compute_widths_triple_team(&self.depth1_, &self.depth2_, &self.depth_mid_, &mut self.width1_, &mut self.width2_, &mut self.width_mid_, self.n_links_, 0, 1);
        }
        if slot_mode {
            for ci in 0..self.n_conduits_ as usize {
                let uj = self.tile_uj_[ci] as usize;
                let yf = self.tile_y_full_[ci];
                let wm = self.tile_w_max_[ci];
                let shape = self.tile_shape_[ci];
                let mut ws = self.get_slot_width(self.depth1_[uj], yf, wm, shape);
                if ws > 0.0 {
                    self.width1_[uj] = ws;
                }
                ws = self.get_slot_width(self.depth2_[uj], yf, wm, shape);
                if ws > 0.0 {
                    self.width2_[uj] = ws;
                }
                ws = self.get_slot_width(self.depth_mid_[uj], yf, wm, shape);
                if ws > 0.0 {
                    self.width_mid_[uj] = ws;
                }
            }
        }

        // STEP C: flow classification + surface area
        let crown_cutoff = self.get_crown_cutoff();
        for ci in 0..self.n_conduits_ as usize {
            let uci = ci;
            let j = self.tile_uj_[uci];
            let uj = j as usize;
            if self.bypassed_[uj] != 0 {
                continue;
            }
            let un1 = self.tile_n1_[uci] as usize;
            let un2 = self.tile_n2_[uci] as usize;

            let yf = self.tile_y_full_[uci];
            let z1_off = self.tile_z1_off_[uci];
            let z2_off = self.tile_z2_off_[uci];
            let beta_j = self.tile_beta_[uci];
            let qmax_j = self.tile_q_max_[uci];
            let inv1 = self.tile_inv1_elev_[uci];
            let inv2 = self.tile_inv2_elev_[uci];
            let mut y1 = self.depth1_[uj];
            let mut y2 = self.depth2_[uj];
            let mut h1 = self.h1_[uj];
            let mut h2 = self.h2_[uj];
            let length = self.tile_length_[uci];

            let mut surf_area1 = 0.0;
            let mut surf_area2 = 0.0;
            let mut fc = FlowClass::Subcritical;
            let mut fasnh = 1.0;

            let mut cached_yn = 0.0;
            let mut cached_yc = 0.0;
            let mut cached_xs: Option<XSectParams> = None;
            macro_rules! xs_ref {
                () => {{
                    if cached_xs.is_none() {
                        cached_xs = Some(build_xsp(ctx, uj));
                    }
                    cached_xs.as_ref().unwrap()
                }};
            }

            let both_full = (y1 >= yf) & (y2 >= yf);
            let both_wet = (y1 > FUDGE) & (y2 > FUDGE);
            let both_dry = (y1 <= FUDGE) & (y2 <= FUDGE);
            let up_dry_only = (y1 <= FUDGE) & (y2 > FUDGE);

            if both_full {
                fc = FlowClass::Subcritical;
            } else if both_dry {
                fc = FlowClass::Dry;
            } else if self.tile_has_offset_[uci] == 0 {
                if up_dry_only {
                    if h2 < inv1 {
                        fc = FlowClass::UpDry;
                    }
                } else if !both_wet {
                    if h1 < inv2 {
                        fc = FlowClass::DnDry;
                    }
                }
            } else {
                let nodes = &ctx.nodes;
                let mut z1_off_eff = z1_off;
                let mut z2_off_eff = z2_off;
                if nodes.type_[un1] == NodeType::Outfall {
                    z1_off_eff = cmax(0.0, z1_off_eff - nodes.depth[un1]);
                }
                if nodes.type_[un2] == NodeType::Outfall {
                    z2_off_eff = cmax(0.0, z2_off_eff - nodes.depth[un2]);
                }
                let q = self.q1_[uj].abs();

                if both_wet {
                    if self.q1_[uj] < 0.0 {
                        if z1_off_eff > 0.0 {
                            let xs = xs_ref!();
                            cached_yn = compute_ynorm(xs, beta_j, qmax_j, q);
                            cached_yc = xsect::get_ycrit(xs, q);
                            let yc_min = cmin(cached_yn, cached_yc);
                            if y1 < yc_min {
                                fc = FlowClass::UpCritical;
                            }
                        }
                    } else if z2_off_eff > 0.0 {
                        let xs = xs_ref!();
                        cached_yn = compute_ynorm(xs, beta_j, qmax_j, q);
                        cached_yc = xsect::get_ycrit(xs, q);
                        let yc_min = cmin(cached_yn, cached_yc);
                        let yc_max = cmax(cached_yn, cached_yc);
                        if y2 < yc_min {
                            fc = FlowClass::DnCritical;
                        } else if y2 < yc_max {
                            if yc_max - yc_min < FUDGE {
                                fasnh = 0.0;
                            } else {
                                fasnh = (yc_max - y2) / (yc_max - yc_min);
                            }
                        }
                    }
                } else if up_dry_only {
                    let z1_elev = inv1 + z1_off;
                    if h2 < z1_elev {
                        fc = FlowClass::UpDry;
                    } else if z1_off_eff > 0.0 {
                        let xs = xs_ref!();
                        cached_yn = compute_ynorm(xs, beta_j, qmax_j, q);
                        cached_yc = xsect::get_ycrit(xs, q);
                        fc = FlowClass::UpCritical;
                    }
                } else {
                    let z2_elev = inv2 + z2_off;
                    if h1 < z2_elev {
                        fc = FlowClass::DnDry;
                    } else if z2_off_eff > 0.0 {
                        let xs = xs_ref!();
                        cached_yn = compute_ynorm(xs, beta_j, qmax_j, q);
                        cached_yc = xsect::get_ycrit(xs, q);
                        fc = FlowClass::DnCritical;
                    }
                }
            }

            match fc {
                FlowClass::Subcritical => {
                    let w1 = self.width1_[uj];
                    let w_m = self.width_mid_[uj];
                    let w2 = self.width2_[uj];
                    surf_area1 = (w1 + w_m) * length / 4.0;
                    surf_area2 = (w_m + w2) * length / 4.0 * fasnh;
                }
                FlowClass::UpCritical => {
                    y1 = cached_yc;
                    if cached_yn < cached_yc {
                        y1 = cached_yn;
                    }
                    y1 = cmax(y1, FUDGE);
                    let z1_elev = inv1 + z1_off;
                    h1 = z1_elev + y1;
                    let y_mid = cmax(0.5 * (y1 + y2), FUDGE);
                    let w2 = self.width2_[uj];
                    let xs = xs_ref!();
                    let w_m = xsect::get_w_of_y(xs, y_mid);
                    let w_slot_m = self.get_slot_width(y_mid, yf, self.tile_w_max_[uci], self.tile_shape_[uci]);
                    let w_msa = if w_slot_m > 0.0 {
                        w_slot_m
                    } else if y_mid / yf >= crown_cutoff && !xsect::is_open(self.tile_xsect_batch_shape_[uci]) {
                        xsect::get_w_of_y(xs, crown_cutoff * yf)
                    } else {
                        w_m
                    };
                    surf_area2 = (w_msa + w2) * length * 0.5;
                    self.depth1_[uj] = y1;
                    self.depth_mid_[uj] = y_mid;
                    self.width_mid_[uj] = w_m;
                    self.h1_[uj] = h1;
                }
                FlowClass::DnCritical => {
                    y2 = cached_yc;
                    if cached_yn < cached_yc {
                        y2 = cached_yn;
                    }
                    y2 = cmax(y2, FUDGE);
                    let z2_elev = inv2 + z2_off;
                    h2 = z2_elev + y2;
                    let y_mid = cmax(0.5 * (y1 + y2), FUDGE);
                    let w1 = self.width1_[uj];
                    let xs = xs_ref!();
                    let w_m = xsect::get_w_of_y(xs, y_mid);
                    let w_slot_m = self.get_slot_width(y_mid, yf, self.tile_w_max_[uci], self.tile_shape_[uci]);
                    let w_msa = if w_slot_m > 0.0 {
                        w_slot_m
                    } else if y_mid / yf >= crown_cutoff && !xsect::is_open(self.tile_xsect_batch_shape_[uci]) {
                        xsect::get_w_of_y(xs, crown_cutoff * yf)
                    } else {
                        w_m
                    };
                    surf_area1 = (w1 + w_msa) * length * 0.5;
                    self.depth2_[uj] = y2;
                    self.depth_mid_[uj] = y_mid;
                    self.width_mid_[uj] = w_m;
                    self.h2_[uj] = h2;
                }
                FlowClass::UpDry => {
                    let y_mid = cmax(0.5 * (FUDGE + y2), FUDGE);
                    let w1 = self.width1_[uj];
                    let w2 = self.width2_[uj];
                    let w_m = self.width_mid_[uj];
                    surf_area2 = (w_m + w2) * length / 4.0;
                    if z1_off <= 0.0 {
                        surf_area1 = (w1 + w_m) * length / 4.0;
                    }
                    self.depth1_[uj] = FUDGE;
                    self.depth_mid_[uj] = y_mid;
                    let xs = xs_ref!();
                    self.width_mid_[uj] = xsect::get_w_of_y(xs, y_mid);
                }
                FlowClass::DnDry => {
                    let y_mid = cmax(0.5 * (y1 + FUDGE), FUDGE);
                    let w1 = self.width1_[uj];
                    let w2 = self.width2_[uj];
                    let w_m = self.width_mid_[uj];
                    surf_area1 = (w_m + w1) * length / 4.0;
                    if z2_off <= 0.0 {
                        surf_area2 = (w2 + w_m) * length / 4.0;
                    }
                    self.depth2_[uj] = FUDGE;
                    self.depth_mid_[uj] = y_mid;
                    let xs = xs_ref!();
                    self.width_mid_[uj] = xsect::get_w_of_y(xs, y_mid);
                }
                FlowClass::Dry => {
                    surf_area1 = FUDGE * length / 2.0;
                    surf_area2 = surf_area1;
                }
                _ => {}
            }
            let _ = (h1, h2);
            let _ = cached_xs;

            self.surf_area1_[uj] = surf_area1;
            self.surf_area2_[uj] = surf_area2;
            self.fasnh_[uj] = fasnh;
            ctx.links.flow_class[uj] = fc;
        }

        // STEP D: areas and hydraulic radii from the (modified) depths
        groups.compute_area_hyd_rad_triple_team(
            &self.depth1_,
            &self.depth2_,
            &self.depth_mid_,
            &mut self.area1_,
            &mut self.area2_,
            &mut self.area_mid_,
            &mut self.hrad1_,
            &mut self.hrad_mid_,
            self.n_links_,
            0,
            1,
        );
        // STEP E (DYNAMIC_SLOT / TPA geometry): refused in execute. The static slot / EXTRAN
        // overrides are in momentum_kernels.
    }

    /// recomputeConduitLossOne: the per-iteration evaporation / seepage rates of one conduit.
    fn recompute_conduit_loss_one(&self, ctx: &mut SimulationContext, dt: f64, ci: usize) {
        let evap = self.evap_rate;
        let uci = ci;
        let uj = self.conduit_idx_[uci];
        let u = uj as usize;
        let ucr = self.tile_crow_[uci] as usize;

        let fc = ctx.links.flow_class[u];
        if fc == FlowClass::Dry || fc == FlowClass::UpDry || fc == FlowClass::DnDry || self.area_mid_[u] <= FUDGE || self.tile_is_closed_[uci] != 0 || ctx.links.setting[u] == 0.0 {
            ctx.link_subtypes.conduits.evap_loss_rate[ucr] = 0.0;
            ctx.link_subtypes.conduits.seep_loss_rate[ucr] = 0.0;
            return;
        }

        let depth = 0.5 * (ctx.links.old_depth[u] + ctx.links.depth[u]);
        let mut evap_loss = 0.0;
        let mut seep_loss = 0.0;

        if depth > FUDGE {
            let cd = &ctx.link_subtypes.conduits;
            let mut length = cd.true_length[ucr];
            if length <= 0.0 {
                length = cd.length[ucr];
            }
            if length <= 0.0 {
                length = cd.mod_length[ucr];
            }
            let shape = ctx.links.xsect_batch_shape[u];

            let gw_two_way = !cd.gw_coupled.is_empty() && cd.gw_coupled[ucr] != 0;
            let k_seep = if gw_two_way { cd.gw_kc[ucr] } else { cd.seep_rate[ucr] };
            let want_evap = xsect::is_open(shape) && evap > 0.0;
            let want_seep = k_seep > 0.0;
            if want_evap || want_seep {
                let xs = build_xsp(ctx, u);
                if want_evap {
                    let top_width = xsect::get_w_of_y(&xs, depth);
                    evap_loss = top_width * length * evap;
                }
                if want_seep {
                    let mut d_seep = depth;
                    if shape != crate::hydraulics::xsect_batch::XSectShape::RectClosed as i32 && d_seep >= xs.yw_max {
                        d_seep = xs.yw_max;
                    }
                    let width = if shape == crate::hydraulics::xsect_batch::XSectShape::RectClosed as i32 { xs.w_max } else { xsect::get_w_of_y(&xs, d_seep) };
                    seep_loss = k_seep * width * length;
                    seep_loss *= ctx.climate_state.infil_factor;
                    if gw_two_way {
                        seep_loss *= aquifer_seep_factor(depth, cd.gw_head_rel[ucr], cd.gw_dc[ucr]);
                        seep_loss = cap_aquifer_gain(seep_loss, cd.gw_gain_max[ucr]);
                    }
                }
            }

            let total = evap_loss + cmax(seep_loss, 0.0);
            if total > 0.0 {
                let q = ctx.links.volume[u] / (dt * cd.barrels[ucr] as f64);
                if total > q {
                    evap_loss = evap_loss * q / total;
                    if seep_loss > 0.0 {
                        seep_loss = seep_loss * q / total;
                    }
                }
            }
        }

        apply_seepage_forcing(ctx, u, &mut seep_loss);

        ctx.link_subtypes.conduits.evap_loss_rate[ucr] = evap_loss;
        ctx.link_subtypes.conduits.seep_loss_rate[ucr] = seep_loss;
    }

    /// vjPrepareIteration: the virtual-junction pair cache for this Picard iteration.
    fn vj_prepare_iteration(&mut self, ctx: &SimulationContext, dt: f64) {
        let links = &ctx.links;
        let full_momentum = ctx.options.virtual_junction_momentum == 1;
        for pi in 0..self.vjunc_.len() {
            let mut p = self.vjunc_[pi];
            p.active = 0;
            p.dirn = 0;
            p.dq4j = 0.0;
            let skip = p.through == 0 || p.up_link < 0 || p.dn_link < 0 || (p.node >= 0 && (p.node as usize) < ctx.nodes.is_inlet.len() && ctx.nodes.is_inlet[p.node as usize] != 0);
            if !skip {
                let uu = p.up_link as usize;
                let ud = p.dn_link as usize;
                let ucu = self.tile_uj_to_ci_[uu] as usize;
                let ucd = self.tile_uj_to_ci_[ud] as usize;
                let aj = self.area2_[uu];
                let wj = self.width2_[uu];
                if aj > FUDGE {
                    let aj_full = links.xsect_a_full[uu];
                    let q_u = links.flow[uu] / self.tile_barrels_d_[ucu];
                    let q_d = links.flow[ud] / self.tile_barrels_d_[ucd];
                    let mut vj = 0.5 * (q_u + q_d) / convey_area(aj, aj_full);
                    if vj.abs() > MAX_VELOCITY {
                        vj = if vj > 0.0 { MAX_VELOCITY } else { -MAX_VELOCITY };
                    }
                    let dh = if wj > FUDGE { aj / wj } else { 0.0 };
                    let frj = if dh > 0.0 { vj.abs() / (constants::GRAVITY * dh).sqrt() } else { 0.0 };
                    p.sigma_j = cmax(0.0, cmin(1.0, 2.0 * (1.0 - frj)));
                    p.a_up_mid = convey_area(self.area_mid_[uu], aj_full);
                    p.r_up_mid = self.hrad_mid_[uu];
                    p.active = 1;
                    if q_u >= 0.0 && q_d >= 0.0 && q_u + q_d > 0.0 {
                        p.dirn = 1;
                    } else if q_u <= 0.0 && q_d <= 0.0 && q_u + q_d < 0.0 {
                        p.dirn = -1;
                    }
                    if full_momentum && p.lambda > 0.0 && p.dirn != 0 {
                        let a_mu = cmax(self.area_mid_[uu], FUDGE);
                        let a_md = cmax(self.area_mid_[ud], FUDGE);
                        let vu = q_u / a_mu;
                        let vd = q_d / a_md;
                        p.dq4j = dt * p.sigma_j * (vd * vd * a_md - vu * vu * a_mu) / p.lambda;
                    }
                }
            }
            self.vjunc_[pi] = p;
        }
    }

    /// vjAccumulateResiduals: the post-step momentum-residual diagnostic (ctx.vj_diag).
    fn vj_accumulate_residuals(&mut self, ctx: &mut SimulationContext) {
        if self.vjunc_.is_empty() {
            return;
        }
        if ctx.vj_diag.node_idx.len() != self.vjunc_.len() {
            let d = &mut ctx.vj_diag;
            d.clear();
            for p in self.vjunc_.iter() {
                d.node_idx.push(p.node);
                d.up_link.push(p.up_link);
                d.dn_link.push(p.dn_link);
                d.resid_max.push(0.0);
                d.resid_sum.push(0.0);
                d.resid_n.push(0);
            }
        }
        for r in 0..self.vjunc_.len() {
            let p = self.vjunc_[r];
            if p.through == 0 || p.up_link < 0 || p.dn_link < 0 {
                continue;
            }
            let uu = p.up_link as usize;
            let ud = p.dn_link as usize;
            let a_u = self.area2_[uu];
            let a_d = self.area1_[ud];
            if a_u <= FUDGE || a_d <= FUDGE {
                continue;
            }
            let ucu = self.tile_uj_to_ci_[uu] as usize;
            let ucd = self.tile_uj_to_ci_[ud] as usize;
            let q_u = ctx.links.flow[uu] / self.tile_barrels_d_[ucu];
            let q_d = ctx.links.flow[ud] / self.tile_barrels_d_[ucd];
            let flux = q_u * q_u / a_u - q_d * q_d / a_d;
            let hyd = constants::GRAVITY * 0.5 * (a_u + a_d) * (self.depth2_[uu] - self.depth1_[ud]);
            let rr = flux + hyd;
            let a_r = rr.abs();
            let p = &mut self.vjunc_[r];
            p.resid_max = cmax(p.resid_max, a_r);
            p.resid_sum += a_r;
            p.resid_n += 1;
            ctx.vj_diag.resid_max[r] = p.resid_max;
            ctx.vj_diag.resid_sum[r] = p.resid_sum;
            ctx.vj_diag.resid_n[r] = p.resid_n;
        }
    }

    /// momentumKernels: the fused per-conduit Picard tail.
    fn momentum_kernels(&mut self, ctx: &mut SimulationContext, dt: f64, step: i32) {
        let static_slot_e = self.surcharge_method != SurchargeMethod::DynamicSlot;
        let do_losses = !self.losses_all_zero_;

        if !self.vjunc_.is_empty() {
            self.vj_prepare_iteration(ctx, dt);
        }
        // The unsteady-friction path runs the same loop serially; this one is serial too.
        for ci in 0..self.n_conduits_ as usize {
            let uci = ci;
            let j = self.tile_uj_[uci];
            let uj = j as usize;

            self.new_flow_[uj] = ctx.links.flow[uj];

            if self.bypassed_[uj] != 0 {
                continue;
            }

            if do_losses {
                self.recompute_conduit_loss_one(ctx, dt, ci);
            }

            self.dqdh_[uj] = 0.0;

            // STEP E: static slot / EXTRAN overrides
            if static_slot_e {
                let yf = self.tile_y_full_[uci];
                let af = self.tile_a_full_[uci];
                let rf = self.tile_r_full_[uci];
                let wm = self.tile_w_max_[uci];
                let shape = self.tile_shape_[uci];

                if self.depth1_[uj] >= yf {
                    self.area1_[uj] = af;
                    self.hrad1_[uj] = rf;
                }
                if self.depth2_[uj] >= yf {
                    self.area2_[uj] = af;
                }
                if self.depth_mid_[uj] >= yf {
                    self.area_mid_[uj] = af;
                    self.hrad_mid_[uj] = rf;
                }

                if self.depth1_[uj] > yf {
                    let w_slot = self.get_slot_width(self.depth1_[uj], yf, wm, shape);
                    self.area1_[uj] = af + (self.depth1_[uj] - yf) * w_slot;
                    self.hrad1_[uj] = rf;
                }
                if self.depth2_[uj] > yf {
                    let w_slot = self.get_slot_width(self.depth2_[uj], yf, wm, shape);
                    self.area2_[uj] = af + (self.depth2_[uj] - yf) * w_slot;
                }
                let y_mid = self.depth_mid_[uj];
                if y_mid > yf {
                    let w_slot = self.get_slot_width(y_mid, yf, wm, shape);
                    self.area_mid_[uj] = af + (y_mid - yf) * w_slot;
                    if w_slot > 0.0 {
                        self.width_mid_[uj] = w_slot;
                    }
                    self.hrad_mid_[uj] = rf;
                }
            }

            // Momentum category
            let fc = ctx.links.flow_class[uj];
            let a_mid = self.area_mid_[uj];
            let yf_c = self.tile_y_full_[uci];
            let is_full = self.depth1_[uj] >= yf_c && self.depth2_[uj] >= yf_c;

            let cat = if fc == FlowClass::Dry || fc == FlowClass::UpDry || fc == FlowClass::DnDry || a_mid <= FUDGE || self.tile_is_closed_[uci] != 0 || ctx.links.setting[uj] == 0.0 {
                MomentumCategory::SkipDry
            } else if self.tile_is_force_main_[uci] != 0 && is_full {
                if ctx.options.force_main_eqn == 1 {
                    MomentumCategory::ForceMainDw
                } else {
                    MomentumCategory::ForceMainHw
                }
            } else if self.tile_is_open_[uci] == 0 && is_full {
                MomentumCategory::ManningClosedFull
            } else if self.tile_is_open_[uci] != 0 {
                MomentumCategory::ManningOpen
            } else {
                MomentumCategory::ManningClosedFs
            };
            self.category_[uj] = cat;

            match cat {
                MomentumCategory::SkipDry => self.process_dry_link(ctx, dt, uj),
                MomentumCategory::ManningOpen | MomentumCategory::ManningClosedFs | MomentumCategory::ManningClosedFull => self.process_manning_link(ctx, dt, step, uj, cat),
                MomentumCategory::ForceMainHw | MomentumCategory::ForceMainDw => self.process_force_main_link(ctx, dt, step, uj, cat),
                _ => {}
            }

            // (storage_area chamber volume: refused in execute)

            ctx.links.flow[uj] = self.new_flow_[uj];
            ctx.links.dqdh[uj] = self.dqdh_[uj];
        }
    }

    /// processDryLink: zero flow, minimal bookkeeping.
    fn process_dry_link(&mut self, ctx: &mut SimulationContext, _dt: f64, uj: usize) {
        let dt_g = self.dt_gravity_;
        let uci = self.tile_uj_to_ci_[uj] as usize;
        let a_mid = self.area_mid_[uj];
        let barrels_d = self.tile_barrels_d_[uci];
        self.area_mid_[uj] = 0.5 * (self.area1_[uj] + self.area2_[uj]);
        self.dqdh_[uj] = dt_g * a_mid / self.tile_length_[uci] * barrels_d;
        self.froude_[uj] = 0.0;
        self.new_flow_[uj] = 0.0;
        self.q1_[uj] = 0.0;
        let yf = self.tile_y_full_[uci];
        ctx.links.depth[uj] = cmin(self.depth_mid_[uj], yf);
        ctx.links.volume[uj] = self.area_mid_[uj] * self.tile_links_length_[uci] * barrels_d;
    }

    /// applyFlowLimits: inlet control, normal-flow limit, relaxation, flap gates, dry nodes,
    /// depth/volume update.
    #[allow(clippy::too_many_arguments)]
    fn apply_flow_limits(&mut self, ctx: &mut SimulationContext, _dt: f64, step: i32, uj: usize, q: &mut f64, q_last: f64, barrels_d: f64, is_full: bool) {
        let normal_flow_ltd = ctx.options.normal_flow_ltd;
        let uci = self.tile_uj_to_ci_[uj] as usize;
        let yf = self.tile_y_full_[uci];
        let h1 = self.h1_[uj];
        let n1 = self.tile_n1_[uci];
        let n2 = self.tile_n2_[uci];
        let un1 = n1 as usize;
        let un2 = n2 as usize;

        let ucr = self.tile_crow_[uci] as usize;
        ctx.link_subtypes.conduits.inlet_control[ucr] = 0;
        ctx.link_subtypes.conduits.normal_flow_limited[ucr] = 0;
        if *q > 0.0 {
            if self.tile_culvert_code_[uci] > 0 && !is_full {
                let y_in = h1 - (self.tile_inv1_elev_[uci] + self.tile_z1_off_[uci]);
                let mut dqdh_culv = 0.0;
                let mut inlet = false;
                let q_inlet = {
                    let xs_c = build_xsp(ctx, uj);
                    culvert::get_inflow(*q, y_in, &xs_c, self.tile_slope_[uci], self.tile_culvert_code_[uci], &mut dqdh_culv, &mut inlet)
                };
                if inlet {
                    *q = q_inlet;
                    self.dqdh_[uj] = dqdh_culv;
                    ctx.link_subtypes.conduits.inlet_control[ucr] = 1;
                }
            } else {
                let fc2 = ctx.links.flow_class[uj];
                let nfl = normal_flow_ltd;
                if nfl != 3 && self.depth1_[uj] < yf && (fc2 == FlowClass::Subcritical || fc2 == FlowClass::Supercritical) {
                    let mut has_outfall = false;
                    if n1 >= 0 && n2 >= 0 {
                        has_outfall = ctx.nodes.type_[un1] == NodeType::Outfall || ctx.nodes.type_[un2] == NodeType::Outfall;
                    }
                    let slope_check = (nfl == 0 || nfl == 2 || has_outfall) && (self.depth1_[uj] < self.depth2_[uj]);
                    let mut froude_check = false;
                    if !slope_check && (nfl == 1 || nfl == 2) && !has_outfall && self.depth1_[uj] > FUDGE && self.depth2_[uj] > FUDGE {
                        let y1d = self.depth1_[uj];
                        let closed_nearfull = self.tile_is_open_[uci] == 0 && (yf - y1d <= FUDGE);
                        if !closed_nearfull {
                            let v1 = *q / self.area1_[uj];
                            let xs = build_xsp(ctx, uj);
                            let w1 = xsect::get_w_of_y(&xs, y1d);
                            let dh1 = self.area1_[uj] / w1;
                            let f1 = v1.abs() / (constants::GRAVITY * dh1).sqrt();
                            froude_check = f1 >= 1.0;
                        }
                    }
                    if slope_check || froude_check {
                        let q_norm = self.tile_beta_[uci] * self.area1_[uj] * self.hrad1_[uj].powf(2.0 / 3.0); // LIBM pow
                        if q_norm < *q {
                            *q = q_norm;
                            ctx.link_subtypes.conduits.normal_flow_limited[ucr] = 1;
                        }
                    }
                }
            }
        }

        // Under-relaxation (after first iteration)
        if step > 0 {
            *q = (1.0 - self.omega) * q_last + self.omega * *q;
            if *q * q_last < 0.0 {
                *q = 0.001 * (if *q > 0.0 { 1.0 } else { -1.0 });
            }
        }

        // Flow limit
        if self.tile_q_limit_[uci] > 0.0 && q.abs() > self.tile_q_limit_[uci] {
            *q = (if *q > 0.0 { 1.0 } else { -1.0 }) * self.tile_q_limit_[uci];
        }

        // Flap gate checks
        if self.tile_has_flap_gate_[uci] != 0 && *q * (self.tile_direction_[uci] as f64) < 0.0 {
            *q = 0.0;
        }
        if *q < 0.0 && n2 >= 0 && ctx.nodes.type_[un2] == NodeType::Outfall {
            let ofr = ctx.node_subtypes.outfall_row(n2);
            if ofr >= 0 && ctx.node_subtypes.outfalls.has_flap_gate[ofr as usize] != 0 {
                *q = 0.0;
            }
        }
        if *q > 0.0 && n1 >= 0 && ctx.nodes.type_[un1] == NodeType::Outfall {
            let ofr = ctx.node_subtypes.outfall_row(n1);
            if ofr >= 0 && ctx.node_subtypes.outfalls.has_flap_gate[ofr as usize] != 0 {
                *q = 0.0;
            }
        }

        // Dry node check
        if *q > FUDGE && lidnode::port_depth(ctx, n1, self.tile_z1_off_[uci]) <= FUDGE {
            *q = FUDGE;
        }
        if *q < -FUDGE && lidnode::port_depth(ctx, n2, self.tile_z2_off_[uci]) <= FUDGE {
            *q = -FUDGE;
        }

        // Save new flow
        self.new_flow_[uj] = *q * barrels_d;
        self.q1_[uj] = *q;

        // Update link depth and volume
        ctx.links.depth[uj] = cmin(self.depth_mid_[uj], yf);
        let a_mid_avg = (self.area1_[uj] + self.area2_[uj]) * 0.5;
        ctx.links.volume[uj] = a_mid_avg * self.tile_links_length_[uci] * barrels_d;

        let a_full = ctx.links.xsect_a_full[uj];
        let up_full = self.area1_[uj] >= a_full;
        let dn_full = self.area2_[uj] >= a_full;
        let mid_full = self.area_mid_[uj] >= a_full;
        ctx.link_subtypes.conduits.full_state[ucr] = ((if up_full { 1 } else { 0 }) | (if dn_full { 2 } else { 0 }) | (if mid_full { 4 } else { 0 })) as i8;
    }

    /// processManningLink: MANNING_OPEN, MANNING_CLOSED_FS and MANNING_CLOSED_FULL.
    fn process_manning_link(&mut self, ctx: &mut SimulationContext, dt: f64, step: i32, uj: usize, cat: MomentumCategory) {
        let inert_damping = ctx.options.inertial_damping;
        let dt_g = self.dt_gravity_;
        let is_closed_full = cat == MomentumCategory::ManningClosedFull;
        let is_open_cat = cat == MomentumCategory::ManningOpen;

        let uci = self.tile_uj_to_ci_[uj] as usize;
        let barrels_d = self.tile_barrels_d_[uci];
        let length = self.tile_length_[uci];
        let a_mid = self.area_mid_[uj];
        let r_mid = self.hrad_mid_[uj];
        let q_last = self.q1_[uj];
        let yf = self.tile_y_full_[uci];
        let is_full = self.depth1_[uj] >= yf && self.depth2_[uj] >= yf;

        let a_full = self.tile_a_full_[uci];
        let a1_conv = convey_area(self.area1_[uj], a_full);
        let a2_conv = convey_area(self.area2_[uj], a_full);
        let a_mid_conv = convey_area(a_mid, a_full);

        // Velocity (clamped)
        let mut v = q_last / a_mid_conv;
        let mut absv = v.abs();
        if absv > MAX_VELOCITY {
            v = if v > 0.0 { MAX_VELOCITY } else { -MAX_VELOCITY };
            absv = MAX_VELOCITY;
        }
        self.velocity_[uj] = v;

        let mut fr = 0.0;
        let mut sig = 0.0;
        if !is_closed_full {
            let mut w_mid = self.width_mid_[uj];
            if ctx.links.flow_class[uj] == FlowClass::Subcritical && self.wcap_dm_[uj] < self.depth_mid_[uj] {
                let xs = build_xsp(ctx, uj);
                w_mid = xsect::get_w_of_y(&xs, self.depth_mid_[uj]);
            }
            let closed_nearfull = self.tile_is_open_[uci] == 0 && (yf - self.depth_mid_[uj] <= FUDGE);
            if self.depth_mid_[uj] > FUDGE && !closed_nearfull {
                let mut a_fr = a_mid;
                if self.depth_mid_[uj] > yf {
                    a_fr = xsect::get_a_of_y(&build_xsp(ctx, uj), self.depth_mid_[uj]);
                }
                let dh = if w_mid > 0.0 { a_fr / w_mid } else { 0.0 };
                fr = if dh > 0.0 { absv / (constants::GRAVITY * dh).sqrt() } else { 0.0 };
            }

            let fc = ctx.links.flow_class[uj];
            if fc == FlowClass::Subcritical && fr > 1.0 {
                ctx.links.flow_class[uj] = FlowClass::Supercritical;
            }

            sig = cmax(0.0, cmin(1.0, 2.0 * (1.0 - fr)));
        }
        self.froude_[uj] = fr;

        let h1 = self.h1_[uj];
        let h2 = self.h2_[uj];

        let r1_val = self.hrad1_[uj];
        let mut rho = 1.0;
        if !is_closed_full && !is_full && q_last > 0.0 && h1 >= h2 {
            rho = sig;
        }
        let a_wtd = a1_conv + (a_mid_conv - a1_conv) * rho;
        let r_wtd = r1_val + (r_mid - r1_val) * rho;

        // Virtual-junction pair coupling (mechanisms 3-4); SWMM_VJ_M3 / SWMM_VJ_M4 unset: on
        let mut a_wtd = a_wtd;
        let mut r_wtd = r_wtd;
        let mut vj_dq4 = 0.0;
        if !self.vjunc_.is_empty() {
            let r1 = self.vj_pair_n1_[uj];
            let r2 = self.vj_pair_n2_[uj];
            if r1 >= 0 {
                let p = self.vjunc_[r1 as usize];
                if p.active != 0 {
                    if p.dirn > 0 && uj as i32 == p.dn_link && !is_closed_full && !is_full && q_last > 0.0 {
                        a_wtd = p.a_up_mid + (a_mid_conv - p.a_up_mid) * rho;
                        r_wtd = p.r_up_mid + (r_mid - p.r_up_mid) * rho;
                    }
                    vj_dq4 += p.dq4j;
                }
            }
            if r2 >= 0 {
                let p = self.vjunc_[r2 as usize];
                if p.active != 0 {
                    vj_dq4 += p.dq4j;
                }
            }
        }

        if !is_closed_full {
            if inert_damping == 0 {
                sig = 1.0;
            } else if inert_damping == 2 {
                sig = 0.0;
            }
            if !is_open_cat && is_full {
                sig = 0.0;
            }
        }
        self.sigma_[uj] = sig;

        let r43 = r_wtd.powf(1.33333); // LIBM pow
        let dq1 = dt * self.tile_rough_factor_[uci] / r43 * absv;

        let dq2 = dt_g * a_wtd * (h2 - h1) / length;

        let a_old = cmax(convey_area(self.area_old_[uj], a_full), FUDGE);
        let mut dq3 = 0.0;
        let mut dq4 = 0.0;
        if sig > 0.0 {
            dq3 = 2.0 * v * (a_mid_conv - a_old) * sig;
            if length > 0.0 {
                dq4 = dt * v * v * (a2_conv - a1_conv) / length * sig;
            }
            dq4 += vj_dq4;
        }

        let mut dq5 = 0.0;
        if self.tile_has_losses_[uci] != 0 {
            let absq = q_last.abs();
            let mut losses = 0.0;
            if a1_conv > FUDGE {
                losses += self.tile_loss_inlet_[uci] * (absq / a1_conv);
            }
            if a2_conv > FUDGE {
                losses += self.tile_loss_outlet_[uci] * (absq / a2_conv);
            }
            if a_mid_conv > FUDGE {
                losses += self.tile_loss_avg_[uci] * (absq / a_mid_conv);
            }
            dq5 = losses / 2.0 / length * dt;
        }

        let mut dq6 = 0.0;
        {
            let conduit_length = self.tile_links_length_[uci];
            let ucr6 = self.tile_crow_[uci] as usize;
            let loss_rate = ctx.link_subtypes.conduits.evap_loss_rate[ucr6] + ctx.link_subtypes.conduits.seep_loss_rate[ucr6];
            if loss_rate > 0.0 && conduit_length > 0.0 {
                dq6 = loss_rate * 2.5 * dt * v / conduit_length;
            }
        }

        let q_old = ctx.links.old_flow[uj] / barrels_d;
        let mut denom = 1.0 + dq1 + dq5;

        // Unsteady friction (issue #156 Phase 3)
        let mut uf_applied = false;
        let mut uf_num = 0.0;
        if self.unsteady_friction != 0 && self.uf_k3 > 0.0 {
            const K_UF_VEL_FLOOR: f64 = 0.01;
            let v_old = if a_old > FUDGE { q_old / a_old } else { 0.0 };
            if v_old.abs() >= K_UF_VEL_FLOOR || absv >= K_UF_VEL_FLOOR {
                let w_mid = self.width_mid_[uj];
                let c_wave = if w_mid > FUDGE && a_mid > FUDGE { (GRAVITY * a_mid / w_mid).sqrt() } else { 0.0 };
                let mut dv_dx = 0.0;
                {
                    let nb_vel = |s: &DWSolver, nbj: i32, sg: i8| -> f64 {
                        let unb = nbj as usize;
                        let uci_nb = s.tile_uj_to_ci_[unb] as usize;
                        let a_nb = cmax(convey_area(s.area_mid_[unb], s.tile_a_full_[uci_nb]), FUDGE);
                        (sg as f64) * (ctx.links.flow[unb] / s.tile_barrels_d_[uci_nb]) / a_nb
                    };
                    let mut v_l = v;
                    let mut v_r = v;
                    let mut dsum = 0.0;
                    if self.uf_nb_up_[uj] >= 0 {
                        let uci_nb = self.tile_uj_to_ci_[self.uf_nb_up_[uj] as usize] as usize;
                        v_l = nb_vel(self, self.uf_nb_up_[uj], self.uf_sg_up_[uj]);
                        dsum += 0.5 * (length + self.tile_length_[uci_nb]);
                    }
                    if self.uf_nb_dn_[uj] >= 0 {
                        let uci_nb = self.tile_uj_to_ci_[self.uf_nb_dn_[uj] as usize] as usize;
                        v_r = nb_vel(self, self.uf_nb_dn_[uj], self.uf_sg_dn_[uj]);
                        dsum += 0.5 * (length + self.tile_length_[uci_nb]);
                    }
                    if dsum > 0.0 {
                        dv_dx = (v_r - v_l).abs() / dsum;
                    } else if length > 0.0 {
                        let absq = q_last.abs();
                        let v1e = if a1_conv > FUDGE { absq / a1_conv } else { 0.0 };
                        let v2e = if a2_conv > FUDGE { absq / a2_conv } else { 0.0 };
                        dv_dx = (v2e - v1e).abs() / length;
                    }
                }
                let sgn = if v > 0.0 {
                    1.0
                } else if v < 0.0 {
                    -1.0
                } else {
                    0.0
                };
                let mut dq_grad = dt * self.uf_k3 * a_wtd * c_wave * sgn * dv_dx;
                let cap = 0.5 * q_old.abs();
                if dq_grad.abs() > cap {
                    dq_grad = if dq_grad > 0.0 { cap } else { -cap };
                }
                denom += self.uf_k3;
                uf_num = self.uf_k3 * a_wtd * v_old - dq_grad;
                uf_applied = true;
            }
        }
        let mut q = if uf_applied { (q_old - dq2 + dq3 + dq4 + dq6 + uf_num) / denom } else { (q_old - dq2 + dq3 + dq4 + dq6) / denom };
        self.dqdh_[uj] = 1.0 / denom * GRAVITY * dt * a_wtd / length * barrels_d;

        self.apply_flow_limits(ctx, dt, step, uj, &mut q, q_last, barrels_d, is_full);
    }

    /// processForceMainLink: FORCE_MAIN_HW and FORCE_MAIN_DW (always full).
    fn process_force_main_link(&mut self, ctx: &mut SimulationContext, dt: f64, step: i32, uj: usize, cat: MomentumCategory) {
        let dt_g = self.dt_gravity_;
        let is_dw = cat == MomentumCategory::ForceMainDw;

        let uci = self.tile_uj_to_ci_[uj] as usize;
        let barrels_d = self.tile_barrels_d_[uci];
        let a_mid = self.area_mid_[uj];
        let r_mid = self.hrad_mid_[uj];
        let q_last = self.q1_[uj];

        let a_full = self.tile_a_full_[uci];
        let a1_conv = convey_area(self.area1_[uj], a_full);
        let a2_conv = convey_area(self.area2_[uj], a_full);
        let a_mid_conv = convey_area(a_mid, a_full);

        let is_full = true;

        let mut v = q_last / a_mid_conv;
        let mut absv = v.abs();
        if absv > MAX_VELOCITY {
            v = if v > 0.0 { MAX_VELOCITY } else { -MAX_VELOCITY };
            absv = MAX_VELOCITY;
        }
        self.velocity_[uj] = v;

        self.froude_[uj] = 0.0;
        self.sigma_[uj] = 0.0;

        let h1 = self.h1_[uj];
        let h2 = self.h2_[uj];

        let a_wtd = a1_conv + (a_mid_conv - a1_conv);

        let sbot = self.tile_fm_sbot_[uci];
        let dq1 = if is_dw {
            let re = 4.0 * r_mid * absv / force_main::VISCOS;
            let f = force_main::get_fric_factor(self.tile_fm_rbot_[uci], r_mid, re);
            dt * (f * sbot * absv / r_mid)
        } else {
            dt * (sbot * absv.powf(0.852) / r_mid.powf(1.1667)) // LIBM pow
        };

        let fm_length = self.tile_length_[uci];
        let dq2 = dt_g * a_wtd * (h2 - h1) / fm_length;

        let mut dq5 = 0.0;
        if self.has_losses_[uj] != 0 {
            let absq = q_last.abs();
            let mut losses = 0.0;
            if a1_conv > FUDGE {
                losses += self.tile_loss_inlet_[uci] * (absq / a1_conv);
            }
            if a2_conv > FUDGE {
                losses += self.tile_loss_outlet_[uci] * (absq / a2_conv);
            }
            if a_mid_conv > FUDGE {
                losses += self.tile_loss_avg_[uci] * (absq / a_mid_conv);
            }
            dq5 = losses / 2.0 / fm_length * dt;
        }

        let mut dq6 = 0.0;
        {
            let conduit_length = self.tile_links_length_[uci];
            let ucr6 = self.tile_crow_[uci] as usize;
            let loss_rate = ctx.link_subtypes.conduits.evap_loss_rate[ucr6] + ctx.link_subtypes.conduits.seep_loss_rate[ucr6];
            if loss_rate > 0.0 && conduit_length > 0.0 {
                dq6 = loss_rate * 2.5 * dt * v / conduit_length;
            }
        }

        let q_old = ctx.links.old_flow[uj] / barrels_d;
        let mut denom = 1.0 + dq1 + dq5;

        let mut uf_applied = false;
        if self.unsteady_friction != 0 && self.uf_k3 > 0.0 {
            const K_UF_VEL_FLOOR: f64 = 0.01;
            let v_old = if a_wtd > FUDGE { q_old / a_wtd } else { 0.0 };
            if v_old.abs() >= K_UF_VEL_FLOOR || v.abs() >= K_UF_VEL_FLOOR {
                denom += self.uf_k3;
                uf_applied = true;
            }
        }
        let mut q = if uf_applied { (q_old - dq2 + dq6 + self.uf_k3 * q_old) / denom } else { (q_old - dq2 + dq6) / denom };
        self.dqdh_[uj] = 1.0 / denom * GRAVITY * dt * a_wtd / fm_length * barrels_d;

        self.apply_flow_limits(ctx, dt, step, uj, &mut q, q_last, barrels_d, is_full);
    }

    /// updateNodeFlows: the serial all-links scatter (only on the no-callback path).
    fn update_node_flows(&mut self, ctx: &mut SimulationContext) {
        for j in 0..self.n_links_ {
            let uj = j as usize;
            let accepted = bound_surface_outfall_flow(ctx, j, self.new_flow_[uj]);
            if accepted != self.new_flow_[uj] {
                self.new_flow_[uj] = accepted;
                self.dqdh_[uj] = 0.0;
                ctx.links.dqdh[uj] = 0.0;
            }
            ctx.links.flow[uj] = self.new_flow_[uj];

            let un1 = ctx.links.node1[uj] as usize;
            let un2 = ctx.links.node2[uj] as usize;
            let q = self.new_flow_[uj];
            let nodes = &mut ctx.nodes;
            if q >= 0.0 {
                nodes.outflow[un1] += q;
                nodes.inflow[un2] += q;
            } else {
                nodes.inflow[un1] -= q;
                nodes.outflow[un2] -= q;
            }
            self.xnode_.sumdqdh[un1] += self.dqdh_[uj];
            self.xnode_.sumdqdh[un2] += self.dqdh_[uj];

            let sa1 = self.surf_area1_[uj];
            let sa2 = self.surf_area2_[uj];

            if ctx.links.type_[uj] == LinkType::Conduit {
                let uci_l = self.tile_uj_to_ci_[uj] as usize;
                let ucr_l = self.tile_crow_[uci_l] as usize;
                let mut conduit_loss = (ctx.link_subtypes.conduits.evap_loss_rate[ucr_l] + ctx.link_subtypes.conduits.seep_loss_rate[ucr_l]) * self.tile_barrels_d_[uci_l];
                let nodes = &mut ctx.nodes;
                if conduit_loss < 0.0 {
                    let mut gain = -conduit_loss;
                    if nodes.type_[un1] != NodeType::Outfall && nodes.type_[un2] != NodeType::Outfall {
                        gain /= 2.0;
                    }
                    if nodes.type_[un1] != NodeType::Outfall {
                        nodes.inflow[un1] += gain;
                    }
                    if nodes.type_[un2] != NodeType::Outfall {
                        nodes.inflow[un2] += gain;
                    }
                } else if conduit_loss > 0.0 {
                    if nodes.type_[un1] != NodeType::Outfall && nodes.type_[un2] != NodeType::Outfall {
                        conduit_loss /= 2.0;
                    }
                    if nodes.type_[un1] != NodeType::Outfall {
                        nodes.outflow[un1] += conduit_loss;
                    }
                    if nodes.type_[un2] != NodeType::Outfall {
                        nodes.outflow[un2] += conduit_loss;
                    }
                }
            }

            let ci_b = self.tile_uj_to_ci_[uj];
            let barrels = if ci_b >= 0 { self.tile_barrels_d_[ci_b as usize] as i32 } else { 1 };
            let b = barrels as f64;
            self.xnode_.new_surf_area[un1] += sa1 * b;
            self.xnode_.new_surf_area[un2] += sa2 * b;
        }
    }

    /// gatherConduitNodeFlows: the node-centric conduit scatter (ascending link order per node).
    fn gather_conduit_node_flows(&mut self, ctx: &mut SimulationContext) {
        let cd = &ctx.link_subtypes.conduits;
        let nodes = &mut ctx.nodes;
        for n in 0..self.n_nodes_ {
            let un = n as usize;
            let b0 = self.csr_row_[un];
            let b1 = self.csr_row_[un + 1];
            if b0 == b1 {
                continue;
            }
            let self_outfall = self.node_tile_[un].is_outfall != 0;
            let mut inflow = nodes.inflow[un];
            let mut outflow = nodes.outflow[un];
            let mut sumdqdh = self.xnode_.sumdqdh[un];
            let mut surf = self.xnode_.new_surf_area[un];

            for k in b0..b1 {
                let uk = k as usize;
                let uj = self.csr_link_[uk] as usize;
                let at_n2 = self.csr_is_n2_[uk] != 0;
                let uci = self.tile_uj_to_ci_[uj] as usize;

                let q = self.new_flow_[uj];
                if q >= 0.0 {
                    if at_n2 {
                        inflow += q;
                    } else {
                        outflow += q;
                    }
                } else if at_n2 {
                    outflow -= q;
                } else {
                    inflow -= q;
                }

                sumdqdh += self.dqdh_[uj];

                let ucr_g = self.tile_crow_[uci] as usize;
                let mut conduit_loss = (cd.evap_loss_rate[ucr_g] + cd.seep_loss_rate[ucr_g]) * self.tile_barrels_d_[uci];
                if conduit_loss < 0.0 {
                    let mut gain = -conduit_loss;
                    if !self_outfall && self.csr_other_outfall_[uk] == 0 {
                        gain /= 2.0;
                    }
                    if !self_outfall {
                        inflow += gain;
                    }
                } else if conduit_loss > 0.0 {
                    if !self_outfall && self.csr_other_outfall_[uk] == 0 {
                        conduit_loss /= 2.0;
                    }
                    if !self_outfall {
                        outflow += conduit_loss;
                    }
                }

                let barrels = self.tile_barrels_d_[uci] as i32;
                surf += (if at_n2 { self.surf_area2_[uj] } else { self.surf_area1_[uj] }) * (barrels as f64);
            }

            nodes.inflow[un] = inflow;
            nodes.outflow[un] = outflow;
            self.xnode_.sumdqdh[un] = sumdqdh;
            self.xnode_.new_surf_area[un] = surf;
        }
    }

    /// limitConduitLosses: budget the incident conduit losses against what each node holds.
    fn limit_conduit_losses(&mut self, ctx: &mut SimulationContext, dt: f64) {
        let adjustable = |ctx: &SimulationContext, j: i32, r: i32| -> bool {
            let modes = &ctx.forcing.link_seepage_mode;
            let cd = &ctx.link_subtypes.conduits;
            ((j as usize) >= modes.len() || (modes[j as usize] != ForcingMode::Override && modes[j as usize] != ForcingMode::Add)) && cd.evap_loss_rate[r as usize] + cd.seep_loss_rate[r as usize] > 0.0
        };
        self.xnode_.conduit_loss.iter_mut().for_each(|x| *x = 0.0);
        for ci in 0..self.n_conduits_ as usize {
            let j = self.tile_uj_[ci];
            let r = self.tile_crow_[ci];
            if !adjustable(ctx, j, r) {
                continue;
            }
            let n1 = self.tile_n1_[ci] as usize;
            let n2 = self.tile_n2_[ci] as usize;
            let cd = &ctx.link_subtypes.conduits;
            let rate = (cd.evap_loss_rate[r as usize] + cd.seep_loss_rate[r as usize]) * self.tile_barrels_d_[ci];
            let share = if self.node_tile_[n1].is_outfall == 0 && self.node_tile_[n2].is_outfall == 0 { 0.5 } else { 1.0 };
            if self.node_tile_[n1].is_outfall == 0 {
                self.xnode_.conduit_loss[n1] += share * rate;
            }
            if self.node_tile_[n2].is_outfall == 0 {
                self.xnode_.conduit_loss[n2] += share * rate;
            }
        }
        for i in 0..self.n_nodes_ as usize {
            self.xnode_.loss_factor[i] = 1.0;
            let rate = self.xnode_.conduit_loss[i];
            if rate <= 0.0 || self.node_tile_[i].is_outfall != 0 {
                continue;
            }
            let nodes = &ctx.nodes;
            let area = if nodes.storage_area[i] > 0.0 { self.xnode_.new_surf_area[i] } else { cmax(self.xnode_.new_surf_area[i], self.min_surf_area_) };
            let storage_rate = cmax(0.0, nodes.old_depth[i] - FUDGE) * area / dt;
            // (storage_area > 0 chamber volume: refused in execute)
            let net = nodes.inflow[i] - nodes.outflow[i] + rate;
            let available = cmin(net + storage_rate, net + nodes.old_net_inflow[i] + 2.0 * storage_rate);
            self.xnode_.loss_factor[i] = cmin(1.0, cmax(0.0, available / rate));
        }
        for ci in 0..self.n_conduits_ as usize {
            let j = self.tile_uj_[ci];
            let r = self.tile_crow_[ci];
            if !adjustable(ctx, j, r) {
                continue;
            }
            let ur = r as usize;
            let n1 = self.tile_n1_[ci] as usize;
            let n2 = self.tile_n2_[ci] as usize;
            let factor = cmin(self.xnode_.loss_factor[n1], self.xnode_.loss_factor[n2]);
            if factor >= 1.0 {
                continue;
            }
            let cd = &mut ctx.link_subtypes.conduits;
            let rate = cd.evap_loss_rate[ur] + cd.seep_loss_rate[ur];
            let positive = cd.evap_loss_rate[ur] + cmax(0.0, cd.seep_loss_rate[ur]);
            let scale = if cd.seep_loss_rate[ur] >= 0.0 { factor } else { (positive - rate * (1.0 - factor)) / positive };
            cd.evap_loss_rate[ur] *= scale;
            if cd.seep_loss_rate[ur] > 0.0 {
                cd.seep_loss_rate[ur] *= scale;
            }
            let reduction = (rate * self.tile_barrels_d_[ci]) * (1.0 - factor);
            let share = if self.node_tile_[n1].is_outfall == 0 && self.node_tile_[n2].is_outfall == 0 { 0.5 } else { 1.0 };
            if self.node_tile_[n1].is_outfall == 0 {
                ctx.nodes.outflow[n1] -= share * reduction;
            }
            if self.node_tile_[n2].is_outfall == 0 {
                ctx.nodes.outflow[n2] -= share * reduction;
            }
        }
    }

    /// computeAASkipFlags: nodes where Anderson acceleration must not mix.
    fn compute_aa_skip_flags(&mut self, ctx: &SimulationContext) {
        if !self.anderson_accel {
            return;
        }
        let links = &ctx.links;
        self.aa_skip_.iter_mut().for_each(|x| *x = 0);

        if self.surcharge_method == SurchargeMethod::Extran && self.node_continuity == NodeContinuity::Explicit {
            for i in 0..self.n_nodes_ as usize {
                if self.xnode_.is_surcharged[i] != 0 {
                    self.aa_skip_[i] = 1;
                }
            }
        }
        // (TPA and DYNAMIC_SLOT walks: refused in execute)
        if self.surcharge_method == SurchargeMethod::Slot {
            for ci in 0..self.n_conduits_ as usize {
                let j = self.conduit_idx_[ci];
                let uj = j as usize;
                if self.is_open_[uj] != 0 {
                    continue;
                }
                let yf = links.xsect_y_full[uj];
                if yf <= 0.0 {
                    continue;
                }
                let ratio = self.depth_mid_[uj] / yf;
                if (0.98..=1.02).contains(&ratio) {
                    self.aa_skip_[links.node1[uj] as usize] = 1;
                    self.aa_skip_[links.node2[uj] as usize] = 1;
                }
            }
        }

        let nodes = &ctx.nodes;
        for j in 0..self.n_links_ {
            let uj = j as usize;
            let lt = links.type_[uj];
            if lt != LinkType::Weir && lt != LinkType::Orifice {
                continue;
            }
            let n1 = links.node1[uj];
            let n2 = links.node2[uj];
            if n1 < 0 || n2 < 0 {
                continue;
            }
            let un1 = n1 as usize;
            let un2 = n2 as usize;
            let y_full = links.xsect_y_full[uj];
            let setting = links.setting[uj];
            if y_full <= 0.0 || setting <= 0.0 {
                continue;
            }
            let hcrown = if lt == LinkType::Weir {
                let wr = ctx.link_subtypes.weir_row(j);
                let crest = if wr >= 0 { ctx.link_subtypes.weirs.crest_height[wr as usize] } else { 0.0 };
                nodes.invert_elev[un1] + crest + y_full
            } else {
                nodes.invert_elev[un1] + links.offset1[uj] + y_full * setting
            };
            let hgl1 = nodes.depth[un1] + nodes.invert_elev[un1];
            let hgl2 = nodes.depth[un2] + nodes.invert_elev[un2];
            let hgl_up = cmax(hgl1, hgl2);
            if hgl_up >= hcrown {
                self.aa_skip_[un1] = 1;
                self.aa_skip_[un2] = 1;
            }
        }

        for j in 0..self.n_links_ as usize {
            if links.type_[j] != LinkType::Pump {
                continue;
            }
            let n1 = links.node1[j];
            let n2 = links.node2[j];
            if n1 >= 0 {
                self.aa_skip_[n1 as usize] = 1;
            }
            if n2 >= 0 {
                self.aa_skip_[n2 as usize] = 1;
            }
        }
    }

    /// updateNodeDepthsTeam: setNodeDepth for every node, Anderson mixing, convergence tally.
    fn update_node_depths_team(&mut self, ctx: &mut SimulationContext, dt: f64, step: i32, unconv_shared: &mut i32) {
        let use_anderson = self.anderson_accel;
        let mut local_unconv = 0;
        for i in 0..self.n_nodes_ {
            let ui = i as usize;
            if self.node_tile_[ui].is_outfall != 0 {
                self.xnode_.converged[ui] = 0;
                continue;
            }

            let y_last = ctx.nodes.depth[ui];
            self.set_node_depth(ctx, i, dt, step);
            let g_k = ctx.nodes.depth[ui];

            if use_anderson && step >= 1 && self.aa_skip_[ui] == 0 && ctx.nodes.storage_area[ui] == 0.0 && !lidnode::active(ctx, i) {
                let r_k = g_k - y_last;
                if r_k.abs() <= 20.0 * self.head_tol {
                    let r_km1 = self.aa_r_prev_[ui];
                    let dr = r_k - r_km1;
                    let dr2 = dr * dr;
                    if dr2 > 1e-30 {
                        let alpha = cmax(0.0, cmin(1.0, r_k * dr / dr2));
                        let y_anderson = (1.0 - alpha) * g_k + alpha * self.aa_g_prev_[ui];
                        if y_anderson >= 0.0 {
                            let d_q = ctx.nodes.inflow[ui] - ctx.nodes.outflow[ui];
                            let d_v = 0.5 * (ctx.nodes.old_net_inflow[ui] + d_q) * dt;
                            self.commit_node_depth_state(ctx, i, y_anderson, d_v, dt);
                        }
                    }
                }
            }

            if use_anderson || ctx.nodes.storage_area[ui] > 0.0 {
                self.aa_y_prev_[ui] = y_last;
                self.aa_g_prev_[ui] = g_k;
                self.aa_r_prev_[ui] = g_k - y_last;
            }

            let raw_residual = (g_k - y_last).abs();
            let accepted_movement = (ctx.nodes.depth[ui] - y_last).abs();
            let tolerance = if ctx.nodes.storage_area[ui] > 0.0 { cmin(self.head_tol, 1.0e-7) } else { self.head_tol };
            let conv_flag: u8 = if raw_residual <= tolerance && accepted_movement <= tolerance { 1 } else { 0 };
            self.xnode_.converged[ui] = conv_flag;
            if conv_flag == 0 {
                local_unconv += 1;
            }
        }
        if local_unconv != 0 {
            *unconv_shared += local_unconv;
        }
    }

    /// setNodeDepth: one node's depth update (EXTRAN surcharge algorithm).
    fn set_node_depth(&mut self, ctx: &mut SimulationContext, node_idx: i32, dt: f64, step: i32) {
        let ui = node_idx as usize;
        let t = self.node_tile_[ui];
        let is_storage = t.is_storage != 0;

        let y_old = ctx.nodes.old_depth[ui];
        let y_last = ctx.nodes.depth[ui];
        let full_depth = t.full_depth;
        let y_crown = t.y_crown;

        let is_coupled = ui < ctx.coupled_node.len() && ctx.coupled_node[ui] != 0;
        let can_pond = (ctx.options.allow_ponding || is_coupled) && (t.ponded_area > 0.0);
        let is_ponded = can_pond && y_last > full_depth;

        ctx.nodes.overflow[ui] = 0.0;
        let mut surf_area = self.xnode_.new_surf_area[ui];
        if t.is_virtual == 0 && ctx.nodes.storage_area[ui] == 0.0 {
            surf_area = cmax(surf_area, self.min_surf_area_);
        } else if t.is_virtual != 0 && ctx.nodes.lat_flow[ui] != 0.0 && ui < self.vj_wet_floor_.len() {
            surf_area = cmax(surf_area, self.vj_wet_floor_[ui]);
        }

        let d_q = ctx.nodes.inflow[ui] - ctx.nodes.outflow[ui];
        let d_v = 0.5 * (ctx.nodes.old_net_inflow[ui] + d_q) * dt;
        // (storage_area > 0 chamber branch and the LID-node volume branch: refused upstream of
        // here, in execute and step_routing)

        let mut is_surcharged = false;
        if self.surcharge_method == SurchargeMethod::Extran || self.surcharge_method == SurchargeMethod::DynamicSlot {
            if is_ponded {
                is_surcharged = false;
            } else if is_storage {
                is_surcharged = t.sur_depth > 0.0 && y_last > full_depth;
            } else {
                is_surcharged = y_crown > 0.0 && y_last > y_crown;
            }
        }
        self.xnode_.is_surcharged[ui] = if is_surcharged { 1 } else { 0 };

        // Virtual junction: zero-storage node on the natural half-link area
        if t.is_virtual != 0 {
            let sum = self.xnode_.sumdqdh[ui];
            let denom_semi = surf_area + 0.5 * dt * sum;
            let semi_implicit = self.node_continuity == NodeContinuity::SemiImplicit && denom_semi > FUDGE;
            let mut y_vj;
            if semi_implicit {
                y_vj = y_old + d_v / denom_semi;
                self.xnode_.old_surf_area[ui] = surf_area;
                if step > 0 {
                    y_vj = (1.0 - self.omega) * y_last + self.omega * y_vj;
                }
            } else if surf_area > FUDGE && !is_surcharged {
                y_vj = y_old + d_v / surf_area;
                self.xnode_.old_surf_area[ui] = surf_area;
                if step > 0 {
                    y_vj = (1.0 - self.omega) * y_last + self.omega * y_vj;
                }
            } else if sum > FUDGE {
                let mut denomv = sum;
                if y_crown > 0.0 && y_last < 1.25 * y_crown {
                    let f = (y_last - y_crown) / y_crown;
                    denomv += (self.xnode_.old_surf_area[ui] / dt - sum) * (-15.0 * f).exp(); // LIBM exp
                }
                y_vj = if denomv != 0.0 { y_last + d_q / denomv } else { y_last };
                if is_surcharged && y_vj < y_crown {
                    y_vj = y_crown - FUDGE;
                }
            } else {
                y_vj = y_last;
            }
            self.commit_node_depth_state(ctx, node_idx, y_vj, d_v, dt);
            return;
        }

        let use_surcharge_dqdh = is_surcharged && self.surcharge_method == SurchargeMethod::Extran;

        let mut y_new;

        if self.node_continuity == NodeContinuity::SemiImplicit {
            let mut denom = surf_area + 0.5 * dt * self.xnode_.sumdqdh[ui];
            if ctx.nodes.storage_area[ui] == 0.0 {
                denom = cmax(denom, self.min_surf_area_);
            }
            let dy = d_v / denom;
            y_new = y_old + dy;
            if !is_ponded {
                self.xnode_.old_surf_area[ui] = surf_area;
            }
            if step > 0 {
                y_new = (1.0 - self.omega) * y_last + self.omega * y_new;
            }
            if is_ponded && y_new < full_depth {
                y_new = full_depth - FUDGE;
            }
        } else if !use_surcharge_dqdh || (is_storage && self.xnode_.sumdqdh[ui] == 0.0) {
            let dy = d_v / surf_area;
            y_new = y_old + dy;
            if !is_ponded {
                self.xnode_.old_surf_area[ui] = surf_area;
            }
            if step > 0 {
                y_new = (1.0 - self.omega) * y_last + self.omega * y_new;
            }
            if is_ponded && y_new < full_depth {
                y_new = full_depth - FUDGE;
            }
        } else {
            let mut corr = 1.0;
            if t.degree < 0 {
                corr = 0.6;
            }
            let mut denom = self.xnode_.sumdqdh[ui];
            if y_last < 1.25 * y_crown && y_crown > 0.0 {
                let f = (y_last - y_crown) / y_crown;
                denom += (self.xnode_.old_surf_area[ui] / dt - self.xnode_.sumdqdh[ui]) * (-15.0 * f).exp(); // LIBM exp
            }
            let mut dy = 0.0;
            if denom != 0.0 {
                dy = corr * d_q / denom;
            }
            y_new = y_last + dy;
            if y_new < y_crown {
                y_new = y_crown - FUDGE;
            }
            if can_pond && y_new > full_depth {
                y_new = full_depth + FUDGE;
            }
        }

        self.commit_node_depth_state(ctx, node_idx, y_new, d_v, dt);
    }

    /// commitNodeDepthState: the canonical commit of an accepted node depth.
    fn commit_node_depth_state(&mut self, ctx: &mut SimulationContext, node_idx: i32, mut y_new: f64, d_v: f64, dt: f64) {
        let ui = node_idx as usize;
        let t = self.node_tile_[ui];

        y_new = cmax(y_new, 0.0);

        let is_inlet_node = ui < ctx.nodes.is_inlet.len() && ctx.nodes.is_inlet[ui] != 0;

        if t.is_virtual != 0 && !is_inlet_node {
            ctx.nodes.overflow[ui] = 0.0;
            ctx.nodes.volume[ui] = 0.0;
            if dt > 0.0 {
                self.xnode_.dYdT[ui] = (y_new - ctx.nodes.old_depth[ui]).abs() / dt;
            }
            ctx.nodes.depth[ui] = y_new;
            ctx.nodes.head[ui] = t.invert_elev + y_new;
            return;
        }

        let is_coupled = ui < ctx.coupled_node.len() && ctx.coupled_node[ui] != 0;
        let can_pond = (ctx.options.allow_ponding || is_coupled) && (t.ponded_area > 0.0);

        let mut y_max = t.full_depth;
        if is_inlet_node && ctx.nodes.rim_depth[ui] > 0.0 {
            y_max = cmax(y_max, ctx.nodes.rim_depth[ui]);
        }
        if !can_pond {
            y_max += t.sur_depth;
        }

        ctx.nodes.overflow[ui] = 0.0;
        if y_new > y_max {
            if !can_pond {
                ctx.nodes.overflow[ui] = d_v / dt;
                ctx.nodes.volume[ui] = t.rpt_full_volume;
                y_new = y_max;
            } else {
                ctx.nodes.volume[ui] = cmax(ctx.nodes.old_volume[ui] + d_v, t.rpt_full_volume);
                ctx.nodes.overflow[ui] = (ctx.nodes.volume[ui] - cmax(ctx.nodes.old_volume[ui], t.rpt_full_volume)) / dt;
            }
            if ctx.nodes.overflow[ui] < FUDGE {
                ctx.nodes.overflow[ui] = 0.0;
            }
        } else if t.is_storage != 0 {
            let v = node::get_volume(&ctx.nodes, node_idx, y_new, Some(&ctx.tables), self.unit_sys_, Some(&ctx.node_subtypes));
            ctx.nodes.volume[ui] = v;
        } else {
            ctx.nodes.volume[ui] = if t.full_depth > 0.0 { t.rpt_full_volume * (y_new / t.full_depth) } else { 0.0 };
        }

        if is_inlet_node {
            ctx.nodes.volume[ui] = 0.0;
        }

        if dt > 0.0 {
            self.xnode_.dYdT[ui] = (y_new - ctx.nodes.old_depth[ui]).abs() / dt;
        }

        ctx.nodes.depth[ui] = y_new;
        ctx.nodes.head[ui] = t.invert_elev + y_new;
    }

    /// getRoutingStep: the CFL-based variable step.
    pub fn get_routing_step(&mut self, ctx: &mut SimulationContext, fixed_step: f64, courant_factor: f64) -> f64 {
        if courant_factor <= 0.0 {
            return fixed_step;
        }
        if fixed_step < MIN_TIMESTEP {
            return fixed_step;
        }
        let min_route_step = cmax(cmin(ctx.options.min_routing_step, fixed_step), MIN_TIMESTEP);

        if self.variable_step_ <= 0.0 {
            self.variable_step_ = min_route_step;
            return self.variable_step_;
        }

        let mut dt_min = fixed_step;
        let mut min_link = -1;
        let mut min_node = -1;

        for j in 0..self.n_links_ {
            let mut t = self.get_link_step(ctx, j);
            if t <= 0.0 || t > 1.0e9 {
                continue;
            }
            t *= courant_factor;
            if t < dt_min {
                dt_min = t;
                min_link = j;
                min_node = -1;
            }
        }

        for i in 0..self.n_nodes_ {
            let ui = i as usize;
            if ctx.nodes.type_[ui] == NodeType::Outfall {
                continue;
            }
            if ctx.nodes.depth[ui] <= FUDGE {
                continue;
            }
            let y_crown = ctx.nodes.crown_elev[ui] - ctx.nodes.invert_elev[ui];
            if ctx.nodes.depth[ui] + FUDGE >= y_crown {
                continue;
            }
            let max_depth = y_crown * 0.25;
            if max_depth < FUDGE {
                continue;
            }
            let d_ydt = self.xnode_.dYdT[ui];
            if d_ydt < FUDGE {
                continue;
            }
            let t = max_depth / d_ydt;
            if t > 0.0 && t < dt_min {
                dt_min = t;
                min_node = i;
                min_link = -1;
            }
        }

        // Virtual-junction pairs (a safety net; none unless the model has them)
        for p in self.vjunc_.iter() {
            if p.through == 0 || p.lambda <= 0.0 {
                continue;
            }
            let uu = p.up_link as usize;
            let aj = self.area2_[uu];
            let wj = self.width2_[uu];
            if aj <= FUDGE || wj <= FUDGE {
                continue;
            }
            let ucu = self.tile_uj_to_ci_[uu] as usize;
            let q_u = ctx.links.flow[uu] / self.tile_barrels_d_[ucu];
            let vmag = q_u.abs() / aj;
            let c = (constants::GRAVITY * (aj / wj)).sqrt();
            if vmag + c <= FUDGE {
                continue;
            }
            let t = p.lambda / (vmag + c) * courant_factor;
            if t > 0.0 && t < dt_min {
                dt_min = t;
                min_node = p.node;
                min_link = -1;
            }
        }

        if min_node >= 0 {
            ctx.nodes.stat_time_courant_critical[min_node as usize] += 1.0;
        } else if min_link >= 0 {
            ctx.links.stat_time_courant_critical[min_link as usize] += 1.0;
        }

        let min_step = min_route_step;
        dt_min = cmax(dt_min, min_step);
        dt_min = (1000.0 * dt_min).floor() / 1000.0;
        dt_min
    }

    /// getLinkStep: one conduit's Courant time.
    pub fn get_link_step(&self, ctx: &SimulationContext, link_idx: i32) -> f64 {
        let uj = link_idx as usize;
        if ctx.links.type_[uj] != LinkType::Conduit {
            return 1.0e10;
        }
        let cd = &ctx.link_subtypes.conduits;
        let ucr = ctx.link_subtypes.conduit_row(link_idx) as usize;

        let barrels = cmax(cd.barrels[ucr], 1);
        let q = ctx.links.flow[uj].abs() / (barrels as f64);
        if q <= FUDGE {
            return 1.0e10;
        }
        let a = self.area_mid_[uj];
        if a <= FUDGE {
            return 1.0e10;
        }
        let mut l = cd.true_length[ucr];
        if l <= 0.0 {
            l = cd.length[ucr];
        }
        let mod_l = cd.mod_length[ucr];

        // (DYNAMIC_SLOT celerity path: that method is refused in execute; its dps_ arrays are
        // never surcharged here)
        if self.surcharge_method == SurchargeMethod::DynamicSlot && (link_idx as usize) < self.tile_uj_to_ci_.len() {
            let ci = self.tile_uj_to_ci_[uj];
            if ci >= 0 {
                let uci = ci as usize;
                if self.dps_.surcharged[uci] != 0 && l > 0.0 {
                    let lscale = if l > 0.0 && mod_l > 0.0 { mod_l / l } else { 1.0 };
                    let v = q / convey_area(a, ctx.links.xsect_a_full[uj]);
                    let p = cmax(self.dps_.P[uci], 1.0);
                    let c_p = self.dps_config_.c_pT / p;
                    let denom = v.abs() + c_p;
                    if denom > 0.0 {
                        return (l * lscale) / denom;
                    }
                }
            }
        }

        let fr = self.froude_[uj];
        if fr <= 0.01 {
            return 1.0e10;
        }
        let vol = ctx.links.volume[uj] / (barrels as f64);
        let mut t = vol / q;
        if l > 0.0 && mod_l > 0.0 {
            t = t * mod_l / l;
        }
        t = t * fr / (1.0 + fr);
        t
    }

    /// accumulateLinkStepStats ([REPORT] LINK_STEPS).
    pub fn accumulate_link_step_stats(&self, ctx: &mut SimulationContext, dt: f64, fixed_step: f64, courant_factor: f64) {
        if dt <= 0.0 {
            return;
        }
        let cf = if courant_factor > 0.0 { courant_factor } else { 1.0 };
        for j in 0..self.n_links_ {
            let uj = j as usize;
            if ctx.links.type_[uj] != LinkType::Conduit {
                continue;
            }
            let n1 = ctx.links.node1[uj] as usize;
            let n2 = ctx.links.node2[uj] as usize;
            let c1 = self.node_tile_[n1].is_outfall != 0 || self.xnode_.converged[n1] != 0;
            let c2 = self.node_tile_[n2].is_outfall != 0 || self.xnode_.converged[n2] != 0;
            ctx.links.stat_conv_total[uj] += 1;
            if c1 && c2 {
                ctx.links.stat_conv_steps[uj] += 1;
            }
            let t = self.get_link_step(ctx, j);
            if t <= 0.0 || t > 1.0e9 {
                continue;
            }
            ctx.links.add_lstep_sample(j, cmin(t * cf, fixed_step), dt, fixed_step);
        }
    }
}

/// node::getPondedArea for a storage node (no LID cells: refused upstream).
fn lidnode_free_ponded_area(ctx: &SimulationContext, idx: i32, depth: f64, unit_sys: i32) -> f64 {
    crate::hydraulics::node_step::get_ponded_area(&ctx.nodes, idx, depth, Some(&ctx.tables), unit_sys, Some(&ctx.node_subtypes))
}
