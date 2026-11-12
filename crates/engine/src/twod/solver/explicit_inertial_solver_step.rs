// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/solver/ExplicitInertialSolver.cpp (advance,
// advanceImpl, syncAndRebuild, refreshDt0, the face and cell firings, the boundary edges,
// settleAccumulators, lazySourcesOnly, the held cell sources, computeLimitedGradientsSwe,
// boundaryFluxSwe, runMacroCycle, runRk2Step, runMomentumStep, applyMomentumPhysics,
// refreshConveyance, refreshDiffusiveSlopes, flushSourceLedgers, finalize, run_stats)
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin,
// licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The step of the serial CPU marcher (phase 2a wave 4a): the LTS macro cycle with its binary
// tiers, the global-step tail, the lazy source tier, the active-set / tier rebuild with the
// FRONT_REBUILD halo, the face laws of LOCAL_INERTIAL, DIFFUSIVE_WAVE and FULL_SWE (first and
// second order), positivity, wet/dry, the boundary edges of every type, the surface-only
// RK2 (RECONSTRUCTION_ORDER 2) and the optional momentum physics (Coriolis, eddy viscosity,
// wall drag), and the run statistics.
//
// What phase 2a leaves out, and why the omission cannot change a result:
//   * The phase-2b owners are never attached here: gw_ (the two-zone kernel), common_stage_,
//     channel_ground_, the source-result recorders, and the RK2 reaction binding are refused by
//     SurfaceRouter2D::initialize; surface species (transport.n_species > 0), mesh heat, mesh
//     snow, the runtime sources and forcings (the swmm_2d_* C API) and the RK2 rain receipts
//     are refused or unreachable too. advance() checks that every one of them is absent before
//     anything changes (stub otherwise), and the code below leaves out the branches they
//     select, keeping each remaining expression as upstream writes it. rk_transport_stage_ is
//     set only with gw_, so it stays false, and with it boundary_stage and preserve_face_water.
//   * The live junction exchange in the tier-0 cell firing (exch_ non-empty: a mesh coupled to
//     1D nodes) was added at the join (12 November 2026, phase 2a wave 4c): computeNodeCouplingQ
//     is twod/coupling/node_coupling.rs (wave 4b), and StepAttachments carries what upstream
//     reaches through state_->node_coupling, state_->nodes_1d, state_->coupling_head_cells and
//     state_->node_pond_area_m2. Its species bookings belong to surface species (phase 2b).
//   * OpenMP. Upstream's loops are `schedule(static)` over `num_threads(opts->num_threads)`.
//     Each iteration writes only its own cell or face, except the per-thread source ledgers
//     (bookSourceLedger / the evaporation slot, folded thread by thread in flushSourceLedgers)
//     and the per-thread partial minima of the step bound (rebuild_dt_partial_). Those two keep
//     upstream's thread partition: libgomp's static schedule gives thread t the contiguous
//     block omp_static_thread() computes, and the serial loop here books each iteration into
//     the slot of the thread upstream would run it on. The team is THREADS capped by
//     OMP_THREAD_LIMIT, as libgomp sizes `num_threads(n)` (core::thread_info::omp_team_size;
//     12 November 2026, wave 5: the CLI runs under OMP_THREAD_LIMIT=1 like the reference's, and
//     a THREADS 2 deck then books every iteration into slot 0 there). compactOrdered (parallel
//     list compaction above 65,536 entries) produces the serial order by construction, so the
//     serial loops stand for it.
//   * The development switches OPENSWMM_2D_MUSCL_OFF and OPENSWMM_2D_MARCHER_CHECK are refused
//     with the crate's stub when set (PLAN.md, the environment-switch policy).
//   * Exceptions upstream throws (std::runtime_error, std::invalid_argument, std::out_of_range
//     from map::at) propagate through SurfaceRouter2D to SWMMEngine::step's catch; here they are
//     Err(Terminated) with exception "throw <type>" (core::engine_step::THROWN_PREFIX).
//   * std::ilogb is ilogb() below (the exponent of a finite positive double, INT_MAX for
//     infinity); std::hypot is f64::hypot (glibc), std::cbrt is core::cmath::glibc_cbrt,
//     std::cos / std::sin are f64::cos / f64::sin (glibc; MomentumPhysics2D::rotate).

#![allow(clippy::all, non_snake_case)]

use crate::core::cmath::{cmax, cmin, cvt_i32, glibc_cbrt};
use crate::core::engine_init_hydraulics::Terminated;
use crate::core::engine_step::stub;
use crate::data::node::NodeData;
use crate::twod::coupling::node_coupling::{compute_node_coupling_q, CouplingAttachments, CouplingPoint};
use crate::twod::data::boundary_data::{BoundaryData, BoundaryType};
use crate::twod::data::mesh_data::MeshData;
use crate::twod::data::solver_options_2d::{CellClosure2D, FaceDepth2D, Momentum2D, SolverOptions2D, Turbulence2D};
use crate::twod::data::surface_state_data::{PondBase, SurfaceStateData};
use crate::twod::solver::diffusive_kernels::{cell_diffusive_dt, face_discharge};
use crate::twod::solver::explicit_inertial_solver::{CellSourceStep, ExplicitInertialSolver, SourceLedger, K_REBUILD_EVERY_CYCLES};
use crate::twod::solver::inertial_kernels::{
    cell_cfl_dt, cell_eta_depth, cell_volume_from_eta, face_depth_from_eta, face_flow_depth, face_flow_depth_vfr, froude_cap, inertial_advection,
    inertial_face_update, momentum_blend_fraction, q_magnitude, K_ETA_DEADBAND, K_GRAVITY,
};
use crate::twod::solver::momentum_physics_2d::MomentumPhysics2D;
use crate::twod::solver::surface_flux_calculator::{compute_boundary_edge_flux, evap_sink, infil_sink};
use crate::twod::solver::swe_kernels::{
    bj_limiter, critical_depth, depth_from_invariant_and_discharge, face_flux, face_flux_recon_bed, friction_update, FaceFlux,
};
use crate::twod::solver::swe_kernels::K_GRAVITY as SWE_GRAVITY;

/// Halo depth (rings) when FRONT_REBUILD is on (kFrontHaloRings).
const K_FRONT_HALO_RINGS: i32 = 5;

/// What the step reads of the router and the 1D network besides the surface state, which
/// upstream's solver reaches through its pointers mesh_, opts_ and state_->boundary,
/// state_->node_coupling, state_->nodes_1d, state_->coupling_head_cells and
/// state_->node_pond_area_m2.
pub struct StepAttachments<'a> {
    /// *mesh_ (the router's mesh_)
    pub mesh: &'a MeshData,
    /// *opts_ (the router's options_)
    pub opts: &'a SolverOptions2D,
    /// state_->boundary (the router's boundary_), when attached
    pub boundary: Option<&'a BoundaryData>,
    /// *state_->node_coupling (the router's node_coupling_points_), when attached
    pub node_coupling: Option<&'a [CouplingPoint]>,
    /// *state_->nodes_1d (the context's node store, the 1D heads and volumes the batch holds),
    /// when attached
    pub nodes_1d: Option<&'a NodeData>,
    /// The targets of state_->coupling_head_cells and state_->node_pond_area_m2 (the router's
    /// coupling_head_cells_ and node_pond_area_m2_), as the exchange law reads them
    pub coupling: CouplingAttachments<'a>,
}

/// A C++ exception on its way to SWMMEngine::step's catch.
pub fn thrown(kind: &'static str, what: String) -> Terminated {
    Terminated { exception: kind, what: Box::leak(what.into_boxed_str()) }
}

fn runtime_error(what: &str) -> Terminated {
    thrown("throw std::runtime_error", what.to_string())
}

/// libgomp's `schedule(static)` without a chunk size: the thread of a team of `t` that runs
/// iteration `k` of `n` (gomp_iter_static_next: q = n / t, r = n % t; the first r threads take
/// q + 1 iterations each, the rest q, in order).
#[inline]
pub fn omp_static_thread(k: usize, n: usize, t: usize) -> usize {
    if t <= 1 {
        return 0;
    }
    let q = n / t;
    let r = n % t;
    let big = (q + 1) * r;
    if k < big {
        k / (q + 1)
    } else {
        r + (k - big) / q
    }
}

/// std::ilogb for the ratios the tiering takes (finite >= 2 or +infinity): the unbiased binary
/// exponent, INT_MAX for an infinity (glibc's FP_ILOGB of a NaN is never asked for here).
#[inline]
fn ilogb(x: f64) -> i32 {
    if x.is_infinite() {
        return i32::MAX;
    }
    if x.is_nan() || x == 0.0 {
        return i32::MIN;
    }
    let bits = x.to_bits();
    let e = ((bits >> 52) & 0x7ff) as i32;
    if e == 0 {
        // subnormal: normalize
        let m = bits & 0x000f_ffff_ffff_ffff;
        return -1023 - (m.leading_zeros() as i32 - 12);
    }
    e - 1023
}

/// The face lists fireFaces / fireCells take, by owner (upstream passes a reference to one of
/// the solver's own vectors).
#[derive(Clone, Copy)]
enum List {
    Active,
    Tier(usize),
}

impl ExplicitInertialSolver {
    fn take_faces(&mut self, l: List) -> Vec<i32> {
        match l {
            List::Active => std::mem::take(&mut self.active_faces_),
            List::Tier(k) => std::mem::take(&mut self.edges_by_tier_[k]),
        }
    }
    fn put_faces(&mut self, l: List, v: Vec<i32>) {
        match l {
            List::Active => self.active_faces_ = v,
            List::Tier(k) => self.edges_by_tier_[k] = v,
        }
    }
    fn take_cells(&mut self, l: List) -> Vec<i32> {
        match l {
            List::Active => std::mem::take(&mut self.active_cells_),
            List::Tier(k) => std::mem::take(&mut self.cells_by_tier_[k]),
        }
    }
    fn put_cells(&mut self, l: List, v: Vec<i32>) {
        match l {
            List::Active => self.active_cells_ = v,
            List::Tier(k) => self.cells_by_tier_[k] = v,
        }
    }

    /// The team `num_threads(max(1, opts->num_threads))` gets: OMP_THREAD_LIMIT caps it.
    fn nthreads(opts: &SolverOptions2D) -> usize {
        cmax(1, crate::core::thread_info::omp_team_size(cmax(1, opts.num_threads))) as usize
    }

    /// The phase-2b owners and the runtime C API stores this translation leaves out (see the
    /// header): Some(what) when one of them is present.
    pub fn step_not_translated(&self, state: &SurfaceStateData) -> Option<&'static str> {
        if self.gw_ {
            return Some("the 2D marcher with the two-zone groundwater kernel attached (phase 2b)");
        }
        if self.common_stage_ || self.channel_ground_ || self.source_history_ || self.engine_source_history_ || self.engine_source_stream_ {
            return Some("the 2D marcher with a common source stage, channel ground or a source-result recorder (phase 2b)");
        }
        if self.rk_reactions_ || state.transport.n_species > 0 || self.species_on_ {
            return Some("2D surface species transport in the marcher (phase 2b)");
        }
        if state.heat.enabled {
            return Some("2D mesh heat in the marcher (phase 2b)");
        }
        if state.snow.active() || !state.snow.pack_.is_empty() {
            return Some("2D mesh snow in the marcher (phase 2b)");
        }
        if !state.runtime_sources.rows.is_empty() || !state.runtime_sources.cells.is_empty() || !state.runtime_forcings.rows.is_empty() {
            return Some("the 2D runtime sources and forcings (the swmm_2d_* C API)");
        }
        if !state.rk2_rain_pending.is_empty() {
            return Some("the 2D RK2 rain receipts (coupled groundwater RK2; phase 2b)");
        }
        if !self.channel_source_ready_.is_empty() || !self.rk_boundary_signals_.is_empty() {
            return Some("the 2D channel-ground sources or RK2 boundary signals (phase 2b)");
        }
        None
    }

    // -----------------------------------------------------------------------------------
    // advance / advanceImpl
    // -----------------------------------------------------------------------------------

    /// advance(t_current, t_target): the coupled-RK2 rollback wrapper (gw_ is never attached in
    /// phase 2a, so it is advanceImpl).
    pub fn advance(&mut self, a: &StepAttachments, st: &mut SurfaceStateData, t_current: f64, t_target: f64) -> Result<f64, Terminated> {
        if let Some(what) = self.step_not_translated(st) {
            return Err(stub(what));
        }
        if std::env::var_os("OPENSWMM_2D_MUSCL_OFF").is_some() {
            return Err(stub("the OPENSWMM_2D_MUSCL_OFF development switch (ExplicitInertialSolver::fireFacesSwe)"));
        }
        if std::env::var_os("OPENSWMM_2D_MARCHER_CHECK").is_some() {
            return Err(stub("the OPENSWMM_2D_MARCHER_CHECK development switch (ExplicitInertialSolver::runMacroCycle)"));
        }
        self.advance_impl(a, st, t_current, t_target)
    }

    /// advanceImpl(t_current, t_target).
    fn advance_impl(&mut self, a: &StepAttachments, st: &mut SurfaceStateData, t_current: f64, t_target: f64) -> Result<f64, Terminated> {
        if !self.initialized_ || t_target <= t_current {
            return Ok(t_target);
        }
        self.refresh_conveyance(a.mesh);
        self.plain_sources_ = !self.second_order_
            && !self.gw_
            && !self.common_stage_
            && !self.channel_ground_
            && st.transport.n_species == 0
            && !st.heat.enabled
            && !st.snow.active()
            && st.runtime_forcings.rows.is_empty()
            && st.runtime_sources.rows.is_empty()
            && st.rk2_rain_pending.is_empty();
        let r = self.advance_body(a, st, t_current, t_target);
        // ResetPlainSources: on every exit
        self.plain_sources_ = false;
        r
    }

    fn advance_body(&mut self, a: &StepAttachments, st: &mut SurfaceStateData, t_current: f64, t_target: f64) -> Result<f64, Terminated> {
        let opts = a.opts;
        let mut t = t_current;
        if self.t_last_sync_ > t_current {
            self.t_last_sync_ = t_current;
        }
        for x in self.bc_accum_.iter_mut() {
            *x = 0.0;
        }
        for x in self.exch_.iter_mut() {
            *x = 0.0;
        }
        for x in self.node_drawn_.iter_mut() {
            *x = 0.0;
        }
        for x in st.transport.exch_mass.iter_mut() {
            *x = 0.0;
        }
        for x in st.transport.exch_spill.iter_mut() {
            *x = 0.0;
        }
        self.exch_tau_ = 0.0;
        self.last_steps_ = 0;
        let mut cycles_since_rebuild = self.cycles_since_rebuild_;

        while t < t_target {
            if cycles_since_rebuild >= K_REBUILD_EVERY_CYCLES || (self.front_breach_ && cycles_since_rebuild > 0) {
                self.sync_and_rebuild(a, st, t)?;
                cycles_since_rebuild = 0;
            } else {
                self.refresh_dt0(a, st);
            }
            if opts.momentum_physics() {
                let dt = self.run_momentum_step(a, st, cmin(self.dt0_, t_target - t), t)?;
                if !(dt > 0.0) || t + dt == t {
                    return Err(runtime_error("Optional momentum physics cannot advance the clock."));
                }
                t = if dt == t_target - t { t_target } else { t + dt };
                self.last_dt_ = dt;
                cycles_since_rebuild += 1;
                continue;
            }
            let k_tiers = self.cells_by_tier_.len() as i32;
            let nsub_full = 1i32 << (k_tiers - 1);
            let target = t_target;
            let remaining = target - t;

            if self.active_cells_.is_empty() {
                // Quiescent: stride the window; the lazy tier keeps accumulating.
                t = target;
                self.last_dt_ = remaining;
                cycles_since_rebuild = K_REBUILD_EVERY_CYCLES;
                if t < t_target {
                    continue;
                }
                break;
            }

            if nsub_full as f64 * self.dt0_ <= remaining {
                let dt0 = self.dt0_;
                if self.second_order_ {
                    self.run_rk2_step(a, st, dt0, t)?;
                } else {
                    self.run_macro_cycle(a, st, dt0, nsub_full, t)?;
                }
                t += nsub_full as f64 * self.dt0_;
                self.last_dt_ = self.dt0_;
                cycles_since_rebuild += 1;
                continue;
            }

            // Tail: split the whole remaining span into nt equal global substeps.
            {
                // (the sub-ulp landing: no groundwater clock, no common stage, no coupled RK2)
                if remaining <= self.dt0_ * 1.0e-9 {
                    t = target;
                    if t < t_target {
                        continue;
                    }
                    break;
                }
                self.settle_accumulators(a, st);
                let nt = cmax(1, cvt_i32((remaining / self.dt0_).ceil()));
                let dt_tail = remaining / nt as f64;
                for s in 0..nt {
                    if self.second_order_ {
                        self.run_rk2_step(a, st, dt_tail, t + s as f64 * dt_tail)?;
                        continue;
                    }
                    self.fire_faces(a, st, List::Active, dt_tail, true)?;
                    self.source_phase_end_ = if s + 1 == nt { target } else { t + (s + 1) as f64 * dt_tail };
                    self.fire_cells(a, st, List::Active, dt_tail, true)?;
                    // captureSourceResults(): no recorder in phase 2a
                    self.substeps_run_ += 1;
                    self.last_steps_ += 1;
                }
                self.accumulators_pending_ = false; // every side just gathered
                self.last_dt_ = dt_tail;
                cycles_since_rebuild += 1;
                t = target; // exact landing
                if t < t_target {
                    continue;
                }
                break;
            }
        }

        self.cycles_since_rebuild_ = cycles_since_rebuild;
        self.settle_accumulators(a, st);
        // Final lazy-source landing.
        if t_target > self.t_last_sync_ {
            if self.cycles_since_rebuild_ >= K_REBUILD_EVERY_CYCLES {
                self.sync_and_rebuild(a, st, t_target)?;
                self.cycles_since_rebuild_ = 0;
            } else {
                self.lazy_sources_only(a, st, t_target);
            }
        }

        // Publish the flux picture the router's output/ledger contract reads.
        for x in st.edge_flux.iter_mut() {
            *x = 0.0;
        }
        let vfr_face = opts.face_reconstruction == FaceDepth2D::VfrFace;
        let ed = &self.edges_;
        for &e in &self.active_faces_ {
            let ue = e as usize;
            let mut qp = 0.0;
            if self.mode_ == Momentum2D::LocalInertial {
                let hf = if vfr_face {
                    face_flow_depth_vfr(st.head[ed.cL[ue] as usize], st.head[ed.cR[ue] as usize], ed.ze_lo[ue], ed.ze_hi[ue])
                } else {
                    face_flow_depth(st.head[ed.cL[ue] as usize], st.head[ed.cR[ue] as usize], ed.zface[ue])
                };
                if hf > opts.dry_depth {
                    qp = froude_cap(self.q_[ue], hf, opts.froude_max) * ed.conveyance[ue];
                }
            } else {
                qp = self.q_[ue];
            }
            let f = qp * ed.xi[ue];
            st.edge_flux[ed.slotL[ue] as usize] = -f;
            st.edge_flux[ed.slotR[ue] as usize] = f;
        }
        let span = t_target - t_current;
        if span > 0.0 {
            for k in 0..self.bc_cell_.len() {
                st.edge_flux[self.bc_slot_[k] as usize] = self.bc_accum_[k] / span;
            }
        }
        // refreshSourceResultBoundary(): no recorder in phase 2a
        self.checkpoint_ready_ = true;
        Ok(t_target)
    }

    // -----------------------------------------------------------------------------------
    // Sources
    // -----------------------------------------------------------------------------------

    /// flushSourceLedgers(): fold every thread's ledger slots into the transport ledgers and
    /// the evaporation total, thread by thread.
    pub fn flush_source_ledgers(&mut self, st: &mut SurfaceStateData) {
        let ns = st.transport.n_species as usize;
        let count = SourceLedger::Count as usize;
        let stride = self.source_ledger_stride_;
        if stride == 0 {
            return;
        }
        let mut base = 0usize;
        while base < self.source_ledgers_.len() {
            for k in 0..count {
                for r in 0..ns {
                    let delta = self.source_ledgers_[base + k * ns + r];
                    let tr = &mut st.transport;
                    let ledger: &mut Vec<f64> = match k {
                        0 => &mut tr.lost_infiltration,
                        1 => &mut tr.lost_coupling,
                        2 => &mut tr.gained_rainfall,
                        3 => &mut tr.gained_coupling,
                        4 => &mut tr.gained_exfiltration,
                        5 => &mut tr.lost_boundary,
                        6 => &mut tr.gained_external,
                        7 => &mut tr.lost_external,
                        _ => &mut tr.gained_bed,
                    };
                    if !ledger.is_empty() {
                        ledger[r] += delta;
                    }
                    self.source_ledgers_[base + k * ns + r] = 0.0;
                }
            }
            let ev = base + count * ns;
            st.evap_loss_total += self.source_ledgers_[ev];
            self.source_ledgers_[ev] = 0.0;
            base += stride;
        }
    }

    /// applyPlainCellSources(i, dt): the water-only held sources (rain, held coupling,
    /// infiltration, evaporation), shared proportionally against the cell's water. `tid` is
    /// the OpenMP thread upstream runs the call on (its ledger slot).
    fn apply_plain_cell_sources(&mut self, a: &StepAttachments, st: &mut SurfaceStateData, i: i32, dt: f64, tid: usize) -> bool {
        let ui = i as usize;
        let opts = a.opts;
        let area = a.mesh.tri_area[ui];
        let area_dt = area * dt;
        let rain = cmax(0.0, st.rainfall[ui]) * area_dt;
        let coupling = st.coupling_flux[ui] * area_dt;
        let requested_infil = cmax(0.0, infil_sink(st.infil_rate[ui], st.depth[ui], opts.dry_depth)) * area_dt;
        let evap_demand = st.evap_rate[ui] * 1.0;
        let requested_evap = cmax(0.0, evap_sink(evap_demand, st.depth[ui], opts.dry_depth)) * area_dt;
        if rain == 0.0 && coupling == 0.0 && requested_infil == 0.0 && requested_evap == 0.0 {
            return false;
        }
        let incoming = rain + cmax(coupling, 0.0) + 0.0;
        let available = cmax(0.0, st.volume[ui] + incoming);
        let requested_out = cmax(-coupling, 0.0);
        let requested = requested_infil + requested_evap + requested_out;
        let scale = if requested > available && requested > 0.0 { available / requested } else { 1.0 };
        let infil = requested_infil * scale;
        let evap = requested_evap * scale;
        let out = requested_out * scale;
        st.infil_applied[ui] += infil / area;
        st.coupling_applied[ui] += cmax(coupling, 0.0) - out;
        self.source_ledgers_[tid * self.source_ledger_stride_] += evap;
        st.volume[ui] = if requested >= available { 0.0 } else { available - requested };
        true
    }

    /// prepareCellSources(i, dt, source) for phase 2a's owners (see the header): None when the
    /// held tuple is empty (upstream's false).
    fn prepare_cell_sources(&mut self, a: &StepAttachments, st: &mut SurfaceStateData, i: i32, dt: f64, tid: usize) -> Option<CellSourceStep> {
        let ui = i as usize;
        let opts = a.opts;
        let area = a.mesh.tri_area[ui];
        let area_dt = area * dt;
        let rain = cmax(0.0, st.liquid_rainfall(i)) * area_dt;
        let coupling = st.coupling_flux[ui] * area_dt;
        let requested_infil = cmax(0.0, infil_sink(st.runtime_forcings.rate(3, i, st.infil_rate[ui]), st.depth[ui], opts.dry_depth)) * area_dt;
        let evap_demand = st.exposed_evaporation(i);
        let requested_evap = cmax(0.0, evap_sink(evap_demand, st.depth[ui], opts.dry_depth)) * area_dt;
        let back = 0.0;
        if rain == 0.0 && coupling == 0.0 && requested_infil == 0.0 && requested_evap == 0.0 && back == 0.0 && !st.heat.enabled {
            return None;
        }
        let incoming = rain + cmax(coupling, 0.0) + back;
        let available = cmax(0.0, st.volume[ui] + incoming);
        let requested_out = cmax(-coupling, 0.0);
        let requested = requested_infil + requested_evap + requested_out;
        let donor_water = available;
        let scale = if requested > donor_water && requested > 0.0 { donor_water / requested } else { 1.0 };
        let infil = requested_infil * scale;
        let evap = requested_evap * scale;
        let out = requested_out * scale;
        st.infil_applied[ui] += infil / area;
        st.coupling_applied[ui] += cmax(coupling, 0.0) - out;
        let ns = st.transport.n_species as usize;
        self.source_ledgers_[tid * self.source_ledger_stride_ + SourceLedger::Count as usize * ns] += evap;
        Some(CellSourceStep {
            area_dt,
            rain,
            coupling,
            available,
            infil,
            evap,
            out,
            final_volume: if requested >= available { 0.0 } else { available - requested },
            requested_infil,
            ledger: true,
        })
    }

    /// applyRuntimeSources(i, dt): the runtime C API sources (none in phase 2a; advance()
    /// refuses a populated store).
    fn apply_runtime_sources(&mut self, st: &SurfaceStateData, i: i32) -> bool {
        let store = &st.runtime_sources;
        if i as usize >= store.cells.len() || store.cells[i as usize].is_empty() {
            return false;
        }
        false
    }

    /// applyCellSources(i, dt, boundary_donor).
    fn apply_cell_sources(&mut self, a: &StepAttachments, st: &mut SurfaceStateData, i: i32, dt: f64, tid: usize) -> bool {
        if self.plain_sources_ {
            return self.apply_plain_cell_sources(a, st, i, dt, tid);
        }
        let source = match self.prepare_cell_sources(a, st, i, dt, tid) {
            None => return self.apply_runtime_sources(st, i),
            Some(s) => s,
        };
        st.volume[i as usize] = source.final_volume;
        self.apply_runtime_sources(st, i);
        true
    }

    /// settleAccumulators(): gather every cell's pending face receipts (and FULL_SWE momentum).
    pub fn settle_accumulators(&mut self, a: &StepAttachments, st: &mut SurfaceStateData) {
        if !self.accumulators_pending_ {
            return;
        }
        self.accumulators_pending_ = false;
        let nt = a.mesh.n_triangles();
        let stride = self.face_accum_stride_;
        let ed = &self.edges_;
        for i in 0..nt {
            let ui = i as usize;
            let mut pending = 0.0;
            for p in ed.cell_ptr[ui] as usize..ed.cell_ptr[ui + 1] as usize {
                let e = ed.cell_edge[p] as usize;
                if ed.cell_sign[p] > 0 {
                    pending += self.face_accum_[stride * e];
                    self.face_accum_[stride * e] = 0.0;
                } else {
                    pending += self.face_accum_[stride * e + 1];
                    self.face_accum_[stride * e + 1] = 0.0;
                }
            }
            if stride == 6 && !self.face_accum_.is_empty() {
                let mut dmx = 0.0;
                let mut dmy = 0.0;
                for p in ed.cell_ptr[ui] as usize..ed.cell_ptr[ui + 1] as usize {
                    let e = ed.cell_edge[p] as usize;
                    if ed.cell_sign[p] > 0 {
                        dmx += self.face_accum_[stride * e + 2];
                        self.face_accum_[stride * e + 2] = 0.0;
                        dmy += self.face_accum_[stride * e + 4];
                        self.face_accum_[stride * e + 4] = 0.0;
                    } else {
                        dmx += self.face_accum_[stride * e + 3];
                        self.face_accum_[stride * e + 3] = 0.0;
                        dmy += self.face_accum_[stride * e + 5];
                        self.face_accum_[stride * e + 5] = 0.0;
                    }
                }
                if dmx != 0.0 || dmy != 0.0 {
                    let area = a.mesh.tri_area[ui];
                    self.qcx_[ui] += dmx / area;
                    self.qcy_[ui] += dmy / area;
                }
            }
            if pending == 0.0 {
                continue;
            }
            let v = st.volume[ui] + pending;
            st.volume[ui] = if v > 0.0 { v } else { 0.0 };
            let (eta, depth) = cell_eta_depth(a.mesh, a.opts, i, st.volume[ui]);
            st.head[ui] = eta;
            st.depth[ui] = depth;
        }
    }

    /// lazySourcesOnly(t): settle the inactive cells' held sources up to t.
    pub fn lazy_sources_only(&mut self, a: &StepAttachments, st: &mut SurfaceStateData, t: f64) {
        let dt_lazy = t - self.t_last_sync_;
        if dt_lazy <= 0.0 {
            return;
        }
        let nt = a.mesh.n_triangles();
        let nthr = Self::nthreads(a.opts);
        for i in 0..nt {
            if self.cell_active_[i as usize] != 0 {
                continue;
            }
            let tid = omp_static_thread(i as usize, nt as usize, nthr);
            if self.apply_cell_sources(a, st, i, dt_lazy, tid) {
                let (eta, depth) = cell_eta_depth(a.mesh, a.opts, i, st.volume[i as usize]);
                st.head[i as usize] = eta;
                st.depth[i as usize] = depth;
            }
        }
        self.flush_source_ledgers(st);
        self.t_last_sync_ = t;
    }

    // -----------------------------------------------------------------------------------
    // The active set and the tiers
    // -----------------------------------------------------------------------------------

    /// syncAndRebuild(t): settle, land the lazy sources, reseed the active set and its halo,
    /// assign the LTS tiers and rebuild the face lists.
    fn sync_and_rebuild(&mut self, a: &StepAttachments, st: &mut SurfaceStateData, t: f64) -> Result<(), Terminated> {
        let mesh = a.mesh;
        let opts = a.opts;
        self.settle_accumulators(a, st);
        let dt_lazy = t - self.t_last_sync_;
        let nt = mesh.n_triangles();
        let unt = nt as usize;
        let nthr = Self::nthreads(opts);

        // 2. Seed.
        let band = cmin(0.001, 0.5 * opts.h_move);
        let h_on = opts.h_move + band;
        self.h_on_front_ = h_on;
        let h_off = cmax(0.0, opts.h_move - band);
        self.rebuild_seed_.clear();
        self.rebuild_seed_.resize(unt, 0u8);
        for i in 0..nt {
            let ui = i as usize;
            if dt_lazy > 0.0 && self.cell_active_[ui] == 0 {
                let tid = omp_static_thread(ui, unt, nthr);
                if self.apply_cell_sources(a, st, i, dt_lazy, tid) {
                    let (eta, depth) = cell_eta_depth(mesh, opts, i, st.volume[ui]);
                    st.head[ui] = eta;
                    st.depth[ui] = depth;
                }
            }
            let thresh = if self.cell_active_[ui] != 0 { h_off } else { h_on };
            if st.depth[ui] >= thresh || st.coupling_flux[ui] != 0.0 || self.pin_t0_[ui] != 0 || st.heat.enabled {
                self.rebuild_seed_[ui] = 1;
            }
        }
        if dt_lazy > 0.0 {
            self.flush_source_ledgers(st);
            self.t_last_sync_ = t;
        }

        // 3. One-ring halo.
        if opts.momentum_physics() {
            for x in self.rebuild_seed_.iter_mut() {
                *x = 1;
            }
        }
        self.cell_active_ = self.rebuild_seed_.clone();
        let n_seed = self.rebuild_seed_.iter().filter(|&&x| x == 1).count() as i32;
        if n_seed != 0 && n_seed != nt {
            let ed = &self.edges_;
            for e in 0..ed.ne as usize {
                let ca = ed.cL[e] as usize;
                let cb = ed.cR[e] as usize;
                if self.rebuild_seed_[ca] != 0 && self.rebuild_seed_[cb] == 0 {
                    self.cell_active_[cb] = 1;
                } else if self.rebuild_seed_[cb] != 0 && self.rebuild_seed_[ca] == 0 {
                    self.cell_active_[ca] = 1;
                }
            }
        }
        if self.front_rebuild_ && n_seed != 0 && n_seed != nt {
            for _ring in 1..K_FRONT_HALO_RINGS {
                self.rebuild_seed_ = self.cell_active_.clone();
                for i in 0..nt {
                    let ui = i as usize;
                    if self.rebuild_seed_[ui] != 0 {
                        continue;
                    }
                    let nvc = mesh.cell_vertex_count(i);
                    for k2 in 0..nvc {
                        let j = mesh.cell_neighbour(i, k2);
                        if j >= 0 && self.rebuild_seed_[j as usize] != 0 {
                            self.cell_active_[ui] = 1;
                            break;
                        }
                    }
                }
            }
        }
        self.active_cells_.clear();
        for i in 0..nt {
            if self.cell_active_[i as usize] != 0 {
                self.active_cells_.push(i);
            }
        }
        if self.front_rebuild_ {
            for i in 0..nt {
                let ui = i as usize;
                let mut f = 0u8;
                if self.cell_active_[ui] != 0 {
                    let nvc = mesh.cell_vertex_count(i);
                    for k2 in 0..nvc {
                        let j = mesh.cell_neighbour(i, k2);
                        if j >= 0 && self.cell_active_[j as usize] == 0 {
                            f = 1;
                            break;
                        }
                    }
                }
                self.frontier_[ui] = f;
            }
            self.front_breach_ = false;
        }

        // 4. Tier assignment from the local CFL step.
        let k_tiers = self.cells_by_tier_.len() as i32;
        let na = self.active_cells_.len();
        if self.mode_ == Momentum2D::DiffusiveWave {
            self.refresh_diffusive_slopes(a, st);
        }
        self.dt0_ = 1.0e30;
        self.rebuild_dt_cell_.resize(na, 0.0);
        self.rebuild_dt_partial_.clear();
        self.rebuild_dt_partial_.resize(cmax(1, opts.num_threads) as usize, 1.0e30); // one per requested thread
        {
            let ed = &self.edges_;
            for k in 0..na {
                let tid = omp_static_thread(k, na, nthr);
                let i = self.active_cells_[k] as usize;
                let h = st.depth[i];
                let mut speed = 0.0;
                if !self.qcx_.is_empty() && h > 1.0e-6 {
                    speed = q_magnitude(self.qcx_[i], self.qcy_[i]) / h;
                }
                let mut dt = if h > opts.dry_depth {
                    cell_cfl_dt(opts.cfl_number, if self.mode_ == Momentum2D::FullSwe { ed.cell_lpos[i] } else { ed.cell_lchar[i] }, h, speed)
                } else {
                    1.0e30
                };
                if self.mode_ == Momentum2D::DiffusiveWave && h > opts.dry_depth {
                    dt = cell_diffusive_dt(opts.cfl_number, ed.cell_lchar[i], h, mesh.mannings_n[i], self.dw_slope_[i], opts.flux_dh_eps / cmax(ed.cell_lchar[i], 1.0e-9));
                }
                dt = cmin(dt, opts.max_timestep);
                self.rebuild_dt_cell_[k] = dt;
                if dt < self.rebuild_dt_partial_[tid] {
                    self.rebuild_dt_partial_[tid] = dt;
                }
            }
        }
        for &v in &self.rebuild_dt_partial_ {
            self.dt0_ = cmin(self.dt0_, v);
        }
        if self.dt0_ >= 1.0e30 {
            self.dt0_ = opts.max_timestep; // fully quiescent
        }

        for v in self.cells_by_tier_.iter_mut() {
            v.clear();
        }
        for v in self.edges_by_tier_.iter_mut() {
            v.clear();
        }
        for k in 0..na {
            let i = self.active_cells_[k] as usize;
            let mut tk = 0i32;
            if k_tiers > 1 {
                let ratio = self.rebuild_dt_cell_[k] / self.dt0_;
                tk = if ratio >= 2.0 { cmin(k_tiers - 1, ilogb(ratio)) } else { 0 };
                if st.coupling_flux[i] != 0.0 || self.pin_t0_[i] != 0 {
                    tk = 0;
                }
            }
            self.tier_[i] = tk as u8;
            if !self.front_rebuild_ {
                self.cells_by_tier_[tk as usize].push(i as i32);
            }
        }
        if self.front_rebuild_ {
            for _ring in 0..K_FRONT_HALO_RINGS {
                let mut changed = false;
                let ed = &self.edges_;
                for e in 0..ed.ne as usize {
                    let ca = ed.cL[e] as usize;
                    let cb = ed.cR[e] as usize;
                    if self.cell_active_[ca] == 0 || self.cell_active_[cb] == 0 {
                        continue;
                    }
                    let dry_a = st.depth[ca] <= opts.dry_depth;
                    let dry_b = st.depth[cb] <= opts.dry_depth;
                    if dry_a && self.tier_[ca] > self.tier_[cb] {
                        self.tier_[ca] = self.tier_[cb];
                        changed = true;
                    }
                    if dry_b && self.tier_[cb] > self.tier_[ca] {
                        self.tier_[cb] = self.tier_[ca];
                        changed = true;
                    }
                }
                if !changed {
                    break;
                }
            }
            for k in 0..na {
                let i = self.active_cells_[k];
                self.cells_by_tier_[self.tier_[i as usize] as usize].push(i);
            }
        }

        // Faces.
        if (na as i32) < nt / 4 && self.active_samples_ > 0 {
            {
                let ed = &self.edges_;
                for &e in &self.active_faces_ {
                    let ue = e as usize;
                    if self.cell_active_[ed.cL[ue] as usize] == 0 || self.cell_active_[ed.cR[ue] as usize] == 0 {
                        self.q_[ue] = 0.0;
                    }
                }
            }
            self.active_faces_.clear();
            {
                let ed = &self.edges_;
                for &i in &self.active_cells_ {
                    let ui = i as usize;
                    for p in ed.cell_ptr[ui] as usize..ed.cell_ptr[ui + 1] as usize {
                        let e = ed.cell_edge[p];
                        if ed.cL[e as usize] == i && self.cell_active_[ed.cR[e as usize] as usize] != 0 {
                            self.active_faces_.push(e);
                        }
                    }
                }
            }
            self.active_faces_.sort();
            let ed = &self.edges_;
            for &e in &self.active_faces_ {
                let ue = e as usize;
                let ft = cmin(self.tier_[ed.cL[ue] as usize], self.tier_[ed.cR[ue] as usize]);
                self.face_tier_[ue] = ft;
                self.edges_by_tier_[ft as usize].push(e);
            }
        } else {
            self.active_faces_.clear();
            let ed = &self.edges_;
            for e in 0..ed.ne as usize {
                let ca = ed.cL[e] as usize;
                let cb = ed.cR[e] as usize;
                if self.cell_active_[ca] != 0 && self.cell_active_[cb] != 0 {
                    let ft = cmin(self.tier_[ca], self.tier_[cb]);
                    self.face_tier_[e] = ft;
                    self.edges_by_tier_[ft as usize].push(e as i32);
                    self.active_faces_.push(e as i32);
                } else {
                    self.q_[e] = 0.0;
                }
            }
        }

        let mut tk = 0usize;
        while tk < self.cells_by_tier_.len() && tk < self.tier_occupancy_.len() {
            self.tier_occupancy_[tk] += self.cells_by_tier_[tk].len() as i64;
            tk += 1;
        }

        let n_active = self.active_cells_.len() as i32;
        if mesh.n_triangles() > 0 {
            let frac = n_active as f64 / mesh.n_triangles() as f64;
            self.active_frac_min_ = cmin(self.active_frac_min_, frac);
            self.active_frac_max_ = cmax(self.active_frac_max_, frac);
            self.active_frac_sum_ += frac;
            self.active_samples_ += 1;
        }
        // telemetry_: OPENSWMM_2D_MARCHER_TELEMETRY is refused at initialize
        Ok(())
    }

    /// refreshDt0(): tighten (never grow) dt0_ between rebuilds.
    fn refresh_dt0(&mut self, a: &StepAttachments, st: &SurfaceStateData) {
        let opts = a.opts;
        let na = self.active_cells_.len();
        if na == 0 {
            return;
        }
        if self.mode_ == Momentum2D::DiffusiveWave {
            self.refresh_diffusive_slopes(a, st);
        }
        let mut fresh = 1.0e30;
        let nthr = Self::nthreads(opts);
        let stride = 128 / std::mem::size_of::<f64>();
        self.rebuild_dt_partial_.clear();
        self.rebuild_dt_partial_.resize(cmax(1, opts.num_threads) as usize * stride, 1.0e30); // one per requested thread
        let ed = &self.edges_;
        for k in 0..na {
            let tid = omp_static_thread(k, na, nthr);
            let i = self.active_cells_[k] as usize;
            let h = st.depth[i];
            if h <= opts.dry_depth {
                continue;
            }
            let mut speed = 0.0;
            if !self.qcx_.is_empty() && h > 1.0e-6 {
                speed = q_magnitude(self.qcx_[i], self.qcy_[i]) / h;
            }
            let mut dt = cell_cfl_dt(opts.cfl_number, if self.mode_ == Momentum2D::FullSwe { ed.cell_lpos[i] } else { ed.cell_lchar[i] }, h, speed);
            if self.mode_ == Momentum2D::DiffusiveWave {
                dt = cell_diffusive_dt(opts.cfl_number, ed.cell_lchar[i], h, a.mesh.mannings_n[i], self.dw_slope_[i], opts.flux_dh_eps / cmax(ed.cell_lchar[i], 1.0e-9));
            }
            dt /= (1i32 << self.tier_[i]) as f64;
            if dt < self.rebuild_dt_partial_[tid * stride] {
                self.rebuild_dt_partial_[tid * stride] = dt;
            }
        }
        for tid in 0..nthr {
            if self.rebuild_dt_partial_[tid * stride] < fresh {
                fresh = self.rebuild_dt_partial_[tid * stride];
            }
        }
        if fresh < self.dt0_ {
            self.dt0_ = fresh;
        }
    }

    /// refreshDiffusiveSlopes(): each active cell's largest |d eta / dx| over its faces.
    fn refresh_diffusive_slopes(&mut self, _a: &StepAttachments, st: &SurfaceStateData) {
        let ed = &self.edges_;
        for k in 0..self.active_cells_.len() {
            let i = self.active_cells_[k] as usize;
            let mut smax = 0.0;
            for p in ed.cell_ptr[i] as usize..ed.cell_ptr[i + 1] as usize {
                let e = ed.cell_edge[p] as usize;
                let sl = ((st.head[ed.cR[e] as usize] - st.head[ed.cL[e] as usize]) * ed.inv_dx_normal[e]).abs();
                if sl > smax {
                    smax = sl;
                }
            }
            self.dw_slope_[i] = smax;
        }
    }

    /// refreshConveyance(): pick up changed face conveyances (LOCAL_INERTIAL).
    fn refresh_conveyance(&mut self, mesh: &MeshData) {
        if self.mode_ != Momentum2D::LocalInertial || !self.edges_.refresh_conveyance(mesh) {
            return;
        }
        let ed = &self.edges_;
        for e in 0..ed.ne as usize {
            if ed.conveyance[e] == 0.0 {
                self.q_[e] = 0.0;
            }
        }
        if self.qcx_.is_empty() {
            return;
        }
        for i in 0..mesh.n_cells() as usize {
            let mut sx = 0.0;
            let mut sy = 0.0;
            for p in ed.cell_ptr[i] as usize..ed.cell_ptr[i + 1] as usize {
                let e = ed.cell_edge[p] as usize;
                let f = ed.cell_sign[p] as f64 * self.q_[e] * ed.xi[e] * ed.conveyance[e];
                sx += f * ed.cell_arm_x[p];
                sy += f * ed.cell_arm_y[p];
            }
            self.qcx_[i] = sx / mesh.tri_area[i];
            self.qcy_[i] = sy / mesh.tri_area[i];
        }
        for k in 0..self.bc_cell_.len() {
            let i = self.bc_cell_[k] as usize;
            let slot = self.bc_slot_[k] as usize;
            let fo = -self.bc_q_[k] * mesh.edge_length[slot];
            let inv_a = 1.0 / mesh.tri_area[i];
            self.qcx_[i] += fo * (mesh.edge_mx[slot] - mesh.tri_cx[i]) * inv_a;
            self.qcy_[i] += fo * (mesh.edge_my[slot] - mesh.tri_cy[i]) * inv_a;
        }
    }

    // -----------------------------------------------------------------------------------
    // Faces
    // -----------------------------------------------------------------------------------

    /// fireFaces(faces, dt_f, global_step).
    fn fire_faces(&mut self, a: &StepAttachments, st: &mut SurfaceStateData, l: List, dt_f: f64, global_step: bool) -> Result<(), Terminated> {
        let faces = self.take_faces(l);
        let r = match self.mode_ {
            Momentum2D::FullSwe => self.fire_faces_swe(a, st, &faces, dt_f, global_step),
            Momentum2D::DiffusiveWave => {
                self.fire_faces_diffusive(a, st, &faces, dt_f, global_step);
                Ok(())
            }
            Momentum2D::LocalInertial => {
                self.fire_faces_inertial(a, st, &faces, dt_f, global_step);
                Ok(())
            }
        };
        self.put_faces(l, faces);
        r
    }

    /// The positivity budget of a face's exporter (shared by the three closures).
    #[inline]
    fn face_budget(&self, mesh: &MeshData, st: &SurfaceStateData, exp_cell: usize, e: usize, global_step: bool, beta3: f64, beta4: f64) -> f64 {
        let refire = if global_step { 1 } else { 1i32 << (self.tier_[exp_cell] as i32 - self.face_tier_[e] as i32) };
        let beta_share = if self.has_quads_ && mesh.cell_vertex_count(exp_cell as i32) == 4 { beta4 } else { beta3 };
        beta_share / refire as f64 * cmax(st.volume[exp_cell], 0.0)
    }

    /// fireFacesInertial(faces, dt_f, global_step): the local-inertial face law.
    fn fire_faces_inertial(&mut self, a: &StepAttachments, st: &SurfaceStateData, faces: &[i32], dt_f: f64, global_step: bool) {
        let opts = a.opts;
        let mesh = a.mesh;
        let na = faces.len();
        let theta = opts.theta;
        let dry = opts.dry_depth;
        let froude_max = opts.froude_max;
        let advection = opts.advection;
        let physical_time_blend = opts.theta_reference_step > 0.0;
        let blend = momentum_blend_fraction(theta, dt_f, opts.theta_reference_step);
        let beta3 = opts.exchange_beta / 3.0;
        let beta4 = opts.exchange_beta / 4.0;
        let vfr_face = opts.face_reconstruction == FaceDepth2D::VfrFace;
        let stride = self.face_accum_stride_;
        for k in 0..na {
            let e = faces[k] as usize;
            let (conveyance, xi, ca, cb) = (self.edges_.conveyance[e], self.edges_.xi[e], self.edges_.cL[e] as usize, self.edges_.cR[e] as usize);
            if conveyance == 0.0 {
                self.q_[e] = 0.0;
                continue;
            }
            let width = xi * conveyance;
            let ed = &self.edges_;
            let hf = if vfr_face {
                face_flow_depth_vfr(st.head[ca], st.head[cb], ed.ze_lo[e], ed.ze_hi[e])
            } else {
                face_flow_depth(st.head[ca], st.head[cb], ed.zface[e])
            };
            if hf <= dry {
                self.q_[e] = 0.0;
                continue;
            }
            let mut qhat = self.q_[e];
            let mut q_mag = self.q_[e].abs();
            if !self.qcx_.is_empty() {
                let qfx = 0.5 * (self.qcx_[ca] + self.qcx_[cb]);
                let qfy = 0.5 * (self.qcy_[ca] + self.qcy_[cb]);
                let qn = qfx * ed.nx[e] + qfy * ed.ny[e];
                qhat = if physical_time_blend { self.q_[e] + blend * (qn - self.q_[e]) } else { theta * self.q_[e] + (1.0 - theta) * qn };
                q_mag = cmax(q_mag, q_magnitude(qfx, qfy));
            }
            let mut deta = st.head[cb] - st.head[ca];
            if deta.abs() < K_ETA_DEADBAND {
                deta = 0.0;
            }
            let slope = deta * ed.inv_dx_normal[e];
            let mut adv = 0.0;
            if advection && !self.qcx_.is_empty() {
                let hl = st.depth[ca];
                let hr = st.depth[cb];
                if hl > dry && hr > dry {
                    let un_l = (self.qcx_[ca] * ed.nx[e] + self.qcy_[ca] * ed.ny[e]) / hl;
                    let un_r = (self.qcx_[cb] * ed.nx[e] + self.qcy_[cb] * ed.ny[e]) / hr;
                    adv = inertial_advection(self.q_[e], un_l, hl, un_r, hr, ed.inv_dx_normal[e]);
                }
            }
            let mut qn1 = qhat - dt_f * (K_GRAVITY * hf * slope + adv);
            if qn1 != 0.0 {
                qn1 = inertial_face_update(self.q_[e], qhat, hf, dt_f, slope, ed.n2_face[e], q_mag, adv);
            }
            qn1 = froude_cap(qn1, hf, froude_max);
            let exp_cell = if qn1 > 0.0 { ca } else { cb };
            let budget = self.face_budget(mesh, st, exp_cell, e, global_step, beta3, beta4);
            let take = qn1.abs() * width * dt_f;
            if take > budget {
                qn1 *= if take > 0.0 { budget / take } else { 0.0 };
            }
            self.q_[e] = qn1;
            let dm = qn1 * width * dt_f;
            self.face_accum_[stride * e] -= dm;
            self.face_accum_[stride * e + 1] += dm;
        }
        self.face_passes_ += na as i64;
        self.accumulators_pending_ = true;
    }

    /// fireFacesSwe(faces, dt_f, global_step): the FULL_SWE Godunov face flux.
    fn fire_faces_swe(&mut self, a: &StepAttachments, st: &SurfaceStateData, faces: &[i32], dt_f: f64, global_step: bool) -> Result<(), Terminated> {
        let opts = a.opts;
        let mesh = a.mesh;
        for &i in &self.active_cells_ {
            let ui = i as usize;
            let h = st.depth[ui];
            let qx = self.qcx_[ui];
            let qy = self.qcy_[ui];
            if !(h >= 0.0)
                || !st.head[ui].is_finite()
                || !qx.is_finite()
                || !qy.is_finite()
                || !(0.5 * SWE_GRAVITY * h * h).is_finite()
                || (h > opts.dry_depth && !((qx / h) * qx + (qy / h) * qy).is_finite())
            {
                return Err(thrown("throw std::runtime_error", format!("FULL_SWE pressure or momentum exceeds finite range at cell {}", i)));
            }
        }
        let na = faces.len();
        let beta3 = opts.exchange_beta / 3.0;
        let beta4 = opts.exchange_beta / 4.0;
        let dry = opts.dry_depth;
        let depth_offset = st.depth.len();
        if self.second_order_ {
            self.compute_limited_gradients_swe(a, st);
        }
        let so = self.second_order_;
        let stride = self.face_accum_stride_;
        for k in 0..na {
            let e = faces[k] as usize;
            let ed = &self.edges_;
            let ca = ed.cL[e] as usize;
            let cb = ed.cR[e] as usize;
            let mut ff = FaceFlux::default();
            let (mut c_lx, mut c_ly, mut c_rx, mut c_ry) = (0.0, 0.0, 0.0, 0.0);
            let wet;
            if so {
                let ha = st.depth[ca];
                let hb = st.depth[cb];
                let za = st.head[ca] - ha;
                let zb = st.head[cb] - hb;
                let axa = ed.mx[e] - mesh.tri_cx[ca];
                let aya = ed.my[e] - mesh.tri_cy[ca];
                let axb = ed.mx[e] - mesh.tri_cx[cb];
                let ayb = ed.my[e] - mesh.tri_cy[cb];
                let ua = if ha > dry { self.qcx_[ca] / ha } else { 0.0 };
                let va = if ha > dry { self.qcy_[ca] / ha } else { 0.0 };
                let ub = if hb > dry { self.qcx_[cb] / hb } else { 0.0 };
                let vb = if hb > dry { self.qcy_[cb] / hb } else { 0.0 };
                let eta_lf = st.head[ca] + self.gex_[ca] * axa + self.gey_[ca] * aya;
                let eta_rf = st.head[cb] + self.gex_[cb] * axb + self.gey_[cb] * ayb;
                let u_lf = ua + self.gux_[ca] * axa + self.guy_[ca] * aya;
                let v_lf = va + self.gvx_[ca] * axa + self.gvy_[ca] * aya;
                let u_rf = ub + self.gux_[cb] * axb + self.guy_[cb] * ayb;
                let v_rf = vb + self.gvx_[cb] * axb + self.gvy_[cb] * ayb;
                let h_lf = cmax(0.0, ha + self.gex_[depth_offset + ca] * axa + self.gey_[depth_offset + ca] * aya);
                let h_rf = cmax(0.0, hb + self.gex_[depth_offset + cb] * axb + self.gey_[depth_offset + cb] * ayb);
                wet = face_flux_recon_bed(
                    eta_lf, u_lf, v_lf, za, ha, eta_rf, u_rf, v_rf, zb, hb, ed.nx[e], ed.ny[e], dry, &mut ff, &mut c_lx, &mut c_ly, &mut c_rx, &mut c_ry,
                    eta_lf - h_lf, eta_rf - h_rf,
                );
            } else {
                wet = face_flux(
                    st.head[ca], st.depth[ca], self.qcx_[ca], self.qcy_[ca], st.head[cb], st.depth[cb], self.qcx_[cb], self.qcy_[cb], ed.nx[e], ed.ny[e], dry,
                    &mut ff, &mut c_lx, &mut c_ly, &mut c_rx, &mut c_ry,
                );
            }
            if !wet {
                self.q_[e] = 0.0;
                continue;
            }
            let conv = mesh.edge_conveyance[ed.slotL[e] as usize];
            let xi = ed.xi[e];
            let mut fh = ff.mass * conv;
            let mut fmx = ff.mx * conv;
            let mut fmy = ff.my * conv;
            if fh != 0.0 {
                let exp_cell = if fh > 0.0 { ca } else { cb };
                let budget = self.face_budget(mesh, st, exp_cell, e, global_step, beta3, beta4);
                let take = fh.abs() * xi * dt_f;
                if take > budget {
                    let lam = if take > 0.0 { budget / take } else { 0.0 };
                    fh *= lam;
                    fmx *= lam;
                    fmy *= lam;
                }
            }
            self.q_[e] = fh;
            let xdt = xi * dt_f;
            let dm = fh * xdt;
            self.face_accum_[stride * e] -= dm;
            self.face_accum_[stride * e + 1] += dm;
            self.face_accum_[stride * e + 2] += (-fmx + c_lx) * xdt;
            self.face_accum_[stride * e + 4] += (-fmy + c_ly) * xdt;
            self.face_accum_[stride * e + 3] += (fmx + c_rx) * xdt;
            self.face_accum_[stride * e + 5] += (fmy + c_ry) * xdt;
        }
        self.face_passes_ += na as i64;
        self.accumulators_pending_ = true;
        Ok(())
    }

    /// fireFacesDiffusive(faces, dt_f, global_step): the Manning quasi-steady face flux.
    fn fire_faces_diffusive(&mut self, a: &StepAttachments, st: &SurfaceStateData, faces: &[i32], dt_f: f64, global_step: bool) {
        let opts = a.opts;
        let mesh = a.mesh;
        let na = faces.len();
        let beta3 = opts.exchange_beta / 3.0;
        let beta4 = opts.exchange_beta / 4.0;
        let vfr_face = opts.face_reconstruction == FaceDepth2D::VfrFace;
        let dry = opts.dry_depth;
        let flux_dh_eps = opts.flux_dh_eps;
        let stride = self.face_accum_stride_;
        for k in 0..na {
            let e = faces[k] as usize;
            let ed = &self.edges_;
            let ca = ed.cL[e] as usize;
            let cb = ed.cR[e] as usize;
            let hf = if vfr_face {
                face_flow_depth_vfr(st.head[ca], st.head[cb], ed.ze_lo[e], ed.ze_hi[e])
            } else {
                face_flow_depth(st.head[ca], st.head[cb], ed.zface[e])
            };
            if hf <= dry {
                self.q_[e] = 0.0;
                continue;
            }
            let mut deta = st.head[cb] - st.head[ca];
            if deta.abs() < K_ETA_DEADBAND {
                deta = 0.0;
            }
            let slope = deta * ed.inv_dx_normal[e];
            let s_eps = flux_dh_eps * ed.inv_dx_normal[e];
            let mut qn1 = face_discharge(hf, slope, ed.n_face[e], s_eps);
            qn1 *= mesh.edge_conveyance[ed.slotL[e] as usize];
            let xi = ed.xi[e];
            let exp_cell = if qn1 > 0.0 { ca } else { cb };
            let budget = self.face_budget(mesh, st, exp_cell, e, global_step, beta3, beta4);
            let take = qn1.abs() * xi * dt_f;
            if take > budget {
                qn1 *= if take > 0.0 { budget / take } else { 0.0 };
            }
            self.q_[e] = qn1;
            let dm = qn1 * xi * dt_f;
            self.face_accum_[stride * e] -= dm;
            self.face_accum_[stride * e + 1] += dm;
        }
        self.face_passes_ += na as i64;
        self.accumulators_pending_ = true;
    }

    /// computeLimitedGradientsSwe(): RECONSTRUCTION_ORDER 2's limited wet-neighbour gradients.
    fn compute_limited_gradients_swe(&mut self, a: &StepAttachments, st: &SurfaceStateData) {
        let mesh = a.mesh;
        let n = st.depth.len();
        let dry = a.opts.dry_depth;
        let nv = mesh.n_vertices() as usize;
        let mut vertex_bounds: Vec<[f64; 4]> = vec![[0.0; 4]; nv];
        for v in 0..nv {
            let mut ulo = 1e300;
            let mut uhi = -1e300;
            let mut vlo = 1e300;
            let mut vhi = -1e300;
            for p in mesh.vert_stencil_ptr[v] as usize..mesh.vert_stencil_ptr[v + 1] as usize {
                let j = mesh.vert_stencil_idx[p] as usize;
                let h = st.depth[j];
                if self.cell_active_[j] == 0 || h <= dry || st.head[j] - mesh.vz[v] <= dry {
                    continue;
                }
                let u = self.qcx_[j] / h;
                let w = self.qcy_[j] / h;
                ulo = cmin(ulo, u);
                uhi = cmax(uhi, u);
                vlo = cmin(vlo, w);
                vhi = cmax(vhi, w);
            }
            vertex_bounds[v] = [ulo, uhi, vlo, vhi];
        }
        let ed = &self.edges_;
        for k in 0..self.active_cells_.len() {
            let i = self.active_cells_[k];
            let ui = i as usize;
            self.gex_[n + ui] = 0.0;
            self.gey_[n + ui] = 0.0;
            self.gex_[ui] = 0.0;
            self.gey_[ui] = 0.0;
            self.gux_[ui] = 0.0;
            self.guy_[ui] = 0.0;
            self.gvx_[ui] = 0.0;
            self.gvy_[ui] = 0.0;
            let hi = st.depth[ui];
            if hi <= 10.0 * dry {
                continue;
            }
            let wi = [st.head[ui], self.qcx_[ui] / hi, self.qcy_[ui] / hi, hi];
            let mut lo = [wi[0], wi[1], wi[2], wi[3]];
            let mut hiw = [wi[0], wi[1], wi[2], wi[3]];
            let (mut xx, mut xy, mut yy) = (0.0, 0.0, 0.0);
            let mut bx = [0.0f64; 4];
            let mut by = [0.0f64; 4];
            let mut connected: u32 = 0;
            let mut count = 0i32;
            let mut shore = false;
            let begin = ed.cell_ptr[ui] as usize;
            let end = ed.cell_ptr[ui + 1] as usize;
            for p in begin..end {
                let e = ed.cell_edge[p] as usize;
                let j = if ed.cL[e] == i { ed.cR[e] } else { ed.cL[e] } as usize;
                let hj = st.depth[j];
                let ej = st.head[j];
                let sill = cmax(wi[0] - hi, ej - hj);
                if hj <= dry || self.cell_active_[j] == 0 || wi[0] - sill <= dry || ej - sill <= dry {
                    shore = true;
                    continue;
                }
                connected |= 1u32 << (p - begin);
                count += 1;
                let dx = mesh.tri_cx[j] - mesh.tri_cx[ui];
                let dy = mesh.tri_cy[j] - mesh.tri_cy[ui];
                let wt = 1.0 / (dx * dx + dy * dy);
                xx += wt * dx * dx;
                xy += wt * dx * dy;
                yy += wt * dy * dy;
                let wj = [ej, self.qcx_[j] / hj, self.qcy_[j] / hj, hj];
                for m in 0..4 {
                    bx[m] += wt * dx * (wj[m] - wi[m]);
                    by[m] += wt * dy * (wj[m] - wi[m]);
                    lo[m] = cmin(lo[m], wj[m]);
                    hiw[m] = cmax(hiw[m], wj[m]);
                }
            }
            if count < 2 {
                continue;
            }
            let det = xx * yy - xy * xy;
            if det <= 1e-2 * (xx + yy) * (xx + yy) {
                continue;
            }
            let mut gx = [0.0f64; 4];
            let mut gy = [0.0f64; 4];
            let mut phi = [1.0f64, 1.0, 1.0, 1.0];
            for m in 0..4 {
                if shore && (m == 1 || m == 2) {
                    continue;
                }
                gx[m] = (yy * bx[m] - xy * by[m]) / det;
                gy[m] = (xx * by[m] - xy * bx[m]) / det;
            }
            for p in begin..end {
                if connected & (1u32 << (p - begin)) == 0 {
                    continue;
                }
                for m in 0..4 {
                    let wf = wi[m] + gx[m] * ed.cell_arm_x[p] + gy[m] * ed.cell_arm_y[p];
                    phi[m] = cmin(phi[m], bj_limiter(wi[m], wf, lo[m], hiw[m]));
                }
            }
            let nvc = mesh.cell_vertex_count(i);
            if !shore {
                for kk in 0..nvc {
                    let v = mesh.cell_vertex(i, kk) as usize;
                    if st.head[ui] - mesh.vz[v] <= dry {
                        phi[1] = 0.0;
                        phi[2] = 0.0;
                        break;
                    }
                    let ax = mesh.vx[v] - mesh.tri_cx[ui];
                    let ay = mesh.vy[v] - mesh.tri_cy[ui];
                    let b = vertex_bounds[v];
                    phi[1] = cmin(phi[1], bj_limiter(wi[1], wi[1] + gx[1] * ax + gy[1] * ay, cmin(wi[1], b[0]), cmax(wi[1], b[1])));
                    phi[2] = cmin(phi[2], bj_limiter(wi[2], wi[2] + gx[2] * ax + gy[2] * ay, cmin(wi[2], b[2]), cmax(wi[2], b[3])));
                }
            }
            if count < nvc {
                let mut positive = 1.0;
                for e in 0..nvc {
                    let va = mesh.cell_vertex(i, e) as usize;
                    let vb = mesh.cell_vertex(i, (e + 1) % nvc) as usize;
                    let ax = 0.5 * (mesh.vx[va] + mesh.vx[vb]) - mesh.tri_cx[ui];
                    let ay = 0.5 * (mesh.vy[va] + mesh.vy[vb]) - mesh.tri_cy[ui];
                    let dh = phi[3] * (gx[3] * ax + gy[3] * ay);
                    if dh < -hi {
                        positive = cmin(positive, -hi / dh);
                    }
                }
                phi[0] *= positive;
                phi[3] *= positive;
            }
            self.gex_[ui] = phi[0] * gx[0];
            self.gey_[ui] = phi[0] * gy[0];
            self.gex_[n + ui] = phi[3] * gx[3];
            self.gey_[n + ui] = phi[3] * gy[3];
            self.gux_[ui] = phi[1] * gx[1];
            self.guy_[ui] = phi[1] * gy[1];
            self.gvx_[ui] = phi[2] * gx[2];
            self.gvy_[ui] = phi[2] * gy[2];
        }
    }

    // -----------------------------------------------------------------------------------
    // Cells
    // -----------------------------------------------------------------------------------

    /// fireCells(cells, dt_c, tier0).
    fn fire_cells(&mut self, a: &StepAttachments, st: &mut SurfaceStateData, l: List, dt_c: f64, tier0: bool) -> Result<(), Terminated> {
        let cells = self.take_cells(l);
        let r = self.fire_cells_impl(a, st, &cells, dt_c, tier0);
        self.put_cells(l, cells);
        r
    }

    /// fireCellsImpl<kSwe, kFuse = false>(cells, dt_c, tier0) (kFuse needs eight or more
    /// species rows, which phase 2a does not have).
    fn fire_cells_impl(&mut self, a: &StepAttachments, st: &mut SurfaceStateData, cells: &[i32], dt_c: f64, tier0: bool) -> Result<(), Terminated> {
        let mesh = a.mesh;
        let opts = a.opts;
        let nc = cells.len();
        let swe = self.mode_ == Momentum2D::FullSwe;
        let perot = !self.qcx_.is_empty() && !swe;
        let front = self.front_rebuild_;
        let physics_boundary = opts.momentum_physics() && !self.bc_cell_.is_empty();
        let physics_volume = if physics_boundary { st.volume.clone() } else { Vec::new() };
        let boundary_qx = if physics_boundary { self.qcx_.clone() } else { Vec::new() };
        let boundary_qy = if physics_boundary { self.qcy_.clone() } else { Vec::new() };
        let dry = opts.dry_depth;
        let stride = self.face_accum_stride_;
        let nthr = Self::nthreads(opts);

        let mut front_breached = 0i32;
        for k in 0..nc {
            let tid = omp_static_thread(k, nc, nthr);
            let i = cells[k];
            let ui = i as usize;
            let wall_h0 = if swe { st.depth[ui] } else { 0.0 };
            let wall_eta0 = if swe { st.head[ui] } else { 0.0 };
            let mut flux_m3 = 0.0;
            let mut sx = 0.0;
            let mut sy = 0.0;
            let mut dmx = 0.0;
            let mut dmy = 0.0;
            {
                let ed = &self.edges_;
                for p in ed.cell_ptr[ui] as usize..ed.cell_ptr[ui + 1] as usize {
                    let e = ed.cell_edge[p] as usize;
                    if ed.cell_sign[p] > 0 {
                        flux_m3 += self.face_accum_[stride * e];
                        self.face_accum_[stride * e] = 0.0;
                        if swe {
                            dmx += self.face_accum_[stride * e + 2];
                            self.face_accum_[stride * e + 2] = 0.0;
                            dmy += self.face_accum_[stride * e + 4];
                            self.face_accum_[stride * e + 4] = 0.0;
                        }
                    } else {
                        flux_m3 += self.face_accum_[stride * e + 1];
                        self.face_accum_[stride * e + 1] = 0.0;
                        if swe {
                            dmx += self.face_accum_[stride * e + 3];
                            self.face_accum_[stride * e + 3] = 0.0;
                            dmy += self.face_accum_[stride * e + 5];
                            self.face_accum_[stride * e + 5] = 0.0;
                        }
                    }
                    if perot {
                        let mut f = ed.cell_sign[p] as f64 * self.q_[e] * ed.xi[e];
                        if self.mode_ == Momentum2D::LocalInertial {
                            f *= ed.conveyance[e];
                        }
                        sx += f * ed.cell_arm_x[p];
                        sy += f * ed.cell_arm_y[p];
                    }
                }
            }
            st.volume[ui] = cmax(0.0, st.volume[ui] + flux_m3);
            self.apply_cell_sources(a, st, i, dt_c, tid);
            let (eta, depth) = cell_eta_depth(mesh, opts, i, st.volume[ui]);
            st.head[ui] = eta;
            st.depth[ui] = depth;
            if front && self.frontier_[ui] != 0 && st.depth[ui] >= self.h_on_front_ {
                front_breached = 1;
            }
            if perot {
                let inv_a = 1.0 / mesh.tri_area[ui];
                self.qcx_[ui] = sx * inv_a;
                self.qcy_[ui] = sy * inv_a;
            }
            if swe {
                let area = mesh.tri_area[ui];
                let h_new = st.depth[ui];
                if wall_h0 > dry {
                    let so = self.second_order_;
                    for w in self.wall_ptr_[ui] as usize..self.wall_ptr_[ui + 1] as usize {
                        let slot = self.wall_slot_[w] as usize;
                        let nx = mesh.edge_nx[slot];
                        let ny = mesh.edge_ny[slot];
                        let mut eta = wall_eta0;
                        let mut hf = wall_h0;
                        let mut ux = self.qcx_[ui] / wall_h0;
                        let mut uy = self.qcy_[ui] / wall_h0;
                        if so {
                            let ax = mesh.edge_mx[slot] - mesh.tri_cx[ui];
                            let ay = mesh.edge_my[slot] - mesh.tri_cy[ui];
                            let n = st.depth.len();
                            eta += self.gex_[ui] * ax + self.gey_[ui] * ay;
                            hf = cmax(0.0, hf + self.gex_[n + ui] * ax + self.gey_[n + ui] * ay);
                            ux += self.gux_[ui] * ax + self.guy_[ui] * ay;
                            uy += self.gvx_[ui] * ax + self.gvy_[ui] * ay;
                        }
                        let un = ux * nx + uy * ny;
                        let ut = -ux * ny + uy * nx;
                        let ug = -un * nx - ut * ny;
                        let vg = -un * ny + ut * nx;
                        let mut ff = FaceFlux::default();
                        let (mut c_lx, mut c_ly, mut c_rx, mut c_ry) = (0.0, 0.0, 0.0, 0.0);
                        let wet = if so {
                            let z = wall_eta0 - wall_h0;
                            let zf = eta - hf;
                            face_flux_recon_bed(eta, ux, uy, z, wall_h0, eta, ug, vg, z, wall_h0, nx, ny, dry, &mut ff, &mut c_lx, &mut c_ly, &mut c_rx, &mut c_ry, zf, zf)
                        } else {
                            face_flux(eta, wall_h0, self.qcx_[ui], self.qcy_[ui], eta, wall_h0, wall_h0 * ug, wall_h0 * vg, nx, ny, dry, &mut ff, &mut c_lx, &mut c_ly, &mut c_rx, &mut c_ry)
                        };
                        if wet {
                            let xdt = mesh.edge_length[slot] * dt_c;
                            dmx += (-ff.mx + c_lx) * xdt;
                            dmy += (-ff.my + c_ly) * xdt;
                        }
                    }
                }
                if h_new > dry {
                    let mut qx = (area * self.qcx_[ui] + dmx) / area;
                    let mut qy = (area * self.qcy_[ui] + dmy) / area;
                    friction_update(&mut qx, &mut qy, h_new, mesh.mannings_n[ui], dt_c);
                    self.qcx_[ui] = qx;
                    self.qcy_[ui] = qy;
                } else {
                    self.qcx_[ui] = 0.0;
                    self.qcy_[ui] = 0.0;
                }
            }
        }

        self.front_breach_ = self.front_breach_ || front_breached != 0;
        self.flush_source_ledgers(st);

        // Boundary edges owned by cells of this firing (serial: perimeter-sized).
        let vfr_face_bc = opts.face_reconstruction == FaceDepth2D::VfrFace;
        let mut k = 0usize;
        while tier0 && k < self.bc_cell_.len() {
            let kk = k;
            k += 1;
            let i = self.bc_cell_[kk];
            let ui = i as usize;
            if self.tier_[ui] != 0 || self.cell_active_[ui] == 0 {
                continue;
            }
            let b = match a.boundary {
                Some(b) => b,
                None => continue, // unreachable: bc_cell_ is built only with a boundary attached
            };
            let idx = self.bc_slot_[kk];
            let uidx = idx as usize;
            let bt_raw = b.edge_bc_type[uidx];
            let bt_stage = bt_raw == BoundaryType::SpecifiedStage as i8;
            let len = mesh.edge_length[uidx];
            let mut f;
            let mut f_requested_swe = 0.0;
            if self.mode_ == Momentum2D::FullSwe {
                self.swe_bc_dqx_ = 0.0;
                self.swe_bc_dqy_ = 0.0;
                let mut held = [0.0f64; 4];
                let use_held = physics_boundary;
                if use_held {
                    let (h0, h1) = cell_eta_depth(mesh, opts, i, physics_volume[ui]);
                    held[0] = h0;
                    held[1] = h1;
                    held[2] = boundary_qx[ui];
                    held[3] = boundary_qy[ui];
                }
                f = self.boundary_flux_swe(a, b, st, kk, i, dt_c, if use_held { Some(&held) } else { None });
                f_requested_swe = f;
            } else if self.mode_ == Momentum2D::DiffusiveWave {
                f = compute_boundary_edge_flux(mesh, st, Some(b), opts, opts.flux_dh_eps, i, idx);
            } else if bt_stage && len > 1.0e-12 {
                let eta_bc = b.edge_bc_head[uidx];
                let hf;
                if vfr_face_bc {
                    let (mut va_e, mut vb_e) = (0, 0);
                    mesh.cell_edge_vertices(i, MeshData::slot_local(idx), &mut va_e, &mut vb_e);
                    let za = mesh.vz[va_e as usize];
                    let zb = mesh.vz[vb_e as usize];
                    hf = face_depth_from_eta(cmax(st.head[ui], eta_bc), cmin(za, zb), cmax(za, zb));
                } else {
                    hf = face_flow_depth(st.head[ui], eta_bc, mesh.tri_cz[ui]);
                }
                if hf <= opts.dry_depth {
                    self.bc_q_[kk] = 0.0;
                    let mut rc = b.receipts.borrow_mut();
                    let receipt = rc.get_mut(&idx).ok_or_else(|| thrown("throw std::out_of_range", "map::at".to_string()))?;
                    receipt.last_flow = 0.0;
                    receipt.last_dt = dt_c;
                    continue;
                }
                let mut deta = st.head[ui] - eta_bc;
                if deta.abs() < K_ETA_DEADBAND {
                    deta = 0.0;
                }
                let slope = if mesh.cell_vertex_count(i) == 3 { deta * (3.0 * len) / (2.0 * mesh.tri_area[ui]) } else { deta / mesh.edge_dist_c[uidx] };
                let n = mesh.mannings_n[ui];
                let mut qn1 = inertial_face_update(self.bc_q_[kk], self.bc_q_[kk], hf, dt_c, slope, n * n, self.bc_q_[kk].abs(), 0.0);
                qn1 = froude_cap(qn1, hf, opts.froude_max);
                f = qn1 * len;
            } else {
                f = compute_boundary_edge_flux(mesh, st, Some(b), opts, opts.flux_dh_eps, i, idx);
            }
            let mut rc = b.receipts.borrow_mut();
            let receipt = rc.get_mut(&idx).ok_or_else(|| thrown("throw std::out_of_range", "map::at".to_string()))?;
            receipt.last_flow = 0.0;
            receipt.last_dt = dt_c;
            if f == 0.0 {
                self.bc_q_[kk] = 0.0;
                continue;
            }
            let requested_volume = f * dt_c;
            let v_old = st.volume[ui];
            let mut v_new = v_old + dt_c * f;
            if bt_stage && self.mode_ != Momentum2D::FullSwe {
                let v_eq = cell_volume_from_eta(mesh, opts, i, b.edge_bc_head[uidx]);
                if f < 0.0 {
                    v_new = cmax(v_new, cmin(v_old, v_eq));
                } else {
                    v_new = cmin(v_new, cmax(v_old, v_eq));
                }
            }
            if v_new < 0.0 {
                v_new = 0.0;
            }
            f = (v_new - v_old) / dt_c;
            if self.mode_ == Momentum2D::FullSwe && f != f_requested_swe && (self.swe_bc_dqx_ != 0.0 || self.swe_bc_dqy_ != 0.0) {
                let lam = if f_requested_swe != 0.0 { f / f_requested_swe } else { 0.0 };
                self.qcx_[ui] += self.swe_bc_dqx_ * (lam - 1.0);
                self.qcy_[ui] += self.swe_bc_dqy_ * (lam - 1.0);
            }
            self.bc_q_[kk] = if len > 1.0e-12 { f / len } else { 0.0 };
            let applied_volume = v_new - v_old;
            receipt.requested += requested_volume;
            receipt.applied += applied_volume;
            receipt.water_in += cmax(0.0, applied_volume);
            receipt.water_out += cmax(0.0, -applied_volume);
            receipt.last_flow = f;
            drop(rc);
            if f != 0.0 {
                st.volume[ui] = v_new;
                self.bc_accum_[kk] += dt_c * f;
                let (eta, depth) = cell_eta_depth(mesh, opts, i, st.volume[ui]);
                st.head[ui] = eta;
                st.depth[ui] = depth;
            }
            if !self.qcx_.is_empty() && self.mode_ == Momentum2D::LocalInertial && self.bc_q_[kk] != 0.0 {
                let (mut va, mut vb) = (0, 0);
                mesh.cell_edge_vertices(i, MeshData::slot_local(idx), &mut va, &mut vb);
                let mxb = 0.5 * (mesh.vx[va as usize] + mesh.vx[vb as usize]);
                let myb = 0.5 * (mesh.vy[va as usize] + mesh.vy[vb as usize]);
                let fo = -f;
                let inv_a = 1.0 / mesh.tri_area[ui];
                self.qcx_[ui] += fo * (mxb - mesh.tri_cx[ui]) * inv_a;
                self.qcy_[ui] += fo * (myb - mesh.tri_cy[ui]) * inv_a;
            }
        }

        // Live junction exchange at tier-0 cadence (windowless coupling): the orifice law
        // against the live 2D heads and the routing step's 1D heads. Drains cap at the
        // exchange-beta share of the source cell; spills cap at the node's ponded volume for the
        // whole advance (the node_drawn_ ledger). The species bookings of a drain or spill
        // (species_on_) are phase 2b, refused in advance().
        if !self.exch_.is_empty() && tier0 && st.node_coupling && st.nodes_1d {
            if let (Some(pts), Some(nodes)) = (a.node_coupling, a.nodes_1d) {
                for k in 0..pts.len() {
                    let cp = &pts[k];
                    let ci = cp.cell_idx;
                    if ci < 0 || self.cell_active_[ci as usize] == 0 {
                        continue;
                    }
                    let uc = ci as usize;
                    let h_off = if k < self.exch_head_slope_.len() { self.exch_head_slope_[k] * self.exch_tau_ } else { 0.0 };
                    let mut q = compute_node_coupling_q(cp, mesh, &*st, &a.coupling, nodes, opts, None, h_off);
                    if q == 0.0 {
                        continue;
                    }
                    if q > 0.0 {
                        // 2D -> 1D drain: the availability share of the cell
                        q = cmin(q, opts.exchange_beta * cmax(st.volume[uc], 0.0) / dt_c);
                    } else {
                        // 1D -> 2D spill: only the node's ponded water, above the 1D router's
                        // rim volume (node_pond_base: rpt_full_volume under DW, else full_volume)
                        let ni = cp.node_idx as usize;
                        let pond_base = match st.node_pond_base {
                            PondBase::RptFullVolume => nodes.rpt_full_volume[ni],
                            PondBase::FullVolume => nodes.full_volume[ni],
                            PondBase::None => nodes.full_volume[ni],
                        };
                        let ponded_1d = cmax(0.0, nodes.volume[ni] - pond_base);
                        let avail = ponded_1d * opts.vol_1d_to_2d - self.node_drawn_[ni];
                        if avail <= 0.0 {
                            continue;
                        }
                        let want = -q * dt_c;
                        let take = cmin(want, avail);
                        self.spill_deficit_ += cmax(0.0, want - take);
                        self.node_drawn_[ni] += take;
                        q = -take / dt_c;
                    }
                    st.coupling_applied[uc] += -q * dt_c;
                    st.volume[uc] -= q * dt_c;
                    if st.volume[uc] < 0.0 {
                        st.volume[uc] = 0.0;
                    }
                    self.exch_[k] += q * dt_c;
                    let (eta, depth) = cell_eta_depth(mesh, opts, ci, st.volume[uc]);
                    st.head[uc] = eta;
                    st.depth[uc] = depth;
                }
            }
        }
        if tier0 {
            self.exch_tau_ += dt_c;
        }
        self.flush_source_ledgers(st);
        Ok(())
    }

    /// boundaryFluxSwe(k, i, dt_c, held): the FULL_SWE ghost-cell Riemann boundary of bc entry
    /// k, inflow-positive m3/s; books its momentum into the cell.
    fn boundary_flux_swe(&mut self, a: &StepAttachments, b: &BoundaryData, st: &SurfaceStateData, k: usize, i: i32, dt_c: f64, held: Option<&[f64; 4]>) -> f64 {
        let mesh = a.mesh;
        let ui = i as usize;
        let idx = self.bc_slot_[k] as usize;
        let bt = b.edge_bc_type[idx];
        let len = mesh.edge_length[idx];
        if len <= 1.0e-12 {
            return 0.0;
        }
        let nx = mesh.edge_nx[idx];
        let ny = mesh.edge_ny[idx];
        let eta_i = match held {
            Some(h) => h[0],
            None => st.head[ui],
        };
        let h_i = match held {
            Some(h) => h[1],
            None => st.depth[ui],
        };
        let z_i = eta_i - h_i;
        let dry = a.opts.dry_depth;
        let qx_i = match held {
            Some(h) => h[2],
            None => self.qcx_[ui],
        };
        let qy_i = match held {
            Some(h) => h[3],
            None => self.qcy_[ui],
        };
        let mut u_n = 0.0;
        let mut u_t = 0.0;
        if h_i > dry {
            let ux = qx_i / h_i;
            let uy = qy_i / h_i;
            u_n = ux * nx + uy * ny;
            u_t = -ux * ny + uy * nx;
        }
        let mut h_g = h_i;
        let mut un_g = u_n;
        let mut ut_g = u_t;
        let mut eta_g = eta_i;
        let mut prescribed_mass = false;
        let mut q_prescribed = 0.0;
        if bt == BoundaryType::Wall as i8 {
            un_g = -u_n;
        } else if bt == BoundaryType::NormalFlow as i8 {
            let s = b.edge_bed_slope[idx];
            let n_man = mesh.mannings_n[ui];
            prescribed_mass = true;
            if s != 0.0 && h_i > 0.0 && n_man > 0.0 {
                let sgn_s = if s > 0.0 { 1.0 } else { -1.0 };
                let h53 = h_i * glibc_cbrt(h_i * h_i); // LIBM (glibc cbrt)
                q_prescribed = sgn_s * h53 * s.abs().sqrt() / n_man;
                un_g = q_prescribed / h_i;
                ut_g = u_t;
            } else {
                q_prescribed = 0.0;
                un_g = -u_n;
            }
        } else if bt == BoundaryType::SpecifiedStage as i8 {
            let eta_bc = b.edge_bc_head[idx];
            let c_i = (SWE_GRAVITY * cmax(h_i, 0.0)).sqrt();
            if !(h_i > dry && u_n >= c_i) {
                eta_g = eta_bc;
                h_g = if eta_bc - z_i > 0.0 { eta_bc - z_i } else { 0.0 };
                if h_i > dry && h_g > 0.0 {
                    un_g = u_n + 2.0 * (c_i - (SWE_GRAVITY * h_g).sqrt());
                }
            }
        } else if bt == BoundaryType::SpecifiedFlow as i8 || bt == BoundaryType::RatingCurve as i8 {
            let q_out = b.edge_bc_flow[idx];
            if h_i > dry {
                h_g = h_i;
                if q_out < 0.0 {
                    let c_i = (SWE_GRAVITY * h_i).sqrt();
                    if u_n.abs() < c_i {
                        h_g = depth_from_invariant_and_discharge(u_n + 2.0 * c_i, q_out, h_i);
                    }
                }
                eta_g = z_i + h_g;
            } else {
                h_g = critical_depth(q_out.abs());
                eta_g = z_i + h_g;
            }
            un_g = if h_g > 0.0 { q_out / h_g } else { 0.0 };
            ut_g = 0.0;
            prescribed_mass = true;
            q_prescribed = q_out;
        }
        let qxg = h_g * (un_g * nx - ut_g * ny);
        let qyg = h_g * (un_g * ny + ut_g * nx);
        let mut ff = FaceFlux::default();
        let (mut c_lx, mut c_ly, mut c_rx, mut c_ry) = (0.0, 0.0, 0.0, 0.0);
        let wet = face_flux(eta_i, h_i, qx_i, qy_i, eta_g, h_g, qxg, qyg, nx, ny, dry, &mut ff, &mut c_lx, &mut c_ly, &mut c_rx, &mut c_ry);
        let conv = mesh.edge_conveyance[idx];
        if !wet {
            return if prescribed_mass { -q_prescribed * conv * len } else { 0.0 };
        }
        let area = mesh.tri_area[ui];
        let fmx = -ff.mx * conv + c_lx;
        let fmy = -ff.my * conv + c_ly;
        self.swe_bc_dqx_ = fmx * len * dt_c / area;
        self.swe_bc_dqy_ = fmy * len * dt_c / area;
        self.qcx_[ui] += self.swe_bc_dqx_;
        self.qcy_[ui] += self.swe_bc_dqy_;
        if prescribed_mass {
            return -q_prescribed * conv * len;
        }
        -ff.mass * conv * len
    }

    // -----------------------------------------------------------------------------------
    // The macro cycle, RK2 and the optional momentum physics
    // -----------------------------------------------------------------------------------

    /// runMacroCycle(dt0, nsub, time): nsub base substeps of the binary LTS ladder.
    fn run_macro_cycle(&mut self, a: &StepAttachments, st: &mut SurfaceStateData, dt0: f64, nsub: i32, time: f64) -> Result<(), Terminated> {
        let k_tiers = self.cells_by_tier_.len() as i32;
        for s in 0..nsub {
            for k in 0..k_tiers {
                if s % (1 << k) != 0 {
                    continue;
                }
                if !self.edges_by_tier_[k as usize].is_empty() {
                    self.fire_faces(a, st, List::Tier(k as usize), (1i32 << k) as f64 * dt0, false)?;
                }
            }
            for k in 0..k_tiers {
                if (s + 1) % (1 << k) != 0 {
                    continue;
                }
                if k == 0 {
                    self.source_phase_end_ = time + (s + 1) as f64 * dt0;
                }
                if !self.cells_by_tier_[k as usize].is_empty() || k == 0 {
                    self.fire_cells(a, st, List::Tier(k as usize), (1i32 << k) as f64 * dt0, k == 0)?;
                }
            }
            // captureSourceResults(): no recorder in phase 2a
            self.substeps_run_ += 1;
            self.last_steps_ += 1;
        }
        Ok(())
    }

    /// runRk2Step(interval, start): the surface-only SSP-RK2 (Heun) step of RECONSTRUCTION_ORDER
    /// 2 (no groundwater, so no forcing knots, no rollback snapshot and one part).
    fn run_rk2_step(&mut self, a: &StepAttachments, st: &mut SurfaceStateData, interval: f64, start: f64) -> Result<(), Terminated> {
        let mesh = a.mesh;
        let opts = a.opts;
        let _time = start;
        let dt = interval;
        let nt = mesh.n_triangles();
        self.rk_v0_ = st.volume.clone();
        let boundary_receipts0 = a.boundary.map(|b| b.receipts.borrow().clone());
        self.rk_qx0_ = self.qcx_.clone();
        self.rk_qy0_ = self.qcy_.clone();
        self.rk_bc0_ = self.bc_accum_.clone();
        self.rk_ex0_ = self.exch_.clone();
        self.rk_inf0_ = st.infil_applied.clone();
        self.rk_cpl0_ = st.coupling_applied.clone();
        self.rk_drawn0_ = self.node_drawn_.clone();
        let heat_evap0 = st.heat.evaporated_m3.clone();
        let evap0 = st.evap_loss_total;
        for _stage in 0..2 {
            self.rk_transport_stage_ = false;
            // resolveRk2BoundarySignals(time + stage * dt): no signals bound in phase 2a
            self.fire_faces(a, st, List::Active, dt, true)?;
            self.fire_cells(a, st, List::Active, dt, true)?;
            self.rk_transport_stage_ = false;
        }
        for i in 0..nt {
            let ui = i as usize;
            if self.cell_active_[ui] != 0 {
                st.volume[ui] = 0.5 * (self.rk_v0_[ui] + st.volume[ui]);
                self.qcx_[ui] = 0.5 * (self.rk_qx0_[ui] + self.qcx_[ui]);
                self.qcy_[ui] = 0.5 * (self.rk_qy0_[ui] + self.qcy_[ui]);
                let (eta, depth) = cell_eta_depth(mesh, opts, i, st.volume[ui]);
                st.head[ui] = eta;
                st.depth[ui] = depth;
                st.infil_applied[ui] = self.rk_inf0_[ui] + 0.5 * (st.infil_applied[ui] - self.rk_inf0_[ui]);
                st.coupling_applied[ui] = self.rk_cpl0_[ui] + 0.5 * (st.coupling_applied[ui] - self.rk_cpl0_[ui]);
            }
            if st.depth[ui] <= opts.dry_depth {
                self.qcx_[ui] = 0.0;
                self.qcy_[ui] = 0.0;
            }
        }
        // The transport members' blend: every one is empty in phase 2a (no species rows).
        for i in 0..st.heat.evaporated_m3.len() {
            st.heat.evaporated_m3[i] = heat_evap0[i] + 0.5 * (st.heat.evaporated_m3[i] - heat_evap0[i]);
        }
        if let (Some(b), Some(before_all)) = (a.boundary, boundary_receipts0) {
            let mut rc = b.receipts.borrow_mut();
            for (key, now) in rc.iter_mut() {
                let before = before_all.get(key).ok_or_else(|| thrown("throw std::out_of_range", "map::at".to_string()))?;
                let blend = |x: f64, y: f64| x + 0.5 * (y - x);
                now.requested = blend(before.requested, now.requested);
                now.applied = blend(before.applied, now.applied);
                now.water_in = blend(before.water_in, now.water_in);
                now.water_out = blend(before.water_out, now.water_out);
                now.heat_requested = blend(before.heat_requested, now.heat_requested);
                now.heat_applied = blend(before.heat_applied, now.heat_applied);
                for r in 0..now.applied_mass.len() {
                    now.requested_mass[r] = blend(before.requested_mass[r], now.requested_mass[r]);
                    now.applied_mass[r] = blend(before.applied_mass[r], now.applied_mass[r]);
                }
            }
        }
        for k in 0..self.bc_accum_.len() {
            self.bc_accum_[k] = self.rk_bc0_[k] + 0.5 * (self.bc_accum_[k] - self.rk_bc0_[k]);
        }
        for k in 0..self.exch_.len() {
            self.exch_[k] = self.rk_ex0_[k] + 0.5 * (self.exch_[k] - self.rk_ex0_[k]);
        }
        for k in 0..self.node_drawn_.len() {
            self.node_drawn_[k] = self.rk_drawn0_[k] + 0.5 * (self.node_drawn_[k] - self.rk_drawn0_[k]);
        }
        st.evap_loss_total = evap0 + 0.5 * (st.evap_loss_total - evap0);
        self.accumulators_pending_ = false;
        self.substeps_run_ += 2;
        self.last_steps_ += 2;
        Ok(())
    }

    /// applyMomentumPhysics(seconds): half wall drag, half rotation, the eddy-viscosity stress,
    /// half rotation, half wall drag.
    fn apply_momentum_physics(&mut self, a: &StepAttachments, st: &SurfaceStateData, seconds: f64) -> Result<(), Terminated> {
        let opts = a.opts;
        let mesh = a.mesh;
        self.momentum_wall(a, st, seconds);
        let f = opts.coriolis_parameter();
        MomentumPhysics2D::rotate(&mut self.qcx_, &mut self.qcy_, 0.5 * seconds * f);
        if opts.turbulence == Turbulence2D::Smagorinsky {
            self.momentum_physics_.advance(mesh, &self.edges_, &st.depth, &st.head, &mut self.qcx_, &mut self.qcy_, opts.dry_depth, opts.smagorinsky_cs, seconds, opts.eddy_viscosity_background, 0.0)?;
        } else if opts.turbulence == Turbulence2D::MixingLength {
            self.momentum_physics_.advance(
                mesh,
                &self.edges_,
                &st.depth,
                &st.head,
                &mut self.qcx_,
                &mut self.qcy_,
                opts.dry_depth,
                0.0,
                seconds,
                opts.eddy_viscosity_background,
                opts.mixing_length_coefficient,
            )?;
        }
        MomentumPhysics2D::rotate(&mut self.qcx_, &mut self.qcy_, 0.5 * seconds * f);
        self.momentum_wall(a, st, seconds);
        Ok(())
    }

    /// applyMomentumPhysics' `wall` lambda.
    fn momentum_wall(&mut self, a: &StepAttachments, st: &SurfaceStateData, seconds: f64) {
        if self.wall_drag_slots_.is_empty() {
            return;
        }
        let b = match a.boundary {
            Some(b) => b,
            None => return, // unreachable: wall drag slots need the boundary (initialize)
        };
        self.wall_drag_energy_ += MomentumPhysics2D::wall_drag(
            a.mesh,
            b,
            &self.wall_drag_slots_,
            &st.depth,
            &st.head,
            &mut self.qcx_,
            &mut self.qcy_,
            a.opts.dry_depth,
            0.5 * seconds,
            a.opts.cell_closure == CellClosure2D::Vfr,
        );
    }

    /// runMomentumStep(limit, time): one globally stepped substep with the optional momentum
    /// physics split around it; returns the step taken.
    fn run_momentum_step(&mut self, a: &StepAttachments, st: &mut SurfaceStateData, limit: f64, time: f64) -> Result<f64, Terminated> {
        let opts = a.opts;
        let mesh = a.mesh;
        if self.gw_
            || self.common_stage_
            || self.channel_ground_
            || st.transport.active()
            || st.heat.enabled
            || st.snow.active()
            || a.node_coupling.map_or(false, |pts| !pts.is_empty())
            || st.runtime_sources.active()
            || !st.runtime_forcings.rows.is_empty()
        {
            return Err(runtime_error("Optional momentum physics has an unqualified coupled source owner."));
        }
        let mut dt = limit;
        let f = opts.coriolis_parameter().abs();
        if f > 0.0 {
            dt = cmin(dt, 0.25 / f);
        }
        if opts.turbulence == Turbulence2D::Smagorinsky && (opts.smagorinsky_cs > 0.0 || opts.eddy_viscosity_background > 0.0) {
            let mut rx = Vec::new();
            let mut ry = Vec::new();
            let b = self.momentum_physics_.evaluate(mesh, &self.edges_, &st.depth, &st.head, &self.qcx_, &self.qcy_, opts.dry_depth, opts.smagorinsky_cs, &mut rx, &mut ry, -1.0, opts.eddy_viscosity_background, 0.0)?;
            dt = cmin(dt, b);
        } else if opts.turbulence == Turbulence2D::MixingLength && (opts.mixing_length_coefficient > 0.0 || opts.eddy_viscosity_background > 0.0) {
            let mut rx = Vec::new();
            let mut ry = Vec::new();
            let b = self.momentum_physics_.evaluate(
                mesh,
                &self.edges_,
                &st.depth,
                &st.head,
                &self.qcx_,
                &self.qcy_,
                opts.dry_depth,
                0.0,
                &mut rx,
                &mut ry,
                -1.0,
                opts.eddy_viscosity_background,
                opts.mixing_length_coefficient,
            )?;
            dt = cmin(dt, b);
        }
        let x0 = self.qcx_.clone();
        let y0 = self.qcy_.clone();
        let wall_energy0 = self.wall_drag_energy_;
        let mut trial = 0;
        loop {
            self.qcx_ = x0.clone();
            self.qcy_ = y0.clone();
            self.wall_drag_energy_ = wall_energy0;
            self.apply_momentum_physics(a, st, 0.5 * dt)?;
            let mut stable = opts.max_timestep;
            for i in 0..mesh.n_cells() as usize {
                if st.depth[i] > opts.dry_depth {
                    stable = cmin(stable, cell_cfl_dt(opts.cfl_number, self.edges_.cell_lpos[i], st.depth[i], self.qcx_[i].hypot(self.qcy_[i]) / st.depth[i]));
                }
            }
            if dt <= stable {
                break;
            }
            if trial == 39 {
                self.qcx_ = x0;
                self.qcy_ = y0;
                self.wall_drag_energy_ = wall_energy0;
                return Err(runtime_error("Optional momentum physics cannot satisfy the wave CFL."));
            }
            dt = cmin(0.5 * dt, stable);
            trial += 1;
        }
        if self.second_order_ {
            self.run_rk2_step(a, st, dt, time)?;
        } else {
            self.run_macro_cycle(a, st, dt, 1, time)?;
        }
        self.apply_momentum_physics(a, st, 0.5 * dt)?;
        Ok(dt)
    }

    /// run_stats(): (nsteps, nrhs, last_h, active fraction min/max/mean, tier occupancy).
    pub fn run_stats(&self) -> RunStats {
        let mut s = RunStats::default();
        s.nsteps = self.substeps_run_;
        s.nrhs = self.face_passes_;
        s.last_h = self.last_dt_;
        if self.active_samples_ > 0 {
            s.active_frac_min = self.active_frac_min_;
            s.active_frac_max = self.active_frac_max_;
            s.active_frac_mean = self.active_frac_sum_ / self.active_samples_ as f64;
        }
        s.n_tiers = cmin(self.cells_by_tier_.len(), self.tier_occupancy_.len()) as i32;
        for k in 0..s.n_tiers as usize {
            s.tier_cells[k] = self.tier_occupancy_[k];
        }
        s
    }
}

/// ISurfaceSolver::RunStats (member order and initializers upstream's).
#[derive(Debug, Clone)]
pub struct RunStats {
    pub nsteps: i64,
    pub nrhs: i64,
    pub last_h: f64,
    pub avg_h: f64,
    pub active_frac_min: f64,
    pub active_frac_mean: f64,
    pub active_frac_max: f64,
    pub tier_cells: [i64; 8],
    pub n_tiers: i32,
}

impl Default for RunStats {
    fn default() -> Self {
        RunStats { nsteps: 0, nrhs: 0, last_h: 0.0, avg_h: 0.0, active_frac_min: -1.0, active_frac_mean: -1.0, active_frac_max: -1.0, tier_cells: [0; 8], n_tiers: 0 }
    }
}
