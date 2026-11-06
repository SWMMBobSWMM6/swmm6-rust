// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/fv/ExplicitFvSolver.cpp,
// ExplicitFvSolver.hpp, INetworkSolver.hpp (FvStepForcing, RunStats, Divergence) and
// CheckedNetworkAdvance.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The CPU explicit FV backend's step: advance / advanceImpl (the substep loop with its
// post-step census and retries, local time stepping with the macro cycle and its
// accumulators, the divergence guard), rebuildActiveLists, censusDt, the stable-step
// functions, assignTiers, fireFaces / fireCells / fireNodes, runMacroCycle,
// settleAccumulators, emitConduitSteps, saveState / restoreState, takeSubstep, reinitialize,
// run_stats, divergence and advanceNetworkChecked. The face fluxes are in explicit_flux.rs,
// the cell and node updates in explicit_nodes.rs.
//
// Differences in form, never in result:
//   * The solver keeps no pointers to the Router's mesh, state and forcing: the step works
//     through `FvStep`, which borrows the solver, the mesh, the state and the forcing for
//     the length of one advance() (upstream: mesh_, state_, forcing_).
//   * FvStepForcing's arrays are owned vectors (upstream: pointers into the Router's
//     vectors). The forcing's refresh hook (a C function pointer and its user pointer) is a
//     closure the Router passes to advance(); it is handed the state the solver holds and
//     the forcing to refill, which is what upstream's hook reads and writes through its
//     pointers. node_withdrawal points at the solver's own withdrawal_rates_ upstream; here
//     it is a flag.
//   * The conduit step sink (a C function pointer and its user pointer) is a closure passed
//     to advance().
//   * Every OpenMP loop is the serial loop: the reference runs on one thread, so each loop
//     runs in order (and none of the step's reductions crosses threads).
//   * perf:: counters and timers, and the OPENSWMM_FV_DT_TRACE print, are left out.
//
// NOT TRANSLATED (Router::fv_not_translated refuses them before the step changes anything):
// FV_PRESSURIZED_IMPLICIT (PressurizedHeadSolver, the nonlinear pressure-storage candidate),
// the TPA pressure closure (updateTpaFlags and the TPA branches), RK2 time integration
// (rkSave / rkAverage), unsteady friction (computeUfGradients), species transport on the FV
// mesh (n_species > 0). No gate deck reaches any of them.

use crate::core::cmath::{cmax, cmin};
use crate::core::engine_step_routing::StepError;
use crate::core::options::fv::TimeIntegration;
use crate::hydraulics::fv::explicit_solver::{ExplicitFvSolver, K_MAX_LTS_TIERS};
use crate::hydraulics::fv::kernels as k;
use crate::hydraulics::fv::{NetworkMeshData, NetworkStateData, K_NODE_OUTFALL, K_NODE_VIRTUAL};

pub const K_STEP_ACCEPT_RATIO: f64 = 0.5;
pub const K_MAX_STEP_RETRIES: i32 = 8;
pub const K_MIN_SUBSTEP: f64 = 1.0e-6;
pub const K_DIVERGED_VELOCITY: f64 = 1000.0;
pub const K_DIVERGED_DEPTH: f64 = 10000.0;
pub const K_REBUILD_INTERVAL: i32 = 8;
pub const K_RETIER_EVERY_CYCLES: i32 = 8;

/// Per-step forcing handed to advance() (INetworkSolver.hpp FvStepForcing).
#[derive(Debug, Clone, Default)]
pub struct FvStepForcing {
    pub node_lateral: Vec<f64>,
    pub node_fixed_head: Vec<f64>,
    /// None where upstream's pointer is null (ctx.surface_outfall_link_limit empty)
    pub node_source_max: Option<Vec<f64>>,
    pub structure_flow: Vec<f64>,
    pub link_q_cap: Vec<f64>,
    pub conduit_loss: Vec<f64>,
    /// a refresh hook is set (FV_STRUCTURE_COUPLING SUBSTEP)
    pub has_refresh: bool,
    pub n_nodes: i32,
    pub n_links: i32,
    /// node_withdrawal points at the solver's withdrawal_rates_ (upstream: non-null)
    pub node_withdrawal: bool,
}

impl FvStepForcing {
    #[inline]
    pub fn fixed(&self, un: usize) -> bool {
        self.node_fixed_head[un].is_finite()
    }
}

/// The forcing's refresh hook: (the state the solver holds, the forcing to refill,
/// seconds since the start of the routing step).
pub type FvRefreshFn<'a> = dyn FnMut(&NetworkStateData, &mut FvStepForcing, f64) -> Result<(), StepError> + 'a;
/// The conduit step sink: (conduit row, step, duration).
pub type FvStepSinkFn<'a> = dyn FnMut(i32, f64, f64) + 'a;

/// INetworkSolver::RunStats.
#[derive(Debug, Clone, Copy)]
pub struct RunStats {
    pub nsteps: i64,
    pub nflux: i64,
    pub last_h: f64,
    pub avg_h: f64,
    pub min_h: f64,
    pub active_frac_min: f64,
    pub active_frac_mean: f64,
    pub active_frac_max: f64,
    pub tier_cells: [i64; 8],
    pub n_tiers: i32,
    pub dt_argmin_pressurized: i64,
    pub dt_argmin_band: i64,
    pub dt_argmin_free: i64,
    pub dt_argmin_node: i64,
    pub n_macro_cycles: i64,
    pub n_macro_rejected: i64,
}

impl Default for RunStats {
    fn default() -> Self {
        RunStats {
            nsteps: 0,
            nflux: 0,
            last_h: 0.0,
            avg_h: 0.0,
            min_h: 0.0,
            active_frac_min: -1.0,
            active_frac_mean: -1.0,
            active_frac_max: -1.0,
            tier_cells: [0; 8],
            n_tiers: 0,
            dt_argmin_pressurized: 0,
            dt_argmin_band: 0,
            dt_argmin_free: 0,
            dt_argmin_node: 0,
            n_macro_cycles: 0,
            n_macro_rejected: 0,
        }
    }
}

/// INetworkSolver::Divergence.
#[derive(Debug, Clone, Copy)]
pub struct Divergence {
    pub link: i32,
    pub value: f64,
    pub what: &'static str,
}

/// One advance(): the solver with the mesh, state and forcing it works on.
pub struct FvStep<'a> {
    pub s: &'a mut ExplicitFvSolver,
    pub m: &'a NetworkMeshData,
    pub st: &'a mut NetworkStateData,
    pub fc: &'a mut FvStepForcing,
}

impl ExplicitFvSolver {
    pub fn last_num_steps(&self) -> i64 {
        self.last_nsteps_
    }
    pub fn last_step_size(&self) -> f64 {
        self.last_h_
    }
    pub fn suggested_step(&self) -> f64 {
        self.suggested_h_
    }

    /// divergence(d): the cell state the last advance() left, when it diverged.
    pub fn divergence(&self) -> Option<Divergence> {
        if self.div_link_ < 0 {
            return None;
        }
        Some(Divergence { link: self.div_link_, value: self.div_value_, what: self.div_what_.unwrap_or("") })
    }

    /// run_stats()
    pub fn run_stats(&self) -> RunStats {
        let mut s = RunStats::default();
        s.nsteps = self.total_steps_;
        s.nflux = self.total_flux_;
        s.last_h = self.last_h_;
        s.min_h = self.min_h_;
        s.avg_h = if self.total_steps_ > 0 { self.sim_time_ / self.total_steps_ as f64 } else { 0.0 };
        if self.active_n_ > 0 {
            s.active_frac_mean = self.active_sum_ / self.active_n_ as f64;
            s.active_frac_min = self.active_min_;
            s.active_frac_max = self.active_max_;
        }
        s.n_tiers = cmin(self.lts_tiers_, K_MAX_LTS_TIERS as i32);
        for kk in 0..s.n_tiers {
            s.tier_cells[kk as usize] = self.tier_occupancy_[kk as usize];
        }
        s.dt_argmin_pressurized = self.dt_argmin_pressurized_;
        s.dt_argmin_band = self.dt_argmin_band_;
        s.dt_argmin_free = self.dt_argmin_free_;
        s.dt_argmin_node = self.dt_argmin_node_;
        s.n_macro_cycles = self.n_macro_cycles_;
        s.n_macro_rejected = self.n_macro_rejected_;
        s
    }

    /// reinitialize(t0)
    pub fn reinitialize(&mut self, mesh: &NetworkMeshData, state: &mut NetworkStateData) {
        if !self.attached_ {
            return;
        }
        self.lts_valid_ = false;
        self.lists_valid_ = false;
        self.census_count_ = 0;
        self.dt_cache_ = 0.0;
        self.dt_census_ = 0.0;
        self.state_clean_ = false;
        self.refresh_depths(mesh, state);
        self.refresh_node_areas(mesh, state);
    }

    /// nodeDepthFromVolume(node, volume)
    pub fn node_depth_from_volume(&self, mesh: &NetworkMeshData, state: &NetworkStateData, node: i32, volume: f64) -> f64 {
        let un = node as usize;
        if volume <= 0.0 {
            return 0.0;
        }
        let full = mesh.node_full_depth[un];
        if mesh.node_can_pond[un] != 0 && full > 0.0 {
            let v_full = self.node_volume_from_depth(mesh, state, node, full);
            if volume > v_full {
                return full + (volume - v_full) / mesh.node_ponded_area[un];
            }
        }
        let off = mesh.node_vol_off[un];
        if off < 0 {
            let a = state.node_surf_area[un];
            return if a > 0.0 { volume / a } else { 0.0 };
        }
        let n = crate::hydraulics::fv::K_NODE_VOL_SAMPLES;
        let dmax = mesh.node_vol_dmax[un];
        let dd = dmax / (n - 1) as f64;
        let vmax = mesh.node_vol_tbl[(off + n - 1) as usize];
        if volume >= vmax {
            let atop = mesh.node_vol_atop[un];
            if atop > 0.0 {
                return dmax + (volume - vmax) / atop;
            }
            return if volume == vmax { dmax } else { full + mesh.node_sur_depth[un] + 1.0 };
        }
        let mut lo = 0i32;
        let mut hi = n - 1;
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if mesh.node_vol_tbl[(off + mid) as usize] <= volume {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let v0 = mesh.node_vol_tbl[(off + lo) as usize];
        let v1 = mesh.node_vol_tbl[(off + lo + 1) as usize];
        let t = if v1 > v0 { (volume - v0) / (v1 - v0) } else { 0.0 };
        (lo as f64 + t) * dd
    }

    /// advance(t_current, t_target, forcing) (the pressure-storage candidate is never
    /// enabled on the run path: enablePressureStorageCandidate has no caller in the engine).
    #[allow(clippy::too_many_arguments)]
    pub fn advance(&mut self, mesh: &NetworkMeshData, state: &mut NetworkStateData, t_current: f64, t_target: f64, forcing: &mut FvStepForcing, refresh: Option<&mut FvRefreshFn>, sink: Option<&mut FvStepSinkFn>) -> Result<f64, StepError> {
        let mut st = FvStep { s: self, m: mesh, st: state, fc: forcing };
        st.advance_impl(t_current, t_target, refresh, sink)
    }
}

/// advanceNetworkChecked(solver, start, target, forcing)
#[allow(clippy::too_many_arguments)]
pub fn advance_network_checked(solver: &mut ExplicitFvSolver, mesh: &NetworkMeshData, state: &mut NetworkStateData, start: f64, target: f64, forcing: &mut FvStepForcing, refresh: Option<&mut FvRefreshFn>, sink: Option<&mut FvStepSinkFn>) -> Result<f64, StepError> {
    if !start.is_finite() || !target.is_finite() || target < start {
        return Err(StepError::Exception("FV routing times must be finite and ordered."));
    }
    let reached = solver.advance(mesh, state, start, target, forcing, refresh, sink)?;
    let tolerance = 16.0 * f64::EPSILON * cmax(cmax(1.0, start.abs()), target.abs());
    if !reached.is_finite() || reached < start || (reached - target).abs() > tolerance {
        return Err(StepError::Exception("FV backend did not reach the requested routing time; partial results were not published."));
    }
    Ok(reached)
}

impl<'a> FvStep<'a> {
    #[inline]
    pub(crate) fn tpa_cell(&self, uc: usize) -> bool {
        self.s.tpa_ && self.st.cell_tpa[uc] != 0
    }
    #[inline]
    pub(crate) fn algebraic_active(&self, n: i32) -> bool {
        self.s.algebraic_active(self.m, self.st, n)
    }
    #[inline]
    pub(crate) fn node_volume_from_depth(&self, n: i32, d: f64) -> f64 {
        self.s.node_volume_from_depth(self.m, self.st, n, d)
    }
    #[inline]
    pub(crate) fn node_depth_from_volume(&self, n: i32, v: f64) -> f64 {
        self.s.node_depth_from_volume(self.m, self.st, n, v)
    }
    #[inline]
    pub(crate) fn refresh_depths(&mut self) {
        self.s.refresh_depths(self.m, self.st);
    }
    pub(crate) fn lts_eligible(&self) -> bool {
        self.s.opts_.lts && self.s.opts_.time_integration != TimeIntegration::Rk2 && self.st.n_species == 0
    }
    #[inline]
    pub(crate) fn face_is_live(&self, f: i32) -> bool {
        self.s.all_faces_live_ || (!self.s.live_stamp_.is_empty() && self.s.live_stamp_[f as usize] == self.s.live_gen_)
    }

    /// anyPressurizedCell()
    pub(crate) fn any_pressurized_cell(&self) -> bool {
        if !self.s.opts_.pressurized_implicit {
            return false;
        }
        let m = self.m;
        for c in 0..m.n_cells() {
            let uc = c as usize;
            if self.s.cell_active_[uc] == 0 {
                continue;
            }
            let g = &m.geom[m.cell_geom[uc] as usize];
            if g.is_open == 0 && self.st.cell_h[uc] >= g.y_crown {
                return true;
            }
        }
        false
    }

    /// advanceImpl(t_current, t_target, forcing)
    fn advance_impl(&mut self, t_current: f64, t_target: f64, mut refresh: Option<&mut FvRefreshFn>, mut sink: Option<&mut FvStepSinkFn>) -> Result<f64, StepError> {
        if !self.s.attached_ || t_target <= t_current {
            return Ok(t_target);
        }
        let m = self.m;

        // Per-routing-step bookkeeping.
        self.s.unmet_withdrawal_.iter_mut().for_each(|x| *x = 0.0);
        if self.s.water_tracking_ {
            self.s.water_.clear();
        }
        self.s.node_exch_.iter_mut().for_each(|x| *x = 0.0);
        self.s.node_in_.iter_mut().for_each(|x| *x = 0.0);
        self.s.node_out_.iter_mut().for_each(|x| *x = 0.0);
        self.s.flood_vol_.iter_mut().for_each(|x| *x = 0.0);
        self.s.inlet_control_.iter_mut().for_each(|x| *x = 0);
        self.s.cell_q_int_.iter_mut().for_each(|x| *x = 0.0);
        self.s.dummy_vol_.iter_mut().for_each(|x| *x = 0.0);

        if !self.s.state_clean_ {
            self.refresh_depths();
        }
        for n in 0..m.n_nodes() {
            let un = n as usize;
            if m.node_kind[un] == K_NODE_VIRTUAL {
                continue;
            }
            if !self.s.node_alg_.is_empty() && self.algebraic_active(n) {
                continue;
            }
            let v = self.node_volume_from_depth(n, cmax(0.0, self.st.node_head[un] - m.node_invert[un]));
            self.st.node_volume[un] = v;
        }
        self.refresh_struct_flows();

        let mut t = t_current;
        let mut steps: i64 = 0;
        self.s.lists_valid_ = false;
        self.s.census_count_ = 0;

        while t < t_target {
            // FV_STRUCTURE_COUPLING SUBSTEP: the caller re-evaluates the head-dependent
            // boundary flows at the top of every cycle but the first.
            if steps > 0 && self.fc.has_refresh && self.s.opts_.structure_coupling == crate::core::options::fv::StructureCoupling::Substep {
                if let Some(r) = refresh.as_mut() {
                    r(&*self.st, &mut *self.fc, t - t_current)?;
                }
                self.refresh_struct_flows();
            }

            if !self.s.lists_valid_ || self.s.since_rebuild_ >= K_REBUILD_INTERVAL {
                self.rebuild_active_lists();
                self.s.census_count_ = 0;
            }

            let press = self.any_pressurized_cell();
            if press && self.s.lts_valid_ {
                self.settle_accumulators();
                self.s.lts_valid_ = false;
            }

            if !press && self.lts_eligible() {
                let now = if self.s.lts_valid_ && self.s.lts_tiers_ > 1 { self.census_dt(false) } else { 0.0 };
                let stiffness_moved = self.s.lts_valid_ && self.s.lts_tiers_ > 1 && (now < 0.9 * self.s.lts_dt0_ || now > 2.0 * self.s.lts_dt0_);
                let mut retiered = false;
                if !self.s.lts_valid_ || self.s.lts_countdown_ <= 0 || stiffness_moved {
                    retiered = true;
                    self.settle_accumulators();
                    let mut dt0 = self.s.lts_dt0_;
                    self.s.lts_tiers_ = self.assign_tiers(&mut dt0, 0);
                    self.s.lts_dt0_ = dt0;
                    self.s.lts_valid_ = true;
                    self.s.lts_countdown_ = K_RETIER_EVERY_CYCLES;
                    self.s.census_count_ = 0;
                }
                let kk = self.s.lts_tiers_;
                let dt0 = self.s.lts_dt0_;
                self.s.lts_countdown_ -= 1;
                let nsub: i32 = 1 << (kk - 1);
                let span = nsub as f64 * dt0;
                if kk > 1 && span <= t_target - t {
                    if self.s.opts_.order < 2 {
                        self.reconstruct_state();
                    }
                    self.save_state();
                    self.run_macro_cycle(dt0, nsub);

                    if self.census_dt(false) < K_STEP_ACCEPT_RATIO * dt0 {
                        self.s.n_macro_rejected_ += 1;
                        self.restore_state();
                        self.s.lts_valid_ = false;
                        self.s.acc_a_.iter_mut().for_each(|x| *x = 0.0);
                        self.s.acc_q_.iter_mut().for_each(|x| *x = 0.0);
                        self.s.acc_nvol_.iter_mut().for_each(|x| *x = 0.0);
                        self.s.acc_dirty_ = false;
                        self.s.census_count_ = 0;
                    } else {
                        self.s.n_macro_cycles_ += 1;
                        if let Some(sk) = sink.as_mut() {
                            self.emit_conduit_steps(dt0, span, true, sk);
                        }
                        t += span;
                        steps += nsub as i64;
                        self.s.since_rebuild_ += nsub;
                        self.s.last_h_ = dt0;
                        self.s.min_h_ = if self.s.min_h_ <= 0.0 { dt0 } else { cmin(self.s.min_h_, dt0) };
                        self.s.dt_cache_ = dt0;
                        self.s.census_count_ = 0;
                        continue;
                    }
                } else if kk > 1 && self.s.lts_fit_ && (t_target - t) >= 2.0 * dt0 {
                    let remaining = t_target - t;
                    let mut kf = 1 + crate::core::cmath::cvt_i32((remaining / dt0).log2().floor()); // LIBM
                    kf = cmin(kf, kk);
                    if kf >= 2 {
                        self.settle_accumulators();
                        let mut dt0f = self.s.lts_dt0_;
                        let kf2 = self.assign_tiers(&mut dt0f, kf);
                        self.s.lts_valid_ = false;
                        self.s.lts_countdown_ = 0;
                        self.s.census_count_ = 0;
                        let nsubf: i32 = 1 << (kf2 - 1);
                        let spanf = nsubf as f64 * dt0f;
                        if kf2 > 1 && spanf <= remaining {
                            if self.s.opts_.order < 2 {
                                self.reconstruct_state();
                            }
                            self.save_state();
                            self.run_macro_cycle(dt0f, nsubf);
                            if self.census_dt(false) < K_STEP_ACCEPT_RATIO * dt0f {
                                self.s.n_macro_rejected_ += 1;
                                self.restore_state();
                                self.s.acc_a_.iter_mut().for_each(|x| *x = 0.0);
                                self.s.acc_q_.iter_mut().for_each(|x| *x = 0.0);
                                self.s.acc_nvol_.iter_mut().for_each(|x| *x = 0.0);
                                self.s.acc_dirty_ = false;
                            } else {
                                self.s.n_macro_cycles_ += 1;
                                if let Some(sk) = sink.as_mut() {
                                    self.emit_conduit_steps(dt0f, spanf, true, sk);
                                }
                                t += spanf;
                                steps += nsubf as i64;
                                self.s.since_rebuild_ += nsubf;
                                self.s.last_h_ = dt0f;
                                self.s.min_h_ = if self.s.min_h_ <= 0.0 { dt0f } else { cmin(self.s.min_h_, dt0f) };
                                self.s.dt_cache_ = dt0f;
                                continue;
                            }
                        }
                    }
                }
                if !retiered && self.s.lts_valid_ && now > 0.0 {
                    self.s.dt_census_ = now;
                    self.s.census_count_ = cmax(1, self.s.opts_.cfl_census_interval);
                }
            }

            if self.s.census_count_ <= 0 {
                self.s.dt_census_ = self.census_dt(press);
                self.s.census_count_ = cmax(1, self.s.opts_.cfl_census_interval);
            }
            self.s.census_count_ -= 1;

            let remaining = t_target - t;
            let mut dt = cmin(self.s.dt_census_, remaining);
            if !(dt > 0.0) {
                dt = remaining;
            }
            if dt < K_MIN_SUBSTEP {
                dt = cmin(K_MIN_SUBSTEP, remaining);
            }

            let mut attempt = 0;
            loop {
                self.save_state();
                // (RK2 is refused before the step; EULER:)
                let _pressure_ok = self.take_substep(dt);
                if attempt >= K_MAX_STEP_RETRIES || dt <= K_MIN_SUBSTEP {
                    break;
                }
                let press2 = self.any_pressurized_cell();
                let dt_post = self.census_dt(press2);
                if dt_post >= K_STEP_ACCEPT_RATIO * dt {
                    break;
                }
                self.restore_state();
                dt = cmax(0.9 * dt_post, K_MIN_SUBSTEP);
                if dt > remaining {
                    dt = remaining;
                }
                attempt += 1;
            }
            self.s.dt_cache_ = dt;

            if dt < self.s.dt_census_ && dt < remaining {
                self.s.dt_census_ = dt;
                self.s.census_count_ = 0;
            }

            if let Some(sk) = sink.as_mut() {
                let d0 = if dt >= remaining { cmax(dt, self.s.dt_census_) } else { dt };
                self.emit_conduit_steps(d0, dt, false, sk);
            }
            t += dt;
            steps += if self.s.opts_.time_integration == TimeIntegration::Rk2 { 2 } else { 1 };
            self.s.since_rebuild_ += 1;
            self.s.last_h_ = dt;
            self.s.min_h_ = if self.s.min_h_ <= 0.0 { dt } else { cmin(self.s.min_h_, dt) };

            if steps > 2000000 {
                break;
            }
        }

        self.settle_accumulators();
        self.s.state_clean_ = !self.s.opts_.pressurized_implicit && !self.s.tpa_;

        for un in 0..self.s.node_pub_.len() {
            if self.s.node_pub_[un] != 0 && !self.algebraic_active(un as i32) {
                self.s.node_pub_[un] = 0;
            }
        }

        // Divergence guard (issue #156 R3).
        self.s.div_link_ = -1;
        {
            let mut worst = 1.0f64;
            let n_cells = m.n_cells();
            for c in 0..n_cells {
                let uc = c as usize;
                let a = self.st.cell_a[uc];
                let q = self.st.cell_q[uc];
                let h = self.st.cell_h[uc];
                let u = self.s.cell_u_[uc];
                let ratio;
                let what: &'static str;
                let mut bad = h;
                if !a.is_finite() || !q.is_finite() || !h.is_finite() {
                    ratio = f64::INFINITY;
                    what = "non-finite state";
                    bad = if !a.is_finite() { a } else if !q.is_finite() { q } else { h };
                } else if u.abs() > K_DIVERGED_VELOCITY * worst {
                    ratio = u.abs() / K_DIVERGED_VELOCITY;
                    what = "velocity (ft/s)";
                } else if h.abs() > K_DIVERGED_DEPTH * worst {
                    ratio = h.abs() / K_DIVERGED_DEPTH;
                    what = "depth magnitude (ft)";
                } else {
                    continue;
                }
                worst = cmin(ratio, f64::MAX);
                self.s.div_link_ = m.conduit_link[m.cell_conduit[uc] as usize];
                self.s.div_value_ = if what.as_bytes()[0] == b'v' { u } else { bad };
                self.s.div_what_ = Some(what);
                if !ratio.is_finite() {
                    break;
                }
            }
        }

        self.s.last_nsteps_ = steps;
        self.s.total_steps_ += steps;
        self.s.sim_time_ += t - t_current;
        self.s.suggested_h_ = if self.s.dt_cache_ > 0.0 && self.s.dt_cache_ < 1.0e29 { self.s.dt_cache_ } else { 0.0 };
        Ok(t)
    }

    /// rebuildActiveLists()
    pub(crate) fn rebuild_active_lists(&mut self) {
        let m = self.m;
        let nc = m.n_cells();
        let nf = m.n_faces();
        self.s.lts_valid_ = false;
        self.s.f_mass_.iter_mut().for_each(|x| *x = 0.0);
        self.s.f_mom_.iter_mut().for_each(|x| *x = 0.0);
        self.s.f_corr_l_.iter_mut().for_each(|x| *x = 0.0);
        self.s.f_corr_r_.iter_mut().for_each(|x| *x = 0.0);
        if !self.s.opts_.compaction {
            self.s.cell_active_.iter_mut().for_each(|x| *x = 1);
            self.s.active_faces_.resize(nf as usize, 0);
            for f in 0..nf {
                self.s.active_faces_[f as usize] = f;
            }
            self.s.lists_valid_ = true;
            self.s.since_rebuild_ = 0;
            return;
        }
        let have_qlat = !self.s.cell_qlat_.is_empty();
        for c in 0..nc {
            let uc = c as usize;
            let on = self.st.cell_h[uc] > k::K_DRY_DEPTH || self.st.cell_q[uc] != 0.0 || (have_qlat && self.s.cell_qlat_[uc] != 0.0) || self.tpa_cell(uc);
            self.s.cell_active_[uc] = if on { 1 } else { 0 };
        }
        for f in 0..nf {
            let uf = f as usize;
            let nd = m.face_node[uf];
            if nd < 0 {
                continue;
            }
            let un = nd as usize;
            let independent_head = self.s.node_alg_[un] == 0 || (!m.node_can_pond.is_empty() && m.node_can_pond[un] != 0);
            if !independent_head && self.st.node_head[un] - m.face_zb[uf] <= k::K_DRY_DEPTH {
                continue;
            }
            let cell = if m.face_cl[uf] >= 0 { m.face_cl[uf] } else { m.face_cr[uf] };
            if cell >= 0 {
                self.s.cell_active_[cell as usize] = 1;
            }
        }
        self.s.halo_prev_.clear();
        self.s.halo_prev_.extend_from_slice(&self.s.cell_active_);
        let halo_levels = K_REBUILD_INTERVAL * if self.s.opts_.time_integration == TimeIntegration::Rk2 { 2 } else { 1 };
        for _level in 0..halo_levels {
            let mut grew = false;
            for f in 0..nf {
                let uf = f as usize;
                let cl = m.face_cl[uf];
                let cr = m.face_cr[uf];
                if cl < 0 || cr < 0 {
                    continue;
                }
                let ul = cl as usize;
                let ur = cr as usize;
                if self.s.halo_prev_[ul] != 0 && self.s.cell_active_[ur] == 0 {
                    self.s.cell_active_[ur] = 1;
                    grew = true;
                }
                if self.s.halo_prev_[ur] != 0 && self.s.cell_active_[ul] == 0 {
                    self.s.cell_active_[ul] = 1;
                    grew = true;
                }
            }
            for un in 0..self.s.node_alg_.len() {
                if self.s.node_alg_[un] == 0 {
                    continue;
                }
                let begin = m.node_face_ptr[un];
                let end = m.node_face_ptr[un + 1];
                let mut connected = false;
                let mut p = begin;
                while p < end && !connected {
                    let f = m.node_face_idx[p as usize] as usize;
                    let c = if m.face_cl[f] >= 0 { m.face_cl[f] } else { m.face_cr[f] };
                    connected = c >= 0 && self.s.halo_prev_[c as usize] != 0;
                    p += 1;
                }
                if !connected {
                    continue;
                }
                for p in begin..end {
                    let f = m.node_face_idx[p as usize] as usize;
                    let c = if m.face_cl[f] >= 0 { m.face_cl[f] } else { m.face_cr[f] };
                    if c >= 0 && self.s.cell_active_[c as usize] == 0 {
                        self.s.cell_active_[c as usize] = 1;
                        grew = true;
                    }
                }
            }
            if !grew {
                break;
            }
            self.s.halo_prev_.clear();
            self.s.halo_prev_.extend_from_slice(&self.s.cell_active_);
        }
        self.s.active_faces_.clear();
        for f in 0..nf {
            let uf = f as usize;
            let cl = m.face_cl[uf];
            let cr = m.face_cr[uf];
            let la = cl >= 0 && self.s.cell_active_[cl as usize] != 0;
            let ra = cr >= 0 && self.s.cell_active_[cr as usize] != 0;
            if la || ra {
                self.s.active_faces_.push(f);
            }
        }
        let frac = if nf > 0 { self.s.active_faces_.len() as f64 / nf as f64 } else { 1.0 };
        self.s.active_sum_ += frac;
        self.s.active_n_ += 1;
        self.s.active_min_ = if self.s.active_min_ < 0.0 { frac } else { cmin(self.s.active_min_, frac) };
        self.s.active_max_ = if self.s.active_max_ < 0.0 { frac } else { cmax(self.s.active_max_, frac) };
        self.s.lists_valid_ = true;
        self.s.since_rebuild_ = 0;
    }

    /// passThroughFarCell(face, nd)
    pub(crate) fn pass_through_far_cell(&self, face: i32, nd: i32) -> i32 {
        let m = self.m;
        if nd < 0 || self.s.node_pass_.is_empty() || self.s.node_pass_[nd as usize] == 0 || !self.algebraic_active(nd) {
            return -1;
        }
        let und = nd as usize;
        for p in m.node_face_ptr[und]..m.node_face_ptr[und + 1] {
            let of = m.node_face_idx[p as usize];
            if of == face {
                continue;
            }
            let uof = of as usize;
            return if m.face_cl[uof] >= 0 { m.face_cl[uof] } else { m.face_cr[uof] };
        }
        -1
    }

    /// censusDt(press_edit). press_edit is true only under FV_PRESSURIZED_IMPLICIT, which is
    /// refused before the step, so the face mode is always PressurizedHeadSolver::kNone.
    pub(crate) fn census_dt(&mut self, press_edit: bool) -> f64 {
        debug_assert!(!press_edit);
        let _ = press_edit;
        let m = self.m;
        let mut dt = f64::MAX;
        let mut trace_face: i32 = -1;
        for ai in 0..self.s.active_faces_.len() {
            let f = self.s.active_faces_[ai];
            let uf = f as usize;
            let cl = m.face_cl[uf];
            let cr = m.face_cr[uf];
            let nd = m.face_node[uf];
            let mut speed = 0.0f64;
            let mut dx_ref = 1.0e30f64;
            let consider_cell = |this: &Self, cell: i32, speed: &mut f64, dx_ref: &mut f64| {
                let uc = cell as usize;
                *dx_ref = cmin(*dx_ref, m.cell_dx[uc]);
                let h = this.st.cell_h[uc];
                let tpa_c = this.tpa_cell(uc);
                if !tpa_c && h <= k::K_DRY_DEPTH {
                    return;
                }
                let g = &m.geom[m.cell_geom[uc] as usize];
                *speed = cmax(*speed, this.s.cell_u_[uc].abs() + k::celerity(this.st.cell_a[uc], if tpa_c { g.t_slot } else { this.s.cell_t_[uc] }));
            };
            if cl >= 0 {
                consider_cell(self, cl, &mut speed, &mut dx_ref);
            }
            if cr >= 0 {
                consider_cell(self, cr, &mut speed, &mut dx_ref);
            }
            if nd >= 0 {
                let other = if cl >= 0 { cl } else { cr };
                let far = if self.s.census_pass_far_ { self.pass_through_far_cell(f, nd) } else { -1 };
                if far >= 0 {
                    let ufar = far as usize;
                    let tpa_f = self.tpa_cell(ufar);
                    if tpa_f || self.st.cell_h[ufar] > k::K_DRY_DEPTH {
                        let gfar = &m.geom[m.cell_geom[ufar] as usize];
                        speed = cmax(speed, self.s.cell_u_[ufar].abs() + k::celerity(self.st.cell_a[ufar], if tpa_f { gfar.t_slot } else { self.s.cell_t_[ufar] }));
                    }
                } else if other >= 0 {
                    let uo = other as usize;
                    let g = &m.geom[m.cell_geom[uo] as usize];
                    let hg = self.st.node_head[nd as usize] - m.face_zb[uf];
                    if hg > k::K_DRY_DEPTH {
                        let ag = k::area_of_depth(g, hg);
                        speed = cmax(speed, self.s.cell_u_[uo].abs() + k::celerity(ag, k::width_of_depth(g, hg)));
                    }
                }
            }
            if speed > 1.0e-12 && dx_ref < 1.0e29 {
                let dt_f = self.s.opts_.cfl * dx_ref / speed;
                if dt_f < dt {
                    dt = dt_f;
                    trace_face = f;
                }
            }
        }
        let mut node_owns_dt = false;
        let nn = m.n_nodes();
        if self.lts_eligible() {
            if dt >= f64::MAX * 0.5 {
                dt = 1.0e30;
            }
            self.count_argmin(node_owns_dt, trace_face);
            return dt;
        }
        for n in 0..nn {
            let un = n as usize;
            if m.node_kind[un] == K_NODE_VIRTUAL {
                continue;
            }
            if !self.s.node_alg_.is_empty() && self.s.node_alg_[un] != 0 {
                if self.s.opts_.node_feedback_dt {
                    // (press_edit && nodePressFolded(n)): press_edit is always false here
                    let dt_n = self.algebraic_node_stable_dt(n);
                    if dt_n < dt {
                        dt = dt_n;
                        node_owns_dt = true;
                    }
                }
                continue;
            }
            if m.node_kind[un] == K_NODE_OUTFALL {
                continue;
            }
            let as_ = self.st.node_surf_area[un];
            if !(as_ > 0.0) {
                continue;
            }
            let b = m.node_face_ptr[un];
            let e = m.node_face_ptr[un + 1];
            for p in b..e {
                let f = m.node_face_idx[p as usize];
                let uf = f as usize;
                let cell = if m.face_cl[uf] >= 0 { m.face_cl[uf] } else { m.face_cr[uf] };
                if cell < 0 {
                    continue;
                }
                let uc = cell as usize;
                if self.s.cell_active_[uc] == 0 {
                    continue;
                }
                let h = self.st.cell_h[uc];
                if h <= k::K_DRY_DEPTH {
                    continue;
                }
                let g = &m.geom[m.cell_geom[uc] as usize];
                let t = k::width_of_depth(g, h);
                if t <= 0.0 {
                    continue;
                }
                let c_wave = k::celerity(self.st.cell_a[uc], t);
                let dt_n = k::face_cfl_dt(self.s.opts_.cfl, as_ / t, self.s.cell_u_[uc], c_wave);
                if dt_n < dt {
                    dt = dt_n;
                    node_owns_dt = true;
                }
            }
        }
        if dt >= f64::MAX * 0.5 {
            dt = 1.0e30;
        }
        self.count_argmin(node_owns_dt, trace_face);
        dt
    }

    /// censusDt's count_argmin lambda (slot program R0 attribution counters).
    fn count_argmin(&mut self, node_owns_dt: bool, trace_face: i32) {
        let m = self.m;
        if node_owns_dt {
            self.s.dt_argmin_node_ += 1;
            return;
        }
        if trace_face < 0 {
            return;
        }
        let uf = trace_face as usize;
        let cl = m.face_cl[uf];
        let cr = m.face_cr[uf];
        let nd = m.face_node[uf];
        let mut h_max = 0.0f64;
        let mut gm: Option<usize> = None;
        let mut consider = |cell: i32, h: f64| {
            if cell < 0 || h <= h_max {
                return;
            }
            h_max = h;
            gm = Some(m.cell_geom[cell as usize] as usize);
        };
        if cl >= 0 {
            consider(cl, self.st.cell_h[cl as usize]);
        }
        if cr >= 0 {
            consider(cr, self.st.cell_h[cr as usize]);
        }
        if nd >= 0 {
            let other = if cl >= 0 { cl } else { cr };
            if other >= 0 {
                consider(other, self.st.node_head[nd as usize] - m.face_zb[uf]);
            }
        }
        let gi = match gm {
            Some(g) => g,
            None => return,
        };
        let g = &m.geom[gi];
        if h_max >= g.y_full {
            self.s.dt_argmin_pressurized_ += 1;
        } else if h_max >= g.y_crown {
            self.s.dt_argmin_band_ += 1;
        } else {
            self.s.dt_argmin_free_ += 1;
        }
    }

    /// saveState()
    pub(crate) fn save_state(&mut self) {
        let s = &mut *self.s;
        s.save_cell_a_.clone_from(&self.st.cell_a);
        s.save_cell_q_.clone_from(&self.st.cell_q);
        s.save_cell_phi_.clone_from(&self.st.cell_phi);
        s.save_node_vol_.clone_from(&self.st.node_volume);
        s.save_node_head_.clone_from(&self.st.node_head);
        s.save_exch_.clone_from(&s.node_exch_);
        s.save_in_.clone_from(&s.node_in_);
        s.save_out_.clone_from(&s.node_out_);
        s.save_flood_.clone_from(&s.flood_vol_);
        s.save_qint_.clone_from(&s.cell_q_int_);
        s.save_carry_.clone_from(&s.node_carry_);
        s.save_dummy_.clone_from(&s.dummy_vol_);
        s.save_unmet_.clone_from(&s.unmet_withdrawal_);
        if s.water_tracking_ {
            s.save_water_.clone_from(&s.water_);
        }
        if s.tpa_ {
            s.save_tpa_.clone_from(&self.st.cell_tpa);
        }
    }

    /// restoreState()
    pub(crate) fn restore_state(&mut self) {
        {
            let s = &mut *self.s;
            self.st.cell_a.clone_from(&s.save_cell_a_);
            self.st.cell_q.clone_from(&s.save_cell_q_);
            self.st.cell_phi.clone_from(&s.save_cell_phi_);
            self.st.node_volume.clone_from(&s.save_node_vol_);
            self.st.node_head.clone_from(&s.save_node_head_);
            s.node_exch_.clone_from(&s.save_exch_);
            s.node_in_.clone_from(&s.save_in_);
            s.node_out_.clone_from(&s.save_out_);
            s.flood_vol_.clone_from(&s.save_flood_);
            s.cell_q_int_.clone_from(&s.save_qint_);
            s.node_carry_.clone_from(&s.save_carry_);
            s.dummy_vol_.clone_from(&s.save_dummy_);
            s.unmet_withdrawal_.clone_from(&s.save_unmet_);
            if s.water_tracking_ {
                s.water_.clone_from(&s.save_water_);
            }
            s.face_lat_spill_.iter_mut().for_each(|x| *x = 0.0);
            if s.tpa_ {
                self.st.cell_tpa.clone_from(&s.save_tpa_);
            }
        }
        self.refresh_depths();
    }

    /// takeSubstep(dt, forcing) on the paths the step translates (no TPA, no pressurized
    /// implicit pass, no unsteady friction, no species: refused before the step).
    pub(crate) fn take_substep(&mut self, dt: f64) -> bool {
        self.reconstruct_state();
        self.compute_fluxes();
        self.s.press_step_ = false;
        self.relax_node_fluxes(dt);
        self.limit_positivity(dt);
        // reconstructScalars(dt): no species (n_species == 0), a no-op
        self.update_cells(dt);
        self.update_nodes(dt);
        // dispersionSolve(dt): no species, a no-op
        true
    }

    /// cellStableDt(c)
    pub(crate) fn cell_stable_dt(&self, c: i32) -> f64 {
        let m = self.m;
        let uc = c as usize;
        let g = &m.geom[m.cell_geom[uc] as usize];
        let dx = m.cell_dx[uc];
        let mut speed = 0.0f64;
        let h = self.st.cell_h[uc];
        let tpa_c = self.tpa_cell(uc);
        if tpa_c || h > k::K_DRY_DEPTH {
            speed = self.s.cell_u_[uc].abs() + k::celerity(self.st.cell_a[uc], if tpa_c { g.t_slot } else { self.s.cell_t_[uc] });
        }
        let faces = [m.cell_face0[uc], m.cell_face1[uc]];
        for &f in faces.iter() {
            let uf = f as usize;
            let nd = m.face_node[uf];
            if nd < 0 {
                continue;
            }
            let far = if self.s.census_pass_far_ { self.pass_through_far_cell(f, nd) } else { -1 };
            if far >= 0 {
                let ufar = far as usize;
                let tpa_f = self.tpa_cell(ufar);
                if tpa_f || self.st.cell_h[ufar] > k::K_DRY_DEPTH {
                    let gfar = &m.geom[m.cell_geom[ufar] as usize];
                    speed = cmax(speed, self.s.cell_u_[ufar].abs() + k::celerity(self.st.cell_a[ufar], if tpa_f { gfar.t_slot } else { self.s.cell_t_[ufar] }));
                }
                continue;
            }
            let hg = self.st.node_head[nd as usize] - m.face_zb[uf];
            if hg <= k::K_DRY_DEPTH {
                continue;
            }
            speed = cmax(speed, self.s.cell_u_[uc].abs() + k::celerity(k::area_of_depth(g, hg), k::width_of_depth(g, hg)));
        }
        if speed > 1.0e-12 { self.s.opts_.cfl * dx / speed } else { 1.0e30 }
    }

    /// algebraicNodeStableDt(n)
    pub(crate) fn algebraic_node_stable_dt(&self, n: i32) -> f64 {
        let m = self.m;
        let un = n as usize;
        let b = m.node_face_ptr[un];
        let e = m.node_face_ptr[un + 1];
        let mut store_min = f64::MAX;
        let mut admit_sum = 0.0f64;
        for p in b..e {
            let f = m.node_face_idx[p as usize];
            let uf = f as usize;
            let cell = if m.face_cl[uf] >= 0 { m.face_cl[uf] } else { m.face_cr[uf] };
            if cell < 0 {
                continue;
            }
            let uc = cell as usize;
            let h = self.st.cell_h[uc];
            let tpa_c = self.tpa_cell(uc);
            if !tpa_c && h <= k::K_DRY_DEPTH {
                continue;
            }
            let g = &m.geom[m.cell_geom[uc] as usize];
            let t = if tpa_c { g.t_slot } else { self.s.cell_t_[uc] };
            if t <= 0.0 {
                continue;
            }
            admit_sum += t * k::celerity(self.st.cell_a[uc], t);
            if tpa_c || h >= g.y_crown {
                let store = 0.5 * m.cell_dx[uc] * t;
                if store < store_min {
                    store_min = store;
                }
            }
        }
        if store_min >= f64::MAX * 0.5 {
            return 1.0e30;
        }
        if !(admit_sum > 0.0) {
            return 1.0e30;
        }
        let tau = store_min / admit_sum;
        if tau > 0.0 { self.s.opts_.cfl * tau } else { 1.0e30 }
    }

    /// nodeStableDt(n)
    pub(crate) fn node_stable_dt(&self, n: i32) -> f64 {
        let m = self.m;
        let un = n as usize;
        if m.node_kind[un] == K_NODE_VIRTUAL {
            return 1.0e30;
        }
        if !self.s.node_alg_.is_empty() && self.s.node_alg_[un] != 0 {
            return 1.0e30;
        }
        if m.node_kind[un] == K_NODE_OUTFALL {
            return 1.0e30;
        }
        let as_ = self.st.node_surf_area[un];
        if !(as_ > 0.0) {
            return 1.0e30;
        }
        let mut dt = 1.0e30f64;
        let b = m.node_face_ptr[un];
        let e = m.node_face_ptr[un + 1];
        for p in b..e {
            let f = m.node_face_idx[p as usize];
            let uf = f as usize;
            let cell = if m.face_cl[uf] >= 0 { m.face_cl[uf] } else { m.face_cr[uf] };
            if cell < 0 {
                continue;
            }
            let uc = cell as usize;
            let h = self.st.cell_h[uc];
            let tpa_c = self.tpa_cell(uc);
            if !tpa_c && h <= k::K_DRY_DEPTH {
                continue;
            }
            let g = &m.geom[m.cell_geom[uc] as usize];
            let t = if tpa_c { g.t_slot } else { self.s.cell_t_[uc] };
            if t <= 0.0 {
                continue;
            }
            dt = cmin(dt, k::face_cfl_dt(self.s.opts_.cfl, as_ / t, self.s.cell_u_[uc], k::celerity(self.st.cell_a[uc], t)));
        }
        dt
    }

    /// assignTiers(dt0, k_cap_override)
    pub(crate) fn assign_tiers(&mut self, dt0: &mut f64, k_cap_override: i32) -> i32 {
        let m = self.m;
        let nc = m.n_cells();
        let nf = m.n_faces();
        let nn = m.n_nodes();
        if self.st.n_species > 0 {
            return 1;
        }
        let mut dt_cell = vec![1.0e30f64; nc as usize];
        let mut dt_node = vec![1.0e30f64; nn as usize];
        let mut dt_min = 1.0e30f64;
        let mut dt_argmin_cell: i32 = -1;
        for c in 0..nc {
            let uc = c as usize;
            dt_cell[uc] = self.cell_stable_dt(c);
            if dt_cell[uc] < dt_min {
                dt_min = dt_cell[uc];
                dt_argmin_cell = c;
            }
        }
        for n in 0..nn {
            let un = n as usize;
            dt_node[un] = self.node_stable_dt(n);
            if dt_node[un] < dt_min {
                dt_min = dt_node[un];
                dt_argmin_cell = -2;
            }
        }
        if self.s.opts_.node_feedback_dt {
            for n in 0..nn {
                let un = n as usize;
                if self.s.node_alg_.is_empty() || self.s.node_alg_[un] == 0 {
                    continue;
                }
                let dtn = self.algebraic_node_stable_dt(n);
                if dtn >= 1.0e29 {
                    continue;
                }
                let b = m.node_face_ptr[un];
                let e = m.node_face_ptr[un + 1];
                for q in b..e {
                    let uf = m.node_face_idx[q as usize] as usize;
                    let cl = m.face_cl[uf];
                    let cr = m.face_cr[uf];
                    if cl >= 0 {
                        dt_cell[cl as usize] = cmin(dt_cell[cl as usize], dtn);
                    }
                    if cr >= 0 {
                        dt_cell[cr as usize] = cmin(dt_cell[cr as usize], dtn);
                    }
                }
                if dtn < dt_min {
                    dt_min = dtn;
                    dt_argmin_cell = -2;
                }
            }
        }
        if !(dt_min > 0.0) || dt_min >= 1.0e29 {
            return 1;
        }
        *dt0 = dt_min;
        if dt_argmin_cell == -2 {
            self.s.dt_argmin_node_ += 1;
        } else if dt_argmin_cell >= 0 {
            let uc = dt_argmin_cell as usize;
            let g = &m.geom[m.cell_geom[uc] as usize];
            let h = self.st.cell_h[uc];
            if h >= g.y_full {
                self.s.dt_argmin_pressurized_ += 1;
            } else if h >= g.y_crown {
                self.s.dt_argmin_band_ += 1;
            } else {
                self.s.dt_argmin_free_ += 1;
            }
        }
        let mut k_cap = cmin(cmax(self.s.opts_.lts_max_tiers, 1), K_MAX_LTS_TIERS as i32);
        if k_cap_override > 0 {
            k_cap = cmin(k_cap, k_cap_override);
        }
        let tier_of = |dt: f64| -> i32 {
            if dt >= 1.0e29 {
                return k_cap - 1;
            }
            let kk = crate::core::cmath::cvt_i32((dt / dt_min).log2().floor()); // LIBM
            cmin(cmax(kk, 0), k_cap - 1)
        };
        let s = &mut *self.s;
        s.cell_tier_.clear();
        s.cell_tier_.resize(nc as usize, 0);
        s.node_tier_.clear();
        s.node_tier_.resize(nn as usize, 0);
        s.face_tier_.clear();
        s.face_tier_.resize(nf as usize, 0);
        for c in 0..nc as usize {
            s.cell_tier_[c] = tier_of(dt_cell[c]) as u8;
        }
        for n in 0..nn as usize {
            s.node_tier_[n] = tier_of(dt_node[n]) as u8;
        }
        for us in 0..m.struct_link.len() {
            s.node_tier_[m.struct_n1[us] as usize] = 0;
            s.node_tier_[m.struct_n2[us] as usize] = 0;
        }
        // Grade: neighbouring volumes differ by at most one tier.
        #[derive(Clone, Copy, PartialEq)]
        enum T {
            Cell(usize),
            Node(usize),
        }
        for _pass in 0..k_cap {
            let mut changed = false;
            for f in 0..nf as usize {
                let cl = m.face_cl[f];
                let cr = m.face_cr[f];
                let nd = m.face_node[f];
                let mut a = if cl >= 0 { Some(T::Cell(cl as usize)) } else { None };
                let mut b = if cr >= 0 { Some(T::Cell(cr as usize)) } else { None };
                if a.is_none() && nd >= 0 {
                    a = Some(T::Node(nd as usize));
                }
                if b.is_none() && nd >= 0 {
                    b = Some(T::Node(nd as usize));
                }
                let (a, b) = match (a, b) {
                    (Some(a), Some(b)) => (a, b),
                    _ => continue,
                };
                let get = |s: &ExplicitFvSolver, t: T| match t {
                    T::Cell(i) => s.cell_tier_[i],
                    T::Node(i) => s.node_tier_[i],
                };
                let set = |s: &mut ExplicitFvSolver, t: T, v: u8| match t {
                    T::Cell(i) => s.cell_tier_[i] = v,
                    T::Node(i) => s.node_tier_[i] = v,
                };
                // integer promotion: *b + 1 is an int
                if get(s, a) as i32 > get(s, b) as i32 + 1 {
                    let v = (get(s, b) as i32 + 1) as u8;
                    set(s, a, v);
                    changed = true;
                }
                if get(s, b) as i32 > get(s, a) as i32 + 1 {
                    let v = (get(s, a) as i32 + 1) as u8;
                    set(s, b, v);
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
        let mut k_max: i32 = 0;
        for c in 0..nc as usize {
            if dt_cell[c] < 1.0e29 {
                k_max = cmax(k_max, s.cell_tier_[c] as i32);
            }
        }
        for n in 0..nn as usize {
            if dt_node[n] < 1.0e29 {
                k_max = cmax(k_max, s.node_tier_[n] as i32);
            }
        }
        let kk = k_max + 1;
        if kk <= 1 {
            return 1;
        }
        for t in s.cell_tier_.iter_mut() {
            *t = cmin(*t, (kk - 1) as u8);
        }
        for t in s.node_tier_.iter_mut() {
            *t = cmin(*t, (kk - 1) as u8);
        }
        for un in 0..nn as usize {
            if s.node_alg_[un] == 0 {
                continue;
            }
            let mut t = (kk - 1) as u8;
            let nb = m.node_face_ptr[un];
            let ne = m.node_face_ptr[un + 1];
            for p in nb..ne {
                let uf = m.node_face_idx[p as usize] as usize;
                let cell = if m.face_cl[uf] >= 0 { m.face_cl[uf] } else { m.face_cr[uf] };
                if cell >= 0 {
                    t = cmin(t, s.cell_tier_[cell as usize]);
                }
            }
            s.node_tier_[un] = t;
        }
        for uf in 0..nf as usize {
            let mut kf = k_cap - 1;
            let cl = m.face_cl[uf];
            let cr = m.face_cr[uf];
            let nd = m.face_node[uf];
            if cl >= 0 {
                kf = cmin(kf, s.cell_tier_[cl as usize] as i32);
            }
            if cr >= 0 {
                kf = cmin(kf, s.cell_tier_[cr as usize] as i32);
            }
            if nd >= 0 {
                kf = cmin(kf, s.node_tier_[nd as usize] as i32);
            }
            s.face_tier_[uf] = kf as u8;
        }
        s.cells_by_tier_ = vec![Vec::new(); kk as usize];
        s.faces_by_tier_ = vec![Vec::new(); kk as usize];
        s.nodes_by_tier_ = vec![Vec::new(); kk as usize];
        for c in 0..nc {
            let uc = c as usize;
            if s.cell_active_[uc] == 0 {
                continue;
            }
            let t = s.cell_tier_[uc] as usize;
            s.cells_by_tier_[t].push(c);
            s.tier_occupancy_[t] += 1;
        }
        for ai in 0..s.active_faces_.len() {
            let f = s.active_faces_[ai];
            let t = s.face_tier_[f as usize] as usize;
            s.faces_by_tier_[t].push(f);
        }
        for n in 0..nn {
            let un = n as usize;
            if m.node_kind[un] == K_NODE_VIRTUAL {
                continue;
            }
            let t = s.node_tier_[un] as usize;
            s.nodes_by_tier_[t].push(n);
        }
        s.due_f_upto_ = vec![Vec::new(); kk as usize];
        s.due_c_upto_ = vec![Vec::new(); kk as usize];
        s.due_n_upto_ = vec![Vec::new(); kk as usize];
        for j in 0..kk as usize {
            if j > 0 {
                s.due_f_upto_[j] = s.due_f_upto_[j - 1].clone();
                s.due_c_upto_[j] = s.due_c_upto_[j - 1].clone();
                s.due_n_upto_[j] = s.due_n_upto_[j - 1].clone();
            }
            let (fb, cb, nb) = (s.faces_by_tier_[j].clone(), s.cells_by_tier_[j].clone(), s.nodes_by_tier_[j].clone());
            s.due_f_upto_[j].extend_from_slice(&fb);
            s.due_c_upto_[j].extend_from_slice(&cb);
            s.due_n_upto_[j].extend_from_slice(&nb);
        }
        s.acc_a_.clear();
        s.acc_a_.resize(nc as usize, 0.0);
        s.acc_q_.clear();
        s.acc_q_.resize(nc as usize, 0.0);
        s.acc_nvol_.clear();
        s.acc_nvol_.resize(nn as usize, 0.0);
        kk
    }

    /// emitConduitSteps(dt0, dur, tiered)
    pub(crate) fn emit_conduit_steps(&self, dt0: f64, dur: f64, tiered: bool, sink: &mut FvStepSinkFn) {
        let m = self.m;
        for r in 0..m.n_conduits() {
            let ur = r as usize;
            let mut dt = dt0;
            if tiered {
                let begin = m.conduit_cell_begin[ur];
                let end = begin + m.conduit_cell_count[ur];
                let mut tmin: i32 = 255;
                for c in begin..end {
                    tmin = cmin(tmin, self.s.cell_tier_[c as usize] as i32);
                }
                if tmin < 255 {
                    dt = (1i32 << tmin) as f64 * dt0;
                }
            }
            sink(r, dt, dur);
        }
    }

    /// runMacroCycle(dt0, nsub, forcing)
    pub(crate) fn run_macro_cycle(&mut self, dt0: f64, nsub: i32) {
        // (tpa_: updateTpaFlags; refused before the step)
        let kk = self.s.cells_by_tier_.len() as i32;
        let trailing_zeros = |v: i32| -> i32 {
            let mut j = 0;
            while ((v >> j) & 1) == 0 {
                j += 1;
            }
            j
        };
        let due_f = std::mem::take(&mut self.s.due_f_upto_);
        let due_c = std::mem::take(&mut self.s.due_c_upto_);
        let due_n = std::mem::take(&mut self.s.due_n_upto_);
        for s in 0..nsub {
            let jf = cmin(if s == 0 { kk - 1 } else { trailing_zeros(s) }, kk - 1) as usize;
            let jv = cmin(trailing_zeros(s + 1), kk - 1) as usize;
            if self.s.opts_.order >= 2 {
                self.reconstruct_state();
            }
            self.fire_faces(&due_f[jf], dt0);
            // (unsteady friction: computeUfGradients; refused before the step)
            self.fire_cells(&due_c[jv], dt0);
            self.fire_nodes(&due_n[jv], dt0);
        }
        self.s.due_f_upto_ = due_f;
        self.s.due_c_upto_ = due_c;
        self.s.due_n_upto_ = due_n;
    }

    /// settleAccumulators()
    pub(crate) fn settle_accumulators(&mut self) {
        let m = self.m;
        if self.s.acc_a_.is_empty() {
            return;
        }
        if !self.s.acc_dirty_ && !self.s.tpa_ {
            return;
        }
        let nc = m.n_cells();
        let nn = m.n_nodes();
        for c in 0..nc {
            let uc = c as usize;
            if self.s.acc_a_[uc] == 0.0 && self.s.acc_q_[uc] == 0.0 {
                continue;
            }
            let inv_dx = 1.0 / m.cell_dx[uc];
            let mut a = self.st.cell_a[uc] + self.s.acc_a_[uc] * inv_dx;
            if a < 0.0 {
                a = 0.0;
            }
            self.st.cell_a[uc] = a;
            self.st.cell_q[uc] += self.s.acc_q_[uc] * inv_dx;
            self.s.acc_a_[uc] = 0.0;
            self.s.acc_q_[uc] = 0.0;
        }
        for n in 0..nn {
            let un = n as usize;
            if self.s.acc_nvol_[un] == 0.0 {
                continue;
            }
            if self.algebraic_active(n) {
                self.s.node_carry_[un] += self.s.acc_nvol_[un];
                self.s.acc_nvol_[un] = 0.0;
                continue;
            }
            let mut vol = self.st.node_volume[un] + self.s.acc_nvol_[un];
            if vol < 0.0 {
                vol = 0.0;
            }
            self.s.acc_nvol_[un] = 0.0;
            self.st.node_volume[un] = vol;
            let d = self.node_depth_from_volume(n, vol);
            self.st.node_head[un] = m.node_invert[un] + d;
        }
        self.s.acc_dirty_ = false;
        // (tpa_: updateTpaFlags; refused before the step)
        self.refresh_depths();
    }

    /// fireFaces(faces, dt0)
    pub(crate) fn fire_faces(&mut self, faces: &[i32], dt0: f64) {
        let m = self.m;
        let n = faces.len();
        if n == 0 {
            return;
        }
        if self.s.live_stamp_.len() != m.n_faces() as usize {
            self.s.live_stamp_.clear();
            self.s.live_stamp_.resize(m.n_faces() as usize, 0);
        }
        self.s.all_faces_live_ = false;
        self.s.live_gen_ = self.s.live_gen_.wrapping_add(1);
        if self.s.live_gen_ == 0 {
            self.s.live_stamp_.iter_mut().for_each(|x| *x = 0);
            self.s.live_gen_ = 1;
        }
        for &f in faces.iter() {
            self.s.live_stamp_[f as usize] = self.s.live_gen_;
        }
        for &f in faces.iter() {
            self.compute_face_flux(f, false);
        }
        self.s.total_flux_ += n as i64;

        // The node solve for the nodes these faces touch (forcing_ is always set here).
        {
            let mut touched = std::mem::take(&mut self.s.ff_touched_);
            let mut touched_nodes = std::mem::take(&mut self.s.ff_touched_nodes_);
            for &node in touched_nodes.iter() {
                touched[node as usize] = 0;
            }
            touched.resize(self.s.acc_nvol_.len(), 0);
            touched_nodes.clear();
            for &f in faces.iter() {
                let nd = m.face_node[f as usize];
                if nd >= 0 && touched[nd as usize] == 0 {
                    touched[nd as usize] = 1;
                    touched_nodes.push(nd);
                }
            }
            self.solve_nodes(&touched_nodes, dt0, true);
            self.s.ff_touched_ = touched;
            self.s.ff_touched_nodes_ = touched_nodes;
        }

        // Positivity over the cycle's outflow (function-static scratch upstream).
        let mut out_cell = std::mem::take(&mut self.s.ff_out_cell_);
        let mut out_node = std::mem::take(&mut self.s.ff_out_node_);
        let mut exporting_cells = std::mem::take(&mut self.s.ff_exporting_cells_);
        let mut exporting_nodes = std::mem::take(&mut self.s.ff_exporting_nodes_);
        let sparse = faces.len() < self.s.acc_a_.len() / 4;
        if sparse && self.s.ff_previous_sparse_ {
            for &c in exporting_cells.iter() {
                out_cell[c] = 0.0;
            }
            for &node in exporting_nodes.iter() {
                out_node[node] = 0.0;
            }
            out_cell.resize(self.s.acc_a_.len(), 0.0);
            out_node.resize(self.s.acc_nvol_.len(), 0.0);
        } else {
            out_cell.clear();
            out_cell.resize(self.s.acc_a_.len(), 0.0);
            out_node.clear();
            out_node.resize(self.s.acc_nvol_.len(), 0.0);
        }
        exporting_cells.clear();
        exporting_nodes.clear();
        self.s.ff_previous_sparse_ = sparse;
        fn record_export(sparse: bool, amounts: &[f64], exporters: &mut Vec<usize>, id: i32, volume: f64) {
            let index = id as usize;
            if sparse && amounts[index] == 0.0 && volume != 0.0 {
                exporters.push(index);
            }
        }
        for &f in faces.iter() {
            let uf = f as usize;
            let fa = self.s.f_mass_[uf];
            if fa == 0.0 {
                continue;
            }
            let vol = fa * ((1i32 << self.s.face_tier_[uf]) as f64 * dt0);
            if vol > 0.0 {
                if m.face_cl[uf] >= 0 {
                    record_export(sparse, &out_cell, &mut exporting_cells, m.face_cl[uf], vol);
                    out_cell[m.face_cl[uf] as usize] += vol;
                } else if m.face_node[uf] >= 0 {
                    record_export(sparse, &out_node, &mut exporting_nodes, m.face_node[uf], vol);
                    out_node[m.face_node[uf] as usize] += vol;
                }
            } else if m.face_cr[uf] >= 0 {
                record_export(sparse, &out_cell, &mut exporting_cells, m.face_cr[uf], vol);
                out_cell[m.face_cr[uf] as usize] -= vol;
            } else if m.face_node[uf] >= 0 {
                record_export(sparse, &out_node, &mut exporting_nodes, m.face_node[uf], vol);
                out_node[m.face_node[uf] as usize] -= vol;
            }
        }
        let n_cells_loop = if sparse { exporting_cells.len() } else { out_cell.len() };
        for i in 0..n_cells_loop {
            let c = if sparse { exporting_cells[i] } else { i };
            if out_cell[c] <= 0.0 {
                out_cell[c] = 1.0;
                continue;
            }
            out_cell[c] = k::positivity_scale(self.st.cell_a[c] * m.cell_dx[c] + self.s.acc_a_[c], out_cell[c], 1.0);
        }
        if m.struct_is_dummy.iter().any(|&v| v != 0) {
            self.refresh_dummy_flows(dt0, Arrived::CellScaled(&out_cell, dt0), false);
        }
        let n_nodes_loop = if sparse { exporting_nodes.len() } else { out_node.len() };
        for i in 0..n_nodes_loop {
            let nn = if sparse { exporting_nodes[i] } else { i };
            if out_node[nn] <= 0.0 {
                out_node[nn] = 1.0;
                continue;
            }
            if self.fc.fixed(nn) {
                out_node[nn] = 1.0;
                continue;
            }
            if self.algebraic_active(nn as i32) {
                let span = (1i32 << self.s.node_tier_[nn]) as f64 * dt0;
                out_node[nn] = if self.s.conservative_node_balance_ {
                    let pending = self.s.acc_nvol_[nn];
                    k::positivity_scale(self.algebraic_available_water(nn, span, pending, &out_cell), out_node[nn], 1.0)
                } else {
                    1.0
                };
                continue;
            }
            out_node[nn] = k::positivity_scale(self.st.node_volume[nn] + self.s.acc_nvol_[nn], out_node[nn], 1.0);
        }

        for &f in faces.iter() {
            let uf = f as usize;
            let cl = m.face_cl[uf];
            let cr = m.face_cr[uf];
            let nd = m.face_node[uf];
            let dt = (1i32 << self.s.face_tier_[uf]) as f64 * dt0;
            let mut fa = self.s.f_mass_[uf];
            if fa != 0.0 && (!sparse || fa * dt != 0.0) {
                let s = if fa > 0.0 {
                    if cl >= 0 {
                        out_cell[cl as usize]
                    } else if nd >= 0 {
                        out_node[nd as usize]
                    } else {
                        1.0
                    }
                } else if cr >= 0 {
                    out_cell[cr as usize]
                } else if nd >= 0 {
                    out_node[nd as usize]
                } else {
                    1.0
                };
                if s < 1.0 {
                    self.s.f_mass_[uf] *= s;
                    self.s.f_mom_[uf] *= s;
                    fa = self.s.f_mass_[uf];
                }
            }
            let fq = self.s.f_mom_[uf];
            if self.s.water_tracking_ {
                self.book_face_water(f, dt);
            }
            self.s.acc_dirty_ = true;
            if cl >= 0 {
                let ul = cl as usize;
                self.s.acc_a_[ul] -= fa * dt;
                self.s.acc_q_[ul] -= m.face_dir_l[uf] as f64 * (fq + self.s.f_corr_l_[uf]) * dt;
                self.s.cell_q_int_[ul] += 0.5 * m.face_dir_l[uf] as f64 * fa * dt;
            }
            if cr >= 0 {
                let ur = cr as usize;
                self.s.acc_a_[ur] += fa * dt;
                self.s.acc_q_[ur] += m.face_dir_r[uf] as f64 * (fq + self.s.f_corr_r_[uf]) * dt;
                self.s.cell_q_int_[ur] += 0.5 * m.face_dir_r[uf] as f64 * fa * dt;
            }
            let vj = m.face_vj_node[uf];
            if vj >= 0 {
                self.book_splice_receipt(vj as usize, fa, dt);
            }
            if nd >= 0 {
                let sign = if cl < 0 { -1.0 } else { 1.0 };
                let un = nd as usize;
                let v = sign * fa * dt;
                self.s.acc_nvol_[un] += v;
                self.s.node_exch_[un] += v;
                if v > 0.0 {
                    self.s.node_in_[un] += v;
                } else {
                    self.s.node_out_[un] -= v;
                }
            }
        }
        self.s.ff_out_cell_ = out_cell;
        self.s.ff_out_node_ = out_node;
        self.s.ff_exporting_cells_ = exporting_cells;
        self.s.ff_exporting_nodes_ = exporting_nodes;
    }

    /// fireCells(cells, dt0, forcing)
    pub(crate) fn fire_cells(&mut self, cells: &[i32], dt0: f64) {
        let m = self.m;
        for &c in cells.iter() {
            let uc = c as usize;
            let dt = (1i32 << self.s.cell_tier_[uc]) as f64 * dt0;
            let g = &m.geom[m.cell_geom[uc] as usize];
            let dx = m.cell_dx[uc];
            let inv_dx = 1.0 / dx;

            let mut dq_src = 0.0;
            if self.s.opts_.order >= 2 && self.s.cell_ho2_[uc] != 0 {
                let half = 0.5 * dx;
                let zp = m.cell_zb[uc] + m.cell_dzdx[uc] * half;
                let zm = m.cell_zb[uc] - m.cell_dzdx[uc] * half;
                let eta_c = self.s.cell_eta_[uc];
                let hp = cmax(0.0, eta_c - zp);
                let hm = cmax(0.0, eta_c - zm);
                let i1p = if hp > 0.0 { k::i1_of_depth(g, hp, k::area_of_depth(g, hp)) } else { 0.0 };
                let i1m = if hm > 0.0 { k::i1_of_depth(g, hm, k::area_of_depth(g, hm)) } else { 0.0 };
                dq_src = k::K_GRAVITY * (i1p - i1m) * dt;
            }

            let mut k_loss = 0.0;
            let faces = [m.cell_face0[uc], m.cell_face1[uc]];
            for &f in faces.iter() {
                let uf = f as usize;
                if m.face_node[uf] >= 0 {
                    k_loss += if m.face_cl[uf] < 0 { m.conduit_loss_inlet[m.cell_conduit[uc] as usize] } else { m.conduit_loss_outlet[m.cell_conduit[uc] as usize] };
                }
            }

            let mut a_new = self.st.cell_a[uc] + self.s.acc_a_[uc] * inv_dx;
            let mut q_new = self.st.cell_q[uc] + (self.s.acc_q_[uc] + dq_src) * inv_dx;
            self.s.acc_a_[uc] = 0.0;
            self.s.acc_q_[uc] = 0.0;

            if !self.s.cell_qlat_.is_empty() && self.s.cell_qlat_[uc] != 0.0 {
                let before_lateral = a_new;
                a_new = self.inject_diverted_lateral(uc, a_new, dt, inv_dx);
                if self.s.water_tracking_ {
                    self.s.water_.cell_lateral[uc] += (cmax(0.0, a_new) - cmax(0.0, before_lateral)) / inv_dx;
                }
            }
            {
                let cr = m.cell_conduit[uc];
                if cr >= 0 {
                    let before_loss = a_new;
                    a_new -= dt * self.fc.conduit_loss[cr as usize];
                    if self.s.water_tracking_ {
                        self.s.water_.cell_loss[uc] += (cmax(0.0, before_loss) - cmax(0.0, a_new)) / inv_dx;
                    }
                }
            }
            if a_new < 0.0 {
                a_new = 0.0;
            }

            let tpa_cell = self.tpa_cell(uc);
            let h_new = if tpa_cell { k::tpa_depth_of_area(g, a_new) } else { k::depth_of_area(g, a_new) };
            self.st.cell_a[uc] = a_new;
            self.st.cell_h[uc] = h_new;
            self.s.cell_eta_[uc] = m.cell_zb[uc] + h_new;
            let radius = self.cache_closure_r(uc, g, h_new);
            if !tpa_cell && h_new <= k::K_DRY_DEPTH {
                self.st.cell_q[uc] = 0.0;
                self.s.cell_u_[uc] = 0.0;
                continue;
            }
            let u = q_new / a_new;
            let cdi = m.cell_conduit[uc] as usize;
            q_new = self.friction_for(g, m.conduit_roughness[cdi], m.conduit_rough_factor[cdi], q_new, u, if tpa_cell { g.y_full } else { h_new }, dt, radius);
            if k_loss > 0.0 {
                q_new = k::local_loss_update(q_new, u, k_loss, dx, dt);
            }
            // (unsteady friction: ufUpdate; refused before the step)
            self.st.cell_q[uc] = q_new;
            self.s.cell_u_[uc] = q_new / a_new;
        }
    }

    /// fireNodes(nodes, dt0, forcing)
    pub(crate) fn fire_nodes(&mut self, nodes: &[i32], dt0: f64) {
        let m = self.m;
        self.refresh_dummy_flows(dt0, Arrived::AccNvol, true);
        if self.s.water_tracking_ {
            self.book_structure_water(dt0);
        }
        self.book_lateral_spills(Some(nodes));

        for &n in nodes.iter() {
            let un = n as usize;
            let dt = (1i32 << self.s.node_tier_[un]) as f64 * dt0;
            let drained = self.s.acc_nvol_[un];
            self.s.acc_nvol_[un] = 0.0;
            self.book_diverted_lateral(un, dt);

            if self.fc.fixed(un) {
                self.st.node_head[un] = self.fc.node_fixed_head[un];
                continue;
            }

            let q_lat = self.node_lateral(un) + self.s.node_qstruct_[un] + self.node_q_dummy(un);

            if self.algebraic_active(n) {
                let bal = self.cap_withdrawal(un, self.s.node_carry_[un] + drained + dt * q_lat, dt);
                self.settle_algebraic_node(n, bal);
                continue;
            }

            let v_prev = self.st.node_volume[un];
            let mut vol = self.cap_withdrawal(un, v_prev + drained + dt * q_lat, dt);
            if vol < 0.0 {
                if self.s.water_tracking_ {
                    self.s.water_.node_deficit[un] -= vol;
                }
                vol = 0.0;
            }
            let mut depth = self.node_depth_from_volume(n, vol);
            self.apply_node_capacity(n, v_prev, &mut vol, &mut depth);
            self.st.node_volume[un] = vol;
            self.st.node_head[un] = m.node_invert[un] + depth;
        }
    }
}

/// refreshDummyFlows' `arrived(node)` callable, one variant per call site.
#[derive(Clone, Copy)]
pub(crate) enum Arrived<'b> {
    /// fireFaces / limitPositivity: the incident face fluxes into the node over dt, each
    /// scaled by its donor cell's outflow scale.
    CellScaled(&'b [f64], f64),
    /// fireNodes: the volume the macro cycle booked into the node (acc_nvol_).
    AccNvol,
    /// updateNodes: the incident face fluxes over dt.
    FaceSum(f64),
}

/// The divergence label from its text (the gate loader's inverse of div_what_).
pub fn div_what_of(b: &[u8]) -> Option<&'static str> {
    match b {
        b"non-finite state" => Some("non-finite state"),
        b"velocity (ft/s)" => Some("velocity (ft/s)"),
        b"depth magnitude (ft)" => Some("depth magnitude (ft)"),
        _ => None,
    }
}
