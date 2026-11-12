// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/SurfaceRouter2D.cpp and SurfaceRouter2D.hpp
// (reportWeight, reportRainfall, reportRainCumulative, reportOldDepth, reportOldHead,
// infilCumulative, outfallExchange, finalize) and src/engine/2d/solver/ExplicitInertialSolver.cpp
// (finalize) (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin,
// licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The 2D router's output side and its end (phase 2a wave 5): the report-time accessors
// SWMMEngine::fillSurfaceSnapshot reads, SurfaceRouter2D::finalize (SWMMEngine::end calls it:
// the partial-batch flush, the solver's run statistics into mass_balance_2d for the "2D Solver
// Statistics" report block, the solver's own finalize).
//
// Differences in form, never in result:
//   * ExplicitInertialSolver::finalize takes the router's mesh, options, boundary, coupling
//     points and the 1D nodes (upstream's solver reaches them through its pointers), as the
//     step does (StepAttachments).
//   * finalize's OPENSWMM_PERF telemetry dump of the subsurface clock is left out: it needs an
//     active two-zone kernel (phase 2b) and only times the run.
// A phase-2b owner finalize would reach (an independent groundwater clock) cannot be active
// here: SurfaceRouter2D::initialize refuses it.

#![allow(clippy::all)]

use crate::cformat;
use crate::core::cmath::cmax;
use crate::core::context::SimulationContext;
use crate::core::engine_init_hydraulics::Terminated;
use crate::core::options::Bytes;
use crate::twod::coupling::node_coupling::CouplingAttachments;
use crate::twod::coupling::outfall_exchange::OutfallExchange;
use crate::twod::data::solver_options_2d::{turbulence_name, SolverOptions2D, Turbulence2D};
use crate::twod::data::surface_state_data::SurfaceStateData;
use crate::twod::solver::explicit_inertial_solver::{eprint_bytes, ExplicitInertialSolver};
use crate::twod::solver::explicit_inertial_solver_step::StepAttachments;
use crate::twod::surface_router_2d::SurfaceRouter2D;

impl SurfaceRouter2D {
    /// reportWeight(report_ms): the clamped interpolation fraction of the report instant inside
    /// the window that crossed it, or 1 without a valid old side or with a nonpositive span.
    pub fn report_weight(&self, report_ms: f64) -> f64 {
        if !self.report_old_valid_ {
            return 1.0;
        }
        let span = self.report_window_new_ms_ - self.report_window_old_ms_;
        if span <= 0.0 {
            return 1.0;
        }
        let mut f = (report_ms - self.report_window_old_ms_) / span;
        if f < 0.0 {
            f = 0.0;
        }
        if f > 1.0 {
            f = 1.0;
        }
        f
    }

    /// reportOldDepth(): the cell depths at the start of the batch that crossed the report
    /// instant.
    pub fn report_old_depth(&self) -> &Vec<f64> {
        &self.report_old_depth_
    }

    /// reportOldHead(): the cell heads at the start of that batch.
    pub fn report_old_head(&self) -> &Vec<f64> {
        &self.report_old_head_
    }

    /// reportRainCumulative(report_ms): the rainfall volume integrated to the report boundary
    /// (exact timestamp match), otherwise the run-time cumulative.
    pub fn report_rain_cumulative(&self, report_ms: f64) -> &Vec<f64> {
        if self.report_rain_ms_ == report_ms {
            &self.report_rain_cum_
        } else {
            &self.rain_cum_
        }
    }

    /// infilCumulative(): the applied infiltration depth per cell.
    pub fn infil_cumulative(&self) -> &Vec<f64> {
        &self.infil_cum_applied_
    }

    /// outfallExchange()
    pub fn outfall_exchange(&self) -> &OutfallExchange {
        &self.outfall_exchange_
    }

    /// reportRainfall(ctx, report_date, out): per-cell rainfall (m/s) as the .out reports
    /// subcatchment rainfall: each gage's series rate at the report date, mapped onto the
    /// cells, then the per-cell rainfall forcings.
    pub fn report_rainfall(&self, ctx: &SimulationContext, report_date: f64, out: &mut Vec<f64>) {
        let n_gages = ctx.n_gages();
        let mut rates = vec![0.0f64; cmax(n_gages, 0) as usize];
        for g in 0..n_gages {
            rates[g as usize] = crate::hydrology::gage_report::get_report_rainfall_from_series(ctx, g, report_date);
        }
        self.gage_rates_to_cells(ctx, &rates, out);
        let st = &self.state_;
        let mut i = 0usize;
        while i < out.len() && i < st.rainfall_forced.len() {
            if st.rainfall_forced[i] == 1 {
                out[i] = st.rainfall_force_val[i];
            } else if st.rainfall_forced[i] == 2 {
                out[i] += st.rainfall_force_val[i];
            }
            i += 1;
        }
    }

    /// finalize(ctx): flush the partial sync batch, publish the solver's cumulative statistics
    /// into ctx.mass_balance_2d, finalize the solver, record the dispersion-limiter advisory,
    /// detach the source owners and mark the router inactive. A C++ exception from the flush
    /// comes back as Err (SWMMEngine::end is noexcept: its caller turns it into a termination).
    pub fn finalize(&mut self, ctx: &mut SimulationContext) -> Result<(), Terminated> {
        // Flush the partial sync batch so the 2D clock ends at the simulation end.
        if self.active_ && self.pending_dt_ > 0.0 {
            let span = self.pending_dt_;
            self.pending_dt_ = 0.0;
            let t = self.last_t_;
            self.co_advance_step(ctx, span, t)?;
        }
        if self.solver_.is_some() {
            // The cumulative solver statistics, before finalize() resets them.
            let s = self.solver_.as_ref().unwrap().run_stats();
            let mb = &mut ctx.mass_balance_2d;
            mb.solver_nsteps = s.nsteps;
            mb.solver_nrhs = s.nrhs;
            mb.solver_last_h = s.last_h;
            mb.solver_avg_h = if s.nsteps > 0 && self.sim_time_ > 0.0 { self.sim_time_ / s.nsteps as f64 } else { 0.0 };
            mb.solver_active_min = s.active_frac_min;
            mb.solver_active_mean = s.active_frac_mean;
            mb.solver_active_max = s.active_frac_max;
            mb.solver_n_tiers = s.n_tiers;
            let mut k = 0;
            while k < s.n_tiers && k < 8 {
                mb.solver_tier_cells[k as usize] = s.tier_cells[k as usize];
                k += 1;
            }
            let mesh_rc = self.mesh_.clone();
            let options_rc = self.options_.clone();
            let boundary_rc = self.boundary_.clone();
            let mesh = mesh_rc.borrow();
            let options = options_rc.borrow();
            let boundary = boundary_rc.borrow();
            let a = StepAttachments {
                mesh: &mesh,
                opts: &options,
                boundary: if self.state_.boundary { Some(&*boundary) } else { None },
                node_coupling: if self.state_.node_coupling { Some(&self.node_coupling_points_[..]) } else { None },
                nodes_1d: if self.state_.nodes_1d { Some(&ctx.nodes) } else { None },
                coupling: CouplingAttachments::of(&self.state_, &self.coupling_head_cells_, &self.node_pond_area_m2_),
            };
            self.solver_.as_mut().unwrap().finalize(&a, &mut self.state_);
        }
        // (perf: the subsurface clock telemetry needs an active two-zone kernel, phase 2b)
        // S2 telemetry: faces whose dispersive exchange was capped.
        if self.state_.transport.dispersion_limiter_binds > 0 {
            let mut buf: Bytes = cformat!(
                b"WARNING: 2D dispersion exchange was limited on %ld face-firing(s) (D\xc2\xb7dt/dx\xc2\xb2 above the explicit limit); the dispersion rate is under-resolved at the marcher's time step there. Reduce DISPERSION or MAX_TIMESTEP, or refine the mesh.",
                self.state_.transport.dispersion_limiter_binds
            );
            buf.truncate(255); // snprintf(buf, 256, ...)
            ctx.warnings.push(buf);
        }
        self.active_ = false;
        self.common_stage_ = false;
        self.source_context_ = false;
        self.joint_hydraulic_ = false;
        self.source_boundary_ready_ = false;
        Ok(())
    }
}

impl ExplicitInertialSolver {
    /// finalize(): the end-of-run notices (wall drag, turbulence), detach the source owners,
    /// settle the pending face receipts and the source ledgers, write the active-set telemetry
    /// file, and release the work arrays and the active-fraction statistics.
    pub fn finalize(&mut self, a: &StepAttachments, st: &mut SurfaceStateData) {
        let opts: Option<&SolverOptions2D> = if self.opts_ { Some(a.opts) } else { None };
        if self.initialized_ && !self.wall_drag_slots_.is_empty() {
            eprint_bytes(&cformat!(
                b"[openswmm 2D] quadratic wall drag: %zu faces, energy removed/rho=%.17g m5/s2\n",
                self.wall_drag_slots_.len(),
                self.wall_drag_energy_
            ));
        }
        if self.initialized_ {
            if let Some(o) = opts {
                if o.turbulence != Turbulence2D::None {
                    eprint_bytes(&cformat!(
                        b"[openswmm 2D] %s mixing: %llu substeps, %llu retries, last diffusion bound %.17g s.\n",
                        turbulence_name(o.turbulence),
                        self.momentum_physics_.substeps,
                        self.momentum_physics_.retries,
                        self.momentum_physics_.last_limit
                    ));
                }
            }
        }
        self.common_stage_ = false;
        self.source_history_ = false;
        self.engine_source_history_ = false;
        self.engine_source_stream_ = false;
        self.checkpoint_ready_ = false;
        self.joint_checkpoint_ = false;
        if self.initialized_ {
            self.settle_accumulators(a, st);
            self.flush_source_ledgers(st);
        }
        if !self.telemetry_path_.is_empty() && !self.telemetry_.is_empty() {
            if let Some(mut f) = crate::core::cprintf::CFile::create(&self.telemetry_path_, false) {
                f.write_bytes(b"t_s,active_cells,active_frac\n");
                let nt = cmax(1, if self.mesh_ { a.mesh.n_triangles() } else { 1 }) as f64;
                for &(t, n) in &self.telemetry_ {
                    f.write_bytes(&cformat!(b"%.3f,%d,%.6f\n", t, n, n as f64 / nt));
                }
                f.flush();
            }
        }
        self.q_.clear();
        self.qcx_.clear();
        self.qcy_.clear();
        self.face_accum_.clear();
        self.dw_slope_.clear();
        self.cell_active_.clear();
        self.active_cells_.clear();
        self.tier_.clear();
        self.face_tier_.clear();
        self.cells_by_tier_.clear();
        self.edges_by_tier_.clear();
        self.bc_cell_.clear();
        self.bc_slot_.clear();
        self.bc_accum_.clear();
        self.bc_q_.clear();
        self.telemetry_.clear();
        self.active_frac_min_ = 1.0e30;
        self.active_frac_max_ = -1.0e30;
        self.active_frac_sum_ = 0.0;
        self.active_samples_ = 0;
        self.initialized_ = false;
    }
}
