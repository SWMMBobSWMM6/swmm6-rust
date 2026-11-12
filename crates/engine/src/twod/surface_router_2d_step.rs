// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/SurfaceRouter2D.cpp and SurfaceRouter2D.hpp (step,
// updateOutfallsPreRouting, updateTrackingPonds, publishEqualizeAreas, refreshDue,
// refreshRenderFields, refreshRenderFieldsIfStale, coAdvanceStep, computeCouplingConductances,
// couplingSyncSpan, prepareLidStep, advancePostRouting, flushPendingBatch, publishNodeRows,
// applyAgingAndReactions, updateRainfall, gageRatesToCells, bookGageRates, resolveBoundaryValues,
// accumulateMassBalance, drainSubcatchRecharge, runtimeActive, nextRuntimeTime, expireRuntime) and
// src/engine/2d/data/SurfaceStateData.hpp (exposedEvaporation, save_state,
// clear_reset_forcings, update_statistics) (HydroCouple/openswmm.engine, branch swmm6_rel,
// commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The 2D router's run-time step on an uncoupled mesh (phase 2a wave 4a): the pre- and
// post-routing hooks SWMMEngine::stepRouting calls, the co-advance batch (rain from the held gage
// rates, evaporation, held infiltration, boundary values, the marcher's advance, the boundary
// and mass-balance ledgers, the statistics and render cadences), and the gage-rate history the
// runoff step books.
//
// The 1D-2D exchange is the join wave's (wave 4b translates it). Every place this step would
// run coupling code refuses with the crate's stub before it changes anything, naming the path:
//   * updateOutfallsPreRouting with outfall or node coupling points, an outfall exchange or
//     tracking ponds (updateOutfallBoundaries, OutfallExchange::recordTailwater / prepare,
//     updateTrackingPonds): the first 2D call of a routing step, so a coupled deck stops there;
//   * advancePostRouting's OutfallExchange::apply, coAdvanceStep's junction ledger (a non-empty
//     last_coupling_exchange), computeCouplingConductances, publishEqualizeAreas with an
//     equalizing point: unreachable once updateOutfallsPreRouting has refused (kept as guards).
// On an uncoupled mesh those calls are no-ops upstream (empty point lists), and the translation
// skips them.
//
// Also left out, each refused at initialize in phase 2a (so the branch cannot run): the common
// source stage, the joint hydraulic pairing, the two-zone subsurface kernel (its producer clock
// setProducerEnd is a member the kernel alone reads), surface species (publishNodeRows and
// applyAgingAndReactions return at once without transport rows), mesh snow and heat
// (prepareMeshHeat and advanceRemainingGroundColumns return at once), the runtime forcing and
// source stores (the swmm_2d_* C API), and the plugin solvers' boundary receipt drain (the
// serial CPU marcher books its receipts itself). The development switches
// OPENSWMM_2D_SYNC_SPAN and OPENSWMM_2D_HEAD_RAMP are refused with the stub when set.
//
// Form: the router's shared members (mesh_, options_, boundary_, infil_) are borrowed for each
// call; the marcher takes them as solver::explicit_inertial_solver_step::StepAttachments.
// Exceptions on their way to SWMMEngine::step's catch are Err(Terminated).

#![allow(clippy::all)]

use crate::core::cmath::{cclamp, cmax};
use crate::core::context::SimulationContext;
use crate::core::datetime;
use crate::core::engine_init_hydraulics::Terminated;
use crate::core::engine_step::stub;
use crate::core::unit_conversion as ucf;
use crate::data::lid_node::LidNodeBottom;
use crate::data::table::{table_lookup_cursor, table_lookup_ex};
use crate::hydrology::infil_bank::infiltration_factors;
use crate::hydrology::lid_node_step as lidnode;
use crate::twod::data::boundary_data::BoundaryType;
use crate::twod::data::mesh_data::MeshData;
use crate::twod::data::solver_options_2d::{Momentum2D, RainfallMode, SolverOptions2D};
use crate::twod::data::surface_state_data::SurfaceStateData;
use crate::twod::mesh::vertex_reconstruction::{reconstruct_vertex_heads, reconstruct_vertex_render_depths};
use crate::twod::solver::explicit_inertial_solver_step::StepAttachments;
use crate::twod::solver::surface_flux_calculator::{compute_cell_continuity, compute_face_velocity, compute_limited_gradients, compute_unlimited_gradients};
use crate::twod::surface_router_2d::SurfaceRouter2D;

/// SurfaceStateData's run-time accessors (SurfaceStateData.hpp).
impl SurfaceStateData {
    /// exposedEvaporation(i): the runtime-resolved evaporation demand over the snow-exposed
    /// surface.
    pub fn exposed_evaporation(&self, i: i32) -> f64 {
        self.runtime_forcings.rate(2, i, self.evap_rate[i as usize]) * if self.snow.owns(i) { self.snow_exposed_fraction(i) } else { 1.0 }
    }

    /// snow.exposedFractions()[i]: snow on the mesh is phase 2b (never owned in 2a).
    fn snow_exposed_fraction(&self, _i: i32) -> f64 {
        1.0
    }

    /// save_state(): the depth and volume snapshot the continuity diagnostic reads.
    pub fn save_state(&mut self) {
        let n = self.depth.len();
        self.old_depth[..n].copy_from_slice(&self.depth);
        let n = self.volume.len();
        self.old_volume[..n].copy_from_slice(&self.volume);
    }

    /// clear_reset_forcings(): expire the one-shot (non-persistent) forcings.
    pub fn clear_reset_forcings(&mut self) {
        if !self.forcing_ever_set {
            return;
        }
        for i in 0..self.rainfall_forced.len() {
            if self.rainfall_persist[i] == 0 && self.rainfall_forced[i] != 0 {
                self.rainfall_forced[i] = 0;
                self.rainfall_force_val[i] = 0.0;
                self.forcing_dirty = true;
            }
            if self.evap_persist[i] == 0 && self.evap_forced[i] != 0 {
                self.evap_forced[i] = 0;
                self.evap_force_val[i] = 0.0;
                self.forcing_dirty = true;
            }
            if self.coupling_persist[i] == 0 && self.coupling_forced[i] != 0 {
                self.coupling_forced[i] = 0;
                self.coupling_force_val[i] = 0.0;
                self.forcing_dirty = true;
            }
        }
    }

    /// update_statistics(tri_area, dt, nthreads): the peak envelopes and the volume-time
    /// integral (each cell writes only its own slots, so the serial loop is the parallel one).
    pub fn update_statistics(&mut self, dt: f64) {
        let n = self.depth.len();
        for i in 0..n {
            if self.depth[i] > self.stat_max_depth[i] {
                self.stat_max_depth[i] = self.depth[i];
            }
            let speed = (self.face_vx[i] * self.face_vx[i] + self.face_vy[i] * self.face_vy[i]).sqrt();
            if speed.is_finite() && (!self.stat_max_velocity[i].is_finite() || speed > self.stat_max_velocity[i]) {
                self.stat_max_velocity[i] = speed;
            }
            let aerr = self.cell_continuity_err[i].abs();
            if aerr > self.stat_max_cont_err[i] {
                self.stat_max_cont_err[i] = aerr;
            }
            self.stat_cum_volume[i] += self.volume[i] * dt;
        }
    }
}

/// refreshOutputGradients(mesh, state, opts): the output gradients when 2D reporting is on.
fn refresh_output_gradients(mesh: &MeshData, state: &mut SurfaceStateData, opts: &SolverOptions2D) {
    if !opts.report_2d {
        return;
    }
    compute_unlimited_gradients(mesh, state, opts.num_threads);
    compute_limited_gradients(mesh, state, opts.limiter_epsilon, opts.num_threads);
}

impl SurfaceRouter2D {
    /// Whether the mesh exchanges water with the 1D network (outfall or node coupling points,
    /// an outfall exchange, tracking ponds): the join wave's paths.
    pub fn coupled_runtime(&self) -> bool {
        !self.coupling_points_.is_empty()
            || !self.node_coupling_points_.is_empty()
            || !self.outfall_exchange_.points_.is_empty()
            || !self.outfall_exchange_.nodes_.is_empty()
            || !self.pond_track_ptr_.is_empty()
    }

    /// runtimeActive(): any runtime surface/aquifer source, boundary override or forcing.
    pub fn runtime_active(&self) -> bool {
        self.state_.runtime_sources.active() || !self.boundary_.borrow().runtime_original.is_empty() || self.state_.runtime_forcings.active()
    }

    /// nextRuntimeTime(now): the earliest runtime source/forcing transition (infinity when there
    /// is none, which is always in phase 2a: only the swmm_2d_* C API installs them).
    pub fn next_runtime_time(&self, _now: f64) -> Result<f64, Terminated> {
        if !self.state_.runtime_sources.rows.is_empty() || !self.state_.runtime_forcings.rows.is_empty() {
            return Err(stub("the 2D runtime sources and forcings (SurfaceRouter2D::nextRuntimeTime; the swmm_2d_* C API)"));
        }
        Ok(f64::INFINITY)
    }

    /// expireRuntime(now): retire the runtime sources and forcings (none in phase 2a).
    pub fn expire_runtime(&mut self, _now: f64) -> Result<(), Terminated> {
        if !self.state_.runtime_sources.rows.is_empty() || !self.state_.runtime_forcings.rows.is_empty() {
            return Err(stub("the 2D runtime sources and forcings (SurfaceRouter2D::expireRuntime; the swmm_2d_* C API)"));
        }
        Ok(())
    }

    /// drainSubcatchRecharge(out): move the pending SUBCATCH_AQUIFER recharge into `out` (m3 per
    /// subcatchment) and zero it.
    pub fn drain_subcatch_recharge(&mut self, out: &mut Vec<f64>) {
        out.clone_from(&self.subcatch_recharge_);
        for x in self.subcatch_recharge_.iter_mut() {
            *x = 0.0;
        }
    }

    /// step(ctx, dt, t): the two hooks around the 1D routing.
    pub fn step(&mut self, ctx: &mut SimulationContext, dt: f64, t: f64) -> Result<(), Terminated> {
        if !self.active_ {
            return Ok(());
        }
        self.update_outfalls_pre_routing(ctx, dt)?;
        self.advance_post_routing(ctx, dt, t)
    }

    /// updateOutfallsPreRouting(ctx, dt): freeze the outfall tailwater and donor reservations
    /// before the 1D routing solve (a coupled mesh only; see the header).
    pub fn update_outfalls_pre_routing(&mut self, _ctx: &mut SimulationContext, _dt: f64) -> Result<(), Terminated> {
        if !self.active_ {
            return Ok(());
        }
        if self.coupled_runtime() {
            return Err(stub(
                "the 1D-2D exchange at run time (SurfaceRouter2D::updateOutfallsPreRouting: updateOutfallBoundaries, OutfallExchange::recordTailwater/prepare, updateTrackingPonds; join wave)",
            ));
        }
        // Uncoupled: updateOutfallBoundaries, recordTailwater, prepare and updateTrackingPonds
        // have no point, node or pond to visit.
        Ok(())
    }

    /// computeCouplingConductances(ctx, out): the node exchange's head sensitivity for the DW
    /// node continuity (nothing on an uncoupled mesh).
    pub fn compute_coupling_conductances(&self, _ctx: &SimulationContext, _out: &mut Vec<(i32, f64)>) -> Result<(), Terminated> {
        if !self.active_ {
            return Ok(());
        }
        if !self.node_coupling_points_.is_empty() {
            return Err(stub("the 1D-2D node exchange conductance (SurfaceRouter2D::computeCouplingConductances; join wave)"));
        }
        Ok(())
    }

    /// couplingSyncSpan(ctx, options): the co-advance batch span.
    fn coupling_sync_span(ctx: &SimulationContext, options: &SolverOptions2D) -> Result<f64, Terminated> {
        if std::env::var_os("OPENSWMM_2D_SYNC_SPAN").is_some() {
            return Err(stub("the OPENSWMM_2D_SYNC_SPAN development switch (SurfaceRouter2D::couplingSyncSpan)"));
        }
        Ok(if options.coupling_sync > 0.0 { cclamp(options.coupling_sync, ctx.options.routing_step, 60.0) } else { ctx.options.routing_step })
    }

    /// prepareLidStep(ctx, dt, evaporation): lidnode::prepareMeshRichardsStep, which is
    /// lidnode::prepareStep unless a Richards node has a MESH_AQUIFER bottom (phase 2b: such a
    /// bottom needs the two-zone kernel).
    pub fn prepare_lid_step(&mut self, ctx: &mut SimulationContext, dt: f64, evaporation: f64) -> Result<(), Terminated> {
        for n in 0..ctx.n_nodes() {
            if lidnode::richards_mode(ctx, n) && lidnode::bottom_mode(ctx, n) == LidNodeBottom::MeshAquifer {
                return Err(stub("a Richards LID node on the 2D mesh aquifer (lidnode::prepareMeshRichardsStep; phase 2b)"));
            }
        }
        lidnode::prepare_step(ctx, dt, evaporation).map_err(Terminated::from)?;
        // subsurface_.lidReceiptsEnabled() is false without the two-zone kernel.
        Ok(())
    }

    /// advancePostRouting(ctx, routing_dt, t): consume the routing step's exchange and run the
    /// pending batch when the sync span is reached.
    pub fn advance_post_routing(&mut self, ctx: &mut SimulationContext, routing_dt: f64, t: f64) -> Result<(), Terminated> {
        if !self.active_ {
            return Ok(());
        }
        // publishNodeRows(ctx): no transport rows. OutfallExchange::apply: no outfall point.
        if self.coupled_runtime() {
            return Err(stub("the 1D-2D exchange at run time (SurfaceRouter2D::advancePostRouting: OutfallExchange::apply; join wave)"));
        }
        self.pending_dt_ += routing_dt;
        // subsurface_.setProducerEnd(sim_time_ + pending_dt_): the kernel's own clock (phase 2b).
        self.last_t_ = t;
        // aquifer_link_shares_ is empty without the two-zone kernel.
        let sync = {
            let options = self.options_.borrow();
            Self::coupling_sync_span(ctx, &options)?
        };
        if self.pending_dt_ + 0.5 * routing_dt >= sync {
            let span = self.pending_dt_;
            self.pending_dt_ = 0.0;
            self.co_advance_step(ctx, span, t)?;
        }
        Ok(())
    }

    /// flushPendingBatch(ctx): run what is pending of the batch now.
    pub fn flush_pending_batch(&mut self, ctx: &mut SimulationContext) -> Result<(), Terminated> {
        if !self.active_ {
            return Ok(());
        }
        if self.pending_dt_ > 0.0 {
            let span = self.pending_dt_;
            self.pending_dt_ = 0.0;
            let t = self.last_t_;
            self.co_advance_step(ctx, span, t)?;
        }
        Ok(())
    }

    /// refreshDue(ctx, dt): the 30 s statistics cadence.
    fn refresh_due(&self, dt: f64) -> bool {
        let refresh_dt = cmax(dt, 30.0);
        self.co_refresh_elapsed_ + dt >= refresh_dt
    }

    /// refreshRenderFields(): face velocity (or FULL_SWE's cell momentum velocity), vertex
    /// heads, output gradients and vertex render depths.
    pub fn refresh_render_fields(&mut self) {
        let mesh_rc = self.mesh_.clone();
        let options_rc = self.options_.clone();
        let boundary_rc = self.boundary_.clone();
        let mesh = mesh_rc.borrow();
        let options = options_rc.borrow();
        let boundary = boundary_rc.borrow();
        if options.momentum == Momentum2D::FullSwe && self.solver_.is_some() {
            self.solver_.as_ref().unwrap().publish_momentum_velocity(&mesh, &mut self.state_, &options);
        } else {
            let b = if self.state_.boundary { Some(&*boundary) } else { None };
            compute_face_velocity(&mesh, &mut self.state_, b, &options);
        }
        reconstruct_vertex_heads(&mesh, &mut self.state_, options.num_threads);
        refresh_output_gradients(&mesh, &mut self.state_, &options);
        reconstruct_vertex_render_depths(&mesh, &mut self.state_, options.dry_depth, options.num_threads, &mut self.render_eta_);
        self.co_render_elapsed_ = 0.0;
        self.render_dirty_ = false;
    }

    /// refreshRenderFieldsIfStale().
    pub fn refresh_render_fields_if_stale(&mut self) {
        if self.render_dirty_ {
            self.refresh_render_fields();
        }
    }

    /// publishEqualizeAreas(ctx, window): no equalizing point on an uncoupled mesh, so only the
    /// window is published.
    fn publish_equalize_areas(&mut self, window: f64) -> Result<(), Terminated> {
        self.state_.node_pond_area_m2 = false;
        self.state_.coupling_window_s = window;
        if self.node_coupling_points_.iter().any(|cp| cp.equalize) {
            return Err(stub("the 1D-2D equalization areas (SurfaceRouter2D::publishEqualizeAreas; join wave)"));
        }
        Ok(())
    }

    /// updateRainfall(ctx): the held per-cell rainfall from this window's gage rates.
    fn update_rainfall(&mut self, ctx: &SimulationContext) {
        let n_gages = cmax(ctx.n_gages(), 0) as usize;
        let rates: Vec<f64> = if self.gage_window_rate_.len() == n_gages { self.gage_window_rate_[..n_gages].to_vec() } else { ctx.gages.rainfall[..n_gages].to_vec() };
        self.rain_gage_last_ = rates;
        let mut out = std::mem::take(&mut self.state_.rainfall);
        let last = std::mem::take(&mut self.rain_gage_last_);
        self.gage_rates_to_cells(ctx, &last, &mut out);
        self.rain_gage_last_ = last;
        self.state_.rainfall = out;
    }

    /// gageRatesToCells(ctx, rates, out): map native gage rates (in/hr or mm/hr) onto the
    /// cells in SI m/s (spatial weights, or the gage mean).
    pub fn gage_rates_to_cells(&self, ctx: &SimulationContext, rates: &[f64], out: &mut Vec<f64>) {
        let n_gages = ctx.n_gages();
        let options = self.options_.borrow();
        let nt = self.mesh_.borrow().n_triangles() as usize;
        out.resize(nt, 0.0);
        if n_gages <= 0 || options.rainfall_mode == RainfallMode::None {
            for x in out.iter_mut() {
                *x = 0.0;
            }
            return;
        }
        let to_ms = if ucf::get_unit_system(ctx.options.flow_units as i32) == 0 { 0.0254 / 3600.0 } else { 0.001 / 3600.0 };
        let mut rain_si = vec![0.0f64; n_gages as usize];
        for g in 0..n_gages as usize {
            rain_si[g] = rates[g] * to_ms;
        }
        if (options.rainfall_mode == RainfallMode::NaturalNeighbour || options.rainfall_mode == RainfallMode::NearestNeighbour) && self.interp_.ready() {
            self.interp_.apply(&rain_si, out);
        } else {
            let mut mean = 0.0;
            for &r in &rain_si {
                mean += r;
            }
            mean /= n_gages as f64;
            for x in out.iter_mut() {
                *x = mean;
            }
        }
    }

    /// bookGageRates(ctx, runoff_start_s): integrate the held gage rates to this runoff
    /// boundary and install the new ones.
    pub fn book_gage_rates(&mut self, ctx: &SimulationContext, runoff_start_s: f64) {
        let n_gages = cmax(ctx.n_gages(), 0) as usize;
        if self.gage_rate_.len() != n_gages {
            self.gage_rate_.clear();
            self.gage_rate_.resize(n_gages, 0.0);
            self.gage_cum_.clear();
            self.gage_cum_.resize(n_gages, 0.0);
            self.gage_window_depth_.clear();
            self.gage_window_depth_.resize(n_gages, 0.0);
            self.gage_report_depth_.clear();
            self.gage_report_depth_.resize(n_gages, 0.0);
            self.gage_since_ = 0.0;
        }
        let report_s = ctx.next_report_ms / 1000.0;
        if self.gage_since_ <= report_s && runoff_start_s > report_s {
            for g in 0..n_gages {
                self.gage_report_depth_[g] = self.gage_cum_[g] + self.gage_rate_[g] * (report_s - self.gage_since_);
            }
        }
        for g in 0..n_gages {
            self.gage_cum_[g] += self.gage_rate_[g] * cmax(0.0, runoff_start_s - self.gage_since_);
            self.gage_rate_[g] = ctx.gages.rainfall[g];
        }
        self.gage_since_ = runoff_start_s;
    }

    /// resolveBoundaryValues(ctx, t): the non-WALL slots' stage and flow from their series and
    /// rating curves at the held routing sample t.
    fn resolve_boundary_values(&mut self, ctx: &mut SimulationContext, t: f64) {
        let boundary_rc = self.boundary_.clone();
        let mut boundary = boundary_rc.borrow_mut();
        let ne = boundary.size();
        if ne == 0 {
            return;
        }
        if self.bc_nonwall_dirty_ {
            self.bc_nonwall_dirty_ = false;
            self.bc_nonwall_slots_.clear();
            for idx in 0..ne {
                if boundary.edge_bc_type[idx as usize] != BoundaryType::Wall as i8 {
                    self.bc_nonwall_slots_.push(idx);
                }
            }
        }
        if !self.boundary_names_resolved_ {
            self.boundary_names_resolved_ = true;
            for idx in 0..ne as usize {
                if boundary.edge_bc_tseries[idx] == -2 {
                    boundary.edge_bc_tseries[idx] = ctx.find_timeseries(&boundary.edge_bc_tseries_name[idx]);
                }
                if boundary.edge_bc_flow_tseries[idx] == -2 {
                    boundary.edge_bc_flow_tseries[idx] = ctx.find_timeseries(&boundary.edge_bc_flow_tseries_name[idx]);
                }
                if boundary.edge_bc_rating_curve[idx] == -2 {
                    boundary.edge_bc_rating_curve[idx] = ctx.find_curve(&boundary.edge_bc_rating_curve_name[idx]);
                }
            }
        }
        // (the coupled RK2 owner's early return needs the two-zone kernel)
        let abs_t = datetime::add_seconds(ctx.options.start_date, t);
        let n_tables = ctx.tables.tables.len() as i32;
        for k in 0..self.bc_nonwall_slots_.len() {
            let idx = self.bc_nonwall_slots_[k] as usize;
            let bt = boundary.edge_bc_type[idx];
            if bt == BoundaryType::SpecifiedStage as i8 {
                let ts = boundary.edge_bc_tseries[idx];
                if ts >= 0 && ts < n_tables {
                    boundary.edge_bc_head[idx] = self.bc_stage_ts_scale_ * table_lookup_cursor(&mut ctx.tables.tables[ts as usize], abs_t);
                }
            } else if bt == BoundaryType::SpecifiedFlow as i8 {
                let ts = boundary.edge_bc_flow_tseries[idx];
                if ts >= 0 && ts < n_tables {
                    boundary.edge_bc_flow[idx] = self.bc_flow_scale_ * table_lookup_cursor(&mut ctx.tables.tables[ts as usize], abs_t);
                }
            } else if bt == BoundaryType::RatingCurve as i8 {
                let cv = boundary.edge_bc_rating_curve[idx];
                if cv >= 0 && cv < n_tables {
                    let i = MeshData::slot_cell(idx as i32) as usize;
                    boundary.edge_bc_flow[idx] = self.bc_flow_scale_ * table_lookup_ex(&ctx.tables.tables[cv as usize], self.state_.head[i] / self.bc_stage_ts_scale_);
                }
            }
        }
    }

    /// accumulateMassBalance(ctx, dt): book one batch's accepted receipts into the 2D ledger.
    fn accumulate_mass_balance(&mut self, ctx: &mut SimulationContext, dt: f64) {
        let mesh_rc = self.mesh_.clone();
        let mesh = mesh_rc.borrow();
        let nt = mesh.n_triangles() as usize;
        let do_infil = self.infil_.borrow().active();
        let track_cum = do_infil && self.infil_cum_applied_.len() == nt;
        let mut rain_vol = 0.0;
        let mut infil_vol = 0.0;
        let mut storage = 0.0;
        let route_aq = do_infil && self.cell_subcatch_.len() == nt && !self.subcatch_recharge_.is_empty();
        for ui in 0..nt {
            let i = ui as i32;
            let area = mesh.tri_area[ui];
            let weather_area = area;
            let rain_rate = self.state_.rainfall[ui];
            rain_vol += (if self.state_.snow.owns(i) { self.state_.liquid_rainfall(i) } else { rain_rate }) * weather_area;
            self.rain_cum_[ui] += rain_rate * weather_area * dt;
            if do_infil {
                let applied = self.state_.infil_applied[ui];
                self.state_.infil_applied[ui] = 0.0;
                infil_vol += applied * area;
                if track_cum {
                    self.infil_cum_applied_[ui] += applied;
                }
                if route_aq && self.cell_subcatch_[ui] >= 0 {
                    let us = self.cell_subcatch_[ui] as usize;
                    if us < self.subcatch_recharge_.len() {
                        self.subcatch_recharge_[us] += applied * area;
                    }
                }
            }
            if self.coupling_cum_applied_.len() == nt {
                self.coupling_cum_applied_[ui] += self.state_.coupling_applied[ui];
                self.state_.coupling_applied[ui] = 0.0;
            }
            storage += self.state_.volume[ui];
        }
        let mb = &mut ctx.mass_balance_2d;
        mb.rainfall_in += rain_vol * dt;
        mb.evap_out = self.state_.evap_loss_total;
        if do_infil {
            mb.infil_out += infil_vol;
        }
        if route_aq {
            let mut pending = 0.0;
            for &v in &self.subcatch_recharge_ {
                pending += v;
            }
            mb.infil_aquifer_pending = pending;
        }
        // (the subsurface's return, aquifer_in, needs the two-zone kernel)
        self.mb_seen_nodes_.clear();
        let vol_1d_to_2d = self.options_.borrow().vol_1d_to_2d;
        for k in 0..self.coupling_points_.len() {
            let cp = &self.coupling_points_[k];
            if cp.is_outfall {
                continue;
            }
            if !self.mb_seen_nodes_.insert(cp.node_idx) {
                continue;
            }
            let vol = ctx.nodes.coupling_volume[cp.node_idx as usize] * vol_1d_to_2d;
            let mb = &mut ctx.mass_balance_2d;
            if vol > 0.0 {
                mb.coupling_2d_to_1d_out += vol;
            } else {
                mb.coupling_1d_to_2d_in += -vol;
            }
        }
        let mut incoming = 0.0;
        let mut outgoing = 0.0;
        {
            let boundary = self.boundary_.borrow();
            let receipts = boundary.receipts.borrow();
            for (_slot, r) in receipts.iter() {
                incoming += r.water_in;
                outgoing += r.water_out;
            }
        }
        let mb = &mut ctx.mass_balance_2d;
        mb.boundary_in += incoming - self.prev_boundary_in_;
        mb.boundary_out += outgoing - self.prev_boundary_out_;
        self.prev_boundary_in_ = incoming;
        self.prev_boundary_out_ = outgoing;
        mb.external_in = 0.0;
        mb.external_out = 0.0;
        for row in &self.state_.runtime_sources.rows {
            mb.external_in += row.receipt.water_in;
            mb.external_out += row.receipt.water_out;
        }
        mb.final_storage = storage;
    }

    /// coAdvanceStep(ctx, dt, t): execute one pending surface batch.
    pub fn co_advance_step(&mut self, ctx: &mut SimulationContext, dt: f64, t: f64) -> Result<(), Terminated> {
        if dt <= 0.0 {
            return Ok(());
        }
        // Report-time blend.
        {
            let window_end_ms = self.report_window_new_ms_ + 1000.0 * dt;
            if window_end_ms >= ctx.next_report_ms {
                let n = self.state_.depth.len();
                self.report_old_depth_[..n].copy_from_slice(&self.state_.depth);
                let n = self.state_.head.len();
                self.report_old_head_[..n].copy_from_slice(&self.state_.head);
                self.report_window_old_ms_ = self.report_window_new_ms_;
                self.report_old_valid_ = true;
            }
            self.report_window_new_ms_ = window_end_ms;
        }
        let refresh_due = self.refresh_due(dt);
        if refresh_due {
            self.state_.save_state();
        }
        if !self.state_.coupling_native.is_empty() {
            self.state_.coupling_flux = self.state_.coupling_native.clone();
        }
        if self.outfall_state_changed_ {
            // The CPU marcher updated its closures on the affected cells.
            self.outfall_state_changed_ = false;
        }
        // publishNodeRows(ctx): no transport rows.

        // Rainfall / forcings on a coarse cadence.
        self.co_forcing_elapsed_ += dt;
        let n_gages = cmax(ctx.n_gages(), 0) as usize;
        let report_s = ctx.next_report_ms / 1000.0;
        let rain_report_due = self.report_window_new_ms_ >= ctx.next_report_ms && report_s >= self.sim_time_ - 1e-8;
        let mut report_gage_depth: Vec<f64> = Vec::new();
        if rain_report_due {
            report_gage_depth = vec![0.0; n_gages];
        }
        if self.gage_rate_.len() == n_gages && n_gages > 0 {
            let t0 = self.sim_time_;
            let t1 = self.sim_time_ + dt;
            self.gage_window_rate_.resize(n_gages, 0.0);
            for g in 0..n_gages {
                let d1 = self.gage_cum_[g] + self.gage_rate_[g] * cmax(0.0, t1 - self.gage_since_);
                self.gage_window_rate_[g] = if self.gage_since_ <= t0 { self.gage_rate_[g] } else { (d1 - self.gage_window_depth_[g]) / dt };
                if rain_report_due {
                    let d_report = if self.gage_since_ <= report_s {
                        self.gage_cum_[g] + self.gage_rate_[g] * (report_s - self.gage_since_)
                    } else {
                        self.gage_report_depth_[g]
                    };
                    report_gage_depth[g] = d_report - self.gage_window_depth_[g];
                }
                self.gage_window_depth_[g] = d1;
            }
        }
        let gage_changed = {
            let rates: &[f64] = if self.gage_window_rate_.len() == n_gages { &self.gage_window_rate_ } else { &ctx.gages.rainfall };
            let mut changed = self.rain_gage_last_.len() != n_gages;
            let mut g = 0;
            while !changed && g < n_gages {
                changed = rates[g] != self.rain_gage_last_[g];
                g += 1;
            }
            changed
        };
        let (evaporation, momentum, report_step, dry_depth, num_threads) = {
            let o = self.options_.borrow();
            (o.evaporation, o.momentum, ctx.options.report_step, o.dry_depth, o.num_threads)
        };
        if self.co_forcing_elapsed_ >= 30.0 || self.co_forcing_first_ || self.state_.forcing_dirty || gage_changed {
            self.state_.forcing_dirty = false;
            self.update_rainfall(ctx);
            let evap_base = if evaporation == 2 { cmax(0.0, ctx.climate_state.evap_rate) * 0.3048 } else { 0.0 };
            for x in self.state_.evap_rate.iter_mut() {
                *x = evap_base;
            }
            let evap_on = evaporation != 0;
            let st = &mut self.state_;
            for i in 0..st.depth.len() {
                if st.rainfall_forced[i] == 1 {
                    st.rainfall[i] = st.rainfall_force_val[i];
                } else if st.rainfall_forced[i] == 2 {
                    st.rainfall[i] += st.rainfall_force_val[i];
                }
                if !evap_on {
                    st.evap_rate[i] = 0.0;
                } else if st.evap_forced[i] == 1 {
                    st.evap_rate[i] = st.evap_force_val[i];
                } else if st.evap_forced[i] == 2 {
                    st.evap_rate[i] += st.evap_force_val[i];
                }
            }
            self.co_forcing_elapsed_ = 0.0;
            self.co_forcing_first_ = false;
        }
        // (mesh snow is a phase-2b refusal)
        if self.state_.forcing_ever_set {
            let st = &mut self.state_;
            st.coupling_native = st.coupling_flux.clone();
            for i in 0..st.coupling_flux.len() {
                if st.coupling_forced[i] == 1 {
                    st.coupling_flux[i] = st.coupling_force_val[i];
                } else if st.coupling_forced[i] == 2 {
                    st.coupling_flux[i] += st.coupling_force_val[i];
                }
            }
        }
        if rain_report_due {
            let mut cum = std::mem::take(&mut self.report_rain_cum_);
            self.gage_rates_to_cells(ctx, &report_gage_depth, &mut cum);
            self.report_rain_cum_ = cum;
            let span = cclamp(report_s - self.sim_time_, 0.0, dt);
            let mesh = self.mesh_.borrow();
            for i in 0..self.report_rain_cum_.len() {
                if self.state_.rainfall_forced[i] == 1 {
                    self.report_rain_cum_[i] = self.state_.rainfall_force_val[i] * span;
                } else if self.state_.rainfall_forced[i] == 2 {
                    self.report_rain_cum_[i] += self.state_.rainfall_force_val[i] * span;
                }
                let area = mesh.tri_area[i];
                self.report_rain_cum_[i] = self.rain_cum_[i] + self.report_rain_cum_[i] * area;
            }
            self.report_rain_ms_ = ctx.next_report_ms;
        }

        // Held infiltration on its own INFIL_STEP cadence.
        {
            let infil_rc = self.infil_.clone();
            let mut infil = infil_rc.borrow_mut();
            if infil.active() {
                self.infil_elapsed_ += dt;
                if self.infil_elapsed_ >= infil.step_seconds() {
                    let mesh = self.mesh_.borrow();
                    infil.update_rates(&mesh, &mut self.state_, self.infil_elapsed_, infiltration_factors(ctx, t), num_threads);
                    self.infil_elapsed_ = 0.0;
                }
            }
        }

        self.resolve_boundary_values(ctx, t);

        if std::env::var_os("OPENSWMM_2D_HEAD_RAMP").is_some() {
            return Err(stub("the OPENSWMM_2D_HEAD_RAMP development switch (SurfaceRouter2D::coAdvanceStep)"));
        }
        // prepareMeshHeat: mesh heat is a phase-2b refusal (returns at once when disabled).
        self.publish_equalize_areas(dt)?;
        {
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
                node_coupling_len: if self.state_.node_coupling { self.node_coupling_points_.len() } else { 0 },
            };
            let t0 = self.sim_time_;
            if let Some(solver) = self.solver_.as_mut() {
                solver.advance(&a, &mut self.state_, t0, t0 + dt)?;
            }
        }
        // (the plugin solvers' boundary receipt drain, the porewater ledger, the ground columns
        // and the ground heat checks belong to owners phase 2a does not have)
        self.sim_time_ += dt;
        // applyAgingAndReactions(ctx, dt): no transport rows.

        // Junction ledger (a coupled mesh only).
        let exch_nonempty = self.solver_.as_ref().map_or(false, |s| !s.exch_.is_empty());
        if exch_nonempty {
            return Err(stub("the 1D-2D junction ledger (SurfaceRouter2D::coAdvanceStep: last_coupling_exchange; join wave)"));
        }

        // Boundary ledger.
        {
            let boundary_rc = self.boundary_.clone();
            let mut boundary = boundary_rc.borrow_mut();
            for k in 0..self.bc_nonwall_slots_.len() {
                let idx = self.bc_nonwall_slots_[k] as usize;
                boundary.edge_bc_cum_flux[idx] += -self.state_.edge_flux[idx] * dt;
            }
        }

        self.accumulate_mass_balance(ctx, dt);
        // (no aquifer node exchange; no junction volumes to queue)
        self.state_.clear_reset_forcings();
        self.reset_window_accumulators();

        self.co_refresh_elapsed_ += dt;
        if refresh_due {
            let mesh_rc = self.mesh_.clone();
            let options_rc = self.options_.clone();
            let boundary_rc = self.boundary_.clone();
            let mesh = mesh_rc.borrow();
            let options = options_rc.borrow();
            compute_cell_continuity(&mesh, &mut self.state_, &options, dt);
            if momentum == Momentum2D::FullSwe && self.solver_.is_some() {
                self.solver_.as_ref().unwrap().publish_momentum_velocity(&mesh, &mut self.state_, &options);
            } else {
                let boundary = boundary_rc.borrow();
                let b = if self.state_.boundary { Some(&*boundary) } else { None };
                compute_face_velocity(&mesh, &mut self.state_, b, &options);
            }
            self.state_.update_statistics(self.co_refresh_elapsed_);
            self.co_refresh_elapsed_ = 0.0;
        }

        self.render_dirty_ = true;
        self.co_render_elapsed_ += dt;
        let render_dt = cmax(dt, if report_step > 0.0 { report_step } else { 30.0 });
        if self.co_render_elapsed_ >= render_dt {
            self.refresh_render_fields();
        }
        let _ = dry_depth;
        Ok(())
    }
}
