// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/SWMMEngine.cpp (SWMMEngine::step, advanceTo,
// runtimeCouplingTime, emit_progress, applyForcings, resetStepMassBalance,
// accumulateNodeRoutingTotals, checkRoutingDiverged, updateStatistics,
// updateRoutingMassBalance, computeFinalStorage, assembleLateralInflows, isBetweenEvents,
// isInSteadyState, ensureXspCache, and the file-local reassertNodeStateOverrides and
// effectiveUserLatFlow), hydraulics/Inlet.cpp (InletSolver::adjustFloodingTotals) and
// transport/components/ReactionModule/MsxLegacyTransport.cpp (bookMsxSystemOutflow), and
// SWMMEngine.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// One routing step (1D build: OPENSWMM_HAS_2D off) and the step's own bookkeeping around the
// solvers. The solver entry points are other files' (see PLAN.md for the signatures):
//   core/engine_step_runoff.rs   step_runoff
//   core/engine_step_routing.rs  step_routing, Router::get_adaptive_step,
//                                Router::accumulate_link_step_stats
//   core/engine_output.rs        post_output_snapshot, accumulate_avg_results,
//                                write_due_hotstart_saves
// assemble_lateral_inflows, is_between_events and is_in_steady_state are stepRouting's
// helpers; they are here, pub(crate), for it to call.
//
// Errors and exceptions. step() is noexcept upstream with a try/catch(std::exception&)
// around the step body. In the translation:
//   * Err(Terminated) from a callee is a termination (an exception escaping a noexcept
//     function: std::terminate) and step() returns it as it is; or, with exception
//     STUB_EXCEPTION, a call into code not translated yet, which the gate harness reports.
//   * A C++ exception still propagating (thrown by a function that is not noexcept, on its
//     way to step()'s catch) is an Err(Terminated) whose exception starts with
//     THROWN_PREFIX; step() catches it as upstream's catch does: checkpoint_ready_ = false,
//     set_error(SWMM_ERR_WRONG_STATE, what), return SWMM_ERR_WRONG_STATE.
//
// Differences of form, never in result:
//   * The StepGuard (step_active_ and runtime_pipeline_ cleared on every return) is the
//     wrapper step() around step_body().
//   * Engine checkpointing (checkpoint_enabled_) is only switched on by untranslated code
//     (prepareNativeHotstartCheckpointing, the checkpoint C API), so its branches are
//     Err(stub) here.
//   * xsp_cache_ keeps only its size (see engine.rs); getVelocity's XSectParams are built
//     from the link store at the use, a verbatim copy of what the cache holds.
//   * The callbacks (progress, step begin and end) are boxed closures; the C callback's
//     engine pointer and user data are their captures.
//   * assembleLateralInflows is assemble_lateral_inflows_sources (every source) followed by
//     inlet_.compute_all, so the gate can run the first part without the routing step's.
//   * The OPENSWMM_FV_SLOT_TRACE environment trace in updateStatistics prints to stdout
//     only; it is left out.

use crate::core::cmath::{cmax, cmin};
use crate::core::context::{EngineState, SimulationContext};
use crate::core::datetime;
use crate::core::engine::{SwmmEngine, CFFI_ERR_BADPARAM, CFFI_ERR_LIFECYCLE, SWMM_ERR_WRONG_STATE, SWMM_OK};
use crate::core::engine_init_hydraulics::Terminated;
use crate::core::engine_init_quality::snow_cover_volume_ft3;
use crate::core::options::{Bytes, RoutingModel};
use crate::data::forcing::ForcingMode;
use crate::data::link::{FlowClass, LinkData, LinkType, XsectShape};
use crate::data::lid_node::LidNodeBottom;
use crate::data::node::NodeType;
use crate::data::pollutant::MassUnits;
use crate::hydraulics::inlet::InletSolver;
use crate::hydraulics::{link, timestep_controller as tc};
use crate::hydrology::lid_node_held as lidnode;

/// The `exception` of a Terminated that marks a call into code not translated yet (a stub
/// file or an untranslated branch), so the gate harness can report where the run stopped.
pub const STUB_EXCEPTION: &str = "swmm6_rust::stub";

/// A Terminated meaning "the translation stops here: `what` is not translated yet".
pub fn stub(what: &'static str) -> Terminated {
    Terminated { exception: STUB_EXCEPTION, what }
}

/// The prefix of a Terminated's `exception` for a C++ exception still propagating to a
/// catch (see the header): "throw std::runtime_error", for instance.
pub const THROWN_PREFIX: &str = "throw ";

/// CffiErrorCode CFFI_ERR_NUMERICAL (core/ErrorCodes.hpp)
pub const CFFI_ERR_NUMERICAL: i32 = 14;

/// reassertNodeStateOverrides: republish node temperature and age overrides.
fn reassert_node_state_overrides(ctx: &mut SimulationContext) {
    let nn = ctx.n_nodes() as usize;
    if ctx.options.heat_transport && ctx.heat_state.node_temp.len() == nn {
        let mut u = 0;
        while u < nn && u < ctx.forcing.node_temperature_mode.len() {
            if ctx.forcing.node_temperature_mode[u] == ForcingMode::Override {
                ctx.heat_state.node_temp[u] = ctx.forcing.node_temperature_value[u];
            }
            u += 1;
        }
    }
    if ctx.options.water_age && ctx.water_age_state.node_age.len() == nn {
        let mut u = 0;
        while u < nn && u < ctx.forcing.node_age_mode.len() {
            if ctx.forcing.node_age_mode[u] == ForcingMode::Override {
                ctx.water_age_state.node_age[u] = ctx.forcing.node_age_value[u] * 3600.0;
            }
            u += 1;
        }
    }
}

/// effectiveUserLatFlow: the persistent user lateral flow with the forcing overlay.
fn effective_user_lat_flow(ctx: &SimulationContext, uj: usize) -> f64 {
    let mut q = ctx.nodes.user_lat_flow[uj];
    if uj < ctx.forcing.node_lat_inflow_mode.len() {
        let m = ctx.forcing.node_lat_inflow_mode[uj];
        if m == ForcingMode::Override {
            q = ctx.forcing.node_lat_inflow_value[uj];
        } else if m == ForcingMode::Add {
            q += ctx.forcing.node_lat_inflow_value[uj];
        }
    }
    q
}

/// transport::bookMsxSystemOutflow: reaction species leaving with a node's system outflow.
fn book_msx_system_outflow(ctx: &mut SimulationContext, node: i32, q: f64, flooding: bool, dt: f64) {
    let rx = &mut ctx.reactions;
    let nsp = rx.n_species();
    if !rx.msx_ledger.armed || nsp <= 0 || !(q > 0.0) || node < 0 {
        return;
    }
    let uns = nsp as usize;
    let base = node as usize * uns;
    if base + uns > rx.msx_node_conc.len() || rx.msx_ledger.outflow.len() != uns {
        return;
    }
    for s in 0..uns {
        if rx.species_is_wall[s] == 0 {
            let v = q * rx.msx_node_conc[base + s] * dt;
            if flooding {
                rx.msx_ledger.flooding[s] += v;
            } else {
                rx.msx_ledger.outflow[s] += v;
            }
        }
    }
}

/// InletSolver::adjustFloodingTotals (legacy inlet_adjustQualOutflows): a capture node's
/// overflow goes back to the street, so it comes back out of the flooding just booked.
fn inlet_adjust_flooding_totals(inlet: &InletSolver, ctx: &mut SimulationContext, dt: f64) {
    let ni = inlet.soa_.count;
    if ni == 0 {
        return;
    }
    let nn = ctx.n_nodes();
    for j in 0..nn {
        let uj = j as usize;
        if uj >= inlet.is_capture_node_.len() || inlet.is_capture_node_[uj] == 0 {
            continue;
        }
        let q = ctx.nodes.overflow[uj];
        if q <= 0.0 {
            continue;
        }
        ctx.mass_balance.routing_flooding -= q * dt;
        ctx.mass_balance.step_flooding -= q;
        let np = ctx.n_pollutants();
        for p in 0..np as usize {
            ctx.mass_balance.qual_routing_flood[p] -= q * ctx.nodes.conc[uj * np as usize + p] * dt;
        }
    }
}

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

impl SwmmEngine {
    /// runtimeCouplingTime: elapsed seconds on the ms recurrence (no 2D source clock here).
    pub fn runtime_coupling_time(&self) -> f64 {
        self.ctx_.elapsed_ms / 1000.0
    }

    /// advanceTo(seconds, actual): step until the coupling clock reaches `seconds`.
    pub fn advance_to(&mut self, seconds: f64, actual: Option<&mut f64>) -> Result<i32, Terminated> {
        let actual = match actual {
            Some(a) => a,
            None => return Ok(CFFI_ERR_BADPARAM),
        };
        if !seconds.is_finite() {
            return Ok(CFFI_ERR_BADPARAM);
        }
        *actual = self.runtime_coupling_time();
        if self.step_active_ || self.ctx_.state != EngineState::Running {
            return Ok(CFFI_ERR_LIFECYCLE);
        }
        if seconds < *actual || seconds > self.ctx_.options.total_duration_ms() / 1000.0 {
            return Ok(CFFI_ERR_BADPARAM);
        }
        self.exchange_target_ = seconds;
        let mut rc = SWMM_OK;
        while self.runtime_coupling_time() < seconds {
            let before = self.ctx_.elapsed_ms;
            let mut days = 0.0;
            rc = match self.step(&mut days) {
                Ok(r) => r,
                Err(t) => {
                    // step() catches every std::exception itself; nothing else reaches here
                    self.exchange_target_ = -1.0;
                    return Err(t);
                }
            };
            if rc != SWMM_OK || self.ctx_.elapsed_ms <= before {
                break;
            }
        }
        self.exchange_target_ = -1.0;
        *actual = self.runtime_coupling_time();
        Ok(rc)
    }

    /// emit_progress: the progress callback with the fraction of the run done.
    pub(crate) fn emit_progress(&mut self) {
        if self.callbacks_.on_progress.is_none() {
            return;
        }
        let total = (self.ctx_.options.end_date - self.ctx_.options.start_date) * tc::SEC_PER_DAY;
        let frac = if total > 0.0 { self.ctx_.current_time / total } else { 0.0 };
        let date = self.ctx_.current_date;
        if let Some(cb) = self.callbacks_.on_progress.as_mut() {
            cb(frac, date);
        }
    }

    /// SWMMEngine::step(elapsed_time).
    pub fn step(&mut self, elapsed_time: &mut f64) -> Result<i32, Terminated> {
        if self.step_active_ {
            return Ok(CFFI_ERR_LIFECYCLE);
        }
        self.step_active_ = true;
        let r = self.step_body(elapsed_time);
        // StepGuard
        self.step_active_ = false;
        self.runtime_pipeline_ = false;
        r
    }

    /// step()'s catch (const std::exception& e): a thrown exception ends the step with
    /// SWMM_ERR_WRONG_STATE; anything else (a termination, a stub) goes on up.
    fn catch_step(&mut self, t: Terminated) -> Result<i32, Terminated> {
        if t.exception.starts_with(THROWN_PREFIX) {
            self.checkpoint_ready_ = false;
            self.set_error(SWMM_ERR_WRONG_STATE, t.what.as_bytes());
            return Ok(SWMM_ERR_WRONG_STATE);
        }
        Err(t)
    }

    fn step_body(&mut self, elapsed_time: &mut f64) -> Result<i32, Terminated> {
        if self.ctx_.state != EngineState::Running {
            *elapsed_time = 0.0;
            self.set_error(SWMM_ERR_WRONG_STATE, b"swmm_engine_step: engine is not running");
            return Ok(SWMM_ERR_WRONG_STATE);
        }
        if tc::simulation_complete(&self.ctx_) {
            *elapsed_time = 0.0;
            self.ctx_.state = EngineState::Ended;
            return Ok(SWMM_OK);
        }

        // ---- the try block ----
        if self.checkpoint_enabled_ {
            return Err(stub("SWMMEngine checkpointing (checkCheckpointBoundary)"));
        }
        let mut dt_next;
        if !self.do_routing_ {
            // IGNORE_ROUTING: MIN(WetStep, ReportStep), bounded by the remaining duration.
            let mut dt = cmin(self.ctx_.options.wet_step, self.ctx_.options.report_step);
            let total_msec = self.ctx_.options.total_duration_ms();
            if self.ctx_.elapsed_ms + 1000.0 * dt > total_msec {
                dt = (total_msec - self.ctx_.elapsed_ms) / 1000.0;
                dt = cmax(dt, 0.001);
            }
            dt_next = dt;
        } else {
            // legacy routing_getRoutingStep between routing events
            let mut dt_between = 0.0;
            if !self.ctx_.events.is_empty() && self.between_events_ && self.next_event_ < self.ctx_.events.len() as i32 {
                let c = &self.ctx_;
                let fixed = c.options.routing_step;
                let ev_start = c.events[self.next_event_ as usize].start;
                let next_ms = cmin(self.new_runoff_ms_, c.next_report_ms);
                let date1 = datetime::add_seconds(c.options.start_date, (c.elapsed_ms + 1.0) / 1000.0);
                let date2 = datetime::add_seconds(c.options.start_date, (next_ms + 1.0) / 1000.0);
                if date2 > date1 && date2 < ev_start {
                    dt_between = (next_ms - c.elapsed_ms) / 1000.0;
                } else {
                    let d1 = datetime::add_seconds(c.options.start_date, (c.elapsed_ms + 1000.0 * fixed + 1.0) / 1000.0);
                    if d1 < ev_start {
                        dt_between = fixed;
                    }
                }
            }
            let mut dt_cfl = self.ctx_.options.routing_step;
            if dt_between > 0.0 {
                dt_cfl = dt_between;
            } else if self.ctx_.options.variable_step > 0.0 {
                let (rs, vs) = (self.ctx_.options.routing_step, self.ctx_.options.variable_step);
                dt_cfl = self.router_.get_adaptive_step(&mut self.ctx_, rs, vs);
                if dt_cfl > 0.0 && dt_cfl < self.ctx_.options.routing_step {
                    let courant = self.ctx_.options.routing_step / dt_cfl;
                    self.ctx_.routing_stats.max_courant = cmax(self.ctx_.routing_stats.max_courant, courant);
                }
            }
            // Control-rule grid alignment (legacy routing.c:196-204)
            if self.ctx_.options.rule_step > 0.0 {
                self.ctx_.dt_controls_remaining = (self.new_rule_time_ms_ + 1000.0 * self.ctx_.options.rule_step - self.ctx_.elapsed_ms) / 1000.0;
            }
            dt_next = tc::compute_next(&self.ctx_, dt_cfl);
            // The between-events jump is not bounded by the routing step.
            if dt_between > self.ctx_.options.routing_step {
                dt_next = dt_between;
                let total_msec = self.ctx_.options.total_duration_ms();
                if self.ctx_.elapsed_ms + 1000.0 * dt_next > total_msec {
                    dt_next = cmax((total_msec - self.ctx_.elapsed_ms) / 1000.0, 0.001);
                }
                if self.ctx_.dt_controls_remaining > 0.0 && self.ctx_.dt_controls_remaining < dt_next {
                    dt_next = self.ctx_.dt_controls_remaining;
                }
            }
        }

        // bound_runtime_step
        let now = self.runtime_coupling_time();
        if self.exchange_target_ >= 0.0 {
            dt_next = cmin(dt_next, self.exchange_target_ - now);
        }
        #[cfg(feature = "twod")]
        {
            dt_next = cmin(dt_next, self.surface_router_.next_runtime_time(now)? - now);
        }
        self.emit_progress();
        if let Some(cb) = self.callbacks_.on_step_begin.as_mut() {
            cb(self.ctx_.current_date, dt_next);
        }
        // A begin callback may install a source expiring during this step.
        let now = self.runtime_coupling_time();
        if self.exchange_target_ >= 0.0 {
            dt_next = cmin(dt_next, self.exchange_target_ - now);
        }
        #[cfg(feature = "twod")]
        {
            dt_next = cmin(dt_next, self.surface_router_.next_runtime_time(now)? - now);
        }
        self.runtime_pipeline_ = true;

        // Snapshot the every-step old state (legacy initSystemInflows' oldLatFlow and the
        // quality old state); the hydraulic old state rolls in stepRouting.
        self.ctx_.save_lat_qual_state();

        // Reset mass balance accumulators (legacy massbal_initTimeStepTotals)
        self.reset_step_mass_balance();

        // ---- user-injected runtime forcings ----
        self.apply_forcings(dt_next)?;

        // ---- the simulation pipeline (legacy swmm_step order) ----
        if let Err(t) = self.step_runoff(dt_next) {
            return self.catch_step(t);
        }
        if self.do_routing_ {
            let fv_window_totals = self.ctx_.options.routing_model == RoutingModel::Fv;
            if !fv_window_totals {
                self.accumulate_node_routing_totals(dt_next / 2.0); // routing.c:227
            }
            if let Err(t) = self.step_routing(dt_next) {
                return self.catch_step(t);
            }
            let drc = self.check_routing_diverged()?;
            if drc != SWMM_OK {
                return Ok(drc);
            }
            let routing_before = self.routing_terms();
            if !self.between_events_ {
                if self.ctx_.n_links() > 0 {
                    let first = self.ctx_.elapsed_ms <= 0.0;
                    let steady = self.last_step_steady_;
                    self.ctx_.routing_stats.record_step(dt_next, steady, first);
                }
                self.update_routing_mass_balance(dt_next);
                self.update_statistics(dt_next);
            }
            // [REPORT] LINK_STEPS
            if self.router_.link_step_stats_enabled() {
                let step_end = datetime::add_seconds(self.ctx_.options.start_date, (self.ctx_.elapsed_ms + 1000.0 * dt_next + 1.0) / 1000.0);
                let in_window = !self.between_events_ && step_end >= self.ctx_.options.report_start;
                self.router_.accumulate_link_step_stats(&mut self.ctx_, dt_next, in_window);
            }
            let after = self.routing_terms();
            let half = dt_next / 2.0;
            for k in 0..11 {
                let rate = if dt_next > 0.0 { (after[k] - routing_before[k]) / dt_next } else { 0.0 };
                let t = &mut self.ctx_.mass_balance.routing_report[k];
                if fv_window_totals {
                    *t += after[k] - routing_before[k];
                } else {
                    *t += self.routing_prev_rates_[k] * half;
                    *t += rate * half;
                }
                self.routing_prev_rates_[k] = rate;
            }
            if !fv_window_totals || !self.between_events_ {
                self.accumulate_node_routing_totals(if fv_window_totals { dt_next } else { dt_next / 2.0 });
            }
        }
        // OPENSWMM_HAS_2D: a runtime prescription or an advanceTo target closes the pending 2D
        // batch now (the checkpointed surface output is off without engine checkpointing).
        #[cfg(feature = "twod")]
        if self.surface_router_.runtime_active() || self.exchange_target_ >= 0.0 {
            if let Err(t) = self.surface_router_.flush_pending_batch(&mut self.ctx_) {
                return self.catch_step(t);
            }
        }
        reassert_node_state_overrides(&mut self.ctx_);
        self.compute_final_storage();
        if !self.ctx_.options.ignore_quality {
            self.compute_final_quality_mass_balance();
        }

        // ---- Clear auto-reset forcings ----
        self.ctx_.forcing.clear_reset_entries();

        // Advance clock (must happen before output_due check)
        tc::advance(&mut self.ctx_, dt_next);

        // Advance the control-rule clock when the step landed on its grid (1 ms tolerance).
        if self.do_routing_ && self.ctx_.options.rule_step > 0.0 && (self.ctx_.elapsed_ms - (self.new_rule_time_ms_ + 1000.0 * self.ctx_.options.rule_step)).abs() < 1.0 {
            self.new_rule_time_ms_ += 1000.0 * self.ctx_.options.rule_step;
        }

        // Time-step averaging with legacy's update ordering (swmm5.c saveResults).
        let avg_out_due = self.ctx_.options.rpt_averages && tc::output_due(&self.ctx_);
        let avg_exact = avg_out_due && (self.ctx_.elapsed_ms == self.ctx_.next_report_ms);
        if self.ctx_.options.rpt_averages && (!avg_out_due || avg_exact) {
            self.accumulate_avg_results()?;
        }

        // Post snapshot after advance so output_due() fires correctly.
        self.post_output_snapshot(dt_next)?;

        if avg_out_due && !avg_exact {
            self.accumulate_avg_results()?;
        }

        // Dated SAVE HOTSTART rows (legacy swmm_step: hotstart_save after saveResults).
        if !self.checkpoint_enabled_ {
            if let Err(t) = self.write_due_hotstart_saves() {
                return self.catch_step(t);
            }
        }

        // Expire a one-step API rainfall prescription.
        self.ctx_.gages.api_rainfall.iter_mut().for_each(|x| *x = -1.0);

        // OPENSWMM_HAS_2D: retire the runtime sources and forcings that expired (the subsurface's
        // independent clock is phase 2b).
        #[cfg(feature = "twod")]
        {
            let completed = self.runtime_coupling_time();
            if let Err(t) = self.surface_router_.expire_runtime(completed) {
                return self.catch_step(t);
            }
        }

        self.runtime_pipeline_ = false;
        if let Some(cb) = self.callbacks_.on_step_end.as_mut() {
            cb(self.ctx_.current_date, dt_next);
        }
        // (checkpoint_enabled_ is false here: see the top of the try block)
        *elapsed_time = self.ctx_.current_time / tc::SEC_PER_DAY;
        Ok(SWMM_OK)
    }

    /// The eleven routing continuity terms, in the order step() reads them.
    fn routing_terms(&self) -> [f64; 11] {
        let m = &self.ctx_.mass_balance;
        [
            m.routing_dry_weather,
            m.routing_wet_weather,
            m.routing_gw_inflow,
            m.routing_link_gw_inflow,
            m.routing_rdii,
            m.routing_external,
            m.routing_flooding,
            m.routing_coupling_out,
            m.routing_outflow,
            m.routing_evap_loss,
            m.routing_seep_loss,
        ]
    }

    /// SWMMEngine::resetStepMassBalance (legacy massbal_initTimeStepTotals).
    pub(crate) fn reset_step_mass_balance(&mut self) {
        let m = &mut self.ctx_.mass_balance;
        m.step_flooding = 0.0;
        m.step_outflow = 0.0;
        m.step_dw_inflow = 0.0;
        m.step_gw_inflow = 0.0;
        m.step_rdii_inflow = 0.0;
        m.step_ext_inflow = 0.0;
        m.step_ext_withdrawal = 0.0;
    }

    /// SWMMEngine::applyForcings(dt): the user-injected runtime forcings of this step.
    pub(crate) fn apply_forcings(&mut self, dt: f64) -> Result<(), Terminated> {
        let c = &mut self.ctx_;
        // ---- Gage rainfall forcing ----
        for g in 0..c.n_gages() {
            let ug = g as usize;
            if c.options.ignore_rainfall && c.forcing.gage_rainfall_mode[ug] != ForcingMode::None {
                c.gages.rainfall[ug] = 0.0;
                continue;
            }
            if c.forcing.gage_rainfall_mode[ug] == ForcingMode::Override {
                c.gages.rainfall[ug] = c.forcing.gage_rainfall_value[ug];
            } else if c.forcing.gage_rainfall_mode[ug] == ForcingMode::Add {
                let r = crate::hydrology::gage_report::get_report_rainfall(c, g as i32, c.current_date);
                c.gages.rainfall[ug] = r + c.forcing.gage_rainfall_value[ug];
            }
        }
        // ---- Link setting forcing ----
        for j in 0..c.n_links() as usize {
            if c.forcing.link_setting_mode[j] == ForcingMode::Override {
                c.links.setting[j] = c.forcing.link_setting_value[j];
                c.links.target_setting[j] = c.forcing.link_setting_value[j];
            } else if c.forcing.link_setting_mode[j] == ForcingMode::Add {
                c.links.setting[j] += c.forcing.link_setting_value[j];
                c.links.target_setting[j] = c.links.setting[j];
            }
        }
        // ---- Link flow forcing ----
        for j in 0..c.n_links() as usize {
            if c.forcing.link_flow_mode[j] == ForcingMode::Override {
                c.links.flow[j] = c.forcing.link_flow_value[j];
            } else if c.forcing.link_flow_mode[j] == ForcingMode::Add {
                c.links.flow[j] += c.forcing.link_flow_value[j];
            }
        }
        // ---- U2: reserved-species node forcing (temperature, age) ----
        let nn = c.n_nodes() as usize;
        if c.options.heat_transport && c.heat_state.node_temp.len() == nn {
            for u in 0..nn {
                if u >= c.forcing.node_temperature_mode.len() {
                    break;
                }
                if c.forcing.node_temperature_mode[u] == ForcingMode::Override {
                    c.heat_state.node_temp[u] = c.forcing.node_temperature_value[u];
                } else if c.forcing.node_temperature_mode[u] == ForcingMode::Add && u < c.heat_state.node_temp_vol_in.len() {
                    c.heat_state.node_temp_vol_in[u] += c.forcing.node_temperature_value[u];
                }
            }
        }
        if c.options.water_age && c.water_age_state.node_age.len() == nn {
            for u in 0..nn {
                if u >= c.forcing.node_age_mode.len() {
                    break;
                }
                if c.forcing.node_age_mode[u] == ForcingMode::Override {
                    c.water_age_state.node_age[u] = c.forcing.node_age_value[u] * 3600.0;
                } else if c.forcing.node_age_mode[u] == ForcingMode::Add && u < c.water_age_state.node_age_vol_in.len() {
                    c.water_age_state.node_age_vol_in[u] += c.forcing.node_age_value[u] * 3600.0;
                }
            }
        }
        // ---- Node quality mass flux forcing ----
        let np = c.n_pollutants();
        if np > 0 {
            for i in 0..nn {
                for p in 0..np as usize {
                    let flat = i * np as usize + p;
                    if c.forcing.node_quality_mode[flat] == ForcingMode::Override {
                        c.nodes.conc[flat] = c.forcing.node_quality_value[flat];
                    } else if c.forcing.node_quality_mode[flat] == ForcingMode::Add {
                        c.nodes.conc[flat] += c.forcing.node_quality_value[flat];
                        c.mass_balance.routing_forcing_qual_inflow[p] += c.forcing.node_quality_value[flat] * dt;
                    }
                }
            }
        }
        Ok(())
    }

    /// SWMMEngine::accumulateNodeRoutingTotals(half_step) (legacy
    /// massbal_updateRoutingTotals, massbal.c:581-595).
    pub(crate) fn accumulate_node_routing_totals(&mut self, half_step: f64) {
        let n = &mut self.ctx_.nodes;
        for uj in 0..self.ctx_.node_names.size() as usize {
            let qin = n.inflow[uj];
            n.stat_total_inflow_vol[uj] += qin * half_step;
            if n.type_[uj] == NodeType::Outfall || (n.degree[uj] == 0 && n.type_[uj] != NodeType::Storage) {
                n.stat_total_outflow_vol[uj] += qin * half_step;
            } else {
                n.stat_total_outflow_vol[uj] += n.outflow[uj] * half_step;
                if n.volume[uj] <= n.full_volume[uj] {
                    n.stat_total_outflow_vol[uj] += n.overflow[uj] * half_step;
                }
            }
        }
    }

    /// SWMMEngine::checkRoutingDiverged: stop at a non-finite or absurd head or flow.
    pub(crate) fn check_routing_diverged(&mut self) -> Result<i32, Terminated> {
        const K_ABSURD: f64 = 1.0e9;
        let diverged = |v: f64| !v.is_finite() || v.abs() > K_ABSURD;
        // FV cell-state guard (issue #156 R3): the FV solver's own check of its cells.
        if let Some(dv) = self.router_.fv_solver_.as_ref().and_then(|s| s.divergence()) {
            let hours = self.ctx_.current_time / 3600.0;
            let name: Vec<u8> = if dv.link >= 0 { crate::core::cstdlib::c_str(self.ctx_.link_names.name_of(dv.link)).to_vec() } else { b"?".to_vec() };
            let what: &[u8] = if dv.what.is_empty() { b"state" } else { dv.what.as_bytes() };
            let buf = cat(&[
                b"ERROR ",
                CFFI_ERR_NUMERICAL.to_string().as_bytes(),
                b": the FV solution diverged in conduit '",
                &name,
                b"' (",
                what,
                b" ",
                &crate::core::cstdio::printf_g(dv.value, 6),
                b") after ",
                format!("{:.4}", hours).as_bytes(),
                b" hours -- a step-size-independent instability the substep retry cannot shrink away; the run cannot continue.",
            ]);
            let buf = truncate255(buf);
            self.ctx_.errors.push(buf.clone());
            self.set_error(CFFI_ERR_NUMERICAL, &buf);
            return Ok(CFFI_ERR_NUMERICAL);
        }
        let hours = self.ctx_.current_time / 3600.0;
        for i in 0..self.ctx_.n_nodes() {
            let ui = i as usize;
            if !diverged(self.ctx_.nodes.head[ui]) {
                continue;
            }
            let buf = cat(&[
                b"ERROR ",
                CFFI_ERR_NUMERICAL.to_string().as_bytes(),
                b": the routing solution diverged at node '",
                crate::core::cstdlib::c_str(self.ctx_.node_names.name_of(i)),
                b"' (head ",
                &crate::core::cstdio::printf_g(self.ctx_.nodes.head[ui], 6),
                b") after ",
                format!("{:.4}", hours).as_bytes(),
                b" hours -- the run cannot continue.",
            ]);
            let buf = truncate255(buf);
            self.ctx_.errors.push(buf.clone());
            self.set_error(CFFI_ERR_NUMERICAL, &buf);
            return Ok(CFFI_ERR_NUMERICAL);
        }
        for j in 0..self.ctx_.n_links() {
            let uj = j as usize;
            if !diverged(self.ctx_.links.flow[uj]) {
                continue;
            }
            let buf = cat(&[
                b"ERROR ",
                CFFI_ERR_NUMERICAL.to_string().as_bytes(),
                b": the routing solution diverged at link '",
                crate::core::cstdlib::c_str(self.ctx_.link_names.name_of(j)),
                b"' (flow ",
                &crate::core::cstdio::printf_g(self.ctx_.links.flow[uj], 6),
                b") after ",
                format!("{:.4}", hours).as_bytes(),
                b" hours -- the run cannot continue.",
            ]);
            let buf = truncate255(buf);
            self.ctx_.errors.push(buf.clone());
            self.set_error(CFFI_ERR_NUMERICAL, &buf);
            return Ok(CFFI_ERR_NUMERICAL);
        }
        Ok(SWMM_OK)
    }

    /// SWMMEngine::ensureXspCache: the cache's size and generation (see engine.rs).
    pub(crate) fn ensure_xsp_cache(&mut self) {
        let n = self.ctx_.n_links() as usize;
        if self.xsp_cache_len_ == n && self.xsp_cache_gen_ == self.ctx_.xsect_generation {
            return;
        }
        self.xsp_cache_len_ = n;
        self.xsp_cache_gen_ = self.ctx_.xsect_generation;
    }

    /// SWMMEngine::updateStatistics(dt_routing): node and link statistics after routing
    /// (legacy stats_updateFlowStats), from REPORT_START on.
    pub(crate) fn update_statistics(&mut self, dt_routing: f64) {
        let np = self.ctx_.n_pollutants();
        let stat_date = datetime::add_seconds(self.ctx_.options.start_date, (self.ctx_.elapsed_ms + 1000.0 * dt_routing + 1.0) / 1000.0);
        if stat_date < self.ctx_.options.report_start {
            return;
        }
        self.ctx_.routing_stats.report_steps += 1;
        self.ctx_.routing_stats.report_time += dt_routing;
        let mut sys_outfall_flow = 0.0;
        let is_dw_model = self.ctx_.options.routing_model == RoutingModel::Dynwave;
        {
            let c = &mut self.ctx_;
            let n = &mut c.nodes;
            for uj in 0..c.node_names.size() as usize {
                let cur_depth = n.depth[uj];
                n.stat_sum_depth[uj] += cur_depth;
                if n.type_[uj] == NodeType::Storage {
                    let volume = n.volume[uj];
                    n.stat_sum_volume[uj] += volume;
                    let capped = cmin(volume, n.full_volume[uj]);
                    if capped > n.stat_storage_max_volume[uj] {
                        n.stat_storage_max_volume[uj] = capped;
                        n.stat_storage_max_volume_date[uj] = stat_date;
                    }
                }
                if cur_depth > n.stat_max_depth[uj] {
                    n.stat_max_depth[uj] = cur_depth;
                    n.stat_max_depth_date[uj] = stat_date;
                }
                if n.overflow[uj] > n.stat_max_overflow[uj] {
                    n.stat_max_overflow[uj] = n.overflow[uj];
                    n.stat_max_overflow_date[uj] = stat_date;
                }
                if n.type_[uj] != NodeType::Outfall {
                    if n.volume[uj] > n.full_volume[uj] || n.overflow[uj] > 0.0 {
                        n.stat_time_flooded[uj] += dt_routing;
                        n.stat_vol_flooded[uj] += n.overflow[uj] * dt_routing;
                    }
                    if is_dw_model && (n.type_[uj] != NodeType::Storage || n.sur_depth[uj] > 0.0) && cur_depth + n.invert_elev[uj] + crate::core::constants::FUDGE >= n.crown_elev[uj] {
                        n.stat_time_surcharged[uj] += dt_routing;
                    }
                }
                let lat = n.lat_flow[uj];
                let total_inflow = n.inflow[uj];
                if lat.abs() > n.stat_max_lat_inflow[uj].abs() {
                    n.stat_max_lat_inflow[uj] = lat;
                }
                if total_inflow > n.stat_max_total_inflow[uj] {
                    n.stat_max_total_inflow[uj] = total_inflow;
                    n.stat_max_inflow_date[uj] = stat_date;
                }
                let lateral_pair = n.old_lat_flow[uj] + lat;
                n.stat_lat_inflow_vol[uj] += lateral_pair * 0.5 * dt_routing;
                if n.type_[uj] == NodeType::Storage && n.outflow[uj] > n.stat_storage_max_outflow[uj] {
                    n.stat_storage_max_outflow[uj] = n.outflow[uj];
                }
                if n.type_[uj] == NodeType::Outfall {
                    let qi = n.inflow[uj];
                    sys_outfall_flow += qi;
                    if qi >= 0.001 {
                        n.stat_outfall_avg_flow[uj] += qi;
                        if qi > n.stat_outfall_max_flow[uj] {
                            n.stat_outfall_max_flow[uj] = qi;
                        }
                        n.stat_outfall_periods[uj] += 1;
                    }
                    if np > 0 {
                        let base = uj * np as usize;
                        for p in 0..np as usize {
                            let idx = base + p;
                            if idx < n.stat_total_load.len() && idx < n.conc.len() {
                                n.stat_total_load[idx] += qi * n.conc[idx] * dt_routing;
                            }
                        }
                    }
                }
            }
        }
        self.ctx_.routing_stats.max_outfall_flow = cmax(self.ctx_.routing_stats.max_outfall_flow, sys_outfall_flow);
        self.ensure_xsp_cache();
        let mut step_slot_vol = 0.0;
        let mut step_stored_vol = 0.0;
        let rm = self.ctx_.options.routing_model;
        let c = &mut self.ctx_;
        for j in 0..c.n_links() {
            let uj = j as usize;
            let q = c.links.flow[uj].abs();
            if q > c.links.stat_max_flow[uj] {
                c.links.stat_max_flow[uj] = q;
                c.links.stat_max_flow_date[uj] = stat_date;
            }
            c.links.stat_vol_flow[uj] += q * dt_routing;
            let mut vel = 0.0;
            if c.links.type_[uj] == LinkType::Conduit {
                let cr = c.link_subtypes.conduit_row(j);
                let nb = if cr >= 0 { c.link_subtypes.conduits.barrels[cr as usize] } else { 1 };
                let xs = link::build_xsect_params(&c.links, uj, Some(&c.transect_tables));
                vel = link::get_velocity(&xs, q, c.links.depth[uj], nb);
            }
            if vel > c.links.stat_max_veloc[uj] {
                c.links.stat_max_veloc[uj] = vel;
            }
            let y_full = c.links.xsect_y_full[uj];
            let d = c.links.depth[uj];
            let filling = if y_full > 0.0 { d / y_full } else { 0.0 };
            if filling > c.links.stat_max_filling[uj] {
                c.links.stat_max_filling[uj] = filling;
            }
            {
                let v = c.links.volume[uj];
                let sv = c.links.slot_volume[uj];
                c.links.stat_vol_dt[uj] += v * dt_routing;
                step_stored_vol += v;
                if sv > 0.0 {
                    c.links.stat_slot_vol_dt[uj] += sv * dt_routing;
                    step_slot_vol += sv;
                    let share = if v > 0.0 { sv / v } else { 0.0 };
                    if share > c.links.stat_peak_slot_share[uj] {
                        c.links.stat_peak_slot_share[uj] = share;
                    }
                    if share > 0.01 {
                        c.links.stat_time_slot_above[uj] += dt_routing;
                    }
                }
            }
            if c.links.type_[uj] == LinkType::Conduit && c.links.xsect_shape[uj] != XsectShape::Dummy {
                let n1 = c.links.node1[uj];
                let n2 = c.links.node2[uj];
                let cr = c.link_subtypes.conduit_row(j);
                let cd = &c.link_subtypes.conduits;
                let fc = c.links.flow_class[uj] as i32;
                if fc >= 0 && fc < LinkData::N_FLOW_CLASSES {
                    let fc_idx = uj * LinkData::N_FLOW_CLASSES as usize + fc as usize;
                    if fc_idx < c.links.stat_flow_class.len() {
                        c.links.stat_flow_class[fc_idx] += dt_routing;
                    }
                }
                if cr >= 0 {
                    let ucr = cr as usize;
                    if q >= cd.q_full[ucr] * cd.barrels[ucr] as f64 {
                        c.links.stat_time_full_both[uj] += dt_routing;
                    }
                }
                let fs: i8 = if cr >= 0 { cd.full_state[cr as usize] } else { 0 };
                let up_full = (fs & 1) != 0;
                let dn_full = (fs & 2) != 0;
                if up_full && dn_full {
                    c.links.stat_time_surcharged[uj] += dt_routing;
                    c.links.stat_time_full_upstream[uj] += dt_routing;
                    c.links.stat_time_full_dnstream[uj] += dt_routing;
                } else if up_full {
                    c.links.stat_time_full_upstream[uj] += dt_routing;
                } else if dn_full {
                    c.links.stat_time_full_dnstream[uj] += dt_routing;
                }
                let dw = rm == RoutingModel::Dynwave;
                if if dw { (fs & 4) != 0 } else { up_full } {
                    let mut cap_ltd = true;
                    if (dw || rm == RoutingModel::Fv) && n1 >= 0 && n2 >= 0 {
                        let h1h = c.nodes.head[n1 as usize];
                        let h2h = c.nodes.head[n2 as usize];
                        let len = if cr >= 0 { cd.length[cr as usize] } else { 0.0 };
                        let slp = (if cr >= 0 { cd.slope[cr as usize] } else { 0.0 }).abs();
                        cap_ltd = (h1h - h2h) > slp * len;
                    }
                    if cap_ltd {
                        c.links.stat_time_capacity_limited[uj] += dt_routing;
                    }
                }
            }
            {
                let crs = c.link_subtypes.conduit_row(j);
                if crs >= 0 {
                    let ucrs = crs as usize;
                    if c.link_subtypes.conduits.normal_flow_limited[ucrs] != 0 {
                        c.links.stat_norm_ltd[uj] += dt_routing;
                    }
                    if c.link_subtypes.conduits.inlet_control[ucrs] != 0 {
                        c.links.stat_inlet_ctrl[uj] += dt_routing;
                    }
                }
            }
            if c.links.type_[uj] == LinkType::Pump {
                const MIN_PUMP_FLOW: f64 = 0.001;
                let is_on = q > MIN_PUMP_FLOW;
                if is_on && c.links.old_flow[uj] < MIN_PUMP_FLOW {
                    c.links.stat_pump_cycles[uj] += 1;
                }
                c.links.stat_pump_was_on[uj] = is_on as _;
                if is_on {
                    c.links.stat_pump_on_time[uj] += dt_routing;
                    c.links.stat_pump_volume[uj] += q * dt_routing;
                    c.links.stat_pump_sum_flow[uj] += q;
                    c.links.stat_pump_periods[uj] += 1;
                    if c.links.flow_class[uj] == FlowClass::DnDry {
                        c.links.stat_pump_off_curve_low[uj] += dt_routing;
                    }
                    if c.links.flow_class[uj] == FlowClass::UpDry {
                        c.links.stat_pump_off_curve_high[uj] += dt_routing;
                    }
                    let pn1 = c.links.node1[uj];
                    let pn2 = c.links.node2[uj];
                    if pn1 >= 0 && pn2 >= 0 {
                        let un1 = pn1 as usize;
                        let un2 = pn2 as usize;
                        let h1 = c.nodes.invert_elev[un1] + c.nodes.depth[un1];
                        let h2 = c.nodes.invert_elev[un2] + c.nodes.depth[un2];
                        let dh = (h1 - h2).abs();
                        const KW_PER_HP: f64 = 0.7457;
                        let power_kw = dh * q / 8.814 * KW_PER_HP;
                        c.links.stat_pump_energy[uj] += power_kw * dt_routing / 3600.0;
                    }
                }
            }
            {
                let dq = c.links.flow[uj] - c.links.old_flow[uj];
                let old_sign = c.links.stat_flow_turn_sign[uj] as i32;
                let new_sign: i32 = if dq > 0.0 { 1 } else if dq < 0.0 { -1 } else { 0 };
                c.links.stat_flow_turn_sign[uj] = new_sign as _;
                if dq.abs() > 0.001 && old_sign * new_sign < 0 {
                    c.links.stat_flow_turns[uj] += 1;
                }
            }
            if np > 0 && q > 0.0 {
                let base = uj * np as usize;
                for p in 0..np as usize {
                    let idx = base + p;
                    if idx < c.links.stat_total_load.len() && idx < c.links.conc.len() {
                        c.links.stat_total_load[idx] += q * c.links.conc[idx] * dt_routing;
                    }
                }
            }
        }
        if step_slot_vol > 0.0 {
            let sys_share = if step_stored_vol > 0.0 { step_slot_vol / step_stored_vol } else { 0.0 };
            if sys_share > c.routing_stats.slot_peak_share {
                c.routing_stats.slot_peak_share = sys_share;
            }
            if sys_share > 0.01 {
                c.routing_stats.slot_time_above_s += dt_routing;
            }
        }
    }

    /// SWMMEngine::updateRoutingMassBalance(dt_routing): routing continuity after routing
    /// (legacy massbal_updateRoutingTotals, removeSystemOutflows, removeStorageLosses).
    pub(crate) fn update_routing_mass_balance(&mut self, dt_routing: f64) {
        let fv_accepted = self.quality_.fv_transport_.accepted_;
        {
            let c = &mut self.ctx_;
            c.mass_balance.step_flooding = 0.0;
            c.mass_balance.step_outflow = c.mass_balance.step_ext_withdrawal;
            c.mass_balance.routing_outflow += c.mass_balance.step_ext_withdrawal * dt_routing;
            for r in 0..c.node_subtypes.storages.lid.len() {
                let st = &c.node_subtypes.storages;
                if st.lid[r].bottom != LidNodeBottom::MeshAquifer || !st.lid_state[r].richards {
                    continue;
                }
                let receipt = &st.lid_state[r].richards_report;
                let (bin, bout) = (receipt.bottom_inflow, receipt.bottom_outflow);
                c.mass_balance.step_ext_inflow += bin / (0.028316846592 * dt_routing);
                c.mass_balance.routing_coupling_out += bout / 0.028316846592;
                c.mass_balance.step_outflow += bout / (0.028316846592 * dt_routing);
            }
            c.mass_balance.routing_dry_weather += c.mass_balance.step_dw_inflow * dt_routing;
            c.mass_balance.routing_gw_inflow += c.mass_balance.step_gw_inflow * dt_routing;
            c.mass_balance.routing_rdii += c.mass_balance.step_rdii_inflow * dt_routing;
            c.mass_balance.routing_external += c.mass_balance.step_ext_inflow * dt_routing;
        }
        let mut step_loss_rate = 0.0;
        let is_dw = self.ctx_.options.routing_model == RoutingModel::Dynwave || self.ctx_.options.routing_model == RoutingModel::Fv;
        let np = self.ctx_.n_pollutants();
        for j in 0..self.ctx_.n_nodes() {
            let uj = j as usize;
            let nt = self.ctx_.nodes.type_[uj];
            let mut quality_outflow = 0.0;
            let mut quality_flooding = false;
            if nt == NodeType::Outfall {
                let q_in = self.ctx_.nodes.inflow[uj];
                let q_out = self.ctx_.nodes.outflow[uj];
                {
                    let ofr = self.ctx_.node_subtypes.outfall_row(j);
                    let sc = if ofr >= 0 { self.ctx_.node_subtypes.outfalls.route_to[ofr as usize] } else { -1 };
                    if sc >= 0 && sc < self.ctx_.n_subcatches() && !self.owns_spatial_source(sc) && q_in > 0.0 {
                        let c = &mut self.ctx_;
                        let usc = sc as usize;
                        c.subcatches.outfall_runon_vol[usc] += q_in * dt_routing;
                        if c.options.water_age && usc < c.water_age_state.subcatch_outfall_age_vol.len() && uj < c.water_age_state.node_age.len() {
                            c.water_age_state.subcatch_outfall_age_vol[usc] += q_in * dt_routing * c.water_age_state.node_age[uj];
                        }
                        if c.options.heat_transport && usc < c.heat_state.subcatch_outfall_temp_vol.len() && uj < c.heat_state.node_temp.len() {
                            c.heat_state.subcatch_outfall_temp_vol[usc] += q_in * dt_routing * c.heat_state.node_temp[uj];
                        }
                    }
                }
                let mut q_sys = 0.0;
                {
                    let c = &mut self.ctx_;
                    if uj < c.surface_outfall_link_limit.len() && c.surface_outfall_link_limit[uj].is_finite() {
                        q_sys = q_in - q_out;
                        if q_sys < 0.0 {
                            c.nodes.inflow[uj] = -q_sys;
                        }
                    } else if c.options.routing_model == RoutingModel::Fv {
                        if q_in > 0.0 && q_out > 0.0 {
                            c.mass_balance.routing_external += q_out * dt_routing;
                            c.mass_balance.step_ext_inflow += q_out;
                            c.outfall_exports.reverse_m3[uj] += q_out * dt_routing * 0.028316846592;
                        }
                        if q_out > 0.0 && q_in <= 0.0 {
                            q_sys = -q_out;
                            c.nodes.inflow[uj] = q_out;
                        } else {
                            q_sys = q_in;
                        }
                    } else if q_out == 0.0 {
                        q_sys = q_in;
                    } else if q_in == 0.0 {
                        q_sys = -q_out;
                        c.nodes.inflow[uj] = q_sys.abs();
                    }
                }
                if q_sys > 0.0 {
                    const CUBIC_METRES_PER_CUBIC_FOOT: f64 = 0.028316846592;
                    {
                        let c = &mut self.ctx_;
                        let volume = q_sys * dt_routing * CUBIC_METRES_PER_CUBIC_FOOT;
                        c.outfall_exports.water_m3[uj] += volume;
                        if c.options.heat_transport && uj < c.heat_state.node_temp.len() {
                            c.outfall_exports.heat_j[uj] += volume * c.heat_state.node_temp[uj] * c.options.water_density * c.options.water_specific_heat;
                        }
                        for p in 0..np as usize {
                            let index = uj * np as usize + p;
                            let units = c.pollutants.units[p];
                            let factor = if units == MassUnits::MgPerL {
                                1.0e-3
                            } else if units == MassUnits::UgPerL {
                                1.0e-6
                            } else {
                                1000.0
                            };
                            c.outfall_exports.pollutant_amounts[index] += volume * c.nodes.conc[index] * factor;
                        }
                    }
                    let row = self.ctx_.node_subtypes.outfall_row(j);
                    let sc = if row >= 0 { self.ctx_.node_subtypes.outfalls.route_to[row as usize] } else { -1 };
                    if self.owns_spatial_source(sc) {
                        self.ctx_.subcatches.outfall_runon_vol[sc as usize] += q_sys * dt_routing;
                    }
                    quality_outflow = q_sys;
                    self.ctx_.mass_balance.routing_outflow += q_sys * dt_routing;
                    self.ctx_.mass_balance.step_outflow += q_sys;
                } else {
                    let c = &mut self.ctx_;
                    c.mass_balance.routing_external += -q_sys * dt_routing;
                    c.outfall_exports.reverse_m3[uj] += -q_sys * dt_routing * 0.028316846592;
                    c.mass_balance.step_ext_inflow += -q_sys;
                }
                if !is_dw {
                    self.ctx_.nodes.overflow[uj] = 0.0;
                    self.ctx_.nodes.volume[uj] = 0.0;
                }
            } else if !is_dw && self.ctx_.nodes.degree[uj] == 0 && nt != NodeType::Storage {
                let c = &mut self.ctx_;
                if c.nodes.outflow[uj] == 0.0 {
                    quality_outflow = c.nodes.inflow[uj];
                    c.mass_balance.routing_outflow += c.nodes.inflow[uj] * dt_routing;
                    c.mass_balance.step_outflow += c.nodes.inflow[uj];
                }
                c.nodes.overflow[uj] = 0.0;
                c.nodes.volume[uj] = 0.0;
            } else {
                let c = &mut self.ctx_;
                let full_vol = if c.options.routing_model == RoutingModel::Fv { c.nodes.full_volume[uj] } else { c.nodes.rpt_full_volume[uj] };
                if c.nodes.overflow[uj] > 0.0 && c.nodes.volume[uj] <= full_vol {
                    quality_outflow = c.nodes.overflow[uj];
                    quality_flooding = true;
                    c.mass_balance.routing_flooding += c.nodes.overflow[uj] * dt_routing;
                    c.mass_balance.step_flooding += c.nodes.overflow[uj];
                }
            }

            let c = &mut self.ctx_;
            if quality_outflow > 0.0 && !fv_accepted {
                for p in 0..np as usize {
                    let v = quality_outflow * c.nodes.conc[uj * np as usize + p] * dt_routing;
                    if quality_flooding {
                        c.mass_balance.qual_routing_flood[p] += v;
                    } else {
                        c.mass_balance.qual_routing_outflow[p] += v;
                    }
                }
                book_msx_system_outflow(c, j, quality_outflow, quality_flooding, dt_routing);
            }
            let spill = cmax(0.0, -c.nodes.coupling_inflow[uj]);
            let lat_withdrawal = spill + cmax(0.0, -(c.nodes.lat_flow[uj] - c.nodes.coupling_inflow[uj]));
            if lat_withdrawal > 0.0 && !fv_accepted {
                for p in 0..np as usize {
                    c.mass_balance.qual_routing_outflow[p] += lat_withdrawal * c.nodes.conc[uj * np as usize + p] * dt_routing;
                }
                book_msx_system_outflow(c, j, lat_withdrawal, false, dt_routing);
            }
            let sr = c.node_subtypes.storage_row(j);
            if sr >= 0 {
                let us = sr as usize;
                c.mass_balance.routing_evap_loss += c.node_subtypes.storages.evap_loss[us];
                c.mass_balance.routing_seep_loss += c.node_subtypes.storages.exfil_loss[us];
            }
            step_loss_rate += c.nodes.losses[uj];
        }

        // A street inlet's capture node hands its overflow back (inlet_adjustQualOutflows).
        inlet_adjust_flooding_totals(&self.inlet_, &mut self.ctx_, dt_routing);

        // Link evaporation and seepage losses
        let c = &mut self.ctx_;
        for j in 0..c.n_links() {
            let uj = j as usize;
            if c.links.type_[uj] == LinkType::Conduit {
                let cr = c.link_subtypes.conduit_row(j);
                let cd = &c.link_subtypes.conduits;
                let ucr = cr as usize;
                let barrels = std::cmp::max(if cr >= 0 { cd.barrels[ucr] } else { 1 }, 1);
                let evap_rate = if cr >= 0 { cd.evap_loss_rate[ucr] } else { 0.0 };
                c.mass_balance.routing_evap_loss += evap_rate * barrels as f64 * dt_routing;
                let seep_rate_s = if cr >= 0 { cd.seep_loss_rate[ucr] } else { 0.0 };
                c.mass_balance.routing_seep_loss += cmax(seep_rate_s, 0.0) * barrels as f64 * dt_routing;
                c.mass_balance.routing_link_gw_inflow += cmax(-seep_rate_s, 0.0) * barrels as f64 * dt_routing;
                step_loss_rate += evap_rate * barrels as f64 + seep_rate_s * barrels as f64;
            }
        }

        // Wet weather inflow from the decomposed runoff_inflow array.
        let step_ww_rate;
        {
            let mut runoff_q = 0.0;
            let mut user_q_total = 0.0;
            let mut coupling_out_q = 0.0;
            for uj in 0..c.node_names.size() as usize {
                if c.nodes.runoff_inflow[uj] > 0.0 {
                    runoff_q += c.nodes.runoff_inflow[uj];
                }
                if uj < self.lid_drain_q_interp_.len() && self.lid_drain_q_interp_[uj] > 0.0 {
                    runoff_q += self.lid_drain_q_interp_[uj];
                }
                let q_user = effective_user_lat_flow(c, uj);
                if q_user > 0.0 {
                    user_q_total += q_user;
                }
                if c.nodes.coupling_inflow[uj] < 0.0 {
                    coupling_out_q += -c.nodes.coupling_inflow[uj];
                }
            }
            if runoff_q > 0.0 {
                c.mass_balance.routing_wet_weather += runoff_q * dt_routing;
            }
            step_ww_rate = runoff_q;
            if user_q_total > 0.0 {
                c.mass_balance.routing_forcing_inflow += user_q_total * dt_routing;
            }
            if coupling_out_q > 0.0 {
                // No surface coupling exists in a build without the 2D feature.
                #[cfg(not(feature = "twod"))]
                let as_flooding = false;
                // ctx.twod_io.options->coupling_in_flooding (the router's [2D_OPTIONS])
                #[cfg(feature = "twod")]
                let as_flooding = self.surface_router_.options().borrow().coupling_in_flooding;
                if as_flooding {
                    c.mass_balance.routing_flooding += coupling_out_q * dt_routing;
                } else {
                    c.mass_balance.routing_coupling_out += coupling_out_q * dt_routing;
                }
                c.mass_balance.step_flooding += coupling_out_q;
            }
        }

        // legacy massbal_getStepFlowError on this step's StepFlowTotals
        {
            let mb = &c.mass_balance;
            let total_in = mb.step_dw_inflow + step_ww_rate + mb.step_gw_inflow + mb.step_rdii_inflow + mb.step_ext_inflow;
            let total_out = mb.step_flooding + mb.step_outflow + step_loss_rate;
            let mut err = 0.0;
            if total_in.abs() > 0.0 {
                err = 1.0 - total_out / total_in;
            } else if total_out.abs() > 0.0 {
                err = total_in / total_out - 1.0;
            }
            c.mass_balance.last_step_flow_error = err;
        }
    }

    /// SWMMEngine::computeFinalStorage: closing runoff, snow, routing and quality storage.
    pub(crate) fn compute_final_storage(&mut self) {
        {
            let soa = &self.runoff_.soa_;
            self.ctx_.mass_balance.runoff_final_store = 0.0;
            for i in 0..self.ctx_.n_subcatches() {
                if self.owns_spatial_source(i) {
                    continue;
                }
                let ui = i as usize;
                let fi = soa.imperv_pct[ui];
                let fp = 1.0 - fi;
                let f0 = soa.frac_imperv0[ui];
                let f1 = soa.frac_imperv1[ui];
                let area = soa.area[ui];
                self.ctx_.mass_balance.runoff_final_store += (soa.depth_imperv0[ui] * f0 + soa.depth_imperv1[ui] * f1 + soa.depth_perv[ui] * fp) * area;
            }
            let lv = self.lid_.stored_volume(self.ordinary_source_owners());
            self.ctx_.mass_balance.runoff_final_store += lv;
            self.ctx_.mass_balance.runoff_final_snow = snow_cover_volume_ft3(&self.ctx_, self.snow_.state());
            self.ctx_.mass_balance.runoff_snowremov = self.snow_.state().removed;
        }

        self.ctx_.mass_balance.routing_final_storage = 0.0;
        for j in 0..self.ctx_.n_nodes() {
            let v = self.reported_node_volume(j);
            self.ctx_.mass_balance.routing_final_storage += v;
        }
        for uj in 0..self.ctx_.n_links() as usize {
            self.ctx_.mass_balance.routing_final_storage += self.ctx_.links.volume[uj];
        }

        // Final stored pollutant mass (legacy massbal_getStoredMass).
        let np = self.ctx_.n_pollutants();
        for p in 0..np {
            let up = p as usize;
            if up >= self.ctx_.mass_balance.qual_routing_final.len() {
                break;
            }
            let c = &self.ctx_;
            let mut m = 0.0;
            for j in 0..c.n_nodes() {
                let idx = j as usize * np as usize + up;
                if idx < c.nodes.conc.len() {
                    m += c.nodes.conc[idx] * c.nodes.volume[j as usize] + lidnode::held_mass(c, j, p);
                }
            }
            for uj in 0..c.n_links() as usize {
                let idx = uj * np as usize + up;
                if idx < c.links.conc.len() {
                    m += c.links.conc[idx] * c.links.volume[uj];
                }
            }
            if up < c.mass_balance.qual_routing_final_dry.len() {
                m += c.mass_balance.qual_routing_final_dry[up];
            }
            self.ctx_.mass_balance.qual_routing_final[up] = m;
        }
    }

    /// SWMMEngine::assembleLateralInflows(dt_routing): every decomposed inflow source into
    /// nodes.lat_flow, in legacy routing_execute's source order, and the step's inflow
    /// totals; then the street inlets' capture transfer. stepRouting calls it.
    #[allow(dead_code)] // stepRouting's helper
    pub(crate) fn assemble_lateral_inflows(&mut self, dt_routing: f64) {
        self.assemble_lateral_inflows_sources(dt_routing);
        // Street inlet capture / backflow (legacy inlet_findCapturedFlows).
        self.inlet_.compute_all(&mut self.ctx_, dt_routing);
    }

    /// assembleLateralInflows up to (not including) its last call, inlet_.computeAll; a
    /// separate fn only so the gate's helper mode can stop before the routing step's code.
    #[allow(dead_code)] // stepRouting's helper
    pub(crate) fn assemble_lateral_inflows_sources(&mut self, dt_routing: f64) {
        let mut sum_dw = 0.0;
        let mut sum_gw = 0.0;
        let mut sum_rdii = 0.0;
        let mut sum_ext = 0.0;
        let mut sum_wd = 0.0;
        let np = self.ctx_.n_pollutants() as usize;
        {
            let c = &mut self.ctx_;
            let remaining = c.coupling_delivery_remaining;
            let drain_queue = |qm: &mut f64| -> f64 {
                let mut rate = 0.0;
                if dt_routing > 0.0 && *qm != 0.0 {
                    if remaining > dt_routing {
                        rate = *qm / remaining;
                        *qm -= rate * dt_routing;
                    } else {
                        rate = *qm / dt_routing;
                        *qm = 0.0;
                    }
                }
                rate
            };
            let n = &mut c.nodes;
            for uj in 0..c.node_names.size() as usize {
                let mut q_couple = 0.0;
                if dt_routing > 0.0 && n.coupling_queue[uj] != 0.0 {
                    if remaining > dt_routing {
                        q_couple = n.coupling_queue[uj] / remaining;
                        n.coupling_queue[uj] -= q_couple * dt_routing;
                    } else {
                        q_couple = n.coupling_queue[uj] / dt_routing;
                        n.coupling_queue[uj] = 0.0;
                    }
                }
                if dt_routing > 0.0 && n.coupling_volume[uj] != 0.0 {
                    q_couple += n.coupling_volume[uj] / dt_routing;
                    n.coupling_volume[uj] = 0.0;
                }
                n.coupling_inflow[uj] = q_couple;
                if !n.coupling_qual_queue.is_empty() {
                    for p in 0..np {
                        let idx = uj * np + p;
                        n.coupling_qual_inflow[idx] = drain_queue(&mut n.coupling_qual_queue[idx]);
                    }
                }
                if uj < n.coupling_age_vol_queue.len() {
                    n.coupling_age_vol_inflow[uj] = drain_queue(&mut n.coupling_age_vol_queue[uj]);
                }
                if uj < n.coupling_temp_vol_queue.len() {
                    n.coupling_temp_vol_inflow[uj] = drain_queue(&mut n.coupling_temp_vol_queue[uj]);
                }

                // legacy routing_execute source order (routing.c:466-473)
                n.lat_flow[uj] = n.ext_inflow[uj] + n.dwf_inflow[uj];
                sum_dw += n.dwf_inflow[uj];
                sum_gw += n.gw_inflow[uj];
                sum_rdii += n.rdii_inflow[uj];
                {
                    let q_ext = n.ext_inflow[uj];
                    if q_ext >= 0.0 {
                        sum_ext += q_ext;
                    } else {
                        sum_wd += -q_ext;
                    }
                }
                sum_ext += n.iface_inflow[uj];
                if n.coupling_inflow[uj] > 0.0 {
                    sum_ext += n.coupling_inflow[uj];
                }
            }
            if c.coupling_delivery_remaining > 0.0 {
                c.coupling_delivery_remaining = cmax(0.0, c.coupling_delivery_remaining - dt_routing);
            }
        }

        // Wet-weather runoff and groundwater, per subcatchment in index order.
        {
            let n_sub = self.ctx_.n_subcatches();
            let n_nodes = self.ctx_.n_nodes();
            if self.wet_q_interp_.len() as i32 == n_sub {
                for ui in 0..n_sub as usize {
                    let out_node = self.ctx_.subcatches.outlet_node[ui];
                    if out_node >= 0 && out_node < n_nodes {
                        self.ctx_.nodes.lat_flow[out_node as usize] += self.wet_q_interp_[ui];
                    }
                }
            }
            if self.gw_q_interp_.len() as i32 == n_sub {
                for ui in 0..n_sub as usize {
                    let gwn = self.gw_q_node_[ui];
                    if gwn >= 0 && gwn < n_nodes {
                        self.ctx_.nodes.lat_flow[gwn as usize] += self.gw_q_interp_[ui];
                    }
                }
            }
        }

        // LID drains per unit, then RDII, iface and the engine-only extensions.
        for dq in &self.lid_drain_q_units_ {
            self.ctx_.nodes.lat_flow[dq.0 as usize] += dq.1;
        }
        for uj in 0..self.ctx_.n_nodes() as usize {
            let c = &mut self.ctx_;
            c.nodes.lat_flow[uj] += c.nodes.rdii_inflow[uj];
            c.nodes.lat_flow[uj] += c.nodes.iface_inflow[uj];
            let q_user = effective_user_lat_flow(c, uj);
            if q_user > 0.0 {
                sum_ext += q_user;
            }
            c.nodes.lat_flow[uj] += q_user + c.nodes.coupling_inflow[uj];
        }

        let mb = &mut self.ctx_.mass_balance;
        mb.step_dw_inflow = sum_dw;
        mb.step_gw_inflow = sum_gw;
        mb.step_rdii_inflow = sum_rdii;
        mb.step_ext_inflow = sum_ext;
        mb.step_ext_withdrawal = sum_wd;
    }

    /// SWMMEngine::isBetweenEvents(current_date) (legacy routing.c).
    #[allow(dead_code)] // stepRouting's helper
    pub(crate) fn is_between_events(&self, current_date: f64) -> bool {
        if self.ctx_.events.is_empty() {
            return false;
        }
        if self.next_event_ >= self.ctx_.events.len() as i32 {
            return true;
        }
        let ev = &self.ctx_.events[self.next_event_ as usize];
        if current_date > ev.end {
            return true;
        }
        if current_date >= ev.start {
            return false;
        }
        true
    }

    /// SWMMEngine::isInSteadyState(action_count) (legacy routing.c isInSteadyState +
    /// inflowHasChanged).
    #[allow(dead_code)] // stepRouting's helper
    pub(crate) fn is_in_steady_state(&self, action_count: i32) -> bool {
        // OPENSWMM_HAS_2D: quiescent 1D inflows do not imply a steady surface.
        #[cfg(feature = "twod")]
        if self.surface_router_.is_active() {
            return false;
        }
        let c = &self.ctx_;
        if !c.options.skip_steady_state {
            return false;
        }
        if c.old_elapsed_ms == 0.0 && c.elapsed_ms == 0.0 {
            return false;
        }
        if c.current_time == 0.0 {
            return false;
        }
        if action_count > 0 {
            return false;
        }
        if c.mass_balance.last_step_flow_error.abs() > c.options.sys_flow_tol {
            return false;
        }
        const TINY: f64 = 1e-6;
        for uj in 0..c.n_nodes() as usize {
            let mut q_old = c.nodes.old_lat_flow[uj];
            let mut q_new = c.nodes.lat_flow[uj];
            let mut diff = if q_old.abs() > TINY {
                (q_new / q_old) - 1.0
            } else if q_new.abs() > TINY {
                1.0
            } else {
                0.0
            };
            if diff.abs() > c.options.lat_flow_tol {
                return false;
            }
            if c.nodes.type_[uj] == NodeType::Outfall || c.nodes.degree[uj] == 0 {
                q_old = c.nodes.old_inflow[uj];
                q_new = c.nodes.inflow[uj];
                diff = if q_old.abs() > TINY {
                    (q_new / q_old) - 1.0
                } else if q_new.abs() > TINY {
                    1.0
                } else {
                    0.0
                };
                if diff.abs() > c.options.lat_flow_tol {
                    return false;
                }
            }
        }
        true
    }
}

/// snprintf into a char[256]: at most 255 bytes.
fn truncate255(mut b: Bytes) -> Bytes {
    b.truncate(255);
    b
}

/// Gate only (cargo feature `gate`): the engine members step()'s helpers read beyond the
/// context, as gate/diff/ctxstep's helper mode writes them (SIDE.eng). Not part of the
/// translation.
#[cfg(feature = "gate")]
#[derive(Debug, Default)]
pub struct GateStepSide {
    pub next_event_: i32,
    pub fv_accepted: bool,
    /// runoff_.soa_: area, imperv_pct, frac_imperv0, frac_imperv1, depth_imperv0,
    /// depth_imperv1, depth_perv
    pub runoff: [Vec<f64>; 7],
    pub snow_enabled: Vec<u8>,
    pub snow_wsnow: Vec<f64>,
    pub snow_fw: Vec<f64>,
    pub snow_farea: Vec<f64>,
    pub snow_removed: f64,
    pub buildup: Vec<f64>,
    pub lid_drain_q_interp_: Vec<f64>,
    pub wet_q_interp_: Vec<f64>,
    pub gw_q_interp_: Vec<f64>,
    pub gw_q_node_: Vec<i32>,
    pub lid_drain_q_units_: Vec<(i32, f64)>,
    /// per LID group: count, surf_depth, soil_moist, stor_depth, pave_depth
    pub lid_groups: Vec<(i32, [Vec<f64>; 4])>,
}

#[cfg(feature = "gate")]
impl SwmmEngine {
    /// Install the side members (see GateStepSide).
    pub fn gate_install_step_side(&mut self, s: GateStepSide) {
        self.next_event_ = s.next_event_;
        self.quality_.fv_transport_.accepted_ = s.fv_accepted;
        let [a, ip, f0, f1, d0, d1, dp] = s.runoff;
        let r = &mut self.runoff_.soa_;
        r.area = a;
        r.imperv_pct = ip;
        r.frac_imperv0 = f0;
        r.frac_imperv1 = f1;
        r.depth_imperv0 = d0;
        r.depth_imperv1 = d1;
        r.depth_perv = dp;
        let sn = &mut self.snow_.soa_;
        sn.enabled = s.snow_enabled;
        sn.wsnow = s.snow_wsnow;
        sn.fw = s.snow_fw;
        sn.fArea = s.snow_farea;
        sn.removed = s.snow_removed;
        self.surface_quality_.buildup = s.buildup;
        self.lid_drain_q_interp_ = s.lid_drain_q_interp_;
        self.wet_q_interp_ = s.wet_q_interp_;
        self.gw_q_interp_ = s.gw_q_interp_;
        self.gw_q_node_ = s.gw_q_node_;
        self.lid_drain_q_units_ = s.lid_drain_q_units_;
        for (t, (count, [sd, sm, st, pd])) in s.lid_groups.into_iter().enumerate() {
            if let Some(g) = self.lid_.groups_.get_mut(t) {
                g.count = count;
                g.surf_depth = sd;
                g.soil_moist = sm;
                g.stor_depth = st;
                g.pave_depth = pd;
            }
        }
    }

    /// Call one of step()'s helpers by the name gate/diff/ctxstep's helper mode uses; the
    /// lines it prints beyond the context dump.
    pub fn gate_step_helper(&mut self, name: &str, dt: f64) -> Result<Bytes, Terminated> {
        let mut o = Bytes::new();
        match name {
            "stats" => self.update_statistics(dt),
            "massbal" => self.update_routing_mass_balance(dt),
            "finalstore" => self.compute_final_storage(),
            "finalqual" => self.compute_final_quality_mass_balance(),
            "nodetotals" => self.accumulate_node_routing_totals(dt / 2.0),
            "diverged" => {
                if self.router_.fv_solver_.is_some() {
                    o.extend_from_slice(b"diverged skipped (FV)\n");
                    return Ok(o);
                }
                let rc = self.check_routing_diverged()?;
                o.extend_from_slice(format!("diverged rc {}\n", rc).as_bytes());
            }
            "events" => {
                let b = self.is_between_events(self.ctx_.current_date);
                let st = self.is_in_steady_state(0);
                o.extend_from_slice(format!("between {} steady {}\n", b as i32, st as i32).as_bytes());
            }
            "lateral" => {
                if self.inlet_.soa_.count > 0 {
                    o.extend_from_slice(b"lateral skipped (inlets)\n");
                    return Ok(o);
                }
                // computeAll returns at once without inlets; its stub would panic, so stop
                // short of it: everything before it is assemble_lateral_inflows' own work
                self.assemble_lateral_inflows_sources(dt);
            }
            "forcings" => {
                self.reset_step_mass_balance();
                self.apply_forcings(dt)?;
            }
            _ => panic!("unknown helper {}", name),
        }
        Ok(o)
    }
}
