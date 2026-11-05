// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/SWMMEngine.cpp (stepRouting, isBetweenEvents,
// isInSteadyState) and SWMMEngine.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit
// 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// SWMMEngine::stepRouting: one routing step (phase B of step()): events, pump hysteresis,
// control rules, setting transitions, external / DWF / RDII / interface inflows, lateral
// inflow assembly, the steady-state skip, the hydraulic routing step through the Router with
// the non-conduit (structure) callback, and the routing statistics' iteration count.
//
// Differences in form, never in result:
//   * stepRouting can throw a std::exception, which step() catches (set_error with
//     SWMM_ERR_WRONG_STATE and e.what()). Here it returns Err(StepError::Exception(what)).
//   * The non-conduit callback, a lambda capturing the engine upstream, is a closure over
//     nothing: the solver hands it the DW node arrays (NcView) and the StructureSolver.
//   * The SWMM_TRACE_RSTEP / SWMM_TRACE_DUMP_STEP trace prints only and is left out.
//   * OPENSWMM_HAS_2D is off in the reference build; its blocks are absent here too.
//
// NOT YET TRANSLATED (StepError::NotTranslated, returned before the hydraulic step changes
// anything it would not also change upstream): the quality stage (QualitySolver::execute,
// the EULERIAN_ARD and LAGRANGIAN engines, inlet quality transfer), species surface washoff
// delivery (msxsurf::deliver for an active surface), [FILES] USE RDII
// (RdiiInterfaceFile::applyFlows), routing interface inflow files (InterfaceManager
// ::readInflows with a file open), engine checkpointing, the FV router, and what the
// solvers report as not translated (KINWAVE / STEADY, DYNAMIC_SLOT, TPA, virtual junctions,
// physical junction areas, node LID storage, storage exfiltration, street inlets, control
// expressions).

use crate::core::context::SimulationContext;
use crate::core::datetime;
use crate::core::engine::SwmmEngine;
use crate::core::options::{Bytes, QualitySolverKind, RoutingModel};
use crate::data::link::LinkType;
use crate::data::node::NodeType;
use crate::data::forcing::ForcingMode;
use crate::hydraulics::dynamic_wave_step::{bound_surface_outfall_flow, NcView};
use crate::hydraulics::hyd_structures::StructureSolver;
use crate::hydrology::lid_node_step as lidnode;
use crate::transport::reaction_legacy_binding::legacy_reactions_active;

/// What stepRouting can end with besides success.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StepError {
    /// A std::exception upstream throws; step() catches it and reports what().
    Exception(Bytes),
    /// A path the translation does not cover yet (never upstream behavior).
    NotTranslated(&'static str),
}

impl SwmmEngine {
    /// isBetweenEvents (legacy routing.c isBetweenEvents)
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

    /// isInSteadyState (legacy routing.c isInSteadyState + inflowHasChanged)
    pub(crate) fn is_in_steady_state(&self, action_count: i32) -> bool {
        let ctx = &self.ctx_;
        if !ctx.options.skip_steady_state {
            return false;
        }
        if ctx.old_elapsed_ms == 0.0 && ctx.elapsed_ms == 0.0 {
            return false;
        }
        if ctx.current_time == 0.0 {
            return false;
        }
        if action_count > 0 {
            return false;
        }
        if ctx.mass_balance.last_step_flow_error.abs() > ctx.options.sys_flow_tol {
            return false;
        }
        const TINY: f64 = 1e-6;
        for j in 0..ctx.n_nodes() as usize {
            let mut q_old = ctx.nodes.old_lat_flow[j];
            let mut q_new = ctx.nodes.lat_flow[j];
            let mut diff = if q_old.abs() > TINY {
                (q_new / q_old) - 1.0
            } else if q_new.abs() > TINY {
                1.0
            } else {
                0.0
            };
            if diff.abs() > ctx.options.lat_flow_tol {
                return false;
            }
            if ctx.nodes.type_[j] == NodeType::Outfall || ctx.nodes.degree[j] == 0 {
                q_old = ctx.nodes.old_inflow[j];
                q_new = ctx.nodes.inflow[j];
                diff = if q_old.abs() > TINY {
                    (q_new / q_old) - 1.0
                } else if q_new.abs() > TINY {
                    1.0
                } else {
                    0.0
                };
                if diff.abs() > ctx.options.lat_flow_tol {
                    return false;
                }
            }
        }
        true
    }

    /// The paths stepRouting does not translate yet, checked before anything changes.
    fn step_routing_not_translated(&self) -> Option<&'static str> {
        let ctx = &self.ctx_;
        if self.checkpoint_enabled_ {
            return Some("engine checkpointing");
        }
        if ctx.options.routing_model == RoutingModel::Fv {
            return Some("FV routing");
        }
        if lidnode::unsupported(ctx) {
            return Some("node LID columns");
        }
        if !ctx.options.ignore_quality && ctx.reactions.surface.active() {
            return Some("species surface washoff delivery");
        }
        if !ctx.options.ignore_rdii && self.rdii_iface_open_ && !self.rdii_iface_writing_ {
            return Some("[FILES] USE RDII");
        }
        if self.iface_.infile_.is_some() {
            return Some("routing interface inflow files");
        }
        if self.exfil_.soa_.count != 0 {
            return Some("storage exfiltration");
        }
        if self.inlet_.soa_.count != 0 {
            return Some("street inlets");
        }
        if (ctx.n_pollutants() > 0 || legacy_reactions_active(ctx) || ctx.options.water_age || ctx.options.heat_transport) && !ctx.options.ignore_quality {
            return Some("quality routing");
        }
        None
    }

    /// SWMMEngine::stepRouting: phase B of step(), one routing step of dt_routing seconds.
    pub fn step_routing(&mut self, dt_routing: f64) -> Result<(), StepError> {
        if let Some(w) = self.step_routing_not_translated() {
            return Err(StepError::NotTranslated(w));
        }
        // quality_.beginFvWindow()
        self.quality_.fv_transport_.accepted_ = false;
        self.quality_.fv_transport_.water.ready = false;
        self.last_step_steady_ = false;

        // The step's date: its start plus legacy's 1 ms, from the ms clock.
        let routing_secs = (self.ctx_.elapsed_ms + 1.0) / 1000.0;
        let routing_date = datetime::add_seconds(self.ctx_.options.start_date, routing_secs);

        // B0a. Between routing events: skip routing
        self.between_events_ = self.is_between_events(routing_date);
        if self.between_events_ {
            self.ctx_.nodes.lat_flow.iter_mut().for_each(|x| *x = 0.0);
            if !self.ctx_.options.ignore_quality {
                self.ctx_.nodes.conc.iter_mut().for_each(|x| *x = 0.0);
                self.ctx_.links.conc.iter_mut().for_each(|x| *x = 0.0);
            }
            if self.next_event_ < self.ctx_.events.len() as i32 && routing_date > self.ctx_.events[self.next_event_ as usize].end {
                self.next_event_ += 1;
            }
            return Ok(());
        }

        // B1a. Pump startup/shutoff hysteresis, once per step
        self.hydstruct_.update_pump_target_settings(&mut self.ctx_);

        // B1b. Control rules
        let rule_time_reached = self.ctx_.options.rule_step <= 0.0 || (self.ctx_.elapsed_ms - self.new_rule_time_ms_).abs() < 1.0;
        self.controls_.evaluate(&mut self.ctx_, routing_date, dt_routing, rule_time_reached)?;
        // Setting transitions toward target_setting
        for j in 0..self.ctx_.n_links() {
            let uj = j as usize;
            let target = self.ctx_.links.target_setting[uj];
            let current = self.ctx_.links.setting[uj];
            if target != current {
                if target * current == 0.0 {
                    self.ctx_.links.time_last_set[uj] = routing_date;
                }
                let orr = self.ctx_.link_subtypes.orifice_row(j);
                let orate_sec = (if orr >= 0 { self.ctx_.link_subtypes.orifices.orate[orr as usize] } else { 0.0 }) * 3600.0;
                if orate_sec > 0.0 && dt_routing > 0.0 {
                    let step = dt_routing / orate_sec;
                    let delta_val = (target - current).abs();
                    if step + 0.001 >= delta_val {
                        self.ctx_.links.setting[uj] = target;
                    } else if target > current {
                        self.ctx_.links.setting[uj] = current + step;
                    } else {
                        self.ctx_.links.setting[uj] = current - step;
                    }
                } else {
                    self.ctx_.links.setting[uj] = target;
                }
            }
        }

        // B2. Inflow sources into the decomposed buffers
        self.ctx_.nodes.clear_inflow_sources();
        self.inflow_.compute_all(&mut self.ctx_, routing_date, dt_routing);
        // (msxsurf::deliver: a no-op without an active species surface, refused above)
        if !self.ctx_.options.ignore_quality {
            drain_msx_coupling_queue(&mut self.ctx_, dt_routing);
        }

        // B2a. RDII inflows
        if !self.ctx_.options.ignore_rdii {
            // (the RDII interface file's applyFlows: refused above)
            apply_rdii_inflows(&self.rdii_, &mut self.ctx_, routing_secs);
        }

        // B2b. Interface file inflows (no file open: refused above otherwise)

        // B2c. Assemble the lateral inflows
        self.assemble_lateral_inflows(dt_routing)?;

        // B2d. Steady state skip
        let action_count = self.controls_.last_action_count();
        let has_lid_nodes = self.ctx_.node_subtypes.storages.lid.iter().any(|c| c.control >= 0);
        if !has_lid_nodes && self.is_in_steady_state(action_count) {
            self.last_step_steady_ = true;
            return Ok(());
        }
        let evap_now = self.ctx_.climate_state.evap_rate;
        lidnode::prepare_step(&mut self.ctx_, dt_routing, evap_now)?;
        self.ctx_.save_hyd_state();

        // B3. Hydraulic routing, with the non-conduit flows computed inside the Picard loop
        let mut non_conduit_fn = |ctx: &mut SimulationContext, dt: f64, step: i32, dw: &mut NcView, hydstruct: &mut StructureSolver| -> Result<(), crate::core::engine_step_routing::StepError> {
            for k in 0..hydstruct.nc_indices_.len() {
                let j = hydstruct.nc_indices_[k];
                let uj = j as usize;
                let q_last = ctx.links.flow[uj];
                let bypassed = dw.is_bypassed(j);

                if !bypassed {
                    let a = ctx.links.node1[uj];
                    let b = ctx.links.node2[uj];
                    if a >= 0 && b >= 0 && (lidnode::active(ctx, a) || lidnode::active(ctx, b)) {
                        return Err(StepError::NotTranslated("node LID ports"));
                    }
                    hydstruct.compute_non_conduit_flow_one(ctx, dt, Some(&mut dw.new_surf_area[..]), j);

                    let mut q_new = ctx.links.flow[uj];
                    if step > 0 && ctx.links.type_[uj] != LinkType::Pump {
                        q_new = (1.0 - 0.5) * q_last + 0.5 * q_new;
                        if q_new * q_last < 0.0 {
                            q_new = 0.001 * (if q_new >= 0.0 { 1.0 } else { -1.0 });
                        }
                    }
                    ctx.links.flow[uj] = q_new;
                } else {
                    hydstruct.scatter_held_surf_area(ctx, Some(&mut dw.new_surf_area[..]), j);
                }

                let before_surface = ctx.links.flow[uj];
                ctx.links.flow[uj] = bound_surface_outfall_flow(ctx, j, before_surface);
                if ctx.links.flow[uj] != before_surface {
                    ctx.links.dqdh[uj] = 0.0;
                }

                let before_port = ctx.links.flow[uj];
                ctx.links.flow[uj] = lidnode::exchange_ports(ctx, j, before_port, dt)?;
                if before_port != 0.0 && ctx.links.flow[uj] != before_port {
                    ctx.links.dqdh[uj] *= ctx.links.flow[uj] / before_port;
                }

                let q = ctx.links.flow[uj];
                let n1 = ctx.links.node1[uj];
                let n2 = ctx.links.node2[uj];
                if n1 < 0 || n2 < 0 {
                    continue;
                }
                let un1 = n1 as usize;
                let un2 = n2 as usize;
                if q >= 0.0 {
                    ctx.nodes.outflow[un1] += q;
                    ctx.nodes.inflow[un2] += q;
                } else {
                    ctx.nodes.inflow[un1] -= q;
                    ctx.nodes.outflow[un2] -= q;
                }

                if dw.initialized {
                    let dqdh = ctx.links.dqdh[uj];
                    let pr_t4 = ctx.link_subtypes.pump_row(j);
                    let is_type4_pump = ctx.links.type_[uj] == LinkType::Pump && pr_t4 >= 0 && ctx.link_subtypes.pumps.curve_type[pr_t4 as usize] == 4;
                    dw.sumdqdh[un1] += dqdh;
                    if !is_type4_pump {
                        dw.sumdqdh[un2] += dqdh;
                    }
                }
            }
            Ok(())
        };

        // B3-pre. Exfiltration (none: refused above otherwise)

        let evap_rate = self.ctx_.climate_state.evap_rate;
        let iters = self.router_.step(&mut self.ctx_, dt_routing, evap_rate, &mut self.hydstruct_, Some(&mut non_conduit_fn))?;
        lidnode::finish_step(&mut self.ctx_)?;
        let converged = self.router_.last_step_converged();
        self.ctx_.routing_stats.update_iterations(iters, converged);

        // B4b / B5. Quality (refused above when it would run)
        let _ = (ForcingMode::None, QualitySolverKind::Legacy);
        Ok(())
    }
}

/// transport::drainMsxCouplingQueue: species mass the 2D surface and the aquifer handed to
/// nodes, at the coupled volume's delivery rate.
fn drain_msx_coupling_queue(ctx: &mut SimulationContext, dt: f64) {
    let nsp = ctx.reactions.n_species();
    if ctx.reactions.msx_coupling_queue.is_empty() || nsp <= 0 || !(dt > 0.0) {
        return;
    }
    let want = (ctx.n_nodes() as usize) * (nsp as usize);
    if ctx.reactions.msx_coupling_queue.len() != want {
        return;
    }
    if ctx.reactions.msx_ext_mass_in.len() != want {
        ctx.reactions.msx_ext_mass_in.clear();
        ctx.reactions.msx_ext_mass_in.resize(want, 0.0);
    }
    let remaining = ctx.coupling_delivery_remaining;
    let rx = &mut ctx.reactions;
    for k in 0..want {
        let qm = rx.msx_coupling_queue[k];
        if qm == 0.0 {
            continue;
        }
        let rate;
        if remaining > dt {
            rate = qm / remaining;
            rx.msx_coupling_queue[k] -= rate * dt;
        } else {
            rate = qm / dt;
            rx.msx_coupling_queue[k] = 0.0;
        }
        rx.msx_ext_mass_in[k] += rate;
    }
}

/// RDIISolver::applyRdiiInflows: the grid row at the step's date into rdii_inflow.
fn apply_rdii_inflows(rdii: &crate::hydrology::rdii::RDIISolver, ctx: &mut SimulationContext, elapsed_sec: f64) {
    if rdii.grid_step_ <= 0.0 || rdii.grid_node_.is_empty() {
        return;
    }
    if elapsed_sec < 0.0 {
        return;
    }
    let row = (elapsed_sec / rdii.grid_step_).floor() as usize;
    let ncol = rdii.grid_node_.len();
    if (row + 1) * ncol > rdii.grid_flows_.len() {
        return;
    }
    let flows = &rdii.grid_flows_[row * ncol..];
    for col in 0..ncol {
        let q = flows[col] as f64;
        if q.abs() < crate::core::constants::FLOW_TOL {
            continue;
        }
        let un = rdii.grid_node_[col] as usize;
        if un < ctx.nodes.rdii_inflow.len() {
            ctx.nodes.rdii_inflow[un] += q;
        }
    }
}
