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
// the non-conduit (structure) callback, the routing statistics' iteration count, and the
// quality stage: the inlet quality transfer (InletSolver::adjustQualInflows, B4b) and the
// LEGACY quality router (QualitySolver::execute, B5; quality/quality_routing.rs).
//
// Differences in form, never in result:
//   * stepRouting can throw a std::exception, which step() catches (set_error with
//     SWMM_ERR_WRONG_STATE and e.what()). Here it returns Err(StepError::Exception(what)).
//   * The non-conduit callback, a lambda capturing the engine upstream, is a closure over
//     nothing: the solver hands it the DW node arrays (NcView) and the StructureSolver.
//   * The SWMM_TRACE_RSTEP / SWMM_TRACE_DUMP_STEP trace prints only and is left out.
//   * OPENSWMM_HAS_2D is off in the reference build; its blocks are absent here too.
//
// NOT YET TRANSLATED (StepError::NotTranslated, returned before the step changes anything:
// step_routing_not_translated): the EULERIAN_ARD and LAGRANGIAN quality engines (ArdEngine,
// LagrangianSolver), the reaction solvers ROS2 and BDF2, EQUIL species, booster sources and
// reaction-species LID controls (msxlid), the heat flux modules (surface, radiative and
// bed exchange), engine checkpointing, FV quality transport of reaction species, heat
// or bed exchange, the FV solver options
// Router::fv_not_translated names,
// node LID Richards columns, and what the solvers report as not translated (DYNAMIC_SLOT,
// TPA, physical junction areas, control-rule expression premises, checkpoint sessions).

use crate::core::context::SimulationContext;
use crate::core::datetime;
use crate::core::engine_init_hydraulics::Terminated;
use crate::core::engine_step::STUB_EXCEPTION;
use crate::core::engine::SwmmEngine;
use crate::core::options::{QualitySolverKind, RoutingModel};
use crate::data::link::LinkType;
use crate::hydraulics::dynamic_wave_step::{bound_surface_outfall_flow, NcView};
use crate::hydraulics::hyd_structures::StructureSolver;
use crate::hydrology::lid_node_step as lidnode;
use crate::transport::reaction_legacy_binding::legacy_reactions_active;

/// What the routing step's solvers can end with besides success (converted to the crate's
/// Terminated at step_routing's boundary).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StepError {
    /// A std::runtime_error upstream throws; step() catches it and reports what().
    Exception(&'static str),
    /// A path the translation does not cover yet (never upstream behavior).
    NotTranslated(&'static str),
}

/// A std::runtime_error on its way to step()'s catch (engine_step's THROWN_PREFIX + the type).
pub const THROWN_RUNTIME_ERROR: &str = "throw std::runtime_error";

impl From<StepError> for Terminated {
    fn from(e: StepError) -> Terminated {
        match e {
            StepError::Exception(what) => Terminated { exception: THROWN_RUNTIME_ERROR, what },
            StepError::NotTranslated(what) => Terminated { exception: STUB_EXCEPTION, what },
        }
    }
}

impl SwmmEngine {
    /// The paths stepRouting does not translate yet, checked before anything changes.
    fn step_routing_not_translated(&self) -> Option<&'static str> {
        let ctx = &self.ctx_;
        if self.checkpoint_enabled_ {
            return Some("engine checkpointing");
        }
        if ctx.options.routing_model == RoutingModel::Fv {
            if let Some(w) = self.router_.fv_not_translated(ctx) {
                return Some(w);
            }
        }
        if lidnode::unsupported(ctx) {
            return Some("Richards node LID columns");
        }
        // B5's gate: the quality stage runs this step
        if (ctx.n_pollutants() > 0 || legacy_reactions_active(ctx) || ctx.options.water_age || ctx.options.heat_transport) && !ctx.options.ignore_quality {
            match ctx.options.quality_solver {
                QualitySolverKind::EulerianArd => return Some("EULERIAN_ARD quality engine"),
                QualitySolverKind::Lagrangian => return Some("LAGRANGIAN quality engine"),
                _ => {}
            }
            if let Some(w) = crate::transport::reaction_legacy_step::reactions_not_translated(ctx) {
                return Some(w);
            }
            if crate::transport::heat_legacy::heat_fluxes_active(ctx) {
                return Some("heat fluxes");
            }
            if ctx.options.routing_model == RoutingModel::Fv && ctx.options.quality_solver == QualitySolverKind::Legacy && (ctx.reactions.n_species() > 0 || ctx.options.heat_transport || legacy_reactions_active(ctx) || crate::transport::heat_legacy::bed_exchange_enabled(ctx)) {
                return Some("FV quality transport of reaction species, heat or bed exchange");
            }
        }
        None
    }

    /// SWMMEngine::stepRouting: phase B of step(), one routing step of dt_routing seconds.
    /// Err(Terminated) carries either a std::runtime_error step() catches (exception
    /// THROWN_RUNTIME_ERROR) or a path not translated yet (STUB_EXCEPTION).
    pub fn step_routing(&mut self, dt_routing: f64) -> Result<(), Terminated> {
        self.step_routing_impl(dt_routing).map_err(Terminated::from)
    }

    fn step_routing_impl(&mut self, dt_routing: f64) -> Result<(), StepError> {
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
        // BW-MSX: this step's species washoff joins the external species loads
        if !self.ctx_.options.ignore_quality {
            crate::transport::msx_surface::deliver(&mut self.ctx_, dt_routing);
        }
        if !self.ctx_.options.ignore_quality {
            drain_msx_coupling_queue(&mut self.ctx_, dt_routing);
        }

        // B2a. RDII inflows
        if !self.ctx_.options.ignore_rdii {
            if self.rdii_iface_file_.is_open() && !self.rdii_iface_file_.is_writing() {
                self.rdii_iface_file_.apply_flows(&mut self.ctx_, routing_date);
            } else {
                self.rdii_.apply_rdii_inflows(&mut self.ctx_, routing_secs);
            }
        }

        // B2b. Interface file inflows (legacy addIfaceInflows(currentDate), the same +1 ms date)
        self.iface_.read_inflows(&mut self.ctx_, routing_date);

        // B2c. Assemble the lateral inflows
        self.assemble_lateral_inflows(dt_routing);

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
        // dw.nodeNewSurfAreaDataMut(): null for the DW solver's empty arrays (FV routing,
        // where that solver is never initialized).
        fn surf_buf(v: &mut Vec<f64>) -> Option<&mut [f64]> {
            if v.is_empty() {
                None
            } else {
                Some(&mut v[..])
            }
        }
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
                        // Each ordinary structure sees the head at its own layer port.
                        let (ua, ub) = (a as usize, b as usize);
                        let da = ctx.nodes.depth[ua];
                        let db = ctx.nodes.depth[ub];
                        let ha = ctx.nodes.head[ua];
                        let hb = ctx.nodes.head[ub];
                        let va = ctx.nodes.old_volume[ua];
                        let vb = ctx.nodes.old_volume[ub];
                        ctx.nodes.depth[ua] = lidnode::port_depth(ctx, a, lidnode::port_offset(ctx, j, a));
                        ctx.nodes.depth[ub] = lidnode::port_depth(ctx, b, lidnode::port_offset(ctx, j, b));
                        ctx.nodes.head[ua] = ctx.nodes.invert_elev[ua] + ctx.nodes.depth[ua];
                        ctx.nodes.head[ub] = ctx.nodes.invert_elev[ub] + ctx.nodes.depth[ub];
                        ctx.nodes.old_volume[ua] += lidnode::held_volume(ctx, a);
                        ctx.nodes.old_volume[ub] += lidnode::held_volume(ctx, b);
                        hydstruct.compute_non_conduit_flow_one(ctx, dt, surf_buf(&mut dw.new_surf_area), j);
                        ctx.nodes.depth[ua] = da;
                        ctx.nodes.depth[ub] = db;
                        ctx.nodes.head[ua] = ha;
                        ctx.nodes.head[ub] = hb;
                        ctx.nodes.old_volume[ua] = va;
                        ctx.nodes.old_volume[ub] = vb;
                    } else {
                        hydstruct.compute_non_conduit_flow_one(ctx, dt, surf_buf(&mut dw.new_surf_area), j);
                    }

                    let mut q_new = ctx.links.flow[uj];
                    if step > 0 && ctx.links.type_[uj] != LinkType::Pump {
                        q_new = (1.0 - 0.5) * q_last + 0.5 * q_new;
                        if q_new * q_last < 0.0 {
                            q_new = 0.001 * (if q_new >= 0.0 { 1.0 } else { -1.0 });
                        }
                    }
                    ctx.links.flow[uj] = q_new;
                } else {
                    hydstruct.scatter_held_surf_area(ctx, surf_buf(&mut dw.new_surf_area), j);
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

        // B3-pre. Exfiltration (storage node Green-Ampt seepage), before the router so the
        //         raw rates reach Router::initNodeFlows' joint evap+exfil cap
        self.exfil_.compute_all(&mut self.ctx_, dt_routing);

        let evap_rate = self.ctx_.climate_state.evap_rate;
        let iters = self.router_.step(&mut self.ctx_, dt_routing, evap_rate, &mut self.hydstruct_, Some(&mut non_conduit_fn))?;
        if self.ctx_.options.routing_model == RoutingModel::Fv && self.ctx_.options.quality_solver == QualitySolverKind::Legacy {
            self.router_.accepted_quality_water(&self.ctx_, dt_routing, &mut self.quality_.fv_transport_.water);
        }
        lidnode::finish_step(&mut self.ctx_)?;
        let converged = self.router_.last_step_converged();
        self.ctx_.routing_stats.update_iterations(iters, converged);

        // B4b. Gap #55: inlet quality mass transfer (bypass <-> capture on the net flow).
        //      (U92: B5's assembleExternalLoads zeroes what it adds.)
        if self.ctx_.n_pollutants() > 0 && !self.ctx_.options.ignore_quality {
            self.inlet_.adjust_qual_inflows(&mut self.ctx_, dt_routing);
        }

        // B5. Water quality routing
        if (self.ctx_.n_pollutants() > 0 || legacy_reactions_active(&self.ctx_) || self.ctx_.options.water_age || self.ctx_.options.heat_transport) && !self.ctx_.options.ignore_quality {
            match self.ctx_.options.quality_solver {
                // refused in step_routing_not_translated
                QualitySolverKind::EulerianArd => return Err(StepError::NotTranslated("EULERIAN_ARD quality engine")),
                QualitySolverKind::Lagrangian => return Err(StepError::NotTranslated("LAGRANGIAN quality engine")),
                _ => self.quality_.execute(&mut self.ctx_, dt_routing)?,
            }
        }
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

