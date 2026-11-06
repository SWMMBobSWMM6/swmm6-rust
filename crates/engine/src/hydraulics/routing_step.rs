// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Routing.cpp and Routing.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The Router's routing step (routing.rs holds its state and init): step/stepImpl,
// getAdaptiveStep, accumulateLinkStepStats, initNodeFlows, computeConduitLosses,
// executeSteadyFlow and updateLinkStates.
//
// Differences in form, never in result:
//   * `structures_` (upstream a pointer to the engine's StructureSolver, set at init) is passed
//     to step() as an argument, and handed on to the non-conduit callback the same way.
//
// NOT YET TRANSLATED (step returns StepError::NotTranslated): the FV router (stepFv,
// acceptedQualityWater, the FV link-step sink), checkpoint sessions (checkCheckpointStep),
// Mesh Richards routing.

use crate::core::cmath::{cmax, cmin};
use crate::core::constants;
use crate::core::context::SimulationContext;
use crate::core::engine_step_routing::StepError;
use crate::core::options::RoutingModel;
use crate::core::unit_conversion as ucf;
use crate::data::link::LinkType;
use crate::data::node::NodeType;
use crate::hydraulics::dynamic_wave::build_xsp;
use crate::hydraulics::dynamic_wave_step::{apply_seepage_forcing, aquifer_seep_factor, cap_aquifer_gain, NonConduitFn};
use crate::hydraulics::hyd_structures::StructureSolver;
use crate::hydraulics::node;
use crate::hydraulics::outfall;
use crate::hydraulics::routing::{RouteModel, Router};
use crate::hydraulics::transect_types::N_TRANSECT_TBL;
use crate::hydraulics::xsect_batch::{XSectParams, XSectShape};
use crate::hydraulics::xsection as xsect;
use crate::hydrology::lid_node_step as lidnode;

impl Router {
    /// Router::step (stepImpl with owned_callback false).
    pub fn step(&mut self, ctx: &mut SimulationContext, dt: f64, evap_rate: f64, structures: &mut StructureSolver, non_conduit_fn: Option<&mut NonConduitFn>) -> Result<i32, StepError> {
        let mesh_richards = lidnode::validate_richards_routing(ctx, dt)?;
        if mesh_richards && self.model_ != RouteModel::Dynwave {
            return Err(StepError::Exception("Joint Mesh Richards routing is qualified only for native dynamic wave."));
        }
        if self.checkpoint_context_attached_ {
            return Err(StepError::NotTranslated("checkpoint sessions"));
        }
        if self.model_ == RouteModel::Fv {
            return Err(StepError::NotTranslated("FV routing"));
        }

        // 2. Init node flows from laterals and losses
        self.init_node_flows(ctx, dt, evap_rate);

        // 2b. Conduit evaporation and seepage (DYNWAVE does it per Picard iteration)
        if self.model_ != RouteModel::Dynwave {
            self.compute_conduit_losses(ctx, dt, evap_rate);
        }

        // 3. Outfall boundary depths at the START of the step
        outfall::set_all_outfall_depths(ctx, 0.0);

        // 4. Dispatch to solver
        let iters = match self.model_ {
            RouteModel::Kinwave => self.kw_solver_.execute(ctx, dt, structures)?,
            RouteModel::Dynwave => {
                self.dw_solver_.evap_rate = evap_rate;
                let Router { dw_solver_, groups_, .. } = self;
                dw_solver_.execute(groups_, ctx, dt, structures, non_conduit_fn)?
            }
            RouteModel::Fv => unreachable!(),
            RouteModel::Steady => self.execute_steady_flow(ctx, dt, structures)?,
        };

        // 6. Update link final states (not under FV)
        if ctx.options.routing_model != RoutingModel::Fv {
            self.update_link_states(ctx);
        }
        // (mesh_richards is false: refused above)
        Ok(iters)
    }

    /// getAdaptiveStep
    pub fn get_adaptive_step(&mut self, ctx: &mut SimulationContext, fixed_step: f64, courant: f64) -> f64 {
        if ctx.n_links() == 0 {
            return fixed_step;
        }
        if self.model_ == RouteModel::Dynwave {
            return self.dw_solver_.get_routing_step(ctx, fixed_step, courant);
        }
        // (FV: suggested_step of the FV backend, not translated; the FV router refuses to step)
        fixed_step
    }

    /// accumulateLinkStepStats ([REPORT] LINK_STEPS). The FV branch (folding the FV solver's
    /// pending samples) is not translated: an FV run never gets here, as step_routing refuses
    /// the FV router before.
    pub fn accumulate_link_step_stats(&mut self, ctx: &mut SimulationContext, dt: f64, in_window: bool) {
        if !self.lstep_enabled_ {
            return;
        }
        if self.model_ == RouteModel::Dynwave && in_window {
            let (rs, vs) = (ctx.options.routing_step, ctx.options.variable_step);
            self.dw_solver_.accumulate_link_step_stats(ctx, dt, rs, vs);
        }
    }

    /// initNodeFlows: storage losses, then inflow = lateral flow, outflow = losses, overflow.
    pub fn init_node_flows(&mut self, ctx: &mut SimulationContext, dt: f64, evap_rate: f64) {
        let n = ctx.n_nodes();
        let unit_sys = ucf::get_unit_system(ctx.options.flow_units as i32);
        for i in 0..n {
            let ui = i as usize;
            let mut loss_rate = 0.0;
            if ctx.nodes.type_[ui] == NodeType::Storage {
                let sr = ctx.node_subtypes.storage_row(i);
                let sru = sr as usize;
                let evap_frac = if sr >= 0 { ctx.node_subtypes.storages.evap_frac[sru] } else { 0.0 };
                let seep_rate = if sr >= 0 { ctx.node_subtypes.storages.seep_rate[sru] } else { 0.0 };
                let stor_evap_rate = evap_rate * evap_frac;
                // (LID cells at a storage node: refused by stepRouting)

                let mut exfil_cfs = if sr >= 0 { ctx.node_subtypes.storages.exfil_rate[sru] } else { 0.0 };

                if stor_evap_rate > 0.0 || seep_rate > 0.0 || exfil_cfs > 0.0 {
                    let depth = ctx.nodes.depth[ui];
                    let area = node::get_surf_area(&ctx.nodes, i, depth, Some(&ctx.tables), unit_sys, Some(&ctx.node_subtypes));
                    let mut evap_cfs = stor_evap_rate;
                    if ctx.nodes.volume[ui] > constants::FUDGE {
                        evap_cfs = area * stor_evap_rate;
                    }
                    let total_loss = (evap_cfs + exfil_cfs) * dt;
                    if total_loss > ctx.nodes.volume[ui] && total_loss > 0.0 {
                        let ratio = ctx.nodes.volume[ui] / total_loss;
                        evap_cfs *= ratio;
                        exfil_cfs *= ratio;
                    }
                    if sr >= 0 {
                        ctx.node_subtypes.storages.evap_loss[sru] = evap_cfs * dt;
                        ctx.node_subtypes.storages.exfil_loss[sru] = exfil_cfs * dt;
                    }
                    loss_rate = evap_cfs + exfil_cfs;
                } else if sr >= 0 {
                    ctx.node_subtypes.storages.evap_loss[sru] = 0.0;
                    ctx.node_subtypes.storages.exfil_loss[sru] = 0.0;
                }
            }
            let nodes = &mut ctx.nodes;
            nodes.losses[ui] = loss_rate;
            nodes.inflow[ui] = nodes.lat_flow[ui];
            nodes.outflow[ui] = nodes.losses[ui];

            let full_vol = if ctx.options.routing_model == RoutingModel::Fv { nodes.full_volume[ui] } else { nodes.rpt_full_volume[ui] };
            if nodes.volume[ui] > full_vol && dt > 0.0 {
                nodes.overflow[ui] = (nodes.volume[ui] - full_vol) / dt;
            } else {
                nodes.overflow[ui] = 0.0;
            }
        }
    }

    /// computeConduitLosses: once-per-step conduit evaporation and seepage (KINWAVE / STEADY).
    pub fn compute_conduit_losses(&mut self, ctx: &mut SimulationContext, dt: f64, evap_rate: f64) {
        let n = ctx.n_links();
        for j in 0..n {
            let uj = j as usize;
            if ctx.links.type_[uj] != LinkType::Conduit {
                continue;
            }
            let ucr = ctx.link_subtypes.conduit_row(j) as usize;

            let depth = 0.5 * (ctx.links.old_depth[uj] + ctx.links.depth[uj]);
            let mut evap_loss = 0.0;
            let mut seep_loss = 0.0;

            if depth > constants::FUDGE {
                let cd = &ctx.link_subtypes.conduits;
                let links = &ctx.links;
                let mut length = cd.true_length[ucr];
                if length <= 0.0 {
                    length = cd.length[ucr];
                }
                if length <= 0.0 {
                    length = cd.mod_length[ucr];
                }
                let batch_shape = links.xsect_batch_shape[uj];

                let gw_two_way = !cd.gw_coupled.is_empty() && cd.gw_coupled[ucr] != 0;
                let k_seep = if gw_two_way { cd.gw_kc[ucr] } else { cd.seep_rate[ucr] };
                let want_evap = xsect::is_open(batch_shape) && evap_rate > 0.0;
                let want_seep = k_seep > 0.0;
                if want_evap || want_seep {
                    let mut xs = XSectParams::default();
                    xs.type_ = batch_shape;
                    xs.y_full = links.xsect_y_full[uj];
                    xs.a_full = links.xsect_a_full[uj];
                    xs.w_max = links.xsect_w_max[uj];
                    xs.yw_max = links.xsect_yw_max[uj];
                    xs.r_full = links.xsect_r_full[uj];
                    xs.s_full = links.xsect_s_full[uj];
                    xs.s_max = links.xsect_s_max[uj];
                    xs.y_bot = links.xsect_y_bot[uj];
                    xs.a_bot = links.xsect_a_bot[uj];
                    xs.s_bot = links.xsect_s_bot[uj];
                    xs.r_bot = links.xsect_r_bot[uj];
                    if batch_shape == XSectShape::Irregular as i32 || batch_shape == XSectShape::Custom as i32 || batch_shape == XSectShape::StreetXsect as i32 {
                        let ci = links.xsect_curve[uj];
                        if ci >= 0 && (ci as usize) < ctx.transect_tables.len() {
                            let td = &ctx.transect_tables[ci as usize];
                            xs.transect = ci;
                            xs.area_tbl = Some(&td.area_tbl[..]);
                            xs.hrad_tbl = Some(&td.hrad_tbl[..]);
                            xs.width_tbl = Some(&td.width_tbl[..]);
                            xs.area_lut = Some(&td.area_lut);
                            xs.transect_tbl_size = N_TRANSECT_TBL;
                        }
                    }

                    if want_evap {
                        let top_width = if batch_shape == XSectShape::Irregular as i32 || batch_shape == XSectShape::Custom as i32 {
                            if xs.y_full > 0.0 {
                                xs.w_max * cmin(depth / xs.y_full, 1.0)
                            } else {
                                0.0
                            }
                        } else {
                            xsect::get_w_of_y(&xs, depth)
                        };
                        evap_loss = top_width * length * evap_rate;
                    }

                    if want_seep {
                        let mut d_seep = depth;
                        let width = if batch_shape == XSectShape::RectClosed as i32 {
                            xs.w_max
                        } else {
                            if d_seep >= xs.yw_max {
                                d_seep = xs.yw_max;
                            }
                            xsect::get_w_of_y(&xs, d_seep)
                        };
                        seep_loss = k_seep * width * length;
                        seep_loss *= ctx.climate_state.infil_factor;
                        if gw_two_way {
                            seep_loss *= aquifer_seep_factor(depth, cd.gw_head_rel[ucr], cd.gw_dc[ucr]);
                            seep_loss = cap_aquifer_gain(seep_loss, cd.gw_gain_max[ucr]);
                        }
                    }
                }

                let total = evap_loss + cmax(seep_loss, 0.0);
                if total > 0.0 && (self.model_ == RouteModel::Dynwave || self.model_ == RouteModel::Fv) {
                    let q_avail = links.volume[uj] / (dt * cd.barrels[ucr] as f64);
                    if total > q_avail && q_avail >= 0.0 {
                        let ratio = q_avail / total;
                        evap_loss *= ratio;
                        if seep_loss > 0.0 {
                            seep_loss *= ratio;
                        }
                    }
                }
            }

            apply_seepage_forcing(ctx, uj, &mut seep_loss);

            ctx.link_subtypes.conduits.evap_loss_rate[ucr] = evap_loss;
            ctx.link_subtypes.conduits.seep_loss_rate[ucr] = seep_loss;
        }
    }

    /// executeSteadyFlow (Gap #33): legacy steadyflow_execute.
    pub fn execute_steady_flow(&mut self, ctx: &mut SimulationContext, dt: f64, structures: &mut StructureSolver) -> Result<i32, StepError> {
        let order: Vec<i32> = if self.steady_sorted_links_.is_empty() { (0..ctx.n_links()).collect() } else { self.steady_sorted_links_.clone() };

        self.steady_storage_updated_.clear();
        self.steady_storage_updated_.resize(ctx.n_nodes() as usize, 0);
        self.steady_y_.clear();
        self.steady_y_.resize(ctx.n_links() as usize, 0.0);

        for idx in 0..order.len() as i32 {
            let j = order[idx as usize];
            let uj = j as usize;

            {
                let n1s = ctx.links.node1[uj];
                if n1s >= 0 && ctx.nodes.type_[n1s as usize] == NodeType::Storage && self.steady_storage_updated_[n1s as usize] == 0 {
                    crate::hydraulics::kinematic_wave_step::update_storage_state(ctx, structures, &order, idx, n1s, dt)?;
                    self.steady_storage_updated_[n1s as usize] = 1;
                }
            }

            if ctx.links.type_[uj] != LinkType::Conduit {
                let n1 = ctx.links.node1[uj];
                let n2 = ctx.links.node2[uj];
                let q = crate::hydraulics::kinematic_wave_step::get_link_inflow(ctx, structures, j, dt)?;
                ctx.links.flow[uj] = q;
                if n1 >= 0 {
                    ctx.nodes.outflow[n1 as usize] += q;
                }
                if n2 >= 0 {
                    ctx.nodes.inflow[n2 as usize] += q;
                }
                continue;
            }
            let ucr = ctx.link_subtypes.conduit_row(j) as usize;

            if ctx.links.xsect_shape[uj] == crate::data::link::XsectShape::Dummy {
                let n1 = ctx.links.node1[uj];
                let n2 = ctx.links.node2[uj];
                let q = crate::hydraulics::kinematic_wave_step::get_link_inflow(ctx, structures, j, dt)?;
                ctx.links.flow[uj] = q;
                if n1 >= 0 {
                    ctx.nodes.outflow[n1 as usize] += q;
                }
                if n2 >= 0 {
                    ctx.nodes.inflow[n2 as usize] += q;
                }
                ctx.links.depth[uj] = 0.0;
                ctx.links.volume[uj] = 0.0;
                continue;
            }

            let n1 = ctx.links.node1[uj];
            let mut qin = 0.0;
            if n1 >= 0 {
                qin = crate::hydraulics::kinematic_wave_step::get_link_inflow(ctx, structures, j, dt)?;
            }

            let barrels = cmax(ctx.link_subtypes.conduits.barrels[ucr], 1) as f64;
            let mut q = qin / barrels;

            let loss_rate = crate::hydraulics::kinematic_wave_step::cap_conduit_loss(ctx, ucr, q.abs());
            q -= loss_rate;
            if q < 0.0 {
                q = 0.0;
            }

            let q_full = ctx.link_subtypes.conduits.q_full[ucr];
            let a_full = ctx.links.xsect_a_full[uj];
            let a;
            {
                let xs = steady_xsp(ctx, uj);
                if q >= q_full {
                    q = q_full;
                    a = a_full;
                    qin = q * barrels;
                } else {
                    let beta = ctx.link_subtypes.conduits.beta[ucr];
                    if beta > 0.0 {
                        let s = q / beta;
                        a = xsect::get_a_of_s(&xs, s);
                    } else {
                        a = 0.0;
                    }
                }
            }

            let qout = q * barrels;
            ctx.links.flow[uj] = qout;

            if n1 >= 0 {
                ctx.nodes.outflow[n1 as usize] += qin;
            }
            let n2 = ctx.links.node2[uj];
            if n2 >= 0 {
                ctx.nodes.inflow[n2 as usize] += qout;
            }

            let y = xsect::get_y_of_a(&steady_xsp(ctx, uj), a);
            let cd = &ctx.link_subtypes.conduits;
            let mut length = cd.true_length[ucr];
            if length <= 0.0 {
                length = cd.length[ucr];
            }

            ctx.links.depth[uj] = y;
            ctx.links.volume[uj] = a * length * barrels;
            self.steady_y_[uj] = y;

            ctx.link_subtypes.conduits.full_state[ucr] = if a_full > 0.0 && a >= a_full { 3 } else { 0 };
        }

        let (su, sy) = (self.steady_storage_updated_.clone(), self.steady_y_.clone());
        crate::hydraulics::kinematic_wave_step::finish_routing(ctx, structures, &order, &su, &sy, &sy, dt)?;
        Ok(1)
    }

    /// updateLinkStates: node heads from depths.
    pub fn update_link_states(&mut self, ctx: &mut SimulationContext) {
        let n = ctx.n_nodes() as usize;
        let nodes = &mut ctx.nodes;
        for i in 0..n {
            nodes.head[i] = nodes.invert_elev[i] + nodes.depth[i];
        }
    }
}

/// Routing.cpp's file-local buildXSP: DynamicWave.cpp's view without yw_max.
pub(crate) fn steady_xsp(ctx: &SimulationContext, uk: usize) -> XSectParams<'_> {
    let mut xs = build_xsp(ctx, uk);
    xs.yw_max = 0.0;
    xs
}
