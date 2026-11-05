// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/SWMMEngine.cpp (assembleLateralInflows,
// effectiveUserLatFlow) (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// SWMMEngine::assembleLateralInflows: the decomposed inflow buffers summed into
// nodes.lat_flow in legacy's source order, the step's inflow mass-balance terms, and the
// street-inlet transfer (InletSolver::computeAll) that follows.
//
// MERGE NOTE: the engine-lifecycle work (branch rs6-init) translates this helper too. This
// copy exists so the routing-step cross test can run stepRouting whole; keep the one whose
// gate is identical (both follow the C++ line for line) and delete the other.
//
// NOT YET TRANSLATED: InletSolver::computeAll for a deck with street inlets
// (StepError::NotTranslated), and the OPENSWMM_HAS_2D source-delivery receipts (off in the
// reference build).

use crate::core::cmath::cmax;
use crate::core::context::SimulationContext;
use crate::core::engine::SwmmEngine;
use crate::core::engine_step_routing::StepError;
use crate::data::forcing::ForcingMode;

/// effectiveUserLatFlow: the runtime-API lateral flow with any forcing overlaid.
pub(crate) fn effective_user_lat_flow(ctx: &SimulationContext, uj: usize) -> f64 {
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

/// The 2D->1D queue drain rule (uniform rate over the remaining span, flushed when
/// remaining <= dt).
fn drain_queue(qm: &mut f64, remaining: f64, dt_routing: f64) -> f64 {
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
}

impl SwmmEngine {
    /// SWMMEngine::assembleLateralInflows
    pub(crate) fn assemble_lateral_inflows(&mut self, dt_routing: f64) -> Result<(), StepError> {
        let ctx = &mut self.ctx_;
        let mut sum_dw = 0.0;
        let mut sum_gw = 0.0;
        let mut sum_rdii = 0.0;
        let mut sum_ext = 0.0;
        let mut sum_wd = 0.0;

        for j in 0..ctx.n_nodes() {
            let uj = j as usize;
            {
                let mut q_couple = 0.0;
                if dt_routing > 0.0 && ctx.nodes.coupling_queue[uj] != 0.0 {
                    if ctx.coupling_delivery_remaining > dt_routing {
                        q_couple = ctx.nodes.coupling_queue[uj] / ctx.coupling_delivery_remaining;
                        ctx.nodes.coupling_queue[uj] -= q_couple * dt_routing;
                    } else {
                        q_couple = ctx.nodes.coupling_queue[uj] / dt_routing;
                        ctx.nodes.coupling_queue[uj] = 0.0;
                    }
                }
                if dt_routing > 0.0 && ctx.nodes.coupling_volume[uj] != 0.0 {
                    q_couple += ctx.nodes.coupling_volume[uj] / dt_routing;
                    ctx.nodes.coupling_volume[uj] = 0.0;
                }
                ctx.nodes.coupling_inflow[uj] = q_couple;

                let remaining = ctx.coupling_delivery_remaining;
                if !ctx.nodes.coupling_qual_queue.is_empty() {
                    let np = ctx.n_pollutants() as usize;
                    for p in 0..np {
                        let idx = uj * np + p;
                        let mut qm = ctx.nodes.coupling_qual_queue[idx];
                        let r = drain_queue(&mut qm, remaining, dt_routing);
                        ctx.nodes.coupling_qual_queue[idx] = qm;
                        ctx.nodes.coupling_qual_inflow[idx] = r;
                    }
                }
                if uj < ctx.nodes.coupling_age_vol_queue.len() {
                    let mut qm = ctx.nodes.coupling_age_vol_queue[uj];
                    let r = drain_queue(&mut qm, remaining, dt_routing);
                    ctx.nodes.coupling_age_vol_queue[uj] = qm;
                    ctx.nodes.coupling_age_vol_inflow[uj] = r;
                }
                if uj < ctx.nodes.coupling_temp_vol_queue.len() {
                    let mut qm = ctx.nodes.coupling_temp_vol_queue[uj];
                    let r = drain_queue(&mut qm, remaining, dt_routing);
                    ctx.nodes.coupling_temp_vol_queue[uj] = qm;
                    ctx.nodes.coupling_temp_vol_inflow[uj] = r;
                }
            }

            ctx.nodes.lat_flow[uj] = ctx.nodes.ext_inflow[uj] + ctx.nodes.dwf_inflow[uj];

            sum_dw += ctx.nodes.dwf_inflow[uj];
            sum_gw += ctx.nodes.gw_inflow[uj];
            sum_rdii += ctx.nodes.rdii_inflow[uj];
            {
                let q_ext = ctx.nodes.ext_inflow[uj];
                if q_ext >= 0.0 {
                    sum_ext += q_ext;
                } else {
                    sum_wd += -q_ext;
                }
            }
            sum_ext += ctx.nodes.iface_inflow[uj];

            if ctx.nodes.coupling_inflow[uj] > 0.0 {
                sum_ext += ctx.nodes.coupling_inflow[uj];
            }
        }

        if ctx.coupling_delivery_remaining > 0.0 {
            ctx.coupling_delivery_remaining = cmax(0.0, ctx.coupling_delivery_remaining - dt_routing);
        }

        {
            let n_sub = ctx.n_subcatches();
            let n_nodes = ctx.n_nodes();
            if self.wet_q_interp_.len() as i32 == n_sub {
                for i in 0..n_sub as usize {
                    let out_node = ctx.subcatches.outlet_node[i];
                    if out_node >= 0 && out_node < n_nodes {
                        ctx.nodes.lat_flow[out_node as usize] += self.wet_q_interp_[i];
                    }
                }
            }
            if self.gw_q_interp_.len() as i32 == n_sub {
                for i in 0..n_sub as usize {
                    let gwn = self.gw_q_node_[i];
                    if gwn >= 0 && gwn < n_nodes {
                        ctx.nodes.lat_flow[gwn as usize] += self.gw_q_interp_[i];
                    }
                }
            }
        }

        for dq in self.lid_drain_q_units_.iter() {
            ctx.nodes.lat_flow[dq.0 as usize] += dq.1;
        }
        for j in 0..ctx.n_nodes() as usize {
            ctx.nodes.lat_flow[j] += ctx.nodes.rdii_inflow[j];
            ctx.nodes.lat_flow[j] += ctx.nodes.iface_inflow[j];
            let q_user = effective_user_lat_flow(ctx, j);
            if q_user > 0.0 {
                sum_ext += q_user;
            }
            ctx.nodes.lat_flow[j] += q_user + ctx.nodes.coupling_inflow[j];
        }

        ctx.mass_balance.step_dw_inflow = sum_dw;
        ctx.mass_balance.step_gw_inflow = sum_gw;
        ctx.mass_balance.step_rdii_inflow = sum_rdii;
        ctx.mass_balance.step_ext_inflow = sum_ext;
        ctx.mass_balance.step_ext_withdrawal = sum_wd;

        // Street inlet capture / backflow (InletSolver::computeAll)
        if self.inlet_.soa_.count != 0 {
            return Err(StepError::NotTranslated("street inlets"));
        }
        Ok(())
    }
}
