// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Exfiltration.cpp (HydroCouple/
// openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed
// under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Storage exfiltration at run time (exfiltration.rs holds the state and init):
// ExfilSolver::computeAll, the raw Green-Ampt bottom and bank seepage rate handed to
// Router::initNodeFlows through storages.exfil_rate.

use crate::core::cmath::cmin;
use crate::core::context::SimulationContext;
use crate::core::unit_conversion as ucf;
use crate::data::lid_node::{GreenAmptState, LidNodeBottom};
use crate::hydraulics::exfiltration::ExfilSolver;
use crate::hydraulics::node;
use crate::hydrology::green_ampt_step::grnampt_get_infil;
use crate::hydrology::infiltration::InfilModel;
use crate::hydrology::lid_node_step as lidnode;

/// The exfiltration rate of one Green-Ampt object (exfil_getLoss's ga_rate lambda).
fn ga_rate(ga: &mut GreenAmptState, d: f64, dt: f64, infil_factor: f64, recovery_factor: f64) -> f64 {
    if ga.IMDmax == 0.0 {
        return ga.Ks * infil_factor;
    }
    grnampt_get_infil(ga, 0.0, d, dt, InfilModel::ModGreenAmpt, infil_factor, recovery_factor)
}

impl ExfilSolver {
    /// ExfilSolver::computeAll: each storage node's raw exfiltration rate (cfs).
    pub fn compute_all(&mut self, ctx: &mut SimulationContext, dt: f64) {
        for k in 0..self.soa_.count as usize {
            let ni = self.soa_.node_idx[k];
            if ni < 0 {
                continue;
            }
            let uni = ni as usize;
            if lidnode::richards_mode(ctx, ni) || lidnode::bottom_mode(ctx, ni) != LidNodeBottom::Current {
                let sr = ctx.node_subtypes.storage_row(ni) as usize;
                ctx.node_subtypes.storages.exfil_rate[sr] = 0.0;
                continue;
            }
            if uni < ctx.nodes.aquifer2d_bed.len() && ctx.nodes.aquifer2d_bed[uni] != 0 {
                let sr0 = ctx.node_subtypes.storage_row(ni);
                if sr0 >= 0 {
                    ctx.node_subtypes.storages.exfil_rate[sr0 as usize] = 0.0;
                }
                continue;
            }

            let depth = ctx.nodes.depth[uni];
            let mut total_loss = 0.0;
            let infil_factor = ctx.climate_state.infil_factor;
            let recovery_factor = ctx.climate_state.recovery_factor;

            let clog = lidnode::bottom_clogging_factor(ctx, ctx.node_subtypes.storage_row(ni));
            total_loss += ga_rate(&mut self.soa_.btm_ga[k], depth, dt, infil_factor, recovery_factor) * self.soa_.btm_area[k] * clog;

            if depth > self.soa_.bank_min_depth[k] {
                let us = ucf::get_unit_system(ctx.options.flow_units as i32);
                let mut area = node::get_surf_area(&ctx.nodes, ni, depth, Some(&ctx.tables), us, Some(&ctx.node_subtypes));
                area = cmin(area, self.soa_.bank_max_area[k]) - self.soa_.btm_area[k];
                if area > 0.0 {
                    let mut bank_depth = depth;
                    if self.soa_.bank_ga[k].IMDmax != 0.0 {
                        if depth > self.soa_.bank_max_depth[k] {
                            bank_depth = depth - self.soa_.bank_max_depth[k] + (self.soa_.bank_max_depth[k] - self.soa_.bank_min_depth[k]) / 2.0;
                        } else {
                            bank_depth = (depth - self.soa_.bank_min_depth[k]) / 2.0;
                        }
                    }
                    total_loss += area * ga_rate(&mut self.soa_.bank_ga[k], bank_depth, dt, infil_factor, recovery_factor);
                }
            }

            let sr = ctx.node_subtypes.storage_row(uni as i32);
            if sr >= 0 {
                ctx.node_subtypes.storages.exfil_rate[sr as usize] = total_loss;
            }
        }
    }
}
