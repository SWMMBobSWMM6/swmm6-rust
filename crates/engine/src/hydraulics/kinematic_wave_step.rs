// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/KinematicWave.cpp and
// KinematicWave.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The kinematic-wave routing step and the helpers the steady-flow router shares with it
// (kinematic_wave.rs holds the KW solver's state and init): KWSolver::execute and
// solveConduit, and kinwave::getLinkInflow, updateStorageState, finishRouting,
// capConduitLoss with the file-local storageConduitOutflow and buildXSP_KW.
//
// Differences in form, never in result:
//   * `hydstruct::StructureSolver* structures` (never null in the engine) is &mut.
//   * The fallback 0..n_links order (a function-local thread_local vector upstream) is built
//     when the sorted order is empty.

use crate::core::cmath::cmax;
use crate::core::constants;
use crate::core::context::SimulationContext;
use crate::core::engine_step_routing::StepError;
use crate::core::unit_conversion as ucf;
use crate::data::link::{LinkType, XsectShape};
use crate::data::node::NodeType;
use crate::data::table::TableType;
use crate::hydraulics::divider;
use crate::hydraulics::hyd_structures::StructureSolver;
use crate::hydraulics::kinematic_wave::{KWSolver, EPSIL, WT, WX};
use crate::hydraulics::node;
use crate::hydraulics::node_step;
use crate::hydraulics::routing_step::steady_xsp;
use crate::hydraulics::xsect_batch::{XSectParams, XSectShape};
use crate::hydraulics::xsection as xsect;
use crate::math::find_root;

const TINY: f64 = 1.0e-6;
const STOR_MAXITER: i32 = 10;
const STOR_OMEGA: f64 = 0.55;
const STOR_STOPTOL: f64 = 0.005;

/// buildXSP_KW: Routing.cpp's view (no yw_max).
fn build_xsp_kw(ctx: &SimulationContext, uk: usize) -> XSectParams<'_> {
    steady_xsp(ctx, uk)
}

/// storageConduitOutflow: normal-depth flow of a conduit draining a storage unit.
fn storage_conduit_outflow(ctx: &SimulationContext, i: i32, j: i32) -> f64 {
    let uj = j as usize;
    let ui = i as usize;
    if ctx.links.type_[uj] != LinkType::Conduit {
        return 0.0;
    }
    let y = ctx.nodes.depth[ui] - ctx.links.offset1[uj];
    if y <= 0.0 {
        return 0.0;
    }
    let cd = &ctx.link_subtypes.conduits;
    let ucr = ctx.link_subtypes.conduit_row(j) as usize;
    if y >= ctx.links.xsect_y_full[uj] {
        return cd.q_full[ucr];
    }
    let xs = build_xsp_kw(ctx, uj);
    let a = xsect::get_a_of_y(&xs, y);
    cd.beta[ucr] * xsect::get_s_of_a(&xs, a)
}

/// kinwave::getLinkInflow (legacy flowrout.c getLinkInflow).
pub fn get_link_inflow(ctx: &mut SimulationContext, structures: &mut StructureSolver, j: i32, dt: f64) -> Result<f64, StepError> {
    let uj = j as usize;
    let n1 = ctx.links.node1[uj];
    if n1 < 0 {
        return Ok(0.0);
    }
    let un1 = n1 as usize;
    let lt = ctx.links.type_[uj];
    let mut q = 0.0;

    if ctx.links.setting[uj] == 0.0 {
        return Ok(0.0);
    }

    if lt == LinkType::Conduit {
        if ctx.nodes.type_[un1] == NodeType::Storage {
            q = storage_conduit_outflow(ctx, n1, j);
        } else if ctx.nodes.type_[un1] == NodeType::Divider {
            q = divider::get_outflow(ctx, n1, j);
        } else {
            q = ctx.nodes.inflow[un1] + ctx.nodes.overflow[un1];
        }
        let qlim = ctx.links.q_limit[uj];
        if qlim > 0.0 && q > qlim {
            q = qlim;
        }
    } else if lt == LinkType::Pump || ctx.nodes.type_[un1] == NodeType::Storage {
        let q_saved = ctx.links.flow[uj];
        structures.compute_non_conduit_flow_one(ctx, dt, None, j);
        q = ctx.links.flow[uj];
        ctx.links.flow[uj] = q_saved;
    }

    Ok(node_step::get_max_outflow(&ctx.nodes, n1, q, dt))
}

/// kinwave::updateStorageState: converge a storage unit's depth (legacy flowrout.c).
pub fn update_storage_state(ctx: &mut SimulationContext, structures: &mut StructureSolver, order: &[i32], pos: i32, i: i32, dt: f64) -> Result<(), StepError> {
    let ui = i as usize;
    let us = ucf::get_unit_system(ctx.options.flow_units as i32);

    let v_fixed = ctx.nodes.old_volume[ui] + 0.5 * (ctx.nodes.old_net_inflow[ui] + ctx.nodes.inflow[ui] - ctx.nodes.outflow[ui]) * dt;
    let full_vol = node::get_volume(&ctx.nodes, i, ctx.nodes.full_depth[ui], Some(&ctx.tables), us, Some(&ctx.node_subtypes));

    let mut d1 = ctx.nodes.depth[ui];
    for _iter in 1..STOR_MAXITER {
        // storageOutflow
        let mut outflow = 0.0;
        let mut k = pos as usize;
        while k < order.len() {
            let m = order[k];
            if ctx.links.node1[m as usize] != i {
                break;
            }
            outflow += get_link_inflow(ctx, structures, m, dt)?;
            k += 1;
        }
        let mut v2 = v_fixed - 0.5 * outflow * dt;
        v2 = cmax(v2, 0.0);

        ctx.nodes.overflow[ui] = 0.0;
        if v2 > full_vol {
            ctx.nodes.overflow[ui] = (v2 - cmax(ctx.nodes.old_volume[ui], full_vol)) / dt;
            if ctx.nodes.overflow[ui] < constants::FUDGE {
                ctx.nodes.overflow[ui] = 0.0;
            }
            if !ctx.options.allow_ponding || ctx.nodes.ponded_area[ui] == 0.0 {
                v2 = full_vol;
            }
        }

        ctx.nodes.volume[ui] = v2;
        let mut d2 = node_step::get_depth(&ctx.nodes, i, v2, Some(&ctx.tables), us, Some(&ctx.node_subtypes));
        ctx.nodes.depth[ui] = d2;

        d2 = (1.0 - STOR_OMEGA) * d1 + STOR_OMEGA * d2;
        let stop = (d2 - d1).abs() <= STOR_STOPTOL;
        d1 = d2;
        if stop {
            break;
        }
    }
    ctx.nodes.head[ui] = ctx.nodes.invert_elev[ui] + ctx.nodes.depth[ui];
    Ok(())
}

impl KWSolver {
    /// solveConduit: the per-conduit kinematic-wave continuity solve (legacy kinwave.c).
    #[allow(clippy::too_many_arguments)]
    fn solve_conduit(&mut self, idx: i32, xs: &XSectParams, q_full: f64, a_full: f64, _s_full: f64, beta: f64, length: f64, dt: f64, loss_rate: f64) -> i32 {
        let ui = idx as usize;
        if q_full <= 0.0 || a_full <= 0.0 {
            self.q_out_[ui] = 0.0;
            self.a_out_[ui] = 0.0;
            return 0;
        }

        let mut q_in_norm = self.q_in_[ui] / q_full;
        let q3 = loss_rate / q_full;
        let beta1 = beta / q_full;

        let prev_q1 = self.q1_[ui] / q_full;
        let prev_a1 = self.a1_[ui] / a_full;
        let prev_q2 = self.q2_[ui] / q_full;
        let prev_a2 = self.a2_[ui] / a_full;

        let a_in_norm = if q_in_norm >= 1.0 {
            1.0
        } else if q_in_norm <= 0.0 {
            0.0
        } else {
            let s_needed = q_in_norm / beta1;
            xsect::get_a_of_s(xs, s_needed) / a_full
        };

        if q_in_norm <= TINY && prev_q2 <= TINY {
            self.a_in_[ui] = a_in_norm * a_full;
            self.a_out_[ui] = 0.0;
            self.q_out_[ui] = 0.0;
            if q_in_norm > 1.0 {
                q_in_norm = 1.0;
            }
            self.q_in_[ui] = q_in_norm * q_full;
            self.q1_[ui] = self.q_in_[ui];
            self.a1_[ui] = self.a_in_[ui];
            self.q2_[ui] = 0.0;
            self.a2_[ui] = 0.0;
            return 1;
        }

        let dxdt = length / dt * a_full / q_full;
        let dq = prev_q2 - prev_q1;

        let c1 = dxdt * WT / WX;
        let mut c2 = (1.0 - WT) * (a_in_norm - prev_a1);
        c2 -= WT * prev_a2;
        c2 = c2 * dxdt / WX;
        c2 = c2 + (1.0 - WX) / WX * dq - q_in_norm;
        c2 += q3 / WX;

        let mut a_hi = 1.0;
        let mut f_hi = 1.0 + c1 + c2;
        let mut a_lo = if xs.type_ == XSectShape::Irregular as i32 || xs.type_ == XSectShape::Custom as i32 { xs.a_bot / a_full } else { xsect::get_amax(xs) };
        let mut f_lo = if a_lo < a_hi { (beta1 * xs.s_max) + (c1 * a_lo) + c2 } else { f_hi };

        if f_hi * f_lo > 0.0 {
            a_hi = a_lo;
            f_hi = f_lo;
            a_lo = 0.0;
            f_lo = c2;
        }

        let mut a = prev_a2;
        let mut result;
        if f_hi * f_lo <= 0.0 {
            if a < a_lo || a > a_hi {
                a = 0.5 * (a_lo + a_hi);
            }
            if f_lo > f_hi {
                std::mem::swap(&mut a_lo, &mut a_hi);
            }
            result = find_root::newton(a_lo, a_hi, &mut a, EPSIL, |an: f64, f: &mut f64, df: &mut f64| {
                *f = (beta1 * xsect::get_s_of_a(xs, an * a_full)) + (c1 * an) + c2;
                *df = (beta1 * a_full * xsect::get_ds_da(xs, an * a_full)) + c1;
            });
            if result <= 0 {
                result = -1;
            }
        } else if f_lo < 0.0 {
            a = if q_in_norm > 1.0 { a_in_norm } else { 1.0 };
            result = -2;
        } else if f_lo > 0.0 {
            a = 0.0;
            result = -3;
        } else {
            result = -1;
        }
        if result <= 0 {
            result = 1;
        }

        let q_out_norm = beta1 * xsect::get_s_of_a(xs, a * a_full);

        self.a_in_[ui] = a_in_norm * a_full;
        self.a_out_[ui] = a * a_full;
        self.q_out_[ui] = q_out_norm * q_full;

        if q_in_norm > 1.0 {
            q_in_norm = 1.0;
        }
        self.q_in_[ui] = q_in_norm * q_full;

        self.q1_[ui] = self.q_in_[ui];
        self.a1_[ui] = self.a_in_[ui];
        self.q2_[ui] = self.q_out_[ui];
        self.a2_[ui] = self.a_out_[ui];

        result
    }

    /// KWSolver::execute: route every link in topological order.
    pub fn execute(&mut self, ctx: &mut SimulationContext, dt: f64, structures: &mut StructureSolver) -> Result<i32, StepError> {
        let mut total_iters = 0;
        let mut n_solved = 0;

        self.storage_updated_.clear();
        self.storage_updated_.resize(ctx.n_nodes() as usize, 0);
        self.y1_.iter_mut().for_each(|x| *x = 0.0);
        self.y2_.iter_mut().for_each(|x| *x = 0.0);

        let order: Vec<i32> = if self.sorted_links_.is_empty() { (0..ctx.n_links()).collect() } else { self.sorted_links_.clone() };

        for idx in 0..order.len() as i32 {
            let j = order[idx as usize];
            let uj = j as usize;

            {
                let n1s = ctx.links.node1[uj];
                if n1s >= 0 && ctx.nodes.type_[n1s as usize] == NodeType::Storage && self.storage_updated_[n1s as usize] == 0 {
                    update_storage_state(ctx, structures, &order, idx, n1s, dt)?;
                    self.storage_updated_[n1s as usize] = 1;
                }
            }

            if ctx.links.type_[uj] != LinkType::Conduit {
                let n1 = ctx.links.node1[uj];
                let n2 = ctx.links.node2[uj];
                let q = get_link_inflow(ctx, structures, j, dt)?;
                ctx.links.flow[uj] = q;
                if n1 >= 0 {
                    ctx.nodes.outflow[n1 as usize] += q;
                }
                if n2 >= 0 {
                    ctx.nodes.inflow[n2 as usize] += q;
                }
                continue;
            }

            if ctx.links.xsect_shape[uj] == XsectShape::Dummy {
                let n1 = ctx.links.node1[uj];
                let n2 = ctx.links.node2[uj];
                let q = get_link_inflow(ctx, structures, j, dt)?;
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
            let ucr = ctx.link_subtypes.conduit_row(j) as usize;

            let n1 = ctx.links.node1[uj];
            let mut qin = 0.0;
            if n1 >= 0 {
                qin = get_link_inflow(ctx, structures, j, dt)?;
            }

            let barrels = cmax(ctx.link_subtypes.conduits.barrels[ucr], 1) as f64;
            let qin_per_barrel = qin / barrels;

            let q_full = ctx.link_subtypes.conduits.q_full[ucr];
            let a_full = ctx.links.xsect_a_full[uj];
            let s_full = ctx.links.xsect_s_full[uj];
            let beta = ctx.link_subtypes.conduits.beta[ucr];
            let mut length = ctx.link_subtypes.conduits.true_length[ucr];
            if length <= 0.0 {
                length = ctx.link_subtypes.conduits.length[ucr];
            }

            let loss_rate = cap_conduit_loss(ctx, ucr, qin_per_barrel.abs());

            self.q_in_[uj] = qin_per_barrel;

            let (y_in, y_out) = {
                let xs = build_xsp_kw(ctx, uj);
                let iters = self.solve_conduit(uj as i32, &xs, q_full, a_full, s_full, beta, length, dt, loss_rate);
                total_iters += iters;
                n_solved += 1;
                (xsect::get_y_of_a(&xs, self.a_in_[uj]), xsect::get_y_of_a(&xs, self.a_out_[uj]))
            };

            let qout = self.q_out_[uj] * barrels;
            qin = self.q_in_[uj] * barrels;
            ctx.links.flow[uj] = qout;
            ctx.links.kw_inflow[uj] = qin;

            if n1 >= 0 {
                ctx.nodes.outflow[n1 as usize] += qin;
            }
            let n2 = ctx.links.node2[uj];
            if n2 >= 0 {
                ctx.nodes.inflow[n2 as usize] += qout;
            }

            ctx.links.depth[uj] = 0.5 * (y_in + y_out);
            ctx.links.volume[uj] = 0.5 * (self.a_in_[uj] + self.a_out_[uj]) * length * barrels;
            self.y1_[uj] = y_in;
            self.y2_[uj] = y_out;

            {
                let mut fs: i8 = 0;
                if a_full > 0.0 {
                    if self.a_in_[uj] >= a_full {
                        fs |= 1;
                    }
                    if self.a_out_[uj] >= a_full {
                        fs |= 2;
                    }
                }
                ctx.link_subtypes.conduits.full_state[ucr] = fs;
            }
        }

        let su = self.storage_updated_.clone();
        finish_routing(ctx, structures, &order, &su, &self.y1_, &self.y2_, dt)?;

        Ok(if n_solved > 0 { total_iters / n_solved } else { 1 })
    }
}

/// kinwave::finishRouting: legacy flowrout_execute's end-of-step passes (KW / steady).
pub fn finish_routing(ctx: &mut SimulationContext, structures: &mut StructureSolver, order: &[i32], storage_updated: &[i8], link_y1: &[f64], link_y2: &[f64], dt: f64) -> Result<(), StepError> {
    let n_nodes = ctx.n_nodes();
    let n_links = ctx.n_links();
    let us = ucf::get_unit_system(ctx.options.flow_units as i32);

    let mut legacy_fv = vec![0.0f64; n_nodes as usize];
    {
        let ucf_vol = ucf::UCF[ucf::VOLUME as usize][us as usize];
        for j in 0..n_links {
            let uj = j as usize;
            if ctx.links.type_[uj] != LinkType::Pump {
                continue;
            }
            let pr = ctx.link_subtypes.pump_row(j);
            let ci = if pr >= 0 { ctx.link_subtypes.pumps.curve[pr as usize] } else { -1 };
            if ci < 0 || ci >= ctx.tables.tables.len() as i32 {
                continue;
            }
            let tbl = &ctx.tables.tables[ci as usize];
            if tbl.type_ != TableType::CurvePump1 {
                continue;
            }
            let n1 = ctx.links.node1[uj];
            if n1 < 0 || ctx.nodes.type_[n1 as usize] == NodeType::Storage {
                continue;
            }
            let mut xmax = tbl.x_max;
            if xmax <= 0.0 && !tbl.x.is_empty() {
                // std::max_element: the first of the largest values
                let mut m = tbl.x[0];
                for &v in tbl.x.iter().skip(1) {
                    if m < v {
                        m = v;
                    }
                }
                xmax = m;
            }
            let un1 = n1 as usize;
            legacy_fv[un1] = cmax(legacy_fv[un1], xmax / ucf_vol);
        }
    }

    let mut degree = vec![0i32; n_nodes as usize];
    for j in 0..n_links as usize {
        let mut n = if ctx.links.direction[j] < 0 { ctx.links.node2[j] } else { ctx.links.node1[j] };
        if n < 0 || n >= n_nodes {
            continue;
        }
        if ctx.nodes.type_[n as usize] == NodeType::Outfall {
            n = if ctx.links.direction[j] < 0 { ctx.links.node1[j] } else { ctx.links.node2[j] };
            if n < 0 || n >= n_nodes {
                continue;
            }
        }
        degree[n as usize] += 1;
    }

    for i in 0..n_nodes {
        let ui = i as usize;
        if ctx.nodes.type_[ui] == NodeType::Storage {
            if ui >= storage_updated.len() || storage_updated[ui] == 0 {
                update_storage_state(ctx, structures, order, order.len() as i32, i, dt)?;
            }
            continue;
        }
        let nodes = &mut ctx.nodes;
        let net = nodes.inflow[ui] - nodes.outflow[ui] - nodes.losses[ui];
        let mut v = nodes.old_volume[ui] + net * dt;
        if v < constants::FUDGE {
            v = 0.0;
        }
        nodes.overflow[ui] = 0.0;
        let can_pond = ctx.options.allow_ponding && nodes.ponded_area[ui] > 0.0;
        if v > legacy_fv[ui] {
            nodes.overflow[ui] = (v - cmax(nodes.old_volume[ui], legacy_fv[ui])) / dt;
            if nodes.overflow[ui] < constants::FUDGE {
                nodes.overflow[ui] = 0.0;
            }
            if !can_pond {
                v = legacy_fv[ui];
            }
        }
        nodes.volume[ui] = v;
        nodes.depth[ui] = 0.0;
        nodes.head[ui] = nodes.invert_elev[ui];
    }

    let raise = |ctx: &mut SimulationContext, ni: i32, mut y: f64| {
        if ni < 0 || ni >= n_nodes {
            return;
        }
        let uni = ni as usize;
        let nodes = &mut ctx.nodes;
        let nt = nodes.type_[uni];
        if nt == NodeType::Storage {
            return;
        }
        if nt != NodeType::Outfall && degree[uni] > 0 && nodes.overflow[uni] > 0.0 {
            y = nodes.full_depth[uni];
        }
        if nodes.depth[uni] < y {
            nodes.depth[uni] = if nodes.full_depth[uni] > 0.0 && y > nodes.full_depth[uni] { nodes.full_depth[uni] } else { y };
            nodes.head[uni] = nodes.invert_elev[uni] + nodes.depth[uni];
        }
    };
    for j in 0..n_links as usize {
        if ctx.links.type_[j] != LinkType::Conduit {
            ctx.links.depth[j] = 0.0;
            ctx.links.volume[j] = 0.0;
            continue;
        }
        let y1 = if j < link_y1.len() { link_y1[j] } else { 0.0 };
        let y2 = if j < link_y2.len() { link_y2[j] } else { 0.0 };
        let (n1, o1, n2, o2) = (ctx.links.node1[j], ctx.links.offset1[j], ctx.links.node2[j], ctx.links.offset2[j]);
        raise(ctx, n1, y1 + o1);
        raise(ctx, n2, y2 + o2);
    }
    Ok(())
}

/// kinwave::capConduitLoss: cap the conduit's loss rates on this solve's inflow.
pub fn cap_conduit_loss(ctx: &mut SimulationContext, ucr: usize, q: f64) -> f64 {
    let cd = &mut ctx.link_subtypes.conduits;
    let mut evap = cd.evap_loss_rate[ucr];
    let mut seep = cd.seep_loss_rate[ucr];
    let losing = evap + cmax(seep, 0.0);
    if losing > q {
        let ratio = q / losing;
        evap *= ratio;
        if seep > 0.0 {
            seep *= ratio;
        }
        cd.evap_loss_rate[ucr] = evap;
        cd.seep_loss_rate[ucr] = seep;
    }
    evap + seep
}
