// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Divider.cpp and Divider.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Flow divider outflow (legacy node.c divider_getOutflow).

use crate::core::context::SimulationContext;
use crate::core::unit_conversion as ucf;
use crate::data::node::DividerType;
use crate::data::table::table_lookup_cursor;

/// Legacy FLOW_TOL (consts.h): minimum significant flow, cfs.
const FLOW_TOL: f64 = 0.00001;

/// divider::getOutflow: the flow a divider sends into one of its outflow links.
pub fn get_outflow(ctx: &mut SimulationContext, node_idx: i32, link_idx: i32) -> f64 {
    let r = ctx.node_subtypes.divider_row(node_idx);
    if r < 0 {
        return 0.0;
    }
    let ur = r as usize;
    let un = node_idx as usize;
    let ucf_flow = ucf::QCF[ctx.options.flow_units as usize];

    let q_in = ctx.nodes.inflow[un] + ctx.nodes.overflow[un];
    let div_link = ctx.node_subtypes.dividers.link[ur];
    let mut q_out;

    match ctx.node_subtypes.dividers.method[ur] {
        DividerType::Cutoff => {
            let cutoff = ctx.node_subtypes.dividers.cutoff[ur];
            q_out = if q_in <= cutoff { 0.0 } else { q_in - cutoff };
        }
        DividerType::OverflowDiv => {
            q_out = if link_idx != div_link { q_in } else { q_in - ctx.nodes.outflow[un] };
            if q_out < FLOW_TOL {
                q_out = 0.0;
            }
            return q_out;
        }
        DividerType::Weir => {
            let divs = &ctx.node_subtypes.dividers;
            let q_min = divs.cutoff[ur];
            let dh_max = divs.max_depth[ur];
            let c_weir = divs.cd[ur];
            let q_max = c_weir * dh_max.powf(1.5) / ucf_flow; // LIBM pow
            if q_in <= q_min {
                q_out = 0.0;
            } else {
                let f = (q_in - q_min) / (q_max - q_min);
                if f > 1.0 {
                    q_out = q_max * f.sqrt();
                } else {
                    q_out = c_weir * (f * dh_max).powf(1.5) / ucf_flow; // LIBM pow
                }
            }
        }
        DividerType::Tabular => {
            let m = ctx.node_subtypes.dividers.curve[ur];
            if m >= 0 && (m as usize) < ctx.tables.tables.len() {
                q_out = table_lookup_cursor(&mut ctx.tables.tables[m as usize], q_in * ucf_flow) / ucf_flow;
            } else {
                q_out = 0.0;
            }
        }
    }

    if q_out > q_in {
        q_out = q_in;
    }
    if link_idx != div_link {
        q_out = q_in - q_out;
    }
    q_out
}
