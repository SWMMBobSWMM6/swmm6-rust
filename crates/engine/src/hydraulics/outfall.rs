// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Outfall.hpp and Outfall.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Outfall boundary conditions. So far only the init-time part is translated: the
// outfall-to-link cache (buildOutfallLinkMap) and its helper legacyOffset.
//
// NOT YET TRANSLATED (run-time step functions, with the routing step): setOutfallDepth,
// setAllOutfallDepths, the tidal/timeseries/normal/critical stage evaluators, the outfall
// 2D feedback and route-to helpers, and everything else Outfall.cpp defines.

use crate::core::context::SimulationContext;
use crate::data::link::LinkType;
use crate::data::node::NodeType;

/// The legacy regulator/conduit offset used at one outfall end (Outfall.cpp legacyOffset).
fn legacy_offset(ctx: &SimulationContext, uk: usize, downstream: bool) -> f64 {
    // A regulator's downstream offset is legacy's un-raised crest (kept in offset2).
    if downstream && ctx.links.type_[uk] != LinkType::Conduit && ctx.links.type_[uk] != LinkType::Pump {
        return ctx.links.offset2[uk];
    }
    match ctx.links.type_[uk] {
        LinkType::Conduit => {
            if downstream {
                ctx.links.offset2[uk]
            } else {
                ctx.links.offset1[uk]
            }
        }
        LinkType::Orifice => ctx.links.offset1[uk],
        LinkType::Weir => {
            let wr = ctx.link_subtypes.weir_row(uk as i32);
            if wr >= 0 {
                ctx.link_subtypes.weirs.crest_height[wr as usize]
            } else {
                0.0
            }
        }
        LinkType::Outlet => {
            let olr = ctx.link_subtypes.outlet_row(uk as i32);
            if olr >= 0 {
                ctx.link_subtypes.outlets.crest_height[olr as usize]
            } else {
                0.0
            }
        }
        _ => 0.0,
    }
}

/// Cache, for every outfall row, the link that sets its depth and that link's offset at
/// the outfall end (legacy link_setOutfallDepth: every link in link order sets ONE end,
/// the downstream outfall if any, else the upstream one; the last such link stands).
pub fn build_outfall_link_map(ctx: &mut SimulationContext) {
    let n_links = ctx.n_links();
    let n_out = ctx.node_subtypes.outfalls.count() as usize;

    // Reset the cached connectivity on every outfall row.
    ctx.node_subtypes.outfalls.link_idx.clear();
    ctx.node_subtypes.outfalls.link_idx.resize(n_out, -1);
    ctx.node_subtypes.outfalls.link_offset.clear();
    ctx.node_subtypes.outfalls.link_offset.resize(n_out, 0.0);

    for k in 0..n_links {
        let uk = k as usize;
        let n1 = ctx.links.node1[uk];
        let n2 = ctx.links.node2[uk];
        let n;
        let z;
        if n2 >= 0 && ctx.nodes.type_[n2 as usize] == NodeType::Outfall {
            n = n2;
            z = legacy_offset(ctx, uk, true);
        } else if n1 >= 0 && ctx.nodes.type_[n1 as usize] == NodeType::Outfall {
            n = n1;
            z = legacy_offset(ctx, uk, false);
        } else {
            continue;
        }
        let r = ctx.node_subtypes.outfall_row(n);
        if r >= 0 {
            ctx.node_subtypes.outfalls.link_idx[r as usize] = k;
            ctx.node_subtypes.outfalls.link_offset[r as usize] = z;
        }
    }
}
