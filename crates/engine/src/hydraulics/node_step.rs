// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Node.cpp (HydroCouple/openswmm.engine,
// branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache
// License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The node functions the routing step calls that hydraulics/node.rs (getVolume, getSurfArea,
// getDepth) does not hold: getPondedArea, getMaxOutflow, getOverflow.

use crate::core::cmath::{cmax, cmin};
use crate::data::node::NodeData;
use crate::data::node_subtypes::NodeSubtypes;
use crate::data::table::TableData;
use crate::hydraulics::node::get_surf_area;

/// node::getDepth is hydraulics::node::get_depth (translated for SWMMEngine::initialize).
pub use crate::hydraulics::node::get_depth;

/// node::getPondedArea
pub fn get_ponded_area(nodes: &NodeData, idx: i32, depth: f64, tables: Option<&TableData>, unit_sys: i32, subs: Option<&NodeSubtypes>) -> f64 {
    let ui = idx as usize;
    if depth <= nodes.full_depth[ui] || nodes.ponded_area[ui] == 0.0 {
        return get_surf_area(nodes, idx, depth, tables, unit_sys, subs);
    }
    let mut a = nodes.ponded_area[ui];
    if a <= 0.0 {
        a = get_surf_area(nodes, idx, nodes.full_depth[ui], tables, unit_sys, subs);
    }
    a
}

/// node::getMaxOutflow: caps only nodes with a legacy full volume.
pub fn get_max_outflow(nodes: &NodeData, idx: i32, mut q: f64, dt: f64) -> f64 {
    let ui = idx as usize;
    if nodes.rpt_full_volume[ui] > 0.0 {
        let q_max = nodes.inflow[ui] + nodes.old_volume[ui] / dt;
        q = cmin(q, q_max);
    }
    cmax(0.0, q)
}

/// node::getOverflow
pub fn get_overflow(new_volume: f64, full_volume: f64, dt: f64) -> f64 {
    if new_volume > full_volume && dt > 0.0 {
        return (new_volume - full_volume) / dt;
    }
    0.0
}
