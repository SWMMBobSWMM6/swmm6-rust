// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/edit/VirtualJunctionOps.cpp and
// VirtualJunctionOps.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The part of namespace openswmm::edit that input/PostParseResolver.cpp reaches: the
// virtual-junction rule check and derived geometry, and the inlet-junction host/shape
// helpers, with the anonymous-namespace helpers they use. The editing operations
// (set/clear virtual, split, fuse, set inlet) are not translated here.
//
// Rule 5 of vj_rule_violation inspects `ctx.twod_io.mesh`, which is non-null only when a
// 2D mesh is loaded. The 1D build never loads one (TwoDModelIO carries no mesh in this
// crate), so the branch cannot be taken and is left out.

use crate::core::context::SimulationContext;
use crate::core::error_codes::{ERR_VJ_LINK_COUNT, ERR_VJ_OFFSET, ERR_VJ_XSECT_MISMATCH};
use crate::data::link::{LinkType, XsectShape};

const VJ_OFFSET_TOL: f64 = 1.0e-9;

// ---- anonymous namespace -------------------------------------------------------

/// Count every link end touching `node_idx`; the first two touching links go to
/// `link_a` / `link_b`.
fn find_attached(ctx: &SimulationContext, node_idx: i32, link_a: &mut i32, link_b: &mut i32) -> i32 {
    *link_a = -1;
    *link_b = -1;
    let mut touches = 0;
    for j in 0..ctx.n_links() {
        let uj = j as usize;
        let ends = [ctx.links.node1[uj], ctx.links.node2[uj]];
        for n in ends {
            if n != node_idx {
                continue;
            }
            if touches == 0 {
                *link_a = j;
            } else if touches == 1 {
                *link_b = j;
            }
            touches += 1;
        }
    }
    touches
}

fn offset_at_node(ctx: &SimulationContext, j: i32, node_idx: i32) -> f64 {
    let uj = j as usize;
    if ctx.links.node1[uj] == node_idx {
        ctx.links.offset1[uj]
    } else {
        ctx.links.offset2[uj]
    }
}

fn ij_conduit_shape_ok(ctx: &SimulationContext, link_idx: i32, for_drop_inlet: bool) -> bool {
    if link_idx < 0 || link_idx >= ctx.n_links() {
        return false;
    }
    let shape = ctx.links.xsect_shape[link_idx as usize];
    if for_drop_inlet {
        return shape == XsectShape::RectOpen || shape == XsectShape::Trapezoidal;
    }
    shape == XsectShape::StreetXsect
}

fn ij_design_is_drop(ctx: &SimulationContext, design_idx: i32) -> bool {
    if design_idx < 0 || design_idx >= ctx.inlets.count() {
        return false;
    }
    let t = &ctx.inlets.inlet_type[design_idx as usize];
    t.as_slice() == b"DROP_GRATE" || t.as_slice() == b"DROP_CURB"
}

// ---- public -------------------------------------------------------------------

/// The first violated virtual-junction rule for node `node_idx` (its error code), or 0.
pub fn vj_rule_violation(ctx: &SimulationContext, node_idx: i32) -> i32 {
    if node_idx < 0 || node_idx >= ctx.n_nodes() {
        return ERR_VJ_LINK_COUNT;
    }

    let mut j1 = -1;
    let mut j2 = -1;
    let touches = find_attached(ctx, node_idx, &mut j1, &mut j2);

    // Rule 1: exactly two attached links, both conduits.
    if touches != 2
        || ctx.links.type_[j1 as usize] != LinkType::Conduit
        || ctx.links.type_[j2 as usize] != LinkType::Conduit
    {
        return ERR_VJ_LINK_COUNT;
    }

    let u1 = j1 as usize;
    let u2 = j2 as usize;

    // Rule 2: identical cross-section (shape, geoms, curve/transect, barrels).
    let cr1 = ctx.link_subtypes.conduit_row(j1);
    let cr2 = ctx.link_subtypes.conduit_row(j2);
    let cd = &ctx.link_subtypes.conduits;
    let b1 = if cr1 >= 0 { cd.barrels[cr1 as usize] } else { 1 };
    let b2 = if cr2 >= 0 { cd.barrels[cr2 as usize] } else { 1 };
    if ctx.links.xsect_shape[u1] != ctx.links.xsect_shape[u2]
        || ctx.links.xsect_geom1[u1] != ctx.links.xsect_geom1[u2]
        || ctx.links.xsect_geom2[u1] != ctx.links.xsect_geom2[u2]
        || ctx.links.xsect_geom3[u1] != ctx.links.xsect_geom3[u2]
        || ctx.links.xsect_geom4[u1] != ctx.links.xsect_geom4[u2]
        || ctx.links.xsect_curve[u1] != ctx.links.xsect_curve[u2]
        || b1 != b2
    {
        return ERR_VJ_XSECT_MISMATCH;
    }

    // Rule 3: zero offsets at the node.
    if offset_at_node(ctx, j1, node_idx).abs() > VJ_OFFSET_TOL || offset_at_node(ctx, j2, node_idx).abs() > VJ_OFFSET_TOL {
        return ERR_VJ_OFFSET;
    }

    // Rule 5 (no 2D surface coupling) applies only with a 2D mesh; see the header.
    0
}

/// Derived geometry of a virtual junction: full depth = the first attached link's crown,
/// no surcharge depth, no ponding, no storage area. rim_depth is left alone.
pub fn vj_apply_derived_geometry(ctx: &mut SimulationContext, node_idx: i32) {
    let mut j1 = -1;
    let mut j2 = -1;
    if find_attached(ctx, node_idx, &mut j1, &mut j2) < 1 || j1 < 0 {
        return;
    }
    let ui = node_idx as usize;
    ctx.nodes.full_depth[ui] = ctx.links.xsect_y_full[j1 as usize];
    ctx.nodes.sur_depth[ui] = 0.0;
    ctx.nodes.ponded_area[ui] = 0.0;
    ctx.nodes.storage_area[ui] = 0.0;
}

/// The conduit that hosts an inlet: for `host_kind` 0 the link itself when it is a
/// conduit; otherwise the attached link whose downstream end is the node, else the first.
pub fn ij_host_conduit(ctx: &SimulationContext, host_kind: i32, host_idx: i32) -> i32 {
    if host_kind == 0 {
        if host_idx < 0 || host_idx >= ctx.n_links() {
            return -1;
        }
        return if ctx.links.type_[host_idx as usize] == LinkType::Conduit { host_idx } else { -1 };
    }
    if host_idx < 0 || host_idx >= ctx.n_nodes() {
        return -1;
    }
    let mut a = -1;
    let mut b = -1;
    if find_attached(ctx, host_idx, &mut a, &mut b) == 0 {
        return -1;
    }
    if a >= 0 && ctx.links.node2[a as usize] == host_idx {
        return a;
    }
    if b >= 0 && ctx.links.node2[b as usize] == host_idx {
        return b;
    }
    a
}

/// Whether inlet design `design_idx` can sit on link `link_idx`'s cross-section.
pub fn ij_usage_shape_ok(ctx: &SimulationContext, design_idx: i32, link_idx: i32) -> bool {
    if design_idx < 0 || design_idx >= ctx.inlets.count() {
        return false;
    }
    if link_idx < 0 || link_idx >= ctx.n_links() {
        return false;
    }
    let t = &ctx.inlets.inlet_type[design_idx as usize];
    if t.as_slice() == b"CUSTOM" {
        return true;
    }
    ij_conduit_shape_ok(ctx, link_idx, ij_design_is_drop(ctx, design_idx))
}
