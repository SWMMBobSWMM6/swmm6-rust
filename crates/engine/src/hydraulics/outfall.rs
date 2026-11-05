// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Outfall.hpp and Outfall.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Outfall boundary conditions: the outfall-to-link cache (buildOutfallLinkMap) and its helper
// legacyOffset (2 November 2026), and the routing step's setAllOutfallDepths with its
// file-local buildXSectParams, getYnorm and outfallStatic (5 November 2026).
//
// NOT YET TRANSLATED: everything else Outfall.cpp defines (the outfall 2D feedback and route-to
// helpers outside the routing step).

use crate::core::cmath::{cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::unit_conversion as ucf;
use crate::data::forcing::ForcingMode;
use crate::data::link::{LinkType, XsectShape};
use crate::data::node::{NodeType, OutfallType};
use crate::data::table::table_lookup_cursor;
use crate::hydraulics::transect_types::N_TRANSECT_TBL;
use crate::hydraulics::xsect_batch::XSectParams;
use crate::hydraulics::xsection as xsect;

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

/// Outfall.cpp buildXSectParams: the batch shape code as stored, no yw_max; tabulated shapes
/// get their transect tables.
fn build_xsect_params(ctx: &SimulationContext, uk: usize) -> XSectParams<'_> {
    let mut xs = XSectParams::default();
    xs.type_ = ctx.links.xsect_batch_shape[uk];
    xs.y_full = ctx.links.xsect_y_full[uk];
    xs.a_full = ctx.links.xsect_a_full[uk];
    xs.w_max = ctx.links.xsect_w_max[uk];
    xs.r_full = ctx.links.xsect_r_full[uk];
    xs.s_full = ctx.links.xsect_s_full[uk];
    xs.s_max = ctx.links.xsect_s_max[uk];
    xs.y_bot = ctx.links.xsect_y_bot[uk];
    xs.a_bot = ctx.links.xsect_a_bot[uk];
    xs.s_bot = ctx.links.xsect_s_bot[uk];
    xs.r_bot = ctx.links.xsect_r_bot[uk];
    let ls = ctx.links.xsect_shape[uk];
    if ls == XsectShape::Irregular || ls == XsectShape::Custom || ls == XsectShape::StreetXsect {
        let ci = ctx.links.xsect_curve[uk];
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
    xs
}

/// Outfall.cpp getYnorm
fn get_ynorm(xs: &XSectParams, beta: f64, q_max: f64, mut q: f64) -> f64 {
    if beta <= 0.0 || q <= 0.0 {
        return 0.0;
    }
    if q > q_max && q_max > 0.0 {
        q = q_max;
    }
    let s = q / beta;
    let a = xsect::get_a_of_s(xs, s);
    xsect::get_y_of_a(xs, a)
}

/// setAllOutfallDepths: every outfall's stage at the END of the step (dt_routing; 0 at the
/// start of a step), from the depth-setting link's normal and critical depths.
pub fn set_all_outfall_depths(ctx: &mut SimulationContext, dt_routing: f64) {
    let unit_sys = ucf::get_unit_system(ctx.options.flow_units as i32);
    let ucf_len = ucf::UCF[ucf::LENGTH as usize][unit_sys as usize];

    let new_routing_ms = ctx.elapsed_ms + 1000.0 * dt_routing;
    const K_MSEC_PER_DAY: f64 = 86400000.0;
    let elapsed_days = new_routing_ms / K_MSEC_PER_DAY;

    let process_outfall = |ctx: &mut SimulationContext, j: i32| {
        let uj = j as usize;
        let mut depth;

        // outfallStatic
        let r = ctx.node_subtypes.outfall_row(j);
        let (bc_type, param, link_idx0, z0) = if r >= 0 {
            let ur = r as usize;
            let o = &ctx.node_subtypes.outfalls;
            (o.bc_type[ur], o.param[ur], o.link_idx[ur], o.link_offset[ur])
        } else {
            (OutfallType::Free, 0.0, -1, 0.0)
        };
        let mut link_idx = link_idx0;
        let mut z = z0;
        if link_idx < 0 && ctx.node_subtypes.outfall_row(j) < 0 {
            for k in 0..ctx.n_links() {
                let uk = k as usize;
                let n2 = ctx.links.node2[uk];
                let n1 = ctx.links.node1[uk];
                if n2 == j {
                    link_idx = k;
                    z = legacy_offset(ctx, uk, true);
                } else if n2 >= 0 && ctx.nodes.type_[n2 as usize] == NodeType::Outfall {
                    continue;
                } else if n1 == j {
                    link_idx = k;
                    z = legacy_offset(ctx, uk, false);
                }
            }
        }
        if link_idx < 0 {
            return;
        }

        let mut y_norm = 0.0;
        let mut y_crit = 0.0;
        if ctx.links.type_[link_idx as usize] == LinkType::Conduit {
            let uk = link_idx as usize;
            let cd = &ctx.link_subtypes.conduits;
            let cr = ctx.link_subtypes.conduit_row(link_idx);
            let barrels = cmax(if cr >= 0 { cd.barrels[cr as usize] } else { 1 }, 1);
            let q = ctx.links.flow[uk].abs() / (barrels as f64);
            let xs = build_xsect_params(ctx, uk);
            let beta_v = if cr >= 0 { cd.beta[cr as usize] } else { 0.0 };
            let qmax_v = if cr >= 0 { cd.q_max[cr as usize] } else { 0.0 };
            y_norm = get_ynorm(&xs, beta_v, qmax_v, q);
            y_crit = xsect::get_ycrit(&xs, q);
        }

        let inv = ctx.nodes.invert_elev[uj];
        match bc_type {
            OutfallType::Free => {
                depth = z + cmin(y_norm, y_crit);
            }
            OutfallType::Normal => {
                depth = z + y_norm;
            }
            OutfallType::Fixed => {
                let stage = param;
                y_crit = cmin(y_crit, y_norm);
                if y_crit + z + inv < stage {
                    depth = stage - inv;
                } else if z > 0.0 {
                    depth = if stage < inv + z { cmax(0.0, stage - inv) } else { z + y_crit };
                } else {
                    depth = y_crit;
                }
            }
            OutfallType::Tidal => {
                let curve_idx = param as i32;
                let mut stage = inv;
                if curve_idx >= 0 && (curve_idx as usize) < ctx.tables.tables.len() {
                    let curve = &mut ctx.tables.tables[curve_idx as usize];
                    let mut x = if curve.x.is_empty() { 0.0 } else { curve.x[0] };
                    x += (elapsed_days - elapsed_days.floor()) * 24.0;
                    stage = table_lookup_cursor(curve, x) / ucf_len;
                }
                y_crit = cmin(y_crit, y_norm);
                if y_crit + z + inv < stage {
                    depth = stage - inv;
                } else if z > 0.0 {
                    depth = if stage < inv + z { cmax(0.0, stage - inv) } else { z + y_crit };
                } else {
                    depth = y_crit;
                }
            }
            OutfallType::Timeseries => {
                let ts_idx = param as i32;
                let mut stage = inv;
                if ts_idx >= 0 && (ts_idx as usize) < ctx.tables.tables.len() {
                    let date = ctx.options.start_date + elapsed_days;
                    stage = table_lookup_cursor(&mut ctx.tables.tables[ts_idx as usize], date) / ucf_len;
                }
                y_crit = cmin(y_crit, y_norm);
                if y_crit + z + inv < stage {
                    depth = stage - inv;
                } else if z > 0.0 {
                    depth = if stage < inv + z { cmax(0.0, stage - inv) } else { z + y_crit };
                } else {
                    depth = y_crit;
                }
            }
        }

        let z_inv = ctx.nodes.invert_elev[uj];
        let mut h_standard = z_inv + depth;

        let orow = ctx.node_subtypes.outfall_row(j);
        let h_2d = if orow >= 0 { ctx.node_subtypes.outfalls.head_2d[orow as usize] } else { -1.0e30 };
        let ramp_2d = if orow >= 0 { ctx.node_subtypes.outfalls.ramp_2d[orow as usize] } else { 0.0 };
        if h_2d > z_inv && ramp_2d > 0.0 {
            let h_blend = h_standard + ramp_2d * cmax(0.0, h_2d - h_standard);
            if h_blend > h_standard {
                h_standard = h_blend;
                depth = h_standard - z_inv;
            }
        }

        if uj < ctx.forcing.node_head_boundary_mode.len() && ctx.forcing.node_head_boundary_mode[uj] == ForcingMode::Override {
            let prescribed_stage = ctx.forcing.node_head_boundary_value[uj];
            depth = cmax(prescribed_stage - z_inv, 0.0);
        }

        ctx.nodes.depth[uj] = depth;
        ctx.nodes.head[uj] = ctx.nodes.invert_elev[uj] + depth;
    };

    if ctx.node_subtypes.subtype_row.len() as i32 == ctx.n_nodes() {
        if ctx.node_subtypes.outfalls.link_idx.len() as i32 != ctx.node_subtypes.outfalls.count() {
            build_outfall_link_map(ctx);
        }
        for r in 0..ctx.node_subtypes.outfalls.count() {
            let n = ctx.node_subtypes.outfalls.node_idx[r as usize];
            process_outfall(ctx, n);
        }
    } else {
        for j in 0..ctx.n_nodes() {
            if ctx.nodes.type_[j as usize] == NodeType::Outfall {
                process_outfall(ctx, j);
            }
        }
    }
}
