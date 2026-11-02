// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/handlers/NodesHandler.cpp and
// NodesHandler.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Section handlers for [JUNCTIONS], [VIRTUAL_JUNCTIONS], [OUTFALLS], [DIVIDERS], [STORAGE]
// and [COORDINATES]. [INLET_JUNCTIONS] is declared in NodesHandler.hpp but implemented in
// InfraHandler.cpp upstream; it lives in infra_handler.rs here. The ponded-area check of
// [JUNCTIONS] calls std::strtod on the token's C string, which is core::cstdlib::strtod on
// the bytes before the first NUL.

use crate::core::cmath::cmax;
use crate::core::context::SimulationContext;
use crate::core::cstdlib::strtod;
use crate::core::error_codes::{format_error, ERR_NUMBER, ERR_VJ_EXTRA_TOKENS};
use crate::core::options::Bytes;
use crate::data::node::{DividerType, NodeType, OutfallType, StorageShape};
use crate::data::storage_geometry::{storage_shape_coeffs, storage_shape_from_keyword};
use crate::input::parse_utils::{add_unique, to_double};
use crate::input::section_parser::parse_section;
use crate::input::tokenizer;

/// Grow base-node arrays through a registered index without resetting prefixes.
fn ensure_node_capacity(ctx: &mut SimulationContext, idx: i32) {
    ctx.nodes.grow_to(idx + 1);
}

/// Extend only node X/Y coordinate columns to an expected count.
fn ensure_spatial_capacity(ctx: &mut SimulationContext, n_nodes: i32) {
    let un = n_nodes as usize;
    if ctx.spatial.node_x.len() < un {
        ctx.spatial.node_x.resize(un, 0.0);
    }
    if ctx.spatial.node_y.len() < un {
        ctx.spatial.node_y.resize(un, 0.0);
    }
}

/// The C string of a token: its bytes up to the first NUL.
fn c_str(s: &[u8]) -> &[u8] {
    match s.iter().position(|&c| c == 0) {
        Some(p) => &s[..p],
        None => s,
    }
}

// ============================================================================
// handle_junctions()
// ============================================================================

pub fn handle_junctions(ctx: &mut SimulationContext, lines: &[Bytes]) {
    // Pre-reserve from the section's row count (capacity only).
    ctx.nodes.reserve_to(ctx.nodes.count() + lines.len() as i32);
    ctx.node_names.reserve(ctx.node_names.size() as usize + lines.len());
    for pl in parse_section(lines) {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.len() < 2 {
            continue;
        }

        let name = &tok[0];

        // Register name; a second definition (any case spelling) is ERR 207
        let idx = add_unique(&mut ctx.node_names, name, &mut ctx.errors);
        if idx < 0 {
            continue; // duplicate ID (ERR 207, legacy input.c parity)
        }

        ensure_node_capacity(ctx, idx);

        // legacy junc_readParams: max depth, initial depth, surcharge depth and ponded
        // area may not be negative (ERROR 211 naming the token, in column order).
        let mut negative_param = false;
        let mut i = 2usize;
        while i < tok.len() && i <= 5 {
            if to_double(&tok[i], 0.0) < 0.0 {
                ctx.errors.push(format_error(ERR_NUMBER, &tok[i]));
                negative_param = true;
                break;
            }
            i += 1;
        }
        if negative_param {
            continue;
        }

        let mut area = 0.0;
        if tok.len() > 6 {
            let cs = c_str(&tok[6]);
            let (v, used) = strtod(cs);
            area = v;
            if used == 0 || used != cs.len() || !area.is_finite() || area < 0.0 {
                ctx.errors.push(format_error(ERR_NUMBER, &tok[6]));
                continue;
            }
        }
        let ui = idx as usize;
        ctx.nodes.storage_area[ui] = area;
        ctx.node_subtypes.set_node_type(&mut ctx.nodes, idx, NodeType::Junction);
        ctx.nodes.invert_elev[ui] = to_double(&tok[1], 0.0); // Elev
        if tok.len() > 2 {
            ctx.nodes.full_depth[ui] = to_double(&tok[2], 0.0); // MaxDepth
        }
        if tok.len() > 3 {
            ctx.nodes.init_depth[ui] = to_double(&tok[3], 0.0); // InitDepth
        }
        if tok.len() > 4 {
            ctx.nodes.sur_depth[ui] = to_double(&tok[4], 0.0); // SurDepth
        }
        if tok.len() > 5 {
            ctx.nodes.ponded_area[ui] = to_double(&tok[5], 0.0); // Aponded
        }
        if !pl.comment.is_empty() {
            ctx.nodes.comments[ui] = pl.comment.clone();
        }
    }
}

// ============================================================================
// handle_virtual_junctions()
// ============================================================================

pub fn handle_virtual_junctions(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for pl in parse_section(lines) {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.len() < 2 {
            continue;
        }

        let name = &tok[0];

        // Name, invert elevation and an OPTIONAL rendering-only max depth; a fourth
        // token is an error.
        if tok.len() > 3 {
            ctx.errors.push(format_error(ERR_VJ_EXTRA_TOKENS, name));
            continue;
        }

        let idx = add_unique(&mut ctx.node_names, name, &mut ctx.errors);
        if idx < 0 {
            continue; // duplicate ID (ERR 207, legacy input.c parity)
        }

        ensure_node_capacity(ctx, idx);

        let ui = idx as usize;
        ctx.node_subtypes.set_node_type(&mut ctx.nodes, idx, NodeType::Junction);
        ctx.nodes.is_virtual[ui] = 1;
        ctx.nodes.invert_elev[ui] = to_double(&tok[1], 0.0);
        // MaxDepth: rim/ground elevation for drawings only. Negatives and unparseable
        // text collapse to 0 = unset.
        if tok.len() > 2 {
            ctx.nodes.rim_depth[ui] = cmax(0.0, to_double(&tok[2], 0.0));
        }
        if !pl.comment.is_empty() {
            ctx.nodes.comments[ui] = pl.comment.clone();
        }
    }
}

// ============================================================================
// handle_outfalls()
// ============================================================================

pub fn handle_outfalls(ctx: &mut SimulationContext, lines: &[Bytes]) {
    // Pre-reserve from the section's row count (capacity only).
    ctx.nodes.reserve_to(ctx.nodes.count() + lines.len() as i32);
    ctx.node_names.reserve(ctx.node_names.size() as usize + lines.len());
    for pl in parse_section(lines) {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.len() < 3 {
            continue;
        }

        let name = &tok[0];
        let idx = add_unique(&mut ctx.node_names, name, &mut ctx.errors);
        if idx < 0 {
            continue; // duplicate ID (ERR 207, legacy input.c parity)
        }

        ensure_node_capacity(ctx, idx);

        let orow = ctx.node_subtypes.set_node_type(&mut ctx.nodes, idx, NodeType::Outfall) as usize;
        ctx.nodes.invert_elev[idx as usize] = to_double(&tok[1], 0.0); // Elev
        let o = &mut ctx.node_subtypes.outfalls;

        // Type: FREE, NORMAL, FIXED, TIDAL, TIMESERIES
        let otype = tokenizer::to_upper(&tok[2]);
        if otype == b"FREE" {
            o.bc_type[orow] = OutfallType::Free;
        } else if otype == b"NORMAL" {
            o.bc_type[orow] = OutfallType::Normal;
        } else if otype == b"FIXED" {
            o.bc_type[orow] = OutfallType::Fixed;
        } else if otype == b"TIDAL" {
            o.bc_type[orow] = OutfallType::Tidal;
        } else if otype == b"TIMESERIES" {
            o.bc_type[orow] = OutfallType::Timeseries;
        }

        // FREE/NORMAL carry no stage-data field, so the gate/routeTo columns shift left
        // by one for those types.
        let has_stage_field = o.bc_type[orow] == OutfallType::Fixed
            || o.bc_type[orow] == OutfallType::Tidal
            || o.bc_type[orow] == OutfallType::Timeseries;
        let mut next = 3usize;

        // FREE/NORMAL: skip a "*" placeholder stage column.
        if !has_stage_field && tok.len() > next && tok[next] == b"*" {
            next += 1;
        }

        if has_stage_field && tok.len() > next {
            if o.bc_type[orow] == OutfallType::Fixed {
                o.param[orow] = to_double(&tok[next], 0.0);
            } else {
                // TIDAL / TIMESERIES: keep the name; -1 marks "unresolved".
                o.param_name[orow] = tok[next].clone();
                o.param[orow] = -1.0;
            }
            next += 1;
        }

        // Gated (YES/NO)
        if tok.len() > next {
            o.has_flap_gate[orow] = tokenizer::parse_boolean(&tok[next]) as u8;
            next += 1;
        }

        // Route-to subcatchment (optional last field)
        if tok.len() > next && !tok[next].is_empty() && tok[next] != b"*" {
            let sc = ctx.subcatch_names.find(&tok[next]);
            ctx.node_subtypes.outfalls.route_to[orow] = sc; // may be -1 if not yet parsed
        }
        if !pl.comment.is_empty() {
            ctx.nodes.comments[idx as usize] = pl.comment.clone();
        }
    }
}

// ============================================================================
// handle_dividers()
// ============================================================================

pub fn handle_dividers(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for pl in parse_section(lines) {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.len() < 4 {
            continue;
        }

        let name = &tok[0];
        let idx = add_unique(&mut ctx.node_names, name, &mut ctx.errors);
        if idx < 0 {
            continue; // duplicate ID (ERR 207, legacy input.c parity)
        }

        ensure_node_capacity(ctx, idx);
        let ui = idx as usize;

        let drow = ctx.node_subtypes.set_node_type(&mut ctx.nodes, idx, NodeType::Divider) as usize;
        ctx.nodes.invert_elev[ui] = to_double(&tok[1], 0.0);

        // tok[2] = diversion link name (resolved in post-parse)
        let div_link_name = &tok[2];
        ctx.node_subtypes.dividers.link_name[drow] = div_link_name.clone();
        let dl = ctx.link_names.find(div_link_name);
        ctx.node_subtypes.dividers.link[drow] = dl; // may be -1 if not yet parsed

        // tok[3] = divider type
        let dtype = tokenizer::to_upper(&tok[3]);
        if dtype == b"CUTOFF" {
            ctx.node_subtypes.dividers.method[drow] = DividerType::Cutoff;
            if tok.len() > 4 {
                ctx.node_subtypes.dividers.cutoff[drow] = to_double(&tok[4], 0.0);
            }
        } else if dtype == b"OVERFLOW" {
            ctx.node_subtypes.dividers.method[drow] = DividerType::OverflowDiv;
        } else if dtype == b"TABULAR" {
            ctx.node_subtypes.dividers.method[drow] = DividerType::Tabular;
            // tok[4] = curve name (deferred)
            if tok.len() > 4 {
                ctx.node_subtypes.dividers.curve_name[drow] = tok[4].clone();
                let ci = ctx.find_curve(&tok[4]);
                ctx.node_subtypes.dividers.curve[drow] = ci; // may be -1
            }
        } else if dtype == b"WEIR" {
            ctx.node_subtypes.dividers.method[drow] = DividerType::Weir;
            // Legacy column order: qMin dhMax cWeir.
            if tok.len() > 4 {
                ctx.node_subtypes.dividers.cutoff[drow] = to_double(&tok[4], 0.0);
            }
            if tok.len() > 5 {
                ctx.node_subtypes.dividers.max_depth[drow] = to_double(&tok[5], 0.0);
            }
            if tok.len() > 6 {
                ctx.node_subtypes.dividers.cd[drow] = to_double(&tok[6], 0.0);
            }
        }

        // The optional tail "maxDepth initDepth surDepth aPond" follows the
        // type-specific parameters: OVERFLOW none, CUTOFF/TABULAR one, WEIR three.
        let mut tail: i32 = 5; // CUTOFF, TABULAR
        if dtype == b"OVERFLOW" {
            tail = 4;
        } else if dtype == b"WEIR" {
            tail = 7;
        }
        let ntoks = tok.len() as i32;
        if ntoks > tail {
            ctx.nodes.full_depth[ui] = to_double(&tok[tail as usize], 0.0);
        }
        if ntoks > tail + 1 {
            ctx.nodes.init_depth[ui] = to_double(&tok[(tail + 1) as usize], 0.0);
        }
        if ntoks > tail + 2 {
            ctx.nodes.sur_depth[ui] = to_double(&tok[(tail + 2) as usize], 0.0);
        }
        if ntoks > tail + 3 {
            ctx.nodes.ponded_area[ui] = to_double(&tok[(tail + 3) as usize], 0.0);
        }
        if !pl.comment.is_empty() {
            ctx.nodes.comments[ui] = pl.comment.clone();
        }
    }
}

// ============================================================================
// handle_storage()
// ============================================================================

pub fn handle_storage(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for pl in parse_section(lines) {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.len() < 4 {
            continue;
        }

        let name = &tok[0];
        let idx = add_unique(&mut ctx.node_names, name, &mut ctx.errors);
        if idx < 0 {
            continue; // duplicate ID (ERR 207, legacy input.c parity)
        }

        ensure_node_capacity(ctx, idx);
        let ui = idx as usize;

        let srow = ctx.node_subtypes.set_node_type(&mut ctx.nodes, idx, NodeType::Storage) as usize;
        ctx.nodes.invert_elev[ui] = to_double(&tok[1], 0.0); // Elev
        ctx.nodes.full_depth[ui] = to_double(&tok[2], 0.0); // MaxDepth
        ctx.nodes.init_depth[ui] = to_double(&tok[3], 0.0); // InitDepth

        if tok.len() < 5 {
            if !pl.comment.is_empty() {
                ctx.nodes.comments[ui] = pl.comment.clone();
            }
            continue;
        }
        let shape = tokenizer::to_upper(&tok[4]);

        // Unknown keyword => leave the row at its FUNCTIONAL default.
        let mut sshape = StorageShape::Functional;
        storage_shape_from_keyword(&shape, &mut sshape);
        let s = &mut ctx.node_subtypes.storages;
        s.shape[srow] = sshape;

        if sshape == StorageShape::Tabular {
            // Next token is curve name — resolved in the post-parse pass
            s.curve[srow] = -1;
            if tok.len() > 5 {
                s.curve_name[srow] = tok[5].clone();
            }
        } else if sshape == StorageShape::Functional {
            // A1, A2, A0
            if tok.len() > 5 {
                s.a[srow] = to_double(&tok[5], 0.0);
            }
            if tok.len() > 6 {
                s.b[srow] = to_double(&tok[6], 0.0);
            }
            if tok.len() > 7 {
                s.c[srow] = to_double(&tok[7], 0.0);
            }
        } else {
            // CYLINDRICAL / CONICAL / PARABOLIC / PYRAMIDAL — three raw dimensions kept
            // verbatim, plus the derived quadratic area coefficients.
            if tok.len() > 5 {
                s.p1[srow] = to_double(&tok[5], 0.0);
            }
            if tok.len() > 6 {
                s.p2[srow] = to_double(&tok[6], 0.0);
            }
            if tok.len() > 7 {
                s.p3[srow] = to_double(&tok[7], 0.0);
            }
            let (mut a, mut b, mut c) = (0.0, 0.0, 0.0);
            if storage_shape_coeffs(sshape, s.p1[srow], s.p2[srow], s.p3[srow], &mut a, &mut b, &mut c) {
                s.a[srow] = a;
                s.b[srow] = b;
                s.c[srow] = c;
            }
            s.curve[srow] = -1;
        }

        // Optional: SurDepth, Fevap, then the exfiltration parameters (one bare Ksat, or
        // Psi Ksat IMD).
        let param_offset: i32 = if sshape == StorageShape::Tabular { 6 } else { 8 };
        if tok.len() as i32 > param_offset {
            ctx.nodes.sur_depth[ui] = to_double(&tok[param_offset as usize], 0.0);
        }
        let s = &mut ctx.node_subtypes.storages;
        if tok.len() as i32 > param_offset + 1 {
            s.evap_frac[srow] = to_double(&tok[(param_offset + 1) as usize], 0.0);
        }
        {
            let n = param_offset + 2;
            let ntoks = tok.len() as i32;
            if ntoks == n + 1 {
                s.exfil_suction[srow] = 0.0;
                s.exfil_ksat[srow] = to_double(&tok[n as usize], 0.0);
                s.exfil_imd[srow] = 0.0;
            } else if ntoks >= n + 3 {
                s.exfil_suction[srow] = to_double(&tok[n as usize], 0.0);
                s.exfil_ksat[srow] = to_double(&tok[(n + 1) as usize], 0.0);
                s.exfil_imd[srow] = to_double(&tok[(n + 2) as usize], 0.0);
            }
        }
        if !pl.comment.is_empty() {
            ctx.nodes.comments[ui] = pl.comment.clone();
        }
    }
}

// ============================================================================
// handle_coordinates()
// ============================================================================

pub fn handle_coordinates(ctx: &mut SimulationContext, lines: &[Bytes]) {
    let nn = ctx.node_names.size();
    ensure_spatial_capacity(ctx, nn);

    // Quoted node names go through the owning tokenizer, as upstream does. (Upstream
    // hoists the token buffer out of the loop to save allocations; here it is per row.)
    for line in lines {
        let quoted: Vec<Bytes>;
        let mut tok: Vec<&[u8]> = Vec::new();
        if !line.contains(&b'"') {
            tokenizer::tokenize_views_into(line, &mut tok);
        } else {
            quoted = tokenizer::tokenize(line);
            tok.extend(quoted.iter().map(|t| t.as_slice()));
        }
        if tok.len() < 3 {
            continue;
        }

        let idx = ctx.node_names.find(tok[0]);
        if idx < 0 {
            continue; // unknown node — silently skip
        }

        ensure_spatial_capacity(ctx, idx + 1);

        let ui = idx as usize;
        ctx.spatial.node_x[ui] = to_double(tok[1], 0.0);
        ctx.spatial.node_y[ui] = to_double(tok[2], 0.0);
        if ctx.spatial.node_has_xy.len() < ctx.spatial.node_x.len() {
            let n = ctx.spatial.node_x.len();
            ctx.spatial.node_has_xy.resize(n, 0);
        }
        ctx.spatial.node_has_xy[ui] = 1; // G-X2
    }
}
