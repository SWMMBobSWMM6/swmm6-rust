// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/handlers/InfraHandler.cpp and
// InfraHandler.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Section handlers for [STREETS], [INLETS], [INLET_USAGE], [INLET_JUNCTIONS] (declared in
// NodesHandler.hpp, implemented here as upstream does), [ADJUSTMENTS] and [EVENTS].

use crate::core::cmath::cmax;
use crate::core::context::{Event, SimulationContext};
use crate::core::error_codes::{
    format_error, ERR_DUP_NAME, ERR_IJ_EXTRA_TOKENS, ERR_ITEMS, ERR_KEYWORD, ERR_NAME, ERR_NUMBER,
};
use crate::core::options::Bytes;
use crate::data::infra::InletUsageStore;
use crate::data::name_index::ieq;
use crate::data::node::NodeType;
use crate::input::parse_utils::{add_unique, parse_date, parse_time_day_fraction, to_double, to_int};
use crate::input::section_parser::parse_section;
use crate::input::tokenizer;

// ============================================================================
// handle_streets()
// ============================================================================

pub fn handle_streets(ctx: &mut SimulationContext, lines: &[Bytes]) {
    // Format: ID  Tcrown  Hcurb  Sx  nRoad  (Hdep  Wg  Sides  Tback  Sback  nBack)
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 5 {
            continue;
        }

        // A second definition in any case spelling is a duplicate ID (ERR 207).
        let mut dup = false;
        for existing in &ctx.streets.names {
            if ieq(existing, &tok[0]) {
                dup = true;
                break;
            }
        }
        if dup {
            ctx.errors.push(format_error(ERR_DUP_NAME, &tok[0]));
            continue;
        }

        let s = &mut ctx.streets;
        s.names.push(tok[0].clone());
        s.t_crown.push(to_double(&tok[1], 0.0));
        s.h_curb.push(to_double(&tok[2], 0.0));
        s.sx.push(to_double(&tok[3], 0.0));
        s.n_road.push(to_double(&tok[4], 0.0));
        s.gutter_depres.push(if tok.len() > 5 { to_double(&tok[5], 0.0) } else { 0.0 });
        s.gutter_width.push(if tok.len() > 6 { to_double(&tok[6], 0.0) } else { 0.0 });
        s.sides.push(if tok.len() > 7 { to_int(&tok[7], 2) } else { 2 });
        s.back_width.push(if tok.len() > 8 { to_double(&tok[8], 0.0) } else { 0.0 });
        s.back_slope.push(if tok.len() > 9 { to_double(&tok[9], 0.0) } else { 0.0 });
        s.back_n.push(if tok.len() > 10 { to_double(&tok[10], 0.0) } else { 0.0 });
    }
}

// ============================================================================
// handle_inlets()
// ============================================================================

/// Grate bar patterns, legacy inlet.c GrateTypeWords order.
const K_GRATE_TYPE_WORDS: [&[u8]; 8] = [
    b"P_BAR-50",
    b"P_BAR-50x100",
    b"P_BAR-30",
    b"CURVED_VANE",
    b"TILT_BAR-45",
    b"TILT_BAR-30",
    b"RETICULINE",
    b"GENERIC",
];

/// Curb-opening throat orientations, legacy inlet.c ThroatAngleWords order.
const K_THROAT_WORDS: [&[u8]; 3] = [b"HORIZONTAL", b"INCLINED", b"VERTICAL"];

/// Placement keywords of the usage tail.
const K_PLACEMENT_WORDS: [&[u8]; 3] = [b"AUTOMATIC", b"ON_GRADE", b"ON_SAG"];

/// First ASCII-insensitive whole-word match in `words`, or -1.
fn find_word(words: &[&[u8]], w: &[u8]) -> i32 {
    for (i, word) in words.iter().enumerate() {
        if ieq(word, w) {
            return i as i32;
        }
    }
    -1
}

/// First ASCII-insensitive inlet-design identifier match, or -1.
fn find_inlet_design(ctx: &SimulationContext, name: &[u8]) -> i32 {
    for i in 0..ctx.inlets.count() {
        if ieq(&ctx.inlets.names[i as usize], name) {
            return i;
        }
    }
    -1
}

/// The optional tail shared by [INLET_USAGE] tokens 4-9 and [INLET_JUNCTIONS] tokens 6-11.
struct UsageTail {
    num_inlets: i32,
    placement: i32, // 0=AUTOMATIC 1=ON_GRADE 2=ON_SAG
    clog_factor: f64,
    flow_limit: f64,
    local_depress: f64,
    local_width: f64,
}

impl Default for UsageTail {
    fn default() -> Self {
        UsageTail {
            num_inlets: 1,
            placement: 0,
            clog_factor: 1.0,
            flow_limit: 0.0,
            local_depress: 0.0,
            local_width: 0.0,
        }
    }
}

/// Parse the optional inlet count, clogging, flow limit, local geometry and placement.
/// Returns false after appending the first refusal; earlier fields stay written.
fn parse_inlet_usage_tail(ctx: &mut SimulationContext, tok: &[Bytes], first: usize, out: &mut UsageTail) -> bool {
    let num = |errors: &mut Vec<Bytes>, k: usize, v: &mut f64, lo: f64, hi: f64| -> bool {
        *v = to_double(&tok[k], -1.0);
        if *v < lo || *v > hi {
            errors.push(format_error(ERR_NUMBER, &tok[k]));
            return false;
        }
        true
    };

    if tok.len() > first {
        out.num_inlets = to_int(&tok[first], 0);
        if out.num_inlets < 1 {
            ctx.errors.push(format_error(ERR_NUMBER, &tok[first]));
            return false;
        }
    }
    if tok.len() > first + 1 {
        let mut pct = 0.0;
        if !num(&mut ctx.errors, first + 1, &mut pct, 0.0, 99.0) {
            return false;
        }
        out.clog_factor = 1.0 - pct / 100.0;
    }
    if tok.len() > first + 2 && !num(&mut ctx.errors, first + 2, &mut out.flow_limit, 0.0, 1.0e30) {
        return false;
    }
    if tok.len() > first + 3 && !num(&mut ctx.errors, first + 3, &mut out.local_depress, 0.0, 1.0e30) {
        return false;
    }
    if tok.len() > first + 4 && !num(&mut ctx.errors, first + 4, &mut out.local_width, 0.0, 1.0e30) {
        return false;
    }
    if tok.len() > first + 5 {
        out.placement = find_word(&K_PLACEMENT_WORDS, &tok[first + 5]);
        if out.placement < 0 {
            ctx.errors.push(format_error(ERR_KEYWORD, &tok[first + 5]));
            return false;
        }
    }
    true
}

/// Copy an already parsed optional tail into an allocated usage row.
fn apply_usage_tail(u: &mut InletUsageStore, row: i32, t: &UsageTail) {
    let ur = row as usize;
    u.num_inlets[ur] = t.num_inlets;
    u.placement[ur] = t.placement;
    u.clog_factor[ur] = t.clog_factor;
    u.flow_limit[ur] = t.flow_limit;
    u.local_depress[ur] = t.local_depress;
    u.local_width[ur] = t.local_width;
}

pub fn handle_inlets(ctx: &mut SimulationContext, lines: &[Bytes]) {
    // A GRATE line and a CURB line sharing one name MERGE into one COMBO row; any other
    // repeat of a name is a duplicate ID.
    for pl in parse_section(lines) {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.len() < 3 {
            continue;
        }

        let type_ = tokenizer::to_upper(&tok[1]);
        let is_grate = type_ == b"GRATE" || type_ == b"DROP_GRATE";
        let is_curb = type_ == b"CURB" || type_ == b"DROP_CURB";
        let is_slotted = type_ == b"SLOTTED";
        let is_custom = type_ == b"CUSTOM";
        if !is_grate && !is_curb && !is_slotted && !is_custom {
            ctx.errors.push(format_error(ERR_KEYWORD, &tok[1]));
            continue;
        }

        // --- Duplicate name: merge a GRATE/CURB pair, reject anything else ---
        let prev = find_inlet_design(ctx, &tok[0]);
        if prev >= 0 {
            let prev_type = &ctx.inlets.inlet_type[prev as usize];
            let combo = (prev_type.as_slice() == b"GRATE" && type_ == b"CURB")
                || (prev_type.as_slice() == b"CURB" && type_ == b"GRATE");
            if !combo {
                ctx.errors.push(format_error(ERR_DUP_NAME, &tok[0]));
                continue;
            }
        }

        // --- Type-specific token count + range checks (legacy read*Params) ---
        let (mut length, mut width, mut open_area, mut splash) = (0.0, 0.0, 0.0, 0.0);
        let mut grate_type: i32 = -1;
        let mut throat: i32 = 2; // default VERTICAL

        if is_grate || is_curb || is_slotted {
            let need: usize = if is_grate { 5 } else { 4 };
            if tok.len() < need {
                ctx.errors.push(format_error(ERR_ITEMS, &tok[0]));
                continue;
            }
            length = to_double(&tok[2], -1.0);
            width = to_double(&tok[3], -1.0); // CURB/DROP_CURB: opening height
            if length <= 0.0 {
                ctx.errors.push(format_error(ERR_NUMBER, &tok[2]));
                continue;
            }
            if width <= 0.0 {
                ctx.errors.push(format_error(ERR_NUMBER, &tok[3]));
                continue;
            }
        }
        if is_grate {
            grate_type = find_word(&K_GRATE_TYPE_WORDS, &tok[4]);
            if grate_type < 0 {
                ctx.errors.push(format_error(ERR_KEYWORD, &tok[4]));
                continue;
            }
            if grate_type == 7 {
                // GENERIC: open area is required
                if tok.len() < 6 {
                    ctx.errors.push(format_error(ERR_ITEMS, &tok[0]));
                    continue;
                }
                open_area = to_double(&tok[5], -1.0);
                if open_area <= 0.0 || open_area > 1.0 {
                    ctx.errors.push(format_error(ERR_NUMBER, &tok[5]));
                    continue;
                }
                if tok.len() > 6 {
                    splash = to_double(&tok[6], -1.0);
                    if splash < 0.0 {
                        ctx.errors.push(format_error(ERR_NUMBER, &tok[6]));
                        continue;
                    }
                }
            }
        }
        if type_ == b"CURB" && tok.len() > 4 {
            // DROP_CURB has no throat token
            throat = find_word(&K_THROAT_WORDS, &tok[4]);
            if throat < 0 {
                ctx.errors.push(format_error(ERR_KEYWORD, &tok[4]));
                continue;
            }
        }

        // --- Write (new row, or merge into the existing half of a COMBO) ---
        let idx = if prev >= 0 { prev } else { ctx.inlets.add_row(&tok[0], &type_) };
        let ui = idx as usize;
        if prev >= 0 {
            ctx.inlets.inlet_type[ui] = b"COMBO".to_vec();
        }

        if is_grate {
            ctx.inlets.length[ui] = length;
            ctx.inlets.width[ui] = width;
            ctx.inlets.grate_type[ui] = K_GRATE_TYPE_WORDS[grate_type as usize].to_vec();
            ctx.inlets.open_area[ui] = open_area;
            ctx.inlets.splash_veloc[ui] = splash;
        } else if is_curb {
            ctx.inlets.curb_length[ui] = length;
            ctx.inlets.curb_height[ui] = width;
            ctx.inlets.curb_throat[ui] = throat;
        } else if is_slotted {
            ctx.inlets.length[ui] = length;
            ctx.inlets.width[ui] = width;
        } else {
            // CUSTOM — the curve is resolved in PostParseResolver
            ctx.inlets.curve_id[ui] = tok[2].clone();
        }
        if !pl.comment.is_empty() && ctx.inlets.comments[ui].is_empty() {
            ctx.inlets.comments[ui] = pl.comment.clone();
        }
    }
}

// ============================================================================
// handle_inlet_usage()
// ============================================================================

pub fn handle_inlet_usage(ctx: &mut SimulationContext, lines: &[Bytes]) {
    // Format: linkID  inletID  nodeID  (#Inlets  %Clog  Qmax  aLocal  wLocal  placement)
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 3 {
            continue;
        }

        let link_idx = ctx.link_names.find(&tok[0]);
        if link_idx < 0 {
            ctx.errors.push(format_error(ERR_NAME, &tok[0]));
            continue;
        }

        let design_idx = find_inlet_design(ctx, &tok[1]);
        if design_idx < 0 {
            ctx.errors.push(format_error(ERR_NAME, &tok[1]));
            continue;
        }

        // The capture node may be declared after this section: deferred by name.
        let node_idx = ctx.node_names.find(&tok[2]);

        let mut tail = UsageTail::default();
        if !parse_inlet_usage_tail(ctx, &tok, 3, &mut tail) {
            continue;
        }

        // street_index stays -1 here; resolved in PostParseResolver.
        let row = ctx.inlet_usages.add_row(link_idx, -1, design_idx, node_idx);
        apply_usage_tail(&mut ctx.inlet_usages, row, &tail);
        if node_idx < 0 {
            ctx.inlet_usages.pending_capture_name[row as usize] = tok[2].clone();
        }
    }
}

// ============================================================================
// handle_inlet_junctions()
// ============================================================================

pub fn handle_inlet_junctions(ctx: &mut SimulationContext, lines: &[Bytes]) {
    // Format: Name Elev MaxDepth Inlet CaptureNode
    //         (#Inlets %Clog Qmax aLocal wLocal Placement)
    for pl in parse_section(lines) {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.len() < 5 {
            continue;
        }

        let name = &tok[0];
        if tok.len() > 11 {
            ctx.errors.push(format_error(ERR_IJ_EXTRA_TOKENS, name));
            continue;
        }

        let mut tail = UsageTail::default();
        if !parse_inlet_usage_tail(ctx, &tok, 5, &mut tail) {
            continue;
        }

        let idx = add_unique(&mut ctx.node_names, name, &mut ctx.errors);
        if idx < 0 {
            continue; // duplicate ID (ERR 207, legacy input.c parity)
        }

        ctx.nodes.grow_to(idx + 1);
        ctx.node_subtypes.set_node_type(&mut ctx.nodes, idx, NodeType::Junction);
        let ui = idx as usize;
        ctx.nodes.is_virtual[ui] = 1;
        ctx.nodes.is_inlet[ui] = 1;
        ctx.nodes.invert_elev[ui] = to_double(&tok[1], 0.0);
        // MaxDepth: flood threshold above the street invert, carried in rim_depth; a
        // missing or negative value collapses to 0 = derive from the street section.
        ctx.nodes.rim_depth[ui] = cmax(0.0, to_double(&tok[2], 0.0));
        if !pl.comment.is_empty() {
            ctx.nodes.comments[ui] = pl.comment.clone();
        }

        // Design and capture node are resolved by name in PostParseResolver.
        let row = ctx.inlet_usages.add_row(-1, idx, -1, -1);
        apply_usage_tail(&mut ctx.inlet_usages, row, &tail);
        let ur = row as usize;
        ctx.inlet_usages.pending_design_name[ur] = tok[3].clone();
        ctx.inlet_usages.pending_capture_name[ur] = tok[4].clone();
    }
}

// ============================================================================
// handle_adjustments()
// ============================================================================

pub fn handle_adjustments(ctx: &mut SimulationContext, lines: &[Bytes]) {
    // Ensure subcatchment pattern vectors are sized
    let ns = ctx.subcatch_names.size() as usize;
    if ctx.subcatch_n_perv_pattern.len() < ns {
        ctx.subcatch_n_perv_pattern.resize(ns, -1);
    }
    if ctx.subcatch_d_store_pattern.len() < ns {
        ctx.subcatch_d_store_pattern.resize(ns, -1);
    }
    if ctx.subcatch_infil_pattern.len() < ns {
        ctx.subcatch_infil_pattern.resize(ns, -1);
    }

    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 2 {
            continue;
        }

        let keyword = tokenizer::to_upper(&tok[0]);
        // legacy climate_readAdjustments uses match(): the keyword is a PREFIX of the token.
        let kw = |k: &[u8]| keyword.len() >= k.len() && keyword.starts_with(k);

        // Monthly arrays: TEMP, EVAP, RAIN, CONDUCT (need 13 tokens total)
        if kw(b"TEMP") && tok.len() >= 13 {
            for i in 0..12usize {
                ctx.adjust_temp[i] = to_double(&tok[i + 1], 0.0);
            }
        } else if kw(b"EVAP") && tok.len() >= 13 {
            for i in 0..12usize {
                ctx.adjust_evap[i] = to_double(&tok[i + 1], 0.0);
            }
        } else if kw(b"RAIN") && tok.len() >= 13 {
            for i in 0..12usize {
                ctx.adjust_rain[i] = to_double(&tok[i + 1], 0.0);
            }
        } else if kw(b"CONDUCT") && tok.len() >= 13 {
            for i in 0..12usize {
                let v = to_double(&tok[i + 1], 0.0);
                ctx.adjust_hydcon[i] = if v <= 0.0 { 1.0 } else { v };
            }
        }
        // Subcatchment pattern assignments: N-PERV, DSTORE, INFIL. The third token names
        // a TIME PATTERN.
        else if kw(b"N-PERV") && tok.len() >= 3 {
            let si = ctx.subcatch_names.find(&tok[1]);
            let pi = ctx.patterns.find(&tok[2]);
            if si >= 0 && pi >= 0 {
                let usi = si as usize;
                if usi >= ctx.subcatch_n_perv_pattern.len() {
                    ctx.subcatch_n_perv_pattern.resize(usi + 1, -1);
                }
                ctx.subcatch_n_perv_pattern[usi] = pi;
            }
        } else if kw(b"DSTORE") && tok.len() >= 3 {
            let si = ctx.subcatch_names.find(&tok[1]);
            let pi = ctx.patterns.find(&tok[2]);
            if si >= 0 && pi >= 0 {
                let usi = si as usize;
                if usi >= ctx.subcatch_d_store_pattern.len() {
                    ctx.subcatch_d_store_pattern.resize(usi + 1, -1);
                }
                ctx.subcatch_d_store_pattern[usi] = pi;
            }
        } else if kw(b"INFIL") && tok.len() >= 3 {
            let si = ctx.subcatch_names.find(&tok[1]);
            let pi = ctx.patterns.find(&tok[2]);
            if si >= 0 && pi >= 0 {
                let usi = si as usize;
                if usi >= ctx.subcatch_infil_pattern.len() {
                    ctx.subcatch_infil_pattern.resize(usi + 1, -1);
                }
                ctx.subcatch_infil_pattern[usi] = pi;
            }
        }
    }
}

// ============================================================================
// handle_events()
// ============================================================================

pub fn handle_events(ctx: &mut SimulationContext, lines: &[Bytes]) {
    // Format: StartDate  StartTime  EndDate  EndTime
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 4 {
            continue;
        }

        let start = parse_date(&tok[0]) + parse_time_day_fraction(&tok[1]);
        let end = parse_date(&tok[2]) + parse_time_day_fraction(&tok[3]);

        if start >= end {
            continue;
        }

        ctx.events.push(Event { start, end });
    }
}
