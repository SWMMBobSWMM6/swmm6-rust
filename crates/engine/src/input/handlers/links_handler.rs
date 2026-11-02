// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/handlers/LinksHandler.cpp and
// LinksHandler.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Section handlers for [CONDUITS], [PUMPS], [ORIFICES], [WEIRS], [OUTLETS], [XSECTIONS],
// [LOSSES] and [TRANSECTS]. `static_cast<int>(double)` on a parsed token (barrels, culvert
// code, outlet type) is core::cmath::cvt_i32, since the token may read as NaN or out of
// range.

use crate::core::cmath::cvt_i32;
use crate::core::constants;
use crate::core::context::SimulationContext;
use crate::core::error_codes::{format_error, format_error_detail, ERR_DUP_NAME, ERR_ITEMS, ERR_KEYWORD, ERR_NUMBER};
use crate::core::options::Bytes;
use crate::data::link::{LinkType, XsectShape};
use crate::data::name_index::ieq;
use crate::input::parse_utils::{add_unique, to_double};
use crate::input::section_parser::parse_section;
use crate::input::tokenizer;

/// Extend base-link arrays through a registered index.
fn ensure_link_capacity(ctx: &mut SimulationContext, idx: i32) {
    ctx.links.grow_to(idx + 1);
}

/// Resolve link end-node names now and retain them for post-parse replay.
fn set_link_nodes(ctx: &mut SimulationContext, idx: i32, n1: &[u8], n2: &[u8]) {
    ctx.links.node1[idx as usize] = ctx.node_names.find(n1);
    ctx.links.node2[idx as usize] = ctx.node_names.find(n2);
    ctx.pending_link_nodes.push((idx, (n1.to_vec(), n2.to_vec())));
}

/// Interpret an authored link offset or elevation sentinel: under LINK_OFFSETS ELEVATION
/// a leading '*' means "at the node invert" and is stored as MISSING.
fn offset_token(ctx: &SimulationContext, tok: &[u8]) -> f64 {
    if ctx.options.link_offsets == 1 && !tok.is_empty() && tok[0] == b'*' {
        return constants::MISSING;
    }
    to_double(tok, 0.0)
}

// ============================================================================
// handle_conduits()
// ============================================================================

pub fn handle_conduits(ctx: &mut SimulationContext, lines: &[Bytes]) {
    // Pre-reserve from the section's row count (capacity only).
    ctx.links.reserve_to(ctx.links.count() + lines.len() as i32);
    ctx.link_names.reserve(ctx.link_names.size() as usize + lines.len());
    for pl in parse_section(lines) {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.len() < 7 {
            continue; // Name Node1 Node2 Length Roughness In Out required
        }

        let name = &tok[0];
        let idx = add_unique(&mut ctx.link_names, name, &mut ctx.errors);
        if idx < 0 {
            continue; // duplicate ID (ERR 207, legacy input.c parity)
        }

        ensure_link_capacity(ctx, idx);

        let cr = ctx.link_subtypes.set_link_type(&mut ctx.links, idx, LinkType::Conduit);
        let ucr = cr as usize;
        let ui = idx as usize;

        set_link_nodes(ctx, idx, &tok[1], &tok[2]);
        ctx.link_subtypes.conduits.length[ucr] = to_double(&tok[3], 0.0);
        ctx.link_subtypes.conduits.roughness[ucr] = to_double(&tok[4], 0.0);
        ctx.links.offset1[ui] = offset_token(ctx, &tok[5]);
        ctx.links.offset2[ui] = offset_token(ctx, &tok[6]);
        if tok.len() > 7 {
            ctx.links.q0[ui] = to_double(&tok[7], 0.0);
        }
        if tok.len() > 8 {
            ctx.links.q_limit[ui] = to_double(&tok[8], 0.0);
        }
        if !pl.comment.is_empty() {
            ctx.links.comments[ui] = pl.comment.clone();
        }
    }
}

// ============================================================================
// handle_pumps()
// ============================================================================

pub fn handle_pumps(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for pl in parse_section(lines) {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.len() < 3 {
            continue;
        }

        let name = &tok[0];
        let idx = add_unique(&mut ctx.link_names, name, &mut ctx.errors);
        if idx < 0 {
            continue; // duplicate ID (ERR 207, legacy input.c parity)
        }

        ensure_link_capacity(ctx, idx);

        let pr = ctx.link_subtypes.set_link_type(&mut ctx.links, idx, LinkType::Pump);
        let upr = pr as usize;
        let ui = idx as usize;
        set_link_nodes(ctx, idx, &tok[1], &tok[2]);
        // tok[3]: pump curve name; "*" is the ideal-pump placeholder (curve stays -1).
        if tok.len() > 3 && tok[3] != b"*" {
            ctx.links.pump_curve_name[ui] = tok[3].clone();
            ctx.link_subtypes.pumps.curve[upr] = ctx.find_curve(&tok[3]);
        }
        // tok[4]: init status (ON/OFF)
        if tok.len() > 4 {
            let on = tokenizer::to_upper(&tok[4]) == b"ON";
            ctx.link_subtypes.pumps.init_state[upr] = if on { 1u8 } else { 0u8 };
            let init_val = if on { 1.0 } else { 0.0 };
            ctx.links.setting[ui] = init_val;
            ctx.links.target_setting[ui] = init_val;
        }
        // tok[5]: startup depth, tok[6]: shutoff depth
        if tok.len() > 5 {
            ctx.link_subtypes.pumps.startup[upr] = to_double(&tok[5], 0.0);
        }
        if tok.len() > 6 {
            ctx.link_subtypes.pumps.shutoff[upr] = to_double(&tok[6], 0.0);
        }
        if !pl.comment.is_empty() {
            ctx.links.comments[ui] = pl.comment.clone();
        }
    }
}

// ============================================================================
// handle_orifices()
// ============================================================================

pub fn handle_orifices(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for pl in parse_section(lines) {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.len() < 3 {
            continue;
        }

        let name = &tok[0];
        let idx = add_unique(&mut ctx.link_names, name, &mut ctx.errors);
        if idx < 0 {
            continue; // duplicate ID (ERR 207, legacy input.c parity)
        }

        ensure_link_capacity(ctx, idx);

        let orr = ctx.link_subtypes.set_link_type(&mut ctx.links, idx, LinkType::Orifice);
        let uorr = orr as usize;
        let ui = idx as usize;
        set_link_nodes(ctx, idx, &tok[1], &tok[2]);
        // tok[3]: SIDE or BOTTOM -> orifice_type (0=BOTTOM, 1=SIDE)
        if tok.len() > 3 {
            ctx.link_subtypes.orifices.orifice_type[uorr] =
                if tokenizer::to_upper(&tok[3]) == b"SIDE" { 1.0 } else { 0.0 };
        }
        // tok[4]: offset (height above invert)
        if tok.len() > 4 {
            ctx.links.offset1[ui] = offset_token(ctx, &tok[4]);
        }
        // tok[5]: discharge coefficient
        if tok.len() > 5 {
            ctx.link_subtypes.orifices.cd[uorr] = to_double(&tok[5], 0.0);
        }
        // tok[6]: flap gate (YES/NO)
        if tok.len() > 6 {
            ctx.links.has_flap_gate[ui] = tokenizer::parse_boolean(&tok[6]) as u8;
        }
        // tok[7]: open/close time (seconds)
        if tok.len() > 7 {
            ctx.link_subtypes.orifices.orate[uorr] = to_double(&tok[7], 0.0);
        }
        if !pl.comment.is_empty() {
            ctx.links.comments[ui] = pl.comment.clone();
        }
    }
}

// ============================================================================
// handle_weirs()
// ============================================================================

pub fn handle_weirs(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for pl in parse_section(lines) {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.len() < 3 {
            continue;
        }

        let name = &tok[0];
        let idx = add_unique(&mut ctx.link_names, name, &mut ctx.errors);
        if idx < 0 {
            continue; // duplicate ID (ERR 207, legacy input.c parity)
        }

        ensure_link_capacity(ctx, idx);

        let wr = ctx.link_subtypes.set_link_type(&mut ctx.links, idx, LinkType::Weir);
        let uwr = wr as usize;
        let ui = idx as usize;
        set_link_nodes(ctx, idx, &tok[1], &tok[2]);
        // tok[3]: weir type (TRANSVERSE=0, SIDEFLOW=1, V-NOTCH=2, TRAPEZOIDAL=3, ROADWAY=4)
        let mut wt = 0.0;
        if tok.len() > 3 {
            let wtype = tokenizer::to_upper(&tok[3]);
            if wtype == b"SIDEFLOW" {
                wt = 1.0;
            } else if wtype == b"V-NOTCH" {
                wt = 2.0;
            } else if wtype == b"TRAPEZOIDAL" {
                wt = 3.0;
            } else if wtype == b"ROADWAY" {
                wt = 4.0;
            }
            ctx.link_subtypes.weirs.weir_type[uwr] = wt;
        }
        // tok[4]: crest height (above invert)
        if tok.len() > 4 {
            ctx.link_subtypes.weirs.crest_height[uwr] = offset_token(ctx, &tok[4]);
        }
        // tok[5]: discharge coefficient
        if tok.len() > 5 {
            ctx.link_subtypes.weirs.cd[uwr] = to_double(&tok[5], 0.0);
        }
        // tok[6]: flap gate (YES/NO)
        if tok.len() > 6 {
            ctx.links.has_flap_gate[ui] = tokenizer::parse_boolean(&tok[6]) as u8;
        }
        // tok[7]: end contractions
        if tok.len() > 7 {
            ctx.link_subtypes.weirs.end_contractions[uwr] = to_double(&tok[7], 0.0);
        }
        // tok[8]: end-section discharge coeff
        if tok.len() > 8 && tok[8] != b"*" {
            ctx.link_subtypes.weirs.cd2[uwr] = to_double(&tok[8], 0.0);
        }
        // tok[9]: can-surcharge flag (default YES)
        if tok.len() > 9 && tok[9] != b"*" {
            ctx.link_subtypes.weirs.can_surcharge[uwr] =
                if tokenizer::to_upper(&tok[9]) == b"YES" { 1u8 } else { 0u8 };
        }
        // tok[10] road width, tok[11] road surface: ROADWAY weir only
        if wt == 4.0 {
            if tok.len() > 10 {
                ctx.link_subtypes.weirs.road_width[uwr] = to_double(&tok[10], 0.0);
            }
            if tok.len() > 11 {
                let surf = tokenizer::to_upper(&tok[11]);
                ctx.link_subtypes.weirs.road_surface[uwr] = if surf == b"PAVED" {
                    1i8
                } else if surf == b"GRAVEL" {
                    2i8
                } else {
                    0i8
                };
            }
        }
        // tok[12]: discharge-coefficient curve, kept on the link's curve-name slot.
        if tok.len() > 12 && tok[12] != b"*" {
            if ui < ctx.links.pump_curve_name.len() {
                ctx.links.pump_curve_name[ui] = tok[12].clone();
            }
        }
        if !pl.comment.is_empty() {
            ctx.links.comments[ui] = pl.comment.clone();
        }
    }
}

// ============================================================================
// handle_outlets()
// ============================================================================

/// `std::string::find(needle) != npos`.
fn contains(hay: &[u8], needle: &[u8]) -> bool {
    hay.windows(needle.len()).any(|w| w == needle)
}

pub fn handle_outlets(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for pl in parse_section(lines) {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.len() < 3 {
            continue;
        }

        let name = &tok[0];
        let idx = add_unique(&mut ctx.link_names, name, &mut ctx.errors);
        if idx < 0 {
            continue; // duplicate ID (ERR 207, legacy input.c parity)
        }

        ensure_link_capacity(ctx, idx);

        let olr = ctx.link_subtypes.set_link_type(&mut ctx.links, idx, LinkType::Outlet);
        let uolr = olr as usize;
        let ui = idx as usize;
        set_link_nodes(ctx, idx, &tok[1], &tok[2]);
        if tok.len() > 3 {
            let mut crest = offset_token(ctx, &tok[3]);
            // legacy outlet_readParams zeroes a negative crest under LINK_OFFSETS DEPTH.
            if ctx.options.link_offsets == 0 && crest < 0.0 {
                crest = 0.0;
            }
            ctx.link_subtypes.outlets.crest_height[uolr] = crest;
        }
        // tok[4]: type string; outlet_type 0=FUNCTIONAL_HEAD, 1=FUNCTIONAL_DEPTH,
        // 2=TABULAR_HEAD, 3=TABULAR_DEPTH (case-sensitive substring tests)
        if tok.len() > 4 {
            let type_str = &tok[4];
            let is_tabular = contains(type_str, b"TABULAR");
            let is_depth_based = contains(type_str, b"DEPTH");
            ctx.link_subtypes.outlets.outlet_type[uolr] = if is_tabular {
                if is_depth_based {
                    3.0
                } else {
                    2.0
                }
            } else if is_depth_based {
                1.0
            } else {
                0.0
            };
        }
        let outlet_type = cvt_i32(ctx.link_subtypes.outlets.outlet_type[uolr]);
        let is_tabular = outlet_type >= 2;
        if tok.len() > 5 {
            if is_tabular {
                // Store curve name in pump_curve_name for PostParseResolver
                if ui < ctx.links.pump_curve_name.len() {
                    ctx.links.pump_curve_name[ui] = tok[5].clone();
                }
            } else {
                ctx.link_subtypes.outlets.coeff[uolr] = to_double(&tok[5], 0.0);
            }
        }
        if tok.len() > 6 && !is_tabular {
            ctx.link_subtypes.outlets.expon[uolr] = to_double(&tok[6], 0.0);
        }
        // Trailing optional Gated column: tok[6] for TABULAR, tok[7] for FUNCTIONAL.
        let gate_tok: usize = if is_tabular { 6 } else { 7 };
        if tok.len() > gate_tok {
            ctx.links.has_flap_gate[ui] = tokenizer::parse_boolean(&tok[gate_tok]) as u8;
        }
        if !pl.comment.is_empty() {
            ctx.links.comments[ui] = pl.comment.clone();
        }
    }
}

// ============================================================================
// handle_xsections()
// ============================================================================

/// Shape keywords in legacy's XsectTypeWords order; the first entry that is a PREFIX of
/// the token wins (legacy findmatch), so "DUMMY_2" reads as DUMMY.
const SHAPE_WORDS: [(&[u8], XsectShape); 26] = [
    (b"DUMMY", XsectShape::Dummy),
    (b"CIRCULAR", XsectShape::Circular),
    (b"FILLED_CIRCULAR", XsectShape::FilledCircular),
    (b"RECT_CLOSED", XsectShape::RectClosed),
    (b"RECT_OPEN", XsectShape::RectOpen),
    (b"TRAPEZOIDAL", XsectShape::Trapezoidal),
    (b"TRIANGULAR", XsectShape::Triangular),
    (b"PARABOLIC", XsectShape::Parabolic),
    (b"POWER", XsectShape::Power),
    (b"RECT_TRIANGULAR", XsectShape::RectTriang),
    (b"RECT_ROUND", XsectShape::RectRound),
    (b"MODBASKETHANDLE", XsectShape::Modbaskethandle),
    (b"HORIZ_ELLIPSE", XsectShape::HorizEllipse),
    (b"VERT_ELLIPSE", XsectShape::VertEllipse),
    (b"ARCH", XsectShape::Arch),
    (b"EGG", XsectShape::Eggshaped),
    (b"HORSESHOE", XsectShape::Horseshoe),
    (b"GOTHIC", XsectShape::Gothic),
    (b"CATENARY", XsectShape::Catenary),
    (b"SEMIELLIPTICAL", XsectShape::Semielliptical),
    (b"BASKETHANDLE", XsectShape::Baskethandle),
    (b"SEMICIRCULAR", XsectShape::Semicircular),
    (b"IRREGULAR", XsectShape::Irregular),
    (b"CUSTOM", XsectShape::Custom),
    (b"FORCE_MAIN", XsectShape::ForceMain),
    (b"STREET", XsectShape::StreetXsect),
];

/// First legacy cross-section keyword that prefixes the (upper-cased) token, or -1.
fn find_shape_word(s: &[u8]) -> i32 {
    for (i, (word, _shape)) in SHAPE_WORDS.iter().enumerate() {
        if s.starts_with(word) {
            return i as i32;
        }
    }
    -1
}

pub fn handle_xsections(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.is_empty() {
            continue;
        }
        // Legacy link.c:182 rejects a row with fewer than 3 items.
        if tok.len() < 3 {
            ctx.errors.push(format_error(ERR_ITEMS, &tok[0]));
            continue;
        }

        let idx = ctx.link_names.find(&tok[0]);
        if idx < 0 {
            // The link is declared further down the file: stash the row for InputReader
            // to re-dispatch once everything is parsed.
            ctx.deferred_section_rows.push((b"XSECTIONS".to_vec(), line.clone()));
            continue;
        }

        ensure_link_capacity(ctx, idx);
        let ui = idx as usize;

        let shape_str = tokenizer::to_upper(&tok[1]);
        let shape_word = find_shape_word(&shape_str);
        if shape_word < 0 {
            // Legacy link.c:190 rejects an unmatched shape keyword outright.
            ctx.errors.push(format_error(ERR_KEYWORD, &tok[1]));
            continue;
        }
        ctx.links.xsect_shape[ui] = SHAPE_WORDS[shape_word as usize].1;
        let shape = ctx.links.xsect_shape[ui];

        // Legacy link.c:221: every shape but IRREGULAR and STREET needs its four
        // geometry items.
        if shape != XsectShape::Irregular && shape != XsectShape::StreetXsect && tok.len() < 6 {
            ctx.errors.push(format_error(ERR_ITEMS, &tok[0]));
            continue;
        }
        // Legacy xsect_setParams refuses Geom1 <= 0 for every shape but DUMMY.
        if shape != XsectShape::Irregular
            && shape != XsectShape::StreetXsect
            && shape != XsectShape::Custom
            && shape != XsectShape::Dummy
            && !(to_double(&tok[2], 0.0) > 0.0)
        {
            let mut detail = b"for Geom1 of link ".to_vec();
            detail.extend_from_slice(&tok[0]);
            detail.extend_from_slice(b" in [XSECTIONS] (must be > 0)");
            ctx.errors.push(format_error_detail(ERR_NUMBER, &tok[2], &detail));
            continue;
        }

        // IRREGULAR / STREET: tok[2] is a transect / street name. CUSTOM: tok[2] = y_full,
        // tok[3] = shape curve name. All resolved later.
        if shape == XsectShape::Irregular || shape == XsectShape::StreetXsect {
            if tok.len() > 2 {
                ctx.links.pump_curve_name[ui] = tok[2].clone(); // Reuse field for transect/street name
                ctx.links.xsect_curve[ui] = -1;
            }
        } else if shape == XsectShape::Custom {
            if tok.len() > 2 {
                ctx.links.xsect_y_full[ui] = to_double(&tok[2], 0.0);
            }
            if tok.len() > 3 {
                ctx.links.pump_curve_name[ui] = tok[3].clone(); // Shape curve name
                ctx.links.xsect_curve[ui] = -1;
            }
        } else {
            // Geom1 = full depth (diameter for circular, etc.)
            if tok.len() > 2 {
                ctx.links.xsect_y_full[ui] = to_double(&tok[2], 0.0);
            }
        }

        // Geom2 = width or second parameter (skip for CUSTOM)
        if shape != XsectShape::Custom {
            if tok.len() > 3 {
                ctx.links.xsect_w_max[ui] = to_double(&tok[3], 0.0);
            }
        }

        // Geom3 = third shape parameter
        if tok.len() > 4 {
            ctx.links.xsect_y_bot[ui] = to_double(&tok[4], 0.0);
        }

        // Geom4 = fourth shape parameter
        if tok.len() > 5 {
            ctx.links.xsect_r_bot[ui] = to_double(&tok[5], 0.0);
        }

        // Retain raw Geom1–Geom4 (display units) for lossless serialization.
        if shape != XsectShape::Irregular && shape != XsectShape::StreetXsect {
            if tok.len() > 2 {
                ctx.links.xsect_geom1[ui] = to_double(&tok[2], 0.0);
            }
            if tok.len() > 3 {
                ctx.links.xsect_geom2[ui] = to_double(&tok[3], 0.0);
            }
            if tok.len() > 4 {
                ctx.links.xsect_geom3[ui] = to_double(&tok[4], 0.0);
            }
            if tok.len() > 5 {
                ctx.links.xsect_geom4[ui] = to_double(&tok[5], 0.0);
            }
        }

        // Barrels: dual-write the conduit side-table row if it exists.
        let cr = ctx.link_subtypes.conduit_row(idx);
        if tok.len() > 6 && cr >= 0 {
            let barrels = cvt_i32(to_double(&tok[6], 0.0));
            if barrels > 0 {
                ctx.link_subtypes.conduits.barrels[cr as usize] = barrels;
            }
        }

        // Culvert code (optional, token 7)
        if tok.len() > 7 && cr >= 0 {
            let cc = cvt_i32(to_double(&tok[7], 0.0));
            if cc > 0 {
                ctx.link_subtypes.conduits.culvert_code[cr as usize] = cc;
            }
        }
    }
}

// ============================================================================
// handle_losses()
// ============================================================================

pub fn handle_losses(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 2 {
            continue; // At minimum: Link Kentry
        }

        let idx = ctx.link_names.find(&tok[0]);
        if idx < 0 {
            // Same deferral as [XSECTIONS].
            ctx.deferred_section_rows.push((b"LOSSES".to_vec(), line.clone()));
            continue;
        }

        ensure_link_capacity(ctx, idx);

        let cr = ctx.link_subtypes.conduit_row(idx);
        let ucr = cr as usize;
        if cr >= 0 {
            if tok.len() > 1 {
                ctx.link_subtypes.conduits.loss_inlet[ucr] = to_double(&tok[1], 0.0);
            }
            if tok.len() > 2 {
                ctx.link_subtypes.conduits.loss_outlet[ucr] = to_double(&tok[2], 0.0);
            }
            if tok.len() > 3 {
                ctx.link_subtypes.conduits.loss_avg[ucr] = to_double(&tok[3], 0.0);
            }
            if tok.len() > 5 {
                ctx.link_subtypes.conduits.seep_rate[ucr] = to_double(&tok[5], 0.0);
            }
        }
        if tok.len() > 4 {
            ctx.links.has_flap_gate[idx as usize] = tokenizer::parse_boolean(&tok[4]) as u8;
        }
    }
}

// ============================================================================
// handle_transects()
// ============================================================================

pub fn handle_transects(ctx: &mut SimulationContext, lines: &[Bytes]) {
    // NC sets Manning's n for the next transect, X1 starts one, GR adds
    // station-elevation pairs to the current one.
    let mut nc_left = 0.0;
    let mut nc_right = 0.0;
    let mut nc_channel = 0.0;

    let mut pending_comment: Bytes = Vec::new();
    for pl in parse_section(lines) {
        if !pl.comment.is_empty() {
            pending_comment = pl.comment.clone();
        }
        let tok = tokenizer::tokenize(&pl.data);
        if tok.is_empty() {
            continue;
        }

        let keyword = tokenizer::to_upper(&tok[0]);

        if keyword == b"NC" {
            // NC nLeft nRight nChannel. A zero component means "unchanged from the
            // preceding NC record"; a negative one is ERR_NUMBER.
            let n_left = if tok.len() > 1 { to_double(&tok[1], 0.0) } else { 0.0 };
            let n_right = if tok.len() > 2 { to_double(&tok[2], 0.0) } else { 0.0 };
            let n_channel = if tok.len() > 3 { to_double(&tok[3], 0.0) } else { 0.0 };
            if n_left < 0.0 {
                ctx.errors.push(format_error(ERR_NUMBER, &tok[1]));
                continue;
            }
            if n_right < 0.0 {
                ctx.errors.push(format_error(ERR_NUMBER, &tok[2]));
                continue;
            }
            if n_channel < 0.0 {
                ctx.errors.push(format_error(ERR_NUMBER, &tok[3]));
                continue;
            }
            if n_left > 0.0 {
                nc_left = n_left;
            }
            if n_right > 0.0 {
                nc_right = n_right;
            }
            if n_channel > 0.0 {
                nc_channel = n_channel;
            }
            // overbank roughness defaults to the channel value, as in legacy
            if nc_left == 0.0 {
                nc_left = nc_channel;
            }
            if nc_right == 0.0 {
                nc_right = nc_channel;
            }
        } else if keyword == b"X1" {
            // X1 Name Nsta Xleft Xright 0 0 Lfactor Xfactor Yfactor; Lfactor/Xfactor of
            // zero mean 1.0.
            if tok.len() < 3 {
                continue;
            }

            let name = &tok[1];

            // A second X1 with the same name is a duplicate ID (ERR 207); keep parsing.
            for existing in &ctx.transects.names {
                if ieq(existing, name) {
                    ctx.errors.push(format_error(ERR_DUP_NAME, name));
                    break;
                }
            }

            ctx.transects.names.push(name.clone());
            ctx.transects.comments.push(pending_comment.clone());
            pending_comment.clear();
            ctx.transects.n_left.push(nc_left);
            ctx.transects.n_right.push(nc_right);
            ctx.transects.n_channel.push(nc_channel);
            ctx.transects.x_left_bank.push(if tok.len() > 3 { to_double(&tok[3], 0.0) } else { 0.0 });
            ctx.transects.x_right_bank.push(if tok.len() > 4 { to_double(&tok[4], 0.0) } else { 0.0 });
            ctx.transects.x_left_encroachment.push(if tok.len() > 5 { to_double(&tok[5], 0.0) } else { 0.0 });
            ctx.transects.x_right_encroachment.push(if tok.len() > 6 { to_double(&tok[6], 0.0) } else { 0.0 });

            let mut l_factor = if tok.len() > 7 { to_double(&tok[7], 0.0) } else { 1.0 };
            let mut x_factor = if tok.len() > 8 { to_double(&tok[8], 0.0) } else { 1.0 };
            let y_factor = if tok.len() > 9 { to_double(&tok[9], 0.0) } else { 0.0 };
            if l_factor == 0.0 {
                l_factor = 1.0;
            }
            if x_factor == 0.0 {
                x_factor = 1.0;
            }
            ctx.transects.length_factor.push(l_factor);
            ctx.transects.x_factor.push(x_factor);
            ctx.transects.y_factor.push(y_factor);

            ctx.transects.stations.push(Vec::new());
            ctx.transects.elevations.push(Vec::new());
        } else if keyword == b"GR" {
            // GR elev station elev station ... (pairs)
            if ctx.transects.names.is_empty() {
                continue; // No active transect
            }

            // (back() on an empty vector is undefined upstream; a panic here.)
            let tr = &mut ctx.transects;
            let sta = tr.stations.last_mut().expect("stations in step with names");
            let elev = tr.elevations.last_mut().expect("elevations in step with names");
            // Pairs start at tok[1]: elev station elev station ...
            let mut i = 1usize;
            while i + 1 < tok.len() {
                elev.push(to_double(&tok[i], 0.0));
                sta.push(to_double(&tok[i + 1], 0.0));
                i += 2;
            }
        }
    }
}
