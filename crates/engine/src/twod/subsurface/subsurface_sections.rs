// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/subsurface/SubsurfaceSections.hpp and
// SubsurfaceSections.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The [2D_AQUIFER*] and [2D_SURFACE_OWNERSHIP] parsers and their registration (what open()
// runs). The resolution functions (gwUnitFactors, resolveSubsurface, mapMeshPolyline,
// resolveLinkSeepage, resolveSubcatchInfiltration, resolveSurfaceOwnership) run from
// SurfaceRouter2D::initialize, and writeSubsurfaceSections from the .inp writer: phase 2b,
// not translated here.
//   * num() is std::from_chars consuming the whole token and finite; inum() narrows with
//     static_cast<int> (core::cmath::cvt_i32).
//   * The handlers reach the router's SubsurfaceConfig and name lists through shared handles.

#![allow(clippy::all)]

use crate::core::charconv::from_chars_double;
use crate::core::cmath::cvt_i32;
use crate::core::context::SimulationContext;
use crate::core::options::Bytes;
use crate::input::section_registry::SectionRegistry;
use crate::input::tokenizer;
use crate::twod::subsurface::subsurface_data::*;
use crate::twod::subsurface::surface_ownership::parse_surface_owner_line;
use crate::twod::Shared;

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut v = Bytes::new();
    for p in parts {
        v.extend_from_slice(p);
    }
    v
}

fn iequals(a: &[u8], b: &[u8]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_ascii_uppercase() == y.to_ascii_uppercase())
}

fn num(s: &[u8], v: &mut f64) -> bool {
    match from_chars_double(s) {
        (Ok(x), used) => {
            *v = x;
            used == s.len() && x.is_finite()
        }
        (Err(_), _) => false,
    }
}

fn inum(s: &[u8], v: &mut i32) -> bool {
    let mut d = 0.0;
    if !num(s, &mut d) {
        return false;
    }
    *v = cvt_i32(d);
    *v as f64 == d
}

/// boolTok(): YES/TRUE/1 or NO/FALSE/0.
fn bool_tok(s: &[u8], v: &mut bool) -> bool {
    if iequals(s, b"YES") || iequals(s, b"TRUE") || s == b"1" {
        *v = true;
        return true;
    }
    if iequals(s, b"NO") || iequals(s, b"FALSE") || s == b"0" {
        *v = false;
        return true;
    }
    false
}

/// parseAquiferOptionsLine() (`tokens` is never empty: the handler skips blank rows)
pub fn parse_aquifer_options_line(tokens: &[Bytes], opts: &mut GwOptions) -> Bytes {
    if tokens.len() < 2 {
        return cat(&[b"expected KEY VALUE: ", &tokens[0]]);
    }
    let k = &tokens[0];
    let v = &tokens[1];

    if iequals(k, b"GW_CLOCK") {
        if !iequals(v, b"SHARED") && !iequals(v, b"INDEPENDENT") {
            return cat(&[b"GW_CLOCK must be SHARED or INDEPENDENT: ", v]);
        }
        opts.independent_clock = iequals(v, b"INDEPENDENT");
        return Bytes::new();
    }
    if iequals(k, b"GW_BASE_STEP") || iequals(k, b"GW_MAX_STEP") {
        let mut seconds = 0.0;
        if !num(v, &mut seconds) || seconds <= 0.0 {
            return cat(&[k, b" must be a finite positive duration in seconds: ", v]);
        }
        if iequals(k, b"GW_BASE_STEP") {
            opts.base_step = seconds;
        } else {
            opts.max_step = seconds;
        }
        return Bytes::new();
    }
    if iequals(k, b"GW_BINDS_SURFACE_DT") {
        if !bool_tok(v, &mut opts.binds_surface_dt) {
            return cat(&[b"GW_BINDS_SURFACE_DT must be YES or NO: ", v]);
        }
        return Bytes::new();
    }
    if iequals(k, b"GW_BUDGET_FRACTION") || iequals(k, b"GW_HEAD_TOLERANCE") {
        let mut value = 0.0;
        let fraction = iequals(k, b"GW_BUDGET_FRACTION");
        if !num(v, &mut value) || value <= 0.0 || value > 1.0 || (!fraction && value == 1.0) {
            return cat(&[k, if fraction { &b" must be in (0,1]: "[..] } else { &b" must be in (0,1): "[..] }, v]);
        }
        if fraction {
            opts.budget_fraction = value;
        } else {
            opts.head_tolerance = value;
        }
        return Bytes::new();
    }
    if iequals(k, b"SOIL_CHAR") {
        if !parse_soil_char(v, &mut opts.soil_char) {
            return cat(&[b"unknown SOIL_CHAR: ", v]);
        }
        return Bytes::new();
    }
    if iequals(k, b"CLOSURE") {
        if !parse_gw_closure(v, &mut opts.closure) {
            return cat(&[b"unknown CLOSURE: ", v]);
        }
        return Bytes::new();
    }
    if iequals(k, b"M_LAYERS") {
        let mut m = 0;
        if !inum(v, &mut m) || m < 2 || m > 128 {
            return cat(&[b"M_LAYERS must be 2..128: ", v]);
        }
        opts.m_layers = m;
        return Bytes::new();
    }
    if iequals(k, b"CAPILLARY_DIFF") {
        if !bool_tok(v, &mut opts.capillary_diff) {
            return cat(&[b"CAPILLARY_DIFF must be YES/NO: ", v]);
        }
        return Bytes::new();
    }
    if iequals(k, b"FORCE_CLOSED_FORM") {
        if !bool_tok(v, &mut opts.force_closed_form) {
            return cat(&[b"FORCE_CLOSED_FORM must be YES/NO: ", v]);
        }
        return Bytes::new();
    }
    if iequals(k, b"DUNNE") {
        if !bool_tok(v, &mut opts.dunne) {
            return cat(&[b"DUNNE must be YES/NO: ", v]);
        }
        return Bytes::new();
    }
    if iequals(k, b"C_GW") || iequals(k, b"C_COL") {
        let mut d = 0.0;
        if !num(v, &mut d) || d <= 0.0 || d > 1.0 {
            return cat(&[if iequals(k, b"C_GW") { &b"C_GW"[..] } else { &b"C_COL"[..] }, b" must be in (0, 1]: ", v]);
        }
        if iequals(k, b"C_GW") {
            opts.c_gw = d;
        } else {
            opts.c_col = d;
        }
        return Bytes::new();
    }
    if iequals(k, b"MODE") {
        if iequals(v, b"MESH") {
            opts.per_subcatch = false;
            return Bytes::new();
        }
        if iequals(v, b"PER_SUBCATCH") {
            opts.per_subcatch = true;
            return Bytes::new();
        }
        return cat(&[b"MODE must be MESH or PER_SUBCATCH: ", v]);
    }
    if iequals(k, b"NODE_ENROLMENT") {
        if iequals(v, b"AUTO") {
            opts.node_auto = true;
            return Bytes::new();
        }
        if iequals(v, b"ROWS") {
            opts.node_auto = false;
            return Bytes::new();
        }
        return cat(&[b"NODE_ENROLMENT must be AUTO or ROWS: ", v]);
    }
    if iequals(k, b"LINK_SEEPAGE") {
        if iequals(v, b"DEFAULT") {
            opts.link_seepage = GwLinkMode::Default;
            return Bytes::new();
        }
        if iequals(v, b"AUTO") || iequals(v, b"ONE_WAY") {
            opts.link_seepage = GwLinkMode::Auto;
            return Bytes::new();
        }
        if iequals(v, b"NONE") {
            opts.link_seepage = GwLinkMode::None;
            return Bytes::new();
        }
        if iequals(v, b"TWO_WAY") {
            opts.link_seepage = GwLinkMode::TwoWay;
            return Bytes::new();
        }
        return cat(&[b"LINK_SEEPAGE must be DEFAULT, ONE_WAY (legacy AUTO), NONE or TWO_WAY: ", v]);
    }
    if iequals(k, b"GW_ET") {
        if iequals(v, b"AUTO") || iequals(v, b"NONE") || iequals(v, b"CAPILLARY_RISE") || iequals(v, b"BOUNDARY_ET") || iequals(v, b"BOTH") {
            opts.gw_et = v.iter().map(|c| c.to_ascii_uppercase()).collect();
            return Bytes::new();
        }
        return cat(&[b"GW_ET must be NONE, CAPILLARY_RISE, BOUNDARY_ET or BOTH: ", v]);
    }
    if iequals(k, b"WILTING_SUCTION") {
        if iequals(v, b"AUTO") {
            opts.wilting_suction = 150.0;
            opts.wilting_suction_set = false;
            return Bytes::new();
        }
        let mut suction = 0.0;
        if !num(v, &mut suction) || !suction.is_finite() || !(suction > 0.0) {
            return b"WILTING_SUCTION must be AUTO or a positive finite project length".to_vec();
        }
        opts.wilting_suction = suction;
        opts.wilting_suction_set = true;
        return Bytes::new();
    }
    cat(&[b"unknown option: ", k])
}

/// parseAquiferLine() (`tokens` is never empty: the handler skips blank rows)
pub fn parse_aquifer_line(tokens: &[Bytes], rows: &mut Vec<GwAquiferRow>) -> Bytes {
    let mut r = GwAquiferRow::default();
    let mut at;

    if tokens[0] == b"*" {
        r.scope = 0;
        at = 1;
    } else if iequals(&tokens[0], b"TAG") {
        if tokens.len() < 2 {
            return b"TAG needs a name".to_vec();
        }
        r.scope = 1;
        r.tag = tokens[1].clone();
        at = 2;
    } else if iequals(&tokens[0], b"CELL") {
        let mut n = 0;
        if tokens.len() < 2 || !inum(&tokens[1], &mut n) || n < 1 {
            return b"CELL needs a 1-based index".to_vec();
        }
        r.scope = 2;
        r.cell = n - 1;
        at = 2;
    } else {
        return cat(&[b"row must start with CELL, TAG or '*': ", &tokens[0]]);
    }

    // Five positional columns: KS ZS THETA_S THETA_R ALPHA.
    const K_POS: [&[u8]; 5] = [b"KS", b"ZS", b"THETA_S", b"THETA_R", b"ALPHA"];
    for i in 0..5 {
        if at >= tokens.len() {
            return cat(&[b"missing ", K_POS[i]]);
        }
        let dst: &mut f64 = match i {
            0 => &mut r.Ks,
            1 => &mut r.zs,
            2 => &mut r.theta_s,
            3 => &mut r.theta_r,
            _ => &mut r.alpha,
        };
        if !num(&tokens[at], dst) {
            return cat(&[b"invalid ", K_POS[i], b": ", &tokens[at]]);
        }
        at += 1;
    }

    // Everything law-specific is a keyword.
    while at < tokens.len() {
        let k = &tokens[at];
        if at + 1 >= tokens.len() {
            return cat(&[b"keyword ", k, b" needs a value"]);
        }
        let v = &tokens[at + 1];
        at += 2;
        if iequals(k, b"PSI_B") {
            if !num(v, &mut r.psi_b) {
                return cat(&[b"invalid PSI_B: ", v]);
            }
        } else if iequals(k, b"LAMBDA") {
            if !num(v, &mut r.lambda) {
                return cat(&[b"invalid LAMBDA: ", v]);
            }
        } else if iequals(k, b"N") {
            if !num(v, &mut r.vg_n) {
                return cat(&[b"invalid N: ", v]);
            }
        } else if iequals(k, b"L") {
            if !num(v, &mut r.vg_L) {
                return cat(&[b"invalid L: ", v]);
            }
        } else if iequals(k, b"C_LOSS") {
            if !num(v, &mut r.c_loss) {
                return cat(&[b"invalid C_LOSS: ", v]);
            }
        } else if iequals(k, b"HG0") {
            if !num(v, &mut r.hg0) {
                return cat(&[b"invalid HG0: ", v]);
            }
        } else if iequals(k, b"SOIL_CHAR") {
            if !parse_soil_char(v, &mut r.soil_char) {
                return cat(&[b"unknown SOIL_CHAR: ", v]);
            }
            r.soil_char_set = true;
        } else if iequals(k, b"CLOSURE") {
            if !parse_gw_closure(v, &mut r.closure) {
                return cat(&[b"unknown CLOSURE: ", v]);
            }
            r.closure_set = true;
        } else if iequals(k, b"M_LAYERS") {
            if !inum(v, &mut r.m_layers) || r.m_layers < 2 || r.m_layers > 128 {
                return cat(&[b"M_LAYERS must be 2..128: ", v]);
            }
        } else {
            return cat(&[b"unknown keyword: ", k]);
        }
    }

    if r.Ks <= 0.0 {
        return b"KS must be > 0".to_vec();
    }
    if r.zs <= 0.0 {
        return b"ZS must be > 0".to_vec();
    }
    if r.theta_s <= 0.0 || r.theta_s > 1.0 {
        return b"THETA_S must be in (0, 1]".to_vec();
    }
    if r.theta_r < 0.0 || r.theta_r >= r.theta_s {
        return b"THETA_R must be in [0, THETA_S)".to_vec();
    }
    if r.alpha <= 0.0 {
        return b"ALPHA must be > 0".to_vec();
    }
    if r.vg_n <= 1.0 {
        return b"N must be > 1".to_vec();
    }
    if r.lambda <= 0.0 {
        return b"LAMBDA must be > 0".to_vec();
    }
    if r.psi_b < 0.0 {
        return b"PSI_B must be >= 0".to_vec();
    }
    if r.c_loss < 0.0 {
        return b"C_LOSS must be >= 0".to_vec();
    }

    rows.push(r);
    Bytes::new()
}

/// parseAquiferNodeLine()
pub fn parse_aquifer_node_line(tokens: &[Bytes], beds: &mut Vec<GwNodeBed>, names: &mut Vec<Bytes>) -> Bytes {
    if tokens.len() < 2 {
        return b"expected NODE CELL|AUTO [KC k DC d] [AREA a] [EXCHANGE YES|NO]".to_vec();
    }
    let mut b = GwNodeBed::default();
    let mut cell = 0;
    if iequals(&tokens[1], b"AUTO") {
        b.locate = true; // G-X2: located at resolve
    } else {
        if !inum(&tokens[1], &mut cell) || cell < 1 {
            return cat(&[b"CELL must be a 1-based index or AUTO: ", &tokens[1]]);
        }
        b.cell = cell - 1;
    }

    let mut at = 2usize;
    while at < tokens.len() {
        let k = &tokens[at];
        if at + 1 >= tokens.len() {
            return cat(&[b"keyword ", k, b" needs a value"]);
        }
        let v = &tokens[at + 1];
        at += 2;
        if iequals(k, b"KC") {
            if !num(v, &mut b.Kc) {
                return cat(&[b"invalid KC: ", v]);
            }
        } else if iequals(k, b"DC") {
            if !num(v, &mut b.dC) {
                return cat(&[b"invalid DC: ", v]);
            }
        } else if iequals(k, b"AREA") {
            if !num(v, &mut b.area) {
                return cat(&[b"invalid AREA: ", v]);
            }
        } else if iequals(k, b"EXCHANGE") {
            if !bool_tok(v, &mut b.exchange) {
                return cat(&[b"EXCHANGE must be YES/NO: ", v]);
            }
        } else {
            return cat(&[b"unknown keyword: ", k]);
        }
    }
    if b.Kc < 0.0 || b.dC < 0.0 || b.area < 0.0 {
        return b"KC, DC and AREA must be >= 0".to_vec();
    }
    if b.Kc > 0.0 && b.dC <= 0.0 {
        return b"KC given without a positive DC (a bed needs a thickness)".to_vec();
    }

    beds.push(b);
    names.push(tokens[0].clone());
    Bytes::new()
}

/// parseAquiferLinkLine(): `LINK [KC k] [DC d] [EXCHANGE YES|NO]`.
pub fn parse_aquifer_link_line(tokens: &[Bytes], rows: &mut Vec<GwLinkRow>, names: &mut Vec<Bytes>) -> Bytes {
    if tokens.is_empty() {
        return b"expected LINK [KC k] [DC d] [EXCHANGE YES|NO]".to_vec();
    }
    let mut r = GwLinkRow::default();
    let mut at = 1usize;
    while at < tokens.len() {
        let k = &tokens[at];
        if at + 1 >= tokens.len() {
            return cat(&[b"keyword ", k, b" needs a value"]);
        }
        let v = &tokens[at + 1];
        at += 2;
        if iequals(k, b"KC") {
            if !num(v, &mut r.Kc) {
                return cat(&[b"invalid KC: ", v]);
            }
        } else if iequals(k, b"DC") {
            if !num(v, &mut r.dC) {
                return cat(&[b"invalid DC: ", v]);
            }
        } else if iequals(k, b"EXCHANGE") {
            if !bool_tok(v, &mut r.exchange) {
                return cat(&[b"EXCHANGE must be YES/NO: ", v]);
            }
        } else {
            return cat(&[b"unknown keyword: ", k]);
        }
    }
    if r.Kc < 0.0 || r.dC < 0.0 {
        return b"KC and DC must be >= 0".to_vec();
    }
    rows.push(r);
    names.push(tokens[0].clone());
    Bytes::new()
}

/// parseAquiferBoundaryLine(): `CELL n EDGE e HEAD|FLOW CONST value|SERIES name`.
pub fn parse_aquifer_boundary_line(t: &[Bytes], rows: &mut Vec<GwBoundaryRow>) -> Bytes {
    if t.len() != 7 || !iequals(&t[0], b"CELL") || !iequals(&t[2], b"EDGE") {
        return b"expected CELL n EDGE e HEAD|FLOW CONST value|SERIES name".to_vec();
    }
    let mut row = GwBoundaryRow::default();
    let mut cell = 0;
    if !inum(&t[1], &mut cell) || cell < 1 || !inum(&t[3], &mut row.edge) || row.edge < 0 || row.edge > 3 {
        return b"CELL must be positive and EDGE must be in [0,3]".to_vec();
    }
    row.cell = cell - 1;
    if iequals(&t[4], b"HEAD") {
        row.head = true;
    } else if !iequals(&t[4], b"FLOW") {
        return b"boundary kind must be HEAD or FLOW".to_vec();
    }
    if iequals(&t[5], b"CONST") {
        if !num(&t[6], &mut row.value) {
            return b"boundary constant must be finite".to_vec();
        }
    } else if iequals(&t[5], b"SERIES") {
        if t[6].is_empty() {
            return b"boundary SERIES requires a name".to_vec();
        }
        row.series = t[6].clone();
    } else {
        return b"boundary forcing must be CONST or SERIES".to_vec();
    }
    for old in rows.iter() {
        if old.cell == row.cell && old.edge == row.edge {
            return b"duplicate hydraulic boundary edge".to_vec();
        }
    }
    rows.push(row);
    Bytes::new()
}

/// makeHandler(): marks the options authored, then parses each non-blank line, stopping at the
/// first diagnostic (error code 5, "[<label>] <err> - line: <raw>").
fn make_handler(
    label: &'static [u8],
    cfg: Shared<SubsurfaceConfig>,
    names: Option<Shared<Vec<Bytes>>>,
    lp: fn(&mut SubsurfaceConfig, Option<&mut Vec<Bytes>>, &[Bytes]) -> Bytes,
) -> impl Fn(&mut SimulationContext, &[Bytes]) + 'static {
    move |ctx: &mut SimulationContext, lines: &[Bytes]| {
        cfg.borrow_mut().options.authored = true;
        for raw in lines {
            let tokens = tokenizer::tokenize(raw);
            if tokens.is_empty() {
                continue;
            }
            let err = {
                let mut c = cfg.borrow_mut();
                match &names {
                    Some(n) => lp(&mut c, Some(&mut n.borrow_mut()), &tokens),
                    None => lp(&mut c, None, &tokens),
                }
            };
            if !err.is_empty() {
                ctx.error_code = 5; // SWMM_ERR_PARSE
                ctx.error_message = cat(&[b"[", label, b"] ", &err, b" - line: ", raw]);
                return;
            }
        }
    }
}

/// registerSubsurfaceSections()
pub fn register_subsurface_sections(cfg: &Shared<SubsurfaceConfig>, node_names: &Shared<Vec<Bytes>>, link_names: &Shared<Vec<Bytes>>, registry: &mut SectionRegistry) {
    {
        let cfg = cfg.clone();
        registry.register_custom(b"2D_SURFACE_OWNERSHIP", move |ctx: &mut SimulationContext, lines: &[Bytes]| {
            for raw in lines {
                let t = tokenizer::tokenize(raw);
                if t.is_empty() {
                    continue;
                }
                let error = parse_surface_owner_line(&t, &mut cfg.borrow_mut().surface_owners);
                if !error.is_empty() {
                    ctx.error_code = 5;
                    ctx.error_message = cat(&[b"[2D_SURFACE_OWNERSHIP] ", &error]);
                    return;
                }
            }
        });
    }
    registry.register_custom(b"2D_AQUIFER_BOUNDARIES", make_handler(b"2D_AQUIFER_BOUNDARIES", cfg.clone(), None, |c, _, t| parse_aquifer_boundary_line(t, &mut c.boundaries)));
    registry.register_custom(b"2D_AQUIFER_OPTIONS", make_handler(b"2D_AQUIFER_OPTIONS", cfg.clone(), None, |c, _, t| parse_aquifer_options_line(t, &mut c.options)));
    registry.register_custom(b"2D_AQUIFER", make_handler(b"2D_AQUIFER", cfg.clone(), None, |c, _, t| parse_aquifer_line(t, &mut c.rows)));
    registry.register_custom(
        b"2D_AQUIFER_NODE",
        make_handler(b"2D_AQUIFER_NODE", cfg.clone(), Some(node_names.clone()), |c, n, t| match n {
            Some(n) => parse_aquifer_node_line(t, &mut c.node_beds, n),
            None => Bytes::new(),
        }),
    );
    registry.register_custom(
        b"2D_AQUIFER_LINKS",
        make_handler(b"2D_AQUIFER_LINKS", cfg.clone(), Some(link_names.clone()), |c, n, t| match n {
            Some(n) => parse_aquifer_link_line(t, &mut c.link_rows, n),
            None => Bytes::new(),
        }),
    );
}
