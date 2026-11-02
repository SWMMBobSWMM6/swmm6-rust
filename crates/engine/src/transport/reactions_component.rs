// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/components/ReactionModule/
// ReactionsComponent.hpp and ReactionsComponent.cpp (HydroCouple/openswmm.engine, branch
// swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache
// License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The `org.hydrocouple.openswmm.reactions` process component (phase R1/R2): parse and
// validate model.rxn [REACTION_*] sections (or the embedded copies kept from the .inp)
// into ctx.reactions, compile every expression into the flat token pool, and commit the
// species to the registry. The save hook (serializeReactionSystem, ReactionsWriter.cpp)
// is not part of open() and is not translated here.
//
// Differences of form only:
//   * The parse helpers take the context and borrow ctx.reactions where upstream holds
//     `auto& rx = ctx.reactions` beside other ctx reads.
//   * parseExpressions picks the form/source vectors by scope at each access instead of
//     binding references to them up front.
//   * recompileReactionSystem's compile_one lambda is a function that names the
//     destination span by slot; the symbol tables are moved out of the reaction data for
//     the compile call and moved back (RxSymbols borrows them; the pool is written).

use crate::core::context::{ProcessComponentSpec, SimulationContext};
use crate::core::cstdlib::{c_str, strtod};
use crate::core::options::{Bytes, QualitySolverKind};
use crate::core::process_components::{ComponentConfigSections, ProcessComponentRegistry};
use crate::data::node::NodeType;
use crate::data::reaction::{
    ReactionAreaUnits, ReactionCoupling, ReactionData, ReactionDataLidLaw, ReactionDataParameterOverride, ReactionDataSource, ReactionDataSourcePattern,
    ReactionDataSubsurfaceLaw, ReactionDataSurfaceLaw, ReactionDataWallArea, ReactionExprForm, ReactionRateUnits, ReactionSolverKind,
};
use crate::data::reaction_tokens::RxExprSpan;
use crate::data::species_registry::SpeciesKind;
use crate::input::tokenizer::{to_upper, tokenize};
use crate::transport::reaction_expression::{compile_reaction_expression, reaction_context_diagnostic, RxSymbols};

const K_COMPONENT_ID: &str = "org.hydrocouple.openswmm.reactions";

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

/// The raw remainder of `line` after its first `n_lead` fields (space/tab/comma delimiters,
/// a field opening with '"' runs to its closing quote); trailing space/tab/CR removed and
/// one enclosing quote pair stripped.
fn expr_payload(line: &[u8], n_lead: usize) -> Bytes {
    let is_delim = |c: u8| c == b' ' || c == b'\t' || c == b',';
    let mut i = 0usize;
    for _ in 0..n_lead {
        while i < line.len() && is_delim(line[i]) {
            i += 1;
        }
        if i < line.len() && line[i] == b'"' {
            // quoted field
            i += 1;
            while i < line.len() && line[i] != b'"' {
                i += 1;
            }
            if i < line.len() {
                i += 1;
            }
        } else {
            while i < line.len() && !is_delim(line[i]) {
                i += 1;
            }
        }
    }
    while i < line.len() && is_delim(line[i]) {
        i += 1;
    }
    let mut e = line.len();
    while e > i && (line[e - 1] == b' ' || line[e - 1] == b'\t' || line[e - 1] == b'\r') {
        e -= 1;
    }
    let mut s = line[i..e].to_vec();
    if s.len() >= 2 && s[0] == b'"' && *s.last().unwrap() == b'"' {
        s = s[1..s.len() - 1].to_vec();
    }
    s
}

/// strtod over the C string with full consumption, and a finite result.
fn to_num(s: &[u8], out: &mut f64) -> bool {
    let (v, n) = strtod(s);
    *out = v;
    n == c_str(s).len() && n != 0 && out.is_finite()
}

fn form_of(u: &[u8]) -> ReactionExprForm {
    if u == b"RATE" {
        ReactionExprForm::Rate
    } else if u == b"EQUIL" {
        ReactionExprForm::Equil
    } else if u == b"FORMULA" {
        ReactionExprForm::Formula
    } else {
        ReactionExprForm::None
    }
}

fn parse_options(ctx: &mut SimulationContext, lines: &[Bytes], errors: &mut Vec<Bytes>) {
    for line in lines {
        let tok = tokenize(line);
        if tok.len() < 2 {
            continue;
        }
        let key = to_upper(&tok[0]);
        let val = to_upper(&tok[1]);
        let mut num = 0.0;
        let rx = &mut ctx.reactions;
        if key == b"SOLVER" {
            if val == b"EUL" {
                rx.solver = ReactionSolverKind::Eul;
            } else if val == b"RK5" {
                rx.solver = ReactionSolverKind::Rk5;
            } else if val == b"ROS2" {
                rx.solver = ReactionSolverKind::Ros2;
            } else if val == b"BDF2" {
                rx.solver = ReactionSolverKind::Bdf2;
            } else {
                errors.push(cat(&[b"[REACTION_OPTIONS] unknown SOLVER '", &tok[1], b"' (EUL|RK5|ROS2|BDF2)."]));
            }
        } else if key == b"COUPLING" {
            if val == b"NONE" {
                rx.coupling = ReactionCoupling::None;
            } else if val == b"FULL" {
                rx.coupling = ReactionCoupling::Full;
            } else {
                errors.push(cat(&[b"[REACTION_OPTIONS] unknown COUPLING '", &tok[1], b"' (NONE|FULL)."]));
            }
        } else if key == b"RATE_UNITS" {
            if val == b"SEC" {
                rx.rate_units = ReactionRateUnits::Sec;
            } else if val == b"MIN" {
                rx.rate_units = ReactionRateUnits::Min;
            } else if val == b"HR" {
                rx.rate_units = ReactionRateUnits::Hr;
            } else if val == b"DAY" {
                rx.rate_units = ReactionRateUnits::Day;
            } else {
                errors.push(cat(&[b"[REACTION_OPTIONS] unknown RATE_UNITS '", &tok[1], b"' (SEC|MIN|HR|DAY)."]));
            }
        } else if key == b"WALL_GEOMETRY" {
            if val == b"YES" {
                rx.wall_geometry = true;
            } else if val == b"NO" {
                rx.wall_geometry = false;
            } else {
                errors.push(b"[REACTION_OPTIONS] WALL_GEOMETRY must be YES or NO.".to_vec());
            }
        } else if key == b"AREA_UNITS" {
            if val == b"FT2" {
                rx.area_units = ReactionAreaUnits::Ft2;
            } else if val == b"M2" {
                rx.area_units = ReactionAreaUnits::M2;
            } else if val == b"CM2" {
                rx.area_units = ReactionAreaUnits::Cm2;
            } else {
                errors.push(cat(&[b"[REACTION_OPTIONS] unknown AREA_UNITS '", &tok[1], b"' (FT2|M2|CM2)."]));
            }
        } else if key == b"TIMESTEP" {
            if to_num(&tok[1], &mut num) && num >= 0.0 {
                rx.timestep = num;
            } else {
                errors.push(cat(&[b"[REACTION_OPTIONS] bad TIMESTEP '", &tok[1], b"'."]));
            }
            // Physics review A3: kept for round trip, read by no binding.
            if ctx.reactions.timestep > 0.0 {
                ctx.warnings.push(cat(&[
                    b"[REACTION_OPTIONS] TIMESTEP ",
                    &tok[1],
                    b" is not used: reactions advance with the quality step and choose their own substeps within it.",
                ]));
            }
        } else if key == b"PATTERN_TIMESTEP" {
            if to_num(&tok[1], &mut num) && num > 0.0 {
                rx.pattern_timestep = num;
            } else {
                errors.push(cat(&[b"[REACTION_OPTIONS] bad PATTERN_TIMESTEP '", &tok[1], b"'."]));
            }
        } else if key == b"ATOL" {
            if to_num(&tok[1], &mut num) && num > 0.0 {
                rx.atol = num;
            } else {
                errors.push(cat(&[b"[REACTION_OPTIONS] bad ATOL '", &tok[1], b"'."]));
            }
        } else if key == b"RTOL" {
            if to_num(&tok[1], &mut num) && num > 0.0 {
                rx.rtol = num;
            } else {
                errors.push(cat(&[b"[REACTION_OPTIONS] bad RTOL '", &tok[1], b"'."]));
            }
        } else if key == b"TEMPERATURE" {
            if to_num(&tok[1], &mut num) {
                rx.default_temp_c = num;
            } else {
                errors.push(cat(&[b"[REACTION_OPTIONS] bad TEMPERATURE '", &tok[1], b"'."]));
            }
        } else {
            errors.push(cat(&[b"[REACTION_OPTIONS] unknown option '", &tok[0], b"'."]));
        }
    }
}

fn parse_species(ctx: &mut SimulationContext, lines: &[Bytes], errors: &mut Vec<Bytes>) {
    for line in lines {
        let tok = tokenize(line);
        if tok.is_empty() {
            continue;
        }
        if tok.len() < 3 {
            errors.push(cat(&[b"[REACTION_SPECIES] row needs BULK|WALL name units [atol rtol]: '", line, b"'."]));
            continue;
        }
        let kind = to_upper(&tok[0]);
        let is_wall = kind == b"WALL";
        if !is_wall && kind != b"BULK" {
            errors.push(cat(&[b"[REACTION_SPECIES] kind must be BULK or WALL: '", &tok[0], b"'."]));
            continue;
        }
        let name = &tok[1];
        if ctx.reactions.find_species(name) >= 0 {
            errors.push(cat(&[b"[REACTION_SPECIES] duplicate species '", name, b"'."]));
            continue;
        }
        // Collision with pollutant (or reserved) names is refused. STAGED ONLY: the registry
        // commit happens after the whole config validates and compiles.
        if ctx.species_registry.find(name) >= 0 {
            errors.push(cat(&[b"[REACTION_SPECIES] '", name, b"' collides with an existing pollutant or species name."]));
            continue;
        }
        let rx = &mut ctx.reactions;
        rx.species_name.push(name.clone());
        rx.species_is_wall.push(if is_wall { 1 } else { 0 });
        rx.species_units.push(tok[2].clone());
        let mut atol = 0.0;
        let mut rtol = 0.0;
        if tok.len() > 3 && !to_num(&tok[3], &mut atol) {
            errors.push(cat(&[b"[REACTION_SPECIES] bad atol for '", name, b"'."]));
        }
        if tok.len() > 4 && !to_num(&tok[4], &mut rtol) {
            errors.push(cat(&[b"[REACTION_SPECIES] bad rtol for '", name, b"'."]));
        }
        rx.species_atol.push(atol);
        rx.species_rtol.push(rtol);
    }
}

fn parse_coefficients(ctx: &mut SimulationContext, lines: &[Bytes], errors: &mut Vec<Bytes>) {
    let rx = &mut ctx.reactions;
    for line in lines {
        let tok = tokenize(line);
        if tok.is_empty() {
            continue;
        }
        if tok.len() < 3 {
            errors.push(cat(&[b"[REACTION_COEFFICIENTS] row needs PARAMETER|CONSTANT name value: '", line, b"'."]));
            continue;
        }
        let kind = to_upper(&tok[0]);
        if kind != b"PARAMETER" && kind != b"CONSTANT" {
            errors.push(cat(&[b"[REACTION_COEFFICIENTS] kind must be PARAMETER or CONSTANT: '", &tok[0], b"'."]));
            continue;
        }
        if rx.find_coef(&tok[1]) >= 0 || rx.find_species(&tok[1]) >= 0 {
            errors.push(cat(&[b"[REACTION_COEFFICIENTS] name '", &tok[1], b"' duplicates a coefficient or species."]));
            continue;
        }
        let mut v = 0.0;
        if !to_num(&tok[2], &mut v) {
            errors.push(cat(&[b"[REACTION_COEFFICIENTS] bad value for '", &tok[1], b"'."]));
            continue;
        }
        rx.coef_name.push(tok[1].clone());
        rx.coef_is_param.push(if kind == b"PARAMETER" { 1 } else { 0 });
        rx.coef_value.push(v);
    }
}

fn parse_terms(ctx: &mut SimulationContext, lines: &[Bytes], errors: &mut Vec<Bytes>) {
    let rx = &mut ctx.reactions;
    for line in lines {
        let tok = tokenize(line);
        if tok.is_empty() {
            continue;
        }
        if tok.len() < 2 {
            errors.push(cat(&[b"[REACTION_TERMS] row needs name expression: '", line, b"'."]));
            continue;
        }
        if rx.find_term(&tok[0]) >= 0 || rx.find_species(&tok[0]) >= 0 || rx.find_coef(&tok[0]) >= 0 {
            errors.push(cat(&[b"[REACTION_TERMS] name '", &tok[0], b"' duplicates a term, species, or coefficient."]));
            continue;
        }
        rx.term_name.push(tok[0].clone());
        rx.term_expr_src.push(expr_payload(line, 1));
    }
}

fn parse_subsurface(ctx: &mut SimulationContext, lines: &[Bytes], errors: &mut Vec<Bytes>) {
    let rx = &mut ctx.reactions;
    for line in lines {
        let tokens = tokenize(line);
        if tokens.is_empty() {
            continue;
        }
        let fail = |errors: &mut Vec<Bytes>, why: &[u8]| {
            errors.push(cat(&[b"[REACTION_SUBSURFACE] ", why, " \u{2014} ".as_bytes(), line]));
        };
        let mut at: usize = 0;
        let mut law = ReactionDataSubsurfaceLaw::default();
        let number = |at: &mut usize, value: &mut i32| -> bool {
            if *at >= tokens.len() {
                return false;
            }
            let mut n = 0.0;
            let ok = to_num(&tokens[*at], &mut n);
            *at += 1;
            if !ok || !n.is_finite() || n < 1.0 || n > 2147483647.0 || n != n.floor() {
                return false;
            }
            *value = n as i32 - 1;
            true
        };
        if tokens[at] == b"*" {
            at += 1;
        } else if to_upper(&tokens[at]) == b"TAG" {
            law.scope = 1;
            at += 1;
            if at >= tokens.len() {
                fail(errors, b"TAG needs a name");
                continue;
            }
            law.tag = tokens[at].clone();
            at += 1;
        } else if to_upper(&tokens[at]) == b"CELL" {
            law.scope = 2;
            at += 1;
            if !number(&mut at, &mut law.cell) {
                fail(errors, b"CELL needs a positive index");
                continue;
            }
        } else {
            fail(errors, b"scope must be *, TAG or CELL");
            continue;
        }
        if at >= tokens.len() {
            fail(errors, b"missing SAT, UNSAT or LAYER");
            continue;
        }
        let zone = to_upper(&tokens[at]);
        at += 1;
        if zone == b"SAT" {
            law.zone = 0;
        } else if zone == b"UNSAT" {
            law.zone = 1;
        } else if zone == b"LAYER" {
            law.zone = 2;
            if !number(&mut at, &mut law.layer) {
                fail(errors, b"LAYER needs a positive index");
                continue;
            }
        } else {
            fail(errors, b"zone must be SAT, UNSAT or LAYER");
            continue;
        }
        if tokens.len() < at + 3 {
            fail(errors, b"needs SPECIES RATE|EQUIL|FORMULA EXPRESSION");
            continue;
        }
        law.species = rx.find_species(&tokens[at]);
        at += 1;
        if law.species < 0 {
            fail(errors, b"species must be a declared MSX BULK species");
            continue;
        }
        if rx.species_is_wall[law.species as usize] != 0 {
            fail(errors, b"groundwater WALL species are unsupported");
            continue;
        }
        let form = to_upper(&tokens[at]);
        at += 1;
        let f = form_of(&form);
        if f == ReactionExprForm::None {
            fail(errors, b"unknown reaction form");
            continue;
        }
        law.form = f;
        law.expression = expr_payload(line, at);
        rx.subsurface_laws.push(law);
    }
}

fn forms_of(rx: &mut ReactionData, scope: i32) -> &mut Vec<ReactionExprForm> {
    if scope == 2 {
        &mut rx.subcatch_form
    } else if scope == 1 {
        &mut rx.tank_form
    } else {
        &mut rx.pipe_form
    }
}

fn srcs_of(rx: &mut ReactionData, scope: i32) -> &mut Vec<Bytes> {
    if scope == 2 {
        &mut rx.subcatch_expr_src
    } else if scope == 1 {
        &mut rx.tank_expr_src
    } else {
        &mut rx.pipe_expr_src
    }
}

/// scope: 0 pipes, 1 tanks, 2 subcatchments.
fn parse_expressions(ctx: &mut SimulationContext, lines: &[Bytes], scope: i32, errors: &mut Vec<Bytes>) {
    let sec: &[u8] = if scope == 2 {
        b"[REACTION_SUBCATCHMENTS]"
    } else if scope == 1 {
        b"[REACTION_TANKS]"
    } else {
        b"[REACTION_PIPES]"
    };
    let ns = ctx.reactions.n_species() as usize;
    forms_of(&mut ctx.reactions, scope).resize(ns, ReactionExprForm::None);
    srcs_of(&mut ctx.reactions, scope).resize(ns, Vec::new());
    for line in lines {
        let tok = tokenize(line);
        if tok.is_empty() {
            continue;
        }
        if tok.len() < 3 {
            errors.push(cat(&[sec, b" row needs RATE|EQUIL|FORMULA species expression: '", line, b"'."]));
            continue;
        }
        let form = to_upper(&tok[0]);
        let f = form_of(&form);
        if f == ReactionExprForm::None {
            errors.push(cat(&[sec, b" unknown form '", &tok[0], b"' (RATE|EQUIL|FORMULA)."]));
            continue;
        }
        let s = ctx.reactions.find_species(&tok[1]);
        if s < 0 {
            if ctx.species_registry.find(&tok[1]) >= 0 {
                errors.push(cat(&[
                    sec,
                    b" '",
                    &tok[1],
                    "' is a pollutant \u{2014} pollutant kinetics (RATE/EQUIL/FORMULA on a pollutant) arrive with plan phase R4b; pollutants may be REFERENCED read-only in MSX expressions today, and kdecay applies under every engine.".as_bytes(),
                ]));
            } else {
                errors.push(cat(&[sec, b" undeclared species '", &tok[1], "' \u{2014} declare it in [REACTION_SPECIES].".as_bytes()]));
            }
            continue;
        }
        let us = s as usize;
        if forms_of(&mut ctx.reactions, scope)[us] != ReactionExprForm::None {
            errors.push(cat(&[sec, b" species '", &tok[1], b"' already has an expression in this scope."]));
            continue;
        }
        forms_of(&mut ctx.reactions, scope)[us] = f;
        srcs_of(&mut ctx.reactions, scope)[us] = expr_payload(line, 2);
    }
}

fn parse_quality(ctx: &mut SimulationContext, lines: &[Bytes], errors: &mut Vec<Bytes>) {
    {
        let rx = &mut ctx.reactions;
        let ns = rx.n_species() as usize;
        rx.init_global.clear();
        rx.init_global.resize(ns, 0.0);
        rx.init_elem_is_link.clear();
        rx.init_elem_idx.clear();
        rx.init_elem_species.clear();
        rx.init_elem_value.clear();
    }
    for line in lines {
        let tok = tokenize(line);
        if tok.is_empty() {
            continue;
        }
        let scope = to_upper(&tok[0]);

        // E-B1: NODE/LINK per-element rows, resolved here (names are final at open).
        if scope == b"NODE" || scope == b"LINK" {
            let link = scope == b"LINK";
            if tok.len() < 4 {
                errors.push(cat(&[b"[REACTION_QUALITY] row needs ", &scope, b" element species value: '", line, b"'."]));
                continue;
            }
            let ei = if link { ctx.link_names.find(&tok[1]) } else { ctx.node_names.find(&tok[1]) };
            if ei < 0 {
                errors.push(cat(&[b"[REACTION_QUALITY] unknown ", if link { b"link" } else { b"node" }, b" '", &tok[1], b"'."]));
                continue;
            }
            let rx = &mut ctx.reactions;
            let s = rx.find_species(&tok[2]);
            if s < 0 {
                errors.push(cat(&[b"[REACTION_QUALITY] undeclared species '", &tok[2], b"'."]));
                continue;
            }
            let mut v = 0.0;
            if !to_num(&tok[3], &mut v) || v < 0.0 {
                errors.push(cat(&[b"[REACTION_QUALITY] bad value for '", &tok[2], b"'."]));
                continue;
            }
            let mut dup = false;
            for k in 0..rx.init_elem_idx.len() {
                if (rx.init_elem_is_link[k] != 0) == link && rx.init_elem_idx[k] == ei && rx.init_elem_species[k] == s {
                    dup = true;
                    break;
                }
            }
            if dup {
                errors.push(cat(&[
                    b"[REACTION_QUALITY] duplicate row for '",
                    &tok[2],
                    b"' at ",
                    if link { b"link" } else { b"node" },
                    b" '",
                    &tok[1],
                    b"'.",
                ]));
                continue;
            }
            rx.init_elem_is_link.push(if link { 1 } else { 0 });
            rx.init_elem_idx.push(ei);
            rx.init_elem_species.push(s);
            rx.init_elem_value.push(v);
            continue;
        }

        if scope != b"GLOBAL" {
            errors.push(cat(&[b"[REACTION_QUALITY] scope '", &tok[0], b"' is not GLOBAL, NODE, or LINK."]));
            continue;
        }
        if tok.len() < 3 {
            errors.push(cat(&[b"[REACTION_QUALITY] row needs GLOBAL species value: '", line, b"'."]));
            continue;
        }
        let rx = &mut ctx.reactions;
        let s = rx.find_species(&tok[1]);
        if s < 0 {
            errors.push(cat(&[b"[REACTION_QUALITY] undeclared species '", &tok[1], b"'."]));
            continue;
        }
        let mut v = 0.0;
        if !to_num(&tok[2], &mut v) || v < 0.0 {
            errors.push(cat(&[b"[REACTION_QUALITY] bad value for '", &tok[1], b"'."]));
            continue;
        }
        rx.init_global[s as usize] = v;
    }
    // U2: [INITIAL_QUALITY] MSX rows live in ctx.initial_quality; bring them back.
    for e in crate::transport::msx_initial_quality::mirror_initial_quality_msx_rows(ctx) {
        errors.push(e);
    }
}

fn parse_patterns(ctx: &mut SimulationContext, lines: &[Bytes], errors: &mut Vec<Bytes>) {
    let rx = &mut ctx.reactions;
    for line in lines {
        let tok = tokenize(line);
        if tok.is_empty() {
            continue;
        }
        if tok.len() < 2 {
            errors.push(cat(&[b"[REACTION_PATTERNS] row needs name factors: '", line, b"'."]));
            continue;
        }
        let mut index: i32 = -1;
        for i in 0..rx.source_patterns.len() as i32 {
            if rx.source_patterns[i as usize].name == tok[0] {
                index = i;
            }
        }
        if index < 0 {
            index = rx.source_patterns.len() as i32;
            rx.source_patterns.push(ReactionDataSourcePattern { name: tok[0].clone(), factors: Vec::new() });
        }
        for i in 1..tok.len() {
            let mut value = 0.0;
            if !to_num(&tok[i], &mut value) || value < 0.0 {
                errors.push(cat(&[b"[REACTION_PATTERNS] factor must be finite and nonnegative: '", &tok[i], b"'."]));
            } else {
                rx.source_patterns[index as usize].factors.push(value);
            }
        }
    }
}

fn parse_sources(ctx: &mut SimulationContext, lines: &[Bytes], errors: &mut Vec<Bytes>) {
    for line in lines {
        let tok = tokenize(line);
        if tok.is_empty() {
            continue;
        }
        let kind = to_upper(&tok[0]);
        if (tok.len() != 4 && tok.len() != 5) || (kind != b"CONC" && kind != b"MASS" && kind != b"FLOWPACED" && kind != b"SETPOINT") {
            errors.push(cat(&[b"[REACTION_SOURCES] row needs CONC|MASS|FLOWPACED|SETPOINT node species value [pattern]: '", line, b"'."]));
            continue;
        }
        let node = ctx.node_names.find(&tok[1]);
        let species = ctx.reactions.find_species(&tok[2]);
        let mut value = 0.0;
        if node < 0 || species < 0 || !to_num(&tok[3], &mut value) || !value.is_finite() || value < 0.0 {
            errors.push(cat(&[b"[REACTION_SOURCES] unknown node/species or invalid value: '", line, b"'."]));
            continue;
        }
        let rx = &mut ctx.reactions;
        if rx.species_is_wall[species as usize] != 0 {
            errors.push(b"[REACTION_SOURCES] WALL species cannot enter with inflow water.".to_vec());
            continue;
        }
        let mut pattern: i32 = -1;
        if tok.len() == 5 && tok[4] != b"*" {
            for i in 0..rx.source_patterns.len() as i32 {
                if rx.source_patterns[i as usize].name == tok[4] {
                    pattern = i;
                }
            }
            if pattern < 0 {
                errors.push(cat(&[b"[REACTION_SOURCES] unknown pattern '", &tok[4], b"'."]));
                continue;
            }
        }
        let mut duplicate = false;
        for row in &rx.sources {
            duplicate |= row.node == node && row.species == species;
        }
        if duplicate {
            errors.push(cat(&[b"[REACTION_SOURCES] duplicate node/species: '", line, b"'."]));
            continue;
        }
        if (kind == b"FLOWPACED" || kind == b"SETPOINT") && node < ctx.nodes.is_virtual.len() as i32 && ctx.nodes.is_virtual[node as usize] != 0 {
            errors.push(b"[REACTION_SOURCES] boosters require a physical mixing node, not a virtual splice.".to_vec());
            continue;
        }
        let booster = if kind == b"FLOWPACED" {
            1
        } else if kind == b"SETPOINT" {
            2
        } else {
            0
        };
        ctx.reactions.sources.push(ReactionDataSource { mass: kind == b"MASS", node, species, pattern, value, booster });
    }
}

fn parse_parameters(ctx: &mut SimulationContext, lines: &[Bytes], errors: &mut Vec<Bytes>) {
    for line in lines {
        let tok = tokenize(line);
        if tok.is_empty() {
            continue;
        }
        let scope = to_upper(&tok[0]);
        if tok.len() != 4 || (scope != b"NODE" && scope != b"LINK") {
            errors.push(cat(&[b"[REACTION_PARAMETERS] row needs NODE|LINK element parameter value: '", line, b"'."]));
            continue;
        }
        let link = scope == b"LINK";
        let element = if link { ctx.link_names.find(&tok[1]) } else { ctx.node_names.find(&tok[1]) };
        let rx = &mut ctx.reactions;
        let coefficient = rx.find_coef(&tok[2]);
        let mut value = 0.0;
        if element < 0 || coefficient < 0 || !to_num(&tok[3], &mut value) {
            errors.push(cat(&[b"[REACTION_PARAMETERS] unknown element/parameter or non-finite value: '", line, b"'."]));
            continue;
        }
        if rx.coef_is_param[coefficient as usize] == 0 {
            errors.push(cat(&[b"[REACTION_PARAMETERS] '", &tok[2], b"' is a CONSTANT, not a PARAMETER."]));
            continue;
        }
        let mut duplicate = false;
        for row in &rx.parameter_overrides {
            duplicate |= row.link == link && row.element == element && row.coefficient == coefficient;
        }
        if duplicate {
            errors.push(cat(&[b"[REACTION_PARAMETERS] duplicate element/parameter: '", line, b"'."]));
            continue;
        }
        rx.parameter_overrides.push(ReactionDataParameterOverride { link, element, coefficient, value });
    }
}

fn parse_report(ctx: &mut SimulationContext, lines: &[Bytes], errors: &mut Vec<Bytes>) {
    let rx = &mut ctx.reactions;
    let ns = rx.n_species() as usize;
    rx.report_enabled.clear();
    rx.report_enabled.resize(ns, 1);
    let mut seen: Vec<bool> = vec![false; rx.report_enabled.len()];
    let mut all_seen = false;
    for line in lines {
        let tok = tokenize(line);
        if tok.is_empty() {
            continue;
        }
        if tok.len() != 3 || to_upper(&tok[0]) != b"SPECIES" || (to_upper(&tok[2]) != b"YES" && to_upper(&tok[2]) != b"NO") {
            errors.push(cat(&[b"[REACTION_REPORT] row needs SPECIES ALL|species YES|NO: '", line, b"'."]));
            continue;
        }
        let value = to_upper(&tok[2]) == b"YES";
        if to_upper(&tok[1]) == b"ALL" {
            if all_seen {
                errors.push(b"[REACTION_REPORT] duplicate ALL row.".to_vec());
            }
            // ALL is a baseline; require it before named overrides.
            let mut named = false;
            for &v in &seen {
                named |= v;
            }
            if named {
                errors.push(b"[REACTION_REPORT] ALL must precede named species rows.".to_vec());
            }
            for x in rx.report_enabled.iter_mut() {
                *x = value as u8;
            }
            all_seen = true;
        } else {
            let m = rx.find_species(&tok[1]);
            if m < 0 {
                errors.push(cat(&[b"[REACTION_REPORT] undeclared species '", &tok[1], b"'."]));
            } else if seen[m as usize] {
                errors.push(cat(&[b"[REACTION_REPORT] duplicate species '", &tok[1], b"'."]));
            } else {
                rx.report_enabled[m as usize] = value as u8;
                seen[m as usize] = true;
            }
        }
    }
}

fn parse_rainfall(ctx: &mut SimulationContext, lines: &[Bytes], errors: &mut Vec<Bytes>) {
    let rx = &mut ctx.reactions;
    let mut seen: Vec<bool> = vec![false; rx.n_species() as usize];
    for line in lines {
        let tok = tokenize(line);
        if tok.is_empty() {
            continue;
        }
        if tok.len() < 2 {
            errors.push(cat(&[b"[REACTION_RAINFALL] row needs species concentration: '", line, b"'."]));
            continue;
        }
        let m = rx.find_species(&tok[0]);
        if m < 0 {
            errors.push(cat(&[b"[REACTION_RAINFALL] undeclared species '", &tok[0], b"'."]));
            continue;
        }
        let mut v = 0.0;
        if !to_num(&tok[1], &mut v) || v < 0.0 {
            errors.push(cat(&[b"[REACTION_RAINFALL] bad concentration '", &tok[1], b"' for '", &tok[0], b"'."]));
            continue;
        }
        if seen[m as usize] {
            errors.push(cat(&[b"[REACTION_RAINFALL] duplicate row for '", &tok[0], b"'."]));
            continue;
        }
        seen[m as usize] = true;
        rx.rain_conc[m as usize] = v;
    }
}

/// reactionSectionTags(): the recognised [REACTION_*] tags (upper-case, no brackets).
pub const REACTION_SECTION_TAGS: [&str; 17] = [
    "REACTION_OPTIONS",
    "REACTION_SPECIES",
    "REACTION_COEFFICIENTS",
    "REACTION_TERMS",
    "REACTION_PIPES",
    "REACTION_TANKS",
    "REACTION_SOURCES",
    "REACTION_QUALITY",
    "REACTION_PARAMETERS",
    "REACTION_PATTERNS",
    "REACTION_REPORT",
    "REACTION_SUBCATCHMENTS",
    "REACTION_RAINFALL",
    "REACTION_SUBSURFACE",
    "REACTION_WALLS",
    "REACTION_LIDS",
    "REACTION_BUILDUP",
];

/// Parse, validate, compile and register one native multispecies reaction system.
/// Any error in `errors` at the final check (including ones the caller brought in)
/// clears the reaction data and registers nothing.
pub fn apply_reaction_sections(ctx: &mut SimulationContext, config: &ComponentConfigSections, errors: &mut Vec<Bytes>) {
    if ctx.reactions.configured {
        errors.push("Reactions system configured twice \u{2014} duplicate registration or embedded sections alongside an external component file.".as_bytes().to_vec());
        return;
    }
    ctx.reactions.clear();

    // Order matters: species before coefficients/terms/expressions so references validate.
    if let Some(s) = config.find(b"REACTION_OPTIONS") {
        parse_options(ctx, s, errors);
    }
    if let Some(s) = config.find(b"REACTION_SPECIES") {
        parse_species(ctx, s, errors);
    }
    if ctx.reactions.n_species() == 0 {
        errors.push("Reactions config declares no [REACTION_SPECIES] \u{2014} at least one species is required.".as_bytes().to_vec());
        return;
    }
    if let Some(s) = config.find(b"REACTION_COEFFICIENTS") {
        parse_coefficients(ctx, s, errors);
    }
    if let Some(s) = config.find(b"REACTION_TERMS") {
        parse_terms(ctx, s, errors);
    }
    let ns = ctx.reactions.n_species() as usize;
    if let Some(s) = config.find(b"REACTION_PIPES") {
        parse_expressions(ctx, s, 0, errors);
    } else {
        ctx.reactions.pipe_form.clear();
        ctx.reactions.pipe_form.resize(ns, ReactionExprForm::None);
        ctx.reactions.pipe_expr_src.resize(ns, Vec::new());
    }
    if let Some(s) = config.find(b"REACTION_TANKS") {
        parse_expressions(ctx, s, 1, errors);
    } else {
        ctx.reactions.tank_form.clear();
        ctx.reactions.tank_form.resize(ns, ReactionExprForm::None);
        ctx.reactions.tank_expr_src.resize(ns, Vec::new());
    }
    ctx.reactions.subcatch_form.clear();
    ctx.reactions.subcatch_form.resize(ns, ReactionExprForm::None);
    ctx.reactions.subcatch_expr_src.resize(ns, Vec::new());
    if let Some(s) = config.find(b"REACTION_SUBCATCHMENTS") {
        let mut global: Vec<Bytes> = Vec::new();
        for line in s {
            let t = tokenize(line);
            if t.is_empty() {
                continue;
            }
            if to_upper(&t[0]) != b"SUBCATCHMENT" {
                global.push(line.clone());
                continue;
            }
            if t.len() < 5 {
                errors.push(b"[REACTION_SUBCATCHMENTS] needs SUBCATCHMENT name form species expression.".to_vec());
                continue;
            }
            let mut law = ReactionDataSurfaceLaw::default();
            law.subcatchment = ctx.subcatch_names.find(&t[1]);
            law.species = ctx.reactions.find_species(&t[3]);
            law.form = form_of(&to_upper(&t[2]));
            if law.subcatchment < 0 || law.species < 0 || law.form == ReactionExprForm::None {
                errors.push(cat(&[b"[REACTION_SUBCATCHMENTS] invalid selector, species or form: ", line]));
                continue;
            }
            if ctx.reactions.wet_laws.iter().any(|l| l.subcatchment == law.subcatchment && l.species == law.species) {
                errors.push(b"[REACTION_SUBCATCHMENTS] duplicate selected law.".to_vec());
                continue;
            }
            law.expression = expr_payload(line, 4);
            ctx.reactions.wet_laws.push(law);
        }
        parse_expressions(ctx, &global, 2, errors);
    }
    if let Some(s) = config.find(b"REACTION_BUILDUP") {
        for line in s {
            let t = tokenize(line);
            if t.is_empty() {
                continue;
            }
            let mut law = ReactionDataSurfaceLaw::default();
            let mut pos: usize = 0;
            let fail = |errors: &mut Vec<Bytes>| {
                errors.push(cat(&[b"[REACTION_BUILDUP] needs [SUBCATCHMENT name] LANDUSE name RATE|EQUIL|FORMULA BULK-species expression: ", line]));
            };
            if to_upper(&t[0]) == b"SUBCATCHMENT" {
                if t.len() < 7 {
                    fail(errors);
                    continue;
                }
                law.subcatchment = ctx.subcatch_names.find(&t[1]);
                if law.subcatchment < 0 {
                    fail(errors);
                    continue;
                }
                pos = 2;
            }
            if t.len() < pos + 5 || to_upper(&t[pos]) != b"LANDUSE" {
                fail(errors);
                continue;
            }
            law.landuse = ctx.landuse_names.find(&t[pos + 1]);
            if law.landuse < 0 {
                fail(errors);
                continue;
            }
            law.form = form_of(&to_upper(&t[pos + 2]));
            law.species = ctx.reactions.find_species(&t[pos + 3]);
            if law.species < 0 || law.form == ReactionExprForm::None {
                fail(errors);
                continue;
            }
            if ctx.reactions.buildup_laws.iter().any(|l| l.subcatchment == law.subcatchment && l.landuse == law.landuse && l.species == law.species) {
                errors.push(b"[REACTION_BUILDUP] duplicate selected law.".to_vec());
                continue;
            }
            law.expression = expr_payload(line, pos + 4);
            ctx.reactions.buildup_laws.push(law);
        }
    }
    if let Some(s) = config.find(b"REACTION_SUBSURFACE") {
        parse_subsurface(ctx, s, errors);
    }
    if let Some(s) = config.find(b"REACTION_QUALITY") {
        parse_quality(ctx, s, errors);
    } else {
        ctx.reactions.init_global.clear();
        ctx.reactions.init_global.resize(ns, 0.0);
    }
    ctx.reactions.rain_conc.clear();
    ctx.reactions.rain_conc.resize(ns, 0.0);
    if let Some(s) = config.find(b"REACTION_RAINFALL") {
        parse_rainfall(ctx, s, errors);
    }

    if let Some(s) = config.find(b"REACTION_PATTERNS") {
        parse_patterns(ctx, s, errors);
    }
    if let Some(s) = config.find(b"REACTION_SOURCES") {
        parse_sources(ctx, s, errors);
    }
    if let Some(s) = config.find(b"REACTION_PARAMETERS") {
        parse_parameters(ctx, s, errors);
    }
    if let Some(s) = config.find(b"REACTION_REPORT") {
        parse_report(ctx, s, errors);
    }
    if let Some(s) = config.find(b"REACTION_WALLS") {
        for line in s {
            let tok = tokenize(line);
            if tok.is_empty() {
                continue;
            }
            let mut area = 0.0;
            let link = to_upper(&tok[0]) == b"LINK";
            let e = if tok.len() == 3 {
                if link {
                    ctx.link_names.find(&tok[1])
                } else {
                    ctx.node_names.find(&tok[1])
                }
            } else {
                -1
            };
            if tok.len() != 3 || (!link && to_upper(&tok[0]) != b"NODE") || e < 0 || !to_num(&tok[2], &mut area) || !area.is_finite() || area < 0.0 {
                errors.push(cat(&[b"[REACTION_WALLS] needs NODE|LINK element nonnegative-area: '", line, b"'."]));
                continue;
            }
            let rows = &mut ctx.reactions.wall_areas;
            if rows.iter().any(|r| r.link == link && r.element == e) {
                errors.push(cat(&[b"[REACTION_WALLS] duplicate contact area: '", line, b"'."]));
            } else {
                rows.push(ReactionDataWallArea { link, element: e, area });
            }
        }
    }

    if let Some(s) = config.find(b"REACTION_LIDS") {
        for line in s {
            let t = tokenize(line);
            if t.is_empty() {
                continue;
            }
            let fail = |errors: &mut Vec<Bytes>, why: &[u8]| {
                errors.push(cat(&[b"[REACTION_LIDS] ", why, b": '", line, b"'."]));
            };
            if t.len() < 3 || to_upper(&t[0]) != b"CONTROL" {
                fail(errors, b"needs CONTROL name TRANSPORT or CONTROL name LAYER number RATE|EQUIL|FORMULA species expression");
                continue;
            }
            let control = ctx.lid_names.find(&t[1]);
            if control < 0 {
                fail(errors, b"unknown LID control");
                continue;
            }
            if !ctx.reactions.lid_enabled(control) {
                ctx.reactions.lid_controls.push(control);
            }
            if t.len() == 3 && to_upper(&t[2]) == b"TRANSPORT" {
                continue;
            }
            let mut number = 0.0;
            if t.len() < 7 || to_upper(&t[2]) != b"LAYER" || !to_num(&t[3], &mut number) || number < 0.0 || number > 10000.0 || number != number.floor() {
                fail(errors, b"invalid layer (0 means every layer)");
                continue;
            }
            let mut law = ReactionDataLidLaw::default();
            law.control = control;
            law.layer = number as i32;
            law.species = ctx.reactions.find_species(&t[5]);
            let f = to_upper(&t[4]);
            if f == b"RATE" {
                law.form = ReactionExprForm::Rate;
            } else if f == b"EQUIL" {
                law.form = ReactionExprForm::Equil;
            } else if f == b"FORMULA" {
                law.form = ReactionExprForm::Formula;
            } else {
                fail(errors, b"unknown form");
                continue;
            }
            if law.species < 0 || ctx.reactions.species_is_wall[law.species as usize] != 0 {
                fail(errors, b"target must be BULK");
                continue;
            }
            if ctx.reactions.lid_laws.iter().any(|r| r.control == control && r.layer == law.layer && r.species == law.species) {
                fail(errors, b"duplicate layer/species law");
                continue;
            }
            law.expression = expr_payload(line, 6);
            ctx.reactions.lid_laws.push(law);
        }
    }
    if !ctx.reactions.wall_geometry && ctx.reactions.buildup_laws.is_empty() && ctx.reactions.area_units != ReactionAreaUnits::Ft2 {
        let u: &[u8] = if ctx.reactions.area_units == ReactionAreaUnits::M2 { b"M2" } else { b"CM2" };
        ctx.warnings.push(cat(&[b"[REACTION_OPTIONS] AREA_UNITS ", u, b" has no effect unless WALL_GEOMETRY YES is selected."]));
    }

    // Unknown REACTION_* sections in the file are typos worth failing on.
    for sect in &config.sections {
        let mut known = false;
        for t in REACTION_SECTION_TAGS {
            if sect.0 == t.as_bytes() {
                known = true;
                break;
            }
        }
        if !known {
            errors.push(cat(&[b"Unknown section [", &sect.0, b"] in reactions config '", &config.source_path, b"'."]));
        }
    }

    if ctx.options.quality_solver == QualitySolverKind::Lagrangian && ctx.reactions.species_is_wall.iter().any(|&wall| wall != 0) {
        errors.push(b"LAGRANGIAN reaction transport does not support attached WALL species; use QUALITY_SOLVER LEGACY.".to_vec());
    }

    if !errors.is_empty() {
        // Transactional: a rejected config leaves NOTHING behind.
        ctx.reactions.clear();
        return;
    }

    // ---- R2: compile every expression into the flat token pool.
    if !recompile_reaction_system(ctx, errors) {
        ctx.reactions.clear();
        return;
    }

    // ---- Commit: registry entries only now, after full success.
    for sidx in 0..ctx.reactions.n_species() {
        let us = sidx as usize;
        let kind = if ctx.reactions.species_is_wall[us] != 0 { SpeciesKind::MsxWall } else { SpeciesKind::MsxBulk };
        let reg = ctx.species_registry.add(ctx.reactions.species_name[us].clone(), kind, ctx.reactions.species_units[us].clone());
        if ctx.reactions.registry_base < 0 {
            ctx.reactions.registry_base = reg;
        }
    }
    ctx.reactions.compiled = true;
    ctx.reactions.configured = true;
}

/// Embedded fallback (called from SWMMEngine::open after component resolution).
pub fn apply_embedded_reaction_sections(ctx: &mut SimulationContext, external_component_registered: bool, errors: &mut Vec<Bytes>) {
    if ctx.embedded_component_sections.is_empty() {
        return;
    }

    if external_component_registered {
        ctx.warnings.push(
            "Embedded [REACTION_*] sections in the .inp are IGNORED because an external reactions component is registered in [PROCESS_COMPONENTS] \u{2014} the external config file wins (TRANSPORT_IO_PLUGIN_CONFIG_PLAN \u{a7}3.2)."
                .as_bytes()
                .to_vec(),
        );
        return;
    }

    ctx.warnings.push(
        "Embedded [REACTION_*] sections found in the .inp. This works, but the clean layout is an external component config file registered via [PROCESS_COMPONENTS] (config=\"model.rxn\") \u{2014} see TRANSPORT_IO_PLUGIN_CONFIG_PLAN.md (style warning, D-UT8)."
            .as_bytes()
            .to_vec(),
    );

    let synth = ComponentConfigSections { source_path: b"(embedded in .inp)".to_vec(), sections: ctx.embedded_component_sections.clone() };
    apply_reaction_sections(ctx, &synth, errors);
}

fn apply_hook(ctx: &mut SimulationContext, _spec: &ProcessComponentSpec, config: &ComponentConfigSections, errors: &mut Vec<Bytes>) {
    apply_reaction_sections(ctx, config, errors);
}

/// Register the reactions component. (Its save hook, serializeReactionSystem, is not
/// part of open() and is not translated here.)
pub fn register_reactions_component(reg: &mut ProcessComponentRegistry) {
    reg.register_component(K_COMPONENT_ID.as_bytes(), b"Multispecies reaction system (EPANET-MSX conventions)", apply_hook);
}

/// The span compile_one writes.
#[derive(Clone, Copy)]
enum SpanSlot {
    Term(usize),
    Pipe(usize),
    Tank(usize),
    Subcatch(usize),
    Subsurface(usize),
    Wet(usize),
    Buildup(usize),
    Lid(usize),
}

fn slot_mut(rx: &mut ReactionData, slot: SpanSlot) -> &mut RxExprSpan {
    match slot {
        SpanSlot::Term(i) => &mut rx.term_expr[i],
        SpanSlot::Pipe(i) => &mut rx.pipe_expr[i],
        SpanSlot::Tank(i) => &mut rx.tank_expr[i],
        SpanSlot::Subcatch(i) => &mut rx.subcatch_expr[i],
        SpanSlot::Subsurface(i) => &mut rx.subsurface_laws[i].compiled,
        SpanSlot::Wet(i) => &mut rx.wet_laws[i].compiled,
        SpanSlot::Buildup(i) => &mut rx.buildup_laws[i].compiled,
        SpanSlot::Lid(i) => &mut rx.lid_laws[i].compiled,
    }
}

/// recompileReactionSystem's compile_one: compile into the pool and the slot's span, then
/// the context diagnostic (physics review A1: refuse what would silently read zero).
#[allow(clippy::too_many_arguments)]
fn compile_one(rx: &mut ReactionData, pollutant_names: &[Bytes], max_term: i32, src_expr: &[u8], slot: SpanSlot, where_: &[u8], context: i32, errors: &mut Vec<Bytes>) {
    let mut col = 0;
    let mut span = *slot_mut(rx, slot);
    let mut err;
    {
        let species = std::mem::take(&mut rx.species_name);
        let coefs = std::mem::take(&mut rx.coef_name);
        let terms = std::mem::take(&mut rx.term_name);
        let sym = RxSymbols { species: Some(&species), coefs: Some(&coefs), terms: Some(&terms), pollutants: Some(pollutant_names), max_term };
        err = compile_reaction_expression(src_expr, &sym, &mut rx.token_pool, &mut span, &mut col);
        rx.species_name = species;
        rx.coef_name = coefs;
        rx.term_name = terms;
    }
    *slot_mut(rx, slot) = span;
    if err.is_empty() {
        err = reaction_context_diagnostic(rx, &rx.token_pool, &span, context);
    }
    if !err.is_empty() {
        errors.push(cat(&[where_, b" (col ", col.to_string().as_bytes(), b"): ", &err, b"."]));
    }
}

/// recompileReactionSystem(): the one compile path. True and rx.compiled on success; on
/// failure false, with the partial pool/spans left in place and compiled false.
pub fn recompile_reaction_system(ctx: &mut SimulationContext, errors: &mut Vec<Bytes>) -> bool {
    let ns = ctx.reactions.n_species() as usize;
    {
        let rx = &mut ctx.reactions;
        rx.compiled = false;
        rx.token_pool.clear();
        let nt = rx.term_name.len();
        rx.term_expr.clear();
        rx.term_expr.resize(nt, RxExprSpan::default());
        rx.pipe_expr.clear();
        rx.pipe_expr.resize(ns, RxExprSpan::default());
        rx.tank_expr.clear();
        rx.tank_expr.resize(ns, RxExprSpan::default());

        rx.subcatch_form.resize(ns, ReactionExprForm::None);
        rx.subcatch_expr_src.resize(ns, Vec::new());
        rx.subcatch_expr.clear();
        rx.subcatch_expr.resize(ns, RxExprSpan::default());
    }
    // R4: pollutants are referencable (read-only), in registry order.
    let mut pollutant_names: Vec<Bytes> = Vec::new();
    for p2 in 0..ctx.species_registry.pollutant_count() {
        pollutant_names.push(ctx.species_registry.name(p2).clone());
    }

    let n_errors_in = errors.len();
    let rx = &mut ctx.reactions;

    for i in 0..rx.term_name.len() {
        let src = rx.term_expr_src[i].clone();
        let where_ = cat(&[b"[REACTION_TERMS] '", &rx.term_name[i], b"'"]);
        // forward-only rule: max_term = i
        compile_one(rx, &pollutant_names, i as i32, &src, SpanSlot::Term(i), &where_, 1, errors);
    }
    let max_term = rx.term_name.len() as i32;
    for sidx in 0..rx.n_species() {
        let us = sidx as usize;
        if rx.pipe_form[us] != ReactionExprForm::None {
            let src = rx.pipe_expr_src[us].clone();
            let where_ = cat(&[b"[REACTION_PIPES] '", &rx.species_name[us], b"'"]);
            compile_one(rx, &pollutant_names, max_term, &src, SpanSlot::Pipe(us), &where_, 0, errors);
        }
        if rx.tank_form[us] != ReactionExprForm::None {
            let src = rx.tank_expr_src[us].clone();
            let where_ = cat(&[b"[REACTION_TANKS] '", &rx.species_name[us], b"'"]);
            compile_one(rx, &pollutant_names, max_term, &src, SpanSlot::Tank(us), &where_, 0, errors);
        }
        if rx.subcatch_form[us] != ReactionExprForm::None {
            if rx.species_is_wall[us] != 0 {
                errors.push(b"[REACTION_SUBCATCHMENTS] requires a BULK species.".to_vec());
            } else {
                let src = rx.subcatch_expr_src[us].clone();
                let where_ = cat(&[b"[REACTION_SUBCATCHMENTS] '", &rx.species_name[us], b"'"]);
                compile_one(rx, &pollutant_names, max_term, &src, SpanSlot::Subcatch(us), &where_, 3, errors);
            }
        }
    }
    for k in 0..rx.subsurface_laws.len() {
        let species = rx.subsurface_laws[k].species;
        if species < 0 || species >= rx.n_species() || rx.species_is_wall[species as usize] != 0 {
            errors.push(b"[REACTION_SUBSURFACE] requires a valid BULK species".to_vec());
            continue;
        }
        let src = rx.subsurface_laws[k].expression.clone();
        let where_ = cat(&[b"[REACTION_SUBSURFACE] '", &rx.species_name[species as usize], b"'"]);
        compile_one(rx, &pollutant_names, max_term, &src, SpanSlot::Subsurface(k), &where_, 2, errors);
    }

    // wet_laws, then buildup_laws
    for which in 0..2 {
        let n = if which == 0 { rx.wet_laws.len() } else { rx.buildup_laws.len() };
        for k in 0..n {
            let (species, src) = {
                let law = if which == 0 { &rx.wet_laws[k] } else { &rx.buildup_laws[k] };
                (law.species, law.expression.clone())
            };
            if species < 0 || species >= rx.n_species() || rx.species_is_wall[species as usize] != 0 {
                errors.push(b"Surface reactions require a BULK species.".to_vec());
                continue;
            }
            let where_ = cat(&[b"Surface reaction '", &rx.species_name[species as usize], b"'"]);
            let slot = if which == 0 { SpanSlot::Wet(k) } else { SpanSlot::Buildup(k) };
            compile_one(rx, &pollutant_names, max_term, &src, slot, &where_, 3, errors);
        }
    }
    for k in 0..rx.lid_laws.len() {
        let species = rx.lid_laws[k].species;
        if species < 0 || species >= rx.n_species() || rx.species_is_wall[species as usize] != 0 {
            errors.push(b"[REACTION_LIDS] target must be BULK.".to_vec());
            continue;
        }
        let src = rx.lid_laws[k].expression.clone();
        let where_ = cat(&[b"[REACTION_LIDS] '", &rx.species_name[species as usize], b"'"]);
        compile_one(rx, &pollutant_names, max_term, &src, SpanSlot::Lid(k), &where_, 4, errors);
        let control = rx.lid_laws[k].control;
        let node = control >= 0 && control < ctx.lid_controls.node_layers.len() as i32 && !ctx.lid_controls.node_layers[control as usize].is_empty();
        let count = if node { ctx.lid_controls.node_layers[control as usize].len() as i32 } else { 4 };
        if rx.lid_laws[k].layer > count {
            errors.push(b"[REACTION_LIDS] layer exceeds the control's stack.".to_vec());
        }
    }
    if !rx.wall_geometry && !rx.wall_areas.is_empty() {
        errors.push(b"[REACTION_WALLS] requires WALL_GEOMETRY YES.".to_vec());
    }
    if rx.wall_geometry {
        for s in 0..rx.n_species() {
            let us = s as usize;
            if rx.species_is_wall[us] != 0 {
                if rx.coupling != ReactionCoupling::Full {
                    errors.push(b"WALL_GEOMETRY requires COUPLING FULL.".to_vec());
                }
                if rx.tank_form[us] != ReactionExprForm::None {
                    for node in 0..ctx.node_names.size() {
                        if ctx.nodes.type_[node as usize] == NodeType::Storage && !rx.wall_areas.iter().any(|r| !r.link && r.element == node) {
                            errors.push(cat(&[
                                b"WALL_GEOMETRY tank reactions require an explicit [REACTION_WALLS] NODE area for '",
                                ctx.node_names.name_of(node),
                                b"'.",
                            ]));
                        }
                    }
                }
                for form in [rx.pipe_form[us], rx.tank_form[us]] {
                    if form == ReactionExprForm::Equil || form == ReactionExprForm::Formula {
                        errors.push(b"WALL_GEOMETRY supports RATE or NONE for WALL species; algebraic WALL projections have no conservative wet/dry transfer.".to_vec());
                    }
                }
            }
        }
    }
    if errors.len() != n_errors_in {
        return false;
    }
    rx.compiled = true;
    true
}
