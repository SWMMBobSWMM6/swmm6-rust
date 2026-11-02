// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/components/EulerianArdComponent/
// ArdConfig.hpp and ArdConfig.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit
// 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The `org.hydrocouple.openswmm.transport.ard` component: its config apply hook
// (model.ard), the FV_DISPERSION bypass warning and the deferred boundary/source row
// resolution that SWMMEngine::open runs after every component has applied. The save hook
// (saveArdConfig) is not part of open() and is not translated here.

use crate::core::context::{ProcessComponentSpec, SimulationContext};
use crate::core::cstdlib::{c_str, strtod};
use crate::core::fs_path::FsPath;
use crate::core::options::fv::{Limiter, ScalarScheme};
use crate::core::options::{Bytes, QualitySolverKind};
use crate::core::process_components::{ComponentConfigSections, ProcessComponentRegistry};
use crate::data::ard_config::{ArdConfigData, ArdDispersionMode, ArdTransportRow, K_LITERS_PER_FT3};

const K_ARD_ID: &str = "org.hydrocouple.openswmm.transport.ard";

/// Split on space and tab only (no quotes, no commas).
pub(crate) fn tokenize_ws(line: &[u8]) -> Vec<Bytes> {
    let mut toks = Vec::new();
    let mut cur: Bytes = Vec::new();
    for &ch in line {
        if ch == b' ' || ch == b'\t' {
            if !cur.is_empty() {
                toks.push(std::mem::take(&mut cur));
            }
        } else {
            cur.push(ch);
        }
    }
    if !cur.is_empty() {
        toks.push(cur);
    }
    toks
}

pub(crate) fn upper(s: &[u8]) -> Bytes {
    s.iter().map(|c| c.to_ascii_uppercase()).collect()
}

/// strtod with full consumption of the C string (`*end == '\0' && end != s`).
pub(crate) fn parse_value(tok: &[u8], out: &mut f64) -> bool {
    let (v, n) = strtod(tok);
    *out = v;
    n != 0 && n == c_str(tok).len()
}

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

fn apply_ard_sections(ctx: &mut SimulationContext, _spec: &ProcessComponentSpec, config: &ComponentConfigSections, errors: &mut Vec<Bytes>) {
    // Reset wholesale so a reopened model never inherits stale rows.
    ctx.ard_config = ArdConfigData::default();

    let errors_before = errors.len();

    for sec in &config.sections {
        let tag: &[u8] = &sec.0;

        if tag == b"TRANSPORT_OPTIONS" {
            for line in &sec.1 {
                let toks = tokenize_ws(line);
                if toks.is_empty() {
                    continue;
                }
                let key = upper(&toks[0]);
                if key == b"DISPERSION" {
                    if toks.len() != 2 {
                        errors.push(cat(&[b"[TRANSPORT_OPTIONS] DISPERSION expects exactly one argument (OFF | FISCHER | value): '", line, b"'."]));
                        continue;
                    }
                    let val = upper(&toks[1]);
                    if val == b"OFF" {
                        ctx.ard_config.dispersion_mode = ArdDispersionMode::Off;
                    } else if val == b"FISCHER" {
                        ctx.ard_config.dispersion_mode = ArdDispersionMode::Fischer;
                    } else {
                        let mut d = 0.0;
                        if !parse_value(&toks[1], &mut d) || d < 0.0 {
                            errors.push(cat(&[
                                b"[TRANSPORT_OPTIONS] DISPERSION '",
                                &toks[1],
                                "' is not OFF, FISCHER, or a non-negative coefficient (display units, len\u{b2}/s).".as_bytes(),
                            ]));
                            continue;
                        }
                        ctx.ard_config.dispersion_mode = ArdDispersionMode::Value;
                        ctx.ard_config.dispersion_value = d;
                    }
                } else if key == b"SCALAR_SCHEME" {
                    if toks.len() != 2 {
                        errors.push(b"[TRANSPORT_OPTIONS] SCALAR_SCHEME expects UPWIND | MUSCL | QUICKEST_ULTIMATE.".to_vec());
                        continue;
                    }
                    let v = upper(&toks[1]);
                    ctx.ard_config.sets_scalar_scheme = true; // IO3c
                    if v == b"UPWIND" {
                        ctx.options.fv.scalar_scheme = ScalarScheme::Upwind;
                    } else if v == b"MUSCL" {
                        ctx.options.fv.scalar_scheme = ScalarScheme::Muscl;
                    } else if v == b"QUICKEST_ULTIMATE" {
                        ctx.options.fv.scalar_scheme = ScalarScheme::QuickestUltimate;
                    } else {
                        errors.push(cat(&[b"[TRANSPORT_OPTIONS] SCALAR_SCHEME '", &toks[1], b"' is not UPWIND, MUSCL, or QUICKEST_ULTIMATE."]));
                    }
                } else if key == b"LIMITER" {
                    if toks.len() != 2 {
                        errors.push(b"[TRANSPORT_OPTIONS] LIMITER expects MINMOD | VANLEER | SUPERBEE.".to_vec());
                        continue;
                    }
                    let v = upper(&toks[1]);
                    ctx.ard_config.sets_limiter = true; // IO3c
                    if v == b"MINMOD" {
                        ctx.options.fv.limiter = Limiter::Minmod;
                    } else if v == b"VANLEER" {
                        ctx.options.fv.limiter = Limiter::Vanleer;
                    } else if v == b"SUPERBEE" {
                        ctx.options.fv.limiter = Limiter::Superbee;
                    } else {
                        errors.push(cat(&[b"[TRANSPORT_OPTIONS] LIMITER '", &toks[1], b"' is not MINMOD, VANLEER, or SUPERBEE."]));
                    }
                } else if key == b"TARGET_DX" {
                    let mut v = 0.0;
                    if toks.len() != 2 || !parse_value(&toks[1], &mut v) || v <= 0.0 {
                        errors.push(b"[TRANSPORT_OPTIONS] TARGET_DX expects a positive cell length (display length units).".to_vec());
                        continue;
                    }
                    ctx.ard_config.target_dx = v;
                } else if key == b"DETAILED_OUTPUT" {
                    if toks.len() != 2 {
                        errors.push(b"[TRANSPORT_OPTIONS] DETAILED_OUTPUT expects one path argument.".to_vec());
                        continue;
                    }
                    // Relative paths resolve against the CONFIG file's own directory.
                    let mut p = FsPath::new(&toks[1]);
                    if !p.is_absolute() {
                        let base = FsPath::new(&config.source_path).parent_path();
                        if !base.empty() {
                            p = base.join(&p);
                        }
                    }
                    ctx.ard_config.detailed_output_path = p.native().to_vec();
                } else {
                    errors.push(cat(&[
                        b"[TRANSPORT_OPTIONS] unknown key '",
                        &toks[0],
                        b"' (E5b recognizes: DISPERSION, SCALAR_SCHEME, LIMITER, TARGET_DX, DETAILED_OUTPUT).",
                    ]));
                }
            }
        } else if tag == b"CONDUIT_DISPERSION" {
            for line in &sec.1 {
                let toks = tokenize_ws(line);
                if toks.len() != 2 {
                    errors.push(cat(&[b"[CONDUIT_DISPERSION] expects '<conduit_name> <value>': '", line, b"'."]));
                    continue;
                }
                let link = ctx.link_names.find(&toks[0]);
                if link < 0 {
                    errors.push(cat(&[b"[CONDUIT_DISPERSION] unknown conduit '", &toks[0], b"'."]));
                    continue;
                }
                if ctx.link_subtypes.conduit_row(link) < 0 {
                    errors.push(cat(&[
                        b"[CONDUIT_DISPERSION] link '",
                        &toks[0],
                        "' is not a conduit \u{2014} dispersion applies to conduits only (structures are zero-volume passthrough).".as_bytes(),
                    ]));
                    continue;
                }
                let mut d = 0.0;
                if !parse_value(&toks[1], &mut d) || d < 0.0 {
                    errors.push(cat(&[
                        b"[CONDUIT_DISPERSION] '",
                        &toks[0],
                        b"': value '",
                        &toks[1],
                        "' is not a non-negative coefficient (display units, len\u{b2}/s).".as_bytes(),
                    ]));
                    continue;
                }
                let mut dup = false;
                for &prev in &ctx.ard_config.conduit_disp_link {
                    if prev == link {
                        dup = true;
                        break;
                    }
                }
                if dup {
                    errors.push(cat(&[b"[CONDUIT_DISPERSION] duplicate row for conduit '", &toks[0], b"'."]));
                    continue;
                }
                ctx.ard_config.conduit_disp_link.push(link);
                ctx.ard_config.conduit_disp_value.push(d);
            }
        } else if tag == b"TRANSPORT_BOUNDARIES" || tag == b"TRANSPORT_SOURCES" {
            let is_bc = tag == b"TRANSPORT_BOUNDARIES";
            for line in &sec.1 {
                let toks = tokenize_ws(line);
                if toks.len() != 4 {
                    errors.push(cat(&[
                        b"[",
                        tag,
                        b"] expects '<",
                        if is_bc { b"node" } else { b"conduit" },
                        b"> <species> VALUE <v> | TIMESERIES <name>': '",
                        line,
                        b"'.",
                    ]));
                    continue;
                }
                let mut row = ArdTransportRow { element: toks[0].clone(), species: toks[1].clone(), ..Default::default() };
                let mode = upper(&toks[2]);
                if mode == b"VALUE" {
                    let mut v = 0.0;
                    if !parse_value(&toks[3], &mut v) || (is_bc && v < 0.0) {
                        errors.push(cat(&[
                            b"[",
                            tag,
                            b"] '",
                            &toks[0],
                            b"': VALUE '",
                            &toks[3],
                            if is_bc { &b"' is not a non-negative number."[..] } else { &b"' is not a number."[..] },
                        ]));
                        continue;
                    }
                    row.value = v;
                } else if mode == b"TIMESERIES" {
                    row.is_ts = true;
                    row.ts_name = toks[3].clone();
                } else {
                    errors.push(cat(&[b"[", tag, b"] '", &toks[0], b"': mode '", &toks[2], b"' is not VALUE or TIMESERIES."]));
                    continue;
                }
                if is_bc {
                    ctx.ard_config.boundary_rows.push(row);
                } else {
                    ctx.ard_config.source_rows.push(row);
                }
            }
        } else if tag == b"STORAGE_MIXING" {
            errors.push("[STORAGE_MIXING] arrives with plan phase E2b (Eulerian ARD plan) \u{2014} not yet implemented.".as_bytes().to_vec());
        } else {
            errors.push(cat(&[
                b"model.ard: unknown section [",
                tag,
                b"] (E5a recognizes: [TRANSPORT_OPTIONS] [CONDUIT_DISPERSION] [TRANSPORT_BOUNDARIES] [TRANSPORT_SOURCES]).",
            ]));
        }
    }

    if errors.len() != errors_before {
        ctx.ard_config = ArdConfigData::default(); // never half-apply
        return;
    }
    ctx.ard_config.configured = true;

    // Silent-bypass enumeration.
    if ctx.ard_config.any_engine_content() {
        if ctx.options.ignore_quality {
            ctx.warnings.push(
                "A transport.ard component is configured but IGNORE_QUALITY is YES \u{2014} the quality stage does not run, so none of its dispersion/boundary/source configuration applies this simulation."
                    .as_bytes()
                    .to_vec(),
            );
        } else if ctx.options.quality_solver != QualitySolverKind::EulerianArd {
            ctx.warnings.push(
                "A transport.ard component is configured but QUALITY_SOLVER is not EULERIAN_ARD \u{2014} its dispersion/boundary/source configuration is inert under the selected engine. Set [OPTIONS] QUALITY_SOLVER EULERIAN_ARD to activate it."
                    .as_bytes()
                    .to_vec(),
            );
        } else if ctx.n_pollutants() == 0 && ctx.ard_config.any_dispersion() {
            let mut has_reactions_component = false;
            for spec in &ctx.process_component_specs {
                if spec.id == b"org.hydrocouple.openswmm.reactions" {
                    has_reactions_component = true;
                }
            }
            if !has_reactions_component {
                ctx.warnings.push(
                    "A transport.ard component configures dispersion but the model has no [POLLUTANTS] and no reactions component \u{2014} there are no transported species, so nothing disperses this simulation."
                        .as_bytes()
                        .to_vec(),
                );
            }
        }
    }
}

pub fn warn_if_fv_dispersion_key_ignored(ctx: &mut SimulationContext) {
    if ctx.options.fv.dispersion == 0.0 {
        return;
    }
    if ctx.options.quality_solver != QualitySolverKind::EulerianArd {
        return;
    }
    ctx.warnings.push(
        "[OPTIONS] FV_DISPERSION is set but does not configure the Eulerian ARD engine \u{2014} it reads dispersion from the transport.ard component (model.ard: [TRANSPORT_OPTIONS] DISPERSION OFF|FISCHER|<value> and [CONDUIT_DISPERSION] rows). No dispersion is applied from this key. Unifying the FV_*/ARD_* option surface arrives with plan phase E5."
            .as_bytes()
            .to_vec(),
    );
}

pub fn ard_boundaries_need_external_volumes(ctx: &SimulationContext) -> bool {
    ctx.options.quality_solver == QualitySolverKind::EulerianArd && !ctx.options.ignore_quality && !ctx.ard_config.bc_node.is_empty()
}

pub fn resolve_ard_transport_rows(ctx: &mut SimulationContext, errors: &mut Vec<Bytes>) {
    {
        let cfg = &mut ctx.ard_config;
        cfg.bc_node.clear();
        cfg.bc_msx.clear();
        cfg.bc_value.clear();
        cfg.bc_ts.clear();
        cfg.src_link.clear();
        cfg.src_msx.clear();
        cfg.src_value.clear();
        cfg.src_ts.clear();
        cfg.transport_rows_resolved = false;
        if !cfg.configured || !cfg.any_transport_rows() {
            cfg.transport_rows_resolved = true;
            return;
        }
    }

    // Species must be MSX.
    let resolve_species = |ctx: &SimulationContext, row: &ArdTransportRow, sec: &[u8], errors: &mut Vec<Bytes>| -> i32 {
        let s = ctx.reactions.find_species(&row.species);
        if s >= 0 {
            return s;
        }
        if ctx.species_registry.find(&row.species) >= 0 {
            errors.push(cat(&[
                b"[",
                sec,
                b"] '",
                &row.element,
                b"': '",
                &row.species,
                "' is a pollutant \u{2014} pollutant loading uses the legacy pathways ([INFLOWS], DWF, washoff); transport boundaries/sources apply to MSX species only (E5a).".as_bytes(),
            ]));
        } else {
            errors.push(cat(&[
                b"[",
                sec,
                b"] '",
                &row.element,
                b"': unknown species '",
                &row.species,
                "' \u{2014} declare it in the reactions component's [REACTION_SPECIES].".as_bytes(),
            ]));
        }
        -1
    };
    let resolve_ts = |ctx: &SimulationContext, row: &ArdTransportRow, sec: &[u8], errors: &mut Vec<Bytes>| -> i32 {
        let t = ctx.tables.find_by_kind(&row.ts_name, true);
        if t < 0 {
            errors.push(cat(&[b"[", sec, b"] '", &row.element, b"': unknown timeseries '", &row.ts_name, b"'."]));
        }
        t
    };

    let boundary_rows = ctx.ard_config.boundary_rows.clone();
    for row in &boundary_rows {
        let nd = ctx.node_names.find(&row.element);
        if nd < 0 {
            errors.push(cat(&[b"[TRANSPORT_BOUNDARIES] unknown node '", &row.element, b"'."]));
            continue;
        }
        let s = resolve_species(ctx, row, b"TRANSPORT_BOUNDARIES", errors);
        if s < 0 {
            continue;
        }
        let mut ts = -1;
        if row.is_ts {
            ts = resolve_ts(ctx, row, b"TRANSPORT_BOUNDARIES", errors);
            if ts < 0 {
                continue;
            }
        }
        let cfg = &mut ctx.ard_config;
        let mut dup = false;
        for i in 0..cfg.bc_node.len() {
            if cfg.bc_node[i] == nd && cfg.bc_msx[i] == s {
                errors.push(cat(&[b"[TRANSPORT_BOUNDARIES] duplicate row for node '", &row.element, b"' species '", &row.species, b"'."]));
                dup = true;
                break;
            }
        }
        if dup {
            continue;
        }
        cfg.bc_node.push(nd);
        cfg.bc_msx.push(s);
        cfg.bc_value.push(row.value);
        cfg.bc_ts.push(if row.is_ts { ts } else { -1 });
    }

    let source_rows = ctx.ard_config.source_rows.clone();
    for row in &source_rows {
        let link = ctx.link_names.find(&row.element);
        if link < 0 || ctx.link_subtypes.conduit_row(link) < 0 {
            errors.push(cat(&[
                b"[TRANSPORT_SOURCES] '",
                &row.element,
                "' is not a conduit \u{2014} distributed sources apply to conduits.".as_bytes(),
            ]));
            continue;
        }
        let s = resolve_species(ctx, row, b"TRANSPORT_SOURCES", errors);
        if s < 0 {
            continue;
        }
        let mut ts = -1;
        if row.is_ts {
            ts = resolve_ts(ctx, row, b"TRANSPORT_SOURCES", errors);
            if ts < 0 {
                continue;
            }
        }
        let mut dup = false;
        for i in 0..ctx.ard_config.src_link.len() {
            if ctx.ard_config.src_link[i] == link && ctx.ard_config.src_msx[i] == s {
                errors.push(cat(&[b"[TRANSPORT_SOURCES] duplicate row for conduit '", &row.element, b"' species '", &row.species, b"'."]));
                dup = true;
                break;
            }
        }
        if dup {
            continue;
        }
        if !row.is_ts && row.value < 0.0 {
            ctx.warnings.push(cat(&[
                b"[TRANSPORT_SOURCES] '",
                &row.element,
                b"' species '",
                &row.species,
                b"': negative rate is treated as EXTRACTION (D-NS1). Removal is clamped to the mass the conduit's cells actually hold.",
            ]));
        }
        let cfg = &mut ctx.ard_config;
        cfg.src_link.push(link);
        cfg.src_msx.push(s);
        // VALUE rates arrive in species mass units per second; divide by L/ft3 once here.
        cfg.src_value.push(row.value / K_LITERS_PER_FT3);
        cfg.src_ts.push(if row.is_ts { ts } else { -1 });
    }

    ctx.ard_config.transport_rows_resolved = true;
}

pub fn register_ard_component(reg: &mut ProcessComponentRegistry) {
    reg.register_component(
        K_ARD_ID.as_bytes(),
        b"Eulerian ARD transport engine configuration (E3: dispersion; boundaries/sources arrive with E5)",
        apply_ard_sections,
    );
}
