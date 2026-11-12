// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/gw/GwTransportSections.hpp and
// GwTransportSections.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The [GW_*] parsers, their registration and resolveGwTransport (run by open() once the mesh,
// species and series are final).
//   * num() is std::from_chars (core::charconv) consuming the whole token and finite; inum()
//     narrows with static_cast<int>, which is core::cmath::cvt_i32 (INT_MIN out of range).
//   * The handlers reach the router's GwTransportData through a shared handle (twod::Shared),
//     where upstream captures a reference.
//   * std::unordered_set is used only for membership tests (HashSet).

#![allow(clippy::all)]

use std::collections::HashSet;

use crate::core::charconv::from_chars_double;
use crate::core::cmath::cvt_i32;
use crate::core::context::SimulationContext;
use crate::core::cstdio::IfStream;
use crate::core::options::Bytes;
use crate::core::path_resolver::{parent_dir, resolve_relative};
use crate::input::section_registry::SectionRegistry;
use crate::input::tokenizer;
use crate::twod::data::mesh_data::MeshData;
use crate::twod::gw::gw_transport_data::*;
use crate::twod::Shared;

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut v = Bytes::new();
    for p in parts {
        v.extend_from_slice(p);
    }
    v
}

/// iequals(): byte case-insensitive equality with a keyword.
fn iequals(a: &[u8], b: &[u8]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_ascii_uppercase() == y.to_ascii_uppercase())
}

/// num(): a whole finite from_chars double. `v` is written whenever from_chars converts.
fn num(s: &[u8], v: &mut f64) -> bool {
    match from_chars_double(s) {
        (Ok(x), used) => {
            *v = x;
            used == s.len() && x.is_finite()
        }
        (Err(_), _) => false,
    }
}

/// inum(): a whole finite double that is an int.
fn inum(s: &[u8], v: &mut i32) -> bool {
    let mut d = 0.0;
    if !num(s, &mut d) {
        return false;
    }
    *v = cvt_i32(d);
    *v as f64 == d
}

/// parseScope(): `* | TAG <t> | CELL <n>` at tokens[0]; `at` set past it.
fn parse_scope(tokens: &[Bytes], sec: &[u8], scope: &mut GwScope, tag: &mut Bytes, cell: &mut i32, at: &mut usize) -> Bytes {
    if tokens.is_empty() {
        return cat(&[sec, b" empty row"]);
    }
    if tokens[0] == b"*" {
        *scope = GwScope::Global;
        *at = 1;
        return Bytes::new();
    }
    if iequals(&tokens[0], b"TAG") {
        if tokens.len() < 2 {
            return cat(&[sec, b" TAG needs a name"]);
        }
        *scope = GwScope::Tag;
        *tag = tokens[1].clone();
        *at = 2;
        return Bytes::new();
    }
    if iequals(&tokens[0], b"CELL") {
        if tokens.len() < 2 {
            return cat(&[sec, b" CELL needs an index"]);
        }
        let mut n = 0;
        if !inum(&tokens[1], &mut n) || n < 1 {
            return cat(&[sec, b" invalid CELL index (1-based): ", &tokens[1]]);
        }
        *scope = GwScope::Cell;
        *cell = n - 1;
        *at = 2;
        return Bytes::new();
    }
    cat(&[sec, b" row must start with CELL, TAG or '*': ", &tokens[0]])
}

/// valueOrSeries(): a finite number fills `v` (series cleared), otherwise the token is a name.
fn value_or_series(token: &[u8], v: &mut f64, ts: &mut Bytes) {
    if num(token, v) {
        ts.clear();
        return;
    }
    *v = 0.0;
    *ts = token.to_vec();
}

/// parseGwZone()
pub fn parse_gw_zone(token: &[u8], z: &mut GwZone) -> bool {
    if iequals(token, b"SAT") {
        *z = GwZone::Sat;
        return true;
    }
    if iequals(token, b"UNSAT") {
        *z = GwZone::Unsat;
        return true;
    }
    if iequals(token, b"LAYER") {
        *z = GwZone::Layer;
        return true;
    }
    false
}

/// parseGwTransportOptionsLine()
pub fn parse_gw_transport_options_line(tokens: &[Bytes], o: &mut GwTransportOptions) -> Bytes {
    if tokens.is_empty() {
        return Bytes::new();
    }
    if tokens.len() < 2 {
        return b"[GW_TRANSPORT_OPTIONS] needs PARAMETER VALUE".to_vec();
    }
    let k = &tokens[0];
    let v = &tokens[1];

    if iequals(k, b"UNSAT_TRANSPORT") {
        if tokens.len() != 2 {
            return b"[GW_TRANSPORT_OPTIONS] UNSAT_TRANSPORT needs one value".to_vec();
        }
        if iequals(v, b"BULK") {
            o.unsat_transport = b"BULK".to_vec();
        } else if iequals(v, b"LAYERED") {
            o.unsat_transport = b"LAYERED".to_vec();
        } else {
            return b"[GW_TRANSPORT_OPTIONS] UNSAT_TRANSPORT expects BULK or LAYERED".to_vec();
        }
        return Bytes::new();
    }

    if iequals(k, b"THERMAL_MODEL") {
        if iequals(v, b"WATER_ONLY") {
            o.thermal_model = b"WATER_ONLY".to_vec();
        } else if iequals(v, b"POROUS_EQUILIBRIUM") {
            o.thermal_model = b"POROUS_EQUILIBRIUM".to_vec();
        } else {
            return b"[GW_TRANSPORT_OPTIONS] THERMAL_MODEL expects WATER_ONLY or POROUS_EQUILIBRIUM".to_vec();
        }
        o.thermal_authored.insert(b"THERMAL_MODEL".to_vec());
        return Bytes::new();
    }
    if iequals(k, b"THERMAL_LAYERS") {
        let mut count = 0.0;
        if !num(v, &mut count) || count < 1.0 || count > 128.0 || count.floor() != count {
            return b"[GW_TRANSPORT_OPTIONS] THERMAL_LAYERS must be an integer in [1,128]".to_vec();
        }
        o.thermal_layers = count as i32;
        o.thermal_authored.insert(b"THERMAL_LAYERS".to_vec());
        return Bytes::new();
    }
    if iequals(k, b"LATERAL_CONDUCTION") {
        if iequals(v, b"YES") || iequals(v, b"ON") {
            o.lateral_conduction = true;
        } else if iequals(v, b"NO") || iequals(v, b"OFF") {
            o.lateral_conduction = false;
        } else {
            return b"[GW_TRANSPORT_OPTIONS] LATERAL_CONDUCTION expects YES or NO".to_vec();
        }
        o.thermal_authored.insert(b"LATERAL_CONDUCTION".to_vec());
        return Bytes::new();
    }
    if iequals(k, b"SOIL_ET_HEAT") {
        if iequals(v, b"LATENT_SURFACE") {
            o.soil_et_heat = b"LATENT_SURFACE".to_vec();
        } else if iequals(v, b"SENSIBLE_ONLY") {
            o.soil_et_heat = b"SENSIBLE_ONLY".to_vec();
        } else {
            return b"[GW_TRANSPORT_OPTIONS] SOIL_ET_HEAT expects LATENT_SURFACE or SENSIBLE_ONLY".to_vec();
        }
        o.thermal_authored.insert(b"SOIL_ET_HEAT".to_vec());
        return Bytes::new();
    }
    for key in [&b"THERMAL_MAX_STEP"[..], b"WATER_THERMAL_CONDUCTIVITY", b"AIR_THERMAL_CONDUCTIVITY", b"SURFACE_THERMAL_RESISTANCE"] {
        if iequals(k, key) {
            let mut value = 0.0;
            if !num(v, &mut value) || !(value > 0.0) {
                return cat(&[b"[GW_TRANSPORT_OPTIONS] ", key, b" must be finite and positive"]);
            }
            if iequals(k, b"THERMAL_MAX_STEP") {
                o.thermal_max_step = value;
            } else if iequals(k, b"WATER_THERMAL_CONDUCTIVITY") {
                o.water_conductivity = value;
            } else if iequals(k, b"AIR_THERMAL_CONDUCTIVITY") {
                o.air_conductivity = value;
            } else {
                o.surface_resistance = value;
            }
            o.thermal_authored.insert(key.to_vec());
            return Bytes::new();
        }
    }

    let yesno = |field: &mut bool| -> Bytes {
        if iequals(v, b"YES") || iequals(v, b"ON") {
            *field = true;
        } else if iequals(v, b"NO") || iequals(v, b"OFF") {
            *field = false;
        } else {
            return cat(&[b"[GW_TRANSPORT_OPTIONS] ", k, b": expected YES|NO, got ", v]);
        }
        Bytes::new()
    };

    if iequals(k, b"TRANSPORT_POLLUTANTS") {
        return yesno(&mut o.transport_pollutants);
    }
    if iequals(k, b"TRANSPORT_MSX") {
        return yesno(&mut o.transport_msx);
    }
    if iequals(k, b"TRANSPORT_AGE") {
        return yesno(&mut o.transport_age);
    }
    if iequals(k, b"TRANSPORT_TEMPERATURE") {
        return yesno(&mut o.transport_temperature);
    }
    if iequals(k, b"DISPERSION") {
        return yesno(&mut o.dispersion);
    }
    if iequals(k, b"CONDUCTION") {
        let error = yesno(&mut o.conduction);
        if error.is_empty() {
            o.thermal_authored.insert(b"CONDUCTION".to_vec());
        }
        return error;
    }

    if iequals(k, b"SURFACE_THERMAL_BC") {
        if iequals(v, b"ADIABATIC") {
            o.surface_thermal_bc = b"ADIABATIC".to_vec();
            o.surface_thermal_arg.clear();
            return Bytes::new();
        }
        if iequals(v, b"AIR") || iequals(v, b"SURFACE_WATER") {
            o.surface_thermal_bc = if iequals(v, b"AIR") { b"AIR".to_vec() } else { b"SURFACE_WATER".to_vec() };
            o.surface_thermal_arg.clear();
        } else if iequals(v, b"FIXED") || iequals(v, b"TIMESERIES") || iequals(v, b"HEAT_FLUX") || iequals(v, b"FLUX_TIMESERIES") {
            if tokens.len() < 3 {
                return cat(&[b"[GW_TRANSPORT_OPTIONS] SURFACE_THERMAL_BC ", v, b" needs a value or series name"]);
            }
            o.surface_thermal_bc = if iequals(v, b"FIXED") {
                b"FIXED".to_vec()
            } else if iequals(v, b"HEAT_FLUX") {
                b"HEAT_FLUX".to_vec()
            } else if iequals(v, b"FLUX_TIMESERIES") {
                b"FLUX_TIMESERIES".to_vec()
            } else {
                b"TIMESERIES".to_vec()
            };
            o.surface_thermal_arg = tokens[2].clone();
            if o.surface_thermal_bc == b"FIXED" || o.surface_thermal_bc == b"HEAT_FLUX" {
                let mut t = 0.0;
                if !num(&tokens[2], &mut t) {
                    return cat(&[b"[GW_TRANSPORT_OPTIONS] SURFACE_THERMAL_BC FIXED needs a temperature, got ", &tokens[2]]);
                }
            }
        } else {
            return cat(&[b"[GW_TRANSPORT_OPTIONS] SURFACE_THERMAL_BC: expected AIR|SURFACE_WATER|FIXED <degC>|TIMESERIES <name>, got ", v]);
        }
        return Bytes::new();
    }

    if iequals(k, b"DEEP_THERMAL_BC") {
        if iequals(v, b"ADIABATIC") {
            o.deep_thermal_bc = b"ADIABATIC".to_vec();
            o.deep_thermal_arg.clear();
            o.deep_depth = 0.0;
            return Bytes::new();
        }
        if !(iequals(v, b"GEOTHERMAL_FLUX") || iequals(v, b"FIXED_TEMP") || iequals(v, b"TIMESERIES") || iequals(v, b"FLUX_TIMESERIES")) {
            return cat(&[b"[GW_TRANSPORT_OPTIONS] DEEP_THERMAL_BC: expected GEOTHERMAL_FLUX|FIXED_TEMP|TIMESERIES, got ", v]);
        }
        if tokens.len() < 3 {
            return b"[GW_TRANSPORT_OPTIONS] DEEP_THERMAL_BC needs a value or series name".to_vec();
        }
        o.deep_thermal_bc = if iequals(v, b"GEOTHERMAL_FLUX") {
            b"GEOTHERMAL_FLUX".to_vec()
        } else if iequals(v, b"FIXED_TEMP") {
            b"FIXED_TEMP".to_vec()
        } else if iequals(v, b"FLUX_TIMESERIES") {
            b"FLUX_TIMESERIES".to_vec()
        } else {
            b"TIMESERIES".to_vec()
        };
        o.deep_thermal_arg = tokens[2].clone();
        if o.deep_thermal_bc != b"TIMESERIES" && o.deep_thermal_bc != b"FLUX_TIMESERIES" {
            let mut d = 0.0;
            if !num(&tokens[2], &mut d) {
                return cat(&[b"[GW_TRANSPORT_OPTIONS] DEEP_THERMAL_BC ", &o.deep_thermal_bc, b" needs a number, got ", &tokens[2]]);
            }
        }
        // Optional trailing `DEPTH <z>`.
        if tokens.len() >= 5 && iequals(&tokens[3], b"DEPTH") {
            let mut z = 0.0;
            if !num(&tokens[4], &mut z) || z < 0.0 {
                return cat(&[b"[GW_TRANSPORT_OPTIONS] DEEP_THERMAL_BC DEPTH must be >= 0, got ", &tokens[4]]);
            }
            o.deep_depth = z;
        }
        return Bytes::new();
    }

    if iequals(k, b"THERMAL_MIXING") {
        if iequals(v, b"ARITHMETIC") {
            o.thermal_mixing = b"ARITHMETIC".to_vec();
        } else if iequals(v, b"GEOMETRIC") {
            o.thermal_mixing = b"GEOMETRIC".to_vec();
        } else {
            return cat(&[b"[GW_TRANSPORT_OPTIONS] THERMAL_MIXING: expected ARITHMETIC|GEOMETRIC, got ", v]);
        }
        o.thermal_authored.insert(b"THERMAL_MIXING".to_vec());
        return Bytes::new();
    }

    if iequals(k, b"C_DIFF") {
        let mut c = 0.0;
        if !num(v, &mut c) || c <= 0.0 || c > 1.0 {
            return cat(&[b"[GW_TRANSPORT_OPTIONS] C_DIFF must be in (0, 1], got ", v]);
        }
        o.c_diff = c;
        return Bytes::new();
    }

    cat(&[b"Unknown GW_TRANSPORT_OPTIONS parameter: ", k])
}

/// parseGwParamsLine()
pub fn parse_gw_params_line(tokens: &[Bytes], rows: &mut Vec<GwParamsRow>) -> Bytes {
    if tokens.is_empty() {
        return Bytes::new();
    }
    let mut r = GwParamsRow::default();
    let mut at = 0usize;
    let e = parse_scope(tokens, b"[GW_TRANSPORT_PARAMS]", &mut r.scope, &mut r.tag, &mut r.cell, &mut at);
    if !e.is_empty() {
        return e;
    }
    const NAMES: [&[u8]; 9] = [b"rho_s", b"c_s", b"lambda_s", b"a_s", b"alpha_L", b"alpha_T", b"D_m", b"D_v", b"geo_flux"];
    const K_SLOTS: usize = 9;
    for k in at..tokens.len() {
        let i = k - at;
        if i >= K_SLOTS {
            return b"[GW_TRANSPORT_PARAMS] too many columns (max 9 after the scope)".to_vec();
        }
        if tokens[k] == b"-" {
            continue; // unset: keep the default
        }
        let mut v = 0.0;
        if !num(&tokens[k], &mut v) {
            return cat(&[b"[GW_TRANSPORT_PARAMS] invalid ", NAMES[i], b": ", &tokens[k]]);
        }
        // Only geo_flux (which may also be a deep temperature) is signed.
        if v < 0.0 && i != 8 {
            return cat(&[b"[GW_TRANSPORT_PARAMS] ", NAMES[i], b" must be >= 0, got ", &tokens[k]]);
        }
        let slot: &mut f64 = match i {
            0 => &mut r.rho_s,
            1 => &mut r.c_s,
            2 => &mut r.lambda_s,
            3 => &mut r.a_s,
            4 => &mut r.alpha_L,
            5 => &mut r.alpha_T,
            6 => &mut r.D_m,
            7 => &mut r.D_v,
            _ => &mut r.geo_flux,
        };
        *slot = v;
    }
    rows.push(r);
    Bytes::new()
}

/// parseGwSorptionLine()
pub fn parse_gw_sorption_line(tokens: &[Bytes], rows: &mut Vec<GwSorptionRow>) -> Bytes {
    if tokens.is_empty() {
        return Bytes::new();
    }
    let mut r = GwSorptionRow::default();
    let mut at = 0usize;
    let e = parse_scope(tokens, b"[GW_SORPTION]", &mut r.scope, &mut r.tag, &mut r.cell, &mut at);
    if !e.is_empty() {
        return e;
    }
    if tokens.len() < at + 2 {
        return b"[GW_SORPTION] needs SCOPE SPECIES Kd [decay]".to_vec();
    }
    r.species = tokens[at].clone();
    if !num(&tokens[at + 1], &mut r.kd) || r.kd < 0.0 {
        return cat(&[b"[GW_SORPTION] Kd must be a finite value >= 0, got ", &tokens[at + 1]]);
    }
    if tokens.len() > at + 2 && tokens[at + 2] != b"-" {
        if !num(&tokens[at + 2], &mut r.decay) || r.decay < 0.0 {
            return cat(&[b"[GW_SORPTION] decay must be a finite value >= 0, got ", &tokens[at + 2]]);
        }
    }
    rows.push(r);
    Bytes::new()
}

/// parseGwInitialQualityLine()
pub fn parse_gw_initial_quality_line(tokens: &[Bytes], rows: &mut Vec<GwInitialQualityRow>, file: &mut Bytes) -> Bytes {
    if tokens.is_empty() {
        return Bytes::new();
    }
    if iequals(&tokens[0], b"FILE") {
        if tokens.len() < 2 {
            return b"[GW_INITIAL_QUALITY] FILE needs a path".to_vec();
        }
        *file = tokens[1].clone();
        return Bytes::new();
    }
    let mut r = GwInitialQualityRow::default();
    let mut at = 0usize;
    let e = parse_scope(tokens, b"[GW_INITIAL_QUALITY]", &mut r.scope, &mut r.tag, &mut r.cell, &mut at);
    if !e.is_empty() {
        return e;
    }
    if tokens.len() < at + 3 {
        return b"[GW_INITIAL_QUALITY] needs SCOPE ZONE SPECIES VALUE".to_vec();
    }
    if !parse_gw_zone(&tokens[at], &mut r.zone) {
        return cat(&[b"[GW_INITIAL_QUALITY] ZONE must be SAT, UNSAT or LAYER <j>, got ", &tokens[at]]);
    }
    at += 1;
    if r.zone == GwZone::Layer {
        if tokens.len() < at + 3 {
            return b"[GW_INITIAL_QUALITY] LAYER needs <j> SPECIES VALUE".to_vec();
        }
        if !inum(&tokens[at], &mut r.layer) || r.layer < 1 {
            return cat(&[b"[GW_INITIAL_QUALITY] LAYER index must be >= 1, got ", &tokens[at]]);
        }
        at += 1;
    }
    r.species = tokens[at].clone();
    if !num(&tokens[at + 1], &mut r.value) {
        return cat(&[b"[GW_INITIAL_QUALITY] VALUE must be a finite number, got ", &tokens[at + 1]]);
    }
    rows.push(r);
    Bytes::new()
}

/// parseGwBoundaryQualityLine()
pub fn parse_gw_boundary_quality_line(tokens: &[Bytes], rows: &mut Vec<GwBoundaryQualityRow>) -> Bytes {
    if tokens.is_empty() {
        return Bytes::new();
    }
    if tokens.len() < 5 {
        return b"[GW_BOUNDARY_QUALITY] needs CELL EDGE SPECIES KIND VALUE".to_vec();
    }
    let mut r = GwBoundaryQualityRow::default();
    let mut cell = 0;
    if !inum(&tokens[0], &mut cell) || cell < 1 {
        return cat(&[b"[GW_BOUNDARY_QUALITY] invalid CELL index (1-based): ", &tokens[0]]);
    }
    r.cell = cell - 1;
    if !inum(&tokens[1], &mut r.edge) || r.edge < 0 {
        return cat(&[b"[GW_BOUNDARY_QUALITY] invalid EDGE (0-based, 0..nv-1): ", &tokens[1]]);
    }
    r.species = tokens[2].clone();
    let kind = &tokens[3];
    if iequals(kind, b"CONC") {
        r.kind = b"CONC".to_vec();
    } else if iequals(kind, b"TS") {
        r.kind = b"TS".to_vec();
    } else if iequals(kind, b"MASSFLUX") {
        r.kind = b"MASSFLUX".to_vec();
    } else if iequals(kind, b"HEATFLUX") {
        r.kind = b"HEATFLUX".to_vec();
    } else {
        return cat(&[b"[GW_BOUNDARY_QUALITY] KIND must be CONC|TS|MASSFLUX|HEATFLUX, got ", kind]);
    }
    if r.kind == b"TS" {
        r.ts_name = tokens[4].clone();
    } else {
        value_or_series(&tokens[4], &mut r.value, &mut r.ts_name);
    }
    rows.push(r);
    Bytes::new()
}

/// parseGwSourcesLine()
pub fn parse_gw_sources_line(tokens: &[Bytes], rows: &mut Vec<GwSourceRow>) -> Bytes {
    if tokens.is_empty() {
        return Bytes::new();
    }
    if tokens.len() < 4 {
        return b"[GW_SOURCES] needs NAME (* | CELL n | TAG t | XY x y) FLOW value".to_vec();
    }
    let mut r = GwSourceRow::default();
    r.name = tokens[0].clone();
    let mut at = 1usize;
    if tokens[at] == b"*" || iequals(&tokens[at], b"GLOBAL") {
        r.scope = GwScope::Global;
        at += 1;
    } else if iequals(&tokens[at], b"CELL") {
        let mut n = 0;
        if !inum(&tokens[at + 1], &mut n) || n < 1 {
            return cat(&[b"[GW_SOURCES] invalid CELL index (1-based): ", &tokens[at + 1]]);
        }
        r.scope = GwScope::Cell;
        r.cell = n - 1;
        at += 2;
    } else if iequals(&tokens[at], b"TAG") {
        r.scope = GwScope::Tag;
        r.tag = tokens[at + 1].clone();
        at += 2;
    } else if iequals(&tokens[at], b"XY") {
        if tokens.len() < at + 3 {
            return b"[GW_SOURCES] XY needs x y".to_vec();
        }
        if !num(&tokens[at + 1], &mut r.x) || !num(&tokens[at + 2], &mut r.y) {
            return b"[GW_SOURCES] XY needs two numbers".to_vec();
        }
        r.scope = GwScope::Cell; // resolved to a cell at resolve
        r.by_xy = true;
        at += 3;
    } else {
        return cat(&[b"[GW_SOURCES] location must be *, CELL, TAG or XY, got ", &tokens[at]]);
    }

    if at >= tokens.len() || !iequals(&tokens[at], b"FLOW") {
        return cat(&[b"[GW_SOURCES] '", &r.name, b"' needs a FLOW value or series"]);
    }
    if at + 1 >= tokens.len() {
        return b"[GW_SOURCES] FLOW needs a value or series name".to_vec();
    }
    value_or_series(&tokens[at + 1], &mut r.flow, &mut r.flow_ts);
    at += 2;

    if at < tokens.len() && iequals(&tokens[at], b"SCALE") {
        if at + 1 >= tokens.len() || !num(&tokens[at + 1], &mut r.scale) || !r.scale.is_finite() || r.scale < 0.0 {
            return b"[GW_SOURCES] SCALE must be a finite nonnegative multiplier".to_vec();
        }
        at += 2;
    }
    // Trailing species terms: SPECIES (CONC|MASS) value|ts, repeated.
    while at < tokens.len() {
        if at + 2 >= tokens.len() {
            return cat(&[b"[GW_SOURCES] '", &r.name, b"' trailing species term needs SPECIES (CONC|MASS) value"]);
        }
        let mut t = GwSourceSpeciesTerm::default();
        t.species = tokens[at].clone();
        if iequals(&tokens[at + 1], b"CONC") {
            t.kind = b"CONC".to_vec();
        } else if iequals(&tokens[at + 1], b"MASS") {
            t.kind = b"MASS".to_vec();
        } else {
            return cat(&[b"[GW_SOURCES] species term kind must be CONC or MASS, got ", &tokens[at + 1]]);
        }
        value_or_series(&tokens[at + 2], &mut t.value, &mut t.ts_name);
        r.species.push(t);
        at += 3;
    }
    rows.push(r);
    Bytes::new()
}

/// makeGwHandler(): marks the options authored, then parses each non-blank line, stopping at
/// the first diagnostic (error code 5, "[GW] <err> — line: <raw>").
fn make_gw_handler(gw: Shared<GwTransportData>, lp: fn(&mut GwTransportData, &[Bytes]) -> Bytes) -> impl Fn(&mut SimulationContext, &[Bytes]) + 'static {
    move |ctx: &mut SimulationContext, lines: &[Bytes]| {
        let mut g = gw.borrow_mut();
        g.options.authored = true;
        for raw in lines {
            let tokens = tokenizer::tokenize(raw);
            if tokens.is_empty() {
                continue;
            }
            let err = lp(&mut *g, &tokens);
            if !err.is_empty() {
                ctx.error_code = 5; // SWMM_ERR_PARSE
                ctx.error_message = cat(&[b"[GW] ", &err, " \u{2014} line: ".as_bytes(), raw]);
                return;
            }
        }
    }
}

/// registerGwTransportSections()
pub fn register_gw_transport_sections(gw: &Shared<GwTransportData>, registry: &mut SectionRegistry) {
    registry.register_custom(b"GW_TRANSPORT_OPTIONS", make_gw_handler(gw.clone(), |g, t| parse_gw_transport_options_line(t, &mut g.options)));
    registry.register_custom(b"GW_TRANSPORT_PARAMS", make_gw_handler(gw.clone(), |g, t| parse_gw_params_line(t, &mut g.params)));
    registry.register_custom(b"GW_SORPTION", make_gw_handler(gw.clone(), |g, t| parse_gw_sorption_line(t, &mut g.sorption)));
    registry.register_custom(
        b"GW_INITIAL_QUALITY",
        make_gw_handler(gw.clone(), |g, t| {
            g.invalidate_authoring_indexes();
            let GwTransportData { initial_quality, initial_quality_file, .. } = g;
            parse_gw_initial_quality_line(t, initial_quality, initial_quality_file)
        }),
    );
    registry.register_custom(b"GW_BOUNDARY_QUALITY", make_gw_handler(gw.clone(), |g, t| parse_gw_boundary_quality_line(t, &mut g.boundary_quality)));
    registry.register_custom(
        b"GW_SOURCES",
        make_gw_handler(gw.clone(), |g, t| {
            g.invalidate_authoring_indexes();
            parse_gw_sources_line(t, &mut g.sources)
        }),
    );
}

/// nearestCell(): the cell whose centroid (from its vertices) is nearest (x, y); the first
/// minimum wins; -1 for a mesh without cells.
pub fn nearest_cell(mesh: &MeshData, x: f64, y: f64) -> i32 {
    let mut best = -1;
    let mut best_d2 = 0.0;
    for c in 0..mesh.n_cells() {
        let mut cx = 0.0;
        let mut cy = 0.0;
        mesh.cell_centroid_xy(c, &mut cx, &mut cy);
        let d2 = (cx - x) * (cx - x) + (cy - y) * (cy - y);
        if best < 0 || d2 < best_d2 {
            best = c;
            best_d2 = d2;
        }
    }
    best
}

/// resolveGwTransport(): validates the [GW_*] rows against the final mesh, species and series
/// (loading a [GW_INITIAL_QUALITY] FILE and placing XY sources); returns every diagnostic.
pub fn resolve_gw_transport(ctx: &SimulationContext, mesh: &MeshData, gw: &mut GwTransportData) -> Vec<Bytes> {
    let mut errs: Vec<Bytes> = Vec::new();
    if gw.empty() {
        return errs;
    }
    gw.invalidate_authoring_indexes();

    let n_cells = mesh.n_cells();

    // `[GW_INITIAL_QUALITY] FILE <csv>`: element,zone,species,value.
    if !gw.initial_quality_file.is_empty() {
        // Drop what a previous resolve loaded (re-resolve is idempotent).
        let mut i = gw.initial_quality.len() as i32 - 1;
        while i >= 0 {
            if gw.initial_quality[i as usize].layer == -2 {
                gw.initial_quality.remove(i as usize);
            }
            i -= 1;
        }
        let dir = parent_dir(&ctx.inp_file_path);
        let path = resolve_relative(&gw.initial_quality_file, &dir);
        match IfStream::open(&path) {
            None => {
                errs.push(cat(&[b"[GW_INITIAL_QUALITY] FILE '", &gw.initial_quality_file, b"' not found or unreadable (", &path, b")."]));
            }
            Some(mut input) => {
                let mut line: Bytes = Bytes::new();
                let mut lineno = 0;
                while input.getline(&mut line) {
                    lineno += 1;
                    for ch in line.iter_mut() {
                        if *ch == b',' || *ch == b';' || *ch == b'\t' {
                            *ch = b' ';
                        }
                    }
                    let tok = tokenizer::tokenize(&line);
                    if tok.is_empty() {
                        continue;
                    }
                    let mut one: Vec<GwInitialQualityRow> = Vec::new();
                    let mut ignored: Bytes = Bytes::new();
                    let e = parse_gw_initial_quality_line(&tok, &mut one, &mut ignored);
                    if !e.is_empty() {
                        if lineno == 1 {
                            continue; // header
                        }
                        errs.push(cat(&[b"[GW_INITIAL_QUALITY] FILE '", &gw.initial_quality_file, b"' line ", lineno.to_string().as_bytes(), b": ", &e]));
                        continue;
                    }
                    for mut r in one {
                        if r.zone != GwZone::Layer {
                            r.layer = -2; // "from the file"
                        }
                        gw.initial_quality.push(r);
                    }
                }
            }
        }
    }

    let check_cell = |errs: &mut Vec<Bytes>, cell: i32, sec: &[u8]| {
        if cell < 0 {
            return; // GLOBAL / TAG rows carry -1
        }
        if cell >= n_cells {
            errs.push(cat(&[sec, b" CELL ", (cell + 1).to_string().as_bytes(), b" is off the mesh (", n_cells.to_string().as_bytes(), b" cells)."]));
        }
    };
    let known_tags: HashSet<&[u8]> = mesh.tri_tag.iter().map(|t| t.as_slice()).collect();
    let check_tag = |errs: &mut Vec<Bytes>, tag: &[u8], sec: &[u8]| {
        if tag.is_empty() || known_tags.contains(tag) {
            return;
        }
        errs.push(cat(&[sec, b" TAG '", tag, b"' matches no cell tag on the mesh."]));
    };
    let check_species = |errs: &mut Vec<Bytes>, name: &[u8], sec: &[u8]| -> bool {
        if name == b"__WATER_AGE__" || name == b"__TEMPERATURE__" {
            return true;
        }
        if ctx.species_registry.find(name) >= 0 {
            return true;
        }
        errs.push(cat(&[
            sec,
            b" unknown species '",
            name,
            "' \u{2014} declare it in [POLLUTANTS] or the reactions component, or name __WATER_AGE__ / __TEMPERATURE__.".as_bytes(),
        ]));
        false
    };
    let check_series = |errs: &mut Vec<Bytes>, name: &[u8], sec: &[u8]| {
        if name.is_empty() {
            return;
        }
        if ctx.find_timeseries(name) < 0 {
            errs.push(cat(&[sec, b" unknown time series '", name, b"'."]));
        }
    };

    for r in &gw.params {
        check_cell(&mut errs, r.cell, b"[GW_TRANSPORT_PARAMS]");
        check_tag(&mut errs, &r.tag, b"[GW_TRANSPORT_PARAMS]");
    }
    for r in &gw.sorption {
        check_cell(&mut errs, r.cell, b"[GW_SORPTION]");
        check_tag(&mut errs, &r.tag, b"[GW_SORPTION]");
        check_species(&mut errs, &r.species, b"[GW_SORPTION]");
    }
    for r in &gw.initial_quality {
        check_cell(&mut errs, r.cell, b"[GW_INITIAL_QUALITY]");
        check_tag(&mut errs, &r.tag, b"[GW_INITIAL_QUALITY]");
        if !check_species(&mut errs, &r.species, b"[GW_INITIAL_QUALITY]") {
            continue;
        }
        let signed_row = r.species == b"__TEMPERATURE__" || r.species == b"__WATER_AGE__";
        if !signed_row && r.value < 0.0 {
            errs.push(cat(&[b"[GW_INITIAL_QUALITY] negative value for '", &r.species, b"'."]));
        }
    }
    for r in &gw.boundary_quality {
        check_cell(&mut errs, r.cell, b"[GW_BOUNDARY_QUALITY]");
        // EDGE 0..nv-1 of ITS cell: a triangle has 3 local edges, a quad 4.
        if r.cell >= 0 && r.cell < n_cells {
            let nv = mesh.cell_nv[r.cell as usize] as i32;
            if r.edge >= nv {
                errs.push(cat(&[
                    b"[GW_BOUNDARY_QUALITY] EDGE ",
                    r.edge.to_string().as_bytes(),
                    b" on CELL ",
                    (r.cell + 1).to_string().as_bytes(),
                    b" (",
                    nv.to_string().as_bytes(),
                    " vertices) \u{2014} EDGE is 0..".as_bytes(),
                    (nv - 1).to_string().as_bytes(),
                    b".",
                ]));
            }
        }
        check_species(&mut errs, &r.species, b"[GW_BOUNDARY_QUALITY]");
        check_series(&mut errs, &r.ts_name, b"[GW_BOUNDARY_QUALITY]");
    }
    for i in 0..gw.sources.len() {
        if gw.sources[i].by_xy {
            let best = nearest_cell(mesh, gw.sources[i].x, gw.sources[i].y);
            if best < 0 {
                errs.push(cat(&[b"[GW_SOURCES] '", &gw.sources[i].name, b"' XY cannot be placed: the mesh has no cells."]));
            }
            gw.sources[i].cell = best;
        }
        let r = &gw.sources[i];
        check_cell(&mut errs, r.cell, b"[GW_SOURCES]");
        check_tag(&mut errs, &r.tag, b"[GW_SOURCES]");
        check_series(&mut errs, &r.flow_ts, b"[GW_SOURCES]");
        for t in &r.species {
            check_species(&mut errs, &t.species, b"[GW_SOURCES]");
            check_series(&mut errs, &t.ts_name, b"[GW_SOURCES]");
        }
    }

    // Duplicate names are an error, not last-wins.
    let mut source_names: HashSet<Bytes> = HashSet::new();
    for source in &gw.sources {
        if !source_names.insert(source.name.clone()) {
            errs.push(cat(&[b"[GW_SOURCES] duplicate source name '", &source.name, b"'."]));
        }
    }

    errs
}

/// gwTransportInertWarning()
pub fn gw_transport_inert_warning(gw: &GwTransportData) -> Bytes {
    if gw.empty() {
        return Bytes::new();
    }
    b"[GW_*] subsurface transport is AUTHORED but INERT this run: the rows configure the two-zone [2D_AQUIFER] kernel, and no [2D_AQUIFER] resolved on this model. The rows are validated, kept and written back unchanged.".to_vec()
}
