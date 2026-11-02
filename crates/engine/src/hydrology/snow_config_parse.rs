// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/surface/UebConfig.cpp,
// src/engine/hydrology/surface/MeshSnowConfig.cpp and (Parameters::validate only)
// src/engine/hydrology/surface/UebSnow.cpp (HydroCouple/openswmm.engine, branch swmm6_rel,
// commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the Apache License,
// Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Only the parsing entry points the [TEMPERATURE], [SNOWPACKS], [UEB_OPTIONS] and
// [2D_SNOW*] section handlers call: setUebOption, setSnowTemperature, setSnowpackUeb,
// setMeshSnowOption and setMeshSnowRow, with Parameters::validate (reached only when
// setUebOption is asked to validate). The writers (uebOptionRows, meshSnowRows, ...) and
// resolveMeshSnowRow are not here. Each function returns the diagnostic as bytes; empty
// means accepted, as upstream's empty std::string does.
//
// Differences in form only: upstream's member-pointer table for the numeric UEB keys is a
// match returning the field (`param_slot`), in the same key order. The two files' private
// `number` helpers differ (UebConfig's rejects an empty token, MeshSnowConfig's does not,
// because parse_double_strict reads "" as 0), so both are kept, as `number_ueb` and
// `number_mesh`. A strict SNOW_STEP clock is split on ':' as std::getline splits it: a
// trailing ':' yields no final empty field.

use std::collections::BTreeSet;

use crate::core::options::Bytes;
use crate::hydrology::ueb_types::{
    ClimateSite, MeshSnowConfig, MeshSnowConfigMode, MeshSnowRow, Parameters, RadiationSourceMode, SnowConfig,
    SnowKind, SnowpackConfig, TemperatureMethod,
};
use crate::input::parse_utils::parse_double_strict;
use crate::input::tokenizer;

// ===========================================================================
// UebSnow.cpp: Parameters::validate
// ===========================================================================

fn finite(v: f64) -> bool {
    v.is_finite()
}
fn nonnegative(v: f64) -> bool {
    finite(v) && v >= 0.0
}
fn fraction(v: f64) -> bool {
    nonnegative(v) && v <= 1.0
}

/// Parameters::validate: empty if every field qualifies, else the first diagnostic.
pub fn parameters_validate(p: &Parameters) -> Bytes {
    if p.temperature != TemperatureMethod::ModifiedForceRestore && p.temperature != TemperatureMethod::Equilibrium {
        return b"Unknown UEB surface temperature method.".to_vec();
    }
    for x in [
        p.rain_temperature,
        p.snow_temperature,
        p.emissivity,
        p.roughness,
        p.measurement_height,
        p.snow_density,
        p.soil_density,
        p.soil_heat_capacity,
        p.active_soil_depth,
        p.liquid_retention,
        p.drainage_conductivity,
        p.equilibrium_conductance,
        p.snow_conductivity,
        p.soil_conductivity,
        p.low_frequency,
        p.amplitude,
        p.stability,
        p.ground_albedo,
        p.visible_albedo,
        p.infrared_albedo,
        p.new_snow_depth,
        p.age_dust,
        p.visible_age,
        p.infrared_age,
        p.max_step,
    ] {
        if !finite(x) {
            return b"Nonfinite UEB parameter.".to_vec();
        }
    }
    if p.rain_temperature <= p.snow_temperature
        || p.snow_temperature <= -273.15
        || !fraction(p.emissivity)
        || p.roughness <= 0.0
        || p.measurement_height <= p.roughness
        || p.snow_density <= 0.0
        || p.snow_density >= 917.0
        || p.soil_density <= 0.0
        || p.soil_heat_capacity <= 0.0
        || p.active_soil_depth <= 0.0
        || p.liquid_retention < 0.0
        || p.liquid_retention >= 1000.0 / p.snow_density - 1000.0 / 917.0
        || p.drainage_conductivity < 0.0
        || p.equilibrium_conductance <= 0.0
        || p.snow_conductivity <= 0.0
        || p.soil_conductivity <= 0.0
        || p.low_frequency <= 0.0
        || p.amplitude <= 0.0
        || !fraction(p.stability)
        || !fraction(p.ground_albedo)
        || !fraction(p.visible_albedo)
        || !fraction(p.infrared_albedo)
        || p.new_snow_depth < 0.0
        || p.age_dust < 0.0
        || !fraction(p.visible_age)
        || !fraction(p.infrared_age)
        || p.max_step < 1.0
        || p.max_step > 21600.0
    {
        return b"Invalid UEB parameter range.".to_vec();
    }
    Vec::new()
}

// ===========================================================================
// UebConfig.cpp
// ===========================================================================

fn upper(s: &[u8]) -> Bytes {
    tokenizer::to_upper(s)
}

/// UebConfig.cpp `number`: nonempty, strictly parsed and finite. `v` may change before a
/// false return (a nonfinite value is stored first), as upstream's may.
fn number_ueb(s: &[u8], v: &mut f64) -> bool {
    !s.is_empty() && parse_double_strict(s, v) && v.is_finite()
}

/// UebConfig.cpp `model`: UEB or DEGREE_DAY, any case; `out` kept on false.
fn model(s: &[u8], out: &mut SnowKind) -> bool {
    let k = upper(s);
    if k == b"UEB" {
        *out = SnowKind::Ueb;
    } else if k == b"DEGREE_DAY" {
        *out = SnowKind::DegreeDay;
    } else {
        return false;
    }
    true
}

/// upstream's `parameters[]` table, key -> member, in its order.
fn param_slot<'a>(p: &'a mut Parameters, key: &[u8]) -> Option<&'a mut f64> {
    Some(match key {
        b"Z0" => &mut p.roughness,
        b"ZM" => &mut p.measurement_height,
        b"KS" => &mut p.equilibrium_conductance,
        b"KSAT" => &mut p.drainage_conductivity,
        b"LC" => &mut p.liquid_retention,
        b"DE" => &mut p.active_soil_depth,
        b"RHO_S" => &mut p.snow_density,
        b"RHO_G" => &mut p.soil_density,
        b"CG" => &mut p.soil_heat_capacity,
        b"LAMBDA_S" => &mut p.snow_conductivity,
        b"LAMBDA_G" => &mut p.soil_conductivity,
        b"OMEGA_LF" => &mut p.low_frequency,
        b"R" => &mut p.amplitude,
        b"FSTAB" => &mut p.stability,
        b"TR" => &mut p.rain_temperature,
        b"TB" => &mut p.snow_temperature,
        b"ABG" => &mut p.ground_albedo,
        b"ALPHA_VO" => &mut p.visible_albedo,
        b"ALPHA_IRO" => &mut p.infrared_albedo,
        b"EPSILON" => &mut p.emissivity,
        b"DNEWS" => &mut p.new_snow_depth,
        b"R3" => &mut p.age_dust,
        b"CV" => &mut p.visible_age,
        b"CIR" => &mut p.infrared_age,
        _ => return None,
    })
}

fn cat(parts: &[&[u8]]) -> Bytes {
    parts.concat()
}

/// setUebOption: parse one UEB_OPTIONS row into a candidate and publish it only when it is
/// accepted (`validate` checks the whole candidate Parameters first).
pub fn set_ueb_option(cfg: &mut SnowConfig, tok: &[Bytes], validate: bool) -> Bytes {
    if tok.len() < 2 {
        return b"UEB option requires a key and value.".to_vec();
    }
    let key = upper(&tok[0]);
    let mut next = cfg.clone();
    if key == b"BC_MONTHLY" {
        let mut month = 0.0f64;
        let mut value = 0.0f64;
        if tok.len() != 3
            || !number_ueb(&tok[1], &mut month)
            || month < 1.0
            || month > 12.0
            || month != month.floor()
            || !number_ueb(&tok[2], &mut value)
            || value <= 0.0
            || value > 100.0
        {
            return "BC_MONTHLY requires month 1–12 and a positive temperature range (°C).".as_bytes().to_vec();
        }
        next.monthly_range[(month as i32 - 1) as usize] = Some(value);
    } else {
        if tok.len() != 2 {
            return b"Unexpected UEB option fields.".to_vec();
        }
        let value = upper(&tok[1]);
        if key == b"SURFACE_TEMP" {
            if value == b"EQUILIBRIUM" {
                next.parameters.temperature = TemperatureMethod::Equilibrium;
            } else if value == b"MODIFIED_FORCE_RESTORE" {
                next.parameters.temperature = TemperatureMethod::ModifiedForceRestore;
            } else {
                return b"Unknown UEB surface temperature scheme.".to_vec();
            }
        } else if key == b"PRECIP_SPLIT" {
            if value != b"DIVT" && value != b"UEB" {
                return b"PRECIP_SPLIT must be UEB or DIVT.".to_vec();
            }
            next.division_temperature = value == b"DIVT";
        } else {
            let mut v = 0.0f64;
            if !number_ueb(&tok[1], &mut v) {
                return b"UEB option must be finite and numeric.".to_vec();
            }
            let mut found = false;
            if let Some(slot) = param_slot(&mut next.parameters, &key) {
                *slot = v;
                found = true;
            }
            if key == b"BC_A" {
                next.bc_a = v;
                found = true;
            }
            if key == b"BC_C" {
                next.bc_c = v;
                found = true;
            }
            if !found {
                return cat(&[b"Unknown UEB option: ", &tok[0]]);
            }
        }
    }
    if validate {
        let error = parameters_validate(&next.parameters);
        if !error.is_empty() {
            return error;
        }
    }
    if next.bc_a <= 0.0 || next.bc_a > 1.0 || next.bc_c <= 0.0 {
        return "Invalid Bristow–Campbell parameters.".as_bytes().to_vec();
    }
    *cfg = next;
    Vec::new()
}

/// std::getline(stream, field, ':') over `s`: the fields, without a final empty one after
/// a trailing ':'.
fn getline_fields(s: &[u8]) -> Vec<&[u8]> {
    let mut out: Vec<&[u8]> = s.split(|&c| c == b':').collect();
    if s.last() == Some(&b':') || s.is_empty() {
        out.pop();
    }
    out
}

/// setSnowTemperature: SITE / SNOWMELT_MODEL / SNOW_STEP / MET_FILE / SHORTWAVE / LONGWAVE.
/// `tok` includes the key.
pub fn set_snow_temperature(cfg: &mut SnowConfig, tok: &[Bytes]) -> Bytes {
    if tok.is_empty() {
        return b"Missing snow/climate key.".to_vec();
    }
    let key = upper(&tok[0]);
    let mut next = cfg.clone();
    if key == b"SITE" {
        let mut site = ClimateSite::default();
        if tok.len() != 5
            || !number_ueb(&tok[1], &mut site.latitude)
            || !number_ueb(&tok[2], &mut site.longitude)
            || !number_ueb(&tok[3], &mut site.timezone)
            || !number_ueb(&tok[4], &mut site.elevation)
            || site.latitude.abs() > 90.0
            || site.longitude.abs() > 180.0
            || site.timezone.abs() > 14.0
            || site.elevation < -500.0
            || site.elevation > 10000.0
        {
            return b"SITE requires latitude, longitude east, UTC offset and elevation in metres.".to_vec();
        }
        next.site = Some(site);
    } else if key == b"SNOWMELT_MODEL" {
        if tok.len() != 2 || !model(&tok[1], &mut next.model) {
            return b"SNOWMELT_MODEL must be DEGREE_DAY or UEB.".to_vec();
        }
        next.model_authored = true;
    } else if key == b"SNOW_STEP" {
        if tok.len() != 2 {
            return b"SNOW_STEP requires seconds or hh:mm[:ss].".to_vec();
        }
        // Strict clock parsing; the legacy permissive parser remains unchanged.
        let s: &[u8] = &tok[1];
        let mut seconds = 0.0f64;
        if !s.contains(&b':') {
            if !number_ueb(s, &mut seconds) {
                return b"Invalid SNOW_STEP.".to_vec();
            }
        } else {
            let mut fields: Vec<f64> = Vec::new();
            for field in getline_fields(s) {
                let mut v = 0.0f64;
                if !number_ueb(field, &mut v) {
                    return b"Invalid SNOW_STEP.".to_vec();
                }
                fields.push(v);
            }
            if s.last() == Some(&b':')
                || fields.len() < 2
                || fields.len() > 3
                || fields[0] < 0.0
                || fields[0] != fields[0].floor()
                || fields[1] < 0.0
                || fields[1] >= 60.0
                || fields[1] != fields[1].floor()
                || (fields.len() == 3 && (fields[2] < 0.0 || fields[2] >= 60.0))
            {
                return b"Invalid SNOW_STEP.".to_vec();
            }
            seconds = fields[0] * 3600.0 + fields[1] * 60.0 + (if fields.len() == 3 { fields[2] } else { 0.0 });
        }
        if seconds < 1.0 || seconds > 21600.0 {
            return b"SNOW_STEP must be between one second and six hours.".to_vec();
        }
        next.step = seconds;
        next.parameters.max_step = seconds;
    } else if key == b"MET_FILE" {
        if tok.len() != 2 || tok[1].is_empty() {
            return b"MET_FILE requires a path.".to_vec();
        }
        next.met_file.set(&tok[1]);
    } else if key == b"SHORTWAVE" || key == b"LONGWAVE" {
        if tok.len() < 2 {
            return b"Radiation source is missing.".to_vec();
        }
        let source = if key == b"SHORTWAVE" { &mut next.shortwave } else { &mut next.longwave };
        let mode = upper(&tok[1]);
        if mode == b"TIMESERIES" && tok.len() == 3 && !tok[2].is_empty() {
            source.mode = RadiationSourceMode::Timeseries;
            source.series = tok[2].clone();
        } else if (mode == b"COMPUTED" || mode == b"FILE") && tok.len() == 2 {
            source.mode = if mode == b"FILE" { RadiationSourceMode::File } else { RadiationSourceMode::Computed };
            source.series.clear();
        } else {
            return b"Radiation source must be COMPUTED, FILE or TIMESERIES name.".to_vec();
        }
        source.authored = true;
    } else {
        return cat(&[b"Unknown snow/climate option: ", &tok[0]]);
    }
    *cfg = next;
    Vec::new()
}

/// setSnowpackUeb: `tok` starts with Fc F Qg aep [U0] [MODEL kind] (the pack name and the
/// UEB keyword already removed).
pub fn set_snowpack_ueb(cfg: &mut SnowConfig, pack: i32, tok: &[Bytes]) -> Bytes {
    if pack < 0 || tok.len() < 4 {
        return b"UEB snowpack requires Fc F Qg aep [U0] [MODEL kind].".to_vec();
    }
    let mut next = SnowpackConfig::default();
    if !number_ueb(&tok[0], &mut next.site.forest)
        || !number_ueb(&tok[1], &mut next.site.drift)
        || !number_ueb(&tok[2], &mut next.site.ground_heat)
        || !number_ueb(&tok[3], &mut next.site.extinction_depth)
    {
        return b"Non-numeric UEB snowpack site value.".to_vec();
    }
    let mut at = 4usize;
    if at < tok.len() && upper(&tok[at]) != b"MODEL" {
        let ok = number_ueb(&tok[at], &mut next.initial_energy);
        at += 1;
        if !ok {
            return b"Invalid UEB initial energy.".to_vec();
        }
    }
    if at < tok.len() {
        if at + 2 != tok.len() || upper(&tok[at]) != b"MODEL" || !model(&tok[at + 1], &mut next.model) {
            return b"Invalid UEB snowpack model override.".to_vec();
        }
    }
    if next.site.forest < 0.0 || next.site.forest > 1.0 || next.site.drift < 0.0 || next.site.extinction_depth <= 0.0 {
        return b"Invalid UEB snowpack site range.".to_vec();
    }
    cfg.packs.insert(pack, next);
    Vec::new()
}

// ===========================================================================
// MeshSnowConfig.cpp
// ===========================================================================

/// MeshSnowConfig.cpp `number`: strictly parsed and finite (an empty token reads as 0).
fn number_mesh(s: &[u8], v: &mut f64) -> bool {
    parse_double_strict(s, v) && v.is_finite()
}

fn has_cr_lf_quote(s: &[u8]) -> bool {
    s.iter().any(|&c| c == b'\r' || c == b'\n' || c == b'"')
}

/// setMeshSnowOption: SNOWMELT AUTO|YES|NO, or SNOW_STEP DEFAULT|<strict step>.
pub fn set_mesh_snow_option(config: &mut MeshSnowConfig, tokens: &[Bytes]) -> Bytes {
    if tokens.len() != 2 {
        return b"Mesh snow option requires one value.".to_vec();
    }
    if upper(&tokens[0]) == b"SNOWMELT" {
        let value = upper(&tokens[1]);
        if value == b"AUTO" {
            config.mode = MeshSnowConfigMode::Auto;
        } else if value == b"YES" {
            config.mode = MeshSnowConfigMode::Yes;
        } else if value == b"NO" {
            config.mode = MeshSnowConfigMode::No;
        } else {
            return b"Mesh SNOWMELT must be AUTO, YES or NO.".to_vec();
        }
    } else if upper(&tokens[0]) == b"SNOW_STEP" {
        if upper(&tokens[1]) == b"DEFAULT" {
            config.step = 0.0;
        } else {
            let mut parsed = SnowConfig::default();
            let error = set_snow_temperature(&mut parsed, tokens);
            if !error.is_empty() {
                return error;
            }
            config.step = parsed.step;
        }
    } else {
        return cat(&[b"Unknown mesh snow option: ", &tokens[0]]);
    }
    Vec::new()
}

/// setMeshSnowRow: [CELL|TAG] key SNOWPACK name [FC v] [F v] [QG v] [AEP v] [SLOPE v ASPECT v].
pub fn set_mesh_snow_row(config: &mut MeshSnowConfig, tokens: &[Bytes], cell: bool) -> Bytes {
    let mut at = 0usize;
    if !tokens.is_empty() && upper(&tokens[0]) == (if cell { &b"CELL"[..] } else { &b"TAG"[..] }) {
        at += 1;
    }
    if tokens.len() < at + 3 {
        return b"Mesh snow row requires a cell/tag and SNOWPACK name.".to_vec();
    }
    let mut key: i32 = -1;
    let tag: Bytes = tokens[at].clone();
    if cell {
        let mut value = 0.0f64;
        if !number_mesh(&tokens[at], &mut value) || value != value.floor() || value < 1.0 || value > i32::MAX as f64 {
            return b"Mesh snow cell must be a positive integer.".to_vec();
        }
        key = value as i32 - 1;
    } else if tag.is_empty() || has_cr_lf_quote(&tag) {
        return b"Invalid mesh snow tag.".to_vec();
    }
    at += 1;
    if upper(&tokens[at]) != b"SNOWPACK" {
        return b"Expected SNOWPACK in mesh snow assignment.".to_vec();
    }
    let mut row = MeshSnowRow::default();
    at += 1;
    row.pack = tokens[at].clone();
    at += 1;
    if row.pack.is_empty() || has_cr_lf_quote(&row.pack) {
        return b"Invalid mesh snowpack name.".to_vec();
    }
    if upper(&row.pack) == b"NONE" {
        row.pack = b"NONE".to_vec();
    }
    let mut seen: BTreeSet<Bytes> = BTreeSet::new();
    while at < tokens.len() {
        let name = upper(&tokens[at]);
        at += 1;
        let mut value = 0.0f64;
        if at >= tokens.len() {
            return b"Mesh snow overrides require unique keys and finite values.".to_vec();
        }
        let ok = number_mesh(&tokens[at], &mut value);
        at += 1;
        if !ok || !seen.insert(name.clone()) {
            return b"Mesh snow overrides require unique keys and finite values.".to_vec();
        }
        if name == b"FC" && value >= 0.0 && value <= 1.0 {
            row.forest = Some(value);
        } else if name == b"F" && value >= 0.0 {
            row.drift = Some(value);
        } else if name == b"QG" {
            row.ground_heat = Some(value);
        } else if name == b"AEP" && value > 0.0 {
            row.extinction_depth = Some(value);
        } else if name == b"SLOPE" && value >= 0.0 && value < 90.0 {
            row.slope = Some(value);
        } else if name == b"ASPECT" && value >= 0.0 && value <= 360.0 {
            row.aspect = Some(value);
        } else {
            return cat(&[b"Invalid mesh snow override: ", &name]);
        }
    }
    if row.slope.is_some() != row.aspect.is_some() {
        return b"Supply both SLOPE and ASPECT, or use cell geometry.".to_vec();
    }
    if cell {
        config.cells.insert(key, row);
    } else {
        config.defaults.insert(tag, row);
    }
    Vec::new()
}
