// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/components/HeatModule/HeatComponent.hpp
// and HeatComponent.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The `org.hydrocouple.openswmm.heat` component's config apply hook (model.heat):
// [HEAT_SOURCES], [HEAT_FLUXES], [RADIATIVE_FLUXES], [SOLAR_RADIATION], [CLOUD_COVER] and
// [SEDIMENT_EXCHANGE], with the per-element override rows (PE2) that
// heat_overrides::resolve_heat_overrides resolves at open. The serializer
// (saveHeatConfig) is not part of open() and is not translated here.
//
// Differences of form only: the `ranged` / `nonneg` / `positive` lambdas take the error
// list and destination as arguments; std::to_string(double) (printf "%f") is
// format!("{:.6}") of the same constants.

use crate::core::context::{ProcessComponentSpec, SimulationContext};
use crate::core::cstdlib::{c_str, strtod};
use crate::core::options::Bytes;
use crate::core::process_components::{ComponentConfigSections, ProcessComponentRegistry};
use crate::data::heat::{DryTempPolicy, HeatConfigData, HeatSource, ShortwaveMode};
use crate::data::heat_override::{is_sediment_attr, HeatAttr, HeatOverrideRow, HeatScope};
use crate::transport::ard_config::{tokenize_ws as tokenize, upper};

const K_HEAT_ID: &str = "org.hydrocouple.openswmm.heat";

const K_MIN_TEMP: f64 = -50.0;
const K_MAX_TEMP: f64 = 100.0;

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

fn parse_finite(tok: &[u8], out: &mut f64) -> bool {
    let (v, n) = strtod(tok);
    *out = v;
    if n != c_str(tok).len() || n == 0 {
        return false;
    }
    out.is_finite()
}

fn parse_celsius(tok: &[u8], out: &mut f64) -> bool {
    let (v, n) = strtod(tok);
    *out = v;
    if n != c_str(tok).len() || n == 0 {
        return false;
    }
    *out >= K_MIN_TEMP && *out <= K_MAX_TEMP
}

fn source_of(u: &[u8]) -> i32 {
    match u {
        b"RAINFALL" => HeatSource::Rainfall as i32,
        b"DWF" => HeatSource::Dwf as i32,
        b"GW" => HeatSource::Gw as i32,
        b"RDII" => HeatSource::Rdii as i32,
        b"EXTERNAL_INFLOW" => HeatSource::ExternalInflow as i32,
        b"IFACE" => HeatSource::Iface as i32,
        b"INITIAL_STATE" => HeatSource::InitialState as i32,
        _ => -1,
    }
}

struct AttrSpec {
    key: &'static str,
    attr: HeatAttr,
    /// refused outside [0,1]
    is_fraction: bool,
    /// refused at or below zero
    positive_only: bool,
}

const K_ATTR_TABLE: [AttrSpec; 16] = [
    AttrSpec { key: "ALBEDO", attr: HeatAttr::Albedo, is_fraction: true, positive_only: false },
    AttrSpec { key: "SHADE_FACTOR", attr: HeatAttr::ShadeFactor, is_fraction: true, positive_only: false },
    AttrSpec { key: "SKY_VIEW", attr: HeatAttr::SkyView, is_fraction: true, positive_only: false },
    AttrSpec { key: "EMISS_WATER", attr: HeatAttr::EmissWater, is_fraction: true, positive_only: false },
    AttrSpec { key: "EMISS_LANDCOVER", attr: HeatAttr::EmissLandcover, is_fraction: true, positive_only: false },
    AttrSpec { key: "LANDCOVER_TEMPERATURE", attr: HeatAttr::LandcoverTemp, is_fraction: false, positive_only: false },
    AttrSpec { key: "THERMAL_DIFFUSIVITY", attr: HeatAttr::SedThermalDiffusivity, is_fraction: false, positive_only: true },
    AttrSpec { key: "SOLUTE_DIFFUSIVITY", attr: HeatAttr::SedSoluteDiffusivity, is_fraction: false, positive_only: false },
    AttrSpec { key: "BED_THICKNESS", attr: HeatAttr::SedBedThickness, is_fraction: false, positive_only: true },
    AttrSpec { key: "GROUND_DEPTH", attr: HeatAttr::SedGroundDepth, is_fraction: false, positive_only: true },
    AttrSpec { key: "GROUND_TEMPERATURE", attr: HeatAttr::SedGroundTemp, is_fraction: false, positive_only: false },
    AttrSpec { key: "HYPORHEIC_VELOCITY", attr: HeatAttr::SedHyporheicVelocity, is_fraction: false, positive_only: false },
    AttrSpec { key: "SEDIMENT_DENSITY", attr: HeatAttr::SedDensity, is_fraction: false, positive_only: true },
    AttrSpec { key: "SEDIMENT_SPECIFIC_HEAT", attr: HeatAttr::SedSpecificHeat, is_fraction: false, positive_only: true },
    AttrSpec { key: "MESH_POROSITY", attr: HeatAttr::SedMeshPorosity, is_fraction: true, positive_only: false },
    AttrSpec { key: "GROUND_SOURCE", attr: HeatAttr::SedGroundSource, is_fraction: false, positive_only: false },
];

fn attr_key(a: HeatAttr) -> &'static str {
    for e in &K_ATTR_TABLE {
        if e.attr == a {
            return e.key;
        }
    }
    "?"
}

struct ScopeRef {
    scope: HeatScope,
    name: Bytes,
    value_at: i32,
    ok: bool,
}

fn parse_scope(toks: &[Bytes], section: &[u8], line: &[u8], errors: &mut Vec<Bytes>) -> ScopeRef {
    let mut r = ScopeRef { scope: HeatScope::Global, name: Vec::new(), value_at: 2, ok: false };
    if toks.len() < 3 {
        errors.push(cat(&[b"[", section, b"] expects '<param> GLOBAL <value>' or '<param> TAG|LINK|NODE <name> <value>': '", line, b"'."]));
        return r;
    }
    let sc = upper(&toks[1]);
    if sc == b"GLOBAL" {
        r.scope = HeatScope::Global;
        r.value_at = 2;
        r.ok = true;
        return r;
    }
    if sc == b"TAG" || sc == b"LINK" || sc == b"NODE" || sc == b"CELL2D" || sc == b"MESH_TAG" {
        if toks.len() < 4 {
            errors.push(cat(&[b"[", section, b"] ", &sc, b" scope needs a name and a value: '", line, b"'."]));
            return r;
        }
        r.scope = if sc == b"TAG" {
            HeatScope::Tag
        } else if sc == b"LINK" {
            HeatScope::Link
        } else if sc == b"NODE" {
            HeatScope::Node
        } else if sc == b"CELL2D" {
            HeatScope::Cell2d
        } else {
            HeatScope::MeshTag
        };
        r.name = toks[2].clone();
        r.value_at = 3;
        r.ok = true;
        return r;
    }
    errors.push(cat(&[b"[", section, b"] unknown scope '", &toks[1], b"' (GLOBAL, TAG, LINK, NODE, CELL2D, MESH_TAG): '", line, b"'."]));
    r
}

fn validate_attr(spec: &AttrSpec, v: f64, section: &[u8], tok: &[u8], errors: &mut Vec<Bytes>) -> bool {
    if spec.is_fraction && (v < 0.0 || v > 1.0) {
        errors.push(cat(&[b"[", section, b"] '", spec.key.as_bytes(), b"' must be a fraction in [0,1], got ", tok, b"."]));
        return false;
    }
    if spec.positive_only && !(v > 0.0) {
        errors.push(cat(&[b"[", section, b"] ", spec.key.as_bytes(), b" must be greater than zero, got ", tok, b"."]));
        return false;
    }
    if spec.attr == HeatAttr::SedHyporheicVelocity && v < 0.0 {
        errors.push(cat(&[
            b"[",
            section,
            b"] HYPORHEIC_VELOCITY cannot be negative: the exchange direction follows the gradient, not this sign. Got ",
            tok,
            b".",
        ]));
        return false;
    }
    if spec.attr == HeatAttr::SedSoluteDiffusivity && v < 0.0 {
        errors.push(cat(&[b"[", section, b"] SOLUTE_DIFFUSIVITY cannot be negative, got ", tok, b"."]));
        return false;
    }
    if spec.attr == HeatAttr::SedGroundSource && v != 0.0 && v != 1.0 {
        errors.push(cat(&[b"[", section, b"] GROUND_SOURCE expects PRESCRIBED or MESH_CONTACTS (0 or 1)."]));
        return false;
    }
    true
}

fn record_override(ctx: &mut SimulationContext, attr: HeatAttr, sr: &ScopeRef, v: f64, section: &[u8], line: &[u8], errors: &mut Vec<Bytes>) -> bool {
    if attr == HeatAttr::SedGroundSource && sr.scope != HeatScope::Tag && sr.scope != HeatScope::Link {
        errors.push(cat(&[b"[", section, b"] GROUND_SOURCE applies only to GLOBAL, network TAG or LINK."]));
        return false;
    }
    let rows = &mut ctx.heat_config.overrides.rows;
    for r in rows.iter() {
        if r.attr == attr && r.scope == sr.scope && r.name == sr.name {
            errors.push(cat(&[
                b"[",
                section,
                b"] duplicate ",
                attr_key(attr).as_bytes(),
                b" row at the same scope for '",
                &sr.name,
                "' \u{2014} one of the two would silently win: '".as_bytes(),
                line,
                b"'.",
            ]));
            return false;
        }
    }
    rows.push(HeatOverrideRow { attr, scope: sr.scope, name: sr.name.clone(), value: v });
    true
}

fn apply_heat_sections(ctx: &mut SimulationContext, _spec: &ProcessComponentSpec, config: &ComponentConfigSections, errors: &mut Vec<Bytes>) {
    ctx.heat_config = HeatConfigData::default(); // reopen hygiene
    let errors_before = errors.len();

    // H6a (D-H6a-3): the three SHORTWAVE spellings are mutually exclusive.
    let mut shortwave_seen = false;
    // Was an actual FRACTION row given, as opposed to coefficients alone?
    let mut cloud_fraction_seen = false;
    // Warnings merge into ctx only if the parse succeeds.
    let mut pending_warnings: Vec<Bytes> = Vec::new();

    for sec in &config.sections {
        let tag: &[u8] = &sec.0;
        if tag == b"RADIATIVE_FLUXES" {
            for line in &sec.1 {
                let toks = tokenize(line);
                let sr = parse_scope(&toks, b"RADIATIVE_FLUXES", line, errors);
                if !sr.ok {
                    continue;
                }
                let k = upper(&toks[0]);

                // PE2 §2.2: SHORTWAVE and the atmosphere keys stay GLOBAL.
                if sr.scope != HeatScope::Global && (k == b"SHORTWAVE" || k == b"ATM_EMISS_COEFF" || k == b"ATM_LW_REFLECTION") {
                    errors.push(cat(&[
                        b"[RADIATIVE_FLUXES] ",
                        &toks[0],
                        " is GLOBAL only \u{2014} it describes the incident resource or the atmosphere, not the reach. Per-element shading is SHADE_FACTOR / SKY_VIEW: '".as_bytes(),
                        line,
                        b"'.",
                    ]));
                    continue;
                }

                // H6a: SHORTWAVE before the numeric parse (two of its spellings are keywords).
                if k == b"SHORTWAVE" {
                    if shortwave_seen {
                        errors.push(
                            "[RADIATIVE_FLUXES] SHORTWAVE is given more than once. The constant, TIMESERIES and COMPUTED spellings are mutually exclusive \u{2014} pick one (plan D-H6a-3)."
                                .as_bytes()
                                .to_vec(),
                        );
                        continue;
                    }
                    let spelling = upper(&toks[2]);

                    if spelling == b"COMPUTED" {
                        if toks.len() != 3 {
                            errors.push(cat(&[b"[RADIATIVE_FLUXES] SHORTWAVE GLOBAL COMPUTED takes no further tokens: '", line, b"'."]));
                            continue;
                        }
                        ctx.heat_config.radiative.sw_mode = ShortwaveMode::Computed;
                        shortwave_seen = true;
                        continue;
                    }

                    if spelling == b"TIMESERIES" {
                        if toks.len() != 4 {
                            errors.push(cat(&[b"[RADIATIVE_FLUXES] SHORTWAVE GLOBAL TIMESERIES takes a series name: '", line, b"'."]));
                            continue;
                        }
                        let ts = ctx.find_timeseries(&toks[3]);
                        if ts < 0 {
                            errors.push(cat(&[b"[RADIATIVE_FLUXES] SHORTWAVE TIMESERIES '", &toks[3], b"' is not a [TIMESERIES] in the model."]));
                            continue;
                        }
                        ctx.heat_config.radiative.sw_mode = ShortwaveMode::Timeseries;
                        ctx.heat_config.radiative.sw_ts_index = ts;
                        shortwave_seen = true;
                        continue;
                    }

                    // Constant — the H3 spelling, still the default.
                    if toks.len() != 3 {
                        errors.push(cat(&[b"[RADIATIVE_FLUXES] malformed SHORTWAVE row: '", line, b"'."]));
                        continue;
                    }
                    let mut sw = 0.0;
                    if !parse_finite(&toks[2], &mut sw) {
                        errors.push(cat(&[b"[RADIATIVE_FLUXES] SHORTWAVE '", &toks[2], b"' is not a finite number, TIMESERIES, or COMPUTED."]));
                        continue;
                    }
                    if sw < 0.0 {
                        errors.push(b"[RADIATIVE_FLUXES] SHORTWAVE must be >= 0 W/m2.".to_vec());
                        continue;
                    }
                    ctx.heat_config.radiative.shortwave_wm2 = sw;
                    ctx.heat_config.radiative.sw_mode = ShortwaveMode::Constant;
                    shortwave_seen = true;
                    continue;
                }

                let vt = sr.value_at as usize;
                let mut v = 0.0;
                if vt >= toks.len() || !parse_finite(&toks[vt], &mut v) {
                    errors.push(cat(&[b"[RADIATIVE_FLUXES] '", &toks[0], b"': missing or non-finite value in '", line, b"'."]));
                    continue;
                }
                if k == b"ATM_EMISS_COEFF" || k == b"ATM_LW_REFLECTION" {
                    if v < 0.0 || v > 1.0 {
                        errors.push(cat(&[b"[RADIATIVE_FLUXES] '", &toks[0], b"' must be a fraction in [0,1], got ", &toks[vt], b"."]));
                        continue;
                    }
                    if k == b"ATM_EMISS_COEFF" {
                        ctx.heat_config.radiative.atm_emiss_coeff = v;
                    } else {
                        ctx.heat_config.radiative.lw_reflection = v;
                    }
                    continue;
                }

                let mut spec: Option<&AttrSpec> = None;
                for e in &K_ATTR_TABLE {
                    if k == e.key.as_bytes() && !is_sediment_attr(e.attr) {
                        spec = Some(e);
                    }
                }
                let spec = match spec {
                    Some(s) => s,
                    None => {
                        errors.push(cat(&[b"[RADIATIVE_FLUXES] unknown parameter '", &toks[0], b"'."]));
                        continue;
                    }
                };
                if !validate_attr(spec, v, b"RADIATIVE_FLUXES", &toks[vt], errors) {
                    continue;
                }

                if sr.scope == HeatScope::Global {
                    let rc = &mut ctx.heat_config.radiative;
                    match spec.attr {
                        HeatAttr::Albedo => rc.albedo = v,
                        HeatAttr::ShadeFactor => rc.shade_factor = v,
                        HeatAttr::SkyView => rc.sky_view = v,
                        HeatAttr::EmissWater => rc.emiss_water = v,
                        HeatAttr::EmissLandcover => rc.emiss_landcover = v,
                        HeatAttr::LandcoverTemp => rc.landcover_temp = v,
                        _ => {}
                    }
                } else {
                    record_override(ctx, spec.attr, &sr, v, b"RADIATIVE_FLUXES", line, errors);
                }
            }
            continue;
        }
        // H6a. Site geometry and Bird atmosphere, consulted only under SHORTWAVE COMPUTED.
        if tag == b"SOLAR_RADIATION" {
            for line in &sec.1 {
                let toks = tokenize(line);
                if toks.len() != 3 || upper(&toks[1]) != b"GLOBAL" {
                    errors.push(cat(&[b"[SOLAR_RADIATION] expects '<param> GLOBAL <value>': '", line, b"'."]));
                    continue;
                }
                let mut v = 0.0;
                if !parse_finite(&toks[2], &mut v) {
                    errors.push(cat(&[b"[SOLAR_RADIATION] '", &toks[0], b"': '", &toks[2], b"' is not a finite number."]));
                    continue;
                }
                let k = upper(&toks[0]);

                // `ranged` returns whether it assigned; the "provided" flags follow it.
                let ranged = |errors: &mut Vec<Bytes>, lo: f64, hi: f64, dst: &mut f64, units: &[u8]| -> bool {
                    if v < lo || v > hi {
                        errors.push(cat(&[
                            b"[SOLAR_RADIATION] '",
                            &toks[0],
                            b"' must be in [",
                            format!("{:.6}", lo).as_bytes(),
                            b", ",
                            format!("{:.6}", hi).as_bytes(),
                            b"] ",
                            units,
                            b", got ",
                            &toks[2],
                            b".",
                        ]));
                        return false;
                    }
                    *dst = v;
                    true
                };
                let nonneg = |errors: &mut Vec<Bytes>, dst: &mut f64| {
                    if v < 0.0 {
                        errors.push(cat(&[b"[SOLAR_RADIATION] '", &toks[0], b"' must be >= 0, got ", &toks[2], b"."]));
                    } else {
                        *dst = v;
                    }
                };

                let sc = &mut ctx.heat_config.solar;
                if k == b"LATITUDE" {
                    if ranged(errors, -90.0, 90.0, &mut sc.latitude_deg, b"degrees") {
                        sc.has_latitude = true;
                    }
                } else if k == b"LONGITUDE" {
                    if ranged(errors, -180.0, 180.0, &mut sc.longitude_deg, b"degrees") {
                        sc.has_longitude = true;
                    }
                } else if k == b"TIMEZONE" {
                    // ±14 covers every real offset (Kiribati is +14).
                    if ranged(errors, -14.0, 14.0, &mut sc.timezone_hours, b"hours from UTC") {
                        sc.has_timezone = true;
                    }
                } else if k == b"ELEVATION" {
                    if ranged(errors, -500.0, 9000.0, &mut sc.elevation_m, b"metres") {
                        sc.has_elevation = true;
                    }
                } else if k == b"TURBIDITY_380" {
                    nonneg(errors, &mut sc.aod380);
                } else if k == b"TURBIDITY_500" {
                    nonneg(errors, &mut sc.aod500);
                } else if k == b"PRECIP_WATER" {
                    nonneg(errors, &mut sc.precip_water_cm);
                } else if k == b"OZONE" {
                    nonneg(errors, &mut sc.ozone_cm);
                } else if k == b"GROUND_ALBEDO" {
                    if v < 0.0 || v > 1.0 {
                        errors.push(cat(&[b"[SOLAR_RADIATION] GROUND_ALBEDO must be a fraction in [0,1], got ", &toks[2], b"."]));
                    } else {
                        sc.ground_albedo = v;
                    }
                } else {
                    errors.push(cat(&[
                        b"[SOLAR_RADIATION] unknown parameter '",
                        &toks[0],
                        b"' (LATITUDE, LONGITUDE, TIMEZONE, ELEVATION, TURBIDITY_380, TURBIDITY_500, PRECIP_WATER, OZONE, GROUND_ALBEDO).",
                    ]));
                }
            }
            continue;
        }

        // H6b — the bed / hyporheic zone.
        if tag == b"SEDIMENT_EXCHANGE" {
            for line in &sec.1 {
                let toks = tokenize(line);
                let sr = parse_scope(&toks, b"SEDIMENT_EXCHANGE", line, errors);
                if !sr.ok {
                    continue;
                }
                let k = upper(&toks[0]);

                // The bed zone is CONDUITS ONLY.
                if sr.scope == HeatScope::Node {
                    errors.push(cat(&[
                        "[SEDIMENT_EXCHANGE] NODE scope is not available \u{2014} the bed zone exists beneath CONDUITS only. Use LINK or TAG: '".as_bytes(),
                        line,
                        b"'.",
                    ]));
                    continue;
                }

                // GROUND_TEMPERATURE is the one key with a TIMESERIES spelling.
                if k == b"GROUND_TEMPERATURE" && sr.scope == HeatScope::Global && upper(&toks[2]) == b"TIMESERIES" {
                    if toks.len() != 4 {
                        errors.push(cat(&[b"[SEDIMENT_EXCHANGE] GROUND_TEMPERATURE GLOBAL TIMESERIES takes a series name: '", line, b"'."]));
                        continue;
                    }
                    let ts = ctx.find_timeseries(&toks[3]);
                    if ts < 0 {
                        errors.push(cat(&[b"[SEDIMENT_EXCHANGE] TIMESERIES '", &toks[3], b"' is not a [TIMESERIES] in the model."]));
                        continue;
                    }
                    ctx.heat_config.sediment.ground_ts_index = ts;
                    ctx.heat_config.sediment.has_ground_temp = true;
                    continue;
                }

                let vt = sr.value_at as usize;
                let mut v = 0.0;
                if k == b"GROUND_SOURCE" {
                    if vt >= toks.len() || toks.len() != vt + 1 {
                        errors.push(b"[SEDIMENT_EXCHANGE] GROUND_SOURCE needs one PRESCRIBED or MESH_CONTACTS value.".to_vec());
                        continue;
                    }
                    let mode = upper(&toks[vt]);
                    if mode == b"PRESCRIBED" {
                        v = 0.0;
                    } else if mode == b"MESH_CONTACTS" {
                        v = 1.0;
                    } else if !parse_finite(&toks[vt], &mut v) || (v != 0.0 && v != 1.0) {
                        errors.push(b"[SEDIMENT_EXCHANGE] GROUND_SOURCE expects PRESCRIBED or MESH_CONTACTS.".to_vec());
                        continue;
                    }
                    if sr.scope == HeatScope::Global {
                        ctx.heat_config.sediment.mesh_ground = v != 0.0;
                        ctx.heat_config.sediment.ground_source_configured = true;
                    } else {
                        record_override(ctx, HeatAttr::SedGroundSource, &sr, v, b"SEDIMENT_EXCHANGE", line, errors);
                    }
                    continue;
                }
                if vt >= toks.len() || !parse_finite(&toks[vt], &mut v) {
                    errors.push(cat(&[b"[SEDIMENT_EXCHANGE] ", &toks[0], b": missing or non-finite value in '", line, b"'."]));
                    continue;
                }

                // PE2: per-element rows go through the shared validator.
                if sr.scope != HeatScope::Global {
                    let mut spec: Option<&AttrSpec> = None;
                    for e in &K_ATTR_TABLE {
                        if k == e.key.as_bytes() && is_sediment_attr(e.attr) {
                            spec = Some(e);
                        }
                    }
                    let spec = match spec {
                        Some(s) => s,
                        None => {
                            errors.push(cat(&[
                                b"[SEDIMENT_EXCHANGE] '",
                                &toks[0],
                                b"' has no per-element form (INITIAL_TEMPERATURE is the bed's own seed and is GLOBAL).",
                            ]));
                            continue;
                        }
                    };
                    if validate_attr(spec, v, b"SEDIMENT_EXCHANGE", &toks[vt], errors) {
                        record_override(ctx, spec.attr, &sr, v, b"SEDIMENT_EXCHANGE", line, errors);
                    }
                    continue;
                }

                // Lengths, diffusivities and material properties are refused at zero or below.
                let positive = |errors: &mut Vec<Bytes>, dst: &mut f64| {
                    if !(v > 0.0) {
                        errors.push(cat(&[b"[SEDIMENT_EXCHANGE] ", &toks[0], b" must be greater than zero, got ", &toks[vt], b"."]));
                        return;
                    }
                    *dst = v;
                };

                let sd = &mut ctx.heat_config.sediment;
                if k == b"THERMAL_DIFFUSIVITY" {
                    positive(errors, &mut sd.thermal_diffusivity);
                } else if k == b"SOLUTE_DIFFUSIVITY" {
                    if v < 0.0 {
                        errors.push(cat(&[b"[SEDIMENT_EXCHANGE] SOLUTE_DIFFUSIVITY cannot be negative, got ", &toks[vt], b"."]));
                    } else {
                        sd.solute_diffusivity = v;
                    }
                } else if k == b"BED_THICKNESS" {
                    positive(errors, &mut sd.bed_thickness);
                } else if k == b"GROUND_DEPTH" {
                    positive(errors, &mut sd.ground_depth);
                } else if k == b"SEDIMENT_DENSITY" {
                    positive(errors, &mut sd.sed_density);
                } else if k == b"SEDIMENT_SPECIFIC_HEAT" {
                    positive(errors, &mut sd.sed_specific_heat);
                } else if k == b"MESH_POROSITY" {
                    if v < 0.0 || v > 1.0 {
                        errors.push(b"[SEDIMENT_EXCHANGE] MESH_POROSITY must be in [0,1].".to_vec());
                    } else {
                        sd.mesh_porosity = v;
                    }
                } else if k == b"HYPORHEIC_VELOCITY" {
                    if v < 0.0 {
                        errors.push(cat(&[
                            b"[SEDIMENT_EXCHANGE] HYPORHEIC_VELOCITY cannot be negative: the exchange direction follows the gradient, not this sign. Got ",
                            &toks[vt],
                            b".",
                        ]));
                    } else {
                        sd.hyporheic_velocity = v;
                    }
                } else if k == b"GROUND_TEMPERATURE" {
                    sd.ground_temp = v;
                    sd.ground_ts_index = -1;
                    sd.has_ground_temp = true;
                } else if k == b"INITIAL_TEMPERATURE" {
                    sd.initial_temp = v;
                    sd.has_initial_temp = true;
                } else {
                    errors.push(cat(&[
                        b"[SEDIMENT_EXCHANGE] unknown parameter '",
                        &toks[0],
                        b"' (THERMAL_DIFFUSIVITY, SOLUTE_DIFFUSIVITY, BED_THICKNESS, GROUND_DEPTH, GROUND_TEMPERATURE, HYPORHEIC_VELOCITY, SEDIMENT_DENSITY, SEDIMENT_SPECIFIC_HEAT, INITIAL_TEMPERATURE).",
                    ]));
                }
            }
            continue;
        }

        if tag == b"CLOUD_COVER" {
            for line in &sec.1 {
                let toks = tokenize(line);
                if toks.len() < 3 || upper(&toks[1]) != b"GLOBAL" {
                    errors.push(cat(&[b"[CLOUD_COVER] expects '<param> GLOBAL <value>': '", line, b"'."]));
                    continue;
                }
                let k = upper(&toks[0]);

                if k == b"FRACTION" {
                    if upper(&toks[2]) == b"TIMESERIES" {
                        if toks.len() != 4 {
                            errors.push(cat(&[b"[CLOUD_COVER] FRACTION GLOBAL TIMESERIES takes a series name: '", line, b"'."]));
                            continue;
                        }
                        let ts = ctx.find_timeseries(&toks[3]);
                        if ts < 0 {
                            errors.push(cat(&[b"[CLOUD_COVER] TIMESERIES '", &toks[3], b"' is not a [TIMESERIES] in the model."]));
                            continue;
                        }
                        let cc = &mut ctx.heat_config.cloud;
                        cc.use_timeseries = true;
                        cc.ts_index = ts;
                        cc.configured = true;
                        cloud_fraction_seen = true;
                        continue;
                    }
                    if toks.len() != 3 {
                        errors.push(cat(&[b"[CLOUD_COVER] malformed FRACTION row: '", line, b"'."]));
                        continue;
                    }
                    let mut c = 0.0;
                    if !parse_finite(&toks[2], &mut c) {
                        errors.push(cat(&[b"[CLOUD_COVER] FRACTION '", &toks[2], b"' is not a finite number or TIMESERIES."]));
                        continue;
                    }
                    if c < 0.0 || c > 1.0 {
                        errors.push(cat(&[b"[CLOUD_COVER] FRACTION must be in [0,1], got ", &toks[2], b" (it is a fraction, not a percent)."]));
                        continue;
                    }
                    let cc = &mut ctx.heat_config.cloud;
                    cc.use_timeseries = false;
                    cc.fraction = c;
                    cc.configured = true;
                    cloud_fraction_seen = true;
                    continue;
                }

                if k != b"SW_ATTEN_K" && k != b"SW_ATTEN_N" && k != b"LW_CLOUD_K" {
                    errors.push(cat(&[b"[CLOUD_COVER] unknown parameter '", &toks[0], b"' (FRACTION, SW_ATTEN_K, SW_ATTEN_N, LW_CLOUD_K)."]));
                    continue;
                }
                if toks.len() != 3 {
                    errors.push(cat(&[b"[CLOUD_COVER] expects '<param> GLOBAL <value>': '", line, b"'."]));
                    continue;
                }
                let mut v = 0.0;
                if !parse_finite(&toks[2], &mut v) {
                    errors.push(cat(&[b"[CLOUD_COVER] '", &toks[0], b"': '", &toks[2], b"' is not a finite number."]));
                    continue;
                }
                if v < 0.0 {
                    errors.push(cat(&[b"[CLOUD_COVER] '", &toks[0], b"' must be >= 0, got ", &toks[2], b"."]));
                    continue;
                }
                let cc = &mut ctx.heat_config.cloud;
                if k == b"SW_ATTEN_K" {
                    cc.sw_atten_k = v;
                } else if k == b"SW_ATTEN_N" {
                    cc.sw_atten_n = v;
                } else {
                    cc.lw_cloud_k = v;
                }
                cc.configured = true;
            }
            continue;
        }

        if tag == b"HEAT_FLUXES" {
            for line in &sec.1 {
                let toks = tokenize(line);
                if toks.len() != 2 {
                    errors.push(cat(&[b"[HEAT_FLUXES] expects '<module> ON|OFF': '", line, b"'."]));
                    continue;
                }
                let module = upper(&toks[0]);
                let val = upper(&toks[1]);

                if module == b"DRY_ELEMENT_TEMPERATURE" {
                    if val == b"HOLD" {
                        ctx.heat_config.dry_temp_policy = DryTempPolicy::Hold;
                    } else if val == b"AIR" {
                        ctx.heat_config.dry_temp_policy = DryTempPolicy::Air;
                    } else if val == b"DEFAULT" {
                        ctx.heat_config.dry_temp_policy = DryTempPolicy::Default;
                    } else {
                        errors.push(cat(&[b"[HEAT_FLUXES] DRY_ELEMENT_TEMPERATURE '", &toks[1], b"' is not HOLD, AIR, or DEFAULT."]));
                    }
                    continue;
                }

                if val != b"ON" && val != b"OFF" && val != b"YES" && val != b"NO" {
                    errors.push(cat(&[b"[HEAT_FLUXES] '", &toks[0], b"': '", &toks[1], b"' is not ON or OFF."]));
                    continue;
                }
                let on = val == b"ON" || val == b"YES";
                if module == b"SURFACE_EXCHANGE" {
                    ctx.heat_config.surface_exchange = on;
                } else if module == b"RADIATIVE_EXCHANGE" {
                    ctx.heat_config.radiative_exchange = on;
                } else if module == b"LAYER_CONDUCTION" {
                    ctx.heat_config.layer_conduction = on;
                } else if module == b"SEDIMENT_EXCHANGE" {
                    ctx.heat_config.sediment_exchange = on; // H6b
                } else {
                    errors.push(cat(&[
                        b"[HEAT_FLUXES] unknown module '",
                        &toks[0],
                        b"' (SURFACE_EXCHANGE in H2, RADIATIVE_EXCHANGE in H3, DRY_ELEMENT_TEMPERATURE in H5a, LAYER_CONDUCTION in H5b, SEDIMENT_EXCHANGE in H6b).",
                    ]));
                }
            }
            continue;
        }
        if tag != b"HEAT_SOURCES" {
            errors.push(cat(&[
                b"model.heat: unknown section [",
                tag,
                b"] (recognized: [HEAT_SOURCES], [HEAT_FLUXES], [RADIATIVE_FLUXES], [SOLAR_RADIATION], [CLOUD_COVER], [SEDIMENT_EXCHANGE]).",
            ]));
            continue;
        }
        for line in &sec.1 {
            let toks = tokenize(line);
            if toks.len() < 3 {
                errors.push(cat(&[b"[HEAT_SOURCES] expects '<source> <scope> [name] <degC>': '", line, b"'."]));
                continue;
            }
            let src = source_of(&upper(&toks[0]));
            if src < 0 {
                errors.push(cat(&[b"[HEAT_SOURCES] unknown source '", &toks[0], b"' (RAINFALL, DWF, GW, RDII, EXTERNAL_INFLOW, IFACE, INITIAL_STATE)."]));
                continue;
            }
            let scope = upper(&toks[1]);

            if scope == b"SUBCATCH" {
                errors.push(b"[HEAT_SOURCES] SUBCATCH scope arrives with plan phase H5 (watershed temperature states).".to_vec());
                continue;
            }
            if scope == b"EDGE_BC" {
                errors.push(b"[HEAT_SOURCES] EDGE_BC scope arrives with phase T6 (2D transport).".to_vec());
                continue;
            }
            let has_name = scope == b"NODE";
            let vpos: usize = if has_name { 3 } else { 2 };
            if toks.len() <= vpos {
                errors.push(cat(&[b"[HEAT_SOURCES] malformed row: '", line, b"'."]));
                continue;
            }
            let vtok = &toks[vpos];
            if upper(vtok).starts_with(b"TIMESERIES") {
                errors.push("[HEAT_SOURCES] TIMESERIES temperatures arrive with a later heat phase \u{2014} H1 takes constant degC.".as_bytes().to_vec());
                continue;
            }
            if toks.len() != vpos + 1 {
                errors.push(cat(&[b"[HEAT_SOURCES] malformed row: '", line, b"'."]));
                continue;
            }
            let mut degc = 0.0;
            if !parse_celsius(vtok, &mut degc) {
                errors.push(cat(&[b"[HEAT_SOURCES] '", &toks[0], b"': '", vtok, b"' is not a temperature in degC between -50 and 100."]));
                continue;
            }

            if scope == b"GLOBAL" {
                ctx.heat_config.global_temp[src as usize] = degc;
                ctx.heat_config.configured_source[src as usize] = true;
            } else if scope == b"NODE" {
                if src != HeatSource::Dwf as i32 && src != HeatSource::ExternalInflow as i32 {
                    errors.push(cat(&[b"[HEAT_SOURCES] NODE scope applies to DWF and EXTERNAL_INFLOW in H1; '", &toks[0], b"' takes GLOBAL."]));
                    continue;
                }
                let nd = ctx.node_names.find(&toks[2]);
                if nd < 0 {
                    errors.push(cat(&[b"[HEAT_SOURCES] unknown node '", &toks[2], b"'."]));
                    continue;
                }
                let hc = &mut ctx.heat_config;
                let mut dup = false;
                for i in 0..hc.node_over_source.len() {
                    if hc.node_over_source[i] == src && hc.node_over_node[i] == nd {
                        errors.push(cat(&[b"[HEAT_SOURCES] duplicate NODE row for '", &toks[0], b"' at '", &toks[2], b"'."]));
                        dup = true;
                        break;
                    }
                }
                if dup {
                    continue;
                }
                hc.node_over_source.push(src);
                hc.node_over_node.push(nd);
                hc.node_over_temp.push(degc);
            } else {
                errors.push(cat(&[b"[HEAT_SOURCES] unknown scope '", &toks[1], b"' (GLOBAL or NODE in H1)."]));
            }
        }
    }

    {
        let rc = &ctx.heat_config.radiative;
        let sc = &ctx.heat_config.solar;
        let site_none = ctx.snowpacks.config.site.is_none();

        if rc.sw_mode == ShortwaveMode::Computed {
            if site_none && (!sc.has_latitude || !sc.has_longitude) {
                errors.push(
                    "[RADIATIVE_FLUXES] SHORTWAVE GLOBAL COMPUTED requires [SOLAR_RADIATION] LATITUDE and LONGITUDE. They are NOT taken from the [TEMPERATURE] SNOWMELT line: that latitude defaults to 0 and its longitude field is a solar-time correction in minutes, not a longitude (plan \u{a7}2.5)."
                        .as_bytes()
                        .to_vec(),
                );
            }
        } else if sc.has_latitude || sc.has_longitude {
            pending_warnings.push(b"[SOLAR_RADIATION] is configured but [RADIATIVE_FLUXES] SHORTWAVE is not COMPUTED, so the site coordinates are unused. Set SHORTWAVE GLOBAL COMPUTED to use them.".to_vec());
        }

        if rc.sw_mode == ShortwaveMode::Timeseries && ctx.heat_config.cloud.configured {
            pending_warnings.push(b"[CLOUD_COVER] is applied on top of a measured SHORTWAVE TIMESERIES. A measured record already contains its clouds, so this attenuates them twice. Drop [CLOUD_COVER] unless the series is explicitly clear-sky. (The longwave cloud correction still applies either way.)".to_vec());
        }

        if ctx.heat_config.cloud.configured && !cloud_fraction_seen {
            pending_warnings.push(
                "[CLOUD_COVER] sets coefficients but no FRACTION, so the cloud fraction is 0 and both cloud corrections are the identity \u{2014} the section has no effect. Add 'FRACTION GLOBAL <C>' or 'FRACTION GLOBAL TIMESERIES <name>'."
                    .as_bytes()
                    .to_vec(),
            );
        }

        if rc.sw_mode == ShortwaveMode::Computed && !sc.has_timezone && site_none {
            pending_warnings.push(
                "[SOLAR_RADIATION] has no TIMEZONE, so local time is taken as UTC. If the deck's clock is local, the computed solar day is shifted by the site's offset \u{2014} up to 12 hours, which is larger than every other error in this module. Add 'TIMEZONE GLOBAL <hours from UTC>'."
                    .as_bytes()
                    .to_vec(),
            );
        }
    }

    if errors.len() != errors_before {
        ctx.heat_config = HeatConfigData::default(); // never half-apply
        return; // ...and never half-warn either
    }
    ctx.heat_config.configured = true;

    for w in pending_warnings {
        ctx.warnings.push(w);
    }

    let heat = &ctx.heat_config;
    let thermal = heat.surface_exchange
        || heat.radiative_exchange
        || heat.layer_conduction
        || heat.sediment_exchange
        || heat.dry_temp_policy != DryTempPolicy::Hold
        || !heat.node_over_source.is_empty()
        || !heat.overrides.rows.is_empty()
        || heat.configured_source.iter().any(|&value| value);
    if !ctx.options.heat_transport && thermal {
        ctx.warnings.push(
            "A heat component is configured but [OPTIONS] HEAT_TRANSPORT is OFF \u{2014} no temperature is tracked this simulation. Set HEAT_TRANSPORT ON."
                .as_bytes()
                .to_vec(),
        );
    }
}

pub fn register_heat_component(reg: &mut ProcessComponentRegistry) {
    reg.register_component(
        K_HEAT_ID.as_bytes(),
        b"Heat transport coordinator ([HEAT_SOURCES] inlet temperatures, [HEAT_FLUXES] module toggles, [RADIATIVE_FLUXES] parameters, and H6a's [SOLAR_RADIATION] + [CLOUD_COVER] shortwave forcing)",
        apply_heat_sections,
    );
}
