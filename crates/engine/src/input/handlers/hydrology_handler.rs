// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/handlers/HydrologyHandler.cpp and
// HydrologyHandler.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Section handlers for [EVAPORATION], [TEMPERATURE], [SNOWPACKS], [UEB_OPTIONS],
// [2D_SNOW_OPTIONS], [2D_SNOW_DEFAULTS], [2D_SNOW], [AQUIFERS], [GROUNDWATER], [GWF],
// [LID_CONTROLS] and [LID_USAGE]. The snow-configuration parsers they call are in
// hydrology::snow_config_parse, gwf_validate in hydrology::groundwater_parse.
//
// Differences in form only: the generic `grow` lambdas are written out per column.
// `std::strtod` (NODE LID layer values) is core::cstdlib::strtod, given the token up to its
// first NUL as c_str() would; `end != c_str` is "consumed > 0" and `*end == '\0'` is
// "consumed everything". Assigning a std::string to a FilePathPair is FilePathPair::set;
// constructing one from a string sets `original` only.

use crate::core::cmath::{cmax, cmin};
use crate::core::constants;
use crate::core::context::SimulationContext;
use crate::core::cstdlib::strtod;
use crate::core::error_codes::{format_error, format_error_detail, ERR_DATETIME, ERR_ITEMS, ERR_KEYWORD, ERR_MATH_EXPR, ERR_NUMBER};
use crate::core::options::{Bytes, FilePathPair};
use crate::data::hydrology::{AquiferStore, LidControlStore, SnowpackStore};
use crate::data::lid_node::{LidNodeLayer, LidNodeLayerKind};
use crate::hydrology::groundwater_parse;
use crate::hydrology::snow_config_parse as surface;
use crate::input::parse_utils::{add_unique, parse_date, parse_double_strict, to_double, to_int};
use crate::input::tokenizer;

fn cat(parts: &[&[u8]]) -> Bytes {
    parts.concat()
}

/// Extend snowpack parameter columns through a registered index (zero-filled).
fn ensure_snowpack_capacity(sp: &mut SnowpackStore, idx: i32) {
    let n = (idx + 1) as usize;
    if sp.plowable.len() < n {
        sp.plowable.resize(n, [0.0; 7]);
    }
    if sp.impervious.len() < n {
        sp.impervious.resize(n, [0.0; 7]);
    }
    if sp.pervious.len() < n {
        sp.pervious.resize(n, [0.0; 7]);
    }
    if sp.removal.len() < n {
        sp.removal.resize(n, [0.0; 6]);
    }
    if sp.removal_subcatch.len() < n {
        sp.removal_subcatch.resize(n, Bytes::new());
    }
}

/// Extend aquifer parameter columns through a registered index, with aquifer defaults.
fn ensure_aquifer_capacity(aq: &mut AquiferStore, idx: i32) {
    let n = (idx + 1) as usize;
    fn grow<T: Clone>(v: &mut Vec<T>, n: usize, def: T) {
        if v.len() < n {
            v.resize(n, def);
        }
    }
    grow(&mut aq.porosity, n, 0.5);
    grow(&mut aq.wilting_point, n, 0.15);
    grow(&mut aq.field_capacity, n, 0.30);
    grow(&mut aq.conductivity, n, 0.1);
    grow(&mut aq.conduct_slope, n, 10.0);
    grow(&mut aq.tension_slope, n, 15.0);
    grow(&mut aq.upper_evap, n, 0.35);
    grow(&mut aq.lower_evap, n, 14.0);
    grow(&mut aq.lower_loss, n, 0.002);
    grow(&mut aq.bottom_elev, n, 0.0);
    grow(&mut aq.water_table_elev, n, 0.0);
    grow(&mut aq.upper_moist, n, 0.30);
    grow(&mut aq.upper_evap_pat, n, Bytes::new());
}

/// Extend LID-control parameter columns through a registered index (zeroed layers, empty
/// NODE layer lists and type tokens); removal/name columns are sized elsewhere.
fn ensure_lid_capacity(lc: &mut LidControlStore, idx: i32) {
    let n = (idx + 1) as usize;
    if lc.node_layers.len() < n {
        lc.node_layers.resize(n, Vec::new());
    }
    if lc.lid_type.len() < n {
        lc.lid_type.resize(n, Bytes::new());
    }
    if lc.surface.len() < n {
        lc.surface.resize(n, [0.0; 5]);
    }
    if lc.soil.len() < n {
        lc.soil.resize(n, [0.0; 7]);
    }
    if lc.pavement.len() < n {
        lc.pavement.resize(n, [0.0; 6]);
    }
    if lc.storage.len() < n {
        lc.storage.resize(n, [0.0; 4]);
    }
    if lc.drain.len() < n {
        lc.drain.resize(n, [0.0; 6]);
    }
    if lc.drainmat.len() < n {
        lc.drainmat.resize(n, [0.0; 3]);
    }
}

/// Extend authored subcatchment groundwater columns through a registered index.
fn ensure_subcatch_gw_capacity(ctx: &mut SimulationContext, idx: i32) {
    let n = (idx + 1) as usize;
    fn grow<T: Clone>(v: &mut Vec<T>, n: usize, def: T) {
        if v.len() < n {
            v.resize(n, def);
        }
    }
    let s = &mut ctx.subcatches;
    grow(&mut s.gw_aquifer, n, -1);
    grow(&mut s.gw_node, n, -1);
    grow(&mut s.gw_surf_elev, n, 0.0);
    grow(&mut s.gw_a1, n, 0.0);
    grow(&mut s.gw_b1, n, 0.0);
    grow(&mut s.gw_a2, n, 0.0);
    grow(&mut s.gw_b2, n, 0.0);
    grow(&mut s.gw_a3, n, 0.0);
    grow(&mut s.gw_tw, n, 0.0);
    grow(&mut s.gw_hstar, n, constants::MISSING);
    grow(&mut s.gw_bot_elev, n, constants::MISSING);
    grow(&mut s.gw_wt_elev, n, constants::MISSING);
    grow(&mut s.gw_upper_moist, n, constants::MISSING);
}

// ============================================================================
// handle_evaporation()
// ============================================================================

/// [EVAPORATION]: CONSTANT / MONTHLY / TIMESERIES / TEMPERATURE / FILE / RECOVERY / DRY_ONLY.
pub fn handle_evaporation(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.is_empty() {
            continue;
        }

        let key = tokenizer::to_upper(&tok[0]);

        if key == b"CONSTANT" && tok.len() >= 2 {
            ctx.options.evap_type = 0;
            let val = to_double(&tok[1], 0.0);
            for i in 0..12 {
                ctx.options.evap_values[i] = val;
            }
        } else if key == b"MONTHLY" && tok.len() >= 13 {
            ctx.options.evap_type = 1;
            for i in 0..12 {
                ctx.options.evap_values[i] = to_double(&tok[1 + i], 0.0);
            }
        } else if key == b"TIMESERIES" && tok.len() >= 2 {
            ctx.options.evap_type = 2;
            ctx.options.evap_ts_name = tok[1].clone();
        } else if key == b"TEMPERATURE" {
            ctx.options.evap_type = 3;
        } else if key == b"FILE" {
            ctx.options.evap_type = 4;
            // Optional pan coefficients: FILE pc1 pc2 ... pc12
            if tok.len() >= 13 {
                for i in 0..12 {
                    ctx.options.pan_coeff[i] = to_double(&tok[1 + i], 0.0);
                }
            }
        } else if key == b"RECOVERY" && tok.len() >= 2 {
            ctx.options.evap_recovery_pat = tok[1].clone();
        } else if key == b"DRY_ONLY" && tok.len() >= 2 {
            ctx.options.evap_dry_only = tokenizer::parse_boolean(&tok[1]);
        }
    }
}

// ============================================================================
// handle_temperature()
// ============================================================================

/// [TEMPERATURE]: snow SITE/model/step/met/radiation keys, TIMESERIES, FILE, WINDSPEED,
/// HUMIDITY, SNOWMELT and ADC.
pub fn handle_temperature(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.is_empty() {
            continue;
        }

        let key = tokenizer::to_upper(&tok[0]);

        if key == b"SITE"
            || key == b"SNOWMELT_MODEL"
            || key == b"SNOW_STEP"
            || key == b"MET_FILE"
            || key == b"SHORTWAVE"
            || key == b"LONGWAVE"
        {
            let error = surface::set_snow_temperature(&mut ctx.snowpacks.config, &tok);
            if !error.is_empty() {
                ctx.errors.push(cat(&[b"[TEMPERATURE] ", &error]));
            }
        } else if key == b"TIMESERIES" && tok.len() >= 2 {
            ctx.options.temp_source = 1;
            ctx.options.temp_ts_name = tok[1].clone();
        } else if key == b"FILE" && tok.len() >= 2 {
            ctx.options.temp_source = 2;
            ctx.options.temp_file.set(&tok[1]);
            // Optional start date (legacy climate.c:565-570): a leading `*` means none;
            // anything else must be a date or it is ERROR 213.
            ctx.options.temp_file_start = 0.0;
            if tok.len() >= 3 && !tok[2].is_empty() && tok[2][0] != b'*' {
                let start = parse_date(&tok[2]);
                if start > 0.0 {
                    ctx.options.temp_file_start = start;
                } else {
                    ctx.errors.push(format_error(ERR_DATETIME, &tok[2]));
                }
            }
            // Optional climate-file temperature units keyword: C10, C or F.
            if tok.len() >= 4 {
                let u = tokenizer::to_upper(&tok[3]);
                if u == b"C10" {
                    ctx.options.temp_units = 0;
                } else if u == b"C" {
                    ctx.options.temp_units = 1;
                } else if u == b"F" {
                    ctx.options.temp_units = 2;
                }
            }
        } else if key == b"WINDSPEED" && tok.len() >= 2 {
            let wtype = tokenizer::to_upper(&tok[1]);
            if wtype == b"MONTHLY" && tok.len() >= 14 {
                ctx.options.wind_type = 0;
                for i in 0..12 {
                    ctx.options.wind_speed[i] = to_double(&tok[2 + i], 0.0);
                }
            } else if wtype == b"FILE" {
                ctx.options.wind_type = 1;
            }
        }
        // H2: relative humidity, %:
        //   HUMIDITY [DEWPOINT] <value> | MONTHLY h1 ... h12 | TIMESERIES <name>
        else if key == b"HUMIDITY" && tok.len() >= 2 {
            let mut p = 1usize;
            ctx.options.humidity_var = 0;
            if tokenizer::to_upper(&tok[p]) == b"DEWPOINT" {
                ctx.options.humidity_var = 1;
                p += 1;
            }
            if p >= tok.len() {
                continue;
            }
            let htype = tokenizer::to_upper(&tok[p]);
            if htype == b"MONTHLY" && tok.len() >= p + 13 {
                ctx.options.humidity_type = 1;
                for i in 0..12 {
                    ctx.options.humidity[i] = to_double(&tok[p + 1 + i], 0.0);
                }
            } else if htype == b"TIMESERIES" && tok.len() >= p + 2 {
                ctx.options.humidity_type = 2;
                ctx.options.humidity_ts_name = tok[p + 1].clone();
            } else if htype != b"MONTHLY" && htype != b"TIMESERIES" {
                ctx.options.humidity_type = 0;
                let h = to_double(&tok[p], 0.0);
                for i in 0..12 {
                    ctx.options.humidity[i] = h;
                }
            }
        } else if key == b"SNOWMELT" && tok.len() >= 7 {
            // Legacy SNOWMELT Stemp ATIwt RNM Elev Lat DTLong [minMelt maxMelt]
            ctx.options.snow_divt = to_double(&tok[1], 0.0);
            ctx.options.snow_ati_wt = to_double(&tok[2], 0.0);
            ctx.options.snow_nrg_ratio = to_double(&tok[3], 0.0);
            ctx.options.snow_elev = to_double(&tok[4], 0.0);
            ctx.options.snow_lat = to_double(&tok[5], 0.0);
            // tok[6] = longitude/solar-time correction (minutes), stored verbatim.
            ctx.options.snow_dtlong = to_double(&tok[6], 0.0);
            if tok.len() >= 9 {
                ctx.options.snow_min_melt = to_double(&tok[7], 0.0);
                ctx.options.snow_max_melt = to_double(&tok[8], 0.0);
            }
        } else if key == b"ADC" && tok.len() >= 12 {
            let surface = tokenizer::to_upper(&tok[1]);
            if surface == b"IMPERVIOUS" {
                for i in 0..10 {
                    ctx.options.adc_imperv[i] = to_double(&tok[2 + i], 0.0);
                }
            } else if surface == b"PERVIOUS" {
                for i in 0..10 {
                    ctx.options.adc_perv[i] = to_double(&tok[2 + i], 0.0);
                }
            }
        }
    }
}

// ============================================================================
// handle_snowpacks()
// ============================================================================

/// [SNOWPACKS]: Name PLOWABLE|IMPERVIOUS|PERVIOUS v1..v7 | REMOVAL v1..v6 [Subcatch] | UEB ...
pub fn handle_snowpacks(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 2 {
            continue;
        }

        let name = &tok[0];
        let mut idx = ctx.snowpack_names.find(name);
        if idx < 0 {
            idx = ctx.snowpack_names.try_add(name); // add(): cannot be a duplicate here
            ctx.snowpacks.names.push(name.clone());
        }
        let u = idx as usize;

        ensure_snowpack_capacity(&mut ctx.snowpacks, idx);

        let surface = tokenizer::to_upper(&tok[1]);

        if surface == b"UEB" {
            let error = surface::set_snowpack_ueb(&mut ctx.snowpacks.config, idx, &tok[2..]);
            if !error.is_empty() {
                ctx.errors.push(cat(&[b"[SNOWPACKS] ", name, b": ", &error]));
            }
        } else if surface == b"PLOWABLE" && tok.len() >= 9 {
            for i in 0..7 {
                ctx.snowpacks.plowable[u][i] = to_double(&tok[2 + i], 0.0);
            }
        } else if surface == b"IMPERVIOUS" && tok.len() >= 9 {
            for i in 0..7 {
                ctx.snowpacks.impervious[u][i] = to_double(&tok[2 + i], 0.0);
            }
        } else if surface == b"PERVIOUS" && tok.len() >= 9 {
            for i in 0..7 {
                ctx.snowpacks.pervious[u][i] = to_double(&tok[2 + i], 0.0);
            }
        } else if surface == b"REMOVAL" && tok.len() >= 8 {
            for i in 0..6 {
                ctx.snowpacks.removal[u][i] = to_double(&tok[2 + i], 0.0);
            }
            if tok.len() >= 9 {
                ctx.snowpacks.removal_subcatch[u] = tok[8].clone();
            }
        }
    }
}

// ============================================================================
// handle_ueb_options()
// ============================================================================

/// [UEB_OPTIONS]: one setUebOption row each, unvalidated across rows (validate = false).
pub fn handle_ueb_options(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.is_empty() {
            continue;
        }
        let error = surface::set_ueb_option(&mut ctx.snowpacks.config, &tok, false);
        if !error.is_empty() {
            ctx.errors.push(cat(&[b"[UEB_OPTIONS] ", &error]));
        }
    }
}

/// meshSnowRows: section 0 options, 2 explicit cell rows, anything else defaults rows.
fn mesh_snow_rows(ctx: &mut SimulationContext, lines: &[Bytes], section: i32) {
    for line in lines {
        let tokens = tokenizer::tokenize(line);
        if tokens.is_empty() {
            continue;
        }
        let error = if section == 0 {
            surface::set_mesh_snow_option(&mut ctx.snowpacks.mesh_config, &tokens)
        } else {
            surface::set_mesh_snow_row(&mut ctx.snowpacks.mesh_config, &tokens, section == 2)
        };
        if !error.is_empty() {
            ctx.errors.push(cat(&[b"[2D_SNOW] ", &error]));
        }
    }
}

/// [2D_SNOW_OPTIONS]
pub fn handle_mesh_snow_options(ctx: &mut SimulationContext, lines: &[Bytes]) {
    mesh_snow_rows(ctx, lines, 0);
}

/// [2D_SNOW_DEFAULTS]
pub fn handle_mesh_snow_defaults(ctx: &mut SimulationContext, lines: &[Bytes]) {
    mesh_snow_rows(ctx, lines, 1);
}

/// [2D_SNOW]
pub fn handle_mesh_snow_cells(ctx: &mut SimulationContext, lines: &[Bytes]) {
    mesh_snow_rows(ctx, lines, 2);
}

// ============================================================================
// handle_aquifers()
// ============================================================================

/// [AQUIFERS]: Name Por WP FC Ksat Kslope Tslope ETu ETs Seep Ebot Egw Umc [ETupat].
pub fn handle_aquifers(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 13 {
            continue;
        }

        let name = &tok[0];
        // Legacy does NOT exempt [AQUIFERS] from duplicate-ID rejection.
        let idx = add_unique(&mut ctx.aquifer_names, name, &mut ctx.errors);
        if idx < 0 {
            continue; // duplicate ID (ERR 207)
        }
        ctx.aquifers.names.push(name.clone());
        let u = idx as usize;

        ensure_aquifer_capacity(&mut ctx.aquifers, idx);

        let a = &mut ctx.aquifers;
        a.porosity[u] = to_double(&tok[1], 0.0);
        a.wilting_point[u] = to_double(&tok[2], 0.0);
        a.field_capacity[u] = to_double(&tok[3], 0.0);
        a.conductivity[u] = to_double(&tok[4], 0.0);
        a.conduct_slope[u] = to_double(&tok[5], 0.0);
        a.tension_slope[u] = to_double(&tok[6], 0.0);
        a.upper_evap[u] = to_double(&tok[7], 0.0);
        a.lower_evap[u] = to_double(&tok[8], 0.0);
        a.lower_loss[u] = to_double(&tok[9], 0.0);
        a.bottom_elev[u] = to_double(&tok[10], 0.0);
        a.water_table_elev[u] = to_double(&tok[11], 0.0);
        a.upper_moist[u] = to_double(&tok[12], 0.0);

        if tok.len() >= 14 {
            a.upper_evap_pat[u] = tok[13].clone();
        }
    }
}

// ============================================================================
// handle_groundwater()
// ============================================================================

/// [GROUNDWATER]: Subcatch Aquifer Node SurfEl A1 B1 A2 B2 A3 Dsw [Egwt Ebot Wgw Umc].
pub fn handle_groundwater(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.is_empty() {
            continue;
        }
        // Legacy gwater_readGroundwaterParams requires >= 11 tokens; fewer is ERR_ITEMS.
        if tok.len() < 11 {
            ctx.errors.push(format_error(ERR_ITEMS, b""));
            continue;
        }

        // A row ahead of [SUBCATCHMENTS] is deferred and replayed at the end.
        let idx = ctx.subcatch_names.find(&tok[0]);
        if idx < 0 {
            ctx.deferred_section_rows.push((b"GROUNDWATER".to_vec(), line.clone()));
            continue;
        }
        let u = idx as usize;

        // Required numerics (legacy x[0..6]): a token getDouble() rejects is ERROR 211.
        let mut req = [0.0f64; 7];
        let mut bad = 0usize;
        let mut k = 0usize;
        while k < 7 && bad == 0 {
            if !parse_double_strict(&tok[k + 3], &mut req[k]) {
                bad = k + 3;
            }
            k += 1;
        }
        // Optional Egwt Ebot Wgw Umc (legacy x[7..10]): absent, or a leading `*`, is MISSING.
        let mut opt = [0.0f64; 4];
        let mut k = 0usize;
        while k < 4 && bad == 0 {
            let t = k + 10;
            opt[k] = constants::MISSING;
            if tok.len() > t && (tok[t].is_empty() || tok[t][0] != b'*') && !parse_double_strict(&tok[t], &mut opt[k]) {
                bad = t;
            }
            k += 1;
        }
        if bad != 0 {
            ctx.errors.push(format_error(ERR_NUMBER, &tok[bad]));
            continue;
        }

        ensure_subcatch_gw_capacity(ctx, idx);

        // Aquifer and node may be defined further down; PostParseResolver resolves them.
        ctx.subcatches.gw_aquifer[u] = ctx.aquifer_names.find(&tok[1]);
        if ctx.subcatches.gw_aquifer[u] < 0 {
            ctx.pending_gw_aquifers.push((idx, tok[1].clone()));
        }
        ctx.subcatches.gw_node[u] = ctx.node_names.find(&tok[2]);
        ctx.pending_gw_nodes.push((idx, tok[2].clone()));
        // Stored in the deck's units; the engine converts at init.
        let s = &mut ctx.subcatches;
        s.gw_surf_elev[u] = req[0];
        s.gw_a1[u] = req[1];
        s.gw_b1[u] = req[2];
        s.gw_a2[u] = req[3];
        s.gw_b2[u] = req[4];
        s.gw_a3[u] = req[5];
        s.gw_tw[u] = req[6];
        s.gw_hstar[u] = opt[0];
        s.gw_bot_elev[u] = opt[1];
        s.gw_wt_elev[u] = opt[2];
        s.gw_upper_moist[u] = opt[3];
    }
}

// ============================================================================
// handle_gwf()
// ============================================================================

/// [GWF]: Subcatch LATERAL|DEEP expression, kept as ext_options["GWF:<subcatch>:<type>"].
pub fn handle_gwf(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 3 {
            continue;
        }

        // Key on the [SUBCATCHMENTS] spelling; a name not defined yet is deferred.
        let subcatch: Bytes = match ctx.subcatch_names.canonical(&tok[0]) {
            Some(c) => c.to_vec(),
            None => {
                ctx.deferred_section_rows.push((b"GWF".to_vec(), line.clone()));
                continue;
            }
        };

        // Legacy gwater.c accepts any "LAT..." spelling for LATERAL.
        let type_tok = tokenizer::to_upper(&tok[1]);
        let type_: &[u8] = if type_tok.starts_with(b"LAT") {
            b"LATERAL"
        } else if type_tok == b"DEEP" {
            b"DEEP"
        } else {
            ctx.errors.push(format_error(ERR_KEYWORD, &tok[1]));
            continue;
        };

        // The expression is free text, taken verbatim from the line (commas kept).
        let mut rest: &[u8] = tokenizer::strip_comment(line);
        for _col in 0..2 {
            // skip Subcatch, Type
            let mut i = 0usize;
            while i < rest.len() && (rest[i] == b' ' || rest[i] == b'\t') {
                i += 1;
            }
            if i < rest.len() && rest[i] == b'"' {
                i += 1;
                while i < rest.len() && rest[i] != b'"' {
                    i += 1;
                }
                if i < rest.len() {
                    i += 1;
                }
            } else {
                while i < rest.len() && rest[i] != b' ' && rest[i] != b'\t' && rest[i] != b',' {
                    i += 1;
                }
            }
            while i < rest.len() && (rest[i] == b' ' || rest[i] == b'\t') {
                i += 1;
            }
            if i < rest.len() && rest[i] == b',' {
                i += 1; // optional CSV comma
            }
            rest = &rest[i..];
        }
        let expr: Bytes = tokenizer::trim(rest).to_vec();

        // mathexpr::parse is lenient, so reject malformed expressions like legacy
        // ERR_MATH_EXPR.
        let mut msg = Bytes::new();
        let mut col = -1i32;
        if groundwater_parse::gwf_validate(&expr, &mut msg, &mut col) != 0 {
            let detail = cat(&[b"in [GWF] ", type_, b" for Subcatchment ", &subcatch, b": ", &msg]);
            ctx.errors.push(format_error_detail(ERR_MATH_EXPR, b"", &detail));
            continue;
        }

        let key = cat(&[b"GWF:", &subcatch, b":", type_]);
        ctx.options.ext_options.insert(key, expr);
    }
}

// ============================================================================
// handle_lid_controls()
// ============================================================================

/// [LID_CONTROLS]: Name Type | Name Layer params... | Name REMOVALS (pollutant %)... and the
/// NODE layer rows.
pub fn handle_lid_controls(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 2 {
            continue;
        }

        let name = &tok[0];
        let mut idx = ctx.lid_names.find(name);
        if idx < 0 {
            idx = ctx.lid_names.try_add(name); // add(): cannot be a duplicate here
            ctx.lid_controls.names.push(name.clone());
        }
        let u = idx as usize;

        ensure_lid_capacity(&mut ctx.lid_controls, idx);

        let layer = tokenizer::to_upper(&tok[1]);

        // First line for a LID sets its type code (e.g., BC, RG, GR, etc.)
        if tok.len() == 2 {
            ctx.lid_controls.lid_type[u] = layer;
            continue;
        }

        if ctx.lid_controls.lid_type[u] == b"NODE" {
            let mut row = LidNodeLayer::default();
            let count: usize;
            if layer == b"SURFACE" {
                row.kind = LidNodeLayerKind::Surface;
                count = 2;
            } else if layer == b"MEDIA" {
                row.kind = LidNodeLayerKind::Media;
                count = 7;
            } else if layer == b"AGGREGATE" {
                row.kind = LidNodeLayerKind::Aggregate;
                count = 3;
            } else if layer == b"BOTTOM" {
                row.kind = LidNodeLayerKind::Bottom;
                count = 2;
            } else {
                ctx.errors.push(cat(&[b"Unknown NODE LID layer: ", &layer]));
                continue;
            }
            let mut valid = tok.len() == count + 2;
            let mut k = 0usize;
            while valid && k < count {
                // std::strtod(tok.c_str(), &end): the C string ends at the first NUL
                let t: &[u8] = &tok[k + 2];
                let cstr = &t[..t.iter().position(|&c| c == 0).unwrap_or(t.len())];
                let (v, end) = strtod(cstr);
                row.params[k] = v;
                valid = end != 0 && end == cstr.len() && row.params[k].is_finite();
                k += 1;
            }
            if !valid {
                ctx.errors.push(cat(&[b"Invalid NODE LID layer: ", line]));
            } else {
                ctx.lid_controls.node_layers[u].push(row);
            }
            continue;
        }

        // Subsequent lines define layer parameters
        let lc = &mut ctx.lid_controls;
        if layer == b"SURFACE" {
            let mut i = 0usize;
            while i < 5 && (i + 2) < tok.len() {
                lc.surface[u][i] = to_double(&tok[2 + i], 0.0);
                i += 1;
            }
        } else if layer == b"SOIL" {
            let mut i = 0usize;
            while i < 7 && (i + 2) < tok.len() {
                lc.soil[u][i] = to_double(&tok[2 + i], 0.0);
                i += 1;
            }
        } else if layer == b"PAVEMENT" {
            let mut i = 0usize;
            while i < 6 && (i + 2) < tok.len() {
                lc.pavement[u][i] = to_double(&tok[2 + i], 0.0);
                i += 1;
            }
        } else if layer == b"STORAGE" {
            let mut i = 0usize;
            while i < 4 && (i + 2) < tok.len() {
                lc.storage[u][i] = to_double(&tok[2 + i], 0.0);
                i += 1;
            }
        } else if layer == b"DRAIN" {
            let mut i = 0usize;
            while i < 6 && (i + 2) < tok.len() {
                lc.drain[u][i] = to_double(&tok[2 + i], 0.0);
                i += 1;
            }
            let nn = lc.names.len();
            lc.drain_curve.resize(nn, Bytes::new());
            if tok.len() > 8 {
                lc.drain_curve[u] = tok[8].clone();
            }
        } else if layer == b"DRAINMAT" {
            let mut i = 0usize;
            while i < 3 && (i + 2) < tok.len() {
                lc.drainmat[u][i] = to_double(&tok[2 + i], 0.0);
                i += 1;
            }
        } else if layer == b"REMOVALS" {
            // LID_ID REMOVALS  PollutName1  %Removal1  PollutName2  %Removal2 ...
            if lc.removals.len() <= u {
                lc.removals.resize(u + 1, Vec::new());
            }
            let mut i = 2usize;
            while i + 1 < tok.len() {
                let pi = ctx.pollutant_names.find(&tok[i]);
                let mut frac = to_double(&tok[i + 1], 0.0) / 100.0;
                frac = cmax(0.0, cmin(1.0, frac));
                let lc = &mut ctx.lid_controls;
                if pi >= 0 {
                    lc.removals[u].push((pi, frac));
                } else {
                    // A reaction species (or a pollutant declared later), bound later.
                    let byname = &mut lc.removals_by_name;
                    if byname.len() <= u {
                        byname.resize(u + 1, Vec::new());
                    }
                    byname[u].push((tok[i].clone(), frac));
                }
                i += 2;
            }
        } else {
            // Unknown layer or it's actually the type code with params
            lc.lid_type[u] = layer;
        }
    }
}

// ============================================================================
// handle_lid_usage()
// ============================================================================

/// [LID_USAGE]: Subcatch LID Number Area Width InitSat FromImp ToPerv [RptFile] [DrainTo]
/// [FromPerv].
pub fn handle_lid_usage(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 8 {
            continue;
        }

        let sc_idx = ctx.subcatch_names.find(&tok[0]);
        let lid_idx = ctx.lid_names.find(&tok[1]);
        if sc_idx < 0 || lid_idx < 0 {
            continue;
        }

        if ctx.lid_controls.lid_type[lid_idx as usize] == b"NODE" {
            ctx.errors.push(cat(&[b"NODE LID controls require [LID_NODES], not [LID_USAGE]: ", &tok[1]]));
            continue;
        }

        // legacy lid_readGroupParams atoi()s the count: negative is ERROR 211, zero adds none.
        let number = to_int(&tok[2], 0);
        if number < 0 {
            ctx.errors.push(format_error(ERR_NUMBER, &tok[2]));
            continue;
        }
        if number == 0 {
            continue;
        }

        let lu = &mut ctx.lid_usage;
        lu.subcatch_index.push(sc_idx);
        lu.lid_index.push(lid_idx);
        lu.number.push(number);
        lu.area.push(to_double(&tok[3], 0.0));
        lu.width.push(to_double(&tok[4], 0.0));
        lu.init_sat.push(to_double(&tok[5], 0.0));
        lu.from_imperv.push(to_double(&tok[6], 0.0));
        lu.to_perv.push(to_int(&tok[7], 0));

        // RptFile `*` is legacy's "no report file" (lid.c:475).
        let rpt: Bytes = if tok.len() > 8 && tok[8] != b"*" { tok[8].clone() } else { Bytes::new() };
        lu.rpt_file.push(FilePathPair { absolute: Bytes::new(), original: rpt });
        lu.drain_to.push(if tok.len() > 9 { tok[9].clone() } else { Bytes::new() });
        lu.from_perv.push(if tok.len() > 10 { to_double(&tok[10], 0.0) } else { 0.0 });
    }
}
