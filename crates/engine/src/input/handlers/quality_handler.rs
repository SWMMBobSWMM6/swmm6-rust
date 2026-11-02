// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/handlers/QualityHandler.cpp and
// QualityHandler.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Section handlers for [POLLUTANTS], [LANDUSES], [COVERAGES], [BUILDUP], [WASHOFF],
// [TREATMENT], [LOADINGS] and [INITIAL_QUALITY]. [TREATMENT] stores each expression as
// text; it is compiled later (quality/Treatment.cpp), not here.
//
// Index arithmetic keeps upstream's widths: `lu_idx * n_pollutants + p_idx` and the like
// are int products converted to size_t afterwards, as in the C++.

use crate::core::charconv::from_chars_double;
use crate::core::constants;
use crate::core::context::SimulationContext;
use crate::core::error_codes::{format_error, ERR_ITEMS, ERR_NAME, ERR_NUMBER};
use crate::core::options::Bytes;
use crate::data::msx_surface::{MsxSurfaceRowsBuildup, MsxSurfaceRowsLoading, MsxSurfaceRowsWashoff};
use crate::data::pollutant::MassUnits;
use crate::input::parse_utils::{add_unique, parse_double_strict, to_double};
use crate::input::section_parser::parse_section;
use crate::input::tokenizer;

fn cat(parts: &[&[u8]]) -> Bytes {
    parts.concat()
}

/// findmatch_prefix: the first nonempty exact-case word that prefixes `token`, or -1.
fn findmatch_prefix(token: &[u8], words: &[&[u8]]) -> i32 {
    let mut i = 0;
    for w in words {
        // The tokenizer has already stripped the leading blanks legacy's match() skips.
        let n = w.len();
        if n > 0 && token.len() >= n && &token[..n] == *w {
            return i;
        }
        i += 1;
    }
    -1
}

// ============================================================================
// handle_pollutants()
// ============================================================================

/// [POLLUTANTS]: Name Units Crain Cgw Crdii Kdecay SnowOnly CoPollut CoFrac Cdwf Cinit [MW].
pub fn handle_pollutants(ctx: &mut SimulationContext, lines: &[Bytes]) {
    let parsed = parse_section(lines);

    // First pass: register all pollutant names so co-pollutant lookup works
    for pl in &parsed {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.len() < 2 {
            continue;
        }
        let name = &tok[0];
        add_unique(&mut ctx.pollutant_names, name, &mut ctx.errors);
        // duplicate ID pushes ERR 207 and is skipped; the second pass repopulates the
        // original.
    }

    // Resize pollutant arrays to accommodate all pollutants
    let n = ctx.pollutant_names.size();
    ctx.pollutants.resize_pollutants(n);

    // Second pass: populate pollutant properties and comments
    for pl in &parsed {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.len() < 2 {
            continue;
        }

        let idx = ctx.pollutant_names.find(&tok[0]);
        if idx < 0 {
            continue;
        }
        let u = idx as usize;
        let p = &mut ctx.pollutants;

        // Units
        let units_str = tokenizer::to_upper(&tok[1]);
        if units_str == b"UG/L" {
            p.units[u] = MassUnits::UgPerL;
        } else if units_str == b"#/L" {
            p.units[u] = MassUnits::CountsPerL;
        } else {
            p.units[u] = MassUnits::MgPerL;
        }

        if tok.len() > 2 {
            p.c_rain[u] = to_double(&tok[2], 0.0);
        }
        if tok.len() > 3 {
            p.c_gw[u] = to_double(&tok[3], 0.0);
        }
        if tok.len() > 4 {
            p.c_rdii[u] = to_double(&tok[4], 0.0);
        }
        // The Kdecay column is 1/day; store 1/sec, as legacy does (KD1).
        if tok.len() > 5 {
            p.k_decay[u] = to_double(&tok[5], 0.0) / constants::SEC_PER_DAY;
        }

        // SnowOnly
        if tok.len() > 6 {
            let snow = tokenizer::to_upper(&tok[6]);
            p.snow_only[u] = snow == b"YES" || snow == b"1";
        }

        // CoPollut
        if tok.len() > 7 {
            if tok[7] == b"*" || tok[7] == b"-" {
                p.co_pollut[u] = -1;
            } else {
                p.co_pollut[u] = ctx.pollutant_names.find(&tok[7]);
            }
        }
        let p = &mut ctx.pollutants;

        // CoFrac
        if tok.len() > 8 {
            p.co_frac[u] = to_double(&tok[8], 0.0);
        }

        // Cdwf — default dry weather sanitary flow concentration
        if tok.len() > 9 {
            p.c_dwf[u] = to_double(&tok[9], 0.0);
        }

        // Native extension after the legacy Cinit column.
        if tok.len() > 11 {
            p.mwt[u] = to_double(&tok[11], 0.0);
        }

        // Cinit
        if tok.len() > 10 {
            p.init_conc[u] = to_double(&tok[10], 0.0);
        }

        if !pl.comment.is_empty() {
            p.comments[u] = pl.comment.clone();
        }
    }
}

// ============================================================================
// handle_landuses()
// ============================================================================

/// [LANDUSES]: Name [SweepInterval] [SweepRemoval] [SweepDays0].
pub fn handle_landuses(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for pl in parse_section(lines) {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.is_empty() {
            continue;
        }

        let name = &tok[0];
        let idx = add_unique(&mut ctx.landuse_names, name, &mut ctx.errors);
        if idx < 0 {
            continue; // duplicate ID (ERR 207, legacy input.c parity)
        }
        let u = idx as usize;

        // Ensure capacity
        let n = (idx + 1) as usize;
        let l = &mut ctx.landuses;
        if l.sweep_interval.len() < n {
            l.sweep_interval.resize(n, 0.0);
            l.sweep_removal.resize(n, 0.0);
            l.last_swept.resize(n, 0.0);
            l.comments.resize(n, Bytes::new());
        }

        if tok.len() > 1 {
            l.sweep_interval[u] = to_double(&tok[1], 0.0);
        }
        if tok.len() > 2 {
            l.sweep_removal[u] = to_double(&tok[2], 0.0);
        }
        if tok.len() > 3 {
            l.last_swept[u] = to_double(&tok[3], 0.0);
        }
        if !pl.comment.is_empty() {
            l.comments[u] = pl.comment.clone();
        }
    }
}

// ============================================================================
// handle_coverages()
// ============================================================================

/// [COVERAGES]: Subcatch LandUse Percent [LandUse Percent ...].
pub fn handle_coverages(ctx: &mut SimulationContext, lines: &[Bytes]) {
    let n_landuses = ctx.landuse_names.size();
    if n_landuses <= 0 {
        return;
    }

    let n_subcatches = ctx.subcatch_names.size();
    if n_subcatches <= 0 {
        return;
    }

    // Allocate flat coverage array if needed
    let total = (n_subcatches * n_landuses) as usize;
    if ctx.subcatches.coverage.len() < total {
        ctx.subcatches.coverage.clear();
        ctx.subcatches.coverage.resize(total, 0.0);
    }
    // The per-(subcatchment, land use) last-swept counters are sized with the coverage
    // matrix, as resize_coverage does.
    if ctx.subcatches.sweep_last_swept.len() < total {
        ctx.subcatches.sweep_last_swept.clear();
        ctx.subcatches.sweep_last_swept.resize(total, 0.0);
    }
    ctx.subcatches.coverage_n_landuses = n_landuses;

    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.is_empty() {
            continue;
        }
        // Legacy subcatch_readLanduseParams: every pair is read, the last value for a
        // repeated land use wins; an unknown name, a missing percent or a bad number is an
        // error.
        if tok.len() < 3 {
            ctx.errors.push(format_error(ERR_ITEMS, b""));
            continue;
        }

        let sub_idx = ctx.subcatch_names.find(&tok[0]);
        if sub_idx < 0 {
            ctx.errors.push(format_error(ERR_NAME, &tok[0]));
            continue;
        }

        let mut k = 1usize;
        while k < tok.len() {
            let lu_idx = ctx.landuse_names.find(&tok[k]);
            if lu_idx < 0 {
                ctx.errors.push(format_error(ERR_NAME, &tok[k]));
                break;
            }
            if k + 1 >= tok.len() {
                ctx.errors.push(format_error(ERR_ITEMS, b""));
                break;
            }
            let mut pct = 0.0f64;
            if !parse_double_strict(&tok[k + 1], &mut pct) {
                ctx.errors.push(format_error(ERR_NUMBER, &tok[k + 1]));
                break;
            }
            ctx.subcatches.coverage[(sub_idx * n_landuses + lu_idx) as usize] = pct;
            k += 2;
        }
    }
}

// ============================================================================
// handle_buildup()
// ============================================================================

const BUILDUP_WORDS: [&[u8]; 5] = [b"NONE", b"POW", b"EXP", b"SAT", b"EXT"];
const WASHOFF_WORDS: [&[u8]; 4] = [b"NONE", b"EXP", b"RC", b"EMC"];

/// [BUILDUP]: Landuse Pollutant Type C1 C2 [C3|TSname] [PerUnit].
pub fn handle_buildup(ctx: &mut SimulationContext, lines: &[Bytes]) {
    let n_landuses = ctx.landuse_names.size();
    let n_pollutants = ctx.pollutant_names.size();
    if n_landuses <= 0 {
        return;
    }

    // Ensure buildup arrays are sized
    if n_pollutants > 0 && (ctx.buildup.n_landuses != n_landuses || ctx.buildup.n_pollutants != n_pollutants) {
        ctx.buildup.resize(n_landuses, n_pollutants);
    }

    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 5 {
            continue;
        }

        let lu_idx = ctx.landuse_names.find(&tok[0]);
        if lu_idx < 0 {
            continue;
        }

        let p_idx = ctx.pollutant_names.find(&tok[1]);
        if p_idx < 0 {
            // BW-MSX: not a pollutant; park the row by name for the reactions component.
            let mut r = MsxSurfaceRowsBuildup::default();
            r.landuse = tok[0].clone();
            r.species = tok[1].clone();
            let ts = tokenizer::to_upper(&tok[2]);
            let mbt = findmatch_prefix(&ts, &BUILDUP_WORDS);
            r.func_type = if mbt >= 0 { mbt } else { 0 };
            r.c1 = to_double(&tok[3], 0.0);
            r.c2 = to_double(&tok[4], 0.0);
            if tok.len() > 5 {
                if r.func_type == 4 {
                    r.ts_name = tok[5].clone();
                } else {
                    r.c3 = to_double(&tok[5], 0.0);
                }
            }
            if tok.len() > 6 {
                r.normalizer = if tokenizer::to_upper(&tok[6]) == b"CURB" { 1 } else { 0 };
            }
            ctx.reactions.surface.rows.buildup.push(r);
            continue;
        }
        if n_pollutants <= 0 {
            continue;
        }

        let flat = (lu_idx * n_pollutants + p_idx) as usize;

        // Type: NONE=0, POW=1, EXP=2, SAT=3, EXT=4, matched by PREFIX.
        let type_str = tokenizer::to_upper(&tok[2]);
        let bt = findmatch_prefix(&type_str, &BUILDUP_WORDS);
        ctx.buildup.func_type[flat] = if bt >= 0 { bt } else { 0 };

        ctx.buildup.coeff1[flat] = to_double(&tok[3], 0.0);
        ctx.buildup.coeff2[flat] = to_double(&tok[4], 0.0);

        if tok.len() > 5 {
            // EXT buildup names a time series in tok[5]; other types give a number.
            if ctx.buildup.func_type[flat] == 4 {
                // EXTERNAL
                let ts_idx = ctx.find_timeseries(&tok[5]);
                ctx.buildup.coeff3[flat] = ts_idx as f64;
            } else {
                ctx.buildup.coeff3[flat] = to_double(&tok[5], 0.0);
            }
        }

        // Normalizer: AREA=0, CURB=1
        if tok.len() > 6 {
            let norm = tokenizer::to_upper(&tok[6]);
            ctx.buildup.normalizer[flat] = if norm == b"CURB" { 1 } else { 0 };
        }
    }
}

// ============================================================================
// handle_washoff()
// ============================================================================

/// [WASHOFF]: Landuse Pollutant Type C1 C2 [SweepEffic] [BMPEffic].
pub fn handle_washoff(ctx: &mut SimulationContext, lines: &[Bytes]) {
    let n_landuses = ctx.landuse_names.size();
    let n_pollutants = ctx.pollutant_names.size();
    if n_landuses <= 0 {
        return;
    }

    // Ensure washoff arrays are sized
    if n_pollutants > 0 && (ctx.washoff.n_landuses != n_landuses || ctx.washoff.n_pollutants != n_pollutants) {
        ctx.washoff.resize(n_landuses, n_pollutants);
    }

    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 5 {
            continue;
        }

        let lu_idx = ctx.landuse_names.find(&tok[0]);
        if lu_idx < 0 {
            continue;
        }

        let p_idx = ctx.pollutant_names.find(&tok[1]);
        if p_idx < 0 {
            // BW-MSX: park by name (see handle_buildup).
            let mut r = MsxSurfaceRowsWashoff::default();
            r.landuse = tok[0].clone();
            r.species = tok[1].clone();
            let ts = tokenizer::to_upper(&tok[2]);
            let mwt = findmatch_prefix(&ts, &WASHOFF_WORDS);
            r.func_type = if mwt >= 0 { mwt } else { 0 };
            r.coeff = to_double(&tok[3], 0.0);
            r.expon = to_double(&tok[4], 0.0);
            if tok.len() > 5 {
                r.sweep_effic = to_double(&tok[5], 0.0);
            }
            if tok.len() > 6 {
                r.bmp_effic = to_double(&tok[6], 0.0);
            }
            ctx.reactions.surface.rows.washoff.push(r);
            continue;
        }
        if n_pollutants <= 0 {
            continue;
        }

        let flat = (lu_idx * n_pollutants + p_idx) as usize;

        // Type: NONE=0, EXP=1, RC=2, EMC=3, matched by PREFIX.
        let type_str = tokenizer::to_upper(&tok[2]);
        let wt = findmatch_prefix(&type_str, &WASHOFF_WORDS);
        ctx.washoff.func_type[flat] = if wt >= 0 { wt } else { 0 };

        ctx.washoff.coeff[flat] = to_double(&tok[3], 0.0);
        ctx.washoff.expon[flat] = to_double(&tok[4], 0.0);

        if tok.len() > 5 {
            ctx.washoff.sweep_effic[flat] = to_double(&tok[5], 0.0);
        }
        if tok.len() > 6 {
            ctx.washoff.bmp_effic[flat] = to_double(&tok[6], 0.0);
        }
    }
}

// ============================================================================
// handle_treatment()
// ============================================================================

/// [TREATMENT]: Node Pollutant Result = Expression, kept as text (tokens rejoined by one
/// space).
pub fn handle_treatment(ctx: &mut SimulationContext, lines: &[Bytes]) {
    let n_pollutants = ctx.pollutant_names.size();
    let n_nodes = ctx.node_names.size();
    if n_pollutants <= 0 || n_nodes <= 0 {
        return;
    }

    // Ensure treatment arrays are sized
    if ctx.treatment.n_nodes != n_nodes || ctx.treatment.n_pollutants != n_pollutants {
        ctx.treatment.resize(n_nodes, n_pollutants);
    }

    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 3 {
            continue;
        }

        let node_idx = ctx.node_names.find(&tok[0]);
        if node_idx < 0 {
            continue;
        }

        let p_idx = ctx.pollutant_names.find(&tok[1]);
        if p_idx < 0 {
            continue;
        }

        let flat = (node_idx * n_pollutants + p_idx) as usize;

        // Reconstruct the expression from remaining tokens (e.g. "R = 0.5 * ...")
        let mut expr = Bytes::new();
        for i in 2..tok.len() {
            if i > 2 {
                expr.push(b' ');
            }
            expr.extend_from_slice(&tok[i]);
        }
        ctx.treatment.expressions[flat] = expr;
    }
}

// ============================================================================
// handle_loadings()
// ============================================================================

/// [LOADINGS]: Subcatch Pollutant InitBuildup (non-pollutant names parked for MSX).
pub fn handle_loadings(ctx: &mut SimulationContext, lines: &[Bytes]) {
    let n_pollutants = ctx.pollutant_names.size();
    // BW-MSX: rows naming a non-pollutant are parked by name, even with no [POLLUTANTS].
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 3 {
            continue;
        }
        if ctx.pollutant_names.find(&tok[1]) >= 0 {
            continue;
        }
        let mut r = MsxSurfaceRowsLoading::default();
        r.subcatch = tok[0].clone();
        r.species = tok[1].clone();
        r.value = to_double(&tok[2], 0.0);
        ctx.reactions.surface.rows.loadings.push(r);
    }
    if n_pollutants <= 0 {
        return;
    }

    // Size the quality arrays now, so the `flat < len()` guard below can admit rows.
    let n_sub = ctx.subcatch_names.size();
    if n_sub > 0
        && (ctx.subcatches.conc_n_pollutants != n_pollutants || ctx.subcatches.conc.len() as i32 != n_sub * n_pollutants)
    {
        ctx.subcatches.resize_quality(n_pollutants);
    }

    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 3 {
            continue;
        }

        let sub_idx = ctx.subcatch_names.find(&tok[0]);
        if sub_idx < 0 {
            continue;
        }

        let p_idx = ctx.pollutant_names.find(&tok[1]);
        if p_idx < 0 {
            continue;
        }

        let flat = (sub_idx * n_pollutants + p_idx) as usize;
        // An initial surface BUILDUP per unit area, not a runoff concentration.
        if flat < ctx.subcatches.init_loading.len() {
            ctx.subcatches.init_loading[flat] = to_double(&tok[2], 0.0);
        }
    }
}

// ============================================================================
// handle_initial_quality()
// ============================================================================

/// [INITIAL_QUALITY]: NODE|LINK element constituent value, or FILE <path>.
pub fn handle_initial_quality(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.is_empty() {
            continue;
        }
        // U2: `FILE <path>` names a CSV sidecar read at open.
        if tokenizer::to_upper(&tok[0]) == b"FILE" {
            if tok.len() < 2 {
                ctx.errors.push(b"[INITIAL_QUALITY] FILE needs a path.".to_vec());
                continue;
            }
            ctx.initial_quality.file = tok[1].clone();
            continue;
        }
        if tok.len() < 4 {
            ctx.errors.push(cat(&[b"[INITIAL_QUALITY] row needs NODE|LINK element constituent value: '", line, b"'."]));
            continue;
        }

        let scope = tokenizer::to_upper(&tok[0]);
        let mut link = false;
        if scope == b"LINK" {
            link = true;
        } else if scope != b"NODE" {
            ctx.errors.push(cat(&[b"[INITIAL_QUALITY] scope must be NODE or LINK: '", &tok[0], b"'."]));
            continue;
        }

        // Value must be fully numeric.
        let vs = &tok[3];
        let (r, used) = from_chars_double(vs);
        let v = match r {
            Ok(v) if used == vs.len() => v,
            _ => {
                ctx.errors.push(cat(&[
                    b"[INITIAL_QUALITY] bad value '",
                    vs,
                    b"' for '",
                    &tok[2],
                    b"' at ",
                    &scope,
                    b" '",
                    &tok[1],
                    b"'.",
                ]));
                continue;
            }
        };

        ctx.initial_quality.add(link, &tok[1], &tok[2], v);
    }
}
