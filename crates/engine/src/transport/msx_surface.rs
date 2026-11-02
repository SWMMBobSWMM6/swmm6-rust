// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/quality/MsxSurfaceQuality.hpp and
// MsxSurfaceQuality.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// BW-MSX: the open-time part of surface buildup / washoff of reactions-component species
// (namespace msxsurf): resolve() binds the [BUILDUP]/[WASHOFF]/[LOADINGS] rows that name
// a species (or a pollutant declared after the section), resolves LID REMOVALS by name,
// and activates the species store for rainfall-only decks. The run-time kernels in the
// same upstream file are not translated here.

use crate::core::cmath::cmax;
use crate::core::context::SimulationContext;
use crate::core::options::{Bytes, SimulationOptions};
use crate::core::unit_conversion::{ucf, MASS};
use crate::data::msx_surface::{MsxSurfaceRows, MsxSurfaceRowsBuildup, MsxSurfaceRowsLoading, MsxSurfaceRowsWashoff};

fn upper(s: &[u8]) -> Bytes {
    s.iter().map(|c| c.to_ascii_uppercase()).collect()
}

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

/// D-A27: concentration-mass -> user-mass factor from the species' units.
pub fn species_unit_factor(units: &[u8], opts: &SimulationOptions) -> f64 {
    let mass_ucf = ucf(MASS, opts);
    let u = upper(units);
    if u.starts_with(b"UG") {
        return mass_ucf / 1000.0;
    }
    if u.starts_with(b"MG") {
        return mass_ucf;
    }
    1.0 // MMOL, counts, dimensionless: the pollutant COUNTS rule
}

pub fn active(ctx: &SimulationContext) -> bool {
    ctx.reactions.surface.active()
}

fn apply_late_pollutant_buildup(ctx: &mut SimulationContext, r: &MsxSurfaceRowsBuildup) -> bool {
    let p = ctx.pollutant_names.find(&r.species);
    let lu = ctx.landuse_names.find(&r.landuse);
    if p < 0 || lu < 0 {
        return p >= 0;
    }
    let np = ctx.n_pollutants();
    let nlu = ctx.n_landuses();
    if ctx.buildup.n_landuses != nlu || ctx.buildup.n_pollutants != np {
        ctx.buildup.resize(nlu, np);
    }
    let k = (lu * np + p) as usize;
    ctx.buildup.func_type[k] = r.func_type;
    ctx.buildup.coeff1[k] = r.c1;
    ctx.buildup.coeff2[k] = r.c2;
    ctx.buildup.coeff3[k] = if r.func_type == 4 { ctx.find_timeseries(&r.ts_name) as f64 } else { r.c3 };
    ctx.buildup.normalizer[k] = r.normalizer;
    true
}

fn apply_late_pollutant_washoff(ctx: &mut SimulationContext, r: &MsxSurfaceRowsWashoff) -> bool {
    let p = ctx.pollutant_names.find(&r.species);
    let lu = ctx.landuse_names.find(&r.landuse);
    if p < 0 || lu < 0 {
        return p >= 0;
    }
    let np = ctx.n_pollutants();
    let nlu = ctx.n_landuses();
    if ctx.washoff.n_landuses != nlu || ctx.washoff.n_pollutants != np {
        ctx.washoff.resize(nlu, np);
    }
    let k = (lu * np + p) as usize;
    ctx.washoff.func_type[k] = r.func_type;
    ctx.washoff.coeff[k] = r.coeff;
    ctx.washoff.expon[k] = r.expon;
    ctx.washoff.sweep_effic[k] = r.sweep_effic;
    ctx.washoff.bmp_effic[k] = r.bmp_effic;
    true
}

fn apply_late_pollutant_loading(ctx: &mut SimulationContext, r: &MsxSurfaceRowsLoading) -> bool {
    let p = ctx.pollutant_names.find(&r.species);
    if p < 0 {
        return false;
    }
    let sc = ctx.subcatch_names.find(&r.subcatch);
    let np = ctx.n_pollutants();
    if sc >= 0 && ctx.subcatches.conc_n_pollutants == np {
        let k = (sc * np + p) as usize;
        if k < ctx.subcatches.conc.len() {
            ctx.subcatches.conc[k] = r.value;
        }
    }
    true
}

/// Multispecies Phase 2 (S3): LID REMOVALS entries kept by name at parse time.
pub fn resolve_lid_removals(ctx: &mut SimulationContext, warnings: &mut Vec<Bytes>) {
    let nm = if ctx.reactions.configured { ctx.reactions.n_species() } else { 0 };
    let nc = ctx.lid_controls.count();
    {
        let s = &mut ctx.reactions.surface;
        s.lid_rmvl.clear();
        s.lid_rmvl.resize(if nm > 0 && nc > 0 { (nc as usize) * (nm as usize) } else { 0 }, 0.0);
    }
    for c in 0..ctx.lid_controls.removals_by_name.len() {
        let rows = std::mem::take(&mut ctx.lid_controls.removals_by_name[c]);
        let mut keep: Vec<(Bytes, f64)> = Vec::new();
        for (name, frac) in &rows {
            let pi = ctx.pollutant_names.find(name);
            if pi >= 0 {
                let lc = &mut ctx.lid_controls;
                if lc.removals.len() <= c {
                    lc.removals.resize(c + 1, Vec::new());
                }
                lc.removals[c].push((pi, *frac));
                continue;
            }
            keep.push((name.clone(), *frac));
            let m = if nm > 0 { ctx.reactions.find_species(name) } else { -1 };
            if m >= 0 && (c as i32) < nc {
                ctx.reactions.surface.lid_rmvl[c * (nm as usize) + (m as usize)] = *frac;
            } else {
                let cname: Bytes = if c < ctx.lid_controls.names.len() { ctx.lid_controls.names[c].clone() } else { b"?".to_vec() };
                warnings.push(cat(&[
                    b"[LID_CONTROLS] ",
                    &cname,
                    b" REMOVALS: '",
                    name,
                    b"' is neither a pollutant nor a reaction species; ignored.",
                ]));
            }
        }
        ctx.lid_controls.removals_by_name[c] = keep;
    }
}

fn resolve_rows(ctx: &mut SimulationContext) -> Vec<Bytes> {
    let mut warnings: Vec<Bytes> = Vec::new();
    ctx.reactions.surface.resolved = false;
    resolve_lid_removals(ctx, &mut warnings);
    if ctx.reactions.surface.rows.is_empty() {
        return warnings;
    }

    // Rows naming a pollutant declared after the section go back to the pollutant tables.
    {
        let rows = std::mem::take(&mut ctx.reactions.surface.rows);
        let mut keep = MsxSurfaceRows::default();
        for r in &rows.buildup {
            if !apply_late_pollutant_buildup(ctx, r) {
                keep.buildup.push(r.clone());
            }
        }
        for r in &rows.washoff {
            if !apply_late_pollutant_washoff(ctx, r) {
                keep.washoff.push(r.clone());
            }
        }
        for r in &rows.loadings {
            if !apply_late_pollutant_loading(ctx, r) {
                keep.loadings.push(r.clone());
            }
        }
        ctx.reactions.surface.rows = keep;
    }
    if ctx.reactions.surface.rows.is_empty() {
        return warnings;
    }

    let nm = if ctx.reactions.configured { ctx.reactions.n_species() } else { 0 };
    let nlu = ctx.n_landuses();
    let nsc = ctx.n_subcatches();
    if nm <= 0 {
        let rows = &ctx.reactions.surface.rows;
        let n = rows.buildup.len() + rows.washoff.len() + rows.loadings.len();
        warnings.push(cat(&[
            b"[BUILDUP]/[WASHOFF]/[LOADINGS]: ",
            n.to_string().as_bytes(),
            b" row(s) name a constituent that is neither a [POLLUTANTS] entry nor a reactions-component species; they are ignored.",
        ]));
        return warnings;
    }
    {
        let s = &mut ctx.reactions.surface;
        s.resize_params(cmax(nlu, 0), nm);
        s.init_loading.clear();
        s.init_loading.resize((cmax(nsc, 0) as usize) * (nm as usize), 0.0);
    }
    for m in 0..nm {
        let um = m as usize;
        let f = species_unit_factor(
            if um < ctx.reactions.species_units.len() { &ctx.reactions.species_units[um] } else { &[] },
            &ctx.options,
        );
        ctx.reactions.surface.mcf[um] = f;
    }

    let species = |ctx: &SimulationContext, name: &[u8], section: &[u8], warnings: &mut Vec<Bytes>| -> i32 {
        let m = ctx.reactions.find_species(name);
        if m < 0 {
            warnings.push(cat(&[b"[", section, b"] '", name, b"' is neither a pollutant nor a reactions-component species; row ignored."]));
            return -1;
        }
        if (m as usize) < ctx.reactions.species_is_wall.len() && ctx.reactions.species_is_wall[m as usize] != 0 {
            warnings.push(cat(&[b"[", section, b"] '", name, b"' is a WALL species and cannot build up on a surface; row ignored."]));
            return -1;
        }
        m
    };
    let landuse = |ctx: &SimulationContext, name: &[u8], section: &[u8], warnings: &mut Vec<Bytes>| -> i32 {
        let lu = ctx.landuse_names.find(name);
        if lu < 0 {
            warnings.push(cat(&[b"[", section, b"] land use '", name, b"' is not defined; row ignored."]));
        }
        lu
    };

    let rows = ctx.reactions.surface.rows.clone();
    let mut any = false;
    for r in &rows.buildup {
        let lu = landuse(ctx, &r.landuse, b"BUILDUP", &mut warnings);
        if lu < 0 {
            continue;
        }
        let m = species(ctx, &r.species, b"BUILDUP", &mut warnings);
        if m < 0 {
            continue;
        }
        let k = ctx.reactions.surface.pidx(lu, m);
        let ts = if r.func_type == 4 { ctx.find_timeseries(&r.ts_name) } else { 0 };
        let s = &mut ctx.reactions.surface;
        s.bu_type[k] = r.func_type;
        s.bu_c1[k] = r.c1;
        s.bu_c2[k] = r.c2;
        if r.func_type == 4 {
            if ts < 0 {
                warnings.push(cat(&[b"[BUILDUP] EXT time series '", &r.ts_name, b"' for '", &r.species, b"' not found; no external buildup."]));
            }
            s.bu_c3[k] = ts as f64;
        } else {
            s.bu_c3[k] = r.c3;
        }
        s.bu_normalizer[k] = r.normalizer;
        any = any || r.func_type != 0;
    }
    for r in &rows.washoff {
        let lu = landuse(ctx, &r.landuse, b"WASHOFF", &mut warnings);
        if lu < 0 {
            continue;
        }
        let m = species(ctx, &r.species, b"WASHOFF", &mut warnings);
        if m < 0 {
            continue;
        }
        let s = &mut ctx.reactions.surface;
        let k = s.pidx(lu, m);
        s.wo_type[k] = r.func_type;
        s.wo_coeff[k] = r.coeff;
        s.wo_expon[k] = r.expon;
        s.wo_sweep_effic[k] = r.sweep_effic;
        s.wo_bmp_effic[k] = r.bmp_effic;
        any = any || r.func_type != 0;
    }
    for r in &rows.loadings {
        let sc = ctx.subcatch_names.find(&r.subcatch);
        if sc < 0 {
            warnings.push(cat(&[b"[LOADINGS] subcatchment '", &r.subcatch, b"' is not defined; row ignored."]));
            continue;
        }
        let m = species(ctx, &r.species, b"LOADINGS", &mut warnings);
        if m < 0 {
            continue;
        }
        let s = &mut ctx.reactions.surface;
        let k = s.sidx(sc, m);
        s.init_loading[k] = r.value;
        any = true;
    }
    ctx.reactions.surface.resolved = any && nlu > 0;
    refresh_derived(ctx);
    warnings
}

/// Multispecies Phase 2 (S4): a species that rains in deposits on every subcatchment.
pub fn activate_for_rain(ctx: &mut SimulationContext) {
    if ctx.reactions.surface.resolved || !ctx.reactions.configured {
        return;
    }
    let rx = &ctx.reactions;
    let nm = rx.n_species();
    let mut rain = false;
    let mut m = 0;
    while m < nm && m < rx.rain_conc.len() as i32 {
        rain = rain || rx.rain_conc[m as usize] > 0.0;
        m += 1;
    }
    if !rain && !rx.has_subcatch_reactions() && rx.buildup_laws.is_empty() && rx.lid_controls.is_empty() {
        return;
    }
    let nlu = cmax(ctx.n_landuses(), 0);
    if ctx.reactions.surface.n_species != nm || ctx.reactions.surface.n_landuses != nlu {
        let nsc = cmax(ctx.n_subcatches(), 0);
        {
            let s = &mut ctx.reactions.surface;
            s.resize_params(nlu, nm);
            s.init_loading.clear();
            s.init_loading.resize((nsc as usize) * (nm as usize), 0.0);
        }
        for m in 0..nm {
            let um = m as usize;
            let f = species_unit_factor(
                if um < ctx.reactions.species_units.len() { &ctx.reactions.species_units[um] } else { &[] },
                &ctx.options,
            );
            ctx.reactions.surface.mcf[um] = f;
        }
        refresh_derived(ctx);
    }
    ctx.reactions.surface.resolved = true;
}

/// Resolve named surface rows, then activate rainfall/reaction-only species storage.
pub fn resolve(ctx: &mut SimulationContext) -> Vec<Bytes> {
    let warnings = resolve_rows(ctx);
    activate_for_rain(ctx);
    warnings
}

/// Recalculate analytic buildup saturation times after a parameter edit.
pub fn refresh_derived(ctx: &mut SimulationContext) {
    let s = &mut ctx.reactions.surface;
    for k in 0..s.bu_type.len() {
        let c0 = s.bu_c1[k];
        let c1 = s.bu_c2[k];
        let c2 = s.bu_c3[k];
        let mut md = 0.0;
        if s.bu_type[k] == 2 && c1 > 0.0 {
            md = -(0.001f64).ln() / c1; // LIBM
        } else if s.bu_type[k] == 1 && c1 > 0.0 && c2 > 0.0 && c0 > 0.0 {
            md = (c0 / c1).powf(1.0 / c2); // LIBM
        } else if s.bu_type[k] == 3 && c2 > 0.0 && c0 > 0.0 {
            md = 999.0 * c2;
        }
        s.bu_max_days[k] = md;
    }
}
