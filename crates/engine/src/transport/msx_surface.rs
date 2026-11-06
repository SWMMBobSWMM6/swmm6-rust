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
// and activates the species store for rainfall-only decks; and the init-time initState
// (with the file-local buildupAt), which SWMMEngine::initQuality calls; and deliver
// (translated 5 November 2026), which SWMMEngine::stepRouting calls each routing step. The
// other run-time kernels in the same upstream file are not translated here.

use crate::core::cmath::{cmax, cmin};
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

/// buildupAt (file-local): the analytic buildup curve at `days` (0 for days <= 0 or NaN, and
/// for NONE / EXTERNAL).
fn buildup_at(type_: i32, c0: f64, c1: f64, c2: f64, days: f64) -> f64 {
    if !(days > 0.0) {
        return 0.0;
    }
    match type_ {
        1 => cmin(c1 * days.powf(c2), c0),        // POW  // LIBM
        2 => c0 * (1.0 - (-c1 * days).exp()),     // EXP  // LIBM
        3 => if c2 + days > 0.0 { c0 * days / (c2 + days) } else { 0.0 }, // SAT
        _ => 0.0,
    }
}

/// initState: reset the species surface state and seed the initial dry buildup from the
/// [LOADINGS] rows or the buildup function at DRY_DAYS (legacy landuse_getInitBuildup).
pub fn init_state(ctx: &mut SimulationContext) {
    activate_for_rain(ctx); // a rain concentration set through the API after open
    if !ctx.reactions.surface.active() {
        return;
    }
    let nsc = ctx.n_subcatches();
    let nlu = ctx.reactions.surface.n_landuses;
    let nm = ctx.reactions.surface.n_species;
    ctx.reactions.surface.resize_state(nsc);
    let dry_days = ctx.options.dry_days;
    for i in 0..nsc {
        let ui = i as usize;
        for lu in 0..nlu {
            let cov = ui.wrapping_mul(nlu as usize).wrapping_add(lu as usize);
            let frac = if cov < ctx.subcatches.coverage.len() { ctx.subcatches.coverage[cov] / 100.0 } else { 0.0 };
            if frac <= 0.0 {
                continue;
            }
            for m in 0..nm {
                let s = &mut ctx.reactions.surface;
                let k = s.pidx(lu, m);
                let norm = if s.bu_normalizer[k] == 0 {
                    frac * ctx.subcatches.area[ui]
                } else {
                    frac * ctx.subcatches.curb_length[ui]
                };
                let mut mass = 0.0;
                let load0 = s.init_loading[s.sidx(i, m)];
                if load0 > 0.0 {
                    if s.bu_normalizer[k] == 0 {
                        mass = load0;
                    } else {
                        mass = if ctx.subcatches.curb_length[ui] > 0.0 {
                            load0 * ctx.subcatches.area[ui] / ctx.subcatches.curb_length[ui]
                        } else {
                            0.0
                        };
                    }
                } else if dry_days > 0.0 && s.bu_type[k] != 0 && s.bu_type[k] != 4 {
                    mass = buildup_at(s.bu_type[k], s.bu_c1[k], s.bu_c2[k], s.bu_c3[k], dry_days);
                }
                if mass <= 0.0 {
                    continue;
                }
                let b = s.bu_idx(i, lu, m);
                s.buildup[b] = mass;
                s.led_init_buildup[m as usize] += mass * norm;
            }
        }
    }
}

// ============================================================================
// msxsurf::remainingBuildup (MsxSurfaceQuality.cpp). Translated 5 November 2026 for
// SWMMEngine::init_modules, which re-books a hot-started species' initial buildup.
// ============================================================================

/// remainingBuildup(ctx, m): species m's mass still on the land (buildup, ponded, LID
/// units and carry, ponds written off), in user mass units.
pub fn remaining_buildup(ctx: &SimulationContext, m: i32) -> f64 {
    let s = &ctx.reactions.surface;
    if !s.active() || m < 0 || m >= s.n_species {
        return 0.0;
    }
    let mut total = 0.0;
    let nlu = s.n_landuses;
    for sc in 0..ctx.n_subcatches() {
        let usc = sc as usize;
        for lu in 0..nlu {
            let cov = usc * nlu as usize + lu as usize;
            let frac = if cov < ctx.subcatches.coverage.len() { ctx.subcatches.coverage[cov] / 100.0 } else { 0.0 };
            if frac <= 0.0 {
                continue;
            }
            let norm = if s.bu_normalizer[s.pidx(lu, m)] == 0 { frac * ctx.subcatches.area[usc] } else { frac * ctx.subcatches.curb_length[usc] };
            let k = s.bu_idx(sc, lu, m);
            if k < s.buildup.len() {
                total += s.buildup[k] * norm;
            }
        }
        let pk = s.sidx(sc, m);
        if pk < s.ponded_qual.len() {
            total += s.ponded_qual[pk] * s.mcf[m as usize];
        }
    }
    let mut k = m as usize;
    while k < ctx.reactions.lid.mass.len() {
        total += ctx.reactions.lid.mass[k] * 28.317 * s.mcf[m as usize];
        k += s.n_species as usize;
    }
    let mut k = m as usize;
    while k < ctx.reactions.lid.surface_carry.len() {
        total += ctx.reactions.lid.surface_carry[k] * 28.317 * s.mcf[m as usize];
        k += s.n_species as usize;
    }
    if (m as usize) < s.led_ponded_lost.len() {
        total += s.led_ponded_lost[m as usize];
    }
    total
}

/// msxsurf::deliver: this routing step's species washoff (and the LID drain species loads)
/// into the external species loads (msx_ext_mass_in) at the outlet nodes, at legacy's
/// interpolation point between the runoff times.
pub fn deliver(ctx: &mut SimulationContext, _dt_routing: f64) {
    if !ctx.reactions.surface.active() {
        return;
    }
    let nm = ctx.reactions.n_species();
    if nm != ctx.reactions.surface.n_species || ctx.n_nodes() <= 0 {
        return;
    }
    let want = ctx.n_nodes() as usize * nm as usize;
    if ctx.reactions.msx_ext_mass_in.len() != want {
        ctx.reactions.msx_ext_mass_in.clear();
        ctx.reactions.msx_ext_mass_in.resize(want, 0.0);
    }
    for i in 0..ctx.n_subcatches() {
        let ui = i as usize;
        let out_node = ctx.subcatches.outlet_node[ui];
        if out_node < 0 || out_node >= ctx.n_nodes() {
            continue;
        }
        let f = ctx.runoff_interp_f;
        let f1 = 1.0 - f;
        let q_old = ctx.subcatches.old_runoff[ui];
        let q_new = ctx.subcatches.runoff[ui];
        let q = f1 * q_old + f * q_new;
        if q <= 0.0 {
            continue;
        }
        let rx = &mut ctx.reactions;
        let s = &rx.surface;
        let has_old = s.washoff_conc_old.len() == s.washoff_conc.len();
        let base = out_node as usize * nm as usize;
        for m in 0..nm {
            let si = s.sidx(i, m);
            let c_new = s.washoff_conc[si];
            let c_old = if has_old { s.washoff_conc_old[si] } else { c_new };
            let mass_rate = f1 * (q_old * c_old) + f * (q_new * c_new);
            if mass_rate > 0.0 {
                rx.msx_ext_mass_in[base + m as usize] += mass_rate;
            }
        }
    }
    if ctx.reactions.surface.lid_drain_load.len() == want {
        let f = ctx.runoff_interp_f;
        let f1 = 1.0 - f;
        let has_old = ctx.reactions.surface.lid_drain_load_old.len() == want;
        for j in 0..ctx.n_nodes() {
            let uj = j as usize;
            let v_new = if uj < ctx.nodes.lid_drain_qual_vol.len() { ctx.nodes.lid_drain_qual_vol[uj] } else { 0.0 };
            let v_old = if uj < ctx.nodes.lid_drain_qual_vol_old.len() { ctx.nodes.lid_drain_qual_vol_old[uj] } else { 0.0 };
            if f1 * v_old + f * v_new <= 0.0 {
                continue;
            }
            let rx = &mut ctx.reactions;
            for m in 0..nm as usize {
                let k = uj * nm as usize + m;
                let w_new = rx.surface.lid_drain_load[k];
                let w_old = if has_old { rx.surface.lid_drain_load_old[k] } else { 0.0 };
                let load = f1 * w_old + f * w_new;
                if load > 0.0 {
                    rx.msx_ext_mass_in[k] += load;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// The runoff-step kernels (6 November 2026): speciesParams, reactPondedPools, reactBuildup,
// step, bookPostLid, addLidDrain, addLidDrainRunon, consumeLidDrainRunon. The species
// Block (speciesBlock) is quality::surface_species_kernel's step_species /
// book_post_lid_species, which the family's MsxSurfaceData is lent to.
//
// Differences of form, never in result:
//   * SurfaceScratchLease's thread_local SurfaceScratch is per-call locals here: every value
//     the kernels read from it is written first in the same call (mixed is assigned, values,
//     before, normalizer, forms and expressions are assigned or filled before use; the
//     workspace is re-initialized by workspace.init at the top of each kernel).
//   * ReactionIntegrator::step's paths that are not translated return StepError (a stub).
// ---------------------------------------------------------------------------

use crate::core::engine_step_routing::StepError;
use crate::data::reaction::ReactionExprForm;
use crate::data::reaction_tokens::RxHydVar;
use crate::quality::landuse::{BuildupParams, BuildupType, WashoffParams, WashoffType};
use crate::quality::surface_species_kernel::{book_post_lid_species, step_species, SpeciesParams, SpeciesReacted};
use crate::transport::reaction_legacy_step::{integrator_step, note_negative_clamp, reaction_area_factor, Ws};

fn species_params(s: &crate::data::msx_surface::MsxSurfaceData) -> SpeciesParams {
    let n = s.bu_type.len();
    let mut p = SpeciesParams { bu: vec![BuildupParams::default(); n], wo: vec![WashoffParams::default(); n] };
    for k in 0..n {
        let bp = &mut p.bu[k];
        bp.type_ = BuildupType::from_i32(s.bu_type[k]);
        bp.coeff[0] = s.bu_c1[k];
        bp.coeff[1] = s.bu_c2[k];
        bp.coeff[2] = s.bu_c3[k];
        bp.normalizer = s.bu_normalizer[k];
        bp.max_days = s.bu_max_days[k];
        let wp = &mut p.wo[k];
        wp.type_ = WashoffType::from_i32(s.wo_type[k]);
        wp.coeff = s.wo_coeff[k];
        wp.expon = s.wo_expon[k];
        wp.sweep_effic = s.wo_sweep_effic[k];
        wp.bmp_effic = s.wo_bmp_effic[k];
    }
    p
}

/// reactPondedPools: react each wet subcatchment's mixed ponded pool; returns the mixed
/// (reacted) masses, sized as ponded_qual.
fn react_ponded_pools(ctx: &mut SimulationContext, dt: f64) -> Result<Vec<f64>, StepError> {
    let ns = ctx.reactions.n_species();
    let mut mixed = vec![0.0f64; ctx.reactions.surface.ponded_qual.len()];
    let mut species = vec![0.0f64; ns as usize];
    let mut before = vec![0.0f64; species.len()];
    let mut workspace = Ws::new(&ctx.reactions);
    let area_factor = ucf(crate::core::unit_conversion::LANDAREA, &ctx.options);
    const LITRES_PER_FT3: f64 = 28.317;
    for sc in 0..ctx.n_subcatches() {
        let i = sc as usize;
        if !ctx.reactions.wet_enabled(sc) {
            continue;
        }
        let mut forms = ctx.reactions.subcatch_form.clone();
        let mut expressions = ctx.reactions.subcatch_expr.clone();
        for law in ctx.reactions.wet_laws.iter() {
            if law.subcatchment == sc {
                forms[law.species as usize] = law.form;
                expressions[law.species as usize] = law.compiled.clone();
            }
        }
        let lid = if i < ctx.subcatches.total_lid_area_ft2.len() { ctx.subcatches.total_lid_area_ft2[i] } else { 0.0 };
        let area = ctx.subcatches.area[i] / area_factor - lid;
        if !(area > 0.0) {
            continue;
        }
        let rain_volume = cmax(0.0, ctx.subcatches.rainfall[i]) * area * dt;
        let volume = ctx.subcatches.old_ponded_depth[i] * area + rain_volume + ctx.subcatches.runon_inflow[i] * dt;
        if !(volume > 0.0) {
            continue;
        }
        let litres = volume * LITRES_PER_FT3;
        for m in 0..ns {
            let k = ctx.reactions.surface.sidx(sc, m);
            let rain = if m < ctx.reactions.rain_conc.len() as i32 { ctx.reactions.rain_conc[m as usize] } else { 0.0 };
            mixed[k] = ctx.reactions.surface.ponded_qual[k] + rain * LITRES_PER_FT3 * rain_volume + ctx.reactions.surface.runon_qual_rate[k] * dt;
        }
        if volume <= crate::quality::quality_routing::ZERO_VOLUME {
            continue;
        }
        for m in 0..ns as usize {
            species[m] = mixed[ctx.reactions.surface.sidx(sc, m as i32)] / litres;
            before[m] = species[m];
        }
        let mut variables = [0.0f64; RxHydVar::Count as usize];
        variables[RxHydVar::Temp as usize] = ctx.reactions.default_temp_c;
        if ctx.options.heat_transport && i < ctx.heat_state.subcatch_runoff_temp.len() && ctx.heat_state.subcatch_runoff_temp[i].is_finite() {
            variables[RxHydVar::Temp as usize] = ctx.heat_state.subcatch_runoff_temp[i];
        }
        let report = {
            let rx = &ctx.reactions;
            integrator_step(rx, &forms, &expressions, dt, &mut species, &mut variables, &mut workspace, None, None, &rx.coef_value)?
        };
        let mut physical = report.ok;
        let mut m = 0;
        while m < ns as usize && physical {
            physical = species[m].is_finite() && (species[m] >= 0.0 || forms[m] == ReactionExprForm::Rate);
            m += 1;
        }
        if !physical {
            if !ctx.reactions.warned_react_failure {
                ctx.reactions.warned_react_failure = true;
                ctx.warnings.push(if report.ok { b"negative algebraic concentration".to_vec() } else { report.error.clone() });
            }
            continue;
        }
        for m in 0..ns as usize {
            let raw = species[m];
            let clamped = cmax(0.0, raw);
            let factor = litres * ctx.reactions.surface.mcf[m];
            ctx.reactions.surface.led_reacted[m] += (before[m] - raw) * factor;
            ctx.reactions.surface.led_reaction_clamped[m] += (clamped - raw) * factor;
            if clamped != raw {
                note_negative_clamp(ctx, b"a wet subcatchment pool");
            }
            mixed[ctx.reactions.surface.sidx(sc, m as i32)] = clamped * litres;
        }
    }
    Ok(mixed)
}

/// reactBuildup: react the dry buildup stores under the [REACTIONS] BUILDUP laws.
fn react_buildup(ctx: &mut SimulationContext, dt: f64) -> Result<(), StepError> {
    if ctx.reactions.buildup_laws.is_empty() {
        return Ok(());
    }
    let ns = ctx.reactions.n_species();
    let uns = ns as usize;
    let mut ws = Ws::new(&ctx.reactions);
    let mut y = vec![0.0f64; uns];
    let mut before = vec![0.0f64; uns];
    let mut normalizer = vec![0.0f64; uns];
    let mut forms = vec![ReactionExprForm::None; uns];
    let mut expr = vec![crate::data::reaction_tokens::RxExprSpan::default(); uns];
    let nlu = ctx.reactions.surface.n_landuses;
    for sc in 0..ctx.n_subcatches() {
        for lu in 0..nlu {
            let cov = sc as usize * nlu as usize + lu as usize;
            let fraction = if cov < ctx.subcatches.coverage.len() { ctx.subcatches.coverage[cov] / 100. } else { 0. };
            let area = fraction * ctx.subcatches.area[sc as usize] / ucf(crate::core::unit_conversion::LANDAREA, &ctx.options) * reaction_area_factor(ctx.reactions.area_units);
            if !(area > 0.0) {
                continue;
            }
            for f in forms.iter_mut() {
                *f = ReactionExprForm::None;
            }
            let mut any = false;
            for selected in [-1, sc] {
                for l in ctx.reactions.buildup_laws.iter() {
                    if l.landuse == lu && l.subcatchment == selected {
                        forms[l.species as usize] = l.form;
                        expr[l.species as usize] = l.compiled.clone();
                        any = true;
                    }
                }
            }
            if !any {
                continue;
            }
            let mut supported = true;
            {
                let s = &ctx.reactions.surface;
                for m in 0..ns {
                    let um = m as usize;
                    let norm = fraction * (if s.bu_normalizer[s.pidx(lu, m)] == 0 { ctx.subcatches.area[sc as usize] } else { ctx.subcatches.curb_length[sc as usize] });
                    normalizer[um] = norm;
                    y[um] = if norm > 0.0 { s.buildup[s.bu_idx(sc, lu, m)] * norm / s.mcf[um] / area } else { 0. };
                    before[um] = y[um];
                    if forms[um] != ReactionExprForm::None && !(norm > 0.0) {
                        supported = false;
                    }
                }
            }
            let mut vars = [0.0f64; RxHydVar::Count as usize];
            vars[RxHydVar::Temp as usize] = ctx.reactions.default_temp_c;
            let report = if supported {
                let rx = &ctx.reactions;
                integrator_step(rx, &forms, &expr, dt, &mut y, &mut vars, &mut ws, None, None, &rx.coef_value)?
            } else {
                crate::transport::reaction_legacy_step::StepReport { ok: true, error: Vec::new() }
            };
            let mut physical = supported && report.ok;
            let mut m = 0;
            while m < uns && physical {
                physical = y[m].is_finite() && (y[m] >= 0.0 || forms[m] == ReactionExprForm::Rate);
                m += 1;
            }
            if !physical {
                if !ctx.reactions.warned_react_failure {
                    ctx.reactions.warned_react_failure = true;
                    ctx.warnings.push(b"Dry buildup reaction failed or has a zero normalizer; original inventory retained.".to_vec());
                }
                continue;
            }
            for m in 0..ns {
                let um = m as usize;
                let bounded = cmax(0., y[um]);
                let factor = area * ctx.reactions.surface.mcf[um];
                ctx.reactions.surface.led_reacted[um] += (before[um] - y[um]) * factor;
                ctx.reactions.surface.led_reaction_clamped[um] += (bounded - y[um]) * factor;
                if y[um] != bounded {
                    note_negative_clamp(ctx, b"a dry buildup store");
                }
                if normalizer[um] > 0.0 {
                    let k = ctx.reactions.surface.bu_idx(sc, lu, m);
                    ctx.reactions.surface.buildup[k] = bounded * factor / normalizer[um];
                }
            }
        }
    }
    Ok(())
}

/// msxsurf::step: dry buildup reactions, wet ponded-pool reactions, then the species
/// family's buildup/washoff/ponded step.
pub fn step(ctx: &mut SimulationContext, dt_runoff: f64, old_runoff_time: f64) -> Result<(), StepError> {
    if !ctx.reactions.surface.active() || !(dt_runoff > 0.0) {
        return Ok(());
    }
    let params = species_params(&ctx.reactions.surface);
    react_buildup(ctx, dt_runoff)?;
    let mut reacted: Vec<f64> = Vec::new();
    let mut wet_mask = vec![0u8; ctx.n_subcatches().max(0) as usize];
    let has = ctx.reactions.has_subcatch_reactions();
    if has {
        reacted = react_ponded_pools(ctx, dt_runoff)?;
        for sc in 0..ctx.n_subcatches() {
            wet_mask[sc as usize] = ctx.reactions.wet_enabled(sc) as u8;
        }
    }
    let mut s = std::mem::take(&mut ctx.reactions.surface);
    let r = SpeciesReacted { ponded_reacted_mass: &reacted, ponded_reactor_mask: &wet_mask };
    step_species(ctx, &mut s, &params, if has { Some(&r) } else { None }, dt_runoff, old_runoff_time);
    ctx.reactions.surface = s;
    Ok(())
}

/// msxsurf::bookPostLid.
pub fn book_post_lid(ctx: &mut SimulationContext, dt_runoff: f64, lid_infil_vol: &[f64]) {
    if !ctx.reactions.surface.active() || ctx.n_subcatches() <= 0 {
        return;
    }
    let mut s = std::mem::take(&mut ctx.reactions.surface);
    book_post_lid_species(ctx, &mut s, dt_runoff, lid_infil_vol);
    ctx.reactions.surface = s;
}

/// msxsurf::addLidDrain: a LID drain's species load to its outlet node (source: the LID
/// unit's own concentrations, or None for the subcatchment's washoff).
#[allow(clippy::too_many_arguments)]
pub fn add_lid_drain(ctx: &mut SimulationContext, control: i32, sc: i32, node: i32, drain_cfs: f64, to_perv: bool, runoff_load: bool, dt_runoff: f64, source: Option<&[f64]>) {
    let n_nodes = ctx.n_nodes();
    let s = &mut ctx.reactions.surface;
    if !s.active() || node < 0 || node >= n_nodes {
        return;
    }
    let nm = s.n_species;
    let want = n_nodes as usize * nm as usize;
    if s.lid_drain_load.len() != want {
        s.lid_drain_load.clear();
        s.lid_drain_load.resize(want, 0.0);
    }
    for m in 0..nm {
        let um = m as usize;
        let rk = (control as usize).wrapping_mul(nm as usize).wrapping_add(um);
        let rmvl = if control >= 0 && rk < s.lid_rmvl.len() { s.lid_rmvl[rk] } else { 0.0 };
        let c_src = match source {
            Some(src) => src[um],
            None => s.washoff_conc[s.sidx(sc, m)],
        };
        s.lid_drain_load[node as usize * nm as usize + um] += drain_cfs * c_src * (1.0 - rmvl);
        if !to_perv {
            let w = drain_cfs * c_src * dt_runoff * 28.317 * s.mcf[um];
            if source.is_some() {
                let k = s.sidx(sc, m);
                s.led_subcatch_load[k] += w;
            }
            if w > 0.0 {
                s.led_bmp_removal[um] += rmvl * w;
                if runoff_load {
                    s.led_runoff_load[um] += w * (1.0 - rmvl);
                }
            }
        }
    }
}

/// msxsurf::addLidDrainRunon: a LID drain's species load routed onto a receiving subcatchment.
#[allow(clippy::too_many_arguments)]
pub fn add_lid_drain_runon(ctx: &mut SimulationContext, control: i32, sc: i32, receiver: i32, drain_cfs: f64, to_perv: bool, runoff_load: bool, dt_runoff: f64, source: Option<&[f64]>) {
    let n_sub = ctx.n_subcatches();
    let s = &mut ctx.reactions.surface;
    if !s.active() || receiver < 0 || receiver >= n_sub {
        return;
    }
    let nm = s.n_species;
    for m in 0..nm {
        let um = m as usize;
        let rk = (control as usize).wrapping_mul(nm as usize).wrapping_add(um);
        let rmvl = if control >= 0 && rk < s.lid_rmvl.len() { s.lid_rmvl[rk] } else { 0.0 };
        let c_src = match source {
            Some(src) => src[um],
            None => s.washoff_conc[s.sidx(sc, m)],
        };
        let k = s.sidx(receiver, m);
        if k < s.lid_drain_runon_load.len() {
            s.lid_drain_runon_load[k] += (drain_cfs * c_src * 28.317) * (1.0 - rmvl);
        }
        if !to_perv {
            let w = drain_cfs * c_src * dt_runoff * 28.317 * s.mcf[um];
            if source.is_some() {
                let k = s.sidx(sc, m);
                s.led_subcatch_load[k] += w;
            }
            if w > 0.0 {
                s.led_bmp_removal[um] += rmvl * w;
                if runoff_load {
                    s.led_runoff_load[um] += w * (1.0 - rmvl);
                }
            }
        }
    }
}

/// msxsurf::consumeLidDrainRunon.
pub fn consume_lid_drain_runon(ctx: &mut SimulationContext, receiver: i32) {
    let s = &mut ctx.reactions.surface;
    if !s.active() {
        return;
    }
    for m in 0..s.n_species {
        let k = s.sidx(receiver, m);
        if k < s.lid_drain_runon_load.len() && k < s.runon_qual_rate.len() {
            s.runon_qual_rate[k] += s.lid_drain_runon_load[k];
            s.lid_drain_runon_load[k] = 0.0;
        }
    }
}
