// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/quality/SurfaceQuality2D.hpp and
// SurfaceQuality2D.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// S7 cell buildup / washoff / sweeping: the [2D_COVERAGES] / [2D_LOADINGS] / [2D_CURB_LENGTH]
// rows (the open-time part: the data, the scope and row parsers, the registration). The
// run-time part (resolve, initState, step, loadParams) belongs to a later 2D wave; its
// members are here with upstream's names (pub, so the gate dumper can read them).
//   * num() is std::stod consuming the whole token (core::cstdlib::stod_pos).
//   * std::to_string(double) is "%f" (core::cprintf::to_string_f64).

#![allow(clippy::all)]

use crate::core::context::SimulationContext;
use crate::core::cprintf::to_string_f64;
use crate::core::cstdlib::stod_pos;
use crate::core::options::Bytes;
use crate::input::section_registry::SectionRegistry;
use crate::input::tokenizer;
use crate::twod::Shared;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i8)]
pub enum SqScope {
    #[default]
    Global = 0,
    Tag = 1,
    Cell = 2,
}

#[derive(Debug, Clone)]
pub struct SqScopeKey {
    pub scope: SqScope,
    pub tag: Bytes,
    pub cell: i32,
}

impl Default for SqScopeKey {
    fn default() -> Self {
        SqScopeKey { scope: SqScope::Global, tag: Bytes::new(), cell: -1 }
    }
}

#[derive(Debug, Clone, Default)]
pub struct CoverageRow2D {
    pub key: SqScopeKey,
    pub uses: Vec<(Bytes, f64)>,
}

#[derive(Debug, Clone, Default)]
pub struct LoadingRow2D {
    pub key: SqScopeKey,
    pub species: Bytes,
    pub value: f64,
}

#[derive(Debug, Clone, Default)]
pub struct CurbRow2D {
    pub key: SqScopeKey,
    pub length: f64,
}

/// SurfaceQuality2D::Params (private upstream)
#[derive(Debug, Clone)]
pub struct SqParams {
    pub bu_type: i32,
    pub c0: f64,
    pub c1: f64,
    pub c2: f64,
    pub ts_idx: i32,
    pub normalizer: i32,
    pub max_days: f64,
    pub wo_type: i32,
    pub coeff: f64,
    pub expon: f64,
    pub sweep_effic: f64,
    pub bmp_effic: f64,
}

impl Default for SqParams {
    fn default() -> Self {
        SqParams {
            bu_type: 0,
            c0: 0.0,
            c1: 0.0,
            c2: 0.0,
            ts_idx: -1,
            normalizer: 0,
            max_days: 0.0,
            wo_type: 0,
            coeff: 0.0,
            expon: 0.0,
            sweep_effic: 0.0,
            bmp_effic: 0.0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct SurfaceQuality2D {
    pub coverage_rows: Vec<CoverageRow2D>,
    pub loading_rows: Vec<LoadingRow2D>,
    pub curb_rows: Vec<CurbRow2D>,
    // ---- run-time members (private upstream)
    pub active_: bool,
    pub n_cells_: i32,
    pub n_lu_: i32,
    pub n_sp_: i32,
    pub si_: bool,
    pub area_to_landarea_: f64,
    pub row_of_: Vec<i32>,
    pub names_: Vec<Bytes>,
    pub mcf_: Vec<f64>,
    pub params_: Vec<SqParams>,
    pub cov_: Vec<f64>,
    pub curb_: Vec<f64>,
    pub area_la_: Vec<f64>,
    pub init_: Vec<f64>,
    pub bu_: Vec<f64>,
    pub last_swept_: Vec<f64>,
    pub runoff_vol_: Vec<f64>,
    pub wet_runoff_vol_: Vec<f64>,
    pub led_init_: Vec<f64>,
    pub led_bu_: Vec<f64>,
    pub led_wo_: Vec<f64>,
    pub led_sw_: Vec<f64>,
    pub led_bmp_: Vec<f64>,
}

impl Default for SurfaceQuality2D {
    fn default() -> Self {
        SurfaceQuality2D {
            coverage_rows: Vec::new(),
            loading_rows: Vec::new(),
            curb_rows: Vec::new(),
            active_: false,
            n_cells_: 0,
            n_lu_: 0,
            n_sp_: 0,
            si_: false,
            area_to_landarea_: 1.0,
            row_of_: Vec::new(),
            names_: Vec::new(),
            mcf_: Vec::new(),
            params_: Vec::new(),
            cov_: Vec::new(),
            curb_: Vec::new(),
            area_la_: Vec::new(),
            init_: Vec::new(),
            bu_: Vec::new(),
            last_swept_: Vec::new(),
            runoff_vol_: Vec::new(),
            wet_runoff_vol_: Vec::new(),
            led_init_: Vec::new(),
            led_bu_: Vec::new(),
            led_wo_: Vec::new(),
            led_sw_: Vec::new(),
            led_bmp_: Vec::new(),
        }
    }
}

impl SurfaceQuality2D {
    pub fn authored(&self) -> bool {
        !self.coverage_rows.is_empty() || !self.loading_rows.is_empty() || !self.curb_rows.is_empty()
    }
    pub fn clear_rows(&mut self) {
        self.coverage_rows.clear();
        self.loading_rows.clear();
        self.curb_rows.clear();
    }
    pub fn active(&self) -> bool {
        self.active_
    }
    pub fn n_land_uses(&self) -> i32 {
        self.n_lu_
    }
    pub fn n_species(&self) -> i32 {
        self.n_sp_
    }
    pub fn n_cells(&self) -> i32 {
        self.n_cells_
    }
}

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

/// num(): std::stod consuming the whole token; `v` keeps its value where stod throws.
fn num(s: &[u8], v: &mut f64) -> bool {
    match stod_pos(s) {
        Ok((x, n)) => {
            *v = x;
            n == s.len()
        }
        Err(_) => false,
    }
}

/// formatSqScope()
pub fn format_sq_scope(k: &SqScopeKey) -> Bytes {
    match k.scope {
        SqScope::Tag => cat(&[b"TAG ", &k.tag]),
        SqScope::Cell => cat(&[b"CELL ", (k.cell + 1).to_string().as_bytes()]),
        _ => b"*".to_vec(),
    }
}

/// parseSqScope()
pub fn parse_sq_scope(tokens: &[Bytes], at: &mut usize, key: &mut SqScopeKey, section: &[u8]) -> Bytes {
    if tokens.is_empty() {
        return cat(&[section, b" empty row"]);
    }
    if tokens[0] == b"*" {
        key.scope = SqScope::Global;
        *at = 1;
        return Bytes::new();
    }
    if iequals(&tokens[0], b"TAG") {
        if tokens.len() < 2 {
            return cat(&[section, b" TAG needs a name"]);
        }
        key.scope = SqScope::Tag;
        key.tag = tokens[1].clone();
        *at = 2;
        return Bytes::new();
    }
    if iequals(&tokens[0], b"CELL") {
        if tokens.len() < 2 {
            return cat(&[section, b" CELL needs an index"]);
        }
        let mut d = 0.0;
        if !num(&tokens[1], &mut d) || d < 1.0 || d != d.floor() {
            return cat(&[section, b" invalid CELL index (1-based): ", &tokens[1]]);
        }
        key.scope = SqScope::Cell;
        key.cell = crate::core::cmath::cvt_i32(d).wrapping_sub(1); // INT_MIN - 1 wraps, as on x86
        *at = 2;
        return Bytes::new();
    }
    cat(&[section, b" row must start with CELL, TAG or '*': ", &tokens[0]])
}

/// parse2DCoveragesLine(): a row replaces any earlier row of the same scope.
pub fn parse_2d_coverages_line(tokens: &[Bytes], rows: &mut Vec<CoverageRow2D>) -> Bytes {
    let mut r = CoverageRow2D::default();
    let mut at = 0usize;
    let err = parse_sq_scope(tokens, &mut at, &mut r.key, b"[2D_COVERAGES]");
    if !err.is_empty() {
        return err;
    }
    if at >= tokens.len() || (tokens.len() - at) % 2 != 0 {
        return b"[2D_COVERAGES] expects LANDUSE PERCENT pairs after the scope".to_vec();
    }
    let mut sum = 0.0;
    while at + 1 < tokens.len() {
        let mut pct = 0.0;
        if !num(&tokens[at + 1], &mut pct) || pct < 0.0 {
            return cat(&[b"[2D_COVERAGES] invalid percent for '", &tokens[at], b"': ", &tokens[at + 1]]);
        }
        sum += pct;
        r.uses.push((tokens[at].clone(), pct));
        at += 2;
    }
    if sum > 100.0 + 1e-9 {
        return cat(&[b"[2D_COVERAGES] percents for scope '", &format_sq_scope(&r.key), b"' sum to ", &to_string_f64(sum), b" (> 100)"]);
    }
    rows.retain(|o| !(o.key.scope == r.key.scope && o.key.tag == r.key.tag && o.key.cell == r.key.cell));
    rows.push(r);
    Bytes::new()
}

/// parse2DLoadingsLine()
pub fn parse_2d_loadings_line(tokens: &[Bytes], rows: &mut Vec<LoadingRow2D>) -> Bytes {
    let mut r = LoadingRow2D::default();
    let mut at = 0usize;
    let err = parse_sq_scope(tokens, &mut at, &mut r.key, b"[2D_LOADINGS]");
    if !err.is_empty() {
        return err;
    }
    if at + 2 != tokens.len() {
        return b"[2D_LOADINGS] expects SPECIES VALUE after the scope".to_vec();
    }
    r.species = tokens[at].clone();
    if !num(&tokens[at + 1], &mut r.value) || r.value < 0.0 {
        return cat(&[b"[2D_LOADINGS] invalid value for '", &r.species, b"': ", &tokens[at + 1]]);
    }
    rows.push(r);
    Bytes::new()
}

/// parse2DCurbLengthLine()
pub fn parse_2d_curb_length_line(tokens: &[Bytes], rows: &mut Vec<CurbRow2D>) -> Bytes {
    let mut r = CurbRow2D::default();
    let mut at = 0usize;
    let err = parse_sq_scope(tokens, &mut at, &mut r.key, b"[2D_CURB_LENGTH]");
    if !err.is_empty() {
        return err;
    }
    if at + 1 != tokens.len() {
        return b"[2D_CURB_LENGTH] expects LENGTH after the scope".to_vec();
    }
    if !num(&tokens[at], &mut r.length) || r.length < 0.0 {
        return cat(&[b"[2D_CURB_LENGTH] invalid length: ", &tokens[at]]);
    }
    rows.push(r);
    Bytes::new()
}

fn make_sq_handler(sq: Shared<SurfaceQuality2D>, parser: fn(&mut SurfaceQuality2D, &[Bytes]) -> Bytes) -> impl Fn(&mut SimulationContext, &[Bytes]) + 'static {
    move |ctx: &mut SimulationContext, lines: &[Bytes]| {
        for raw in lines {
            let tokens = tokenizer::tokenize(raw);
            if tokens.is_empty() {
                continue;
            }
            let err = parser(&mut *sq.borrow_mut(), &tokens);
            if !err.is_empty() {
                ctx.error_code = 5; // SWMM_ERR_PARSE
                ctx.error_message = cat(&[b"[2D] ", &err, " \u{2014} line: ".as_bytes(), raw]);
                return;
            }
        }
    }
}

/// registerSurfaceQualitySections()
pub fn register_surface_quality_sections(sq: &Shared<SurfaceQuality2D>, registry: &mut SectionRegistry) {
    registry.register_custom(b"2D_COVERAGES", make_sq_handler(sq.clone(), |s, t| parse_2d_coverages_line(t, &mut s.coverage_rows)));
    registry.register_custom(b"2D_LOADINGS", make_sq_handler(sq.clone(), |s, t| parse_2d_loadings_line(t, &mut s.loading_rows)));
    registry.register_custom(b"2D_CURB_LENGTH", make_sq_handler(sq.clone(), |s, t| parse_2d_curb_length_line(t, &mut s.curb_rows)));
}
