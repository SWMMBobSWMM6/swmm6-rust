// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/infil/Infil2D.hpp and Infil2D.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Per-cell 2D infiltration: the authored [2D_INFILTRATION_OPTIONS] / [2D_INFILTRATION_DEFAULTS]
// / [2D_INFILTRATION] rows (filled at open) and the members resolve() fills: the data, the token
// helpers (parseInfil2DMethod, infil2DMethodToken, parseInfil2DDest, infil2DDestToken,
// infil2DParamCount) and reset() (wave 1); resolve/resolveImpl, validateRow and updateRates
// (wave 3, 12 November 2026), which SurfaceRouter2D::initialize runs. updateRates' OpenMP loop
// (n >= 512 cells, more than one thread) writes only each cell's own rate and advances only
// each cell's own kernel, so the serial loop gives the same bits.

#![allow(clippy::all)]

use crate::core::options::{Bytes, SimulationOptions};
use crate::hydrology::infil_bank::{Factors, InfilBank, Owner};
use crate::hydrology::infiltration::InfilModel;
use crate::twod::data::mesh_data::MeshData;
use crate::twod::data::surface_state_data::SurfaceStateData;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i32)]
pub enum Infil2DDest {
    #[default]
    Lost = 0,
    SubcatchAquifer = 1,
    Aquifer2D = 2,
}

pub const K_INFIL_2D_MAX_PARAMS: i32 = 5;

#[derive(Debug, Clone)]
pub struct Infil2DRow {
    pub has_method: bool,
    pub method: InfilModel,
    pub p: [f64; K_INFIL_2D_MAX_PARAMS as usize],
    pub dest: Infil2DDest,
    pub dest_explicit: bool,
}

impl Default for Infil2DRow {
    fn default() -> Self {
        Infil2DRow { has_method: false, method: InfilModel::Horton, p: [0.0; 5], dest: Infil2DDest::Lost, dest_explicit: false }
    }
}

#[derive(Debug, Clone, Default)]
pub struct Infil2DDefault {
    pub tag: Bytes,
    pub row: Infil2DRow,
}

#[derive(Debug, Clone)]
pub struct Infil2DOverride {
    pub tri: i32,
    pub row: Infil2DRow,
}

impl Default for Infil2DOverride {
    fn default() -> Self {
        Infil2DOverride { tri: -1, row: Infil2DRow::default() }
    }
}

#[derive(Debug, Clone, Default)]
pub struct Infil2DOptions {
    pub infil_step: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(u8)]
pub enum Infil2DProvenance {
    #[default]
    None = 0,
    Star = 1,
    Tag = 2,
    Override = 3,
}

/// Infil2D. Upstream's members are private; here they are pub (as the other solver structs
/// of the crate) so the gate dumper can read them. Names and order are upstream's.
#[derive(Debug, Clone, Default)]
pub struct Infil2D {
    pub aquifer_owners_: Vec<u8>,
    pub ownership_messages_: Vec<Bytes>,
    pub defaults_: Vec<Infil2DDefault>,
    pub overrides_: Vec<Infil2DOverride>,
    pub options_: Infil2DOptions,
    pub active_: bool,
    pub step_seconds_: f64,
    pub aquifer_2d_available_: bool,
    pub resolved_: Vec<Infil2DRow>,
    pub prov_: Vec<Infil2DProvenance>,
    pub bank_: InfilBank,
}

impl Infil2D {
    pub fn defaults(&self) -> &Vec<Infil2DDefault> {
        &self.defaults_
    }
    pub fn defaults_mut(&mut self) -> &mut Vec<Infil2DDefault> {
        &mut self.defaults_
    }
    pub fn overrides(&self) -> &Vec<Infil2DOverride> {
        &self.overrides_
    }
    pub fn overrides_mut(&mut self) -> &mut Vec<Infil2DOverride> {
        &mut self.overrides_
    }
    pub fn options(&self) -> &Infil2DOptions {
        &self.options_
    }
    pub fn options_mut(&mut self) -> &mut Infil2DOptions {
        &mut self.options_
    }
    pub fn set_aquifer_owners(&mut self, owners: Vec<u8>) {
        self.aquifer_owners_ = owners;
    }
    pub fn ownership_messages(&self) -> &Vec<Bytes> {
        &self.ownership_messages_
    }
    pub fn set_aquifer_2d_available(&mut self, on: bool) {
        self.aquifer_2d_available_ = on;
    }
    pub fn aquifer_2d_available(&self) -> bool {
        self.aquifer_2d_available_
    }
    pub fn active(&self) -> bool {
        self.active_
    }
    pub fn deactivate(&mut self) {
        self.active_ = false;
    }
    pub fn resolved_rows(&self) -> &Vec<Infil2DRow> {
        &self.resolved_
    }
    pub fn provenance(&self) -> &Vec<Infil2DProvenance> {
        &self.prov_
    }
    pub fn step_seconds(&self) -> f64 {
        self.step_seconds_
    }

    /// reset(): back to an unauthored, inactive owner.
    pub fn cumulative(&self) -> &Vec<f64> {
        self.bank_.cumulative()
    }

    /// resolve(mesh, opts, err): resolve into a candidate and adopt it only on success.
    pub fn resolve(&mut self, mesh: &MeshData, opts: &SimulationOptions, err: &mut Bytes) -> bool {
        let mut candidate = self.clone();
        if !candidate.resolve_impl(mesh, opts, err) {
            return false;
        }
        *self = candidate;
        true
    }

    fn resolve_impl(&mut self, mesh: &MeshData, opts: &SimulationOptions, err: &mut Bytes) -> bool {
        let nt = mesh.n_triangles();
        let nt_u = nt as usize;
        self.resolved_.clear();
        self.resolved_.resize(nt_u, Infil2DRow::default());
        self.prov_.clear();
        self.prov_.resize(nt_u, Infil2DProvenance::None);
        self.bank_.init(nt);
        self.active_ = false;
        self.ownership_messages_.clear();
        if self.aquifer_owners_.is_empty() {
            let v = if self.aquifer_2d_available_ { 1 } else { 0 };
            self.aquifer_owners_.resize(nt_u, v);
        }
        if self.aquifer_owners_.len() != nt_u {
            *err = b"2D infiltration: invalid aquifer coverage size".to_vec();
            return false;
        }
        let owners = self.aquifer_owners_.clone();
        let owned = |i: usize| owners[i] != 0;
        for d in &self.defaults_ {
            let mut applied = 0;
            let mut skipped = 0;
            for i in 0..nt_u {
                if d.tag.as_slice() != b"*" && (i >= mesh.tri_tag.len() || mesh.tri_tag[i] != d.tag) {
                    continue;
                }
                if owned(i) {
                    skipped += 1;
                } else {
                    applied += 1;
                }
            }
            let mut m = b"2D infiltration default '".to_vec();
            m.extend_from_slice(&d.tag);
            m.extend_from_slice(format!("': {} applied, {} skipped (aquifer-owned)", applied, skipped).as_bytes());
            self.ownership_messages_.push(m);
        }
        // D-I1: the cadence is INFIL_STEP, falling back to the project WET_STEP (seconds).
        self.step_seconds_ = if self.options_.infil_step > 0.0 { self.options_.infil_step } else { opts.wet_step };

        // --- D-I3 precedence, least specific first: '*' -> tag -> override ---
        let mut star: Option<usize> = None;
        for (k, d) in self.defaults_.iter().enumerate() {
            let mut who = b"tag '".to_vec();
            who.extend_from_slice(&d.tag);
            who.push(b'\'');
            if !validate_row(&d.row, &who, err) {
                return false;
            }
            if d.tag.as_slice() == b"*" {
                star = Some(k);
            }
        }
        if let Some(sk) = star {
            if self.defaults_[sk].row.has_method {
                let row = self.defaults_[sk].row.clone();
                for i in 0..nt_u {
                    if owned(i) {
                        continue;
                    }
                    self.resolved_[i] = row.clone();
                    self.prov_[i] = Infil2DProvenance::Star;
                }
            }
        }
        // A tag row spelled NONE deliberately clears the '*' default for its cells.
        let n_tagged = std::cmp::min(nt_u, mesh.tri_tag.len());
        for d in &self.defaults_ {
            if d.tag.as_slice() == b"*" {
                continue;
            }
            for i in 0..n_tagged {
                if mesh.tri_tag[i] != d.tag || owned(i) {
                    continue;
                }
                if d.row.has_method {
                    self.resolved_[i] = d.row.clone();
                    self.prov_[i] = Infil2DProvenance::Tag;
                } else {
                    self.resolved_[i] = Infil2DRow::default();
                    self.prov_[i] = Infil2DProvenance::None;
                }
            }
        }
        for o in &self.overrides_ {
            if o.tri < 0 || o.tri >= nt {
                *err = format!("2D infiltration cell {}: triangle index out of range (mesh has {} triangles)", o.tri + 1, nt).into_bytes();
                return false;
            }
            let who = format!("cell {}", o.tri + 1).into_bytes();
            if !validate_row(&o.row, &who, err) {
                return false;
            }
            let ui = o.tri as usize;
            if owned(ui) && o.row.has_method {
                *err = format!("2D infiltration cell {}: explicit surface method conflicts with aquifer ownership; remove this row", o.tri + 1).into_bytes();
                return false;
            }
            if o.row.has_method {
                self.resolved_[ui] = o.row.clone();
                self.prov_[ui] = Infil2DProvenance::Override;
            } else {
                self.resolved_[ui] = Infil2DRow::default();
                self.prov_[ui] = Infil2DProvenance::None;
            }
        }

        // --- kernel state, from the PROJECT-UNIT parameters ---
        for i in 0..nt_u {
            if owned(i) {
                self.bank_.set_owner(i as i32, Owner::External);
                self.active_ = true;
                continue;
            }
            if !self.resolved_[i].has_method {
                continue;
            }
            self.active_ = true;
            let (method, p) = (self.resolved_[i].method, self.resolved_[i].p);
            self.bank_.set_method(i as i32, method as i32, &p, opts);
        }
        // Nothing resolved: drop back to the unconfigured fast path.
        if !self.active_ {
            self.resolved_.clear();
            self.prov_.clear();
            self.bank_.clear();
        }
        true
    }

    /// updateRates(mesh, state, dt, factors, num_threads): advance every model-carrying
    /// cell's kernel by dt and publish its held rate into state.infil_rate.
    pub fn update_rates(&mut self, _mesh: &MeshData, state: &mut SurfaceStateData, dt: f64, factors: Factors, _num_threads: i32) {
        if !self.active_ || dt <= 0.0 {
            return;
        }
        let n = std::cmp::min(self.resolved_.len(), state.infil_rate.len());
        for i in 0..n {
            if !self.resolved_[i].has_method {
                continue;
            }
            if self.bank_.owner(i as i32) == Owner::External {
                continue;
            }
            let precip = state.liquid_rainfall(i as i32);
            let depth = state.depth[i];
            state.infil_rate[i] = self.bank_.rate(i as i32, precip, 0.0, depth, dt, factors);
        }
    }

    pub fn reset(&mut self) {
        self.defaults_.clear();
        self.overrides_.clear();
        self.aquifer_owners_.clear();
        self.ownership_messages_.clear();
        self.options_ = Infil2DOptions::default();
        self.active_ = false;
        self.step_seconds_ = 0.0;
        self.aquifer_2d_available_ = false;
        self.resolved_.clear();
        self.prov_.clear();
        self.bank_.clear();
    }
}

/// iequals() of Infil2D.cpp: equal lengths and equal bytes under C-locale toupper.
/// validateRow(row, who, aquifer_2d, err).
fn validate_row(row: &Infil2DRow, who: &[u8], err: &mut Bytes) -> bool {
    let fail = |err: &mut Bytes, tail: &str| {
        let mut m = b"2D infiltration ".to_vec();
        m.extend_from_slice(who);
        m.extend_from_slice(tail.as_bytes());
        *err = m;
        false
    };
    if row.dest == Infil2DDest::Aquifer2D {
        return fail(err, ": destination AQUIFER_2D is obsolete; remove the surface infiltration row on aquifer-owned cells. The aquifer now computes its own receiving capacity.");
    }
    if !row.has_method {
        return true;
    }
    match row.method {
        InfilModel::Horton | InfilModel::ModHorton => {
            if row.p[0] < 0.0 || row.p[1] < 0.0 {
                return fail(err, ": Horton rates f0 and fmin must be >= 0");
            }
        }
        InfilModel::GreenAmpt | InfilModel::ModGreenAmpt => {
            if row.p[1] < 0.0 {
                return fail(err, ": Green-Ampt conductivity Ks must be >= 0");
            }
            if row.p[2] < 0.0 || row.p[2] > 1.0 {
                return fail(err, ": Green-Ampt initial moisture deficit IMD must be in 0..1");
            }
        }
        InfilModel::CurveNum => {
            if row.p[0] < 1.0 || row.p[0] > 100.0 {
                return fail(err, ": curve number must be in 1..100");
            }
        }
        InfilModel::Constant => {
            if row.p[0] < 0.0 {
                return fail(err, ": constant rate must be >= 0");
            }
        }
    }
    true
}

fn iequals(a: &[u8], b: &[u8]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_ascii_uppercase() == y.to_ascii_uppercase())
}

/// parseInfil2DMethod(): NONE clears has_method; a known method sets it; an unknown token
/// clears it and returns false (`method` keeps its value).
pub fn parse_infil_2d_method(token: &[u8], method: &mut InfilModel, has_method: &mut bool) -> bool {
    if iequals(token, b"NONE") {
        *method = InfilModel::Horton; // meaningless when has_method is false
        *has_method = false;
        return true;
    }
    *has_method = true;
    if iequals(token, b"HORTON") {
        *method = InfilModel::Horton;
        return true;
    }
    if iequals(token, b"MODIFIED_HORTON") {
        *method = InfilModel::ModHorton;
        return true;
    }
    if iequals(token, b"GREEN_AMPT") {
        *method = InfilModel::GreenAmpt;
        return true;
    }
    if iequals(token, b"MODIFIED_GREEN_AMPT") {
        *method = InfilModel::ModGreenAmpt;
        return true;
    }
    if iequals(token, b"CURVE_NUMBER") {
        *method = InfilModel::CurveNum;
        return true;
    }
    if iequals(token, b"CONSTANT") {
        *method = InfilModel::Constant;
        return true;
    }
    *has_method = false;
    false
}

/// infil2DMethodToken()
pub fn infil_2d_method_token(row: &Infil2DRow) -> &'static [u8] {
    if !row.has_method {
        return b"NONE";
    }
    match row.method {
        InfilModel::Horton => b"HORTON",
        InfilModel::ModHorton => b"MODIFIED_HORTON",
        InfilModel::GreenAmpt => b"GREEN_AMPT",
        InfilModel::ModGreenAmpt => b"MODIFIED_GREEN_AMPT",
        InfilModel::CurveNum => b"CURVE_NUMBER",
        InfilModel::Constant => b"CONSTANT",
    }
}

/// parseInfil2DDest()
pub fn parse_infil_2d_dest(token: &[u8], dest: &mut Infil2DDest) -> bool {
    if iequals(token, b"LOST") {
        *dest = Infil2DDest::Lost;
        return true;
    }
    if iequals(token, b"SUBCATCH_AQUIFER") {
        *dest = Infil2DDest::SubcatchAquifer;
        return true;
    }
    if iequals(token, b"AQUIFER_2D") {
        *dest = Infil2DDest::Aquifer2D;
        return true;
    }
    false
}

/// infil2DDestToken()
pub fn infil_2d_dest_token(dest: Infil2DDest) -> &'static [u8] {
    match dest {
        Infil2DDest::Lost => b"LOST",
        Infil2DDest::SubcatchAquifer => b"SUBCATCH_AQUIFER",
        Infil2DDest::Aquifer2D => b"AQUIFER_2D",
    }
}

/// infil2DParamCount()
pub fn infil_2d_param_count(method: InfilModel) -> i32 {
    match method {
        InfilModel::Horton | InfilModel::ModHorton => 5,
        InfilModel::GreenAmpt | InfilModel::ModGreenAmpt => 3,
        InfilModel::CurveNum => 3,
        InfilModel::Constant => 1,
    }
}
