// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/infil/Infil2D.hpp and Infil2D.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Per-cell 2D infiltration: the authored [2D_INFILTRATION_OPTIONS] / [2D_INFILTRATION_DEFAULTS]
// / [2D_INFILTRATION] rows (filled at open) and the members resolve() fills. Translated so far:
// the data, the token helpers (parseInfil2DMethod, infil2DMethodToken, parseInfil2DDest,
// infil2DDestToken, infil2DParamCount) and reset(). Not yet: resolve/resolveImpl and
// updateRates, which run from SurfaceRouter2D::initialize and the step (later waves).

#![allow(clippy::all)]

use crate::core::options::Bytes;
use crate::hydrology::infil_bank::InfilBank;
use crate::hydrology::infiltration::InfilModel;

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
