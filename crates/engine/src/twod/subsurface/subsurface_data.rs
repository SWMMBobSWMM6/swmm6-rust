// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/subsurface/SubsurfaceData.hpp and
// SubsurfaceData.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright
// 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The [2D_AQUIFER*] authoring store (SubsurfaceConfig) and its row types, in the user's own
// units, as open() fills them, with the enum token helpers. The kernel's run-time state
// (SubsurfaceState, GwLinkShare, GwSubcatchShare, GwUnitFactors) belongs to phase 2b and is not
// translated here.

#![allow(clippy::all)]

use crate::core::options::Bytes;
use crate::twod::subsurface::surface_ownership::SurfaceOwnerRecord;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i8)]
pub enum SoilChar {
    #[default]
    Russo = 0,
    Gardner = 1,
    BrooksCorey = 2,
    VanGenuchten = 3,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i8)]
pub enum GwClosure {
    #[default]
    Auto = -1,
    ClosedForm = 0,
    Enslaved = 1,
    Sigma = 2,
}

#[derive(Debug, Clone)]
#[allow(non_snake_case)]
pub struct GwNodeBed {
    pub node: i32,
    pub cell: i32,
    pub Kc: f64,
    pub dC: f64,
    pub area: f64,
    pub locate: bool,
    pub exchange: bool,
    pub automatic: bool,
}

impl Default for GwNodeBed {
    fn default() -> Self {
        GwNodeBed { node: -1, cell: -1, Kc: 0.0, dC: 0.0, area: 0.0, locate: false, exchange: true, automatic: false }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i8)]
pub enum GwLinkMode {
    None = 0,
    Auto = 1,
    TwoWay = 2,
    #[default]
    Default = 3,
}

#[derive(Debug, Clone)]
#[allow(non_snake_case)]
pub struct GwLinkRow {
    pub link: i32,
    pub Kc: f64,
    pub dC: f64,
    pub exchange: bool,
}

impl Default for GwLinkRow {
    fn default() -> Self {
        GwLinkRow { link: -1, Kc: 0.0, dC: 0.0, exchange: true }
    }
}

#[derive(Debug, Clone)]
pub struct GwOptions {
    pub independent_clock: bool,
    pub binds_surface_dt: bool,
    pub base_step: f64,
    pub max_step: f64,
    pub budget_fraction: f64,
    pub head_tolerance: f64,
    pub soil_char: SoilChar,
    pub closure: GwClosure,
    pub m_layers: i32,
    pub capillary_diff: bool,
    pub c_gw: f64,
    pub c_col: f64,
    pub force_closed_form: bool,
    pub per_subcatch: bool,
    pub dunne: bool,
    pub gw_et: Bytes,
    pub wilting_suction: f64,
    pub wilting_suction_set: bool,
    pub node_auto: bool,
    pub link_seepage: GwLinkMode,
    pub authored: bool,
}

impl Default for GwOptions {
    fn default() -> Self {
        GwOptions {
            independent_clock: false,
            binds_surface_dt: false,
            base_step: 1.0,
            max_step: 3600.0,
            budget_fraction: 0.5,
            head_tolerance: 0.05,
            soil_char: SoilChar::Russo,
            closure: GwClosure::Auto,
            m_layers: 8,
            capillary_diff: false,
            c_gw: 0.5,
            c_col: 0.9,
            force_closed_form: false,
            per_subcatch: false,
            dunne: true,
            gw_et: b"AUTO".to_vec(),
            wilting_suction: 150.0,
            wilting_suction_set: false,
            node_auto: true,
            link_seepage: GwLinkMode::Default,
            authored: false,
        }
    }
}

#[derive(Debug, Clone)]
#[allow(non_snake_case)]
pub struct GwAquiferRow {
    pub scope: i32,
    pub tag: Bytes,
    pub cell: i32,
    pub Ks: f64,
    pub zs: f64,
    pub theta_s: f64,
    pub theta_r: f64,
    pub alpha: f64,
    pub psi_b: f64,
    pub lambda: f64,
    pub vg_n: f64,
    pub vg_L: f64,
    pub c_loss: f64,
    pub hg0: f64,
    pub soil_char: SoilChar,
    pub closure: GwClosure,
    pub m_layers: i32,
    pub soil_char_set: bool,
    pub closure_set: bool,
}

impl Default for GwAquiferRow {
    fn default() -> Self {
        GwAquiferRow {
            scope: 0,
            tag: Bytes::new(),
            cell: -1,
            Ks: 1.0e-5,
            zs: 5.0,
            theta_s: 0.45,
            theta_r: 0.10,
            alpha: 2.0,
            psi_b: 0.20,
            lambda: 0.40,
            vg_n: 1.6,
            vg_L: 0.5,
            c_loss: 0.0,
            hg0: -1.0,
            soil_char: SoilChar::Russo,
            closure: GwClosure::Auto,
            m_layers: -1,
            soil_char_set: false,
            closure_set: false,
        }
    }
}

#[derive(Debug, Clone)]
pub struct GwBoundaryRow {
    pub cell: i32,
    pub edge: i32,
    pub head: bool,
    pub value: f64,
    pub series: Bytes,
}

impl Default for GwBoundaryRow {
    fn default() -> Self {
        GwBoundaryRow { cell: -1, edge: -1, head: false, value: 0.0, series: Bytes::new() }
    }
}

#[derive(Debug, Clone, Default)]
pub struct SubsurfaceConfig {
    pub options: GwOptions,
    pub rows: Vec<GwAquiferRow>,
    pub node_beds: Vec<GwNodeBed>,
    pub boundaries: Vec<GwBoundaryRow>,
    pub link_rows: Vec<GwLinkRow>,
    pub surface_owners: Vec<SurfaceOwnerRecord>,
}

impl SubsurfaceConfig {
    pub fn empty(&self) -> bool {
        !self.options.authored && self.rows.is_empty() && self.node_beds.is_empty() && self.link_rows.is_empty() && self.boundaries.is_empty()
    }
    pub fn clear(&mut self) {
        *self = SubsurfaceConfig::default();
    }
}

/// upperTrim(): every non-space byte, upper-cased (white space anywhere is dropped).
fn upper_trim(t: &[u8]) -> Bytes {
    t.iter().filter(|&&c| !crate::core::cstdlib::c_isspace(c)).map(|c| c.to_ascii_uppercase()).collect()
}

/// soilCharToken()
pub fn soil_char_token(s: SoilChar) -> &'static [u8] {
    match s {
        SoilChar::Russo => b"RUSSO",
        SoilChar::Gardner => b"GARDNER",
        SoilChar::BrooksCorey => b"BROOKS_COREY",
        SoilChar::VanGenuchten => b"VAN_GENUCHTEN",
    }
}

/// parseSoilChar()
pub fn parse_soil_char(t: &[u8], s: &mut SoilChar) -> bool {
    let u = upper_trim(t);
    let u = u.as_slice();
    if u == b"RUSSO" {
        *s = SoilChar::Russo;
        return true;
    }
    if u == b"GARDNER" {
        *s = SoilChar::Gardner;
        return true;
    }
    if u == b"BROOKS_COREY" || u == b"BROOKSCOREY" || u == b"BROOKS-COREY" || u == b"BC" {
        *s = SoilChar::BrooksCorey;
        return true;
    }
    if u == b"VAN_GENUCHTEN" || u == b"VANGENUCHTEN" || u == b"VAN-GENUCHTEN" || u == b"VG" {
        *s = SoilChar::VanGenuchten;
        return true;
    }
    false
}

/// gwClosureToken()
pub fn gw_closure_token(c: GwClosure) -> &'static [u8] {
    match c {
        GwClosure::Auto => b"AUTO",
        GwClosure::ClosedForm => b"CLOSED_FORM",
        GwClosure::Enslaved => b"ENSLAVED",
        GwClosure::Sigma => b"SIGMA",
    }
}

/// parseGwClosure()
pub fn parse_gw_closure(t: &[u8], c: &mut GwClosure) -> bool {
    let u = upper_trim(t);
    let u = u.as_slice();
    if u == b"AUTO" {
        *c = GwClosure::Auto;
        return true;
    }
    if u == b"CLOSED_FORM" || u == b"CLOSEDFORM" || u == b"CLOSED-FORM" || u == b"A" {
        *c = GwClosure::ClosedForm;
        return true;
    }
    if u == b"ENSLAVED" {
        *c = GwClosure::Enslaved;
        return true;
    }
    if u == b"SIGMA" || u == b"COLUMN" || u == b"B" {
        *c = GwClosure::Sigma;
        return true;
    }
    false
}
