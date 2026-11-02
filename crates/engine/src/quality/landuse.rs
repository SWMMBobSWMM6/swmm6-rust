// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/quality/Landuse.hpp and Landuse.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Land-use buildup/washoff (namespace landuse), init part: the function-type enums, the
// per (land use x pollutant) parameter structs, the per-subcatchment surface quality store
// (SurfaceQualitySoA::resize, bu_idx) and LanduseSolver::init.
//
// Differences of form:
//   * Upstream fills BuildupParams::type / WashoffParams::type with
//     `static_cast<BuildupType>(int)`, which keeps any int. `BuildupType::from_i32` and
//     `WashoffType::from_i32` accept the enumerators only and panic on anything else; the
//     [BUILDUP] / [WASHOFF] handlers store only 0..4 and 0..3 (an unmatched keyword is 0).
//   * `double coeff[3] = {}` is `[f64; 3]`.
//
// Run-time functions, translated with the quality phase (not here): getBuildupMass,
// getBuildupDays (file-local), LanduseSolver::{computeBuildup, computeWashoff,
// applyCoPollutant}.

/// Buildup function types (`enum class BuildupType : int`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(i32)]
pub enum BuildupType {
    #[default]
    None = 0,
    Power = 1,
    Expon = 2,
    Satur = 3,
    External = 4,
}

impl BuildupType {
    /// `static_cast<BuildupType>(v)` for the values the input handlers store.
    pub fn from_i32(v: i32) -> BuildupType {
        match v {
            0 => BuildupType::None,
            1 => BuildupType::Power,
            2 => BuildupType::Expon,
            3 => BuildupType::Satur,
            4 => BuildupType::External,
            _ => panic!("BuildupType {} has no enumerator (the [BUILDUP] handler stores 0..4)", v),
        }
    }
}

/// Washoff function types (`enum class WashoffType : int`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(i32)]
pub enum WashoffType {
    #[default]
    None = 0,
    Expon = 1,
    Rating = 2,
    Emc = 3,
}

impl WashoffType {
    /// `static_cast<WashoffType>(v)` for the values the input handlers store.
    pub fn from_i32(v: i32) -> WashoffType {
        match v {
            0 => WashoffType::None,
            1 => WashoffType::Expon,
            2 => WashoffType::Rating,
            3 => WashoffType::Emc,
            _ => panic!("WashoffType {} has no enumerator (the [WASHOFF] handler stores 0..3)", v),
        }
    }
}

// ============================================================================
// Per (landuse x pollutant) buildup/washoff parameters
// ============================================================================

/// All member initializers are Rust's defaults (NONE is 0), so Default is derived.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct BuildupParams {
    pub type_: BuildupType,
    /// c0 (max), c1 (rate), c2 (exponent)
    pub coeff: [f64; 3],
    /// Time to reach max buildup
    pub max_days: f64,
    /// 0=PER_AREA, 1=PER_CURB
    pub normalizer: i32,
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct WashoffParams {
    pub type_: WashoffType,
    /// Washoff coefficient
    pub coeff: f64,
    /// Washoff exponent
    pub expon: f64,
    /// Street sweep removal [0-1]
    pub sweep_effic: f64,
    /// BMP removal [0-1]
    pub bmp_effic: f64,
}

// ============================================================================
// Surface quality SoA (per subcatchment x pollutant)
// ============================================================================

#[derive(Debug, Clone, Default, PartialEq)]
pub struct SurfaceQualitySoA {
    pub n_subcatch: i32,
    pub n_landuses: i32,
    pub n_pollutants: i32,
    /// Per-landuse buildup: [subcatch * n_landuses * n_pollutants + lu * n_pollutants + p]
    pub buildup: Vec<f64>,
    /// Washoff concentration [subcatch * n_pollutants + pollutant] (mass/vol)
    pub washoff_conc: Vec<f64>,
}

impl SurfaceQualitySoA {
    /// size_t(sc) * size_t(n_landuses * n_pollutants) + size_t(lu) * size_t(n_pollutants)
    /// + size_t(p); the inner product is taken in int, as upstream does.
    pub fn bu_idx(&self, sc: i32, lu: i32, p: i32) -> usize {
        (sc as usize)
            .wrapping_mul((self.n_landuses * self.n_pollutants) as usize)
            .wrapping_add((lu as usize).wrapping_mul(self.n_pollutants as usize))
            .wrapping_add(p as usize)
    }

    /// The products are taken in int, then converted (`static_cast<size_t>(n_sc*n_lu*n_poll)`).
    pub fn resize(&mut self, n_sc: i32, n_lu: i32, n_poll: i32) {
        self.n_subcatch = n_sc;
        self.n_landuses = n_lu;
        self.n_pollutants = n_poll;
        self.buildup.clear();
        self.buildup.resize((n_sc * n_lu * n_poll) as usize, 0.0);
        self.washoff_conc.clear();
        self.washoff_conc.resize((n_sc * n_poll) as usize, 0.0);
    }
}

// ============================================================================
// Landuse quality solver
// ============================================================================

#[derive(Debug, Clone, Default, PartialEq)]
pub struct LanduseSolver {
    /// Per (landuse x pollutant) parameters. Index: [lu * n_pollutants + p]
    pub buildup_params: Vec<BuildupParams>,
    pub washoff_params: Vec<WashoffParams>,
    pub n_landuses_: i32,
    pub n_pollutants_: i32,
}

impl LanduseSolver {
    /// Sets the counts and resizes the parameter tables; `resize` keeps the overlapping old
    /// parameters and default-initializes added entries (not a full reset).
    pub fn init(&mut self, n_landuses: i32, n_pollutants: i32) {
        self.n_landuses_ = n_landuses;
        self.n_pollutants_ = n_pollutants;
        let total = (n_landuses * n_pollutants) as usize;
        self.buildup_params.resize(total, BuildupParams::default());
        self.washoff_params.resize(total, WashoffParams::default());
    }
}
