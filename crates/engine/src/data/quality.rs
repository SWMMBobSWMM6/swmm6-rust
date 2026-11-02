// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/QualityData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson (first pass by tools/soa2rs.py, finished by hand). Not an
// OpenSWMM release.
//
// The land-use, buildup, washoff and treatment stores, one Vec per attribute, sized and
// cleared exactly as upstream's resize() does. The size arithmetic keeps upstream's
// casts: Buildup/Washoff convert each count to size_t before multiplying (unsigned, so
// wrapping_mul), while TreatmentData multiplies the two ints first and converts after.
// `compiled.resize(total)` keeps the existing prefix, as std::vector::resize does.
// `hasAny` is `has_any`.

use crate::core::options::Bytes;
use crate::data::treatment_types::TreatExpr;

// ============================================================================
// Land use definitions
// ============================================================================

#[derive(Debug, Clone, Default)]
pub struct LanduseData {
    /// Days between sweeps.
    pub sweep_interval: Vec<f64>,
    /// Max removal fraction (0-100).
    pub sweep_removal: Vec<f64>,
    /// Days since last swept.
    pub last_swept: Vec<f64>,
    /// Object comment from the INP file, joined by literal "\\n"; empty means none.
    pub comments: Vec<Bytes>,
}

impl LanduseData {
    pub fn count(&self) -> i32 {
        self.sweep_interval.len() as i32
    }

    pub fn resize(&mut self, n: i32) {
        let un = n as usize;
        self.sweep_interval.clear();
        self.sweep_interval.resize(un, 0.0);
        self.sweep_removal.clear();
        self.sweep_removal.resize(un, 0.0);
        self.last_swept.clear();
        self.last_swept.resize(un, 0.0);
        self.comments.clear();
        self.comments.resize(un, Vec::new());
    }

    pub fn shrink_to_fit(&mut self) {
        self.sweep_interval.shrink_to_fit();
        self.sweep_removal.shrink_to_fit();
        self.last_swept.shrink_to_fit();
        self.comments.shrink_to_fit();
    }
}

// ============================================================================
// Buildup function per (landuse x pollutant)
// ============================================================================

/// Index: [landuse * n_pollutants + pollutant].
#[derive(Debug, Clone, Default)]
pub struct BuildupData {
    /// 0=NONE,1=POW,2=EXP,3=SAT,4=EXT
    pub func_type: Vec<i32>,
    pub coeff1: Vec<f64>,
    pub coeff2: Vec<f64>,
    pub coeff3: Vec<f64>,
    /// 0=PER_AREA, 1=PER_CURB
    pub normalizer: Vec<i32>,
    pub n_landuses: i32,
    pub n_pollutants: i32,
}

// The member initializers (n_landuses = n_pollutants = 0) are Rust's defaults, so Default
// is derived, here and for WashoffData and TreatmentData.

impl BuildupData {
    pub fn resize(&mut self, nlu: i32, npoll: i32) {
        self.n_landuses = nlu;
        self.n_pollutants = npoll;
        let total = (nlu as usize).wrapping_mul(npoll as usize);
        self.func_type.clear();
        self.func_type.resize(total, 0);
        self.coeff1.clear();
        self.coeff1.resize(total, 0.0);
        self.coeff2.clear();
        self.coeff2.resize(total, 0.0);
        self.coeff3.clear();
        self.coeff3.resize(total, 0.0);
        self.normalizer.clear();
        self.normalizer.resize(total, 0);
    }

    pub fn shrink_to_fit(&mut self) {
        self.func_type.shrink_to_fit();
        self.coeff1.shrink_to_fit();
        self.coeff2.shrink_to_fit();
        self.coeff3.shrink_to_fit();
        self.normalizer.shrink_to_fit();
    }
}

// ============================================================================
// Washoff function per (landuse x pollutant)
// ============================================================================

/// Index: [landuse * n_pollutants + pollutant].
#[derive(Debug, Clone, Default)]
pub struct WashoffData {
    /// 0=NONE,1=EXP,2=RC,3=EMC
    pub func_type: Vec<i32>,
    pub coeff: Vec<f64>,
    pub expon: Vec<f64>,
    /// 0-100
    pub sweep_effic: Vec<f64>,
    /// 0-100
    pub bmp_effic: Vec<f64>,
    pub n_landuses: i32,
    pub n_pollutants: i32,
}

impl WashoffData {
    pub fn resize(&mut self, nlu: i32, npoll: i32) {
        self.n_landuses = nlu;
        self.n_pollutants = npoll;
        let total = (nlu as usize).wrapping_mul(npoll as usize);
        self.func_type.clear();
        self.func_type.resize(total, 0);
        self.coeff.clear();
        self.coeff.resize(total, 0.0);
        self.expon.clear();
        self.expon.resize(total, 0.0);
        self.sweep_effic.clear();
        self.sweep_effic.resize(total, 0.0);
        self.bmp_effic.clear();
        self.bmp_effic.resize(total, 0.0);
    }

    pub fn shrink_to_fit(&mut self) {
        self.func_type.shrink_to_fit();
        self.coeff.shrink_to_fit();
        self.expon.shrink_to_fit();
        self.sweep_effic.shrink_to_fit();
        self.bmp_effic.shrink_to_fit();
    }
}

// ============================================================================
// Treatment expression per (node x pollutant)
// ============================================================================

#[derive(Debug, Clone, Default)]
pub struct TreatmentData {
    /// Index: [node * n_pollutants + pollutant], e.g. "R = 0.5 * exp(-0.1 * DT)".
    pub expressions: Vec<Bytes>,
    pub n_nodes: i32,
    pub n_pollutants: i32,
    /// Compiled treatment expressions (same indexing as expressions).
    pub compiled: Vec<TreatExpr>,
    /// Per-node flag: true if any pollutant has a treatment expression.
    pub has_treatment: Vec<bool>,
    /// Per-node inflow concentrations (size = n_pollutants, reused each timestep).
    pub cin: Vec<f64>,
    /// Per-node removal fractions (size = n_pollutants, reused each timestep):
    /// -1.0 = not computed, 0-1 = computed, > 1.0 = being evaluated (cycle detect).
    pub removal: Vec<f64>,
}

impl TreatmentData {
    pub fn resize(&mut self, nn: i32, npoll: i32) {
        self.n_nodes = nn;
        self.n_pollutants = npoll;
        // static_cast<std::size_t>(nn * npoll): the product is taken in int.
        let total = (nn * npoll) as usize;
        self.expressions.clear();
        self.expressions.resize(total, Vec::new());
        self.compiled.resize(total, TreatExpr::default());
        self.has_treatment.clear();
        self.has_treatment.resize(nn as usize, false);
        self.cin.clear();
        self.cin.resize(npoll as usize, 0.0);
        self.removal.clear();
        self.removal.resize(npoll as usize, -1.0);
    }

    pub fn shrink_to_fit(&mut self) {
        self.expressions.shrink_to_fit();
        self.compiled.shrink_to_fit();
    }

    /// True if any authored expression is nonempty (whitespace alone counts).
    pub fn has_any(&self) -> bool {
        for e in &self.expressions {
            if !e.is_empty() {
                return true;
            }
        }
        false
    }
}
