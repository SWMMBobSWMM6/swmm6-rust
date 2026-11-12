// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/data/SurfaceTransportState.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The species state carried on the 2D surface (transport rows, per-cell mass and the ledgers).
// Phase 2a sizes it (SurfaceRouter2D::initialize resizes it with zero rows when the model
// carries none) and refuses a model whose surface carries transported rows (phase 2b), so the
// rows here stay empty in every run the translation completes. The std::map member is a
// BTreeMap (the same key order).

#![allow(clippy::all)]

use std::collections::BTreeMap;

use crate::core::options::Bytes;

/// SurfaceTransportState::BoundaryQualityRow.
#[derive(Debug, Clone)]
pub struct BoundaryQualityRow {
    pub slot: i32,
    pub species: i32,
    pub conc: f64,
}

impl Default for BoundaryQualityRow {
    fn default() -> Self {
        BoundaryQualityRow { slot: -1, species: -1, conc: 0.0 }
    }
}

#[derive(Debug, Clone)]
pub struct SurfaceTransportState {
    pub n_species: i32,
    pub n_cells: i32,
    pub n_pollut: i32,
    pub n_msx: i32,
    pub age_row: i32,
    pub temp_row: i32,
    pub row_names: Vec<Bytes>,
    pub cell_mass: Vec<f64>,
    pub lost_infiltration: Vec<f64>,
    pub lost_boundary: Vec<f64>,
    pub gained_external: Vec<f64>,
    pub lost_external: Vec<f64>,
    pub gained_bed: Vec<f64>,
    pub lost_coupling: Vec<f64>,
    pub exch_mass: Vec<f64>,
    pub exch_spill: Vec<f64>,
    pub rain_conc: Vec<f64>,
    pub bc_conc: Vec<f64>,
    pub runtime_bc_conc: BTreeMap<i32, Vec<f64>>,
    pub bc_quality_rows: Vec<BoundaryQualityRow>,
    pub dispersion_limiter_binds: i64,
    pub gained_rainfall: Vec<f64>,
    pub gained_boundary: Vec<f64>,
    pub gained_coupling: Vec<f64>,
    pub gained_exfiltration: Vec<f64>,
    pub generated_age: Vec<f64>,
    pub gained_washoff: Vec<f64>,
    pub cell_runoff_vol: Vec<f64>,
    pub coupling_src: Vec<f64>,
}

impl Default for SurfaceTransportState {
    fn default() -> Self {
        SurfaceTransportState {
            n_species: 0,
            n_cells: 0,
            n_pollut: 0,
            n_msx: 0,
            age_row: -1,
            temp_row: -1,
            row_names: Vec::new(),
            cell_mass: Vec::new(),
            lost_infiltration: Vec::new(),
            lost_boundary: Vec::new(),
            gained_external: Vec::new(),
            lost_external: Vec::new(),
            gained_bed: Vec::new(),
            lost_coupling: Vec::new(),
            exch_mass: Vec::new(),
            exch_spill: Vec::new(),
            rain_conc: Vec::new(),
            bc_conc: Vec::new(),
            runtime_bc_conc: BTreeMap::new(),
            bc_quality_rows: Vec::new(),
            dispersion_limiter_binds: 0,
            gained_rainfall: Vec::new(),
            gained_boundary: Vec::new(),
            gained_coupling: Vec::new(),
            gained_exfiltration: Vec::new(),
            generated_age: Vec::new(),
            gained_washoff: Vec::new(),
            cell_runoff_vol: Vec::new(),
            coupling_src: Vec::new(),
        }
    }
}

fn assign(v: &mut Vec<f64>, n: usize, x: f64) {
    v.clear();
    v.resize(n, x);
}

impl SurfaceTransportState {
    pub fn signed_row(&self, s: i32) -> bool {
        s == self.temp_row
    }

    pub fn row_index(&self, name: &[u8]) -> i32 {
        for (i, r) in self.row_names.iter().enumerate() {
            if r.as_slice() == name {
                return i as i32;
            }
        }
        -1
    }

    pub fn active(&self) -> bool {
        self.n_species > 0 && self.n_cells > 0
    }

    pub fn resize(&mut self, n_spec: i32, n_tri: i32, n_coupling_points: i32) {
        self.n_species = if n_spec > 0 { n_spec } else { 0 };
        self.n_cells = if n_tri > 0 { n_tri } else { 0 };
        self.n_pollut = self.n_species;
        self.n_msx = 0;
        self.age_row = -1;
        self.temp_row = -1; // S4: router refines
        self.row_names.clear();
        let ns = self.n_species as usize;
        let ncp = if n_coupling_points > 0 { n_coupling_points as usize } else { 0 };
        assign(&mut self.cell_mass, ns * self.n_cells as usize, 0.0);
        assign(&mut self.lost_infiltration, ns, 0.0);
        assign(&mut self.lost_boundary, ns, 0.0);
        assign(&mut self.gained_external, ns, 0.0);
        assign(&mut self.lost_external, ns, 0.0);
        self.gained_bed.clear();
        assign(&mut self.lost_coupling, ns, 0.0);
        assign(&mut self.exch_mass, ns * ncp, 0.0);
        assign(&mut self.exch_spill, ncp, 0.0);
        self.rain_conc.clear(); // S2: router fills when pollutants carry one
        self.runtime_bc_conc.clear();
        self.bc_conc.clear(); // S2: solver sizes at initialize
        self.bc_quality_rows.clear();
        self.dispersion_limiter_binds = 0;
        assign(&mut self.gained_rainfall, ns, 0.0);
        assign(&mut self.gained_boundary, ns, 0.0);
        assign(&mut self.gained_coupling, ns, 0.0);
        assign(&mut self.gained_exfiltration, ns, 0.0); // T7.1
        self.generated_age.clear();
        self.coupling_src.clear(); // S3: router sizes when outfalls exist
    }

    pub fn clear(&mut self) {
        *self = SurfaceTransportState::default();
    }

    pub fn idx(&self, s: i32, c: i32) -> usize {
        s as usize * self.n_cells as usize + c as usize
    }
}
