// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/gw/GwTransportData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The [GW_*] subsurface-transport authoring rows, owned by SurfaceRouter2D.
//   * std::set<std::string> is BTreeSet<Bytes> (the same byte order).
//   * The two lookup indexes are unordered_maps upstream; they are only searched, never
//     iterated, so HashMap keeps their behavior. The private members are pub here (with
//     upstream's names) so the gate dumper can read them.

#![allow(clippy::all)]

use std::collections::{BTreeSet, HashMap};

use crate::core::options::Bytes;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(i8)]
pub enum GwZone {
    #[default]
    Sat = 0,
    Unsat = 1,
    Layer = 2,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
#[repr(i8)]
pub enum GwScope {
    #[default]
    Global = 0,
    Tag = 1,
    Cell = 2,
}

#[derive(Debug, Clone)]
pub struct GwTransportOptions {
    pub unsat_transport: Bytes,
    pub thermal_model: Bytes,
    pub thermal_layers: i32,
    pub lateral_conduction: bool,
    pub thermal_max_step: f64,
    pub water_conductivity: f64,
    pub air_conductivity: f64,
    pub surface_resistance: f64,
    pub soil_et_heat: Bytes,
    pub thermal_authored: BTreeSet<Bytes>,
    pub transport_pollutants: bool,
    pub transport_msx: bool,
    pub transport_age: bool,
    pub transport_temperature: bool,
    pub dispersion: bool,
    pub conduction: bool,
    pub surface_thermal_bc: Bytes,
    pub surface_thermal_arg: Bytes,
    pub deep_thermal_bc: Bytes,
    pub deep_thermal_arg: Bytes,
    pub deep_depth: f64,
    pub thermal_mixing: Bytes,
    pub c_diff: f64,
    pub authored: bool,
}

impl Default for GwTransportOptions {
    fn default() -> Self {
        GwTransportOptions {
            unsat_transport: b"BULK".to_vec(),
            thermal_model: b"WATER_ONLY".to_vec(),
            thermal_layers: 16,
            lateral_conduction: false,
            thermal_max_step: 60.0,
            water_conductivity: 0.606,
            air_conductivity: 0.025,
            surface_resistance: 0.1,
            soil_et_heat: b"LATENT_SURFACE".to_vec(),
            thermal_authored: BTreeSet::new(),
            transport_pollutants: true,
            transport_msx: true,
            transport_age: true,
            transport_temperature: true,
            dispersion: true,
            conduction: true,
            surface_thermal_bc: b"AIR".to_vec(),
            surface_thermal_arg: Bytes::new(),
            deep_thermal_bc: b"GEOTHERMAL_FLUX".to_vec(),
            deep_thermal_arg: Bytes::new(),
            deep_depth: 0.0,
            thermal_mixing: b"ARITHMETIC".to_vec(),
            c_diff: 0.4,
            authored: false,
        }
    }
}

#[derive(Debug, Clone)]
#[allow(non_snake_case)]
pub struct GwParamsRow {
    pub scope: GwScope,
    pub tag: Bytes,
    pub cell: i32,
    pub rho_s: f64,
    pub c_s: f64,
    pub lambda_s: f64,
    pub a_s: f64,
    pub alpha_L: f64,
    pub alpha_T: f64,
    pub D_m: f64,
    pub D_v: f64,
    pub geo_flux: f64,
}

impl Default for GwParamsRow {
    fn default() -> Self {
        GwParamsRow {
            scope: GwScope::Global,
            tag: Bytes::new(),
            cell: -1,
            rho_s: 2650.0,
            c_s: 880.0,
            lambda_s: 2.0,
            a_s: 0.0,
            alpha_L: 1.0,
            alpha_T: 0.1,
            D_m: 1.0e-9,
            D_v: 1.0e-9,
            geo_flux: 0.065,
        }
    }
}

#[derive(Debug, Clone)]
pub struct GwSorptionRow {
    pub scope: GwScope,
    pub tag: Bytes,
    pub cell: i32,
    pub species: Bytes,
    pub kd: f64,
    pub decay: f64,
}

impl Default for GwSorptionRow {
    fn default() -> Self {
        GwSorptionRow { scope: GwScope::Global, tag: Bytes::new(), cell: -1, species: Bytes::new(), kd: 0.0, decay: -1.0 }
    }
}

#[derive(Debug, Clone)]
pub struct GwInitialQualityRow {
    pub scope: GwScope,
    pub tag: Bytes,
    pub cell: i32,
    pub zone: GwZone,
    pub layer: i32,
    pub species: Bytes,
    pub value: f64,
}

impl Default for GwInitialQualityRow {
    fn default() -> Self {
        GwInitialQualityRow {
            scope: GwScope::Global,
            tag: Bytes::new(),
            cell: -1,
            zone: GwZone::Sat,
            layer: -1,
            species: Bytes::new(),
            value: 0.0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct GwBoundaryQualityRow {
    pub cell: i32,
    pub edge: i32,
    pub species: Bytes,
    pub kind: Bytes,
    pub value: f64,
    pub ts_name: Bytes,
}

impl Default for GwBoundaryQualityRow {
    fn default() -> Self {
        GwBoundaryQualityRow { cell: -1, edge: -1, species: Bytes::new(), kind: b"CONC".to_vec(), value: 0.0, ts_name: Bytes::new() }
    }
}

#[derive(Debug, Clone)]
pub struct GwSourceSpeciesTerm {
    pub species: Bytes,
    pub kind: Bytes,
    pub value: f64,
    pub ts_name: Bytes,
}

impl Default for GwSourceSpeciesTerm {
    fn default() -> Self {
        GwSourceSpeciesTerm { species: Bytes::new(), kind: b"CONC".to_vec(), value: 0.0, ts_name: Bytes::new() }
    }
}

#[derive(Debug, Clone)]
pub struct GwSourceRow {
    pub name: Bytes,
    pub scope: GwScope,
    pub tag: Bytes,
    pub cell: i32,
    pub x: f64,
    pub y: f64,
    pub by_xy: bool,
    pub flow: f64,
    pub flow_ts: Bytes,
    pub scale: f64,
    pub species: Vec<GwSourceSpeciesTerm>,
}

impl Default for GwSourceRow {
    fn default() -> Self {
        GwSourceRow {
            name: Bytes::new(),
            scope: GwScope::Cell,
            tag: Bytes::new(),
            cell: -1,
            x: 0.0,
            y: 0.0,
            by_xy: false,
            flow: 0.0,
            flow_ts: Bytes::new(),
            scale: 1.0,
            species: Vec::new(),
        }
    }
}

/// GwTransportData::QualityKey: the identity of an initial-quality row.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct QualityKey {
    pub scope: GwScope,
    pub tag: Bytes,
    pub cell: i32,
    pub zone: GwZone,
    pub layer: i32,
    pub species: Bytes,
}

impl QualityKey {
    pub fn new(r: &GwInitialQualityRow) -> Self {
        QualityKey { scope: r.scope, tag: r.tag.clone(), cell: r.cell, zone: r.zone, layer: r.layer, species: r.species.clone() }
    }
}

#[derive(Debug, Clone)]
pub struct GwTransportData {
    pub options: GwTransportOptions,
    pub params: Vec<GwParamsRow>,
    pub sorption: Vec<GwSorptionRow>,
    pub initial_quality: Vec<GwInitialQualityRow>,
    pub initial_quality_file: Bytes,
    pub boundary_quality: Vec<GwBoundaryQualityRow>,
    pub sources: Vec<GwSourceRow>,
    pub source_index_: HashMap<Bytes, usize>,
    pub quality_index_: HashMap<QualityKey, usize>,
    pub source_index_size_: usize,
    pub quality_index_size_: usize,
}

impl Default for GwTransportData {
    fn default() -> Self {
        GwTransportData {
            options: GwTransportOptions::default(),
            params: Vec::new(),
            sorption: Vec::new(),
            initial_quality: Vec::new(),
            initial_quality_file: Bytes::new(),
            boundary_quality: Vec::new(),
            sources: Vec::new(),
            source_index_: HashMap::new(),
            quality_index_: HashMap::new(),
            source_index_size_: GwTransportData::NO_ROW,
            quality_index_size_: GwTransportData::NO_ROW,
        }
    }
}

impl GwTransportData {
    pub const NO_ROW: usize = usize::MAX;

    pub fn empty(&self) -> bool {
        !self.options.authored
            && self.params.is_empty()
            && self.sorption.is_empty()
            && self.initial_quality.is_empty()
            && self.initial_quality_file.is_empty()
            && self.boundary_quality.is_empty()
            && self.sources.is_empty()
    }

    pub fn clear(&mut self) {
        *self = GwTransportData::default();
    }

    pub fn invalidate_authoring_indexes(&mut self) {
        self.source_index_size_ = Self::NO_ROW;
        self.quality_index_size_ = Self::NO_ROW;
    }

    pub fn find_source(&mut self, name: &[u8]) -> usize {
        if self.source_index_size_ != self.sources.len() {
            self.source_index_size_ = Self::NO_ROW;
            self.source_index_.clear();
            self.source_index_.reserve(self.sources.len());
            for i in 0..self.sources.len() {
                self.source_index_.entry(self.sources[i].name.clone()).or_insert(i); // emplace: the first wins
            }
            self.source_index_size_ = self.sources.len();
        }
        match self.source_index_.get(name) {
            Some(&i) => i,
            None => Self::NO_ROW,
        }
    }

    pub fn find_initial_quality(&mut self, row: &GwInitialQualityRow) -> usize {
        if self.quality_index_size_ != self.initial_quality.len() {
            self.quality_index_size_ = Self::NO_ROW;
            self.quality_index_.clear();
            self.quality_index_.reserve(self.initial_quality.len());
            for i in 0..self.initial_quality.len() {
                self.quality_index_.entry(QualityKey::new(&self.initial_quality[i])).or_insert(i);
            }
            self.quality_index_size_ = self.initial_quality.len();
        }
        match self.quality_index_.get(&QualityKey::new(row)) {
            Some(&i) => i,
            None => Self::NO_ROW,
        }
    }

    pub fn source_appended(&mut self) {
        if self.source_index_size_ != Self::NO_ROW && self.source_index_size_ + 1 == self.sources.len() {
            let n = self.sources.len();
            self.source_index_.entry(self.sources[n - 1].name.clone()).or_insert(n - 1);
            self.source_index_size_ = n;
        }
    }

    pub fn initial_quality_appended(&mut self) {
        if self.quality_index_size_ != Self::NO_ROW && self.quality_index_size_ + 1 == self.initial_quality.len() {
            let n = self.initial_quality.len();
            self.quality_index_.entry(QualityKey::new(&self.initial_quality[n - 1])).or_insert(n - 1);
            self.quality_index_size_ = n;
        }
    }

    pub fn remove_source(&mut self, index: usize) {
        if self.source_index_size_ == self.sources.len() && index + 1 == self.sources.len() {
            let key = self.sources[index].name.clone();
            if self.source_index_.get(&key) == Some(&index) {
                self.source_index_.remove(&key);
            }
            self.source_index_size_ -= 1;
        } else {
            self.source_index_size_ = Self::NO_ROW;
        }
        self.sources.remove(index);
    }

    pub fn remove_initial_quality(&mut self, index: usize) {
        if self.quality_index_size_ == self.initial_quality.len() && index + 1 == self.initial_quality.len() {
            let key = QualityKey::new(&self.initial_quality[index]);
            if self.quality_index_.get(&key) == Some(&index) {
                self.quality_index_.remove(&key);
            }
            self.quality_index_size_ -= 1;
        } else {
            self.quality_index_size_ = Self::NO_ROW;
        }
        self.initial_quality.remove(index);
    }
}

/// gwZoneToken()
pub fn gw_zone_token(z: GwZone) -> &'static [u8] {
    match z {
        GwZone::Sat => b"SAT",
        GwZone::Unsat => b"UNSAT",
        GwZone::Layer => b"LAYER",
    }
}

/// gwScopeToken()
pub fn gw_scope_token(s: GwScope) -> &'static [u8] {
    match s {
        GwScope::Global => b"*",
        GwScope::Tag => b"TAG",
        GwScope::Cell => b"CELL",
    }
}
