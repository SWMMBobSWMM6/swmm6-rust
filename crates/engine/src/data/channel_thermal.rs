// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/ChannelThermalData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Solver-neutral finite channel water regions and their thermal / species increments.
// Region IDs belong to the actual transport engine (ARD regions are cells, LEGACY/LARD
// regions are links). Defaults are upstream's member initializers; row names are bytes.

use crate::core::options::Bytes;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChannelThermalRegion {
    pub link: i32,
    pub region: i32,
    pub begin: f64,
    pub end: f64,
    /// SI m3
    pub volume: f64,
    /// degC
    pub temperature: f64,
    /// SI m
    pub length: f64,
    /// SI m2
    pub surface_area: f64,
    /// SI m2
    pub bed_area: f64,
    /// SI m
    pub depth: f64,
}

impl Default for ChannelThermalRegion {
    fn default() -> Self {
        ChannelThermalRegion {
            link: -1,
            region: -1,
            begin: 0.0,
            end: 1.0,
            volume: 0.0,
            temperature: 0.0,
            length: 0.0,
            surface_area: 0.0,
            bed_area: 0.0,
            depth: 0.0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChannelThermalIncrement {
    pub link: i32,
    pub region: i32,
    pub celsius: f64,
}

impl Default for ChannelThermalIncrement {
    fn default() -> Self {
        ChannelThermalIncrement { link: -1, region: -1, celsius: 0.0 }
    }
}

/// Finite channel region and signed concentration update for reciprocal 2D porewater
/// contacts. Rows use the channel transport owner's pollutant, age, temperature ordering.
#[derive(Debug, Clone, PartialEq)]
pub struct ChannelSpeciesRegion {
    pub link: i32,
    pub region: i32,
    pub begin: f64,
    pub end: f64,
    pub volume: f64,
    pub row_names: Vec<Bytes>,
    pub concentration: Vec<f64>,
}

impl Default for ChannelSpeciesRegion {
    fn default() -> Self {
        ChannelSpeciesRegion {
            link: -1,
            region: -1,
            begin: 0.0,
            end: 1.0,
            volume: 0.0,
            row_names: Vec::new(),
            concentration: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ChannelSpeciesIncrement {
    pub link: i32,
    pub region: i32,
    pub concentration_delta: Vec<f64>,
}

impl Default for ChannelSpeciesIncrement {
    fn default() -> Self {
        ChannelSpeciesIncrement { link: -1, region: -1, concentration_delta: Vec::new() }
    }
}
