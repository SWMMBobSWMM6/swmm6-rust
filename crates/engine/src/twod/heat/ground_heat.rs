// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/heat/GroundHeat.hpp (HydroCouple/openswmm.engine,
// branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the Apache
// License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// GroundHeatState as far as phase 2a reaches it: the porous ground-column heat owner is phase
// 2b, so its thermal cells are never allocated here and active() is false. The members the 2a
// run path touches (cells, lid_owners) and the scalar configuration keep upstream's names and
// defaults; the borrowed pointers (context, snow_owner, config) are flags (set or not).
// GroundThermalCell and GroundLidThermalOwner carry no fields here: no 2a path creates one.

#![allow(clippy::all)]

/// A ground thermal cell (phase 2b; never created in 2a).
#[derive(Debug, Clone, Default)]
pub struct GroundThermalCell;

/// A Mesh Richards LID's thermal footprint on a soil cell (phase 2b; never created in 2a).
#[derive(Debug, Clone, Default)]
pub struct GroundLidThermalOwner;

#[derive(Debug, Clone)]
pub struct GroundHeatState {
    pub context: bool,
    pub snow_owner: bool,
    pub config: bool,
    pub cells: Vec<GroundThermalCell>,
    pub lid_owners: Vec<Vec<GroundLidThermalOwner>>,
    pub top_now: f64,
    pub deep_now: f64,
    pub density: f64,
    pub water_capacity: f64,
    pub top_series: i32,
    pub deep_series: i32,
}

impl Default for GroundHeatState {
    fn default() -> Self {
        GroundHeatState {
            context: false,
            snow_owner: false,
            config: false,
            cells: Vec::new(),
            lid_owners: Vec::new(),
            top_now: 0.0,
            deep_now: 0.0,
            density: 1000.0,
            water_capacity: 4184000.0,
            top_series: -1,
            deep_series: -1,
        }
    }
}

impl GroundHeatState {
    pub fn active(&self) -> bool {
        !self.cells.is_empty()
    }
    pub fn clear(&mut self) {
        *self = GroundHeatState::default();
    }
}
