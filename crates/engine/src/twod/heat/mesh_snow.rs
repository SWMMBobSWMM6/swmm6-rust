// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/surface/MeshSnow.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// surface::MeshSnow as far as phase 2a reaches it: snow on the mesh is phase 2b
// (SurfaceRouter2D::initialize refuses it before MeshSnow::initialize), so the object is only
// ever default-constructed (SurfaceStateData::resize assigns `snow = {}`). The three members the
// 2a queries read (active_cells_ for active(), pack_ for owns(), liquid_rate_ for
// liquidRates()) keep upstream's names; the others are not modelled.

#![allow(clippy::all)]

#[derive(Debug, Clone, Default)]
pub struct MeshSnow {
    pub active_cells_: Vec<i32>,
    pub pack_: Vec<i32>,
    pub liquid_rate_: Vec<f64>,
}

impl MeshSnow {
    pub fn active(&self) -> bool {
        !self.active_cells_.is_empty()
    }
    pub fn liquid_rates(&self) -> &Vec<f64> {
        &self.liquid_rate_
    }
    pub fn owns(&self, cell: i32) -> bool {
        cell >= 0 && (cell as usize) < self.pack_.len() && self.pack_[cell as usize] >= 0
    }
}
