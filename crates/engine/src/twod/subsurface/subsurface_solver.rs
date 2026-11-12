// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/subsurface/SubsurfaceSolver.hpp and
// SubsurfaceClock.cpp (prepareClock) (HydroCouple/openswmm.engine, branch swmm6_rel, commit
// 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The two-zone subsurface (groundwater) kernel is phase 2b: SurfaceRouter2D::initialize
// refuses an aquifer that would run before SubsurfaceSolver::initialize. What phase 2a reads of
// the router's (inactive) kernel is modelled here: active() (state_.active), the coverage mask,
// the transported-row layout, the ground heat owner, the node beds and the clock queries, each
// at the value a never-initialized kernel has upstream. prepareClock returns at once for a
// kernel without an independent clock, which an inactive kernel never has.

#![allow(clippy::all)]

use crate::core::options::Bytes;
use crate::twod::heat::ground_heat::GroundHeatState;
use crate::twod::subsurface::subsurface_data::GwNodeBed;

/// The part of SubsurfaceTransportState phase 2a reads (row layout and the RK2 age ledger).
#[derive(Debug, Clone, Default)]
pub struct SubsurfaceTransportState {
    pub n_species: i32,
    pub n_layers: i32,
    pub row_names: Vec<Bytes>,
    pub generated_age: Vec<f64>,
}

#[derive(Debug, Clone, Default)]
pub struct SubsurfaceSolver {
    /// state_.active (SubsurfaceState::active).
    pub state_active: bool,
    pub covered_: Vec<u8>,
    pub tr_: SubsurfaceTransportState,
    pub heat_: GroundHeatState,
    pub node_beds_: Vec<GwNodeBed>,
    /// options_.independent_clock.
    pub independent_clock: bool,
}

impl SubsurfaceSolver {
    pub fn active(&self) -> bool {
        self.state_active
    }
    pub fn coverage(&self) -> &Vec<u8> {
        &self.covered_
    }
    pub fn covers(&self, cell: i32) -> bool {
        cell >= 0 && (cell as usize) < self.covered_.len() && self.covered_[cell as usize] != 0
    }
    pub fn transport(&self) -> &SubsurfaceTransportState {
        &self.tr_
    }
    pub fn heat(&self) -> &GroundHeatState {
        &self.heat_
    }
    pub fn heat_mut(&mut self) -> &mut GroundHeatState {
        &mut self.heat_
    }
    pub fn node_beds(&self) -> &Vec<GwNodeBed> {
        &self.node_beds_
    }
    pub fn independent_clock(&self) -> bool {
        self.active() && self.independent_clock
    }
    /// prepareClock(epoch): nothing to prepare without an independent clock.
    pub fn prepare_clock(&mut self) {
        if !self.independent_clock() {
            return;
        }
        // An independent clock needs an active kernel, which phase 2a never has.
    }
}
