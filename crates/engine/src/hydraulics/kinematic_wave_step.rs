// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/KinematicWave.cpp and
// KinematicWave.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The kinematic-wave and steady-flow routing step (kinematic_wave.rs holds the KW solver's
// state and init).
//
// NOT YET TRANSLATED: KWSolver::execute and the kinwave:: helpers (getLinkInflow,
// updateStorageState, capConduitLoss, finishRouting); they return StepError::NotTranslated.

use crate::core::context::SimulationContext;
use crate::core::engine_step_routing::StepError;
use crate::hydraulics::hyd_structures::StructureSolver;
use crate::hydraulics::kinematic_wave::KWSolver;

impl KWSolver {
    /// KWSolver::execute
    pub fn execute(&mut self, _ctx: &mut SimulationContext, _dt: f64, _structures: &mut StructureSolver) -> Result<i32, StepError> {
        Err(StepError::NotTranslated("KINWAVE routing"))
    }
}

/// kinwave::getLinkInflow
pub fn get_link_inflow(_ctx: &mut SimulationContext, _structures: &mut StructureSolver, _j: i32, _dt: f64) -> Result<f64, StepError> {
    Err(StepError::NotTranslated("STEADY routing"))
}

/// kinwave::updateStorageState
pub fn update_storage_state(_ctx: &mut SimulationContext, _structures: &mut StructureSolver, _order: &[i32], _idx: i32, _n: i32, _dt: f64) -> Result<(), StepError> {
    Err(StepError::NotTranslated("STEADY routing"))
}

/// kinwave::capConduitLoss
pub fn cap_conduit_loss(_ctx: &mut SimulationContext, _ucr: usize, _q: f64) -> f64 {
    0.0
}

/// kinwave::finishRouting
pub fn finish_routing(_ctx: &mut SimulationContext, _structures: &mut StructureSolver, _order: &[i32], _updated: &[i8], _y1: &[f64], _y2: &[f64], _dt: f64) -> Result<(), StepError> {
    Err(StepError::NotTranslated("STEADY routing"))
}
