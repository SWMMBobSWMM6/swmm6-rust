// SPDX-License-Identifier: Apache-2.0
//
// STUB, written 5 November 2026 by Claude (Anthropic) at the request of Robert Dickinson, so
// that SWMMEngine::step (core/engine_step.rs) compiles and the step gate can report where it
// stops. The translation of SWMMEngine::stepRunoff, stepGroundwater, stepSurfaceQuality and
// the runoff helpers (OpenSWMM Engine src/engine/core/SWMMEngine.cpp,
// HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6, Copyright 2026 Caleb Buahin,
// Apache-2.0) REPLACES THIS WHOLE FILE.
//
// The signature step() calls (see PLAN.md):
//   SwmmEngine::step_runoff(&mut self, dt_routing: f64) -> Result<(), Terminated>

use crate::core::engine::SwmmEngine;
use crate::core::engine_init_hydraulics::Terminated;
use crate::core::engine_step::stub;

impl SwmmEngine {
    /// STUB of SWMMEngine::stepRunoff.
    pub(crate) fn step_runoff(&mut self, _dt_routing: f64) -> Result<(), Terminated> {
        Err(stub("SWMMEngine::stepRunoff"))
    }
    /// STUB of SWMMEngine::stepGroundwater (stepRunoff calls it).
    #[allow(dead_code)]
    pub(crate) fn step_groundwater(&mut self, _dt_runoff: f64) -> Result<(), Terminated> {
        Err(stub("SWMMEngine::stepGroundwater"))
    }
    /// STUB of SWMMEngine::stepSurfaceQuality (stepRunoff calls it).
    #[allow(dead_code)]
    pub(crate) fn step_surface_quality(&mut self, _dt_runoff: f64) -> Result<(), Terminated> {
        Err(stub("SWMMEngine::stepSurfaceQuality"))
    }
}
