// SPDX-License-Identifier: Apache-2.0
//
// STUB, written 5 November 2026 by Claude (Anthropic) at the request of Robert Dickinson, so
// that SWMMEngine::step (core/engine_step.rs) compiles and the step gate can report where it
// stops. The routing agent's translation of SWMMEngine::stepRouting (OpenSWMM Engine
// src/engine/core/SWMMEngine.cpp, HydroCouple/openswmm.engine, branch swmm6_rel, commit
// 88329d6, Copyright 2026 Caleb Buahin, Apache-2.0) REPLACES THIS WHOLE FILE.
//
// The signatures step() calls (see PLAN.md):
//   SwmmEngine::step_routing(&mut self, dt_routing: f64) -> Result<(), Terminated>
//   Router::get_adaptive_step(&mut self, ctx: &mut SimulationContext, fixed_step: f64, courant: f64) -> f64
//   Router::accumulate_link_step_stats(&mut self, ctx: &mut SimulationContext, dt: f64, in_window: bool)
//   InletSolver::compute_all(&mut self, ctx: &mut SimulationContext, dt: f64)   (from
//     SwmmEngine::assemble_lateral_inflows, core/engine_step.rs)
// The Router and InletSolver methods belong in hydraulics/routing.rs with the rest of Routing.cpp; they
// (and hydraulics/inlet.rs) are stubbed here only so that this file is the one thing the
// merge replaces. They return no Result upstream, so their stubs panic with "STUB: ..." (the
// step harness catches it).

use crate::core::context::SimulationContext;
use crate::core::engine::SwmmEngine;
use crate::core::engine_init_hydraulics::Terminated;
use crate::core::engine_step::stub;
use crate::hydraulics::inlet::InletSolver;
use crate::hydraulics::routing::Router;

impl SwmmEngine {
    /// STUB of SWMMEngine::stepRouting.
    pub(crate) fn step_routing(&mut self, _dt_routing: f64) -> Result<(), Terminated> {
        Err(stub("SWMMEngine::stepRouting"))
    }
}

impl Router {
    /// STUB of Router::getAdaptiveStep.
    pub fn get_adaptive_step(&mut self, _ctx: &mut SimulationContext, _fixed_step: f64, _courant: f64) -> f64 {
        panic!("STUB: Router::getAdaptiveStep")
    }
    /// STUB of Router::accumulateLinkStepStats.
    pub fn accumulate_link_step_stats(&mut self, _ctx: &mut SimulationContext, _dt: f64, _in_window: bool) {
        panic!("STUB: Router::accumulateLinkStepStats")
    }
}

impl InletSolver {
    /// STUB of InletSolver::computeAll.
    pub fn compute_all(&mut self, _ctx: &mut SimulationContext, _dt: f64) {
        panic!("STUB: InletSolver::computeAll")
    }
}
