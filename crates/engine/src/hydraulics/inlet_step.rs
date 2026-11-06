// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Inlet.cpp (HydroCouple/openswmm.engine,
// branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache
// License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The street-inlet solver's routing-step call (inlet.rs holds its state and init).
//
// NOT YET TRANSLATED: InletSolver::computeAll for a deck with inlets. stepRouting refuses
// such a deck (StepError::NotTranslated) before assembleLateralInflows reaches this call, so
// only the no-inlet early return is ever taken here.

use crate::core::context::SimulationContext;
use crate::hydraulics::inlet::InletSolver;

impl InletSolver {
    /// InletSolver::computeAll (the no-inlet case; see the file notes).
    pub fn compute_all(&mut self, _ctx: &mut SimulationContext, _dt: f64) {
        if self.soa_.count == 0 {
            return;
        }
        // unreachable through step_routing (refused there)
    }
}
