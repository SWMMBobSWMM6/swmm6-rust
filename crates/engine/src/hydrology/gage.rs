// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/Gage.hpp and Gage.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The part of the rain-gage module that SWMMEngine::open reaches: gageIsUsed. The rest
// (rainfall state machine, conversion) is translated with hydrology (phase 1c).
//
// 1D build: the 2D-mesh clause (every gage is used under a rain-on-grid RAINFALL_MODE)
// tests ctx.twod_io.mesh, which a 1D build never sets, so it is always false here.

use crate::core::context::SimulationContext;

/// True when a subcatchment or an RDII unit hydrograph reads the gage.
pub fn gage_is_used(ctx: &SimulationContext, gage_idx: i32) -> bool {
    for i in 0..ctx.n_subcatches() {
        if ctx.subcatches.gage[i as usize] == gage_idx {
            return true;
        }
    }
    for gname in &ctx.unit_hyds.gage_names {
        if ctx.gage_names.find(gname) == gage_idx {
            return true;
        }
    }
    false
}
