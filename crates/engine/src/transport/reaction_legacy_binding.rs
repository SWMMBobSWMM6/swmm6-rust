// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/components/ReactionModule/
// ReactionLegacyBinding.hpp and ReactionLegacyBinding.cpp (HydroCouple/openswmm.engine,
// branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache
// License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The open-time part of the reactions' LEGACY binding: legacyReactionsActive and the
// IGNORE_QUALITY bypass warning. The run-time binding (species seeding, decay, kinetics)
// is not translated here.

use crate::core::context::SimulationContext;

/// A configured and compiled reactions component.
pub fn legacy_reactions_active(ctx: &SimulationContext) -> bool {
    ctx.reactions.configured && ctx.reactions.compiled
}

/// Warn when configured compiled reactions are bypassed by IGNORE_QUALITY.
pub fn warn_if_legacy_binding_bypassed(ctx: &mut SimulationContext) {
    if !legacy_reactions_active(ctx) {
        return;
    }
    if ctx.options.ignore_quality {
        ctx.warnings.push(
            "A reactions component is configured but IGNORE_QUALITY is YES \u{2014} the quality stage does not run, so no species react and no pollutant decays this simulation."
                .as_bytes()
                .to_vec(),
        );
    }
}
