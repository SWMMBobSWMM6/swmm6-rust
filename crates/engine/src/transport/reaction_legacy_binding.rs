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
// IGNORE_QUALITY bypass warning; and the init-time ensureMsxState (SWMMEngine::initQuality
// seeds the element species state with it). The run-time binding (decay, kinetics,
// noteNegativeClamp, clampRateSpecies, the scratch workspace) is not translated here.

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

/// ensureMsxState: lazily size and seed the reaction species element state. Only the node
/// vector's shape triggers the reset; a matching node shape with a wrong link shape is not
/// repaired.
pub fn ensure_msx_state(ctx: &mut SimulationContext) {
    let ns = ctx.reactions.n_species() as usize;
    let nn = ctx.n_nodes() as usize;
    let nl = ctx.n_links() as usize;
    let want_nodes = nn * ns;
    let want_links = nl * ns;
    let rx = &mut ctx.reactions;
    if rx.msx_node_conc.len() != want_nodes {
        rx.msx_node_conc.clear();
        rx.msx_node_conc.resize(want_nodes, 0.0);
        rx.msx_link_conc.clear();
        rx.msx_link_conc.resize(want_links, 0.0);
        for e in 0..nn {
            for s in 0..ns {
                rx.msx_node_conc[e * ns + s] = rx.init_global[s];
            }
        }
        for e in 0..nl {
            for s in 0..ns {
                rx.msx_link_conc[e * ns + s] = rx.init_global[s];
            }
        }
        // E-B2: [REACTION_QUALITY] NODE/LINK rows override the GLOBAL fill.
        for k in 0..rx.init_elem_idx.len() {
            let e = rx.init_elem_idx[k] as usize;
            let s = rx.init_elem_species[k] as usize;
            if s >= ns {
                continue;
            }
            let v = rx.init_elem_value[k];
            let arr = if rx.init_elem_is_link[k] != 0 { &mut rx.msx_link_conc } else { &mut rx.msx_node_conc };
            let i = e.wrapping_mul(ns).wrapping_add(s);
            if i < arr.len() {
                arr[i] = v;
            }
        }
    }
}
