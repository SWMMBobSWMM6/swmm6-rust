// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/LidNode.cpp (active, heldVolume) and
// LidNodeTreatment.cpp (heldMass) (HydroCouple/openswmm.engine, branch swmm6_rel, commit
// 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The retained-water and retained-pollutant reads of LID storage nodes (namespace lidnode)
// that SWMMEngine::reportedNodeVolume and SWMMEngine::initMassBalance make. The rest of
// LidNode.cpp / LidNodeTreatment.cpp is in hydrology::lid_node (readers, validate) or not yet
// translated (see that file's list).

use crate::core::context::SimulationContext;

/// lidnode::active: the node has a storage row whose LID state holds cells.
pub fn active(ctx: &SimulationContext, node: i32) -> bool {
    let r = ctx.node_subtypes.storage_row(node);
    r >= 0 && !ctx.node_subtypes.storages.lid_state[r as usize].cells.is_empty()
}

/// lidnode::heldVolume: retained moisture of the node's LID cells (ft3), 0 without a row.
pub fn held_volume(ctx: &SimulationContext, node: i32) -> f64 {
    let r = ctx.node_subtypes.storage_row(node);
    if r >= 0 {
        ctx.node_subtypes.storages.lid_state[r as usize].held_volume
    } else {
        0.0
    }
}

/// lidnode::heldMass (LidNodeTreatment.cpp): pollutant `p` held in the node's LID cells, or,
/// without per-cell quality mass, the held volume at the node concentration.
pub fn held_mass(ctx: &SimulationContext, node: i32, p: i32) -> f64 {
    let r = ctx.node_subtypes.storage_row(node);
    if r < 0 {
        return 0.0;
    }
    let state = &ctx.node_subtypes.storages.lid_state[r as usize];
    let mut mass = 0.0;
    let np = ctx.n_pollutants();
    if np <= 0 || p < 0 || p >= np {
        return 0.0;
    }
    let mut i = p as usize;
    while i < state.quality_mass.len() {
        mass += state.quality_mass[i];
        i += np as usize;
    }
    if state.quality_mass.is_empty() && np > 0 && node * np + p < ctx.nodes.conc.len() as i32 {
        mass = held_volume(ctx, node) * ctx.nodes.conc[(node * np + p) as usize];
    }
    mass
}
