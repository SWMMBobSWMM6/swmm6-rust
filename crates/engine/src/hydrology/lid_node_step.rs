// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/LidNode.cpp and LidNodeRichards.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The node-LID (storage LID columns, SWMM6-only) calls the routing step makes: active,
// heldVolume, portOffset, portDepth, exchangePorts, resetPorts, prepareStep, finishStep,
// validateRichardsRouting, richardsMode.
//
// NOT YET TRANSLATED: a storage node with LID cells or a Richards column. `unsupported`
// reports one, and SWMMEngine::stepRouting (core/engine_step_routing.rs) returns
// StepError::NotTranslated for such a deck before any state changes. Every function here is
// therefore exact for the remaining decks, where no storage row has cells: active() is false,
// portDepth() is the node depth, exchangePorts() returns the flow, prepareStep/finishStep do
// nothing, and resetPorts() still clears every row's port state as upstream does.

use crate::core::context::SimulationContext;
use crate::core::engine_step_routing::StepError;

/// A storage row with LID cells or a Richards column (not translated yet).
pub fn unsupported(ctx: &SimulationContext) -> bool {
    ctx.node_subtypes.storages.lid_state.iter().any(|s| !s.cells.is_empty() || s.richards)
}

/// lidnode::active
pub fn active(ctx: &SimulationContext, node: i32) -> bool {
    let r = ctx.node_subtypes.storage_row(node);
    r >= 0 && !ctx.node_subtypes.storages.lid_state[r as usize].cells.is_empty()
}

/// lidnode::richardsMode
pub fn richards_mode(ctx: &SimulationContext, node: i32) -> bool {
    let r = ctx.node_subtypes.storage_row(node);
    r >= 0 && ctx.node_subtypes.storages.lid_state[r as usize].richards
}

/// lidnode::heldVolume
pub fn held_volume(ctx: &SimulationContext, node: i32) -> f64 {
    let r = ctx.node_subtypes.storage_row(node);
    if r >= 0 {
        ctx.node_subtypes.storages.lid_state[r as usize].held_volume
    } else {
        0.0
    }
}

/// lidnode::portOffset: the port elevation of a link at one of its nodes.
pub fn port_offset(ctx: &SimulationContext, link: i32, node: i32) -> f64 {
    let ul = link as usize;
    let mut offset = ctx.links.offset1[ul];
    let wr = ctx.link_subtypes.weir_row(link);
    let out = ctx.link_subtypes.outlet_row(link);
    if wr >= 0 {
        offset = ctx.link_subtypes.weirs.crest_height[wr as usize];
    } else if out >= 0 {
        offset = ctx.link_subtypes.outlets.crest_height[out as usize];
    } else {
        return if ctx.links.node1[ul] == node { offset } else { ctx.links.offset2[ul] };
    }
    offset + (ctx.nodes.invert_elev[ctx.links.node1[ul] as usize] - ctx.nodes.invert_elev[node as usize])
}

/// lidnode::portDepth. With no LID cells at the node (portCell gives -1) it is the node depth.
pub fn port_depth(ctx: &SimulationContext, node: i32, _offset: f64) -> f64 {
    ctx.nodes.depth[node as usize]
}

/// lidnode::exchangePorts. With no LID node at either end it returns the flow unchanged.
pub fn exchange_ports(ctx: &mut SimulationContext, link: i32, flow: f64, dt: f64) -> Result<f64, StepError> {
    if !(dt > 0.0) || flow == 0.0 {
        return Ok(flow);
    }
    let n1 = ctx.links.node1[link as usize];
    let n2 = ctx.links.node2[link as usize];
    if richards_mode(ctx, n1) || richards_mode(ctx, n2) || active(ctx, n1) || active(ctx, n2) {
        return Err(StepError::NotTranslated("node LID ports"));
    }
    Ok(flow)
}

/// lidnode::resetPorts
pub fn reset_ports(ctx: &mut SimulationContext) -> Result<(), StepError> {
    if ctx.node_subtypes.storages.lid_state.iter().any(|s| s.richards_mesh_routed) {
        return Err(StepError::Exception("Mesh Richards final routing ports must finish before reset."));
    }
    for s in ctx.node_subtypes.storages.lid_state.iter_mut() {
        s.port_delta.iter_mut().for_each(|x| *x = 0.0);
        s.mobile_delta = 0.0;
        s.quality_ports.clear();
    }
    Ok(())
}

/// lidnode::prepareStep (a no-op without LID cells or Richards columns).
pub fn prepare_step(ctx: &mut SimulationContext, _dt: f64, _evaporation: f64) -> Result<(), StepError> {
    if unsupported(ctx) {
        return Err(StepError::NotTranslated("node LID columns"));
    }
    Ok(())
}

/// lidnode::finishStep (a no-op without LID cells or Richards columns).
pub fn finish_step(ctx: &mut SimulationContext) -> Result<(), StepError> {
    if unsupported(ctx) {
        return Err(StepError::NotTranslated("node LID columns"));
    }
    Ok(())
}

/// lidnode::validateRichardsRouting: false without mesh-aquifer Richards columns.
pub fn validate_richards_routing(ctx: &SimulationContext, _dt: f64) -> Result<bool, StepError> {
    if ctx.node_subtypes.storages.lid_state.iter().any(|s| s.richards) {
        return Err(StepError::NotTranslated("Richards node LID columns"));
    }
    Ok(false)
}
