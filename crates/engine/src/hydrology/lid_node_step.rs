// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/LidNode.cpp and LidNodeRichards.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The node-LID (storage LID columns, SWMM6-only) calls the routing step makes: active,
// richardsMode, heldVolume, portOffset, portCell, portDepth, exchangePorts, resetPorts,
// prepareStep (prepareStepImpl), finishStep (finishRows), validateRichardsRouting,
// bottomMode, bottomCloggingFactor and the file-local infiltrationDeficit and
// reconcileInfiltration.
//
// Differences in form, never in result:
//   * perf::ProcessTimer (timing only) is left out.
//   * std::max over an initializer list is a left fold of cmax, which is what max_element's
//     `largest < *it` test computes.
//
// NOT YET TRANSLATED: Richards columns (LidNodeRichards.cpp: prepareRichardsStep's column
// advance, exchangeRichardsPorts, refreshRichardsState, the Mesh-aquifer owner checks).
// `unsupported` reports a Richards row, and SWMMEngine::stepRouting
// (core/engine_step_routing.rs) returns StepError::NotTranslated for such a deck before any
// state changes, so every Richards branch below is unreachable and returns NotTranslated.

use crate::core::cmath::{cclamp, cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::engine_step_routing::StepError;
use crate::core::unit_conversion as ucf;
use crate::data::lid_node::{LidNodeBottom, LidNodeLayerKind, LidNodeState, LidPortTransfer, LidWaterTransfer};
use crate::hydraulics::node;
use crate::hydraulics::node_step;
use crate::hydrology::green_ampt_step::grnampt_get_infil;
use crate::hydrology::infiltration::InfilModel;
use crate::hydrology::lid_node::layers;

/// A storage row with a Richards column (not translated yet).
pub fn unsupported(ctx: &SimulationContext) -> bool {
    ctx.node_subtypes.storages.lid_state.iter().any(|s| s.richards)
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

/// infiltrationDeficit: the exposed finite Green-Ampt upper-zone moisture deficit.
fn infiltration_deficit(state: &LidNodeState, head: f64) -> f64 {
    if state.infiltration_cell < 0 {
        return 0.0;
    }
    let entry = &state.cells[state.infiltration_cell as usize];
    let bottom = entry.top - state.surface_infil.Lu;
    let mut empty = 0.0;
    let mut geometric = 0.0;
    for cell in state.cells.iter() {
        if cell.layer != entry.layer {
            continue;
        }
        let lo = cmax(bottom, cell.bottom);
        let hi = cmin(entry.top, cell.top);
        let dz = cmax(0.0, hi - lo);
        geometric += cell.area * dz;
        empty += cell.area * cmax(0.0, hi - cmax(lo, head)) * cmax(0.0, cell.porosity - cell.theta);
    }
    if geometric > 0.0 {
        empty / geometric
    } else {
        0.0
    }
}

/// reconcileInfiltration: Green-Ampt memory after backwater/port wetting and recession.
fn reconcile_infiltration(state: &mut LidNodeState, head: f64, port_wetting: bool) {
    if state.infiltration_cell < 0 {
        return;
    }
    let entry_top = state.cells[state.infiltration_cell as usize].top;
    let tolerance = 32.0 * f64::EPSILON * cmax(1.0, entry_top.abs());
    let submerged = head >= entry_top - tolerance;
    let previously_submerged = state.infiltration_head >= entry_top - tolerance;
    let deficit = if submerged { 0.0 } else { infiltration_deficit(state, head) };
    let ga = &mut state.surface_infil;
    let wet = cclamp((ga.IMDmax - deficit) * ga.Lu, 0.0, ga.Fumax);
    if head > state.infiltration_head + 1.0e-12 || port_wetting {
        ga.IMD = cmin(ga.IMD, deficit);
        ga.Fu = cmax(ga.Fu, wet);
        if submerged {
            ga.IMD = 0.0;
            ga.Fu = ga.Fumax;
            ga.F = 0.0;
            ga.saturated = true;
        }
    }
    if !submerged && (previously_submerged || (head < state.infiltration_head - 1.0e-12 && ga.F == 0.0)) {
        ga.IMD = deficit;
        ga.Fu = wet;
        ga.F = 0.0;
        ga.saturated = false;
    }
    state.infiltration_head = head;
}

/// The capacity and head refresh prepareStepImpl and finishRows both end with.
fn refresh_node(ctx: &mut SimulationContext, n: usize, us: i32) {
    ctx.nodes.full_volume[n] = 0.0;
    let fv = node::get_volume(&ctx.nodes, n as i32, ctx.nodes.full_depth[n], Some(&ctx.tables), us, Some(&ctx.node_subtypes));
    ctx.nodes.full_volume[n] = fv;
    ctx.nodes.rpt_full_volume[n] = ctx.nodes.full_volume[n];
    ctx.nodes.depth[n] = node_step::get_depth(&ctx.nodes, n as i32, ctx.nodes.volume[n], Some(&ctx.tables), us, Some(&ctx.node_subtypes));
    if ctx.options.allow_ponding && ctx.nodes.ponded_area[n] > 0.0 && ctx.nodes.volume[n] > ctx.nodes.full_volume[n] {
        ctx.nodes.depth[n] = ctx.nodes.full_depth[n] + (ctx.nodes.volume[n] - ctx.nodes.full_volume[n]) / ctx.nodes.ponded_area[n];
    }
    ctx.nodes.head[n] = ctx.nodes.invert_elev[n] + ctx.nodes.depth[n];
}

/// lidnode::prepareStep (prepareStepImpl with mesh_prepared false).
pub fn prepare_step(ctx: &mut SimulationContext, dt: f64, evaporation: f64) -> Result<(), StepError> {
    if !(dt > 0.0) {
        return Ok(());
    }
    // prepareRichardsStep: with no Richards row it collects no column and returns
    if unsupported(ctx) {
        return Err(StepError::NotTranslated("Richards node LID columns"));
    }
    let us = ucf::get_unit_system(ctx.options.flow_units as i32);
    let count = ctx.node_subtypes.storages.count();
    for r in 0..count as usize {
        {
            let state = &ctx.node_subtypes.storages.lid_state[r];
            if state.cells.is_empty() || state.richards {
                continue;
            }
        }
        let n = ctx.node_subtypes.storages.node_idx[r] as usize;
        let evap_frac = ctx.node_subtypes.storages.evap_frac[r];
        let infil_factor = ctx.climate_state.infil_factor;
        let recovery_factor = ctx.climate_state.recovery_factor;
        let nodes = &mut ctx.nodes;
        let state = &mut ctx.node_subtypes.storages.lid_state[r];
        let before = state.held_volume;
        state.quality_old_mobile = nodes.volume[n];
        state.quality_old_water.clear();
        state.quality_transfers.clear();
        for i in 0..state.cells.len() {
            let c = &state.cells[i];
            let w = c.theta * c.geometric_volume;
            state.quality_old_water.push(w);
        }
        state.captured_flow = 0.0;
        state.evap_volume = 0.0;
        reconcile_infiltration(state, nodes.depth[n], false);
        for i in 0..state.cells.len() {
            let c = &mut state.cells[i];
            if c.top <= nodes.depth[n] {
                let v = (c.field_capacity - c.theta) * c.geometric_volume;
                if v > 0.0 {
                    state.quality_transfers.push(LidWaterTransfer { from: -1, to: i as i32, volume: v });
                } else if v < 0.0 {
                    state.quality_transfers.push(LidWaterTransfer { from: i as i32, to: -1, volume: -v });
                }
                c.theta = c.field_capacity;
            }
        }

        // Lateral inflows enter the top of the column.
        let capture;
        {
            let first = &mut state.cells[0];
            let capacity = cmax(0.0, (first.porosity - first.theta) * first.geometric_volume);
            capture = cmin(cmax(0.0, nodes.lat_flow[n]) * dt, capacity);
            if first.geometric_volume > 0.0 {
                first.theta += capture / first.geometric_volume;
            }
        }
        state.captured_flow = capture / dt;
        if capture > 0.0 {
            state.quality_transfers.push(LidWaterTransfer { from: -2, to: 0, volume: capture });
        }
        state.treated_volume += capture;
        nodes.lat_flow[n] -= state.captured_flow;
        // Evaporation draws from top to bottom, never below wilting point.
        let mut evap = cmax(0.0, evaporation * evap_frac * state.cells[0].area * dt);
        for i in 0..state.cells.len() {
            let c = &mut state.cells[i];
            let take = cmin(evap, cmax(0.0, (c.theta - c.wilting_point) * c.geometric_volume));
            if c.geometric_volume > 0.0 {
                c.theta -= take / c.geometric_volume;
            }
            state.evap_volume += take;
            evap -= take;
            if take > 0.0 {
                state.quality_transfers.push(LidWaterTransfer { from: i as i32, to: -3, volume: take });
            }
        }
        // Bounded forward-Euler method of lines, in substeps of at most 1 s.
        let ncell = state.cells.len();
        let mut remaining = dt;
        while remaining > 0.0 {
            let step = cmin(remaining, 1.0);
            let mut transfer = vec![0.0f64; ncell];
            for i in 0..ncell {
                let c = state.cells[i];
                if c.top <= nodes.depth[n] {
                    continue;
                }
                let mut q = 0.0;
                if c.kind == LidNodeLayerKind::Surface && i + 1 < ncell {
                    let below = state.cells[i + 1];
                    if state.infiltration_cell == (i + 1) as i32 {
                        let mut trial = state.surface_infil;
                        let factor = infil_factor;
                        let head = c.theta / c.porosity * (c.top - c.bottom);
                        let wet_from_below = nodes.depth[n] > below.top - trial.Lu;
                        if head <= 0.0 && !wet_from_below {
                            trial.saturated = false;
                        }
                        if factor > 0.0 && (head > 0.0 || !wet_from_below) {
                            q = grnampt_get_infil(&mut trial, 0.0, head, step, InfilModel::ModGreenAmpt, factor, cmax(recovery_factor, 1.0e-12)) * c.area;
                        }
                        state.surface_infil.T = trial.T;
                        state.surface_infil.saturated = trial.saturated;
                        if head <= 0.0 && !wet_from_below {
                            let deficit = infiltration_deficit(state, nodes.depth[n]);
                            trial.Fu = cmax(trial.Fu, (trial.IMDmax - deficit) * trial.Lu);
                            trial.IMD = cmin(trial.IMD, deficit);
                            trial.F = cmax(0.0, trial.F);
                            state.surface_infil = trial;
                        }
                    } else if below.kind == LidNodeLayerKind::Aggregate {
                        q = below.conductivity * c.area;
                    }
                } else if c.theta > c.field_capacity {
                    // Legacy soil conductivity slope, not a power exponent.
                    q = c.conductivity
                        * (if c.kind == LidNodeLayerKind::Media {
                            (-c.conductivity_slope * cmax(0.0, c.porosity - c.theta)).exp() // LIBM exp
                        } else {
                            1.0
                        })
                        * c.area;
                }
                let mut v = cmin(q * step, cmax(0.0, (c.theta - c.field_capacity) * c.geometric_volume));
                if i + 1 < ncell {
                    let below = &state.cells[i + 1];
                    if below.bottom >= nodes.depth[n] {
                        v = cmin(v, cmax(0.0, (below.porosity - below.theta) * below.geometric_volume));
                    }
                }
                transfer[i] = v;
                if c.kind == LidNodeLayerKind::Surface && state.infiltration_cell == (i + 1) as i32 {
                    let ga = &mut state.surface_infil;
                    let accepted = v / c.area;
                    ga.F += accepted;
                    ga.Fu = cmin(ga.Fumax, ga.Fu + accepted);
                }
            }
            for i in 0..ncell {
                let dest = if i + 1 < ncell && state.cells[i + 1].bottom >= nodes.depth[n] { (i + 1) as i32 } else { -1 };
                if transfer[i] > 0.0 {
                    state.quality_transfers.push(LidWaterTransfer { from: i as i32, to: dest, volume: transfer[i] });
                }
                {
                    let c = &mut state.cells[i];
                    if c.geometric_volume > 0.0 {
                        c.theta -= transfer[i] / c.geometric_volume;
                    }
                }
                if i + 1 < ncell {
                    let below = &mut state.cells[i + 1];
                    if below.bottom >= nodes.depth[n] && below.geometric_volume > 0.0 {
                        below.theta += transfer[i] / below.geometric_volume;
                    }
                }
            }
            remaining -= step;
        }
        state.held_volume = 0.0;
        for i in 0..ncell {
            let c = &state.cells[i];
            state.held_volume += c.theta * c.geometric_volume;
        }
        nodes.volume[n] += before + capture - state.evap_volume - state.held_volume;
        // Independent retained/mobile sums can leave an empty mobile store a few ULPs below zero.
        let mut m = 1.0;
        for x in [state.quality_old_mobile.abs(), before, capture, state.evap_volume, state.held_volume] {
            m = cmax(m, x);
        }
        let roundoff = 16.0 * f64::EPSILON * m;
        if nodes.volume[n] < 0.0 && nodes.volume[n] >= -roundoff {
            nodes.volume[n] = 0.0;
        }
        // An internal redistribution before the routing step: refresh capacity and head.
        refresh_node(ctx, n, us);
    }
    Ok(())
}

/// finishRows over every row (no Richards row reaches here).
fn finish_rows(ctx: &mut SimulationContext, rows: &[i32]) {
    let us = ucf::get_unit_system(ctx.options.flow_units as i32);
    for &r in rows {
        let ur = r as usize;
        if ctx.node_subtypes.storages.lid_state[ur].cells.is_empty() {
            continue;
        }
        let n = ctx.node_subtypes.storages.node_idx[ur] as usize;
        {
            let state = &mut ctx.node_subtypes.storages.lid_state[ur];
            // (richards_mesh_pending_finish: Richards only)
            for i in 0..state.cells.len() {
                let pd = state.port_delta[i];
                let cell = &mut state.cells[i];
                if cell.geometric_volume > 0.0 {
                    cell.theta += pd / cell.geometric_volume;
                }
            }
            state.held_volume = 0.0;
            for i in 0..state.cells.len() {
                let cell = &state.cells[i];
                state.held_volume += cell.theta * cell.geometric_volume;
            }
        }
        refresh_node(ctx, n, us);
        let depth = ctx.nodes.depth[n];
        let st = &mut ctx.node_subtypes.storages;
        let node_r = st.node_idx[ur] as usize;
        let state = &mut st.lid_state[ur];
        let port_wetting = state.infiltration_cell >= 0
            && state.quality_ports.iter().any(|port| port.volume > 0.0 && state.cells[port.cell as usize].layer == state.cells[state.infiltration_cell as usize].layer);
        reconcile_infiltration(state, depth, port_wetting);
        ctx.nodes.lat_flow[node_r] += state.captured_flow;
        st.evap_loss[ur] += st.lid_state[ur].evap_volume;
        st.exfil_loss[ur] += st.lid_state[ur].richards_bottom_loss;
        st.lid_state[ur].captured_flow = 0.0;
    }
}

/// lidnode::finishStep
pub fn finish_step(ctx: &mut SimulationContext) -> Result<(), StepError> {
    let count = ctx.node_subtypes.storages.count();
    let mut rows: Vec<i32> = Vec::new();
    for r in 0..count as usize {
        rows.push(r as i32);
        let st = &ctx.node_subtypes.storages;
        if !st.lid_state[r].richards || st.lid[r].bottom != LidNodeBottom::MeshAquifer {
            continue;
        }
        return Err(StepError::NotTranslated("Richards node LID columns"));
    }
    finish_rows(ctx, &rows);
    Ok(())
}

/// lidnode::validateRichardsRouting: false without mesh-aquifer Richards columns.
pub fn validate_richards_routing(ctx: &SimulationContext, _dt: f64) -> Result<bool, StepError> {
    let st = &ctx.node_subtypes.storages;
    for r in 0..st.count() as usize {
        if !st.lid_state[r].richards || st.lid[r].bottom != LidNodeBottom::MeshAquifer {
            continue;
        }
        return Err(StepError::NotTranslated("Richards node LID columns"));
    }
    Ok(false)
}

/// lidnode::bottomMode
pub fn bottom_mode(c: &SimulationContext, node: i32) -> LidNodeBottom {
    if node < 0 || node >= c.n_nodes() {
        return LidNodeBottom::Current;
    }
    let r = c.node_subtypes.storage_row(node);
    if r >= 0 && c.node_subtypes.storages.lid[r as usize].control >= 0 {
        c.node_subtypes.storages.lid[r as usize].bottom
    } else {
        LidNodeBottom::Current
    }
}

/// lidnode::bottomCloggingFactor
pub fn bottom_clogging_factor(c: &SimulationContext, r: i32) -> f64 {
    if r < 0 || c.node_subtypes.storages.lid[r as usize].control < 0 {
        return 1.0;
    }
    let ur = r as usize;
    let st = &c.node_subtypes.storages;
    let stack = layers(c, st.lid[ur].control);
    match stack.last() {
        Some(b) if b.kind == LidNodeLayerKind::Bottom && b.params[1] > 0.0 => {
            let mut pores = 0.0;
            for cell in st.lid_state[ur].cells.iter() {
                if cell.kind != LidNodeLayerKind::Surface {
                    pores += cell.porosity * cell.geometric_volume;
                }
            }
            if pores > 0.0 {
                cmax(0.0, 1.0 - st.lid_state[ur].treated_volume / (b.params[1] * pores))
            } else {
                1.0
            }
        }
        _ => 1.0,
    }
}

/// lidnode::portCell: the cell a port at `offset` belongs to (an interface belongs to the cell
/// above it).
pub fn port_cell(state: &LidNodeState, offset: f64) -> i32 {
    for i in 0..state.cells.len() {
        let cell = &state.cells[i];
        let mut scale = 1.0;
        for x in [offset.abs(), cell.bottom.abs(), cell.top.abs()] {
            scale = cmax(scale, x);
        }
        let tolerance = 32.0 * f64::EPSILON * scale;
        if offset >= cell.bottom - tolerance && offset < cell.top - tolerance {
            return i as i32;
        }
    }
    if state.cells.is_empty() {
        -1
    } else {
        0
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

/// lidnode::portDepth: the head a structure or conduit end sees at its port.
pub fn port_depth(ctx: &SimulationContext, node: i32, offset: f64) -> f64 {
    let head = ctx.nodes.depth[node as usize];
    let r = ctx.node_subtypes.storage_row(node);
    if r < 0 {
        return head;
    }
    let state = &ctx.node_subtypes.storages.lid_state[r as usize];
    let i = port_cell(state, offset);
    if state.richards && i > 0 {
        let ui = i as usize;
        return 0.5 * (state.cells[ui].bottom + state.cells[ui].top) + state.richards_pressure[ui];
    }
    if i < 0 || state.cells[i as usize].kind != LidNodeLayerKind::Surface || state.cells[i as usize].theta <= 0.0 {
        return head;
    }
    let c = &state.cells[i as usize];
    cmax(head, c.bottom + c.theta / c.porosity * (c.top - c.bottom))
}

/// lidnode::exchangePorts: moves a link's flow through the retained/mobile stores of LID
/// nodes at its ends; returns the accepted flow.
pub fn exchange_ports(ctx: &mut SimulationContext, link: i32, flow: f64, dt: f64) -> Result<f64, StepError> {
    if !(dt > 0.0) || flow == 0.0 {
        return Ok(flow);
    }
    let ul = link as usize;
    let n1 = ctx.links.node1[ul];
    let n2 = ctx.links.node2[ul];
    if richards_mode(ctx, n1) || richards_mode(ctx, n2) {
        return Err(StepError::NotTranslated("Richards node LID ports"));
    }
    if !active(ctx, n1) && !active(ctx, n2) {
        return Ok(flow);
    }
    let source = if flow > 0.0 { n1 } else { n2 };
    let dest = if flow > 0.0 { n2 } else { n1 };
    let source_offset = port_offset(ctx, link, source);
    let dest_offset = port_offset(ctx, link, dest);
    let us_ = source as usize;
    let ud = dest as usize;
    let mut volume = flow.abs() * dt;
    let mut sr = ctx.node_subtypes.storage_row(source);
    if active(ctx, source) && source_offset < ctx.nodes.depth[us_] {
        let state = &mut ctx.node_subtypes.storages.lid_state[sr as usize];
        let available = cmax(0.0, ctx.nodes.old_volume[us_] + (ctx.nodes.lat_flow[us_] - ctx.nodes.losses[us_]) * dt + state.mobile_delta);
        volume = cmin(volume, available);
        state.mobile_delta -= volume;
    }

    if sr >= 0 && source_offset >= ctx.nodes.depth[us_] {
        let s = &mut ctx.node_subtypes.storages.lid_state[sr as usize];
        let i = port_cell(s, source_offset);
        if i >= 0 {
            let ui = i as usize;
            let c = s.cells[ui];
            let available = cmax(0.0, (c.theta - c.wilting_point) * c.geometric_volume + s.port_delta[ui]);
            let throttle = if c.kind == LidNodeLayerKind::Media { cclamp(c.theta / c.porosity, 0.0, 1.0) } else { 1.0 };
            volume = cmin(volume * throttle, available);
            s.port_delta[ui] -= volume;
            s.quality_ports.push(LidPortTransfer { link, cell: i, volume: -volume });
            // The routing equation handles only mobile water: credit the same flux internally.
            ctx.nodes.inflow[us_] += volume / dt;
        }
    }
    sr = ctx.node_subtypes.storage_row(dest);
    let mut mobile_received = volume;
    if sr >= 0 && dest_offset >= ctx.nodes.depth[ud] {
        let depth_d = ctx.nodes.depth[ud];
        let s = &mut ctx.node_subtypes.storages.lid_state[sr as usize];
        let i = port_cell(s, dest_offset);
        if i >= 0 && s.cells[i as usize].bottom >= depth_d {
            let ui = i as usize;
            let c = s.cells[ui];
            let capacity = cmax(0.0, (c.porosity - c.theta) * c.geometric_volume - s.port_delta[ui]);
            let held = cmin(volume, capacity);
            s.port_delta[ui] += held;
            s.quality_ports.push(LidPortTransfer { link, cell: i, volume: held });
            mobile_received -= held;
            ctx.nodes.outflow[ud] += held / dt;
        }
    }
    if active(ctx, dest) {
        ctx.node_subtypes.storages.lid_state[sr as usize].mobile_delta += mobile_received;
    }
    Ok((volume / dt).copysign(flow))
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
