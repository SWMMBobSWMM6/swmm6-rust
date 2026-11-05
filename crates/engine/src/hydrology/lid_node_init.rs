// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/LidNode.cpp (initialize,
// initializeInfiltration, the file-local infiltrationDeficit), LidNodeRichards.cpp
// (refreshRichardsState with its file-local refresh) and RichardsColumn.cpp (richards::theta,
// storage, pressure) (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// lidnode::initialize: SWMMEngine::initialize calls it, after lidnode::validate, to build every
// LID node's cells (top to bottom), seed their moisture and the Richards column, and set up
// the Green-Ampt surface infiltration of an ordinary stack.
//
// Differences of form, never in result:
//   * refresh() throws std::runtime_error on a non-finite pressure. Nothing between it and
//     SWMMEngine::initialize (noexcept) catches it, so upstream terminates; here initialize
//     returns Err(Terminated) with that exception and what().

use crate::core::cmath::{cclamp, cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::engine_init_hydraulics::Terminated;
use crate::core::options::Bytes;
use crate::core::unit_conversion as ucf;
use crate::data::lid_node::{richards, LidNodeBottom, LidNodeCell, LidNodeLayer, LidNodeLayerKind, LidNodeState};
use crate::hydraulics::node;
use crate::hydrology::infiltration::grnampt_init;
use crate::hydrology::lid_node::{layers, thickness, validate_stack};

const METRES: f64 = 0.3048;

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

/// richards::theta: van Genuchten water content at pressure h (m).
pub fn richards_theta(p: &richards::Material, ts: f64, h: f64) -> f64 {
    if h >= 0.0 {
        return ts;
    }
    let m = 1.0 - 1.0 / p.n;
    p.theta_r + (ts - p.theta_r) * (1.0 + (-p.alpha * h).powf(p.n)).powf(-m) // LIBM
}

/// richards::storage: water content plus specific storage above zero pressure.
pub fn richards_storage(p: &richards::Material, ts: f64, h: f64) -> f64 {
    richards_theta(p, ts, h) + p.specific_storage * cmax(0.0, h)
}

/// richards::pressure: the inverse of storage.
pub fn richards_pressure(p: &richards::Material, ts: f64, w: f64) -> f64 {
    if w >= ts {
        return (w - ts) / p.specific_storage;
    }
    if w <= p.theta_r {
        return f64::NEG_INFINITY;
    }
    let se = (w - p.theta_r) / (ts - p.theta_r);
    let m = 1.0 - 1.0 / p.n;
    -(-se.ln() / m).exp_m1().powf(1.0 / p.n) / p.alpha // LIBM (log, expm1, pow)
}

/// rainDepth(): display rain-depth units per native ft.
fn rain_depth(ctx: &SimulationContext) -> f64 {
    ucf::UCF[ucf::RAINDEPTH as usize][ucf::get_unit_system(ctx.options.flow_units as i32) as usize]
}

/// LidNodeRichards.cpp refresh: pressure and theta of every porous cell from its water.
fn refresh(s: &mut LidNodeState, stack: &[LidNodeLayer], name: &[u8]) -> Result<(), Terminated> {
    s.held_volume = 0.0;
    for i in 0..s.cells.len() {
        if i == 0 {
            s.cells[i].theta = 0.0;
            continue;
        }
        let p = stack[(s.cells[i].layer - 1) as usize].retention;
        let porosity = s.cells[i].porosity;
        let h = richards_pressure(&p, porosity, s.richards_water[i] / s.cells[i].geometric_volume);
        if !h.is_finite() {
            let what = cat(&[b"Richards non-finite pressure at LID node ", name]);
            return Err(Terminated { exception: "std::runtime_error", what: Box::leak(String::from_utf8_lossy(&what).into_owned().into_boxed_str()) });
        }
        s.richards_pressure[i] = h / METRES;
        s.cells[i].theta = richards_theta(&p, porosity, h);
        s.held_volume += s.richards_water[i];
    }
    Ok(())
}

/// lidnode::refreshRichardsState(ctx, row).
pub fn refresh_richards_state(ctx: &mut SimulationContext, r: i32) -> Result<(), Terminated> {
    let ur = r as usize;
    let stack = layers(ctx, ctx.node_subtypes.storages.lid[ur].control);
    let name = ctx.node_names.name_of(ctx.node_subtypes.storages.node_idx[ur]).to_vec();
    refresh(&mut ctx.node_subtypes.storages.lid_state[ur], &stack, &name)
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
    for cell in &state.cells {
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

/// lidnode::initializeInfiltration(ctx, node).
pub fn initialize_infiltration(ctx: &mut SimulationContext, node: i32) {
    let r = ctx.node_subtypes.storage_row(node) as usize;
    let depth = ctx.nodes.depth[node as usize];
    let us = ucf::get_unit_system(ctx.options.flow_units as i32);
    let rd = rain_depth(ctx);
    let opts = &ctx.options;
    let state = &mut ctx.node_subtypes.storages.lid_state[r];
    state.surface_infil = Default::default();
    state.infiltration_cell = -1;
    state.infiltration_head = depth;
    if state.cells.len() < 2
        || state.cells[0].kind != LidNodeLayerKind::Surface
        || state.cells[1].kind != LidNodeLayerKind::Media
        || state.cells[1].conductivity <= 0.0
    {
        return;
    }
    state.infiltration_cell = 1;
    let entry = state.cells[1];
    grnampt_init(
        &mut state.surface_infil,
        entry.suction * rd,
        entry.conductivity * ucf::UCF[ucf::RAINFALL as usize][us as usize],
        entry.porosity - entry.wilting_point,
        opts,
    );
    // The empirical recovery zone cannot extend past this authored media.
    let mut bottom = entry.bottom;
    for cell in &state.cells {
        if cell.layer == entry.layer {
            bottom = cmin(bottom, cell.bottom);
        }
    }
    let ga = &mut state.surface_infil;
    ga.Lu = cmin(ga.Lu, entry.top - bottom);
    ga.Fumax = ga.IMDmax * ga.Lu;
    let imd = infiltration_deficit(state, depth);
    let ga = &mut state.surface_infil;
    ga.IMD = imd;
    ga.Fu = cclamp((ga.IMDmax - ga.IMD) * ga.Lu, 0.0, ga.Fumax);
    ga.saturated = depth >= entry.top;
}

/// lidnode::initialize(ctx): build every LID node's cells and seed its moisture.
pub fn initialize(ctx: &mut SimulationContext) -> Result<(), Terminated> {
    // without OPENSWMM_HAS_2D a MESH_AQUIFER bottom is refused
    for r in 0..ctx.node_subtypes.storages.count() {
        let ur = r as usize;
        if ctx.node_subtypes.storages.lid[ur].control >= 0 && ctx.node_subtypes.storages.lid[ur].bottom == LidNodeBottom::MeshAquifer {
            let name = ctx.node_names.name_of(ctx.node_subtypes.storages.node_idx[ur]).to_vec();
            ctx.errors.push(cat(&[b"LID node ", &name, b": MESH_AQUIFER bottom requires a build with the 2D groundwater solver"]));
            return Ok(());
        }
    }
    let us = ucf::get_unit_system(ctx.options.flow_units as i32);
    let rd = rain_depth(ctx);
    let rf = ucf::UCF[ucf::RAINFALL as usize][us as usize];
    for r in 0..ctx.node_subtypes.storages.count() {
        let ur = r as usize;
        ctx.node_subtypes.storages.lid_state[ur] = LidNodeState::default();
        if ctx.node_subtypes.storages.lid[ur].control < 0 {
            continue;
        }
        let n = ctx.node_subtypes.storages.node_idx[ur];
        let un = n as usize;
        ctx.nodes.full_volume[un] = 0.0;
        let stack = layers(ctx, ctx.node_subtypes.storages.lid[ur].control);
        if !validate_stack(&stack).is_empty() {
            continue;
        }
        let initial_saturation = ctx.node_subtypes.storages.lid[ur].initial_saturation;
        let mut top = thickness(&stack) / rd;
        let mut cells: Vec<LidNodeCell> = Vec::new();
        let mut layer = 0;
        for l in &stack {
            if l.kind == LidNodeLayerKind::Bottom {
                continue;
            }
            layer += 1;
            let count = if stack[0].flow.enabled && (l.kind == LidNodeLayerKind::Media || l.kind == LidNodeLayerKind::Aggregate) {
                stack[0].flow.cells_per_layer
            } else if l.kind == LidNodeLayerKind::Media {
                5
            } else {
                1
            };
            let dz = l.params[0] / rd / count as f64;
            for _ in 0..count {
                let mut c = LidNodeCell::default();
                c.top = top;
                c.bottom = cmax(0.0, top - dz);
                c.kind = l.kind;
                c.layer = layer;
                // No runtime cells have been installed yet, so these calls evaluate the
                // authored geometric storage relation.
                let vtop = node::get_volume(&ctx.nodes, n, c.top, Some(&ctx.tables), us, Some(&ctx.node_subtypes));
                let vbottom = node::get_volume(&ctx.nodes, n, c.bottom, Some(&ctx.tables), us, Some(&ctx.node_subtypes));
                c.geometric_volume = cmax(0.0, vtop - vbottom);
                c.area = c.geometric_volume / dz;
                c.porosity = if l.kind == LidNodeLayerKind::Surface { 1.0 - l.params[1] } else { l.params[1] };
                if l.kind == LidNodeLayerKind::Media {
                    c.field_capacity = l.params[2];
                    c.wilting_point = l.params[3];
                    c.conductivity = l.params[4] / rf;
                    c.conductivity_slope = l.params[5];
                    c.suction = l.params[6] / rd;
                    c.theta = c.wilting_point + initial_saturation / 100.0 * (c.porosity - c.wilting_point);
                } else if l.kind == LidNodeLayerKind::Aggregate {
                    c.conductivity = l.params[2] / rf;
                    c.theta = initial_saturation / 100.0 * c.porosity;
                }
                if c.area <= 0.0 {
                    let name = ctx.node_names.name_of(n).to_vec();
                    ctx.errors.push(cat(&[b"LID node ", &name, b": layers require positive geometric storage area"]));
                }
                cells.push(c);
                top = c.bottom;
            }
        }
        {
            let state = &mut ctx.node_subtypes.storages.lid_state[ur];
            state.cells = cells;
            state.port_delta.clear();
            state.port_delta.resize(state.cells.len(), 0.0);
        }
        if stack[0].flow.enabled {
            let initial_head = ctx.nodes.depth[un];
            {
                let state = &mut ctx.node_subtypes.storages.lid_state[ur];
                state.richards = true;
                state.richards_water.clear();
                state.richards_water.resize(state.cells.len(), 0.0);
                state.richards_pressure.clear();
                state.richards_pressure.resize(state.cells.len(), 0.0);
                let surface = &state.cells[0];
                ctx.nodes.volume[un] = cmax(0.0, initial_head - surface.bottom) * surface.area * surface.porosity;
                for i in 1..state.cells.len() {
                    let cell = state.cells[i];
                    let p = stack[(cell.layer - 1) as usize].retention;
                    let base = if cell.kind == LidNodeLayerKind::Media { cell.wilting_point } else { p.theta_r };
                    let initial = base + cmax(initial_saturation / 100.0, 1.0e-8) * (cell.porosity - base);
                    let mut h = richards_pressure(&p, cell.porosity, initial);
                    let center = 0.5 * (cell.top + cell.bottom);
                    if initial_head > center {
                        h = (initial_head - center) * 0.3048;
                    }
                    state.richards_water[i] = richards_storage(&p, cell.porosity, h) * cell.geometric_volume;
                }
            }
            refresh_richards_state(ctx, r)?;
            let surface = ctx.node_subtypes.storages.lid_state[ur].cells[0];
            ctx.nodes.depth[un] = surface.bottom + ctx.nodes.volume[un] / (surface.area * surface.porosity);
            ctx.nodes.head[un] = ctx.nodes.invert_elev[un] + ctx.nodes.depth[un];
            ctx.nodes.full_volume[un] = 0.0;
            continue;
        }
        // Fully saturated cells connected to the base belong to the mobile hydraulic store
        // above field capacity, including InitSat=100%.
        let water_table = node::get_depth(&ctx.nodes, n, 0.0, Some(&ctx.tables), us, Some(&ctx.node_subtypes));
        {
            let state = &mut ctx.node_subtypes.storages.lid_state[ur];
            let mut held = state.held_volume;
            for c in state.cells.iter_mut() {
                if c.top <= water_table && c.theta > c.field_capacity {
                    ctx.nodes.volume[un] += (c.theta - c.field_capacity) * c.geometric_volume;
                    c.theta = c.field_capacity;
                }
                held += c.theta * c.geometric_volume;
            }
            state.held_volume = held;
        }
        ctx.nodes.full_volume[un] = 0.0;
        initialize_infiltration(ctx, n);
    }
    Ok(())
}
