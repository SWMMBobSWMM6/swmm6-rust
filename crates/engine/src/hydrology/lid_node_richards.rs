// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/LidNodeRichards.cpp and
// LidNodeRichardsState.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 6 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The routing step's Richards node-LID columns: exchangeRichardsPorts, validateRichardsRouting,
// prepareRichardsStep and the detail namespace (checkMeshInterval, column, draft, publish,
// routed). readRichards is in hydrology/lid_node.rs, refresh and refreshRichardsState in
// hydrology/lid_node_init.rs, richardsMode in hydrology/lid_node_step.rs, the column kernel
// in hydrology/richards_column.rs.
//
// Differences in form, never in result:
//   * A std::runtime_error is StepError::Exception with upstream's what(); a message built
//     at run time (a node name in it) is leaked to a &'static str, as engine_step_runoff.rs
//     does.
//   * RichardsNodeDraft owns a copy of the row's LidNodeState; publish swaps it in.
//   * std::max over an initializer list is a left fold of cmax (max_element's
//     `largest < *it` test).
//
// The Mesh-aquifer branches (bottom MESH_AQUIFER) are translated as they stand. In a 1D
// build nothing sets richards_mesh_pending_finish (only LidNodeAquifer.cpp's 2D owner does),
// so a Mesh-aquifer Richards row stops at prepareRichardsStep's runtime_error, as upstream's
// does.

use crate::core::cmath::{cclamp, cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::engine_init_hydraulics::Terminated;
use crate::core::engine_step_routing::StepError;
use crate::core::unit_conversion as ucf;
use crate::data::lid_node::richards::Report;
use crate::data::lid_node::{LidNodeBottom, LidNodeLayerKind, LidNodeState, LidPortTransfer, LidWaterTransfer};
use crate::hydraulics::node;
use crate::hydrology::lid_node::layers;
use crate::hydrology::lid_node_init::refresh;
use crate::hydrology::lid_node_step::{active, bottom_clogging_factor, bottom_mode, port_cell, port_depth, port_offset, richards_mode};
use crate::hydrology::richards_column::{self as richards, Bottom, Cell, Column};

const METRES: f64 = 0.3048;
const CUBIC_METRES: f64 = 0.028316846592;

/// A std::runtime_error with a message built at run time.
fn thrown(what: Vec<u8>) -> StepError {
    StepError::Exception(Box::leak(String::from_utf8_lossy(&what).into_owned().into_boxed_str()))
}

/// refresh's runtime_error ("Richards non-finite pressure at LID node ...") on the step side.
fn from_refresh(t: Terminated) -> StepError {
    StepError::Exception(t.what)
}

fn max4(a: f64, b: f64, c: f64, d: f64) -> f64 {
    cmax(cmax(cmax(a, b), c), d)
}

/// lidnode::detail::checkMeshInterval
pub fn check_mesh_interval(s: &LidNodeState, start: f64, end: f64) -> Result<(), StepError> {
    if !s.richards
        || !s.richards_mesh_pending_finish
        || !start.is_finite()
        || start < 0.0
        || !end.is_finite()
        || end <= start
        || !(end - start).is_finite()
        || s.richards_mesh_start != start
        || s.richards_completed_time != end
    {
        return Err(StepError::Exception("Mesh Richards routing requires the matching accepted owner interval."));
    }
    Ok(())
}

/// lidnode::exchangeRichardsPorts: a link's flow through the Richards (and ordinary) stores
/// of the LID nodes at its ends; returns the accepted flow.
pub fn exchange_richards_ports(ctx: &mut SimulationContext, link: i32, flow: f64, dt: f64) -> Result<f64, StepError> {
    let ul = link as usize;
    let mut mesh: Option<usize> = None;
    for n in [ctx.links.node1[ul], ctx.links.node2[ul]] {
        let r = ctx.node_subtypes.storage_row(n);
        if r < 0 || !richards_mode(ctx, n) || bottom_mode(ctx, n) != LidNodeBottom::MeshAquifer {
            continue;
        }
        let st = &ctx.node_subtypes.storages;
        let s = &st.lid_state[r as usize];
        let (start, end) = match mesh {
            Some(m) => (st.lid_state[m].richards_mesh_start, st.lid_state[m].richards_completed_time),
            None => (s.richards_mesh_start, s.richards_completed_time),
        };
        check_mesh_interval(s, start, end)?;
        if s.richards_mesh_routed {
            return Err(StepError::Exception("Mesh Richards final routing ports await finish."));
        }
        let interval = s.richards_completed_time - s.richards_mesh_start;
        if !dt.is_finite() || dt <= 0.0 || (dt - interval).abs() > 8.0 * f64::EPSILON * max4(1., dt, s.richards_mesh_start, s.richards_completed_time) {
            return Err(StepError::Exception("Mesh Richards port duration differs from its accepted owner interval."));
        }
        mesh = Some(r as usize);
    }
    if mesh.is_some() && !flow.is_finite() {
        return Err(StepError::Exception("Nonfinite Mesh Richards port flow."));
    }
    let source = if flow > 0.0 { ctx.links.node1[ul] } else { ctx.links.node2[ul] };
    let dest = if flow > 0.0 { ctx.links.node2[ul] } else { ctx.links.node1[ul] };
    let sr = ctx.node_subtypes.storage_row(source);
    let dr = ctx.node_subtypes.storage_row(dest);
    let so = port_offset(ctx, link, source);
    let doff = port_offset(ctx, link, dest);
    let source_active = active(ctx, source);
    let dest_active = active(ctx, dest);
    let si = if source_active { port_cell(&ctx.node_subtypes.storages.lid_state[sr as usize], so) } else { -1 };
    let di = if dest_active { port_cell(&ctx.node_subtypes.storages.lid_state[dr as usize], doff) } else { -1 };
    let source_held = source_active && if richards_mode(ctx, source) { si > 0 } else { so >= ctx.nodes.depth[source as usize] };
    let dest_held = dest_active && if richards_mode(ctx, dest) { di > 0 } else { doff >= ctx.nodes.depth[dest as usize] };
    let mut volume = flow.abs() * dt;
    if source_active {
        let s = &ctx.node_subtypes.storages.lid_state[sr as usize];
        if source_held {
            let usi = si as usize;
            let c = s.cells[usi];
            let available;
            if s.richards {
                let p = layers(ctx, ctx.node_subtypes.storages.lid[sr as usize].control)[(c.layer - 1) as usize].retention;
                let cutoff = (so - 0.5 * (c.top + c.bottom)) * METRES;
                available = s.richards_water[usi] + s.port_delta[usi] - richards::storage(&p, c.porosity, cutoff) * c.geometric_volume;
            } else {
                available = (c.theta - c.wilting_point) * c.geometric_volume + s.port_delta[usi];
                if c.kind == LidNodeLayerKind::Media {
                    volume *= cclamp(c.theta / c.porosity, 0., 1.);
                }
            }
            volume = cmin(volume, cmax(0., available));
        } else {
            let us_ = source as usize;
            volume = cmin(volume, cmax(0., ctx.nodes.old_volume[us_] + (ctx.nodes.lat_flow[us_] - ctx.nodes.losses[us_]) * dt + s.mobile_delta));
        }
    }
    if dest_held && richards_mode(ctx, dest) {
        let s = &ctx.node_subtypes.storages.lid_state[dr as usize];
        let udi = di as usize;
        let c = s.cells[udi];
        let p = layers(ctx, ctx.node_subtypes.storages.lid[dr as usize].control)[(c.layer - 1) as usize].retention;
        let h = ctx.nodes.invert_elev[source as usize] + port_depth(ctx, source, so) - ctx.nodes.invert_elev[dest as usize];
        let target = (h - 0.5 * (c.top + c.bottom)) * METRES;
        let capacity = richards::storage(&p, c.porosity, target) * c.geometric_volume - s.richards_water[udi] - s.port_delta[udi];
        volume = cmin(volume, cmax(0., capacity));
    }
    if source_active {
        let s = &mut ctx.node_subtypes.storages.lid_state[sr as usize];
        if source_held {
            s.port_delta[si as usize] -= volume;
            s.quality_ports.push(LidPortTransfer { link, cell: si, volume: -volume });
            ctx.nodes.inflow[source as usize] += volume / dt;
        } else {
            s.mobile_delta -= volume;
        }
    }
    if dest_active {
        let s = &mut ctx.node_subtypes.storages.lid_state[dr as usize];
        let mut held = 0.0;
        if dest_held {
            let udi = di as usize;
            let c = s.cells[udi];
            held = if s.richards { volume } else { cmin(volume, cmax(0., (c.porosity - c.theta) * c.geometric_volume - s.port_delta[udi])) };
            s.port_delta[udi] += held;
            s.quality_ports.push(LidPortTransfer { link, cell: di, volume: held });
            ctx.nodes.outflow[dest as usize] += held / dt;
        }
        s.mobile_delta += volume - held;
    }
    Ok((volume / dt).copysign(flow))
}

/// lidnode::detail::column: one initialized native node stack as an owned SI column.
pub fn column(ctx: &SimulationContext, r: i32, dt: f64, evaporation: f64, continuous: bool) -> Column {
    let ur = r as usize;
    let st = &ctx.node_subtypes.storages;
    let s = &st.lid_state[ur];
    let node = st.node_idx[ur] as usize;
    let stack = layers(ctx, st.lid[ur].control);
    let mut col = Column { options: stack[0].flow, ..Column::default() };
    let capture = cmax(0., ctx.nodes.lat_flow[node]) * dt;
    for i in 0..s.cells.len() {
        let c = &s.cells[i];
        let mut cell = Cell::default();
        cell.volume = c.geometric_volume * CUBIC_METRES;
        cell.area = c.area * METRES * METRES;
        cell.z = (if i != 0 { 0.5 * (c.top + c.bottom) } else { c.bottom }) * METRES;
        cell.dz = (c.top - c.bottom) * METRES;
        cell.theta_s = c.porosity;
        cell.ks = c.conductivity * METRES * cmax(0., ctx.climate_state.infil_factor);
        cell.wilting = c.wilting_point;
        cell.material = stack[(c.layer - 1) as usize].retention;
        if i == 0 {
            cell.volume *= c.porosity;
        }
        col.cells.push(cell);
        col.water.push((if i != 0 { s.richards_water[i] } else { ctx.nodes.volume[node] + if continuous { 0. } else { capture } }) * CUBIC_METRES);
    }
    if continuous {
        col.inflow = cmax(0., ctx.nodes.lat_flow[node]) * CUBIC_METRES;
    }
    let demand = cmax(0., evaporation * st.evap_frac[ur] * s.cells[0].area) * CUBIC_METRES;
    if col.water[0] > 0.0 {
        col.cells[0].potential_et = demand;
    } else {
        let mut thickness = 0.0;
        for i in 1..col.cells.len() {
            thickness += col.cells[i].dz;
        }
        for i in 1..col.cells.len() {
            col.cells[i].potential_et = demand * col.cells[i].dz / thickness;
        }
    }
    let us = ucf::get_unit_system(ctx.options.flow_units as i32);
    let back = stack.last().unwrap();
    if st.lid[ur].bottom == LidNodeBottom::Current && back.kind == LidNodeLayerKind::Bottom && back.params[0] > 0.0 {
        col.bottom = Bottom::FreeDrainage;
        col.bottom_k = back.params[0] / ucf::UCF[ucf::RAINFALL as usize][us as usize] * METRES * bottom_clogging_factor(ctx, r);
    }
    col
}

/// lidnode::detail::RichardsNodeDraft
pub struct RichardsNodeDraft {
    pub row: i32,
    pub node: i32,
    pub state: LidNodeState,
    pub volume: f64,
    pub depth: f64,
    pub head: f64,
    pub full_volume: f64,
    pub lat_flow: f64,
}

/// lidnode::detail::draft: an owned native publication from an accepted SI column/report.
pub fn draft(ctx: &SimulationContext, r: i32, col: &Column, report: &Report, dt: f64, completed: f64, mesh: bool) -> Result<RichardsNodeDraft, StepError> {
    let ur = r as usize;
    let st = &ctx.node_subtypes.storages;
    let node = st.node_idx[ur];
    let un = node as usize;
    let mut s = st.lid_state[ur].clone();
    s.richards_report = report.clone();
    s.quality_old_mobile = ctx.nodes.volume[un];
    s.quality_old_water = s.richards_water.clone();
    s.quality_transfers.clear();
    let capture = if mesh { report.inflow / CUBIC_METRES } else { cmax(0., ctx.nodes.lat_flow[un]) * dt };
    s.captured_flow = capture / dt;
    s.treated_volume += capture;
    if !mesh && capture > 0.0 {
        s.quality_transfers.push(LidWaterTransfer { from: -2, to: -1, volume: capture });
    }
    let lat_flow = ctx.nodes.lat_flow[un] - s.captured_flow;
    let volume = col.water[0] / CUBIC_METRES;
    for i in 1..col.water.len() {
        s.richards_water[i] = col.water[i] / CUBIC_METRES;
    }
    for k in 0..s.richards_report.transfers.len() {
        let t = s.richards_report.transfers[k];
        s.quality_transfers.push(LidWaterTransfer { from: if t.from == 0 { -1 } else { t.from }, to: if t.to == 0 { -1 } else { t.to }, volume: t.volume / CUBIC_METRES });
    }
    s.evap_volume = s.richards_report.evaporation / CUBIC_METRES;
    s.richards_bottom_loss = if mesh { 0. } else { cmax(0., s.richards_report.bottom_volume) / CUBIC_METRES };
    refresh(&mut s, &layers(ctx, st.lid[ur].control), ctx.node_names.name_of(node)).map_err(from_refresh)?;
    s.richards_completed_time = completed;
    s.richards_mesh_pending_finish = mesh;
    s.richards_mesh_routed = false;
    // Richards geometry depends only on the surface cell; accepted porous water does not
    // change these curves, so the current context is read-only.
    let us = ucf::get_unit_system(ctx.options.flow_units as i32);
    let full_volume = node::get_volume(&ctx.nodes, node, ctx.nodes.full_depth[un], None, us, Some(&ctx.node_subtypes));
    let depth = node::get_depth(&ctx.nodes, node, volume, None, us, Some(&ctx.node_subtypes));
    let head = ctx.nodes.invert_elev[un] + depth;
    Ok(RichardsNodeDraft { row: r, node, state: s, volume, depth, head, full_volume, lat_flow })
}

/// lidnode::detail::publish
pub fn publish(ctx: &mut SimulationContext, mut d: RichardsNodeDraft) {
    std::mem::swap(&mut ctx.node_subtypes.storages.lid_state[d.row as usize], &mut d.state);
    let n = d.node as usize;
    ctx.nodes.lat_flow[n] = d.lat_flow;
    ctx.nodes.volume[n] = d.volume;
    ctx.nodes.full_volume[n] = d.full_volume;
    ctx.nodes.rpt_full_volume[n] = d.full_volume;
    ctx.nodes.depth[n] = d.depth;
    ctx.nodes.head[n] = d.head;
}

/// lidnode::detail::routed: a Mesh state with its accepted porous port water.
pub fn routed(ctx: &SimulationContext, r: i32) -> Result<LidNodeState, StepError> {
    let ur = r as usize;
    let st = &ctx.node_subtypes.storages;
    let mut s = st.lid_state[ur].clone();
    if s.richards_water.len() != s.cells.len() || s.port_delta.len() != s.cells.len() {
        return Err(StepError::Exception("Invalid Mesh Richards accepted port layout."));
    }
    for i in 0..s.cells.len() {
        if !s.port_delta[i].is_finite() {
            return Err(StepError::Exception("Nonfinite Mesh Richards accepted port water."));
        }
        if i == 0 && s.port_delta[i] != 0.0 {
            return Err(StepError::Exception("Mesh Richards surface port water belongs to the mobile owner."));
        }
        s.richards_water[i] += s.port_delta[i];
    }
    refresh(&mut s, &layers(ctx, st.lid[ur].control), ctx.node_names.name_of(st.node_idx[ur])).map_err(from_refresh)?;
    Ok(s)
}

/// lidnode::validateRichardsRouting: true when Mesh-aquifer Richards owners accepted this
/// interval (none in a 1D build).
pub fn validate_richards_routing(ctx: &SimulationContext, dt: f64) -> Result<bool, StepError> {
    let st = &ctx.node_subtypes.storages;
    let mut first: Option<usize> = None;
    for r in 0..st.count() as usize {
        let s = &st.lid_state[r];
        if !s.richards || st.lid[r].bottom != LidNodeBottom::MeshAquifer {
            continue;
        }
        let end = match first {
            Some(f) => st.lid_state[f].richards_completed_time,
            None => s.richards_completed_time,
        };
        check_mesh_interval(s, ctx.current_time, end)?;
        if s.richards_mesh_routed {
            return Err(StepError::Exception("Mesh Richards owner interval was already routed."));
        }
        let interval = s.richards_completed_time - s.richards_mesh_start;
        if !dt.is_finite() || dt <= 0.0 || (dt - interval).abs() > 8.0 * f64::EPSILON * max4(1., dt, s.richards_mesh_start, s.richards_completed_time) {
            return Err(StepError::Exception("Mesh Richards global routing duration differs from its accepted owner interval."));
        }
        first = Some(r);
    }
    Ok(first.is_some())
}

/// lidnode::prepareRichardsStep: advance every standalone Richards column over the routing
/// interval and publish the accepted states (all or none).
pub fn prepare_richards_step(ctx: &mut SimulationContext, dt: f64, evaporation: f64, mesh_prepared: bool) -> Result<(), StepError> {
    if mesh_prepared {
        validate_richards_routing(ctx, dt)?;
    }
    let count = ctx.node_subtypes.storages.count() as usize;
    // Refuse the whole batch before capturing forcing or changing any owner.
    for r in 0..count {
        let st = &ctx.node_subtypes.storages;
        if !mesh_prepared && st.lid_state[r].richards && st.lid[r].bottom == LidNodeBottom::MeshAquifer {
            return Err(StepError::Exception("Mesh Richards nodes require joint aquifer owner preparation."));
        }
    }
    let mut columns: Vec<Column> = Vec::with_capacity(count);
    let mut rows: Vec<i32> = Vec::with_capacity(count);
    for r in 0..count {
        let st = &ctx.node_subtypes.storages;
        if st.lid_state[r].richards && st.lid[r].bottom != LidNodeBottom::MeshAquifer {
            rows.push(r as i32);
            columns.push(column(ctx, r as i32, dt, evaporation, false));
        }
    }
    if columns.is_empty() {
        return Ok(());
    }
    let reports = richards::advance(&mut columns, dt);
    for a in 0..rows.len() {
        if !reports[a].ok {
            let st = &ctx.node_subtypes.storages;
            let mut what: Vec<u8> = b"LID Richards at ".to_vec();
            what.extend_from_slice(ctx.node_names.name_of(st.node_idx[rows[a] as usize]));
            what.extend_from_slice(b": ");
            what.extend_from_slice(&reports[a].error);
            return Err(thrown(what));
        }
    }
    let mut drafts: Vec<RichardsNodeDraft> = Vec::with_capacity(rows.len());
    for a in 0..rows.len() {
        drafts.push(draft(ctx, rows[a], &columns[a], &reports[a], dt, ctx.current_time + dt, false)?);
    }
    for d in drafts {
        publish(ctx, d);
    }
    Ok(())
}
