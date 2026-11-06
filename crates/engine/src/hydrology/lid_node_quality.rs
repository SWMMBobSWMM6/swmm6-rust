// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/LidNodeTreatment.cpp and LidNode.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The node-LID (storage LID columns, SWMM6-only) pollutant calls the LEGACY quality router
// makes each routing step: prepareQuality (retained-cell transfers, layer treatment and decay),
// prepareOutletQuality (treated mobile outlets), receiveQuality, outletQuality, and the
// file-local layerVolume and treatedConcentration. (heldMass is in lid_node_held.rs;
// validTreatment and readTreatment, the input side, are in lid_node.rs.)
//
// Differences in form, never in result:
//   * A storage row's LidNodeState is moved out of the context while its row is prepared
//     (std::mem::take) and put back after, so the row can be changed while the rest of the
//     context is read; nothing these functions call reads lid_state.
//   * treatment::parse can throw (std::stod) inside treatedConcentration. These functions are
//     not noexcept upstream, so the exception reaches step()'s catch: here
//     Err(StepError::Exception("stod")), with the row's state put back first.
//   * The recursive co-treatment lambda is a function taking the removal scratch.

use crate::core::cmath::{cclamp, cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::engine_step_routing::StepError;
use crate::core::options::SimulationOptions;
use crate::core::unit_conversion as ucf;
use crate::data::lid_node::{LidNodeCell, LidNodeLayer, LidNodeLayerKind, LidNodeState};
use crate::data::name_index::NameIndex;
use crate::data::treatment_types::{TreatExpr, TreatVar};
use crate::hydrology::lid_node::layers;
use crate::hydrology::lid_node::treatment as tparse;
use crate::hydrology::lid_node_step::{port_cell, port_offset};
use crate::quality::treatment;

/// layerVolume: the water of every cell in `cell`'s authored layer, ft3.
fn layer_volume(state: &LidNodeState, water: &[f64], cell: i32) -> f64 {
    let mut total = 0.0;
    for i in 0..water.len() {
        if state.cells[i].layer == state.cells[cell as usize].layer {
            total += water[i];
        }
    }
    total
}

/// What treatedConcentration reads from the context.
struct Tctx<'a> {
    names: &'a NameIndex,
    opts: &'a SimulationOptions,
    np: i32,
}

/// treatedConcentration's recursive `evaluate` lambda: the removal of `pollutant` at this layer.
#[allow(clippy::too_many_arguments)]
fn eval_removal(t: &Tctx, cell: &LidNodeCell, layer: &LidNodeLayer, pollutant: i32, concentrations: &[f64], volume: f64, flow: f64, dt: f64, removals: &mut Vec<f64>) -> Result<f64, StepError> {
    let np = t.np;
    let up = pollutant as usize;
    if removals[up] >= 0.0 {
        return Ok(if removals[up] <= 1.0 { removals[up] } else { 0.0 });
    }
    removals[up] = 2.0; // break cyclic co-treatment references
    let it = layer.treatment.iter().find(|rule| t.names.find(&rule.pollutant) == pollutant);
    let rule = match it {
        None => {
            removals[up] = 0.0;
            return Ok(0.0);
        }
        Some(r) => r,
    };
    let mut r = rule.removal;
    if !rule.expression.is_empty() {
        let mut expr = TreatExpr::default();
        let rc = tparse::parse_names(&rule.expression, &mut expr, t.names.names()).map_err(|_| StepError::Exception("stod"))?;
        if rc == 0 {
            for token in expr.tokens.iter() {
                if token.var == TreatVar::RPollut && token.pollut_ref >= 0 && token.pollut_ref < np {
                    eval_removal(t, cell, layer, token.pollut_ref, concentrations, volume, flow, dt, removals)?;
                }
            }
            let mut safe = removals.clone();
            for value in safe.iter_mut() {
                if *value < 0.0 || *value > 1.0 {
                    *value = 0.0;
                }
            }
            let c = concentrations[up] * (1.0 - r);
            let length = ucf::ucf(ucf::LENGTH, t.opts);
            let result = treatment::evaluate(
                &expr,
                c,
                dt,
                if flow > 0.0 { volume / flow / 3600.0 } else { dt / 3600.0 },
                flow * ucf::ucf(ucf::FLOW, t.opts),
                volume,
                layer.params[0] / ucf::ucf(ucf::RAINDEPTH, t.opts) * length,
                Some(concentrations),
                Some(&safe),
                np,
                cell.area * length * length,
                Some(concentrations),
            );
            if result.is_finite() {
                let out = if expr.is_removal { c * (1.0 - cclamp(result, 0.0, 1.0)) } else { cclamp(result, 0.0, c) };
                r = if concentrations[up] > 0.0 { 1.0 - out / concentrations[up] } else { 0.0 };
            }
        }
    }
    removals[up] = r;
    Ok(r)
}

/// treatedConcentration: a layer outlet's fixed removal and optional R=/C= expression applied
/// to `concentration` (co-treatment resolved recursively, cycles broken as zero).
#[allow(clippy::too_many_arguments)]
fn treated_concentration(t: &Tctx, cell: &LidNodeCell, layer: &LidNodeLayer, p: i32, concentration: f64, concentrations: &[f64], volume: f64, flow: f64, dt: f64) -> Result<f64, StepError> {
    let mut removals = vec![-1.0f64; t.np as usize];
    let r = eval_removal(t, cell, layer, p, concentrations, volume, flow, dt, &mut removals)?;
    Ok(concentration * (1.0 - r))
}

/// lidnode::prepareQuality: move retained-cell pollutant mass through this step's accepted
/// water transfers, treat what leaves a layer, debit retained outlets, decay resident mass and
/// publish the mobile mixture as the node's old concentration.
pub fn prepare_quality(ctx: &mut SimulationContext, dt: f64) -> Result<(), StepError> {
    let np = ctx.n_pollutants();
    if np == 0 || dt <= 0.0 {
        return Ok(());
    }
    for r in 0..ctx.node_subtypes.storages.count() {
        let ur = r as usize;
        {
            let state = &ctx.node_subtypes.storages.lid_state[ur];
            if state.cells.is_empty() || state.quality_old_water.len() != state.cells.len() {
                continue;
            }
        }
        let mut state = std::mem::take(&mut ctx.node_subtypes.storages.lid_state[ur]);
        let res = prepare_row(ctx, &mut state, ur, dt);
        ctx.node_subtypes.storages.lid_state[ur] = state;
        res?;
    }
    Ok(())
}

fn prepare_row(ctx: &mut SimulationContext, state: &mut LidNodeState, r: usize, dt: f64) -> Result<(), StepError> {
    let np = ctx.n_pollutants();
    let unp = np as usize;
    let node = ctx.node_subtypes.storages.node_idx[r];
    let un = node as usize;
    let base = (node * np) as usize;
    let mut water = state.quality_old_water.clone();
    let mut mobile = state.quality_old_mobile;
    let mut mobile_mass = vec![0.0f64; unp];
    let mut lateral_c = vec![0.0f64; unp];
    if state.quality_mass.len() != water.len() * unp {
        state.quality_mass.resize(water.len() * unp, 0.0);
        for i in 0..water.len() {
            for p in 0..unp {
                state.quality_mass[i * unp + p] = water[i] * ctx.nodes.conc_old[base + p];
            }
        }
    }
    for p in 0..unp {
        mobile_mass[p] = ctx.nodes.conc_old[base + p] * mobile;
        if state.richards {
            mobile_mass[p] += state.quality_mass[p]; // dry pond solute redissolves on supply
            state.quality_mass[p] = 0.0;
        }
        lateral_c[p] = if ctx.nodes.qual_vol_in[un] > 0.0 { cmax(0.0, ctx.nodes.qual_mass_in[base + p] * dt / ctx.nodes.qual_vol_in[un]) } else { 0.0 };
    }
    let stack = layers(ctx, ctx.node_subtypes.storages.lid[r].control);
    let transfers = state.quality_transfers.clone();
    for transfer in transfers.iter() {
        let from = transfer.from;
        let to = transfer.to;
        let v = transfer.volume;
        if to == -3 {
            if from < 0 {
                mobile = cmax(0., mobile - v);
            } else {
                water[from as usize] = cmax(0., water[from as usize] - v);
            }
            continue; // evaporation leaves solute
        }
        let available = if from == -2 {
            v
        } else if from < 0 {
            mobile
        } else {
            water[from as usize]
        };
        let mut c = vec![0.0f64; unp];
        for p in 0..unp {
            c[p] = if from == -2 {
                lateral_c[p]
            } else if available > 0.0 {
                (if from < 0 { mobile_mass[p] } else { state.quality_mass[from as usize * unp + p] }) / available
            } else {
                0.0
            };
        }
        let exits = from >= 0 && (to < 0 || state.cells[from as usize].layer != state.cells[to as usize].layer);
        for p in 0..unp {
            let resident = if from == -2 {
                cmax(0.0, v * c[p])
            } else {
                cmax(0.0, if from < 0 { mobile_mass[p] } else { state.quality_mass[from as usize * unp + p] })
            };
            let mass = if state.richards { cclamp(cmin(v, available) * c[p], 0.0, resident) } else { cmin(v, available) * c[p] };
            let t = Tctx { names: &ctx.pollutant_names, opts: &ctx.options, np };
            let removal = if state.richards && exits && c[p] > 0.0 {
                let fc = from as usize;
                cclamp(treated_concentration(&t, &state.cells[fc], &stack[(state.cells[fc].layer - 1) as usize], p as i32, c[p], &c, layer_volume(state, &water, from), v / dt, dt)? / c[p], 0.0, 1.0)
            } else {
                1.0
            };
            let out = if state.richards {
                mass * removal
            } else if exits {
                let fc = from as usize;
                cmin(v, available) * treated_concentration(&t, &state.cells[fc], &stack[(state.cells[fc].layer - 1) as usize], p as i32, c[p], &c, layer_volume(state, &water, from), v / dt, dt)?
            } else {
                mass
            };
            if from == -2 {
                let incoming = &mut ctx.nodes.qual_mass_in[base + p];
                let prior = *incoming;
                *incoming -= mass / dt;
                // A fully captured lateral load can cancel to a tiny
                // negative rate. It is roundoff, not an extraction source.
                if *incoming < 0.0 && -*incoming <= 64.0 * f64::EPSILON * cmax(1.0, prior.abs()) {
                    *incoming = 0.0;
                }
            } else if from < 0 {
                mobile_mass[p] -= mass;
            } else {
                state.quality_mass[from as usize * unp + p] -= mass;
            }
            if to == -4 {
                ctx.mass_balance.qual_routing_seep[p] += out;
            } else if to < 0 {
                mobile_mass[p] += out;
            } else {
                state.quality_mass[to as usize * unp + p] += out;
            }
            ctx.mass_balance.qual_routing_reacted[p] += mass - out;
        }
        if from == -2 {
            ctx.nodes.qual_vol_in[un] = cmax(0.0, ctx.nodes.qual_vol_in[un] - v);
        } else if from < 0 {
            mobile = cmax(0.0, mobile - v);
        } else {
            water[from as usize] = cmax(0.0, water[from as usize] - v);
        }
        if to != -4 {
            if to < 0 {
                mobile += v;
            } else {
                water[to as usize] += v;
            }
        }
    }
    state.quality_outlet_conc.clear();
    state.quality_outlet_conc.resize((ctx.n_links() * np) as usize, f64::NAN);
    let ports = state.quality_ports.clone();
    for port in ports.iter() {
        if port.volume < 0.0 {
            let i = port.cell as usize;
            let v = -port.volume;
            let mut c = vec![0.0f64; unp];
            for p in 0..unp {
                c[p] = if water[i] > 0.0 { state.quality_mass[i * unp + p] / water[i] } else { 0.0 };
            }
            for p in 0..unp {
                let t = Tctx { names: &ctx.pollutant_names, opts: &ctx.options, np };
                let out = treated_concentration(&t, &state.cells[i], &stack[(state.cells[i].layer - 1) as usize], p as i32, c[p], &c, layer_volume(state, &water, i as i32), v / dt, dt)?;
                let mass = if state.richards {
                    cclamp(cmin(v, water[i]) * c[p], 0.0, cmax(0.0, state.quality_mass[i * unp + p]))
                } else {
                    cmin(v, water[i]) * c[p]
                };
                let accepted_out = if state.richards {
                    if c[p] > 0.0 {
                        mass * cclamp(out / c[p], 0.0, 1.0)
                    } else {
                        0.0
                    }
                } else {
                    cmin(v, water[i]) * out
                };
                state.quality_mass[i * unp + p] -= mass;
                state.quality_outlet_conc[(port.link * np) as usize + p] = if v > 0.0 { accepted_out / v } else { 0.0 };
                ctx.mass_balance.qual_routing_reacted[p] += mass - accepted_out;
            }
            water[i] = cmax(0.0, water[i] - v);
        }
    }
    for i in 0..water.len() {
        for p in 0..unp {
            let rules = &stack[(state.cells[i].layer - 1) as usize].treatment;
            let rule = rules.iter().find(|value| ctx.pollutant_names.find(&value.pollutant) == p as i32);
            // Layer decay supplements the pollutant's background decay.
            let rate = ctx.pollutants.k_decay[p] + (match rule {
                Some(rr) => rr.decay / 86400.0,
                None => 0.0,
            });
            let mass = &mut state.quality_mass[i * unp + p];
            let after = *mass * (-rate * dt).exp(); // LIBM
            ctx.mass_balance.qual_routing_reacted[p] += *mass - after;
            *mass = after;
        }
    }
    // Saturated water shares the storage-node reactor. Weight each layer's
    // decay by its share of that mobile water, so saturated aggregate and
    // media remain active treatment volumes as the water table rises.
    let mobile_volume = ctx.nodes.old_volume[un];
    let mut mobile_rates = vec![0.0f64; unp];
    for cell in state.cells.iter() {
        if state.richards && cell.kind != LidNodeLayerKind::Surface {
            continue;
        }
        let fraction = cclamp((ctx.nodes.old_depth[un] - cell.bottom) / (cell.top - cell.bottom), 0.0, 1.0);
        let v = cell.geometric_volume * cmax(0.0, cell.porosity - cell.theta) * fraction;
        for rule in stack[(cell.layer - 1) as usize].treatment.iter() {
            let p = ctx.pollutant_names.find(&rule.pollutant);
            if p >= 0 && mobile_volume > 0.0 {
                mobile_rates[p as usize] += rule.decay * cmin(1.0, v / mobile_volume);
            }
        }
    }
    for p in 0..unp {
        let after = mobile_mass[p] * (-mobile_rates[p] * dt / 86400.0).exp(); // LIBM
        ctx.mass_balance.qual_routing_reacted[p] += mobile_mass[p] - after;
        mobile_mass[p] = after;
        if mobile_volume > (if state.richards { 1.0e-12 } else { 0.0 }) {
            ctx.nodes.conc_old[base + p] = after / mobile_volume;
        } else {
            // Preserve residual solute when the last mobile water moves
            // into retained pores or evaporates; it redissolves there.
            state.quality_mass[p] += after;
            ctx.nodes.conc_old[base + p] = 0.0;
        }
    }
    let _ = mobile;
    Ok(())
}

/// lidnode::prepareOutletQuality: refresh the treated concentration of each mobile outlet
/// against the node's current mixture (retained outlets were debited in prepareQuality).
pub fn prepare_outlet_quality(ctx: &mut SimulationContext, dt: f64, book_reaction: bool) -> Result<(), StepError> {
    let np = ctx.n_pollutants();
    if np <= 0 {
        return Ok(());
    }
    for r in 0..ctx.node_subtypes.storages.count() {
        let ur = r as usize;
        if ctx.node_subtypes.storages.lid_state[ur].cells.is_empty() {
            continue;
        }
        let mut state = std::mem::take(&mut ctx.node_subtypes.storages.lid_state[ur]);
        let res = outlet_row(ctx, &mut state, ur, dt, book_reaction);
        ctx.node_subtypes.storages.lid_state[ur] = state;
        res?;
    }
    Ok(())
}

fn outlet_row(ctx: &mut SimulationContext, state: &mut LidNodeState, r: usize, dt: f64, book_reaction: bool) -> Result<(), StepError> {
    let np = ctx.n_pollutants();
    let unp = np as usize;
    let node = ctx.node_subtypes.storages.node_idx[r];
    let un = node as usize;
    let base = (node * np) as usize;
    let stack = layers(ctx, ctx.node_subtypes.storages.lid[r].control);
    let mut mobile_layers = vec![0.0f64; stack.len()];
    for cell in state.cells.iter() {
        let fraction = cclamp((ctx.nodes.old_depth[un] - cell.bottom) / (cell.top - cell.bottom), 0.0, 1.0);
        mobile_layers[(cell.layer - 1) as usize] += cell.geometric_volume * cmax(0.0, cell.porosity - cell.theta) * fraction;
    }
    let mut c = vec![0.0f64; unp];
    for p in 0..unp {
        c[p] = ctx.nodes.conc[base + p];
    }
    for link in 0..ctx.n_links() {
        let ul = link as usize;
        let q = ctx.links.flow[ul];
        let source = if q >= 0.0 { ctx.links.node1[ul] } else { ctx.links.node2[ul] };
        if source != node || q == 0.0 {
            continue;
        }
        if state.quality_ports.iter().any(|port| port.link == link && port.volume < 0.0) {
            continue;
        }
        let offset = port_offset(ctx, link, node);
        let cell_index = port_cell(state, offset);
        if cell_index < 0 {
            continue;
        }
        let cell = &state.cells[cell_index as usize];
        let lay = &stack[(cell.layer - 1) as usize];
        if lay.treatment.is_empty() {
            continue;
        }
        for p in 0..unp {
            let t = Tctx { names: &ctx.pollutant_names, opts: &ctx.options, np };
            let out = treated_concentration(&t, cell, lay, p as i32, c[p], &c, mobile_layers[(cell.layer - 1) as usize], q.abs(), dt)?;
            state.quality_outlet_conc[ul * unp + p] = out;
            if book_reaction {
                ctx.mass_balance.qual_routing_reacted[p] += (c[p] - out) * q.abs() * dt;
            }
        }
    }
    Ok(())
}

/// lidnode::receiveQuality: the share of a link's delivered pollutant mass that enters a
/// retained cell through a port (moved out of the node's mobile load); true when a port took it.
pub fn receive_quality(ctx: &mut SimulationContext, node: i32, link: i32, volume: f64, p: i32, mass: f64, dt: f64) -> bool {
    let r = ctx.node_subtypes.storage_row(node);
    if r < 0 || volume <= 0.0 {
        return false;
    }
    let np = ctx.n_pollutants();
    let state = &mut ctx.node_subtypes.storages.lid_state[r as usize];
    for port in state.quality_ports.iter() {
        if port.link == link && port.volume > 0.0 {
            let held = mass * cclamp(port.volume / volume, 0.0, 1.0);
            let k = port.cell * np + p;
            if state.quality_mass.len() > k as usize {
                state.quality_mass[k as usize] += held * dt;
            }
            ctx.nodes.qual_mass_in[(node * np + p) as usize] -= held;
            return true;
        }
    }
    false
}

/// lidnode::outletQuality: a cached treated outlet concentration, or `fallback`.
pub fn outlet_quality(ctx: &SimulationContext, node: i32, link: i32, p: i32, fallback: f64) -> f64 {
    let r = ctx.node_subtypes.storage_row(node);
    if r < 0 {
        return fallback;
    }
    let values = &ctx.node_subtypes.storages.lid_state[r as usize].quality_outlet_conc;
    let i = link * ctx.n_pollutants() + p;
    if i >= 0 && (i as usize) < values.len() && values[i as usize].is_finite() {
        values[i as usize]
    } else {
        fallback
    }
}
