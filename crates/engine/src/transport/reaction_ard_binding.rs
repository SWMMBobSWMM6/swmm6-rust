// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine
// src/engine/transport/components/ReactionModule/ReactionArdBinding.cpp and
// ReactionArdBinding.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The EULERIAN_ARD engine's reaction stage (reactArdStage): exact first-order pollutant decay on
// the mesh cells and node stores, then the reactions component's species per cell (pipe scope)
// and per node store (tank scope) through the shared integrator.
//
// Differences in form, never in result:
//   * Upstream's function-local static ArdScratch (integrator workspace, gather blocks) is a
//     per-call local here. Every gather block entry is written before the integrator reads it,
//     and the integrator writes every workspace entry it reads (the legacy binding's
//     react_elements builds its workspace per call the same way, gated identical).
//   * The cell and store arrays are slices instead of raw pointers.

use crate::core::context::SimulationContext;
use crate::core::engine_step_routing::StepError;
use crate::data::reaction_tokens::RxHydVar;
use crate::transport::reaction_legacy_step::{clamp_rate_species, integrator_step, note_negative_clamp, reaction_area_volume, Ws};

/// ardReactionsActive
pub fn ard_reactions_active(ctx: &SimulationContext) -> bool {
    ctx.reactions.configured && ctx.reactions.compiled
}

/// ardHasWallSpecies
pub fn ard_has_wall_species(ctx: &SimulationContext) -> bool {
    ctx.reactions.species_is_wall.iter().any(|&w| w != 0)
}

/// containFailure: one warning per run for a failed integration.
fn contain_failure(ctx: &mut SimulationContext, what: &str, index: i32, error: &[u8]) {
    if ctx.reactions.warned_react_failure {
        return;
    }
    ctx.reactions.warned_react_failure = true;
    let mut w = format!("Reaction step failed at ARD {} index {} (element state left unchanged): ", what, index).into_bytes();
    w.extend_from_slice(error);
    ctx.warnings.push(w);
}

/// The mesh arrays reactArdStage reads and writes.
pub struct ArdArrays<'a> {
    pub cell_phi: &'a mut [f64],
    pub cell_a: &'a [f64],
    pub cell_dx: &'a [f64],
    pub n_cells: i32,
    pub node_mass: &'a mut [f64],
    pub node_vol: &'a [f64],
    pub n_nodes: i32,
    pub cell_conduit: &'a [i32],
    pub conduit_link: &'a [i32],
}

/// reactArdStage
pub fn react_ard_stage(ctx: &mut SimulationContext, dt: f64, a: ArdArrays, n_pollut: i32, ns_total: i32, min_store_vol: f64, temp_row: i32) -> Result<(), StepError> {
    if dt <= 0.0 {
        return Ok(());
    }
    let ArdArrays { cell_phi, cell_a, cell_dx, n_cells, node_mass, node_vol, n_nodes, cell_conduit, conduit_link } = a;
    let unc = n_cells as usize;
    let uns = ns_total as usize;
    // 1. Pollutant kdecay: exact exponential, booked into qual_routing_reacted
    for p in 0..n_pollut {
        let k = ctx.pollutants.k_decay[p as usize];
        if k == 0.0 {
            continue;
        }
        let f = (-k * dt).exp(); // LIBM
        let loss_frac = 1.0 - f;
        let mut removed = 0.0;
        let row = p as usize * unc;
        for c in 0..n_cells {
            let ucc = c as usize;
            let uc = row + ucc;
            if cell_phi[uc] > 0.0 {
                removed += cell_a[ucc] * cell_dx[ucc] * cell_phi[uc] * loss_frac;
            }
            cell_phi[uc] *= f;
            if cell_phi[uc] < 0.0 {
                cell_phi[uc] = 0.0;
            }
        }
        for nd in 0..n_nodes {
            let ui = nd as usize * uns + p as usize;
            if node_mass[ui] > 0.0 {
                removed += node_mass[ui] * loss_frac;
            }
            node_mass[ui] *= f;
            if node_mass[ui] < 0.0 {
                node_mass[ui] = 0.0;
            }
        }
        if (p as usize) < ctx.mass_balance.qual_routing_reacted.len() {
            ctx.mass_balance.qual_routing_reacted[p as usize] += removed;
        }
    }
    // 2. MSX species: per-cell pipe scope, per-store tank scope.
    if !ard_reactions_active(ctx) {
        return Ok(());
    }
    let nm = ctx.reactions.n_species();
    if nm == 0 || ns_total < n_pollut + nm {
        return Ok(());
    }
    let unm = nm as usize;
    let mut ws = Ws::new(&ctx.reactions);
    let mut coef_scratch: Vec<f64> = Vec::new();
    let mut msx = vec![0.0f64; unm];
    let mut pollut = vec![0.0f64; if n_pollut > 0 { n_pollut as usize } else { 1 }];
    let mut hydvar = [0.0f64; RxHydVar::Count as usize];

    // Cells (pipe scope)
    for c in 0..n_cells {
        let ucell = c as usize;
        for m in 0..unm {
            msx[m] = cell_phi[(n_pollut as usize + m) * unc + ucell];
        }
        for p in 0..n_pollut as usize {
            pollut[p] = cell_phi[p * unc + ucell];
        }
        hydvar[RxHydVar::Hrt as usize] = 0.0;
        let link = if !cell_conduit.is_empty() && !conduit_link.is_empty() { conduit_link[cell_conduit[ucell] as usize] } else { -1 };
        hydvar[RxHydVar::Av as usize] = if ctx.reactions.wall_geometry && link >= 0 { reaction_area_volume(ctx, true, link, ctx.links.volume[link as usize]) } else { 0.0 };
        hydvar[RxHydVar::Temp as usize] = if temp_row >= 0 { cell_phi[temp_row as usize * unc + ucell] } else { ctx.reactions.default_temp_c };
        let rep = {
            let rx = &ctx.reactions;
            let coefs: Vec<f64> = rx.coefficients(true, link, &mut coef_scratch).to_vec();
            let pl: Option<&[f64]> = if n_pollut > 0 { Some(&pollut) } else { None };
            integrator_step(rx, &rx.pipe_form, &rx.pipe_expr, dt, &mut msx, &mut hydvar, &mut ws, pl, None, &coefs)?
        };
        if rep.ok {
            if clamp_rate_species(&ctx.reactions, false, &mut msx, None) > 0.0 {
                note_negative_clamp(ctx, b"a transport cell");
            }
            for m in 0..unm {
                cell_phi[(n_pollut as usize + m) * unc + ucell] = msx[m];
            }
        } else {
            contain_failure(ctx, "cell", c, &rep.error);
        }
    }
    // Node stores (tank scope): integrate concentrations, write back mass.
    for nd in 0..n_nodes {
        let und = nd as usize;
        let vol = node_vol[und];
        if vol <= min_store_vol {
            continue;
        }
        hydvar[RxHydVar::Av as usize] = if ctx.reactions.wall_geometry { reaction_area_volume(ctx, false, nd, vol) } else { 0. };
        for m in 0..unm {
            msx[m] = node_mass[und * uns + n_pollut as usize + m] / vol;
        }
        for p in 0..n_pollut as usize {
            pollut[p] = node_mass[und * uns + p] / vol;
        }
        hydvar[RxHydVar::Hrt as usize] = if und < ctx.nodes.hrt.len() { ctx.nodes.hrt[und] } else { 0.0 };
        hydvar[RxHydVar::Temp as usize] = if temp_row >= 0 { node_mass[und * uns + temp_row as usize] / vol } else { ctx.reactions.default_temp_c };
        let rep = {
            let rx = &ctx.reactions;
            let coefs: Vec<f64> = rx.coefficients(false, nd, &mut coef_scratch).to_vec();
            let pl: Option<&[f64]> = if n_pollut > 0 { Some(&pollut) } else { None };
            integrator_step(rx, &rx.tank_form, &rx.tank_expr, dt, &mut msx, &mut hydvar, &mut ws, pl, None, &coefs)?
        };
        if rep.ok {
            if clamp_rate_species(&ctx.reactions, true, &mut msx, None) > 0.0 {
                note_negative_clamp(ctx, b"a node store");
            }
            for m in 0..unm {
                node_mass[und * uns + n_pollut as usize + m] = msx[m] * vol;
            }
        } else {
            contain_failure(ctx, "node store", nd, &rep.error);
        }
    }
    Ok(())
}
