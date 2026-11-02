// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/components/HeatModule/HeatLid.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Phase H5b, temperature through the LID layer stack: initLidLayerTemperature, which
// SWMMEngine::initHydrology calls right after LIDSolver::init.
//
// NOT YET TRANSLATED (run-time step functions): setLidInflowTemperature,
// routeLidLayerTemperature and HeatLid.cpp's file-static helpers.

use crate::core::context::SimulationContext;
use crate::data::heat::HeatSource;
use crate::data::lid_layer_species::{LidLayer, LidSpecies};
use crate::hydrology::lid::LIDSolver;
use crate::transport::lid_layer_common::{build_offsets, layer_volumes, K_NL};

const K_TEMP: i32 = LidSpecies::Temperature as i32;

/// `static_cast<LidLayer>(k)` for a layer index 0..3.
fn lid_layer_of(k: usize) -> LidLayer {
    match k {
        0 => LidLayer::Surface,
        1 => LidLayer::Pavement,
        2 => LidLayer::Soil,
        3 => LidLayer::Storage,
        _ => LidLayer::Count,
    }
}

/// Size the shared layer scalar state and seed every conventional-LID layer temperature
/// (wet or dry) with the INITIAL_STATE temperature.
pub fn init_lid_layer_temperature(ctx: &mut SimulationContext, solver: &LIDSolver) {
    if !ctx.options.heat_transport {
        return;
    }
    let offsets = build_offsets(solver);
    let units = if offsets.is_empty() { 0 } else { *offsets.last().unwrap() };
    if units <= 0 {
        return;
    }

    let t_init = ctx.heat_config.global_temp[HeatSource::InitialState as usize];
    let st = &mut ctx.lid_layer_state;
    st.ensure_sized(units, LidSpecies::Count as i32, &offsets);

    for t in 0..solver.num_groups() {
        let g = solver.group(t);
        for u in 0..g.count {
            let ui = u as usize;
            let flat = offsets[t as usize] + u;
            let mut v = [0.0f64; K_NL];
            layer_volumes(g, ui, &mut v);
            for k in 0..K_NL {
                let idx = (flat as usize) * K_NL + k;
                st.vol_prev[idx] = v[k];
                let li = st.layer_index(flat, lid_layer_of(k), K_TEMP);
                st.value[li] = t_init;
            }
        }
    }
}
