// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/components/WaterAgeModule/WaterAgeLid.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Phase A4, water age through the LID layer stack: initLidLayerAge, which
// SWMMEngine::initHydrology calls right after LIDSolver::init.
//
// NOT YET TRANSLATED (run-time step functions): setLidInflowAge, routeLidLayerAge and
// the rest of WaterAgeLid.cpp.

use crate::core::context::SimulationContext;
use crate::data::lid_layer_species::{LidLayer, LidSpecies};
use crate::data::water_age::WaterAgeSource;
use crate::hydrology::lid::LIDSolver;
use crate::transport::lid_layer_common::{build_offsets, layer_volumes, K_NL};

const K_AGE: i32 = LidSpecies::Age as i32;

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

/// ft of water per unit area
const K_TINY_VOL: f64 = 1.0e-12;

/// Size the shared conventional-LID scalar storage and seed retained-water layer ages.
pub fn init_lid_layer_age(ctx: &mut SimulationContext, solver: &LIDSolver) {
    if !ctx.options.water_age {
        return;
    }
    let offsets = build_offsets(solver);
    let units = if offsets.is_empty() { 0 } else { *offsets.last().unwrap() };
    if units <= 0 {
        return;
    }

    let a_init = ctx.water_age_config.global_age[WaterAgeSource::InitialState as usize];
    let st = &mut ctx.lid_layer_state;
    // H5b: ensureSized, not resize (the block is shared with the heat track).
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
                if v[k] > K_TINY_VOL {
                    let li = st.layer_index(flat, lid_layer_of(k), K_AGE);
                    st.value[li] = a_init;
                }
            }
        }
    }
}
