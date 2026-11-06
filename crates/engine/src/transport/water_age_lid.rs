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
// The run-time step functions setLidInflowAge and routeLidLayerAge were added on 6 November
// 2026.

use crate::core::context::SimulationContext;
use crate::data::lid_layer_species::{LidLayer, LidSpecies};
use crate::data::water_age::WaterAgeSource;
use crate::core::cmath::cmax;
use crate::hydrology::lid::{LIDSolver, LIDType};
use crate::transport::lid_layer_common::{build_offsets, donors_for, layer_volumes, K_ABSENT, K_EXTERNAL, K_NL};

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

/// setLidInflowAge: the age of a unit's inflow (rain at the RAINFALL age, the captured surface
/// water at the subcatchment's runoff age), weighted by the four inflow rates.
#[allow(clippy::too_many_arguments)]
pub fn set_lid_inflow_age(ctx: &mut SimulationContext, type_index: i32, unit: i32, subcatch: i32, rain_rate: f64, q_imperv: f64, q_perv: f64, q_runon: f64, lid_area: f64) {
    if !ctx.options.water_age || !ctx.lid_layer_state.active() {
        return;
    }
    let st = &ctx.lid_layer_state;
    let ut = type_index as usize;
    if ut + 1 >= st.group_offset.len() {
        return;
    }
    let flat = st.group_offset[ut] + unit;
    if flat < 0 || flat >= st.n_units {
        return;
    }
    let q_rain = rain_rate * lid_area;
    let a_rain = ctx.water_age_config.global_age[WaterAgeSource::Rainfall as usize];
    let ws = &ctx.water_age_state;
    let sc = subcatch;
    let a_surface = if sc >= 0 && (sc as usize) < ws.subcatch_runoff_age.len() && ws.subcatch_runoff_age[sc as usize] > 0.0 { ws.subcatch_runoff_age[sc as usize] } else { a_rain };
    let q_surface = q_imperv + q_perv + q_runon;
    let q_tot = q_rain + q_surface;
    let age = if q_tot > 0.0 { (q_rain * a_rain + q_surface * a_surface) / q_tot } else { a_rain };
    let k = flat as usize * st.n_species as usize + K_AGE as usize;
    ctx.lid_layer_state.inflow_value[k] = age;
}

/// routeLidLayerAge: age and complete-mix each layer on the published inflows (donors read
/// from the old state); the drain leaves at the storage (roof disconnection: surface) age.
pub fn route_lid_layer_age(ctx: &mut SimulationContext, solver: &LIDSolver, dt: f64) {
    if !ctx.options.water_age || !ctx.lid_layer_state.active() || !(dt > 0.0) {
        return;
    }
    let st = &mut ctx.lid_layer_state;
    for t in 0..solver.num_groups() {
        let g = solver.group(t);
        if g.count == 0 {
            continue;
        }
        let ut = t as usize;
        if ut + 1 >= st.group_offset.len() {
            continue;
        }
        for u in 0..g.count {
            let ui = u as usize;
            let flat = st.group_offset[ut] + u;
            if flat < 0 || flat >= st.n_units {
                continue;
            }
            let uf = flat as usize;
            let d = donors_for(g.type_, g.soil_thick[ui] > 0.0);
            let donor = [d.surface, d.pavement, d.soil, d.storage];
            let in_rate = [g.in_surf[ui], g.in_pave[ui], g.in_soil[ui], g.in_stor[ui]];
            let mut v_new = [0.0f64; K_NL];
            layer_volumes(g, ui, &mut v_new);
            let a_ext = st.inflow_value[uf * st.n_species as usize + K_AGE as usize];
            let mut a_old = [0.0f64; K_NL];
            for k in 0..K_NL {
                a_old[k] = st.value[st.layer_index(flat, lid_layer_of(k), K_AGE)];
            }
            for k in 0..K_NL {
                let vidx = uf * K_NL + k;
                let v_old = st.vol_prev[vidx];
                let v_in = cmax(0.0, in_rate[k]) * dt;
                if donor[k] == K_ABSENT || (v_old <= K_TINY_VOL && v_in <= K_TINY_VOL) {
                    let li = st.layer_index(flat, lid_layer_of(k), K_AGE);
                    st.value[li] = 0.0;
                    st.vol_prev[vidx] = v_new[k];
                    continue;
                }
                let mut a = a_old[k] + dt;
                if v_in > K_TINY_VOL {
                    let a_in = if donor[k] == K_EXTERNAL {
                        a_ext
                    } else if donor[k] >= 0 {
                        a_old[donor[k] as usize] + dt
                    } else {
                        a_ext
                    };
                    let den = v_old + v_in;
                    a = if den > K_TINY_VOL { (a * v_old + a_in * v_in) / den } else { a_in };
                }
                let li = st.layer_index(flat, lid_layer_of(k), K_AGE);
                st.value[li] = cmax(a, 0.0);
                st.vol_prev[vidx] = v_new[k];
            }
            let from_surface = g.type_ == LIDType::RoofDiscon;
            let src = if from_surface { LidLayer::Surface } else { LidLayer::Storage };
            let di = uf * st.n_species as usize + K_AGE as usize;
            st.drain_value[di] = st.value[st.layer_index(flat, src, K_AGE)];
        }
    }
}
