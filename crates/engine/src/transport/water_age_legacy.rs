// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/components/WaterAgeModule/
// WaterAgeLegacy.cpp, src/engine/transport/InitialQualitySeeds.hpp
// (applyInitialAgeOverrides) and src/engine/quality/NegativeSources.hpp
// (bookNegativeAgeClamp) (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The LEGACY water-age mirror (WATER_AGE ON under QUALITY_SOLVER LEGACY): routeLegacyAge,
// run last in QualitySolver::execute on the fully accumulated qual_vol_in.
//
// Differences in form, never in result:
//   * Upstream keeps its old-state snapshots and the age-mass accumulator in a process-wide
//     static scratch (AgeScratch, resized on a count change). Every element of it is written
//     before it is read on each call, so here it is a local allocated per call.

use crate::core::cmath::{cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::options::RoutingModel;
use crate::data::initial_quality::InitialQualityData;
use crate::data::node::NodeType;
use crate::data::water_age::WaterAgeSource;
use crate::quality::quality_routing::{conduit_mixing_inflow, link_is_dry, link_takes_upstream_value, node_is_dry, node_is_reactor, LEGACY_ZERO, ZERO_DEPTH};

/// transport::applyInitialAgeOverrides: [INITIAL_QUALITY] __WATER_AGE__ rows (hours) onto the
/// node/link age arrays (seconds); a no-op after a hotstart load.
pub fn apply_initial_age_overrides(ctx: &mut SimulationContext) {
    let ws = &mut ctx.water_age_state;
    if ws.hotstart_loaded {
        return;
    }
    let iq = &ctx.initial_quality;
    for r in 0..iq.count() {
        let ur = r as usize;
        if iq.kind[ur] != InitialQualityData::K_KIND_WATER_AGE {
            continue;
        }
        let ei = iq.elem_idx[ur];
        if ei < 0 {
            continue;
        }
        let arr = if iq.is_link[ur] != 0 { &mut ws.link_age } else { &mut ws.node_age };
        if (ei as usize) < arr.len() {
            arr[ei as usize] = iq.value[ur] * 3600.0;
        }
    }
}

/// quality::bookNegativeAgeClamp: count an age extraction clamp (not ledgered).
pub fn book_negative_age_clamp(ctx: &mut SimulationContext, node: i32) {
    let st = &mut ctx.negsrc;
    st.age_clamp_events += 1;
    if !st.first_clamp_recorded {
        st.first_clamp_recorded = true;
        st.first_node = node;
    }
}

/// transport::routeLegacyAge: age the network by dt, then mix it on the quality path's
/// kernels (no evaporation factor, no reaction).
pub fn route_legacy_age(ctx: &mut SimulationContext, dt: f64) {
    if !ctx.options.water_age || dt <= 0.0 {
        return;
    }
    let nn = ctx.n_nodes();
    let nl = ctx.n_links();
    if ctx.water_age_state.node_age.len() != nn as usize {
        let ns = ctx.n_subcatches();
        ctx.water_age_state.resize(nn, nl, ns); // A3: keep watershed rows
    }

    let mut node_old = vec![0.0f64; nn as usize];
    let mut link_old = vec![0.0f64; nl as usize];
    let mut age_in = vec![0.0f64; nn as usize];

    // ---- 0. INITIAL_STATE seeding, once
    if !ctx.water_age_state.legacy_seeded {
        let a0 = ctx.water_age_config.global_age[WaterAgeSource::InitialState as usize];
        if a0 > 0.0 {
            let ws = &mut ctx.water_age_state;
            ws.node_age.iter_mut().for_each(|x| *x = a0);
            ws.link_age.iter_mut().for_each(|x| *x = a0);
        }
        apply_initial_age_overrides(ctx);
        ctx.water_age_state.legacy_seeded = true;
    }

    // ---- 1. Aging: +dt, then the aged values are this step's "old" state
    for i in 0..nn as usize {
        node_old[i] = ctx.water_age_state.node_age[i] + dt;
    }
    for j in 0..nl as usize {
        link_old[j] = ctx.water_age_state.link_age[j] + dt;
    }

    // ---- 2. Age-mass accumulation + the loaders' per-source rates
    age_in.iter_mut().for_each(|x| *x = 0.0);
    for j in 0..nl {
        let uj = j as usize;
        let q = ctx.links.flow[uj].abs();
        if q <= 0.0 {
            continue;
        }
        let downstream = if ctx.links.flow[uj] >= 0.0 { ctx.links.node2[uj] } else { ctx.links.node1[uj] };
        if downstream < 0 || downstream >= nn {
            continue;
        }
        age_in[downstream as usize] += q * link_old[uj];
    }
    for i in 0..nn as usize {
        if i < ctx.water_age_state.node_age_vol_in.len() {
            age_in[i] += ctx.water_age_state.node_age_vol_in[i];
        }
    }

    // ---- 3. Node mixing
    for i in 0..nn {
        let ui = i as usize;
        let v_old = ctx.nodes.old_volume[ui];
        let v_in = ctx.nodes.qual_vol_in[ui];
        let q_in = if dt > 0.0 { v_in / dt } else { 0.0 };
        let a_old = node_old[ui];

        if !node_is_reactor(ctx, i) {
            if q_in > LEGACY_ZERO {
                let mut mass_in = age_in[ui] * dt;
                if mass_in < 0.0 && a_old * v_old + mass_in < 0.0 {
                    book_negative_age_clamp(ctx, i);
                    mass_in = -(a_old * v_old);
                }
                ctx.water_age_state.node_age[ui] = cmax(mass_in / v_in, 0.0);
            } else if ctx.options.outfall_backflow_zero && ctx.nodes.type_[ui] == NodeType::Outfall {
                ctx.water_age_state.node_age[ui] = 0.0;
            } else {
                ctx.water_age_state.node_age[ui] = if ctx.nodes.depth[ui] > ZERO_DEPTH { a_old } else { 0.0 };
            }
            continue;
        }

        let mut mass_in = age_in[ui] * dt;
        if mass_in < 0.0 && a_old * v_old + mass_in < 0.0 {
            book_negative_age_clamp(ctx, i);
            mass_in = -(a_old * v_old);
        }
        let mut a_new;
        if q_in <= LEGACY_ZERO {
            a_new = a_old;
        } else {
            let a_in = mass_in / v_in;
            let a_max = cmax(a_old, a_in);
            a_new = (a_old * v_old + mass_in) / (v_old + v_in);
            a_new = cmin(a_new, a_max);
        }
        if node_is_dry(ctx, i, q_in) {
            a_new = 0.0;
        }
        ctx.water_age_state.node_age[ui] = cmax(a_new, 0.0);
    }

    // ---- 4. Link update
    let is_steady = ctx.options.routing_model == RoutingModel::Steady;
    for j in 0..nl {
        let uj = j as usize;
        let upstream = if is_steady {
            ctx.links.node1[uj]
        } else if ctx.links.flow[uj] >= 0.0 {
            ctx.links.node1[uj]
        } else {
            ctx.links.node2[uj]
        };
        if upstream < 0 || upstream >= nn {
            continue;
        }
        let a_up = ctx.water_age_state.node_age[upstream as usize];

        if link_takes_upstream_value(ctx, j) {
            ctx.water_age_state.link_age[uj] = cmax(a_up, 0.0);
            continue;
        }

        let v_old = ctx.links.old_volume[uj];
        let a_old = link_old[uj];

        let mut a_new;
        if is_steady {
            a_new = a_up;
        } else {
            let q_in = conduit_mixing_inflow(ctx, j, dt);
            if q_in <= LEGACY_ZERO {
                a_new = a_old;
            } else {
                let v_in = q_in * dt;
                let a_max = cmax(a_old, a_up);
                a_new = (a_old * v_old + a_up * v_in) / (v_old + v_in);
                a_new = cmin(a_new, a_max);
            }
            if link_is_dry(ctx, j) {
                a_new = 0.0;
            }
        }
        ctx.water_age_state.link_age[uj] = cmax(a_new, 0.0);
    }
}
