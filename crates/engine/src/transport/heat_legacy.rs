// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/components/HeatModule/HeatLegacy.cpp,
// src/engine/transport/InitialQualitySeeds.hpp (applyInitialTempOverrides) and the guard of
// src/engine/transport/components/HeatFluxModules/HeatFluxes.cpp (applyHeatFluxes) and
// BedExchange.cpp (bedExchangeEnabled) (HydroCouple/openswmm.engine, branch swmm6_rel,
// commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The LEGACY temperature mirror (HEAT_TRANSPORT ON under QUALITY_SOLVER LEGACY):
// routeLegacyHeat, run in QualitySolver::execute after the age mirror.
//
// Differences in form, never in result:
//   * Upstream's process-wide static HeatScratch (old-state snapshots, the temperature-volume
//     accumulator) is written in full before it is read on every call; here it is a local.
//
// NOT YET TRANSLATED: the heat flux modules (heat::applyHeatFluxes past its guard: surface
// and radiative exchange, the bed, solar forcing). With any of them on, applyHeatFluxes
// returns StepError::NotTranslated; SWMMEngine::stepRouting refuses such a deck first.

use crate::core::cmath::{cclamp, cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::engine_step_routing::StepError;
use crate::core::options::RoutingModel;
use crate::data::heat::HeatSource;
use crate::data::initial_quality::InitialQualityData;
use crate::quality::quality_routing::{conduit_mixing_inflow, link_is_dry, link_takes_upstream_value, node_is_reactor, LEGACY_ZERO};

/// heat::bedExchangeEnabled
pub fn bed_exchange_enabled(ctx: &SimulationContext) -> bool {
    ctx.options.heat_transport && ctx.heat_config.sediment_exchange
}

/// Whether heat::applyHeatFluxes would act past its guard (not translated yet).
pub fn heat_fluxes_active(ctx: &SimulationContext) -> bool {
    ctx.options.heat_transport && (ctx.heat_config.surface_exchange || ctx.heat_config.radiative_exchange || bed_exchange_enabled(ctx))
}

/// heat::applyHeatFluxes: its guard; the flux modules themselves are not translated.
fn apply_heat_fluxes(ctx: &mut SimulationContext, dt: f64) -> Result<(), StepError> {
    if !ctx.options.heat_transport {
        return Ok(());
    }
    let bed_on = bed_exchange_enabled(ctx);
    if !ctx.heat_config.surface_exchange && !ctx.heat_config.radiative_exchange && !bed_on {
        return Ok(());
    }
    if !(dt > 0.0) {
        return Ok(());
    }
    Err(StepError::NotTranslated("heat fluxes"))
}

/// transport::applyInitialTempOverrides: [INITIAL_QUALITY] __TEMPERATURE__ rows (degC).
pub fn apply_initial_temp_overrides(ctx: &mut SimulationContext) {
    let hs = &mut ctx.heat_state;
    let iq = &ctx.initial_quality;
    for r in 0..iq.count() {
        let ur = r as usize;
        if iq.kind[ur] != InitialQualityData::K_KIND_TEMPERATURE {
            continue;
        }
        let ei = iq.elem_idx[ur];
        if ei < 0 {
            continue;
        }
        let arr = if iq.is_link[ur] != 0 { &mut hs.link_temp } else { &mut hs.node_temp };
        if (ei as usize) < arr.len() {
            arr[ei as usize] = iq.value[ur];
        }
    }
}

/// transport::routeLegacyHeat: the flux modules, then advective mixing of temperature on the
/// quality path's kernels (two-sided clamps; a dry node or link keeps or takes a temperature
/// rather than zero).
pub fn route_legacy_heat(ctx: &mut SimulationContext, dt: f64) -> Result<(), StepError> {
    if !ctx.options.heat_transport || dt <= 0.0 {
        return Ok(());
    }
    let nn = ctx.n_nodes();
    let nl = ctx.n_links();
    let t_init = ctx.heat_config.global_temp[HeatSource::InitialState as usize];
    if ctx.heat_state.node_temp.len() != nn as usize {
        ctx.heat_state.resize(nn, nl, t_init);
    }

    let mut node_old = vec![0.0f64; nn as usize];
    let mut link_old = vec![0.0f64; nl as usize];
    let mut temp_in = vec![0.0f64; nn as usize];

    // ---- 0. INITIAL_STATE seeding, once
    if !ctx.heat_state.legacy_seeded {
        let hs = &mut ctx.heat_state;
        hs.node_temp.iter_mut().for_each(|x| *x = t_init);
        hs.link_temp.iter_mut().for_each(|x| *x = t_init);
        apply_initial_temp_overrides(ctx);
        ctx.heat_state.legacy_seeded = true;
    }

    // ---- 1. Source terms, then old state
    apply_heat_fluxes(ctx, dt)?;

    for i in 0..nn as usize {
        node_old[i] = ctx.heat_state.node_temp[i];
    }
    for j in 0..nl as usize {
        link_old[j] = ctx.heat_state.link_temp[j];
    }

    // ---- 2. Temperature-volume accumulation + the loaders' per-source rates
    temp_in.iter_mut().for_each(|x| *x = 0.0);
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
        temp_in[downstream as usize] += q * link_old[uj];
    }
    for i in 0..nn as usize {
        if i < ctx.heat_state.node_temp_vol_in.len() {
            temp_in[i] += ctx.heat_state.node_temp_vol_in[i];
        }
    }

    // ---- 3. Node mixing
    for i in 0..nn {
        let ui = i as usize;
        let v_old = ctx.nodes.old_volume[ui];
        let v_in = ctx.nodes.qual_vol_in[ui];
        let q_in = if dt > 0.0 { v_in / dt } else { 0.0 };
        let t_old = node_old[ui];

        if !node_is_reactor(ctx, i) {
            ctx.heat_state.node_temp[ui] = if q_in > LEGACY_ZERO { (temp_in[ui] * dt) / v_in } else { t_old };
            continue;
        }

        if q_in <= LEGACY_ZERO {
            ctx.heat_state.node_temp[ui] = t_old;
            continue;
        }
        let mass_in = temp_in[ui] * dt;
        let t_in = mass_in / v_in;
        let mut t_new = (t_old * v_old + mass_in) / (v_old + v_in);
        t_new = cclamp(t_new, cmin(t_old, t_in), cmax(t_old, t_in));
        ctx.heat_state.node_temp[ui] = t_new;
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
        let t_up = ctx.heat_state.node_temp[upstream as usize];

        if link_takes_upstream_value(ctx, j) {
            ctx.heat_state.link_temp[uj] = t_up;
            continue;
        }

        let v_old = ctx.links.old_volume[uj];
        let t_old = link_old[uj];

        let t_new;
        if is_steady {
            t_new = t_up;
        } else if link_is_dry(ctx, j) {
            t_new = t_up;
        } else {
            let q_in = conduit_mixing_inflow(ctx, j, dt);
            if q_in <= LEGACY_ZERO {
                t_new = t_old;
            } else {
                let v_in = q_in * dt;
                let t = (t_old * v_old + t_up * v_in) / (v_old + v_in);
                t_new = cclamp(t, cmin(t_old, t_up), cmax(t_old, t_up));
            }
        }
        ctx.heat_state.link_temp[uj] = t_new;
    }
    Ok(())
}
