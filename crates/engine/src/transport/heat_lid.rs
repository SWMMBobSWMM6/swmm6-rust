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
// The run-time step functions setLidInflowTemperature and routeLidLayerTemperature and
// HeatLid.cpp's file-static helpers were added on 6 November 2026. lidLayerHeatCapacity (the
// C API's LID heat ledger) is not translated.

use crate::core::cmath::{cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::engine_step_routing::StepError;
use crate::data::heat::{ConductionConfig, DryTempPolicy, HeatConfigData, HeatSource};
use crate::data::heat_override::{HeatElemKind, HeatElement};
use crate::data::lid_layer_species::{LidLayer, LidSpecies};
use crate::hydrology::lid::{LIDGroupSoA, LIDSolver, LIDType};
use crate::transport::heat_fluxes::{air_temp_celsius, net_flux_out, update_solar_forcing, K_PROBE_C};
use crate::transport::lid_layer_common::{build_offsets, donors_for, layer_volumes, K_ABSENT, K_EXTERNAL, K_NL, K_PAVE, K_SOIL, K_STOR, K_SURF};

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

const K_TINY_VOL: f64 = 1.0e-12;
const K_FT_TO_M: f64 = 0.3048;

fn dry_temperature(ctx: &SimulationContext, held: f64, t_air_c: f64) -> f64 {
    match ctx.heat_config.dry_temp_policy {
        DryTempPolicy::Air => t_air_c,
        DryTempPolicy::Default => HeatConfigData::K_DEFAULT_TEMP,
        DryTempPolicy::Hold => held,
    }
}

/// layerThickness: each layer's thickness, m (surface: its storage depth).
fn layer_thickness(g: &LIDGroupSoA, ui: usize, d: &mut [f64; K_NL]) {
    d[K_SURF as usize] = g.surf_store[ui] * K_FT_TO_M;
    d[K_PAVE as usize] = g.pave_thick[ui] * K_FT_TO_M;
    d[K_SOIL as usize] = g.soil_thick[ui] * K_FT_TO_M;
    d[K_STOR as usize] = g.stor_thick[ui] * K_FT_TO_M;
}

/// heatCapacity: volumetric heat capacity of a layer (water plus solid matrix), J/m3/K.
fn heat_capacity(c: &ConductionConfig, layer: i32, water_m: f64, thick_m: f64, rho_w: f64, cp_w: f64) -> f64 {
    let water = rho_w * cp_w * water_m;
    if layer == K_SURF {
        return if thick_m > 0.0 { water / thick_m } else { 0.0 };
    }
    let solid_m = if thick_m > water_m { thick_m - water_m } else { 0.0 };
    let total = water + c.sed_density * c.sed_specific_heat * solid_m;
    if thick_m > 0.0 {
        total / thick_m
    } else {
        0.0
    }
}

/// conductivity: water/solid blend by the layer's water fraction, W/m/K.
fn conductivity(c: &ConductionConfig, layer: i32, water_m: f64, thick_m: f64) -> f64 {
    if layer == K_SURF {
        return c.water_conductivity;
    }
    if !(thick_m > 0.0) {
        return c.water_conductivity;
    }
    let theta = cmin(1.0, water_m / thick_m);
    theta * c.water_conductivity + (1.0 - theta) * c.sed_conductivity
}

/// solveTridiagonal (Thomas); a zero pivot returns the right-hand side unchanged.
fn solve_tridiagonal(n: usize, a: &[f64; K_NL], b: &[f64; K_NL], c: &[f64; K_NL], d: &[f64; K_NL], x: &mut [f64; K_NL]) {
    let mut cp = [0.0f64; K_NL];
    let mut dp = [0.0f64; K_NL];
    if !(b[0].abs() > 0.0) {
        x[..n].copy_from_slice(&d[..n]);
        return;
    }
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];
    for i in 1..n {
        let m = b[i] - a[i] * cp[i - 1];
        if !(m.abs() > 0.0) {
            x[..n].copy_from_slice(&d[..n]);
            return;
        }
        cp[i] = c[i] / m;
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m;
    }
    x[n - 1] = dp[n - 1];
    let mut i = n as i32 - 2;
    while i >= 0 {
        let ui = i as usize;
        x[ui] = dp[ui] - cp[ui] * x[ui + 1];
        i -= 1;
    }
}

/// setLidInflowTemperature: the temperature of a unit's inflow (rain at the RAINFALL
/// temperature, captured surface water at the subcatchment's runoff temperature).
#[allow(clippy::too_many_arguments)]
pub fn set_lid_inflow_temperature(ctx: &mut SimulationContext, type_index: i32, unit: i32, subcatch: i32, rain_rate: f64, q_imperv: f64, q_perv: f64, q_runon: f64, lid_area: f64) {
    if !ctx.options.heat_transport || !ctx.lid_layer_state.active() {
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
    let t_rain = ctx.heat_config.global_temp[HeatSource::Rainfall as usize];
    let hs = &ctx.heat_state;
    let t_surface = if subcatch >= 0 && (subcatch as usize) < hs.subcatch_runoff_temp.len() { hs.subcatch_runoff_temp[subcatch as usize] } else { t_rain };
    let q_surface = q_imperv + q_perv + q_runon;
    let q_tot = q_rain + q_surface;
    let temp = if q_tot > 0.0 { (q_rain * t_rain + q_surface * t_surface) / q_tot } else { t_rain };
    let k = flat as usize * st.n_species as usize + K_TEMP as usize;
    ctx.lid_layer_state.inflow_value[k] = temp;
}

/// routeLidLayerTemperature: advection (complete mix on the published inflows), then one
/// coupled backward-Euler thermal solve over the column (conduction and the atmospheric flux on
/// the topmost layer); the drain leaves at the storage (roof disconnection: surface)
/// temperature.
pub fn route_lid_layer_temperature(ctx: &mut SimulationContext, solver: &LIDSolver, dt: f64) -> Result<(), StepError> {
    if !ctx.options.heat_transport || !ctx.lid_layer_state.active() || !(dt > 0.0) {
        return Ok(());
    }
    let t_air = air_temp_celsius(ctx);
    let rho_w = ctx.options.water_density;
    let cp_w = ctx.options.water_specific_heat;
    let cc = ctx.heat_config.conduction.clone();
    let do_cond = ctx.heat_config.layer_conduction;
    let do_flux = ctx.heat_config.surface_exchange || ctx.heat_config.radiative_exchange;
    update_solar_forcing(ctx, dt)?;
    for t in 0..solver.num_groups() {
        let g = solver.group(t);
        if g.count == 0 {
            continue;
        }
        let ut = t as usize;
        if ut + 1 >= ctx.lid_layer_state.group_offset.len() {
            continue;
        }
        for u in 0..g.count {
            let ui = u as usize;
            let flat = ctx.lid_layer_state.group_offset[ut] + u;
            if flat < 0 || flat >= ctx.lid_layer_state.n_units {
                continue;
            }
            let uf = flat as usize;
            let d = donors_for(g.type_, g.soil_thick[ui] > 0.0);
            let donor = [d.surface, d.pavement, d.soil, d.storage];
            let in_rate = [g.in_surf[ui], g.in_pave[ui], g.in_soil[ui], g.in_stor[ui]];
            let mut v_new = [0.0f64; K_NL];
            let mut thick = [0.0f64; K_NL];
            layer_volumes(g, ui, &mut v_new);
            layer_thickness(g, ui, &mut thick);
            let nsp = ctx.lid_layer_state.n_species as usize;
            let t_ext = ctx.lid_layer_state.inflow_value[uf * nsp + K_TEMP as usize];
            let mut t_old = [0.0f64; K_NL];
            for k in 0..K_NL {
                let st = &ctx.lid_layer_state;
                t_old[k] = st.value[st.layer_index(flat, lid_layer_of(k), K_TEMP)];
            }
            // 1. Advection
            let mut tk = [0.0f64; K_NL];
            let mut live = [false; K_NL];
            let mut mass = [false; K_NL];
            for k in 0..K_NL {
                let vidx = uf * K_NL + k;
                let v_old = ctx.lid_layer_state.vol_prev[vidx];
                let v_in = cmax(0.0, in_rate[k]) * dt;
                live[k] = donor[k] != K_ABSENT && (v_old > K_TINY_VOL || v_in > K_TINY_VOL);
                mass[k] = thick[k] > 0.0 && (if k as i32 == K_SURF { live[k] } else { donor[k] != K_ABSENT });
                if !live[k] {
                    tk[k] = if mass[k] { t_old[k] } else { dry_temperature(ctx, t_old[k], t_air) };
                    ctx.lid_layer_state.vol_prev[vidx] = v_new[k];
                    continue;
                }
                let mut tt = t_old[k];
                if v_in > K_TINY_VOL {
                    let t_in = if donor[k] == K_EXTERNAL {
                        t_ext
                    } else if donor[k] >= 0 {
                        t_old[donor[k] as usize]
                    } else {
                        t_ext
                    };
                    let den = v_old + v_in;
                    tt = if den > K_TINY_VOL { (tt * v_old + t_in * v_in) / den } else { t_in };
                }
                tk[k] = tt;
                ctx.lid_layer_state.vol_prev[vidx] = v_new[k];
            }
            if !do_cond && !do_flux {
                for k in 0..K_NL {
                    let li = ctx.lid_layer_state.layer_index(flat, lid_layer_of(k), K_TEMP);
                    ctx.lid_layer_state.value[li] = tk[k];
                }
            } else {
                // 2. One coupled thermal solve over the layers with thermal mass
                let mut idx = [0usize; K_NL];
                let mut n = 0usize;
                for k in 0..K_NL {
                    if mass[k] {
                        idx[n] = k;
                        n += 1;
                    }
                }
                if n == 0 {
                    for k in 0..K_NL {
                        let li = ctx.lid_layer_state.layer_index(flat, lid_layer_of(k), K_TEMP);
                        ctx.lid_layer_state.value[li] = tk[k];
                    }
                    continue;
                }
                let mut cap = [0.0f64; K_NL];
                let mut cond = [0.0f64; K_NL];
                for i in 0..n {
                    let k = idx[i];
                    let water_m = v_new[k] * K_FT_TO_M;
                    cap[i] = heat_capacity(&cc, k as i32, water_m, thick[k], rho_w, cp_w) * thick[k];
                    cond[i] = conductivity(&cc, k as i32, water_m, thick[k]);
                }
                let mut a = [0.0f64; K_NL];
                let mut b = [0.0f64; K_NL];
                let mut c = [0.0f64; K_NL];
                let mut rhs = [0.0f64; K_NL];
                let mut out = [0.0f64; K_NL];
                for i in 0..n {
                    b[i] = if cap[i] > 0.0 { cap[i] } else { 1.0 };
                    rhs[i] = b[i] * tk[idx[i]];
                }
                if do_cond {
                    let mut i = 0usize;
                    while i + 1 < n {
                        let r = 0.5 * thick[idx[i]] / cmax(cond[i], 1.0e-12) + 0.5 * thick[idx[i + 1]] / cmax(cond[i + 1], 1.0e-12);
                        if r > 0.0 {
                            let h = dt / r;
                            b[i] += h;
                            c[i] = -h;
                            b[i + 1] += h;
                            a[i + 1] = -h;
                        }
                        i += 1;
                    }
                }
                if do_flux {
                    let top = idx[0];
                    let el = HeatElement { kind: HeatElemKind::Lid, index: flat };
                    let j0 = net_flux_out(ctx, &el, tk[top]);
                    let j1 = net_flux_out(ctx, &el, tk[top] + K_PROBE_C);
                    let djdt = (j1 - j0) / K_PROBE_C;
                    if djdt > 0.0 {
                        b[0] += dt * djdt;
                        rhs[0] += dt * (djdt * tk[top] - j0);
                    } else {
                        rhs[0] -= dt * j0;
                    }
                }
                solve_tridiagonal(n, &a, &b, &c, &rhs, &mut out);
                for k in 0..K_NL {
                    let li = ctx.lid_layer_state.layer_index(flat, lid_layer_of(k), K_TEMP);
                    ctx.lid_layer_state.value[li] = tk[k];
                }
                for i in 0..n {
                    if out[i].is_finite() {
                        let li = ctx.lid_layer_state.layer_index(flat, lid_layer_of(idx[i]), K_TEMP);
                        ctx.lid_layer_state.value[li] = out[i];
                    }
                }
            }
            let from_surface = g.type_ == LIDType::RoofDiscon;
            let src = if from_surface { LidLayer::Surface } else { LidLayer::Storage };
            let st = &mut ctx.lid_layer_state;
            let di = uf * st.n_species as usize + K_TEMP as usize;
            st.drain_value[di] = st.value[st.layer_index(flat, src, K_TEMP)];
        }
    }
    Ok(())
}
