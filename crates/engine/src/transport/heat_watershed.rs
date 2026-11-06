// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/components/HeatModule/HeatWatershed.cpp and
// HeatWatershed.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright
// 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 6 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// Phase H5a, temperature on the runoff clock: each subarea's ponded water exchanges heat with
// the atmosphere (the shared flux evaluators), mixes with what arrived (precipitation and
// run-on, flow-weighted), and the runoff leaves at the volume-weighted mean of the wet
// subareas; addRunonTemperature(At) books the run-on's temperature beside its flow.
//
// Differences in form, never in result:
//   * routeSubcatchmentTemperature returns Result: updateSolarForcing can throw or meet code
//     not translated yet.

use crate::core::context::SimulationContext;
use crate::core::engine_step_routing::StepError;
use crate::data::heat::{DryTempPolicy, HeatConfigData, HeatSource, HeatState};
use crate::data::heat_override::HeatElement;
use crate::hydrology::runoff::RunoffSoA;
use crate::transport::heat_fluxes::{air_temp_celsius, net_flux_out, relax_t, update_solar_forcing, K_CU_FT_TO_CU_M, K_PROBE_C, K_SQ_FT_TO_SQ_M};
use crate::transport::water_age_watershed::{arriving_precip_rate, arriving_precip_temperature};

const K_NSUB: usize = HeatState::K_N_SUB_AREA as usize;
const K_TINY_VOL: f64 = 1.0e-12;

fn dry_temperature(ctx: &SimulationContext, held: f64, t_air_c: f64) -> f64 {
    match ctx.heat_config.dry_temp_policy {
        DryTempPolicy::Air => t_air_c,
        DryTempPolicy::Default => HeatConfigData::K_DEFAULT_TEMP,
        DryTempPolicy::Hold => held,
    }
}

/// addRunonTemperature: the donor subcatchment's runoff temperature, at rate q.
pub fn add_runon_temperature(ctx: &mut SimulationContext, donor_sc: i32, receiver_sc: i32, q: f64) {
    if !ctx.options.heat_transport {
        return;
    }
    let hs = &ctx.heat_state;
    let ud = donor_sc as usize;
    let ur = receiver_sc as usize;
    if ud >= hs.subcatch_runoff_temp.len() || ur >= hs.subcatch_runon_temp_vol_in.len() || ur >= hs.subcatch_runon_temp_rate.len() {
        return;
    }
    let t = hs.subcatch_runoff_temp[ud];
    add_runon_temperature_at(ctx, receiver_sc, q, t);
}

/// addRunonTemperatureAt: q * temp into the numerator and q into the rate, together.
pub fn add_runon_temperature_at(ctx: &mut SimulationContext, receiver_sc: i32, q: f64, temp_c: f64) {
    if !ctx.options.heat_transport || !(q > 0.0) {
        return;
    }
    let hs = &mut ctx.heat_state;
    let ur = receiver_sc as usize;
    if ur >= hs.subcatch_runon_temp_vol_in.len() || ur >= hs.subcatch_runon_temp_rate.len() {
        return;
    }
    hs.subcatch_runon_temp_vol_in[ur] += q * temp_c;
    hs.subcatch_runon_temp_rate[ur] += q;
}

/// routeSubcatchmentTemperature
pub fn route_subcatchment_temperature(ctx: &mut SimulationContext, soa: &RunoffSoA, dt: f64) -> Result<(), StepError> {
    if !ctx.options.heat_transport || dt <= 0.0 {
        return Ok(());
    }
    let nsc = ctx.n_subcatches();
    if nsc <= 0 {
        return Ok(());
    }
    if !ctx.heat_state.watershed_sized(nsc) {
        let t0 = ctx.heat_config.global_temp[HeatSource::InitialState as usize];
        ctx.heat_state.resize_watershed(nsc, t0);
    }
    let t_rain = ctx.heat_config.global_temp[HeatSource::Rainfall as usize];
    let t_air = air_temp_celsius(ctx);
    let rho = ctx.options.water_density;
    let cp = ctx.options.water_specific_heat;
    let do_surface = ctx.heat_config.surface_exchange;
    let do_radiative = ctx.heat_config.radiative_exchange;
    update_solar_forcing(ctx, dt)?;

    let mut i = 0;
    while i < nsc && (i as usize) < soa.area.len() {
        let ui = i as usize;
        let area = soa.area[ui];
        let fi = soa.imperv_pct[ui];
        let fp = 1.0 - fi;
        let f0 = soa.frac_imperv0[ui];
        let f1 = soa.frac_imperv1[ui];
        let frac = [f0, f1, fp];
        let depth = [soa.depth_imperv0[ui], soa.depth_imperv1[ui], soa.depth_perv[ui]];
        let runon_rate = if ui < ctx.subcatches.runon_inflow.len() { ctx.subcatches.runon_inflow[ui] } else { 0.0 };
        let hs = &ctx.heat_state;
        let known_rate = if ui < hs.subcatch_runon_temp_rate.len() { hs.subcatch_runon_temp_rate[ui] } else { 0.0 };
        let runon_temp = if known_rate > K_TINY_VOL && ui < hs.subcatch_runon_temp_vol_in.len() { hs.subcatch_runon_temp_vol_in[ui] / known_rate } else { t_rain };
        let runon_depth_rate = if area > 0.0 { runon_rate / area } else { 0.0 };
        let mut out_num = 0.0;
        let mut out_den = 0.0;
        for k in 0..K_NSUB {
            let idx = ui * K_NSUB + k;
            if idx >= ctx.heat_state.subarea_temp.len() {
                break;
            }
            let rain_rate = arriving_precip_rate(ctx, ui, k as i32);
            let in_rate = rain_rate + runon_depth_rate;
            let t_precip = arriving_precip_temperature(ctx, ui, k as i32);
            let t_in = if in_rate > 0.0 { (rain_rate * t_precip + runon_depth_rate * runon_temp) / in_rate } else { t_precip };
            let v_new = depth[k] * frac[k] * area;
            let v_old = ctx.heat_state.subarea_vol_prev[idx];
            let v_in = in_rate * frac[k] * area * dt;
            let a_ft2 = frac[k] * area;
            let mut t = ctx.heat_state.subarea_temp[idx];
            if v_old <= K_TINY_VOL && v_in <= K_TINY_VOL {
                ctx.heat_state.subarea_temp[idx] = dry_temperature(ctx, t, t_air);
                ctx.heat_state.subarea_vol_prev[idx] = v_new;
                continue;
            }
            if v_old > K_TINY_VOL && (do_surface || do_radiative) {
                let er = HeatElement::subcatch(i);
                t += relax_t(net_flux_out(ctx, &er, t), net_flux_out(ctx, &er, t + K_PROBE_C), K_PROBE_C, a_ft2 * K_SQ_FT_TO_SQ_M, v_old * K_CU_FT_TO_CU_M, dt, rho, cp);
            }
            if v_in > K_TINY_VOL {
                let denom = v_old + v_in;
                t = if denom > K_TINY_VOL { (t * v_old + t_in * v_in) / denom } else { t_in };
            }
            ctx.heat_state.subarea_temp[idx] = t;
            ctx.heat_state.subarea_vol_prev[idx] = v_new;
            if v_new > K_TINY_VOL {
                out_num += t * v_new;
                out_den += v_new;
            }
        }
        if ui < ctx.heat_state.subcatch_runoff_temp.len() {
            ctx.heat_state.subcatch_runoff_temp[ui] = if out_den > K_TINY_VOL { out_num / out_den } else { t_rain };
        }
        i += 1;
    }
    Ok(())
}
