// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/components/WaterAgeModule/
// WaterAgeWatershed.cpp and src/engine/transport/components/WatershedCommon.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Subcatchment surface water age (Phase A3): addRunonAge and routeSubcatchmentAge, with
// the precipitation helpers they read (arrivingPrecipRate, arrivingMeltFraction,
// arrivingPrecipAge). arrivingPrecipTemperature (heat) is not translated here.

use crate::core::cmath::cmax;
use crate::core::context::SimulationContext;
use crate::data::water_age::WaterAgeSource;
use crate::hydrology::runoff::RunoffSoA;

const K_SUB_PERV: i32 = 2;
const K_NSUB: usize = 3;
const K_TINY_VOL: f64 = 1.0e-12;

/// arrivingPrecipRate: the water reaching one subarea (ft/s): the gage rate, or under a
/// snow pack the snow-modified net precipitation.
pub fn arriving_precip_rate(ctx: &SimulationContext, ui: usize, subarea: i32) -> f64 {
    let gage = if ui < ctx.subcatches.rainfall.len() { ctx.subcatches.rainfall[ui] } else { 0.0 };
    if ctx.options.ignore_snow_melt || ui >= ctx.subcatches.snowpack.len() || ctx.subcatches.snowpack[ui] < 0 {
        return gage;
    }
    let src = if subarea == K_SUB_PERV { &ctx.subcatches.snow_net_perv } else { &ctx.subcatches.snow_net_imperv };
    if ui >= src.len() {
        return gage;
    }
    let net = src[ui];
    if net >= 0.0 {
        net
    } else {
        gage
    }
}

/// arrivingMeltFraction: the meltwater share of the arriving water.
pub fn arriving_melt_fraction(ctx: &SimulationContext, ui: usize, subarea: i32) -> f64 {
    if ctx.options.ignore_snow_melt || ui >= ctx.subcatches.snowpack.len() || ctx.subcatches.snowpack[ui] < 0 {
        return 0.0;
    }
    let perv = subarea == K_SUB_PERV;
    let net_v = if perv { &ctx.subcatches.snow_net_perv } else { &ctx.subcatches.snow_net_imperv };
    let melt_v = if perv { &ctx.subcatches.snow_melt_perv } else { &ctx.subcatches.snow_melt_imperv };
    if ui >= net_v.len() || ui >= melt_v.len() {
        return 0.0;
    }
    let net = net_v[ui];
    let melt = melt_v[ui];
    if !(net > 0.0) || !(melt >= 0.0) {
        return 0.0;
    }
    let f = melt / net;
    if f < 0.0 {
        0.0
    } else if f > 1.0 {
        1.0
    } else {
        f
    }
}

/// arrivingPrecipAge: the age of the arriving water (rain at the configured source age,
/// meltwater at the pack's age).
pub fn arriving_precip_age(ctx: &SimulationContext, ui: usize, subarea: i32) -> f64 {
    let a_rain = ctx.water_age_config.global_age[WaterAgeSource::Rainfall as usize];
    let f = arriving_melt_fraction(ctx, ui, subarea);
    if f <= 0.0 {
        return a_rain;
    }
    let age_v = if subarea == K_SUB_PERV { &ctx.subcatches.snow_melt_age_perv } else { &ctx.subcatches.snow_melt_age_imperv };
    if ui >= age_v.len() {
        return a_rain;
    }
    let a_pack = age_v[ui];
    if !(a_pack >= 0.0) {
        return a_rain;
    }
    (1.0 - f) * a_rain + f * a_pack
}

/// kMeltwaterTempC (WatershedCommon.hpp)
pub const K_MELTWATER_TEMP_C: f64 = 0.0;

/// arrivingPrecipTemperature: the configured rainfall temperature blended with meltwater (0
/// degC) by the arriving melt fraction.
pub fn arriving_precip_temperature(ctx: &SimulationContext, ui: usize, subarea: i32) -> f64 {
    let t_rain = ctx.heat_config.global_temp[crate::data::heat::HeatSource::Rainfall as usize];
    let f = arriving_melt_fraction(ctx, ui, subarea);
    if f <= 0.0 {
        return t_rain;
    }
    (1.0 - f) * t_rain + f * K_MELTWATER_TEMP_C
}

/// addRunonAge: add the donor's runoff-age rate (q * age) to the receiver's run-on.
pub fn add_runon_age(ctx: &mut SimulationContext, donor_sc: i32, receiver_sc: i32, q: f64) {
    if !ctx.options.water_age {
        return;
    }
    let ws = &mut ctx.water_age_state;
    let ud = donor_sc as usize;
    let ur = receiver_sc as usize;
    if ud >= ws.subcatch_runoff_age.len() || ur >= ws.subcatch_runon_age_vol_in.len() {
        return;
    }
    ws.subcatch_runon_age_vol_in[ur] += q * ws.subcatch_runoff_age[ud];
}

/// routeSubcatchmentAge: age and mix each subarea's water store; publish the runoff age.
pub fn route_subcatchment_age(ctx: &mut SimulationContext, soa: &RunoffSoA, dt: f64) {
    if !ctx.options.water_age || dt <= 0.0 {
        return;
    }
    let nsc = ctx.n_subcatches();
    if nsc <= 0 {
        return;
    }
    let want = nsc as usize * K_NSUB;
    if ctx.water_age_state.subarea_age.len() != want {
        let (nn, nl) = (ctx.n_nodes(), ctx.n_links());
        ctx.water_age_state.resize(nn, nl, nsc);
    }
    let a_rain = ctx.water_age_config.global_age[WaterAgeSource::Rainfall as usize];

    let mut i = 0usize;
    while (i as i32) < nsc && i < soa.area.len() {
        let ui = i;
        let area = soa.area[ui];
        let fi = soa.imperv_pct[ui];
        let fp = 1.0 - fi;
        let f0 = soa.frac_imperv0[ui];
        let f1 = soa.frac_imperv1[ui];
        let frac = [f0, f1, fp];
        let depth = [soa.depth_imperv0[ui], soa.depth_imperv1[ui], soa.depth_perv[ui]];

        let ws_len = ctx.water_age_state.subcatch_runon_age_vol_in.len();
        let runon_rate = if ui < ws_len && ui < ctx.subcatches.runon_inflow.len() { ctx.subcatches.runon_inflow[ui] } else { 0.0 };
        let runon_age = if runon_rate > 0.0 && ui < ws_len { ctx.water_age_state.subcatch_runon_age_vol_in[ui] / runon_rate } else { 0.0 };
        let runon_depth_rate = if area > 0.0 { runon_rate / area } else { 0.0 };

        let mut out_num = 0.0;
        let mut out_den = 0.0;
        for k in 0..K_NSUB {
            let idx = ui * K_NSUB + k;
            let v_new = depth[k] * frac[k] * area;
            let v_old = ctx.water_age_state.subarea_vol_prev[idx];
            let rain_rate = arriving_precip_rate(ctx, ui, k as i32);
            let a_precip = arriving_precip_age(ctx, ui, k as i32);
            let in_rate = rain_rate + runon_depth_rate;
            let a_in = if in_rate > 0.0 { (rain_rate * a_precip + runon_depth_rate * runon_age) / in_rate } else { a_precip };

            let mut a = ctx.water_age_state.subarea_age[idx] + dt;
            let v_in = in_rate * frac[k] * area * dt;
            if v_in > K_TINY_VOL {
                let denom = v_old + v_in;
                a = if denom > K_TINY_VOL { (a * v_old + a_in * v_in) / denom } else { a_in };
            }
            ctx.water_age_state.subarea_age[idx] = cmax(a, 0.0);
            ctx.water_age_state.subarea_vol_prev[idx] = v_new;
            if v_new > K_TINY_VOL {
                out_num += a * v_new;
                out_den += v_new;
            }
        }
        if ui < ctx.water_age_state.subcatch_runoff_age.len() {
            ctx.water_age_state.subcatch_runoff_age[ui] = if out_den > K_TINY_VOL { out_num / out_den } else { a_rain };
        }
        i += 1;
    }
}
