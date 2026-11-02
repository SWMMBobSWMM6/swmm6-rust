// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/surface/UebSubcatchments.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// UEB snow on subcatchments, the part initHydrology reaches: initializeSubcatchments
// and the file-static helpers it uses (curveFor, area, publish). NativeState is in
// ueb_types.rs. Errors are Bytes, empty on success.
//
// NOT YET TRANSLATED (run-time step functions, phase 1c): nextSnowBoundary,
// advanceSubcatchments, and the file-static transfer and plow.

use crate::core::cmath::cmax;
use crate::core::context::SimulationContext;
use crate::core::options::Bytes;
use crate::core::unit_conversion as ucf;
use crate::hydrology::snow::{SnowSoA, SNOW_PERV};
use crate::hydrology::ueb_lumped::initialize_lumped;
use crate::hydrology::ueb_types::{HeldState, LumpedFluxes, NativeState, SnowKind};

/// The native depletion curve for one snow surface kind (SNOW_PERV selects pervious,
/// every other kind impervious).
pub(crate) fn curve_for(s: &SnowSoA, surface: i32) -> [f64; 10] {
    if surface == SNOW_PERV {
        s.adc_perv
    } else {
        s.adc_imperv
    }
}

/// A native snow column's non-LID horizontal footprint, m2.
pub(crate) fn area(ctx: &SimulationContext, s: &SnowSoA, column: i32) -> f64 {
    let sc = column / 3;
    let lid = if (sc as i64) < ctx.subcatches.total_lid_area_ft2.len() as i64 { ctx.subcatches.total_lid_area_ft2[sc as usize] } else { 0.0 };
    cmax(0.0, ctx.subcatches.area[sc as usize] / ucf::ucf(ucf::LANDAREA, &ctx.options) - lid) * s.fArea[column as usize] * 0.3048 * 0.3048
}

/// Mirror accepted UEB storage and the latest flux into the native report arrays.
pub(crate) fn publish(s: &mut SnowSoA, state: &NativeState, i: i32, dt: f64) {
    let ui = i as usize;
    let column = &state.columns[ui].accepted;
    s.wsnow[ui] = column.water / 0.3048;
    s.fw[ui] = 0.0;
    s.asc[ui] = column.cover;
    s.coldc[ui] = cmax(0.0, -column.point.energy * column.cover) / (333500.0 * 0.3048);
    if dt > 0.0 {
        s.imelt[ui] = state.last[ui].liquid / 0.3048 / dt;
        s.external_net[ui] = (state.last[ui].liquid + state.last[ui].throughfall) / 0.3048 / dt;
    }
}

pub fn initialize_subcatchments(ctx: &mut SimulationContext, s: &mut SnowSoA) -> Bytes {
    let mut next = NativeState::default();
    next.columns.resize(s.wsnow.len(), HeldState::default());
    next.last.resize(s.wsnow.len(), LumpedFluxes::default());
    next.sample_end.resize(s.wsnow.len(), 0.0);
    next.sample_seconds.resize(s.wsnow.len(), 0.0);
    next.rain.clear();
    next.rain.resize(s.n_subcatch as usize, -1.0);
    next.snow.clear();
    next.snow.resize(s.n_subcatch as usize, -1.0);
    let mut external = s.external.clone();
    {
        let config = &ctx.snowpacks.config;
        for sc in 0..s.n_subcatch {
            let usc = sc as usize;
            let pack = ctx.subcatches.snowpack[usc];
            external[usc] = pack >= 0 && config.model_for(pack) == SnowKind::Ueb;
            if !external[usc] {
                continue;
            }
            next.active = true;
            let mut fractions = 0.0;
            for k in 0..5 {
                let fraction = s.sfrac[(sc * 5 + k) as usize];
                if !fraction.is_finite() || fraction < 0.0 || fraction > 1.0 {
                    return b"Invalid UEB snow removal fraction.".to_vec();
                }
                fractions += fraction;
            }
            if fractions > 1.0 + 1e-12 {
                return b"UEB snow removal fractions exceed one.".to_vec();
            }
            let initial_energy = match config.packs.get(&pack) {
                None => 0.0,
                Some(found) => found.initial_energy,
            };
            for k in 0..3 {
                let i = sc * 3 + k;
                let ui = i as usize;
                if s.fArea[ui] <= 0.0 {
                    continue;
                }
                let curve = curve_for(s, k);
                let surfaces = if k == 0 {
                    &ctx.snowpacks.plowable
                } else if k == 1 {
                    &ctx.snowpacks.impervious
                } else {
                    &ctx.snowpacks.pervious
                };
                let free_water = if (pack as i64) < surfaces.len() as i64 { surfaces[pack as usize][5] / ucf::ucf(ucf::RAINDEPTH, &ctx.options) * 0.3048 } else { s.fw[ui] * 0.3048 };
                let water = s.wsnow[ui] * 0.3048 + free_water;
                let energy = initial_energy + free_water * 333500.0;
                let error = initialize_lumped(&mut next.columns[ui].accepted, water, energy, s.si[ui] * 0.3048, &curve);
                if !error.is_empty() {
                    let mut m = b"Snowpack ".to_vec();
                    m.extend_from_slice(&ctx.snowpacks.names[pack as usize]);
                    m.extend_from_slice(b": ");
                    m.extend_from_slice(&error);
                    return m;
                }
                next.initial_water += water * area(ctx, s, i);
                let l = &next.columns[ui].accepted;
                let e = (if l.cover > 0.0 { l.cover } else { 1.0 }) * l.point.energy * area(ctx, s, i);
                next.initial_energy += e;
            }
        }
    }
    ctx.snowpacks.runtime = next;
    s.external = external;
    for sc in 0..s.n_subcatch {
        if s.external[sc as usize] {
            for k in 0..3 {
                publish(s, &ctx.snowpacks.runtime, sc * 3 + k, 0.0);
            }
        }
    }
    Bytes::new()
}
