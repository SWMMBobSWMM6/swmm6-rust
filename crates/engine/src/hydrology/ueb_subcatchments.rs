// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/surface/UebSubcatchments.cpp and
// SnowCheckpoint.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// UEB snow on subcatchments: initializeSubcatchments and the file-static helpers it uses
// (curveFor, area, publish); and (6 November 2026) the run-time step: nextSnowBoundary,
// advanceSubcatchments, the file-static transfer and plow, and SnowCheckpoint.cpp's
// snowRuntimeForcing and snowTransientForcing (the only parts of that file the step uses).
// NativeState is in ueb_types.rs. Errors are Bytes, empty on success.
//
// Differences of form, never in result:
//   * advance_subcatchments returns Result<Bytes, String>: Ok(error) is upstream's returned
//     std::string; Err(what) is the std::invalid_argument surface::solarIntegral throws
//     inside climate::resolveRadiation, which nothing in this file catches (the caller
//     turns both into the same throw, with or without the "UEB subcatchment N: " text).
//   * The perf::ProcessTimer is left out (timing only).
//   * advanceSubcatchments reads `held` (next.columns[i]) and `flux` (next.last[i]) through
//     references; here they are indexed in place.

use crate::core::cmath::{cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::options::Bytes;
use crate::core::unit_conversion as ucf;
use crate::data::forcing::{ForcingMode, ForcingPersist};
use crate::hydrology::gage::gage_next_rate_change;
use crate::hydrology::snow::{SnowSoA, SNOW_PERV};
use crate::hydrology::ueb_forcing::resolve_ueb_forcing;
use crate::hydrology::ueb_held::{close_held, consume_held, prepare_held};
use crate::hydrology::ueb_lumped::{depletion_cover, initialize_lumped};
use crate::hydrology::ueb_snow::bulk_temperature;
use crate::hydrology::ueb_types::{Forcing, HeldState, LumpedFluxes, NativeState, Site, SnowKind};

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

/// transfer: add (or take) plowed water and its heat to a column's accepted inventory and
/// abandon any held interval.
fn transfer(held: &mut HeldState, water: f64, heat: f64, curve: &[f64; 10]) -> Bytes {
    let mut next = held.accepted.clone();
    let total_energy = (if next.cover > 0.0 { next.cover } else { 1.0 }) * next.point.energy + heat;
    next.water += water;
    if !next.water.is_finite() || next.water < -1e-12 || !total_energy.is_finite() {
        return b"Invalid plowed snow transfer.".to_vec();
    }
    next.water = cmax(0.0, next.water);
    next.maximum = if next.water > 0.0 { cmax(next.maximum, next.water) } else { 0.0 };
    next.cover = if next.water > 0.0 { depletion_cover(curve, next.water / next.maximum) } else { 0.0 };
    if next.water > 0.0 && next.cover <= 0.0 {
        return b"Plowed snow has no receiving covered area.".to_vec();
    }
    next.point.water = if next.cover > 0.0 { next.water / next.cover } else { 0.0 };
    next.point.energy = total_energy / (if next.cover > 0.0 { next.cover } else { 1.0 });
    held.accepted = next;
    close_held(held);
    Bytes::new()
}

/// plow: move plowable snow above weplow to the removal destinations.
fn plow(ctx: &SimulationContext, state: &mut NativeState, soa: &mut SnowSoA, dt: f64) -> Bytes {
    let p = &ctx.snowpacks.config.parameters;
    for sc in 0..soa.n_subcatch {
        let usc = sc as usize;
        if soa.enabled[usc] == 0 || !soa.external[usc] {
            continue;
        }
        let donor = sc * 3;
        let ud = donor as usize;
        let a = area(ctx, soa, donor);
        let inventory = state.columns[ud].accepted.clone();
        if a <= 0.0 || inventory.water < soa.weplow[usc] * 0.3048 {
            continue;
        }
        let age = soa.age[ud];
        let temperature = bulk_temperature(&inventory.point, p);
        let specific = if inventory.point.energy <= 0.0 {
            1000.0 * 2.09 * temperature
        } else if inventory.point.water > 0.0 {
            cmin(333500.0, inventory.point.energy / inventory.point.water)
        } else {
            0.0
        };
        let mut moved = 0.0;
        for destination in 0..5 {
            let depth = soa.sfrac[(sc * 5 + destination) as usize] * inventory.water;
            if depth <= 0.0 {
                continue;
            }
            let volume = depth * a;
            let heat = volume * specific;
            if destination == 0 {
                soa.removed += volume / (0.3048 * 0.3048 * 0.3048);
                state.plowed_energy += heat;
                state.removed += volume;
            } else if destination == 3 {
                let flux = &mut state.last[ud];
                if soa.track_age {
                    soa.out_age[ud] = (flux.liquid * soa.out_age[ud] + depth * age) / (flux.liquid + depth);
                }
                flux.liquid += depth;
                state.liquid += volume;
                state.plowed_energy += heat;
            } else {
                let receiver = if destination == 1 {
                    donor + 1
                } else if destination == 2 {
                    donor + 2
                } else {
                    soa.to_subcatch[usc] * 3 + 2
                };
                if receiver < 0 || receiver >= state.columns.len() as i32 || receiver == donor {
                    continue;
                }
                let ur = receiver as usize;
                let receiving_area = area(ctx, soa, receiver);
                if receiving_area <= 0.0 {
                    continue;
                }
                let added = volume / receiving_area;
                let previous = soa.wsnow[ur] * 0.3048 + soa.fw[ur] * 0.3048;
                if soa.track_age {
                    soa.age[ur] = (previous * soa.age[ur] + added * age) / (previous + added);
                }
                if soa.external[(receiver / 3) as usize] {
                    let curve = curve_for(soa, receiver % 3);
                    let error = transfer(&mut state.columns[ur], added, heat / receiving_area, &curve);
                    if !error.is_empty() {
                        return error;
                    }
                    publish(soa, state, receiver, dt);
                } else {
                    soa.wsnow[ur] += added / 0.3048;
                    soa.coldc[ur] += cmax(0.0, -heat / receiving_area) / (333500.0 * 0.3048);
                    state.plowed_energy += heat;
                    state.transferred_water -= volume;
                }
            }
            moved += depth;
        }
        if moved > 0.0 {
            let curve = curve_for(soa, 0);
            let error = transfer(&mut state.columns[ud], -moved, -moved * specific, &curve);
            if !error.is_empty() {
                return error;
            }
            if state.columns[ud].accepted.water <= 0.0 {
                soa.age[ud] = 0.0;
            }
            publish(soa, state, donor, dt);
        }
    }
    Bytes::new()
}

/// snowRuntimeForcing (SnowCheckpoint.cpp): the modes and values of the runtime climate
/// overrides, the key a held interval was prepared under.
pub fn snow_runtime_forcing(ctx: &SimulationContext) -> [f64; 4] {
    let f = &ctx.forcing;
    [
        f.climate_temperature_mode as i32 as f64,
        if f.climate_temperature_mode == ForcingMode::None { 0.0 } else { f.climate_temperature_value },
        f.climate_wind_mode as i32 as f64,
        if f.climate_wind_mode == ForcingMode::None { 0.0 } else { f.climate_wind_value },
    ]
}

/// snowTransientForcing (SnowCheckpoint.cpp): a climate override that resets after a step.
pub fn snow_transient_forcing(ctx: &SimulationContext) -> bool {
    let f = &ctx.forcing;
    (f.climate_temperature_mode != ForcingMode::None && f.climate_temperature_persist == ForcingPersist::Reset)
        || (f.climate_wind_mode != ForcingMode::None && f.climate_wind_persist == ForcingPersist::Reset)
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

/// nextSnowBoundary: the runoff step shortened to the next UEB snow step (and hourly
/// meteorology hour) boundary.
pub fn next_snow_boundary(ctx: &SimulationContext, elapsed: f64, mut proposed: f64) -> f64 {
    if !ctx.snowpacks.runtime.active || ctx.options.ignore_snow_melt {
        return proposed;
    }
    let step = ctx.snowpacks.config.step;
    let next = (((elapsed + 1e-7) / step).floor() + 1.0) * step;
    proposed = cmin(proposed, next - elapsed);
    if !ctx.snowpacks.meteorology.is_empty() {
        let start = ctx.options.start_date;
        let hour = (start - start.floor()) * 86400. + elapsed;
        let boundary = (((hour + 1e-5) / 3600.).floor() + 1.0) * 3600.;
        proposed = cmin(proposed, boundary - hour);
    }
    proposed
}

/// advanceSubcatchments: advance every UEB subcatchment column over one runoff step (staged;
/// committed only when every column and the plowing succeed).
pub fn advance_subcatchments(ctx: &mut SimulationContext, soa: &mut SnowSoA, date: f64, elapsed: f64, seconds: f64, rain: &[f64], snow: &[f64]) -> Result<Bytes, String> {
    if !ctx.snowpacks.runtime.active {
        return Ok(Bytes::new());
    }
    let mut next = ctx.snowpacks.runtime.clone();
    let mut staged = soa.clone();
    let forcing_key = snow_runtime_forcing(ctx);
    if forcing_key != next.forcing {
        for column in next.columns.iter_mut() {
            close_held(column);
        }
    }
    next.forcing = forcing_key;
    let config_step = ctx.snowpacks.config.step;
    let parameters = ctx.snowpacks.config.parameters;
    let mut evaporation_ft3 = 0.0;
    for sc in 0..soa.n_subcatch {
        let usc = sc as usize;
        if soa.enabled[usc] == 0 || !soa.external[usc] {
            continue;
        }
        let pack = ctx.subcatches.snowpack[usc];
        let mut site = Site::default();
        if let Some(found) = ctx.snowpacks.config.packs.get(&pack) {
            site = found.site;
        }
        site.drift = 1.0; // already applied by the shared gage partition
        let mut horizon = next_snow_boundary(ctx, elapsed, config_step);
        let gage = ctx.subcatches.gage[usc];
        if gage >= 0 {
            let limit = date + horizon / 86400.0;
            let change = gage_next_rate_change(ctx, gage, date, limit);
            if change < limit && change > date + 1e-10 {
                horizon = cmin(horizon, (change - date) * 86400.0);
            }
        }
        horizon = cmin(horizon, ctx.options.total_duration_ms() / 1000.0 - elapsed);
        if snow_transient_forcing(ctx) {
            horizon = cmin(horizon, seconds);
        }
        if horizon < seconds && seconds - horizon < 1e-5 {
            horizon = seconds;
        }
        if horizon < seconds {
            return Ok(b"Runoff step crosses a UEB forcing boundary.".to_vec());
        }
        let mut forcing = Forcing::default();
        let mut resolved = false;
        for k in 0..3 {
            let i = sc * 3 + k;
            let ui = i as usize;
            next.last[ui] = LumpedFluxes::default();
            if soa.fArea[ui] <= 0.0 {
                continue;
            }
            let incoming = (soa.wsnow[ui] + soa.fw[ui]) * 0.3048 - next.columns[ui].accepted.water;
            if incoming.abs() > 1e-12 {
                let heat = staged.ueb_transfer_energy[ui];
                let error = transfer(&mut next.columns[ui], incoming, heat, &curve_for(soa, k));
                if !error.is_empty() {
                    return Ok(error);
                }
                staged.ueb_transfer_energy[ui] = 0.0;
                next.transferred_water += incoming * area(ctx, soa, i);
                next.plowed_energy -= heat * area(ctx, soa, i);
            }
            if next.rain[usc] != rain[usc] || next.snow[usc] != snow[usc] {
                close_held(&mut next.columns[ui]);
            }
            if !next.columns[ui].pending() {
                if !resolved {
                    let error = resolve_ueb_forcing(ctx, date, horizon, 0.0, 0.0, &mut forcing, &mut site.pressure)?;
                    if !error.is_empty() {
                        return Ok(error);
                    }
                    forcing.precipitation = (rain[usc] + snow[usc]) * 0.3048 * horizon;
                    forcing.snow_fraction = if rain[usc] + snow[usc] > 0.0 { snow[usc] / (rain[usc] + snow[usc]) } else { 0.0 };
                    forcing.snow_catch = 1.0;
                    resolved = true;
                }
                let error = prepare_held(&mut next.columns[ui], &parameters, &site, &forcing, &curve_for(soa, k));
                if !error.is_empty() {
                    let mut m = b"UEB subcatchment ".to_vec();
                    m.extend_from_slice(sc.to_string().as_bytes());
                    m.extend_from_slice(b": ");
                    m.extend_from_slice(&error);
                    return Ok(m);
                }
            }
            {
                let NativeState { columns, last, .. } = &mut next;
                let error = consume_held(&mut columns[ui], seconds, &curve_for(soa, k), &mut last[ui]);
                if !error.is_empty() {
                    return Ok(error);
                }
            }
            next.sample_end[ui] = elapsed + seconds;
            next.sample_seconds[ui] = seconds;
            let a = area(ctx, soa, i);
            let flux = next.last[ui];
            next.precipitation += (flux.rain + flux.snow) * a;
            next.liquid += flux.liquid * a;
            next.throughfall += flux.throughfall * a;
            next.sublimation += flux.sublimation * a;
            next.net_energy += flux.net_energy * a;
            next.geometry_energy += flux.geometry_energy * a;
            evaporation_ft3 += flux.sublimation * a / (0.3048 * 0.3048 * 0.3048);
            if staged.track_age {
                let initial = soa.wsnow[ui] + soa.fw[ui];
                let arrived = (flux.rain + flux.snow - flux.throughfall) / 0.3048;
                let volume = initial + arrived;
                let age = if volume > 0.0 { (initial * (soa.age[ui] + seconds) + arrived * soa.precip_age) / volume } else { 0.0 };
                staged.out_age[ui] = if flux.liquid > 0.0 { age } else { 0.0 };
                staged.age[ui] = if next.columns[ui].accepted.water > 0.0 { age } else { 0.0 };
            }
            publish(&mut staged, &next, i, seconds);
        }
        next.rain[usc] = rain[usc];
        next.snow[usc] = snow[usc];
    }
    let error = plow(ctx, &mut next, &mut staged, seconds);
    if !error.is_empty() {
        return Ok(error);
    }
    next.last_seconds = seconds;
    next.completed_seconds = elapsed + seconds;
    ctx.snowpacks.runtime = next;
    *soa = staged;
    if evaporation_ft3 >= 0.0 {
        ctx.mass_balance.runoff_evap += evaporation_ft3;
    } else {
        ctx.mass_balance.runoff_rainfall -= evaporation_ft3; // atmospheric deposition is an input
    }
    Ok(Bytes::new())
}
