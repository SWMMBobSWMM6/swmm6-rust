// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/surface/UebSnow.cpp and UebSnow.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin,
// licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 6 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The UEB point snow model's driver: advance (validation, the max_step partition), boundedStep
// (halving on a failed step, to depth 16), step (the reference kernel call and its
// conservation checks), the 24-hour thermal history, bulkTemperature. Errors are Bytes, empty
// on success (std::string upstream).
//
// NOT YET TRANSLATED: the EQUILIBRIUM surface temperature method's kernel
// (ueb_reference.rs); ueb::advanceSubcatchments refuses a deck that selects it before any
// state changes, so `step` never reaches it.

use crate::core::cmath::cmax;
use crate::core::options::Bytes;
use crate::hydrology::ueb_reference::force_restore_step;
use crate::hydrology::ueb_types::{Fluxes, Forcing, Parameters, Site, State, TemperatureMethod, ThermalSample};

const WATER_DENSITY: f64 = 1000.0;
const LATENT_HEAT: f64 = 333.5;
const ICE_HEAT: f64 = 2.09;

fn finite(v: f64) -> bool {
    v.is_finite()
}
fn nonnegative(v: f64) -> bool {
    finite(v) && v >= 0.0
}
fn fraction(v: f64) -> bool {
    nonnegative(v) && v <= 1.0
}

fn mean(state: &State, surface: bool, fallback: f64) -> f64 {
    let mut sum = 0.0;
    let mut time = 0.0;
    for h in state.history.iter() {
        time += h.seconds;
        sum += h.seconds * (if surface { h.surface } else { h.bulk });
    }
    if time > 0.0 {
        sum / time
    } else {
        fallback
    }
}

fn remember(state: &mut State, dt: f64, ts: f64, tb: f64) {
    state.history.push(ThermalSample { seconds: dt, surface: ts, bulk: tb });
    let mut duration = 0.0;
    for h in state.history.iter() {
        duration += h.seconds;
    }
    let mut excess = cmax(0.0, duration - 86400.0);
    let mut erase = 0usize;
    while erase < state.history.len() && excess >= state.history[erase].seconds {
        excess -= state.history[erase].seconds;
        erase += 1;
    }
    state.history.drain(0..erase);
    if !state.history.is_empty() {
        state.history[0].seconds -= excess;
    }
}

fn state_error(s: &State) -> Bytes {
    if !finite(s.energy) || !nonnegative(s.water) || !nonnegative(s.age) || !nonnegative(s.refreezing_depth) || !nonnegative(s.total_refreezing_depth) {
        return b"Invalid UEB snow inventory.".to_vec();
    }
    let mut duration = 0.0;
    for h in s.history.iter() {
        if !finite(h.seconds) || h.seconds <= 0.0 || !finite(h.surface) || h.surface <= -273.15 || !finite(h.bulk) || h.bulk <= -273.15 {
            return b"Invalid UEB thermal history.".to_vec();
        }
        duration += h.seconds;
    }
    if !finite(duration) || duration > 86400.0 + 1e-7 {
        return b"UEB thermal history exceeds 24 hours.".to_vec();
    }
    Vec::new()
}

/// step: one kernel step; Ok(false) for a failed finite-state or conservation check.
fn step(state: &mut State, p: &Parameters, site: &Site, f: &Forcing, flux: &mut Fluxes) -> Result<bool, String> {
    let hours = f.seconds / 3600.0;
    let mut input = [0.0f64; 15];
    input[8] = f.temperature;
    input[9] = f.precipitation / hours;
    input[10] = f.wind;
    input[11] = f.humidity;
    input[12] = f.shortwave * 3.6;
    input[13] = f.longwave * 3.6;
    input[14] = f.cosine;
    let location = [0.0, site.forest, site.drift * f.snow_catch, site.pressure, site.ground_heat * 3.6, site.extinction_depth];
    let mut water = [0.0, state.energy, state.water, state.age, state.refreezing_depth, state.total_refreezing_depth];
    let mut output = [0.0f64; 15];
    let diagnostics;
    if p.temperature == TemperatureMethod::Equilibrium {
        // refused by ueb::advanceSubcatchments before any state changes
        panic!("STUB: the UEB EQUILIBRIUM kernel (ueb_reference::equilibriumStep)");
    } else {
        let parameters: [f64; 23] = [0.0, p.rain_temperature, p.snow_temperature, p.emissivity, p.soil_heat_capacity, p.measurement_height, p.roughness, p.snow_density, p.soil_density, p.liquid_retention, p.drainage_conductivity, p.active_soil_depth, p.ground_albedo, p.visible_albedo, p.infrared_albedo, p.snow_conductivity, p.soil_conductivity, p.low_frequency, p.amplitude, p.stability, 0.0, p.new_snow_depth, 0.0];
        let initial_bulk = bulk_temperature(state, p);
        let old_surface = match state.history.last() {
            None => crate::core::cmath::cmin(0.0, initial_bulk),
            Some(h) => h.surface,
        };
        let old_bulk = match state.history.last() {
            None => initial_bulk,
            Some(h) => h.bulk,
        };
        let mut surface = old_surface;
        let mut average = old_bulk;
        diagnostics = force_restore_step(hours, &input, &location, &mut water, &mut surface, &mut average, &parameters, &mut output, p.age_dust, mean(state, true, old_surface), mean(state, false, old_bulk), p.visible_age, p.infrared_age, f.snow_fraction)?;
    }
    let mut next = Fluxes::default();
    next.rain = output[1] * hours;
    next.snow = output[2] * hours;
    next.albedo = output[3];
    next.sublimation = output[6] * hours;
    next.liquid = output[7] * hours;
    next.outgoing_melt_energy = output[8] * hours;
    next.net_energy = output[9] * hours;
    next.bulk_temperature = output[11];
    next.surface_temperature = output[12];
    next.diagnostics = diagnostics;
    next.water_residual = state.water + next.rain + next.snow - next.sublimation - next.liquid - water[2];
    next.energy_residual = state.energy + next.net_energy - water[1];
    for x in water.iter() {
        if !finite(*x) {
            return Ok(false);
        }
    }
    for x in output.iter() {
        if !finite(*x) {
            return Ok(false);
        }
    }
    if water[2] < 0.0
        || water[3] < 0.0
        || water[4] < 0.0
        || water[5] < 0.0
        || next.liquid < 0.0
        || !fraction(next.albedo)
        || next.water_residual.abs() > 1e-10 * cmax(1.0, state.water + next.rain + next.snow)
        || next.energy_residual.abs() > 1e-10 * cmax(cmax(1.0, state.energy.abs()), next.net_energy.abs())
        || next.surface_temperature <= -273.15
        || next.bulk_temperature <= -273.15
        || diagnostics > 0
    {
        return Ok(false);
    }
    state.energy = water[1];
    state.water = water[2];
    state.age = water[3];
    state.refreezing_depth = water[4];
    state.total_refreezing_depth = water[5];
    remember(state, f.seconds, next.surface_temperature, next.bulk_temperature);
    next.steps = 1;
    *flux = next;
    Ok(true)
}

fn add(sum: &mut Fluxes, x: &Fluxes) {
    sum.rain += x.rain;
    sum.snow += x.snow;
    sum.liquid += x.liquid;
    sum.sublimation += x.sublimation;
    sum.net_energy += x.net_energy;
    sum.outgoing_melt_energy += x.outgoing_melt_energy;
    sum.water_residual += x.water_residual;
    sum.energy_residual += x.energy_residual;
    sum.steps += x.steps;
    sum.retries += x.retries;
    sum.diagnostics += x.diagnostics;
    sum.surface_temperature = x.surface_temperature;
    sum.bulk_temperature = x.bulk_temperature;
    sum.albedo = x.albedo;
}

fn bounded_step(state: &mut State, p: &Parameters, site: &Site, forcing: &Forcing, flux: &mut Fluxes, depth: i32) -> Result<bool, String> {
    let mut candidate = state.clone();
    if step(&mut candidate, p, site, forcing, flux)? {
        *state = candidate;
        return Ok(true);
    }
    if depth == 16 || forcing.seconds <= 0.001 {
        return Ok(false);
    }
    let mut half = *forcing;
    half.seconds *= 0.5;
    half.precipitation *= 0.5;
    let mut a = Fluxes::default();
    let mut b = Fluxes::default();
    if !bounded_step(&mut candidate, p, site, &half, &mut a, depth + 1)? || !bounded_step(&mut candidate, p, site, &half, &mut b, depth + 1)? {
        return Ok(false);
    }
    add(&mut a, &b);
    a.retries += 1;
    *flux = a;
    *state = candidate;
    Ok(true)
}

/// bulkTemperature: the bulk (snow and active soil) temperature, C.
pub fn bulk_temperature(s: &State, p: &Parameters) -> f64 {
    let capacity = WATER_DENSITY * s.water * ICE_HEAT + p.soil_density * p.active_soil_depth * p.soil_heat_capacity;
    if s.energy <= 0.0 {
        s.energy / capacity
    } else {
        cmax(0.0, s.energy - s.water * WATER_DENSITY * LATENT_HEAT) / capacity
    }
}

/// advance: one forcing interval, partitioned at max_step.
pub fn advance(state: &mut State, p: &Parameters, site: &Site, f: &Forcing, output: &mut Fluxes) -> Bytes {
    let error = p.validate();
    if !error.is_empty() {
        return error;
    }
    let error = state_error(state);
    if !error.is_empty() {
        return error;
    }
    if !fraction(site.forest) || !nonnegative(site.drift) || !finite(site.pressure) || site.pressure <= 0.0 || !finite(site.ground_heat) || !finite(site.extinction_depth) || site.extinction_depth <= 0.0 {
        return b"Invalid UEB site parameters.".to_vec();
    }
    if !finite(f.seconds) || f.seconds <= 0.0 || !nonnegative(f.precipitation) || !finite(f.temperature) || f.temperature <= -100.0 || f.temperature > 70.0 || !nonnegative(f.wind) || !fraction(f.humidity) || !nonnegative(f.shortwave) || !nonnegative(f.longwave) || !fraction(f.cosine) || !nonnegative(f.snow_catch) || (f.snow_fraction != -1.0 && !fraction(f.snow_fraction)) {
        return b"Invalid UEB meteorological forcing.".to_vec();
    }
    let mut candidate = state.clone();
    let mut result = Fluxes::default();
    let count = (f.seconds / p.max_step).ceil();
    if count > 100000.0 {
        return b"UEB interval exceeds the kernel step limit.".to_vec();
    }
    let steps = count as i32;
    let mut elapsed = 0.0;
    let mut precipitation = 0.0;
    for i in 0..steps {
        let end = if i + 1 == steps { f.seconds } else { crate::core::cmath::cmin(f.seconds, (i + 1) as f64 * p.max_step) };
        let mut part = *f;
        part.seconds = end - elapsed;
        part.precipitation = if i + 1 == steps { f.precipitation - precipitation } else { f.precipitation * (part.seconds / f.seconds) };
        let mut value = Fluxes::default();
        match bounded_step(&mut candidate, p, site, &part, &mut value, 0) {
            Ok(true) => {}
            Ok(false) => return b"UEB step failed its finite-state or conservation checks.".to_vec(),
            Err(what) => return what.into_bytes(),
        }
        add(&mut result, &value);
        elapsed = end;
        precipitation += part.precipitation;
    }
    *state = candidate;
    *output = result;
    Vec::new()
}
