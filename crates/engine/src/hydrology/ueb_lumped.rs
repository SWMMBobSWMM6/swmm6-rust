// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/surface/UebLumped.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The lumped UEB snow inventory's setup: validateDepletion, depletionCover and
// initializeLumped, which ueb::initializeSubcatchments calls. The types (LumpedState,
// LumpedFluxes) are in ueb_types.rs. Errors are std::string upstream (empty = success);
// here they are Bytes, empty on success.
//
// NOT YET TRANSLATED (run-time step functions, phase 1c): advanceLumped and the
// file-static helpers it uses (soilCapacity, representedEnergy, remap).

use crate::core::cmath::{cmax, cmin};
use crate::core::options::Bytes;
use crate::hydrology::ueb_types::{Forcing, LumpedFluxes, LumpedState, Parameters, Site, State};

pub fn validate_depletion(curve: &[f64; 10]) -> Bytes {
    let mut previous = 0.0;
    for &value in curve.iter() {
        if !value.is_finite() || value < previous || value > 1.0 {
            return b"Snow depletion ordinates must be finite, nondecreasing fractions.".to_vec();
        }
        previous = value;
    }
    Bytes::new()
}

pub fn depletion_cover(curve: &[f64; 10], water: f64) -> f64 {
    if water <= 0.0 {
        return 0.0;
    }
    if water >= 1.0 {
        return 1.0;
    }
    let x = water * 10.0;
    // std::min(9, int(x)); x is in (0, 10) here, so the conversion is in range.
    let i = cmin(9, x as i32);
    let fi = i as f64;
    curve[i as usize] + (x - fi) * ((if i == 9 { 1.0 } else { curve[i as usize + 1] }) - curve[i as usize])
}

pub fn initialize_lumped(state: &mut LumpedState, water: f64, energy: f64, maximum: f64, curve: &[f64; 10]) -> Bytes {
    let error = validate_depletion(curve);
    if !error.is_empty() {
        return error;
    }
    if !water.is_finite() || water < 0.0 || !energy.is_finite() || !maximum.is_finite() || maximum < 0.0 {
        return b"Invalid initial lumped snow inventory.".to_vec();
    }
    let mut next = LumpedState::default();
    next.water = water;
    next.maximum = if water > 0.0 { cmax(water, maximum) } else { 0.0 };
    next.cover = if water > 0.0 { depletion_cover(curve, water / next.maximum) } else { 0.0 };
    if water > 0.0 && next.cover <= 0.0 {
        return b"Snow depletion curve gives no covered area for a nonempty pack.".to_vec();
    }
    next.point.water = if next.cover > 0.0 { water / next.cover } else { 0.0 };
    next.point.energy = energy;
    *state = next;
    Bytes::new()
}

// ---------------------------------------------------------------------------
// The run-time step (6 November 2026): advanceLumped and its file-static helpers.
// ---------------------------------------------------------------------------

const RHO: f64 = 1000.0;
const FUSION: f64 = 333.5;
const ICE_HEAT: f64 = 2.09;

fn soil_capacity(p: &Parameters) -> f64 {
    p.soil_density * p.active_soil_depth * p.soil_heat_capacity
}

fn represented_energy(s: &LumpedState) -> f64 {
    (if s.cover > 0.0 { s.cover } else { 1.0 }) * s.point.energy
}

fn remap(point: &mut State, water: f64, p: &Parameters) {
    let old_water = point.water;
    if point.energy <= 0.0 {
        let temperature = point.energy / (soil_capacity(p) + RHO * old_water * ICE_HEAT);
        point.energy = temperature * (soil_capacity(p) + RHO * water * ICE_HEAT);
    } else if old_water > 0.0 {
        let liquid = crate::core::cmath::cmin(1.0, point.energy / (RHO * FUSION * old_water));
        let warm = crate::core::cmath::cmax(0.0, point.energy - RHO * FUSION * old_water) / (soil_capacity(p) + RHO * old_water * ICE_HEAT);
        point.energy = liquid * RHO * FUSION * water + warm * (soil_capacity(p) + RHO * water * ICE_HEAT);
    }
    point.water = water;
}

/// advanceLumped: the point model on the covered fraction, the depletion geometry after it.
pub fn advance_lumped(state: &mut LumpedState, p: &Parameters, site: &Site, f: &Forcing, curve: &[f64; 10], output: &mut LumpedFluxes) -> Bytes {
    let error = validate_depletion(curve);
    if !error.is_empty() {
        return error;
    }
    if !state.water.is_finite() || state.water < 0.0 || !state.maximum.is_finite() || state.maximum < state.water || !state.cover.is_finite() || state.cover < 0.0 || state.cover > 1.0 || (state.water > 0.0 && state.cover <= 0.0) || (state.point.water * state.cover - state.water).abs() > 1e-12 * crate::core::cmath::cmax(1.0, state.water) {
        return b"Inconsistent lumped snow inventory.".to_vec();
    }
    let mut next = state.clone();
    let mut result = LumpedFluxes::default();
    let error = crate::hydrology::ueb_snow::advance(&mut next.point, p, site, f, &mut result.point);
    if !error.is_empty() {
        return error;
    }
    let change = next.point.water - state.point.water;
    let fresh_growth = result.point.snow > 0.0 && change > 0.0;
    let fraction = if fresh_growth { 1.0 } else { state.cover };
    result.active_fraction = fraction;
    result.rain = result.point.rain;
    result.snow = result.point.snow;
    result.liquid = fraction * result.point.liquid;
    result.sublimation = fraction * result.point.sublimation;
    result.throughfall = (1.0 - fraction) * (result.rain + result.snow);
    next.water = if state.cover == 1.0 { next.point.water } else { state.water + fraction * change };
    if next.water < -1e-12 {
        return b"Lumped snow step overdraws its water inventory.".to_vec();
    }
    next.water = crate::core::cmath::cmax(0.0, next.water);
    next.maximum = if next.water > 0.0 { crate::core::cmath::cmax(state.maximum, next.water) } else { 0.0 };
    next.cover = if next.water > 0.0 { depletion_cover(curve, next.water / next.maximum) } else { 0.0 };
    if next.water > 0.0 && next.cover <= 0.0 {
        return b"Snow depletion curve gives no covered area for a nonempty pack.".to_vec();
    }
    let target_water = if next.cover > 0.0 { next.water / next.cover } else { 0.0 };
    let energy_fraction = if state.cover > 0.0 { fraction } else { 1.0 };
    let advanced_energy = represented_energy(state) + energy_fraction * result.point.net_energy;
    if target_water != next.point.water {
        remap(&mut next.point, target_water, p);
    }
    result.net_energy = energy_fraction * result.point.net_energy;
    result.geometry_energy = represented_energy(&next) - advanced_energy;
    result.energy_residual = represented_energy(state) + result.net_energy + result.geometry_energy - represented_energy(&next);
    result.water_residual = state.water + result.rain + result.snow - result.liquid - result.throughfall - result.sublimation - next.water;
    if !result.geometry_energy.is_finite() || result.water_residual.abs() > 1e-10 * crate::core::cmath::cmax(1.0, state.water) {
        return b"Lumped snow step failed its conservation checks.".to_vec();
    }
    *state = next;
    *output = result;
    Vec::new()
}
