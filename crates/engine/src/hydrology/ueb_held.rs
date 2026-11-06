// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/surface/UebHeld.cpp and UebHeld.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin,
// licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 6 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// A held UEB interval: prepared once over the snow step, consumed by the runoff steps inside
// it (interpolated inventory, prorated fluxes). Errors are Bytes, empty on success.

use crate::core::cmath::cmax;
use crate::core::options::Bytes;
use crate::hydrology::ueb_lumped::{advance_lumped, depletion_cover, validate_depletion};
use crate::hydrology::ueb_types::{Fluxes, Forcing, HeldState, LumpedFluxes, LumpedState, Parameters, Site, State, ThermalSample};

fn energy(s: &LumpedState) -> f64 {
    (if s.cover > 0.0 { s.cover } else { 1.0 }) * s.point.energy
}

fn blend(a: f64, b: f64, fraction: f64) -> f64 {
    a + fraction * (b - a)
}

fn history_at(state: &mut State, start: &State, elapsed: f64, f: &Fluxes) {
    state.history = start.history.clone();
    if elapsed > 0.0 {
        state.history.push(ThermalSample { seconds: elapsed, surface: f.surface_temperature, bulk: f.bulk_temperature });
    }
    let mut total = 0.0;
    for h in state.history.iter() {
        total += h.seconds;
    }
    let mut excess = cmax(0.0, total - 86400.0);
    let mut remove = 0usize;
    while remove < state.history.len() && excess >= state.history[remove].seconds {
        excess -= state.history[remove].seconds;
        remove += 1;
    }
    state.history.drain(0..remove);
    if !state.history.is_empty() {
        state.history[0].seconds -= excess;
    }
}

/// prepareHeld: advance the accepted inventory over the whole interval, held for consumption.
pub fn prepare_held(s: &mut HeldState, p: &Parameters, site: &Site, f: &Forcing, curve: &[f64; 10]) -> Bytes {
    if s.pending() {
        return b"The previous snow interval has not been consumed.".to_vec();
    }
    let mut end = s.accepted.clone();
    let mut flux = LumpedFluxes::default();
    let error = advance_lumped(&mut end, p, site, f, curve, &mut flux);
    if !error.is_empty() {
        return error;
    }
    s.begin = s.accepted.clone();
    s.end = end;
    s.flux = flux;
    s.duration = f.seconds;
    s.elapsed = 0.0;
    Vec::new()
}

/// consumeHeld: `seconds` of the prepared interval.
pub fn consume_held(s: &mut HeldState, seconds: f64, curve: &[f64; 10], output: &mut LumpedFluxes) -> Bytes {
    const CLOCK_TOLERANCE: f64 = 1e-5;
    if !seconds.is_finite() || seconds <= 0.0 || !s.pending() || seconds > s.duration - s.elapsed + CLOCK_TOLERANCE {
        return b"Snow consumption must lie inside its prepared interval.".to_vec();
    }
    let error = validate_depletion(curve);
    if !error.is_empty() {
        return error;
    }
    let end = if (s.duration - s.elapsed - seconds).abs() <= CLOCK_TOLERANCE { s.duration } else { s.elapsed + seconds };
    let share = (end - s.elapsed) / s.duration;
    let progress = end / s.duration;
    let mut next = s.accepted.clone();
    let mut result = s.flux;
    result.rain *= share;
    result.snow *= share;
    result.liquid *= share;
    result.throughfall *= share;
    result.sublimation *= share;
    result.net_energy *= share;
    if end == s.duration {
        next = s.end.clone();
    } else {
        next.water = blend(s.begin.water, s.end.water, progress);
        next.maximum = if next.water > 0.0 { cmax(s.begin.maximum, next.water) } else { 0.0 };
        next.cover = if next.water > 0.0 { depletion_cover(curve, next.water / next.maximum) } else { 0.0 };
        if next.water > 0.0 && next.cover <= 0.0 {
            return b"Snow depletion curve gives no covered area for a nonempty pack.".to_vec();
        }
        next.point.water = if next.cover > 0.0 { next.water / next.cover } else { 0.0 };
        next.point.energy = blend(energy(&s.begin), energy(&s.end), progress) / (if next.cover > 0.0 { next.cover } else { 1.0 });
        next.point.age = blend(s.begin.point.age, s.end.point.age, progress);
        next.point.refreezing_depth = blend(s.begin.point.refreezing_depth, s.end.point.refreezing_depth, progress);
        next.point.total_refreezing_depth = blend(s.begin.point.total_refreezing_depth, s.end.point.total_refreezing_depth, progress);
        let begin_point = s.begin.point.clone();
        history_at(&mut next.point, &begin_point, end, &s.flux.point);
    }
    result.geometry_energy = energy(&next) - energy(&s.accepted) - result.net_energy;
    result.energy_residual = energy(&s.accepted) + result.net_energy + result.geometry_energy - energy(&next);
    result.water_residual = s.accepted.water + result.rain + result.snow - result.liquid - result.throughfall - result.sublimation - next.water;
    if !next.point.energy.is_finite() || result.water_residual.abs() > 1e-10 * cmax(1.0, next.water) {
        return b"Held snow consumption failed its conservation checks.".to_vec();
    }
    if s.elapsed > 0.0 {
        result.point.steps = 0;
        result.point.retries = 0;
        result.point.diagnostics = 0;
    }
    s.accepted = next;
    s.elapsed = end;
    *output = result;
    Vec::new()
}

/// closeHeld: abandon the prepared interval.
pub fn close_held(s: &mut HeldState) {
    s.duration = 0.0;
    s.elapsed = 0.0;
}
