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
use crate::hydrology::ueb_types::LumpedState;

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
