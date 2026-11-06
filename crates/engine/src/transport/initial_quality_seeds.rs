// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/InitialQualitySeeds.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin,
// licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// Per-element [INITIAL_QUALITY] seeds for the reserved age and temperature rows, read by the
// EULERIAN_ARD and LAGRANGIAN engines' init (applyInitialAgeOverrides and
// applyInitialTempOverrides are in water_age_legacy.rs and heat_legacy.rs).

use crate::core::context::SimulationContext;
use crate::data::initial_quality::InitialQualityData;

/// initialAgeSecondsFor: the last matching __WATER_AGE__ row (hours, returned in seconds), or
/// `fallback`.
pub fn initial_age_seconds_for(ctx: &SimulationContext, is_link: bool, idx: i32, fallback: f64) -> f64 {
    let mut v = fallback;
    let iq = &ctx.initial_quality;
    for r in 0..iq.count() {
        let ur = r as usize;
        if iq.kind[ur] != InitialQualityData::K_KIND_WATER_AGE {
            continue;
        }
        if (iq.is_link[ur] != 0) != is_link {
            continue;
        }
        if iq.elem_idx[ur] != idx {
            continue;
        }
        v = iq.value[ur] * 3600.0;
    }
    v
}

/// initialTempFor: the last matching __TEMPERATURE__ row (degC), or `fallback`.
pub fn initial_temp_for(ctx: &SimulationContext, is_link: bool, idx: i32, fallback: f64) -> f64 {
    let mut v = fallback;
    let iq = &ctx.initial_quality;
    for r in 0..iq.count() {
        let ur = r as usize;
        if iq.kind[ur] != InitialQualityData::K_KIND_TEMPERATURE {
            continue;
        }
        if (iq.is_link[ur] != 0) != is_link {
            continue;
        }
        if iq.elem_idx[ur] != idx {
            continue;
        }
        v = iq.value[ur];
    }
    v
}
