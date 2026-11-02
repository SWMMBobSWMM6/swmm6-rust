// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/Gage.hpp and Gage.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The parts of the rain-gage module translated so far: gageIsUsed (SWMMEngine::open
// reaches it), and the constants, RainType and GageState, which SWMMEngine::initHydrology
// sizes (gage_states_). The rest (rainfall state machine, conversion: convertRainfall,
// gageNextRateChange, ...) is translated with the hydrology step functions (phase 1c).
// `double past_rain[MAXPASTRAIN]` is a fixed array; member initializers carry into a
// hand-written Default for GageState.
//
// 1D build: the 2D-mesh clause (every gage is used under a rain-on-grid RAINFALL_MODE)
// tests ctx.twod_io.mesh, which a 1D build never sets, so it is always false here.

use crate::core::context::SimulationContext;

/// One second in days
pub const ONE_SECOND: f64 = 1.1574074e-5;
/// Max past hours tracked
pub const MAXPASTRAIN: usize = 48;

/// Rain type codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum RainType {
    Intensity = 0,
    Volume = 1,
    Cumulative = 2,
}

/// Per-gage state.
#[derive(Debug, Clone, Copy)]
pub struct GageState {
    /// Current intensity (project in/hr or mm/hr)
    pub rainfall: f64,
    /// Current snow intensity (project in/hr or mm/hr)
    pub snowfall: f64,
    /// rainfall + snowfall
    pub total_precip: f64,
    /// Cumulative rain accumulator (for CUMULATIVE type)
    pub rain_accum: f64,
    /// API-overridden rainfall (-1 = no override)
    pub api_rainfall: f64,
    /// Snow catch factor
    pub snow_factor: f64,
    /// Rainfall scaling factor (1.0 = no scaling)
    pub scale_factor: f64,
    /// Unit conversion factor
    pub units_factor: f64,
    /// Monthly/seasonal adjustment factor
    pub adjust_factor: f64,
    /// Rain recording interval (seconds)
    pub rain_interval: f64,
    pub rain_type: RainType,

    /// Past hourly rainfall totals
    pub past_rain: [f64; MAXPASTRAIN],
    /// Accumulator for current hour
    pub past_rain_accum: f64,
    /// Seconds in the current partial hour (0..3600)
    pub past_rain_interval: f64,
}

impl Default for GageState {
    fn default() -> Self {
        GageState {
            rainfall: 0.0,
            snowfall: 0.0,
            total_precip: 0.0,
            rain_accum: 0.0,
            api_rainfall: -1.0,
            snow_factor: 1.0,
            scale_factor: 1.0,
            units_factor: 1.0,
            adjust_factor: 1.0,
            rain_interval: 0.0,
            rain_type: RainType::Intensity,
            past_rain: [0.0; MAXPASTRAIN],
            past_rain_accum: 0.0,
            past_rain_interval: 0.0,
        }
    }
}

/// True when a subcatchment or an RDII unit hydrograph reads the gage.
pub fn gage_is_used(ctx: &SimulationContext, gage_idx: i32) -> bool {
    for i in 0..ctx.n_subcatches() {
        if ctx.subcatches.gage[i as usize] == gage_idx {
            return true;
        }
    }
    for gname in &ctx.unit_hyds.gage_names {
        if ctx.gage_names.find(gname) == gage_idx {
            return true;
        }
    }
    false
}
