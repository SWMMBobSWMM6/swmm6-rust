// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/Climate.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The climate types of namespace openswmm::climate: the evaporation method, the
// seven-day temperature history and the daily climate state. The header defines no
// inline functions: MovingAvg7's push/avg_temp/avg_range and the free functions
// (hargreaves ... batchDistributeEvap) are declared only, named in comments below and
// translated with Climate.cpp. C arrays (`double ta[7]`) are fixed Rust arrays. Member
// initializers carry into a hand-written Default for ClimateState; MovingAvg7's are all
// zero, so its Default is derived.

#![allow(clippy::all)]
use super::radiation::RadiationState;

// ============================================================================
// Constants
// ============================================================================

pub const MM_PER_INCH: f64 = 25.40;

// ============================================================================
// Evaporation method
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum EvapMethod {
    Constant = 0,
    Monthly = 1,
    Timeseries = 2,
    /// Hargreaves method
    Temperature = 3,
    Pan = 4,
}

// ============================================================================
// 7-day temperature moving average (for Hargreaves ET)
// ============================================================================

#[derive(Debug, Clone, Copy, Default)]
pub struct MovingAvg7 {
    /// Daily average temps (deg F)
    pub ta: [f64; 7],
    /// Daily temp ranges (deg F)
    pub tr: [f64; 7],
    /// Circular buffer write index
    pub front: i32,
    /// Number of values stored (max 7)
    pub count: i32,
}

// Declared in Climate.hpp as members of MovingAvg7, defined in Climate.cpp (translated
// later):
//   void MovingAvg7::push(double t_avg, double t_range);
//   double MovingAvg7::avg_temp() const;
//   double MovingAvg7::avg_range() const;

// ============================================================================
// Daily climate state (scalar — broadcast to all subcatchments)
// ============================================================================

#[derive(Debug, Clone, Copy)]
pub struct ClimateState {
    pub radiation: RadiationState,
    /// Air temperature (deg F) — effective (post-forcing)
    pub temperature: f64,
    /// Daily temperature range (deg F)
    pub temp_range: f64,
    /// Evaporation rate (ft/sec)
    pub evap_rate: f64,
    /// Wind speed (mph) — effective (post-forcing)
    pub wind_speed: f64,
    /// Relative humidity (%)
    pub humidity: f64,

    // Source/default bases (pre-forcing). The per-step source lookup writes these;
    // forcing then resolves effective values on top.
    /// Source/default air temperature (deg F)
    pub temperature_src: f64,
    /// Last daily min read from a climate file (deg F)
    pub file_tmin: f64,
    /// Last daily max read from a climate file (deg F)
    pub file_tmax: f64,
    /// Source/default wind speed (mph)
    pub wind_speed_src: f64,

    // Derived values
    /// Psychrometric constant
    pub gamma: f64,
    /// Saturation vapor pressure
    pub ea: f64,

    // Hargreaves parameters
    /// Latitude (degrees)
    pub latitude: f64,

    // Site elevation for psychrometric constant (matching legacy Temp.elev)
    /// Site elevation above sea level (ft)
    pub elev: f64,
    /// Rain/snow dividing temperature (degF), legacy Snow.snotmp
    pub snow_divt: f64,

    // Monthly evaporation table (for MONTHLY method)
    pub monthly_evap: [f64; 12],

    // EVAPRATE unit-conversion factor (in/day US or mm/day SI → ft/sec) for the
    // project's unit system. Set at init; defaults to the US factor.
    pub evaprate_ucf: f64,

    // Monthly adjustment factors
    /// ADDED to the rate (ft/s), legacy Adjust.evap /= UCF(EVAPRATE)
    pub adjust_evap: [f64; 12],
    // Legacy climate.c NextEvapDate / NextEvapRate: the date the evaporation source next
    // changes and the rate that entry carries. A time-series evaporation is STEP-WISE.
    pub next_evap_date: f64,
    /// user units
    pub next_evap_rate: f64,
    /// index of the last series entry read (legacy table_getNextEntry cursor)
    pub evap_ts_pos: i32,
    pub adjust_temp: [f64; 12],
    pub adjust_rain: [f64; 12],
    /// Infiltration conductivity multipliers
    pub adjust_hydcon: [f64; 12],

    pub evap_method: EvapMethod,

    // Monthly adjustment for infiltration capacity (scales f0/fmin/Ks)
    pub infil_factor: f64,

    // 7-day moving average for Hargreaves ET (matching legacy TMovAve Tma)
    pub temp_ma: MovingAvg7,

    // Recovery factor for infiltration models
    pub recovery_factor: f64,

    // Resolved timeseries/pattern indices (set at init, -1 = none)
    /// Temperature timeseries table index
    pub temp_ts_index: i32,
    /// Evaporation timeseries table index
    pub evap_ts_index: i32,
    /// Humidity / dew-point timeseries table index
    pub humidity_ts_index: i32,
    /// Recovery pattern index in ctx.patterns
    pub recovery_pat_index: i32,

    // Sub-daily sinusoidal temperature interpolation (Gap #9). Valid only when
    // temperature comes from a daily min/max source.
    /// Current day's minimum temperature (deg F)
    pub tmin_daily: f64,
    /// Current day's maximum temperature (deg F)
    pub tmax_daily: f64,
    /// Previous day's maximum temperature (deg F)
    pub prev_tmax: f64,
    /// Sunrise hour (time of min temp), 0-24
    pub hrsr: f64,
    /// Sunset-3 hour (time of max temp), 0-24
    pub hrss: f64,
    /// Mid-hour between hrsr and hrss
    pub hrday: f64,
    /// hrsr - hrss (negative)
    pub dhrdy: f64,
    /// 24 + hrsr - hrss (hours from max to next min)
    pub dydif: f64,
    /// Longitude correction (hrs); 0 = solar time
    pub dtlong: f64,
    /// True when tmin_daily/tmax_daily are valid
    pub has_minmax: bool,
    /// Day-of-year when hrsr/hrss were last computed
    pub last_temp_doy: i32,
}

impl Default for ClimateState {
    fn default() -> Self {
        ClimateState {
            radiation: RadiationState::default(),
            temperature: 70.0,
            temp_range: 0.0,
            evap_rate: 0.0,
            wind_speed: 0.0,
            humidity: 50.0,
            temperature_src: 70.0,
            file_tmin: 70.0,
            file_tmax: 70.0,
            wind_speed_src: 0.0,
            gamma: 0.0,
            ea: 0.0,
            latitude: 0.0,
            elev: 0.0,
            snow_divt: 34.0,
            monthly_evap: [0.0; 12],
            evaprate_ucf: 1036800.0,
            adjust_evap: [0.0; 12],
            next_evap_date: 0.0,
            next_evap_rate: 0.0,
            evap_ts_pos: -1,
            adjust_temp: [0.0; 12],
            adjust_rain: [1.0; 12],
            adjust_hydcon: [1.0; 12],
            evap_method: EvapMethod::Constant,
            infil_factor: 1.0,
            temp_ma: MovingAvg7::default(),
            recovery_factor: 1.0,
            temp_ts_index: -1,
            evap_ts_index: -1,
            humidity_ts_index: -1,
            recovery_pat_index: -1,
            tmin_daily: 0.0,
            tmax_daily: 0.0,
            prev_tmax: 0.0,
            hrsr: 6.0,
            hrss: 15.0,
            hrday: 10.5,
            dhrdy: -9.0,
            dydif: 15.0,
            dtlong: 0.0,
            has_minmax: false,
            last_temp_doy: -1,
        }
    }
}

// ============================================================================
// Functions
// ============================================================================

// Declared in Climate.hpp, defined in Climate.cpp (translated later). SimulationContext
// is forward-declared there:
//   double hargreaves(double latitude, int day_of_year, double t_avg, double t_range);
//   void updateDailyClimate(ClimateState& state, int day_of_year, int month);
//   void updateTempTimes(ClimateState& state, int doy);
//   double getSubdailyTemp(const ClimateState& state, double hour);
//   void batchDistributeEvap(double evap_rate, const double* ponded_depth, double* evap_out,
//                            int n, double dt);
