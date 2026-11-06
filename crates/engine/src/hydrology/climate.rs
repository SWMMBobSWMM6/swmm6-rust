// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/Climate.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The climate functions the runoff step calls: MovingAvg7 push/avg_temp/avg_range,
// hargreaves, updateDailyClimate, updateTempTimes, getSubdailyTemp. (The types are in
// climate_types.rs.) batchDistributeEvap is not on the run path; NOT YET TRANSLATED.

use crate::core::cmath::{cmax, cmin};
use crate::core::constants::PI;
use crate::core::unit_conversion as ucf;
use crate::hydrology::climate_types::{ClimateState, EvapMethod, MovingAvg7};

impl MovingAvg7 {
    pub fn push(&mut self, t_avg: f64, t_range: f64) {
        self.ta[self.front as usize] = t_avg;
        self.tr[self.front as usize] = t_range;
        self.front = (self.front + 1) % 7;
        if self.count < 7 {
            self.count += 1;
        }
    }

    pub fn avg_temp(&self) -> f64 {
        if self.count == 0 {
            return 0.0;
        }
        let mut sum = 0.0;
        for i in 0..self.count as usize {
            sum += self.ta[i];
        }
        sum / self.count as f64
    }

    pub fn avg_range(&self) -> f64 {
        if self.count == 0 {
            return 0.0;
        }
        let mut sum = 0.0;
        for i in 0..self.count as usize {
            sum += self.tr[i];
        }
        sum / self.count as f64
    }
}

/// Hargreaves evapotranspiration (mm/day).
pub fn hargreaves(latitude: f64, day_of_year: i32, t_avg: f64, t_range: f64) -> f64 {
    let a = 2.0 * PI / 365.0;
    let doy = day_of_year as f64;

    // Convert to Celsius
    let ta = (t_avg - 32.0) * 5.0 / 9.0;
    let mut tr = t_range * 5.0 / 9.0;
    tr = cmax(tr, 0.0);

    // Latent heat of vaporization
    let mut lambda = 2.50 - 0.002361 * ta;
    lambda = cmax(lambda, 0.01);

    // Relative earth-sun distance
    let dr = 1.0 + 0.033 * (a * doy).cos(); // LIBM

    // Solar declination
    let delta = 0.4093 * (a * (284.0 + doy)).sin(); // LIBM

    // Latitude in radians
    let phi = latitude * 2.0 * PI / 360.0;

    // Sunset hour angle
    let mut x = -phi.tan() * delta.tan(); // LIBM
    x = cmax(-1.0, cmin(1.0, x));
    let omega = x.acos(); // LIBM

    // Extraterrestrial radiation (MJ/m2/day)
    let mut ra = 37.6 * dr * (omega * phi.sin() * delta.sin() + phi.cos() * delta.cos() * omega.sin()); // LIBM
    ra = cmax(ra, 0.0);

    // Hargreaves evapotranspiration (mm/day)
    let mut e = 0.0023 * ra / lambda * tr.sqrt() * (ta + 17.8);
    e = cmax(e, 0.0);
    e
}

/// updateDailyClimate: the evaporation rate of the climate source, its monthly adjustment,
/// the saturation vapour pressure and the psychrometric constant.
pub fn update_daily_climate(state: &mut ClimateState, day_of_year: i32, month: i32) {
    let m = month as usize;
    match state.evap_method {
        EvapMethod::Constant => {
            state.evap_rate = state.monthly_evap[0] / state.evaprate_ucf;
        }
        EvapMethod::Monthly => {
            state.evap_rate = state.monthly_evap[m] / state.evaprate_ucf;
        }
        EvapMethod::Temperature => {
            state.temp_ma.push(state.temperature, state.temp_range);
            let e_mm = hargreaves(state.latitude, day_of_year, state.temp_ma.avg_temp(), state.temp_ma.avg_range());
            state.evap_rate = e_mm / ucf::UCF[ucf::EVAPRATE as usize][1];
        }
        EvapMethod::Timeseries | EvapMethod::Pan => {}
    }

    if state.evap_method != EvapMethod::Timeseries && state.evap_method != EvapMethod::Pan {
        state.evap_rate += state.adjust_evap[m];
    }

    // Saturation vapor pressure (for snowmelt rain-on-snow)
    let ta = state.temperature;
    state.ea = 8.1175e6 * (-7701.544 / (ta + 405.0265)).exp(); // LIBM

    // Psychrometric constant (legacy climate.c)
    {
        let z = state.elev / 1000.0;
        let pa = if z <= 0.0 { 29.9 } else { 29.9 - 1.02 * z + 0.0032 * z.powf(2.4) }; // LIBM
        state.gamma = 0.000359 * pa;
    }
}

/// updateTempTimes: sunrise/sunset hours for the sub-daily temperature curve (legacy
/// climate.c updateTempTimes).
pub fn update_temp_times(state: &mut ClimateState, doy: i32) {
    let decl = 0.40928 * (0.017202 * (172.0 - doy as f64)).cos(); // LIBM

    let lat_rad = state.latitude * PI / 180.0;
    let tan_lat = lat_rad.tan(); // LIBM
    let arg = -decl.tan() * tan_lat; // LIBM

    let hrsr;
    let hrss;
    if arg >= 1.0 {
        hrsr = 12.0;
        hrss = 12.0;
    } else if arg <= -1.0 {
        hrsr = 0.0;
        hrss = 24.0;
    } else {
        let hrang = 3.8197 * arg.acos(); // LIBM
        hrsr = 12.0 - hrang + state.dtlong;
        hrss = 12.0 + hrang + state.dtlong - 3.0;
    }

    state.hrsr = hrsr;
    state.hrss = hrss;
    state.hrday = (hrsr + hrss) / 2.0;
    state.dhrdy = hrsr - hrss;
    state.dydif = 24.0 + hrsr - hrss;
}

/// getSubdailyTemp: the sinusoidal temperature at `hour` (legacy climate.c setTemp).
pub fn get_subdaily_temp(state: &ClimateState, hour: f64) -> f64 {
    let tmin = state.tmin_daily;
    let tmax = state.tmax_daily;
    let tave = (tmin + tmax) / 2.0;
    let trng = (tmax - tmin) / 2.0;
    let trng1 = state.prev_tmax - tmin;

    let hrsr = state.hrsr;
    let hrss = state.hrss;
    let hrday = state.hrday;
    let dhrdy = state.dhrdy;
    let dydif = state.dydif;

    if hour < hrsr {
        let denom = if dydif.abs() > 1e-10 { dydif } else { 1.0 };
        tmin + (trng1 / 2.0) * (PI / denom * (hrsr - hour)).sin() // LIBM
    } else if hour <= hrss {
        let denom = if dhrdy.abs() > 1e-10 { dhrdy } else { -1.0 };
        tave + trng * (PI / denom * (hrday - hour)).sin() // LIBM
    } else {
        let denom = if dydif.abs() > 1e-10 { dydif } else { 1.0 };
        tmax - trng * (PI / denom * (hour - hrss)).sin() // LIBM
    }
}
