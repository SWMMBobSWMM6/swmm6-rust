// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/HourlyMeteorology.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The hourly meteorology record, the daily snow-weather cache row and the hourly cache
// class of namespace openswmm::climate. The class keeps its one private member
// (`records_`, here `records`) private. Its two inline accessors are here; `empty()` is
// `is_empty()`, as elsewhere in this crate. read/load/average are declared only and are
// translated with HourlyMeteorology.cpp.

#![allow(clippy::all)]

/// Date/time marks the START of the hourly interval. Temperatures are °C, humidity is
/// percent, wind m/s and radiation W/m². Every initializer is zero (`values{}`), so
/// Default is derived.
#[derive(Debug, Clone, Copy, Default)]
pub struct HourlyMetRecord {
    pub date: f64,
    /// temperature, RH, wind, Qsi, Qli
    pub values: [f64; 5],
}

/// Immutable daily source cache for independent surface clocks.
#[derive(Debug, Clone, Copy)]
pub struct SnowDailyWeather {
    pub date: f64,
    pub tmin: f64,
    pub tmax: f64,
    pub previous_max: f64,
    pub wind: f64,
}

impl Default for SnowDailyWeather {
    fn default() -> Self {
        SnowDailyWeather {
            date: 0.0,
            tmin: 70.0,
            tmax: 70.0,
            previous_max: 70.0,
            wind: 0.0,
        }
    }
}

/// Owned consecutive hourly records with interval-weighted means.
#[derive(Debug, Clone, Default)]
pub struct HourlyMeteorology {
    records: Vec<HourlyMetRecord>,
}

impl HourlyMeteorology {
    // Declared in HourlyMeteorology.hpp, defined in HourlyMeteorology.cpp (translated
    // later):
    //   std::string read(std::istream& stream);
    //   std::string load(const std::string& path);
    //   std::string average(double start, double seconds,
    //                       std::array<double, 5>& result) const;

    /// True when no validated records are stored (upstream empty()).
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Borrow the immutable hourly starts and their complete-hour weather values.
    #[inline]
    pub fn records(&self) -> &Vec<HourlyMetRecord> {
        &self.records
    }
}
