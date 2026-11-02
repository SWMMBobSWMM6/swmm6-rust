// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/surface/UebSnow.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Only Parameters::validate() and the three predicates it uses (anonymous namespace):
// the part of UebSnow.cpp that input/PostParseResolver.cpp reaches. The rest of the UEB
// point model is translated with the snow solver.

use crate::core::options::Bytes;
use crate::hydrology::ueb_types::{Parameters, TemperatureMethod};

fn finite(v: f64) -> bool {
    v.is_finite()
}

fn nonnegative(v: f64) -> bool {
    finite(v) && v >= 0.0
}

fn fraction(v: f64) -> bool {
    nonnegative(v) && v <= 1.0
}

impl Parameters {
    /// Empty when every parameter is finite and in range, else the reason.
    pub fn validate(&self) -> Bytes {
        if self.temperature != TemperatureMethod::ModifiedForceRestore && self.temperature != TemperatureMethod::Equilibrium {
            return b"Unknown UEB surface temperature method.".to_vec();
        }
        for x in [
            self.rain_temperature,
            self.snow_temperature,
            self.emissivity,
            self.roughness,
            self.measurement_height,
            self.snow_density,
            self.soil_density,
            self.soil_heat_capacity,
            self.active_soil_depth,
            self.liquid_retention,
            self.drainage_conductivity,
            self.equilibrium_conductance,
            self.snow_conductivity,
            self.soil_conductivity,
            self.low_frequency,
            self.amplitude,
            self.stability,
            self.ground_albedo,
            self.visible_albedo,
            self.infrared_albedo,
            self.new_snow_depth,
            self.age_dust,
            self.visible_age,
            self.infrared_age,
            self.max_step,
        ] {
            if !finite(x) {
                return b"Nonfinite UEB parameter.".to_vec();
            }
        }
        if self.rain_temperature <= self.snow_temperature
            || self.snow_temperature <= -273.15
            || !fraction(self.emissivity)
            || self.roughness <= 0.0
            || self.measurement_height <= self.roughness
            || self.snow_density <= 0.0
            || self.snow_density >= 917.0
            || self.soil_density <= 0.0
            || self.soil_heat_capacity <= 0.0
            || self.active_soil_depth <= 0.0
            || self.liquid_retention < 0.0
            || self.liquid_retention >= 1000.0 / self.snow_density - 1000.0 / 917.0
            || self.drainage_conductivity < 0.0
            || self.equilibrium_conductance <= 0.0
            || self.snow_conductivity <= 0.0
            || self.soil_conductivity <= 0.0
            || self.low_frequency <= 0.0
            || self.amplitude <= 0.0
            || !fraction(self.stability)
            || !fraction(self.ground_albedo)
            || !fraction(self.visible_albedo)
            || !fraction(self.infrared_albedo)
            || self.new_snow_depth < 0.0
            || self.age_dust < 0.0
            || !fraction(self.visible_age)
            || !fraction(self.infrared_age)
            || self.max_step < 1.0
            || self.max_step > 21600.0
        {
            return b"Invalid UEB parameter range.".to_vec();
        }
        Bytes::new()
    }
}
