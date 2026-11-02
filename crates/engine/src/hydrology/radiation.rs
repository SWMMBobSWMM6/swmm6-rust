// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/Radiation.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The radiation types of namespace openswmm::climate: the resolved interval radiation
// record, its origin markers, the two Bird-Hulstrom constants and the solar-position
// result. The header defines no inline functions; every function it declares
// (resolveRadiation ... updateSolarForcing) is named in a comment below and translated
// with Radiation.cpp. Upstream's `inline constexpr` constants are Rust consts
// (kSolarConstant is K_SOLAR_CONSTANT, kStdPressureMb K_STD_PRESSURE_MB). The scoped
// enum has no underlying type, so it is `int` in C++ and #[repr(i32)] here. Member
// initializers carry into hand-written Default impls.

#![allow(clippy::all)]

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum RadiationOrigin {
    Computed = 0,
    Timeseries = 1,
    MetFile = 2,
    LegacyHeat = 3,
}

#[derive(Debug, Clone, Copy)]
pub struct RadiationState {
    /// W/m², horizontal
    pub shortwave: f64,
    pub longwave: f64,
    pub cosine: f64,
    pub horizontal: f64,
    pub transmittance: f64,
    pub cloud: f64,
    /// °C
    pub temperature: f64,
    /// fraction
    pub humidity: f64,
    /// m/s
    pub wind: f64,
    pub pressure: f64,
    pub declination: f64,
    pub hour_angle: f64,
    pub latitude: f64,
    pub date: f64,
    pub seconds: f64,
    pub shortwave_source: RadiationOrigin,
    pub longwave_source: RadiationOrigin,
    pub shared_heat_forcing: bool,
}

impl Default for RadiationState {
    fn default() -> Self {
        RadiationState {
            shortwave: 0.0,
            longwave: 0.0,
            cosine: 0.0,
            horizontal: 0.0,
            transmittance: 0.0,
            cloud: 0.0,
            temperature: 0.0,
            humidity: 0.5,
            wind: 0.0,
            pressure: 101325.0,
            declination: 0.0,
            hour_angle: 0.0,
            latitude: 0.0,
            date: 0.0,
            seconds: 0.0,
            shortwave_source: RadiationOrigin::Computed,
            longwave_source: RadiationOrigin::Computed,
            shared_heat_forcing: false,
        }
    }
}

// Declared in Radiation.hpp, defined in Radiation.cpp (translated later).
// SimulationContext is forward-declared; surface::ueb::Forcing is
// crate::hydrology::ueb_types::Forcing:
//   std::string resolveRadiation(SimulationContext& ctx, double date, double seconds,
//                                RadiationState& output);
//   std::string initializeRadiation(SimulationContext& ctx);
//   std::string resolveUebForcing(SimulationContext& ctx, double date, double seconds,
//                                 double slope, double aspect, surface::ueb::Forcing& forcing,
//                                 double& pressure);

/// Solar constant, W/m² (Bird & Hulstrom 1981; WMO 1982).
pub const K_SOLAR_CONSTANT: f64 = 1367.0;

/// Sea-level standard pressure, millibars — Bird's pressure reference.
pub const K_STD_PRESSURE_MB: f64 = 1013.25;

/// Sun geometry at one instant. Angles in DEGREES. `cos_zenith` is carried alongside
/// `zenith_deg` because every downstream consumer wants the cosine.
#[derive(Debug, Clone, Copy)]
pub struct SolarPosition {
    /// Solar declination, +N
    pub declination_deg: f64,
    /// Equation of time, minutes
    pub eq_time_min: f64,
    /// Solar hour angle, 0 at solar noon
    pub hour_angle_deg: f64,
    /// Solar zenith angle
    pub zenith_deg: f64,
    /// Clamped at 0 below the horizon
    pub cos_zenith: f64,
    /// 90 − zenith (may be negative)
    pub elevation_deg: f64,
    /// Earth–sun distance correction `E0 = (r0/r)²`, dimensionless.
    pub eccentricity: f64,
}

impl Default for SolarPosition {
    fn default() -> Self {
        SolarPosition {
            declination_deg: 0.0,
            eq_time_min: 0.0,
            hour_angle_deg: 0.0,
            zenith_deg: 90.0,
            cos_zenith: 0.0,
            elevation_deg: 0.0,
            eccentricity: 1.0,
        }
    }
}

// Declared in Radiation.hpp, defined in Radiation.cpp (translated later). SolarConfig is
// forward-declared there and lives in crate::data::heat::SolarConfig:
//   SolarPosition solarPosition(int day_of_year, double hour_local, double lat_deg,
//                               double lon_deg, double tz_hours) noexcept;
//   double airMass(double zenith_deg) noexcept;
//   double pressureFromElevation(double elevation_m) noexcept;
//   double birdClearSkyGHI(const SolarPosition& pos, double pressure_mb,
//                          const SolarConfig& cfg) noexcept;
//   double cloudShortwaveFactor(double cloud_fraction, double k, double n) noexcept;
//   double cloudLongwaveFactor(double cloud_fraction, double k_lw) noexcept;
//   void updateSolarForcing(SimulationContext& ctx, double seconds = 1);  // default argument
//   void updateSolarForcing(SimulationContext& ctx, double seconds, bool mesh_enabled);
// The one-argument call updateSolarForcing(ctx) means seconds = 1.0 and, through the
// two-argument overload, mesh_enabled = false.
