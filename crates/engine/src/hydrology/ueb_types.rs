// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/surface/UebSnow.hpp, UebLumped.hpp,
// UebHeld.hpp, UebConfig.hpp, MeshSnowConfig.hpp and UebSubcatchments.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The Utah Energy Balance (UEB) snow types: point parameters/site/forcing/state/fluxes,
// the lumped and held wrappers, the snow and mesh-snow configuration, and the native
// per-subcatchment state. Upstream spreads these over two namespaces,
// openswmm::surface::ueb (UebSnow/UebLumped/UebHeld/UebSubcatchments) and
// openswmm::surface (UebConfig/MeshSnowConfig); no names collide, so they share this one
// module and keep upstream's type names. Nested types take the outer name as a prefix
// (RadiationSource::Mode is RadiationSourceMode, MeshSnowConfig::Mode is
// MeshSnowConfigMode). Scoped enums without an underlying type are `int` in C++, hence
// #[repr(i32)]. Member initializers carry into Default, by hand wherever one differs
// from Rust's default; SnowKind has no Default because its two owners default it
// differently (SnowpackConfig to Ueb, SnowConfig to DegreeDay).
//
// std::map becomes BTreeMap. std::map<std::string, _> orders keys by
// std::char_traits<char>::compare, which compares as unsigned char, the same order as
// Vec<u8>'s Ord, so BTreeMap<Bytes, _> iterates in upstream's order. std::optional<double>
// is Option<f64>; std::array<double, N> is [f64; N].
//
// Only the inline functions are here (HeldState::pending, SnowConfig::model_for,
// MeshSnowConfig::assigned). Every function the headers declare without defining is
// named in a comment where upstream declares it; they are translated with the .cpp files.

#![allow(clippy::all)]
use crate::core::options::{Bytes, FilePathPair};
use std::collections::BTreeMap;

// ===========================================================================
// UebSnow.hpp (namespace openswmm::surface::ueb)
// ===========================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum TemperatureMethod {
    ModifiedForceRestore = 0,
    Equilibrium = 1,
}

/// Point UEB material, radiation and integration parameters. Temperatures are Celsius,
/// lengths m, densities kg/m3, heat capacities kJ/(kg K), conductivities kJ/(m K h),
/// drainage/conductance m/h, max_step seconds.
#[derive(Debug, Clone, Copy)]
pub struct Parameters {
    pub temperature: TemperatureMethod,
    pub rain_temperature: f64,
    pub snow_temperature: f64,
    pub emissivity: f64,
    pub roughness: f64,
    pub measurement_height: f64,
    pub snow_density: f64,
    pub soil_density: f64,
    pub soil_heat_capacity: f64,
    pub active_soil_depth: f64,
    /// reference parameter units
    pub liquid_retention: f64,
    /// m/h, reference parameter units
    pub drainage_conductivity: f64,
    /// m/h
    pub equilibrium_conductance: f64,
    /// kJ/(m K h)
    pub snow_conductivity: f64,
    /// kJ/(m K h)
    pub soil_conductivity: f64,
    pub low_frequency: f64,
    pub amplitude: f64,
    pub stability: f64,
    pub ground_albedo: f64,
    pub visible_albedo: f64,
    pub infrared_albedo: f64,
    pub new_snow_depth: f64,
    pub age_dust: f64,
    pub visible_age: f64,
    pub infrared_age: f64,
    /// seconds; kernel subdivides longer forcing intervals
    pub max_step: f64,
}

impl Default for Parameters {
    fn default() -> Self {
        Parameters {
            temperature: TemperatureMethod::ModifiedForceRestore,
            rain_temperature: 3.0,
            snow_temperature: -1.0,
            emissivity: 0.99,
            roughness: 0.01,
            measurement_height: 2.0,
            snow_density: 200.0,
            soil_density: 1700.0,
            soil_heat_capacity: 2.09,
            active_soil_depth: 0.1,
            liquid_retention: 0.02,
            drainage_conductivity: 200.0,
            equilibrium_conductance: 0.02,
            snow_conductivity: 0.33,
            soil_conductivity: 6.5,
            low_frequency: 0.0654,
            amplitude: 1.0,
            stability: 1.0,
            ground_albedo: 0.25,
            visible_albedo: 0.85,
            infrared_albedo: 0.65,
            new_snow_depth: 0.002,
            age_dust: 0.03,
            visible_age: 0.2,
            infrared_age: 0.5,
            max_step: 3600.0,
        }
    }
}

// Declared in UebSnow.hpp, defined in UebSnow.cpp (translated later):
//   std::string Parameters::validate() const;

/// Per-column environmental settings: forest canopy fraction, drift snow multiplier,
/// pressure Pa, ground_heat W/m2 into the snow system, extinction_depth m.
#[derive(Debug, Clone, Copy)]
pub struct Site {
    pub forest: f64,
    pub drift: f64,
    pub pressure: f64,
    pub ground_heat: f64,
    pub extinction_depth: f64,
}

impl Default for Site {
    fn default() -> Self {
        Site {
            forest: 0.0,
            drift: 1.0,
            pressure: 101325.0,
            ground_heat: 0.0,
            extinction_depth: 0.1,
        }
    }
}

/// Held point meteorology over one requested interval. snow_fraction = -1 selects the
/// UEB temperature partition; [0,1] prescribes it.
#[derive(Debug, Clone, Copy)]
pub struct Forcing {
    pub seconds: f64,
    pub precipitation: f64,
    pub temperature: f64,
    pub wind: f64,
    pub humidity: f64,
    /// W/m2
    pub shortwave: f64,
    /// W/m2
    pub longwave: f64,
    /// cosine of incident sun angle
    pub cosine: f64,
    pub snow_catch: f64,
    /// -1 uses UEB temperature partition; [0,1] prescribes it
    pub snow_fraction: f64,
}

impl Default for Forcing {
    fn default() -> Self {
        Forcing {
            seconds: 0.0,
            precipitation: 0.0,
            temperature: 0.0,
            wind: 0.0,
            humidity: 0.5,
            shortwave: 0.0,
            longwave: 0.0,
            cosine: 0.0,
            snow_catch: 1.0,
            snow_fraction: -1.0,
        }
    }
}

/// One duration-weighted preceding thermal sample (seconds and Celsius). Every
/// initializer is zero, so Default is derived.
#[derive(Debug, Clone, Copy, Default)]
pub struct ThermalSample {
    pub seconds: f64,
    pub surface: f64,
    pub bulk: f64,
}

/// One owned point snow/active-soil energy and thermal-memory inventory. Every
/// initializer is zero (history starts empty), so Default is derived.
#[derive(Debug, Clone, Default)]
pub struct State {
    pub energy: f64,
    pub water: f64,
    pub age: f64,
    pub refreezing_depth: f64,
    pub total_refreezing_depth: f64,
    /// at most the preceding 24 hours, duration weighted
    pub history: Vec<ThermalSample>,
}

/// Integrated accepted point exchanges and final thermal diagnostics. Every initializer
/// is zero, so Default is derived.
#[derive(Debug, Clone, Copy, Default)]
pub struct Fluxes {
    pub rain: f64,
    pub snow: f64,
    pub liquid: f64,
    /// integrated m; sublimation may be negative
    pub sublimation: f64,
    /// integrated kJ/m2
    pub net_energy: f64,
    /// integrated kJ/m2
    pub outgoing_melt_energy: f64,
    pub surface_temperature: f64,
    pub bulk_temperature: f64,
    pub albedo: f64,
    pub water_residual: f64,
    pub energy_residual: f64,
    pub steps: i32,
    pub retries: i32,
    pub diagnostics: i32,
}

// Declared in UebSnow.hpp, defined in UebSnow.cpp (translated later):
//   std::string advance(State& state, const Parameters& parameters, const Site& site,
//                       const Forcing& forcing, Fluxes& output);
//   double bulkTemperature(const State& state, const Parameters& parameters) noexcept;

// ===========================================================================
// UebLumped.hpp (namespace openswmm::surface::ueb)
// ===========================================================================

/// The point state describes the snow-covered part; water is the average over the
/// entire element. Every initializer is zero, so Default is derived.
#[derive(Debug, Clone, Default)]
pub struct LumpedState {
    pub point: State,
    pub water: f64,
    pub maximum: f64,
    pub cover: f64,
}

/// Every initializer is zero, so Default is derived.
#[derive(Debug, Clone, Copy, Default)]
pub struct LumpedFluxes {
    pub point: Fluxes,
    pub rain: f64,
    pub snow: f64,
    pub liquid: f64,
    pub throughfall: f64,
    pub sublimation: f64,
    pub active_fraction: f64,
    pub water_residual: f64,
    // UEB's depletion adjustment preserves temperature/liquid fraction. Its change in
    // represented energy is reported separately from physical forcing; it must never be
    // disguised as numerical energy closure.
    pub net_energy: f64,
    pub geometry_energy: f64,
    pub energy_residual: f64,
}

// Declared in UebLumped.hpp, defined in UebLumped.cpp (translated later); `curve` is
// std::array<double, 10>, i.e. &[f64; 10]:
//   std::string validateDepletion(const std::array<double, 10>& curve);
//   double depletionCover(const std::array<double, 10>& curve, double relative_water) noexcept;
//   std::string initializeLumped(LumpedState& state, double water, double point_energy,
//                                double maximum, const std::array<double, 10>& curve);
//   std::string advanceLumped(LumpedState& state, const Parameters& parameters,
//                             const Site& site, const Forcing& forcing,
//                             const std::array<double, 10>& curve, LumpedFluxes& output);

// ===========================================================================
// UebHeld.hpp (namespace openswmm::surface::ueb)
// ===========================================================================

/// A snow solve spanning many runoff/marcher steps; the end state is private until
/// consumed. Every initializer is zero, so Default is derived.
#[derive(Debug, Clone, Default)]
pub struct HeldState {
    pub accepted: LumpedState,
    pub begin: LumpedState,
    pub end: LumpedState,
    pub flux: LumpedFluxes,
    pub duration: f64,
    pub elapsed: f64,
}

impl HeldState {
    /// `elapsed < duration` exactly (false when either is NaN, as in C++).
    #[inline]
    pub fn pending(&self) -> bool {
        self.elapsed < self.duration
    }
}

// Declared in UebHeld.hpp, defined in UebHeld.cpp (translated later):
//   std::string prepareHeld(HeldState& state, const Parameters& parameters, const Site& site,
//                           const Forcing& forcing, const std::array<double, 10>& curve);
//   std::string consumeHeld(HeldState& state, double seconds,
//                           const std::array<double, 10>& curve, LumpedFluxes& output);
//   void closeHeld(HeldState& state) noexcept;

// ===========================================================================
// UebConfig.hpp (namespace openswmm::surface)
// ===========================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum SnowKind {
    DegreeDay = 0,
    Ueb = 1,
}

#[derive(Debug, Clone, Copy)]
pub struct SnowpackConfig {
    pub model: SnowKind,
    pub site: Site,
    pub initial_energy: f64,
}

impl Default for SnowpackConfig {
    fn default() -> Self {
        SnowpackConfig {
            model: SnowKind::Ueb,
            site: Site::default(),
            initial_energy: 0.0,
        }
    }
}

/// RadiationSource::Mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum RadiationSourceMode {
    Computed = 0,
    Timeseries = 1,
    File = 2,
}

#[derive(Debug, Clone)]
pub struct RadiationSource {
    pub mode: RadiationSourceMode,
    pub series: Bytes,
    pub authored: bool,
}

impl Default for RadiationSource {
    fn default() -> Self {
        RadiationSource {
            mode: RadiationSourceMode::Computed,
            series: Bytes::new(),
            authored: false,
        }
    }
}

/// Every initializer is zero, so Default is derived.
#[derive(Debug, Clone, Copy, Default)]
pub struct ClimateSite {
    pub latitude: f64,
    pub longitude: f64,
    pub timezone: f64,
    pub elevation: f64,
}

#[derive(Debug, Clone)]
pub struct SnowConfig {
    pub model: SnowKind,
    pub model_authored: bool,
    pub division_temperature: bool,
    pub step: f64,
    pub parameters: Parameters,
    pub bc_a: f64,
    pub bc_c: f64,
    // NaN/unspecified is represented explicitly: never invent a temperature range.
    pub monthly_range: [Option<f64>; 12],
    pub site: Option<ClimateSite>,
    pub shortwave: RadiationSource,
    pub longwave: RadiationSource,
    pub met_file: FilePathPair,
    pub packs: BTreeMap<i32, SnowpackConfig>,
}

impl Default for SnowConfig {
    fn default() -> Self {
        SnowConfig {
            model: SnowKind::DegreeDay,
            model_authored: false,
            division_temperature: false,
            step: 3600.0,
            parameters: Parameters::default(),
            bc_a: 0.8,
            bc_c: 2.4,
            monthly_range: [None; 12],
            site: None,
            shortwave: RadiationSource::default(),
            longwave: RadiationSource::default(),
            met_file: FilePathPair::default(),
            packs: BTreeMap::new(),
        }
    }
}

impl SnowConfig {
    /// Resolve a per-pack model override or the global snow-model fallback
    /// (upstream modelFor). No pack-existence/range validation.
    pub fn model_for(&self, pack: i32) -> SnowKind {
        match self.packs.get(&pack) {
            None => self.model,
            Some(it) => it.model,
        }
    }

    // Declared in UebConfig.hpp, defined in UebConfig.cpp (translated later):
    //   bool SnowConfig::hasTemperatureRows() const;
}

// Declared in UebConfig.hpp, defined in UebConfig.cpp (translated later):
//   std::string setUebOption(SnowConfig& config, const std::vector<std::string>& tokens,
//                            bool validate = true);   // note the default argument
//   std::vector<std::string> uebOptionRows(const SnowConfig& config);
//   std::optional<std::string> uebOptionValue(const SnowConfig& config, const std::string& key);
//   std::string setSnowTemperature(SnowConfig& config, const std::vector<std::string>& tokens);
//   std::string setSnowpackUeb(SnowConfig& config, int pack,
//                              const std::vector<std::string>& tokens);

// ===========================================================================
// MeshSnowConfig.hpp (namespace openswmm::surface)
// ===========================================================================

/// One authored pack selection and optional site/geometry override. Missing optionals
/// inherit lower-precedence scope values; explicit zero is a real override. pack is an
/// unresolved model name or normalized NONE token. Default is derived (empty pack, all
/// options absent), as upstream's initializer-less members give.
#[derive(Debug, Clone, Default)]
pub struct MeshSnowRow {
    pub pack: Bytes,
    pub forest: Option<f64>,
    pub drift: Option<f64>,
    pub ground_heat: Option<f64>,
    pub extinction_depth: Option<f64>,
    pub slope: Option<f64>,
    pub aspect: Option<f64>,
}

/// MeshSnowConfig::Mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum MeshSnowConfigMode {
    Auto = 0,
    Yes = 1,
    No = 2,
}

#[derive(Debug, Clone)]
pub struct MeshSnowConfig {
    pub mode: MeshSnowConfigMode,
    /// inherit the shared SNOW_STEP
    pub step: f64,
    /// * whole mesh, otherwise mesh tag
    pub defaults: BTreeMap<Bytes, MeshSnowRow>,
    /// zero-based cell index
    pub cells: BTreeMap<i32, MeshSnowRow>,
}

impl Default for MeshSnowConfig {
    fn default() -> Self {
        MeshSnowConfig {
            mode: MeshSnowConfigMode::Auto,
            step: 0.0,
            defaults: BTreeMap::new(),
            cells: BTreeMap::new(),
        }
    }
}

impl MeshSnowConfig {
    /// True for any stored tag/default/cell assignment row.
    pub fn assigned(&self) -> bool {
        !self.defaults.is_empty() || !self.cells.is_empty()
    }
}

// Declared in MeshSnowConfig.hpp, defined in MeshSnowConfig.cpp (translated later):
//   std::string setMeshSnowOption(MeshSnowConfig& config, const std::vector<std::string>& tokens);
//   std::string setMeshSnowRow(MeshSnowConfig& config, const std::vector<std::string>& tokens,
//                              bool cell);
//   std::vector<std::string> meshSnowOptionRows(const MeshSnowConfig& config);
//   std::vector<std::string> meshSnowRows(const MeshSnowConfig& config, bool cell);
//   MeshSnowRow resolveMeshSnowRow(const MeshSnowConfig& config, int cell, const std::string& tag);

// ===========================================================================
// UebSubcatchments.hpp (namespace openswmm::surface::ueb)
// ===========================================================================

/// Owned per-surface UEB predictions, accepted outputs and whole-model ledgers. Three
/// surface columns per subcatchment; rain/snow caches are per subcatchment. Every
/// initializer is zero (`forcing{}` value-initializes to zeros), so Default is derived.
#[derive(Debug, Clone, Default)]
pub struct NativeState {
    /// modes and values of runtime climate overrides
    pub forcing: [f64; 4],
    pub columns: Vec<HeldState>,
    pub last: Vec<LumpedFluxes>,
    pub rain: Vec<f64>,
    pub snow: Vec<f64>,
    /// per column: independent owner clocks
    pub sample_end: Vec<f64>,
    pub sample_seconds: Vec<f64>,
    pub active: bool,
    pub sublimation: f64,
    pub precipitation: f64,
    pub liquid: f64,
    pub throughfall: f64,
    pub net_energy: f64,
    pub geometry_energy: f64,
    pub plowed_energy: f64,
    pub initial_energy: f64,
    pub initial_water: f64,
    pub removed: f64,
    pub transferred_water: f64,
    pub last_seconds: f64,
    pub completed_seconds: f64,
}

// Declared in UebSubcatchments.hpp, defined in UebSubcatchments.cpp (translated later).
// SimulationContext and snow::SnowSoA are forward-declared there:
//   std::string initializeSubcatchments(SimulationContext& context, snow::SnowSoA& snow);
//   std::string advanceSubcatchments(SimulationContext& context, snow::SnowSoA& native,
//                                    double date, double elapsed, double seconds,
//                                    const double* rain, const double* snow);
//   double nextSnowBoundary(const SimulationContext& context, double elapsed,
//                           double proposed) noexcept;
