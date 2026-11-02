// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/HeatData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Heat-transport configuration and runtime state (heat plan H1-H6, PE). Defaults are
// upstream's member initializers value for value; RadiativeConfig::landcover_temp is
// quiet_NaN, whose bits (0x7ff8000000000000) f64::NAN shares. Upstream's static
// constexpr members are associated consts (kDefaultTemp is K_DEFAULT_TEMP, kNSubArea
// K_N_SUB_AREA); camelCase methods are snake_case; HeatOverrideData::empty() is
// is_empty(). Sizes come from `static_cast<std::size_t>(int)` unchecked, which `as usize`
// matches bit for bit (a negative count becomes huge, and resize then fails, as in C++).

use crate::data::bed_zone::SedimentConfig;
use crate::data::heat_override::HeatOverrideRow;

/// Source pathways with a configurable inlet temperature ([HEAT_SOURCES]). Order is the
/// storage index of HeatConfigData::global_temp, matching WaterAgeSource one-for-one.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum HeatSource {
    /// washoff runoff
    Rainfall = 0,
    Dwf = 1,
    Gw = 2,
    Rdii = 3,
    ExternalInflow = 4,
    Iface = 5,
    /// water in the network at t = 0
    InitialState = 6,
    Count = 7,
}

/// H5a: which subarea of a subcatchment a ponded temperature belongs to. Deliberately a
/// separate enum from water_age::SubArea.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum HeatSubArea {
    Imperv0 = 0,
    Imperv1 = 1,
    Perv = 2,
    Count = 3,
}

/// What a dry or absent element reports (plan D-H5c).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum DryTempPolicy {
    /// Freeze the last wet value; rewetting mixes against it.
    Hold = 0,
    /// Track air temperature.
    Air = 1,
    /// Fall to `HeatConfigData::K_DEFAULT_TEMP`.
    Default = 2,
}

/// Vertical conduction between LID layers (plan §6.1 D-H5b, H5b). GWComponent's
/// constructor values, not its dead in-class initializers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConductionConfig {
    /// W/m/K
    pub water_conductivity: f64,
    /// W/m/K
    pub sed_conductivity: f64,
    /// kg/m³
    pub sed_density: f64,
    /// J/kg/K
    pub sed_specific_heat: f64,
}

impl Default for ConductionConfig {
    fn default() -> Self {
        ConductionConfig {
            water_conductivity: 0.606,
            sed_conductivity: 2.6,
            sed_density: 1970.0,
            sed_specific_heat: 2758.0,
        }
    }
}

/// Where incoming shortwave `Jin` comes from (plan §2.5, phase H6a).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum ShortwaveMode {
    /// `SHORTWAVE GLOBAL <W/m²>`
    Constant = 0,
    /// `SHORTWAVE GLOBAL TIMESERIES <name>`
    Timeseries = 1,
    /// `SHORTWAVE GLOBAL COMPUTED`
    Computed = 2,
}

/// `[SOLAR_RADIATION]` — site geometry and Bird atmosphere (H6a).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolarConfig {
    /// +N. REQUIRED under COMPUTED.
    pub latitude_deg: f64,
    /// +E. REQUIRED under COMPUTED.
    pub longitude_deg: f64,
    /// Offset from UTC, +E (e.g. MST = -7).
    pub timezone_hours: f64,
    pub has_latitude: bool,
    pub has_longitude: bool,
    pub has_timezone: bool,
    /// Site elevation, metres, for the Bird pressure term.
    pub elevation_m: f64,
    pub has_elevation: bool,
    /// Aerosol optical depth at 380 nm
    pub aod380: f64,
    /// Aerosol optical depth at 500 nm
    pub aod500: f64,
    /// Precipitable water vapour, cm
    pub precip_water_cm: f64,
    /// Ozone column, cm (NTP)
    pub ozone_cm: f64,
    /// Surface albedo for the sky-ground multiple-reflection term (not the water's).
    pub ground_albedo: f64,
}

impl Default for SolarConfig {
    fn default() -> Self {
        SolarConfig {
            latitude_deg: 0.0,
            longitude_deg: 0.0,
            timezone_hours: 0.0,
            has_latitude: false,
            has_longitude: false,
            has_timezone: false,
            elevation_m: 0.0,
            has_elevation: false,
            aod380: 0.30,
            aod500: 0.20,
            precip_water_cm: 1.42,
            ozone_cm: 0.34,
            ground_albedo: 0.20,
        }
    }
}

/// `[CLOUD_COVER]` — one fraction driving both shortwave and longwave (H6a, D-H6a-2).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CloudConfig {
    /// No `[CLOUD_COVER]` section → clear sky.
    pub configured: bool,
    pub use_timeseries: bool,
    /// Index into `ctx.tables`, -1 = none.
    pub ts_index: i32,
    /// C ∈ [0,1], constant spelling.
    pub fraction: f64,
    /// Kasten–Czeplak k
    pub sw_atten_k: f64,
    /// Kasten–Czeplak n
    pub sw_atten_n: f64,
    /// Bolz k_lw
    pub lw_cloud_k: f64,
}

impl Default for CloudConfig {
    fn default() -> Self {
        CloudConfig {
            configured: false,
            use_timeseries: false,
            ts_index: -1,
            fraction: 0.0,
            sw_atten_k: 0.75,
            sw_atten_n: 3.4,
            lw_cloud_k: 0.17,
        }
    }
}

/// `[RADIATIVE_FLUXES]` parameters (heat plan §2.2, phase H3).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RadiativeConfig {
    /// Incoming solar Jin, W/m² (0 = night)
    pub shortwave_wm2: f64,
    /// Rs, reference default
    pub albedo: f64,
    /// fs, 0 = unshaded
    pub shade_factor: f64,
    /// fsky, 1 = open sky
    pub sky_view: f64,
    /// εw
    pub emiss_water: f64,
    /// εlc
    pub emiss_landcover: f64,
    /// Brunt Aa
    pub atm_emiss_coeff: f64,
    /// RL
    pub lw_reflection: f64,
    /// Land-cover radiating temperature, °C. NaN ⇒ use air temperature (PE2).
    pub landcover_temp: f64,
    pub sw_mode: ShortwaveMode,
    /// Index into `ctx.tables` under TIMESERIES; -1 otherwise.
    pub sw_ts_index: i32,
}

impl Default for RadiativeConfig {
    fn default() -> Self {
        RadiativeConfig {
            shortwave_wm2: 0.0,
            albedo: 0.0,
            shade_factor: 0.0,
            sky_view: 1.0,
            emiss_water: 0.97,
            emiss_landcover: 0.97,
            atm_emiss_coeff: 0.5,
            lw_reflection: 0.03,
            landcover_temp: f64::NAN, // std::numeric_limits<double>::quiet_NaN()
            sw_mode: ShortwaveMode::Constant,
            sw_ts_index: -1,
        }
    }
}

/// Dense per-element attribute storage (PE2, D-PE2). Sized ONLY when at least one
/// override row targets the family; empty means "use the global".
#[derive(Debug, Clone, Default)]
pub struct HeatOverrideData {
    /// as parsed; the serializer's input
    pub rows: Vec<HeatOverrideRow>,
    /// [link], empty ⇒ use global
    pub rad_link: Vec<RadiativeConfig>,
    /// [node], empty ⇒ use global
    pub rad_node: Vec<RadiativeConfig>,
    /// [link], empty ⇒ use global
    pub sed_link: Vec<SedimentConfig>,
    /// [surface mesh cell]
    pub rad_cell: Vec<RadiativeConfig>,
    /// [surface mesh cell], not ARD cells
    pub sed_cell: Vec<SedimentConfig>,
    pub resolved: bool,
}

// resolved = false is Rust's default, so Default is derived.

impl HeatOverrideData {
    /// rows.empty(); does not inspect the resolved dense vectors or the resolved flag.
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }
    pub fn clear(&mut self) {
        *self = HeatOverrideData::default();
    }
}

/// Parsed `model.heat` state (heat component, phase H1).
#[derive(Debug, Clone)]
pub struct HeatConfigData {
    pub configured: bool,
    /// `[HEAT_FLUXES] RADIATIVE_EXCHANGE ON` (H3).
    pub radiative_exchange: bool,
    pub radiative: RadiativeConfig,
    /// `[HEAT_FLUXES] SURFACE_EXCHANGE ON` (H2).
    pub surface_exchange: bool,
    /// `[HEAT_FLUXES] DRY_ELEMENT_TEMPERATURE HOLD|AIR|DEFAULT` (D-H5c).
    pub dry_temp_policy: DryTempPolicy,
    /// `[HEAT_FLUXES] LAYER_CONDUCTION ON` (H5b).
    pub layer_conduction: bool,
    pub conduction: ConductionConfig,
    /// `[HEAT_FLUXES] SEDIMENT_EXCHANGE ON` (H6b).
    pub sediment_exchange: bool,
    /// `[SEDIMENT_EXCHANGE]`
    pub sediment: SedimentConfig,
    /// PE2 — per-element overrides for `radiative` and `sediment`.
    pub overrides: HeatOverrideData,
    /// `[SOLAR_RADIATION]`
    pub solar: SolarConfig,
    /// `[CLOUD_COVER]`
    pub cloud: CloudConfig,
    /// GLOBAL inlet temperature per source, °C.
    pub global_temp: [f64; HeatSource::Count as usize],
    /// Whether the user set this source explicitly (vs. taking the default).
    pub configured_source: [bool; HeatSource::Count as usize],
    /// HeatSource as int (NODE-scope overrides, parallel arrays)
    pub node_over_source: Vec<i32>,
    /// node index
    pub node_over_node: Vec<i32>,
    /// °C
    pub node_over_temp: Vec<f64>,
}

impl HeatConfigData {
    /// Default inlet temperature when a source has no row (°C).
    pub const K_DEFAULT_TEMP: f64 = 20.0;

    /// First matching node/source override, otherwise global_temp[s].
    pub fn source_temp(&self, s: HeatSource, node: i32) -> f64 {
        for i in 0..self.node_over_source.len() {
            if self.node_over_source[i] == s as i32 && self.node_over_node[i] == node {
                return self.node_over_temp[i];
            }
        }
        self.global_temp[s as i32 as usize]
    }
}

impl Default for HeatConfigData {
    fn default() -> Self {
        HeatConfigData {
            configured: false,
            radiative_exchange: false,
            radiative: RadiativeConfig::default(),
            surface_exchange: false,
            dry_temp_policy: DryTempPolicy::Hold,
            layer_conduction: false,
            conduction: ConductionConfig::default(),
            sediment_exchange: false,
            sediment: SedimentConfig::default(),
            overrides: HeatOverrideData::default(),
            solar: SolarConfig::default(),
            cloud: CloudConfig::default(),
            global_temp: [Self::K_DEFAULT_TEMP; HeatSource::Count as usize],
            configured_source: [false; HeatSource::Count as usize],
            node_over_source: Vec::new(),
            node_over_node: Vec::new(),
            node_over_temp: Vec::new(),
        }
    }
}

/// Runtime heat state shared by the engines (phase H1).
#[derive(Debug, Clone)]
pub struct HeatState {
    /// Empty in prescribed compatibility mode.
    pub mesh_contact_link: Vec<u8>,
    /// [node], °C·ft³/s
    pub node_temp_vol_in: Vec<f64>,
    /// [node], °C·ft³/s arriving through LID underdrains
    pub node_lid_drain_temp_vol_in: Vec<f64>,
    /// [node], °C
    pub node_temp: Vec<f64>,
    /// [link], °C
    pub link_temp: Vec<f64>,
    /// The LEGACY mirror seeds INITIAL_STATE on its first step.
    pub legacy_seeded: bool,
    /// Incoming shortwave at the current step, W/m². NEGATIVE means "not yet resolved".
    pub shortwave_now: f64,
    /// Cloud fraction at the current step, C ∈ [0,1]. 0 = clear.
    pub cloud_now: f64,
    /// [subcatch * K_N_SUB_AREA + subarea], °C
    pub subarea_temp: Vec<f64>,
    /// [subcatch * K_N_SUB_AREA + subarea], ft³ at the END of the previous step
    pub subarea_vol_prev: Vec<f64>,
    /// [subcatch], °C
    pub subcatch_runoff_temp: Vec<f64>,
    /// [subcatch], °C·ft³/s
    pub subcatch_runon_temp_vol_in: Vec<f64>,
    /// [subcatch], ft³/s — the run-on rate whose temperature is KNOWN
    pub subcatch_runon_temp_rate: Vec<f64>,
    /// [subcatch], °C·ft³
    pub subcatch_outfall_temp_vol: Vec<f64>,
    /// LID underdrain water sent to a subcatchment: cfs x temperature
    pub subcatch_lid_drain_temp_cfs: Vec<f64>,
}

impl Default for HeatState {
    fn default() -> Self {
        HeatState {
            mesh_contact_link: Vec::new(),
            node_temp_vol_in: Vec::new(),
            node_lid_drain_temp_vol_in: Vec::new(),
            node_temp: Vec::new(),
            link_temp: Vec::new(),
            legacy_seeded: false,
            shortwave_now: -1.0,
            cloud_now: 0.0,
            subarea_temp: Vec::new(),
            subarea_vol_prev: Vec::new(),
            subcatch_runoff_temp: Vec::new(),
            subcatch_runon_temp_vol_in: Vec::new(),
            subcatch_runon_temp_rate: Vec::new(),
            subcatch_outfall_temp_vol: Vec::new(),
            subcatch_lid_drain_temp_cfs: Vec::new(),
        }
    }
}

impl HeatState {
    pub const K_N_SUB_AREA: i32 = HeatSubArea::Count as i32;

    /// Replace network heat arrays and reset current forcing caches. Preserves
    /// mesh_contact_link and all watershed arrays.
    pub fn resize(&mut self, n_nodes: i32, n_links: i32, initial_temp: f64) {
        self.node_temp_vol_in.clear();
        self.node_temp_vol_in.resize(n_nodes as usize, 0.0);
        self.node_lid_drain_temp_vol_in.clear();
        self.node_lid_drain_temp_vol_in.resize(n_nodes as usize, 0.0);
        self.node_temp.clear();
        self.node_temp.resize(n_nodes as usize, initial_temp);
        self.link_temp.clear();
        self.link_temp.resize(n_links as usize, initial_temp);
        self.legacy_seeded = false;
        self.shortwave_now = -1.0;
        self.cloud_now = 0.0;
    }

    /// Replace subcatchment heat arrays and their three-subarea histories.
    pub fn resize_watershed(&mut self, n_subcatch: i32, initial_temp: f64) {
        let n = n_subcatch as usize;
        self.subarea_temp.clear();
        self.subarea_temp.resize(n * Self::K_N_SUB_AREA as usize, initial_temp);
        self.subarea_vol_prev.clear();
        self.subarea_vol_prev.resize(n * Self::K_N_SUB_AREA as usize, 0.0);
        self.subcatch_runoff_temp.clear();
        self.subcatch_runoff_temp.resize(n, initial_temp);
        self.subcatch_runon_temp_vol_in.clear();
        self.subcatch_runon_temp_vol_in.resize(n, 0.0);
        self.subcatch_runon_temp_rate.clear();
        self.subcatch_runon_temp_rate.resize(n, 0.0);
        self.subcatch_outfall_temp_vol.clear();
        self.subcatch_outfall_temp_vol.resize(n, 0.0);
        self.subcatch_lid_drain_temp_cfs.clear();
        self.subcatch_lid_drain_temp_cfs.resize(n, 0.0);
    }

    /// True if subarea_temp has n_subcatch*K_N_SUB_AREA entries; other arrays unchecked.
    pub fn watershed_sized(&self, n_subcatch: i32) -> bool {
        self.subarea_temp.len() == (n_subcatch as usize) * (Self::K_N_SUB_AREA as usize)
    }

    pub fn clear(&mut self) {
        *self = HeatState::default();
    }
}
