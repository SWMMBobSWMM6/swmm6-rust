// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/ForcingData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson (resize() body first pass by tools/soa2rs.py, the rest by
// hand). Not an OpenSWMM release.
//
// Per-element runtime forcing state, one Vec per attribute. Upstream's nested
// ForcingData::ElemClimateChannel is ForcingDataElemClimateChannel (Rust has no nested
// type declarations); its C array members `ElemClimateChannel x[2]` are `[_; 2]`.
// camelCase methods are snake_case (elemSlot is elem_slot, resetPerStep reset_per_step).
// The set_none / clear_resets lambdas are the free functions below. Sizes come from an
// unchecked `static_cast<std::size_t>(int)`, which `as usize` matches bit for bit.

use crate::data::heat_override::{HeatElemKind, HeatElement};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i8)]
pub enum ForcingMode {
    /// Use model-computed value (no forcing)
    None = 0,
    /// Replace computed value with user value
    Override = 1,
    /// Add user value to computed value
    Add = 2,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i8)]
pub enum ForcingPersist {
    /// Auto-clear after each timestep
    Reset = 0,
    /// Keep until explicitly cleared
    Persist = 1,
}

/// ForcingData::ElemClimateChannel — PE4 per-element climate, sized lazily on first use.
#[derive(Debug, Clone, Default)]
pub struct ForcingDataElemClimateChannel {
    pub mode: Vec<ForcingMode>,
    pub value: Vec<f64>,
    pub persist: Vec<ForcingPersist>,
}

impl ForcingDataElemClimateChannel {
    /// Allocate a channel shape, or keep one whose mode length already equals n.
    pub fn ensure(&mut self, n: usize) {
        if self.mode.len() == n {
            return;
        }
        self.mode.clear();
        self.mode.resize(n, ForcingMode::None);
        self.value.clear();
        self.value.resize(n, 0.0);
        self.persist.clear();
        self.persist.resize(n, ForcingPersist::Reset);
    }

    /// base for absent/NONE; value for OVERRIDE; base+value for any other mode.
    pub fn apply(&self, i: i32, base: f64) -> f64 {
        if i < 0 {
            return base;
        }
        let u = i as usize;
        if u >= self.mode.len() || self.mode[u] == ForcingMode::None {
            return base;
        }
        if self.mode[u] == ForcingMode::Override {
            self.value[u]
        } else {
            base + self.value[u]
        }
    }

    /// Disable only the entries marked RESET.
    pub fn reset_per_step(&mut self) {
        for i in 0..self.mode.len() {
            if self.persist[i] == ForcingPersist::Reset {
                self.mode[i] = ForcingMode::None;
            }
        }
    }

    pub fn clear(&mut self) {
        *self = ForcingDataElemClimateChannel::default();
    }
}

#[derive(Debug, Clone)]
pub struct ForcingData {
    // ------ Node forcing (sized to n_nodes)
    pub node_lat_inflow_mode: Vec<ForcingMode>,
    pub node_lat_inflow_value: Vec<f64>,
    pub node_lat_inflow_persist: Vec<ForcingPersist>,

    pub node_head_boundary_mode: Vec<ForcingMode>,
    pub node_head_boundary_value: Vec<f64>,
    pub node_head_boundary_persist: Vec<ForcingPersist>,

    /// flat 2D [node_idx * n_pollutants + pollutant_idx]
    pub node_quality_mode: Vec<ForcingMode>,
    /// mass rate (mass/sec)
    pub node_quality_value: Vec<f64>,
    pub node_quality_persist: Vec<ForcingPersist>,

    pub node_temperature_mode: Vec<ForcingMode>,
    pub node_temperature_value: Vec<f64>,
    pub node_temperature_persist: Vec<ForcingPersist>,

    pub node_age_mode: Vec<ForcingMode>,
    pub node_age_value: Vec<f64>,
    pub node_age_persist: Vec<ForcingPersist>,

    // ------ Link forcing (sized to n_links)
    pub link_flow_mode: Vec<ForcingMode>,
    pub link_flow_value: Vec<f64>,
    pub link_flow_persist: Vec<ForcingPersist>,

    pub link_setting_mode: Vec<ForcingMode>,
    /// 0.0–1.0 for pump/orifice/weir
    pub link_setting_value: Vec<f64>,
    pub link_setting_persist: Vec<ForcingPersist>,

    /// G-X4: conduit ⇄ groundwater exchange rate (`+` OUT of the conduit)
    pub link_seepage_mode: Vec<ForcingMode>,
    pub link_seepage_value: Vec<f64>,
    pub link_seepage_persist: Vec<ForcingPersist>,

    /// flattened link × pollutant
    pub link_quality_mode: Vec<ForcingMode>,
    /// OVERRIDE: concentration; ADD: mass rate (mass/sec)
    pub link_quality_value: Vec<f64>,
    pub link_quality_persist: Vec<ForcingPersist>,

    // ------ Subcatchment forcing (sized to n_subcatches)
    pub subcatch_rainfall_mode: Vec<ForcingMode>,
    /// user units (in/hr or mm/hr)
    pub subcatch_rainfall_value: Vec<f64>,
    pub subcatch_rainfall_persist: Vec<ForcingPersist>,

    pub subcatch_evap_mode: Vec<ForcingMode>,
    /// prescribed PET rate, ft/sec (internal units)
    pub subcatch_evap_value: Vec<f64>,
    pub subcatch_evap_persist: Vec<ForcingPersist>,

    pub subcatch_snowfall_mode: Vec<ForcingMode>,
    /// ft/sec (internal)
    pub subcatch_snowfall_value: Vec<f64>,
    pub subcatch_snowfall_persist: Vec<ForcingPersist>,

    // ------ Gage forcing (sized to n_gages)
    pub gage_rainfall_mode: Vec<ForcingMode>,
    /// user units (in/hr or mm/hr)
    pub gage_rainfall_value: Vec<f64>,
    pub gage_rainfall_persist: Vec<ForcingPersist>,

    // ------ Climate forcing (scalar — system-wide)
    pub climate_temperature_mode: ForcingMode,
    /// deg F (internal)
    pub climate_temperature_value: f64,
    pub climate_temperature_persist: ForcingPersist,

    pub climate_wind_mode: ForcingMode,
    /// mph (internal)
    pub climate_wind_value: f64,
    pub climate_wind_persist: ForcingPersist,

    pub climate_evap_mode: ForcingMode,
    /// ft/sec (internal)
    pub climate_evap_value: f64,
    pub climate_evap_persist: ForcingPersist,

    // ------ PE4: per-ELEMENT climate (API only), [0] = LINK, [1] = NODE
    /// deg F
    pub elem_air_temp: [ForcingDataElemClimateChannel; 2],
    /// %
    pub elem_humidity: [ForcingDataElemClimateChannel; 2],
    /// mph
    pub elem_wind: [ForcingDataElemClimateChannel; 2],
    /// W/m2
    pub elem_shortwave: [ForcingDataElemClimateChannel; 2],

    // ------ Counts (for iteration)
    pub n_nodes_: i32,
    pub n_links_: i32,
    pub n_subcatches_: i32,
    pub n_gages_: i32,
    pub n_pollutants_: i32,
}

impl Default for ForcingData {
    fn default() -> Self {
        ForcingData {
            node_lat_inflow_mode: Vec::new(),
            node_lat_inflow_value: Vec::new(),
            node_lat_inflow_persist: Vec::new(),
            node_head_boundary_mode: Vec::new(),
            node_head_boundary_value: Vec::new(),
            node_head_boundary_persist: Vec::new(),
            node_quality_mode: Vec::new(),
            node_quality_value: Vec::new(),
            node_quality_persist: Vec::new(),
            node_temperature_mode: Vec::new(),
            node_temperature_value: Vec::new(),
            node_temperature_persist: Vec::new(),
            node_age_mode: Vec::new(),
            node_age_value: Vec::new(),
            node_age_persist: Vec::new(),
            link_flow_mode: Vec::new(),
            link_flow_value: Vec::new(),
            link_flow_persist: Vec::new(),
            link_setting_mode: Vec::new(),
            link_setting_value: Vec::new(),
            link_setting_persist: Vec::new(),
            link_seepage_mode: Vec::new(),
            link_seepage_value: Vec::new(),
            link_seepage_persist: Vec::new(),
            link_quality_mode: Vec::new(),
            link_quality_value: Vec::new(),
            link_quality_persist: Vec::new(),
            subcatch_rainfall_mode: Vec::new(),
            subcatch_rainfall_value: Vec::new(),
            subcatch_rainfall_persist: Vec::new(),
            subcatch_evap_mode: Vec::new(),
            subcatch_evap_value: Vec::new(),
            subcatch_evap_persist: Vec::new(),
            subcatch_snowfall_mode: Vec::new(),
            subcatch_snowfall_value: Vec::new(),
            subcatch_snowfall_persist: Vec::new(),
            gage_rainfall_mode: Vec::new(),
            gage_rainfall_value: Vec::new(),
            gage_rainfall_persist: Vec::new(),
            climate_temperature_mode: ForcingMode::None,
            climate_temperature_value: 0.0,
            climate_temperature_persist: ForcingPersist::Reset,
            climate_wind_mode: ForcingMode::None,
            climate_wind_value: 0.0,
            climate_wind_persist: ForcingPersist::Reset,
            climate_evap_mode: ForcingMode::None,
            climate_evap_value: 0.0,
            climate_evap_persist: ForcingPersist::Reset,
            elem_air_temp: Default::default(),
            elem_humidity: Default::default(),
            elem_wind: Default::default(),
            elem_shortwave: Default::default(),
            n_nodes_: 0,
            n_links_: 0,
            n_subcatches_: 0,
            n_gages_: 0,
            n_pollutants_: 0,
        }
    }
}

/// The `set_none` lambda of clear_all().
fn set_none(mode_vec: &mut [ForcingMode]) {
    for m in mode_vec.iter_mut() {
        *m = ForcingMode::None;
    }
}

/// The `clear_resets` lambda of clear_reset_entries().
fn clear_resets(mode_vec: &mut [ForcingMode], persist_vec: &[ForcingPersist]) {
    for i in 0..mode_vec.len() {
        if persist_vec[i] == ForcingPersist::Reset {
            mode_vec[i] = ForcingMode::None;
        }
    }
}

impl ForcingData {
    /// 0 for LINK, 1 for NODE, otherwise -1.
    pub fn elem_slot(k: HeatElemKind) -> i32 {
        if k == HeatElemKind::Link {
            return 0;
        }
        if k == HeatElemKind::Node {
            return 1;
        }
        -1
    }

    /// Effective air temperature for one heat element, deg F.
    pub fn element_air_temp_f(&self, e: &HeatElement, base: f64) -> f64 {
        let s = Self::elem_slot(e.kind);
        if s < 0 {
            base
        } else {
            self.elem_air_temp[s as usize].apply(e.index, base)
        }
    }

    /// Effective relative humidity for one heat element, %.
    pub fn element_humidity(&self, e: &HeatElement, base: f64) -> f64 {
        let s = Self::elem_slot(e.kind);
        if s < 0 {
            base
        } else {
            self.elem_humidity[s as usize].apply(e.index, base)
        }
    }

    /// Effective wind speed for one heat element, mph.
    pub fn element_wind_mph(&self, e: &HeatElement, base: f64) -> f64 {
        let s = Self::elem_slot(e.kind);
        if s < 0 {
            base
        } else {
            self.elem_wind[s as usize].apply(e.index, base)
        }
    }

    /// Effective incoming shortwave for one heat element, W/m2.
    pub fn element_shortwave(&self, e: &HeatElement, base: f64) -> f64 {
        let s = Self::elem_slot(e.kind);
        if s < 0 {
            base
        } else {
            self.elem_shortwave[s as usize].apply(e.index, base)
        }
    }

    /// Replace object-channel arrays with NONE/zero/RESET defaults. Scalar climate and the
    /// lazily allocated element-climate channels are unchanged.
    pub fn resize(
        &mut self,
        n_nodes: i32,
        n_links: i32,
        n_subcatches: i32,
        n_gages: i32,
        n_pollutants: i32,
    ) {
        self.n_nodes_ = n_nodes;
        self.n_links_ = n_links;
        self.n_subcatches_ = n_subcatches;
        self.n_gages_ = n_gages;
        self.n_pollutants_ = n_pollutants;

        let un = n_nodes as usize;
        let ul = n_links as usize;
        let us = n_subcatches as usize;
        let ug = n_gages as usize;
        let unp = (n_nodes as usize) * (n_pollutants as usize);

        self.node_lat_inflow_mode.clear();
        self.node_lat_inflow_mode.resize(un, ForcingMode::None);
        self.node_lat_inflow_value.clear();
        self.node_lat_inflow_value.resize(un, 0.0);
        self.node_lat_inflow_persist.clear();
        self.node_lat_inflow_persist.resize(un, ForcingPersist::Reset);
        self.node_head_boundary_mode.clear();
        self.node_head_boundary_mode.resize(un, ForcingMode::None);
        self.node_head_boundary_value.clear();
        self.node_head_boundary_value.resize(un, 0.0);
        self.node_head_boundary_persist.clear();
        self.node_head_boundary_persist.resize(un, ForcingPersist::Reset);
        self.node_quality_mode.clear();
        self.node_quality_mode.resize(unp, ForcingMode::None);
        self.node_quality_value.clear();
        self.node_quality_value.resize(unp, 0.0);
        self.node_quality_persist.clear();
        self.node_quality_persist.resize(unp, ForcingPersist::Reset);
        self.node_temperature_mode.clear();
        self.node_temperature_mode.resize(un, ForcingMode::None);
        self.node_temperature_value.clear();
        self.node_temperature_value.resize(un, 0.0);
        self.node_temperature_persist.clear();
        self.node_temperature_persist.resize(un, ForcingPersist::Reset);
        self.node_age_mode.clear();
        self.node_age_mode.resize(un, ForcingMode::None);
        self.node_age_value.clear();
        self.node_age_value.resize(un, 0.0);
        self.node_age_persist.clear();
        self.node_age_persist.resize(un, ForcingPersist::Reset);
        self.link_flow_mode.clear();
        self.link_flow_mode.resize(ul, ForcingMode::None);
        self.link_flow_value.clear();
        self.link_flow_value.resize(ul, 0.0);
        self.link_flow_persist.clear();
        self.link_flow_persist.resize(ul, ForcingPersist::Reset);
        self.link_setting_mode.clear();
        self.link_setting_mode.resize(ul, ForcingMode::None);
        self.link_setting_value.clear();
        self.link_setting_value.resize(ul, 0.0);
        self.link_setting_persist.clear();
        self.link_setting_persist.resize(ul, ForcingPersist::Reset);
        self.link_seepage_mode.clear();
        self.link_seepage_mode.resize(ul, ForcingMode::None);
        self.link_seepage_value.clear();
        self.link_seepage_value.resize(ul, 0.0);
        self.link_seepage_persist.clear();
        self.link_seepage_persist.resize(ul, ForcingPersist::Reset);
        let ulp = (n_links as usize) * (n_pollutants as usize);
        self.link_quality_mode.clear();
        self.link_quality_mode.resize(ulp, ForcingMode::None);
        self.link_quality_value.clear();
        self.link_quality_value.resize(ulp, 0.0);
        self.link_quality_persist.clear();
        self.link_quality_persist.resize(ulp, ForcingPersist::Reset);
        self.subcatch_rainfall_mode.clear();
        self.subcatch_rainfall_mode.resize(us, ForcingMode::None);
        self.subcatch_rainfall_value.clear();
        self.subcatch_rainfall_value.resize(us, 0.0);
        self.subcatch_rainfall_persist.clear();
        self.subcatch_rainfall_persist.resize(us, ForcingPersist::Reset);
        self.subcatch_evap_mode.clear();
        self.subcatch_evap_mode.resize(us, ForcingMode::None);
        self.subcatch_evap_value.clear();
        self.subcatch_evap_value.resize(us, 0.0);
        self.subcatch_evap_persist.clear();
        self.subcatch_evap_persist.resize(us, ForcingPersist::Reset);
        self.subcatch_snowfall_mode.clear();
        self.subcatch_snowfall_mode.resize(us, ForcingMode::None);
        self.subcatch_snowfall_value.clear();
        self.subcatch_snowfall_value.resize(us, 0.0);
        self.subcatch_snowfall_persist.clear();
        self.subcatch_snowfall_persist.resize(us, ForcingPersist::Reset);
        self.gage_rainfall_mode.clear();
        self.gage_rainfall_mode.resize(ug, ForcingMode::None);
        self.gage_rainfall_value.clear();
        self.gage_rainfall_value.resize(ug, 0.0);
        self.gage_rainfall_persist.clear();
        self.gage_rainfall_persist.resize(ug, ForcingPersist::Reset);
    }

    /// Disable object and global-climate forcing modes (values, persistence, dimensions
    /// and the per-element channels are untouched).
    pub fn clear_all(&mut self) {
        set_none(&mut self.node_lat_inflow_mode);
        set_none(&mut self.node_head_boundary_mode);
        set_none(&mut self.node_quality_mode);
        set_none(&mut self.node_temperature_mode);
        set_none(&mut self.node_age_mode);
        set_none(&mut self.link_flow_mode);
        set_none(&mut self.link_setting_mode);
        set_none(&mut self.link_seepage_mode); // G-X4
        set_none(&mut self.link_quality_mode);
        set_none(&mut self.subcatch_rainfall_mode);
        set_none(&mut self.subcatch_evap_mode);
        set_none(&mut self.subcatch_snowfall_mode);
        set_none(&mut self.gage_rainfall_mode);
        self.climate_temperature_mode = ForcingMode::None;
        self.climate_wind_mode = ForcingMode::None;
        self.climate_evap_mode = ForcingMode::None;
    }

    /// Clear only RESET-persistence entries (called at end of each step), including the
    /// per-element climate channels.
    pub fn clear_reset_entries(&mut self) {
        clear_resets(&mut self.node_lat_inflow_mode, &self.node_lat_inflow_persist);
        clear_resets(&mut self.node_head_boundary_mode, &self.node_head_boundary_persist);
        clear_resets(&mut self.node_quality_mode, &self.node_quality_persist);
        clear_resets(&mut self.node_temperature_mode, &self.node_temperature_persist);
        clear_resets(&mut self.node_age_mode, &self.node_age_persist);
        clear_resets(&mut self.link_flow_mode, &self.link_flow_persist);
        clear_resets(&mut self.link_setting_mode, &self.link_setting_persist);
        clear_resets(&mut self.link_seepage_mode, &self.link_seepage_persist); // G-X4
        clear_resets(&mut self.link_quality_mode, &self.link_quality_persist);
        clear_resets(&mut self.subcatch_rainfall_mode, &self.subcatch_rainfall_persist);
        clear_resets(&mut self.subcatch_evap_mode, &self.subcatch_evap_persist);
        clear_resets(&mut self.subcatch_snowfall_mode, &self.subcatch_snowfall_persist);
        clear_resets(&mut self.gage_rainfall_mode, &self.gage_rainfall_persist);
        if self.climate_temperature_persist == ForcingPersist::Reset {
            self.climate_temperature_mode = ForcingMode::None;
        }
        if self.climate_wind_persist == ForcingPersist::Reset {
            self.climate_wind_mode = ForcingMode::None;
        }
        if self.climate_evap_persist == ForcingPersist::Reset {
            self.climate_evap_mode = ForcingMode::None;
        }
        // PE4: the per-element channels join the SAME sweep.
        for k in 0..2 {
            self.elem_air_temp[k].reset_per_step();
            self.elem_humidity[k].reset_per_step();
            self.elem_wind[k].reset_per_step();
            self.elem_shortwave[k].reset_per_step();
        }
    }

    /// Effective air temperature (deg F internal).
    pub fn effective_temperature(&self, broadcast: f64) -> f64 {
        match self.climate_temperature_mode {
            ForcingMode::Override => self.climate_temperature_value,
            ForcingMode::Add => broadcast + self.climate_temperature_value,
            _ => broadcast,
        }
    }

    /// Effective wind speed (mph internal).
    pub fn effective_wind(&self, broadcast: f64) -> f64 {
        match self.climate_wind_mode {
            ForcingMode::Override => self.climate_wind_value,
            ForcingMode::Add => broadcast + self.climate_wind_value,
            _ => broadcast,
        }
    }

    /// Effective rainfall for a subcatchment, in the caller's user units.
    pub fn effective_rainfall(&self, ui: usize, gage_rainfall: f64) -> f64 {
        if ui >= self.subcatch_rainfall_mode.len() {
            return gage_rainfall;
        }
        match self.subcatch_rainfall_mode[ui] {
            ForcingMode::Override => self.subcatch_rainfall_value[ui],
            ForcingMode::Add => gage_rainfall + self.subcatch_rainfall_value[ui],
            _ => gage_rainfall,
        }
    }

    /// Effective system-wide evaporation rate (ft/sec).
    pub fn effective_climate_evap(&self, broadcast: f64) -> f64 {
        match self.climate_evap_mode {
            ForcingMode::Override => self.climate_evap_value,
            ForcingMode::Add => broadcast + self.climate_evap_value,
            _ => broadcast,
        }
    }

    /// Effective snowfall for a subcatchment (ft/sec).
    pub fn effective_snowfall(&self, ui: usize, gage_snowfall: f64) -> f64 {
        if ui >= self.subcatch_snowfall_mode.len() {
            return gage_snowfall;
        }
        match self.subcatch_snowfall_mode[ui] {
            ForcingMode::Override => self.subcatch_snowfall_value[ui],
            ForcingMode::Add => gage_snowfall + self.subcatch_snowfall_value[ui],
            _ => gage_snowfall,
        }
    }

    /// Effective potential evaporation rate for a subcatchment (ft/sec).
    pub fn effective_evap_rate(&self, ui: usize, broadcast_rate: f64) -> f64 {
        // Total over unallocated forcing (hand-built / partially-initialized contexts).
        if ui >= self.subcatch_evap_mode.len() {
            return broadcast_rate;
        }
        match self.subcatch_evap_mode[ui] {
            ForcingMode::Override => self.subcatch_evap_value[ui],
            ForcingMode::Add => broadcast_rate + self.subcatch_evap_value[ui],
            _ => broadcast_rate,
        }
    }
}
