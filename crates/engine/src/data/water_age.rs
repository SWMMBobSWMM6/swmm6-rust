// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/WaterAgeData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Water-age configuration and runtime state (water age plan A1-A4, Z1). Every member
// initializer is Rust's default, so Default is derived. node_ext_inflow_age is filled
// with quiet_NaN, whose bits (0x7ff8000000000000) f64::NAN shares. Sizes come from an
// unchecked `static_cast<std::size_t>(int)`, which `as usize` matches bit for bit.

/// Source pathways with configurable initial age ([WATER_AGE_SOURCES]). Order is the
/// storage index of WaterAgeConfigData::global_age.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum WaterAgeSource {
    /// washoff runoff (and LID drains until A4)
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

/// Parsed `model.age` state (waterage component, phase A1a).
#[derive(Debug, Clone, Default)]
pub struct WaterAgeConfigData {
    pub configured: bool,
    /// GLOBAL initial age per source, SECONDS (default 0 = fresh water).
    pub global_age: [f64; WaterAgeSource::Count as usize],
    /// WaterAgeSource as int (NODE-scope overrides, parallel arrays)
    pub node_over_source: Vec<i32>,
    /// node index
    pub node_over_node: Vec<i32>,
    /// seconds
    pub node_over_age: Vec<f64>,
}

impl WaterAgeConfigData {
    /// First matching node/source override, otherwise global_age[s]. Seconds.
    pub fn source_age(&self, s: WaterAgeSource, node: i32) -> f64 {
        for i in 0..self.node_over_source.len() {
            if self.node_over_source[i] == s as i32 && self.node_over_node[i] == node {
                return self.node_over_age[i];
            }
        }
        self.global_age[s as i32 as usize]
    }
}

/// A3: subarea index within a subcatchment, matching the RunoffSolver's three depths.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum SubArea {
    /// impervious, no depression storage
    Imperv0 = 0,
    /// impervious with depression storage
    Imperv1 = 1,
    Perv = 2,
    Count = 3,
}

/// Runtime age state shared by the engines (phase A1a).
#[derive(Debug, Clone, Default)]
pub struct WaterAgeState {
    /// [node], age·ft³/s
    pub node_age_vol_in: Vec<f64>,
    /// [node], seconds
    pub node_age: Vec<f64>,
    /// [link], seconds
    pub link_age: Vec<f64>,
    /// A1b: the LEGACY mirror seeds INITIAL_STATE on its first step.
    pub legacy_seeded: bool,
    /// [subcatch*3 + subarea] mean age of the water ponded on each subarea, SECONDS.
    pub subarea_age: Vec<f64>,
    /// [subcatch*3 + subarea] the previous step's stored volume (ft³).
    pub subarea_vol_prev: Vec<f64>,
    /// [subcatch] age of the water LEAVING as runoff, SECONDS.
    pub subcatch_runoff_age: Vec<f64>,
    /// [subcatch] age·ft³/s arriving as RUN-ON from upstream subcatchments.
    pub subcatch_runon_age_vol_in: Vec<f64>,
    /// A4: [node] age·ft³/s arriving through LID underdrains.
    pub node_lid_drain_age_vol_in: Vec<f64>,
    /// A4: [subcatch] age·ft³ (LID drain) waiting for assembleRunon.
    pub subcatch_lid_drain_age_cfs: Vec<f64>,
    /// A4: [subcatch] age·ft³ (outfall) waiting for assembleRunon.
    pub subcatch_outfall_age_vol: Vec<f64>,
    /// Z1: [node] age (SECONDS) prescribed by an [INFLOWS] row; quiet_NaN where none.
    pub node_ext_inflow_age: Vec<f64>,
    /// A2a: set by HotStartManager::apply when a V3 file restored ages.
    pub hotstart_loaded: bool,
}

impl WaterAgeState {
    /// Replace network and watershed age state, preserving node_ext_inflow_age when its
    /// length already matches n_nodes.
    pub fn resize(&mut self, n_nodes: i32, n_links: i32, n_subcatch: i32) {
        self.node_age_vol_in.clear();
        self.node_age_vol_in.resize(n_nodes as usize, 0.0);
        self.node_age.clear();
        self.node_age.resize(n_nodes as usize, 0.0);
        self.link_age.clear();
        self.link_age.resize(n_links as usize, 0.0);
        self.node_lid_drain_age_vol_in.clear();
        self.node_lid_drain_age_vol_in.resize(n_nodes as usize, 0.0);
        // Z1: PRESERVE when already sized.
        if self.node_ext_inflow_age.len() != n_nodes as usize {
            self.node_ext_inflow_age.clear();
            self.node_ext_inflow_age.resize(n_nodes as usize, f64::NAN); // quiet_NaN
        }
        let ns3 = (n_subcatch as usize) * (SubArea::Count as i32 as usize);
        self.subarea_age.clear();
        self.subarea_age.resize(ns3, 0.0);
        self.subarea_vol_prev.clear();
        self.subarea_vol_prev.resize(ns3, 0.0);
        self.subcatch_runoff_age.clear();
        self.subcatch_runoff_age.resize(n_subcatch as usize, 0.0);
        self.subcatch_runon_age_vol_in.clear();
        self.subcatch_runon_age_vol_in.resize(n_subcatch as usize, 0.0);
        self.subcatch_lid_drain_age_cfs.clear();
        self.subcatch_lid_drain_age_cfs.resize(n_subcatch as usize, 0.0);
        self.subcatch_outfall_age_vol.clear();
        self.subcatch_outfall_age_vol.resize(n_subcatch as usize, 0.0);
        self.legacy_seeded = false;
        self.hotstart_loaded = false;
    }

    pub fn clear(&mut self) {
        *self = WaterAgeState::default();
    }
}
