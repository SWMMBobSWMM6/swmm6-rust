// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/HydrologyData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The structure-of-arrays stores for snowpacks, aquifers, LID controls and LID usage,
// one Vec per attribute in upstream's field order (including SnowpackStore's, whose
// runtime/config members precede its names column). The stores have no member
// initializers of their own; each member's Default carries its type's upstream values,
// so Default is derived. std::pair<A, B> is (A, B); std::string is Bytes. count() is
// `static_cast<int>(size())`, which `len() as i32` matches bit for bit.
// resize_wb takes `static_cast<std::size_t>(n)` unchecked, as `n as usize` does: a
// negative n becomes huge, and assign/resize then fails (std::length_error or
// bad_alloc in C++, a capacity-overflow panic or allocation abort here).

#![allow(clippy::all)]
use crate::core::options::{Bytes, FilePathPair};
use crate::data::lid_node::LidNodeLayer;
use crate::hydrology::hourly_meteorology::{HourlyMeteorology, SnowDailyWeather};
use crate::hydrology::ueb_types::{MeshSnowConfig, NativeState, SnowConfig};

// ============================================================================
// Snowpack definitions (from [SNOWPACKS] section)
// ============================================================================

/// SoA storage for snowpack parameter sets. Each snowpack has PLOWABLE, IMPERVIOUS,
/// PERVIOUS surface params (7 values each) and optional REMOVAL params.
#[derive(Debug, Clone, Default)]
pub struct SnowpackStore {
    // Derived at report time from accepted native and completed source owners.
    // Not restart state; it is rebuilt after every run/resume.
    pub report_values: Vec<(Bytes, f64)>,
    pub runtime: NativeState,
    pub config: SnowConfig,
    pub mesh_config: MeshSnowConfig,
    pub meteorology: HourlyMeteorology,
    pub daily_weather: Vec<SnowDailyWeather>,

    pub names: Vec<Bytes>,

    /// PLOWABLE surface: cmin, cmax, tbase, fwfrac, sd0, fw0, snn0
    pub plowable: Vec<[f64; 7]>,

    /// IMPERVIOUS surface: cmin, cmax, tbase, fwfrac, sd0, fw0, snn0
    pub impervious: Vec<[f64; 7]>,

    /// PERVIOUS surface: cmin, cmax, tbase, fwfrac, sd0, fw0, snn0
    pub pervious: Vec<[f64; 7]>,

    /// REMOVAL: dSnow, fOut, fImp, fPerv, fImelt, fToSubcatch, subcatch_name
    pub removal: Vec<[f64; 6]>,
    pub removal_subcatch: Vec<Bytes>,
}

impl SnowpackStore {
    /// names.size() converted to int.
    pub fn count(&self) -> i32 {
        self.names.len() as i32
    }
}

// ============================================================================
// Aquifer definitions (from [AQUIFERS] section)
// ============================================================================

/// SoA storage for aquifer parameter sets.
#[derive(Debug, Clone, Default)]
pub struct AquiferStore {
    pub names: Vec<Bytes>,
    pub porosity: Vec<f64>,
    pub wilting_point: Vec<f64>,
    pub field_capacity: Vec<f64>,
    pub conductivity: Vec<f64>,
    pub conduct_slope: Vec<f64>,
    pub tension_slope: Vec<f64>,
    pub upper_evap: Vec<f64>,
    pub lower_evap: Vec<f64>,
    pub lower_loss: Vec<f64>,
    pub bottom_elev: Vec<f64>,
    pub water_table_elev: Vec<f64>,
    pub upper_moist: Vec<f64>,
    /// ET pattern name (optional)
    pub upper_evap_pat: Vec<Bytes>,
}

impl AquiferStore {
    /// names.size() converted to int.
    pub fn count(&self) -> i32 {
        self.names.len() as i32
    }
}

// ============================================================================
// LID control definitions (from [LID_CONTROLS] section)
// ============================================================================

/// SoA storage for LID control type definitions.
#[derive(Debug, Clone, Default)]
pub struct LidControlStore {
    pub node_layers: Vec<Vec<LidNodeLayer>>,

    pub names: Vec<Bytes>,

    /// LID type code string: BC, RG, GR, IT, PP, RB, RD, VS.
    pub lid_type: Vec<Bytes>,

    /// SURFACE layer params (up to 5 values).
    pub surface: Vec<[f64; 5]>,

    /// SOIL layer params (up to 7 values).
    pub soil: Vec<[f64; 7]>,

    /// PAVEMENT layer params (up to 6 values).
    pub pavement: Vec<[f64; 6]>,

    /// STORAGE layer params (up to 4 values).
    pub storage: Vec<[f64; 4]>,

    /// DRAIN layer params (up to 6 values).
    pub drain: Vec<[f64; 6]>,
    pub drain_curve: Vec<Bytes>,

    /// DRAINMAT layer params (up to 3 values).
    pub drainmat: Vec<[f64; 3]>,

    /// REMOVALS: removals[lid_index] = vector of {pollutant_index, fraction} pairs.
    pub removals: Vec<Vec<(i32, f64)>>,

    /// REMOVALS entries naming something that was not a pollutant when the row was read,
    /// kept by name and bound by msxsurf::resolve.
    pub removals_by_name: Vec<Vec<(Bytes, f64)>>,
}

impl LidControlStore {
    /// names.size() converted to int.
    pub fn count(&self) -> i32 {
        self.names.len() as i32
    }
}

// ============================================================================
// LID usage assignments (from [LID_USAGE] section)
// ============================================================================

/// SoA storage for LID usage assignments to subcatchments.
#[derive(Debug, Clone, Default)]
pub struct LidUsageStore {
    /// Subcatchment index
    pub subcatch_index: Vec<i32>,
    /// Index into LidControlStore
    pub lid_index: Vec<i32>,
    /// Number of replicate units
    pub number: Vec<i32>,
    /// Area of each unit (ft2 or m2)
    pub area: Vec<f64>,
    /// Top width of overland flow (ft or m)
    pub width: Vec<f64>,
    /// Initial saturation (0–1)
    pub init_sat: Vec<f64>,
    /// % of impervious area routed to LID
    pub from_imperv: Vec<f64>,
    /// 1 = route outflow to pervious area
    pub to_perv: Vec<i32>,
    /// Per-unit LID report file (optional).
    pub rpt_file: Vec<FilePathPair>,
    /// Drain-to subcatchment name (optional)
    pub drain_to: Vec<Bytes>,
    /// % of pervious area routed to LID
    pub from_perv: Vec<f64>,

    // Water balance results (ft depth, filled at end of simulation from LIDGroupSoA)
    /// Total inflow depth (ft)
    pub wb_inflow: Vec<f64>,
    /// Total evaporation depth (ft)
    pub wb_evap: Vec<f64>,
    /// Total exfiltration depth (ft)
    pub wb_infil: Vec<f64>,
    /// Total surface outflow depth (ft)
    pub wb_surf_flow: Vec<f64>,
    /// Total drain outflow depth (ft)
    pub wb_drain_flow: Vec<f64>,
    /// Initial stored depth (ft)
    pub wb_init_vol: Vec<f64>,
    /// Final stored depth (ft)
    pub wb_final_vol: Vec<f64>,
}

impl LidUsageStore {
    /// subcatch_index.size() converted to int.
    pub fn count(&self) -> i32 {
        self.subcatch_index.len() as i32
    }

    /// Replace all seven water-balance columns with n zero totals.
    pub fn resize_wb(&mut self, n: i32) {
        let un = n as usize;
        self.wb_inflow.clear();
        self.wb_inflow.resize(un, 0.0);
        self.wb_evap.clear();
        self.wb_evap.resize(un, 0.0);
        self.wb_infil.clear();
        self.wb_infil.resize(un, 0.0);
        self.wb_surf_flow.clear();
        self.wb_surf_flow.resize(un, 0.0);
        self.wb_drain_flow.clear();
        self.wb_drain_flow.resize(un, 0.0);
        self.wb_init_vol.clear();
        self.wb_init_vol.resize(un, 0.0);
        self.wb_final_vol.clear();
        self.wb_final_vol.resize(un, 0.0);
    }
}
