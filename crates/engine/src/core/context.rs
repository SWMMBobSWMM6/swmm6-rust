// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/SimulationContext.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The one context passed to every solver function: options, every object store, every
// name index, mass balance and run statistics. Member for member as upstream, in
// upstream's order. Differences of form only:
//   * Nested structs are top-level types (MassBalance, RoutingStepStats, ...).
//   * StateAccessors' std::function members are Option<Box<dyn Fn...>>.
//   * TwoDModelIO holds raw pointers to 2D objects that a 1D build never creates; it is
//     an empty struct here until 2D is translated.
//   * The C arrays of MaxStats are [MaxStats; MAX_STATS].

#![allow(clippy::all)]

use std::collections::BTreeMap;

use crate::core::cmath::{cmax, cmin};
use crate::core::options::{Bytes, FilePathPair, RoutingModel, SimulationOptions};
use crate::core::spatial_frame::SpatialFrame;
use crate::core::user_flags::UserFlags;
use crate::data::ard_config::ArdConfigData;
use crate::data::bed_zone::BedZoneState;
use crate::data::channel_thermal::{ChannelSpeciesIncrement, ChannelSpeciesRegion, ChannelThermalIncrement, ChannelThermalRegion};
use crate::data::forcing::ForcingData;
use crate::data::gage::GageData;
use crate::data::heat::{HeatConfigData, HeatState};
use crate::data::hydrology::{AquiferStore, LidControlStore, LidUsageStore, SnowpackStore};
use crate::data::inflow::{DwfData, ExtInflowData, PatternData, RDIIAssignData, RDIIDecayData, UnitHydData};
use crate::data::infra::{ControlRuleStore, InletStore, InletUsageStore, StreetStore, TransectStore};
use crate::data::initial_quality::InitialQualityData;
use crate::data::lid_layer_species::LidLayerSpeciesState;
use crate::data::lid_node::LidNodeOutlet;
use crate::data::link::LinkData;
use crate::data::link_subtypes::LinkSubtypes;
use crate::data::name_index::NameIndex;
use crate::data::node::NodeData;
use crate::data::node_coupling_2d::NodeCoupling2DRows;
use crate::data::node_subtypes::NodeSubtypes;
use crate::data::pollutant::PollutantData;
use crate::data::quality::{BuildupData, LanduseData, TreatmentData, WashoffData};
use crate::data::reaction::ReactionData;
use crate::data::species_registry::SpeciesRegistry;
use crate::data::subcatch::SubcatchData;
use crate::data::table::TableData;
use crate::data::water_age::{WaterAgeConfigData, WaterAgeState};
use crate::hydraulics::transect_types::TransectData;
use crate::hydrology::climate_types::ClimateState;

#[derive(Debug, Clone, Default)]
pub struct PluginSpec {
    pub path: Bytes,
    pub init_args: Vec<Bytes>,
}

#[derive(Debug, Clone, Default)]
pub struct ProcessComponentSpec {
    pub id: Bytes,
    pub config_path: Bytes,
    pub args: Vec<(Bytes, Bytes)>,
    pub resolved_config_path: Bytes,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum FileMode {
    #[default]
    None,
    Save,
    Use,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i32)]
pub enum HotstartFormat {
    #[default]
    Auto = 0,
    Legacy = 1,
    Native = 2,
}

#[derive(Debug, Clone, Default)]
pub struct HotstartSaveEntry {
    pub path: FilePathPair,
    pub datetime: f64,
    pub format: HotstartFormat,
}

#[derive(Debug, Clone, Default)]
pub struct FilesSpec {
    pub rainfall_mode: FileMode,
    pub rainfall_path: FilePathPair,
    pub runoff_mode: FileMode,
    pub runoff_path: FilePathPair,
    pub rdii_mode: FileMode,
    pub rdii_path: FilePathPair,
    pub inflows_path: FilePathPair,
    pub outflows_path: FilePathPair,
    pub hotstart_use_path: FilePathPair,
    pub hotstart_saves: Vec<HotstartSaveEntry>,
}

impl FilesSpec {
    pub fn has_any(&self) -> bool {
        self.rainfall_mode != FileMode::None
            || self.runoff_mode != FileMode::None
            || self.rdii_mode != FileMode::None
            || !self.inflows_path.is_empty()
            || !self.outflows_path.is_empty()
            || !self.hotstart_use_path.is_empty()
            || !self.hotstart_saves.is_empty()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i32)]
pub enum EngineState {
    #[default]
    Created = 0,
    Opened = 1,
    Initialized = 2,
    Running = 3,
    Paused = 4,
    Ended = 5,
    Reported = 6,
    Closed = 7,
    ErrorState = 8,
    Building = 9,
}

/// Run-time hooks other subsystems install so hotstart and the API can reach their state.
#[derive(Default)]
pub struct StateAccessors {
    pub has_fv_retained_dry_mass: Option<Box<dyn Fn() -> bool>>,
    pub get_fv_dry_state: Option<Box<dyn Fn() -> Bytes>>,
    pub set_fv_dry_state: Option<Box<dyn Fn(&[u8], bool) -> Bytes>>,
    pub get_fv_retained_pollutant_mass: Option<Box<dyn Fn() -> Vec<f64>>>,
    pub get_fv_retained_species_mass: Option<Box<dyn Fn() -> Vec<f64>>>,
    pub get_channel_thermal_regions: Option<Box<dyn Fn() -> Vec<ChannelThermalRegion>>>,
    pub apply_channel_thermal_increments: Option<Box<dyn Fn(&[ChannelThermalIncrement]) -> Bytes>>,
    pub get_channel_species_regions: Option<Box<dyn Fn() -> Vec<ChannelSpeciesRegion>>>,
    pub apply_channel_species_increments: Option<Box<dyn Fn(&[ChannelSpeciesIncrement]) -> Bytes>>,
    pub get_snow_state: Option<Box<dyn Fn() -> Bytes>>,
    pub set_snow_state: Option<Box<dyn Fn(&[u8], bool) -> Bytes>>,
    pub get_surface_state: Option<Box<dyn Fn() -> Bytes>>,
    pub set_surface_state: Option<Box<dyn Fn(&[u8], bool) -> Bytes>>,
    /// (subcatch_index, out_model, out_infil) -> found
    pub get_infil_state: Option<Box<dyn Fn(i32, &mut i32, &mut [f64]) -> bool>>,
    pub set_infil_state: Option<Box<dyn Fn(i32, i32, &[f64]) -> bool>>,
    pub get_gw_state: Option<Box<dyn Fn(i32, &mut f64, &mut f64) -> bool>>,
    pub set_gw_state: Option<Box<dyn Fn(i32, f64, f64) -> bool>>,
}

impl StateAccessors {
    pub fn can_read(&self) -> bool {
        self.get_infil_state.is_some() || self.get_gw_state.is_some()
    }
    pub fn can_write(&self) -> bool {
        self.set_infil_state.is_some() || self.set_gw_state.is_some()
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct Event {
    pub start: f64,
    pub end: f64,
}

/// Upstream: raw pointers to the 2D model's objects, all null in a 1D build.
#[derive(Debug, Clone, Default)]
pub struct TwoDModelIO {}

#[derive(Debug, Clone, Default)]
pub struct OutfallExports {
    pub water_m3: Vec<f64>,
    pub reverse_m3: Vec<f64>,
    pub heat_j: Vec<f64>,
    pub pollutant_amounts: Vec<f64>,
}

#[derive(Debug, Clone, Default)]
pub struct MassBalance {
    pub runoff_rainfall: f64,
    pub runoff_runon: f64,
    pub runoff_evap: f64,
    pub runoff_infil: f64,
    pub runoff_runoff: f64,
    pub runoff_lid_drain: f64,
    pub runoff_snowremov: f64,
    pub runoff_init_store: f64,
    pub runoff_final_store: f64,
    pub runoff_init_snow: f64,
    pub runoff_final_snow: f64,
    pub rdii_rainfall: f64,
    pub rdii_produced: f64,
    pub routing_dry_weather: f64,
    pub routing_wet_weather: f64,
    pub routing_gw_inflow: f64,
    pub routing_link_gw_inflow: f64,
    pub routing_rdii: f64,
    pub routing_external: f64,
    pub routing_flooding: f64,
    pub routing_coupling_out: f64,
    pub routing_outflow: f64,
    pub routing_evap_loss: f64,
    pub routing_seep_loss: f64,
    pub routing_init_storage: f64,
    pub routing_final_storage: f64,
    pub routing_report: [f64; 11],
    pub routing_forcing_inflow: f64,
    pub routing_forcing_qual_inflow: Vec<f64>,
    pub gw_infil: f64,
    pub gw_upper_evap: f64,
    pub gw_lower_evap: f64,
    pub gw_lower_perc: f64,
    pub gw_lateral_flow: f64,
    pub gw_init_storage: f64,
    pub gw_final_storage: f64,
    pub gw_infil_2d_recharge: f64,
    pub step_flooding: f64,
    pub step_outflow: f64,
    pub step_dw_inflow: f64,
    pub step_gw_inflow: f64,
    pub step_rdii_inflow: f64,
    pub step_ext_inflow: f64,
    pub last_step_flow_error: f64,
    pub step_ext_withdrawal: f64,
    pub qual_init_buildup: Vec<f64>,
    pub qual_final_buildup: Vec<f64>,
    pub qual_ponded_lost: Vec<f64>,
    pub qual_surface_buildup: Vec<f64>,
    pub qual_wet_deposition: Vec<f64>,
    pub qual_sweeping: Vec<f64>,
    pub qual_bmp_removal: Vec<f64>,
    pub qual_infil_loss: Vec<f64>,
    pub qual_runoff_load: Vec<f64>,
    pub qual_routing_wet: Vec<f64>,
    pub qual_routing_outflow: Vec<f64>,
    pub qual_routing_flood: Vec<f64>,
    pub qual_routing_init: Vec<f64>,
    pub qual_routing_final: Vec<f64>,
    pub qual_routing_reacted: Vec<f64>,
    pub qual_routing_ii_in: Vec<f64>,
    pub qual_routing_dw_in: Vec<f64>,
    pub qual_routing_gw_in: Vec<f64>,
    pub qual_routing_ex_in: Vec<f64>,
    pub qual_routing_seep: Vec<f64>,
    pub qual_routing_evap: Vec<f64>,
    pub qual_routing_final_dry: Vec<f64>,
}

impl MassBalance {
    pub fn resize_quality(&mut self, n_pollutants: i32) {
        let np = n_pollutants as usize;
        for v in [
            &mut self.qual_init_buildup,
            &mut self.qual_final_buildup,
            &mut self.qual_ponded_lost,
            &mut self.qual_surface_buildup,
            &mut self.qual_wet_deposition,
            &mut self.qual_sweeping,
            &mut self.qual_bmp_removal,
            &mut self.qual_infil_loss,
            &mut self.qual_runoff_load,
            &mut self.qual_routing_wet,
            &mut self.qual_routing_outflow,
            &mut self.qual_routing_flood,
            &mut self.qual_routing_init,
            &mut self.qual_routing_final,
            &mut self.qual_routing_reacted,
            &mut self.qual_routing_ii_in,
            &mut self.qual_routing_dw_in,
            &mut self.qual_routing_gw_in,
            &mut self.qual_routing_ex_in,
            &mut self.qual_routing_seep,
            &mut self.qual_routing_evap,
            &mut self.qual_routing_final_dry,
            &mut self.routing_forcing_qual_inflow,
        ] {
            v.clear();
            v.resize(np, 0.0);
        }
    }

    /// Every scalar back to zero; the per-pollutant vectors keep their length and are
    /// zeroed, except qual_init_buildup, which upstream leaves holding its values.
    pub fn reset(&mut self) {
        let keep = MassBalance {
            qual_init_buildup: std::mem::take(&mut self.qual_init_buildup),
            qual_final_buildup: std::mem::take(&mut self.qual_final_buildup),
            qual_ponded_lost: std::mem::take(&mut self.qual_ponded_lost),
            qual_surface_buildup: std::mem::take(&mut self.qual_surface_buildup),
            qual_wet_deposition: std::mem::take(&mut self.qual_wet_deposition),
            qual_sweeping: std::mem::take(&mut self.qual_sweeping),
            qual_bmp_removal: std::mem::take(&mut self.qual_bmp_removal),
            qual_infil_loss: std::mem::take(&mut self.qual_infil_loss),
            qual_runoff_load: std::mem::take(&mut self.qual_runoff_load),
            qual_routing_wet: std::mem::take(&mut self.qual_routing_wet),
            qual_routing_outflow: std::mem::take(&mut self.qual_routing_outflow),
            qual_routing_flood: std::mem::take(&mut self.qual_routing_flood),
            qual_routing_init: std::mem::take(&mut self.qual_routing_init),
            qual_routing_final: std::mem::take(&mut self.qual_routing_final),
            qual_routing_reacted: std::mem::take(&mut self.qual_routing_reacted),
            qual_routing_ii_in: std::mem::take(&mut self.qual_routing_ii_in),
            qual_routing_dw_in: std::mem::take(&mut self.qual_routing_dw_in),
            qual_routing_gw_in: std::mem::take(&mut self.qual_routing_gw_in),
            qual_routing_ex_in: std::mem::take(&mut self.qual_routing_ex_in),
            qual_routing_seep: std::mem::take(&mut self.qual_routing_seep),
            qual_routing_evap: std::mem::take(&mut self.qual_routing_evap),
            qual_routing_final_dry: std::mem::take(&mut self.qual_routing_final_dry),
            routing_forcing_qual_inflow: std::mem::take(&mut self.routing_forcing_qual_inflow),
            ..MassBalance::default()
        };
        *self = keep;
        for v in [
            &mut self.qual_final_buildup,
            &mut self.qual_ponded_lost,
            &mut self.qual_surface_buildup,
            &mut self.qual_wet_deposition,
            &mut self.qual_sweeping,
            &mut self.qual_bmp_removal,
            &mut self.qual_infil_loss,
            &mut self.qual_runoff_load,
            &mut self.qual_routing_wet,
            &mut self.qual_routing_outflow,
            &mut self.qual_routing_flood,
            &mut self.qual_routing_init,
            &mut self.qual_routing_final,
            &mut self.qual_routing_reacted,
            &mut self.qual_routing_ii_in,
            &mut self.qual_routing_dw_in,
            &mut self.qual_routing_gw_in,
            &mut self.qual_routing_ex_in,
            &mut self.qual_routing_seep,
            &mut self.qual_routing_evap,
            &mut self.qual_routing_final_dry,
            &mut self.routing_forcing_qual_inflow,
        ] {
            v.fill(0.0);
        }
    }

    pub fn runoff_error(&self) -> f64 {
        let total_in = self.runoff_rainfall + self.runoff_runon + self.runoff_init_store + self.runoff_init_snow;
        let total_out = self.runoff_evap
            + self.runoff_infil
            + self.runoff_runoff
            + self.runoff_lid_drain
            + self.runoff_final_store
            + self.runoff_snowremov
            + self.runoff_final_snow;
        if total_in > 0.0 {
            (total_in - total_out) / total_in
        } else {
            0.0
        }
    }

    pub fn routing_error(&self) -> f64 {
        let total_in = self.routing_dry_weather
            + self.routing_wet_weather
            + self.routing_gw_inflow
            + self.routing_link_gw_inflow
            + self.routing_rdii
            + self.routing_external
            + self.routing_init_storage;
        let total_out = self.routing_flooding
            + self.routing_coupling_out
            + self.routing_outflow
            + self.routing_evap_loss
            + self.routing_seep_loss
            + self.routing_final_storage;
        if total_in > 0.0 {
            (total_in - total_out) / total_in
        } else {
            0.0
        }
    }

    pub fn gw_error(&self) -> f64 {
        let total_in = self.gw_infil + self.gw_init_storage;
        let total_out = self.gw_upper_evap + self.gw_lower_evap + self.gw_lower_perc + self.gw_lateral_flow + self.gw_final_storage;
        if total_in > 0.0 {
            (total_in - total_out) / total_in
        } else {
            0.0
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct NegativeSourceStats {
    pub clamp_events: i64,
    pub shortfall_mass: f64,
    pub age_clamp_events: i64,
    pub first_node: i32,
    pub first_clamp_recorded: bool,
    pub api_warned: bool,
}

impl Default for NegativeSourceStats {
    fn default() -> Self {
        NegativeSourceStats {
            clamp_events: 0,
            shortfall_mass: 0.0,
            age_clamp_events: 0,
            first_node: -1,
            first_clamp_recorded: false,
            api_warned: false,
        }
    }
}

impl NegativeSourceStats {
    pub fn reset(&mut self) {
        *self = NegativeSourceStats::default();
    }
}

#[derive(Debug, Clone, Default)]
pub struct QualityOnlyInflow {
    pub vol: Vec<f64>,
    pub mass: Vec<f64>,
    pub msx: Vec<f64>,
    pub age_vol: Vec<f64>,
    pub temp_vol: Vec<f64>,
    pub any: bool,
}

impl QualityOnlyInflow {
    pub fn ensure(&mut self, nn: i32, np: i32, ns: i32) {
        let un = if nn > 0 { nn } else { 0 } as usize;
        if self.vol.len() != un {
            for v in [&mut self.vol, &mut self.age_vol, &mut self.temp_vol] {
                v.clear();
                v.resize(un, 0.0);
            }
        }
        let npu = if np > 0 { np } else { 0 } as usize;
        if self.mass.len() != un * npu {
            self.mass.clear();
            self.mass.resize(un * npu, 0.0);
        }
        let nsu = if ns > 0 { ns } else { 0 } as usize;
        if self.msx.len() != un * nsu {
            self.msx.clear();
            self.msx.resize(un * nsu, 0.0);
        }
    }
    pub fn clear(&mut self) {
        for v in [&mut self.vol, &mut self.mass, &mut self.msx, &mut self.age_vol, &mut self.temp_vol] {
            v.fill(0.0);
        }
        self.any = false;
    }
}

#[derive(Debug, Clone, Copy)]
pub struct MassBalance2D {
    pub init_storage: f64,
    pub final_storage: f64,
    pub external_in: f64,
    pub external_out: f64,
    pub rainfall_in: f64,
    pub coupling_1d_to_2d_in: f64,
    pub coupling_2d_to_1d_out: f64,
    pub outfall_in: f64,
    pub outfall_out: f64,
    pub boundary_in: f64,
    pub boundary_out: f64,
    pub evap_out: f64,
    pub infil_out: f64,
    pub infil_to_aquifer: f64,
    pub infil_aquifer_pending: f64,
    pub aquifer_in: f64,
    pub active: bool,
    pub solver_nsteps: i64,
    pub solver_nrhs: i64,
    pub solver_avg_h: f64,
    pub solver_last_h: f64,
    pub solver_active_min: f64,
    pub solver_active_mean: f64,
    pub solver_active_max: f64,
    pub solver_tier_cells: [i64; 8],
    pub solver_n_tiers: i32,
}

impl Default for MassBalance2D {
    fn default() -> Self {
        MassBalance2D {
            init_storage: 0.0,
            final_storage: 0.0,
            external_in: 0.0,
            external_out: 0.0,
            rainfall_in: 0.0,
            coupling_1d_to_2d_in: 0.0,
            coupling_2d_to_1d_out: 0.0,
            outfall_in: 0.0,
            outfall_out: 0.0,
            boundary_in: 0.0,
            boundary_out: 0.0,
            evap_out: 0.0,
            infil_out: 0.0,
            infil_to_aquifer: 0.0,
            infil_aquifer_pending: 0.0,
            aquifer_in: 0.0,
            active: false,
            solver_nsteps: -1,
            solver_nrhs: 0,
            solver_avg_h: 0.0,
            solver_last_h: 0.0,
            solver_active_min: -1.0,
            solver_active_mean: -1.0,
            solver_active_max: -1.0,
            solver_tier_cells: [0; 8],
            solver_n_tiers: 0,
        }
    }
}

impl MassBalance2D {
    pub fn error(&self) -> f64 {
        let total_in = self.rainfall_in
            + self.coupling_1d_to_2d_in
            + self.outfall_in
            + self.boundary_in
            + self.aquifer_in
            + self.init_storage
            + self.external_in;
        let total_out = self.coupling_2d_to_1d_out
            + self.outfall_out
            + self.boundary_out
            + self.evap_out
            + self.infil_out
            + self.final_storage
            + self.external_out;
        if total_in > 0.0 {
            (total_in - total_out) / total_in
        } else {
            0.0
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaxStats {
    pub obj_type: i32,
    pub index: i32,
    pub value: f64,
}

impl Default for MaxStats {
    fn default() -> Self {
        MaxStats { obj_type: -1, index: -1, value: 0.0 }
    }
}

pub const MAX_STATS: usize = 5;
pub const N_TIME_BINS: usize = 5;

#[derive(Debug, Clone, Copy)]
pub struct RoutingStepStats {
    pub min_step: f64,
    pub max_step: f64,
    pub sum_step: f64,
    pub n_steps: i64,
    pub steady_pct: f64,
    pub steady_time: f64,
    pub report_steps: i64,
    pub report_time: f64,
    pub max_outfall_flow: f64,
    pub step_counts: [i64; N_TIME_BINS + 1],
    pub step_intervals: [f64; N_TIME_BINS + 1],
    pub n_non_converged: i64,
    pub sum_iterations: f64,
    pub max_courant: f64,
    pub slot_peak_share: f64,
    pub slot_time_above_s: f64,
    pub fv_nsteps: i64,
    pub fv_nflux: i64,
    pub fv_avg_h: f64,
    pub fv_last_h: f64,
    pub fv_min_h: f64,
    pub fv_active_min: f64,
    pub fv_active_mean: f64,
    pub fv_active_max: f64,
    pub fv_tier_cells: [i64; 8],
    pub fv_n_tiers: i32,
    pub fv_dt_argmin_pressurized: i64,
    pub fv_dt_argmin_band: i64,
    pub fv_dt_argmin_free: i64,
    pub fv_dt_argmin_node: i64,
    pub fv_macro_cycles: i64,
    pub fv_macro_rejected: i64,
    pub lstep_mode: i32,
}

impl Default for RoutingStepStats {
    fn default() -> Self {
        RoutingStepStats {
            min_step: 1.0e30,
            max_step: 0.0,
            sum_step: 0.0,
            n_steps: 0,
            steady_pct: 0.0,
            steady_time: 0.0,
            report_steps: 0,
            report_time: 0.0,
            max_outfall_flow: 0.0,
            step_counts: [0; N_TIME_BINS + 1],
            step_intervals: [0.0; N_TIME_BINS + 1],
            n_non_converged: 0,
            sum_iterations: 0.0,
            max_courant: 0.0,
            slot_peak_share: 0.0,
            slot_time_above_s: 0.0,
            fv_nsteps: -1,
            fv_nflux: 0,
            fv_avg_h: 0.0,
            fv_last_h: 0.0,
            fv_min_h: 0.0,
            fv_active_min: -1.0,
            fv_active_mean: -1.0,
            fv_active_max: -1.0,
            fv_tier_cells: [0; 8],
            fv_n_tiers: 0,
            fv_dt_argmin_pressurized: 0,
            fv_dt_argmin_band: 0,
            fv_dt_argmin_free: 0,
            fv_dt_argmin_node: 0,
            fv_macro_cycles: 0,
            fv_macro_rejected: 0,
            lstep_mode: 0,
        }
    }
}

impl RoutingStepStats {
    pub fn update(&mut self, dt: f64) {
        self.min_step = cmin(self.min_step, dt);
        self.max_step = cmax(self.max_step, dt);
        self.sum_step += dt;
        self.n_steps += 1;
    }

    pub fn record_step(&mut self, dt: f64, steady: bool, first_step: bool) {
        if steady {
            self.steady_time += dt;
            return;
        }
        if !first_step {
            self.min_step = cmin(self.min_step, dt);
            for i in 0..N_TIME_BINS {
                if dt >= self.step_intervals[i + 1] {
                    self.step_counts[i] += 1;
                    break;
                }
            }
        }
        self.max_step = cmax(self.max_step, dt);
        self.sum_step += dt;
        self.n_steps += 1;
    }

    pub fn update_iterations(&mut self, iters: i32, converged: bool) {
        self.sum_iterations += iters as f64;
        if !converged {
            self.n_non_converged += 1;
        }
    }

    pub fn init_histogram(&mut self, route_step: f64, mut min_route_step: f64) {
        if route_step <= 0.0 {
            return;
        }
        if min_route_step <= 0.0 {
            min_route_step = route_step;
        }
        let log_hi = route_step.log10(); // LIBM
        let log_lo = min_route_step.log10(); // LIBM
        let delta = (log_hi - log_lo) / N_TIME_BINS as f64;
        self.step_intervals[0] = route_step;
        for i in 1..=N_TIME_BINS {
            self.step_intervals[i] = 10.0f64.powf(log_hi - i as f64 * delta); // LIBM
        }
        self.step_intervals[N_TIME_BINS] = min_route_step;
    }

    pub fn build_histogram(&mut self) {
        if self.n_steps == 0 || self.max_step <= 0.0 {
            return;
        }
        let hi = self.max_step;
        let mut lo = if self.min_step < 1.0e30 { self.min_step } else { 0.0 };
        if lo <= 0.0 {
            lo = hi;
        }
        self.init_histogram(hi, lo);
    }

    pub fn record_step_bin(&mut self, dt: f64) {
        for i in 0..N_TIME_BINS {
            if dt >= self.step_intervals[i + 1] {
                self.step_counts[i] += 1;
                return;
            }
        }
        self.step_counts[N_TIME_BINS - 1] += 1;
    }

    pub fn avg_step(&self) -> f64 {
        if self.n_steps > 0 {
            self.sum_step / self.n_steps as f64
        } else {
            0.0
        }
    }

    pub fn pct_non_converged(&self) -> f64 {
        if self.n_steps > 0 {
            100.0 * self.n_non_converged as f64 / self.n_steps as f64
        } else {
            0.0
        }
    }

    pub fn computed_avg_iterations(&self) -> f64 {
        if self.n_steps > 0 {
            self.sum_iterations / self.n_steps as f64
        } else {
            0.0
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct VJDiag {
    pub node_idx: Vec<i32>,
    pub up_link: Vec<i32>,
    pub dn_link: Vec<i32>,
    pub resid_max: Vec<f64>,
    pub resid_sum: Vec<f64>,
    pub resid_n: Vec<i64>,
}

impl VJDiag {
    pub fn clear(&mut self) {
        self.node_idx.clear();
        self.up_link.clear();
        self.dn_link.clear();
        self.resid_max.clear();
        self.resid_sum.clear();
        self.resid_n.clear();
    }
}

#[derive(Debug, Clone, Default)]
pub struct InletDiag {
    pub host_node: Vec<i32>,
    pub up_link: Vec<i32>,
    pub is_sag: Vec<u8>,
    pub num_inlets: Vec<i32>,
    pub flow_periods: Vec<i32>,
    pub capture_periods: Vec<i32>,
    pub backflow_periods: Vec<i32>,
    pub peak_flow: Vec<f64>,
    pub peak_flow_capture: Vec<f64>,
    pub avg_flow_capture: Vec<f64>,
    pub bypass_freq: Vec<f64>,
}

impl InletDiag {
    pub fn count(&self) -> i32 {
        self.host_node.len() as i32
    }
    pub fn resize(&mut self, n: i32) {
        let un = n as usize;
        macro_rules! assign {
            ($v:expr, $d:expr) => {{
                $v.clear();
                $v.resize(un, $d);
            }};
        }
        assign!(self.host_node, -1);
        assign!(self.up_link, -1);
        assign!(self.is_sag, 0u8);
        assign!(self.num_inlets, 1);
        assign!(self.flow_periods, 0);
        assign!(self.capture_periods, 0);
        assign!(self.backflow_periods, 0);
        assign!(self.peak_flow, 0.0);
        assign!(self.peak_flow_capture, 0.0);
        assign!(self.avg_flow_capture, 0.0);
        assign!(self.bypass_freq, 0.0);
    }
    pub fn clear(&mut self) {
        self.host_node.clear();
        self.up_link.clear();
        self.is_sag.clear();
        self.num_inlets.clear();
        self.flow_periods.clear();
        self.capture_periods.clear();
        self.backflow_periods.clear();
        self.peak_flow.clear();
        self.peak_flow_capture.clear();
        self.avg_flow_capture.clear();
        self.bypass_freq.clear();
    }
}

#[derive(Debug, Clone, Copy, Default)]
pub struct ControlLogEntry {
    pub link_idx: i32,
    pub rule_idx: i32,
    pub new_setting: f64,
    pub date: f64,
}

pub const K_MAX_CONTROL_LOG: usize = 1_000_000;

pub struct SimulationContext {
    pub state: EngineState,
    pub wall_start: i64,
    pub title_notes: Vec<Bytes>,
    pub passthrough_sections: BTreeMap<Bytes, Vec<Bytes>>,
    pub options: SimulationOptions,
    pub climate_state: ClimateState,
    pub current_time: f64,
    pub current_date: f64,
    pub elapsed_ms: f64,
    pub old_elapsed_ms: f64,
    pub next_report_ms: f64,
    pub runoff_interp_f: f64,
    pub dt_controls_remaining: f64,
    pub nodes: NodeData,
    pub node_subtypes: NodeSubtypes,
    pub links: LinkData,
    pub link_subtypes: LinkSubtypes,
    pub subcatches: SubcatchData,
    pub gages: GageData,
    pub pollutants: PollutantData,
    pub species_registry: SpeciesRegistry,
    pub reactions: ReactionData,
    pub ard_config: ArdConfigData,
    pub water_age_config: WaterAgeConfigData,
    pub water_age_state: WaterAgeState,
    pub lid_layer_state: LidLayerSpeciesState,
    pub heat_config: HeatConfigData,
    pub heat_state: HeatState,
    pub bed_state: BedZoneState,
    pub reported_species_names: Vec<Bytes>,
    pub tables: TableData,
    pub node_names: NameIndex,
    pub link_names: NameIndex,
    pub subcatch_names: NameIndex,
    pub gage_names: NameIndex,
    pub pollutant_names: NameIndex,
    pub landuse_names: NameIndex,
    pub landuses: LanduseData,
    pub buildup: BuildupData,
    pub washoff: WashoffData,
    pub treatment: TreatmentData,
    pub initial_quality: InitialQualityData,
    pub ext_inflows: ExtInflowData,
    pub dwf_inflows: DwfData,
    pub rdii_assigns: RDIIAssignData,
    pub unit_hyds: UnitHydData,
    pub rdii_decay: RDIIDecayData,
    pub patterns: PatternData,
    pub transects: TransectStore,
    pub transect_tables: Vec<TransectData>,
    pub xsect_generation: u64,
    pub gpkg_units_internal: bool,
    pub streets: StreetStore,
    pub inlets: InletStore,
    pub inlet_usages: InletUsageStore,
    pub control_rules: ControlRuleStore,
    pub events: Vec<Event>,
    pub adjust_temp: [f64; 12],
    pub adjust_evap: [f64; 12],
    pub adjust_rain: [f64; 12],
    pub adjust_hydcon: [f64; 12],
    pub subcatch_n_perv_pattern: Vec<i32>,
    pub subcatch_d_store_pattern: Vec<i32>,
    pub subcatch_infil_pattern: Vec<i32>,
    pub base_n_perv: Vec<f64>,
    pub base_ds_perv: Vec<f64>,
    pub has_subcatch_adj_patterns: bool,
    pub snowpacks: SnowpackStore,
    pub snowpack_names: NameIndex,
    pub aquifers: AquiferStore,
    pub aquifer_names: NameIndex,
    pub lid_controls: LidControlStore,
    pub lid_node_outlets: Vec<LidNodeOutlet>,
    pub lid_names: NameIndex,
    pub lid_usage: LidUsageStore,
    pub pending_gw_nodes: Vec<(i32, Bytes)>,
    pub pending_gw_aquifers: Vec<(i32, Bytes)>,
    pub pending_link_nodes: Vec<(i32, (Bytes, Bytes))>,
    pub spatial: SpatialFrame,
    pub user_flags: UserFlags,
    pub forcing: ForcingData,
    pub plugin_specs: Vec<PluginSpec>,
    pub process_component_specs: Vec<ProcessComponentSpec>,
    pub embedded_component_sections: Vec<(Bytes, Vec<Bytes>)>,
    pub files: FilesSpec,
    pub state_accessors: StateAccessors,
    pub twod_io: TwoDModelIO,
    pub deferred_section_rows: Vec<(Bytes, Bytes)>,
    pub unknown_section_headers: Vec<(Bytes, i64)>,
    pub error_code: i32,
    pub warning_code: i32,
    pub error_message: Bytes,
    pub warnings: Vec<Bytes>,
    pub errors: Vec<Bytes>,
    pub node_coupling_2d: NodeCoupling2DRows,
    pub coupled_node: Vec<u8>,
    pub coupled_pond_area: Vec<f64>,
    pub coupled_pond_band: Vec<f64>,
    pub coupling_delivery_remaining: f64,
    pub outfall_exports: OutfallExports,
    pub mass_balance: MassBalance,
    pub negsrc: NegativeSourceStats,
    pub surface_outfall_link_limit: Vec<f64>,
    pub surface_outfall_quality: QualityOnlyInflow,
    pub link_gain_quality: QualityOnlyInflow,
    pub runtime_coupling_used: bool,
    pub mass_balance_2d: MassBalance2D,
    pub max_courant_crit: [MaxStats; MAX_STATS],
    pub max_flow_turns: [MaxStats; MAX_STATS],
    pub max_non_converged: [MaxStats; MAX_STATS],
    pub max_mass_bal_errs: [MaxStats; MAX_STATS],
    pub routing_stats: RoutingStepStats,
    pub vj_diag: VJDiag,
    pub inlet_diag: InletDiag,
    pub control_log: Vec<ControlLogEntry>,
    pub control_rule_names: Vec<Bytes>,
    pub control_log_dropped: usize,
    pub inp_file_path: Bytes,
}

impl Default for SimulationContext {
    fn default() -> Self {
        SimulationContext {
            state: EngineState::Created,
            wall_start: 0,
            title_notes: Vec::new(),
            passthrough_sections: BTreeMap::new(),
            options: SimulationOptions::default(),
            climate_state: ClimateState::default(),
            current_time: 0.0,
            current_date: 0.0,
            elapsed_ms: 0.0,
            old_elapsed_ms: 0.0,
            next_report_ms: 0.0,
            runoff_interp_f: 1.0,
            dt_controls_remaining: 0.0,
            nodes: NodeData::default(),
            node_subtypes: NodeSubtypes::default(),
            links: LinkData::default(),
            link_subtypes: LinkSubtypes::default(),
            subcatches: SubcatchData::default(),
            gages: GageData::default(),
            pollutants: PollutantData::default(),
            species_registry: SpeciesRegistry::default(),
            reactions: ReactionData::default(),
            ard_config: ArdConfigData::default(),
            water_age_config: WaterAgeConfigData::default(),
            water_age_state: WaterAgeState::default(),
            lid_layer_state: LidLayerSpeciesState::default(),
            heat_config: HeatConfigData::default(),
            heat_state: HeatState::default(),
            bed_state: BedZoneState::default(),
            reported_species_names: Vec::new(),
            tables: TableData::default(),
            node_names: NameIndex::default(),
            link_names: NameIndex::default(),
            subcatch_names: NameIndex::default(),
            gage_names: NameIndex::default(),
            pollutant_names: NameIndex::default(),
            landuse_names: NameIndex::default(),
            landuses: LanduseData::default(),
            buildup: BuildupData::default(),
            washoff: WashoffData::default(),
            treatment: TreatmentData::default(),
            initial_quality: InitialQualityData::default(),
            ext_inflows: ExtInflowData::default(),
            dwf_inflows: DwfData::default(),
            rdii_assigns: RDIIAssignData::default(),
            unit_hyds: UnitHydData::default(),
            rdii_decay: RDIIDecayData::default(),
            patterns: PatternData::default(),
            transects: TransectStore::default(),
            transect_tables: Vec::new(),
            xsect_generation: 0,
            gpkg_units_internal: false,
            streets: StreetStore::default(),
            inlets: InletStore::default(),
            inlet_usages: InletUsageStore::default(),
            control_rules: ControlRuleStore::default(),
            events: Vec::new(),
            adjust_temp: [0.0; 12],
            adjust_evap: [0.0; 12],
            adjust_rain: [1.0; 12],
            adjust_hydcon: [1.0; 12],
            subcatch_n_perv_pattern: Vec::new(),
            subcatch_d_store_pattern: Vec::new(),
            subcatch_infil_pattern: Vec::new(),
            base_n_perv: Vec::new(),
            base_ds_perv: Vec::new(),
            has_subcatch_adj_patterns: false,
            snowpacks: SnowpackStore::default(),
            snowpack_names: NameIndex::default(),
            aquifers: AquiferStore::default(),
            aquifer_names: NameIndex::default(),
            lid_controls: LidControlStore::default(),
            lid_node_outlets: Vec::new(),
            lid_names: NameIndex::default(),
            lid_usage: LidUsageStore::default(),
            pending_gw_nodes: Vec::new(),
            pending_gw_aquifers: Vec::new(),
            pending_link_nodes: Vec::new(),
            spatial: SpatialFrame::default(),
            user_flags: UserFlags::default(),
            forcing: ForcingData::default(),
            plugin_specs: Vec::new(),
            process_component_specs: Vec::new(),
            embedded_component_sections: Vec::new(),
            files: FilesSpec::default(),
            state_accessors: StateAccessors::default(),
            twod_io: TwoDModelIO::default(),
            deferred_section_rows: Vec::new(),
            unknown_section_headers: Vec::new(),
            error_code: 0,
            warning_code: 0,
            error_message: Vec::new(),
            warnings: Vec::new(),
            errors: Vec::new(),
            node_coupling_2d: NodeCoupling2DRows::default(),
            coupled_node: Vec::new(),
            coupled_pond_area: Vec::new(),
            coupled_pond_band: Vec::new(),
            coupling_delivery_remaining: 0.0,
            outfall_exports: OutfallExports::default(),
            mass_balance: MassBalance::default(),
            negsrc: NegativeSourceStats::default(),
            surface_outfall_link_limit: Vec::new(),
            surface_outfall_quality: QualityOnlyInflow::default(),
            link_gain_quality: QualityOnlyInflow::default(),
            runtime_coupling_used: false,
            mass_balance_2d: MassBalance2D::default(),
            max_courant_crit: [MaxStats::default(); MAX_STATS],
            max_flow_turns: [MaxStats::default(); MAX_STATS],
            max_non_converged: [MaxStats::default(); MAX_STATS],
            max_mass_bal_errs: [MaxStats::default(); MAX_STATS],
            routing_stats: RoutingStepStats::default(),
            vj_diag: VJDiag::default(),
            inlet_diag: InletDiag::default(),
            control_log: Vec::new(),
            control_rule_names: Vec::new(),
            control_log_dropped: 0,
            inp_file_path: Vec::new(),
        }
    }
}

impl SimulationContext {
    pub fn refresh_reported_species(&mut self) {
        self.reported_species_names.clear();
        for p in 0..self.n_pollutants() {
            self.reported_species_names.push(self.pollutant_names.name_of(p).to_vec());
        }
        if self.options.water_age {
            self.reported_species_names.push(b"__WATER_AGE__".to_vec());
        }
        if self.options.heat_transport {
            self.reported_species_names.push(b"__TEMPERATURE__".to_vec());
        }
        self.reactions.reported_species.clear();
        for m in 0..self.reactions.n_species() {
            if !self.reactions.report_enabled.is_empty() && self.reactions.report_enabled[m as usize] == 0 {
                continue;
            }
            self.reactions.reported_species.push(m);
            let name = self.reactions.species_name[m as usize].clone();
            self.reported_species_names.push(name);
        }
    }

    pub fn n_reported_species(&self) -> i32 {
        self.reported_species_names.len() as i32
    }

    pub fn log_control_action(&mut self, e: ControlLogEntry) {
        if self.control_log.len() >= K_MAX_CONTROL_LOG {
            self.control_log_dropped += 1;
            return;
        }
        self.control_log.push(e);
    }

    /// Insert `value` into a top-five list by magnitude, pushing smaller entries down.
    pub fn update_max_stats(arr: &mut [MaxStats; MAX_STATS], obj_type: i32, idx: i32, value: f64) {
        let mut candidate = MaxStats { obj_type, index: idx, value };
        for k in 0..MAX_STATS {
            if candidate.value.abs() > arr[k].value.abs() {
                std::mem::swap(&mut arr[k], &mut candidate);
            }
        }
    }

    pub fn finalize_max_stats(&mut self) {
        for k in 0..MAX_STATS {
            self.max_courant_crit[k] = MaxStats { obj_type: 0, index: -1, value: -1.0 };
            self.max_flow_turns[k] = MaxStats { obj_type: 1, index: -1, value: -1.0 };
            self.max_non_converged[k] = MaxStats { obj_type: 0, index: -1, value: 0.0 };
        }
        let rpt_steps = self.routing_stats.report_steps;
        if rpt_steps > 2 {
            let z = 100.0 / (2.0 / 3.0 * (rpt_steps - 2) as f64);
            for j in 0..self.n_links() {
                let x = self.links.stat_flow_turns[j as usize] as f64 * z;
                Self::update_max_stats(&mut self.max_flow_turns, 1, j, x);
            }
        }
        let step_count = self.routing_stats.n_steps;
        let dw = self.options.routing_model == RoutingModel::Dynwave;
        if dw {
            for j in 0..self.n_nodes() {
                let v = self.nodes.stat_non_converged_count[j as usize] as f64 / step_count as f64;
                Self::update_max_stats(&mut self.max_non_converged, 0, j, v);
            }
        }
        if !dw || self.options.variable_step == 0.0 || step_count == 0 {
            return;
        }
        let inv_steps = 1.0 / step_count as f64;
        for j in 0..self.n_nodes() {
            let x = self.nodes.stat_time_courant_critical[j as usize] * inv_steps;
            Self::update_max_stats(&mut self.max_courant_crit, 0, j, 100.0 * x);
        }
        for j in 0..self.n_links() {
            let x = self.links.stat_time_courant_critical[j as usize] * inv_steps;
            Self::update_max_stats(&mut self.max_courant_crit, 1, j, 100.0 * x);
        }
    }

    /// Upstream's reset(): only the members it names, to the values it names.
    pub fn reset(&mut self) {
        self.outfall_exports = OutfallExports::default();
        self.runtime_coupling_used = false;
        self.surface_outfall_link_limit.clear();
        self.surface_outfall_quality = QualityOnlyInflow::default();
        self.link_gain_quality = QualityOnlyInflow::default();
        self.state = EngineState::Created;
        self.control_log.clear();
        self.control_rule_names.clear();
        self.control_log_dropped = 0;
        self.current_time = 0.0;
        self.current_date = 0.0;
        self.dt_controls_remaining = 0.0;
        self.elapsed_ms = 0.0;
        self.old_elapsed_ms = 0.0;
        self.next_report_ms = 0.0;
        self.error_code = 0;
        self.warning_code = 0;
        self.error_message.clear();
        self.warnings.clear();
        self.errors.clear();
        self.title_notes.clear();
        self.deferred_section_rows.clear();
        self.unknown_section_headers.clear();
        self.pending_gw_nodes.clear();
        self.pending_gw_aquifers.clear();
        self.pending_link_nodes.clear();
        self.nodes = NodeData::default();
        self.links = LinkData::default();
        self.node_subtypes.clear();
        self.link_subtypes.clear();
        self.node_coupling_2d.clear();
        self.control_rules = ControlRuleStore::default();
        self.subcatches = SubcatchData::default();
        self.gages = GageData::default();
        self.pollutants = PollutantData::default();
        self.tables = TableData::default();
        self.snowpacks = SnowpackStore::default();
        self.node_names.clear();
        self.link_names.clear();
        self.subcatch_names.clear();
        self.gage_names.clear();
        self.pollutant_names.clear();
        self.snowpack_names.clear();
        self.rdii_decay = RDIIDecayData::default();
        self.ard_config = ArdConfigData::default();
        self.water_age_config = WaterAgeConfigData::default();
        self.water_age_state.clear();
        self.heat_config = HeatConfigData::default();
        self.heat_state.clear();
        self.lid_layer_state.clear();
        self.bed_state.clear();
        self.vj_diag.clear();
        self.inlet_diag.clear();
        self.climate_state = ClimateState::default();
        self.spatial = SpatialFrame::default();
        self.user_flags.clear();
        self.events.clear();
        self.adjust_temp.fill(0.0);
        self.adjust_evap.fill(0.0);
        self.adjust_rain.fill(1.0);
        self.adjust_hydcon.fill(1.0);
        self.subcatch_n_perv_pattern.clear();
        self.subcatch_d_store_pattern.clear();
        self.subcatch_infil_pattern.clear();
        self.base_n_perv.clear();
        self.base_ds_perv.clear();
        self.has_subcatch_adj_patterns = false;
        self.forcing = ForcingData::default();
    }

    pub fn save_state(&mut self) {
        self.nodes.save_state();
        self.links.save_state();
    }

    pub fn save_lat_qual_state(&mut self) {
        self.nodes.save_lat_qual_state();
        self.links.save_qual_state();
    }

    pub fn save_hyd_state(&mut self) {
        self.nodes.save_hyd_state();
        self.links.save_hyd_state();
    }

    pub fn reset_state(&mut self) {
        self.runtime_coupling_used = false;
        self.nodes.reset_state();
        self.links.reset_state();
        self.subcatches.reset_state();
        self.gages.reset_state();
        self.tables.reset_cursors();
        self.current_time = 0.0;
        self.current_date = self.options.start_date;
        self.dt_controls_remaining = 0.0;
        self.elapsed_ms = 0.0;
        self.old_elapsed_ms = 0.0;
        self.next_report_ms = 1000.0 * self.options.report_step;
    }

    pub fn allocate_objects(&mut self) {
        self.nodes.resize(self.node_names.size());
        self.links.resize(self.link_names.size());
        self.subcatches.resize(self.subcatch_names.size());
        self.gages.resize(self.gage_names.size());
        self.pollutants.resize_pollutants(self.pollutant_names.size());
        let np = self.pollutant_names.size();
        self.nodes.resize_quality(np);
        self.links.resize_quality(np);
        self.subcatches.resize_quality(np);
        let nn = self.node_names.size() as usize;
        let nl = self.link_names.size() as usize;
        let ns = self.subcatch_names.size() as usize;
        let ng = self.gage_names.size() as usize;
        let sp = &mut self.spatial;
        sp.node_x.clear();
        sp.node_x.resize(nn, 0.0);
        sp.node_y.clear();
        sp.node_y.resize(nn, 0.0);
        sp.node_has_xy.clear();
        sp.node_has_xy.resize(nn, 0);
        sp.link_x.clear();
        sp.link_x.resize(nl, 0.0);
        sp.link_y.clear();
        sp.link_y.resize(nl, 0.0);
        sp.subcatch_x.clear();
        sp.subcatch_x.resize(ns, 0.0);
        sp.subcatch_y.clear();
        sp.subcatch_y.resize(ns, 0.0);
        sp.link_vertices_x.resize(nl, Vec::new());
        sp.link_vertices_y.resize(nl, Vec::new());
        sp.subcatch_polygon_x.resize(ns, Vec::new());
        sp.subcatch_polygon_y.resize(ns, Vec::new());
        sp.gage_x.clear();
        sp.gage_x.resize(ng, 0.0);
        sp.gage_y.clear();
        sp.gage_y.resize(ng, 0.0);
    }

    pub fn shrink_all_to_fit(&mut self) {
        self.nodes.shrink_to_fit();
        self.links.shrink_to_fit();
        self.subcatches.shrink_to_fit();
        self.gages.shrink_to_fit();
        self.pollutants.shrink_to_fit();
        self.landuses.shrink_to_fit();
        self.buildup.shrink_to_fit();
        self.washoff.shrink_to_fit();
        self.treatment.shrink_to_fit();
        self.spatial.shrink_to_fit();
    }

    pub fn n_nodes(&self) -> i32 {
        self.node_names.size()
    }
    pub fn n_links(&self) -> i32 {
        self.link_names.size()
    }
    pub fn n_subcatches(&self) -> i32 {
        self.subcatch_names.size()
    }
    pub fn n_gages(&self) -> i32 {
        self.gage_names.size()
    }
    pub fn n_pollutants(&self) -> i32 {
        self.pollutant_names.size()
    }
    pub fn n_landuses(&self) -> i32 {
        self.landuse_names.size()
    }
    pub fn n_tables(&self) -> i32 {
        self.tables.count() as i32
    }
    pub fn find_timeseries(&self, name: &[u8]) -> i32 {
        self.tables.find_by_kind(name, true)
    }
    pub fn find_curve(&self, name: &[u8]) -> i32 {
        self.tables.find_by_kind(name, false)
    }
    pub fn find_table_any(&self, name: &[u8]) -> i32 {
        let ts = self.find_timeseries(name);
        if ts >= 0 {
            return ts;
        }
        self.find_curve(name)
    }
}
