// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/LID.hpp and LID.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// LID controls, batch by type group: the LID type enum, the per-group store
// (LIDGroupSoA, resize), the solver's members and LIDSolver::init, which
// SWMMEngine::initHydrology calls first (it fills subcatches.total_lid_area_ft2, which
// RunoffSolver::init subtracts). Also storedVolume and completedUnitStorage, which read
// only what init builds.
//
// Differences of form:
//   * LidRptFile's std::ofstream is an Option<std::fs::File>; the files are opened
//     lazily by the step functions (none at init).
//   * The solver's members are public so the gate dumper can read them; upstream's are
//     private.
//   * NOT YET TRANSLATED (run-time step functions, phase 1c): waterTrial,
//     executeCompleted, checkpoint, validateCheckpoint, totalStoredVolume,
//     totalInitVolume, totalInfilVolume, totalEvapVolume, execute, anyWet,
//     saveLidReport, the eight batch*Flux kernels and LID.cpp's file-static process
//     code (LidProcView, runUnitLegacy, ...).

use std::collections::BTreeMap;

use crate::core::cmath::cmax;
use crate::core::context::SimulationContext;
use crate::core::unit_conversion as ucf;
use crate::hydrology::infiltration::{self as infil, GreenAmptState};

/// LID type codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum LIDType {
    BioCell = 0,
    RainGarden = 1,
    GreenRoof = 2,
    InfilTrench = 3,
    PermPavement = 4,
    RainBarrel = 5,
    VegSwale = 6,
    RoofDiscon = 7,
}

pub fn has_native_bottom(type_: LIDType, storage_thickness: f64, storage_ksat: f64) -> bool {
    if type_ == LIDType::VegSwale {
        return true;
    }
    if type_ != LIDType::BioCell && type_ != LIDType::RainGarden && type_ != LIDType::InfilTrench && type_ != LIDType::PermPavement {
        return false;
    }
    storage_thickness <= 0.0 || storage_ksat > 0.0
}

/// The string overload of hasNativeBottom.
pub fn has_native_bottom_code(type_: &[u8], thickness: f64, ksat: f64) -> bool {
    if type_ == b"VS" {
        return true;
    }
    (type_ == b"BC" || type_ == b"RG" || type_ == b"IT" || type_ == b"PP") && (thickness <= 0.0 || ksat > 0.0)
}

pub const SURF: i32 = 0;
pub const SOIL: i32 = 1;
pub const STOR: i32 = 2;
pub const PAVE: i32 = 3;
pub const N_LAYERS: i32 = 4;

/// LID unit parameters, one group per LID type.
#[derive(Debug, Clone)]
pub struct LIDGroupSoA {
    pub type_: LIDType,
    pub count: i32,

    /// ft -> in|mm (RAINDEPTH)
    pub ucf_raindepth: f64,
    /// ft/sec -> in/hr|mm/hr (RAINFALL)
    pub ucf_rainfall: f64,

    pub subcatch_idx: Vec<i32>,
    pub control_idx: Vec<i32>,
    /// Unit area (ft2)
    pub area: Vec<f64>,
    pub from_imperv: Vec<f64>,
    pub from_perv: Vec<f64>,
    pub to_perv: Vec<i32>,
    pub drain_node: Vec<i32>,
    pub drain_subcatch: Vec<i32>,
    pub inflow: Vec<f64>,
    pub evap_rate_unit: Vec<f64>,
    pub subcatch_rain: Vec<f64>,
    pub rpt_path: Vec<Vec<u8>>,
    pub rpt_surf_infil: Vec<f64>,
    pub rpt_pave_perc: Vec<f64>,
    pub rpt_soil_perc: Vec<f64>,

    // Surface layer
    pub surf_store: Vec<f64>,
    pub surf_rough: Vec<f64>,
    pub surf_slope: Vec<f64>,

    // Soil layer
    pub soil_thick: Vec<f64>,
    pub soil_poros: Vec<f64>,
    pub soil_fc: Vec<f64>,
    pub soil_wp: Vec<f64>,
    pub soil_ksat: Vec<f64>,
    pub soil_kslope: Vec<f64>,
    pub soil_suction: Vec<f64>,

    // Storage layer
    pub stor_thick: Vec<f64>,
    pub stor_void: Vec<f64>,
    pub stor_ksat: Vec<f64>,
    pub stor_clog: Vec<f64>,
    pub stor_covered: Vec<i32>,

    // Drain
    pub drain_coeff: Vec<f64>,
    pub drain_expon: Vec<f64>,
    pub drain_offset: Vec<f64>,
    pub drain_delay: Vec<f64>,
    pub drain_hopen: Vec<f64>,
    pub drain_hclose: Vec<f64>,
    pub drain_open: Vec<i32>,

    // Pavement layer (PERM_PAVEMENT)
    pub pave_thick: Vec<f64>,
    pub pave_void: Vec<f64>,
    pub pave_imperv_frac: Vec<f64>,
    pub pave_ksat: Vec<f64>,
    pub pave_clog_factor: Vec<f64>,
    pub pave_regen_days: Vec<f64>,
    pub pave_regen_deg: Vec<f64>,
    pub next_regen_day: Vec<f64>,

    // Drainage mat layer (GREEN_ROOF)
    pub drainmat_thick: Vec<f64>,
    pub drainmat_void: Vec<f64>,
    pub drainmat_rough: Vec<f64>,

    // Surface geometry
    pub surf_void_frac: Vec<f64>,
    pub surf_alpha: Vec<f64>,
    pub surf_side_slope: Vec<f64>,
    pub full_width: Vec<f64>,
    pub unit_area: Vec<f64>,
    pub unit_width: Vec<f64>,
    pub dry_time: Vec<f64>,

    // State variables
    pub surf_depth: Vec<f64>,
    pub soil_moist: Vec<f64>,
    pub stor_depth: Vec<f64>,
    pub pave_depth: Vec<f64>,

    pub soil_infil: Vec<GreenAmptState>,
    pub old_drain_flow: Vec<f64>,
    pub is_wet: Vec<u8>,
    pub can_overflow: Vec<u8>,
    pub drainmat_alpha: Vec<f64>,

    // Outputs (per unit)
    pub surface_runoff: Vec<f64>,
    pub drain_flow: Vec<f64>,
    pub evap_loss: Vec<f64>,
    pub infil_loss: Vec<f64>,

    /// drain_rmvl[unit * n_pollutants + pollutant]
    pub drain_rmvl: Vec<f64>,
    pub n_pollutants: i32,

    // A4: per-layer inflow rates
    pub in_surf: Vec<f64>,
    pub in_pave: Vec<f64>,
    pub in_soil: Vec<f64>,
    pub in_stor: Vec<f64>,
    pub evap_surf: Vec<f64>,
    pub evap_pave: Vec<f64>,
    pub evap_soil: Vec<f64>,
    pub evap_stor: Vec<f64>,

    // Previous flux rates (Modified Puls time weighting)
    pub f_old_surf: Vec<f64>,
    pub f_old_soil: Vec<f64>,
    pub f_old_stor: Vec<f64>,
    pub f_old_pave: Vec<f64>,

    // Water balance tracking (cumulative per unit)
    pub wb_inflow: Vec<f64>,
    pub wb_evap: Vec<f64>,
    pub spatial_infil_loss: Vec<f64>,
    pub wb_spatial_infil: Vec<f64>,
    pub wb_infil: Vec<f64>,
    pub wb_surf_flow: Vec<f64>,
    pub wb_drain_flow: Vec<f64>,
    pub wb_init_vol: Vec<f64>,
    pub wb_final_vol: Vec<f64>,
    pub vol_treated: Vec<f64>,
}

impl Default for LIDGroupSoA {
    fn default() -> Self {
        LIDGroupSoA {
            type_: LIDType::BioCell,
            count: 0,
            ucf_raindepth: 12.0,
            ucf_rainfall: 43200.0,
            subcatch_idx: Vec::new(),
            control_idx: Vec::new(),
            area: Vec::new(),
            from_imperv: Vec::new(),
            from_perv: Vec::new(),
            to_perv: Vec::new(),
            drain_node: Vec::new(),
            drain_subcatch: Vec::new(),
            inflow: Vec::new(),
            evap_rate_unit: Vec::new(),
            subcatch_rain: Vec::new(),
            rpt_path: Vec::new(),
            rpt_surf_infil: Vec::new(),
            rpt_pave_perc: Vec::new(),
            rpt_soil_perc: Vec::new(),
            surf_store: Vec::new(),
            surf_rough: Vec::new(),
            surf_slope: Vec::new(),
            soil_thick: Vec::new(),
            soil_poros: Vec::new(),
            soil_fc: Vec::new(),
            soil_wp: Vec::new(),
            soil_ksat: Vec::new(),
            soil_kslope: Vec::new(),
            soil_suction: Vec::new(),
            stor_thick: Vec::new(),
            stor_void: Vec::new(),
            stor_ksat: Vec::new(),
            stor_clog: Vec::new(),
            stor_covered: Vec::new(),
            drain_coeff: Vec::new(),
            drain_expon: Vec::new(),
            drain_offset: Vec::new(),
            drain_delay: Vec::new(),
            drain_hopen: Vec::new(),
            drain_hclose: Vec::new(),
            drain_open: Vec::new(),
            pave_thick: Vec::new(),
            pave_void: Vec::new(),
            pave_imperv_frac: Vec::new(),
            pave_ksat: Vec::new(),
            pave_clog_factor: Vec::new(),
            pave_regen_days: Vec::new(),
            pave_regen_deg: Vec::new(),
            next_regen_day: Vec::new(),
            drainmat_thick: Vec::new(),
            drainmat_void: Vec::new(),
            drainmat_rough: Vec::new(),
            surf_void_frac: Vec::new(),
            surf_alpha: Vec::new(),
            surf_side_slope: Vec::new(),
            full_width: Vec::new(),
            unit_area: Vec::new(),
            unit_width: Vec::new(),
            dry_time: Vec::new(),
            surf_depth: Vec::new(),
            soil_moist: Vec::new(),
            stor_depth: Vec::new(),
            pave_depth: Vec::new(),
            soil_infil: Vec::new(),
            old_drain_flow: Vec::new(),
            is_wet: Vec::new(),
            can_overflow: Vec::new(),
            drainmat_alpha: Vec::new(),
            surface_runoff: Vec::new(),
            drain_flow: Vec::new(),
            evap_loss: Vec::new(),
            infil_loss: Vec::new(),
            drain_rmvl: Vec::new(),
            n_pollutants: 0,
            in_surf: Vec::new(),
            in_pave: Vec::new(),
            in_soil: Vec::new(),
            in_stor: Vec::new(),
            evap_surf: Vec::new(),
            evap_pave: Vec::new(),
            evap_soil: Vec::new(),
            evap_stor: Vec::new(),
            f_old_surf: Vec::new(),
            f_old_soil: Vec::new(),
            f_old_stor: Vec::new(),
            f_old_pave: Vec::new(),
            wb_inflow: Vec::new(),
            wb_evap: Vec::new(),
            spatial_infil_loss: Vec::new(),
            wb_spatial_infil: Vec::new(),
            wb_infil: Vec::new(),
            wb_surf_flow: Vec::new(),
            wb_drain_flow: Vec::new(),
            wb_init_vol: Vec::new(),
            wb_final_vol: Vec::new(),
            vol_treated: Vec::new(),
        }
    }
}

fn assign<T: Clone>(v: &mut Vec<T>, n: usize, d: T) {
    v.clear();
    v.resize(n, d);
}

impl LIDGroupSoA {
    /// Sizes and clears exactly the columns upstream's resize does (drain_rmvl is cleared;
    /// evap_surf/pave/soil/stor are assigned; ucf_*, type and n_pollutants are kept).
    pub fn resize(&mut self, n: i32) {
        self.count = n;
        let un = n as usize;

        assign(&mut self.subcatch_idx, un, -1);
        assign(&mut self.control_idx, un, -1);
        assign(&mut self.area, un, 0.0);
        assign(&mut self.from_imperv, un, 0.0);
        assign(&mut self.from_perv, un, 0.0);
        assign(&mut self.to_perv, un, 0);
        assign(&mut self.drain_node, un, -1);
        assign(&mut self.drain_subcatch, un, -1);
        assign(&mut self.inflow, un, 0.0);
        assign(&mut self.evap_rate_unit, un, 0.0);

        assign(&mut self.surf_store, un, 0.0);
        assign(&mut self.surf_rough, un, 0.01);
        assign(&mut self.surf_slope, un, 0.01);

        assign(&mut self.soil_thick, un, 0.0);
        assign(&mut self.soil_poros, un, 0.4);
        assign(&mut self.soil_fc, un, 0.2);
        assign(&mut self.soil_wp, un, 0.1);
        assign(&mut self.soil_ksat, un, 0.0);
        assign(&mut self.soil_kslope, un, 0.0);
        assign(&mut self.soil_suction, un, 0.0);

        assign(&mut self.stor_thick, un, 0.0);
        assign(&mut self.stor_void, un, 0.5);
        assign(&mut self.stor_ksat, un, 0.0);
        assign(&mut self.stor_clog, un, 0.0);
        assign(&mut self.stor_covered, un, 0);

        assign(&mut self.drain_coeff, un, 0.0);
        assign(&mut self.drain_expon, un, 0.5);
        assign(&mut self.drain_offset, un, 0.0);
        assign(&mut self.drain_delay, un, 0.0);
        assign(&mut self.drain_hopen, un, 0.0);
        assign(&mut self.drain_hclose, un, 0.0);
        assign(&mut self.drain_open, un, 1);

        assign(&mut self.pave_thick, un, 0.0);
        assign(&mut self.pave_void, un, 0.15);
        assign(&mut self.pave_imperv_frac, un, 0.0);
        assign(&mut self.pave_ksat, un, 0.0);
        assign(&mut self.pave_clog_factor, un, 0.0);
        assign(&mut self.pave_regen_days, un, 0.0);
        assign(&mut self.pave_regen_deg, un, 0.0);
        assign(&mut self.next_regen_day, un, 0.0);

        assign(&mut self.drainmat_thick, un, 0.0);
        assign(&mut self.drainmat_void, un, 0.5);
        assign(&mut self.drainmat_rough, un, 0.1);

        assign(&mut self.surf_void_frac, un, 1.0);
        assign(&mut self.surf_alpha, un, 0.0);
        assign(&mut self.surf_side_slope, un, 0.0);
        assign(&mut self.full_width, un, 0.0);
        assign(&mut self.unit_area, un, 0.0);
        assign(&mut self.unit_width, un, 0.0);
        assign(&mut self.dry_time, un, 0.0);
        assign(&mut self.subcatch_rain, un, 0.0);

        assign(&mut self.surf_depth, un, 0.0);
        assign(&mut self.soil_moist, un, 0.2);
        assign(&mut self.stor_depth, un, 0.0);
        assign(&mut self.pave_depth, un, 0.0);

        // drain_rmvl sized separately in model builder (needs n_pollutants)
        self.drain_rmvl.clear();

        assign(&mut self.in_surf, un, 0.0);
        assign(&mut self.in_pave, un, 0.0);
        assign(&mut self.in_soil, un, 0.0);
        assign(&mut self.in_stor, un, 0.0);
        assign(&mut self.evap_surf, un, 0.0);
        assign(&mut self.evap_pave, un, 0.0);
        assign(&mut self.evap_soil, un, 0.0);
        assign(&mut self.evap_stor, un, 0.0);

        assign(&mut self.f_old_surf, un, 0.0);
        assign(&mut self.f_old_soil, un, 0.0);
        assign(&mut self.f_old_stor, un, 0.0);
        assign(&mut self.f_old_pave, un, 0.0);

        assign(&mut self.soil_infil, un, GreenAmptState::default());
        assign(&mut self.old_drain_flow, un, 0.0);
        assign(&mut self.is_wet, un, 0);
        assign(&mut self.rpt_path, un, Vec::new());
        assign(&mut self.rpt_surf_infil, un, 0.0);
        assign(&mut self.rpt_pave_perc, un, 0.0);
        assign(&mut self.rpt_soil_perc, un, 0.0);
        assign(&mut self.can_overflow, un, 1);
        assign(&mut self.drainmat_alpha, un, 0.0);
        assign(&mut self.surface_runoff, un, 0.0);
        assign(&mut self.drain_flow, un, 0.0);
        assign(&mut self.evap_loss, un, 0.0);
        assign(&mut self.infil_loss, un, 0.0);

        assign(&mut self.wb_inflow, un, 0.0);
        assign(&mut self.wb_evap, un, 0.0);
        assign(&mut self.wb_infil, un, 0.0);
        assign(&mut self.spatial_infil_loss, un, 0.0);
        assign(&mut self.wb_spatial_infil, un, 0.0);
        assign(&mut self.wb_surf_flow, un, 0.0);
        assign(&mut self.wb_drain_flow, un, 0.0);
        assign(&mut self.wb_init_vol, un, 0.0);
        assign(&mut self.wb_final_vol, un, 0.0);
        assign(&mut self.vol_treated, un, 0.0);
    }
}

/// One unit's stored water (ft3) as completed-source accounting measures it.
pub fn completed_unit_storage(g: &LIDGroupSoA, u: i32) -> f64 {
    let u = u as usize;
    if g.type_ == LIDType::VegSwale {
        let top = cmax(g.unit_width[u], 0.5);
        let height = g.surf_store[u];
        let mut side = g.surf_side_slope[u];
        let mut bottom = top - 2.0 * side * height;
        if bottom < 0.5 {
            bottom = 0.5;
            side = 0.5 * (top - 0.5) / height;
        }
        let depth = g.surf_depth[u];
        return g.area[u] / top * depth * (bottom + side * depth) * g.surf_void_frac[u];
    }
    (g.surf_depth[u] * g.surf_void_frac[u] + g.soil_moist[u] * g.soil_thick[u] + g.stor_depth[u] * g.stor_void[u] + g.pave_depth[u] * g.pave_void[u] * (1.0 - g.pave_imperv_frac[u])) * g.area[u]
}

/// Decode an exact, case-sensitive authored LID type code: BC/RG/GR/IT/PP/RB/VS/RD give
/// 0..7, anything else -1.
fn type_code_to_index(code: &[u8]) -> i32 {
    match code {
        b"BC" => 0,
        b"RG" => 1,
        b"GR" => 2,
        b"IT" => 3,
        b"PP" => 4,
        b"RB" => 5,
        b"VS" => 6,
        b"RD" => 7,
        _ => -1,
    }
}

/// One legacy TLidRptFile per reporting unit, opened lazily on the unit's first step.
#[derive(Debug, Default)]
pub struct LidRptFile {
    pub f: Option<std::fs::File>,
    pub was_dry: i32,
    pub held: Vec<u8>,
    pub failed: bool,
}

/// LIDSolver::CompletedUnitInput.
#[derive(Debug, Clone, Copy)]
pub struct CompletedUnitInput {
    pub type_: i32,
    pub unit: i32,
    pub inflow: f64,
    pub rain: f64,
    pub pet: f64,
    pub native_infil: f64,
    pub max_native_infil: f64,
    pub infil_factor: f64,
    pub spatial_fraction: f64,
    pub spatial_max_infil: f64,
}

impl Default for CompletedUnitInput {
    fn default() -> Self {
        CompletedUnitInput { type_: -1, unit: -1, inflow: 0.0, rain: 0.0, pet: 0.0, native_infil: 0.0, max_native_infil: 1.0e10, infil_factor: 1.0, spatial_fraction: -1.0, spatial_max_infil: 1.0e10 }
    }
}

#[derive(Debug)]
pub struct LIDSolver {
    pub groups_: Vec<LIDGroupSoA>,
    pub usage_order_: Vec<(i32, i32)>,
    pub old_runoff_sec_: f64,
    pub rpt_files_: BTreeMap<(i32, i32), Box<LidRptFile>>,
    pub recovery_factor_: f64,
    pub native_infil_: Vec<f64>,
    pub max_native_infil_: Vec<f64>,
    pub infil_factor_: Vec<f64>,
}

impl Default for LIDSolver {
    fn default() -> Self {
        LIDSolver { groups_: Vec::new(), usage_order_: Vec::new(), old_runoff_sec_: 0.0, rpt_files_: BTreeMap::new(), recovery_factor_: 1.0, native_infil_: Vec::new(), max_native_infil_: Vec::new(), infil_factor_: Vec::new() }
    }
}

impl LIDSolver {
    pub fn group(&self, type_index: i32) -> &LIDGroupSoA {
        &self.groups_[type_index as usize]
    }
    pub fn group_mut(&mut self, type_index: i32) -> &mut LIDGroupSoA {
        &mut self.groups_[type_index as usize]
    }
    pub fn num_groups(&self) -> i32 {
        self.groups_.len() as i32
    }
    pub fn usage_order(&self) -> &Vec<(i32, i32)> {
        &self.usage_order_
    }
    pub fn set_runoff_time(&mut self, old_runoff_sec: f64) {
        self.old_runoff_sec_ = old_runoff_sec;
    }
    pub fn set_native_infil(&mut self, native: Vec<f64>, max_native: Vec<f64>) {
        self.native_infil_ = native;
        self.max_native_infil_ = max_native;
    }
    pub fn set_infil_factors(&mut self, infil_factor: Vec<f64>, recovery_factor: f64) {
        self.infil_factor_ = infil_factor;
        self.recovery_factor_ = recovery_factor;
    }

    pub fn init(&mut self, ctx: &mut SimulationContext) {
        const N_LID_TYPES: i32 = 8;
        self.groups_.resize(N_LID_TYPES as usize, LIDGroupSoA::default());

        // Assign type codes
        self.groups_[0].type_ = LIDType::BioCell;
        self.groups_[1].type_ = LIDType::RainGarden;
        self.groups_[2].type_ = LIDType::GreenRoof;
        self.groups_[3].type_ = LIDType::InfilTrench;
        self.groups_[4].type_ = LIDType::PermPavement;
        self.groups_[5].type_ = LIDType::RainBarrel;
        self.groups_[6].type_ = LIDType::VegSwale;
        self.groups_[7].type_ = LIDType::RoofDiscon;

        // Unit-conversion factors for the model's flow units (issue #102).
        let ucf_rain_depth = ucf::ucf(ucf::RAINDEPTH, &ctx.options);
        let ucf_rainfall = ucf::ucf(ucf::RAINFALL, &ctx.options);
        let ucf_length = ucf::ucf(ucf::LENGTH, &ctx.options);
        let ucf_length2 = ucf_length * ucf_length;

        // If no LID usage data, leave all groups empty
        let n_usage = ctx.lid_usage.count();
        ctx.lid_usage.resize_wb(n_usage);
        if n_usage == 0 {
            for g in self.groups_.iter_mut() {
                g.resize(0);
            }
            return;
        }

        // 1. Count units per type from usage entries
        let mut type_counts = [0i32; 8];
        for j in 0..n_usage {
            let uj = j as usize;
            let li = ctx.lid_usage.lid_index[uj];
            if li < 0 || li >= ctx.lid_controls.count() {
                continue;
            }
            let ti = type_code_to_index(&ctx.lid_controls.lid_type[li as usize]);
            if ti < 0 || ti >= N_LID_TYPES {
                continue;
            }
            type_counts[ti as usize] += 1;
        }

        // 2. Resize groups
        for t in 0..N_LID_TYPES {
            self.groups_[t as usize].resize(type_counts[t as usize]);
        }

        // The group's unit-conversion pair feeds the drain law itself.
        for g in self.groups_.iter_mut() {
            g.ucf_raindepth = ucf_rain_depth;
            g.ucf_rainfall = ucf_rainfall;
        }

        // 3. Populate per-unit parameters from LidControlStore + LidUsageStore
        let mut type_cursor = [0i32; 8];
        self.usage_order_.clear();
        self.usage_order_.reserve(n_usage as usize);
        for j in 0..n_usage {
            let uj = j as usize;
            let li = ctx.lid_usage.lid_index[uj];
            if li < 0 || li >= ctx.lid_controls.count() {
                continue;
            }
            let uli = li as usize;
            let ti = type_code_to_index(&ctx.lid_controls.lid_type[uli]);
            if ti < 0 || ti >= N_LID_TYPES {
                continue;
            }

            let g = &mut self.groups_[ti as usize];
            let slot = type_cursor[ti as usize];
            type_cursor[ti as usize] += 1;
            let us = slot as usize;
            self.usage_order_.push((ti, slot)); // parse order; reversed per sub below

            let usage = &ctx.lid_usage;
            let controls = &ctx.lid_controls;
            // Usage-level fields: area and full_width scaled by the replicate count,
            // converted from display units one unit first (lid.c:551-552), then times
            // the count (lid.c:1168/1688).
            let n_units = usage.number[uj] as f64;
            g.subcatch_idx[us] = usage.subcatch_index[uj];
            g.control_idx[us] = li;
            g.area[us] = (usage.area[uj] / ucf_length2) * n_units;
            g.full_width[us] = (usage.width[uj] / ucf_length) * n_units;
            g.unit_area[us] = usage.area[uj] / ucf_length2;
            g.unit_width[us] = usage.width[uj] / ucf_length;
            g.from_imperv[us] = usage.from_imperv[uj] / 100.0; // % -> fraction
            g.from_perv[us] = if uj < usage.from_perv.len() { usage.from_perv[uj] / 100.0 } else { 0.0 };
            g.to_perv[us] = usage.to_perv[uj];
            if uj < usage.rpt_file.len() {
                let fp = &usage.rpt_file[uj];
                g.rpt_path[us] = if fp.absolute.is_empty() { fp.original.clone() } else { fp.absolute.clone() };
            }
            // legacy validateLidGroup: no pervious area to return to on a (near-)fully
            // impervious subcatchment.
            {
                let sc0 = usage.subcatch_index[uj];
                if sc0 >= 0 && (sc0 as usize) < ctx.subcatches.frac_imperv.len() && ctx.subcatches.frac_imperv[sc0 as usize] >= 0.999 {
                    g.to_perv[us] = 0;
                }
            }

            // Resolve drain-to target
            if uj < usage.drain_to.len() && !usage.drain_to[uj].is_empty() {
                let dt_name = &usage.drain_to[uj];
                let ni = ctx.node_names.find(dt_name);
                if ni >= 0 {
                    g.drain_node[us] = ni;
                } else {
                    let si = ctx.subcatch_names.find(dt_name);
                    if si >= 0 {
                        g.drain_subcatch[us] = si;
                    }
                }
            }

            // SURFACE layer: [0]=StorHt, [1]=VegVolFrac, [2]=Roughness, [3]=SurfSlope, [4]=SideSlope
            if uli < controls.surface.len() {
                let p = &controls.surface[uli];
                let veg = if p[0] == 0.0 { 0.0 } else { p[1] };
                g.surf_store[us] = p[0] / ucf_rain_depth;
                g.surf_void_frac[us] = 1.0 - veg;
                g.surf_rough[us] = p[2];
                g.surf_slope[us] = p[3] / 100.0;
                g.surf_side_slope[us] = p[4];
            }

            // SOIL layer: [0]=Thick, [1]=Poros, [2]=FC, [3]=WP, [4]=Ksat, [5]=Kslope, [6]=Suction
            if uli < controls.soil.len() {
                let p = &controls.soil[uli];
                g.soil_thick[us] = p[0] / ucf_rain_depth;
                g.soil_poros[us] = p[1];
                g.soil_fc[us] = p[2];
                g.soil_wp[us] = p[3];
                g.soil_ksat[us] = p[4] / ucf_rainfall;
                g.soil_kslope[us] = p[5];
                g.soil_suction[us] = p[6] / ucf_rain_depth;
            }

            // STORAGE layer: [0]=Thick, [1]=VoidRatio, [2]=Ksat, [3]=ClogFactor
            if uli < controls.storage.len() {
                let p = &controls.storage[uli];
                g.stor_thick[us] = p[0] / ucf_rain_depth;
                g.stor_void[us] = if p[1] > 0.0 { p[1] / (p[1] + 1.0) } else { 0.0 };
                g.stor_ksat[us] = p[2] / ucf_rainfall;
                g.stor_clog[us] = p[3];
            }

            // DRAIN layer: [0]=Coeff, [1]=Expon, [2]=Offset, [3]=Delay, [4]=hOpen, [5]=hClose
            if uli < controls.drain.len() {
                let p = &controls.drain[uli];
                g.drain_coeff[us] = p[0];
                g.drain_expon[us] = p[1];
                g.drain_offset[us] = p[2] / ucf_rain_depth;
                g.drain_delay[us] = p[3] * 3600.0;
                g.drain_hopen[us] = p[4] / ucf_rain_depth;
                g.drain_hclose[us] = p[5] / ucf_rain_depth;
            }

            // PAVEMENT layer: [0]=Thick, [1]=VoidRatio, [2]=FracImperv, [3]=Ksat, [4]=ClogFactor, [5]=RegenDays
            if uli < controls.pavement.len() {
                let p = &controls.pavement[uli];
                g.pave_thick[us] = p[0] / ucf_rain_depth;
                g.pave_void[us] = if p[1] > 0.0 { p[1] / (p[1] + 1.0) } else { 0.0 };
                g.pave_imperv_frac[us] = p[2];
                g.pave_ksat[us] = p[3] / ucf_rainfall;
                g.pave_clog_factor[us] = p[4];
                if p[5] > 0.0 {
                    g.pave_regen_days[us] = p[5];
                    // upstream: `(uli < pavement.size() && p[5] > 0.0) ? 0.0 : 0.0`
                    g.pave_regen_deg[us] = 0.0;
                    g.next_regen_day[us] = p[5]; // first regen after p[5] days
                }
            }

            // DRAINMAT layer: [0]=Thick, [1]=VoidRatio, [2]=Roughness
            if uli < controls.drainmat.len() {
                let p = &controls.drainmat[uli];
                g.drainmat_thick[us] = p[0] / ucf_rain_depth;
                g.drainmat_void[us] = p[1];
                g.drainmat_rough[us] = p[2];
            }

            // Legacy validateLidProc (lid.c): the derived parameters, in its order.
            {
                const PHI: f64 = 1.486;
                if g.type_ == LIDType::VegSwale {
                    g.surf_alpha[us] = if g.surf_rough[us] * g.surf_slope[us] > 0.0 && g.surf_store[us] != 0.0 { PHI * g.surf_slope[us].sqrt() / g.surf_rough[us] } else { 0.0 };
                } else {
                    g.surf_alpha[us] = if g.surf_rough[us] > 0.0 { PHI / g.surf_rough[us] * g.surf_slope[us].sqrt() } else { 0.0 };
                }
                g.drainmat_alpha[us] = if g.drainmat_rough[us] > 0.0 { PHI / g.drainmat_rough[us] * g.surf_slope[us].sqrt() } else { 0.0 };
                // Clogging factors become the treated volume (ft) at which the layer is
                // fully clogged: void volume times the user's factor.
                if g.pave_thick[us] > 0.0 {
                    g.pave_clog_factor[us] *= g.pave_thick[us] * g.pave_void[us] * (1.0 - g.pave_imperv_frac[us]);
                }
                if g.stor_thick[us] > 0.0 {
                    g.stor_clog[us] *= g.stor_thick[us] * g.stor_void[us];
                } else {
                    g.stor_clog[us] = 0.0;
                    g.stor_void[us] = 1.0; // no storage layer
                    g.drain_offset[us] = 0.0;
                }
                let mut can_overflow = true;
                match g.type_ {
                    LIDType::RoofDiscon => can_overflow = false,
                    LIDType::InfilTrench | LIDType::PermPavement | LIDType::BioCell | LIDType::RainGarden | LIDType::GreenRoof => {
                        if g.surf_alpha[us] > 0.0 {
                            can_overflow = false;
                        }
                    }
                    _ => {}
                }
                g.can_overflow[us] = if can_overflow { 1 } else { 0 };
                if g.type_ == LIDType::RainBarrel {
                    g.stor_void[us] = 1.0;
                    g.stor_ksat[us] = 0.0;
                }
                if g.type_ == LIDType::GreenRoof {
                    // The drainage mat IS the storage layer.
                    g.stor_thick[us] = g.drainmat_thick[us];
                    g.stor_void[us] = g.drainmat_void[us];
                    g.stor_clog[us] = 0.0;
                    g.stor_ksat[us] = 0.0;
                }
            }

            // Initial saturation (percent -> fraction, legacy x[2]/100)
            let init_sat = usage.init_sat[uj] / 100.0;
            if g.soil_thick[us] > 0.0 {
                g.soil_moist[us] = g.soil_wp[us] + init_sat * (g.soil_poros[us] - g.soil_wp[us]);
            }
            if g.stor_thick[us] > 0.0 {
                g.stor_depth[us] = init_sat * g.stor_thick[us];
            }

            // The unit's own Green-Ampt state for surface-to-soil infiltration (legacy
            // validateLidGroup lid.c:1173-1185 + lid_initState), in the deck's units.
            g.soil_infil[us] = GreenAmptState::default();
            if g.soil_thick[us] > 0.0 {
                let (s, k, imd) = (g.soil_suction[us] * ucf_rain_depth, g.soil_ksat[us] * ucf_rainfall, (g.soil_poros[us] - g.soil_wp[us]) * (1.0 - init_sat));
                infil::grnampt_init(&mut g.soil_infil[us], s, k, imd, &ctx.options);
            }
            if g.type_ == LIDType::VegSwale {
                let sc = g.subcatch_idx[us];
                let usc = sc as usize;
                if sc >= 0 && usc < ctx.subcatches.infil_model.len() && (ctx.subcatches.infil_model[usc] == 2 || ctx.subcatches.infil_model[usc] == 3) {
                    infil::grnampt_init(&mut g.soil_infil[us], ctx.subcatches.infil_p1[usc], ctx.subcatches.infil_p2[usc], ctx.subcatches.infil_p3[usc], &ctx.options);
                }
            }
            // legacy lid_initState: dryTime starts at DRY_DAYS
            g.dry_time[us] = ctx.options.dry_days * 86400.0;
            g.old_drain_flow[us] = 0.0;
        }

        // 4. Initialize water balance initial volumes
        for g in self.groups_.iter_mut() {
            for i in 0..g.count {
                let ui = i as usize;
                let init_vol = g.surf_depth[ui] * g.surf_void_frac[ui] + g.soil_moist[ui] * g.soil_thick[ui] + g.stor_depth[ui] * g.stor_void[ui] + g.pave_depth[ui] * g.pave_void[ui] * (1.0 - g.pave_imperv_frac[ui]);
                g.wb_init_vol[ui] = init_vol;
                g.wb_final_vol[ui] = init_vol;
            }
        }

        // 5. Populate drain pollutant removal fractions
        let np = ctx.n_pollutants();
        if np > 0 {
            let mut cursor2 = [0i32; 8];
            for j in 0..n_usage {
                let uj = j as usize;
                let li = ctx.lid_usage.lid_index[uj];
                if li < 0 || li >= ctx.lid_controls.count() {
                    continue;
                }
                let uli = li as usize;
                let ti = type_code_to_index(&ctx.lid_controls.lid_type[uli]);
                if ti < 0 || ti >= N_LID_TYPES {
                    continue;
                }

                let g = &mut self.groups_[ti as usize];
                let slot = cursor2[ti as usize];
                cursor2[ti as usize] += 1;

                // Ensure drain_rmvl is sized
                if g.n_pollutants != np {
                    g.n_pollutants = np;
                    g.drain_rmvl.clear();
                    g.drain_rmvl.resize((g.count * np) as usize, 0.0);
                }

                // Copy removal fractions from LidControlStore
                if uli < ctx.lid_controls.removals.len() {
                    for pr in &ctx.lid_controls.removals[uli] {
                        let pi = pr.0;
                        let frac = pr.1;
                        if pi >= 0 && pi < np {
                            g.drain_rmvl[(slot * np + pi) as usize] = frac;
                        }
                    }
                }
            }
        }

        // Legacy's per-subcatchment lidList is built by prepending each parsed row
        // (lid.c:544-545), so units are walked in reverse parse order: subcatchments
        // ascending, each one's rows reversed.
        {
            let mut by_sub: BTreeMap<i32, Vec<(i32, i32)>> = BTreeMap::new();
            for tu in &self.usage_order_ {
                let gg = &self.groups_[tu.0 as usize];
                by_sub.entry(gg.subcatch_idx[tu.1 as usize]).or_default().push(*tu);
            }
            self.usage_order_.clear();
            for (_, v) in by_sub.iter() {
                for tu in v.iter().rev() {
                    self.usage_order_.push(*tu);
                }
            }
        }

        // Gap #60: accumulate total LID area (ft2) per subcatchment, per [LID_USAGE] row
        // in parse order (lid.c addLidUnit).
        let n_sc = ctx.n_subcatches();
        if n_sc > 0 {
            ctx.subcatches.total_lid_area_ft2.clear();
            ctx.subcatches.total_lid_area_ft2.resize(n_sc as usize, 0.0);
            for tu in &self.usage_order_ {
                let grp = &self.groups_[tu.0 as usize];
                let uu = tu.1 as usize;
                let sc = grp.subcatch_idx[uu];
                if sc >= 0 && sc < n_sc {
                    ctx.subcatches.total_lid_area_ft2[sc as usize] += grp.area[uu];
                }
            }
            // Legacy lid_validate() (lid.c:1234): snap the LID total to the full
            // subcatchment area when within 0.1% (issue #131).
            let ucf_area = ucf::ucf(ucf::LANDAREA, &ctx.options);
            for sc in 0..n_sc {
                let usc = sc as usize;
                let full_ft2 = ctx.subcatches.area[usc] / ucf_area;
                if ctx.subcatches.total_lid_area_ft2[usc] > 0.999 * full_ft2 {
                    ctx.subcatches.total_lid_area_ft2[usc] = full_ft2;
                }
            }
        }
    }

    /// Void-weighted layer depths times the unit footprint, summed (legacy
    /// lid_getStoredVolume). Units whose source has a nonnegative owner are skipped.
    pub fn stored_volume(&self, source_owners: Option<&Vec<i32>>) -> f64 {
        let mut total = 0.0;
        for g in &self.groups_ {
            for i in 0..g.count {
                let ui = i as usize;
                if let Some(owners) = source_owners {
                    if g.subcatch_idx[ui] >= 0 && owners[g.subcatch_idx[ui] as usize] >= 0 {
                        continue;
                    }
                }
                total += (g.surf_depth[ui] * g.surf_void_frac[ui] + g.soil_moist[ui] * g.soil_thick[ui] + g.stor_depth[ui] * g.stor_void[ui] + g.pave_depth[ui] * g.pave_void[ui] * (1.0 - g.pave_imperv_frac[ui])) * g.area[ui];
            }
        }
        total
    }
}
