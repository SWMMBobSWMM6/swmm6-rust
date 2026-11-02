// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/NodeData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson (first pass by tools/soa2rs.py, finished by hand). Not an
// OpenSWMM release.
//
// The node store, one Vec per attribute (structure of arrays), sized and cleared exactly
// as upstream sizes and clears each one: resize() sets every field it names, grow_to()
// only the ones it names, and quality arrays wait for resize_quality(). `type` is
// `type_` (a Rust keyword). std::copy into a longer destination copies the prefix, as
// here; into a shorter one it is undefined in C++ and a panic here.

#![allow(clippy::all)]
use crate::core::options::Bytes;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i8)]
pub enum NodeType {
    Junction = 0,
    Outfall = 1,
    Divider = 2,
    Storage = 3,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i8)]
pub enum OutfallType {
    Free = 0,
    Normal = 1,
    Fixed = 2,
    Tidal = 3,
    Timeseries = 4,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i8)]
pub enum DividerType {
    Cutoff = 0,
    OverflowDiv = 1,
    Tabular = 2,
    Weir = 3,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i8)]
pub enum StorageShape {
    Tabular = 0,
    Functional = 1,
    Cylindrical = 2,
    Conical = 3,
    Paraboloid = 4,
    Pyramidal = 5,
}

#[derive(Debug, Clone, Default)]
pub struct NodeData {
    pub type_: Vec<NodeType>,
    pub invert_elev: Vec<f64>,
    pub full_depth: Vec<f64>,
    pub init_depth: Vec<f64>,
    pub sur_depth: Vec<f64>,
    pub ponded_area: Vec<f64>,
    pub storage_area: Vec<f64>,
    pub is_virtual: Vec<u8>,
    pub is_inlet: Vec<u8>,
    pub rim_depth: Vec<f64>,
    pub depth: Vec<f64>,
    pub head: Vec<f64>,
    pub volume: Vec<f64>,
    pub lat_flow: Vec<f64>,
    pub user_lat_flow: Vec<f64>,
    pub runoff_inflow: Vec<f64>,
    pub gw_inflow: Vec<f64>,
    pub ext_inflow: Vec<f64>,
    pub dwf_inflow: Vec<f64>,
    pub rdii_inflow: Vec<f64>,
    pub iface_inflow: Vec<f64>,
    pub coupling_inflow: Vec<f64>,
    pub coupling_volume: Vec<f64>,
    pub aquifer2d_bed: Vec<u8>,
    pub coupling_queue: Vec<f64>,
    pub coupling_qual_queue: Vec<f64>,
    pub coupling_qual_inflow: Vec<f64>,
    pub coupling_age_vol_queue: Vec<f64>,
    pub coupling_temp_vol_queue: Vec<f64>,
    pub coupling_age_vol_inflow: Vec<f64>,
    pub coupling_temp_vol_inflow: Vec<f64>,
    pub coupling_tuple_age: bool,
    pub coupling_tuple_temp: bool,
    pub qual_mass_in: Vec<f64>,
    pub qual_vol_in: Vec<f64>,
    pub iface_qual_mass: Vec<f64>,
    pub ext_qual_mass: Vec<f64>,
    pub dwf_qual_mass: Vec<f64>,
    pub lid_drain_qual_load: Vec<f64>,
    pub lid_drain_qual_vol: Vec<f64>,
    pub lid_drain_qual_load_old: Vec<f64>,
    pub lid_drain_qual_vol_old: Vec<f64>,
    pub lid_drain_inflow: Vec<f64>,
    pub conc: Vec<f64>,
    pub conc_old: Vec<f64>,
    pub hrt: Vec<f64>,
    pub user_conc_mass_flux: Vec<f64>,
    pub conc_n_pollutants: i32,
    pub inflow: Vec<f64>,
    pub outflow: Vec<f64>,
    pub overflow: Vec<f64>,
    pub losses: Vec<f64>,
    pub crown_elev: Vec<f64>,
    pub degree: Vec<i32>,
    pub old_net_inflow: Vec<f64>,
    pub full_volume: Vec<f64>,
    pub rpt_full_volume: Vec<f64>,
    pub old_depth: Vec<f64>,
    pub old_volume: Vec<f64>,
    pub old_lat_flow: Vec<f64>,
    pub old_inflow: Vec<f64>,
    pub comments: Vec<Bytes>,
    pub tags: Vec<Bytes>,
    pub rpt_flag: Vec<i8>,
    pub stat_vol_flooded: Vec<f64>,
    pub stat_time_flooded: Vec<f64>,
    pub stat_max_depth: Vec<f64>,
    pub stat_max_overflow: Vec<f64>,
    pub stat_max_overflow_date: Vec<f64>,
    pub stat_sum_depth: Vec<f64>,
    pub stat_sum_volume: Vec<f64>,
    pub stat_storage_max_volume: Vec<f64>,
    pub stat_storage_max_volume_date: Vec<f64>,
    pub stat_max_depth_date: Vec<f64>,
    pub stat_max_rpt_depth: Vec<f64>,
    pub stat_max_inflow_date: Vec<f64>,
    pub stat_time_surcharged: Vec<f64>,
    pub stat_max_surcharge_height: Vec<f64>,
    pub stat_outfall_avg_flow: Vec<f64>,
    pub stat_max_lat_inflow: Vec<f64>,
    pub stat_max_total_inflow: Vec<f64>,
    pub stat_storage_max_outflow: Vec<f64>,
    pub stat_lat_inflow_vol: Vec<f64>,
    pub stat_total_inflow_vol: Vec<f64>,
    pub stat_total_outflow_vol: Vec<f64>,
    pub stat_outfall_max_flow: Vec<f64>,
    pub stat_outfall_periods: Vec<i64>,
    pub stat_non_converged_count: Vec<i32>,
    pub stat_time_courant_critical: Vec<f64>,
    pub stat_total_load: Vec<f64>,
    pub stat_n_pollutants: i32,
}

// The member initializers (coupling_tuple_age/temp = false, conc_n_pollutants =
// stat_n_pollutants = 0) are Rust's defaults, so Default is derived.

impl NodeData {
    pub fn count(&self) -> i32 {
        self.type_.len() as i32
    }

    pub fn resize(&mut self, n: i32) {
        let un = usize::try_from(n).unwrap_or(0);
        self.type_.clear();
        self.type_.resize(un, NodeType::Junction);
        self.invert_elev.clear();
        self.invert_elev.resize(un, 0.0);
        self.full_depth.clear();
        self.full_depth.resize(un, 0.0);
        self.init_depth.clear();
        self.init_depth.resize(un, 0.0);
        self.sur_depth.clear();
        self.sur_depth.resize(un, 0.0);
        self.ponded_area.clear();
        self.ponded_area.resize(un, 0.0);
        self.storage_area.clear();
        self.storage_area.resize(un, 0.0);
        self.is_virtual.clear();
        self.is_virtual.resize(un, 0);
        self.is_inlet.clear();
        self.is_inlet.resize(un, 0);
        self.rim_depth.clear();
        self.rim_depth.resize(un, 0.0);
        self.depth.clear();
        self.depth.resize(un, 0.0);
        self.head.clear();
        self.head.resize(un, 0.0);
        self.volume.clear();
        self.volume.resize(un, 0.0);
        self.lat_flow.clear();
        self.lat_flow.resize(un, 0.0);
        self.user_lat_flow.clear();
        self.user_lat_flow.resize(un, 0.0);
        self.runoff_inflow.clear();
        self.runoff_inflow.resize(un, 0.0);
        self.gw_inflow.clear();
        self.gw_inflow.resize(un, 0.0);
        self.ext_inflow.clear();
        self.ext_inflow.resize(un, 0.0);
        self.dwf_inflow.clear();
        self.dwf_inflow.resize(un, 0.0);
        self.rdii_inflow.clear();
        self.rdii_inflow.resize(un, 0.0);
        self.iface_inflow.clear();
        self.iface_inflow.resize(un, 0.0);
        self.coupling_inflow.clear();
        self.coupling_inflow.resize(un, 0.0);
        self.coupling_volume.clear();
        self.coupling_volume.resize(un, 0.0);
        self.coupling_queue.clear();
        self.coupling_queue.resize(un, 0.0);
        self.coupling_age_vol_queue.clear();
        self.coupling_age_vol_queue.resize(un, 0.0);
        self.coupling_temp_vol_queue.clear();
        self.coupling_temp_vol_queue.resize(un, 0.0);
        self.coupling_age_vol_inflow.clear();
        self.coupling_age_vol_inflow.resize(un, 0.0);
        self.coupling_temp_vol_inflow.clear();
        self.coupling_temp_vol_inflow.resize(un, 0.0);
        self.qual_mass_in.clear();
        self.iface_qual_mass.clear();
        self.coupling_qual_queue.clear();
        self.coupling_qual_inflow.clear();
        self.ext_qual_mass.clear();
        self.dwf_qual_mass.clear();
        self.qual_vol_in.clear();
        self.qual_vol_in.resize(un, 0.0);
        self.lid_drain_qual_load.clear();
        self.lid_drain_qual_vol.clear();
        self.lid_drain_qual_vol.resize(un, 0.0);
        self.lid_drain_qual_load_old.clear();
        self.lid_drain_qual_vol_old.clear();
        self.lid_drain_qual_vol_old.resize(un, 0.0);
        self.lid_drain_inflow.clear();
        self.lid_drain_inflow.resize(un, 0.0);
        self.inflow.clear();
        self.inflow.resize(un, 0.0);
        self.outflow.clear();
        self.outflow.resize(un, 0.0);
        self.overflow.clear();
        self.overflow.resize(un, 0.0);
        self.losses.clear();
        self.losses.resize(un, 0.0);
        self.crown_elev.clear();
        self.crown_elev.resize(un, 0.0);
        self.degree.clear();
        self.degree.resize(un, 0);
        self.old_net_inflow.clear();
        self.old_net_inflow.resize(un, 0.0);
        self.full_volume.clear();
        self.full_volume.resize(un, 0.0);
        self.rpt_full_volume.clear();
        self.rpt_full_volume.resize(un, 0.0);
        self.old_depth.clear();
        self.old_depth.resize(un, 0.0);
        self.old_volume.clear();
        self.old_volume.resize(un, 0.0);
        self.old_lat_flow.clear();
        self.old_lat_flow.resize(un, 0.0);
        self.old_inflow.clear();
        self.old_inflow.resize(un, 0.0);
        self.comments.clear();
        self.comments.resize(un, Vec::new());
        self.tags.clear();
        self.tags.resize(un, Vec::new());
        self.rpt_flag.clear();
        self.rpt_flag.resize(un, 0);
        self.stat_vol_flooded.clear();
        self.stat_vol_flooded.resize(un, 0.0);
        self.stat_time_flooded.clear();
        self.stat_time_flooded.resize(un, 0.0);
        self.stat_max_depth.clear();
        self.stat_max_depth.resize(un, 0.0);
        self.stat_max_overflow.clear();
        self.stat_max_overflow.resize(un, 0.0);
        self.stat_max_overflow_date.clear();
        self.stat_max_overflow_date.resize(un, 0.0);
        self.stat_sum_depth.clear();
        self.stat_sum_depth.resize(un, 0.0);
        self.stat_sum_volume.clear();
        self.stat_sum_volume.resize(un, 0.0);
        self.stat_storage_max_volume_date.clear();
        self.stat_storage_max_volume_date.resize(un, 0.0);
        self.stat_storage_max_volume.clear();
        self.stat_storage_max_volume.resize(un, 0.0);
        self.stat_max_depth_date.clear();
        self.stat_max_depth_date.resize(un, 0.0);
        self.stat_max_rpt_depth.clear();
        self.stat_max_rpt_depth.resize(un, 0.0);
        self.stat_max_inflow_date.clear();
        self.stat_max_inflow_date.resize(un, 0.0);
        self.stat_time_surcharged.clear();
        self.stat_time_surcharged.resize(un, 0.0);
        self.stat_max_surcharge_height.clear();
        self.stat_max_surcharge_height.resize(un, 0.0);
        self.stat_max_lat_inflow.clear();
        self.stat_max_lat_inflow.resize(un, 0.0);
        self.stat_max_total_inflow.clear();
        self.stat_max_total_inflow.resize(un, 0.0);
        self.stat_storage_max_outflow.clear();
        self.stat_storage_max_outflow.resize(un, 0.0);
        self.stat_lat_inflow_vol.clear();
        self.stat_lat_inflow_vol.resize(un, 0.0);
        self.stat_total_inflow_vol.clear();
        self.stat_total_inflow_vol.resize(un, 0.0);
        self.stat_total_outflow_vol.clear();
        self.stat_total_outflow_vol.resize(un, 0.0);
        self.stat_outfall_avg_flow.clear();
        self.stat_outfall_avg_flow.resize(un, 0.0);
        self.stat_outfall_max_flow.clear();
        self.stat_outfall_max_flow.resize(un, 0.0);
        self.stat_outfall_periods.clear();
        self.stat_outfall_periods.resize(un, 0);
        self.stat_non_converged_count.clear();
        self.stat_non_converged_count.resize(un, 0);
        self.stat_time_courant_critical.clear();
        self.stat_time_courant_critical.resize(un, 0.0);
    }

    pub fn grow_to(&mut self, n: i32) {
        if n <= self.count() {
            return;
        }
        let un = usize::try_from(n).unwrap_or(0);
        self.type_.resize(un, NodeType::Junction);
        self.invert_elev.resize(un, 0.0);
        self.full_depth.resize(un, 0.0);
        self.init_depth.resize(un, 0.0);
        self.sur_depth.resize(un, 0.0);
        self.ponded_area.resize(un, 0.0);
        self.storage_area.resize(un, 0.0);
        self.is_virtual.resize(un, 0);
        self.is_inlet.resize(un, 0);
        self.rim_depth.resize(un, 0.0);
        self.depth.resize(un, 0.0);
        self.head.resize(un, 0.0);
        self.volume.resize(un, 0.0);
        self.lat_flow.resize(un, 0.0);
        self.user_lat_flow.resize(un, 0.0);
        self.runoff_inflow.resize(un, 0.0);
        self.gw_inflow.resize(un, 0.0);
        self.ext_inflow.resize(un, 0.0);
        self.dwf_inflow.resize(un, 0.0);
        self.rdii_inflow.resize(un, 0.0);
        self.iface_inflow.resize(un, 0.0);
        self.coupling_inflow.resize(un, 0.0);
        self.coupling_volume.resize(un, 0.0);
        self.coupling_queue.resize(un, 0.0);
        self.coupling_age_vol_queue.resize(un, 0.0);
        self.coupling_temp_vol_queue.resize(un, 0.0);
        self.coupling_age_vol_inflow.resize(un, 0.0);
        self.coupling_temp_vol_inflow.resize(un, 0.0);
        self.qual_vol_in.resize(un, 0.0);
        self.lid_drain_qual_vol.resize(un, 0.0);
        self.lid_drain_qual_vol_old.resize(un, 0.0);
        self.lid_drain_inflow.resize(un, 0.0);
        self.inflow.resize(un, 0.0);
        self.outflow.resize(un, 0.0);
        self.overflow.resize(un, 0.0);
        self.losses.resize(un, 0.0);
        self.crown_elev.resize(un, 0.0);
        self.degree.resize(un, 0);
        self.old_net_inflow.resize(un, 0.0);
        self.full_volume.resize(un, 0.0);
        self.rpt_full_volume.resize(un, 0.0);
        self.old_depth.resize(un, 0.0);
        self.old_volume.resize(un, 0.0);
        self.old_lat_flow.resize(un, 0.0);
        self.old_inflow.resize(un, 0.0);
        self.comments.resize(un, Vec::new());
        self.tags.resize(un, Vec::new());
        self.rpt_flag.resize(un, 0);
        self.stat_vol_flooded.resize(un, 0.0);
        self.stat_time_flooded.resize(un, 0.0);
        self.stat_max_depth.resize(un, 0.0);
        self.stat_max_overflow.resize(un, 0.0);
        self.stat_max_overflow_date.resize(un, 0.0);
        self.stat_sum_depth.resize(un, 0.0);
        self.stat_sum_volume.resize(un, 0.0);
        self.stat_storage_max_volume_date.resize(un, 0.0);
        self.stat_storage_max_volume.resize(un, 0.0);
        self.stat_max_depth_date.resize(un, 0.0);
        self.stat_max_rpt_depth.resize(un, 0.0);
        self.stat_max_inflow_date.resize(un, 0.0);
        self.stat_time_surcharged.resize(un, 0.0);
        self.stat_max_surcharge_height.resize(un, 0.0);
        self.stat_max_lat_inflow.resize(un, 0.0);
        self.stat_max_total_inflow.resize(un, 0.0);
        self.stat_storage_max_outflow.resize(un, 0.0);
        self.stat_lat_inflow_vol.resize(un, 0.0);
        self.stat_total_inflow_vol.resize(un, 0.0);
        self.stat_total_outflow_vol.resize(un, 0.0);
        self.stat_outfall_avg_flow.resize(un, 0.0);
        self.stat_outfall_max_flow.resize(un, 0.0);
        self.stat_outfall_periods.resize(un, 0);
        self.stat_non_converged_count.resize(un, 0);
        self.stat_time_courant_critical.resize(un, 0.0);
    }

    pub fn reserve_to(&mut self, n: i32) {
        if n <= 0 {
            return;
        }
        let un = usize::try_from(n).unwrap_or(0);
        if self.type_.capacity() >= un {
            return;
        }
        self.type_.reserve(un.saturating_sub(self.type_.len()));
        self.invert_elev.reserve(un.saturating_sub(self.invert_elev.len()));
        self.full_depth.reserve(un.saturating_sub(self.full_depth.len()));
        self.init_depth.reserve(un.saturating_sub(self.init_depth.len()));
        self.sur_depth.reserve(un.saturating_sub(self.sur_depth.len()));
        self.ponded_area.reserve(un.saturating_sub(self.ponded_area.len()));
        self.storage_area.reserve(un.saturating_sub(self.storage_area.len()));
        self.is_virtual.reserve(un.saturating_sub(self.is_virtual.len()));
        self.is_inlet.reserve(un.saturating_sub(self.is_inlet.len()));
        self.rim_depth.reserve(un.saturating_sub(self.rim_depth.len()));
        self.depth.reserve(un.saturating_sub(self.depth.len()));
        self.head.reserve(un.saturating_sub(self.head.len()));
        self.volume.reserve(un.saturating_sub(self.volume.len()));
        self.lat_flow.reserve(un.saturating_sub(self.lat_flow.len()));
        self.user_lat_flow.reserve(un.saturating_sub(self.user_lat_flow.len()));
        self.runoff_inflow.reserve(un.saturating_sub(self.runoff_inflow.len()));
        self.gw_inflow.reserve(un.saturating_sub(self.gw_inflow.len()));
        self.ext_inflow.reserve(un.saturating_sub(self.ext_inflow.len()));
        self.dwf_inflow.reserve(un.saturating_sub(self.dwf_inflow.len()));
        self.rdii_inflow.reserve(un.saturating_sub(self.rdii_inflow.len()));
        self.iface_inflow.reserve(un.saturating_sub(self.iface_inflow.len()));
        self.coupling_inflow.reserve(un.saturating_sub(self.coupling_inflow.len()));
        self.coupling_volume.reserve(un.saturating_sub(self.coupling_volume.len()));
        self.coupling_queue.reserve(un.saturating_sub(self.coupling_queue.len()));
        self.inflow.reserve(un.saturating_sub(self.inflow.len()));
        self.outflow.reserve(un.saturating_sub(self.outflow.len()));
        self.overflow.reserve(un.saturating_sub(self.overflow.len()));
        self.coupling_age_vol_queue.reserve(un.saturating_sub(self.coupling_age_vol_queue.len()));
        self.coupling_temp_vol_queue.reserve(un.saturating_sub(self.coupling_temp_vol_queue.len()));
        self.coupling_age_vol_inflow.reserve(un.saturating_sub(self.coupling_age_vol_inflow.len()));
        self.coupling_temp_vol_inflow.reserve(un.saturating_sub(self.coupling_temp_vol_inflow.len()));
        self.losses.reserve(un.saturating_sub(self.losses.len()));
        self.crown_elev.reserve(un.saturating_sub(self.crown_elev.len()));
        self.degree.reserve(un.saturating_sub(self.degree.len()));
        self.old_net_inflow.reserve(un.saturating_sub(self.old_net_inflow.len()));
        self.full_volume.reserve(un.saturating_sub(self.full_volume.len()));
        self.rpt_full_volume.reserve(un.saturating_sub(self.rpt_full_volume.len()));
        self.old_depth.reserve(un.saturating_sub(self.old_depth.len()));
        self.old_volume.reserve(un.saturating_sub(self.old_volume.len()));
        self.old_lat_flow.reserve(un.saturating_sub(self.old_lat_flow.len()));
        self.old_inflow.reserve(un.saturating_sub(self.old_inflow.len()));
        self.rpt_flag.reserve(un.saturating_sub(self.rpt_flag.len()));
        self.stat_vol_flooded.reserve(un.saturating_sub(self.stat_vol_flooded.len()));
        self.stat_time_flooded.reserve(un.saturating_sub(self.stat_time_flooded.len()));
        self.stat_max_depth.reserve(un.saturating_sub(self.stat_max_depth.len()));
        self.stat_max_overflow.reserve(un.saturating_sub(self.stat_max_overflow.len()));
        self.stat_max_overflow_date.reserve(un.saturating_sub(self.stat_max_overflow_date.len()));
        self.stat_sum_depth.reserve(un.saturating_sub(self.stat_sum_depth.len()));
        self.stat_sum_volume.reserve(un.saturating_sub(self.stat_sum_volume.len()));
        self.stat_storage_max_volume_date.reserve(un.saturating_sub(self.stat_storage_max_volume_date.len()));
        self.stat_storage_max_volume.reserve(un.saturating_sub(self.stat_storage_max_volume.len()));
        self.stat_max_depth_date.reserve(un.saturating_sub(self.stat_max_depth_date.len()));
        self.stat_max_rpt_depth.reserve(un.saturating_sub(self.stat_max_rpt_depth.len()));
        self.stat_max_inflow_date.reserve(un.saturating_sub(self.stat_max_inflow_date.len()));
        self.stat_time_surcharged.reserve(un.saturating_sub(self.stat_time_surcharged.len()));
        self.stat_max_surcharge_height.reserve(un.saturating_sub(self.stat_max_surcharge_height.len()));
        self.stat_max_lat_inflow.reserve(un.saturating_sub(self.stat_max_lat_inflow.len()));
        self.stat_max_total_inflow.reserve(un.saturating_sub(self.stat_max_total_inflow.len()));
        self.stat_storage_max_outflow.reserve(un.saturating_sub(self.stat_storage_max_outflow.len()));
        self.stat_lat_inflow_vol.reserve(un.saturating_sub(self.stat_lat_inflow_vol.len()));
        self.stat_total_inflow_vol.reserve(un.saturating_sub(self.stat_total_inflow_vol.len()));
        self.stat_total_outflow_vol.reserve(un.saturating_sub(self.stat_total_outflow_vol.len()));
        self.stat_outfall_avg_flow.reserve(un.saturating_sub(self.stat_outfall_avg_flow.len()));
        self.stat_outfall_max_flow.reserve(un.saturating_sub(self.stat_outfall_max_flow.len()));
        self.stat_outfall_periods.reserve(un.saturating_sub(self.stat_outfall_periods.len()));
        self.stat_non_converged_count.reserve(un.saturating_sub(self.stat_non_converged_count.len()));
        self.stat_time_courant_critical.reserve(un.saturating_sub(self.stat_time_courant_critical.len()));
        self.qual_vol_in.reserve(un.saturating_sub(self.qual_vol_in.len()));
        self.lid_drain_qual_vol.reserve(un.saturating_sub(self.lid_drain_qual_vol.len()));
        self.lid_drain_inflow.reserve(un.saturating_sub(self.lid_drain_inflow.len()));
        self.comments.reserve(un.saturating_sub(self.comments.len()));
        self.tags.reserve(un.saturating_sub(self.tags.len()));
    }

    pub fn erase_at(&mut self, idx: i32) {
        let ui = idx as usize;
        if ui < self.type_.len() { self.type_.remove(ui); }
        if ui < self.invert_elev.len() { self.invert_elev.remove(ui); }
        if ui < self.full_depth.len() { self.full_depth.remove(ui); }
        if ui < self.init_depth.len() { self.init_depth.remove(ui); }
        if ui < self.sur_depth.len() { self.sur_depth.remove(ui); }
        if ui < self.ponded_area.len() { self.ponded_area.remove(ui); }
        if ui < self.storage_area.len() { self.storage_area.remove(ui); }
        if ui < self.is_virtual.len() { self.is_virtual.remove(ui); }
        if ui < self.is_inlet.len() { self.is_inlet.remove(ui); }
        if ui < self.rim_depth.len() { self.rim_depth.remove(ui); }
        if ui < self.depth.len() { self.depth.remove(ui); }
        if ui < self.head.len() { self.head.remove(ui); }
        if ui < self.volume.len() { self.volume.remove(ui); }
        if ui < self.lat_flow.len() { self.lat_flow.remove(ui); }
        if ui < self.user_lat_flow.len() { self.user_lat_flow.remove(ui); }
        if ui < self.runoff_inflow.len() { self.runoff_inflow.remove(ui); }
        if ui < self.gw_inflow.len() { self.gw_inflow.remove(ui); }
        if ui < self.ext_inflow.len() { self.ext_inflow.remove(ui); }
        if ui < self.dwf_inflow.len() { self.dwf_inflow.remove(ui); }
        if ui < self.rdii_inflow.len() { self.rdii_inflow.remove(ui); }
        if ui < self.iface_inflow.len() { self.iface_inflow.remove(ui); }
        if ui < self.coupling_inflow.len() { self.coupling_inflow.remove(ui); }
        if ui < self.coupling_volume.len() { self.coupling_volume.remove(ui); }
        if ui < self.coupling_queue.len() { self.coupling_queue.remove(ui); }
        if ui < self.coupling_age_vol_queue.len() { self.coupling_age_vol_queue.remove(ui); }
        if ui < self.coupling_temp_vol_queue.len() { self.coupling_temp_vol_queue.remove(ui); }
        if ui < self.coupling_age_vol_inflow.len() { self.coupling_age_vol_inflow.remove(ui); }
        if ui < self.coupling_temp_vol_inflow.len() { self.coupling_temp_vol_inflow.remove(ui); }
        if ui < self.qual_vol_in.len() { self.qual_vol_in.remove(ui); }
        if ui < self.lid_drain_qual_vol.len() { self.lid_drain_qual_vol.remove(ui); }
        if ui < self.lid_drain_inflow.len() { self.lid_drain_inflow.remove(ui); }
        if ui < self.inflow.len() { self.inflow.remove(ui); }
        if ui < self.outflow.len() { self.outflow.remove(ui); }
        if ui < self.overflow.len() { self.overflow.remove(ui); }
        if ui < self.losses.len() { self.losses.remove(ui); }
        if ui < self.crown_elev.len() { self.crown_elev.remove(ui); }
        if ui < self.degree.len() { self.degree.remove(ui); }
        if ui < self.old_net_inflow.len() { self.old_net_inflow.remove(ui); }
        if ui < self.full_volume.len() { self.full_volume.remove(ui); }
        if ui < self.rpt_full_volume.len() { self.rpt_full_volume.remove(ui); }
        if ui < self.old_depth.len() { self.old_depth.remove(ui); }
        if ui < self.old_volume.len() { self.old_volume.remove(ui); }
        if ui < self.old_lat_flow.len() { self.old_lat_flow.remove(ui); }
        if ui < self.old_inflow.len() { self.old_inflow.remove(ui); }
        if ui < self.comments.len() { self.comments.remove(ui); }
        if ui < self.tags.len() { self.tags.remove(ui); }
        if ui < self.rpt_flag.len() { self.rpt_flag.remove(ui); }
        if ui < self.stat_vol_flooded.len() { self.stat_vol_flooded.remove(ui); }
        if ui < self.stat_time_flooded.len() { self.stat_time_flooded.remove(ui); }
        if ui < self.stat_max_depth.len() { self.stat_max_depth.remove(ui); }
        if ui < self.stat_max_overflow.len() { self.stat_max_overflow.remove(ui); }
        if ui < self.stat_max_overflow_date.len() { self.stat_max_overflow_date.remove(ui); }
        if ui < self.stat_sum_depth.len() { self.stat_sum_depth.remove(ui); }
        if ui < self.stat_sum_volume.len() { self.stat_sum_volume.remove(ui); }
        if ui < self.stat_storage_max_volume_date.len() { self.stat_storage_max_volume_date.remove(ui); }
        if ui < self.stat_storage_max_volume.len() { self.stat_storage_max_volume.remove(ui); }
        if ui < self.stat_max_depth_date.len() { self.stat_max_depth_date.remove(ui); }
        if ui < self.stat_max_rpt_depth.len() { self.stat_max_rpt_depth.remove(ui); }
        if ui < self.stat_max_inflow_date.len() { self.stat_max_inflow_date.remove(ui); }
        if ui < self.stat_time_surcharged.len() { self.stat_time_surcharged.remove(ui); }
        if ui < self.stat_max_surcharge_height.len() { self.stat_max_surcharge_height.remove(ui); }
        if ui < self.stat_outfall_avg_flow.len() { self.stat_outfall_avg_flow.remove(ui); }
        if ui < self.stat_max_lat_inflow.len() { self.stat_max_lat_inflow.remove(ui); }
        if ui < self.stat_max_total_inflow.len() { self.stat_max_total_inflow.remove(ui); }
        if ui < self.stat_storage_max_outflow.len() { self.stat_storage_max_outflow.remove(ui); }
        if ui < self.stat_lat_inflow_vol.len() { self.stat_lat_inflow_vol.remove(ui); }
        if ui < self.stat_total_inflow_vol.len() { self.stat_total_inflow_vol.remove(ui); }
        if ui < self.stat_total_outflow_vol.len() { self.stat_total_outflow_vol.remove(ui); }
        if ui < self.stat_outfall_max_flow.len() { self.stat_outfall_max_flow.remove(ui); }
        if ui < self.stat_outfall_periods.len() { self.stat_outfall_periods.remove(ui); }
        if ui < self.stat_non_converged_count.len() { self.stat_non_converged_count.remove(ui); }
        if ui < self.stat_time_courant_critical.len() { self.stat_time_courant_critical.remove(ui); }
        if self.conc_n_pollutants > 0 {
            let np = self.conc_n_pollutants as usize;
            let base = ui * np;
            for v in [
                &mut self.conc,
                &mut self.conc_old,
                &mut self.qual_mass_in,
                &mut self.iface_qual_mass,
                &mut self.coupling_qual_queue,
                &mut self.coupling_qual_inflow,
                &mut self.ext_qual_mass,
                &mut self.dwf_qual_mass,
                &mut self.lid_drain_qual_load,
                &mut self.user_conc_mass_flux,
            ] {
                if base + np <= v.len() {
                    v.drain(base..base + np);
                }
            }
            if ui < self.hrt.len() {
                self.hrt.remove(ui);
            }
        }
        if self.stat_n_pollutants > 0 {
            let np = self.stat_n_pollutants as usize;
            let base = ui * np;
            if base + np <= self.stat_total_load.len() {
                self.stat_total_load.drain(base..base + np);
            }
        }
    }

    pub fn resize_loads(&mut self, n_pollutants: i32) {
        self.stat_n_pollutants = n_pollutants;
        if n_pollutants > 0 {
            let total = self.count() as usize * n_pollutants as usize;
            self.stat_total_load.clear();
            self.stat_total_load.resize(total, 0.0);
        }
    }

    pub fn resize_quality(&mut self, n_pollutants: i32) {
        self.conc_n_pollutants = n_pollutants;
        if n_pollutants > 0 {
        let total = self.count() as usize * n_pollutants as usize;
        self.conc.clear();
        self.conc.resize(total, 0.0);
        self.conc_old.clear();
        self.conc_old.resize(total, 0.0);
        let cnt = self.count() as usize;
        self.hrt.clear();
        self.hrt.resize(cnt, 0.0);
        self.user_conc_mass_flux.clear();
        self.user_conc_mass_flux.resize(total, 0.0);
        self.qual_mass_in.clear();
        self.qual_mass_in.resize(total, 0.0);
        self.iface_qual_mass.clear();
        self.iface_qual_mass.resize(total, 0.0);
        self.coupling_qual_queue.clear();
        self.coupling_qual_queue.resize(total, 0.0);
        self.coupling_qual_inflow.clear();
        self.coupling_qual_inflow.resize(total, 0.0);
        self.ext_qual_mass.clear();
        self.ext_qual_mass.resize(total, 0.0);
        self.dwf_qual_mass.clear();
        self.dwf_qual_mass.resize(total, 0.0);
        self.lid_drain_qual_load.clear();
        self.lid_drain_qual_load.resize(total, 0.0);
        }
    }

    pub fn shrink_to_fit(&mut self) {
        self.type_.shrink_to_fit();
        self.invert_elev.shrink_to_fit();
        self.full_depth.shrink_to_fit();
        self.init_depth.shrink_to_fit();
        self.sur_depth.shrink_to_fit();
        self.ponded_area.shrink_to_fit();
        self.storage_area.shrink_to_fit();
        self.is_virtual.shrink_to_fit();
        self.is_inlet.shrink_to_fit();
        self.rim_depth.shrink_to_fit();
        self.depth.shrink_to_fit();
        self.head.shrink_to_fit();
        self.volume.shrink_to_fit();
        self.lat_flow.shrink_to_fit();
        self.user_lat_flow.shrink_to_fit();
        self.runoff_inflow.shrink_to_fit();
        self.gw_inflow.shrink_to_fit();
        self.ext_inflow.shrink_to_fit();
        self.dwf_inflow.shrink_to_fit();
        self.rdii_inflow.shrink_to_fit();
        self.iface_inflow.shrink_to_fit();
        self.coupling_inflow.shrink_to_fit();
        self.coupling_volume.shrink_to_fit();
        self.coupling_queue.shrink_to_fit();
        self.coupling_age_vol_queue.shrink_to_fit();
        self.coupling_temp_vol_queue.shrink_to_fit();
        self.coupling_age_vol_inflow.shrink_to_fit();
        self.coupling_temp_vol_inflow.shrink_to_fit();
        self.qual_mass_in.shrink_to_fit();
        self.iface_qual_mass.shrink_to_fit();
        self.coupling_qual_queue.shrink_to_fit();
        self.coupling_qual_inflow.shrink_to_fit();
        self.ext_qual_mass.shrink_to_fit();
        self.dwf_qual_mass.shrink_to_fit();
        self.qual_vol_in.shrink_to_fit();
        self.conc.shrink_to_fit();
        self.conc_old.shrink_to_fit();
        self.user_conc_mass_flux.shrink_to_fit();
        self.inflow.shrink_to_fit();
        self.outflow.shrink_to_fit();
        self.overflow.shrink_to_fit();
        self.losses.shrink_to_fit();
        self.crown_elev.shrink_to_fit();
        self.degree.shrink_to_fit();
        self.old_net_inflow.shrink_to_fit();
        self.full_volume.shrink_to_fit();
        self.rpt_full_volume.shrink_to_fit();
        self.old_depth.shrink_to_fit();
        self.old_volume.shrink_to_fit();
        self.old_lat_flow.shrink_to_fit();
        self.old_inflow.shrink_to_fit();
        self.comments.shrink_to_fit();
        self.tags.shrink_to_fit();
        self.rpt_flag.shrink_to_fit();
        self.stat_vol_flooded.shrink_to_fit();
        self.stat_time_flooded.shrink_to_fit();
        self.stat_max_depth.shrink_to_fit();
        self.stat_max_overflow.shrink_to_fit();
        self.stat_max_overflow_date.shrink_to_fit();
        self.stat_sum_depth.shrink_to_fit();
        self.stat_max_depth_date.shrink_to_fit();
        self.stat_max_rpt_depth.shrink_to_fit();
        self.stat_max_inflow_date.shrink_to_fit();
        self.stat_time_surcharged.shrink_to_fit();
        self.stat_max_surcharge_height.shrink_to_fit();
        self.stat_max_lat_inflow.shrink_to_fit();
        self.stat_max_total_inflow.shrink_to_fit();
        self.stat_storage_max_outflow.shrink_to_fit();
        self.stat_lat_inflow_vol.shrink_to_fit();
        self.stat_total_inflow_vol.shrink_to_fit();
        self.stat_total_outflow_vol.shrink_to_fit();
        self.stat_outfall_avg_flow.shrink_to_fit();
        self.stat_outfall_max_flow.shrink_to_fit();
        self.stat_outfall_periods.shrink_to_fit();
        self.stat_total_load.shrink_to_fit();
    }

    pub fn save_state(&mut self) {
        self.save_hyd_state();
        self.save_lat_qual_state();
    }

    pub fn save_hyd_state(&mut self) {
        { let n = self.depth.len(); self.old_depth[..n].copy_from_slice(&self.depth); }
        { let n = self.volume.len(); self.old_volume[..n].copy_from_slice(&self.volume); }
        { let n = self.inflow.len(); self.old_inflow[..n].copy_from_slice(&self.inflow); }
        for i in 0..self.inflow.len() {
            self.old_net_inflow[i] = self.inflow[i] - self.outflow[i];
        }
    }

    pub fn save_lat_qual_state(&mut self) {
        { let n = self.lat_flow.len(); self.old_lat_flow[..n].copy_from_slice(&self.lat_flow); }
        { let n = self.conc.len(); self.conc_old[..n].copy_from_slice(&self.conc); }
    }

    pub fn reset_state(&mut self) {
        let n = self.depth.len();
        for i in 0..n {
            self.depth[i] = self.init_depth[i];
            self.old_depth[i] = self.init_depth[i];
            self.head[i] = self.invert_elev[i] + self.init_depth[i];
        }
        self.volume.fill(0.0);
        self.old_volume.fill(0.0);
        self.lat_flow.fill(0.0);
        self.inflow.fill(0.0);
        self.outflow.fill(0.0);
        self.overflow.fill(0.0);
        self.losses.fill(0.0);
        self.old_net_inflow.fill(0.0);
        self.old_lat_flow.fill(0.0);
        self.old_inflow.fill(0.0);
        self.coupling_inflow.fill(0.0);
        self.coupling_volume.fill(0.0);
        self.coupling_queue.fill(0.0);
        self.coupling_qual_queue.fill(0.0);
        self.coupling_qual_inflow.fill(0.0);
        self.coupling_age_vol_queue.fill(0.0);
        self.coupling_temp_vol_queue.fill(0.0);
        self.coupling_age_vol_inflow.fill(0.0);
        self.coupling_temp_vol_inflow.fill(0.0);
        self.clear_inflow_sources();
        self.conc.fill(0.0);
        self.conc_old.fill(0.0);
    }

    pub fn clear_inflow_sources(&mut self) {
        self.ext_inflow.fill(0.0);
        self.dwf_inflow.fill(0.0);
        self.rdii_inflow.fill(0.0);
        self.iface_inflow.fill(0.0);
        self.iface_qual_mass.fill(0.0);
        self.ext_qual_mass.fill(0.0);
        self.dwf_qual_mass.fill(0.0);
        self.qual_mass_in.fill(0.0);
        self.qual_vol_in.fill(0.0);
    }

}

