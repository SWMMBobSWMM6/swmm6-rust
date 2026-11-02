// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/SubcatchData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson (first pass by tools/soa2rs.py, finished by hand). Not an
// OpenSWMM release.
//
// The subcatchment store, one Vec per attribute (structure of arrays), sized and cleared
// exactly as upstream sizes and clears each one. resize() assigns every base column
// except the three parse-time name columns (outlet_name, snowpack_name, gage_name), which
// upstream resizes and so keeps their prefixes; grow_to() extends the same set; the
// per-pollutant matrices wait for resize_quality()/resize_washoff_load()/
// resize_total_load(), and the land-use matrices for resize_coverage(). sweep_now and
// sweep_last_date are sized by no method here, as upstream. `std::vector<char>` is
// Vec<i8>, `std::vector<unsigned char>` Vec<u8>, `std::vector<long>` Vec<i64>.
// std::numeric_limits<double>::max() is f64::MAX. Products of two casts to size_t use
// wrapping_mul, as size_t arithmetic wraps. std::copy into a longer destination copies
// the prefix, as here; into a shorter one it is undefined in C++ and a panic here.

#![allow(clippy::all)]
use crate::core::options::Bytes;

#[derive(Debug, Clone, Default)]
pub struct SubcatchData {
    pub outlet_node: Vec<i32>,
    pub outlet_subcatch: Vec<i32>,
    pub outlet_name: Vec<Bytes>,
    pub snowpack_name: Vec<Bytes>,
    pub gage_name: Vec<Bytes>,
    pub gage: Vec<i32>,
    pub area: Vec<f64>,
    pub width: Vec<f64>,
    pub slope: Vec<f64>,
    pub curb_length: Vec<f64>,
    pub rain_scale_factor: Vec<f64>,
    pub snow_scale_factor: Vec<f64>,
    pub frac_imperv: Vec<f64>,
    pub frac_imperv_no_store: Vec<f64>,
    pub pct_zero: Vec<f64>,
    pub n_imperv: Vec<f64>,
    pub n_perv: Vec<f64>,
    pub ds_imperv: Vec<f64>,
    pub ds_perv: Vec<f64>,
    pub subarea_routing: Vec<i32>,
    pub pct_routed: Vec<f64>,
    pub infil_model: Vec<i32>,
    pub infil_p1: Vec<f64>,
    pub infil_p2: Vec<f64>,
    pub infil_p3: Vec<f64>,
    pub infil_p4: Vec<f64>,
    pub infil_p5: Vec<f64>,
    pub runoff: Vec<f64>,
    pub rainfall: Vec<f64>,
    pub evap_loss: Vec<f64>,
    pub infil_loss: Vec<f64>,
    pub ponded_depth: Vec<f64>,
    pub old_ponded_depth: Vec<f64>,
    pub washoff_vout1: Vec<f64>,
    pub gw_flow: Vec<f64>,
    pub snow_depth: Vec<f64>,
    pub lid_drain_flow: Vec<f64>,
    pub old_runoff: Vec<f64>,
    pub old_gw_flow: Vec<f64>,
    pub old_snow_depth: Vec<f64>,
    pub old_lid_drain_flow: Vec<f64>,
    pub runon_inflow: Vec<f64>,
    pub old_runon_inflow: Vec<f64>,
    pub runon_rate: Vec<f64>,
    pub outfall_runon_vol: Vec<f64>,
    pub gw_sw_head: Vec<f64>,
    pub gw_node_avail_flow: Vec<f64>,
    pub gw_max_infil_vol: Vec<f64>,
    pub washoff_load: Vec<f64>,
    pub conc: Vec<f64>,
    pub conc_old: Vec<f64>,
    pub ponded_qual: Vec<f64>,
    pub runon_qual_rate: Vec<f64>,
    pub lid_drain_runon_qual: Vec<f64>,
    pub subarea_runoff_rate: Vec<f64>,
    pub init_loading: Vec<f64>,
    pub conc_n_pollutants: i32,
    pub comments: Vec<Bytes>,
    pub tags: Vec<Bytes>,
    pub rpt_flag: Vec<i8>,
    pub stat_precip_vol: Vec<f64>,
    pub stat_runon_vol: Vec<f64>,
    pub stat_gw_time: Vec<f64>,
    pub stat_evap_vol: Vec<f64>,
    pub stat_infil_vol: Vec<f64>,
    pub stat_imperv_vol: Vec<f64>,
    pub stat_perv_vol: Vec<f64>,
    pub stat_runoff_vol: Vec<f64>,
    pub stat_max_runoff: Vec<f64>,
    pub stat_gw_infil_vol: Vec<f64>,
    pub stat_gw_upper_evap_vol: Vec<f64>,
    pub stat_gw_lower_evap_vol: Vec<f64>,
    pub stat_gw_deep_perc_vol: Vec<f64>,
    pub stat_gw_flow_vol: Vec<f64>,
    pub stat_gw_max_flow: Vec<f64>,
    pub stat_gw_sum_theta: Vec<f64>,
    pub stat_gw_sum_depth: Vec<f64>,
    pub stat_gw_final_theta: Vec<f64>,
    pub stat_gw_final_depth: Vec<f64>,
    pub stat_gw_steps: Vec<i64>,
    pub gw_aquifer: Vec<i32>,
    pub gw_node: Vec<i32>,
    pub gw_surf_elev: Vec<f64>,
    pub gw_a1: Vec<f64>,
    pub gw_b1: Vec<f64>,
    pub gw_a2: Vec<f64>,
    pub gw_b2: Vec<f64>,
    pub gw_a3: Vec<f64>,
    pub gw_tw: Vec<f64>,
    pub gw_hstar: Vec<f64>,
    pub gw_bot_elev: Vec<f64>,
    pub gw_wt_elev: Vec<f64>,
    pub gw_upper_moist: Vec<f64>,
    pub snowpack: Vec<i32>,
    pub snow_net_imperv: Vec<f64>,
    pub snow_net_perv: Vec<f64>,
    pub snow_melt_imperv: Vec<f64>,
    pub snow_melt_perv: Vec<f64>,
    pub snow_melt_age_imperv: Vec<f64>,
    pub snow_melt_age_perv: Vec<f64>,
    pub total_lid_area_ft2: Vec<f64>,
    pub lid_return_to_perv_cfs: Vec<f64>,
    pub lid_drain_runon_cfs: Vec<f64>,
    pub total_load: Vec<f64>,
    pub total_load_n_pollutants: i32,
    pub coverage: Vec<f64>,
    pub coverage_n_landuses: i32,
    pub sweep_last_swept: Vec<f64>,
    pub sweep_now: Vec<u8>,
    pub sweep_last_date: Vec<f64>,
}

// The member initializers (conc_n_pollutants = total_load_n_pollutants =
// coverage_n_landuses = 0) are Rust's defaults, so Default is derived.

impl SubcatchData {
    pub fn resize_total_load(&mut self, n_sc: i32, n_poll: i32) {
        self.total_load_n_pollutants = n_poll;
        self.total_load.clear();
        self.total_load.resize((n_sc as usize).wrapping_mul(n_poll as usize), 0.0);
    }

    pub fn resize_coverage(&mut self, n_sc: i32, n_lu: i32) {
        self.coverage_n_landuses = n_lu;
        let sz = (n_sc as usize).wrapping_mul(n_lu as usize);
        self.coverage.clear();
        self.coverage.resize(sz, 0.0);
        self.sweep_last_swept.clear();
        self.sweep_last_swept.resize(sz, 0.0);
    }

    pub fn count(&self) -> i32 {
        self.area.len() as i32
    }

    pub fn resize(&mut self, n: i32) {
        let un = usize::try_from(n).unwrap_or(0);
        self.outlet_node.clear();
        self.outlet_node.resize(un, -1);
        self.outlet_subcatch.clear();
        self.outlet_subcatch.resize(un, -1);
        self.outlet_name.resize(un, Vec::new());
        self.snowpack_name.resize(un, Vec::new());
        self.gage_name.resize(un, Vec::new());
        self.gage.clear();
        self.gage.resize(un, -1);
        self.area.clear();
        self.area.resize(un, 0.0);
        self.width.clear();
        self.width.resize(un, 0.0);
        self.slope.clear();
        self.slope.resize(un, 0.0);
        self.curb_length.clear();
        self.curb_length.resize(un, 0.0);
        self.rain_scale_factor.clear();
        self.rain_scale_factor.resize(un, 1.0);
        self.snow_scale_factor.clear();
        self.snow_scale_factor.resize(un, 1.0);
        self.frac_imperv.clear();
        self.frac_imperv.resize(un, 0.0);
        self.frac_imperv_no_store.clear();
        self.frac_imperv_no_store.resize(un, 0.0);
        self.pct_zero.clear();
        self.pct_zero.resize(un, 0.0);
        self.subarea_runoff_rate.clear();
        self.subarea_runoff_rate.resize(un, 0.0);
        self.n_imperv.clear();
        self.n_imperv.resize(un, 0.013);
        self.n_perv.clear();
        self.n_perv.resize(un, 0.1);
        self.ds_imperv.clear();
        self.ds_imperv.resize(un, 0.0);
        self.ds_perv.clear();
        self.ds_perv.resize(un, 0.0);
        self.subarea_routing.clear();
        self.subarea_routing.resize(un, 0);
        self.pct_routed.clear();
        self.pct_routed.resize(un, 0.0);
        self.infil_model.clear();
        self.infil_model.resize(un, 0);
        self.infil_p1.clear();
        self.infil_p1.resize(un, 0.0);
        self.infil_p2.clear();
        self.infil_p2.resize(un, 0.0);
        self.infil_p3.clear();
        self.infil_p3.resize(un, 0.0);
        self.infil_p4.clear();
        self.infil_p4.resize(un, 0.0);
        self.infil_p5.clear();
        self.infil_p5.resize(un, 0.0);
        self.runoff.clear();
        self.runoff.resize(un, 0.0);
        self.rainfall.clear();
        self.rainfall.resize(un, 0.0);
        self.evap_loss.clear();
        self.evap_loss.resize(un, 0.0);
        self.infil_loss.clear();
        self.infil_loss.resize(un, 0.0);
        self.ponded_depth.clear();
        self.ponded_depth.resize(un, 0.0);
        self.old_ponded_depth.clear();
        self.old_ponded_depth.resize(un, 0.0);
        self.washoff_vout1.clear();
        self.washoff_vout1.resize(un, 0.0);
        self.gw_flow.clear();
        self.gw_flow.resize(un, 0.0);
        self.snow_depth.clear();
        self.snow_depth.resize(un, 0.0);
        self.lid_drain_flow.clear();
        self.lid_drain_flow.resize(un, 0.0);
        self.old_runoff.clear();
        self.old_runoff.resize(un, 0.0);
        self.old_gw_flow.clear();
        self.old_gw_flow.resize(un, 0.0);
        self.old_snow_depth.clear();
        self.old_snow_depth.resize(un, 0.0);
        self.old_lid_drain_flow.clear();
        self.old_lid_drain_flow.resize(un, 0.0);
        self.runon_inflow.clear();
        self.runon_inflow.resize(un, 0.0);
        self.old_runon_inflow.clear();
        self.old_runon_inflow.resize(un, 0.0);
        self.runon_rate.clear();
        self.runon_rate.resize(un, 0.0);
        self.outfall_runon_vol.clear();
        self.outfall_runon_vol.resize(un, 0.0);
        self.gw_sw_head.clear();
        self.gw_sw_head.resize(un, 0.0);
        self.gw_node_avail_flow.clear();
        self.gw_node_avail_flow.resize(un, 0.0);
        self.gw_max_infil_vol.clear();
        self.gw_max_infil_vol.resize(un, f64::MAX);
        self.comments.clear();
        self.comments.resize(un, Vec::new());
        self.tags.clear();
        self.tags.resize(un, Vec::new());
        self.rpt_flag.clear();
        self.rpt_flag.resize(un, 0);
        self.stat_precip_vol.clear();
        self.stat_precip_vol.resize(un, 0.0);
        self.stat_gw_time.clear();
        self.stat_gw_time.resize(un, 0.0);
        self.stat_runon_vol.clear();
        self.stat_runon_vol.resize(un, 0.0);
        self.stat_evap_vol.clear();
        self.stat_evap_vol.resize(un, 0.0);
        self.stat_infil_vol.clear();
        self.stat_infil_vol.resize(un, 0.0);
        self.stat_imperv_vol.clear();
        self.stat_imperv_vol.resize(un, 0.0);
        self.stat_perv_vol.clear();
        self.stat_perv_vol.resize(un, 0.0);
        self.stat_runoff_vol.clear();
        self.stat_runoff_vol.resize(un, 0.0);
        self.stat_max_runoff.clear();
        self.stat_max_runoff.resize(un, 0.0);
        self.stat_gw_infil_vol.clear();
        self.stat_gw_infil_vol.resize(un, 0.0);
        self.stat_gw_upper_evap_vol.clear();
        self.stat_gw_upper_evap_vol.resize(un, 0.0);
        self.stat_gw_lower_evap_vol.clear();
        self.stat_gw_lower_evap_vol.resize(un, 0.0);
        self.stat_gw_deep_perc_vol.clear();
        self.stat_gw_deep_perc_vol.resize(un, 0.0);
        self.stat_gw_flow_vol.clear();
        self.stat_gw_flow_vol.resize(un, 0.0);
        self.stat_gw_max_flow.clear();
        self.stat_gw_max_flow.resize(un, 0.0);
        self.stat_gw_sum_theta.clear();
        self.stat_gw_sum_theta.resize(un, 0.0);
        self.stat_gw_sum_depth.clear();
        self.stat_gw_sum_depth.resize(un, 0.0);
        self.stat_gw_final_theta.clear();
        self.stat_gw_final_theta.resize(un, 0.0);
        self.stat_gw_final_depth.clear();
        self.stat_gw_final_depth.resize(un, 0.0);
        self.stat_gw_steps.clear();
        self.stat_gw_steps.resize(un, 0);
        self.gw_aquifer.clear();
        self.gw_aquifer.resize(un, -1);
        self.gw_node.clear();
        self.gw_node.resize(un, -1);
        self.gw_surf_elev.clear();
        self.gw_surf_elev.resize(un, 0.0);
        self.gw_a1.clear();
        self.gw_a1.resize(un, 0.0);
        self.gw_b1.clear();
        self.gw_b1.resize(un, 0.0);
        self.gw_a2.clear();
        self.gw_a2.resize(un, 0.0);
        self.gw_b2.clear();
        self.gw_b2.resize(un, 0.0);
        self.gw_a3.clear();
        self.gw_a3.resize(un, 0.0);
        self.gw_tw.clear();
        self.gw_tw.resize(un, 0.0);
        self.gw_hstar.clear();
        self.gw_hstar.resize(un, -1.0e10);
        self.gw_bot_elev.clear();
        self.gw_bot_elev.resize(un, -1.0e10);
        self.gw_wt_elev.clear();
        self.gw_wt_elev.resize(un, -1.0e10);
        self.gw_upper_moist.clear();
        self.gw_upper_moist.resize(un, -1.0e10);
        self.snowpack.clear();
        self.snowpack.resize(un, -1);
        self.snow_net_imperv.clear();
        self.snow_net_imperv.resize(un, -1.0);
        self.snow_net_perv.clear();
        self.snow_net_perv.resize(un, -1.0);
        self.snow_melt_imperv.clear();
        self.snow_melt_imperv.resize(un, -1.0);
        self.snow_melt_perv.clear();
        self.snow_melt_perv.resize(un, -1.0);
        self.snow_melt_age_imperv.clear();
        self.snow_melt_age_imperv.resize(un, -1.0);
        self.snow_melt_age_perv.clear();
        self.snow_melt_age_perv.resize(un, -1.0);
        self.total_lid_area_ft2.clear();
        self.total_lid_area_ft2.resize(un, 0.0);
        self.lid_return_to_perv_cfs.clear();
        self.lid_return_to_perv_cfs.resize(un, 0.0);
        self.lid_drain_runon_cfs.clear();
        self.lid_drain_runon_cfs.resize(un, 0.0);
    }

    pub fn grow_to(&mut self, n: i32) {
        if n <= self.count() {
            return;
        }
        let un = usize::try_from(n).unwrap_or(0);
        self.outlet_node.resize(un, -1);
        self.outlet_subcatch.resize(un, -1);
        self.outlet_name.resize(un, Vec::new());
        self.snowpack_name.resize(un, Vec::new());
        self.gage_name.resize(un, Vec::new());
        self.gage.resize(un, -1);
        self.area.resize(un, 0.0);
        self.width.resize(un, 0.0);
        self.slope.resize(un, 0.0);
        self.curb_length.resize(un, 0.0);
        self.rain_scale_factor.resize(un, 1.0);
        self.snow_scale_factor.resize(un, 1.0);
        self.frac_imperv.resize(un, 0.0);
        self.frac_imperv_no_store.resize(un, 0.0);
        self.pct_zero.resize(un, 0.0);
        self.subarea_runoff_rate.resize(un, 0.0);
        self.n_imperv.resize(un, 0.013);
        self.n_perv.resize(un, 0.1);
        self.ds_imperv.resize(un, 0.0);
        self.ds_perv.resize(un, 0.0);
        self.subarea_routing.resize(un, 0);
        self.pct_routed.resize(un, 0.0);
        self.infil_model.resize(un, 0);
        self.infil_p1.resize(un, 0.0);
        self.infil_p2.resize(un, 0.0);
        self.infil_p3.resize(un, 0.0);
        self.infil_p4.resize(un, 0.0);
        self.infil_p5.resize(un, 0.0);
        self.runoff.resize(un, 0.0);
        self.rainfall.resize(un, 0.0);
        self.evap_loss.resize(un, 0.0);
        self.infil_loss.resize(un, 0.0);
        self.ponded_depth.resize(un, 0.0);
        self.old_ponded_depth.resize(un, 0.0);
        self.washoff_vout1.resize(un, 0.0);
        self.gw_flow.resize(un, 0.0);
        self.snow_depth.resize(un, 0.0);
        self.lid_drain_flow.resize(un, 0.0);
        self.old_runoff.resize(un, 0.0);
        self.old_gw_flow.resize(un, 0.0);
        self.old_snow_depth.resize(un, 0.0);
        self.old_lid_drain_flow.resize(un, 0.0);
        self.runon_inflow.resize(un, 0.0);
        self.old_runon_inflow.resize(un, 0.0);
        self.runon_rate.resize(un, 0.0);
        self.gw_sw_head.resize(un, 0.0);
        self.gw_node_avail_flow.resize(un, 0.0);
        self.gw_max_infil_vol.resize(un, f64::MAX);
        self.outfall_runon_vol.resize(un, 0.0);
        self.comments.resize(un, Vec::new());
        self.tags.resize(un, Vec::new());
        self.rpt_flag.resize(un, 0);
        self.stat_precip_vol.resize(un, 0.0);
        self.stat_gw_time.resize(un, 0.0);
        self.stat_runon_vol.resize(un, 0.0);
        self.stat_evap_vol.resize(un, 0.0);
        self.stat_infil_vol.resize(un, 0.0);
        self.stat_imperv_vol.resize(un, 0.0);
        self.stat_perv_vol.resize(un, 0.0);
        self.stat_runoff_vol.resize(un, 0.0);
        self.stat_max_runoff.resize(un, 0.0);
        self.stat_gw_infil_vol.resize(un, 0.0);
        self.stat_gw_upper_evap_vol.resize(un, 0.0);
        self.stat_gw_lower_evap_vol.resize(un, 0.0);
        self.stat_gw_deep_perc_vol.resize(un, 0.0);
        self.stat_gw_flow_vol.resize(un, 0.0);
        self.stat_gw_max_flow.resize(un, 0.0);
        self.stat_gw_sum_theta.resize(un, 0.0);
        self.stat_gw_sum_depth.resize(un, 0.0);
        self.stat_gw_final_theta.resize(un, 0.0);
        self.stat_gw_final_depth.resize(un, 0.0);
        self.stat_gw_steps.resize(un, 0);
        self.gw_aquifer.resize(un, -1);
        self.gw_node.resize(un, -1);
        self.gw_surf_elev.resize(un, 0.0);
        self.gw_a1.resize(un, 0.0);
        self.gw_b1.resize(un, 0.0);
        self.gw_a2.resize(un, 0.0);
        self.gw_b2.resize(un, 0.0);
        self.gw_a3.resize(un, 0.0);
        self.gw_tw.resize(un, 0.0);
        self.gw_hstar.resize(un, -1.0e10);
        self.gw_bot_elev.resize(un, -1.0e10);
        self.gw_wt_elev.resize(un, -1.0e10);
        self.gw_upper_moist.resize(un, -1.0e10);
        self.snowpack.resize(un, -1);
        self.snow_net_imperv.resize(un, -1.0);
        self.snow_net_perv.resize(un, -1.0);
        self.snow_melt_imperv.resize(un, -1.0);
        self.snow_melt_perv.resize(un, -1.0);
        self.snow_melt_age_imperv.resize(un, -1.0);
        self.snow_melt_age_perv.resize(un, -1.0);
        self.total_lid_area_ft2.resize(un, 0.0);
        self.lid_return_to_perv_cfs.resize(un, 0.0);
        self.lid_drain_runon_cfs.resize(un, 0.0);
    }

    pub fn reserve_to(&mut self, n: i32) {
        if n <= 0 {
            return;
        }
        let un = usize::try_from(n).unwrap_or(0);
        if self.outlet_node.capacity() >= un {
            return;
        }
        self.outlet_node.reserve(un.saturating_sub(self.outlet_node.len()));
        self.outlet_subcatch.reserve(un.saturating_sub(self.outlet_subcatch.len()));
        self.gage.reserve(un.saturating_sub(self.gage.len()));
        self.area.reserve(un.saturating_sub(self.area.len()));
        self.width.reserve(un.saturating_sub(self.width.len()));
        self.slope.reserve(un.saturating_sub(self.slope.len()));
        self.curb_length.reserve(un.saturating_sub(self.curb_length.len()));
        self.rain_scale_factor.reserve(un.saturating_sub(self.rain_scale_factor.len()));
        self.snow_scale_factor.reserve(un.saturating_sub(self.snow_scale_factor.len()));
        self.frac_imperv.reserve(un.saturating_sub(self.frac_imperv.len()));
        self.frac_imperv_no_store.reserve(un.saturating_sub(self.frac_imperv_no_store.len()));
        self.pct_zero.reserve(un.saturating_sub(self.pct_zero.len()));
        self.subarea_runoff_rate.reserve(un.saturating_sub(self.subarea_runoff_rate.len()));
        self.n_imperv.reserve(un.saturating_sub(self.n_imperv.len()));
        self.n_perv.reserve(un.saturating_sub(self.n_perv.len()));
        self.ds_imperv.reserve(un.saturating_sub(self.ds_imperv.len()));
        self.ds_perv.reserve(un.saturating_sub(self.ds_perv.len()));
        self.subarea_routing.reserve(un.saturating_sub(self.subarea_routing.len()));
        self.pct_routed.reserve(un.saturating_sub(self.pct_routed.len()));
        self.infil_model.reserve(un.saturating_sub(self.infil_model.len()));
        self.infil_p1.reserve(un.saturating_sub(self.infil_p1.len()));
        self.infil_p2.reserve(un.saturating_sub(self.infil_p2.len()));
        self.infil_p3.reserve(un.saturating_sub(self.infil_p3.len()));
        self.infil_p4.reserve(un.saturating_sub(self.infil_p4.len()));
        self.infil_p5.reserve(un.saturating_sub(self.infil_p5.len()));
        self.runoff.reserve(un.saturating_sub(self.runoff.len()));
        self.rainfall.reserve(un.saturating_sub(self.rainfall.len()));
        self.evap_loss.reserve(un.saturating_sub(self.evap_loss.len()));
        self.infil_loss.reserve(un.saturating_sub(self.infil_loss.len()));
        self.ponded_depth.reserve(un.saturating_sub(self.ponded_depth.len()));
        self.gw_flow.reserve(un.saturating_sub(self.gw_flow.len()));
        self.snow_depth.reserve(un.saturating_sub(self.snow_depth.len()));
        self.lid_drain_flow.reserve(un.saturating_sub(self.lid_drain_flow.len()));
        self.old_runoff.reserve(un.saturating_sub(self.old_runoff.len()));
        self.old_gw_flow.reserve(un.saturating_sub(self.old_gw_flow.len()));
        self.old_snow_depth.reserve(un.saturating_sub(self.old_snow_depth.len()));
        self.old_lid_drain_flow.reserve(un.saturating_sub(self.old_lid_drain_flow.len()));
        self.runon_inflow.reserve(un.saturating_sub(self.runon_inflow.len()));
        self.old_runon_inflow.reserve(un.saturating_sub(self.old_runon_inflow.len()));
        self.runon_rate.reserve(un.saturating_sub(self.runon_rate.len()));
        self.gw_sw_head.reserve(un.saturating_sub(self.gw_sw_head.len()));
        self.gw_node_avail_flow.reserve(un.saturating_sub(self.gw_node_avail_flow.len()));
        self.gw_max_infil_vol.reserve(un.saturating_sub(self.gw_max_infil_vol.len()));
        self.outfall_runon_vol.reserve(un.saturating_sub(self.outfall_runon_vol.len()));
        self.rpt_flag.reserve(un.saturating_sub(self.rpt_flag.len()));
        self.stat_precip_vol.reserve(un.saturating_sub(self.stat_precip_vol.len()));
        self.stat_gw_time.reserve(un.saturating_sub(self.stat_gw_time.len()));
        self.stat_runon_vol.reserve(un.saturating_sub(self.stat_runon_vol.len()));
        self.stat_evap_vol.reserve(un.saturating_sub(self.stat_evap_vol.len()));
        self.stat_infil_vol.reserve(un.saturating_sub(self.stat_infil_vol.len()));
        self.stat_imperv_vol.reserve(un.saturating_sub(self.stat_imperv_vol.len()));
        self.stat_perv_vol.reserve(un.saturating_sub(self.stat_perv_vol.len()));
        self.stat_runoff_vol.reserve(un.saturating_sub(self.stat_runoff_vol.len()));
        self.stat_max_runoff.reserve(un.saturating_sub(self.stat_max_runoff.len()));
        self.stat_gw_infil_vol.reserve(un.saturating_sub(self.stat_gw_infil_vol.len()));
        self.stat_gw_upper_evap_vol.reserve(un.saturating_sub(self.stat_gw_upper_evap_vol.len()));
        self.stat_gw_lower_evap_vol.reserve(un.saturating_sub(self.stat_gw_lower_evap_vol.len()));
        self.stat_gw_deep_perc_vol.reserve(un.saturating_sub(self.stat_gw_deep_perc_vol.len()));
        self.stat_gw_flow_vol.reserve(un.saturating_sub(self.stat_gw_flow_vol.len()));
        self.stat_gw_max_flow.reserve(un.saturating_sub(self.stat_gw_max_flow.len()));
        self.stat_gw_sum_theta.reserve(un.saturating_sub(self.stat_gw_sum_theta.len()));
        self.stat_gw_sum_depth.reserve(un.saturating_sub(self.stat_gw_sum_depth.len()));
        self.stat_gw_final_theta.reserve(un.saturating_sub(self.stat_gw_final_theta.len()));
        self.stat_gw_final_depth.reserve(un.saturating_sub(self.stat_gw_final_depth.len()));
        self.gw_aquifer.reserve(un.saturating_sub(self.gw_aquifer.len()));
        self.gw_node.reserve(un.saturating_sub(self.gw_node.len()));
        self.gw_surf_elev.reserve(un.saturating_sub(self.gw_surf_elev.len()));
        self.gw_a1.reserve(un.saturating_sub(self.gw_a1.len()));
        self.gw_b1.reserve(un.saturating_sub(self.gw_b1.len()));
        self.gw_a2.reserve(un.saturating_sub(self.gw_a2.len()));
        self.gw_b2.reserve(un.saturating_sub(self.gw_b2.len()));
        self.gw_a3.reserve(un.saturating_sub(self.gw_a3.len()));
        self.gw_tw.reserve(un.saturating_sub(self.gw_tw.len()));
        self.gw_hstar.reserve(un.saturating_sub(self.gw_hstar.len()));
        self.gw_bot_elev.reserve(un.saturating_sub(self.gw_bot_elev.len()));
        self.gw_wt_elev.reserve(un.saturating_sub(self.gw_wt_elev.len()));
        self.gw_upper_moist.reserve(un.saturating_sub(self.gw_upper_moist.len()));
        self.snowpack.reserve(un.saturating_sub(self.snowpack.len()));
        self.snow_net_imperv.reserve(un.saturating_sub(self.snow_net_imperv.len()));
        self.snow_net_perv.reserve(un.saturating_sub(self.snow_net_perv.len()));
        self.snow_melt_imperv.reserve(un.saturating_sub(self.snow_melt_imperv.len()));
        self.snow_melt_perv.reserve(un.saturating_sub(self.snow_melt_perv.len()));
        self.snow_melt_age_imperv.reserve(un.saturating_sub(self.snow_melt_age_imperv.len()));
        self.snow_melt_age_perv.reserve(un.saturating_sub(self.snow_melt_age_perv.len()));
        self.total_lid_area_ft2.reserve(un.saturating_sub(self.total_lid_area_ft2.len()));
        self.lid_return_to_perv_cfs.reserve(un.saturating_sub(self.lid_return_to_perv_cfs.len()));
        self.lid_drain_runon_cfs.reserve(un.saturating_sub(self.lid_drain_runon_cfs.len()));
        self.outlet_name.reserve(un.saturating_sub(self.outlet_name.len()));
        self.gage_name.reserve(un.saturating_sub(self.gage_name.len()));
        self.comments.reserve(un.saturating_sub(self.comments.len()));
        self.tags.reserve(un.saturating_sub(self.tags.len()));
        self.stat_gw_steps.reserve(un.saturating_sub(self.stat_gw_steps.len()));
    }

    pub fn erase_at(&mut self, idx: i32) {
        let ui = idx as usize;
        if ui < self.outlet_node.len() { self.outlet_node.remove(ui); }
        if ui < self.outlet_subcatch.len() { self.outlet_subcatch.remove(ui); }
        if ui < self.outlet_name.len() { self.outlet_name.remove(ui); }
        if ui < self.gage.len() { self.gage.remove(ui); }
        if ui < self.area.len() { self.area.remove(ui); }
        if ui < self.width.len() { self.width.remove(ui); }
        if ui < self.slope.len() { self.slope.remove(ui); }
        if ui < self.curb_length.len() { self.curb_length.remove(ui); }
        if ui < self.rain_scale_factor.len() { self.rain_scale_factor.remove(ui); }
        if ui < self.snow_scale_factor.len() { self.snow_scale_factor.remove(ui); }
        if ui < self.frac_imperv.len() { self.frac_imperv.remove(ui); }
        if ui < self.frac_imperv_no_store.len() { self.frac_imperv_no_store.remove(ui); }
        if ui < self.pct_zero.len() { self.pct_zero.remove(ui); }
        if ui < self.subarea_runoff_rate.len() { self.subarea_runoff_rate.remove(ui); }
        if ui < self.n_imperv.len() { self.n_imperv.remove(ui); }
        if ui < self.n_perv.len() { self.n_perv.remove(ui); }
        if ui < self.ds_imperv.len() { self.ds_imperv.remove(ui); }
        if ui < self.ds_perv.len() { self.ds_perv.remove(ui); }
        if ui < self.subarea_routing.len() { self.subarea_routing.remove(ui); }
        if ui < self.pct_routed.len() { self.pct_routed.remove(ui); }
        if ui < self.infil_model.len() { self.infil_model.remove(ui); }
        if ui < self.infil_p1.len() { self.infil_p1.remove(ui); }
        if ui < self.infil_p2.len() { self.infil_p2.remove(ui); }
        if ui < self.infil_p3.len() { self.infil_p3.remove(ui); }
        if ui < self.infil_p4.len() { self.infil_p4.remove(ui); }
        if ui < self.infil_p5.len() { self.infil_p5.remove(ui); }
        if ui < self.runoff.len() { self.runoff.remove(ui); }
        if ui < self.rainfall.len() { self.rainfall.remove(ui); }
        if ui < self.evap_loss.len() { self.evap_loss.remove(ui); }
        if ui < self.infil_loss.len() { self.infil_loss.remove(ui); }
        if ui < self.ponded_depth.len() { self.ponded_depth.remove(ui); }
        if ui < self.gw_flow.len() { self.gw_flow.remove(ui); }
        if ui < self.snow_depth.len() { self.snow_depth.remove(ui); }
        if ui < self.lid_drain_flow.len() { self.lid_drain_flow.remove(ui); }
        if ui < self.old_runoff.len() { self.old_runoff.remove(ui); }
        if ui < self.old_gw_flow.len() { self.old_gw_flow.remove(ui); }
        if ui < self.old_snow_depth.len() { self.old_snow_depth.remove(ui); }
        if ui < self.old_lid_drain_flow.len() { self.old_lid_drain_flow.remove(ui); }
        if ui < self.runon_inflow.len() { self.runon_inflow.remove(ui); }
        if ui < self.old_runon_inflow.len() { self.old_runon_inflow.remove(ui); }
        if ui < self.runon_rate.len() { self.runon_rate.remove(ui); }
        if ui < self.outfall_runon_vol.len() { self.outfall_runon_vol.remove(ui); }
        if ui < self.gw_sw_head.len() { self.gw_sw_head.remove(ui); }
        if ui < self.gw_node_avail_flow.len() { self.gw_node_avail_flow.remove(ui); }
        if ui < self.gw_max_infil_vol.len() { self.gw_max_infil_vol.remove(ui); }
        if ui < self.comments.len() { self.comments.remove(ui); }
        if ui < self.tags.len() { self.tags.remove(ui); }
        if ui < self.rpt_flag.len() { self.rpt_flag.remove(ui); }
        if ui < self.stat_precip_vol.len() { self.stat_precip_vol.remove(ui); }
        if ui < self.stat_gw_time.len() { self.stat_gw_time.remove(ui); }
        if ui < self.stat_runon_vol.len() { self.stat_runon_vol.remove(ui); }
        if ui < self.stat_evap_vol.len() { self.stat_evap_vol.remove(ui); }
        if ui < self.stat_infil_vol.len() { self.stat_infil_vol.remove(ui); }
        if ui < self.stat_imperv_vol.len() { self.stat_imperv_vol.remove(ui); }
        if ui < self.stat_perv_vol.len() { self.stat_perv_vol.remove(ui); }
        if ui < self.stat_runoff_vol.len() { self.stat_runoff_vol.remove(ui); }
        if ui < self.stat_max_runoff.len() { self.stat_max_runoff.remove(ui); }
        if ui < self.stat_gw_infil_vol.len() { self.stat_gw_infil_vol.remove(ui); }
        if ui < self.stat_gw_upper_evap_vol.len() { self.stat_gw_upper_evap_vol.remove(ui); }
        if ui < self.stat_gw_lower_evap_vol.len() { self.stat_gw_lower_evap_vol.remove(ui); }
        if ui < self.stat_gw_deep_perc_vol.len() { self.stat_gw_deep_perc_vol.remove(ui); }
        if ui < self.stat_gw_flow_vol.len() { self.stat_gw_flow_vol.remove(ui); }
        if ui < self.stat_gw_max_flow.len() { self.stat_gw_max_flow.remove(ui); }
        if ui < self.stat_gw_sum_theta.len() { self.stat_gw_sum_theta.remove(ui); }
        if ui < self.stat_gw_sum_depth.len() { self.stat_gw_sum_depth.remove(ui); }
        if ui < self.stat_gw_final_theta.len() { self.stat_gw_final_theta.remove(ui); }
        if ui < self.stat_gw_final_depth.len() { self.stat_gw_final_depth.remove(ui); }
        if ui < self.stat_gw_steps.len() { self.stat_gw_steps.remove(ui); }
        if ui < self.gw_aquifer.len() { self.gw_aquifer.remove(ui); }
        if ui < self.gw_node.len() { self.gw_node.remove(ui); }
        if ui < self.gw_surf_elev.len() { self.gw_surf_elev.remove(ui); }
        if ui < self.gw_a1.len() { self.gw_a1.remove(ui); }
        if ui < self.gw_b1.len() { self.gw_b1.remove(ui); }
        if ui < self.gw_a2.len() { self.gw_a2.remove(ui); }
        if ui < self.gw_b2.len() { self.gw_b2.remove(ui); }
        if ui < self.gw_a3.len() { self.gw_a3.remove(ui); }
        if ui < self.gw_tw.len() { self.gw_tw.remove(ui); }
        if ui < self.gw_hstar.len() { self.gw_hstar.remove(ui); }
        if ui < self.gw_bot_elev.len() { self.gw_bot_elev.remove(ui); }
        if ui < self.gw_wt_elev.len() { self.gw_wt_elev.remove(ui); }
        if ui < self.gw_upper_moist.len() { self.gw_upper_moist.remove(ui); }
        if ui < self.snowpack.len() { self.snowpack.remove(ui); }
        if ui < self.snow_net_imperv.len() { self.snow_net_imperv.remove(ui); }
        if ui < self.snow_net_perv.len() { self.snow_net_perv.remove(ui); }
        if ui < self.snow_melt_imperv.len() { self.snow_melt_imperv.remove(ui); }
        if ui < self.snow_melt_perv.len() { self.snow_melt_perv.remove(ui); }
        if ui < self.snow_melt_age_imperv.len() { self.snow_melt_age_imperv.remove(ui); }
        if ui < self.snow_melt_age_perv.len() { self.snow_melt_age_perv.remove(ui); }
        if ui < self.total_lid_area_ft2.len() { self.total_lid_area_ft2.remove(ui); }
        if ui < self.lid_return_to_perv_cfs.len() { self.lid_return_to_perv_cfs.remove(ui); }
        if ui < self.lid_drain_runon_cfs.len() { self.lid_drain_runon_cfs.remove(ui); }

        // Flat 2D quality arrays: [sc * np + p]
        if self.conc_n_pollutants > 0 {
            let np = self.conc_n_pollutants as usize;
            let base = ui * np;
            for v in [
                &mut self.conc,
                &mut self.conc_old,
                &mut self.init_loading,
                &mut self.ponded_qual,
                &mut self.runon_qual_rate,
                &mut self.lid_drain_runon_qual,
                &mut self.washoff_load,
            ] {
                if base + np <= v.len() {
                    v.drain(base..base + np);
                }
            }
        }

        // Flat 2D total load: [sc * np + p]
        if self.total_load_n_pollutants > 0 {
            let np = self.total_load_n_pollutants as usize;
            let base = ui * np;
            if base + np <= self.total_load.len() {
                self.total_load.drain(base..base + np);
            }
        }

        // Flat 2D coverage/sweep matrices: [sc * n_lu + lu]
        if self.coverage_n_landuses > 0 {
            let nlu = self.coverage_n_landuses as usize;
            let base = ui * nlu;
            for v in [&mut self.coverage, &mut self.sweep_last_swept] {
                if base + nlu <= v.len() {
                    v.drain(base..base + nlu);
                }
            }
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
            self.init_loading.clear();
            self.init_loading.resize(total, 0.0);
            self.ponded_qual.clear();
            self.ponded_qual.resize(total, 0.0);
            self.runon_qual_rate.clear();
            self.runon_qual_rate.resize(total, 0.0);
            self.lid_drain_runon_qual.clear();
            self.lid_drain_runon_qual.resize(total, 0.0);
            self.washoff_load.clear();
            self.washoff_load.resize(total, 0.0);
        }
    }

    pub fn resize_washoff_load(&mut self, n_pollutants: i32) {
        if n_pollutants > 0 {
            self.washoff_load.clear();
            self.washoff_load.resize(self.count() as usize * n_pollutants as usize, 0.0);
        }
    }

    pub fn shrink_to_fit(&mut self) {
        self.outlet_node.shrink_to_fit();
        self.outlet_subcatch.shrink_to_fit();
        self.outlet_name.shrink_to_fit();
        self.gage.shrink_to_fit();
        self.area.shrink_to_fit();
        self.width.shrink_to_fit();
        self.slope.shrink_to_fit();
        self.curb_length.shrink_to_fit();
        self.rain_scale_factor.shrink_to_fit();
        self.snow_scale_factor.shrink_to_fit();
        self.frac_imperv.shrink_to_fit();
        self.frac_imperv_no_store.shrink_to_fit();
        self.pct_zero.shrink_to_fit();
        self.subarea_runoff_rate.shrink_to_fit();
        self.n_imperv.shrink_to_fit();
        self.n_perv.shrink_to_fit();
        self.ds_imperv.shrink_to_fit();
        self.ds_perv.shrink_to_fit();
        self.subarea_routing.shrink_to_fit();
        self.pct_routed.shrink_to_fit();
        self.infil_model.shrink_to_fit();
        self.infil_p1.shrink_to_fit();
        self.infil_p2.shrink_to_fit();
        self.infil_p3.shrink_to_fit();
        self.infil_p4.shrink_to_fit();
        self.infil_p5.shrink_to_fit();
        self.runoff.shrink_to_fit();
        self.rainfall.shrink_to_fit();
        self.evap_loss.shrink_to_fit();
        self.infil_loss.shrink_to_fit();
        self.ponded_depth.shrink_to_fit();
        self.gw_flow.shrink_to_fit();
        self.snow_depth.shrink_to_fit();
        self.lid_drain_flow.shrink_to_fit();
        self.old_runoff.shrink_to_fit();
        self.old_gw_flow.shrink_to_fit();
        self.old_snow_depth.shrink_to_fit();
        self.old_lid_drain_flow.shrink_to_fit();
        self.runon_inflow.shrink_to_fit();
        self.old_runon_inflow.shrink_to_fit();
        self.runon_rate.shrink_to_fit();
        self.outfall_runon_vol.shrink_to_fit();
        self.gw_sw_head.shrink_to_fit();
        self.gw_node_avail_flow.shrink_to_fit();
        self.gw_max_infil_vol.shrink_to_fit();
        self.comments.shrink_to_fit();
        self.tags.shrink_to_fit();
        self.rpt_flag.shrink_to_fit();
        self.stat_precip_vol.shrink_to_fit();
        self.stat_gw_time.shrink_to_fit();
        self.stat_runon_vol.shrink_to_fit();
        self.stat_evap_vol.shrink_to_fit();
        self.stat_infil_vol.shrink_to_fit();
        self.stat_imperv_vol.shrink_to_fit();
        self.stat_perv_vol.shrink_to_fit();
        self.stat_runoff_vol.shrink_to_fit();
        self.stat_max_runoff.shrink_to_fit();
        self.gw_aquifer.shrink_to_fit();
        self.gw_node.shrink_to_fit();
        self.gw_surf_elev.shrink_to_fit();
        self.gw_a1.shrink_to_fit();
        self.gw_b1.shrink_to_fit();
        self.gw_a2.shrink_to_fit();
        self.gw_b2.shrink_to_fit();
        self.gw_a3.shrink_to_fit();
        self.gw_tw.shrink_to_fit();
        self.gw_hstar.shrink_to_fit();
        self.gw_bot_elev.shrink_to_fit();
        self.gw_wt_elev.shrink_to_fit();
        self.gw_upper_moist.shrink_to_fit();
        self.snowpack.shrink_to_fit();
        self.snow_net_imperv.shrink_to_fit();
        self.snow_net_perv.shrink_to_fit();
        self.snow_melt_imperv.shrink_to_fit();
        self.snow_melt_perv.shrink_to_fit();
        self.snow_melt_age_imperv.shrink_to_fit();
        self.snow_melt_age_perv.shrink_to_fit();
        self.conc.shrink_to_fit();
        self.conc_old.shrink_to_fit();
        self.ponded_qual.shrink_to_fit();
        self.runon_qual_rate.shrink_to_fit();
        self.lid_drain_runon_qual.shrink_to_fit();
        self.washoff_load.shrink_to_fit();
        self.total_load.shrink_to_fit();
        self.coverage.shrink_to_fit();
        self.sweep_last_swept.shrink_to_fit();
    }

    pub fn save_state(&mut self) {
        { let n = self.runoff.len(); self.old_runoff[..n].copy_from_slice(&self.runoff); }
        { let n = self.runon_inflow.len(); self.old_runon_inflow[..n].copy_from_slice(&self.runon_inflow); }
        { let n = self.conc.len(); self.conc_old[..n].copy_from_slice(&self.conc); }
        { let n = self.ponded_depth.len(); self.old_ponded_depth[..n].copy_from_slice(&self.ponded_depth); }
    }

    pub fn reset_state(&mut self) {
        self.runoff.fill(0.0);
        self.rainfall.fill(0.0);
        self.snow_net_imperv.fill(-1.0);
        self.snow_net_perv.fill(-1.0);
        self.snow_melt_imperv.fill(-1.0);
        self.snow_melt_perv.fill(-1.0);
        self.snow_melt_age_imperv.fill(-1.0);
        self.snow_melt_age_perv.fill(-1.0);
        self.evap_loss.fill(0.0);
        self.infil_loss.fill(0.0);
        self.ponded_depth.fill(0.0);
        self.old_runoff.fill(0.0);
        self.old_gw_flow.fill(0.0);
        self.snow_depth.fill(0.0);
        self.old_snow_depth.fill(0.0);
        self.lid_drain_flow.fill(0.0);
        self.old_lid_drain_flow.fill(0.0);
        self.runon_inflow.fill(0.0);
        self.old_runon_inflow.fill(0.0);
        self.runon_rate.fill(0.0);
        self.gw_sw_head.fill(0.0);
        self.gw_node_avail_flow.fill(0.0);
        self.washoff_load.fill(0.0);
        self.conc.fill(0.0);
        self.conc_old.fill(0.0);
    }
}
