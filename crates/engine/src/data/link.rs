// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/LinkData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson (first pass by tools/soa2rs.py, finished by hand). Not an
// OpenSWMM release.
//
// The base link store, one Vec per attribute (structure of arrays), sized and cleared
// exactly as upstream sizes and clears each one. Upstream's methods do not all visit the
// same fields, and that is kept: erase_at() leaves time_last_set, slot_volume and the
// four slot statistics unerased; reserve_to() and shrink_to_fit() each skip their own
// sets. resize() keeps pump_curve_name's prefix (upstream calls resize, not assign).
// `type` is `type_` (a Rust keyword); the static constexpr counts are associated
// consts. In erase_at() the flat-array offsets use wrapping usize arithmetic, as size_t
// does; an index that would make C++ erase an invalid range panics here. std::copy into
// a shorter destination is undefined in C++ and a panic here.

#![allow(clippy::all)]
use crate::core::options::Bytes;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i8)]
pub enum LinkType {
    Conduit = 0,
    Pump = 1,
    Orifice = 2,
    Weir = 3,
    Outlet = 4,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i16)]
pub enum XsectShape {
    Circular = 0,
    FilledCircular = 1,
    RectClosed = 2,
    RectOpen = 3,
    Trapezoidal = 4,
    Triangular = 5,
    Parabolic = 6,
    Power = 7,
    Modbaskethandle = 8,
    Eggshaped = 9,
    Horseshoe = 10,
    Gothic = 11,
    Catenary = 12,
    Semielliptical = 13,
    Baskethandle = 14,
    Semicircular = 15,
    RectTriang = 16,
    RectRound = 17,
    HorizEllipse = 18,
    VertEllipse = 19,
    Arch = 20,
    Irregular = 21,
    Custom = 22,
    ForceMain = 23,
    StreetXsect = 24,
    Dummy = 25,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i8)]
pub enum FlowClass {
    Dry = 0,
    UpDry = 1,
    DnDry = 2,
    Subcritical = 3,
    Supercritical = 4,
    UpCritical = 5,
    DnCritical = 6,
}

#[derive(Debug, Clone, Default)]
pub struct LinkData {
    pub type_: Vec<LinkType>,
    pub node1: Vec<i32>,
    pub node2: Vec<i32>,
    pub offset1: Vec<f64>,
    pub offset2: Vec<f64>,
    pub q0: Vec<f64>,
    pub q_limit: Vec<f64>,
    pub xsect_shape: Vec<XsectShape>,
    pub xsect_y_full: Vec<f64>,
    pub xsect_a_full: Vec<f64>,
    pub xsect_w_max: Vec<f64>,
    pub xsect_geom1: Vec<f64>,
    pub xsect_geom2: Vec<f64>,
    pub xsect_geom3: Vec<f64>,
    pub xsect_geom4: Vec<f64>,
    pub xsect_curve: Vec<i32>,
    pub xsect_r_full: Vec<f64>,
    pub xsect_s_full: Vec<f64>,
    pub xsect_s_max: Vec<f64>,
    pub xsect_y_bot: Vec<f64>,
    pub xsect_a_bot: Vec<f64>,
    pub xsect_s_bot: Vec<f64>,
    pub xsect_r_bot: Vec<f64>,
    pub xsect_yw_max: Vec<f64>,
    pub xsect_batch_shape: Vec<i32>,
    pub setting: Vec<f64>,
    pub target_setting: Vec<f64>,
    pub time_last_set: Vec<f64>,
    pub direction: Vec<i32>,
    pub pump_curve_name: Vec<Bytes>,
    pub has_flap_gate: Vec<u8>,
    pub dqdh: Vec<f64>,
    pub flow: Vec<f64>,
    pub kw_inflow: Vec<f64>,
    pub depth: Vec<f64>,
    pub volume: Vec<f64>,
    pub slot_volume: Vec<f64>,
    pub stat_slot_vol_dt: Vec<f64>,
    pub stat_vol_dt: Vec<f64>,
    pub stat_peak_slot_share: Vec<f64>,
    pub stat_time_slot_above: Vec<f64>,
    pub froude: Vec<f64>,
    pub flow_class: Vec<FlowClass>,
    pub is_closed: Vec<u8>,
    pub old_flow: Vec<f64>,
    pub old_depth: Vec<f64>,
    pub old_volume: Vec<f64>,
    pub conc: Vec<f64>,
    pub conc_old: Vec<f64>,
    pub conc_n_pollutants: i32,
    pub comments: Vec<Bytes>,
    pub tags: Vec<Bytes>,
    pub rpt_flag: Vec<i8>,
    pub stat_vol_flow: Vec<f64>,
    pub stat_max_flow: Vec<f64>,
    pub stat_max_veloc: Vec<f64>,
    pub stat_max_filling: Vec<f64>,
    pub stat_time_surcharged: Vec<f64>,
    pub stat_flow_class: Vec<f64>,
    pub stat_norm_ltd: Vec<f64>,
    pub stat_inlet_ctrl: Vec<f64>,
    pub stat_max_flow_date: Vec<f64>,
    pub stat_time_full_upstream: Vec<f64>,
    pub stat_time_full_dnstream: Vec<f64>,
    pub stat_time_full_both: Vec<f64>,
    pub stat_time_capacity_limited: Vec<f64>,
    pub stat_total_load: Vec<f64>,
    pub stat_n_pollutants: i32,
    pub stat_pump_cycles: Vec<i32>,
    pub stat_pump_on_time: Vec<f64>,
    pub stat_pump_volume: Vec<f64>,
    pub stat_pump_energy: Vec<f64>,
    pub stat_pump_sum_flow: Vec<f64>,
    pub stat_pump_periods: Vec<i64>,
    pub stat_pump_off_curve_low: Vec<f64>,
    pub stat_pump_off_curve_high: Vec<f64>,
    pub stat_pump_was_on: Vec<bool>,
    pub stat_flow_turns: Vec<i64>,
    pub stat_flow_turn_sign: Vec<i32>,
    pub stat_time_courant_critical: Vec<f64>,
    pub stat_lstep_min: Vec<f64>,
    pub stat_lstep_max: Vec<f64>,
    pub stat_lstep_dt_time: Vec<f64>,
    pub stat_lstep_time: Vec<f64>,
    pub stat_lstep_bin_time: Vec<f64>,
    pub stat_conv_steps: Vec<i64>,
    pub stat_conv_total: Vec<i64>,
}

// The member initializers (conc_n_pollutants = stat_n_pollutants = 0) are Rust's
// defaults, so Default is derived.

impl LinkData {
    /// Flow classes per link in `stat_flow_class` (FlowClass DRY..DN_CRITICAL).
    pub const N_FLOW_CLASSES: i32 = 7;
    /// Log2 step bins per link in `stat_lstep_bin_time`.
    pub const N_LSTEP_BINS: i32 = 5;

    pub fn lstep_bin(dt: f64, rs: f64) -> i32 {
        let mut b: i32 = 0;
        let mut edge: f64 = 0.5 * rs;
        while b < Self::N_LSTEP_BINS - 1 && dt <= edge {
            b += 1;
            edge *= 0.5;
        }
        b
    }

    pub fn add_lstep_sample(&mut self, j: i32, dt: f64, dur: f64, rs: f64) {
        let u = j as usize;
        if self.stat_lstep_time[u] <= 0.0 || dt < self.stat_lstep_min[u] {
            self.stat_lstep_min[u] = dt;
        }
        if dt > self.stat_lstep_max[u] {
            self.stat_lstep_max[u] = dt;
        }
        self.stat_lstep_dt_time[u] += dt * dur;
        self.stat_lstep_time[u] += dur;
        self.stat_lstep_bin_time[u * Self::N_LSTEP_BINS as usize + Self::lstep_bin(dt, rs) as usize] += dur;
    }

    pub fn count(&self) -> i32 {
        self.type_.len() as i32
    }

    pub fn resize(&mut self, n: i32) {
        let un = usize::try_from(n).unwrap_or(0);
        self.type_.clear();
        self.type_.resize(un, LinkType::Conduit);
        self.node1.clear();
        self.node1.resize(un, -1);
        self.node2.clear();
        self.node2.resize(un, -1);
        self.offset1.clear();
        self.offset1.resize(un, 0.0);
        self.offset2.clear();
        self.offset2.resize(un, 0.0);
        self.q0.clear();
        self.q0.resize(un, 0.0);
        self.q_limit.clear();
        self.q_limit.resize(un, 0.0);
        self.xsect_shape.clear();
        self.xsect_shape.resize(un, XsectShape::Circular);
        self.xsect_y_full.clear();
        self.xsect_y_full.resize(un, 0.0);
        self.xsect_a_full.clear();
        self.xsect_a_full.resize(un, 0.0);
        self.xsect_w_max.clear();
        self.xsect_w_max.resize(un, 0.0);
        self.xsect_geom1.clear();
        self.xsect_geom1.resize(un, 0.0);
        self.xsect_geom2.clear();
        self.xsect_geom2.resize(un, 0.0);
        self.xsect_geom3.clear();
        self.xsect_geom3.resize(un, 0.0);
        self.xsect_geom4.clear();
        self.xsect_geom4.resize(un, 0.0);
        self.xsect_curve.clear();
        self.xsect_curve.resize(un, -1);
        self.xsect_r_full.clear();
        self.xsect_r_full.resize(un, 0.0);
        self.xsect_s_full.clear();
        self.xsect_s_full.resize(un, 0.0);
        self.xsect_s_max.clear();
        self.xsect_s_max.resize(un, 0.0);
        self.xsect_y_bot.clear();
        self.xsect_y_bot.resize(un, 0.0);
        self.xsect_a_bot.clear();
        self.xsect_a_bot.resize(un, 0.0);
        self.xsect_s_bot.clear();
        self.xsect_s_bot.resize(un, 0.0);
        self.xsect_r_bot.clear();
        self.xsect_r_bot.resize(un, 0.0);
        self.xsect_yw_max.clear();
        self.xsect_yw_max.resize(un, 0.0);
        self.xsect_batch_shape.clear();
        self.xsect_batch_shape.resize(un, 0);
        self.setting.clear();
        self.setting.resize(un, 1.0);
        self.target_setting.clear();
        self.target_setting.resize(un, 1.0);
        self.time_last_set.clear();
        self.time_last_set.resize(un, 0.0);
        self.direction.clear();
        self.direction.resize(un, 1);
        self.has_flap_gate.clear();
        self.has_flap_gate.resize(un, 0);
        self.dqdh.clear();
        self.dqdh.resize(un, 0.0);
        self.pump_curve_name.resize(un, Vec::new());
        self.flow.clear();
        self.flow.resize(un, 0.0);
        self.kw_inflow.clear();
        self.kw_inflow.resize(un, 0.0);
        self.depth.clear();
        self.depth.resize(un, 0.0);
        self.volume.clear();
        self.volume.resize(un, 0.0);
        self.slot_volume.clear();
        self.slot_volume.resize(un, 0.0);
        self.stat_slot_vol_dt.clear();
        self.stat_slot_vol_dt.resize(un, 0.0);
        self.stat_vol_dt.clear();
        self.stat_vol_dt.resize(un, 0.0);
        self.stat_peak_slot_share.clear();
        self.stat_peak_slot_share.resize(un, 0.0);
        self.stat_time_slot_above.clear();
        self.stat_time_slot_above.resize(un, 0.0);
        self.froude.clear();
        self.froude.resize(un, 0.0);
        self.flow_class.clear();
        self.flow_class.resize(un, FlowClass::Dry);
        self.is_closed.clear();
        self.is_closed.resize(un, 0);
        self.old_flow.clear();
        self.old_flow.resize(un, 0.0);
        self.old_depth.clear();
        self.old_depth.resize(un, 0.0);
        self.old_volume.clear();
        self.old_volume.resize(un, 0.0);
        self.comments.clear();
        self.comments.resize(un, Vec::new());
        self.tags.clear();
        self.tags.resize(un, Vec::new());
        self.rpt_flag.clear();
        self.rpt_flag.resize(un, 0);
        self.stat_vol_flow.clear();
        self.stat_vol_flow.resize(un, 0.0);
        self.stat_max_flow.clear();
        self.stat_max_flow.resize(un, 0.0);
        self.stat_max_veloc.clear();
        self.stat_max_veloc.resize(un, 0.0);
        self.stat_max_filling.clear();
        self.stat_max_filling.resize(un, 0.0);
        self.stat_time_surcharged.clear();
        self.stat_time_surcharged.resize(un, 0.0);
        self.stat_flow_class.clear();
        self.stat_flow_class.resize(un * Self::N_FLOW_CLASSES as usize, 0.0);
        self.stat_norm_ltd.clear();
        self.stat_norm_ltd.resize(un, 0.0);
        self.stat_inlet_ctrl.clear();
        self.stat_inlet_ctrl.resize(un, 0.0);
        self.stat_max_flow_date.clear();
        self.stat_max_flow_date.resize(un, 0.0);
        self.stat_time_full_upstream.clear();
        self.stat_time_full_upstream.resize(un, 0.0);
        self.stat_time_full_dnstream.clear();
        self.stat_time_full_dnstream.resize(un, 0.0);
        self.stat_time_full_both.clear();
        self.stat_time_full_both.resize(un, 0.0);
        self.stat_time_capacity_limited.clear();
        self.stat_time_capacity_limited.resize(un, 0.0);
        self.stat_pump_cycles.clear();
        self.stat_pump_cycles.resize(un, 0);
        self.stat_pump_on_time.clear();
        self.stat_pump_on_time.resize(un, 0.0);
        self.stat_pump_volume.clear();
        self.stat_pump_volume.resize(un, 0.0);
        self.stat_pump_energy.clear();
        self.stat_pump_energy.resize(un, 0.0);
        self.stat_pump_off_curve_high.clear();
        self.stat_pump_off_curve_high.resize(un, 0.0);
        self.stat_pump_off_curve_low.clear();
        self.stat_pump_off_curve_low.resize(un, 0.0);
        self.stat_pump_periods.clear();
        self.stat_pump_periods.resize(un, 0); // assign(un, 0.0) on vector<long>
        self.stat_pump_sum_flow.clear();
        self.stat_pump_sum_flow.resize(un, 0.0);
        self.stat_pump_was_on.clear();
        self.stat_pump_was_on.resize(un, false);
        self.stat_flow_turns.clear();
        self.stat_flow_turns.resize(un, 0);
        self.stat_flow_turn_sign.clear();
        self.stat_flow_turn_sign.resize(un, 0);
        self.stat_time_courant_critical.clear();
        self.stat_time_courant_critical.resize(un, 0.0);
        self.stat_lstep_min.clear();
        self.stat_lstep_min.resize(un, 0.0);
        self.stat_lstep_max.clear();
        self.stat_lstep_max.resize(un, 0.0);
        self.stat_lstep_dt_time.clear();
        self.stat_lstep_dt_time.resize(un, 0.0);
        self.stat_lstep_time.clear();
        self.stat_lstep_time.resize(un, 0.0);
        self.stat_lstep_bin_time.clear();
        self.stat_lstep_bin_time.resize(un * Self::N_LSTEP_BINS as usize, 0.0);
        self.stat_conv_steps.clear();
        self.stat_conv_steps.resize(un, 0);
        self.stat_conv_total.clear();
        self.stat_conv_total.resize(un, 0);
    }

    pub fn grow_to(&mut self, n: i32) {
        if n <= self.count() {
            return;
        }
        let un = usize::try_from(n).unwrap_or(0);
        self.type_.resize(un, LinkType::Conduit);
        self.node1.resize(un, -1);
        self.node2.resize(un, -1);
        self.offset1.resize(un, 0.0);
        self.offset2.resize(un, 0.0);
        self.q0.resize(un, 0.0);
        self.q_limit.resize(un, 0.0);
        self.xsect_shape.resize(un, XsectShape::Circular);
        self.xsect_y_full.resize(un, 0.0);
        self.xsect_a_full.resize(un, 0.0);
        self.xsect_w_max.resize(un, 0.0);
        self.xsect_geom1.resize(un, 0.0);
        self.xsect_geom2.resize(un, 0.0);
        self.xsect_geom3.resize(un, 0.0);
        self.xsect_geom4.resize(un, 0.0);
        self.xsect_curve.resize(un, -1);
        self.xsect_r_full.resize(un, 0.0);
        self.xsect_s_full.resize(un, 0.0);
        self.xsect_s_max.resize(un, 0.0);
        self.xsect_y_bot.resize(un, 0.0);
        self.xsect_a_bot.resize(un, 0.0);
        self.xsect_s_bot.resize(un, 0.0);
        self.xsect_r_bot.resize(un, 0.0);
        self.xsect_yw_max.resize(un, 0.0);
        self.xsect_batch_shape.resize(un, 0);
        self.setting.resize(un, 1.0);
        self.target_setting.resize(un, 1.0);
        self.time_last_set.resize(un, 0.0);
        self.direction.resize(un, 1);
        self.has_flap_gate.resize(un, 0);
        self.dqdh.resize(un, 0.0);
        self.pump_curve_name.resize(un, Vec::new());
        self.flow.resize(un, 0.0);
        self.depth.resize(un, 0.0);
        self.volume.resize(un, 0.0);
        self.slot_volume.resize(un, 0.0);
        self.kw_inflow.resize(un, 0.0);
        self.stat_slot_vol_dt.resize(un, 0.0);
        self.stat_vol_dt.resize(un, 0.0);
        self.stat_peak_slot_share.resize(un, 0.0);
        self.stat_time_slot_above.resize(un, 0.0);
        self.froude.resize(un, 0.0);
        self.flow_class.resize(un, FlowClass::Dry);
        self.is_closed.resize(un, 0);
        self.old_flow.resize(un, 0.0);
        self.old_depth.resize(un, 0.0);
        self.old_volume.resize(un, 0.0);
        self.comments.resize(un, Vec::new());
        self.tags.resize(un, Vec::new());
        self.rpt_flag.resize(un, 0);
        self.stat_vol_flow.resize(un, 0.0);
        self.stat_max_flow.resize(un, 0.0);
        self.stat_max_veloc.resize(un, 0.0);
        self.stat_max_filling.resize(un, 0.0);
        self.stat_time_surcharged.resize(un, 0.0);
        self.stat_max_flow_date.resize(un, 0.0);
        self.stat_time_full_upstream.resize(un, 0.0);
        self.stat_time_full_dnstream.resize(un, 0.0);
        self.stat_time_full_both.resize(un, 0.0);
        self.stat_time_capacity_limited.resize(un, 0.0);
        self.stat_pump_cycles.resize(un, 0);
        self.stat_pump_on_time.resize(un, 0.0);
        self.stat_pump_volume.resize(un, 0.0);
        self.stat_pump_energy.resize(un, 0.0);
        self.stat_pump_off_curve_high.resize(un, 0.0);
        self.stat_pump_off_curve_low.resize(un, 0.0);
        self.stat_pump_periods.resize(un, 0); // assign(un, 0.0) on vector<long>
        self.stat_pump_sum_flow.resize(un, 0.0);
        self.stat_pump_was_on.resize(un, false);
        self.stat_flow_turns.resize(un, 0);
        self.stat_flow_turn_sign.resize(un, 0);
        self.stat_time_courant_critical.resize(un, 0.0);
        self.stat_lstep_min.resize(un, 0.0);
        self.stat_lstep_max.resize(un, 0.0);
        self.stat_lstep_dt_time.resize(un, 0.0);
        self.stat_lstep_time.resize(un, 0.0);
        self.stat_conv_steps.resize(un, 0);
        self.stat_conv_total.resize(un, 0);
        self.stat_lstep_bin_time.resize(un * Self::N_LSTEP_BINS as usize, 0.0);
        self.stat_norm_ltd.resize(un, 0.0);
        self.stat_inlet_ctrl.resize(un, 0.0);
        // stat_flow_class is flat 2D [n * N_FLOW_CLASSES]
        self.stat_flow_class.resize(un * Self::N_FLOW_CLASSES as usize, 0.0);
        // Note: conc, conc_old handled by resize_quality()
        // Note: stat_total_load handled by resize_loads()
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
        self.node1.reserve(un.saturating_sub(self.node1.len()));
        self.node2.reserve(un.saturating_sub(self.node2.len()));
        self.offset1.reserve(un.saturating_sub(self.offset1.len()));
        self.offset2.reserve(un.saturating_sub(self.offset2.len()));
        self.q0.reserve(un.saturating_sub(self.q0.len()));
        self.q_limit.reserve(un.saturating_sub(self.q_limit.len()));
        self.xsect_shape.reserve(un.saturating_sub(self.xsect_shape.len()));
        self.xsect_y_full.reserve(un.saturating_sub(self.xsect_y_full.len()));
        self.xsect_a_full.reserve(un.saturating_sub(self.xsect_a_full.len()));
        self.xsect_w_max.reserve(un.saturating_sub(self.xsect_w_max.len()));
        self.xsect_geom1.reserve(un.saturating_sub(self.xsect_geom1.len()));
        self.xsect_geom2.reserve(un.saturating_sub(self.xsect_geom2.len()));
        self.xsect_geom3.reserve(un.saturating_sub(self.xsect_geom3.len()));
        self.xsect_geom4.reserve(un.saturating_sub(self.xsect_geom4.len()));
        self.xsect_curve.reserve(un.saturating_sub(self.xsect_curve.len()));
        self.xsect_r_full.reserve(un.saturating_sub(self.xsect_r_full.len()));
        self.xsect_s_full.reserve(un.saturating_sub(self.xsect_s_full.len()));
        self.xsect_s_max.reserve(un.saturating_sub(self.xsect_s_max.len()));
        self.xsect_y_bot.reserve(un.saturating_sub(self.xsect_y_bot.len()));
        self.xsect_a_bot.reserve(un.saturating_sub(self.xsect_a_bot.len()));
        self.xsect_s_bot.reserve(un.saturating_sub(self.xsect_s_bot.len()));
        self.xsect_r_bot.reserve(un.saturating_sub(self.xsect_r_bot.len()));
        self.xsect_yw_max.reserve(un.saturating_sub(self.xsect_yw_max.len()));
        self.xsect_batch_shape.reserve(un.saturating_sub(self.xsect_batch_shape.len()));
        self.setting.reserve(un.saturating_sub(self.setting.len()));
        self.target_setting.reserve(un.saturating_sub(self.target_setting.len()));
        self.time_last_set.reserve(un.saturating_sub(self.time_last_set.len()));
        self.direction.reserve(un.saturating_sub(self.direction.len()));
        self.has_flap_gate.reserve(un.saturating_sub(self.has_flap_gate.len()));
        self.dqdh.reserve(un.saturating_sub(self.dqdh.len()));
        self.flow.reserve(un.saturating_sub(self.flow.len()));
        self.kw_inflow.reserve(un.saturating_sub(self.kw_inflow.len()));
        self.depth.reserve(un.saturating_sub(self.depth.len()));
        self.volume.reserve(un.saturating_sub(self.volume.len()));
        self.froude.reserve(un.saturating_sub(self.froude.len()));
        self.flow_class.reserve(un.saturating_sub(self.flow_class.len()));
        self.is_closed.reserve(un.saturating_sub(self.is_closed.len()));
        self.old_flow.reserve(un.saturating_sub(self.old_flow.len()));
        self.old_depth.reserve(un.saturating_sub(self.old_depth.len()));
        self.old_volume.reserve(un.saturating_sub(self.old_volume.len()));
        self.rpt_flag.reserve(un.saturating_sub(self.rpt_flag.len()));
        self.stat_vol_flow.reserve(un.saturating_sub(self.stat_vol_flow.len()));
        self.stat_max_flow.reserve(un.saturating_sub(self.stat_max_flow.len()));
        self.stat_max_veloc.reserve(un.saturating_sub(self.stat_max_veloc.len()));
        self.stat_max_filling.reserve(un.saturating_sub(self.stat_max_filling.len()));
        self.stat_time_surcharged.reserve(un.saturating_sub(self.stat_time_surcharged.len()));
        self.stat_max_flow_date.reserve(un.saturating_sub(self.stat_max_flow_date.len()));
        self.stat_time_full_upstream.reserve(un.saturating_sub(self.stat_time_full_upstream.len()));
        self.stat_time_full_dnstream.reserve(un.saturating_sub(self.stat_time_full_dnstream.len()));
        self.stat_time_full_both.reserve(un.saturating_sub(self.stat_time_full_both.len()));
        self.stat_time_capacity_limited.reserve(un.saturating_sub(self.stat_time_capacity_limited.len()));
        self.stat_pump_cycles.reserve(un.saturating_sub(self.stat_pump_cycles.len()));
        self.stat_pump_on_time.reserve(un.saturating_sub(self.stat_pump_on_time.len()));
        self.stat_pump_volume.reserve(un.saturating_sub(self.stat_pump_volume.len()));
        self.stat_pump_energy.reserve(un.saturating_sub(self.stat_pump_energy.len()));
        self.stat_pump_off_curve_high.reserve(un.saturating_sub(self.stat_pump_off_curve_high.len()));
        self.stat_pump_off_curve_low.reserve(un.saturating_sub(self.stat_pump_off_curve_low.len()));
        self.stat_pump_periods.reserve(un.saturating_sub(self.stat_pump_periods.len()));
        self.stat_pump_sum_flow.reserve(un.saturating_sub(self.stat_pump_sum_flow.len()));
        self.stat_pump_was_on.reserve(un.saturating_sub(self.stat_pump_was_on.len()));
        self.stat_flow_turns.reserve(un.saturating_sub(self.stat_flow_turns.len()));
        self.stat_flow_turn_sign.reserve(un.saturating_sub(self.stat_flow_turn_sign.len()));
        self.stat_time_courant_critical.reserve(un.saturating_sub(self.stat_time_courant_critical.len()));
        self.stat_norm_ltd.reserve(un.saturating_sub(self.stat_norm_ltd.len()));
        self.stat_inlet_ctrl.reserve(un.saturating_sub(self.stat_inlet_ctrl.len()));
        self.pump_curve_name.reserve(un.saturating_sub(self.pump_curve_name.len()));
        self.comments.reserve(un.saturating_sub(self.comments.len()));
        self.tags.reserve(un.saturating_sub(self.tags.len()));
        self.stat_flow_class.reserve(un.saturating_sub(self.stat_flow_class.len()));
        self.stat_lstep_min.reserve(un.saturating_sub(self.stat_lstep_min.len()));
        self.stat_lstep_max.reserve(un.saturating_sub(self.stat_lstep_max.len()));
        self.stat_lstep_dt_time.reserve(un.saturating_sub(self.stat_lstep_dt_time.len()));
        self.stat_lstep_time.reserve(un.saturating_sub(self.stat_lstep_time.len()));
        self.stat_conv_steps.reserve(un.saturating_sub(self.stat_conv_steps.len()));
        self.stat_conv_total.reserve(un.saturating_sub(self.stat_conv_total.len()));
        {
            let t = un * Self::N_LSTEP_BINS as usize;
            self.stat_lstep_bin_time.reserve(t.saturating_sub(self.stat_lstep_bin_time.len()));
        }
    }

    pub fn erase_at(&mut self, idx: i32) {
        let ui = idx as usize;
        if ui < self.type_.len() { self.type_.remove(ui); }
        if ui < self.node1.len() { self.node1.remove(ui); }
        if ui < self.node2.len() { self.node2.remove(ui); }
        if ui < self.offset1.len() { self.offset1.remove(ui); }
        if ui < self.offset2.len() { self.offset2.remove(ui); }
        if ui < self.q0.len() { self.q0.remove(ui); }
        if ui < self.q_limit.len() { self.q_limit.remove(ui); }
        if ui < self.xsect_shape.len() { self.xsect_shape.remove(ui); }
        if ui < self.xsect_y_full.len() { self.xsect_y_full.remove(ui); }
        if ui < self.xsect_a_full.len() { self.xsect_a_full.remove(ui); }
        if ui < self.xsect_w_max.len() { self.xsect_w_max.remove(ui); }
        if ui < self.xsect_curve.len() { self.xsect_curve.remove(ui); }
        if ui < self.xsect_geom1.len() { self.xsect_geom1.remove(ui); }
        if ui < self.xsect_geom2.len() { self.xsect_geom2.remove(ui); }
        if ui < self.xsect_geom3.len() { self.xsect_geom3.remove(ui); }
        if ui < self.xsect_geom4.len() { self.xsect_geom4.remove(ui); }
        if ui < self.xsect_r_full.len() { self.xsect_r_full.remove(ui); }
        if ui < self.xsect_s_full.len() { self.xsect_s_full.remove(ui); }
        if ui < self.xsect_s_max.len() { self.xsect_s_max.remove(ui); }
        if ui < self.xsect_y_bot.len() { self.xsect_y_bot.remove(ui); }
        if ui < self.xsect_a_bot.len() { self.xsect_a_bot.remove(ui); }
        if ui < self.xsect_s_bot.len() { self.xsect_s_bot.remove(ui); }
        if ui < self.xsect_r_bot.len() { self.xsect_r_bot.remove(ui); }
        if ui < self.xsect_yw_max.len() { self.xsect_yw_max.remove(ui); }
        if ui < self.xsect_batch_shape.len() { self.xsect_batch_shape.remove(ui); }
        if ui < self.setting.len() { self.setting.remove(ui); }
        if ui < self.target_setting.len() { self.target_setting.remove(ui); }
        if ui < self.direction.len() { self.direction.remove(ui); }
        if ui < self.pump_curve_name.len() { self.pump_curve_name.remove(ui); }
        if ui < self.has_flap_gate.len() { self.has_flap_gate.remove(ui); }
        if ui < self.dqdh.len() { self.dqdh.remove(ui); }
        if ui < self.flow.len() { self.flow.remove(ui); }
        if ui < self.depth.len() { self.depth.remove(ui); }
        if ui < self.volume.len() { self.volume.remove(ui); }
        if ui < self.froude.len() { self.froude.remove(ui); }
        if ui < self.flow_class.len() { self.flow_class.remove(ui); }
        if ui < self.is_closed.len() { self.is_closed.remove(ui); }
        if ui < self.kw_inflow.len() { self.kw_inflow.remove(ui); }
        if ui < self.old_flow.len() { self.old_flow.remove(ui); }
        if ui < self.old_depth.len() { self.old_depth.remove(ui); }
        if ui < self.old_volume.len() { self.old_volume.remove(ui); }
        if ui < self.comments.len() { self.comments.remove(ui); }
        if ui < self.tags.len() { self.tags.remove(ui); }
        if ui < self.rpt_flag.len() { self.rpt_flag.remove(ui); }
        if ui < self.stat_vol_flow.len() { self.stat_vol_flow.remove(ui); }
        if ui < self.stat_max_flow.len() { self.stat_max_flow.remove(ui); }
        if ui < self.stat_max_veloc.len() { self.stat_max_veloc.remove(ui); }
        if ui < self.stat_max_filling.len() { self.stat_max_filling.remove(ui); }
        if ui < self.stat_time_surcharged.len() { self.stat_time_surcharged.remove(ui); }
        if ui < self.stat_norm_ltd.len() { self.stat_norm_ltd.remove(ui); }
        if ui < self.stat_inlet_ctrl.len() { self.stat_inlet_ctrl.remove(ui); }
        if ui < self.stat_max_flow_date.len() { self.stat_max_flow_date.remove(ui); }
        if ui < self.stat_time_full_upstream.len() { self.stat_time_full_upstream.remove(ui); }
        if ui < self.stat_time_full_dnstream.len() { self.stat_time_full_dnstream.remove(ui); }
        if ui < self.stat_time_full_both.len() { self.stat_time_full_both.remove(ui); }
        if ui < self.stat_time_capacity_limited.len() { self.stat_time_capacity_limited.remove(ui); }
        if ui < self.stat_pump_cycles.len() { self.stat_pump_cycles.remove(ui); }
        if ui < self.stat_pump_on_time.len() { self.stat_pump_on_time.remove(ui); }
        if ui < self.stat_pump_volume.len() { self.stat_pump_volume.remove(ui); }
        if ui < self.stat_pump_energy.len() { self.stat_pump_energy.remove(ui); }
        if ui < self.stat_pump_off_curve_high.len() { self.stat_pump_off_curve_high.remove(ui); }
        if ui < self.stat_pump_off_curve_low.len() { self.stat_pump_off_curve_low.remove(ui); }
        if ui < self.stat_pump_periods.len() { self.stat_pump_periods.remove(ui); }
        if ui < self.stat_pump_sum_flow.len() { self.stat_pump_sum_flow.remove(ui); }
        if ui < self.stat_pump_was_on.len() { self.stat_pump_was_on.remove(ui); }
        if ui < self.stat_flow_turns.len() { self.stat_flow_turns.remove(ui); }
        if ui < self.stat_flow_turn_sign.len() { self.stat_flow_turn_sign.remove(ui); }
        if ui < self.stat_time_courant_critical.len() { self.stat_time_courant_critical.remove(ui); }
        if ui < self.stat_lstep_min.len() { self.stat_lstep_min.remove(ui); }
        if ui < self.stat_lstep_max.len() { self.stat_lstep_max.remove(ui); }
        if ui < self.stat_lstep_dt_time.len() { self.stat_lstep_dt_time.remove(ui); }
        if ui < self.stat_lstep_time.len() { self.stat_lstep_time.remove(ui); }
        if ui < self.stat_conv_steps.len() { self.stat_conv_steps.remove(ui); }
        if ui < self.stat_conv_total.len() { self.stat_conv_total.remove(ui); }
        {
            let base = ui.wrapping_mul(Self::N_LSTEP_BINS as usize);
            let end = base.wrapping_add(Self::N_LSTEP_BINS as usize);
            if end <= self.stat_lstep_bin_time.len() {
                self.stat_lstep_bin_time.drain(base..end);
            }
        }

        // Flat 2D: stat_flow_class [link * N_FLOW_CLASSES + class]
        {
            let base = ui.wrapping_mul(Self::N_FLOW_CLASSES as usize);
            let end = base.wrapping_add(Self::N_FLOW_CLASSES as usize);
            if end <= self.stat_flow_class.len() {
                self.stat_flow_class.drain(base..end);
            }
        }

        // Flat 2D quality arrays: [link * np + p]
        if self.conc_n_pollutants > 0 {
            let np = self.conc_n_pollutants as usize;
            let base = ui.wrapping_mul(np);
            for v in [&mut self.conc, &mut self.conc_old] {
                if base.wrapping_add(np) <= v.len() {
                    v.drain(base..base.wrapping_add(np));
                }
            }
        }

        // Flat 2D stat load: [link * np + p]
        if self.stat_n_pollutants > 0 {
            let np = self.stat_n_pollutants as usize;
            let base = ui.wrapping_mul(np);
            if base.wrapping_add(np) <= self.stat_total_load.len() {
                self.stat_total_load.drain(base..base.wrapping_add(np));
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
        }
    }

    pub fn shrink_to_fit(&mut self) {
        self.type_.shrink_to_fit();
        self.node1.shrink_to_fit();
        self.node2.shrink_to_fit();
        self.offset1.shrink_to_fit();
        self.offset2.shrink_to_fit();
        self.q0.shrink_to_fit();
        self.q_limit.shrink_to_fit();
        self.xsect_shape.shrink_to_fit();
        self.xsect_y_full.shrink_to_fit();
        self.xsect_a_full.shrink_to_fit();
        self.xsect_w_max.shrink_to_fit();
        self.xsect_curve.shrink_to_fit();
        self.xsect_r_full.shrink_to_fit();
        self.xsect_s_full.shrink_to_fit();
        self.xsect_s_max.shrink_to_fit();
        self.xsect_y_bot.shrink_to_fit();
        self.xsect_a_bot.shrink_to_fit();
        self.xsect_s_bot.shrink_to_fit();
        self.xsect_r_bot.shrink_to_fit();
        self.xsect_yw_max.shrink_to_fit();
        self.xsect_batch_shape.shrink_to_fit();
        self.setting.shrink_to_fit();
        self.target_setting.shrink_to_fit();
        self.direction.shrink_to_fit();
        self.has_flap_gate.shrink_to_fit();
        self.dqdh.shrink_to_fit();
        self.pump_curve_name.shrink_to_fit();
        self.flow.shrink_to_fit();
        self.kw_inflow.shrink_to_fit();
        self.depth.shrink_to_fit();
        self.volume.shrink_to_fit();
        self.froude.shrink_to_fit();
        self.flow_class.shrink_to_fit();
        self.is_closed.shrink_to_fit();
        self.old_flow.shrink_to_fit();
        self.old_depth.shrink_to_fit();
        self.old_volume.shrink_to_fit();
        self.conc.shrink_to_fit();
        self.conc_old.shrink_to_fit();
        self.comments.shrink_to_fit();
        self.tags.shrink_to_fit();
        self.rpt_flag.shrink_to_fit();
        self.stat_vol_flow.shrink_to_fit();
        self.stat_max_flow.shrink_to_fit();
        self.stat_max_veloc.shrink_to_fit();
        self.stat_max_filling.shrink_to_fit();
        self.stat_time_surcharged.shrink_to_fit();
        self.stat_flow_class.shrink_to_fit();
        self.stat_norm_ltd.shrink_to_fit();
        self.stat_inlet_ctrl.shrink_to_fit();
        self.stat_max_flow_date.shrink_to_fit();
        self.stat_time_full_upstream.shrink_to_fit();
        self.stat_time_full_dnstream.shrink_to_fit();
        self.stat_time_full_both.shrink_to_fit();
        self.stat_time_capacity_limited.shrink_to_fit();
        self.stat_pump_cycles.shrink_to_fit();
        self.stat_pump_on_time.shrink_to_fit();
        self.stat_pump_volume.shrink_to_fit();
        self.stat_pump_energy.shrink_to_fit();
        self.stat_pump_off_curve_high.shrink_to_fit();
        self.stat_pump_off_curve_low.shrink_to_fit();
        self.stat_pump_periods.shrink_to_fit();
        self.stat_pump_sum_flow.shrink_to_fit();
        self.stat_pump_was_on.shrink_to_fit();
        self.stat_total_load.shrink_to_fit();
    }

    pub fn save_state(&mut self) {
        self.save_hyd_state();
        self.save_qual_state();
    }

    pub fn save_hyd_state(&mut self) {
        { let n = self.flow.len(); self.old_flow[..n].copy_from_slice(&self.flow); }
        { let n = self.depth.len(); self.old_depth[..n].copy_from_slice(&self.depth); }
        { let n = self.volume.len(); self.old_volume[..n].copy_from_slice(&self.volume); }
    }

    pub fn save_qual_state(&mut self) {
        { let n = self.conc.len(); self.conc_old[..n].copy_from_slice(&self.conc); }
    }

    pub fn reset_state(&mut self) {
        self.flow.fill(0.0);
        self.kw_inflow.fill(0.0);
        self.depth.fill(0.0);
        self.volume.fill(0.0);
        self.froude.fill(0.0);
        self.flow_class.fill(FlowClass::Dry);
        self.old_flow.fill(0.0);
        self.old_depth.fill(0.0);
        self.old_volume.fill(0.0);
        self.conc.fill(0.0);
        self.conc_old.fill(0.0);
    }
}
