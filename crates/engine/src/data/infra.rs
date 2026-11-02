// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/InfraData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson (first pass by tools/soa2rs.py, finished by hand). Not an
// OpenSWMM release.
//
// Authoring stores for [TRANSECTS], [STREETS], [INLETS], [INLET_USAGE] and [CONTROLS],
// one Vec per column in upstream's field order (InletUsageStore's statistic columns are
// declared after its methods upstream and so come last here too). Strings are bytes. A
// parameter named `type` is `type_`. InletStore::erase_row is unchecked upstream; a bad
// index is undefined there and a panic here. InletUsageStore::erase_row converts the
// index to size_t, so a negative one wraps and skips every column, as `idx as usize`
// does here. resize_stats treats a negative count as zero, as data/node.rs does.

#![allow(clippy::all)]
use crate::core::options::Bytes;

// ============================================================================
// Transect definitions (from [TRANSECTS] section)
// ============================================================================

#[derive(Debug, Clone, Default)]
pub struct TransectStore {
    pub names: Vec<Bytes>,
    /// Free-form description per transect (DA-ENG-09)
    pub comments: Vec<Bytes>,
    pub n_left: Vec<f64>,
    pub n_right: Vec<f64>,
    pub n_channel: Vec<f64>,
    pub x_left_bank: Vec<f64>,
    pub x_right_bank: Vec<f64>,
    /// Encroachment stations — independent of bank stations (BQ-TR-02).
    /// Default to 0.0 at add-time.
    pub x_left_encroachment: Vec<f64>,
    pub x_right_encroachment: Vec<f64>,
    pub x_factor: Vec<f64>,
    pub y_factor: Vec<f64>,
    /// Meander factor (channel/floodplain length ratio); default 1.0
    pub length_factor: Vec<f64>,
    /// Station-elevation pairs per transect
    pub stations: Vec<Vec<f64>>,
    pub elevations: Vec<Vec<f64>>,
}

impl TransectStore {
    pub fn count(&self) -> i32 {
        self.names.len() as i32
    }
}

// ============================================================================
// Street definitions (from [STREETS] section)
// ============================================================================

#[derive(Debug, Clone, Default)]
pub struct StreetStore {
    pub names: Vec<Bytes>,
    pub t_crown: Vec<f64>,
    pub h_curb: Vec<f64>,
    /// Cross slope (%)
    pub sx: Vec<f64>,
    pub n_road: Vec<f64>,
    pub gutter_depres: Vec<f64>,
    pub gutter_width: Vec<f64>,
    pub sides: Vec<i32>,
    pub back_width: Vec<f64>,
    pub back_slope: Vec<f64>,
    pub back_n: Vec<f64>,
}

impl StreetStore {
    pub fn count(&self) -> i32 {
        self.names.len() as i32
    }
}

// ============================================================================
// Inlet definitions (from [INLETS] section)
// ============================================================================

/// Inlet design store. All lengths/velocities are kept in USER (display) units exactly
/// as read from the .inp; the runtime InletSolver converts to internal units at init.
///
/// Column usage by type:
///   GRATE / DROP_GRATE : length, width, grate_type, open_area, splash_veloc
///   SLOTTED            : length, width
///   CURB / DROP_CURB   : curb_length, curb_height, curb_throat
///   COMBO              : grate columns + curb columns
///   CUSTOM             : curve_id, curve_index, curve_kind
#[derive(Debug, Clone, Default)]
pub struct InletStore {
    pub names: Vec<Bytes>,
    /// GRATE/CURB/COMBO/SLOTTED/DROP_GRATE/DROP_CURB/CUSTOM
    pub inlet_type: Vec<Bytes>,
    /// grate or slotted length
    pub length: Vec<f64>,
    /// grate or slotted width
    pub width: Vec<f64>,
    /// P_BAR-50, P_BAR-50x100, P_BAR-30, CURVED_VANE, TILT_BAR-45, TILT_BAR-30,
    /// RETICULINE, GENERIC
    pub grate_type: Vec<Bytes>,
    /// GENERIC only: open-area fraction (0,1]
    pub open_area: Vec<f64>,
    /// GENERIC only: splash-over velocity
    pub splash_veloc: Vec<f64>,
    /// CURB/DROP_CURB/COMBO
    pub curb_length: Vec<f64>,
    /// CURB/DROP_CURB/COMBO
    pub curb_height: Vec<f64>,
    /// 0=HORIZONTAL 1=INCLINED 2=VERTICAL (legacy ThroatAngleWords order); default 2
    pub curb_throat: Vec<i32>,
    /// CUSTOM: capture curve name
    pub curve_id: Vec<Bytes>,
    /// CUSTOM: resolved curve index (-1 until PostParseResolver)
    pub curve_index: Vec<i32>,
    /// CUSTOM: 0=unresolved, 1=DIVERSION (captured vs approach flow), 2=RATING
    /// (captured vs depth)
    pub curve_kind: Vec<i32>,
    /// optional description (round-tripped as ';' comment)
    pub comments: Vec<Bytes>,
}

impl InletStore {
    pub fn count(&self) -> i32 {
        self.names.len() as i32
    }

    /// Append a copied design name/type and default parameter columns; returns the new
    /// row index.
    pub fn add_row(&mut self, name: &[u8], type_: &[u8]) -> i32 {
        self.names.push(name.to_vec());
        self.inlet_type.push(type_.to_vec());
        self.length.push(0.0);
        self.width.push(0.0);
        self.grate_type.push(Bytes::new());
        self.open_area.push(0.0);
        self.splash_veloc.push(0.0);
        self.curb_length.push(0.0);
        self.curb_height.push(0.0);
        self.curb_throat.push(2);
        self.curve_id.push(Bytes::new());
        self.curve_index.push(-1);
        self.curve_kind.push(0);
        self.comments.push(Bytes::new());
        self.count() - 1
    }

    /// Erase one inlet design from every owned column (unchecked index).
    pub fn erase_row(&mut self, idx: i32) {
        let u = idx as usize;
        self.names.remove(u);
        self.inlet_type.remove(u);
        self.length.remove(u);
        self.width.remove(u);
        self.grate_type.remove(u);
        self.open_area.remove(u);
        self.splash_veloc.remove(u);
        self.curb_length.remove(u);
        self.curb_height.remove(u);
        self.curb_throat.remove(u);
        self.curve_id.remove(u);
        self.curve_index.remove(u);
        self.curve_kind.remove(u);
        self.comments.remove(u);
    }
}

// ============================================================================
// Inlet usage (from [INLET_USAGE] section)
// ============================================================================

/// One row per inlet placement. Two host kinds share the store:
///   * conduit-attribute usage ([INLET_USAGE]): link_index >= 0, node_host == -1
///   * inlet junction ([INLET_JUNCTIONS]):      link_index == -1, node_host >= 0
/// flow_limit / local_depress / local_width are kept in USER units (as read).
#[derive(Debug, Clone, Default)]
pub struct InletUsageStore {
    /// Host conduit link index (-1 for an inlet junction)
    pub link_index: Vec<i32>,
    /// Host inlet-junction node index (-1 for a conduit usage)
    pub node_host: Vec<i32>,
    /// Index into InletStore
    pub design_index: Vec<i32>,
    /// Capture (receiving / underdrain) node index
    pub node_index: Vec<i32>,
    /// Number of inlets per side
    pub num_inlets: Vec<i32>,
    /// 0=auto, 1=on_grade, 2=on_sag
    pub placement: Vec<i32>,
    /// 1.0 - pctClogged/100
    pub clog_factor: Vec<f64>,
    /// Max capture flow per inlet (user flow units), 0=unlimited
    pub flow_limit: Vec<f64>,
    /// Local gutter depression (user length units)
    pub local_depress: Vec<f64>,
    /// Local depression width (user length units)
    pub local_width: Vec<f64>,
    /// Index into StreetStore, resolved by PostParseResolver (-1 if host xsect is not
    /// STREET)
    pub street_index: Vec<i32>,
    /// Design / capture-node names of an [INLET_JUNCTIONS] row, kept until
    /// PostParseResolver turns them into design_index / node_index.
    pub pending_design_name: Vec<Bytes>,
    pub pending_capture_name: Vec<Bytes>,

    // Stats (populated by InletSolver::gatherStats() before reporting)
    /// Total captured volume (ft³)
    pub stat_capture_vol: Vec<f64>,
    /// Total bypassed volume (ft³)
    pub stat_bypass_vol: Vec<f64>,
    /// Total backflow volume (ft³)
    pub stat_backflow_vol: Vec<f64>,
    /// Peak captured flow (cfs)
    pub stat_peak_flow: Vec<f64>,
}

/// erase_row's lambda `e`: erase element `u` if the column holds it.
fn erase_if_present<T>(v: &mut Vec<T>, u: usize) {
    if u < v.len() {
        v.remove(u);
    }
}

impl InletUsageStore {
    pub fn count(&self) -> i32 {
        self.design_index.len() as i32
    }

    /// Append an unresolved/default inlet usage with explicit host references; returns
    /// the new row index. Does not grow the statistic columns.
    pub fn add_row(
        &mut self,
        link_idx: i32,
        node_host_idx: i32,
        design_idx: i32,
        capture_node_idx: i32,
    ) -> i32 {
        self.link_index.push(link_idx);
        self.node_host.push(node_host_idx);
        self.design_index.push(design_idx);
        self.node_index.push(capture_node_idx);
        self.num_inlets.push(1);
        self.placement.push(0);
        self.clog_factor.push(1.0);
        self.flow_limit.push(0.0);
        self.local_depress.push(0.0);
        self.local_width.push(0.0);
        self.street_index.push(-1);
        self.pending_design_name.push(Bytes::new());
        self.pending_capture_name.push(Bytes::new());
        self.count() - 1
    }

    /// Erase a usage row from each column that contains it, statistic columns included.
    pub fn erase_row(&mut self, idx: i32) {
        let u = idx as usize;
        erase_if_present(&mut self.link_index, u);
        erase_if_present(&mut self.node_host, u);
        erase_if_present(&mut self.design_index, u);
        erase_if_present(&mut self.node_index, u);
        erase_if_present(&mut self.num_inlets, u);
        erase_if_present(&mut self.placement, u);
        erase_if_present(&mut self.clog_factor, u);
        erase_if_present(&mut self.flow_limit, u);
        erase_if_present(&mut self.local_depress, u);
        erase_if_present(&mut self.local_width, u);
        erase_if_present(&mut self.street_index, u);
        erase_if_present(&mut self.pending_design_name, u);
        erase_if_present(&mut self.pending_capture_name, u);
        erase_if_present(&mut self.stat_capture_vol, u);
        erase_if_present(&mut self.stat_bypass_vol, u);
        erase_if_present(&mut self.stat_backflow_vol, u);
        erase_if_present(&mut self.stat_peak_flow, u);
    }

    /// First usage whose conduit-host key equals `link_idx`, or -1.
    pub fn find_by_link(&self, link_idx: i32) -> i32 {
        for i in 0..self.count() {
            if self.link_index[i as usize] == link_idx {
                return i;
            }
        }
        -1
    }

    /// First usage whose inlet-junction host key equals `node_idx`, or -1.
    pub fn find_by_node_host(&self, node_idx: i32) -> i32 {
        for i in 0..self.count() {
            if self.node_host[i as usize] == node_idx {
                return i;
            }
        }
        -1
    }

    /// Replace the statistic columns with `n` zeros each.
    pub fn resize_stats(&mut self, n: i32) {
        let un = usize::try_from(n).unwrap_or(0);
        self.stat_capture_vol.clear();
        self.stat_capture_vol.resize(un, 0.0);
        self.stat_bypass_vol.clear();
        self.stat_bypass_vol.resize(un, 0.0);
        self.stat_backflow_vol.clear();
        self.stat_backflow_vol.resize(un, 0.0);
        self.stat_peak_flow.clear();
        self.stat_peak_flow.resize(un, 0.0);
    }
}

// ============================================================================
// Control rule text (from [CONTROLS] section)
// ============================================================================

#[derive(Debug, Clone, Default)]
pub struct ControlRuleStore {
    /// Each element is the full multi-line text block for one rule:
    /// "RULE name\nIF ...\nTHEN ...\nELSE ...\nPRIORITY ..."
    pub rule_text: Vec<Bytes>,
}

impl ControlRuleStore {
    pub fn count(&self) -> i32 {
        self.rule_text.len() as i32
    }
}
