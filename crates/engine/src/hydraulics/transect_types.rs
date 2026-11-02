// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Transect.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Irregular transect cross-section geometry (HEC-2 style): the header only, namely the
// table size and the TransectData record. The header has no inline functions; the three
// functions it declares are defined in Transect.cpp and are translated with it (see the
// notes at the end of this file). Upstream's forward declaration
// `struct TransectStore;` names crate::data::infra::TransectStore.
//
// Differences in form: `std::string name` is bytes; the fixed `double[N_TRANSECT_TBL]`
// tables are `[f64; N_TRANSECT_TBL]`; Default is written by hand because
// length_factor starts at 1.0 and Rust has no derived Default for 51-element arrays.

use crate::core::options::Bytes;
use crate::hydraulics::xsect_lookup::LocateLut;

pub const N_TRANSECT_TBL: i32 = 51;

#[derive(Debug, Clone)]
pub struct TransectData {
    pub name: Bytes,
    /// Left overbank Manning's n
    pub n_left: f64,
    /// Right overbank Manning's n
    pub n_right: f64,
    /// Main channel Manning's n
    pub n_channel: f64,
    /// Main-channel / flood-plain length ratio (Lfactor)
    pub length_factor: f64,

    // Station-elevation pairs
    /// x positions (ft)
    pub stations: Vec<f64>,
    /// y elevations (ft)
    pub elevations: Vec<f64>,

    /// Left bank station
    pub x_left_bank: f64,
    /// Right bank station
    pub x_right_bank: f64,

    // Computed geometry tables (normalized)
    pub y_full: f64,
    pub w_max: f64,
    pub a_full: f64,
    pub r_full: f64,
    // CUSTOM shapes only (legacy TShape.sMax/aMax, shape.c getSmax): max section factor
    // and its area for a shape of UNIT height, computed on the unnormalized tables.
    // Scaled by yFull in the caller (xsect.c:685-686).
    pub s_max: f64,
    pub a_max: f64,

    pub area_tbl: [f64; N_TRANSECT_TBL as usize],
    pub hrad_tbl: [f64; N_TRANSECT_TBL as usize],
    pub width_tbl: [f64; N_TRANSECT_TBL as usize],

    /// Bucket map over `area_tbl`, used by getYofA's invLookup (plan XSECT_LOOKUP_ACCEL
    /// A1). Rebuilt by the two table builders; stays disabled (and callers bisect as
    /// before) if either bails out early.
    pub area_lut: LocateLut,
}

impl Default for TransectData {
    fn default() -> Self {
        TransectData {
            name: Bytes::new(),
            n_left: 0.0,
            n_right: 0.0,
            n_channel: 0.0,
            length_factor: 1.0,
            stations: Vec::new(),
            elevations: Vec::new(),
            x_left_bank: 0.0,
            x_right_bank: 0.0,
            y_full: 0.0,
            w_max: 0.0,
            a_full: 0.0,
            r_full: 0.0,
            s_max: 0.0,
            a_max: 0.0,
            area_tbl: [0.0; N_TRANSECT_TBL as usize],
            hrad_tbl: [0.0; N_TRANSECT_TBL as usize],
            width_tbl: [0.0; N_TRANSECT_TBL as usize],
            area_lut: LocateLut::default(),
        }
    }
}

// Declared in Transect.hpp, defined in Transect.cpp; translated with Transect.cpp in
// hydraulics::transect (build_tables, build_custom_tables, build_from_store):
//
//   void buildTables(TransectData& td, bool add_end_walls = true);
//       Build tabulated geometry (area, hydraulic radius, width at N_TRANSECT_TBL equal
//       depth increments, normalized) from station-elevation data; add_end_walls closes
//       the section with vertical walls as legacy transect_validate does (false for
//       STREET sections).
//
//   void buildCustomTables(TransectData& td, double y_full, const double* curve_x,
//                          const double* curve_y, int n_pts);
//       Build tabulated geometry from a CUSTOM shape curve (y/yFull vs w/wMax).
//
//   bool buildFromStore(const TransectStore& ts, int index, double ucf_length,
//                       TransectData& td);
//       Copy one [TRANSECTS] store row into td, apply the legacy parse-time transforms
//       to the copy, then buildTables(); false when the channel Manning's n is invalid.
