// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/ArdConfigData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Parsed configuration of the Eulerian ARD transport component (`model.ard`): the
// dispersion model, per-conduit overrides, raw and resolved boundary/source rows, and the
// E5b options. Names and paths are bytes. kLitersPerFt3 is K_LITERS_PER_FT3.

use crate::core::options::Bytes;

/// One raw [TRANSPORT_BOUNDARIES] or [TRANSPORT_SOURCES] row (E5a), stored unresolved.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ArdTransportRow {
    /// node name (boundary) / conduit name (source)
    pub element: Bytes,
    /// MSX species name
    pub species: Bytes,
    pub is_ts: bool,
    /// VALUE mode: conc (BC, species units) or mass rate (source, species mass units/s)
    pub value: f64,
    /// TIMESERIES mode
    pub ts_name: Bytes,
}

// The member initializers (is_ts = false, value = 0.0) are Rust's defaults.

/// Global dispersion coefficient model ([TRANSPORT_OPTIONS] DISPERSION).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum ArdDispersionMode {
    /// no global dispersion (per-conduit overrides still apply)
    Off = 0,
    /// Fischer et al. (1979): D = 0.011 v²B²/(Y·U*), U* = √(gYS)
    Fischer = 1,
    /// uniform user value (display units, len²/s)
    Value = 2,
}

/// Parsed `model.ard` state consumed by ArdEngine::init (phase E3).
#[derive(Debug, Clone)]
pub struct ArdConfigData {
    /// True once a [PROCESS_COMPONENTS] transport.ard registration applied.
    pub configured: bool,
    pub dispersion_mode: ArdDispersionMode,
    /// Global coefficient for VALUE mode (display units, len²/s).
    pub dispersion_value: f64,
    /// Per-conduit overrides (parallel arrays): link index
    pub conduit_disp_link: Vec<i32>,
    /// Per-conduit overrides: coefficient in display units
    pub conduit_disp_value: Vec<f64>,
    // ---- E5a: raw boundary/source rows
    pub boundary_rows: Vec<ArdTransportRow>,
    pub source_rows: Vec<ArdTransportRow>,
    // ---- E5a: resolved rows
    /// node index
    pub bc_node: Vec<i32>,
    /// ReactionData species index
    pub bc_msx: Vec<i32>,
    /// VALUE-mode concentration
    pub bc_value: Vec<f64>,
    /// timeseries table index, -1 ⇒ VALUE
    pub bc_ts: Vec<i32>,
    /// link index (conduit)
    pub src_link: Vec<i32>,
    pub src_msx: Vec<i32>,
    /// VALUE-mode rate, internal conc·ft³/s
    pub src_value: Vec<f64>,
    /// ts index, -1 ⇒ VALUE (ts values are species-mass/s; engine converts)
    pub src_ts: Vec<i32>,
    pub transport_rows_resolved: bool,
    // ---- E5b
    /// [TRANSPORT_OPTIONS] TARGET_DX (display length units; 0 ⇒ FV_CELL_LENGTH rules).
    pub target_dx: f64,
    /// IO3c provenance: this component's file carried the SCALAR_SCHEME alias.
    pub sets_scalar_scheme: bool,
    /// IO3c provenance: this component's file carried the LIMITER alias.
    pub sets_limiter: bool,
    /// [TRANSPORT_OPTIONS] DETAILED_OUTPUT: per-cell CSV sidecar path (empty ⇒ off).
    pub detailed_output_path: Bytes,
}

impl Default for ArdConfigData {
    fn default() -> Self {
        ArdConfigData {
            configured: false,
            dispersion_mode: ArdDispersionMode::Off,
            dispersion_value: 0.0,
            conduit_disp_link: Vec::new(),
            conduit_disp_value: Vec::new(),
            boundary_rows: Vec::new(),
            source_rows: Vec::new(),
            bc_node: Vec::new(),
            bc_msx: Vec::new(),
            bc_value: Vec::new(),
            bc_ts: Vec::new(),
            src_link: Vec::new(),
            src_msx: Vec::new(),
            src_value: Vec::new(),
            src_ts: Vec::new(),
            transport_rows_resolved: false,
            target_dx: 0.0,
            sets_scalar_scheme: false,
            sets_limiter: false,
            detailed_output_path: Vec::new(),
        }
    }
}

impl ArdConfigData {
    /// True for mode != OFF or any conduit_disp_link row.
    pub fn any_dispersion(&self) -> bool {
        self.dispersion_mode != ArdDispersionMode::Off || !self.conduit_disp_link.is_empty()
    }
    /// True if boundary_rows or source_rows is nonempty.
    pub fn any_transport_rows(&self) -> bool {
        !self.boundary_rows.is_empty() || !self.source_rows.is_empty()
    }
    /// True for dispersion, raw transport rows, a detailed-output path or target_dx > 0.
    pub fn any_engine_content(&self) -> bool {
        self.any_dispersion()
            || self.any_transport_rows()
            || !self.detailed_output_path.is_empty()
            || self.target_dx > 0.0
    }
}

/// Litres per cubic foot — the species-mass/s → internal conc·ft³/s source conversion.
pub const K_LITERS_PER_FT3: f64 = 28.316846592;
