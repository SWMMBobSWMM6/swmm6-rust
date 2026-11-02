// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/InflowData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson (first pass by tools/soa2rs.py, finished by hand). Not an
// OpenSWMM release.
//
// Authoring stores for [INFLOWS], [DWF], [RDII], [HYDROGRAPHS], [RDII_DECAY] and
// [PATTERNS], one Vec per column in upstream's field order. Strings are bytes. Upstream's
// trailing `name = {}` default argument has no Rust form: callers pass `b""`. A
// parameter named `type` is `type_`. UnitHydEntry has no member initializers upstream,
// so a default-initialized one is indeterminate in C++; Default here gives zeros, which
// is what value-initialization (`UnitHydEntry{}`) gives. RDIIDecayEntry keeps upstream's
// field names (k_T, T_ref, ...), so non_snake_case is allowed on it. PatternData::find
// uses the shared ASCII-only `ieq`.

#![allow(clippy::all)]
use crate::core::options::Bytes;
use crate::data::name_index::ieq;

// ============================================================================
// External inflow definitions (from [INFLOWS] section)
// ============================================================================

/// External inflow authoring rows stored in aligned columns.
#[derive(Debug, Clone, Default)]
pub struct ExtInflowData {
    /// Target node index
    pub node_idx: Vec<i32>,
    /// Raw node name from input (for post-parse re-resolution)
    pub node_name: Vec<Bytes>,
    /// "FLOW" or pollutant name
    pub constituent: Vec<Bytes>,
    /// Timeseries name ("" if none)
    pub ts_name: Vec<Bytes>,
    /// "FLOW","CONCEN","MASS"
    pub inflow_type: Vec<Bytes>,
    /// Multiplier factor
    pub m_factor: Vec<f64>,
    /// Scaling factor
    pub s_factor: Vec<f64>,
    /// Baseline value
    pub baseline: Vec<f64>,
    /// Baseline pattern name
    pub pattern_name: Vec<Bytes>,
}

impl ExtInflowData {
    pub fn count(&self) -> i32 {
        self.node_idx.len() as i32
    }

    /// Append a copied authoring row. `name` is upstream's optional raw node name
    /// (default empty).
    pub fn add(
        &mut self,
        ni: i32,
        cons: &[u8],
        ts: &[u8],
        type_: &[u8],
        mf: f64,
        sf: f64,
        base: f64,
        pat: &[u8],
        name: &[u8],
    ) {
        self.node_idx.push(ni);
        self.node_name.push(name.to_vec());
        self.constituent.push(cons.to_vec());
        self.ts_name.push(ts.to_vec());
        self.inflow_type.push(type_.to_vec());
        self.m_factor.push(mf);
        self.s_factor.push(sf);
        self.baseline.push(base);
        self.pattern_name.push(pat.to_vec());
    }

    /// Erase an aligned row; invalid indices are a no-op.
    pub fn erase(&mut self, idx: i32) {
        if idx < 0 || idx >= self.count() {
            return;
        }
        let u = idx as usize;
        self.node_idx.remove(u);
        self.node_name.remove(u);
        self.constituent.remove(u);
        self.ts_name.remove(u);
        self.inflow_type.remove(u);
        self.m_factor.remove(u);
        self.s_factor.remove(u);
        self.baseline.remove(u);
        self.pattern_name.remove(u);
    }
}

// ============================================================================
// Dry weather flow definitions (from [DWF] section)
// ============================================================================

/// Dry-weather-flow authoring rows stored in aligned columns.
#[derive(Debug, Clone, Default)]
pub struct DwfData {
    /// Target node
    pub node_idx: Vec<i32>,
    /// Raw node name from input (for post-parse re-resolution)
    pub node_name: Vec<Bytes>,
    /// "FLOW" or pollutant name
    pub constituent: Vec<Bytes>,
    /// Average value
    pub avg_value: Vec<f64>,
    /// Monthly pattern name
    pub pat1: Vec<Bytes>,
    /// Daily pattern name
    pub pat2: Vec<Bytes>,
    /// Hourly pattern name
    pub pat3: Vec<Bytes>,
    /// Weekend pattern name
    pub pat4: Vec<Bytes>,
}

impl DwfData {
    pub fn count(&self) -> i32 {
        self.node_idx.len() as i32
    }

    /// Append a copied authoring row. `name` is upstream's optional raw node name
    /// (default empty).
    pub fn add(
        &mut self,
        ni: i32,
        cons: &[u8],
        avg: f64,
        p1: &[u8],
        p2: &[u8],
        p3: &[u8],
        p4: &[u8],
        name: &[u8],
    ) {
        self.node_idx.push(ni);
        self.node_name.push(name.to_vec());
        self.constituent.push(cons.to_vec());
        self.avg_value.push(avg);
        self.pat1.push(p1.to_vec());
        self.pat2.push(p2.to_vec());
        self.pat3.push(p3.to_vec());
        self.pat4.push(p4.to_vec());
    }

    /// Erase an aligned row; invalid indices are a no-op.
    pub fn erase(&mut self, idx: i32) {
        if idx < 0 || idx >= self.count() {
            return;
        }
        let u = idx as usize;
        self.node_idx.remove(u);
        self.node_name.remove(u);
        self.constituent.remove(u);
        self.avg_value.remove(u);
        self.pat1.remove(u);
        self.pat2.remove(u);
        self.pat3.remove(u);
        self.pat4.remove(u);
    }
}

// ============================================================================
// RDII assignments (from [RDII] section)
// ============================================================================

/// Node-to-hydrograph authoring assignments stored in aligned columns.
#[derive(Debug, Clone, Default)]
pub struct RDIIAssignData {
    /// Target node
    pub node_idx: Vec<i32>,
    /// Raw node name from input (for post-parse re-resolution)
    pub node_name: Vec<Bytes>,
    /// Unit hydrograph name
    pub uh_name: Vec<Bytes>,
    /// Tributary sewer area
    pub sewer_area: Vec<f64>,
}

impl RDIIAssignData {
    pub fn count(&self) -> i32 {
        self.node_idx.len() as i32
    }

    /// Append a copied authoring row. `name` is upstream's optional raw node name
    /// (default empty).
    pub fn add(&mut self, ni: i32, uh: &[u8], area: f64, name: &[u8]) {
        self.node_idx.push(ni);
        self.node_name.push(name.to_vec());
        self.uh_name.push(uh.to_vec());
        self.sewer_area.push(area);
    }

    /// Erase an aligned row; invalid indices are a no-op.
    pub fn erase(&mut self, idx: i32) {
        if idx < 0 || idx >= self.count() {
            return;
        }
        let u = idx as usize;
        self.node_idx.remove(u);
        self.node_name.remove(u);
        self.uh_name.remove(u);
        self.sewer_area.remove(u);
    }
}

// ============================================================================
// Unit Hydrograph data (from [HYDROGRAPHS] section)
// ============================================================================

/// No member initializers upstream; Default gives value-initialization's zeros.
#[derive(Debug, Clone, Default)]
pub struct UnitHydEntry {
    /// UH group name
    pub name: Bytes,
    /// Associated rain gage name
    pub gage_name: Bytes,
    /// Month index (0-11, or -1 for ALL)
    pub month: i32,
    /// 0=SHORT, 1=MEDIUM, 2=LONG
    pub response: i32,
    /// Fraction of rainfall volume
    pub r: f64,
    /// Time to peak (hours)
    pub t: f64,
    /// Recession-limb-to-peak-time ratio (tBase = t*(1+k); k >= 0)
    pub k: f64,
    /// Max initial abstraction depth
    pub dmax: f64,
    /// IA recovery rate
    pub drecov: f64,
    /// Initial IA used
    pub dinit: f64,
}

/// Hydrograph parameter rows and independent group-to-gage assignments.
#[derive(Debug, Clone, Default)]
pub struct UnitHydData {
    pub entries: Vec<UnitHydEntry>,

    /// Rain gage names associated with each UH group (name → gage name)
    /// UH group names
    pub gage_assignments: Vec<Bytes>,
    /// Assigned rain gage names
    pub gage_names: Vec<Bytes>,
}

impl UnitHydData {
    pub fn count(&self) -> i32 {
        self.entries.len() as i32
    }

    /// Append a copied group-to-gage assignment without deduplication.
    pub fn add_gage(&mut self, uh_name: &[u8], gage: &[u8]) {
        self.gage_assignments.push(uh_name.to_vec());
        self.gage_names.push(gage.to_vec());
    }

    /// Append a copied parameter row.
    pub fn add(&mut self, e: &UnitHydEntry) {
        self.entries.push(e.clone());
    }
}

// ============================================================================
// RDII exponential-decay parameters (from [RDII_DECAY] section)
// ============================================================================
//
// One row per (UH group, response). Granularity is per-response, NOT per-month.
// A UH group with no row here uses the legacy linear IA recovery for every
// response. A group with one row falls back to linear for the two unspecified
// responses, so adoption is incremental.

#[allow(non_snake_case)]
#[derive(Debug, Clone)]
pub struct RDIIDecayEntry {
    /// Matches UnitHydEntry::name
    pub uh_name: Bytes,
    /// 0=SHORT, 1=MEDIUM, 2=LONG
    pub response: i32,
    /// Depletion rate (1/project rain-depth unit: 1/in or 1/mm) — temperature-independent
    pub k_dep: f64,
    /// Base recovery rate (1/hr)
    pub k_0: f64,
    /// Thermal recovery rate at T_ref (1/hr)
    pub k_T: f64,
    /// Reference temperature (deg C)
    pub T_ref: f64,
    /// Temperature sensitivity (1/deg C)
    pub theta_rec: f64,
    /// Recovery suppressed when T < T_freeze (deg C)
    pub T_freeze: f64,

    // Optional degree-day snow model (SNOW keyword + 2 extra columns).
    // When on, precipitation at T <= snow_T accumulates as SWE (no liquid
    // input); at T > snow_T, melt = min(SWE, snow_ddf*(T - snow_T)*dt_days)
    // is added to rainfall before the IA update (rain-on-snow adds).
    /// Degree-day snow model enabled for this row
    pub snow_on: bool,
    /// Rain/snow partition threshold & melt base (deg C)
    pub snow_T: f64,
    /// Degree-day melt factor (project rain-depth unit/degC/day)
    pub snow_ddf: f64,
}

impl Default for RDIIDecayEntry {
    fn default() -> Self {
        RDIIDecayEntry {
            uh_name: Bytes::new(),
            response: -1,
            k_dep: 0.0,
            k_0: 0.0,
            k_T: 0.0,
            T_ref: 10.0,
            theta_rec: 0.0,
            T_freeze: 0.0,
            snow_on: false,
            snow_T: 1.0,
            snow_ddf: 0.0,
        }
    }
}

/// Exponential IA-decay authoring rows keyed by group and response.
#[derive(Debug, Clone, Default)]
pub struct RDIIDecayData {
    pub entries: Vec<RDIIDecayEntry>,
}

impl RDIIDecayData {
    pub fn count(&self) -> i32 {
        self.entries.len() as i32
    }

    /// Append a copied decay row.
    pub fn add(&mut self, e: &RDIIDecayEntry) {
        self.entries.push(e.clone());
    }
}

// ============================================================================
// Time patterns (from [PATTERNS] section)
// ============================================================================

/// Time-pattern authoring names, types and multiplier vectors.
#[derive(Debug, Clone, Default)]
pub struct PatternData {
    /// Pattern name
    pub names: Vec<Bytes>,
    /// 0=MONTHLY,1=DAILY,2=HOURLY,3=WEEKEND
    pub types: Vec<i32>,
    /// Up to 24 multiplier values
    pub factors: Vec<Vec<f64>>,
}

impl PatternData {
    pub fn count(&self) -> i32 {
        self.names.len() as i32
    }

    /// Append a copied pattern row.
    pub fn add(&mut self, name: &[u8], type_: i32, facs: &[f64]) {
        self.names.push(name.to_vec());
        self.types.push(type_);
        self.factors.push(facs.to_vec());
    }

    /// First ASCII-insensitive name match (first entry wins among duplicates), or -1.
    pub fn find(&self, name: &[u8]) -> i32 {
        for i in 0..self.count() {
            if ieq(&self.names[i as usize], name) {
                return i;
            }
        }
        -1
    }
}
