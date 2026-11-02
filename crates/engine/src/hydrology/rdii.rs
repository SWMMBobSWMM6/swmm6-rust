// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/RDII.hpp and RDII.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// RDII by unit-hydrograph convolution: the parameter and state types, RDIIGroupSoA,
// and RDIISolver::init with what it calls (addUnitHydParams, findUnitHyd,
// getRainInterval, getMaxPeriods, validateExpDecay). Upstream calls rdii_.init from
// SWMMEngine::initHydraulics (step 10a), before initHydrology.
//
// Differences of form:
//   * UhNameMap is `std::unordered_map<std::string, int, CiHash, CiEqual>`. Here it is a
//     Vec of (key, index) in insertion order with ASCII case-insensitive lookup; a key
//     keeps its first spelling, as operator[] keeps it. HAZARD: upstream iterates this
//     map in SWMMEngine::initHydraulics (Gap #80, the UH validation), and that order is
//     libstdc++'s bucket order for CiHash (FNV-1a over upper-cased bytes); whoever
//     translates that loop must reproduce it. Nothing in init iterates it.
//   * The init-local maps (params_seen, uh_gage_map) are only looked up.
//   * `static_cast<long>(double)` is cvt_i64 below (the reference build gives LONG_MIN
//     for NaN or out-of-range values, where `as i64` saturates).
//   * The members are public so the gate dumper can read them; upstream's groups_,
//     uh_name_to_idx_, node_rdii_flow_ and the grid members are private.
//   * NOT YET TRANSLATED (run-time step functions, phase 1c): getRecoveryRate,
//     updateIA_exp, uhOrdinate, advance, emitTick, applyRdiiInflows and the file-static
//     helpers (updateIA_linear, fToC, updateDryPeriod, ...).

#![allow(non_snake_case)]

use crate::core::cmath::{cmax, cmin, cvt_i32};
use crate::core::context::SimulationContext;
use crate::core::datetime;
use crate::core::error_codes::{format_warning, WARN_UH_PARAMS_REPEATED};
use crate::core::options::Bytes;
use crate::data::name_index::ieq;

#[derive(Debug, Clone, Copy, Default)]
pub struct UnitHydParams {
    /// Rainfall fraction per month x response
    pub r: [[f64; 3]; 12],
    /// Time to peak (sec) per month x response
    pub tPeak: [[f64; 3]; 12],
    /// Base time (sec) per month x response
    pub tBase: [[f64; 3]; 12],
    /// Maximum abstraction, project rain-depth units
    pub iaMax: [[f64; 3]; 12],
    /// Linear recovery, project rain-depth units/day
    pub iaRecov: [[f64; 3]; 12],
    /// Initial used abstraction, project rain-depth units
    pub iaInit: [[f64; 3]; 12],
}

/// Exponential IA decay parameters for one RDII response (SHORT/MEDIUM/LONG).
#[derive(Debug, Clone, Copy)]
pub struct ExpDecayParams {
    pub active: bool,
    pub k_dep: f64,
    pub k_0: f64,
    pub k_T: f64,
    pub T_ref: f64,
    pub theta_rec: f64,
    pub T_freeze: f64,
    pub snow_on: bool,
    pub snow_T: f64,
    pub snow_ddf: f64,
}

impl Default for ExpDecayParams {
    fn default() -> Self {
        ExpDecayParams { active: false, k_dep: 0.0, k_0: 0.0, k_T: 0.0, T_ref: 10.0, theta_rec: 0.0, T_freeze: 0.0, snow_on: false, snow_T: 1.0, snow_ddf: 0.0 }
    }
}

/// Per-response (SHORT/MEDIUM/LONG) unit hydrograph data (legacy TUHData).
#[derive(Debug, Clone, Default)]
pub struct UHResponseData {
    /// circular buffer of past rainfall depths
    pub past_rain: Vec<f64>,
    /// month for each past rainfall entry
    pub past_month: Vec<i32>,
    /// current buffer write position
    pub period: i32,
    /// buffer capacity
    pub max_periods: i32,
    /// true if any non-zero past rain
    pub has_past_rain: i32,
    /// initial abstraction used so far
    pub ia_used: f64,
    /// snow water equivalent (project depth; degree-day snow model)
    pub swe: f64,
    /// seconds since last non-zero rainfall
    pub dry_seconds: i64,
}

impl UHResponseData {
    pub fn allocate(&mut self, n: i32) {
        self.max_periods = n;
        self.past_rain.clear();
        self.past_rain.resize(n as usize, 0.0);
        self.past_month.clear();
        self.past_month.resize(n as usize, 0);
        self.period = 0;
        self.has_past_rain = 0;
        self.ia_used = 0.0;
        self.swe = 0.0;
        self.dry_seconds = (n as i64) * 300 + 1; // start dry
    }
}

#[derive(Debug, Clone, Default)]
pub struct RDIIGroupSoA {
    pub count: i32,
    /// Target node index
    pub node_idx: Vec<i32>,
    /// Unit hydrograph parameter index
    pub uh_idx: Vec<i32>,
    /// Rain gage index per UH group (legacy: UnitHyd[j].rainGage)
    pub gage_idx: Vec<i32>,
    /// Contributing area, project land-area units (acres/ha)
    pub area: Vec<f64>,
    /// Per-response data: [group * 3 + response]
    pub uh_data: Vec<UHResponseData>,
    /// Rain processing interval (sec) per group
    pub rain_interval: Vec<i32>,
    /// Chunk cursor (sec, legacy UHGroup.gageDate)
    pub gage_elapsed: Vec<f64>,
    /// Chunk starts with recorded rates (sec)
    pub rate_recorded_until: Vec<f64>,
    /// FIFO of recorded chunk rates
    pub pending_rates: Vec<Vec<f64>>,
}

fn assign<T: Clone>(v: &mut Vec<T>, n: usize, d: T) {
    v.clear();
    v.resize(n, d);
}

impl RDIIGroupSoA {
    pub fn resize(&mut self, n: i32) {
        self.count = n;
        let un = n as usize;
        assign(&mut self.node_idx, un, -1);
        assign(&mut self.uh_idx, un, -1);
        assign(&mut self.gage_idx, un, -1);
        assign(&mut self.area, un, 0.0);
        self.uh_data.resize(un * 3, UHResponseData::default()); // 3 responses per group
        assign(&mut self.rain_interval, un, 300);
        assign(&mut self.gage_elapsed, un, 0.0);
        assign(&mut self.rate_recorded_until, un, 0.0);
        assign(&mut self.pending_rates, un, Vec::new());
    }
}

/// `static_cast<long>(x)` as the reference build computes it: LONG_MIN for NaN and
/// out-of-range values.
fn cvt_i64(x: f64) -> i64 {
    if x.is_nan() || x >= 9.223372036854775808e18 || x < -9.223372036854775808e18 {
        i64::MIN
    } else {
        x as i64
    }
}

/// A case-insensitive name map (CiHash/CiEqual), insertion ordered.
fn ci_find(map: &[(Bytes, i32)], name: &[u8]) -> Option<usize> {
    map.iter().position(|(k, _)| ieq(k, name))
}

#[derive(Debug, Clone)]
pub struct RDIISolver {
    pub uh_params: Vec<UnitHydParams>,
    /// Exponential decay params per UH group x response ([SHORT, MEDIUM, LONG]).
    pub decay_params: Vec<[ExpDecayParams; 3]>,

    pub groups_: RDIIGroupSoA,
    /// UH group name -> index (case-insensitive), insertion order
    pub uh_name_to_idx_: Vec<(Bytes, i32)>,
    /// Latest tick's per-node RDII flow (CFS)
    pub node_rdii_flow_: Vec<f64>,

    /// RdiiStep (sec) = wet_step
    pub grid_step_: f64,
    /// Next un-emitted tick (elapsed sec)
    pub next_tick_: f64,
    /// Column -> node index (sorted unique)
    pub grid_node_: Vec<i32>,
    /// Emitted rows x grid_node_.size(), REAL4
    pub grid_flows_: Vec<f32>,
}

impl Default for RDIISolver {
    fn default() -> Self {
        RDIISolver { uh_params: Vec::new(), decay_params: Vec::new(), groups_: RDIIGroupSoA::default(), uh_name_to_idx_: Vec::new(), node_rdii_flow_: Vec::new(), grid_step_: 0.0, next_tick_: 0.0, grid_node_: Vec::new(), grid_flows_: Vec::new() }
    }
}

const RESPONSE_NAMES: [&[u8]; 3] = [b"SHORT", b"MEDIUM", b"LONG"];

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

impl RDIISolver {
    pub fn add_unit_hyd_params(&mut self, name: &[u8], params: &UnitHydParams) -> i32 {
        if let Some(p) = ci_find(&self.uh_name_to_idx_, name) {
            let idx = self.uh_name_to_idx_[p].1;
            self.uh_params[idx as usize] = *params;
            return idx;
        }
        let idx = self.uh_params.len() as i32;
        self.uh_params.push(*params);
        self.uh_name_to_idx_.push((name.to_vec(), idx));
        idx
    }

    pub fn find_unit_hyd(&self, name: &[u8]) -> i32 {
        match ci_find(&self.uh_name_to_idx_, name) {
            None => -1,
            Some(p) => self.uh_name_to_idx_[p].1,
        }
    }

    /// The UH name map (see the module notes on its iteration order).
    pub fn uh_name_index(&self) -> &Vec<(Bytes, i32)> {
        &self.uh_name_to_idx_
    }

    pub fn node_flows(&self) -> &Vec<f64> {
        &self.node_rdii_flow_
    }

    pub fn rdii_node_list(&self) -> Vec<i32> {
        self.grid_node_.clone()
    }

    pub fn get_rain_interval(uh: &UnitHydParams, wet_step: f64) -> i32 {
        let mut ri = cvt_i32(wet_step);
        if ri <= 0 {
            ri = 300;
        }
        for m in 0..12 {
            for k in 0..3 {
                let tp = uh.tPeak[m][k];
                if tp > 0.0 {
                    let t_limb = cvt_i32(tp);
                    if t_limb > 0 {
                        ri = cmin(ri, t_limb);
                    }
                    let t_fall = cvt_i32(uh.tBase[m][k] - tp);
                    if t_fall > 0 {
                        ri = cmin(ri, t_fall);
                    }
                }
            }
        }
        cmax(ri, 1)
    }

    pub fn get_max_periods(uh: &UnitHydParams, response: i32, rain_interval: i32) -> i32 {
        if rain_interval <= 0 {
            return 0;
        }
        let k = (response % 3) as usize;
        let mut n_max = 0;
        for m in 0..12 {
            let n = cvt_i32(uh.tBase[m][k] / rain_interval as f64).wrapping_add(1);
            n_max = cmax(n, n_max);
        }
        n_max
    }

    /// Warn when exponential decay is configured without a temperature source, and on
    /// degenerate per-row configurations.
    fn validate_exp_decay(&self, ctx: &mut SimulationContext) {
        let mut any_active = false;
        for triple in &self.decay_params {
            for dp in triple.iter() {
                if dp.active {
                    any_active = true;
                    break;
                }
            }
            if any_active {
                break;
            }
        }
        if any_active && ctx.options.temp_source == 0 {
            ctx.warnings.push(
                b"WARNING: [RDII_DECAY] exponential IA model is active but no temperature source is configured. Recovery rate will be evaluated at T_ref for every group; no seasonal variation will be produced."
                    .to_vec(),
            );
        }

        let entries = ctx.rdii_decay.entries.clone();
        for e in &entries {
            if self.find_unit_hyd(&e.uh_name) < 0 {
                continue;
            }
            if e.response < 0 || e.response > 2 {
                continue;
            }
            let resp = RESPONSE_NAMES[e.response as usize];
            if e.k_dep == 0.0 {
                ctx.warnings.push(cat(&[b"WARNING: [RDII_DECAY] k_dep = 0 for group '", &e.uh_name, b"' (", resp, b"): initial abstraction is disabled; excess equals rainfall and the recovery parameters have no effect."]));
            }
            if e.T_freeze >= e.T_ref {
                ctx.warnings.push(cat(&[b"WARNING: [RDII_DECAY] T_freeze >= T_ref for group '", &e.uh_name, b"' (", resp, b"): recovery is suppressed at the reference temperature; check the freeze threshold."]));
            }
            if e.snow_on && ctx.options.temp_source == 0 {
                if e.T_ref <= e.snow_T {
                    ctx.warnings.push(cat(&[
                        b"WARNING: [RDII_DECAY] snow model is on for group '",
                        &e.uh_name,
                        b"' (",
                        resp,
                        b") but no temperature source is configured and T_ref <= snow_T: all rainfall will accumulate as snow and never melt, producing zero RDII for this response.",
                    ]));
                } else {
                    ctx.warnings.push(cat(&[
                        b"WARNING: [RDII_DECAY] snow model is on for group '",
                        &e.uh_name,
                        b"' (",
                        resp,
                        b") but no temperature source is configured: temperature is fixed at T_ref, so no snow will ever accumulate and the snow parameters have no effect.",
                    ]));
                }
            }
        }
    }

    /// Populate UH params from parsed data and allocate per-response buffers.
    pub fn init(&mut self, ctx: &mut SimulationContext) {
        // Which (month, response) slots each UH group has already been assigned
        // (legacy WARN13 from rdii_readUnitHydParams).
        let mut params_seen: std::collections::HashMap<i32, [[bool; 3]; 12]> = std::collections::HashMap::new();

        // Populate UH params from parsed [HYDROGRAPHS] data
        let entries = ctx.unit_hyds.entries.clone();
        for entry in &entries {
            let mut idx = self.find_unit_hyd(&entry.name);
            if idx < 0 {
                let params = UnitHydParams::default();
                idx = self.add_unit_hyd_params(&entry.name, &params);
            }

            // PARITY rdii.c:424-426 — base time formed in HOURS, both times truncated to
            // whole seconds.
            let t_base_hrs = entry.t * (1.0 + entry.k);
            let t_peak_sec = cvt_i64(entry.t * 3600.0) as f64;
            let t_base_sec = cvt_i64(t_base_hrs * 3600.0) as f64;

            let m_start = if entry.month < 0 { 0 } else { entry.month };
            let m_end = if entry.month < 0 { 11 } else { entry.month };
            let k = entry.response as usize;

            // Warn once per entry that re-assigns an already-assigned month.
            let seen = params_seen.entry(idx).or_insert([[false; 3]; 12]);
            let mut repeated = false;
            for m in m_start..=m_end {
                if seen[m as usize][k] {
                    repeated = true;
                }
                seen[m as usize][k] = true;
            }
            if repeated {
                ctx.warnings.push(format_warning(WARN_UH_PARAMS_REPEATED, &entry.name));
            }

            let uh = &mut self.uh_params[idx as usize];
            for m in m_start..=m_end {
                let m = m as usize;
                uh.r[m][k] = entry.r;
                uh.tPeak[m][k] = t_peak_sec;
                uh.tBase[m][k] = t_base_sec;
                uh.iaMax[m][k] = entry.dmax;
                uh.iaRecov[m][k] = entry.drecov;
                uh.iaInit[m][k] = entry.dinit;
            }
        }

        // A gage-only UH is a valid zero-response group in legacy.
        let gage_assignments = ctx.unit_hyds.gage_assignments.clone();
        for name in &gage_assignments {
            if self.find_unit_hyd(name) < 0 {
                self.add_unit_hyd_params(name, &UnitHydParams::default());
            }
        }

        // decay_params parallel to uh_params; defaults are inactive (linear IA).
        self.decay_params.clear();
        self.decay_params.resize(self.uh_params.len(), [ExpDecayParams::default(); 3]);

        // Populate decay_params from parsed [RDII_DECAY] data
        for d in &ctx.rdii_decay.entries {
            let idx = self.find_unit_hyd(&d.uh_name);
            if idx < 0 {
                continue; // unknown UH group — silently skip
            }
            if d.response < 0 || d.response > 2 {
                continue;
            }
            let dp = &mut self.decay_params[idx as usize][d.response as usize];
            dp.active = true;
            dp.k_dep = d.k_dep;
            dp.k_0 = d.k_0;
            dp.k_T = d.k_T;
            dp.T_ref = d.T_ref;
            dp.theta_rec = d.theta_rec;
            dp.T_freeze = d.T_freeze;
            dp.snow_on = d.snow_on;
            dp.snow_T = d.snow_T;
            dp.snow_ddf = d.snow_ddf;
        }
        self.validate_exp_decay(ctx);

        // UH name -> gage index from the [HYDROGRAPHS] gage lines (case-insensitive).
        let mut uh_gage_map: Vec<(Bytes, i32)> = Vec::new();
        for gi in 0..ctx.unit_hyds.gage_assignments.len() {
            let uh_name = &ctx.unit_hyds.gage_assignments[gi];
            let gage_name = &ctx.unit_hyds.gage_names[gi];
            let gidx = ctx.gage_names.find(gage_name);
            if gidx >= 0 {
                match ci_find(&uh_gage_map, uh_name) {
                    Some(p) => uh_gage_map[p].1 = gidx,
                    None => uh_gage_map.push((uh_name.clone(), gidx)),
                }
            }
        }

        let assigns = &ctx.rdii_assigns;
        let n_assigns = assigns.count();

        // Strict-grid state — reset for re-initialization.
        self.grid_step_ = 0.0;
        self.next_tick_ = 0.0;
        self.grid_node_.clear();
        self.grid_flows_.clear();

        if n_assigns == 0 {
            self.groups_.count = 0;
            return;
        }

        self.groups_.resize(n_assigns);
        self.node_rdii_flow_.clear();
        self.node_rdii_flow_.resize(ctx.n_nodes() as usize, 0.0);
        let wet_step = ctx.options.wet_step;
        self.grid_step_ = wet_step;

        for i in 0..n_assigns {
            let ui = i as usize;
            self.groups_.node_idx[ui] = assigns.node_idx[ui];
            self.groups_.area[ui] = assigns.sewer_area[ui];

            let uh_name = &assigns.uh_name[ui];
            let uh_i = self.find_unit_hyd(uh_name);
            self.groups_.uh_idx[ui] = uh_i;

            // Resolve per-UH rain gage (legacy: UnitHyd[j].rainGage)
            self.groups_.gage_idx[ui] = match ci_find(&uh_gage_map, uh_name) {
                Some(p) => uh_gage_map[p].1,
                None => -1,
            };

            if uh_i < 0 || uh_i >= self.uh_params.len() as i32 {
                self.groups_.rain_interval[ui] = 300;
                for k in 0..3 {
                    self.groups_.uh_data[ui * 3 + k].max_periods = 0;
                }
                continue;
            }

            let uh = self.uh_params[uh_i as usize];
            let rain_interval = Self::get_rain_interval(&uh, wet_step);
            self.groups_.rain_interval[ui] = rain_interval;

            // Starting month for iaInit initialization
            let start_month = (datetime::month_of_year(ctx.options.start_date) - 1) as usize;

            for k in 0..3 {
                let max_per = Self::get_max_periods(&uh, k as i32, rain_interval);
                let rd = &mut self.groups_.uh_data[ui * 3 + k];
                if max_per > 0 {
                    rd.allocate(max_per);
                    // legacy rdii.c 1079-1085
                    rd.dry_seconds = (max_per as i64) * (rain_interval as i64) + 1;
                    rd.period = max_per + 1;
                    rd.ia_used = uh.iaInit[start_month][k];
                }
            }
        }

        // Grid columns: the RDII node set (legacy RdiiNodeIndex file header).
        for &ni in &self.groups_.node_idx {
            if ni >= 0 {
                self.grid_node_.push(ni);
            }
        }
        self.grid_node_.sort();
        self.grid_node_.dedup();
    }
}
