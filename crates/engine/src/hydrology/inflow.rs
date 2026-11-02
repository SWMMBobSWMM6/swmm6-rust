// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/Inflow.hpp and Inflow.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// External and dry-weather inflows: the runtime SoA blocks and InflowSolver::init with
// refreshPatterns.
//
// NOT YET TRANSLATED (run-time step functions, with the routing step):
// InflowSolver::computeAll and getPatternFactor.
//
// Differences in form, never in result:
//   * The case-insensitive `unordered_map<string,int,CiHash,CiEqual>` pattern map and the
//     `unordered_map<long long,int>` supersede map are only looked up, never iterated, so
//     they are HashMaps keyed by the ASCII-uppercased name and by the key value.
//   * Members keep upstream's names and are `pub` (upstream's private ones too) so the
//     gate's state dumper can read them.

use std::collections::HashMap;

use crate::core::context::SimulationContext;
use crate::core::options::Bytes;
use crate::core::unit_conversion as ucf;
use crate::data::water_age::WaterAgeSource;

pub const MONTHLY_PATTERN: usize = 0;
pub const DAILY_PATTERN: usize = 1;
pub const HOURLY_PATTERN: usize = 2;
pub const WEEKEND_PATTERN: usize = 3;

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExtInflowKind {
    Flow = 0,
    Concen = 1,
    Mass = 2,
    Age = 3,
    MsxConcen = 4,
    MsxMass = 5,
}

#[derive(Debug, Clone, Default)]
pub struct ExtInflowSoA {
    pub count: i32,
    pub node_idx: Vec<i32>,
    pub ts_idx: Vec<i32>,
    pub base_pat_idx: Vec<i32>,
    pub baseline: Vec<f64>,
    pub scale_factor: Vec<f64>,
    pub conv_factor: Vec<f64>,
    pub kind: Vec<i32>,
    pub pollut_idx: Vec<i32>,
}

fn assign<T: Clone>(v: &mut Vec<T>, n: usize, d: T) {
    v.clear();
    v.resize(n, d);
}

impl ExtInflowSoA {
    pub fn resize(&mut self, n: i32) {
        self.count = n;
        let un = n as usize;
        assign(&mut self.node_idx, un, -1);
        assign(&mut self.ts_idx, un, -1);
        assign(&mut self.base_pat_idx, un, -1);
        assign(&mut self.baseline, un, 0.0);
        assign(&mut self.scale_factor, un, 1.0);
        assign(&mut self.conv_factor, un, 1.0);
        assign(&mut self.kind, un, ExtInflowKind::Flow as i32);
        assign(&mut self.pollut_idx, un, -1);
    }
}

#[derive(Debug, Clone, Default)]
pub struct DwfInflowSoA {
    pub count: i32,
    pub node_idx: Vec<i32>,
    pub avg_value: Vec<f64>,
    pub pat_monthly: Vec<i32>,
    pub pat_daily: Vec<i32>,
    pub pat_hourly: Vec<i32>,
    pub pat_weekend: Vec<i32>,
    pub is_flow: Vec<u8>,
    pub pollut_idx: Vec<i32>,
    pub msx_idx: Vec<i32>,
}

impl DwfInflowSoA {
    pub fn resize(&mut self, n: i32) {
        self.count = n;
        let un = n as usize;
        assign(&mut self.node_idx, un, -1);
        assign(&mut self.avg_value, un, 0.0);
        assign(&mut self.pat_monthly, un, -1);
        assign(&mut self.pat_daily, un, -1);
        assign(&mut self.pat_hourly, un, -1);
        assign(&mut self.pat_weekend, un, -1);
        assign(&mut self.is_flow, un, 0);
        assign(&mut self.pollut_idx, un, -1);
        assign(&mut self.msx_idx, un, -1);
    }
}

/// 0=monthly, 1=daily, 2=hourly, 3=weekend; up to 24 factors.
#[derive(Debug, Clone, Copy, Default)]
pub struct TimePattern {
    pub type_: i32,
    pub factors: [f64; 24],
}

#[derive(Debug, Clone, Default)]
pub struct InflowSolver {
    pub ext_inflows_: ExtInflowSoA,
    pub dwf_inflows_: DwfInflowSoA,
    pub patterns_: Vec<TimePattern>,
}

/// The legacy case-insensitive FLOW keyword prefix: leading ASCII spaces (only), then
/// FLOW as a four-byte prefix.
fn matches_flow_keyword(s: &[u8]) -> bool {
    const K_FLOW: &[u8; 4] = b"FLOW";
    let mut i = match s.iter().position(|&c| c != b' ') {
        None => return false,
        Some(i) => i,
    };
    for &k in K_FLOW.iter() {
        if i >= s.len() {
            return false;
        }
        if s[i].to_ascii_uppercase() != k {
            return false;
        }
        i += 1;
    }
    true
}

fn ci_key(s: &[u8]) -> Bytes {
    s.iter().map(|c| c.to_ascii_uppercase()).collect()
}

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

impl InflowSolver {
    pub fn refresh_patterns(&mut self, ctx: &SimulationContext) {
        let np = ctx.patterns.count();
        self.patterns_.resize(np as usize, TimePattern::default());
        for i in 0..np {
            let ui = i as usize;
            self.patterns_[ui].type_ = ctx.patterns.types[ui];
            let facs = &ctx.patterns.factors[ui];
            for k in 0..24 {
                self.patterns_[ui].factors[k] = 1.0;
            }
            let mut k = 0;
            while k < facs.len() && k < 24 {
                self.patterns_[ui].factors[k] = facs[k];
                k += 1;
            }
        }
    }

    pub fn init(&mut self, ctx: &mut SimulationContext) {
        // Pattern name -> index (case-insensitive, last duplicate wins).
        let mut pattern_map: HashMap<Bytes, i32> = HashMap::new();
        let np = ctx.patterns.count();
        for i in 0..np {
            pattern_map.insert(ci_key(&ctx.patterns.names[i as usize]), i);
        }

        self.refresh_patterns(ctx);

        if ctx.options.water_age {
            let n = ctx.n_nodes() as usize;
            assign(&mut ctx.water_age_state.node_ext_inflow_age, n, f64::NAN); // quiet_NaN
        }

        // ---- Populate external inflows with name resolution ----
        let ne = ctx.ext_inflows.count();
        self.ext_inflows_.resize(ne);
        for i in 0..ne {
            let ui = i as usize;
            self.ext_inflows_.node_idx[ui] = ctx.ext_inflows.node_idx[ui];
            self.ext_inflows_.baseline[ui] = ctx.ext_inflows.baseline[ui];
            self.ext_inflows_.scale_factor[ui] = ctx.ext_inflows.s_factor[ui];

            let mut cf = ctx.ext_inflows.m_factor[ui];
            let cons = ctx.ext_inflows.constituent[ui].clone();
            let upper = |s: &[u8]| -> Bytes { s.iter().map(|c| c.to_ascii_uppercase()).collect() };
            // Warnings dedupe against ctx.warnings.
            let warn_once = |ctx: &mut SimulationContext, msg: Bytes| {
                for w in &ctx.warnings {
                    if *w == msg {
                        return;
                    }
                }
                ctx.warnings.push(msg);
            };
            if ctx.pollutant_names.find(&cons) < 0 && matches_flow_keyword(&cons) {
                self.ext_inflows_.kind[ui] = ExtInflowKind::Flow as i32;
                self.ext_inflows_.pollut_idx[ui] = -1;
                let fu = ctx.options.flow_units as i32;
                if (0..6).contains(&fu) {
                    cf /= ucf::QCF[fu as usize];
                }
            } else if cons == b"__WATER_AGE__" {
                self.ext_inflows_.kind[ui] = ExtInflowKind::Age as i32;
                self.ext_inflows_.pollut_idx[ui] = -1;
                cf *= 3600.0; // hours -> seconds
                let type_u = upper(&ctx.ext_inflows.inflow_type[ui]);
                if type_u == b"MASS" {
                    let m = cat(&[
                        b"[INFLOWS] __WATER_AGE__ at node '",
                        &ctx.ext_inflows.node_name[ui],
                        "': MASS has no meaning for the age species \u{2014} the value is taken as the inflow's age in hours (CONCEN).".as_bytes(),
                    ]);
                    warn_once(ctx, m);
                }
                if !ctx.options.water_age {
                    warn_once(ctx, "[INFLOWS] __WATER_AGE__ rows are present but [OPTIONS] WATER_AGE is OFF \u{2014} the rows are inert this simulation.".as_bytes().to_vec());
                }
                for k in 0..ctx.water_age_config.node_over_source.len() {
                    if ctx.water_age_config.node_over_source[k] == WaterAgeSource::ExternalInflow as i32 && ctx.water_age_config.node_over_node[k] == self.ext_inflows_.node_idx[ui] {
                        let m = cat(&[
                            b"[INFLOWS] __WATER_AGE__ at node '",
                            &ctx.ext_inflows.node_name[ui],
                            "' overrides the [WATER_AGE_SOURCES] EXTERNAL_INFLOW row for the same node \u{2014} the inflow row wins (amendment D-Y4).".as_bytes(),
                        ]);
                        warn_once(ctx, m);
                    }
                }
            } else if cons == b"__TEMPERATURE__" {
                self.ext_inflows_.kind[ui] = ExtInflowKind::Concen as i32;
                self.ext_inflows_.pollut_idx[ui] = -1;
                let m = cat(&[
                    "[INFLOWS] __TEMPERATURE__ rows arrive with heat's round \u{2014} the row at node '".as_bytes(),
                    &ctx.ext_inflows.node_name[ui],
                    b"' is ignored.",
                ]);
                warn_once(ctx, m);
            } else {
                let type_u = upper(&ctx.ext_inflows.inflow_type[ui]);
                let is_mass = type_u == b"MASS";
                self.ext_inflows_.kind[ui] = if is_mass { ExtInflowKind::Mass } else { ExtInflowKind::Concen } as i32;
                self.ext_inflows_.pollut_idx[ui] = ctx.pollutant_names.find(&cons);
                if self.ext_inflows_.pollut_idx[ui] < 0 && ctx.reactions.configured {
                    let m = ctx.reactions.find_species(&cons);
                    if m >= 0 {
                        self.ext_inflows_.kind[ui] = if is_mass { ExtInflowKind::MsxMass } else { ExtInflowKind::MsxConcen } as i32;
                        self.ext_inflows_.pollut_idx[ui] = m;
                        if (m as usize) < ctx.reactions.species_is_wall.len() && ctx.reactions.species_is_wall[m as usize] != 0 {
                            let msg = cat(&[
                                b"[INFLOWS] '",
                                &cons,
                                b"' at node '",
                                &ctx.ext_inflows.node_name[ui],
                                "' is a WALL species \u{2014} inflow water carries no wall-bound mass; the row is ignored.".as_bytes(),
                            ]);
                            warn_once(ctx, msg);
                        }
                    }
                }
                if self.ext_inflows_.pollut_idx[ui] < 0 {
                    let msg = cat(&[
                        b"[INFLOWS] constituent '",
                        &cons,
                        b"' at node '",
                        &ctx.ext_inflows.node_name[ui],
                        "' matches no pollutant, reactions species or reserved species \u{2014} the row is ignored.".as_bytes(),
                    ]);
                    warn_once(ctx, msg);
                }
                const L_PER_FT3: f64 = 28.317; // legacy consts.h LperFT3
                if is_mass {
                    cf /= L_PER_FT3;
                }
            }
            self.ext_inflows_.conv_factor[ui] = cf;

            // Resolve timeseries name -> table index
            let ts_name = &ctx.ext_inflows.ts_name[ui];
            self.ext_inflows_.ts_idx[ui] = if !ts_name.is_empty() { ctx.find_timeseries(ts_name) } else { -1 };

            // Resolve baseline pattern name -> pattern index
            let pat_name = &ctx.ext_inflows.pattern_name[ui];
            self.ext_inflows_.base_pat_idx[ui] = if !pat_name.is_empty() { pattern_map.get(&ci_key(pat_name)).copied().unwrap_or(-1) } else { -1 };
        }

        // ---- Populate DWF inflows with name resolution + pattern sorting ----
        let nd = ctx.dwf_inflows.count();
        self.dwf_inflows_.resize(nd);
        for i in 0..nd {
            let ui = i as usize;
            self.dwf_inflows_.node_idx[ui] = ctx.dwf_inflows.node_idx[ui];

            let mut avg_val = ctx.dwf_inflows.avg_value[ui];
            let constituent = &ctx.dwf_inflows.constituent[ui];
            let dwf_pollut = ctx.pollutant_names.find(constituent);
            if dwf_pollut < 0 && matches_flow_keyword(constituent) {
                let fu = ctx.options.flow_units as usize;
                avg_val /= ucf::QCF[fu];
                self.dwf_inflows_.is_flow[ui] = 1;
            } else {
                self.dwf_inflows_.pollut_idx[ui] = dwf_pollut;
                if self.dwf_inflows_.pollut_idx[ui] < 0 && ctx.reactions.configured {
                    self.dwf_inflows_.msx_idx[ui] = ctx.reactions.find_species(constituent);
                }
            }
            self.dwf_inflows_.avg_value[ui] = avg_val;

            let mut tmp_pats = [-1i32; 4];
            let pat_fields = [&ctx.dwf_inflows.pat1[ui], &ctx.dwf_inflows.pat2[ui], &ctx.dwf_inflows.pat3[ui], &ctx.dwf_inflows.pat4[ui]];
            for p in pat_fields {
                if p.is_empty() {
                    continue;
                }
                let pat_idx = match pattern_map.get(&ci_key(p)) {
                    None => continue,
                    Some(&v) => v,
                };
                let pat_type = ctx.patterns.types[pat_idx as usize];
                if (0..4).contains(&pat_type) {
                    tmp_pats[pat_type as usize] = pat_idx;
                }
            }

            self.dwf_inflows_.pat_monthly[ui] = tmp_pats[MONTHLY_PATTERN];
            self.dwf_inflows_.pat_daily[ui] = tmp_pats[DAILY_PATTERN];
            self.dwf_inflows_.pat_hourly[ui] = tmp_pats[HOURLY_PATTERN];
            self.dwf_inflows_.pat_weekend[ui] = tmp_pats[WEEKEND_PATTERN];
        }

        // ---- A repeated (node, constituent) row REPLACES the earlier one ----
        {
            const K_NONE: i64 = i64::MIN;
            const K_MSX_BASE: i64 = 1i64 << 32;
            const K_AGE_KEY: i64 = 1i64 << 33;
            let ext_param = |kind: i32, idx: i32| -> i64 {
                match kind {
                    k if k == ExtInflowKind::Flow as i32 => -1,
                    k if k == ExtInflowKind::Age as i32 => K_AGE_KEY,
                    k if k == ExtInflowKind::MsxConcen as i32 || k == ExtInflowKind::MsxMass as i32 => {
                        if idx >= 0 {
                            K_MSX_BASE + idx as i64
                        } else {
                            K_NONE
                        }
                    }
                    _ => {
                        if idx >= 0 {
                            idx as i64
                        } else {
                            K_NONE
                        }
                    }
                }
            };
            let mut last: HashMap<i64, i32> = HashMap::new();
            let supersede = |last: &mut HashMap<i64, i32>, node_idx: &mut Vec<i32>, row: i32, param: i64| {
                if param == K_NONE {
                    return;
                }
                // (node << 34) + param + 3, in long long arithmetic
                let key = ((node_idx[row as usize] as i64) << 34).wrapping_add(param).wrapping_add(3);
                if let Some(&prev) = last.get(&key) {
                    node_idx[prev as usize] = -1;
                }
                last.insert(key, row);
            };

            for i in 0..ne {
                let ui = i as usize;
                if self.ext_inflows_.node_idx[ui] < 0 {
                    continue;
                }
                let p = ext_param(self.ext_inflows_.kind[ui], self.ext_inflows_.pollut_idx[ui]);
                supersede(&mut last, &mut self.ext_inflows_.node_idx, i, p);
            }

            last.clear();
            for i in 0..nd {
                let ui = i as usize;
                if self.dwf_inflows_.node_idx[ui] < 0 {
                    continue;
                }
                let param = if self.dwf_inflows_.is_flow[ui] != 0 {
                    -1
                } else if self.dwf_inflows_.msx_idx[ui] >= 0 {
                    K_MSX_BASE + self.dwf_inflows_.msx_idx[ui] as i64
                } else if self.dwf_inflows_.pollut_idx[ui] >= 0 {
                    self.dwf_inflows_.pollut_idx[ui] as i64
                } else {
                    K_NONE
                };
                supersede(&mut last, &mut self.dwf_inflows_.node_idx, i, param);
            }
        }
    }
}
