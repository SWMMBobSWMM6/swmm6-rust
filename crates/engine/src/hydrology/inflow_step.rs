// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/Inflow.cpp and Inflow.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The routing step's external and dry-weather inflows (inflow.rs holds the solver's state
// and init): InflowSolver::computeAll and getPatternFactor.

use crate::core::cmath::{cmax, cmin};
use crate::core::constants;
use crate::core::context::SimulationContext;
use crate::core::datetime;
use crate::data::table::table_tseries_lookup_cursor;
use crate::hydrology::inflow::{ExtInflowKind, InflowSolver};

impl InflowSolver {
    /// getPatternFactor
    pub fn get_pattern_factor(&self, pat_idx: i32, month: i32, day: i32, hour: i32) -> f64 {
        if pat_idx < 0 || pat_idx >= self.patterns_.len() as i32 {
            return 1.0;
        }
        let pat = &self.patterns_[pat_idx as usize];
        match pat.type_ {
            0 => pat.factors[(month % 12) as usize],
            1 => pat.factors[(day % 7) as usize],
            2 => pat.factors[(hour % 24) as usize],
            3 => pat.factors[(hour % 24) as usize],
            _ => 1.0,
        }
    }

    /// Value of one [INFLOWS] row at the date: cf * (tsv + blv).
    fn row_value(&self, ctx: &mut SimulationContext, ui: usize, current_date: f64, month: i32, day: i32, hour: i32) -> f64 {
        let mut base = self.ext_inflows_.baseline[ui];
        let bp = self.ext_inflows_.base_pat_idx[ui];
        if bp >= 0 {
            base *= self.get_pattern_factor(bp, month, day, hour);
        }
        let mut ts_val = 0.0;
        let ts = self.ext_inflows_.ts_idx[ui];
        if ts >= 0 && ts < ctx.tables.count() as i32 {
            ts_val = table_tseries_lookup_cursor(&mut ctx.tables.tables[ts as usize], current_date);
            ts_val *= self.ext_inflows_.scale_factor[ui];
        }
        self.ext_inflows_.conv_factor[ui] * (ts_val + base)
    }

    /// Pattern factor of one DWF row: monthly * daily * (hourly | weekend).
    fn dwf_factor(&self, ui: usize, month: i32, day: i32, hour: i32) -> f64 {
        let mut factor = 1.0;
        let pm = self.dwf_inflows_.pat_monthly[ui];
        if pm >= 0 {
            factor *= self.get_pattern_factor(pm, month, day, hour);
        }
        let pd = self.dwf_inflows_.pat_daily[ui];
        if pd >= 0 {
            factor *= self.get_pattern_factor(pd, month, day, hour);
        }
        let ph = self.dwf_inflows_.pat_hourly[ui];
        let pw = self.dwf_inflows_.pat_weekend[ui];
        if pw >= 0 {
            if day == 0 || day == 6 {
                factor *= self.get_pattern_factor(pw, month, day, hour);
            } else if ph >= 0 {
                factor *= self.get_pattern_factor(ph, month, day, hour);
            }
        } else if ph >= 0 {
            factor *= self.get_pattern_factor(ph, month, day, hour);
        }
        factor
    }

    /// InflowSolver::computeAll: this step's external and dry-weather inflows into the
    /// decomposed node buffers (and the pollutant / species mass rates).
    pub fn compute_all(&mut self, ctx: &mut SimulationContext, current_date: f64, _dt: f64) {
        let (h_tmp, _m, _s) = datetime::decode_time(current_date);
        let hour = h_tmp;
        let total_days = current_date.floor() as i32;
        let day = (total_days + constants::DATE_DELTA) % 7;
        let month = datetime::month_of_year(current_date) - 1;

        let np = ctx.n_pollutants();

        // Pass 1 — FLOW rows
        for i in 0..self.ext_inflows_.count {
            let ui = i as usize;
            if self.ext_inflows_.kind[ui] != ExtInflowKind::Flow as i32 {
                continue;
            }
            let mut q = self.row_value(ctx, ui, current_date, month, day, hour);
            if q.abs() < constants::FLOW_TOL {
                q = 0.0;
            }
            let ni = self.ext_inflows_.node_idx[ui];
            if ni >= 0 && ni < ctx.nodes.ext_inflow.len() as i32 {
                ctx.nodes.ext_inflow[ni as usize] += q;
            }
        }

        // Pass 2 — pollutant rows
        if np > 0 && !ctx.nodes.ext_qual_mass.is_empty() {
            for i in 0..self.ext_inflows_.count {
                let ui = i as usize;
                let kind = self.ext_inflows_.kind[ui];
                if kind != ExtInflowKind::Concen as i32 && kind != ExtInflowKind::Mass as i32 {
                    continue;
                }
                let p = self.ext_inflows_.pollut_idx[ui];
                let ni = self.ext_inflows_.node_idx[ui];
                if p < 0 || p >= np || ni < 0 || ni >= ctx.n_nodes() {
                    continue;
                }
                let mut w = self.row_value(ctx, ui, current_date, month, day, hour);
                if kind == ExtInflowKind::Concen as i32 {
                    w *= ctx.nodes.ext_inflow[ni as usize];
                }
                let idx = (ni as usize) * (np as usize) + (p as usize);
                if idx < ctx.nodes.ext_qual_mass.len() {
                    ctx.nodes.ext_qual_mass[idx] += w;
                }
            }
        }

        // Pass 2b — reaction-species rows into msx_ext_mass_in
        {
            let nsp = if ctx.reactions.configured { ctx.reactions.n_species() } else { 0 };
            if nsp > 0 && ctx.n_nodes() > 0 {
                let mut any_msx = false;
                let mut i = 0;
                while i < self.ext_inflows_.count && !any_msx {
                    let kind = self.ext_inflows_.kind[i as usize];
                    any_msx = kind == ExtInflowKind::MsxConcen as i32 || kind == ExtInflowKind::MsxMass as i32;
                    i += 1;
                }
                let want = (ctx.n_nodes() as usize) * (nsp as usize);
                if any_msx || ctx.reactions.msx_ext_mass_in.len() == want {
                    if ctx.reactions.msx_ext_mass_in.len() != want {
                        ctx.reactions.msx_ext_mass_in.clear();
                        ctx.reactions.msx_ext_mass_in.resize(want, 0.0);
                    } else {
                        ctx.reactions.msx_ext_mass_in.iter_mut().for_each(|x| *x = 0.0);
                    }
                }
                if any_msx {
                    for i in 0..self.ext_inflows_.count {
                        let ui = i as usize;
                        let kind = self.ext_inflows_.kind[ui];
                        if kind != ExtInflowKind::MsxConcen as i32 && kind != ExtInflowKind::MsxMass as i32 {
                            continue;
                        }
                        let m = self.ext_inflows_.pollut_idx[ui];
                        let ni = self.ext_inflows_.node_idx[ui];
                        if m < 0 || m >= nsp || ni < 0 || ni >= ctx.n_nodes() {
                            continue;
                        }
                        if (m as usize) < ctx.reactions.species_is_wall.len() && ctx.reactions.species_is_wall[m as usize] != 0 {
                            continue;
                        }
                        let mut w = self.row_value(ctx, ui, current_date, month, day, hour);
                        if kind == ExtInflowKind::MsxConcen as i32 {
                            w *= ctx.nodes.ext_inflow[ni as usize];
                        }
                        ctx.reactions.msx_ext_mass_in[(ni as usize) * (nsp as usize) + (m as usize)] += w;
                    }
                }
            }
        }

        // Pass 3 — __WATER_AGE__ rows
        if ctx.options.water_age && !ctx.water_age_state.node_ext_inflow_age.is_empty() {
            for i in 0..self.ext_inflows_.count {
                let ui = i as usize;
                if self.ext_inflows_.kind[ui] != ExtInflowKind::Age as i32 {
                    continue;
                }
                let ni = self.ext_inflows_.node_idx[ui];
                if ni >= 0 && (ni as usize) < ctx.water_age_state.node_ext_inflow_age.len() {
                    let v = self.row_value(ctx, ui, current_date, month, day, hour);
                    ctx.water_age_state.node_ext_inflow_age[ni as usize] = v;
                }
            }
        }

        // DWF Pass 1 — FLOW rows
        for i in 0..self.dwf_inflows_.count {
            let ui = i as usize;
            if self.dwf_inflows_.is_flow[ui] == 0 {
                continue;
            }
            let mut q = self.dwf_factor(ui, month, day, hour) * self.dwf_inflows_.avg_value[ui];
            if q.abs() < constants::FLOW_TOL {
                q = 0.0;
            }
            let ni = self.dwf_inflows_.node_idx[ui];
            if ni >= 0 && ni < ctx.nodes.dwf_inflow.len() as i32 {
                ctx.nodes.dwf_inflow[ni as usize] += q;
            }
        }

        // DWF Pass 2 — pollutant rows
        if np > 0 && !ctx.nodes.dwf_qual_mass.is_empty() {
            for i in 0..self.dwf_inflows_.count {
                let ui = i as usize;
                let p = self.dwf_inflows_.pollut_idx[ui];
                if p < 0 || p >= np {
                    continue;
                }
                let ni = self.dwf_inflows_.node_idx[ui];
                if ni < 0 || ni >= ctx.n_nodes() {
                    continue;
                }
                let q = ctx.nodes.dwf_inflow[ni as usize];
                if q <= 0.0 {
                    continue;
                }
                let mut w = q * (self.dwf_factor(ui, month, day, hour) * self.dwf_inflows_.avg_value[ui]);
                let c_glob = ctx.pollutants.c_dwf[p as usize];
                if c_glob > 0.0 {
                    w -= q * c_glob;
                }
                let idx = (ni as usize) * (np as usize) + (p as usize);
                if idx < ctx.nodes.dwf_qual_mass.len() {
                    ctx.nodes.dwf_qual_mass[idx] += w;
                }
            }
        }

        // DWF Pass 2b — species DWF rows
        {
            let nsp = if ctx.reactions.configured { ctx.reactions.n_species() } else { 0 };
            if nsp > 0 && ctx.n_nodes() > 0 {
                let mut any = false;
                let mut i = 0;
                while i < self.dwf_inflows_.count && !any {
                    any = self.dwf_inflows_.msx_idx[i as usize] >= 0;
                    i += 1;
                }
                if any {
                    let want = (ctx.n_nodes() as usize) * (nsp as usize);
                    if ctx.reactions.msx_ext_mass_in.len() != want {
                        ctx.reactions.msx_ext_mass_in.clear();
                        ctx.reactions.msx_ext_mass_in.resize(want, 0.0);
                    }
                    for i in 0..self.dwf_inflows_.count {
                        let ui = i as usize;
                        let m = self.dwf_inflows_.msx_idx[ui];
                        if m < 0 || m >= nsp {
                            continue;
                        }
                        if (m as usize) < ctx.reactions.species_is_wall.len() && ctx.reactions.species_is_wall[m as usize] != 0 {
                            continue;
                        }
                        let ni = self.dwf_inflows_.node_idx[ui];
                        if ni < 0 || ni >= ctx.n_nodes() {
                            continue;
                        }
                        let q = ctx.nodes.dwf_inflow[ni as usize];
                        if q <= 0.0 {
                            continue;
                        }
                        ctx.reactions.msx_ext_mass_in[(ni as usize) * (nsp as usize) + (m as usize)] += q * (self.dwf_factor(ui, month, day, hour) * self.dwf_inflows_.avg_value[ui]);
                    }
                }
            }
        }

        // Reaction sources
        let n_nodes = ctx.n_nodes() as usize;
        let start_date = ctx.options.start_date;
        let rx = &mut ctx.reactions;
        if rx.configured && !rx.sources.is_empty() {
            let ns = rx.n_species() as usize;
            let want = n_nodes * ns;
            if rx.msx_ext_mass_in.len() != want {
                rx.msx_ext_mass_in.clear();
                rx.msx_ext_mass_in.resize(want, 0.0);
            }
            let seconds = cmax(0.0, (current_date - start_date) * 86400.0);
            for si in 0..rx.sources.len() {
                let source = rx.sources[si].clone();
                if source.booster != 0 {
                    continue;
                }
                let mut factor = 1.0;
                if source.pattern >= 0 {
                    let factors = &rx.source_patterns[source.pattern as usize].factors;
                    let count = factors.len() as f64;
                    let period = rx.pattern_timestep * count;
                    let phase = if period.is_finite() { seconds % period } else { seconds }; // LIBM fmod
                    let epsilon = cmin(1.0e-6, rx.pattern_timestep * 1.0e-6);
                    let step = ((phase + epsilon) / rx.pattern_timestep).floor();
                    factor = factors[(step % count) as usize]; // LIBM fmod
                }
                let node = source.node as usize;
                let water = cmax(0.0, ctx.nodes.ext_inflow[node]) + cmax(0.0, ctx.nodes.dwf_inflow[node]);
                let rate = source.value * factor * (if source.mass { 1.0 / 28.317 } else { water });
                rx.msx_ext_mass_in[node * ns + source.species as usize] += rate;
            }
        }
    }
}
