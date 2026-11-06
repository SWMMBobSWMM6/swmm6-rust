// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/surface/InfilBank.hpp and
// InfilBank.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// One infiltration state owner for subcatchments (and, in 2D, mesh cells): init, clear,
// setMethod, setOwner/owner and cumulative, which RunoffSolver::init uses.
//
// Differences of form:
//   * methods_ holds the raw int of each InfilModel (`std::vector<InfilModel>` is a vector
//     of ints): setMethod stores whatever `static_cast<InfilModel>(int)` gave, an unknown
//     value included, and dispatches nothing for it, as upstream's switch does.
//   * The members are public so the gate dumper (an example crate) can read them;
//     upstream's are private.
//   * infiltrationFactors, rate, rateFeet and evaluate were translated on 5 November 2026
//     for the runoff step. NOT YET TRANSLATED: checkpoint, validateCheckpoint (engine
//     checkpoints).

use crate::core::cmath::cmax;
use crate::core::context::SimulationContext;
use crate::core::datetime;
use crate::core::options::SimulationOptions;
use crate::core::unit_conversion as ucf;
use crate::hydrology::infiltration::{self as infil, CurveNumState, GreenAmptState, HortonState, InfilModel};

/// InfilBank::Owner.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum Owner {
    Bank = 0,
    External = 1,
}

/// InfilBank::Factors.
#[derive(Debug, Clone, Copy)]
pub struct Factors {
    pub infiltration: f64,
    pub recovery: f64,
}

impl Default for Factors {
    fn default() -> Self {
        Factors { infiltration: 1.0, recovery: 1.0 }
    }
}

#[derive(Debug, Clone, Default)]
pub struct InfilBank {
    /// InfilModel as int, per entry
    pub methods_: Vec<i32>,
    pub owners_: Vec<Owner>,
    pub horton_: Vec<HortonState>,
    pub grnampt_: Vec<GreenAmptState>,
    pub curvenum_: Vec<CurveNumState>,
    pub constant_: Vec<f64>,
    pub cumulative_: Vec<f64>,
}

impl InfilBank {
    /// Clear all ownership/history and allocate n default method entries (Horton, BANK,
    /// zero constants and cumulative); per-model state arrays stay empty until set_method.
    pub fn init(&mut self, n: i32) {
        self.clear();
        let un = n as usize;
        self.methods_.clear();
        self.methods_.resize(un, InfilModel::Horton as i32);
        self.owners_.clear();
        self.owners_.resize(un, Owner::Bank);
        self.constant_.clear();
        self.constant_.resize(un, 0.0);
        self.cumulative_.clear();
        self.cumulative_.resize(un, 0.0);
    }

    pub fn clear(&mut self) {
        self.methods_.clear();
        self.owners_.clear();
        self.horton_.clear();
        self.grnampt_.clear();
        self.curvenum_.clear();
        self.constant_.clear();
        self.cumulative_.clear();
    }

    /// Select a method and convert its authored parameters for entry `e`. A per-model
    /// state bank is sized to the whole entry count on first use.
    pub fn set_method(&mut self, e: i32, method: i32, p: &[f64; 5], opts: &SimulationOptions) {
        let ue = e as usize;
        self.methods_[ue] = method;
        if method == InfilModel::Horton as i32 || method == InfilModel::ModHorton as i32 {
            if self.horton_.is_empty() {
                self.horton_.resize(self.methods_.len(), HortonState::default());
            }
            infil::horton_init(&mut self.horton_[ue], p[0], p[1], p[2], p[3], p[4], opts);
        } else if method == InfilModel::GreenAmpt as i32 || method == InfilModel::ModGreenAmpt as i32 {
            if self.grnampt_.is_empty() {
                self.grnampt_.resize(self.methods_.len(), GreenAmptState::default());
            }
            infil::grnampt_init(&mut self.grnampt_[ue], p[0], p[1], p[2], opts);
        } else if method == InfilModel::CurveNum as i32 {
            if self.curvenum_.is_empty() {
                self.curvenum_.resize(self.methods_.len(), CurveNumState::default());
            }
            infil::curvenum_init(&mut self.curvenum_[ue], p[0], p[2]);
        } else if method == InfilModel::Constant as i32 {
            // Dedicated stateless rate, preserving the original init operation order.
            self.constant_[ue] = p[0] / ucf::ucf(ucf::RAINFALL, opts);
        }
    }

    pub fn set_owner(&mut self, e: i32, owner: Owner) {
        self.owners_[e as usize] = owner;
    }

    pub fn owner(&self, e: i32) -> Owner {
        self.owners_[e as usize]
    }

    /// Per-entry integrated infiltration depth in metres.
    pub fn cumulative(&self) -> &Vec<f64> {
        &self.cumulative_
    }
}

// InfilBank::pack / unpack (translated 5 November 2026 for the legacy hot-start runoff
// section, SWMMEngine::applyLegacyRunoffSection).
impl InfilBank {
    /// pack(i, model, state): the method code and its six state slots (zeros for an
    /// external owner, a CONSTANT method or an index out of range, where model is 0).
    pub fn pack(&self, i: i32, model: &mut i32, state: &mut [f64; 6]) {
        *state = [0.0; 6];
        if i < 0 || i as usize >= self.methods_.len() {
            *model = 0;
            return;
        }
        let ui = i as usize;
        *model = self.methods_[ui];
        if self.owners_[ui] == Owner::External {
            return;
        }
        match self.methods_[ui] {
            m if m == InfilModel::Horton as i32 || m == InfilModel::ModHorton as i32 => {
                let h = &self.horton_[ui];
                state[0] = h.tp;
                state[1] = h.Fe;
                state[2] = h.Fmh;
            }
            m if m == InfilModel::GreenAmpt as i32 || m == InfilModel::ModGreenAmpt as i32 => {
                let g = &self.grnampt_[ui];
                state[0] = g.IMD;
                state[1] = g.F;
                state[2] = g.Fu;
                state[3] = g.T;
                state[4] = if g.saturated { 1.0 } else { 0.0 };
            }
            m if m == InfilModel::CurveNum as i32 => {
                let c = &self.curvenum_[ui];
                state[0] = c.S;
                state[1] = c.Se;
                state[2] = c.P;
                state[3] = c.F;
                state[4] = c.f;
                state[5] = c.T;
            }
            _ => {}
        }
    }

    /// unpack(i, model, state): restore the six slots when `model` is the entry's method.
    pub fn unpack(&mut self, i: i32, model: i32, state: &[f64; 6]) {
        if i < 0 || i as usize >= self.methods_.len() {
            return;
        }
        let ui = i as usize;
        if self.owners_[ui] == Owner::External {
            return;
        }
        if model != self.methods_[ui] {
            return;
        }
        match self.methods_[ui] {
            m if m == InfilModel::Horton as i32 || m == InfilModel::ModHorton as i32 => {
                let h = &mut self.horton_[ui];
                h.tp = state[0];
                h.Fe = state[1];
                h.Fmh = state[2];
            }
            m if m == InfilModel::GreenAmpt as i32 || m == InfilModel::ModGreenAmpt as i32 => {
                let g = &mut self.grnampt_[ui];
                g.IMD = state[0];
                g.F = state[1];
                g.Fu = state[2];
                g.T = state[3];
                g.saturated = state[4] != 0.0;
            }
            m if m == InfilModel::CurveNum as i32 => {
                let c = &mut self.curvenum_[ui];
                c.S = state[0];
                c.Se = state[1];
                c.P = state[2];
                c.F = state[3];
                c.f = state[4];
                c.T = state[5];
            }
            _ => {}
        }
    }
}

/// infiltrationFactors: the global monthly factors at a surface clock (elapsed seconds).
pub fn infiltration_factors(ctx: &SimulationContext, elapsed_s: f64) -> Factors {
    let month = datetime::month_of_year(datetime::add_seconds(ctx.options.start_date, elapsed_s)) - 1;
    let mut result = Factors { infiltration: ctx.adjust_hydcon[month as usize], recovery: 1.0 };
    let pattern = ctx.climate_state.recovery_pat_index;
    if pattern >= 0 && (pattern as usize) < ctx.patterns.factors.len() {
        let values = &ctx.patterns.factors[pattern as usize];
        if (month as usize) < values.len() {
            result.recovery = values[month as usize];
        }
    }
    result
}

// InfilBank::rateFeet / rate / evaluate (translated 5 November 2026 for the runoff step).
impl InfilBank {
    /// rateFeet: advance one bank-owned kernel in feet; zero for an EXTERNAL entry or dt <= 0.
    pub fn rate_feet(&mut self, e: i32, precip: f64, runon: f64, depth: f64, dt: f64, factors: Factors) -> f64 {
        let ue = e as usize;
        if self.owners_[ue] == Owner::External || dt <= 0.0 {
            return 0.0;
        }
        let f = self.evaluate(e, precip, runon, depth, dt, factors);
        self.cumulative_[ue] += cmax(0.0, f * 0.3048) * dt;
        f
    }

    /// rate: the SI boundary (m/s in and out).
    pub fn rate(&mut self, e: i32, precip: f64, runon: f64, depth: f64, dt: f64, factors: Factors) -> f64 {
        const FEET_PER_METRE: f64 = 3.280839895013123;
        cmax(0.0, self.rate_feet(e, precip * FEET_PER_METRE, runon * FEET_PER_METRE, depth * FEET_PER_METRE, dt, factors) * 0.3048)
    }

    fn evaluate(&mut self, i: i32, precip: f64, runon: f64, depth: f64, dt: f64, factors: Factors) -> f64 {
        let ui = i as usize;
        let mut infil = 0.0;
        let im = self.methods_[ui];
        if im == InfilModel::Horton as i32 || im == InfilModel::ModHorton as i32 {
            let hs = &mut self.horton_[ui];
            let (save_f0, save_fmin, save_regen) = (hs.f0, hs.fmin, hs.regen);
            hs.f0 *= factors.infiltration;
            hs.fmin *= factors.infiltration;
            hs.regen *= factors.recovery;
            infil = if im == InfilModel::Horton as i32 { infil::horton_get_infil(hs, precip + runon, depth, dt) } else { infil::mod_horton_get_infil(hs, precip + runon, depth, dt) };
            hs.f0 = save_f0;
            hs.fmin = save_fmin;
            hs.regen = save_regen;
        } else if im == InfilModel::GreenAmpt as i32 || im == InfilModel::ModGreenAmpt as i32 {
            infil = infil::grnampt_get_infil(&mut self.grnampt_[ui], precip + runon, depth, dt, im, factors.infiltration, factors.recovery);
        } else if im == InfilModel::CurveNum as i32 {
            // CN treats runon as ponded depth, not as a rainfall rate.
            let cn_depth = depth + runon * dt;
            infil = infil::curvenum_get_infil(&mut self.curvenum_[ui], precip, cn_depth, dt, factors.recovery);
        } else if im == InfilModel::Constant as i32 {
            infil = infil::constant_get_infil(self.constant_[ui] * factors.infiltration, precip + runon, depth, dt);
        }
        infil
    }
}
