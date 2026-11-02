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
//   * NOT YET TRANSLATED (run-time step functions, phase 1c): infiltrationFactors, rate,
//     rateFeet, evaluate, pack, unpack, checkpoint, validateCheckpoint.

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
