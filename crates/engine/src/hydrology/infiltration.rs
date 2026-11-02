// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/Infiltration.hpp and Infiltration.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The infiltration model enum, the three per-entry state types and their init functions
// (horton_init, grnampt_init, curvenum_init), which RunoffSolver::init (through
// InfilBank::setMethod) and LIDSolver::init call. CONSTANT has no state of its own; its
// rate is set by InfilBank::setMethod.
//
// Differences of form:
//   * GreenAmptState was first translated in data/lid_node.rs (LidNodeState holds one);
//     it is re-exported here so there is one type, as upstream has one.
//   * NOT YET TRANSLATED (run-time step functions, phase 1c): horton_getInfil,
//     modHorton_getInfil, grnampt_getInfil (and its static helpers grnampt_getF2,
//     grnampt_getUnsatInfil, grnampt_getSatInfil), curvenum_getInfil, constant_getInfil.

#![allow(non_snake_case)]
use crate::core::cmath::{cmax, cmin};
use crate::core::options::SimulationOptions;
use crate::core::unit_conversion as ucf;

pub use crate::data::lid_node::GreenAmptState;

/// Infiltration model type (`enum class InfilModel : int`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum InfilModel {
    Horton = 0,
    /// Modified Horton: linear decay fp = f0 - kd * Fe
    ModHorton = 1,
    GreenAmpt = 2,
    /// Modified Green-Ampt: does not reset F between events
    ModGreenAmpt = 3,
    CurveNum = 4,
    /// Constant rate, capacity-bounded. 2D-only (no legacy [INFILTRATION] token).
    Constant = 5,
}

/// Horton infiltration state (per subcatchment). Every initializer is zero.
#[derive(Debug, Clone, Copy, Default)]
pub struct HortonState {
    /// Initial infiltration rate (ft/sec)
    pub f0: f64,
    /// Minimum (ultimate) rate (ft/sec)
    pub fmin: f64,
    /// Decay constant k (1/sec)
    pub decay: f64,
    /// Regeneration constant kr (1/sec)
    pub regen: f64,
    /// Max cumulative infiltration (ft)
    pub Fmax: f64,
    /// Cumulative infiltration time (sec)
    pub tp: f64,
    /// Cumulative excess infiltration (ft)
    pub Fe: f64,
    /// Cumulative infiltration for Modified Horton (ft)
    pub Fmh: f64,
}

/// SCS Curve Number state (per subcatchment). Every initializer is zero.
#[derive(Debug, Clone, Copy, Default)]
pub struct CurveNumState {
    /// Max retention S = (1000/CN - 10)/12 (ft)
    pub Smax: f64,
    /// Current retention (ft)
    pub S: f64,
    /// Effective retention at event start (ft)
    pub Se: f64,
    /// Cumulative precipitation (ft)
    pub P: f64,
    /// Cumulative infiltration (ft)
    pub F: f64,
    /// Previous infiltration rate (ft/sec)
    pub f: f64,
    /// Regeneration rate (1/sec)
    pub regen: f64,
    /// Time since last rainfall (sec)
    pub T: f64,
    /// Inter-event time (sec)
    pub Tmax: f64,
}

/// Infiltration.cpp's file constant (legacy ZERO = 1.E-10); used by the step functions.
#[allow(dead_code)]
const TINY: f64 = 1.0e-10;

// ============================================================================
// Horton
// ============================================================================

/// Convert Horton parameters from project units to internal ft and seconds. Fmh is not
/// touched, as upstream leaves it.
pub fn horton_init(state: &mut HortonState, f0: f64, fmin: f64, decay: f64, regen: f64, Fmax: f64, opts: &SimulationOptions) {
    // Convert from project rain units (in/hr or mm/hr) to ft/sec via UCF
    let ucf_rain = ucf::ucf(ucf::RAINFALL, opts);
    let ucf_depth = ucf::ucf(ucf::RAINDEPTH, opts);
    state.f0 = f0 / ucf_rain;
    state.fmin = fmin / ucf_rain;
    state.decay = decay / 3600.0; // 1/hr -> 1/sec
    // regen parameter is drying time in days (legacy infil.c line 350):
    // regen = -log(1.0 - 0.98) / p[3] / SECperDAY
    state.regen = if regen > 0.0 { -(1.0f64 - 0.98).ln() / regen / ucf::SEC_PER_DAY } else { 0.0 }; // LIBM
    state.Fmax = if Fmax > 0.0 { Fmax / ucf_depth } else { 0.0 }; // in/mm -> ft (0 = unlimited)
    state.tp = 0.0;
    state.Fe = 0.0;
}

// ============================================================================
// Green-Ampt
// ============================================================================

pub fn grnampt_init(state: &mut GreenAmptState, S: f64, Ks: f64, IMD: f64, opts: &SimulationOptions) {
    let ucf_rain = ucf::ucf(ucf::RAINFALL, opts);
    let ucf_depth = ucf::ucf(ucf::RAINDEPTH, opts);
    state.S = S / ucf_depth; // in/mm -> ft
    state.Ks = Ks / ucf_rain; // in/hr or mm/hr -> ft/sec
    state.IMDmax = IMD;
    state.IMD = IMD;
    state.F = 0.0;
    // Legacy grnampt_setParams (infil.c) derives the upper-zone depth from Ks in in/hr
    // whatever the unit system: ksat = Ks * 12. * 3600.; Lu = 4.0 * sqrt(ksat) / 12.
    let ksat = state.Ks * 12.0 * 3600.0;
    state.Lu = 4.0 * ksat.sqrt() / 12.0;
    state.Fumax = state.IMDmax * state.Lu;
    state.Fu = 0.0; // Legacy: starts dry (empty upper zone = max storage)
    state.saturated = false;
}

// ============================================================================
// SCS Curve Number
// ============================================================================

pub fn curvenum_init(state: &mut CurveNumState, CN: f64, regen_days: f64) {
    // SCS formula: S in inches, then convert to feet (always /12)
    let mut CN = cmax(CN, 10.0);
    CN = cmin(CN, 99.0);
    state.Smax = (1000.0 / CN - 10.0) / 12.0;
    state.S = state.Smax;
    state.Se = state.Smax;
    state.P = 0.0;
    state.F = 0.0;
    state.f = 0.0;
    // Legacy: regen = 1.0 / (dryingTime_days * SECperDAY)
    const SEC_PER_DAY: f64 = 86400.0;
    state.regen = if regen_days > 0.0 { 1.0 / (regen_days * SEC_PER_DAY) } else { 0.0 };
    state.T = 0.0;
    state.Tmax = if state.regen > 0.0 { 0.06 / state.regen } else { 1.0e10 };
}
