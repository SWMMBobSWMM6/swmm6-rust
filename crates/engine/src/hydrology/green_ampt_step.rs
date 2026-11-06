// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/Infiltration.cpp (HydroCouple/
// openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed
// under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The Green-Ampt rate (grnampt_getInfil and its file-local grnampt_getF2,
// grnampt_getSatInfil, grnampt_getUnsatInfil), which the routing step reaches through storage
// exfiltration (ExfilSolver::computeAll). Op for op as upstream (legacy infil.c).
//
// MERGE NOTE: the runoff step needs the same functions; if another branch translates them
// into hydrology/infiltration.rs, keep one copy.

use crate::core::cmath::{cmax, cmin};
use crate::data::lid_node::GreenAmptState;
use crate::hydrology::infiltration::InfilModel;

/// legacy ZERO (consts.h)
const TINY: f64 = 1.0e-10;

/// grnampt_getF2: cumulative saturated infiltration after `ts` seconds.
fn grnampt_get_f2(f1: f64, c1: f64, ks: f64, ts: f64) -> f64 {
    let mut f2 = f1;
    let f2min = f1 + ks * ts;
    if c1 == 0.0 {
        return f2min;
    }
    if ts < 10.0 && f1 > 0.01 * c1 {
        f2 = f1 + ks * (1.0 + c1 / f1) * ts;
        return cmax(f2, f2min);
    }
    let c2 = c1 * (f1 + c1).ln() - ks * ts; // LIBM log
    for _i in 1..=20 {
        let df2 = (f2 - f1 - c1 * (f2 + c1).ln() + c2) / (1.0 - c1 / (f2 + c1)); // LIBM log
        if df2.abs() < 0.00001 {
            return cmax(f2, f2min);
        }
        f2 -= df2;
    }
    f2min
}

/// grnampt_getSatInfil
#[allow(clippy::too_many_arguments)]
fn grnampt_get_sat_infil(s: &mut GreenAmptState, tstep: f64, irate: f64, depth: f64, ks: f64, lu: f64, fumax: f64, recovery_factor: f64) -> f64 {
    let ia = irate + depth / tstep;
    if ia < TINY {
        return 0.0;
    }
    s.T = 5400.0 / lu / recovery_factor;
    let c1 = (s.S + depth) * s.IMD;
    let f2 = grnampt_get_f2(s.F, c1, ks, tstep);
    let mut d_f = f2 - s.F;
    if d_f > ia * tstep {
        d_f = ia * tstep;
        s.saturated = false;
    }
    s.F += d_f;
    s.Fu += d_f;
    s.Fu = cmin(s.Fu, fumax);
    d_f / tstep
}

/// grnampt_getUnsatInfil
#[allow(clippy::too_many_arguments)]
fn grnampt_get_unsat_infil(s: &mut GreenAmptState, tstep: f64, irate: f64, depth: f64, model_type: InfilModel, ks: f64, lu: f64, fumax: f64, recovery_factor: f64) -> f64 {
    let mut ia = irate + depth / tstep;
    if ia < TINY {
        ia = 0.0;
    }
    if ia == 0.0 {
        if s.Fu <= 0.0 {
            return 0.0;
        }
        let kr = lu / 90000.0 * recovery_factor;
        let d_f = kr * fumax * tstep;
        s.F -= d_f;
        s.Fu -= d_f;
        if s.Fu <= 0.0 {
            s.Fu = 0.0;
            s.F = 0.0;
            s.IMD = s.IMDmax;
            return 0.0;
        }
        if s.T <= 0.0 {
            s.IMD = (fumax - s.Fu) / lu;
            s.F = 0.0;
        }
        return 0.0;
    }
    if ia <= ks {
        let d_f = ia * tstep;
        s.F += d_f;
        s.Fu += d_f;
        s.Fu = cmin(s.Fu, fumax);
        if model_type == InfilModel::GreenAmpt && s.T <= 0.0 {
            s.IMD = (fumax - s.Fu) / lu;
            s.F = 0.0;
        }
        return ia;
    }
    s.T = 5400.0 / lu / recovery_factor;
    let fs = ks * (s.S + depth) * s.IMD / (ia - ks);
    if s.F > fs {
        s.saturated = true;
        return grnampt_get_sat_infil(s, tstep, irate, depth, ks, lu, fumax, recovery_factor);
    }
    if s.F + ia * tstep < fs {
        let d_f = ia * tstep;
        s.F += d_f;
        s.Fu += d_f;
        s.Fu = cmin(s.Fu, fumax);
        return ia;
    }
    let mut ts = tstep - (fs - s.F) / ia;
    if ts <= 0.0 {
        ts = 0.0;
    }
    let c1 = (s.S + depth) * s.IMD;
    let mut f2 = grnampt_get_f2(fs, c1, ks, ts);
    if f2 > fs + ia * ts {
        f2 = fs + ia * ts;
    }
    let d_f = f2 - s.F;
    s.F = f2;
    s.Fu += d_f;
    s.Fu = cmin(s.Fu, fumax);
    s.saturated = true;
    d_f / tstep
}

/// grnampt_getInfil: the Green-Ampt infiltration rate (ft/s) over `dt` seconds.
pub fn grnampt_get_infil(state: &mut GreenAmptState, precip: f64, depth: f64, dt: f64, model_type: InfilModel, infil_factor: f64, recovery_factor: f64) -> f64 {
    let sqrt_if = infil_factor.sqrt();
    let fumax = state.IMDmax * state.Lu * sqrt_if;
    let ks = state.Ks * infil_factor;
    let lu = state.Lu * sqrt_if;
    state.T -= dt;
    if state.saturated {
        return grnampt_get_sat_infil(state, dt, precip, depth, ks, lu, fumax, recovery_factor);
    }
    grnampt_get_unsat_infil(state, dt, precip, depth, model_type, ks, lu, fumax, recovery_factor)
}
