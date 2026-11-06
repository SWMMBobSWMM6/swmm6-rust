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
//   * The run-time kernels (horton_getInfil, modHorton_getInfil, grnampt_getInfil with its
//     static helpers grnampt_getF2, grnampt_getUnsatInfil and grnampt_getSatInfil,
//     curvenum_getInfil, constant_getInfil) were translated on 5 November 2026 for the
//     runoff step, op for op.

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

// ============================================================================
// Run-time kernels (Infiltration.cpp)
// ============================================================================

/// horton_getInfil: the legacy infil.c Horton rate over one step (ft/s), the dry-period
/// regeneration included.
pub fn horton_get_infil(state: &mut HortonState, precip: f64, depth: f64, dt: f64) -> f64 {
    let ia = precip + depth / dt; // available water rate (ft/sec)
    let f0 = state.f0;
    let fmin = state.fmin;
    let df = f0 - fmin;
    let kd = state.decay;
    let Fmax = state.Fmax;
    let mut tp = state.tp;
    let kr = state.regen;

    // Special cases: no infiltration or constant infiltration
    if df < 0.0 || kd < 0.0 || kr < 0.0 {
        return 0.0;
    }
    if df == 0.0 || kd == 0.0 {
        let mut fp = f0;
        fp = cmin(fp, ia);
        return cmax(0.0, fp);
    }

    let mut fp = 0.0;
    let tlim = 16.0 / kd;

    if ia > TINY {
        // Compute average infiltration rate over timestep
        let t1 = tp + dt;
        let Fp;
        let mut F1;
        if tp >= tlim {
            Fp = fmin * tp + df / kd;
            F1 = Fp + fmin * dt;
        } else {
            Fp = fmin * tp + df / kd * (1.0 - (-kd * tp).exp()); // LIBM
            F1 = fmin * t1 + df / kd * (1.0 - (-kd * t1).exp()); // LIBM
        }
        fp = (F1 - Fp) / dt;
        // Legacy MAX(fp, fmin): fmin when fp is NaN
        fp = if fp >= fmin { fp } else { fmin };

        // Limit infiltration rate to available water
        fp = cmin(fp, ia);

        // Update cumulative time tp
        if t1 > tlim {
            tp = t1;
        } else if fp < ia {
            tp = t1;
        } else {
            // Infiltration limited by available water: Newton-Raphson for tp
            F1 = Fp + fp * dt;
            tp = tp + dt / 2.0;
            for _iter in 0..20 {
                let kt = if 60.0 <= kd * tp { 60.0 } else { kd * tp };
                let ex = (-kt).exp(); // LIBM
                let FF = fmin * tp + df / kd * (1.0 - ex) - F1;
                let FF1 = fmin + df * ex;
                let r = FF / FF1;
                tp -= r;
                if r.abs() <= 0.001 * dt {
                    break;
                }
            }
        }

        // Limit cumulative infiltration to Fmax
        if Fmax > 0.0 {
            if state.Fe + fp * dt > Fmax {
                fp = (Fmax - state.Fe) / dt;
            }
            fp = cmax(fp, 0.0);
            state.Fe += fp * dt;
        }
    }
    // Dry period with regeneration
    else if kr > 0.0 {
        let r = (-kr * dt).exp(); // LIBM
        tp = 1.0 - (-kd * tp).exp(); // LIBM
        tp = -(1.0 - r * tp).ln() / kd; // LIBM

        if Fmax > 0.0 {
            state.Fe = fmin * tp + (df / kd) * (1.0 - (-kd * tp).exp()); // LIBM
        }
    }

    state.tp = tp;
    fp
}

/// modHorton_getInfil: the modified Horton rate (linear decay fp = f0 - kd * Fe).
pub fn mod_horton_get_infil(state: &mut HortonState, precip: f64, depth: f64, dt: f64) -> f64 {
    let mut f = 0.0;
    let f0 = state.f0;
    let fmin = state.fmin;
    let df = f0 - fmin;
    let kd = state.decay;
    let kr = state.regen;

    if df < 0.0 || kd < 0.0 || kr < 0.0 {
        return 0.0;
    }
    if df == 0.0 || kd == 0.0 {
        let mut fp = f0;
        let fa = precip + depth / dt;
        fp = cmin(fp, fa);
        return cmax(0.0, fp);
    }

    let fa = precip + depth / dt;

    if fa > TINY {
        if state.Fmax > 0.0 && state.Fe >= state.Fmax {
            return 0.0;
        }
        let mut fp = f0 - kd * state.Fe;
        fp = cmax(fp, fmin);
        f = cmin(fa, fp);
        if state.Fmax > 0.0 {
            if state.Fmh + f * dt > state.Fmax {
                f = (state.Fmax - state.Fmh) / dt;
            }
            f = cmax(f, 0.0);
            state.Fmh += f * dt;
        }
        state.Fe += cmax(f - fmin, 0.0) * dt;
        if state.Fmax > 0.0 {
            state.Fe = cmin(state.Fe, state.Fmax);
        }
    } else if kr > 0.0 {
        let decay_factor = (-kr * dt).exp(); // LIBM
        state.Fe *= decay_factor;
        state.Fe = cmax(state.Fe, 0.0);
        state.Fmh *= decay_factor;
        state.Fmh = cmax(state.Fmh, 0.0);
    }

    f
}

/// grnampt_getF2: cumulative saturated infiltration after a duration (static upstream).
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

    let c2 = c1 * (f1 + c1).ln() - ks * ts; // LIBM
    for _i in 1..=20 {
        let df2 = (f2 - f1 - c1 * (f2 + c1).ln() + c2) / (1.0 - c1 / (f2 + c1)); // LIBM
        if df2.abs() < 0.00001 {
            return cmax(f2, f2min);
        }
        f2 -= df2;
    }
    f2min
}

/// grnampt_getSatInfil (static upstream).
#[allow(clippy::too_many_arguments)]
fn grnampt_get_sat_infil(s: &mut GreenAmptState, tstep: f64, irate: f64, depth: f64, ks: f64, lu: f64, Fumax: f64, recovery_factor: f64) -> f64 {
    let ia = irate + depth / tstep;
    if ia < TINY {
        return 0.0;
    }

    s.T = 5400.0 / lu / recovery_factor;

    let c1 = (s.S + depth) * s.IMD;
    let F2 = grnampt_get_f2(s.F, c1, ks, tstep);
    let mut dF = F2 - s.F;

    if dF > ia * tstep {
        dF = ia * tstep;
        s.saturated = false;
    }

    s.F += dF;
    s.Fu += dF;
    s.Fu = cmin(s.Fu, Fumax);
    dF / tstep
}

/// grnampt_getUnsatInfil (static upstream).
#[allow(clippy::too_many_arguments)]
fn grnampt_get_unsat_infil(s: &mut GreenAmptState, tstep: f64, irate: f64, depth: f64, model_type: i32, ks: f64, lu: f64, Fumax: f64, recovery_factor: f64) -> f64 {
    let mut ia = irate + depth / tstep;
    if ia < TINY {
        ia = 0.0;
    }

    if ia == 0.0 {
        if s.Fu <= 0.0 {
            return 0.0;
        }
        let kr = lu / 90000.0 * recovery_factor;
        let dF = kr * Fumax * tstep;
        s.F -= dF;
        s.Fu -= dF;
        if s.Fu <= 0.0 {
            s.Fu = 0.0;
            s.F = 0.0;
            s.IMD = s.IMDmax;
            return 0.0;
        }

        if s.T <= 0.0 {
            s.IMD = (Fumax - s.Fu) / lu;
            s.F = 0.0;
        }
        return 0.0;
    }

    if ia <= ks {
        let dF = ia * tstep;
        s.F += dF;
        s.Fu += dF;
        s.Fu = cmin(s.Fu, Fumax);
        if model_type == InfilModel::GreenAmpt as i32 && s.T <= 0.0 {
            s.IMD = (Fumax - s.Fu) / lu;
            s.F = 0.0;
        }
        return ia;
    }

    s.T = 5400.0 / lu / recovery_factor;

    let Fs = ks * (s.S + depth) * s.IMD / (ia - ks);

    if s.F > Fs {
        s.saturated = true;
        return grnampt_get_sat_infil(s, tstep, irate, depth, ks, lu, Fumax, recovery_factor);
    }

    if s.F + ia * tstep < Fs {
        let dF = ia * tstep;
        s.F += dF;
        s.Fu += dF;
        s.Fu = cmin(s.Fu, Fumax);
        return ia;
    }

    let mut ts = tstep - (Fs - s.F) / ia;
    if ts <= 0.0 {
        ts = 0.0;
    }

    let c1 = (s.S + depth) * s.IMD;
    let mut F2 = grnampt_get_f2(Fs, c1, ks, ts);
    if F2 > Fs + ia * ts {
        F2 = Fs + ia * ts;
    }

    let dF = F2 - s.F;
    s.F = F2;
    s.Fu += dF;
    s.Fu = cmin(s.Fu, Fumax);
    s.saturated = true;
    dF / tstep
}

/// grnampt_getInfil: the factor-scaled locals (legacy's per-call scaling), then dispatch on
/// the saturation flag. `model_type` is the InfilModel as int.
pub fn grnampt_get_infil(state: &mut GreenAmptState, precip: f64, depth: f64, dt: f64, model_type: i32, infil_factor: f64, recovery_factor: f64) -> f64 {
    let sqrt_if = infil_factor.sqrt();
    let Fumax = state.IMDmax * state.Lu * sqrt_if;
    let ks = state.Ks * infil_factor;
    let lu = state.Lu * sqrt_if;
    state.T -= dt;
    if state.saturated {
        return grnampt_get_sat_infil(state, dt, precip, depth, ks, lu, Fumax, recovery_factor);
    }
    grnampt_get_unsat_infil(state, dt, precip, depth, model_type, ks, lu, Fumax, recovery_factor)
}

/// curvenum_getInfil: legacy infil.c's curve number rate.
pub fn curvenum_get_infil(state: &mut CurveNumState, precip: f64, depth: f64, dt: f64, recovery_factor: f64) -> f64 {
    const MIN_TOTAL_DEPTH: f64 = 0.004167; // ft
    let fa = precip + depth / dt;
    let mut f1 = 0.0;

    if precip > TINY {
        if state.T >= state.Tmax {
            state.P = 0.0;
            state.F = 0.0;
            state.f = 0.0;
            state.Se = state.S;
        }
        state.T = 0.0;

        state.P += precip * dt;

        let F1 = state.P * (1.0 - state.P / (state.P + state.Se));

        f1 = (F1 - state.F) / dt;
        if f1 < 0.0 || state.S <= 0.0 {
            f1 = 0.0;
        }
    } else if depth > MIN_TOTAL_DEPTH && state.S > 0.0 {
        f1 = state.f;
        if f1 * dt > state.S {
            f1 = state.S / dt;
        }
    } else {
        state.T += dt;
    }

    if f1 > 0.0 {
        f1 = cmin(f1, fa);
        f1 = cmax(f1, 0.0);

        state.F += f1 * dt;

        if state.regen > 0.0 {
            state.S -= f1 * dt;
            state.S = cmax(state.S, 0.0);
        }
    } else {
        state.S += state.regen * state.Smax * dt * recovery_factor;
        state.S = cmin(state.S, state.Smax);
    }

    state.f = f1;
    f1
}

/// constant_getInfil: a capacity-bounded constant rate (2D only upstream).
pub fn constant_get_infil(rate_ftsec: f64, precip: f64, depth: f64, dt: f64) -> f64 {
    let ia = precip + depth / dt;
    cmax(0.0, cmin(rate_ftsec, ia))
}
