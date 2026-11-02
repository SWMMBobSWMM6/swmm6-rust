// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/HydClosureKernels.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Portable scalar closures shared by the host solvers and the device backend: FHWA HEC-5
// culvert inlet control and force-main friction (upstream namespace openswmm::hydkernels).
//
// Differences in form, never in result:
//   * OPENSWMM_KERNEL_FN (plain `inline` on host builds) is dropped; these are ordinary
//     functions.
//   * culvertInflow's `double& dqdh` out-parameter is `&mut f64`.
//   * Upstream pins kViscosity/kGravity to constants::VISCOS/GRAVITY with static_asserts in
//     ForceMain.cpp; the same check is the const assertion at the end of this file.

#![allow(clippy::needless_late_init)]

use crate::core::constants;

// ===========================================================================
// Force-main friction
// ===========================================================================

/// Kinematic viscosity of water, ft²/s; must equal constants::VISCOS.
pub const K_VISCOSITY: f64 = 1.1e-5;

/// Gravity in internal units, matching constants::GRAVITY.
pub const K_GRAVITY: f64 = 32.2;

/// Hazen-Williams friction slope: (|v| / (1.318 C R^0.63))^(1/0.54); 0 when C or R <= 0.
pub fn fric_slope_hw(velocity: f64, hyd_rad: f64, c_hw: f64) -> f64 {
    if c_hw <= 0.0 || hyd_rad <= 0.0 {
        return 0.0;
    }
    let v_abs = if velocity < 0.0 { -velocity } else { velocity };
    (v_abs / (1.318 * c_hw * hyd_rad.powf(0.63))).powf(1.0 / 0.54) // LIBM (both)
}

/// Darcy-Weisbach friction slope with a laminar or Swamee-Jain friction factor.
pub fn fric_slope_dw(velocity: f64, hyd_rad: f64, roughness: f64) -> f64 {
    if hyd_rad <= 0.0 {
        return 0.0;
    }
    let v_abs = if velocity < 0.0 { -velocity } else { velocity };
    let diameter = 4.0 * hyd_rad;

    let re = v_abs * diameter / K_VISCOSITY;
    if re <= 0.0 {
        return 0.0;
    }

    let f: f64;
    if re <= 2000.0 {
        f = 64.0 / re; // laminar
    } else {
        let e_over_d = roughness / diameter;
        let arg = e_over_d / 3.7 + 5.74 / re.powf(0.9); // LIBM
        if arg <= 0.0 {
            return 0.0;
        }
        let logarg = arg.log10(); // LIBM
        f = 0.25 / (logarg * logarg);
    }
    f * v_abs * v_abs / (2.0 * K_GRAVITY * diameter)
}

// ===========================================================================
// Culvert inlet control (FHWA HEC-5)
// ===========================================================================

/// Inlet-control curve coefficients for one culvert type code.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[allow(non_snake_case)]
pub struct CulvertCurve {
    pub K: f64,
    pub M: f64,
    pub C: f64,
    pub Y: f64,
}

/// Limit a proposed culvert flow by the three-regime inlet-control curve. `dqdh` is set to
/// 0 first, then to the curve's derivative wherever the curve is evaluated.
#[allow(clippy::too_many_arguments, non_snake_case)]
pub fn culvert_inflow(
    q_proposed: f64,
    head: f64,
    y_full: f64,
    a_full: f64,
    slope: f64,
    cc: &CulvertCurve,
    mitered: bool,
    dqdh: &mut f64,
) -> f64 {
    *dqdh = 0.0;
    if head <= 0.0 || y_full <= 0.0 || a_full <= 0.0 {
        return q_proposed;
    }
    if cc.K == 0.0 {
        return q_proposed;
    }

    let AD = a_full * y_full.sqrt();
    let y_norm = head / y_full;
    let scf = if mitered { -7.0 * slope } else { 0.5 * slope };

    let y1_norm = 0.95;
    let mut y2_norm = 16.0 * cc.C + cc.Y - scf;
    if y2_norm < y1_norm {
        y2_norm = y1_norm + 0.01;
    }

    let q_inlet: f64;
    if y_norm <= y1_norm {
        let arg = y_norm / cc.K;
        if arg <= 0.0 {
            return q_proposed;
        }
        q_inlet = AD * arg.powf(1.0 / cc.M); // LIBM
        *dqdh = q_inlet / (head * cc.M);
    } else if y_norm >= y2_norm {
        let arg = (y_norm - cc.Y + scf) / cc.C;
        if arg <= 0.0 {
            return q_proposed;
        }
        q_inlet = AD * arg.sqrt();
        *dqdh = 0.5 * q_inlet / (arg * y_full * cc.C);
    } else {
        let arg1 = y1_norm / cc.K;
        let q1 = AD * arg1.powf(1.0 / cc.M); // LIBM
        let arg2 = (y2_norm - cc.Y + scf) / cc.C;
        let q2 = if arg2 > 0.0 { AD * arg2.sqrt() } else { q1 };

        let frac = (y_norm - y1_norm) / (y2_norm - y1_norm);
        q_inlet = q1 + frac * (q2 - q1);
        *dqdh = (q2 - q1) / ((y2_norm - y1_norm) * y_full);
    }
    if q_inlet < q_proposed {
        q_inlet
    } else {
        q_proposed
    }
}

// static_assert(hydkernels::kViscosity == constants::VISCOS) and
// static_assert(hydkernels::kGravity == constants::GRAVITY), from ForceMain.cpp.
const _: () = assert!(K_VISCOSITY == constants::VISCOS);
const _: () = assert!(K_GRAVITY == constants::GRAVITY);
