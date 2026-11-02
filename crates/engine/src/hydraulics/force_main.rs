// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/ForceMain.hpp and
// src/engine/hydraulics/ForceMain.cpp (HydroCouple/openswmm.engine, branch swmm6_rel,
// commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the Apache License,
// Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Force main friction: Hazen-Williams and Darcy-Weisbach (legacy forcmain.c). The two
// friction-slope laws live in hyd_closure_kernels; these are the public entry points.
//
// Differences in form, never in result:
//   * Names are snake_case (getFricSlope_HW -> get_fric_slope_hw, getEquivN ->
//     get_equiv_n, ...).
//   * FrictionModel is a closed Rust enum, so the `return n_raw;` / `return 0.0;` that
//     follow upstream's switches (reached only by an out-of-range model value) have no
//     counterpart.
//   * batchFricSlope's raw arrays are slices plus upstream's int count.

use crate::hydraulics::hyd_closure_kernels as hydkernels;

/// Kinematic viscosity @ 20C (ft2/sec)
pub const VISCOS: f64 = 1.1e-5;

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum FrictionModel {
    HazenWilliams = 0,
    DarcyWeisbach = 1,
}

/// Friction slope for Hazen-Williams.
pub fn get_fric_slope_hw(velocity: f64, hyd_rad: f64, c_hw: f64) -> f64 {
    hydkernels::fric_slope_hw(velocity, hyd_rad, c_hw)
}

/// Friction slope for Darcy-Weisbach.
pub fn get_fric_slope_dw(velocity: f64, hyd_rad: f64, roughness: f64) -> f64 {
    hydkernels::fric_slope_dw(velocity, hyd_rad, roughness)
}

// ============================================================================
// getFricFactor — Darcy-Weisbach friction factor (matching legacy)
// ============================================================================

/// Darcy-Weisbach friction factor f (legacy forcemain_getFricFactor, bit-exact).
pub fn get_fric_factor(r_bot: f64, hrad: f64, re: f64) -> f64 {
    let mut re = re;
    if hrad <= 0.0 {
        return 0.0;
    }
    // Legacy forcmain.c forcemain_getFricFactor, bit-exact (Swamee-Jain
    // approximation to Colebrook-White). The three arms — laminar, the
    // transitional 2000<Re<4000 blend, and turbulent — plus the Re<10 floor
    // and the conditional 5.74 term must all be present. Grouping matches
    // legacy: e/3.7/(4·hrad) and 0.25/f/f, left to right.
    let mut f: f64;
    if re < 10.0 {
        re = 10.0;
    }
    if re <= 2000.0 {
        f = 64.0 / re;
    } else if re < 4000.0 {
        f = get_fric_factor(r_bot, hrad, 4000.0);
        f = 0.032 + (f - 0.032) * (re - 2000.0) / 2000.0;
    } else {
        f = r_bot / 3.7 / (4.0 * hrad);
        if re < 1.0e10 {
            f += 5.74 / re.powf(0.9); // LIBM
        }
        f = f.log10(); // LIBM
        f = 0.25 / f / f;
    }
    f
}

// ============================================================================
// getEquivN — equivalent Manning's n for DW force main (Gap #22)
// ============================================================================

/// Equivalent Manning's n of a force main (legacy forcemain_getEquivN).
pub fn get_equiv_n(model: FrictionModel, r_bot: f64, y_full: f64, slope: f64, n_raw: f64) -> f64 {
    if slope <= 0.0 || y_full <= 0.0 {
        return n_raw;
    }
    match model {
        FrictionModel::HazenWilliams => {
            // legacy: 1.067 / rBot * pow(d/slope, 0.04)
            1.067 / r_bot * (y_full / slope).powf(0.04) // LIBM
        }
        FrictionModel::DarcyWeisbach => {
            // legacy: sqrt(f/185.0) * pow(d, 1/6)
            let hrad = y_full / 4.0;
            let f = get_fric_factor(r_bot, hrad, 1.0e12);
            (f / 185.0).sqrt() * y_full.powf(1.0 / 6.0) // LIBM
        }
    }
}

// ============================================================================
// getRoughFactor — roughness adjustment for artificially-lengthened force main
// ============================================================================

/// Roughness factor stored in xsect.sBot (legacy forcemain_getRoughFactor).
pub fn get_rough_factor(model: FrictionModel, r_bot: f64, length_factor: f64) -> f64 {
    let mut length_factor = length_factor;
    if length_factor <= 0.0 {
        length_factor = 1.0;
    }
    const G: f64 = 32.2;
    match model {
        FrictionModel::HazenWilliams => {
            // legacy: GRAVITY / pow(1.318 * rBot * pow(lengthFactor, 0.54), 1.852)
            let denom = 1.318 * r_bot * length_factor.powf(0.54); // LIBM
            if denom > 0.0 {
                G / denom.powf(1.852) // LIBM
            } else {
                0.0
            }
        }
        // legacy: 1.0 / 8.0 / lengthFactor
        FrictionModel::DarcyWeisbach => 1.0 / (8.0 * length_factor),
    }
}

// ============================================================================
// Batch friction slope
// ============================================================================

/// Friction slope of `count` force mains under one model.
pub fn batch_fric_slope(
    velocity: &[f64],
    hyd_rad: &[f64],
    param: &[f64],
    fric_slope: &mut [f64],
    model: FrictionModel,
    count: i32,
) {
    if model == FrictionModel::HazenWilliams {
        for i in 0..count.max(0) as usize {
            fric_slope[i] = get_fric_slope_hw(velocity[i], hyd_rad[i], param[i]);
        }
    } else {
        for i in 0..count.max(0) as usize {
            fric_slope[i] = get_fric_slope_dw(velocity[i], hyd_rad[i], param[i]);
        }
    }
}
