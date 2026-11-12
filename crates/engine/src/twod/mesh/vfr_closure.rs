// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/mesh/VfrClosure.hpp (HydroCouple/openswmm.engine,
// branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache
// License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The volume/free-surface relation (VFR) of a planar-bed triangular cell after Begnudelli &
// Sanders (2006, 2007), with its wet-area-fraction regularization (upstream namespace
// openswmm::twoD). With sorted vertex beds z1 <= z2 <= z3, zbar their mean and hbar = V/A:
//
//   z1 < eta <= z2 : hbar = (eta-z1)^3 / (3(z2-z1)(z3-z1))
//   z2 < eta <= z3 : hbar = (eta-zbar) + (z3-eta)^3 / (3(z3-z1)(z3-z2))
//   eta >= z3      : hbar = eta - zbar
//
// Differences in form, never in result:
//   * OPENSWMM_KERNEL_FN (plain `inline` on host builds) is dropped.
//   * vfrSort3's three `double&` are `&mut f64`.
//   * std::cbrt is core::cmath::glibc_cbrt (glibc's s_cbrt.c; Rust's f64::cbrt is not glibc's).
//
// Checked bit for bit against upstream's header by gate/diff/twod_vfr.cpp and
// examples/twod_vfr.rs.

// The negated comparisons (`!(x > 0.0)`) are upstream's and decide what a NaN does; the
// swaps, argument lists, clamps and late initializations follow upstream's code as written.
#![allow(
    clippy::neg_cmp_op_on_partial_ord,
    clippy::too_many_arguments,
    clippy::manual_swap,
    clippy::needless_late_init,
    clippy::manual_range_contains,
    clippy::manual_clamp
)]

use crate::core::cmath::glibc_cbrt;

/// Relief below which a cell is treated as flat (metres).
pub const K_VFR_FLAT_RELIEF: f64 = 1.0e-9;

/// Sort three vertex elevations in place so z1 <= z2 <= z3 (vfrSort3).
#[inline]
pub fn vfr_sort3(z1: &mut f64, z2: &mut f64, z3: &mut f64) {
    if *z1 > *z2 {
        let t = *z1;
        *z1 = *z2;
        *z2 = t;
    }
    if *z2 > *z3 {
        let t = *z2;
        *z2 = *z3;
        *z3 = t;
    }
    if *z1 > *z2 {
        let t = *z1;
        *z1 = *z2;
        *z2 = t;
    }
}

/// Wetted-area fraction A_wet/A of the planar-bed cell at stage `eta` (vfrWetFraction).
/// Inputs must be sorted.
#[inline]
pub fn vfr_wet_fraction(z1: f64, z2: f64, z3: f64, eta: f64) -> f64 {
    let relief = z3 - z1;
    if relief < K_VFR_FLAT_RELIEF {
        return if eta > z1 { 1.0 } else { 0.0 };
    }
    if eta <= z1 {
        return 0.0;
    }
    if eta >= z3 {
        return 1.0;
    }
    if eta <= z2 {
        let d21 = z2 - z1;
        if d21 < K_VFR_FLAT_RELIEF {
            return 0.0; // measure-zero band z1==z2
        }
        let t = eta - z1;
        return t * t / (d21 * relief);
    }
    // z2 < eta < z3 (z3 > z2 strictly here)
    let d32 = z3 - z2;
    let t = z3 - eta;
    1.0 - t * t / (relief * d32)
}

/// EXACT cell-mean depth hbar(eta) of the planar-bed cell (vfrMeanDepthFromEtaExact).
/// Inputs must be sorted.
#[inline]
pub fn vfr_mean_depth_from_eta_exact(z1: f64, z2: f64, z3: f64, eta: f64) -> f64 {
    let zbar = (z1 + z2 + z3) / 3.0;
    let relief = z3 - z1;
    if relief < K_VFR_FLAT_RELIEF {
        return if eta > zbar { eta - zbar } else { 0.0 };
    }
    if eta <= z1 {
        return 0.0;
    }
    if eta >= z3 {
        return eta - zbar;
    }
    if eta <= z2 {
        let d21 = z2 - z1;
        if d21 < K_VFR_FLAT_RELIEF {
            return 0.0; // measure-zero band z1==z2
        }
        let t = eta - z1;
        return t * t * t / (3.0 * d21 * relief);
    }
    let d32 = z3 - z2;
    let t = z3 - eta;
    (eta - zbar) + t * t * t / (3.0 * relief * d32)
}

/// Stage at which the wetted-area fraction equals `eps` (vfrStageAtWetFraction). Inputs must
/// be sorted and non-flat.
#[inline]
pub fn vfr_stage_at_wet_fraction(z1: f64, z2: f64, z3: f64, eps: f64) -> f64 {
    let relief = z3 - z1;
    let d21 = z2 - z1;
    // Lower branch reaches wet fraction (z2-z1)/(z3-z1) at eta = z2.
    if eps <= d21 / relief {
        return z1 + (eps * d21 * relief).sqrt();
    }
    let d32 = z3 - z2;
    z3 - ((1.0 - eps) * relief * d32).sqrt()
}

/// Regularized cell-mean depth hbar(eta) (vfrMeanDepthFromEta): exact for wet fraction >= eps,
/// tangent-line continuation below, clamped at 0. Inputs must be sorted.
#[inline]
pub fn vfr_mean_depth_from_eta(z1: f64, z2: f64, z3: f64, eta: f64, eps: f64) -> f64 {
    let relief = z3 - z1;
    if eps <= 0.0 || relief < K_VFR_FLAT_RELIEF {
        return vfr_mean_depth_from_eta_exact(z1, z2, z3, eta);
    }
    let eta_s = vfr_stage_at_wet_fraction(z1, z2, z3, eps);
    if eta >= eta_s {
        return vfr_mean_depth_from_eta_exact(z1, z2, z3, eta);
    }
    let h_s = vfr_mean_depth_from_eta_exact(z1, z2, z3, eta_s);
    let h = h_s - eps * (eta_s - eta);
    if h > 0.0 {
        h
    } else {
        0.0
    }
}

/// Regularized free-surface elevation eta(hbar), the solver closure (vfrEtaFromMeanDepth).
/// Inputs must be sorted.
pub fn vfr_eta_from_mean_depth(z1: f64, z2: f64, z3: f64, mean_depth: f64, eps: f64) -> f64 {
    let zbar = (z1 + z2 + z3) / 3.0;
    let relief = z3 - z1;

    // Flat (or degenerate) cell: the flat closure is exact.
    if relief < K_VFR_FLAT_RELIEF {
        return zbar + if mean_depth > 0.0 { mean_depth } else { 0.0 };
    }

    // Fully wet: flat closure is exact (checked first; see upstream's comment).
    if mean_depth >= z3 - zbar {
        return zbar + mean_depth;
    }

    // Regularized tail: below the switch depth h_s the closure is linear. (Upstream declares
    // eta_s = z1, h_s = 0.0 before the test; nothing reads them outside this branch.)
    if eps > 0.0 {
        let eta_s = vfr_stage_at_wet_fraction(z1, z2, z3, eps);
        let h_s = vfr_mean_depth_from_eta_exact(z1, z2, z3, eta_s);
        let h = if mean_depth > 0.0 { mean_depth } else { 0.0 };
        if h <= h_s {
            return eta_s - (h_s - h) / eps;
        }
    } else if !(mean_depth > 0.0) {
        return z1; // exact relation, dry limit
    }

    // Lower branch (z1 < eta <= z2): closed-form cube root.
    let d21 = z2 - z1;
    let h_at_z2 = d21 * d21 / (3.0 * relief);
    if mean_depth <= h_at_z2 {
        return z1 + glibc_cbrt(3.0 * mean_depth * d21 * relief); // LIBM (glibc cbrt)
    }

    // z2 == z3 (to rounding): the upper branch is an empty interval.
    if z3 - z2 < K_VFR_FLAT_RELIEF {
        return zbar + mean_depth;
    }

    // Upper branch (z2 < eta < z3): safeguarded Newton on the bracket [z2, z3].
    let d32 = z3 - z2;
    let denom = 3.0 * relief * d32;
    let mut lo = z2;
    let mut hi = z3;
    let mut eta = zbar + mean_depth; // flat-closure initial guess
    if eta <= lo || eta >= hi {
        eta = 0.5 * (lo + hi);
    }
    for _it in 0..64 {
        let dz3 = z3 - eta;
        let f = (eta - zbar) + dz3 * dz3 * dz3 / denom - mean_depth;
        if f > 0.0 {
            hi = eta;
        } else {
            lo = eta;
        }
        let df = 1.0 - dz3 * dz3 / (relief * d32); // A_wet/A
        let mut next = if df > 1.0e-12 { eta - f / df } else { 0.5 * (lo + hi) };
        if next <= lo || next >= hi {
            next = 0.5 * (lo + hi); // safeguard
        }
        if (next - eta).abs() < 1.0e-12 * (1.0 + relief) {
            return next;
        }
        eta = next;
    }
    eta
}

/// Dry-cell free surface of the regularized closure, eta(V = 0) (vfrDryEta). Inputs must be
/// sorted.
#[inline]
pub fn vfr_dry_eta(z1: f64, z2: f64, z3: f64, eps: f64) -> f64 {
    vfr_eta_from_mean_depth(z1, z2, z3, 0.0, eps)
}

/// d eta / d hbar of the regularized closure at stage `eta`: 1/max(w(eta), eps)
/// (vfrDEtaDMeanDepth). Inputs must be sorted.
#[inline]
pub fn vfr_d_eta_d_mean_depth(z1: f64, z2: f64, z3: f64, eta: f64, eps: f64) -> f64 {
    let mut w = vfr_wet_fraction(z1, z2, z3, eta);
    if w < eps {
        w = eps;
    }
    if w < 1.0e-12 {
        w = 1.0e-12; // eps == 0 safety (exact relation, dry cell)
    }
    1.0 / w
}
