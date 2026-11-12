// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/mesh/QuadVfr.hpp (HydroCouple/openswmm.engine,
// branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache
// License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The volume/free-surface relation of a quadrilateral cell after Begnudelli & Sanders (2007):
// the quad is split along a diagonal chosen from its elevation ordering into two planar
// sub-triangles, and its mean depth is the area-weighted sum of the two triangle closures
// (vfr_closure). Upstream namespace openswmm::twoD.
//
// Differences in form, never in result:
//   * OPENSWMM_KERNEL_FN is dropped.
//   * `const double* zs` (six sorted elevations) is `&[f64]`; quadVfrPrecompute's x/y/z
//     pointers are slices of four and its `double* zs` a `&mut [f64]` of six; its `double&`
//     areas are `&mut f64`.
//
// Checked bit for bit against upstream's header by gate/diff/twod_quadvfr.cpp and
// examples/twod_quadvfr.rs.

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

use super::vfr_closure::{
    vfr_dry_eta, vfr_mean_depth_from_eta, vfr_sort3, vfr_stage_at_wet_fraction, vfr_wet_fraction,
    K_VFR_FLAT_RELIEF,
};

/// Number of doubles of precomputed elevation data per quad cell (two sorted triples).
pub const K_QUAD_VFR_Z: i32 = 6;

/// Wetted-area fraction of the quad at stage eta (quadWetFraction).
#[inline]
pub fn quad_wet_fraction(zs: &[f64], a1: f64, a2: f64, eta: f64) -> f64 {
    let w1 = vfr_wet_fraction(zs[0], zs[1], zs[2], eta);
    let w2 = vfr_wet_fraction(zs[3], zs[4], zs[5], eta);
    (a1 * w1 + a2 * w2) / (a1 + a2)
}

/// Regularized cell-mean depth hbar(eta) of the quad (quadMeanDepthFromEta).
#[inline]
pub fn quad_mean_depth_from_eta(zs: &[f64], a1: f64, a2: f64, eta: f64, eps: f64) -> f64 {
    let d1 = vfr_mean_depth_from_eta(zs[0], zs[1], zs[2], eta, eps);
    let d2 = vfr_mean_depth_from_eta(zs[3], zs[4], zs[5], eta, eps);
    (a1 * d1 + a2 * d2) / (a1 + a2)
}

/// Slope d hbar / d eta of the regularized sum (quadMeanDepthSlope).
pub fn quad_mean_depth_slope(zs: &[f64], a1: f64, a2: f64, eta: f64, eps: f64) -> f64 {
    let mut s = 0.0;
    for k in 0..2usize {
        let z1 = zs[3 * k];
        let z2 = zs[3 * k + 1];
        let z3 = zs[3 * k + 2];
        let a = if k == 0 { a1 } else { a2 };
        let relief = z3 - z1;
        let w;
        if eps <= 0.0 || relief < K_VFR_FLAT_RELIEF {
            w = vfr_wet_fraction(z1, z2, z3, eta);
        } else {
            let eta_s = vfr_stage_at_wet_fraction(z1, z2, z3, eps);
            if eta >= eta_s {
                w = vfr_wet_fraction(z1, z2, z3, eta);
            } else {
                let d = vfr_mean_depth_from_eta(z1, z2, z3, eta, eps);
                w = if d > 0.0 { eps } else { 0.0 };
            }
        }
        s += a * w;
    }
    s / (a1 + a2)
}

/// Regularized free-surface elevation eta(hbar) of the quad, the solver closure
/// (quadEtaFromMeanDepth).
pub fn quad_eta_from_mean_depth(zs: &[f64], a1: f64, a2: f64, mean_depth: f64, eps: f64) -> f64 {
    let a = a1 + a2;
    let zbar1 = (zs[0] + zs[1] + zs[2]) / 3.0;
    let zbar2 = (zs[3] + zs[4] + zs[5]) / 3.0;
    let zw = (a1 * zbar1 + a2 * zbar2) / a;
    let ztop = if zs[2] > zs[5] { zs[2] } else { zs[5] };
    let zlow = if zs[0] < zs[3] { zs[0] } else { zs[3] };
    let relief = ztop - zlow;

    // Flat (or degenerate) cell: the flat closure is exact.
    if relief < K_VFR_FLAT_RELIEF {
        return zw + if mean_depth > 0.0 { mean_depth } else { 0.0 };
    }

    // Fully wet: both sub-triangles submerged - flat closure exact.
    if mean_depth >= ztop - zw {
        return zw + mean_depth;
    }

    // Dry limit: the lower of the two sub-triangle dry stages.
    let dry1 = vfr_dry_eta(zs[0], zs[1], zs[2], eps);
    let dry2 = vfr_dry_eta(zs[3], zs[4], zs[5], eps);
    let lo0 = if dry1 < dry2 { dry1 } else { dry2 };
    if !(mean_depth > 0.0) {
        return lo0;
    }

    // Safeguarded Newton on the monotone sum over [lo0, ztop].
    let mut lo = lo0;
    let mut hi = ztop;
    let mut eta = zw + mean_depth; // flat-closure guess
    if eta <= lo || eta >= hi {
        eta = 0.5 * (lo + hi);
    }
    for _it in 0..80 {
        let f = quad_mean_depth_from_eta(zs, a1, a2, eta, eps) - mean_depth;
        if f > 0.0 {
            hi = eta;
        } else {
            lo = eta;
        }
        let df = quad_mean_depth_slope(zs, a1, a2, eta, eps);
        let mut next = if df > 1.0e-12 { eta - f / df } else { 0.5 * (lo + hi) };
        if next <= lo || next >= hi {
            next = 0.5 * (lo + hi); // safeguard
        }
        if (next - eta).abs() < 1.0e-13 * (1.0 + relief) {
            return next;
        }
        eta = next;
    }
    eta
}

/// d eta / d hbar of the regularized quad closure at stage eta (quadDEtaDMeanDepth).
#[inline]
pub fn quad_d_eta_d_mean_depth(zs: &[f64], a1: f64, a2: f64, eta: f64, eps: f64) -> f64 {
    let mut w = quad_mean_depth_slope(zs, a1, a2, eta, eps);
    if w < eps {
        w = eps;
    }
    if w < 1.0e-12 {
        w = 1.0e-12;
    }
    1.0 / w
}

/// Choose the B&S 2007 diagonal for a quad and emit its precomputed VFR data
/// (quadVfrPrecompute). `x`, `y`, `z` hold the four vertices in cyclic order; `zs` receives the
/// two sorted sub-triangle triples, `a1`/`a2` their planimetric areas. Returns the case (1, 2
/// or 3).
pub fn quad_vfr_precompute(x: &[f64], y: &[f64], z: &[f64], zs: &mut [f64], a1: &mut f64, a2: &mut f64) -> i32 {
    // Sort the four cyclic positions by elevation (stable insertion sort).
    let mut p: [i32; 4] = [0, 1, 2, 3];
    for i in 1..4i32 {
        let key = p[i as usize];
        let mut j = i - 1;
        while j >= 0 && z[p[j as usize] as usize] > z[key as usize] {
            p[(j + 1) as usize] = p[j as usize];
            j -= 1;
        }
        p[(j + 1) as usize] = key;
    }
    let n1 = p[0];
    let n2 = p[1];
    let n3 = p[2];
    let n4 = p[3];
    let adjacent = |a: i32, b: i32| {
        let d = (a - b + 4) % 4;
        d == 1 || d == 3
    };

    let t1: [i32; 3];
    let t2: [i32; 3];
    let kase;
    if !adjacent(n1, n4) {
        // Case 1: n1-n4 is a diagonal
        kase = 1;
        t1 = [n1, n2, n4];
        t2 = [n1, n3, n4];
    } else if adjacent(n2, n1) {
        // Case 2: diagonal n2-n4
        kase = 2;
        t1 = [n1, n2, n4];
        t2 = [n2, n3, n4];
    } else {
        // Case 3: diagonal n3-n4
        kase = 3;
        t1 = [n1, n3, n4];
        t2 = [n2, n3, n4];
    }
    let tri_area = |t: &[i32; 3]| {
        let (t0, t1_, t2_) = (t[0] as usize, t[1] as usize, t[2] as usize);
        let dx1 = x[t1_] - x[t0];
        let dy1 = y[t1_] - y[t0];
        let dx2 = x[t2_] - x[t0];
        let dy2 = y[t2_] - y[t0];
        0.5 * (dx1 * dy2 - dx2 * dy1).abs()
    };
    *a1 = tri_area(&t1);
    *a2 = tri_area(&t2);
    zs[0] = z[t1[0] as usize];
    zs[1] = z[t1[1] as usize];
    zs[2] = z[t1[2] as usize];
    zs[3] = z[t2[0] as usize];
    zs[4] = z[t2[1] as usize];
    zs[5] = z[t2[2] as usize];
    for k in [0usize, 3] {
        // vfrSort3(zs[k], zs[k + 1], zs[k + 2]) on copies written back (no aliasing in Rust)
        let (mut za, mut zb, mut zc) = (zs[k], zs[k + 1], zs[k + 2]);
        vfr_sort3(&mut za, &mut zb, &mut zc);
        zs[k] = za;
        zs[k + 1] = zb;
        zs[k + 2] = zc;
    }
    // Degenerate sub-triangle (collinear): upstream's comment says "fold everything into the
    // other", but the code only zeroes a non-positive (or NaN) area (UPSTREAM_NOTES U350).
    if !(*a1 > 0.0) {
        *a1 = 0.0;
    }
    if !(*a2 > 0.0) {
        *a2 = 0.0;
    }
    kase
}
