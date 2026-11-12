// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/coupling/PondModel.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0 (the file's @license tag; it carries
// no SPDX line).
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The above-rim pond of a 2D-coupled junction (upstream namespace openswmm::twoD::pond): a C1
// surface-area ramp over the band [0, delta] above the rim, its closed-form volume and inverse,
// and the head-equalization cap. DynamicWave.cpp's coupled-junction surface area calls `area`.
//
//   area(x)   = A S(x/delta),               S(t) = 3t^2 - 2t^3
//   volume(x) = A delta (t^3 - t^4/2), t = x/delta   (x <= delta)
//             = A (x - delta/2)                      (x > delta)
//
// Differences in form, never in result: OPENSWMM_POND_HD is dropped; std::cbrt is
// core::cmath::glibc_cbrt; HUGE_VAL is f64::INFINITY.
//
// Checked bit for bit against upstream's header by gate/diff/twod_pond.cpp and
// examples/twod_pond.rs.

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

/// Pond surface area at depth x above the rim (area).
#[inline]
pub fn area(x: f64, a_eq: f64, band: f64) -> f64 {
    if !(x > 0.0) {
        return 0.0;
    }
    if !(band > 0.0) || x >= band {
        return a_eq;
    }
    let t = x / band;
    a_eq * t * t * (3.0 - 2.0 * t)
}

/// Pond volume held between the rim and depth x above it (volume).
#[inline]
pub fn volume(x: f64, a_eq: f64, band: f64) -> f64 {
    if !(x > 0.0) {
        return 0.0;
    }
    if !(band > 0.0) {
        return a_eq * x;
    }
    if x >= band {
        return a_eq * (x - 0.5 * band);
    }
    let t = x / band;
    a_eq * band * t * t * t * (1.0 - 0.5 * t)
}

/// Depth above the rim at which the pond holds v, the inverse of `volume` (depthForVolume).
pub fn depth_for_volume(v: f64, a_eq: f64, band: f64) -> f64 {
    if !(v > 0.0) || !(a_eq > 0.0) {
        return 0.0;
    }
    if !(band > 0.0) {
        return v / a_eq;
    }
    let full = 0.5 * a_eq * band; // volume at x = band
    if v >= full {
        return v / a_eq + 0.5 * band;
    }
    // t^3 - t^4/2 = r on [0, 1], monotone: safeguarded Newton from t ~ cbrt(r).
    let r = v / (a_eq * band);
    let mut lo = 0.0;
    let mut hi = 1.0;
    let mut t = glibc_cbrt(r); // LIBM (glibc cbrt)
    if t > 1.0 {
        t = 1.0;
    }
    for _it in 0..60 {
        let f = t * t * t * (1.0 - 0.5 * t) - r;
        if f > 0.0 {
            hi = t;
        } else {
            lo = t;
        }
        let df = t * t * (3.0 - 2.0 * t);
        let mut tn = if df > 0.0 { t - f / df } else { 0.5 * (lo + hi) };
        if !(tn > lo && tn < hi) {
            tn = 0.5 * (lo + hi);
        }
        if (tn - t).abs() <= 1e-15 {
            t = tn;
            break;
        }
        t = tn;
    }
    t * band
}

/// Head-equalization cap on one exchange, |Q| <= |dh| A1 / T (equalizeCap); +inf when the
/// window or area is unknown.
#[inline]
pub fn equalize_cap(dh: f64, window: f64, a1: f64) -> f64 {
    if !(window > 0.0) || !(a1 > 0.0) {
        return f64::INFINITY; // HUGE_VAL
    }
    dh.abs() * a1 / window
}
