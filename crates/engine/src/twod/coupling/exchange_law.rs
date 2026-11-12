// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/coupling/ExchangeLaw.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0 (the file's @license tag; it carries
// no SPDX line).
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The regime-aware 1D <-> 2D exchange law, COUPLING_LAW WEIR_ORIFICE (upstream namespace
// openswmm::twoD::exchange): free and submerged weir and orifice, surface -> sewer and
// sewer -> surface, blended over the submergence band with the Hermite smoothstep and combined
// through the p = 16 norm minimum, with the C1-regularized square root phi. SI throughout;
// Q > 0 drains 2D -> 1D, Q < 0 spills 1D -> 2D.
//
// Differences in form, never in result:
//   * OPENSWMM_EXCHANGE_HD (plain `inline` on host builds) is dropped.
//   * softmin's `double& wa, double& wb` outputs are returned with the result as (q, wa, wb).
//   * The struct `Result` is `ExchangeResult` (Rust's prelude already has a `Result`).
//
// Checked bit for bit against upstream's header by gate/diff/twod_exchange.cpp and
// examples/twod_exchange.rs.

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

/// m/s^2
pub const K_GRAVITY: f64 = 9.80665;
/// Submergence blend width.
pub const K_SIGMA_BAND: f64 = 0.1;

/// A coupling point's structure, SI (Geometry).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Geometry {
    /// Crest / rim elevation (m, 2D frame).
    pub z_crest: f64,
    /// Weir crest length (m).
    pub l_weir: f64,
    /// Orifice / open area (m2).
    pub a_orif: f64,
}

/// Discharge coefficients and the phi regularization head (Coefficients).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Coefficients {
    /// Free weir.
    pub c_w: f64,
    /// Submerged weir.
    pub c_sw: f64,
    /// Orifice, surface -> sewer.
    pub c_oi: f64,
    /// Orifice, sewer -> surface.
    pub c_oo: f64,
    /// phi regularization head (m).
    pub h_eps: f64,
}

impl Default for Coefficients {
    fn default() -> Self {
        Self { c_w: 0.540, c_sw: 0.360, c_oi: 0.167, c_oo: 0.160, h_eps: 0.02 }
    }
}

/// C1-regularized sqrt(x): sqrt(x) for x >= eps, (3/(2 sqrt eps)) x - x^2/(2 eps^1.5) below
/// (phi).
#[inline]
pub fn phi(x: f64, eps: f64) -> f64 {
    if !(x > 0.0) {
        return 0.0;
    }
    if x >= eps {
        return x.sqrt();
    }
    let inv_sqrt_e = 1.0 / eps.sqrt();
    (1.5 * inv_sqrt_e) * x - (0.5 * inv_sqrt_e / eps) * x * x
}

/// d phi / dx, finite at 0 (dphi).
#[inline]
pub fn dphi(x: f64, eps: f64) -> f64 {
    if !(x > 0.0) {
        return 1.5 / eps.sqrt();
    }
    if x >= eps {
        return 0.5 / x.sqrt();
    }
    let inv_sqrt_e = 1.0 / eps.sqrt();
    1.5 * inv_sqrt_e - (inv_sqrt_e / eps) * x
}

/// Hermite smoothstep of t clamped to [0, 1] (smoothstep).
#[inline]
pub fn smoothstep(t: f64) -> f64 {
    let t = if t < 0.0 {
        0.0
    } else if t > 1.0 {
        1.0
    } else {
        t
    };
    t * t * (3.0 - 2.0 * t)
}

/// p-norm minimum, p = 16 (softmin), returned as (q, wa, wb): (a^-16 + b^-16)^(-1/16) for
/// a, b > 0, else 0, with the weights dQ/da, dQ/db.
#[inline]
pub fn softmin(a: f64, b: f64) -> (f64, f64, f64) {
    let mut wa = 0.0;
    let mut wb = 0.0;
    if !(a > 0.0) || !(b > 0.0) {
        return (0.0, wa, wb);
    }
    let a_min = a <= b;
    let lo = if a_min { a } else { b };
    let r = lo / if a_min { b } else { a };
    let mut r16 = r * r;
    r16 *= r16;
    r16 *= r16;
    r16 *= r16;
    let root = (1.0 + r16).sqrt().sqrt().sqrt().sqrt();
    let s = 1.0 / root; // Q / min
    let q = lo * s;
    // dQ/dx = (Q/x)^17 for each argument x.
    let mut s17 = s * s;
    s17 *= s17;
    s17 *= s17;
    s17 *= s17;
    s17 *= s;
    let rs = r * s;
    let mut rs17 = rs * rs;
    rs17 *= rs17;
    rs17 *= rs17;
    rs17 *= rs17;
    rs17 *= rs;
    if a_min {
        wa = s17;
        wb = rs17;
    } else {
        wb = s17;
        wa = rs17;
    }
    (q, wa, wb)
}

/// Exchange Q (m3/s) and its head sensitivity G = -dQ/dh1 >= 0 (Result).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ExchangeResult {
    pub q: f64,
    pub g: f64,
}

/// Evaluate the bidirectional regularized weir/orifice exchange and its damping (exchange).
pub fn exchange(geo: &Geometry, c: &Coefficients, h2: f64, h1: f64) -> ExchangeResult {
    let mut r = ExchangeResult::default();
    let s2g = (2.0 * K_GRAVITY).sqrt();
    let z = geo.z_crest;
    let eps = c.h_eps;

    if h2 > h1 {
        let hs = h2 - z;
        if !(hs > 0.0) {
            return r;
        }
        let dh = h2 - h1;
        let s_ = if h1 > z { smoothstep((h1 - z) / (K_SIGMA_BAND * hs)) } else { 0.0 };

        let weir_free = (2.0 / 3.0) * c.c_w * geo.l_weir * s2g * hs * phi(hs, eps);
        let weir_sub = c.c_sw * geo.l_weir * s2g * hs * phi(dh, eps);
        let orif_free = c.c_oi * geo.a_orif * s2g * phi(hs, eps);
        let orif_sub = c.c_oi * geo.a_orif * s2g * phi(dh, eps);
        let qw = (1.0 - s_) * weir_free + s_ * weir_sub;
        let qo = (1.0 - s_) * orif_free + s_ * orif_sub;

        let (q, ww, wo) = softmin(qw, qo);
        r.q = q;
        // Only the submerged forms depend on h1.
        let dqw = s_ * c.c_sw * geo.l_weir * s2g * hs * dphi(dh, eps);
        let dqo = s_ * c.c_oi * geo.a_orif * s2g * dphi(dh, eps);
        r.g = ww * dqw + wo * dqo;
        return r;
    }
    if h1 > h2 && h1 > z {
        let head_dry = h1 - z;
        let s_ = if h2 > z { smoothstep((h2 - z) / (K_SIGMA_BAND * head_dry)) } else { 0.0 };
        let k = c.c_oo * geo.a_orif * s2g;
        r.q = -k * ((1.0 - s_) * phi(head_dry, eps) + s_ * phi(h1 - h2, eps));
        r.g = k * ((1.0 - s_) * dphi(head_dry, eps) + s_ * dphi(h1 - h2, eps));
        return r;
    }
    r
}
