// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/solver/InertialKernels.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The scalar kernels of the explicit local-inertial marcher (upstream namespace
// openswmm::twoD::inertial): the de Almeida & Bates (2013) face update with semi-implicit
// Manning friction, the convective term, the Froude cap, the CFL bound, the positivity scale,
// the face flow depths (MEAN and the Begnudelli & Sanders Eq. 14 wetted-edge form) and the
// scalar cores of the FLAT/VFR volume <-> free-surface closure.
//
// NOT TRANSLATED here (they gather from MeshData and SolverOptions2D, the 2D data structs of
// wave 1, and belong with the mesh in wave 3): cellEtaDepthQuad, cellEtaDepth, cellWetArea,
// cellVolumeFromEta. Each is a gather that forwards to the scalar cores below.
//
// Differences in form, never in result:
//   * OPENSWMM_KERNEL_FN is dropped.
//   * The scalar closures' `double& eta, double& depth` outputs are a returned (eta, depth).
//   * `const double* zs` is `&[f64]`.
//   * inertialFaceUpdate's defaulted `adv = 0.0` is an explicit argument.
//   * std::cbrt is core::cmath::glibc_cbrt; std::log and std::expm1 are f64::ln and
//     f64::exp_m1, which call glibc's log and expm1 (both checked by the differential).
//   * froudeCap's std::clamp is written out as libstdc++ 13 defines it, min(max(v, lo), hi),
//     not core::cmath::cclamp (they differ when hi < lo; see froude_cap).
//
// Checked bit for bit against upstream's header by gate/diff/twod_inertial.cpp and
// examples/twod_inertial.rs.

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

use crate::core::cmath::{cmax, cmin, glibc_cbrt};
use crate::twod::mesh::quad_vfr::{quad_eta_from_mean_depth, quad_mean_depth_from_eta};
use crate::twod::mesh::vfr_closure::{vfr_eta_from_mean_depth, vfr_mean_depth_from_eta, vfr_sort3};

/// Gravity, m/s^2 (matches the existing kernels).
pub const K_GRAVITY: f64 = 9.80665;

/// Free-surface difference below which a face slope is numerical zero (m).
pub const K_ETA_DEADBAND: f64 = 1.0e-12;

/// Fraction of reconstructed momentum mixed over one physical face step
/// (momentumBlendFraction).
#[inline]
pub fn momentum_blend_fraction(theta: f64, dt: f64, reference_step: f64) -> f64 {
    if !(reference_step > 0.0) {
        return 1.0 - theta;
    }
    if !(dt > 0.0) || theta == 1.0 {
        return 0.0;
    }
    -(theta.ln() * (dt / reference_step)).exp_m1() // LIBM (log, expm1)
}

/// |(x, y)| for the marcher's discharge vectors (qMagnitude): a plain sqrt, not hypot.
#[inline]
pub fn q_magnitude(x: f64, y: f64) -> f64 {
    (x * x + y * y).sqrt()
}

/// Scalar core of the V -> (eta, depth) closure (etaDepthScalar), returned as (eta, depth).
/// The vertex beds are sorted on local copies.
#[inline]
pub fn eta_depth_scalar(
    area: f64,
    cz: f64,
    z1: f64,
    z2: f64,
    z3: f64,
    vfr: bool,
    vfr_min_wet_frac: f64,
    v_in: f64,
) -> (f64, f64) {
    let v = if v_in > 0.0 { v_in } else { 0.0 };
    let depth = if area > 1.0e-30 { v / area } else { 0.0 };
    let eta;
    if vfr {
        let (mut z1, mut z2, mut z3) = (z1, z2, z3);
        vfr_sort3(&mut z1, &mut z2, &mut z3);
        eta = vfr_eta_from_mean_depth(z1, z2, z3, depth, vfr_min_wet_frac);
    } else {
        eta = cz + depth;
    }
    (eta, depth)
}

/// Scalar core of the eta -> V inverse (volumeFromEtaScalar).
#[inline]
pub fn volume_from_eta_scalar(
    area: f64,
    cz: f64,
    z1: f64,
    z2: f64,
    z3: f64,
    vfr: bool,
    vfr_min_wet_frac: f64,
    eta: f64,
) -> f64 {
    if vfr {
        let (mut z1, mut z2, mut z3) = (z1, z2, z3);
        vfr_sort3(&mut z1, &mut z2, &mut z3);
        return area * vfr_mean_depth_from_eta(z1, z2, z3, eta, vfr_min_wet_frac);
    }
    let d = eta - cz;
    if d > 0.0 {
        area * d
    } else {
        0.0
    }
}

/// Quad (B&S 2007 two-plane) V -> (eta, depth) closure core (etaDepthQuadScalar), returned as
/// (eta, depth).
#[inline]
pub fn eta_depth_quad_scalar(
    area: f64,
    cz: f64,
    zs: &[f64],
    a1: f64,
    a2: f64,
    vfr: bool,
    vfr_min_wet_frac: f64,
    v_in: f64,
) -> (f64, f64) {
    let v = if v_in > 0.0 { v_in } else { 0.0 };
    let depth = if area > 1.0e-30 { v / area } else { 0.0 };
    let eta = if vfr { quad_eta_from_mean_depth(zs, a1, a2, depth, vfr_min_wet_frac) } else { cz + depth };
    (eta, depth)
}

/// Quad eta -> V inverse (volumeFromEtaQuadScalar).
#[inline]
pub fn volume_from_eta_quad_scalar(
    area: f64,
    cz: f64,
    zs: &[f64],
    a1: f64,
    a2: f64,
    vfr: bool,
    vfr_min_wet_frac: f64,
    eta: f64,
) -> f64 {
    if vfr {
        return area * quad_mean_depth_from_eta(zs, a1, a2, eta, vfr_min_wet_frac);
    }
    let d = eta - cz;
    if d > 0.0 {
        area * d
    } else {
        0.0
    }
}

/// Flow depth at a face, FACE_RECONSTRUCTION = MEAN (faceFlowDepth).
#[inline]
pub fn face_flow_depth(eta_l: f64, eta_r: f64, zface: f64) -> f64 {
    cmax(eta_l, eta_r) - zface
}

/// Wetted-edge flow depth over an edge with endpoint beds z_lo <= z_hi, Begnudelli & Sanders
/// (2007) Eq. 14 (faceDepthFromEta).
#[inline]
pub fn face_depth_from_eta(eta: f64, z_lo: f64, z_hi: f64) -> f64 {
    if eta <= z_lo {
        return 0.0;
    }
    let dz = z_hi - z_lo;
    if dz < 1.0e-9 {
        return eta - z_lo; // level edge
    }
    if eta <= z_hi {
        let t = eta - z_lo;
        return t * t / (2.0 * dz);
    }
    eta - 0.5 * (z_lo + z_hi)
}

/// Flow depth at an interior face under FACE_RECONSTRUCTION = VFR_FACE (faceFlowDepthVfr).
#[inline]
pub fn face_flow_depth_vfr(eta_l: f64, eta_r: f64, z_lo: f64, z_hi: f64) -> f64 {
    face_depth_from_eta(cmax(eta_l, eta_r), z_lo, z_hi)
}

/// One local-inertial face update over dt (inertialFaceUpdate). `_q` is unused upstream too.
#[inline]
pub fn inertial_face_update(
    _q: f64,
    qhat: f64,
    hf: f64,
    dt: f64,
    slope: f64,
    n2: f64,
    q_mag: f64,
    adv: f64,
) -> f64 {
    let h73 = hf * hf * glibc_cbrt(hf); // LIBM (glibc cbrt)
    let num = qhat - dt * (K_GRAVITY * hf * slope + adv);
    let den = 1.0 + K_GRAVITY * dt * n2 * q_mag / h73;
    num / den
}

/// Convective momentum flux difference at an interior face, ADVECTION (inertialAdvection).
#[inline]
pub fn inertial_advection(q_f: f64, un_l: f64, h_l: f64, un_r: f64, h_r: f64, inv_dx_normal: f64) -> f64 {
    let fl = un_l * if un_l > 0.0 { h_l * un_l } else { q_f };
    let fr = un_r * if un_r > 0.0 { q_f } else { h_r * un_r };
    (fr - fl) * inv_dx_normal
}

/// Froude-number clamp |q| <= Fr_max h_f sqrt(g h_f) (froudeCap).
#[inline]
pub fn froude_cap(q: f64, hf: f64, fr_max: f64) -> f64 {
    let qcap = fr_max * hf * (K_GRAVITY * hf).sqrt();
    // std::clamp(q, -qcap, qcap) as the reference build's libstdc++ 13 defines it,
    // std::min(std::max(v, lo), hi). core::cmath::cclamp is the ternary form, which returns lo
    // instead of hi when hi < lo and v < lo (here: a negative fr_max); PLAN.md, hazards.
    cmin(cmax(q, -qcap), qcap)
}

/// Per-cell CFL stable step (cellCflDt).
#[inline]
pub fn cell_cfl_dt(alpha: f64, lchar: f64, h: f64, speed: f64) -> f64 {
    let c = (K_GRAVITY * h).sqrt() + speed;
    if c > 1.0e-12 {
        alpha * lchar / c
    } else {
        1.0e30
    }
}

/// Positivity scale factor for a cell about to export volume (positivityScale).
#[inline]
pub fn positivity_scale(v: f64, outflow_m3s: f64, dt: f64, beta: f64) -> f64 {
    if outflow_m3s <= 0.0 || dt <= 0.0 {
        return 1.0;
    }
    let avail = beta * cmax(v, 0.0);
    let take = outflow_m3s * dt;
    if take <= avail {
        1.0
    } else {
        avail / take
    }
}
