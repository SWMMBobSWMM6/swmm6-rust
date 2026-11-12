// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/solver/SweKernels.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// Face kernels for MOMENTUM_EQUATION FULL_SWE (upstream namespace openswmm::twoD::swe): the
// Audusse hydrostatic reconstruction, the HLLC Riemann flux (Toro 2001) in the face-normal
// frame with its wave-speed estimates, the bed-slope corrections, the Barth-Jespersen limiter,
// the semi-implicit Manning friction on a cell's discharge, the critical depth and the
// subcritical prescribed-discharge ghost depth.
//
// Differences in form, never in result:
//   * OPENSWMM_KERNEL_FN is dropped.
//   * physicalFlux, waveSpeeds and hllcFlux return their `double&` outputs as tuples
//     ((fh, fn, ft), (sl, sr), (fh, fn, ft, sstar)); `fn` is a Rust keyword, so it is `fn_`.
//   * faceFlux, faceFluxReconBed and faceFluxRecon keep their `FaceFlux&` and `double&`
//     outputs as `&mut` and return the bool.
//   * frictionUpdate's in/out `double& qx, double& qy` are `&mut f64`.
//   * std::cbrt is core::cmath::glibc_cbrt.
//
// Checked bit for bit against upstream's header by gate/diff/twod_swe.cpp and
// examples/twod_swe.rs.

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

pub const K_GRAVITY: f64 = 9.80665;

/// One side of a face in the face-normal frame after hydrostatic reconstruction.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct FaceSide {
    pub h: f64,
    pub un: f64,
    pub ut: f64,
}

/// Flux through the face, positive L->R: mass (m2/s), momentum in x/y (m3/s2), and the HLLC
/// contact speed.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct FaceFlux {
    pub mass: f64,
    pub mx: f64,
    pub my: f64,
    pub sstar: f64,
}

/// Physical flux of the SWE in the normal frame (physicalFlux), returned as (fh, fn, ft).
#[inline]
pub fn physical_flux(s: &FaceSide) -> (f64, f64, f64) {
    let fh = s.h * s.un;
    let fn_ = s.h * s.un * s.un + 0.5 * K_GRAVITY * s.h * s.h;
    let ft = s.h * s.un * s.ut;
    (fh, fn_, ft)
}

/// Wave-speed estimates S_L <= S_R (waveSpeeds), returned as (sl, sr).
#[inline]
pub fn wave_speeds(l: &FaceSide, r: &FaceSide) -> (f64, f64) {
    let cl = (K_GRAVITY * l.h).sqrt();
    let cr = (K_GRAVITY * r.h).sqrt();
    let sl;
    let sr;
    if l.h > 0.0 && r.h > 0.0 {
        let ustar = 0.5 * (l.un + r.un) + cl - cr;
        let cstar = 0.5 * (cl + cr) + 0.25 * (l.un - r.un);
        sl = cmin(l.un - cl, ustar - cstar);
        sr = cmax(r.un + cr, ustar + cstar);
    } else if l.h > 0.0 {
        // dry right
        sl = l.un - cl;
        sr = l.un + 2.0 * cl;
    } else if r.h > 0.0 {
        // dry left
        sl = r.un - 2.0 * cr;
        sr = r.un + cr;
    } else {
        sl = 0.0;
        sr = 0.0;
    }
    (sl, sr)
}

/// HLLC flux in the face-normal frame (hllcFlux), returned as (fh, fn, ft, sstar).
#[inline]
pub fn hllc_flux(l: &FaceSide, r: &FaceSide) -> (f64, f64, f64, f64) {
    if l.h <= 0.0 && r.h <= 0.0 {
        return (0.0, 0.0, 0.0, 0.0);
    }
    let (sl, sr) = wave_speeds(l, r);
    // Contact speed (Toro Eq. 10.70); the mean normal velocity when the denominator degenerates.
    let dl = l.h * (l.un - sl);
    let dr = r.h * (r.un - sr);
    let den = dr - dl;
    let sstar = if den.abs() > 1.0e-14 { (sl * dr - sr * dl) / den } else { 0.5 * (l.un + r.un) };

    let fh;
    let fn_;
    let ft;
    if sl >= 0.0 {
        let (fhl, fnl, ftl) = physical_flux(l);
        fh = fhl;
        fn_ = fnl;
        ft = ftl;
    } else if sr <= 0.0 {
        let (fhr, fnr, ftr) = physical_flux(r);
        fh = fhr;
        fn_ = fnr;
        ft = ftr;
    } else {
        // HLL middle state for mass and normal momentum.
        let (fhl, fnl, _ftl) = physical_flux(l);
        let (fhr, fnr, _ftr) = physical_flux(r);
        let inv = 1.0 / (sr - sl);
        fh = (sr * fhl - sl * fhr + sl * sr * (r.h - l.h)) * inv;
        fn_ = (sr * fnl - sl * fnr + sl * sr * (r.h * r.un - l.h * l.un)) * inv;
        // Tangential momentum: passive scalar upwinded on the contact.
        ft = fh * if sstar >= 0.0 { l.ut } else { r.ut };
    }
    (fh, fn_, ft, sstar)
}

/// Full face evaluation (faceFlux): hydrostatic reconstruction, rotation, HLLC, rotation back
/// and the per-side bed-slope corrections. Returns true when the face carries flux or a
/// bed-pressure correction.
#[inline]
pub fn face_flux(
    eta_l: f64,
    h_l: f64,
    qx_l: f64,
    qy_l: f64,
    eta_r: f64,
    h_r: f64,
    qx_r: f64,
    qy_r: f64,
    nx: f64,
    ny: f64,
    h_dry: f64,
    out: &mut FaceFlux,
    corr_l_x: &mut f64,
    corr_l_y: &mut f64,
    corr_r_x: &mut f64,
    corr_r_y: &mut f64,
) -> bool {
    // Beds consistent with the storage: z = eta - hbar.
    let z_l = eta_l - h_l;
    let z_r = eta_r - h_r;
    let zf = if z_l > z_r { z_l } else { z_r };
    let mut l = FaceSide::default();
    let mut r = FaceSide::default();
    l.h = if eta_l - zf > 0.0 { eta_l - zf } else { 0.0 };
    r.h = if eta_r - zf > 0.0 { eta_r - zf } else { 0.0 };
    let tx = -ny;
    let ty = nx;
    if h_l > h_dry && l.h > 0.0 {
        let ux = qx_l / h_l;
        let uy = qy_l / h_l;
        l.un = ux * nx + uy * ny;
        l.ut = ux * tx + uy * ty;
    }
    if h_r > h_dry && r.h > 0.0 {
        let ux = qx_r / h_r;
        let uy = qy_r / h_r;
        r.un = ux * nx + uy * ny;
        r.ut = ux * tx + uy * ty;
    }
    // Hydrostatic pressure remains on the wet side of a blocked face.
    let c_l = 0.5 * K_GRAVITY * (l.h * l.h - h_l * h_l);
    let c_r = 0.5 * K_GRAVITY * (r.h * r.h - h_r * h_r);
    *corr_l_x = c_l * nx;
    *corr_l_y = c_l * ny;
    *corr_r_x = -c_r * nx;
    *corr_r_y = -c_r * ny;
    if l.h <= 0.0 && r.h <= 0.0 {
        *out = FaceFlux::default();
        return c_l != 0.0 || c_r != 0.0;
    }
    let (fh, fn_, ft, sstar) = hllc_flux(&l, &r);
    out.mass = fh;
    out.mx = fn_ * nx + ft * tx;
    out.my = fn_ * ny + ft * ty;
    out.sstar = sstar;
    true
}

/// Second-order face flux with independently reconstructed depth/bed (faceFluxReconBed).
#[inline]
pub fn face_flux_recon_bed(
    eta_lf: f64,
    ux_lf: f64,
    uy_lf: f64,
    z_l: f64,
    h_l_cell: f64,
    eta_rf: f64,
    ux_rf: f64,
    uy_rf: f64,
    z_r: f64,
    h_r_cell: f64,
    nx: f64,
    ny: f64,
    h_dry: f64,
    out: &mut FaceFlux,
    corr_l_x: &mut f64,
    corr_l_y: &mut f64,
    corr_r_x: &mut f64,
    corr_r_y: &mut f64,
    z_lf: f64,
    z_rf: f64,
) -> bool {
    let zf = if z_lf > z_rf { z_lf } else { z_rf };
    let mut l = FaceSide::default();
    let mut r = FaceSide::default();
    l.h = if eta_lf - zf > 0.0 { eta_lf - zf } else { 0.0 };
    r.h = if eta_rf - zf > 0.0 { eta_rf - zf } else { 0.0 };
    let tx = -ny;
    let ty = nx;
    if h_l_cell > h_dry && l.h > 0.0 {
        l.un = ux_lf * nx + uy_lf * ny;
        l.ut = ux_lf * tx + uy_lf * ty;
    }
    if h_r_cell > h_dry && r.h > 0.0 {
        r.un = ux_rf * nx + uy_rf * ny;
        r.ut = ux_rf * tx + uy_rf * ty;
    }
    let h_lf = cmax(0.0, eta_lf - z_lf);
    let h_rf = cmax(0.0, eta_rf - z_rf);
    // Hydrostatic face correction plus a centered within-cell bed source.
    let c_l = 0.5 * K_GRAVITY * (l.h * l.h - h_lf * h_lf - (h_l_cell + h_lf) * (z_lf - z_l));
    let c_r = 0.5 * K_GRAVITY * (r.h * r.h - h_rf * h_rf - (h_r_cell + h_rf) * (z_rf - z_r));
    *corr_l_x = c_l * nx;
    *corr_l_y = c_l * ny;
    *corr_r_x = -c_r * nx;
    *corr_r_y = -c_r * ny;
    if l.h <= 0.0 && r.h <= 0.0 {
        *out = FaceFlux::default();
        return c_l != 0.0 || c_r != 0.0;
    }
    let (fh, fn_, ft, sstar) = hllc_flux(&l, &r);
    out.mass = fh;
    out.mx = fn_ * nx + ft * tx;
    out.my = fn_ * ny + ft * ty;
    out.sstar = sstar;
    true
}

/// Piecewise-constant-bed compatibility entry point (faceFluxRecon): cell beds as face beds.
#[inline]
pub fn face_flux_recon(
    eta_lf: f64,
    ux_lf: f64,
    uy_lf: f64,
    z_l: f64,
    h_l_cell: f64,
    eta_rf: f64,
    ux_rf: f64,
    uy_rf: f64,
    z_r: f64,
    h_r_cell: f64,
    nx: f64,
    ny: f64,
    h_dry: f64,
    out: &mut FaceFlux,
    corr_l_x: &mut f64,
    corr_l_y: &mut f64,
    corr_r_x: &mut f64,
    corr_r_y: &mut f64,
) -> bool {
    face_flux_recon_bed(
        eta_lf, ux_lf, uy_lf, z_l, h_l_cell, eta_rf, ux_rf, uy_rf, z_r, h_r_cell, nx, ny, h_dry, out, corr_l_x,
        corr_l_y, corr_r_x, corr_r_y, z_l, z_r,
    )
}

/// Barth-Jespersen limiter factor (bjLimiter).
#[inline]
pub fn bj_limiter(w: f64, wf: f64, wmin: f64, wmax: f64) -> f64 {
    let d = wf - w;
    if d > 1.0e-14 {
        cmin(1.0, (wmax - w) / d)
    } else if d < -1.0e-14 {
        cmin(1.0, (wmin - w) / d)
    } else {
        1.0
    }
}

/// Semi-implicit Manning friction on a cell's unit discharge over dt (frictionUpdate).
#[inline]
pub fn friction_update(qx: &mut f64, qy: &mut f64, h: f64, n: f64, dt: f64) {
    if !(h > 0.0) {
        *qx = 0.0;
        *qy = 0.0;
        return;
    }
    let qm = (*qx * *qx + *qy * *qy).sqrt();
    if qm == 0.0 {
        return;
    }
    let h73 = h * h * glibc_cbrt(h); // LIBM (glibc cbrt)
    let den = 1.0 + K_GRAVITY * dt * n * n * qm / h73;
    *qx /= den;
    *qy /= den;
}

/// Critical depth (q^2/g)^(1/3) for a per-metre discharge q (criticalDepth).
#[inline]
pub fn critical_depth(q: f64) -> f64 {
    glibc_cbrt(q * q / K_GRAVITY) // LIBM (glibc cbrt)
}

/// Ghost depth of a SUBCRITICAL prescribed-discharge boundary
/// (depthFromInvariantAndDischarge): Newton on s = sqrt(h) for q/h + 2 sqrt(g h) = r.
pub fn depth_from_invariant_and_discharge(r: f64, q: f64, h_interior: f64) -> f64 {
    let rg = K_GRAVITY.sqrt();
    let mut s = (if h_interior > 0.0 { h_interior } else { 1.0e-6 }).sqrt();
    for _it in 0..40 {
        let f = 2.0 * rg * s * s * s - r * s * s + q;
        let df = 6.0 * rg * s * s - 2.0 * r * s;
        if !(df.abs() > 1.0e-12) {
            break;
        }
        let s1 = s - f / df;
        if !(s1 > 0.0) || !s1.is_finite() {
            break;
        }
        let step = (s1 - s).abs();
        s = s1;
        if step < 1.0e-12 * (1.0 + s) {
            let h = s * s;
            return if h > 0.0 && h.is_finite() { h } else { h_interior };
        }
    }
    let h = s * s;
    // Accept only a converged, physical root; otherwise keep the interior depth.
    if h > 0.0 && h.is_finite() && (2.0 * rg * s * s * s - r * s * s + q).abs() < 1.0e-6 * (1.0 + q.abs()) {
        h
    } else {
        h_interior
    }
}
