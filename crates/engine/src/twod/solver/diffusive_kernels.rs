// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/solver/DiffusiveKernels.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// Face law and step bound for MOMENTUM_EQUATION DIFFUSIVE_WAVE (upstream namespace
// openswmm::twoD::diffusive): the Manning quasi-steady face discharge of Hunter et al. (2005)
// with a regularized slope, and the explicit-diffusion step bound.
//
// Differences in form, never in result: OPENSWMM_KERNEL_FN is dropped; std::cbrt is
// core::cmath::glibc_cbrt.
//
// Checked bit for bit against upstream's header by gate/diff/twod_diffusive.cpp and
// examples/twod_diffusive.rs.

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

/// Manning quasi-steady face discharge, m2/s positive L->R (faceDischarge).
#[inline]
pub fn face_discharge(hf: f64, slope: f64, nf: f64, s_eps: f64) -> f64 {
    if !(hf > 0.0) || !(nf > 0.0) || slope == 0.0 {
        return 0.0;
    }
    let h53 = hf * glibc_cbrt(hf * hf); // LIBM (glibc cbrt)
    let sa = slope.abs();
    let sr = (if sa > s_eps { sa } else { s_eps }).sqrt();
    -h53 * slope / (nf * sr)
}

/// Explicit-diffusion step bound for a cell (cellDiffusiveDt); 1e30 s for a dry cell.
#[inline]
pub fn cell_diffusive_dt(alpha: f64, lchar: f64, h: f64, n: f64, slope_abs: f64, s_eps: f64) -> f64 {
    if !(h > 0.0) {
        return 1.0e30;
    }
    let h53 = h * glibc_cbrt(h * h); // LIBM (glibc cbrt)
    let s = if slope_abs > s_eps { slope_abs } else { s_eps };
    alpha * lchar * lchar * n * s.sqrt() / (4.0 * h53)
}
