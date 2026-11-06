// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/fv/NativeGmres.cpp and NativeGmres.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 6 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Restarted GMRES with modified Gram-Schmidt (one reorthogonalization pass), Givens
// rotations and an optional right preconditioner; the conservative quality mixer
// (quality::ConservativeMixing) is its one caller on the run path.
//
// Differences in form, never in result:
//   * `dot` accumulates in x87 extended precision upstream (`long double`); here it is the
//     exact software emulation in core::ext80 (gate/diff/ext80dot checks it against g++).
//   * The workspace bookkeeping (capacity counting against workspace_limit_bytes, the
//     bad_alloc catch, workspace_growths / reused_workspace) is reduced to the size checks:
//     the vectors are sized exactly as upstream sizes them, and a workspace over the limit
//     returns MemoryLimit as upstream does. Capacity history cannot change a solution.
//   * The matrix and preconditioner callables (std::function<bool(const double*, double*)>)
//     are closures over slices; an empty preconditioner is None.

use crate::core::cmath::{cmax, cmin};
use crate::core::ext80::dot_ext80;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GmresStatus {
    Converged,
    InvalidInput,
    MemoryLimit,
    NonFinite,
    OperatorFailure,
    PreconditionerFailure,
    Breakdown,
    Stagnated,
    IterationLimit,
}

#[derive(Debug, Clone, Copy)]
pub struct GmresControls {
    pub restart: i32,
    pub max_iterations: i32,
    pub stagnant_restarts: i32,
    pub absolute_tolerance: f64,
    pub relative_tolerance: f64,
    pub workspace_limit_bytes: usize,
}

impl Default for GmresControls {
    fn default() -> Self {
        GmresControls { restart: 40, max_iterations: 200, stagnant_restarts: 3, absolute_tolerance: 1.0e-12, relative_tolerance: 1.0e-10, workspace_limit_bytes: 64 * 1024 * 1024 }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct GmresResult {
    pub status: GmresStatus,
    pub iterations: i32,
    pub restarts: i32,
    pub matvecs: i32,
    pub preconditioner_applications: i32,
    pub reorthogonalizations: i32,
    pub true_residual: f64,
    pub estimated_residual: f64,
    pub tolerance: f64,
    pub happy_breakdown: bool,
}

impl Default for GmresResult {
    fn default() -> Self {
        GmresResult { status: GmresStatus::InvalidInput, iterations: 0, restarts: 0, matvecs: 0, preconditioner_applications: 0, reorthogonalizations: 0, true_residual: 0.0, estimated_residual: 0.0, tolerance: 0.0, happy_breakdown: false }
    }
}

impl GmresResult {
    pub fn converged(&self) -> bool {
        self.status == GmresStatus::Converged
    }
}

/// The GMRES workspace (NativeGmres).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct NativeGmres {
    basis_: Vec<f64>,
    hessenberg_: Vec<f64>,
    cosine_: Vec<f64>,
    sine_: Vec<f64>,
    g_: Vec<f64>,
    y_: Vec<f64>,
    x_: Vec<f64>,
    r_: Vec<f64>,
    w_: Vec<f64>,
    z_: Vec<f64>,
    update_: Vec<f64>,
    image_: Vec<f64>,
}

fn finite(x: &[f64]) -> bool {
    x.iter().all(|v| v.is_finite())
}

/// `norm`: hypot accumulated in order.
fn norm(x: &[f64]) -> f64 {
    let mut value = 0.0f64;
    for &v in x {
        value = value.hypot(v); // LIBM
    }
    value
}

pub type GmresApply<'a> = dyn FnMut(&[f64], &mut [f64]) -> bool + 'a;

impl NativeGmres {
    fn release(&mut self) {
        *self = NativeGmres::default();
    }

    /// solve(rhs, matrix, inverse, controls, output)
    pub fn solve(&mut self, b: &[f64], matrix: &mut GmresApply, mut inverse: Option<&mut GmresApply>, controls: &GmresControls, output: &mut Vec<f64>) -> GmresResult {
        let mut result = GmresResult::default();
        let n = b.len();
        if n == 0 || controls.restart <= 0 || controls.max_iterations < 0 || controls.stagnant_restarts <= 0 || !controls.absolute_tolerance.is_finite() || controls.absolute_tolerance < 0.0 || !controls.relative_tolerance.is_finite() || controls.relative_tolerance < 0.0 || (controls.absolute_tolerance == 0.0 && controls.relative_tolerance == 0.0) || !finite(b) {
            return result;
        }
        let m = cmin(n, cmin(controls.restart, cmax(1, controls.max_iterations)) as usize);
        let limit = controls.workspace_limit_bytes / std::mem::size_of::<f64>();
        let mut count: usize = 0;
        let mut add_product = |a: usize, c: usize| -> bool {
            if a != 0 && c > (limit - count) / a {
                return false;
            }
            count += a * c;
            true
        };
        if !add_product(m + 1, n) || !add_product(m + 1, m) || !add_product(6, n) || !add_product(4, m) || !add_product(1, 1) {
            result.status = GmresStatus::MemoryLimit;
            return result;
        }
        let total = (m + 1) * n + (m + 1) * m + 4 * m + 1 + 6 * n;
        if total * std::mem::size_of::<f64>() > controls.workspace_limit_bytes {
            self.release();
            result.status = GmresStatus::MemoryLimit;
            return result;
        }
        self.basis_.resize((m + 1) * n, 0.0);
        self.hessenberg_.resize((m + 1) * m, 0.0);
        self.cosine_.resize(m, 0.0);
        self.sine_.resize(m, 0.0);
        self.g_.resize(m + 1, 0.0);
        self.y_.resize(m, 0.0);
        for v in [&mut self.x_, &mut self.r_, &mut self.w_, &mut self.z_, &mut self.update_, &mut self.image_] {
            v.resize(n, 0.0);
        }

        self.x_.iter_mut().for_each(|v| *v = 0.0);
        let bnorm = norm(b);
        result.tolerance = controls.absolute_tolerance + controls.relative_tolerance * bnorm;
        if !bnorm.is_finite() || !result.tolerance.is_finite() {
            result.status = GmresStatus::NonFinite;
            return result;
        }

        // apply(in, out) and precondition(in, out) as functions of the workspace slices
        fn apply(matrix: &mut GmresApply, result: &mut GmresResult, inp: &[f64], out: &mut [f64]) -> bool {
            result.matvecs += 1;
            if !matrix(inp, out) {
                result.status = GmresStatus::OperatorFailure;
                return false;
            }
            if !finite(out) {
                result.status = GmresStatus::NonFinite;
                return false;
            }
            true
        }
        fn precondition(inverse: &mut Option<&mut GmresApply>, result: &mut GmresResult, inp: &[f64], out: &mut [f64]) -> bool {
            if let Some(inv) = inverse.as_mut() {
                result.preconditioner_applications += 1;
                if !inv(inp, out) {
                    result.status = GmresStatus::PreconditionerFailure;
                    return false;
                }
            } else {
                out.copy_from_slice(inp);
            }
            if !finite(out) {
                result.status = GmresStatus::NonFinite;
                return false;
            }
            true
        }
        // trueResidual()
        let true_residual = |this: &mut NativeGmres, matrix: &mut GmresApply, result: &mut GmresResult| -> bool {
            if !apply(matrix, result, &this.x_, &mut this.image_) {
                return false;
            }
            for i in 0..n {
                this.r_[i] = b[i] - this.image_[i];
            }
            result.true_residual = norm(&this.r_);
            if !result.true_residual.is_finite() {
                result.status = GmresStatus::NonFinite;
                return false;
            }
            true
        };

        if !true_residual(self, matrix, &mut result) {
            return result;
        }
        if result.true_residual <= result.tolerance {
            output.clone_from(&self.x_);
            result.status = GmresStatus::Converged;
            return result;
        }
        let mut stagnant = 0;
        let ld = m + 1;
        while result.iterations < controls.max_iterations {
            let old_residual = result.true_residual;
            self.hessenberg_.iter_mut().for_each(|v| *v = 0.0);
            self.g_.iter_mut().for_each(|v| *v = 0.0);
            self.g_[0] = old_residual;
            for i in 0..n {
                self.basis_[i] = self.r_[i] / old_residual;
            }
            let mut used = 0usize;
            let mut breakdown = false;
            let mut j = 0usize;
            while j < m && result.iterations < controls.max_iterations {
                if !precondition(&mut inverse, &mut result, &self.basis_[j * n..(j + 1) * n], &mut self.z_) || !apply(matrix, &mut result, &self.z_, &mut self.w_) {
                    return result;
                }
                let before = norm(&self.w_);
                if !before.is_finite() {
                    result.status = GmresStatus::NonFinite;
                    return result;
                }
                let hb = j * ld;
                for kk in 0..=j {
                    self.hessenberg_[hb + kk] = dot_ext80(&self.basis_[kk * n..(kk + 1) * n], &self.w_);
                    let hk = self.hessenberg_[hb + kk];
                    for i in 0..n {
                        self.w_[i] -= hk * self.basis_[kk * n + i];
                    }
                }
                let mut after = norm(&self.w_);
                if after < 0.7 * before {
                    result.reorthogonalizations += 1;
                    for kk in 0..=j {
                        let correction = dot_ext80(&self.basis_[kk * n..(kk + 1) * n], &self.w_);
                        self.hessenberg_[hb + kk] += correction;
                        for i in 0..n {
                            self.w_[i] -= correction * self.basis_[kk * n + i];
                        }
                    }
                    after = norm(&self.w_);
                }
                if !after.is_finite() {
                    result.status = GmresStatus::NonFinite;
                    return result;
                }
                breakdown = after <= 64.0 * f64::EPSILON * before;
                self.hessenberg_[hb + j + 1] = after;
                if !breakdown {
                    for i in 0..n {
                        self.basis_[(j + 1) * n + i] = self.w_[i] / after;
                    }
                }
                for kk in 0..j {
                    let a = self.hessenberg_[hb + kk];
                    let bnext = self.hessenberg_[hb + kk + 1];
                    self.hessenberg_[hb + kk] = self.cosine_[kk] * a + self.sine_[kk] * bnext;
                    self.hessenberg_[hb + kk + 1] = -self.sine_[kk] * a + self.cosine_[kk] * bnext;
                }
                let rotation = self.hessenberg_[hb + j].hypot(self.hessenberg_[hb + j + 1]); // LIBM
                result.iterations += 1;
                if !rotation.is_finite() || rotation == 0.0 {
                    result.status = GmresStatus::Breakdown;
                    return result;
                }
                self.cosine_[j] = self.hessenberg_[hb + j] / rotation;
                self.sine_[j] = self.hessenberg_[hb + j + 1] / rotation;
                self.hessenberg_[hb + j] = rotation;
                self.hessenberg_[hb + j + 1] = 0.0;
                let a = self.g_[j];
                self.g_[j] = self.cosine_[j] * a;
                self.g_[j + 1] = -self.sine_[j] * a;
                result.estimated_residual = self.g_[j + 1].abs();
                used = j + 1;
                if breakdown || result.estimated_residual <= result.tolerance {
                    break;
                }
                j += 1;
            }
            let mut i = used;
            while i > 0 {
                i -= 1;
                let mut value = self.g_[i];
                for jj in i + 1..used {
                    value -= self.hessenberg_[jj * ld + i] * self.y_[jj];
                }
                self.y_[i] = value / self.hessenberg_[i * ld + i];
                if !self.y_[i].is_finite() {
                    result.status = GmresStatus::Breakdown;
                    return result;
                }
            }
            self.update_.iter_mut().for_each(|v| *v = 0.0);
            for jj in 0..used {
                for i in 0..n {
                    self.update_[i] += self.basis_[jj * n + i] * self.y_[jj];
                }
            }
            if !finite(&self.update_) || !precondition(&mut inverse, &mut result, &self.update_, &mut self.z_) {
                if result.status == GmresStatus::InvalidInput {
                    result.status = GmresStatus::NonFinite;
                }
                return result;
            }
            for i in 0..n {
                self.x_[i] += self.z_[i];
            }
            if !finite(&self.x_) {
                result.status = GmresStatus::NonFinite;
                return result;
            }
            if !true_residual(self, matrix, &mut result) {
                return result;
            }
            if result.true_residual <= result.tolerance {
                output.clone_from(&self.x_);
                result.status = GmresStatus::Converged;
                result.happy_breakdown = breakdown;
                return result;
            }
            if breakdown {
                result.status = GmresStatus::Breakdown;
                return result;
            }
            if result.true_residual >= old_residual * (1.0 - 1.0e-3) {
                stagnant += 1;
            } else {
                stagnant = 0;
            }
            if stagnant >= controls.stagnant_restarts {
                result.status = GmresStatus::Stagnated;
                return result;
            }
            if result.iterations < controls.max_iterations {
                result.restarts += 1;
            }
        }
        result.status = GmresStatus::IterationLimit;
        result
    }
}
