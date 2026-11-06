// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/math/OdeSolver.hpp, OdeSolverDetail.hpp and
// OdeSolver.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The RK45 Cash-Karp integrator (legacy odesolve.c): integrate, rkqs, rkck.
//
// Differences of form, never in result:
//   * Upstream keeps a thread_local workspace of calloc'ed arrays that grows with n. Every
//     element integrate() reads is written first in the same call (y from ystart, dydx and
//     yscal at each step, the stage arrays in rkck), so the workspace carries nothing from
//     one call to the next. Here it is a stack workspace sized by the const generic N (the
//     number of equations), and ODE_TOO_MANY cannot arise.
//   * The derivative callback is a closure `FnMut(x, &[f64; N], &mut [f64; N])`.
//   * integrate_batch_2eq (one explicit Euler step) is not on the 1D run path the runoff
//     step reaches; NOT YET TRANSLATED.

/// Maximum integration steps
pub const MAXSTP: i32 = 10000;
/// Underflow protection
pub const TINY: f64 = 1.0e-30;
/// Step adjustment safety factor
pub const SAFETY: f64 = 0.9;
/// Exponent for step increase
pub const PGROW: f64 = -0.2;
/// Exponent for step decrease
pub const PSHRNK: f64 = -0.25;
/// = (5/SAFETY)^(1/PGROW)
pub const ERRCON: f64 = 1.89e-4;

pub const ODE_OK: i32 = 0;
pub const ODE_TOO_MANY: i32 = 1;
pub const ODE_UNDERFLOW: i32 = 2;
pub const ODE_MAX_STEPS: i32 = 3;

// Cash-Karp RK45 coefficients
const A2: f64 = 0.2;
const A3: f64 = 0.3;
const A4: f64 = 0.6;
const A5: f64 = 1.0;
const A6: f64 = 0.875;
const B21: f64 = 0.2;
const B31: f64 = 3.0 / 40.0;
const B32: f64 = 9.0 / 40.0;
const B41: f64 = 0.3;
const B42: f64 = -0.9;
const B43: f64 = 1.2;
const B51: f64 = -11.0 / 54.0;
const B52: f64 = 2.5;
const B53: f64 = -70.0 / 27.0;
const B54: f64 = 35.0 / 27.0;
const B61: f64 = 1631.0 / 55296.0;
const B62: f64 = 175.0 / 512.0;
const B63: f64 = 575.0 / 13824.0;
const B64: f64 = 44275.0 / 110592.0;
const B65: f64 = 253.0 / 4096.0;
const C1: f64 = 37.0 / 378.0;
const C3: f64 = 250.0 / 621.0;
const C4: f64 = 125.0 / 594.0;
const C6: f64 = 512.0 / 1771.0;
const DC1: f64 = C1 - 2825.0 / 27648.0;
const DC3: f64 = C3 - 18575.0 / 48384.0;
const DC4: f64 = C4 - 13525.0 / 55296.0;
const DC5: f64 = -277.0 / 14336.0;
const DC6: f64 = C6 - 0.25;

struct Ws<const N: usize> {
    y: [f64; N],
    yscal: [f64; N],
    dydx: [f64; N],
    yerr: [f64; N],
    ytemp: [f64; N],
    ak2: [f64; N],
    ak3: [f64; N],
    ak4: [f64; N],
    ak5: [f64; N],
    ak6: [f64; N],
}

fn rkck<const N: usize, F: FnMut(f64, &[f64; N], &mut [f64; N])>(ws: &mut Ws<N>, x: f64, h: f64, derivs: &mut F) {
    for i in 0..N {
        ws.ytemp[i] = ws.y[i] + B21 * h * ws.dydx[i];
    }
    derivs(x + A2 * h, &ws.ytemp, &mut ws.ak2);

    for i in 0..N {
        ws.ytemp[i] = ws.y[i] + h * (B31 * ws.dydx[i] + B32 * ws.ak2[i]);
    }
    derivs(x + A3 * h, &ws.ytemp, &mut ws.ak3);

    for i in 0..N {
        ws.ytemp[i] = ws.y[i] + h * (B41 * ws.dydx[i] + B42 * ws.ak2[i] + B43 * ws.ak3[i]);
    }
    derivs(x + A4 * h, &ws.ytemp, &mut ws.ak4);

    for i in 0..N {
        ws.ytemp[i] = ws.y[i] + h * (B51 * ws.dydx[i] + B52 * ws.ak2[i] + B53 * ws.ak3[i] + B54 * ws.ak4[i]);
    }
    derivs(x + A5 * h, &ws.ytemp, &mut ws.ak5);

    for i in 0..N {
        ws.ytemp[i] = ws.y[i] + h * (B61 * ws.dydx[i] + B62 * ws.ak2[i] + B63 * ws.ak3[i] + B64 * ws.ak4[i] + B65 * ws.ak5[i]);
    }
    derivs(x + A6 * h, &ws.ytemp, &mut ws.ak6);

    for i in 0..N {
        ws.ytemp[i] = ws.y[i] + h * (C1 * ws.dydx[i] + C3 * ws.ak3[i] + C4 * ws.ak4[i] + C6 * ws.ak6[i]);
    }

    for i in 0..N {
        ws.yerr[i] = h * (DC1 * ws.dydx[i] + DC3 * ws.ak3[i] + DC4 * ws.ak4[i] + DC5 * ws.ak5[i] + DC6 * ws.ak6[i]);
    }
}

#[allow(clippy::too_many_arguments)]
fn rkqs<const N: usize, F: FnMut(f64, &[f64; N], &mut [f64; N])>(ws: &mut Ws<N>, x: &mut f64, htry: f64, eps: f64, hdid: &mut f64, hnext: &mut f64, derivs: &mut F) -> i32 {
    let xold = *x;
    let mut h = htry;
    loop {
        rkck(ws, xold, h, derivs);

        let mut errmax = 0.0;
        for i in 0..N {
            let err = (ws.yerr[i] / ws.yscal[i]).abs();
            if err > errmax {
                errmax = err;
            }
        }
        errmax /= eps;

        if errmax > 1.0 {
            let htemp = SAFETY * h * errmax.powf(PSHRNK); // LIBM
            if h >= 0.0 {
                if htemp > 0.1 * h {
                    h = htemp;
                } else {
                    h = 0.1 * h;
                }
            } else if htemp < 0.1 * h {
                h = htemp;
            } else {
                h = 0.1 * h;
            }
            let xnew = xold + h;
            if xnew == xold {
                return ODE_UNDERFLOW;
            }
            continue;
        } else {
            if errmax > ERRCON {
                *hnext = SAFETY * h * errmax.powf(PGROW); // LIBM
            } else {
                *hnext = 5.0 * h;
            }
            *hdid = h;
            *x += h;
            ws.y = ws.ytemp;
            return ODE_OK;
        }
    }
}

/// ode::integrate (detail::integrate): integrate the N equations `ystart` from x1 to x2
/// with tolerance eps and first step h1. ystart is written only on success (ODE_OK).
pub fn integrate<const N: usize, F: FnMut(f64, &[f64; N], &mut [f64; N])>(ystart: &mut [f64; N], x1: f64, x2: f64, eps: f64, h1: f64, mut derivs: F) -> i32 {
    let mut ws = Ws::<N> { y: [0.0; N], yscal: [0.0; N], dydx: [0.0; N], yerr: [0.0; N], ytemp: [0.0; N], ak2: [0.0; N], ak3: [0.0; N], ak4: [0.0; N], ak5: [0.0; N], ak6: [0.0; N] };
    let mut hdid = 0.0;
    let mut hnext = 0.0;
    let mut x = x1;
    let mut h = h1;

    ws.y = *ystart;

    for _nstp in 1..=MAXSTP {
        derivs(x, &ws.y, &mut ws.dydx);
        for i in 0..N {
            ws.yscal[i] = ws.y[i].abs() + (ws.dydx[i] * h).abs() + TINY;
        }
        if (x + h - x2) * (x + h - x1) > 0.0 {
            h = x2 - x;
        }

        let errcode = rkqs(&mut ws, &mut x, h, eps, &mut hdid, &mut hnext, &mut derivs);
        if errcode != 0 {
            break;
        }

        if (x - x2) * (x2 - x1) >= 0.0 {
            *ystart = ws.y;
            return ODE_OK;
        }
        if hnext.abs() <= 0.0 {
            return ODE_UNDERFLOW;
        }
        h = hnext;
    }
    ODE_MAX_STEPS
}
