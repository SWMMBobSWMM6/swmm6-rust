// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/math/FindRoot.hpp and src/engine/math/FindRoot.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Newton-Raphson and Ridder root finders (legacy findroot.c, op for op). The evaluation
// SEQUENCE is part of the result: the culvert's Form-1 flow reads the last evaluation of its
// equation, not the root.
//
// Differences in form, never in result:
//   * `std::function` callbacks are generic `FnMut` closures. Upstream's callbacks are
//     const-callable lambdas that mutate through captured pointers; FnMut says the same.
//   * newton's `double* rts` is `&mut f64`, and its callback's two out-pointers are
//     `&mut f64`. Upstream leaves f and df uninitialized before the first call (the
//     callback must fill both); here they start at 0.0.

/// Maximum iterations
pub const MAXIT: i32 = 60;

/// Find a root with Newton steps and a bisection fallback. f(x1) < 0 < f(x2) is the
/// caller's job; returns the number of function evaluations, 0 past MAXIT. `rts` holds
/// the initial guess and receives the final estimate.
pub fn newton<F>(x1: f64, x2: f64, rts: &mut f64, xacc: f64, mut func: F) -> i32
where
    F: FnMut(f64, &mut f64, &mut f64),
{
    // legacy findroot_Newton, op for op: the caller brackets the root with
    // f(x1) < 0 < f(x2) (switching x1/x2 when f(x1) > f(x2), as kinwave.c
    // solveContinuity does) — no orienting pre-evaluations here, so the
    // function is called exactly in legacy's sequence and count; the
    // return is the number of evaluations, 0 past MAXIT.
    let mut n: i32 = 0;
    let mut f: f64 = 0.0;
    let mut df: f64 = 0.0;
    let mut dx: f64;
    let mut dxold: f64;
    let mut temp: f64;
    let mut x = *rts;
    let mut xlo = x1;
    let mut xhi = x2;
    dxold = (x2 - x1).abs();
    dx = dxold;
    func(x, &mut f, &mut df);
    n += 1;

    let mut j = 1;
    while j <= MAXIT {
        // Bisect if Newton out of range or not decreasing fast enough.
        if (((x - xhi) * df - f) * ((x - xlo) * df - f) >= 0.0)
            || ((2.0 * f).abs() > (dxold * df).abs())
        {
            dxold = dx;
            dx = 0.5 * (xhi - xlo);
            x = xlo + dx;
            if xlo == x {
                break;
            }
        }
        // Newton step acceptable. Take it.
        else {
            dxold = dx;
            dx = f / df;
            temp = x;
            x -= dx;
            if temp == x {
                break;
            }
        }

        // Convergence criterion.
        if dx.abs() < xacc {
            break;
        }

        // Evaluate function. Maintain bracket on the root.
        func(x, &mut f, &mut df);
        n += 1;
        if f < 0.0 {
            xlo = x;
        } else {
            xhi = x;
        }
        j += 1;
    }
    *rts = x;
    if n <= MAXIT {
        n
    } else {
        0
    }
}

/// Estimate a bracketed root with derivative-free Ridder iterations; -1.0e20 when the
/// endpoints do not straddle a root.
pub fn ridder<F>(x1: f64, x2: f64, xacc: f64, mut func: F) -> f64
where
    F: FnMut(f64) -> f64,
{
    // legacy findroot_Ridder, op for op — including the evaluation SEQUENCE,
    // which the culvert's Form-1 flow depends on (it reads the last
    // form1Eqn evaluation, not the root), the initial estimate 0.5*(x1+x2),
    // the return of the PREVIOUS estimate when the new one is within xacc
    // (`break` then `return ans`), the SIGN-based bracket updates and the
    // -1e20 no-bracket return.
    let sign = |x: f64, y: f64| -> f64 {
        if y >= 0.0 {
            x.abs()
        } else {
            -x.abs()
        }
    };
    let mut flo = func(x1);
    let mut fhi = func(x2);
    if flo == 0.0 {
        return x1;
    }
    if fhi == 0.0 {
        return x2;
    }
    let mut ans = 0.5 * (x1 + x2);
    if (flo > 0.0 && fhi < 0.0) || (flo < 0.0 && fhi > 0.0) {
        let mut xlo = x1;
        let mut xhi = x2;
        let mut j = 1;
        while j <= MAXIT {
            let xm = 0.5 * (xlo + xhi);
            let fm = func(xm);
            let s = (fm * fm - flo * fhi).sqrt();
            if s == 0.0 {
                return ans;
            }
            let xnew = xm + (xm - xlo) * ((if flo >= fhi { 1.0 } else { -1.0 }) * fm / s);
            if (xnew - ans).abs() <= xacc {
                break;
            }
            ans = xnew;
            let fnew = func(ans);
            if sign(fm, fnew) != fm {
                xlo = xm;
                flo = fm;
                xhi = ans;
                fhi = fnew;
            } else if sign(flo, fnew) != flo {
                xhi = ans;
                fhi = fnew;
            } else if sign(fhi, fnew) != fhi {
                xlo = ans;
                flo = fnew;
            } else {
                return ans;
            }
            if (xhi - xlo).abs() <= xacc {
                return ans;
            }
            j += 1;
        }
        return ans;
    }
    -1.0e20
}
