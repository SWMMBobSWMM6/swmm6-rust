// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/fv/FvClosureKernels.hpp,
// FvKernels.hpp and SectionGeometry.hpp (HydroCouple/openswmm.engine, branch swmm6_rel,
// commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the Apache License,
// Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The FV section closure (Hermite area/width/I1 tables and their inverse), the slot
// taper, and the exact section geometry (secgeom) the closure tables are sampled from.
//
// The step kernels (5 November 2026, with the FV step): closureEval with the radius output,
// closureI1, i1OfDepth, celerity, waveSpeeds, physicalFlux, riemannMassFlux, riemannFlux,
// frictionUpdate, localLossUpdate, ufUpdate, the tpa* kernels, faceCflDt and positivityScale
// (FaceState and FaceFlux are in explicit_solver). NOT TRANSLATED: speciesFlux (species
// transport on the FV mesh; no gate deck carries FV species).
//
// Differences in form, never in result:
//   * `static_cast<int>(double)` is cvt_i32 (INT_MIN for NaN / out of range, as x86).
//   * closureHermite's three `double&` outputs are a returned tuple.
//   * secgeom's templated adaptive Simpson takes the integrand as a closure; the
//     recursion and every operation keep upstream's order.

use crate::core::cmath::cvt_i32;
use crate::hydraulics::fv::{FvClosure, FvGeometry, K_CLOSURE_PANELS, K_CLOSURE_POLYNOMIAL};
use crate::hydraulics::xsect_batch::{XSectParams, XSectShape};
use crate::hydraulics::xsection as xsect;

pub const K_GRAVITY: f64 = 32.2;
pub const K_DRY_DEPTH: f64 = 1.0e-7;
pub const K_DRY_AREA: f64 = 1.0e-12;
pub const K_ETA_DEADBAND: f64 = 1.0e-12;

const NP: i32 = K_CLOSURE_PANELS as i32;

// ===========================================================================
// Closure kernels (FvClosureKernels.hpp)
// ===========================================================================

/// A, T and I1 at one depth, from one panel locate.
#[derive(Debug, Clone, Copy, Default)]
pub struct ClosureEval {
    pub a: f64,
    pub t: f64,
    pub i1: f64,
}

pub fn closure_locate(c: &FvClosure, h: f64) -> (i32, f64) {
    let x = h * c.inv_dh;
    let mut i = cvt_i32(x);
    if i < 0 {
        i = 0;
    }
    if i > NP - 1 {
        i = NP - 1;
    }
    let mut t = x - i as f64;
    if t < 0.0 {
        t = 0.0;
    }
    if t > 1.0 {
        t = 1.0;
    }
    (i, t)
}

/// (a, dadh, int_dh) of panel i at panel coordinate t.
pub fn closure_hermite(c: &FvClosure, i: i32, t: f64) -> (f64, f64, f64) {
    let ui = i as usize;
    let a0 = c.A[ui];
    let a1 = c.A[ui + 1];
    let m0 = c.M[ui] * c.dh;
    let m1 = c.M[ui + 1] * c.dh;
    let t2 = t * t;
    let t3 = t2 * t;
    let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h10 = t3 - 2.0 * t2 + t;
    let h01 = -2.0 * t3 + 3.0 * t2;
    let h11 = t3 - t2;
    let a = h00 * a0 + h10 * m0 + h01 * a1 + h11 * m1;
    let d00 = 6.0 * t2 - 6.0 * t;
    let d10 = 3.0 * t2 - 4.0 * t + 1.0;
    let d01 = -6.0 * t2 + 6.0 * t;
    let d11 = 3.0 * t2 - 2.0 * t;
    let dadh = (d00 * a0 + d10 * m0 + d01 * a1 + d11 * m1) * c.inv_dh;
    let t4 = t2 * t2;
    let i00 = 0.5 * t4 - t3 + t;
    let i10 = 0.25 * t4 - (2.0 / 3.0) * t3 + 0.5 * t2;
    let i01 = -0.5 * t4 + t3;
    let i11 = 0.25 * t4 - t3 / 3.0;
    let int_dh = (i00 * a0 + i10 * m0 + i01 * a1 + i11 * m1) * c.dh;
    (a, dadh, int_dh)
}

pub fn closure_eval(c: &FvClosure, h: f64) -> ClosureEval {
    let mut e = ClosureEval::default();
    if h <= 0.0 {
        return e;
    }
    if h >= c.y_full {
        let d = h - c.y_full;
        e.a = c.a_crown + c.t_slot * d;
        e.t = c.t_slot;
        e.i1 = c.i1_crown + c.a_crown * d + 0.5 * c.t_slot * d * d;
        return e;
    }
    if c.kind == K_CLOSURE_POLYNOMIAL {
        e.a = (c.c2 * h + c.c1) * h;
        e.t = 2.0 * c.c2 * h + c.c1;
        e.i1 = (c.c2 * h / 3.0 + 0.5 * c.c1) * h * h;
        return e;
    }
    let (i, t) = closure_locate(c, h);
    let (a, dadh, idh) = closure_hermite(c, i, t);
    e.a = a;
    e.t = if dadh > 0.0 { dadh } else { 0.0 };
    e.i1 = c.I1[i as usize] + idh;
    e
}

/// closureEval(c, h, radius): the same evaluation, with the hydraulic radius as well.
pub fn closure_eval_r(c: &FvClosure, h: f64) -> (ClosureEval, f64) {
    let mut e = ClosureEval::default();
    let mut radius = 0.0;
    if h <= 0.0 {
        return (e, radius);
    }
    if h >= c.y_full {
        radius = c.r_full;
        let d = h - c.y_full;
        e.a = c.a_crown + c.t_slot * d;
        e.t = c.t_slot;
        e.i1 = c.i1_crown + c.a_crown * d + 0.5 * c.t_slot * d * d;
        return (e, radius);
    }
    if c.kind == K_CLOSURE_POLYNOMIAL {
        e.a = (c.c2 * h + c.c1) * h;
        let p = c.p0 + c.p1 * h;
        radius = if p > 0.0 { (c.c2 * h + c.c1) * h / p } else { 0.0 };
        e.t = 2.0 * c.c2 * h + c.c1;
        e.i1 = (c.c2 * h / 3.0 + 0.5 * c.c1) * h * h;
        return (e, radius);
    }
    let (i, t) = closure_locate(c, h);
    let ui = i as usize;
    radius = c.R[ui] + t * (c.R[ui + 1] - c.R[ui]);
    let (a, dadh, idh) = closure_hermite(c, i, t);
    e.a = a;
    e.t = if dadh > 0.0 { dadh } else { 0.0 };
    e.i1 = c.I1[ui] + idh;
    (e, radius)
}

pub fn closure_i1(c: &FvClosure, h: f64) -> f64 {
    closure_eval(c, h).i1
}

pub fn closure_area(c: &FvClosure, h: f64) -> f64 {
    if h <= 0.0 {
        return 0.0;
    }
    if h >= c.y_full {
        return c.a_crown + c.t_slot * (h - c.y_full);
    }
    if c.kind == K_CLOSURE_POLYNOMIAL {
        return (c.c2 * h + c.c1) * h;
    }
    let (i, t) = closure_locate(c, h);
    closure_hermite(c, i, t).0
}

pub fn closure_width(c: &FvClosure, h: f64) -> f64 {
    if h <= 0.0 {
        return 0.0;
    }
    if h >= c.y_full {
        return c.t_slot;
    }
    if c.kind == K_CLOSURE_POLYNOMIAL {
        return 2.0 * c.c2 * h + c.c1;
    }
    let (i, t) = closure_locate(c, h);
    let dadh = closure_hermite(c, i, t).1;
    if dadh > 0.0 {
        dadh
    } else {
        0.0
    }
}

pub fn closure_hyd_rad(c: &FvClosure, h: f64) -> f64 {
    if h <= 0.0 {
        return 0.0;
    }
    if h >= c.y_full {
        return c.r_full;
    }
    if c.kind == K_CLOSURE_POLYNOMIAL {
        let p = c.p0 + c.p1 * h;
        return if p > 0.0 { (c.c2 * h + c.c1) * h / p } else { 0.0 };
    }
    let (i, t) = closure_locate(c, h);
    let ui = i as usize;
    c.R[ui] + t * (c.R[ui + 1] - c.R[ui])
}

pub fn closure_depth_of_area(c: &FvClosure, a: f64) -> f64 {
    if a <= 0.0 {
        return 0.0;
    }
    if a >= c.a_crown {
        return c.y_full + (a - c.a_crown) / c.t_slot;
    }
    if c.kind == K_CLOSURE_POLYNOMIAL {
        if c.c2 <= 0.0 {
            return if c.c1 > 0.0 { a / c.c1 } else { 0.0 };
        }
        return 2.0 * a / (c.c1 + (c.c1 * c.c1 + 4.0 * c.c2 * a).sqrt());
    }

    let mut k = cvt_i32(a * c.inv_a_crown_n);
    if k < 0 {
        k = 0;
    }
    if k > NP {
        k = NP;
    }
    let mut i = c.j_of_a[k as usize] as i32;
    while i > 0 && c.A[i as usize] > a {
        i -= 1;
    }
    while i < NP - 1 && c.A[(i + 1) as usize] < a {
        i += 1;
    }

    let a0 = c.A[i as usize];
    let a1 = c.A[(i + 1) as usize];
    if a <= a0 {
        return i as f64 * c.dh;
    }
    if a >= a1 {
        return (i + 1) as f64 * c.dh;
    }

    let tol = 1.0e-15 * c.y_full * c.inv_dh; // in units of t
    let mut lo = 0.0;
    let mut hi = 1.0;
    let frac = (a - a0) / (a1 - a0);
    let mut t = if c.M[i as usize] > 0.0 { frac } else { frac.sqrt() };
    for _it in 0..12 {
        let (mut f, dadh, _idh) = closure_hermite(c, i, t);
        f -= a;
        if f == 0.0 {
            break;
        }
        if f > 0.0 {
            hi = t;
        } else {
            lo = t;
        }
        let dadt = dadh * c.dh;
        let mut step = if dadt > 0.0 { -f / dadt } else { 0.0 };
        if dadt > 0.0 && step.abs() <= tol {
            break;
        }
        let mut tn = t + step;
        if !(dadt > 0.0) || !(tn > lo && tn < hi) {
            tn = 0.5 * (lo + hi);
        }
        step = tn - t;
        t = tn;
        if step.abs() <= tol {
            break;
        }
    }
    i as f64 * c.dh + t * c.dh
}

// ===========================================================================
// FvKernels.hpp (the closure entry points and the slot taper)
// ===========================================================================

pub fn slot_ramp(s: f64) -> f64 {
    if s <= 0.0 {
        return 0.0;
    }
    if s >= 1.0 {
        return 1.0;
    }
    s * s * (3.0 - 2.0 * s)
}

pub fn slot_ramp_integral(s: f64) -> f64 {
    if s <= 0.0 {
        return 0.0;
    }
    if s >= 1.0 {
        return s - 0.5;
    }
    let s3 = s * s * s;
    s3 - 0.5 * s3 * s
}

pub fn area_of_depth(g: &FvGeometry, h: f64) -> f64 {
    closure_area(&g.closure_tbl, h)
}

pub fn width_of_depth(g: &FvGeometry, h: f64) -> f64 {
    closure_width(&g.closure_tbl, h)
}

pub fn hyd_rad_of_depth(g: &FvGeometry, h: f64) -> f64 {
    closure_hyd_rad(&g.closure_tbl, h)
}

pub fn depth_of_area(g: &FvGeometry, a: f64) -> f64 {
    closure_depth_of_area(&g.closure_tbl, a)
}

pub fn i1_of_depth(g: &FvGeometry, h: f64, _area_at_h: f64) -> f64 {
    closure_i1(&g.closure_tbl, h)
}

pub fn celerity(a: f64, t: f64) -> f64 {
    if a <= K_DRY_AREA || t <= 0.0 {
        return 0.0;
    }
    (K_GRAVITY * a / t).sqrt()
}

use crate::core::cmath::{cmax, cmin};
use crate::hydraulics::fv::explicit_solver::{FaceFlux, FaceState};

/// waveSpeeds(L, R, sl, sr), returned as (sl, sr).
pub fn wave_speeds(l: &FaceState, r: &FaceState) -> (f64, f64) {
    let wet_l = l.a > K_DRY_AREA;
    let wet_r = r.a > K_DRY_AREA;
    let (mut sl, mut sr);
    if wet_l && wet_r {
        sl = cmin(l.u - l.c, r.u - r.c);
        sr = cmax(l.u + l.c, r.u + r.c);
        if l.press != r.press && l.press < 2 && r.press < 2 {
            if l.press != 0 {
                let da = l.a - r.a;
                let w = if da > 0.0 { (l.q - r.q) / da } else { sr };
                sr = cmax(r.u + r.c, cmin(sr, w));
            } else {
                let da = r.a - l.a;
                let w = if da > 0.0 { (r.q - l.q) / da } else { sl };
                sl = cmin(l.u - l.c, cmax(sl, w));
            }
        }
    } else if wet_r {
        sl = r.u - 2.0 * r.c;
        sr = r.u + r.c;
    } else if wet_l {
        sl = l.u - l.c;
        sr = l.u + 2.0 * l.c;
    } else {
        sl = 0.0;
        sr = 0.0;
    }
    (sl, sr)
}

/// physicalFlux(S, fa, fq), returned as (fa, fq).
pub fn physical_flux(s: &FaceState) -> (f64, f64) {
    (s.q, s.q * s.u + K_GRAVITY * s.i1)
}

pub fn riemann_mass_flux(l: &FaceState, r: &FaceState) -> f64 {
    if l.a <= K_DRY_AREA && r.a <= K_DRY_AREA {
        return 0.0;
    }
    let (sl, sr) = wave_speeds(l, r);
    let (fal, _fql) = physical_flux(l);
    let (far, _fqr) = physical_flux(r);
    if sl >= 0.0 {
        return fal;
    }
    if sr <= 0.0 {
        return far;
    }
    (sr * fal - sl * far + sl * sr * (r.a - l.a)) / (sr - sl)
}

pub fn riemann_flux(l: &FaceState, r: &FaceState) -> FaceFlux {
    let mut out = FaceFlux::default();
    if l.a <= K_DRY_AREA && r.a <= K_DRY_AREA {
        return out;
    }
    let (sl, sr) = wave_speeds(l, r);
    out.sl = sl;
    out.sr = sr;
    let (fal, fql) = physical_flux(l);
    let (far, fqr) = physical_flux(r);
    if sl >= 0.0 {
        out.mass = fal;
        out.mom = fql;
        out.sstar = sl;
        return out;
    }
    if sr <= 0.0 {
        out.mass = far;
        out.mom = fqr;
        out.sstar = sr;
        return out;
    }
    let dsr = sr - sl;
    out.mass = (sr * fal - sl * far + sl * sr * (r.a - l.a)) / dsr;
    out.mom = (sr * fql - sl * fqr + sl * sr * (r.q - l.q)) / dsr;
    let den = r.a * (r.u - sr) - l.a * (l.u - sl);
    if den.abs() > 1.0e-14 {
        out.sstar = (sl * r.a * (r.u - sr) - sr * l.a * (l.u - sl)) / den;
    } else {
        let a_hll = (sr * r.a - sl * l.a - (far - fal)) / dsr;
        out.sstar = if a_hll > K_DRY_AREA { out.mass / a_hll } else { 0.0 };
    }
    out
}

pub fn friction_update(q: f64, u: f64, r: f64, dt: f64, rough_factor: f64) -> f64 {
    if r <= 0.0 || rough_factor <= 0.0 {
        return q;
    }
    let r43 = r * crate::core::cmath::glibc_cbrt(r); // LIBM (std::cbrt: glibc's, not Rust's f64::cbrt)
    let den = 1.0 + dt * rough_factor * u.abs() / r43;
    q / den
}

pub fn local_loss_update(q: f64, u: f64, k: f64, dx: f64, dt: f64) -> f64 {
    if k <= 0.0 || dx <= 0.0 {
        return q;
    }
    q / (1.0 + dt * k * u.abs() / (2.0 * dx))
}

pub fn uf_update(q: f64, a: f64, u_old: f64, k3: f64, grad_term: f64, dt: f64) -> f64 {
    if k3 <= 0.0 || a <= 0.0 {
        return q;
    }
    const K_UF_VEL_FLOOR: f64 = 0.01;
    if u_old.abs() < K_UF_VEL_FLOOR && (q / a).abs() < K_UF_VEL_FLOOR {
        return q;
    }
    let mut qn = (q + k3 * a * u_old) / (1.0 + k3);
    qn -= dt * k3 * a * grad_term;
    let dq = qn - q;
    let cap = 0.5 * q.abs();
    if dq.abs() > cap {
        qn = q + if dq > 0.0 { cap } else { -cap };
    }
    qn
}

pub fn tpa_depth_of_area(g: &FvGeometry, a: f64) -> f64 {
    g.y_full + (a - g.a_crown) / g.t_slot
}

pub fn tpa_area_of_depth(g: &FvGeometry, h: f64) -> f64 {
    g.a_crown + g.t_slot * (h - g.y_full)
}

pub fn tpa_i1_of_depth(g: &FvGeometry, h: f64) -> f64 {
    g.i1_crown + g.a_crown * (h - g.y_full)
}

pub fn face_cfl_dt(cfl: f64, dx: f64, u: f64, c: f64) -> f64 {
    let s = u.abs() + c;
    if s > 1.0e-12 { cfl * dx / s } else { 1.0e30 }
}

pub fn positivity_scale(vol: f64, outflow: f64, dt: f64) -> f64 {
    if outflow <= 0.0 || dt <= 0.0 {
        return 1.0;
    }
    let take = outflow * dt;
    let avail = if vol > 0.0 { vol } else { 0.0 };
    if take <= avail { 1.0 } else { avail / take }
}

// ===========================================================================
// SectionGeometry.hpp (secgeom)
// ===========================================================================

pub mod secgeom {
    use super::*;

    pub const K_PI: f64 = std::f64::consts::PI; // 3.14159265358979323846
    pub const K_RECT_ALF_MAX: f64 = 0.97;

    fn circ_theta(r: f64, y: f64) -> f64 {
        let mut c = 1.0 - y / r;
        if c > 1.0 {
            c = 1.0;
        }
        if c < -1.0 {
            c = -1.0;
        }
        2.0 * c.acos() // LIBM
    }
    fn circ_area(r: f64, y: f64) -> f64 {
        if y <= 0.0 {
            return 0.0;
        }
        if y >= 2.0 * r {
            return K_PI * r * r;
        }
        let th = circ_theta(r, y);
        0.5 * r * r * (th - th.sin()) // LIBM
    }
    fn circ_width(r: f64, y: f64) -> f64 {
        if y <= 0.0 || y >= 2.0 * r {
            return 0.0;
        }
        let s = y * (2.0 * r - y);
        if s > 0.0 {
            2.0 * s.sqrt()
        } else {
            0.0
        }
    }
    fn circ_perim(r: f64, y: f64) -> f64 {
        if y <= 0.0 {
            return 0.0;
        }
        if y >= 2.0 * r {
            return 2.0 * K_PI * r;
        }
        r * circ_theta(r, y)
    }

    #[allow(clippy::too_many_arguments)]
    fn simpson_step(f: &dyn Fn(f64) -> f64, a: f64, b: f64, fa: f64, fm: f64, fb: f64, whole: f64, eps: f64, depth: i32) -> f64 {
        let m = 0.5 * (a + b);
        let lm = 0.5 * (a + m);
        let rm = 0.5 * (m + b);
        let flm = f(lm);
        let frm = f(rm);
        let left = (m - a) / 6.0 * (fa + 4.0 * flm + fm);
        let right = (b - m) / 6.0 * (fm + 4.0 * frm + fb);
        let delta = left + right - whole;
        if depth <= 0 || delta.abs() <= 15.0 * eps {
            return left + right + delta / 15.0;
        }
        simpson_step(f, a, m, fa, flm, fm, left, 0.5 * eps, depth - 1) + simpson_step(f, m, b, fm, frm, fb, right, 0.5 * eps, depth - 1)
    }

    fn integrate(f: &dyn Fn(f64) -> f64, a: f64, b: f64, eps: f64, depth: i32) -> f64 {
        if !(b > a) {
            return 0.0;
        }
        let fa = f(a);
        let fb = f(b);
        let fm = f(0.5 * (a + b));
        let whole = (b - a) / 6.0 * (fa + 4.0 * fm + fb);
        simpson_step(f, a, b, fa, fm, fb, whole, eps, depth)
    }

    #[allow(dead_code)]
    fn ellipse_perim(ax: f64, by: f64, y: f64) -> f64 {
        if y <= 0.0 {
            return 0.0;
        }
        let mut c = 1.0 - y / by;
        if c < -1.0 {
            c = -1.0;
        }
        if c > 1.0 {
            c = 1.0;
        }
        let phi_max = c.acos(); // LIBM
        let ds = |phi: f64| -> f64 {
            let cx = ax * phi.cos(); // LIBM
            let sy = by * phi.sin(); // LIBM
            (cx * cx + sy * sy).sqrt()
        };
        2.0 * integrate(&ds, 0.0, phi_max, 1.0e-12 * (ax + by), 40)
    }

    fn shape(t: i32) -> Option<XSectShape> {
        XSectShape::from_code(t)
    }

    pub fn has_closed_form(type_: i32) -> bool {
        matches!(
            shape(type_),
            Some(XSectShape::Circular)
                | Some(XSectShape::ForceMain)
                | Some(XSectShape::FilledCircular)
                | Some(XSectShape::RectClosed)
                | Some(XSectShape::RectOpen)
                | Some(XSectShape::Trapezoidal)
                | Some(XSectShape::Triangular)
                | Some(XSectShape::Parabolic)
                | Some(XSectShape::Powerfunc)
        )
    }

    pub fn is_polynomial_open(type_: i32) -> bool {
        matches!(shape(type_), Some(XSectShape::RectOpen) | Some(XSectShape::Trapezoidal) | Some(XSectShape::Triangular))
    }

    pub fn area_of_depth(xs: &XSectParams, h: f64) -> f64 {
        let mut h = h;
        if h <= 0.0 {
            return 0.0;
        }
        if h > xs.y_full {
            h = xs.y_full;
        }
        match shape(xs.type_) {
            Some(XSectShape::Circular) | Some(XSectShape::ForceMain) => circ_area(0.5 * xs.y_full, h),
            Some(XSectShape::FilledCircular) => {
                let r = 0.5 * (xs.y_full + xs.y_bot);
                circ_area(r, xs.y_bot + h) - circ_area(r, xs.y_bot)
            }
            Some(XSectShape::RectClosed) | Some(XSectShape::RectOpen) => h * xs.w_max,
            Some(XSectShape::Trapezoidal) => (xs.y_bot + xs.s_bot * h) * h,
            Some(XSectShape::Triangular) => xs.s_bot * h * h,
            Some(XSectShape::Parabolic) => (4.0 / 3.0) * xs.r_bot * h * h.sqrt(),
            Some(XSectShape::Powerfunc) => xs.r_bot * h.powf(xs.s_bot + 1.0), // LIBM
            _ => xsect::get_a_of_y(xs, h),
        }
    }

    pub fn width_of_depth(xs: &XSectParams, h: f64) -> f64 {
        let mut h = h;
        if h < 0.0 {
            return 0.0;
        }
        if h > xs.y_full {
            h = xs.y_full;
        }
        match shape(xs.type_) {
            Some(XSectShape::Circular) | Some(XSectShape::ForceMain) => circ_width(0.5 * xs.y_full, h),
            Some(XSectShape::FilledCircular) => circ_width(0.5 * (xs.y_full + xs.y_bot), xs.y_bot + h),
            Some(XSectShape::RectClosed) | Some(XSectShape::RectOpen) => xs.w_max,
            Some(XSectShape::Trapezoidal) => xs.y_bot + 2.0 * xs.s_bot * h,
            Some(XSectShape::Triangular) => 2.0 * xs.s_bot * h,
            Some(XSectShape::Parabolic) => 2.0 * xs.r_bot * h.sqrt(),
            Some(XSectShape::Powerfunc) => {
                if h > 0.0 {
                    (xs.s_bot + 1.0) * xs.r_bot * h.powf(xs.s_bot) // LIBM
                } else {
                    0.0
                }
            }
            _ => xsect::get_w_of_y(xs, h),
        }
    }

    pub fn perimeter_of_depth(xs: &XSectParams, h: f64) -> f64 {
        let mut h = h;
        if h <= 0.0 {
            return 0.0;
        }
        if h > xs.y_full {
            h = xs.y_full;
        }
        match shape(xs.type_) {
            Some(XSectShape::Circular) | Some(XSectShape::ForceMain) => circ_perim(0.5 * xs.y_full, h),
            Some(XSectShape::FilledCircular) => {
                let r = 0.5 * (xs.y_full + xs.y_bot);
                circ_perim(r, xs.y_bot + h) - circ_perim(r, xs.y_bot) + circ_width(r, xs.y_bot)
            }
            Some(XSectShape::RectClosed) => {
                let mut p = xs.w_max + 2.0 * h;
                let alf = if xs.a_full > 0.0 { (h * xs.w_max) / xs.a_full } else { 0.0 };
                if alf > K_RECT_ALF_MAX {
                    p += (alf - K_RECT_ALF_MAX) / (1.0 - K_RECT_ALF_MAX) * xs.w_max;
                }
                p
            }
            Some(XSectShape::RectOpen) => xs.w_max + (2.0 - xs.s_bot) * h,
            Some(XSectShape::Trapezoidal) => xs.y_bot + h * xs.r_bot,
            Some(XSectShape::Triangular) => 2.0 * h * xs.r_bot,
            Some(XSectShape::Parabolic) => {
                let x = 2.0 * h.sqrt() / xs.r_bot;
                let t = (1.0 + x * x).sqrt();
                0.5 * xs.r_bot * xs.r_bot * (x * t + (x + t).ln()) // LIBM
            }
            Some(XSectShape::Powerfunc) => {
                let hc = (xs.s_bot + 1.0) * xs.r_bot / 2.0;
                let m = xs.s_bot;
                let ds = |y: f64| -> f64 {
                    let dxdy = if y > 0.0 { hc * m * y.powf(m - 1.0) } else { 0.0 }; // LIBM
                    (1.0 + dxdy * dxdy).sqrt()
                };
                let y0 = 1.0e-9 * h;
                2.0 * (integrate(&ds, y0, h, 1.0e-12 * h, 40) + y0.hypot(hc * y0.powf(m))) // LIBM
            }
            _ => {
                let r = xsect::get_r_of_y(xs, h);
                if r > 0.0 {
                    xsect::get_a_of_y(xs, h) / r
                } else {
                    0.0
                }
            }
        }
    }

    pub fn hyd_rad_of_depth(xs: &XSectParams, h: f64) -> f64 {
        if h <= 0.0 {
            return 0.0;
        }
        if !has_closed_form(xs.type_) {
            return xsect::get_r_of_y(xs, crate::core::cmath::cmin(h, xs.y_full));
        }
        let p = perimeter_of_depth(xs, h);
        if p > 0.0 {
            area_of_depth(xs, h) / p
        } else {
            0.0
        }
    }
}
