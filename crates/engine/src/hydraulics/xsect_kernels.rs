// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/XSectKernels.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The cross-section geometry layer: the analytic shape formulas (`shape`), the table
// binding (`XsectTables`) and the evaluator (`XsectEval`) whose methods hold every
// per-shape formula of legacy xsect.c. Every expression keeps upstream's operand order
// and grouping; `a * b * c` is `(a * b) * c` in both languages.
//
// Differences in form, never in result on upstream's defined domain:
//   * Names are snake_case: getAofY -> get_a_of_y, getdSdA -> get_ds_da,
//     getYcircular -> get_ycircular, circ_getAofY -> circ_get_a_of_y,
//     findroot_Newton -> findroot_newton, shape::rectClosedRofA -> rect_closed_r_of_a.
//     XsectTables keeps upstream's field names (A_Circ, N_A_Circ, ...).
//   * XsectTables binds `&'static [f64]` slices instead of raw pointers (the only binding
//     in a Rust build is the host one, see xsection::host_tables), and `luts` is an
//     `Option<&'static [LocateLut]>`. Default gives empty slices and zero counts where
//     upstream gives null pointers.
//   * `switch (static_cast<XSectShape>(xs.type))` matches on XSectShape::from_code, so an
//     unknown code takes the `default` arm as it does upstream.
//   * `static_cast<int>` of a value that can be non-finite is core::cmath::cvt_i32
//     (x86's INT_MIN). Where upstream would then index a table with a negative index
//     (tabular_getdSdA with alpha <= -delta or NaN), or loop on INT_MIN (getYcrit's
//     enumeration when y0 is NaN), upstream is undefined; Rust panics on the index and
//     loops in the same way. `i1 + 1`/`i1 - 1` there are wrapping, as x86 executes them.
//   * findroot_Newton/findroot_Ridder take a closure; upstream's template takes a functor
//     by value. The `double* rts` in/out argument is `&mut f64`.
//   * Upstream leaves df/f uninitialized before the first callback, which assigns both;
//     they start at 0.0 here and are overwritten by the same callback.
//   * Recursion between get_s_of_a and get_r_of_a for DUMMY or an unknown code with a
//     positive (or NaN) area does not terminate upstream either (outside the documented
//     domain); it is translated as written.
//
// Transcendental calls (pow, sin, cos, acos, asin, log) go to the platform libm through
// f64::powf/sin/cos/acos/asin/ln and are marked `// LIBM`. sqrt, fabs and copysign are
// exact.
#![allow(non_snake_case, clippy::too_many_arguments, clippy::neg_cmp_op_on_partial_ord)]
#![allow(
    clippy::manual_range_contains,
    clippy::approx_constant,
    clippy::excessive_precision
)]
#![allow(clippy::manual_clamp, clippy::assign_op_pattern)]

use crate::core::cmath::{cmin, cvt_i32};
use crate::hydraulics::xsect_batch::{XSectParams, XSectShape};
use crate::hydraulics::xsect_lookup::{locate_bisect, locate_maybe_lut, lookup_exact, LocateLut};

// Constants matching legacy consts.h / xsect.c exactly.
pub const TINY: f64 = 1.0e-6;
pub const PI: f64 = 3.141592654;
pub const GRAVITY: f64 = 32.2;

pub const RECT_ALFMAX: f64 = 0.97;
pub const RECT_TRIANG_ALFMAX: f64 = 0.98;
pub const RECT_ROUND_ALFMAX: f64 = 0.98;

/// The analytic shape formulas, shared by the per-element evaluator and the batch
/// kernels (upstream `xsect::shape`).
pub mod shape {
    use super::RECT_ALFMAX;

    // ---- Area from depth ----

    /// Rectangular area: depth times width, no clamp to full depth.
    #[inline]
    pub fn rect_a_of_y(y: f64, w_max: f64) -> f64 {
        y * w_max
    }

    /// Trapezoid area: (bottom width + side slope * depth) * depth.
    #[inline]
    pub fn trapez_a_of_y(y: f64, y_bot: f64, s_bot: f64) -> f64 {
        (y_bot + s_bot * y) * y
    }

    /// Triangular area: y*y*s_bot.
    #[inline]
    pub fn triang_a_of_y(y: f64, s_bot: f64) -> f64 {
        y * y * s_bot
    }

    /// Parabolic area: (4/3)*r_bot*y*sqrt(y).
    #[inline]
    pub fn parab_a_of_y(y: f64, r_bot: f64) -> f64 {
        (4.0 / 3.0) * r_bot * y * y.sqrt()
    }

    /// Power-function area: r_bot*y^(s_bot+1).
    #[inline]
    pub fn powerfunc_a_of_y(y: f64, s_bot: f64, r_bot: f64) -> f64 {
        r_bot * y.powf(s_bot + 1.0) // LIBM
    }

    // ---- Hydraulic radius ----

    /// Closed rectangular radius with the near-crown perimeter correction; 0 when a <= 0.
    #[inline]
    pub fn rect_closed_r_of_a(a: f64, w_max: f64, a_full: f64) -> f64 {
        if a <= 0.0 {
            return 0.0;
        }
        let mut p = w_max + 2.0 * a / w_max;
        if a / a_full > RECT_ALFMAX {
            p += (a / a_full - RECT_ALFMAX) / (1.0 - RECT_ALFMAX) * w_max;
        }
        a / p
    }

    /// Open rectangular radius excluding s_bot banks; 0 when a <= 0.
    #[inline]
    pub fn rect_open_r_of_a(a: f64, w_max: f64, s_bot: f64) -> f64 {
        if a <= 0.0 {
            return 0.0;
        }
        a / (w_max + (2.0 - s_bot) * a / w_max)
    }

    /// Trapezoid radius from depth; exactly 0 at y == 0.
    #[inline]
    pub fn trapez_r_of_y(y: f64, y_bot: f64, s_bot: f64, r_bot: f64) -> f64 {
        if y == 0.0 {
            return 0.0;
        }
        trapez_a_of_y(y, y_bot, s_bot) / (y_bot + y * r_bot)
    }

    /// Triangular radius: y*s_bot/(2*r_bot).
    #[inline]
    pub fn triang_r_of_y(y: f64, s_bot: f64, r_bot: f64) -> f64 {
        (y * s_bot) / (2.0 * r_bot)
    }

    // ---- Top width ----

    /// Closed rectangle: 0 exactly at y_norm == 1, otherwise w_max.
    #[inline]
    pub fn rect_closed_w_of_y(y_norm: f64, w_max: f64) -> f64 {
        if y_norm == 1.0 {
            return 0.0;
        }
        w_max
    }

    /// Trapezoid top width: y_bot + 2*y*s_bot.
    #[inline]
    pub fn trapez_w_of_y(y: f64, y_bot: f64, s_bot: f64) -> f64 {
        y_bot + 2.0 * y * s_bot
    }

    /// Triangular top width: 2*s_bot*y.
    #[inline]
    pub fn triang_w_of_y(y: f64, s_bot: f64) -> f64 {
        2.0 * s_bot * y
    }

    /// Parabolic top width: 2*r_bot*sqrt(y).
    #[inline]
    pub fn parab_w_of_y(y: f64, r_bot: f64) -> f64 {
        2.0 * r_bot * y.sqrt()
    }

    /// Power-function top width: (s_bot+1)*r_bot*y^s_bot.
    #[inline]
    pub fn powerfunc_w_of_y(y: f64, s_bot: f64, r_bot: f64) -> f64 {
        (s_bot + 1.0) * r_bot * y.powf(s_bot) // LIBM
    }
}

/// The tables `inv_lookup()` inverts, and so the ones that carry a LocateLut.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LutId {
    YGothic = 0,
    YCatenary,
    YSemiEllip,
    YSemiCirc,
    AHorizEllipse,
    AVertEllipse,
    AArch,
    SCirc,
    SEgg,
    SHorseshoe,
    SGothic,
    SCatenary,
    SSemiEllip,
    SBasketHandle,
    SSemiCirc,
    Count,
}

/// Bindings of the shared geometry tables (upstream `XsectTables`).
#[derive(Debug, Clone, Copy, Default)]
pub struct XsectTables {
    pub A_Arch: &'static [f64],
    pub A_Baskethandle: &'static [f64],
    pub A_Circ: &'static [f64],
    pub A_Egg: &'static [f64],
    pub A_HorizEllipse: &'static [f64],
    pub A_Horseshoe: &'static [f64],
    pub A_VertEllipse: &'static [f64],
    pub R_Arch: &'static [f64],
    pub R_Baskethandle: &'static [f64],
    pub R_Circ: &'static [f64],
    pub R_Egg: &'static [f64],
    pub R_HorizEllipse: &'static [f64],
    pub R_Horseshoe: &'static [f64],
    pub R_VertEllipse: &'static [f64],
    pub S_BasketHandle: &'static [f64],
    pub S_Catenary: &'static [f64],
    pub S_Circ: &'static [f64],
    pub S_Egg: &'static [f64],
    pub S_Gothic: &'static [f64],
    pub S_Horseshoe: &'static [f64],
    pub S_SemiCirc: &'static [f64],
    pub S_SemiEllip: &'static [f64],
    pub W_Arch: &'static [f64],
    pub W_BasketHandle: &'static [f64],
    pub W_Catenary: &'static [f64],
    pub W_Circ: &'static [f64],
    pub W_Egg: &'static [f64],
    pub W_Gothic: &'static [f64],
    pub W_HorizEllipse: &'static [f64],
    pub W_Horseshoe: &'static [f64],
    pub W_SemiCirc: &'static [f64],
    pub W_SemiEllip: &'static [f64],
    pub W_VertEllipse: &'static [f64],
    pub Y_BasketHandle: &'static [f64],
    pub Y_Catenary: &'static [f64],
    pub Y_Circ: &'static [f64],
    pub Y_Egg: &'static [f64],
    pub Y_Gothic: &'static [f64],
    pub Y_Horseshoe: &'static [f64],
    pub Y_SemiCirc: &'static [f64],
    pub Y_SemiEllip: &'static [f64],
    pub Amax: &'static [f64],
    /// Bucket maps for the inverted tables, indexed by LutId; None == bisect.
    pub luts: Option<&'static [LocateLut]>,
    pub N_A_Arch: i32,
    pub N_A_Baskethandle: i32,
    pub N_A_Circ: i32,
    pub N_A_Egg: i32,
    pub N_A_HorizEllipse: i32,
    pub N_A_Horseshoe: i32,
    pub N_A_VertEllipse: i32,
    pub N_R_Arch: i32,
    pub N_R_Baskethandle: i32,
    pub N_R_Circ: i32,
    pub N_R_Egg: i32,
    pub N_R_HorizEllipse: i32,
    pub N_R_Horseshoe: i32,
    pub N_R_VertEllipse: i32,
    pub N_S_BasketHandle: i32,
    pub N_S_Catenary: i32,
    pub N_S_Circ: i32,
    pub N_S_Egg: i32,
    pub N_S_Gothic: i32,
    pub N_S_Horseshoe: i32,
    pub N_S_SemiCirc: i32,
    pub N_S_SemiEllip: i32,
    pub N_W_Arch: i32,
    pub N_W_BasketHandle: i32,
    pub N_W_Catenary: i32,
    pub N_W_Circ: i32,
    pub N_W_Egg: i32,
    pub N_W_Gothic: i32,
    pub N_W_HorizEllipse: i32,
    pub N_W_Horseshoe: i32,
    pub N_W_SemiCirc: i32,
    pub N_W_SemiEllip: i32,
    pub N_W_VertEllipse: i32,
    pub N_Y_BasketHandle: i32,
    pub N_Y_Catenary: i32,
    pub N_Y_Circ: i32,
    pub N_Y_Egg: i32,
    pub N_Y_Gothic: i32,
    pub N_Y_Horseshoe: i32,
    pub N_Y_SemiCirc: i32,
    pub N_Y_SemiEllip: i32,
}

impl XsectTables {
    /// The optional bucket map for one table; None when no maps are bound.
    #[inline]
    pub fn lut(&self, id: LutId) -> Option<&'static LocateLut> {
        self.luts.map(|l| &l[id as usize])
    }
}

/// The geometry layer, parameterized on its table binding (upstream `XsectEval`).
#[derive(Debug, Clone, Copy, Default)]
pub struct XsectEval {
    pub tbl: XsectTables,
}

type S = XSectShape;

impl XsectEval {
    // ============================================================================
    // Root finders — ports of legacy findroot.c (Numerical Recipes).
    // ============================================================================

    /// Newton-Raphson + bisection; func(x, &f, &df). Root bracketed in [x1, x2].
    /// Returns the evaluation count when <= 60, otherwise 0; `rts` holds the final x.
    pub fn findroot_newton<F: FnMut(f64, &mut f64, &mut f64)>(
        &self,
        x1: f64,
        x2: f64,
        rts: &mut f64,
        xacc: f64,
        mut func: F,
    ) -> i32 {
        const MAXIT: i32 = 60;
        let mut n: i32 = 0;
        let mut df: f64 = 0.0;
        let mut f: f64 = 0.0;
        let mut dx: f64;
        let mut dxold: f64;
        let mut x: f64;
        let mut temp: f64;
        let mut xhi: f64;
        let mut xlo: f64;

        x = *rts;
        xlo = x1;
        xhi = x2;
        dxold = (x2 - x1).abs();
        dx = dxold;
        func(x, &mut f, &mut df);
        n += 1;
        for _j in 1..=MAXIT {
            if ((x - xhi) * df - f) * ((x - xlo) * df - f) >= 0.0 || ((2.0 * f).abs() > (dxold * df).abs()) {
                dxold = dx;
                dx = 0.5 * (xhi - xlo);
                x = xlo + dx;
                if xlo == x {
                    break;
                }
            } else {
                dxold = dx;
                dx = f / df;
                temp = x;
                x -= dx;
                if temp == x {
                    break;
                }
            }
            if dx.abs() < xacc {
                break;
            }
            func(x, &mut f, &mut df);
            n += 1;
            if f < 0.0 {
                xlo = x;
            } else {
                xhi = x;
            }
        }
        *rts = x;
        if n <= MAXIT {
            n
        } else {
            0
        }
    }

    /// Ridder's method; func(x) -> residual. Root bracketed in [x1, x2]; -1e20 when the
    /// endpoint residuals do not bracket a root.
    pub fn findroot_ridder<F: FnMut(f64) -> f64>(&self, x1: f64, x2: f64, xacc: f64, mut func: F) -> f64 {
        const MAXIT: i32 = 60;
        let sign = |a: f64, b: f64| if b >= 0.0 { a.abs() } else { -a.abs() };

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
            for _j in 1..=MAXIT {
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
            }
            return ans;
        }
        -1.0e20
    }

    // ============================================================================
    // Lookup helpers — legacy lookup()/invLookup()/locate()
    // ============================================================================

    /// Lower interpolation bracket by legacy bisection.
    #[inline]
    pub fn locate(&self, y: f64, table: &[f64], j_last: i32) -> i32 {
        locate_bisect(y, table, j_last)
    }

    /// Interpolate a normalized table: NaN/nonpositive x gives table[0], x > 1 the final
    /// entry, otherwise lookup_exact (x == 1 included, as legacy).
    #[inline]
    pub fn lookup(&self, x: f64, table: &[f64], n_items: i32) -> f64 {
        if !(x > 0.0) {
            return table[0]; // x <= 0 or NaN
        }
        if x > 1.0 {
            return table[(n_items - 1) as usize]; // x > 1 or +inf
        }
        lookup_exact(x, table, n_items)
    }

    /// Invert a normalized table, including the legacy declining-top branch.
    pub fn inv_lookup(&self, y: f64, table: &[f64], n_items: i32, lut: Option<&LocateLut>) -> f64 {
        let dx = 1.0 / ((n_items - 1) as f64);
        let mut n = n_items;

        // Truncate item count if last 2 entries are decreasing (section-factor tables)
        if table[(n - 3) as usize] > table[(n - 1) as usize] {
            n -= 2;
        }

        let i: i32;
        if n < n_items && y > table[(n_items - 1) as usize] {
            if y >= table[(n_items - 3) as usize] {
                return ((n - 1) as f64) * dx;
            }
            if y <= table[(n_items - 2) as usize] {
                i = n_items - 2;
            } else {
                i = n_items - 3;
            }
        } else {
            i = locate_maybe_lut(y, table, n - 1, lut);
        }
        if i >= n - 1 {
            return ((n - 1) as f64) * dx;
        }

        let iu = i as usize;
        let x0 = (i as f64) * dx;
        let dy = table[iu + 1] - table[iu];
        let mut x = if dy == 0.0 {
            x0
        } else {
            x0 + (y - table[iu]) * dx / dy
        };
        if x < 0.0 {
            x = 0.0;
        }
        if x > 1.0 {
            x = 1.0;
        }
        x
    }

    // ============================================================================
    // Circular small-area special functions
    // ============================================================================

    /// Circular central angle from normalized area (Newton, up to 40 steps).
    pub fn get_theta_of_alpha(&self, alpha: f64) -> f64 {
        let mut theta: f64;
        if alpha > 0.04 {
            theta = 1.2 + 5.08 * (alpha - 0.04) / 0.96;
        } else {
            theta = 0.031715 - 12.79384 * alpha + 8.28479 * alpha.sqrt();
        }
        let theta1 = theta;
        let ap = (2.0 * PI) * alpha;
        for _k in 1..=40 {
            let mut d = -(ap - theta + theta.sin()) / (1.0 - theta.cos()); // LIBM
            if d > 1.0 {
                d = 1.0f64.copysign(d);
            }
            theta -= d;
            if d.abs() <= 0.0001 {
                return theta;
            }
        }
        theta1
    }

    /// Circular central angle from normalized section factor (Newton, up to 40 steps).
    pub fn get_theta_of_psi(&self, psi: f64) -> f64 {
        let mut theta: f64;
        if psi > 0.90 {
            theta = 4.17 + 1.12 * (psi - 0.90) / 0.176;
        } else if psi > 0.5 {
            theta = 3.14 + 1.03 * (psi - 0.5) / 0.4;
        } else if psi > 0.015 {
            theta = 1.2 + 1.94 * (psi - 0.015) / 0.485;
        } else {
            theta = 0.12103 - 55.5075 * psi + 15.62254 * psi.sqrt();
        }
        let theta1 = theta;
        let ap = (2.0 * PI) * psi;
        for _k in 1..=40 {
            theta = theta.abs();
            let tt = theta - theta.sin(); // LIBM
            let tt23 = tt.powf(2.0 / 3.0); // LIBM
            let t3 = theta.powf(1.0 / 3.0); // LIBM
            let mut d = ap * theta / t3 - tt * tt23;
            d = d / (ap * (2.0 / 3.0) / t3 - (5.0 / 3.0) * tt23 * (1.0 - theta.cos())); // LIBM
            theta -= d;
            if d.abs() <= 0.0001 {
                return theta;
            }
        }
        theta1
    }

    /// Circular normalized depth from area ratio.
    pub fn get_ycircular(&self, alpha: f64) -> f64 {
        if alpha >= 1.0 {
            return 1.0;
        }
        if alpha <= 0.0 {
            return 0.0;
        }
        if alpha <= 1.0e-5 {
            let theta = (37.6911 * alpha).powf(1.0 / 3.0); // LIBM
            return theta * theta / 16.0;
        }
        let theta = self.get_theta_of_alpha(alpha);
        (1.0 - (theta / 2.0).cos()) / 2.0 // LIBM
    }

    /// Circular normalized section factor from area ratio.
    pub fn get_scircular(&self, alpha: f64) -> f64 {
        if alpha >= 1.0 {
            return 1.0;
        }
        if alpha <= 0.0 {
            return 0.0;
        }
        if alpha <= 1.0e-5 {
            let theta = (37.6911 * alpha).powf(1.0 / 3.0); // LIBM
            return theta.powf(13.0 / 3.0) / 124.4797; // LIBM
        }
        let theta = self.get_theta_of_alpha(alpha);
        (theta - theta.sin()).powf(5.0 / 3.0) / (2.0 * PI) / theta.powf(2.0 / 3.0)
        // LIBM
    }

    /// Circular normalized area from section-factor ratio.
    pub fn get_acircular(&self, psi: f64) -> f64 {
        if psi >= 1.0 {
            return 1.0;
        }
        if psi <= 0.0 {
            return 0.0;
        }
        if psi <= 1.0e-6 {
            let theta = (124.4797 * psi).powf(3.0 / 13.0); // LIBM
            return theta * theta * theta / 37.6911;
        }
        let theta = self.get_theta_of_psi(psi);
        (theta - theta.sin()) / (2.0 * PI) // LIBM
    }

    // ============================================================================
    // Generic / tabular section-factor derivatives
    // ============================================================================

    /// Centered normalized-area difference of the section factor (alpha +- 0.001).
    pub fn generic_get_ds_da(&self, xs: &XSectParams, a: f64) -> f64 {
        let alpha = a / xs.a_full;
        let mut alpha1 = alpha - 0.001;
        let alpha2 = alpha + 0.001;
        if alpha1 < 0.0 {
            alpha1 = 0.0;
        }
        let a1 = alpha1 * xs.a_full;
        let a2 = alpha2 * xs.a_full;
        (self.get_s_of_a(xs, a2) - self.get_s_of_a(xs, a1)) / (a2 - a1)
    }

    /// Slope of the local normalized S-table segment, scaled by s_full/a_full.
    pub fn tabular_get_ds_da(&self, xs: &XSectParams, a: f64, table: &[f64], n_items: i32) -> f64 {
        let alpha = a / xs.a_full;
        let delta = 1.0 / ((n_items - 1) as f64);
        let mut i = cvt_i32(alpha / delta);
        if i >= n_items - 1 {
            i = n_items - 2;
        }
        let iu = i as usize; // negative i: upstream reads before the table (UB); panics here
        let d_s_d_a = (table[iu + 1] - table[iu]) / delta;
        d_s_d_a * xs.s_full / xs.a_full
    }

    // ============================================================================
    // Circular
    // ============================================================================

    pub fn circ_get_a_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        let y_norm = y / xs.y_full;
        xs.a_full * self.lookup(y_norm, self.tbl.A_Circ, self.tbl.N_A_Circ)
    }

    pub fn circ_get_y_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        let alpha = a / xs.a_full;
        if alpha < 0.04 {
            return xs.y_full * self.get_ycircular(alpha);
        }
        xs.y_full * self.lookup(alpha, self.tbl.Y_Circ, self.tbl.N_Y_Circ)
    }

    pub fn circ_get_s_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        let alpha = a / xs.a_full;
        if alpha < 0.04 {
            return xs.s_full * self.get_scircular(alpha);
        }
        xs.s_full * self.lookup(alpha, self.tbl.S_Circ, self.tbl.N_S_Circ)
    }

    pub fn circ_get_a_of_s(&self, xs: &XSectParams, s: f64) -> f64 {
        let psi = s / xs.s_full;
        if psi == 0.0 {
            return 0.0;
        }
        if psi >= 1.0 {
            return xs.a_full;
        }
        if psi <= 0.015 {
            return xs.a_full * self.get_acircular(psi);
        }
        xs.a_full
            * self.inv_lookup(
                psi,
                self.tbl.S_Circ,
                self.tbl.N_S_Circ,
                self.tbl.lut(LutId::SCirc),
            )
    }

    pub fn circ_get_ds_da(&self, xs: &XSectParams, a: f64) -> f64 {
        let alpha = a / xs.a_full;
        if alpha <= 1.0e-30 {
            return 1.0e-30;
        }
        if alpha < 0.04 {
            let theta = self.get_theta_of_alpha(alpha);
            let p = theta * xs.y_full / 2.0;
            let r = a / p;
            let d_p_d_a = 4.0 / xs.y_full / (1.0 - theta.cos()); // LIBM
            return (5.0 / 3.0 - (2.0 / 3.0) * d_p_d_a * r) * r.powf(2.0 / 3.0); // LIBM
        }
        self.tabular_get_ds_da(xs, a, self.tbl.S_Circ, self.tbl.N_S_Circ)
    }

    // ============================================================================
    // Filled circular (yBot/aBot/sBot/rBot = filled bottom: depth/area/width/perim)
    // ============================================================================

    pub fn filled_circ_get_a_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        // Temporarily restore the unfilled full circle (legacy expands yFull/aFull).
        let mut t = *xs;
        t.y_full += xs.y_bot;
        t.a_full += xs.a_bot;
        let a = self.circ_get_a_of_y(&t, y + xs.y_bot);
        a - xs.a_bot
    }

    pub fn filled_circ_get_y_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        let mut t = *xs;
        t.y_full += xs.y_bot;
        t.a_full += xs.a_bot;
        let y = self.circ_get_y_of_a(&t, a + xs.a_bot);
        y - xs.y_bot
    }

    pub fn filled_circ_get_r_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        let mut t = *xs;
        t.y_full += xs.y_bot;
        t.a_full += xs.a_bot;
        let yy = y + xs.y_bot;
        let mut a = self.circ_get_a_of_y(&t, yy);
        let r = 0.25 * t.y_full * self.lookup(yy / t.y_full, self.tbl.R_Circ, self.tbl.N_R_Circ);
        let mut p = a / r;
        a -= xs.a_bot;
        p = p - xs.r_bot + xs.s_bot; // rBot = filled perimeter, sBot = filled width
        a / p
    }

    // ============================================================================
    // Rectangular closed / open
    // ============================================================================

    pub fn rect_closed_get_r_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        shape::rect_closed_r_of_a(a, xs.w_max, xs.a_full)
    }

    pub fn rect_closed_get_s_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        if a / xs.a_full > RECT_ALFMAX {
            return xs.s_max + (xs.s_full - xs.s_max) * (a / xs.a_full - RECT_ALFMAX) / (1.0 - RECT_ALFMAX);
        }
        a * self.rect_closed_get_r_of_a(xs, a).powf(2.0 / 3.0) // LIBM
    }

    pub fn rect_closed_get_ds_da(&self, xs: &XSectParams, a: f64) -> f64 {
        let alpha = a / xs.a_full;
        if alpha > RECT_ALFMAX {
            return (xs.s_full - xs.s_max) / ((1.0 - RECT_ALFMAX) * xs.a_full);
        }
        if alpha <= 1.0e-30 {
            return self.generic_get_ds_da(xs, a);
        }
        let r = self.rect_closed_get_r_of_a(xs, a);
        (5.0 / 3.0 - (2.0 / 3.0) * (2.0 / xs.w_max) * r) * r.powf(2.0 / 3.0) // LIBM
    }

    pub fn rect_open_get_s_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        let y = a / xs.w_max;
        let r = a / ((2.0 - xs.s_bot) * y + xs.w_max);
        a * r.powf(2.0 / 3.0) // LIBM
    }

    pub fn rect_open_get_ds_da(&self, xs: &XSectParams, a: f64) -> f64 {
        if a / xs.a_full <= 1.0e-30 {
            return self.generic_get_ds_da(xs, a);
        }
        let r = self.get_r_of_a(xs, a);
        let d_p_d_a = (2.0 - xs.s_bot) / xs.w_max;
        (5.0 / 3.0 - (2.0 / 3.0) * d_p_d_a * r) * r.powf(2.0 / 3.0) // LIBM
    }

    // ============================================================================
    // Rect-triangular (triangular bottom + rectangular top)
    // ============================================================================

    pub fn rect_triang_get_a_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        if y <= xs.y_bot {
            return y * y * xs.s_bot;
        }
        xs.a_bot + (y - xs.y_bot) * xs.w_max
    }

    pub fn rect_triang_get_w_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        if y <= xs.y_bot {
            return 2.0 * xs.s_bot * y;
        }
        xs.w_max
    }

    pub fn rect_triang_get_y_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        if a <= xs.a_bot {
            return (a / xs.s_bot).sqrt();
        }
        xs.y_bot + (a - xs.a_bot) / xs.w_max
    }

    pub fn rect_triang_get_r_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        if a <= 0.0 {
            return 0.0;
        }
        let y = self.rect_triang_get_y_of_a(xs, a);
        if y <= xs.y_bot {
            return a / (2.0 * y * xs.r_bot);
        }
        let mut p = 2.0 * xs.y_bot * xs.r_bot + 2.0 * (y - xs.y_bot);
        let alf = (a / xs.a_full) - RECT_TRIANG_ALFMAX;
        if alf > 0.0 {
            p += alf / (1.0 - RECT_TRIANG_ALFMAX) * xs.w_max;
        }
        a / p
    }

    pub fn rect_triang_get_r_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        if y <= xs.y_bot {
            return y * xs.s_bot / (2.0 * xs.r_bot);
        }
        let a = xs.a_bot + (y - xs.y_bot) * xs.w_max;
        let mut p = 2.0 * xs.y_bot * xs.r_bot + 2.0 * (y - xs.y_bot);
        let alf = (a / xs.a_full) - RECT_TRIANG_ALFMAX;
        if alf > 0.0 {
            p += alf / (1.0 - RECT_TRIANG_ALFMAX) * xs.w_max;
        }
        a / p
    }

    pub fn rect_triang_get_s_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        let alf_max = RECT_TRIANG_ALFMAX;
        if a / xs.a_full > alf_max {
            return xs.s_max + (xs.s_full - xs.s_max) * (a / xs.a_full - alf_max) / (1.0 - alf_max);
        }
        a * self.rect_triang_get_r_of_a(xs, a).powf(2.0 / 3.0) // LIBM
    }

    pub fn rect_triang_get_ds_da(&self, xs: &XSectParams, a: f64) -> f64 {
        let alf_max = RECT_TRIANG_ALFMAX;
        let alpha = a / xs.a_full;
        if alpha > alf_max {
            return (xs.s_full - xs.s_max) / ((1.0 - alf_max) * xs.a_full);
        }
        if alpha <= 1.0e-30 {
            return self.generic_get_ds_da(xs, a);
        }
        let d_p_d_a = if a > xs.a_bot {
            2.0 / xs.w_max
        } else {
            xs.r_bot / (a * xs.s_bot).sqrt()
        };
        let r = self.rect_triang_get_r_of_a(xs, a);
        (5.0 / 3.0 - (2.0 / 3.0) * d_p_d_a * r) * r.powf(2.0 / 3.0) // LIBM
    }

    // ============================================================================
    // Rect-round (circular invert + rectangular top)
    // ============================================================================

    pub fn rect_round_get_y_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        if a > xs.a_bot {
            return xs.y_bot + (a - xs.a_bot) / xs.w_max;
        }
        let alpha = a / (PI * xs.r_bot * xs.r_bot);
        if alpha < 0.04 {
            return (2.0 * xs.r_bot) * self.get_ycircular(alpha);
        }
        (2.0 * xs.r_bot) * self.lookup(alpha, self.tbl.Y_Circ, self.tbl.N_Y_Circ)
    }

    pub fn rect_round_get_a_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        if y > xs.y_bot {
            return xs.a_bot + (y - xs.y_bot) * xs.w_max;
        }
        let theta1 = 2.0 * (1.0 - y / xs.r_bot).acos(); // LIBM
        0.5 * xs.r_bot * xs.r_bot * (theta1 - theta1.sin()) // LIBM
    }

    pub fn rect_round_get_w_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        if y > xs.y_bot {
            return xs.w_max;
        }
        2.0 * (y * (2.0 * xs.r_bot - y)).sqrt()
    }

    pub fn rect_round_get_r_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        if a <= 0.0 {
            return 0.0;
        }
        if a > xs.a_bot {
            let y1 = (a - xs.a_bot) / xs.w_max;
            let theta1 = 2.0 * (xs.w_max / 2.0 / xs.r_bot).asin(); // LIBM
            let mut p = xs.r_bot * theta1 + 2.0 * y1;
            let arg = (a / xs.a_full) - RECT_ROUND_ALFMAX;
            if arg > 0.0 {
                p += arg / (1.0 - RECT_ROUND_ALFMAX) * xs.w_max;
            }
            return a / p;
        }
        let y1 = self.rect_round_get_y_of_a(xs, a);
        let theta1 = 2.0 * (1.0 - y1 / xs.r_bot).acos(); // LIBM
        let p = xs.r_bot * theta1;
        a / p
    }

    pub fn rect_round_get_r_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        if y <= 0.0 {
            return 0.0;
        }
        if y > xs.y_bot {
            return self.rect_round_get_r_of_a(xs, self.rect_round_get_a_of_y(xs, y));
        }
        let theta1 = 2.0 * (1.0 - y / xs.r_bot).acos(); // LIBM
        0.5 * xs.r_bot * (1.0 - theta1.sin()) / theta1 // LIBM
    }

    pub fn rect_round_get_s_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        let alf_max = RECT_ROUND_ALFMAX;
        if a / xs.a_full > alf_max {
            return xs.s_max + (xs.s_full - xs.s_max) * (a / xs.a_full - alf_max) / (1.0 - alf_max);
        }
        if a > xs.a_bot {
            return a * self.rect_round_get_r_of_a(xs, a).powf(2.0 / 3.0); // LIBM
        }
        let a_full = PI * xs.r_bot * xs.r_bot;
        let alpha = a / a_full;
        let s_full = xs.s_bot;
        if alpha < 0.04 {
            return s_full * self.get_scircular(alpha);
        }
        s_full * self.lookup(alpha, self.tbl.S_Circ, self.tbl.N_S_Circ)
    }

    pub fn rect_round_get_ds_da(&self, xs: &XSectParams, a: f64) -> f64 {
        let alf_max = RECT_ROUND_ALFMAX;
        if a / xs.a_full > alf_max {
            return (xs.s_full - xs.s_max) / ((1.0 - alf_max) * xs.a_full);
        }
        if a > xs.a_bot {
            let r = self.rect_round_get_r_of_a(xs, a);
            let d_p_d_a = 2.0 / xs.w_max;
            return (5.0 / 3.0 - (2.0 / 3.0) * d_p_d_a * r) * r.powf(2.0 / 3.0); // LIBM
        }
        self.generic_get_ds_da(xs, a)
    }

    // ============================================================================
    // Modified baskethandle (rectangular bottom + circular-arc top).
    // rBot/yBot/aBot/sBot refer to the CIRCULAR TOP portion.
    // ============================================================================

    pub fn mod_basket_get_y_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        if a <= xs.a_full - xs.a_bot {
            return a / xs.w_max;
        }
        let alpha = (xs.a_full - a) / (PI * xs.r_bot * xs.r_bot);
        let mut y1 = if alpha < 0.04 {
            self.get_ycircular(alpha)
        } else {
            self.lookup(alpha, self.tbl.Y_Circ, self.tbl.N_Y_Circ)
        };
        y1 = 2.0 * xs.r_bot * y1;
        xs.y_full - y1
    }

    pub fn mod_basket_get_a_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        if y <= xs.y_full - xs.y_bot {
            return y * xs.w_max;
        }
        let y1 = xs.y_full - y;
        let theta1 = 2.0 * (1.0 - y1 / xs.r_bot).acos(); // LIBM
        let a1 = 0.5 * xs.r_bot * xs.r_bot * (theta1 - theta1.sin()); // LIBM
        xs.a_full - a1
    }

    pub fn mod_basket_get_w_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        if y <= 0.0 {
            return 0.0;
        }
        if y <= xs.y_full - xs.y_bot {
            return xs.w_max;
        }
        let y1 = xs.y_full - y;
        2.0 * (y1 * (2.0 * xs.r_bot - y1)).sqrt()
    }

    pub fn mod_basket_get_r_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        if a <= xs.a_full - xs.a_bot {
            return a / (xs.w_max + 2.0 * a / xs.w_max);
        }
        let mut y1 = xs.y_full - self.mod_basket_get_y_of_a(xs, a);
        let theta1 = 2.0 * (1.0 - y1 / xs.r_bot).acos(); // LIBM
        let mut p = (xs.s_bot - theta1) * xs.r_bot; // sBot = full circular opening angle
        y1 = xs.y_full - xs.y_bot;
        p = p + 2.0 * y1 + xs.w_max;
        a / p
    }

    pub fn mod_basket_get_ds_da(&self, xs: &XSectParams, a: f64) -> f64 {
        if a <= xs.a_full - xs.a_bot && a / xs.a_full > 1.0e-30 {
            let r = a / (xs.w_max + 2.0 * a / xs.w_max);
            let d_p_d_a = 2.0 / xs.w_max;
            return (5.0 / 3.0 - (2.0 / 3.0) * d_p_d_a * r) * r.powf(2.0 / 3.0); // LIBM
        }
        self.generic_get_ds_da(xs, a)
    }

    // ============================================================================
    // Trapezoidal (yBot = bottom width, sBot = avg side slope, rBot = side length/depth)
    // ============================================================================

    pub fn trapez_get_a_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        shape::trapez_a_of_y(y, xs.y_bot, xs.s_bot)
    }

    pub fn trapez_get_w_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        shape::trapez_w_of_y(y, xs.y_bot, xs.s_bot)
    }

    pub fn trapez_get_y_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        if xs.s_bot == 0.0 {
            return a / xs.y_bot;
        }
        ((xs.y_bot * xs.y_bot + 4.0 * xs.s_bot * a).sqrt() - xs.y_bot) / (2.0 * xs.s_bot)
    }

    pub fn trapez_get_r_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        a / (xs.y_bot + self.trapez_get_y_of_a(xs, a) * xs.r_bot)
    }

    pub fn trapez_get_r_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        shape::trapez_r_of_y(y, xs.y_bot, xs.s_bot, xs.r_bot)
    }

    pub fn trapez_get_ds_da(&self, xs: &XSectParams, a: f64) -> f64 {
        if a / xs.a_full <= 1.0e-30 {
            return self.generic_get_ds_da(xs, a);
        }
        let r = self.trapez_get_r_of_a(xs, a);
        let d_p_d_a = xs.r_bot / (xs.y_bot * xs.y_bot + 4.0 * xs.s_bot * a).sqrt();
        (5.0 / 3.0 - (2.0 / 3.0) * d_p_d_a * r) * r.powf(2.0 / 3.0) // LIBM
    }

    // ============================================================================
    // Triangular
    // ============================================================================

    pub fn triang_get_a_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        shape::triang_a_of_y(y, xs.s_bot)
    }

    pub fn triang_get_w_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        shape::triang_w_of_y(y, xs.s_bot)
    }

    pub fn triang_get_y_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        (a / xs.s_bot).sqrt()
    }

    pub fn triang_get_r_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        a / (2.0 * self.triang_get_y_of_a(xs, a) * xs.r_bot)
    }

    pub fn triang_get_r_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        shape::triang_r_of_y(y, xs.s_bot, xs.r_bot)
    }

    pub fn triang_get_ds_da(&self, xs: &XSectParams, a: f64) -> f64 {
        if a / xs.a_full <= 1.0e-30 {
            return self.generic_get_ds_da(xs, a);
        }
        let r = self.triang_get_r_of_a(xs, a);
        let d_p_d_a = xs.r_bot / (a * xs.s_bot).sqrt();
        (5.0 / 3.0 - (2.0 / 3.0) * d_p_d_a * r) * r.powf(2.0 / 3.0) // LIBM
    }

    // ============================================================================
    // Parabolic (rBot = 1/sqrt(c) where y = c*x^2)
    // ============================================================================

    pub fn parab_get_a_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        shape::parab_a_of_y(y, xs.r_bot)
    }

    pub fn parab_get_w_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        shape::parab_w_of_y(y, xs.r_bot)
    }

    pub fn parab_get_y_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        ((3.0 / 4.0) * a / xs.r_bot).powf(2.0 / 3.0) // LIBM
    }

    /// Analytical wetted perimeter (legacy parab_getPofY).
    pub fn parab_get_p_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        let x = 2.0 * y.sqrt() / xs.r_bot;
        let t = (1.0 + x * x).sqrt();
        0.5 * xs.r_bot * xs.r_bot * (x * t + (x + t).ln()) // LIBM
    }

    pub fn parab_get_r_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        if y <= 0.0 {
            return 0.0;
        }
        self.parab_get_a_of_y(xs, y) / self.parab_get_p_of_y(xs, y)
    }

    pub fn parab_get_r_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        if a <= 0.0 {
            return 0.0;
        }
        a / self.parab_get_p_of_y(xs, self.parab_get_y_of_a(xs, a))
    }

    // ============================================================================
    // Power function (sBot = 1/exponent, rBot = coefficient)
    // ============================================================================

    pub fn powerfunc_get_a_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        shape::powerfunc_a_of_y(y, xs.s_bot, xs.r_bot)
    }

    pub fn powerfunc_get_w_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        shape::powerfunc_w_of_y(y, xs.s_bot, xs.r_bot)
    }

    pub fn powerfunc_get_y_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        (a / xs.r_bot).powf(1.0 / (xs.s_bot + 1.0)) // LIBM
    }

    /// Numerical wetted perimeter (legacy powerfunc_getPofY, stepwise integration in
    /// depth increments of 0.02*y_full; does not terminate unless that step is positive).
    pub fn powerfunc_get_p_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        let dy1 = 0.02 * xs.y_full;
        let h = (xs.s_bot + 1.0) * xs.r_bot / 2.0;
        let m = xs.s_bot;
        let mut p = 0.0;
        let mut y1 = 0.0;
        let mut x1 = 0.0;
        let mut y2: f64;
        loop {
            y2 = y1 + dy1;
            if y2 > y {
                y2 = y;
            }
            let x2 = h * y2.powf(m); // LIBM
            let dx = x2 - x1;
            let dy = y2 - y1;
            p += (dx * dx + dy * dy).sqrt();
            x1 = x2;
            y1 = y2;
            if !(y2 < y) {
                break;
            }
        }
        2.0 * p
    }

    pub fn powerfunc_get_r_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        if y <= 0.0 {
            return 0.0;
        }
        self.powerfunc_get_a_of_y(xs, y) / self.powerfunc_get_p_of_y(xs, y)
    }

    pub fn powerfunc_get_r_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        if a <= 0.0 {
            return 0.0;
        }
        a / self.powerfunc_get_p_of_y(xs, self.powerfunc_get_y_of_a(xs, a))
    }

    // ============================================================================
    // Main dispatch: getAofY
    // ============================================================================

    /// Wetted area from depth; 0 for y <= 0, nonpositive full depth, a missing optional
    /// area table or an unsupported code.
    pub fn get_a_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        if y <= 0.0 {
            return 0.0;
        }
        // A section with no (or not-yet-set) full depth has no area.
        if xs.y_full <= 0.0 {
            return 0.0;
        }
        let y_norm = y / xs.y_full;
        let t = &self.tbl;

        match XSectShape::from_code(xs.type_) {
            Some(S::ForceMain) | Some(S::Circular) => xs.a_full * self.lookup(y_norm, t.A_Circ, t.N_A_Circ),
            Some(S::FilledCircular) => self.filled_circ_get_a_of_y(xs, y),
            Some(S::Eggshaped) => xs.a_full * self.lookup(y_norm, t.A_Egg, t.N_A_Egg),
            Some(S::Horseshoe) => xs.a_full * self.lookup(y_norm, t.A_Horseshoe, t.N_A_Horseshoe),
            Some(S::Gothic) => {
                xs.a_full * self.inv_lookup(y_norm, t.Y_Gothic, t.N_Y_Gothic, t.lut(LutId::YGothic))
            }
            Some(S::Catenary) => {
                xs.a_full * self.inv_lookup(y_norm, t.Y_Catenary, t.N_Y_Catenary, t.lut(LutId::YCatenary))
            }
            Some(S::Semielliptical) => {
                xs.a_full * self.inv_lookup(y_norm, t.Y_SemiEllip, t.N_Y_SemiEllip, t.lut(LutId::YSemiEllip))
            }
            Some(S::Baskethandle) => xs.a_full * self.lookup(y_norm, t.A_Baskethandle, t.N_A_Baskethandle),
            Some(S::Semicircular) => {
                xs.a_full * self.inv_lookup(y_norm, t.Y_SemiCirc, t.N_Y_SemiCirc, t.lut(LutId::YSemiCirc))
            }
            Some(S::HorizEllipse) => xs.a_full * self.lookup(y_norm, t.A_HorizEllipse, t.N_A_HorizEllipse),
            Some(S::VertEllipse) => xs.a_full * self.lookup(y_norm, t.A_VertEllipse, t.N_A_VertEllipse),
            Some(S::Arch) => xs.a_full * self.lookup(y_norm, t.A_Arch, t.N_A_Arch),
            Some(S::RectClosed) | Some(S::RectOpen) => shape::rect_a_of_y(y, xs.w_max),
            Some(S::RectTriang) => self.rect_triang_get_a_of_y(xs, y),
            Some(S::RectRound) => self.rect_round_get_a_of_y(xs, y),
            Some(S::ModBasket) => self.mod_basket_get_a_of_y(xs, y),
            Some(S::Trapezoidal) => self.trapez_get_a_of_y(xs, y),
            Some(S::Triangular) => self.triang_get_a_of_y(xs, y),
            Some(S::Parabolic) => self.parab_get_a_of_y(xs, y),
            Some(S::Powerfunc) => self.powerfunc_get_a_of_y(xs, y),
            Some(S::Irregular) | Some(S::Custom) | Some(S::StreetXsect) => {
                // Tabulated transect (legacy xsect_getAofY IRREGULAR case).
                if let Some(tb) = xs.area_tbl {
                    return xs.a_full * self.lookup(y_norm, tb, xs.transect_tbl_size);
                }
                0.0
            }
            _ => 0.0,
        }
    }

    // ============================================================================
    // Main dispatch: getWofY
    // ============================================================================

    /// Top width from depth. No universal nonpositive-depth guard.
    pub fn get_w_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        let y_norm = y / xs.y_full;
        let t = &self.tbl;

        match XSectShape::from_code(xs.type_) {
            Some(S::ForceMain) | Some(S::Circular) => xs.w_max * self.lookup(y_norm, t.W_Circ, t.N_W_Circ),
            Some(S::FilledCircular) => {
                let yn = (y + xs.y_bot) / (xs.y_full + xs.y_bot);
                xs.w_max * self.lookup(yn, t.W_Circ, t.N_W_Circ)
            }
            Some(S::Eggshaped) => xs.w_max * self.lookup(y_norm, t.W_Egg, t.N_W_Egg),
            Some(S::Horseshoe) => xs.w_max * self.lookup(y_norm, t.W_Horseshoe, t.N_W_Horseshoe),
            Some(S::Gothic) => xs.w_max * self.lookup(y_norm, t.W_Gothic, t.N_W_Gothic),
            Some(S::Catenary) => xs.w_max * self.lookup(y_norm, t.W_Catenary, t.N_W_Catenary),
            Some(S::Semielliptical) => xs.w_max * self.lookup(y_norm, t.W_SemiEllip, t.N_W_SemiEllip),
            Some(S::Baskethandle) => xs.w_max * self.lookup(y_norm, t.W_BasketHandle, t.N_W_BasketHandle),
            Some(S::Semicircular) => xs.w_max * self.lookup(y_norm, t.W_SemiCirc, t.N_W_SemiCirc),
            Some(S::HorizEllipse) => xs.w_max * self.lookup(y_norm, t.W_HorizEllipse, t.N_W_HorizEllipse),
            Some(S::VertEllipse) => xs.w_max * self.lookup(y_norm, t.W_VertEllipse, t.N_W_VertEllipse),
            Some(S::Arch) => xs.w_max * self.lookup(y_norm, t.W_Arch, t.N_W_Arch),
            Some(S::RectClosed) => shape::rect_closed_w_of_y(y_norm, xs.w_max),
            Some(S::RectOpen) => xs.w_max,
            Some(S::RectTriang) => self.rect_triang_get_w_of_y(xs, y),
            Some(S::RectRound) => self.rect_round_get_w_of_y(xs, y),
            Some(S::ModBasket) => self.mod_basket_get_w_of_y(xs, y),
            Some(S::Trapezoidal) => self.trapez_get_w_of_y(xs, y),
            Some(S::Triangular) => self.triang_get_w_of_y(xs, y),
            Some(S::Parabolic) => self.parab_get_w_of_y(xs, y),
            Some(S::Powerfunc) => self.powerfunc_get_w_of_y(xs, y),
            Some(S::Irregular) | Some(S::Custom) | Some(S::StreetXsect) => {
                // Tabulated transect (legacy xsect_getWofY IRREGULAR case).
                if let Some(tb) = xs.width_tbl {
                    return xs.w_max * self.lookup(y_norm, tb, xs.transect_tbl_size);
                }
                0.0
            }
            _ => 0.0,
        }
    }

    // ============================================================================
    // Main dispatch: getRofY
    // ============================================================================

    /// Hydraulic radius from depth; RECT_CLOSED, RECT_OPEN, MOD_BASKET and the S-table
    /// shapes go through get_r_of_a(get_a_of_y).
    pub fn get_r_of_y(&self, xs: &XSectParams, y: f64) -> f64 {
        let y_norm = y / xs.y_full;
        let t = &self.tbl;

        match XSectShape::from_code(xs.type_) {
            Some(S::ForceMain) | Some(S::Circular) => xs.r_full * self.lookup(y_norm, t.R_Circ, t.N_R_Circ),
            Some(S::FilledCircular) => {
                if xs.y_bot == 0.0 {
                    return xs.r_full * self.lookup(y_norm, t.R_Circ, t.N_R_Circ);
                }
                self.filled_circ_get_r_of_y(xs, y)
            }
            Some(S::Eggshaped) => xs.r_full * self.lookup(y_norm, t.R_Egg, t.N_R_Egg),
            Some(S::Horseshoe) => xs.r_full * self.lookup(y_norm, t.R_Horseshoe, t.N_R_Horseshoe),
            Some(S::Baskethandle) => xs.r_full * self.lookup(y_norm, t.R_Baskethandle, t.N_R_Baskethandle),
            Some(S::HorizEllipse) => xs.r_full * self.lookup(y_norm, t.R_HorizEllipse, t.N_R_HorizEllipse),
            Some(S::VertEllipse) => xs.r_full * self.lookup(y_norm, t.R_VertEllipse, t.N_R_VertEllipse),
            Some(S::Arch) => xs.r_full * self.lookup(y_norm, t.R_Arch, t.N_R_Arch),
            Some(S::RectTriang) => self.rect_triang_get_r_of_y(xs, y),
            Some(S::RectRound) => self.rect_round_get_r_of_y(xs, y),
            Some(S::Trapezoidal) => self.trapez_get_r_of_y(xs, y),
            Some(S::Triangular) => self.triang_get_r_of_y(xs, y),
            Some(S::Parabolic) => self.parab_get_r_of_y(xs, y),
            Some(S::Powerfunc) => self.powerfunc_get_r_of_y(xs, y),
            Some(S::Irregular) | Some(S::Custom) | Some(S::StreetXsect) => {
                // Tabulated transect (legacy xsect_getRofY IRREGULAR case). Must be
                // explicit: the generic default would recurse forever for these shapes.
                if let Some(tb) = xs.hrad_tbl {
                    return xs.r_full * self.lookup(y_norm, tb, xs.transect_tbl_size);
                }
                0.0
            }
            // RECT_CLOSED, RECT_OPEN, MOD_BASKET, tabulated S-only shapes
            _ => self.get_r_of_a(xs, self.get_a_of_y(xs, y)),
        }
    }

    // ============================================================================
    // Main dispatch: getYofA
    // ============================================================================

    /// Depth from wetted area; 0 for a <= 0, a missing optional area table or an
    /// unsupported code.
    pub fn get_y_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        if a <= 0.0 {
            return 0.0;
        }
        let alpha = a / xs.a_full;
        let t = &self.tbl;

        match XSectShape::from_code(xs.type_) {
            Some(S::ForceMain) | Some(S::Circular) => self.circ_get_y_of_a(xs, a),
            Some(S::FilledCircular) => self.filled_circ_get_y_of_a(xs, a),
            Some(S::Eggshaped) => xs.y_full * self.lookup(alpha, t.Y_Egg, t.N_Y_Egg),
            Some(S::Horseshoe) => xs.y_full * self.lookup(alpha, t.Y_Horseshoe, t.N_Y_Horseshoe),
            Some(S::Gothic) => xs.y_full * self.lookup(alpha, t.Y_Gothic, t.N_Y_Gothic),
            Some(S::Catenary) => xs.y_full * self.lookup(alpha, t.Y_Catenary, t.N_Y_Catenary),
            Some(S::Semielliptical) => xs.y_full * self.lookup(alpha, t.Y_SemiEllip, t.N_Y_SemiEllip),
            Some(S::Baskethandle) => xs.y_full * self.lookup(alpha, t.Y_BasketHandle, t.N_Y_BasketHandle),
            Some(S::Semicircular) => xs.y_full * self.lookup(alpha, t.Y_SemiCirc, t.N_Y_SemiCirc),
            Some(S::HorizEllipse) => {
                xs.y_full
                    * self.inv_lookup(
                        alpha,
                        t.A_HorizEllipse,
                        t.N_A_HorizEllipse,
                        t.lut(LutId::AHorizEllipse),
                    )
            }
            Some(S::VertEllipse) => {
                xs.y_full
                    * self.inv_lookup(
                        alpha,
                        t.A_VertEllipse,
                        t.N_A_VertEllipse,
                        t.lut(LutId::AVertEllipse),
                    )
            }
            Some(S::Arch) => xs.y_full * self.inv_lookup(alpha, t.A_Arch, t.N_A_Arch, t.lut(LutId::AArch)),
            Some(S::RectClosed) | Some(S::RectOpen) => a / xs.w_max,
            Some(S::RectTriang) => self.rect_triang_get_y_of_a(xs, a),
            Some(S::RectRound) => self.rect_round_get_y_of_a(xs, a),
            Some(S::ModBasket) => self.mod_basket_get_y_of_a(xs, a),
            Some(S::Trapezoidal) => self.trapez_get_y_of_a(xs, a),
            Some(S::Triangular) => self.triang_get_y_of_a(xs, a),
            Some(S::Parabolic) => self.parab_get_y_of_a(xs, a),
            Some(S::Powerfunc) => self.powerfunc_get_y_of_a(xs, a),
            Some(S::Irregular) | Some(S::Custom) | Some(S::StreetXsect) => {
                // Invert the normalized area table (legacy xsect_getYofA IRREGULAR).
                if let Some(tb) = xs.area_tbl {
                    return xs.y_full * self.inv_lookup(alpha, tb, xs.transect_tbl_size, xs.area_lut);
                }
                0.0
            }
            _ => 0.0,
        }
    }

    // ============================================================================
    // Main dispatch: getSofA
    // ============================================================================

    /// Section factor from area: S tables, crown corrections, or A*R^(2/3).
    pub fn get_s_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        let alpha = a / xs.a_full;
        let t = &self.tbl;

        match XSectShape::from_code(xs.type_) {
            Some(S::ForceMain) | Some(S::Circular) => self.circ_get_s_of_a(xs, a),
            Some(S::Eggshaped) => xs.s_full * self.lookup(alpha, t.S_Egg, t.N_S_Egg),
            Some(S::Horseshoe) => xs.s_full * self.lookup(alpha, t.S_Horseshoe, t.N_S_Horseshoe),
            Some(S::Gothic) => xs.s_full * self.lookup(alpha, t.S_Gothic, t.N_S_Gothic),
            Some(S::Catenary) => xs.s_full * self.lookup(alpha, t.S_Catenary, t.N_S_Catenary),
            Some(S::Semielliptical) => xs.s_full * self.lookup(alpha, t.S_SemiEllip, t.N_S_SemiEllip),
            Some(S::Baskethandle) => xs.s_full * self.lookup(alpha, t.S_BasketHandle, t.N_S_BasketHandle),
            Some(S::Semicircular) => xs.s_full * self.lookup(alpha, t.S_SemiCirc, t.N_S_SemiCirc),
            Some(S::RectClosed) => self.rect_closed_get_s_of_a(xs, a),
            Some(S::RectOpen) => self.rect_open_get_s_of_a(xs, a),
            Some(S::RectTriang) => self.rect_triang_get_s_of_a(xs, a),
            Some(S::RectRound) => self.rect_round_get_s_of_a(xs, a),
            _ => {
                if a == 0.0 {
                    return 0.0;
                }
                let r = self.get_r_of_a(xs, a);
                if r < TINY {
                    return 0.0;
                }
                a * r.powf(2.0 / 3.0) // LIBM
            }
        }
    }

    // ============================================================================
    // getRofA — hydraulic radius from area
    // ============================================================================

    /// Hydraulic radius from area; a <= 0 gives 0. The S-table shapes invert (S/A)^(3/2).
    pub fn get_r_of_a(&self, xs: &XSectParams, a: f64) -> f64 {
        if a <= 0.0 {
            return 0.0;
        }
        match XSectShape::from_code(xs.type_) {
            Some(S::HorizEllipse)
            | Some(S::VertEllipse)
            | Some(S::Arch)
            | Some(S::Irregular)
            | Some(S::FilledCircular)
            | Some(S::Custom)
            | Some(S::StreetXsect) => self.get_r_of_y(xs, self.get_y_of_a(xs, a)),
            Some(S::RectClosed) => self.rect_closed_get_r_of_a(xs, a),
            Some(S::RectOpen) => shape::rect_open_r_of_a(a, xs.w_max, xs.s_bot),
            Some(S::RectTriang) => self.rect_triang_get_r_of_a(xs, a),
            Some(S::RectRound) => self.rect_round_get_r_of_a(xs, a),
            Some(S::ModBasket) => self.mod_basket_get_r_of_a(xs, a),
            Some(S::Trapezoidal) => self.trapez_get_r_of_a(xs, a),
            Some(S::Triangular) => self.triang_get_r_of_a(xs, a),
            Some(S::Parabolic) => self.parab_get_r_of_a(xs, a),
            Some(S::Powerfunc) => self.powerfunc_get_r_of_a(xs, a),
            _ => {
                let s = self.get_s_of_a(xs, a);
                if s < TINY || a < TINY {
                    return 0.0;
                }
                (s / a).powf(3.0 / 2.0) // LIBM
            }
        }
    }

    // ============================================================================
    // getdSdA — derivative of section factor w.r.t. area
    // ============================================================================

    /// Section-factor derivative: circular, table segment, analytic, else finite difference.
    pub fn get_ds_da(&self, xs: &XSectParams, a: f64) -> f64 {
        let t = &self.tbl;
        match XSectShape::from_code(xs.type_) {
            Some(S::ForceMain) | Some(S::Circular) => self.circ_get_ds_da(xs, a),
            Some(S::Eggshaped) => self.tabular_get_ds_da(xs, a, t.S_Egg, t.N_S_Egg),
            Some(S::Horseshoe) => self.tabular_get_ds_da(xs, a, t.S_Horseshoe, t.N_S_Horseshoe),
            Some(S::Gothic) => self.tabular_get_ds_da(xs, a, t.S_Gothic, t.N_S_Gothic),
            Some(S::Catenary) => self.tabular_get_ds_da(xs, a, t.S_Catenary, t.N_S_Catenary),
            Some(S::Semielliptical) => self.tabular_get_ds_da(xs, a, t.S_SemiEllip, t.N_S_SemiEllip),
            Some(S::Baskethandle) => self.tabular_get_ds_da(xs, a, t.S_BasketHandle, t.N_S_BasketHandle),
            Some(S::Semicircular) => self.tabular_get_ds_da(xs, a, t.S_SemiCirc, t.N_S_SemiCirc),
            Some(S::RectClosed) => self.rect_closed_get_ds_da(xs, a),
            Some(S::RectOpen) => self.rect_open_get_ds_da(xs, a),
            Some(S::RectTriang) => self.rect_triang_get_ds_da(xs, a),
            Some(S::RectRound) => self.rect_round_get_ds_da(xs, a),
            Some(S::ModBasket) => self.mod_basket_get_ds_da(xs, a),
            Some(S::Trapezoidal) => self.trapez_get_ds_da(xs, a),
            Some(S::Triangular) => self.triang_get_ds_da(xs, a),
            _ => self.generic_get_ds_da(xs, a),
        }
    }

    // ============================================================================
    // getAofS — area from section factor
    // ============================================================================

    /// Area from section factor. psi is computed before the local s_max cap (legacy).
    pub fn get_a_of_s(&self, xs: &XSectParams, s: f64) -> f64 {
        let psi = s / xs.s_full;
        if s <= 0.0 {
            return 0.0;
        }
        let s = if s > xs.s_max { xs.s_max } else { s };
        let t = &self.tbl;

        match XSectShape::from_code(xs.type_) {
            Some(S::Dummy) => 0.0,
            Some(S::ForceMain) | Some(S::Circular) => self.circ_get_a_of_s(xs, s),
            Some(S::Eggshaped) => xs.a_full * self.inv_lookup(psi, t.S_Egg, t.N_S_Egg, t.lut(LutId::SEgg)),
            Some(S::Horseshoe) => {
                xs.a_full * self.inv_lookup(psi, t.S_Horseshoe, t.N_S_Horseshoe, t.lut(LutId::SHorseshoe))
            }
            Some(S::Gothic) => {
                xs.a_full * self.inv_lookup(psi, t.S_Gothic, t.N_S_Gothic, t.lut(LutId::SGothic))
            }
            Some(S::Catenary) => {
                xs.a_full * self.inv_lookup(psi, t.S_Catenary, t.N_S_Catenary, t.lut(LutId::SCatenary))
            }
            Some(S::Semielliptical) => {
                xs.a_full * self.inv_lookup(psi, t.S_SemiEllip, t.N_S_SemiEllip, t.lut(LutId::SSemiEllip))
            }
            Some(S::Baskethandle) => {
                xs.a_full
                    * self.inv_lookup(
                        psi,
                        t.S_BasketHandle,
                        t.N_S_BasketHandle,
                        t.lut(LutId::SBasketHandle),
                    )
            }
            Some(S::Semicircular) => {
                xs.a_full * self.inv_lookup(psi, t.S_SemiCirc, t.N_S_SemiCirc, t.lut(LutId::SSemiCirc))
            }
            _ => {
                // Newton-Raphson on S(a) = s, bracketed in [a1, a2] (legacy generic_getAofS).
                // PARITY: legacy xsect_getAmax returns aBot for IRREGULAR and CUSTOM.
                let sh = xs.type_;
                let a_max = if sh == S::Custom as i32 || sh == S::Irregular as i32 {
                    xs.a_bot
                } else {
                    xs.a_full * self.get_amax(xs)
                };
                let (a1, a2) = if (s <= xs.s_max && s >= xs.s_full) && xs.s_max != xs.s_full {
                    (xs.a_full, a_max) // sFull < sMax: root lies between aFull and aMax
                } else {
                    (0.0, a_max)
                };
                let mut a = 0.5 * (a1 + a2);
                let tol = 0.0001 * xs.a_full;
                self.findroot_newton(a1, a2, &mut a, tol, |aa, f, df| {
                    *f = self.get_s_of_a(xs, aa) - s;
                    *df = self.get_ds_da(xs, aa);
                });
                a
            }
        }
    }

    // ============================================================================
    // getAmax — ratio of area at max flow to full area
    // ============================================================================

    /// Built-in Amax ratio for codes 0..25, otherwise 1.
    pub fn get_amax(&self, xs: &XSectParams) -> f64 {
        if xs.type_ >= 0 && xs.type_ <= 25 {
            return self.tbl.Amax[xs.type_ as usize];
        }
        1.0
    }

    // ============================================================================
    // getYcrit — critical depth for a given flow rate (legacy xsect_getYcrit)
    // ============================================================================

    /// Critical depth for discharge q, capped by y_full.
    pub fn get_ycrit(&self, xs: &XSectParams, q: f64) -> f64 {
        if q <= 0.0 {
            return 0.0;
        }
        let q2g = q * q / GRAVITY;
        if q2g == 0.0 {
            return 0.0;
        }

        let y: f64;
        match XSectShape::from_code(xs.type_) {
            Some(S::Dummy) => return 0.0,
            Some(S::RectOpen) | Some(S::RectClosed) => {
                y = (q2g / (xs.w_max * xs.w_max)).powf(1.0 / 3.0); // LIBM
            }
            Some(S::Triangular) => {
                y = (2.0 * q2g / (xs.s_bot * xs.s_bot)).powf(1.0 / 5.0); // LIBM
            }
            Some(S::Parabolic) => {
                y = (27.0 / 32.0 * q2g / (xs.r_bot * xs.r_bot)).powf(1.0 / 4.0);
                // LIBM
            }
            Some(S::Powerfunc) => {
                let e = 1.0 / (2.0 * xs.s_bot + 3.0);
                y = (q2g * (xs.s_bot + 1.0) / (xs.r_bot * xs.r_bot)).powf(e); // LIBM
            }
            _ => {
                // Critical flow function Q_c(yc) - qTarget (legacy getQcritical).
                let q_critical = |yc: f64, q_target: f64| -> f64 {
                    let a = self.get_a_of_y(xs, yc);
                    let w = self.get_w_of_y(xs, yc);
                    if w > 0.0 {
                        return a * (GRAVITY * a / w).sqrt() - q_target;
                    }
                    -q_target
                };

                // Initial estimate from equivalent circular conduit.
                let mut y0 = 1.01 * (q2g / xs.y_full).powf(0.25); // LIBM
                if y0 >= xs.y_full {
                    y0 = 0.97 * xs.y_full;
                }

                // Ratio of conduit area to equivalent circular area.
                let r = xs.a_full / (PI / 4.0 * xs.y_full * xs.y_full);

                if r >= 0.5 && r <= 2.0 {
                    // --- interval enumeration (legacy getYcritEnum), 25 increments
                    const N_INC: i32 = 25;
                    let dy = xs.y_full / (N_INC as f64);
                    let i1 = cvt_i32(y0 / dy);
                    let mut q0 = q_critical((i1 as f64) * dy, 0.0);
                    let mut yy: f64;
                    if q0 < q {
                        yy = xs.y_full;
                        let mut i = i1.wrapping_add(1);
                        while i <= N_INC {
                            let qc = q_critical((i as f64) * dy, 0.0);
                            if qc >= q {
                                yy = ((q - q0) / (qc - q0) + ((i - 1) as f64)) * dy;
                                break;
                            }
                            q0 = qc;
                            i += 1;
                        }
                    } else {
                        yy = 0.0;
                        let mut i = i1.wrapping_sub(1);
                        while i >= 0 {
                            let qc = q_critical((i as f64) * dy, 0.0);
                            if qc < q {
                                yy = ((q - qc) / (q0 - qc) + (i as f64)) * dy;
                                break;
                            }
                            q0 = qc;
                            i -= 1;
                        }
                    }
                    y = yy;
                } else {
                    // --- Ridder's method (legacy getYcritRidder)
                    let mut y1 = 0.0;
                    let mut y2 = 0.99 * xs.y_full;
                    let q2 = q_critical(y2, 0.0);
                    if q2 < q {
                        y = xs.y_full;
                    } else {
                        let q0 = q_critical(y0, 0.0);
                        let q1 = q_critical(0.5 * xs.y_full, 0.0);
                        if q0 > q {
                            y2 = y0;
                            if q1 < q {
                                y1 = 0.5 * xs.y_full;
                            }
                        } else {
                            y1 = y0;
                            if q1 > q {
                                y2 = 0.5 * xs.y_full;
                            }
                        }
                        y = self.findroot_ridder(y1, y2, 0.001, |yc| q_critical(yc, q));
                    }
                }
            }
        }
        cmin(y, xs.y_full)
    }

    // ============================================================================
    // isOpen — returns true for open shapes
    // ============================================================================

    /// Legacy xsect_isOpen() (Amax[type] >= 1.0): DUMMY, RECT_OPEN, TRAPEZOIDAL,
    /// TRIANGULAR, PARABOLIC, POWERFUNC, IRREGULAR, STREET_XSECT.
    pub fn is_open(&self, type_: i32) -> bool {
        matches!(
            XSectShape::from_code(type_),
            Some(S::Dummy)
                | Some(S::RectOpen)
                | Some(S::Trapezoidal)
                | Some(S::Triangular)
                | Some(S::Parabolic)
                | Some(S::Powerfunc)
                | Some(S::Irregular)
                | Some(S::StreetXsect)
        )
    }
}
