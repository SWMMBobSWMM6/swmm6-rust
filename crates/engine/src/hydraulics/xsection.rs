// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/XSection.cpp and the per-element
// declarations of src/engine/hydraulics/XSectBatch.hpp (HydroCouple/openswmm.engine,
// branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the
// Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The host table binding, the bucket maps for the inverted built-in tables, the public
// per-element geometry functions (one-line forwarders onto XsectEval) and setParams, the
// port of legacy xsect_setParams.
//
// Differences in form, never in result:
//   * Function-local statics built on first use (builtinLuts, hostTables, hostEval)
//     are std::sync::OnceLock values; like upstream's, each is built once and only read
//     afterwards.
//   * Names are snake_case (getAofY -> get_a_of_y, setParams -> set_params, ...).
//   * `const double p[]` is `&[f64]`; set_params reads exactly the slots upstream reads.
//   * `static_cast<int>(std::floor(code)) - 1` is `cvt_i32(code.floor()).wrapping_sub(1)`:
//     a code at or above 2^31 gives INT_MIN on x86, and INT_MIN - 1 (undefined in C++)
//     wraps to INT_MAX, which the range test rejects, so both return -1.
//   * The nullable LocateLut pointer is an Option.

#![allow(clippy::manual_range_contains, clippy::excessive_precision)]

use std::sync::OnceLock;

use crate::core::cmath::{cmax, cvt_i32};
use crate::hydraulics::xsect_batch::{XSectParams, XSectShape};
use crate::hydraulics::xsect_kernels::{
    LutId, XsectEval, XsectTables, PI, RECT_ALFMAX, RECT_ROUND_ALFMAX, RECT_TRIANG_ALFMAX,
};
use crate::hydraulics::xsect_lookup::{build_invlookup_lut, LocateLut};
use crate::hydraulics::xsect_tables as xt;

// ============================================================================
// Bucket maps for the inverted built-in tables (plan XSECT_LOOKUP_ACCEL A1).
// One block, built on first use and never mutated afterwards.
// ============================================================================

fn builtin_luts() -> &'static [LocateLut] {
    static LUTS: OnceLock<[LocateLut; LutId::Count as usize]> = OnceLock::new();
    LUTS.get_or_init(|| {
        let mut a = [LocateLut::default(); LutId::Count as usize];
        let mut set = |id: LutId, t: &[f64], n: i32| {
            build_invlookup_lut(&mut a[id as usize], Some(t), n);
        };
        set(LutId::YGothic, &xt::Y_Gothic, xt::N_Y_Gothic);
        set(LutId::YCatenary, &xt::Y_Catenary, xt::N_Y_Catenary);
        set(LutId::YSemiEllip, &xt::Y_SemiEllip, xt::N_Y_SemiEllip);
        set(LutId::YSemiCirc, &xt::Y_SemiCirc, xt::N_Y_SemiCirc);
        set(LutId::AHorizEllipse, &xt::A_HorizEllipse, xt::N_A_HorizEllipse);
        set(LutId::AVertEllipse, &xt::A_VertEllipse, xt::N_A_VertEllipse);
        set(LutId::AArch, &xt::A_Arch, xt::N_A_Arch);
        set(LutId::SCirc, &xt::S_Circ, xt::N_S_Circ);
        set(LutId::SEgg, &xt::S_Egg, xt::N_S_Egg);
        set(LutId::SHorseshoe, &xt::S_Horseshoe, xt::N_S_Horseshoe);
        set(LutId::SGothic, &xt::S_Gothic, xt::N_S_Gothic);
        set(LutId::SCatenary, &xt::S_Catenary, xt::N_S_Catenary);
        set(LutId::SSemiEllip, &xt::S_SemiEllip, xt::N_S_SemiEllip);
        set(LutId::SBasketHandle, &xt::S_BasketHandle, xt::N_S_BasketHandle);
        set(LutId::SSemiCirc, &xt::S_SemiCirc, xt::N_S_SemiCirc);
        a
    })
}

/// The host binding of every built-in geometry table, its length and the inverse maps.
pub fn host_tables() -> &'static XsectTables {
    static T: OnceLock<XsectTables> = OnceLock::new();
    T.get_or_init(|| XsectTables {
        A_Arch: &xt::A_Arch,
        A_Baskethandle: &xt::A_Baskethandle,
        A_Circ: &xt::A_Circ,
        A_Egg: &xt::A_Egg,
        A_HorizEllipse: &xt::A_HorizEllipse,
        A_Horseshoe: &xt::A_Horseshoe,
        A_VertEllipse: &xt::A_VertEllipse,
        R_Arch: &xt::R_Arch,
        R_Baskethandle: &xt::R_Baskethandle,
        R_Circ: &xt::R_Circ,
        R_Egg: &xt::R_Egg,
        R_HorizEllipse: &xt::R_HorizEllipse,
        R_Horseshoe: &xt::R_Horseshoe,
        R_VertEllipse: &xt::R_VertEllipse,
        S_BasketHandle: &xt::S_BasketHandle,
        S_Catenary: &xt::S_Catenary,
        S_Circ: &xt::S_Circ,
        S_Egg: &xt::S_Egg,
        S_Gothic: &xt::S_Gothic,
        S_Horseshoe: &xt::S_Horseshoe,
        S_SemiCirc: &xt::S_SemiCirc,
        S_SemiEllip: &xt::S_SemiEllip,
        W_Arch: &xt::W_Arch,
        W_BasketHandle: &xt::W_BasketHandle,
        W_Catenary: &xt::W_Catenary,
        W_Circ: &xt::W_Circ,
        W_Egg: &xt::W_Egg,
        W_Gothic: &xt::W_Gothic,
        W_HorizEllipse: &xt::W_HorizEllipse,
        W_Horseshoe: &xt::W_Horseshoe,
        W_SemiCirc: &xt::W_SemiCirc,
        W_SemiEllip: &xt::W_SemiEllip,
        W_VertEllipse: &xt::W_VertEllipse,
        Y_BasketHandle: &xt::Y_BasketHandle,
        Y_Catenary: &xt::Y_Catenary,
        Y_Circ: &xt::Y_Circ,
        Y_Egg: &xt::Y_Egg,
        Y_Gothic: &xt::Y_Gothic,
        Y_Horseshoe: &xt::Y_Horseshoe,
        Y_SemiCirc: &xt::Y_SemiCirc,
        Y_SemiEllip: &xt::Y_SemiEllip,
        Amax: &xt::Amax,
        luts: Some(builtin_luts()),
        N_A_Arch: xt::N_A_Arch,
        N_A_Baskethandle: xt::N_A_Baskethandle,
        N_A_Circ: xt::N_A_Circ,
        N_A_Egg: xt::N_A_Egg,
        N_A_HorizEllipse: xt::N_A_HorizEllipse,
        N_A_Horseshoe: xt::N_A_Horseshoe,
        N_A_VertEllipse: xt::N_A_VertEllipse,
        N_R_Arch: xt::N_R_Arch,
        N_R_Baskethandle: xt::N_R_Baskethandle,
        N_R_Circ: xt::N_R_Circ,
        N_R_Egg: xt::N_R_Egg,
        N_R_HorizEllipse: xt::N_R_HorizEllipse,
        N_R_Horseshoe: xt::N_R_Horseshoe,
        N_R_VertEllipse: xt::N_R_VertEllipse,
        N_S_BasketHandle: xt::N_S_BasketHandle,
        N_S_Catenary: xt::N_S_Catenary,
        N_S_Circ: xt::N_S_Circ,
        N_S_Egg: xt::N_S_Egg,
        N_S_Gothic: xt::N_S_Gothic,
        N_S_Horseshoe: xt::N_S_Horseshoe,
        N_S_SemiCirc: xt::N_S_SemiCirc,
        N_S_SemiEllip: xt::N_S_SemiEllip,
        N_W_Arch: xt::N_W_Arch,
        N_W_BasketHandle: xt::N_W_BasketHandle,
        N_W_Catenary: xt::N_W_Catenary,
        N_W_Circ: xt::N_W_Circ,
        N_W_Egg: xt::N_W_Egg,
        N_W_Gothic: xt::N_W_Gothic,
        N_W_HorizEllipse: xt::N_W_HorizEllipse,
        N_W_Horseshoe: xt::N_W_Horseshoe,
        N_W_SemiCirc: xt::N_W_SemiCirc,
        N_W_SemiEllip: xt::N_W_SemiEllip,
        N_W_VertEllipse: xt::N_W_VertEllipse,
        N_Y_BasketHandle: xt::N_Y_BasketHandle,
        N_Y_Catenary: xt::N_Y_Catenary,
        N_Y_Circ: xt::N_Y_Circ,
        N_Y_Egg: xt::N_Y_Egg,
        N_Y_Gothic: xt::N_Y_Gothic,
        N_Y_Horseshoe: xt::N_Y_Horseshoe,
        N_Y_SemiCirc: xt::N_Y_SemiCirc,
        N_Y_SemiEllip: xt::N_Y_SemiEllip,
    })
}

/// The evaluator bound to the host tables.
pub fn host_eval() -> &'static XsectEval {
    static E: OnceLock<XsectEval> = OnceLock::new();
    E.get_or_init(|| XsectEval { tbl: *host_tables() })
}

pub fn locate(y: f64, table: &[f64], j_last: i32) -> i32 {
    host_eval().locate(y, table, j_last)
}
pub fn lookup(x: f64, table: &[f64], n: i32) -> f64 {
    host_eval().lookup(x, table, n)
}
pub fn inv_lookup(y: f64, table: &[f64], n: i32, lut: Option<&LocateLut>) -> f64 {
    host_eval().inv_lookup(y, table, n, lut)
}
pub fn get_ycircular(alpha: f64) -> f64 {
    host_eval().get_ycircular(alpha)
}
pub fn get_scircular(alpha: f64) -> f64 {
    host_eval().get_scircular(alpha)
}

pub fn get_a_of_y(xs: &XSectParams, y: f64) -> f64 {
    host_eval().get_a_of_y(xs, y)
}
pub fn get_w_of_y(xs: &XSectParams, y: f64) -> f64 {
    host_eval().get_w_of_y(xs, y)
}
pub fn get_r_of_y(xs: &XSectParams, y: f64) -> f64 {
    host_eval().get_r_of_y(xs, y)
}
pub fn get_y_of_a(xs: &XSectParams, a: f64) -> f64 {
    host_eval().get_y_of_a(xs, a)
}
pub fn get_s_of_a(xs: &XSectParams, a: f64) -> f64 {
    host_eval().get_s_of_a(xs, a)
}
pub fn get_r_of_a(xs: &XSectParams, a: f64) -> f64 {
    host_eval().get_r_of_a(xs, a)
}
pub fn get_ds_da(xs: &XSectParams, a: f64) -> f64 {
    host_eval().get_ds_da(xs, a)
}
pub fn get_a_of_s(xs: &XSectParams, s: f64) -> f64 {
    host_eval().get_a_of_s(xs, s)
}
pub fn get_amax(xs: &XSectParams) -> f64 {
    host_eval().get_amax(xs)
}
pub fn get_ycrit(xs: &XSectParams, q: f64) -> f64 {
    host_eval().get_ycrit(xs, q)
}
pub fn is_open(type_: i32) -> bool {
    host_eval().is_open(type_)
}

// ============================================================================
// setParams — compute full-depth properties from raw geometry inputs.
//
// Port of legacy xsect_setParams. Returns 0 on success, -1 on invalid geometry (some
// fields may already have changed). IRREGULAR / CUSTOM / STREET_XSECT / DUMMY and
// unknown codes change only `type_`.
// ============================================================================

/// Derive full-flow geometry from authored dimensions `p` (lengths divided by `ucf`).
pub fn set_params(xs: &mut XSectParams, type_: i32, p: &[f64], ucf: f64) -> i32 {
    xs.type_ = type_;
    let ev = host_eval();

    match XSectShape::from_code(type_) {
        Some(XSectShape::Circular) => {
            let d = p[0] / ucf;
            xs.y_full = d;
            xs.w_max = d;
            xs.yw_max = d / 2.0;
            xs.a_full = PI / 4.0 * d * d;
            xs.r_full = d / 4.0;
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            xs.s_max = 1.08 * xs.s_full;
            0
        }
        Some(XSectShape::ForceMain) => {
            let d = p[0] / ucf;
            xs.y_full = d;
            xs.w_max = d;
            xs.yw_max = d / 2.0;
            xs.a_full = PI / 4.0 * d * d;
            xs.r_full = d / 4.0;
            xs.s_full = xs.a_full * xs.r_full.powf(0.63); // LIBM; Hazen-Williams
            xs.s_max = 1.06949 * xs.s_full;
            xs.r_bot = p[1]; // C-factor / roughness stored in rBot
            0
        }
        Some(XSectShape::RectClosed) => {
            let y = p[0] / ucf;
            let w = p[1] / ucf;
            xs.y_full = y;
            xs.w_max = w;
            xs.yw_max = y;
            xs.a_full = y * w;
            xs.r_full = xs.a_full / (2.0 * (y + w));
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            let a_max = RECT_ALFMAX * xs.a_full;
            xs.s_max = a_max * ev.rect_closed_get_r_of_a(xs, a_max).powf(2.0 / 3.0); // LIBM
            0
        }
        Some(XSectShape::RectOpen) => {
            let y = p[0] / ucf;
            let w = p[1] / ucf;
            xs.y_full = y;
            xs.w_max = w;
            xs.yw_max = y;
            xs.s_bot = p[2]; // # of sides to ignore (0, 1, or 2)
            xs.a_full = y * w;
            xs.r_full = xs.a_full / ((2.0 - xs.s_bot) * y + w);
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            xs.s_max = xs.s_full;
            0
        }
        Some(XSectShape::Trapezoidal) => {
            let y = p[0] / ucf;
            let yb = p[1] / ucf; // bottom width
            let m1 = p[2];
            let m2 = p[3]; // side slopes
            xs.y_full = y;
            xs.yw_max = y;
            xs.y_bot = yb;
            xs.s_bot = (m1 + m2) / 2.0;
            xs.r_bot = (1.0 + m1 * m1).sqrt() + (1.0 + m2 * m2).sqrt();
            xs.w_max = yb + y * (m1 + m2);
            xs.a_full = (yb + xs.s_bot * y) * y;
            xs.r_full = xs.a_full / (yb + y * xs.r_bot);
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            xs.s_max = xs.s_full;
            0
        }
        Some(XSectShape::Triangular) => {
            let y = p[0] / ucf;
            let w = p[1] / ucf;
            xs.y_full = y;
            xs.w_max = w;
            xs.yw_max = y;
            xs.s_bot = w / y / 2.0;
            xs.r_bot = (1.0 + xs.s_bot * xs.s_bot).sqrt();
            xs.a_full = y * y * xs.s_bot;
            xs.r_full = xs.a_full / (2.0 * y * xs.r_bot);
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            xs.s_max = xs.s_full;
            0
        }
        Some(XSectShape::Parabolic) => {
            let y = p[0] / ucf;
            let w = p[1] / ucf;
            xs.y_full = y;
            xs.w_max = w;
            xs.yw_max = y;
            xs.r_bot = w / 2.0 / y.sqrt();
            xs.a_full = (2.0 / 3.0) * y * w;
            xs.r_full = ev.parab_get_r_of_y(xs, y);
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            xs.s_max = xs.s_full;
            0
        }
        Some(XSectShape::Powerfunc) => {
            let y = p[0] / ucf;
            let w = p[1] / ucf;
            xs.y_full = y;
            xs.w_max = w;
            xs.yw_max = y;
            xs.s_bot = 1.0 / p[2];
            xs.r_bot = w / (xs.s_bot + 1.0) / y.powf(xs.s_bot); // LIBM
            xs.a_full = y * w / (xs.s_bot + 1.0);
            xs.r_full = ev.powerfunc_get_r_of_y(xs, y);
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            xs.s_max = xs.s_full;
            0
        }
        Some(XSectShape::FilledCircular) => {
            if p[1] >= p[0] {
                return -1; // fill must be below crown
            }
            xs.y_full = p[0] / ucf;
            xs.w_max = xs.y_full;
            xs.a_full = PI / 4.0 * xs.y_full * xs.y_full;
            xs.r_full = 0.25 * xs.y_full;
            xs.y_bot = p[1] / ucf;
            xs.a_bot = ev.circ_get_a_of_y(xs, xs.y_bot);
            xs.s_bot = get_w_of_y(xs, xs.y_bot); // type==FILLED here (legacy)
            xs.r_bot = xs.a_bot / (xs.r_full * lookup(xs.y_bot / xs.y_full, &xt::R_Circ, xt::N_R_Circ));
            xs.a_full -= xs.a_bot;
            xs.r_full = xs.a_full / (PI * xs.y_full - xs.r_bot + xs.s_bot);
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            xs.s_max = 1.08 * xs.s_full;
            xs.y_full -= xs.y_bot;
            xs.yw_max = 0.5 * xs.y_full;
            0
        }

        // --- Fixed-geometry tabulated shapes (full props from constants) ---
        Some(XSectShape::Eggshaped) => {
            xs.y_full = p[0] / ucf;
            xs.a_full = 0.5105 * xs.y_full * xs.y_full;
            xs.r_full = 0.1931 * xs.y_full;
            xs.w_max = 2.0 / 3.0 * xs.y_full;
            xs.yw_max = 0.64 * xs.y_full;
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            xs.s_max = 1.065 * xs.s_full;
            0
        }
        Some(XSectShape::Horseshoe) => {
            xs.y_full = p[0] / ucf;
            xs.a_full = 0.8293 * xs.y_full * xs.y_full;
            xs.r_full = 0.2538 * xs.y_full;
            xs.w_max = 1.0 * xs.y_full;
            xs.yw_max = 0.5 * xs.y_full;
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            xs.s_max = 1.077 * xs.s_full;
            0
        }
        Some(XSectShape::Gothic) => {
            xs.y_full = p[0] / ucf;
            xs.a_full = 0.6554 * xs.y_full * xs.y_full;
            xs.r_full = 0.2269 * xs.y_full;
            xs.w_max = 0.84 * xs.y_full;
            xs.yw_max = 0.45 * xs.y_full;
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            xs.s_max = 1.065 * xs.s_full;
            0
        }
        Some(XSectShape::Catenary) => {
            xs.y_full = p[0] / ucf;
            xs.a_full = 0.70277 * xs.y_full * xs.y_full;
            xs.r_full = 0.23172 * xs.y_full;
            xs.w_max = 0.9 * xs.y_full;
            xs.yw_max = 0.25 * xs.y_full;
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            xs.s_max = 1.05 * xs.s_full;
            0
        }
        Some(XSectShape::Semielliptical) => {
            xs.y_full = p[0] / ucf;
            xs.a_full = 0.785 * xs.y_full * xs.y_full;
            xs.r_full = 0.242 * xs.y_full;
            xs.w_max = 1.0 * xs.y_full;
            xs.yw_max = 0.15 * xs.y_full;
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            xs.s_max = 1.045 * xs.s_full;
            0
        }
        Some(XSectShape::Baskethandle) => {
            xs.y_full = p[0] / ucf;
            xs.a_full = 0.7862 * xs.y_full * xs.y_full;
            xs.r_full = 0.2464 * xs.y_full;
            xs.w_max = 0.944 * xs.y_full;
            xs.yw_max = 0.2 * xs.y_full;
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            xs.s_max = 1.06078 * xs.s_full;
            0
        }
        Some(XSectShape::Semicircular) => {
            xs.y_full = p[0] / ucf;
            xs.a_full = 1.2697 * xs.y_full * xs.y_full;
            xs.r_full = 0.2946 * xs.y_full;
            xs.w_max = 1.64 * xs.y_full;
            xs.yw_max = 0.15 * xs.y_full;
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            xs.s_max = 1.06637 * xs.s_full;
            0
        }

        Some(XSectShape::RectTriang) => {
            if p[1] <= 0.0 || p[2] <= 0.0 {
                return -1;
            }
            xs.y_full = p[0] / ucf;
            xs.w_max = p[1] / ucf;
            xs.y_bot = p[2] / ucf; // triangle height
            xs.yw_max = xs.y_full;
            xs.a_bot = xs.y_bot * xs.w_max / 2.0;
            xs.s_bot = xs.w_max / xs.y_bot / 2.0;
            xs.r_bot = (1.0 + xs.s_bot * xs.s_bot).sqrt();
            xs.a_full = xs.w_max * (xs.y_full - xs.y_bot / 2.0);
            xs.r_full = xs.a_full / (2.0 * xs.y_bot * xs.r_bot + 2.0 * (xs.y_full - xs.y_bot) + xs.w_max);
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            let a_max = RECT_TRIANG_ALFMAX * xs.a_full;
            xs.s_max = a_max * ev.rect_triang_get_r_of_a(xs, a_max).powf(2.0 / 3.0); // LIBM
            0
        }
        Some(XSectShape::RectRound) => {
            if p[1] <= 0.0 {
                return -1;
            }
            xs.y_full = p[0] / ucf;
            xs.w_max = p[1] / ucf;
            xs.r_bot = cmax(p[2], p[1] / 2.0) / ucf; // invert radius
            let theta = 2.0 * (xs.w_max / 2.0 / xs.r_bot).asin(); // LIBM
            xs.a_bot = xs.r_bot * xs.r_bot / 2.0 * (theta - theta.sin()); // LIBM
            xs.s_bot = PI * xs.r_bot * xs.r_bot * (xs.r_bot / 2.0).powf(2.0 / 3.0); // LIBM
            xs.y_bot = xs.r_bot * (1.0 - (theta / 2.0).cos()); // LIBM
            if xs.y_bot > xs.y_full {
                return -1;
            }
            xs.yw_max = xs.y_full;
            xs.a_full = xs.w_max * (xs.y_full - xs.y_bot) + xs.a_bot;
            xs.r_full = xs.a_full / (xs.r_bot * theta + 2.0 * (xs.y_full - xs.y_bot) + xs.w_max);
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            let a_max = RECT_ROUND_ALFMAX * xs.a_full;
            xs.s_max = a_max * ev.rect_round_get_r_of_a(xs, a_max).powf(2.0 / 3.0); // LIBM
            0
        }
        Some(XSectShape::ModBasket) => {
            if p[1] <= 0.0 {
                return -1;
            }
            xs.y_full = p[0] / ucf;
            xs.w_max = p[1] / ucf;
            xs.r_bot = cmax(p[2], p[1] / 2.0) / ucf; // arc radius
            let theta = 2.0 * (xs.w_max / 2.0 / xs.r_bot).asin(); // LIBM
            xs.s_bot = theta;
            xs.y_bot = xs.r_bot * (1.0 - (theta / 2.0).cos()); // LIBM
            if xs.y_bot > xs.y_full {
                return -1;
            }
            xs.yw_max = xs.y_full - xs.y_bot;
            xs.a_bot = xs.r_bot * xs.r_bot / 2.0 * (theta - theta.sin()); // LIBM
            xs.a_full = (xs.y_full - xs.y_bot) * xs.w_max + xs.a_bot;
            xs.r_full = xs.a_full / (xs.r_bot * theta + 2.0 * (xs.y_full - xs.y_bot) + xs.w_max);
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            xs.s_max = get_s_of_a(xs, get_amax(xs) * xs.a_full);
            0
        }

        Some(XSectShape::HorizEllipse) => {
            let mut code = p[2];
            if p[1] == 0.0 {
                code = p[0];
            }
            if code > 0.0 {
                // standard ellipse pipe
                let idx = cvt_i32(code.floor()).wrapping_sub(1);
                if idx < 0 || idx >= xt::NumCodesEllipse {
                    return -1;
                }
                let iu = idx as usize;
                xs.y_full = xt::MinorAxis_Ellipse[iu] / 12.0;
                xs.w_max = xt::MajorAxis_Ellipse[iu] / 12.0;
                xs.a_full = xt::Afull_Ellipse[iu];
                xs.r_full = xt::Rfull_Ellipse[iu];
            } else {
                if p[1] < 0.0 {
                    return -1;
                }
                xs.y_full = p[0] / ucf;
                xs.w_max = p[1] / ucf;
                xs.a_full = 1.2692 * xs.y_full * xs.y_full;
                xs.r_full = 0.3061 * xs.y_full;
            }
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            xs.s_max = xs.s_full;
            xs.yw_max = 0.48 * xs.y_full;
            0
        }
        Some(XSectShape::VertEllipse) => {
            let mut code = p[2];
            if p[1] == 0.0 {
                code = p[0];
            }
            if code > 0.0 {
                // standard ellipse pipe
                let idx = cvt_i32(code.floor()).wrapping_sub(1);
                if idx < 0 || idx >= xt::NumCodesEllipse {
                    return -1;
                }
                let iu = idx as usize;
                xs.y_full = xt::MajorAxis_Ellipse[iu] / 12.0;
                xs.w_max = xt::MinorAxis_Ellipse[iu] / 12.0;
                xs.a_full = xt::Afull_Ellipse[iu];
                xs.r_full = xt::Rfull_Ellipse[iu];
            } else {
                if p[1] < 0.0 {
                    return -1;
                }
                xs.y_full = p[0] / ucf;
                xs.w_max = p[1] / ucf;
                xs.a_full = 1.2692 * xs.w_max * xs.w_max;
                xs.r_full = 0.3061 * xs.w_max;
            }
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            xs.s_max = xs.s_full;
            xs.yw_max = 0.48 * xs.y_full;
            0
        }
        Some(XSectShape::Arch) => {
            let mut code = p[2];
            if p[1] == 0.0 {
                code = p[0];
            }
            if code > 0.0 {
                // standard arch pipe
                let idx = cvt_i32(code.floor()).wrapping_sub(1);
                if idx < 0 || idx >= xt::NumCodesArch {
                    return -1;
                }
                let iu = idx as usize;
                xs.y_full = xt::Yfull_Arch[iu] / 12.0; // inches
                xs.w_max = xt::Wmax_Arch[iu] / 12.0;
                xs.a_full = xt::Afull_Arch[iu];
                xs.r_full = xt::Rfull_Arch[iu];
            } else {
                if p[1] < 0.0 {
                    return -1;
                }
                xs.y_full = p[0] / ucf;
                xs.w_max = p[1] / ucf;
                xs.a_full = 0.7879 * xs.y_full * xs.w_max;
                xs.r_full = 0.2991 * xs.y_full;
            }
            xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
            xs.s_max = xs.s_full;
            xs.yw_max = 0.28 * xs.y_full;
            0
        }

        // IRREGULAR / CUSTOM / STREET_XSECT / DUMMY — geometry comes from
        // transect/curve/street tables resolved externally; leave fields.
        _ => 0,
    }
}
