// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Link.hpp and
// src/engine/hydraulics/Link.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit
// 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Link hydraulics: velocity, Froude number, Manning conveyance, normal depth, and the
// bridge from the link store and the transect tables to XSectParams (legacy link.c and
// xsect.c).
//
// Differences in form, never in result:
//   * Names are snake_case (getVelocity -> get_velocity, translateShape -> translate_shape,
//     ...). Default arguments become explicit ones (barrels = 1, q_max = -1.0).
//   * computeConveyance's three `double&` outputs are `&mut f64`.
//   * buildXSectParams' nullable `const std::vector<TransectData>*` is an Option of a
//     slice, and the returned XSectParams borrows the tables it attaches.
//   * translateShape's `default:` arm (an XsectShape value outside the enum) cannot occur
//     with Rust's closed enum, so the match is exhaustive without it.

#![allow(clippy::field_reassign_with_default)]

use crate::core::cmath::cmax;
use crate::core::constants;
use crate::data::link::{LinkData, LinkType, XsectShape};
use crate::hydraulics::transect_types::{TransectData, N_TRANSECT_TBL};
use crate::hydraulics::xsect_batch::{XSectParams, XSectShape};
use crate::hydraulics::xsection as xsect;

/// Map a link-store section type (LinkData's XsectShape) to the kernel numbering
/// (XSectShape); the two enums order the shapes differently.
pub fn translate_shape(s: XsectShape) -> i32 {
    match s {
        XsectShape::Circular => XSectShape::Circular as i32,
        XsectShape::FilledCircular => XSectShape::FilledCircular as i32,
        XsectShape::RectClosed => XSectShape::RectClosed as i32,
        XsectShape::RectOpen => XSectShape::RectOpen as i32,
        XsectShape::Trapezoidal => XSectShape::Trapezoidal as i32,
        XsectShape::Triangular => XSectShape::Triangular as i32,
        XsectShape::Parabolic => XSectShape::Parabolic as i32,
        XsectShape::Power => XSectShape::Powerfunc as i32,
        XsectShape::Modbaskethandle => XSectShape::ModBasket as i32,
        XsectShape::Eggshaped => XSectShape::Eggshaped as i32,
        XsectShape::Horseshoe => XSectShape::Horseshoe as i32,
        XsectShape::Gothic => XSectShape::Gothic as i32,
        XsectShape::Catenary => XSectShape::Catenary as i32,
        XsectShape::Semielliptical => XSectShape::Semielliptical as i32,
        XsectShape::Baskethandle => XSectShape::Baskethandle as i32,
        XsectShape::Semicircular => XSectShape::Semicircular as i32,
        XsectShape::RectTriang => XSectShape::RectTriang as i32,
        XsectShape::RectRound => XSectShape::RectRound as i32,
        XsectShape::HorizEllipse => XSectShape::HorizEllipse as i32,
        XsectShape::VertEllipse => XSectShape::VertEllipse as i32,
        XsectShape::Arch => XSectShape::Arch as i32,
        XsectShape::Irregular => XSectShape::Irregular as i32,
        XsectShape::Custom => XSectShape::Custom as i32,
        XsectShape::ForceMain => XSectShape::ForceMain as i32,
        XsectShape::StreetXsect => XSectShape::StreetXsect as i32,
        XsectShape::Dummy => XSectShape::Dummy as i32,
    }
}

// ============================================================================
// Per-element: getVelocity
// ============================================================================

/// V = Q / (barrels * A(y)); 0 when depth <= 0.01 ft or area <= FUDGE.
pub fn get_velocity(xs: &XSectParams, flow: f64, depth: f64, barrels: i32) -> f64 {
    if depth <= 0.01 {
        return 0.0;
    }

    let q = flow / cmax(barrels, 1) as f64;
    let area = xsect::get_a_of_y(xs, depth);

    if area > constants::FUDGE {
        return q / area;
    }
    0.0
}

// ============================================================================
// Per-element: getFroude
// ============================================================================

/// Fr = |V| / sqrt(g * A/T); 0 for a dry section or a full closed conduit.
pub fn get_froude(xs: &XSectParams, v: f64, y: f64) -> f64 {
    if y <= constants::FUDGE {
        return 0.0;
    }

    // For closed conduits: return 0 if full
    if !xsect::is_open(xs.type_) && (xs.y_full - y) <= constants::FUDGE {
        return 0.0;
    }

    // Hydraulic depth D_h = A / T
    let area = xsect::get_a_of_y(xs, y);
    let width = xsect::get_w_of_y(xs, y);
    if width <= constants::FUDGE {
        return 0.0;
    }

    let d_h = area / width;
    if d_h <= 0.0 {
        return 0.0;
    }

    v.abs() / (constants::GRAVITY * d_h).sqrt()
}

// ============================================================================
// Per-element: computeConveyance
// ============================================================================

/// beta = PHI*sqrt(|slope|)/n, rough_factor = g*(n/PHI)^2, q_full = s_full*beta.
pub fn compute_conveyance(
    roughness: f64,
    slope: f64,
    s_full: f64,
    beta: &mut f64,
    rough_factor: &mut f64,
    q_full: &mut f64,
) {
    // beta = PHI * sqrt(|slope|) / n
    *beta = constants::PHI * slope.abs().sqrt() / roughness;

    // rough_factor = g * (n / PHI)^2
    let n_over_phi = roughness / constants::PHI;
    *rough_factor = constants::GRAVITY * n_over_phi * n_over_phi;

    // q_full = section_factor_at_full * beta
    *q_full = s_full * *beta;
}

// ============================================================================
// Per-element: getDepthFromFlow
// ============================================================================

/// Normal depth for flow q: A = getAofS(q/beta), y = getYofA(A). q is first capped at
/// q_max (negative disables the cap), as legacy link_getYnorm does.
pub fn get_depth_from_flow(xs: &XSectParams, beta: f64, q: f64, q_max: f64) -> f64 {
    let mut q = q;
    if beta <= 0.0 {
        return 0.0;
    }
    // legacy link_getYnorm (link.c:800): cap the flow at the conduit's qMax
    // BEFORE forming the section factor — the tabular shapes derive their
    // inverse-lookup argument from the uncapped psi inside getAofS.
    if q_max >= 0.0 && q > q_max {
        q = q_max;
    }
    if q <= 0.0 {
        return 0.0;
    }

    let s = q / beta; // section factor needed
    let a = xsect::get_a_of_s(xs, s); // area from section factor
    xsect::get_y_of_a(xs, a) // depth from area
}

// ============================================================================
// Per-element: getCapacity
// ============================================================================

/// A(depth)/a_full, or 0 when y_full <= 0 or a_full <= 0.
pub fn get_capacity(xs: &XSectParams, depth: f64) -> f64 {
    if xs.y_full <= 0.0 {
        return 0.0;
    }
    let area = xsect::get_a_of_y(xs, depth);
    if xs.a_full <= 0.0 {
        return 0.0;
    }
    area / xs.a_full
}

// ============================================================================
// Per-element: getHydPower
// ============================================================================

/// P = gamma * |Q| * |hL| in ft-lb/s (gamma = 62.4 lb/ft3).
pub fn get_hyd_power(flow: f64, head_upstream: f64, head_downstream: f64) -> f64 {
    const GAMMA: f64 = 62.4; // lb/ft³
    let h_l = head_upstream - head_downstream;
    GAMMA * flow.abs() * h_l.abs()
}

// ============================================================================
// Per-element: buildXSectParams
// ============================================================================

/// XSectParams of link `uj` from the link store; IRREGULAR, CUSTOM and STREET_XSECT links
/// get their transect tables attached when `transects` is given.
pub fn build_xsect_params<'a>(
    links: &LinkData,
    uj: usize,
    transects: Option<&'a [TransectData]>,
) -> XSectParams<'a> {
    let mut xs = XSectParams::default();
    xs.type_ = links.xsect_batch_shape[uj];
    xs.y_full = links.xsect_y_full[uj];
    xs.a_full = links.xsect_a_full[uj];
    xs.w_max = links.xsect_w_max[uj];
    xs.r_full = links.xsect_r_full[uj];
    xs.s_full = links.xsect_s_full[uj];
    xs.s_max = links.xsect_s_max[uj];
    xs.y_bot = links.xsect_y_bot[uj];
    xs.a_bot = links.xsect_a_bot[uj];
    xs.s_bot = links.xsect_s_bot[uj];
    xs.r_bot = links.xsect_r_bot[uj];

    // Attach the per-link transect tables for tabulated shapes. xsect_curve is
    // the unified index into ctx.transect_tables for IRREGULAR (transect),
    // CUSTOM (shape curve), and STREET_XSECT (street transect).
    let shape = xs.type_;
    if let Some(transects) = transects {
        if shape == XSectShape::Irregular as i32
            || shape == XSectShape::Custom as i32
            || shape == XSectShape::StreetXsect as i32
        {
            let ti = links.xsect_curve[uj];
            if ti >= 0 && (ti as usize) < transects.len() {
                let td = &transects[ti as usize];
                xs.transect = ti;
                xs.area_tbl = Some(&td.area_tbl[..]);
                xs.hrad_tbl = Some(&td.hrad_tbl[..]);
                xs.width_tbl = Some(&td.width_tbl[..]);
                xs.area_lut = Some(&td.area_lut);
                xs.transect_tbl_size = N_TRANSECT_TBL;
            }
        }
    }
    xs
}

/// Sediment-depth adjustment of a resolved FILLED_CIRCULAR conduit (its xsect_y_bot), else
/// 0.
#[inline]
pub fn filled_circular_offset_bump(links: &LinkData, j: usize, resolved: bool) -> f64 {
    if !resolved {
        return 0.0;
    }
    if links.type_[j] != LinkType::Conduit || links.xsect_shape[j] != XsectShape::FilledCircular {
        return 0.0;
    }
    links.xsect_y_bot[j]
}

// ============================================================================
// Per-element: applyTabulatedXSectParams
// ============================================================================

// Height at the LOWEST widest point, from the normalized width table
// (xsect.c:688-699 for CUSTOM, xsect.c:1347-1358 for IRREGULAR). Scan up while
// the width is still increasing; the first decrease ends it.
fn widest_point_height(td: &TransectData, y_full: f64) -> f64 {
    let mut i_max: i32 = 0;
    let mut wt_max = td.width_tbl[0];
    for i in 1..N_TRANSECT_TBL as usize {
        if td.width_tbl[i] < wt_max {
            break;
        }
        wt_max = td.width_tbl[i];
        i_max = i as i32;
    }
    y_full * i_max as f64 / (N_TRANSECT_TBL - 1) as f64
}

/// Full-flow properties of a tabulated section (IRREGULAR, CUSTOM, STREET_XSECT) from its
/// built table, with the table attached; `xs.type_` must already hold the kernel code.
pub fn apply_tabulated_xsect_params<'a>(
    xs: &mut XSectParams<'a>,
    td: &'a TransectData,
    y_full: f64,
    transect_idx: i32,
) {
    if xs.type_ == XSectShape::Custom as i32 {
        // PARITY: legacy xsect_setCustomXsectParams (xsect.c:668-700). The
        // shape curve is built at UNIT height, so the extensive properties
        // scale by y_full here rather than in buildCustomTables.
        xs.y_full = y_full;
        xs.a_full = td.a_full;
        xs.r_full = td.r_full;
        xs.w_max = td.w_max;
        xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM; xsect.c:684 (no clamp)
        xs.s_max = td.s_max * y_full * y_full * y_full.powf(2.0 / 3.0); // LIBM; xsect.c:685-686
        xs.a_bot = td.a_max * y_full * y_full;
        xs.yw_max = widest_point_height(td, xs.y_full);
    } else if xs.type_ == XSectShape::Irregular as i32 {
        // PARITY: legacy getTransectParams (xsect.c:1339-1358). The table is
        // built at true scale, so these are absolute.
        xs.y_full = td.y_full;
        xs.a_full = td.a_full;
        xs.r_full = td.r_full;
        xs.w_max = td.w_max;
        xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM; xsect.c:1343
        xs.s_max = td.s_max; // xsect.c:1344 (absolute)
        xs.a_bot = td.a_max; // xsect.c:1345
        xs.yw_max = widest_point_height(td, xs.y_full);
    } else {
        // STREET_XSECT. street::buildTransect does not compute a peak section
        // factor, so s_max stays at s_full and a_bot is left alone; yw_max is
        // pinned at full depth rather than scanned.
        xs.y_full = td.y_full;
        xs.a_full = td.a_full;
        xs.r_full = td.r_full;
        xs.w_max = td.w_max;
        xs.s_full = xs.a_full * xs.r_full.powf(2.0 / 3.0); // LIBM
        xs.s_max = xs.s_full;
        xs.yw_max = xs.y_full;
    }

    xs.transect = transect_idx;
    xs.area_tbl = Some(&td.area_tbl[..]);
    xs.hrad_tbl = Some(&td.hrad_tbl[..]);
    xs.width_tbl = Some(&td.width_tbl[..]);
    xs.area_lut = Some(&td.area_lut);
    xs.transect_tbl_size = N_TRANSECT_TBL;
}
