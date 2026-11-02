// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Culvert.hpp and
// src/engine/hydraulics/Culvert.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit
// 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Culvert inlet control, FHWA HEC-5 equations (legacy culvert.c, op for op). The
// unsubmerged Form-1 flow comes out of Ridder's method on form1Eqn: the flow is whatever
// the LAST evaluation of form1Eqn left in qc (the root itself is discarded), so the root
// finder calls the equation in legacy's exact sequence (math::find_root::ridder).
//
// Differences in form, never in result:
//   * getCoeffs/getInflow are get_coeffs/get_inflow; getInflow's `double& dqdh` and
//     `bool& controls` are `&mut`.
//   * The TCulvert workspace holds a reference to the section where upstream holds a
//     pointer, and is passed as `&mut` where upstream passes `TCulvert*`.
//   * COEFFS keeps upstream's initializer order {form, K, M, C, Y} (upstream's comment
//     above the table says "{K, M, C, Y}"; the struct's first member is `form`).

use crate::core::constants;
use crate::hydraulics::xsect_batch::XSectParams;
use crate::hydraulics::xsection as xsect;
use crate::math::find_root;

/// Culvert curve coefficients per type code — legacy culvert.c Params[]: the FHWA equation
/// FORM (1 or 2) and K, M, C, Y.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
#[allow(non_snake_case)]
pub struct CulvertCoeffs {
    pub form: f64,
    pub K: f64,
    pub M: f64,
    pub C: f64,
    pub Y: f64,
}

const fn cc(form: f64, k: f64, m: f64, c: f64, y: f64) -> CulvertCoeffs {
    CulvertCoeffs { form, K: k, M: m, C: c, Y: y }
}

/// Complete 58-entry FHWA HEC-5 culvert coefficient table.
/// Legacy reference: Params[58][5] in src/legacy/engine/culvert.c.
#[rustfmt::skip]
static COEFFS: [CulvertCoeffs; 58] = [
    // 0: unused
    cc(0.0, 0.0, 0.0, 0.0, 0.00),

    // Circular concrete
    cc(1.0, 0.0098, 2.00, 0.0398, 0.67), //  1: Square edge w/headwall
    cc(1.0, 0.0018, 2.00, 0.0292, 0.74), //  2: Groove end w/headwall
    cc(1.0, 0.0045, 2.00, 0.0317, 0.69), //  3: Groove end projecting

    // Circular Corrugated Metal Pipe
    cc(1.0, 0.0078, 2.00, 0.0379, 0.69), //  4: Headwall
    cc(1.0, 0.0210, 1.33, 0.0463, 0.75), //  5: Mitered to slope
    cc(1.0, 0.0340, 1.50, 0.0553, 0.54), //  6: Projecting

    // Circular Pipe, Beveled Ring Entrance
    cc(1.0, 0.0018, 2.50, 0.0300, 0.74), //  7: Beveled ring, 45 deg bevels
    cc(1.0, 0.0018, 2.50, 0.0243, 0.83), //  8: Beveled ring, 33.7 deg bevels

    // Rectangular Box with Flared Wingwalls
    cc(1.0, 0.026, 1.0, 0.0347, 0.81),  //  9: 30-75 deg wingwall flares
    cc(1.0, 0.061, 0.75, 0.0400, 0.80), // 10: 90 or 15 deg wingwall flares
    cc(1.0, 0.061, 0.75, 0.0423, 0.82), // 11: 0 deg wingwall flares (straight sides)

    // Rectangular Box with Flared Wingwalls & Top Edge Bevel
    cc(2.0, 0.510, 0.667, 0.0309, 0.80), // 12: 45 deg flare; 0.43D top edge bevel
    cc(2.0, 0.486, 0.667, 0.0249, 0.83), // 13: 18-33.7 deg flare; 0.083D top edge bevel

    // Rectangular Box; 90-deg Headwall; Chamfered or Beveled Inlet Edges
    cc(2.0, 0.515, 0.667, 0.0375, 0.79),  // 14: chamfered 3/4-in
    cc(2.0, 0.495, 0.667, 0.0314, 0.82),  // 15: beveled 1/2-in/ft at 45 deg (1:1)
    cc(2.0, 0.486, 0.667, 0.0252, 0.865), // 16: beveled 1-in/ft at 33.7 deg (1:1.5)

    // Rectangular Box; Skewed Headwall; Chamfered or Beveled Inlet Edges
    cc(2.0, 0.545, 0.667, 0.04505, 0.73), // 17: 3/4" chamfered edge, 45 deg skewed headwall
    cc(2.0, 0.533, 0.667, 0.0425, 0.705), // 18: 3/4" chamfered edge, 30 deg skewed headwall
    cc(2.0, 0.522, 0.667, 0.0402, 0.68),  // 19: 3/4" chamfered edge, 15 deg skewed headwall
    cc(2.0, 0.498, 0.667, 0.0327, 0.75),  // 20: 45 deg beveled edge, 10-45 deg skewed headwall

    // Rectangular box, Non-offset Flared Wingwalls; 3/4" Chamfer at Top of Inlet
    cc(2.0, 0.497, 0.667, 0.0339, 0.803), // 21: 45 deg (1:1) wingwall flare
    cc(2.0, 0.493, 0.667, 0.0361, 0.806), // 22: 18.4 deg (3:1) wingwall flare
    cc(2.0, 0.495, 0.667, 0.0386, 0.71),  // 23: 18.4 deg (3:1) wingwall flare, 30 deg inlet skew

    // Rectangular box, Offset Flared Wingwalls, Beveled Edge at Inlet Top
    cc(2.0, 0.497, 0.667, 0.0302, 0.835), // 24: 45 deg (1:1) flare, 0.042D top edge bevel
    cc(2.0, 0.495, 0.667, 0.0252, 0.881), // 25: 33.7 deg (1.5:1) flare, 0.083D top edge bevel
    cc(2.0, 0.493, 0.667, 0.0227, 0.887), // 26: 18.4 deg (3:1) flare, 0.083D top edge bevel

    // Corrugated Metal Box
    cc(1.0, 0.0083, 2.00, 0.0379, 0.69), // 27: 90 deg headwall
    cc(1.0, 0.0145, 1.75, 0.0419, 0.64), // 28: Thick wall projecting
    cc(1.0, 0.0340, 1.50, 0.0496, 0.57), // 29: Thin wall projecting

    // Horizontal Ellipse Concrete
    cc(1.0, 0.0100, 2.00, 0.0398, 0.67), // 30: Square edge w/headwall
    cc(1.0, 0.0018, 2.50, 0.0292, 0.74), // 31: Grooved end w/headwall
    cc(1.0, 0.0045, 2.00, 0.0317, 0.69), // 32: Grooved end projecting

    // Vertical Ellipse Concrete
    cc(1.0, 0.0100, 2.00, 0.0398, 0.67), // 33: Square edge w/headwall
    cc(1.0, 0.0018, 2.50, 0.0292, 0.74), // 34: Grooved end w/headwall
    cc(1.0, 0.0095, 2.00, 0.0317, 0.69), // 35: Grooved end projecting

    // Pipe Arch, 18" Corner Radius, Corrugated Metal
    cc(1.0, 0.0083, 2.00, 0.0379, 0.69), // 36: 90 deg headwall
    cc(1.0, 0.0300, 1.00, 0.0463, 0.75), // 37: Mitered to slope
    cc(1.0, 0.0340, 1.50, 0.0496, 0.57), // 38: Projecting

    // Pipe Arch, 18" Corner Radius, Corrugated Metal (cont.)
    cc(1.0, 0.0300, 1.50, 0.0496, 0.57), // 39: Projecting
    cc(1.0, 0.0088, 2.00, 0.0368, 0.68), // 40: No bevels
    cc(1.0, 0.0030, 2.00, 0.0269, 0.77), // 41: 33.7 deg bevels

    // Pipe Arch, 31" Corner Radius, Corrugated Metal
    cc(1.0, 0.0300, 1.50, 0.0496, 0.57), // 42: Projecting
    cc(1.0, 0.0088, 2.00, 0.0368, 0.68), // 43: No bevels
    cc(1.0, 0.0030, 2.00, 0.0269, 0.77), // 44: 33.7 deg bevels

    // Arch, Corrugated Metal
    cc(1.0, 0.0083, 2.00, 0.0379, 0.69), // 45: 90 deg headwall
    cc(1.0, 0.0300, 1.00, 0.0473, 0.75), // 46: Mitered to slope
    cc(1.0, 0.0340, 1.50, 0.0496, 0.57), // 47: Thin wall projecting

    // Circular Culvert
    cc(2.0, 0.534, 0.555, 0.0196, 0.90), // 48: Smooth tapered inlet throat
    cc(2.0, 0.519, 0.640, 0.0210, 0.90), // 49: Rough tapered inlet throat

    // Elliptical Inlet Face
    cc(2.0, 0.536, 0.622, 0.0368, 0.83),  // 50: Tapered inlet, beveled edges
    cc(2.0, 0.5035, 0.719, 0.0478, 0.80), // 51: Tapered inlet, square edges
    cc(2.0, 0.547, 0.800, 0.0598, 0.75),  // 52: Tapered inlet, thin edge projecting

    // Rectangular
    cc(2.0, 0.475, 0.667, 0.0179, 0.97), // 53: Tapered inlet throat

    // Rectangular Concrete
    cc(2.0, 0.560, 0.667, 0.0446, 0.85), // 54: Side tapered, less favorable edges
    cc(2.0, 0.560, 0.667, 0.0378, 0.87), // 55: Side tapered, more favorable edges

    // Rectangular Concrete
    cc(2.0, 0.500, 0.667, 0.0446, 0.65), // 56: Slope tapered, less favorable edges
    cc(2.0, 0.500, 0.667, 0.0378, 0.71), // 57: Slope tapered, more favorable edges
];
const MAX_CULVERT_CODE: i32 = 57;

/// Copy of the FHWA coefficient record for codes 1..57; all zeros for any other code.
pub fn get_coeffs(code: i32) -> CulvertCoeffs {
    if (1..=MAX_CULVERT_CODE).contains(&code) {
        return COEFFS[code as usize];
    }
    COEFFS[0]
}

// ---------------------------------------------------------------------------
// legacy culvert.c, op for op.
// ---------------------------------------------------------------------------

const BIG: f64 = 1.0e10; // legacy consts.h

struct TCulvert<'x, 'a> {
    y_full: f64, // full depth of culvert (ft)
    scf: f64,    // slope correction factor
    dqdh: f64,   // derivative of flow w.r.t. head
    qc: f64,     // unsubmerged critical flow
    kk: f64,     // coeffs. for unsubmerged flow
    mm: f64,
    ad: f64,
    h_plus: f64, // intermediate terms
    xsect: &'x XSectParams<'a>,
}

// h/yFull + 0.5*s = yc/yFull + yh/2/yFull + K[ac/aFull*sqrt(g*yh/yFull)]^M
fn form1_eqn(yc: f64, culvert: &mut TCulvert) -> f64 {
    let ac = xsect::get_a_of_y(culvert.xsect, yc);
    let wc = xsect::get_w_of_y(culvert.xsect, yc);
    let yh = ac / wc;
    culvert.qc = ac * (constants::GRAVITY * yh).sqrt();
    culvert.h_plus - yc / culvert.y_full - yh / 2.0 / culvert.y_full
        - culvert.kk * (culvert.qc / culvert.ad).powf(culvert.mm) // LIBM
}

fn get_form1_flow(h: f64, culvert: &mut TCulvert) -> f64 {
    culvert.h_plus = h / culvert.y_full + culvert.scf;
    // Ridder's method between 0.01h and h; the result is culvert->qc
    find_root::ridder(0.01 * h, h, 0.001, |yc| form1_eqn(yc, culvert));
    culvert.qc
}

fn get_unsubmerged_flow(p: &CulvertCoeffs, h: f64, culvert: &mut TCulvert) -> f64 {
    culvert.kk = p.K;
    culvert.mm = p.M;
    let arg = h / culvert.y_full / culvert.kk;
    let q = if p.form == 1.0 {
        get_form1_flow(h, culvert)
    } else {
        culvert.ad * arg.powf(1.0 / culvert.mm) // LIBM
    };
    culvert.dqdh = q / h / culvert.mm;
    q
}

fn get_submerged_flow(p: &CulvertCoeffs, h: f64, culvert: &mut TCulvert) -> f64 {
    let cc = p.C;
    let yy = p.Y;
    let arg = (h / culvert.y_full - yy + culvert.scf) / cc;
    if arg <= 0.0 {
        culvert.dqdh = 0.0;
        return BIG;
    }
    let q = arg.sqrt() * culvert.ad;
    culvert.dqdh = 0.5 * q / arg / culvert.y_full / cc;
    q
}

fn get_transition_flow(p: &CulvertCoeffs, h: f64, h1: f64, h2: f64, culvert: &mut TCulvert) -> f64 {
    let q1 = get_unsubmerged_flow(p, h1, culvert);
    let q2 = get_submerged_flow(p, h2, culvert);
    let q = q1 + (q2 - q1) * (h - h1) / (h2 - h1);
    culvert.dqdh = (q2 - q1) / (h2 - h1);
    q
}

/// Apply the FHWA inlet-capacity limit to a proposed per-barrel flow q0 at head y. Returns
/// the inlet capacity, with `controls` true and `dqdh` replaced, when it is below q0;
/// otherwise q0 (and `dqdh` untouched). Codes outside 1..57 return q0.
pub fn get_inflow(
    q0: f64,
    y: f64,
    xs: &XSectParams,
    slope: f64,
    code: i32,
    dqdh: &mut f64,
    controls: &mut bool,
) -> f64 {
    *controls = false;
    if code <= 0 || code > MAX_CULVERT_CODE {
        return q0;
    }
    let p = COEFFS[code as usize];

    let mut culvert = TCulvert {
        y_full: 0.0,
        scf: 0.0,
        dqdh: 0.0,
        qc: 0.0,
        kk: 0.0,
        mm: 0.0,
        ad: 0.0,
        h_plus: 0.0,
        xsect: xs,
    };
    culvert.y_full = xs.y_full;
    culvert.ad = xs.a_full * culvert.y_full.sqrt();

    // slope correction factor (-7 for mitered inlets, 0.5 for others)
    culvert.scf = match code {
        5 | 37 | 46 => -7.0 * slope,
        _ => 0.5 * slope,
    };

    let q: f64;
    // submerged flow (FHWA criterion Q/AD > 4)
    let y2 = culvert.y_full * (16.0 * p.C + p.Y - culvert.scf);
    if y >= y2 {
        q = get_submerged_flow(&p, y, &mut culvert);
    } else {
        // unsubmerged below the arbitrary 0.95-full limit, transition between
        let y1 = 0.95 * culvert.y_full;
        if y <= y1 {
            q = get_unsubmerged_flow(&p, y, &mut culvert);
        } else {
            q = get_transition_flow(&p, y, y1, y2, &mut culvert);
        }
    }

    if q < q0 {
        *controls = true;
        *dqdh = culvert.dqdh;
        return q;
    }
    q0
}
