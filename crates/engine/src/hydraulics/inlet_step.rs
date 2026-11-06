// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Inlet.cpp and Inlet.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Street and channel inlet capture at run time (inlet.rs holds the state, init and the
// backflow ratios): InletSolver::computeAll (called at the end of assembleLateralInflows),
// updateStats, geomOf, and the HEC-22 kernels getEo, getFlowSpread, getGutterFlowRatio,
// getGutterAreaRatio, getSplashOverVelocity, getGrateInletCapture, getCurbInletCapture,
// getOnGradeInletCapture, getOnGradeCapturedFlow, findOnSagGrateFlows, getCurbOrificeFlow,
// findOnSagCurbFlows, getOnSagSlottedFlow, getOnSagInletCapture, getOnSagCapturedFlow,
// getCustomCapturedFlow, with SPLASH_COEFFS.
//
// Differences in form, never in result:
//   * Geom's `const XSectParams* xs` (null for street inlets) is an Option<XSectParams> held
//     by value in the Geom.
//   * The out-parameters of findOnSagGrateFlows / findOnSagCurbFlows are &mut f64.
//
// NOT YET TRANSLATED here: adjustQualInflows (quality), adjustFloodingTotals and gatherStats
// (the step's bookkeeping, outside stepRouting).

use crate::core::cmath::{cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::options::RoutingModel;
use crate::core::unit_conversion as ucf;
use crate::data::table::{table_lookup_ex, Table};
use crate::hydraulics::inlet::{CurveKind, Design, GrateType, InletSolver, InletType, ThroatAngle, UsageParams, INLET_BIG, INLET_FUDGE, MIN_RUNOFF_FLOW};
use crate::hydraulics::link;
use crate::hydraulics::xsect_batch::{XSectParams, XSectShape};
use crate::hydraulics::xsection as xsect;

/// Cubic fits to the splash-over velocity curves of HEC-22 Chart 5B (seven rows for eight
/// grate types; GENERIC is guarded in getSplashOverVelocity).
const SPLASH_COEFFS: [[f64; 4]; 7] = [
    [2.22, 4.03, 0.65, 0.06],
    [0.74, 2.44, 0.27, 0.02],
    [1.76, 3.12, 0.45, 0.03],
    [0.30, 4.85, 1.31, 0.15],
    [0.99, 2.64, 0.36, 0.03],
    [0.51, 2.34, 0.20, 0.01],
    [0.28, 2.28, 0.18, 0.01],
];

/// Grate opening ratios (HEC-22 Chart 9B), as in inlet.rs.
const GRATE_OPEN_RATIOS: [f64; 8] = [0.90, 0.80, 0.60, 0.35, 0.17, 0.34, 0.80, 1.00];

/// The street / channel geometry one inlet sees (inlet::Geom).
#[derive(Debug, Clone, Copy)]
pub struct Geom<'a> {
    pub sx: f64,
    pub sl: f64,
    pub sw: f64,
    pub a: f64,
    pub w: f64,
    pub t: f64,
    pub n: f64,
    pub nsides: i32,
    pub t_crown: f64,
    pub beta: f64,
    pub qfactor: f64,
    pub xs: Option<XSectParams<'a>>,
}

fn grate_open_ratio(d: &Design) -> f64 {
    let gt = d.grate_type;
    if gt < 0 || gt >= GrateType::Generic as i32 {
        return d.frac_open_area;
    }
    GRATE_OPEN_RATIOS[gt as usize]
}

/// getEo: HEC-22 Eq(4-4)
pub fn get_eo(sr: f64, ts: f64, w: f64) -> f64 {
    let mut x = sr / (ts / w);
    x = (1.0 + x).powf(2.67) - 1.0; // LIBM pow
    x = 1.0 + sr / x;
    1.0 / x
}

/// getFlowSpread (inlet.c:1200)
pub fn get_flow_spread(g: &Geom, q: f64) -> f64 {
    let f = g.qfactor;
    let mut ts1;
    if g.a == 0.0 {
        ts1 = (q / f).powf(0.375); // LIBM pow
    } else {
        let f1 = f * ((g.a / g.w) / g.sx).powf(1.67); // LIBM pow
        let tw = (q / f1).powf(0.375); // LIBM pow
        if tw <= g.w {
            ts1 = tw;
        } else {
            let sr = (g.sx + g.a / g.w) / g.sx;
            let mut iter = 1;
            let mut ts2 = 0.0;
            ts1 = (q / f).powf(0.375) - g.w; // LIBM pow
            if ts1 <= 0.0 {
                ts1 = tw - g.w;
            }
            while iter < 11 {
                let eo = get_eo(sr, ts1, g.w);
                let qs = (1.0 - eo) * q;
                ts2 = (qs / f).powf(0.375); // LIBM pow
                if (ts2 - ts1).abs() < 0.01 {
                    break;
                }
                ts1 = ts2;
                iter += 1;
            }
            ts1 = ts2 + g.w;
        }
    }
    cmin(ts1, g.t_crown)
}

/// getGutterFlowRatio (inlet.c:1518)
pub fn get_gutter_flow_ratio(g: &Geom, w: f64) -> f64 {
    if g.t <= w {
        1.0
    } else if g.a > 0.0 {
        get_eo(g.sw / g.sx, g.t - w, w)
    } else {
        1.0 - (1.0 - w / g.t).powf(2.67) // LIBM pow
    }
}

/// getGutterAreaRatio (inlet.c:1535)
pub fn get_gutter_area_ratio(g: &Geom, wg: f64, a: f64) -> f64 {
    if wg >= g.w {
        return 1.0;
    }
    if g.t <= wg {
        return 1.0;
    }
    if g.t <= g.w {
        return wg / g.t;
    }
    let as_ = 0.5 * (g.t - g.w) * (g.t - g.w) * g.sx;
    let ag = wg * ((g.t * g.sx) + g.a - (wg * g.sw / 2.0));
    ag / (a - as_)
}

/// getSplashOverVelocity (inlet.c:1557)
pub fn get_splash_over_velocity(grate_type: i32, l: f64) -> f64 {
    if !(0..7).contains(&grate_type) {
        return 0.0;
    }
    let c = SPLASH_COEFFS[grate_type as usize];
    c[0] + c[1] * l - c[2] * l * l + c[3] * l * l * l
}

/// getGrateInletCapture (inlet.c:1387)
pub fn get_grate_inlet_capture(d: &Design, g: &mut Geom, q: f64) -> f64 {
    let lg = d.grate_length;
    let mut wg = d.grate_width;
    let a;
    let mut eo;
    let mut qo = q;
    let mut rf = 1.0;
    let mut rs = 0.0;

    let channel = match &g.xs {
        Some(xs) => xs.type_ == XSectShape::Trapezoidal as i32 || xs.type_ == XSectShape::RectOpen as i32,
        None => false,
    };

    if channel {
        let xs = g.xs.unwrap();
        a = xsect::get_a_of_s(&xs, q / g.beta);
        let y = xsect::get_y_of_a(&xs, a);
        g.t = xsect::get_w_of_y(&xs, y);
        eo = g.beta * (y * wg).powf(1.67) / (wg + 2.0 * y).powf(0.67) / q; // LIBM pow
        if wg > 0.99 * xs.y_bot && xs.type_ == XSectShape::Trapezoidal as i32 && xs.s_bot > 0.0 {
            wg = xs.y_bot; // dead upstream too: the channel branch reads Wg no further
            g.sx = 1.0 / xs.s_bot;
            let _ = wg;
        }
    } else if g.a == 0.0 {
        a = g.t * g.t * g.sx / 2.0;
        eo = get_gutter_flow_ratio(g, wg);
        if g.t >= g.t_crown {
            qo = g.qfactor * g.t_crown.powf(2.67); // LIBM pow
        }
    } else {
        if g.t <= g.w {
            a = g.t * g.t * g.sw / 2.0;
        } else {
            a = (g.t * g.t * g.sx + g.a * g.w) / 2.0;
        }
        eo = get_gutter_flow_ratio(g, g.w);
        if eo < 1.0 {
            if g.t >= g.t_crown {
                qo = g.qfactor * g.t_crown.powf(2.67) / (1.0 - eo); // LIBM pow
            }
            eo *= get_gutter_area_ratio(g, wg, a);
        }
    }

    let v = qo / a;
    let vo = if d.grate_type < 0 || d.grate_type == GrateType::Generic as i32 { d.splash_veloc } else { get_splash_over_velocity(d.grate_type, lg) };

    if v > vo {
        rf = 1.0 - 0.09 * (v - vo);
    }
    if eo < 1.0 {
        rs = 1.0 / (1.0 + (0.15 * v.powf(1.8) / g.sx / lg.powf(2.3))); // LIBM pow
    }
    q * (rf * eo + rs * (1.0 - eo))
}

/// getCurbInletCapture (inlet.c:1477)
pub fn get_curb_inlet_capture(g: &Geom, q: f64, l: f64) -> f64 {
    let mut se = g.sx;
    let mut e = 1.0;
    if g.a > 0.0 {
        let sr = g.sw / g.sx;
        let eo = get_eo(sr, g.t - g.w, g.w);
        se = g.sx + (g.a / g.w) * eo;
    }
    let lt = 0.6 * q.powf(0.42) * g.sl.powf(0.3) * (1.0 / (g.n * se)).powf(0.6); // LIBM pow
    if l < lt {
        e = 1.0 - (l / lt);
        e = 1.0 - e.powf(1.8); // LIBM pow
    }
    e = cmin(e, 1.0);
    e = cmax(e, 0.0);
    e * q
}

/// getOnGradeInletCapture (inlet.c:1325)
pub fn get_on_grade_inlet_capture(d: &Design, g: &mut Geom, q: f64, depth: f64) -> f64 {
    let mut q1 = q;
    let mut qc = 0.0;
    if d.type_ == InletType::DropCurb as i32 {
        qc = get_on_sag_inlet_capture(d, g, depth);
        return cmin(qc, q);
    }
    if d.type_ == InletType::DropGrate as i32 {
        qc = get_grate_inlet_capture(d, g, q);
        return cmin(qc, q);
    }
    g.t = get_flow_spread(g, q);
    if d.type_ == InletType::Slotted as i32 {
        qc = get_curb_inlet_capture(g, q, d.slotted_length);
        return cmin(qc, q);
    }
    let lcurb = d.curb_length;
    let lgrate = d.grate_length;
    if lcurb > 0.0 {
        let lsweep = lcurb - lgrate;
        if lsweep > 0.0 {
            qc = get_curb_inlet_capture(g, q1, lsweep);
            q1 -= qc;
        }
    }
    if lgrate > 0.0 && q1 > 0.0 {
        if q1 != q {
            g.t = get_flow_spread(g, q1);
        }
        qc += get_grate_inlet_capture(d, g, q1);
    }
    qc
}

/// getOnGradeCapturedFlow (inlet.c:1267): inlets in series.
pub fn get_on_grade_captured_flow(d: &Design, u: &UsageParams, g: &mut Geom, q: f64, depth: f64) -> f64 {
    if u.num_inlets == 0 {
        return 0.0;
    }
    let mut q_approach = q;
    if q_approach < MIN_RUNOFF_FLOW {
        return 0.0;
    }
    q_approach /= g.nsides as f64;
    let mut q_bypassed = q_approach;
    let mut q_captured = 0.0;
    let mut q_max = INLET_BIG;
    if u.flow_limit > 0.0 {
        q_max = u.flow_limit;
    }
    for _i in 1..=u.num_inlets {
        let mut qc = get_on_grade_inlet_capture(d, g, q_bypassed, depth) * u.clog_factor;
        qc = cmin(qc, q_max);
        qc = cmin(qc, q_bypassed);
        q_captured += qc;
        q_bypassed -= qc;
        if q_bypassed < MIN_RUNOFF_FLOW {
            break;
        }
    }
    q_captured * g.nsides as f64
}

/// findOnSagGrateFlows (inlet.c:1644)
pub fn find_on_sag_grate_flows(d: &Design, g: &Geom, depth: f64, qw: &mut f64, qo: &mut f64) {
    let lg = d.grate_length;
    let mut wg = d.grate_width;
    let p;
    let di;
    if d.type_ == InletType::DropGrate as i32 {
        di = depth;
        p = 2.0 * (lg + wg);
    } else {
        if depth <= wg * g.sw {
            wg = depth / g.sw;
        }
        di = depth - (wg / 2.0) * g.sw;
        p = lg + 2.0 * wg;
    }
    let ao = lg * wg * grate_open_ratio(d);
    if depth <= 1.79 * ao / p {
        *qw = 3.0 * p * di.powf(1.5); // LIBM pow
    } else {
        *qo = 0.67 * ao * (2.0 * 32.16 * di).sqrt();
    }
}

/// getCurbOrificeFlow (inlet.c:1764)
pub fn get_curb_orifice_flow(di: f64, h: f64, l: f64, throat_angle: i32) -> f64 {
    let mut d = di;
    if throat_angle == ThroatAngle::Horizontal as i32 {
        d = di - h / 2.0;
    } else if throat_angle == ThroatAngle::Inclined as i32 {
        d = di - (h / 2.0) * 0.7071;
    }
    0.67 * h * l * (2.0 * 32.16 * d).sqrt()
}

/// findOnSagCurbFlows (inlet.c:1703)
pub fn find_on_sag_curb_flows(d: &Design, g: &Geom, depth: f64, mut l: f64, qw: &mut f64, qo: &mut f64) {
    let throat = d.curb_throat;
    let h = d.curb_height;
    let qweir;
    let dweir;
    if l <= 0.0 {
        return;
    }
    if d.type_ == InletType::DropCurb as i32 {
        l *= 4.0;
    }
    let dorif = 1.4 * h;
    if depth > dorif {
        *qo = get_curb_orifice_flow(depth, h, l, throat);
        return;
    }
    if g.a == 0.0 || l > 12.0 {
        dweir = h;
        if depth < dweir {
            *qw = 3.0 * l * depth.powf(1.5); // LIBM pow
            return;
        }
        qweir = 3.0 * l * dweir.powf(1.5); // LIBM pow
    } else {
        let p = l + 1.8 * g.w;
        dweir = h + g.a;
        if depth < dweir {
            *qw = 2.3 * p * depth.powf(1.5); // LIBM pow
            return;
        }
        qweir = 2.3 * p * dweir.powf(1.5); // LIBM pow
    }
    let qorif = get_curb_orifice_flow(dorif, h, l, throat);
    let r = (depth - dweir) / (dorif - dweir);
    *qw = (1.0 - r) * qweir;
    *qo = r * qorif;
}

/// getOnSagSlottedFlow (inlet.c:1785)
pub fn get_on_sag_slotted_flow(d: &Design, depth: f64) -> f64 {
    let l = d.slotted_length;
    let w = d.slotted_width;
    if depth <= 2.587 * w {
        return 2.48 * l * depth.powf(1.5); // LIBM pow
    }
    0.8 * l * w * (64.32 * depth).sqrt()
}

/// getOnSagInletCapture (inlet.c:1610)
pub fn get_on_sag_inlet_capture(d: &Design, g: &Geom, depth: f64) -> f64 {
    let mut qsw = 0.0;
    let mut qso = 0.0;
    let mut qgw = 0.0;
    let mut qgo = 0.0;
    let mut qcw = 0.0;
    let mut qco = 0.0;
    if d.slotted_length > 0.0 {
        return get_on_sag_slotted_flow(d, depth);
    }
    let lgrate = d.grate_length;
    if lgrate > 0.0 {
        find_on_sag_grate_flows(d, g, depth, &mut qgw, &mut qgo);
    }
    let lcurb = d.curb_length;
    if lcurb > 0.0 {
        let lsweep = lcurb - lgrate;
        if lsweep > 0.0 {
            find_on_sag_curb_flows(d, g, depth, lsweep, &mut qsw, &mut qso);
        }
        if qgo > 0.0 {
            find_on_sag_curb_flows(d, g, depth, lgrate, &mut qcw, &mut qco);
        }
    }
    let _ = qcw;
    qgw + qgo + qsw + qso + qco
}

/// getOnSagCapturedFlow (inlet.c:1574): `nsides_prev` is the previous inlet's street side
/// count, as legacy's module-level Nsides stood on entry.
pub fn get_on_sag_captured_flow(d: &Design, u: &UsageParams, g: &Geom, depth: f64, nsides_prev: i32) -> f64 {
    if u.num_inlets == 0 {
        return 0.0;
    }
    let total_inlets = nsides_prev * u.num_inlets;
    let mut q_max = INLET_BIG;
    if u.flow_limit > 0.0 {
        q_max = u.flow_limit;
    }
    let mut q_captured = get_on_sag_inlet_capture(d, g, depth.abs());
    q_captured *= u.clog_factor;
    q_captured = cmin(q_captured, q_max);
    q_captured *= total_inlets as f64;
    q_captured
}

/// getCustomCapturedFlow (inlet.c:1902)
#[allow(clippy::too_many_arguments)]
pub fn get_custom_captured_flow(u: &UsageParams, g: &Geom, curve: &Table, curve_kind: i32, q: f64, depth: f64, ucf_flow: f64, ucf_len: f64) -> f64 {
    if u.num_inlets == 0 {
        return 0.0;
    }
    let mut q_max = INLET_BIG;
    if u.flow_limit > 0.0 {
        q_max = u.flow_limit;
    }
    let sides = g.nsides;
    let mut q_bypassed = q / sides as f64;
    let mut q_captured = 0.0;
    if curve_kind == CurveKind::Diversion as i32 {
        for _j in 1..=u.num_inlets {
            let mut q_increment = u.clog_factor * table_lookup_ex(curve, q_bypassed * ucf_flow) / ucf_flow;
            q_increment = cmin(q_increment, q_max);
            q_increment = cmin(q_increment, q_bypassed);
            q_captured += q_increment;
            q_bypassed -= q_increment;
            if q_bypassed < MIN_RUNOFF_FLOW {
                break;
            }
        }
    } else if curve_kind == CurveKind::Rating as i32 {
        q_captured = u.num_inlets as f64 * u.clog_factor * table_lookup_ex(curve, depth * ucf_len) / ucf_flow;
    }
    q_captured * sides as f64
}

impl InletSolver {
    /// geomOf: the street geometry of one usage row.
    pub fn geom_of<'a>(&self, _ctx: &SimulationContext, ii: i32, xs: Option<XSectParams<'a>>) -> Geom<'a> {
        let ui = ii as usize;
        let s = &self.soa_;
        let mut g = Geom { sx: s.sx[ui], sl: s.slope[ui], sw: 0.01, a: s.gutter_depression[ui], w: s.gutter_width[ui], t: 0.0, n: s.road_roughness[ui], nsides: s.n_sides[ui], t_crown: s.t_crown[ui], beta: 0.0, qfactor: s.flow_factor[ui], xs };
        if s.local_depress[ui] * s.local_width[ui] > 0.0 {
            g.a += s.local_depress[ui];
            g.w = s.local_width[ui];
        }
        g.sw = if g.w * g.a > 0.0 { g.sx + g.a / g.w } else { g.sx };
        g.beta = if g.n > 0.0 && g.sl > 0.0 { 1.486 * g.sl.sqrt() / g.n } else { 0.0 };
        g
    }

    /// computeAll: capture, backflow and the lateral-flow coupling of every inlet.
    pub fn compute_all(&mut self, ctx: &mut SimulationContext, dt: f64) {
        let ni = self.soa_.count;
        if ni == 0 {
            return;
        }
        let nn = ctx.n_nodes();
        let nl = ctx.n_links();
        let is_dw = ctx.options.routing_model == RoutingModel::Dynwave;
        let report_on = ctx.current_date > ctx.options.report_start;

        let us = ucf::get_unit_system(ctx.options.flow_units as i32);
        let ucf_len = ucf::UCF[ucf::LENGTH as usize][us as usize];
        let ucf_flow = ucf::QCF[ctx.options.flow_units as usize];

        if (self.node_inlet_flow_.len() as i32) < nn {
            self.node_inlet_flow_.clear();
            self.node_inlet_flow_.resize(nn as usize, 0.0);
            self.node_inflow_.clear();
            self.node_inflow_.resize(nn as usize, 0.0);
        }
        self.node_inlet_flow_.iter_mut().for_each(|x| *x = 0.0);

        if !is_dw {
            for j in 0..nn as usize {
                self.node_inflow_[j] = cmax(0.0, ctx.nodes.lat_flow[j]);
            }
            for j in 0..nl as usize {
                let n2 = ctx.links.node2[j];
                if n2 >= 0 && n2 < nn {
                    self.node_inflow_[n2 as usize] += cmax(0.0, ctx.links.flow[j]);
                }
            }
        }

        // first pass: capture + backflow, in reverse row order (legacy's prepended list)
        for ii in (0..ni).rev() {
            let ui = ii as usize;
            self.soa_.flow_capture[ui] = 0.0;
            self.soa_.q_approach[ui] = 0.0;
            self.soa_.backflow[ui] = 0.0;

            let host = self.soa_.bypass_node[ui];
            let capture = self.soa_.node_idx[ui];
            if host < 0 || host >= nn {
                continue;
            }
            if capture < 0 || capture >= nn {
                continue;
            }
            let up = self.soa_.up_link[ui];
            let up2 = self.soa_.up_link2[ui];

            let mut q = 0.0;
            if self.soa_.host_node[ui] >= 0 {
                if up >= 0 && ctx.links.node2[up as usize] == host {
                    q += cmax(0.0, ctx.links.flow[up as usize]);
                }
                if self.soa_.two_inbound[ui] != 0 && up2 >= 0 {
                    q += cmax(0.0, ctx.links.flow[up2 as usize]);
                }
            } else if up >= 0 {
                q = ctx.links.flow[up as usize].abs();
            }
            self.soa_.q_approach[ui] = q;

            let depth = ctx.nodes.depth[host as usize];
            let d = self.design_of(ii);
            let u = self.usage_of(ii);

            let xsp = if up >= 0 && (d.type_ == InletType::DropGrate as i32 || d.type_ == InletType::DropCurb as i32) { Some(link::build_xsect_params(&ctx.links, up as usize, Some(&ctx.transect_tables[..]))) } else { None };
            let mut g = self.geom_of(ctx, ii, xsp);

            let mut qc = 0.0;
            if d.type_ == InletType::Custom as i32 {
                let ci = self.soa_.curve_index[ui];
                if ci >= 0 && (ci as usize) < ctx.tables.tables.len() {
                    qc = get_custom_captured_flow(&u, &g, &ctx.tables.tables[ci as usize], self.soa_.curve_kind[ui], q, depth, ucf_flow, ucf_len);
                }
            } else if self.soa_.is_sag[ui] == 0 {
                qc = get_on_grade_captured_flow(&d, &u, &mut g, q, depth);
                if u.num_inlets != 0 && q >= MIN_RUNOFF_FLOW {
                    self.nsides_prev_ = g.nsides;
                }
            } else {
                qc = get_on_sag_captured_flow(&d, &u, &g, depth, self.nsides_prev_);
                if u.num_inlets != 0 {
                    self.nsides_prev_ = g.nsides;
                }
            }
            if qc.abs() < INLET_FUDGE {
                qc = 0.0;
            }
            self.soa_.flow_capture[ui] = qc;
            self.node_inlet_flow_[host as usize] += qc;

            let mut qbf = ctx.nodes.overflow[capture as usize] * self.soa_.backflow_ratio[ui];
            if qbf.abs() < INLET_FUDGE {
                qbf = 0.0;
            }
            self.soa_.backflow[ui] = qbf;
        }

        // second pass: throttle, couple, accumulate
        for ii in (0..ni).rev() {
            let ui = ii as usize;
            let host = self.soa_.bypass_node[ui];
            let capture = self.soa_.node_idx[ui];
            if host < 0 || host >= nn {
                continue;
            }
            if capture < 0 || capture >= nn {
                continue;
            }
            let uh = host as usize;

            if !is_dw && self.soa_.is_sag[ui] != 0 && self.soa_.host_node[ui] < 0 {
                let mut qlim = ctx.nodes.volume[uh] / dt;
                qlim += cmax(self.node_inflow_[uh], 0.0);
                if self.node_inlet_flow_[uh] > qlim {
                    self.soa_.flow_capture[ui] *= qlim / self.node_inlet_flow_[uh];
                }
            }
            if self.soa_.host_node[ui] >= 0 {
                let qlim = self.soa_.q_approach[ui] + ctx.nodes.volume[uh] / dt;
                self.soa_.flow_capture[ui] = cmin(self.soa_.flow_capture[ui], qlim);
            }

            let qcap = self.soa_.flow_capture[ui];
            let qbf = self.soa_.backflow[ui];

            ctx.nodes.lat_flow[uh] -= qcap - qbf;
            ctx.nodes.lat_flow[capture as usize] += qcap;

            let qnet = qcap - qbf;
            self.soa_.stat_capture_vol[ui] += qcap * dt;
            self.soa_.stat_bypass_vol[ui] += cmax(0.0, self.soa_.q_approach[ui] - qnet) * dt;
            self.soa_.stat_backflow_vol[ui] += qbf * dt;
            if qcap > self.soa_.stat_peak_flow[ui] {
                self.soa_.stat_peak_flow[ui] = qcap;
            }
            if report_on {
                let qa = self.soa_.q_approach[ui];
                self.update_stats(ii, qa);
            }
        }
    }

    /// updateStats (legacy updateInletStats)
    pub fn update_stats(&mut self, ii: i32, q: f64) {
        let ui = ii as usize;
        let s = &mut self.soa_;
        let q_capture = s.flow_capture[ui];
        let q_backflow = s.backflow[ui];
        let q_net = q_capture - q_backflow;
        let q_bypass = q - q_net;
        let mut f_capture = 0.0;
        if q < MIN_RUNOFF_FLOW && q_backflow <= 0.0 {
            return;
        }
        s.flow_periods[ui] += 1;
        if q_net > 0.0 {
            s.capture_periods[ui] += 1;
            f_capture = q_net / q;
            f_capture = cmin(f_capture, 1.0);
            s.avg_flow_capture[ui] += f_capture;
            if q_bypass > MIN_RUNOFF_FLOW {
                s.bypass_freq[ui] += 1.0;
            }
        } else {
            s.backflow_periods[ui] += 1;
        }
        if q > s.peak_flow[ui] {
            s.peak_flow[ui] = q;
            s.peak_flow_capture[ui] = f_capture * 100.0;
        }
    }
}
