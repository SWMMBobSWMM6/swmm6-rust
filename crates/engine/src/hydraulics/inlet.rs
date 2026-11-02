// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Inlet.hpp and Inlet.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Street and channel inlets (HEC-22 capture): the enums, the per-usage SoA, the design /
// usage records, getInletArea and grateOpenRatio, and InletSolver::init with
// computeBackflowRatios and the per-row accessors designOf / usageOf.
//
// NOT YET TRANSLATED (run-time step functions, with the routing step): getEo,
// getFlowSpread, getGutterFlowRatio, getGutterAreaRatio, getSplashOverVelocity,
// getGrateInletCapture, getCurbInletCapture, getOnGradeInletCapture,
// getOnGradeCapturedFlow, findOnSagGrateFlows, findOnSagCurbFlows, getCurbOrificeFlow,
// getOnSagSlottedFlow, getOnSagInletCapture, getOnSagCapturedFlow,
// getCustomCapturedFlow, InletSolver::geomOf, computeAll, adjustQualInflows,
// adjustFloodingTotals, gatherStats, updateStats, and SPLASH_COEFFS.
//
// Differences in form, never in result:
//   * Members keep upstream's names and are `pub` (upstream's private ones too) so the
//     gate's state dumper can read them.

use crate::core::constants;
use crate::core::context::SimulationContext;
use crate::core::unit_conversion as ucf;
use crate::data::link::LinkType;
use crate::data::node::NodeType;

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InletType {
    Grate = 0,
    Curb = 1,
    Combo = 2,
    Slotted = 3,
    DropGrate = 4,
    DropCurb = 5,
    Custom = 6,
}

#[allow(non_camel_case_types)]
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GrateType {
    PBar50 = 0,
    PBar50x100 = 1,
    PBar30 = 2,
    CurvedVane = 3,
    TiltBar45 = 4,
    TiltBar30 = 5,
    Reticuline = 6,
    Generic = 7,
}

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Placement {
    Automatic = 0,
    OnGrade = 1,
    OnSag = 2,
}

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThroatAngle {
    Horizontal = 0,
    Inclined = 1,
    Vertical = 2,
}

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveKind {
    Unresolved = 0,
    Diversion = 1,
    Rating = 2,
}

pub const MIN_RUNOFF_FLOW: f64 = 0.001;
pub const INLET_FUDGE: f64 = 0.0001;
pub const INLET_BIG: f64 = 1.0e10;

/// Grate opening ratios (HEC-22 Chart 9B, inlet.c:162-170)
const GRATE_OPEN_RATIOS: [f64; 8] = [
    0.90, // P_BAR-50
    0.80, // P_BAR-50x100
    0.60, // P_BAR-30
    0.35, // CURVED_VANE
    0.17, // TILT_BAR-45 (assumed)
    0.34, // TILT_BAR-30
    0.80, // RETICULINE
    1.00, // GENERIC
];

/// One inlet design, in internal units.
#[derive(Debug, Clone, Copy)]
pub struct Design {
    pub type_: i32,
    pub grate_type: i32,
    pub grate_length: f64,
    pub grate_width: f64,
    pub frac_open_area: f64,
    pub splash_veloc: f64,
    pub curb_length: f64,
    pub curb_height: f64,
    pub curb_throat: i32,
    pub slotted_length: f64,
    pub slotted_width: f64,
}

impl Default for Design {
    fn default() -> Self {
        Design {
            type_: 0,
            grate_type: GrateType::Generic as i32,
            grate_length: 0.0,
            grate_width: 0.0,
            frac_open_area: 0.0,
            splash_veloc: 0.0,
            curb_length: 0.0,
            curb_height: 0.0,
            curb_throat: ThroatAngle::Vertical as i32,
            slotted_length: 0.0,
            slotted_width: 0.0,
        }
    }
}

/// Per-usage multipliers.
#[derive(Debug, Clone, Copy)]
pub struct UsageParams {
    pub num_inlets: i32,
    pub clog_factor: f64,
    pub flow_limit: f64,
}

impl Default for UsageParams {
    fn default() -> Self {
        UsageParams { num_inlets: 1, clog_factor: 1.0, flow_limit: 0.0 }
    }
}

/// Open-area ratio of a grate design.
fn grate_open_ratio(d: &Design) -> f64 {
    let gt = d.grate_type;
    if gt < 0 || gt >= GrateType::Generic as i32 {
        return d.frac_open_area;
    }
    GRATE_OPEN_RATIOS[gt as usize]
}

/// Total effective opening area of one usage (inlet.c:1870).
pub fn get_inlet_area(d: &Design, u: &UsageParams) -> f64 {
    let mut area = 0.0;
    if d.grate_length > 0.0 {
        area = d.grate_length * d.grate_width * grate_open_ratio(d);
    }

    let curb_length = d.curb_length - d.grate_length;
    if curb_length > 0.0 {
        area += curb_length * d.curb_height;
    }

    if d.slotted_length > 0.0 {
        area = d.slotted_length * d.slotted_width;
    }

    area * u.num_inlets as f64 * u.clog_factor
}

#[derive(Debug, Clone, Default)]
pub struct InletSoA {
    pub count: i32,
    pub link_idx: Vec<i32>,
    pub host_node: Vec<i32>,
    pub node_idx: Vec<i32>,
    pub bypass_node: Vec<i32>,
    pub up_link: Vec<i32>,
    pub up_link2: Vec<i32>,
    pub two_inbound: Vec<u8>,
    pub inlet_type: Vec<i32>,
    pub grate_type: Vec<i32>,
    pub grate_length: Vec<f64>,
    pub grate_width: Vec<f64>,
    pub open_area: Vec<f64>,
    pub splash_veloc: Vec<f64>,
    pub curb_length: Vec<f64>,
    pub curb_height: Vec<f64>,
    pub curb_throat: Vec<i32>,
    pub slotted_length: Vec<f64>,
    pub slotted_width: Vec<f64>,
    pub curve_index: Vec<i32>,
    pub curve_kind: Vec<i32>,
    pub num_inlets: Vec<i32>,
    pub clog_factor: Vec<f64>,
    pub flow_limit: Vec<f64>,
    pub local_depress: Vec<f64>,
    pub local_width: Vec<f64>,
    pub placement: Vec<i32>,
    pub is_sag: Vec<u8>,
    pub sx: Vec<f64>,
    pub gutter_depression: Vec<f64>,
    pub gutter_width: Vec<f64>,
    pub road_roughness: Vec<f64>,
    pub t_crown: Vec<f64>,
    pub n_sides: Vec<i32>,
    pub slope: Vec<f64>,
    pub flow_factor: Vec<f64>,
    pub q_approach: Vec<f64>,
    pub flow_capture: Vec<f64>,
    pub backflow: Vec<f64>,
    pub backflow_ratio: Vec<f64>,
    pub flow_periods: Vec<i32>,
    pub capture_periods: Vec<i32>,
    pub backflow_periods: Vec<i32>,
    pub peak_flow: Vec<f64>,
    pub peak_flow_capture: Vec<f64>,
    pub avg_flow_capture: Vec<f64>,
    pub bypass_freq: Vec<f64>,
    pub stat_capture_vol: Vec<f64>,
    pub stat_bypass_vol: Vec<f64>,
    pub stat_backflow_vol: Vec<f64>,
    pub stat_peak_flow: Vec<f64>,
}

impl InletSoA {
    pub fn resize(&mut self, n: i32) {
        self.count = n;
        let un = n as usize;
        fn a<T: Clone>(v: &mut Vec<T>, n: usize, d: T) {
            v.clear();
            v.resize(n, d);
        }
        a(&mut self.link_idx, un, -1);
        a(&mut self.host_node, un, -1);
        a(&mut self.node_idx, un, -1);
        a(&mut self.bypass_node, un, -1);
        a(&mut self.up_link, un, -1);
        a(&mut self.up_link2, un, -1);
        a(&mut self.two_inbound, un, 0);

        a(&mut self.inlet_type, un, InletType::Custom as i32);
        a(&mut self.grate_type, un, GrateType::Generic as i32);
        a(&mut self.grate_length, un, 0.0);
        a(&mut self.grate_width, un, 0.0);
        a(&mut self.open_area, un, 0.0);
        a(&mut self.splash_veloc, un, 0.0);
        a(&mut self.curb_length, un, 0.0);
        a(&mut self.curb_height, un, 0.0);
        a(&mut self.curb_throat, un, ThroatAngle::Vertical as i32);
        a(&mut self.slotted_length, un, 0.0);
        a(&mut self.slotted_width, un, 0.0);
        a(&mut self.curve_index, un, -1);
        a(&mut self.curve_kind, un, CurveKind::Unresolved as i32);

        a(&mut self.num_inlets, un, 1);
        a(&mut self.clog_factor, un, 1.0);
        a(&mut self.flow_limit, un, 0.0);
        a(&mut self.local_depress, un, 0.0);
        a(&mut self.local_width, un, 0.0);
        a(&mut self.placement, un, Placement::Automatic as i32);
        a(&mut self.is_sag, un, 0);

        a(&mut self.sx, un, 0.01);
        a(&mut self.gutter_depression, un, 0.0);
        a(&mut self.gutter_width, un, 0.0);
        a(&mut self.road_roughness, un, 0.016);
        a(&mut self.t_crown, un, 100.0);
        a(&mut self.n_sides, un, 1);
        a(&mut self.slope, un, 0.0);
        a(&mut self.flow_factor, un, 0.0);

        a(&mut self.q_approach, un, 0.0);
        a(&mut self.flow_capture, un, 0.0);
        a(&mut self.backflow, un, 0.0);
        a(&mut self.backflow_ratio, un, 0.0);

        a(&mut self.flow_periods, un, 0);
        a(&mut self.capture_periods, un, 0);
        a(&mut self.backflow_periods, un, 0);
        a(&mut self.peak_flow, un, 0.0);
        a(&mut self.peak_flow_capture, un, 0.0);
        a(&mut self.avg_flow_capture, un, 0.0);
        a(&mut self.bypass_freq, un, 0.0);

        a(&mut self.stat_capture_vol, un, 0.0);
        a(&mut self.stat_bypass_vol, un, 0.0);
        a(&mut self.stat_backflow_vol, un, 0.0);
        a(&mut self.stat_peak_flow, un, 0.0);
    }
}

#[derive(Debug, Clone, Default)]
pub struct InletSolver {
    pub soa_: InletSoA,
    pub node_inlet_flow_: Vec<f64>,
    pub node_inflow_: Vec<f64>,
    pub nsides_prev_: i32,
    pub is_capture_node_: Vec<u8>,
}

impl InletSolver {
    pub fn design_of(&self, ii: i32) -> Design {
        let ui = ii as usize;
        let s = &self.soa_;
        Design {
            type_: s.inlet_type[ui],
            grate_type: s.grate_type[ui],
            grate_length: s.grate_length[ui],
            grate_width: s.grate_width[ui],
            frac_open_area: s.open_area[ui],
            splash_veloc: s.splash_veloc[ui],
            curb_length: s.curb_length[ui],
            curb_height: s.curb_height[ui],
            curb_throat: s.curb_throat[ui],
            slotted_length: s.slotted_length[ui],
            slotted_width: s.slotted_width[ui],
        }
    }

    pub fn usage_of(&self, ii: i32) -> UsageParams {
        let ui = ii as usize;
        UsageParams { num_inlets: self.soa_.num_inlets[ui], clog_factor: self.soa_.clog_factor[ui], flow_limit: self.soa_.flow_limit[ui] }
    }

    pub fn soa(&self) -> &InletSoA {
        &self.soa_
    }

    pub fn has_inlets(&self) -> bool {
        self.soa_.count > 0
    }

    /// Resolve designs, geometry, host topology and placement.
    pub fn init(&mut self, ctx: &mut SimulationContext) {
        let usages = &ctx.inlet_usages;
        let n = usages.count();
        self.soa_.resize(n);
        if n == 0 {
            return;
        }

        let nn = ctx.n_nodes();
        let nl = ctx.n_links();
        self.node_inlet_flow_.clear();
        self.node_inlet_flow_.resize(nn as usize, 0.0);
        self.node_inflow_.clear();
        self.node_inflow_.resize(nn as usize, 0.0);
        self.is_capture_node_.clear();
        self.is_capture_node_.resize(nn as usize, 0);

        let us = ucf::get_unit_system(ctx.options.flow_units as i32);
        let inv_len = ucf::UCF_INV[ucf::LENGTH as usize][us as usize];
        let inv_flow = ucf::QCF_INV[ctx.options.flow_units as usize];

        // Legacy Node[].degree = outflow-link count, rebuilt here.
        let mut out_degree = vec![0i32; nn as usize];
        for j in 0..nl {
            let uj = j as usize;
            let mut m = if ctx.links.direction[uj] < 0 { ctx.links.node2[uj] } else { ctx.links.node1[uj] };
            if m < 0 || m >= nn {
                continue;
            }
            if ctx.nodes.type_[m as usize] == NodeType::Outfall {
                m = if ctx.links.direction[uj] < 0 { ctx.links.node1[uj] } else { ctx.links.node2[uj] };
                if m < 0 || m >= nn {
                    continue;
                }
            }
            out_degree[m as usize] += 1;
        }

        for i in 0..n {
            let ui = i as usize;
            let usages = &ctx.inlet_usages;
            let li = usages.link_index[ui];
            let nh = usages.node_host[ui];
            let di = usages.design_index[ui];
            let s = &mut self.soa_;

            s.link_idx[ui] = li;
            s.host_node[ui] = nh;
            s.node_idx[ui] = usages.node_index[ui];
            if s.node_idx[ui] >= 0 && s.node_idx[ui] < nn {
                self.is_capture_node_[s.node_idx[ui] as usize] = 1;
            }
            s.num_inlets[ui] = usages.num_inlets[ui];
            s.clog_factor[ui] = usages.clog_factor[ui];
            s.flow_limit[ui] = usages.flow_limit[ui] * inv_flow;
            s.local_depress[ui] = usages.local_depress[ui] * inv_len;
            s.local_width[ui] = usages.local_width[ui] * inv_len;
            s.placement[ui] = usages.placement[ui];

            // Design (InletStore rows are in user units)
            if di >= 0 && di < ctx.inlets.count() {
                let ud = di as usize;
                let inl = &ctx.inlets;
                let type_str = inl.inlet_type[ud].as_slice();
                let itype = match type_str {
                    b"GRATE" => InletType::Grate,
                    b"CURB" => InletType::Curb,
                    b"COMBO" => InletType::Combo,
                    b"SLOTTED" => InletType::Slotted,
                    b"DROP_GRATE" => InletType::DropGrate,
                    b"DROP_CURB" => InletType::DropCurb,
                    _ => InletType::Custom,
                } as i32;
                s.inlet_type[ui] = itype;

                let gt = match inl.grate_type[ud].as_slice() {
                    b"P_BAR-50" => GrateType::PBar50,
                    b"P_BAR-50x100" => GrateType::PBar50x100,
                    b"P_BAR-30" => GrateType::PBar30,
                    b"CURVED_VANE" => GrateType::CurvedVane,
                    b"TILT_BAR-45" => GrateType::TiltBar45,
                    b"TILT_BAR-30" => GrateType::TiltBar30,
                    b"RETICULINE" => GrateType::Reticuline,
                    _ => GrateType::Generic,
                } as i32;
                s.grate_type[ui] = gt;

                if itype == InletType::Slotted as i32 {
                    s.slotted_length[ui] = inl.length[ud] * inv_len;
                    s.slotted_width[ui] = inl.width[ud] * inv_len;
                } else if itype == InletType::Grate as i32 || itype == InletType::DropGrate as i32 || itype == InletType::Combo as i32 {
                    s.grate_length[ui] = inl.length[ud] * inv_len;
                    s.grate_width[ui] = inl.width[ud] * inv_len;
                }
                s.open_area[ui] = inl.open_area[ud];
                s.splash_veloc[ui] = inl.splash_veloc[ud] * inv_len;

                if itype == InletType::Curb as i32 || itype == InletType::DropCurb as i32 || itype == InletType::Combo as i32 {
                    s.curb_length[ui] = inl.curb_length[ud] * inv_len;
                    s.curb_height[ui] = inl.curb_height[ud] * inv_len;
                    s.curb_throat[ui] = inl.curb_throat[ud];
                }

                s.curve_index[ui] = inl.curve_index[ud];
                s.curve_kind[ui] = inl.curve_kind[ud];
            }

            // Host topology: bypass node, approach conduit(s)
            let mut bypass = -1i32;
            let mut up = -1i32;
            let mut up2 = -1i32;
            if nh >= 0 && nh < nn {
                bypass = nh;
                let mut j = 0;
                while j < nl && up2 < 0 {
                    let uj = j as usize;
                    j += 1;
                    if ctx.links.type_[uj] != LinkType::Conduit {
                        continue;
                    }
                    if ctx.links.node2[uj] != nh {
                        continue;
                    }
                    if up < 0 {
                        up = uj as i32;
                    } else {
                        up2 = uj as i32;
                    }
                }
                if up < 0 {
                    for j in 0..nl {
                        let uj = j as usize;
                        if ctx.links.type_[uj] == LinkType::Conduit && ctx.links.node1[uj] == nh {
                            up = j;
                            break;
                        }
                    }
                }
                if up2 >= 0 {
                    s.two_inbound[ui] = 1;
                    let ca = ctx.link_subtypes.conduit_row(up);
                    let cb = ctx.link_subtypes.conduit_row(up2);
                    let sa = if ca >= 0 { ctx.link_subtypes.conduits.slope[ca as usize] } else { 0.0 };
                    let sb = if cb >= 0 { ctx.link_subtypes.conduits.slope[cb as usize] } else { 0.0 };
                    if sb > sa {
                        std::mem::swap(&mut up, &mut up2);
                    }
                }
            } else if li >= 0 && li < nl {
                bypass = ctx.links.node2[li as usize];
                up = li;
            }
            s.bypass_node[ui] = bypass;
            s.up_link[ui] = up;
            s.up_link2[ui] = up2;

            // Street / channel geometry (inlet.c:1147-1196)
            let si = usages.street_index[ui];
            if si >= 0 && si < ctx.streets.count() {
                let su = si as usize;
                s.sx[ui] = ctx.streets.sx[su] / 100.0; // % -> fraction
                s.gutter_depression[ui] = ctx.streets.gutter_depres[su] * inv_len;
                s.gutter_width[ui] = ctx.streets.gutter_width[su] * inv_len;
                s.road_roughness[ui] = ctx.streets.n_road[su];
                s.n_sides[ui] = ctx.streets.sides[su];
                s.t_crown[ui] = ctx.streets.t_crown[su] * inv_len;
            } else {
                let cr = if up >= 0 { ctx.link_subtypes.conduit_row(up) } else { -1 };
                s.sx[ui] = 0.01;
                s.gutter_depression[ui] = 0.0;
                s.gutter_width[ui] = 0.0;
                s.road_roughness[ui] = if cr >= 0 { ctx.link_subtypes.conduits.roughness[cr as usize] } else { 0.01 };
                s.n_sides[ui] = 1;
                s.t_crown[ui] = 100.0;
            }

            // Longitudinal slope of the approach conduit + Izzard's flow factor.
            {
                let cr = if up >= 0 { ctx.link_subtypes.conduit_row(up) } else { -1 };
                let sl = if cr >= 0 { ctx.link_subtypes.conduits.slope[cr as usize] } else { 0.0 };
                s.slope[ui] = sl;
                // inlet.c:520 -- pow(Sx, 1.67), NOT the exact 5/3 exponent.
                s.flow_factor[ui] = if s.road_roughness[ui] > 0.0 && sl > 0.0 {
                    (0.56 / s.road_roughness[ui]) * sl.powf(0.5) * s.sx[ui].powf(1.67) // LIBM
                } else {
                    0.0
                };
            }

            // Placement
            let mut sag: u8 = 0;
            if s.placement[ui] == Placement::OnSag as i32 {
                sag = 1;
            } else if s.placement[ui] == Placement::OnGrade as i32 {
                sag = 0;
            } else if nh >= 0 {
                let z = ctx.nodes.invert_elev[nh as usize];
                let mut n_toward = 0;
                let mut n_attached = 0;
                for j in 0..nl {
                    let uj = j as usize;
                    if ctx.links.type_[uj] != LinkType::Conduit {
                        continue;
                    }
                    let a = ctx.links.node1[uj];
                    let b = ctx.links.node2[uj];
                    let far;
                    if a == nh {
                        far = b;
                    } else if b == nh {
                        far = a;
                    } else {
                        continue;
                    }
                    n_attached += 1;
                    if far >= 0 && far < nn && ctx.nodes.invert_elev[far as usize] > z + constants::FUDGE {
                        n_toward += 1;
                    }
                }
                sag = if n_attached > 0 && n_toward == n_attached { 1 } else { 0 };
            } else if bypass >= 0 && bypass < nn {
                sag = if out_degree[bypass as usize] > 0 { 0 } else { 1 };
            }
            s.is_sag[ui] = sag;
        }

        // Seed legacy's carried-over Nsides (the FIRST [INLET_USAGE] row).
        let mut ii = self.soa_.count - 1;
        while ii >= 0 {
            self.nsides_prev_ = self.soa_.n_sides[ii as usize];
            ii -= 1;
        }

        self.compute_backflow_ratios();
    }

    /// Legacy getBackflowRatios (inlet.c:1805), including its single aggregate bucket.
    pub fn compute_backflow_ratios(&mut self) {
        let ni = self.soa_.count;
        if ni == 0 {
            return;
        }

        let mut n_links = 0i32;
        let mut n_std_links = 0i32;
        let mut n_custom = 0i32;
        let mut total_area = 0.0f64;

        let mut ii = ni - 1;
        while ii >= 0 {
            let m = self.soa_.node_idx[ii as usize];
            if m >= 0 {
                n_links += 1;
                let area = get_inlet_area(&self.design_of(ii), &self.usage_of(ii));
                if area > 0.0 {
                    n_std_links += 1;
                    total_area += area;
                } else {
                    n_custom += self.soa_.num_inlets[ii as usize];
                }
            }
            ii -= 1;
        }
        if n_links == 0 {
            return;
        }

        let f = n_std_links as f64 / n_links as f64;

        let mut ii = ni - 1;
        while ii >= 0 {
            let ui = ii as usize;
            self.soa_.backflow_ratio[ui] = 0.0;
            if self.soa_.node_idx[ui] >= 0 {
                let area = get_inlet_area(&self.design_of(ii), &self.usage_of(ii));
                if area > 0.0 {
                    if total_area > 0.0 {
                        self.soa_.backflow_ratio[ui] = area / total_area * f;
                    }
                } else if n_custom > 0 {
                    self.soa_.backflow_ratio[ui] = self.soa_.num_inlets[ui] as f64 / n_custom as f64 * (1.0 - f);
                }
            }
            ii -= 1;
        }
    }
}
