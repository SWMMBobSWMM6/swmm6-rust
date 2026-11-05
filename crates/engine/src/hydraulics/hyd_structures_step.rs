// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/HydStructures.cpp and
// HydStructures.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The structures' run-time half (hyd_structures.rs holds the groups and init): tailwaterHead,
// flapGateBlocks, updatePumpTargetSettings, computePumpFlowK, scatterOrificeSurfArea,
// computeOrificeFlowK, computeWeirFlowK, computeOutletFlowK, computeAllFlows,
// computeNonConduitFlowOne, scatterHeldSurfArea and the file-local buildXSP.
//
// Differences in form, never in result:
//   * `double* node_new_surf_area` (nullable upstream) is Option<&mut [f64]>.
//   * The SWMM_TRACE_ORIF trace prints only and is left out.

use crate::core::cmath::{cmax, cmin};
use crate::core::constants;
use crate::core::context::SimulationContext;
use crate::core::options::RoutingModel;
use crate::core::unit_conversion as ucf;
use crate::data::link::{FlowClass, LinkType, XsectShape};
use crate::data::node::NodeType;
use crate::data::table::{table_get_slope, table_interval_lookup, table_lookup_cursor};
use crate::hydraulics::divider;
use crate::hydraulics::hyd_structures::StructureSolver;
use crate::hydraulics::roadway;
use crate::hydraulics::routing_step::steady_xsp;
use crate::hydraulics::xsection as xsect;

/// tailwaterHead: the upstream invert under the tree routers, else the downstream head.
fn tailwater_head(ctx: &SimulationContext, un1: usize, un2: usize) -> f64 {
    let rm = ctx.options.routing_model;
    if rm == RoutingModel::Kinwave || rm == RoutingModel::Steady {
        return ctx.nodes.invert_elev[un1];
    }
    ctx.nodes.depth[un2] + ctx.nodes.invert_elev[un2]
}

/// flapGateBlocks (legacy link_setFlapGate): the link's own gate against reverse flow, or a
/// gated outfall on the inflow end.
pub fn flap_gate_blocks(ctx: &SimulationContext, j: i32, n1: i32, n2: i32, q: f64) -> bool {
    let uj = j as usize;
    if ctx.links.has_flap_gate[uj] != 0 && q * (ctx.links.direction[uj] as f64) < 0.0 {
        return true;
    }
    let mut n = -1;
    if q < 0.0 {
        n = n2;
    }
    if q > 0.0 {
        n = n1;
    }
    if n < 0 {
        return false;
    }
    if ctx.nodes.type_[n as usize] != NodeType::Outfall {
        return false;
    }
    let ofr = ctx.node_subtypes.outfall_row(n);
    ofr >= 0 && ctx.node_subtypes.outfalls.has_flap_gate[ofr as usize] != 0
}

impl StructureSolver {
    /// updatePumpTargetSettings: startup/shutoff depth hysteresis, once per step.
    pub fn update_pump_target_settings(&mut self, ctx: &mut SimulationContext) {
        for k in 0..self.pumps_.count as usize {
            let j = self.pumps_.link_idx[k];
            let uj = j as usize;
            let n1 = ctx.links.node1[uj];
            if n1 < 0 {
                continue;
            }
            let links = &mut ctx.links;
            links.target_setting[uj] = links.setting[uj];
            let depth = ctx.nodes.depth[n1 as usize];
            let y_off = self.pumps_.y_off[k];
            let y_on = self.pumps_.y_on[k];
            if y_off > 0.0 && links.setting[uj] > 0.0 && depth < y_off {
                links.target_setting[uj] = 0.0;
            }
            if y_on > 0.0 && links.setting[uj] == 0.0 && depth > y_on {
                links.target_setting[uj] = 1.0;
            }
        }
    }

    /// computePumpFlowK
    pub fn compute_pump_flow_k(&mut self, ctx: &mut SimulationContext, dt: f64, node_new_surf_area: Option<&[f64]>, k: i32) {
        let fu = ctx.options.flow_units as i32;
        let ucf_flow = ucf::QCF[fu as usize];
        let unit_sys = ucf::get_unit_system(fu);
        let ucf_len = ucf::UCF[ucf::LENGTH as usize][unit_sys as usize];
        let ucf_vol = ucf::UCF[ucf::VOLUME as usize][unit_sys as usize];

        let uk = k as usize;
        let j = self.pumps_.link_idx[uk];
        let uj = j as usize;
        let n1 = ctx.links.node1[uj];
        let n2 = ctx.links.node2[uj];
        if n1 < 0 || n2 < 0 {
            ctx.links.flow[uj] = 0.0;
            return;
        }
        let un1 = n1 as usize;
        let un2 = n2 as usize;

        let depth = ctx.nodes.depth[un1];
        let head = (ctx.nodes.depth[un2] + ctx.nodes.invert_elev[un2]) - (ctx.nodes.depth[un1] + ctx.nodes.invert_elev[un1]);

        ctx.links.setting[uj] = ctx.links.target_setting[uj];
        ctx.links.flow_class[uj] = FlowClass::Dry;

        if ctx.links.setting[uj] == 0.0 {
            ctx.links.flow[uj] = 0.0;
            return;
        }

        let mut q = 0.0;
        let ct = self.pumps_.curve_type[uk];
        let ci = self.pumps_.curve_idx[uk];
        let uci = ci as usize;
        let have_curve = ci >= 0 && uci < ctx.tables.tables.len();

        match ct {
            1 => {
                if have_curve {
                    let vol = ctx.nodes.volume[un1] * ucf_vol;
                    q = table_interval_lookup(&ctx.tables.tables[uci], vol);
                }
            }
            2 => {
                if have_curve {
                    q = table_interval_lookup(&ctx.tables.tables[uci], depth * ucf_len);
                }
            }
            3 | 5 => {
                let s = if ct == 5 { ctx.links.setting[uj] } else { 1.0 };
                let h = if s > 0.0 { cmax(head / (s * s), 0.0) } else { 0.0 };
                if have_curve {
                    let curve = &mut ctx.tables.tables[uci];
                    q = table_lookup_cursor(curve, h * ucf_len);
                    if s > 0.0 {
                        let slope = table_get_slope(curve, h * ucf_len);
                        ctx.links.dqdh[uj] = -slope * ucf_len / ucf_flow / s;
                    }
                }
            }
            4 => {
                if have_curve {
                    q = table_lookup_cursor(&mut ctx.tables.tables[uci], depth * ucf_len);
                    const DH: f64 = 0.001;
                    let q1 = table_lookup_cursor(&mut ctx.tables.tables[uci], (depth + DH) * ucf_len);
                    ctx.links.dqdh[uj] = (q1 - q) / (DH * ucf_flow);
                }
            }
            6 => {
                q = ctx.nodes.inflow[un1] + ctx.nodes.overflow[un1];
            }
            _ => {}
        }

        // Off-curve classification (report only)
        if have_curve && (1..=5).contains(&ct) {
            let x = &ctx.tables.tables[uci].x;
            if !x.is_empty() {
                let mut argument = depth * ucf_len;
                if ct == 1 {
                    argument = ctx.nodes.volume[un1] * ucf_vol;
                }
                if ct == 3 || ct == 5 {
                    let speed = if ct == 5 { ctx.links.setting[uj] } else { 1.0 };
                    argument = cmax(head / (speed * speed), 0.0) * ucf_len;
                }
                if argument < x[0] || argument > x[x.len() - 1] {
                    ctx.links.flow_class[uj] = if ct == 4 && argument < x[0] { FlowClass::DnDry } else { FlowClass::UpDry };
                }
            }
        }

        q = cmax(q, 0.0);
        if ct != 6 {
            q /= ucf_flow;
        }
        q *= ctx.links.setting[uj];

        if q > 0.0 {
            let nodes = &ctx.nodes;
            if nodes.type_[un1] == NodeType::Storage || nodes.storage_area[un1] > 0.0 || self.pumps_.curve_type[uk] == 1 {
                let full_vol = if ctx.options.routing_model == RoutingModel::Dynwave { nodes.rpt_full_volume[un1] } else { nodes.full_volume[un1] };
                if full_vol > 0.0 {
                    let max_q = nodes.inflow[un1] + nodes.old_volume[un1] / dt;
                    if q > max_q {
                        q = max_q;
                    }
                }
                if q < 0.0 {
                    q = 0.0;
                }
            } else if ct == 2 || ct == 3 || ct == 4 {
                let net_inflow = nodes.inflow[un1] - nodes.outflow[un1] - q;
                let net_vol = 0.5 * (nodes.old_net_inflow[un1] + net_inflow) * dt;
                let surf = match node_new_surf_area {
                    Some(a) => a[un1],
                    None => constants::MIN_SURFAREA,
                };
                let y_new = nodes.old_depth[un1] + net_vol / surf;
                if y_new <= 0.0 {
                    q = cmax(nodes.inflow[un1], 0.0);
                }
            }
        }

        ctx.links.flow[uj] = q;
    }

    /// scatterOrificeSurfArea (legacy findNonConduitSurfArea for an orifice).
    fn scatter_orifice_surf_area(&self, ctx: &SimulationContext, node_new_surf_area: Option<&mut [f64]>, uk: usize, uj: usize, un1: usize, un2: usize) {
        let a = match node_new_surf_area {
            Some(a) => a,
            None => return,
        };
        let mut sa1 = self.orifices_.surf_area[uk] * 0.5;
        let mut sa2 = sa1;
        let fc = ctx.links.flow_class[uj];
        if fc == FlowClass::UpCritical || ctx.nodes.type_[un1] == NodeType::Storage {
            sa1 = 0.0;
        }
        if fc == FlowClass::DnCritical || ctx.nodes.type_[un2] == NodeType::Storage {
            sa2 = 0.0;
        }
        a[un1] += sa1;
        a[un2] += sa2;
    }

    /// computeOrificeFlowK
    pub fn compute_orifice_flow_k(&mut self, ctx: &mut SimulationContext, mut node_new_surf_area: Option<&mut [f64]>, k: i32) {
        use constants::GRAVITY;
        const FUDGE_ORI: f64 = 0.0001;

        let uk = k as usize;
        let j = self.orifices_.link_idx[uk];
        let uj = j as usize;
        let n1 = ctx.links.node1[uj];
        let n2 = ctx.links.node2[uj];
        if n1 < 0 || n2 < 0 {
            ctx.links.flow[uj] = 0.0;
            return;
        }
        let un1 = n1 as usize;
        let un2 = n2 as usize;

        let setting = ctx.links.setting[uj];
        let y_full = ctx.links.xsect_y_full[uj];
        let orow = ctx.link_subtypes.orifice_row(j) as usize;
        let cd_val = ctx.link_subtypes.orifices.cd[orow];
        let h_open = setting * y_full;
        if setting == 0.0 {
            ctx.links.flow[uj] = 0.0;
            self.scatter_orifice_surf_area(ctx, node_new_surf_area, uk, uj, un1, un2);
            return;
        }

        let (a_eff, xs_w_at, xs_a_at) = {
            let xs = steady_xsp(ctx, uj);
            (xsect::get_a_of_y(&xs, h_open), (), ())
        };
        let _ = (xs_w_at, xs_a_at);
        let f_area = a_eff * (2.0 * GRAVITY).sqrt();
        let c_orif = cd_val * f_area;

        let is_side = ctx.link_subtypes.orifices.orifice_type[orow] > 0.5;
        let mut h_crit;
        let c_weir;
        if !is_side {
            let a_over_l = if ctx.links.xsect_shape[uj] == XsectShape::Circular {
                h_open / 4.0
            } else {
                let w = ctx.links.xsect_w_max[uj];
                if w > 0.0 && h_open > 0.0 {
                    (h_open * w) / (2.0 * (h_open + w))
                } else {
                    h_open / 4.0
                }
            };
            h_crit = (cd_val / 0.414) * a_over_l;
            c_weir = cd_val * h_crit.sqrt() * f_area;
        } else {
            h_crit = h_open;
            c_weir = cd_val * (h_open / 2.0).sqrt() * f_area;
        }
        h_crit = cmax(h_crit, FUDGE_ORI);

        let mut h1 = ctx.nodes.depth[un1] + ctx.nodes.invert_elev[un1];
        let mut h2 = tailwater_head(ctx, un1, un2);
        let dir = if h1 >= h2 { 1.0 } else { -1.0 };

        let mut y1 = ctx.nodes.depth[un1];
        if dir < 0.0 {
            std::mem::swap(&mut h1, &mut h2);
            y1 = ctx.nodes.depth[un2];
        }

        let hcrest = ctx.nodes.invert_elev[un1] + ctx.links.offset1[uj];

        let mut head;
        let mut f;
        if !is_side {
            if h1 < hcrest {
                head = 0.0;
            } else if h2 > hcrest {
                head = h1 - h2;
            } else {
                head = h1 - hcrest;
            }
            f = cmin(head / h_crit, 1.0);
        } else {
            let hcrown = hcrest + y_full * setting;
            let hmidpt = (hcrest + hcrown) / 2.0;
            if h1 < hcrown && hcrown > hcrest {
                f = (h1 - hcrest) / (hcrown - hcrest);
            } else {
                f = 1.0;
            }
            if f < 1.0 {
                head = h1 - hcrest;
            } else if h2 < hmidpt {
                head = h1 - hmidpt;
            } else {
                head = h1 - h2;
            }
        }

        if head <= FUDGE_ORI || y1 <= FUDGE_ORI {
            ctx.links.flow[uj] = 0.0;
            ctx.links.depth[uj] = 0.0;
            ctx.links.dqdh[uj] = 0.0;
            ctx.links.flow_class[uj] = FlowClass::Dry;
            self.orifices_.surf_area[uk] = FUDGE_ORI * self.orifices_.length_eff[uk];
            self.scatter_orifice_surf_area(ctx, node_new_surf_area, uk, uj, un1, un2);
            return;
        }

        if flap_gate_blocks(ctx, j, n1, n2, dir) {
            ctx.links.flow[uj] = 0.0;
            ctx.links.depth[uj] = 0.0;
            ctx.links.dqdh[uj] = 0.0;
            ctx.links.flow_class[uj] = FlowClass::Dry;
            self.orifices_.surf_area[uk] = FUDGE_ORI * self.orifices_.length_eff[uk];
            self.scatter_orifice_surf_area(ctx, node_new_surf_area, uk, uj, un1, un2);
            return;
        }

        if hcrest > h2 {
            ctx.links.flow_class[uj] = if dir > 0.0 { FlowClass::DnCritical } else { FlowClass::UpCritical };
        } else {
            ctx.links.flow_class[uj] = FlowClass::Subcritical;
        }

        let y_link = y_full * setting;
        if is_side {
            let link_depth = y_link * cmax(f, 0.0);
            ctx.links.depth[uj] = link_depth;
            let width_at_depth = xsect::get_w_of_y(&steady_xsp(ctx, uj), link_depth);
            self.orifices_.surf_area[uk] = width_at_depth * self.orifices_.length_eff[uk];
        } else {
            ctx.links.depth[uj] = y_link;
            self.orifices_.surf_area[uk] = xsect::get_a_of_y(&steady_xsp(ctx, uj), y_link);
        }

        let mut q;
        let mut dqdh = 0.0;
        if f <= 0.0 {
            q = 0.0;
        } else if f < 1.0 {
            q = c_weir * f.powf(1.5); // LIBM pow
            dqdh = 1.5 * q / (f * h_crit);
        } else {
            q = c_orif * head.sqrt();
            dqdh = q / (2.0 * head);
        }

        if ctx.links.has_flap_gate[uj] != 0 && q > 0.0 && a_eff > FUDGE_ORI {
            let veloc = q / a_eff;
            let h_loss = (4.0 / GRAVITY) * veloc * veloc * (-1.15 * veloc / head.sqrt()).exp(); // LIBM exp
            if f < 1.0 {
                f -= h_loss / h_crit;
                f = cmax(f, 0.0);
            } else {
                head -= h_loss;
                head = cmax(head, 0.0);
            }
            if f <= 0.0 || head <= 0.0 {
                q = 0.0;
                dqdh = 0.0;
            } else if f < 1.0 {
                q = c_weir * f.powf(1.5); // LIBM pow
                dqdh = 1.5 * q / (f * h_crit);
            } else {
                q = c_orif * head.sqrt();
                dqdh = q / (2.0 * head);
            }
        }

        if f < 1.0 && h2 > hcrest {
            let ratio = (h2 - hcrest) / (h1 - hcrest);
            q *= (1.0 - ratio.powf(1.5)).powf(0.385); // LIBM pow
        }

        q = cmax(q, 0.0);
        ctx.links.flow[uj] = q * dir;
        ctx.links.dqdh[uj] = dqdh;

        self.scatter_orifice_surf_area(ctx, node_new_surf_area.as_deref_mut(), uk, uj, un1, un2);
    }

    /// computeWeirFlowK
    pub fn compute_weir_flow_k(&mut self, ctx: &mut SimulationContext, _node_new_surf_area: Option<&mut [f64]>, k: i32) {
        const FUDGE_W: f64 = 0.0001;
        let uk = k as usize;
        let j = self.weirs_.link_idx[uk];
        let uj = j as usize;
        let n1 = ctx.links.node1[uj];
        let n2 = ctx.links.node2[uj];
        if n1 < 0 || n2 < 0 {
            ctx.links.flow[uj] = 0.0;
            ctx.links.depth[uj] = 0.0;
            self.weirs_.surf_area[uk] = 0.0;
            return;
        }
        let un1 = n1 as usize;
        let un2 = n2 as usize;

        if ctx.links.setting[uj] == 0.0 {
            ctx.links.flow[uj] = 0.0;
            return;
        }

        let mut hgl1 = ctx.nodes.depth[un1] + ctx.nodes.invert_elev[un1];
        let mut hgl2 = tailwater_head(ctx, un1, un2);
        let mut dir = if hgl1 >= hgl2 { 1.0 } else { -1.0 };

        let wrow = ctx.link_subtypes.weir_row(j) as usize;

        if self.weirs_.weir_type[uk] == 4 {
            dir = if hgl1 > hgl2 { 1.0 } else { -1.0 };
            if dir < 0.0 {
                std::mem::swap(&mut hgl1, &mut hgl2);
            }
            let h_road = ctx.nodes.invert_elev[un1] + ctx.link_subtypes.weirs.crest_height[wrow];
            let si = ucf::get_unit_system(ctx.options.flow_units as i32) == 1;
            let mut dqdh_r = 0.0;
            let mut depth_r = 0.0;
            ctx.links.flow[uj] = roadway::get_inflow(dir, h_road, hgl1, hgl2, self.weirs_.c_disch1[uk], si, self.weirs_.road_width[uk], self.weirs_.road_surface[uk], ctx.links.xsect_w_max[uj], &mut dqdh_r, &mut depth_r);
            ctx.links.dqdh[uj] = dqdh_r;
            ctx.links.depth[uj] = depth_r;
            self.weirs_.surf_area[uk] = 0.0;
            return;
        }

        if flap_gate_blocks(ctx, j, n1, n2, dir) {
            ctx.links.flow[uj] = 0.0;
            ctx.links.depth[uj] = 0.0;
            ctx.links.dqdh[uj] = 0.0;
            self.weirs_.surf_area[uk] = 0.0;
            return;
        }

        if dir < 0.0 {
            std::mem::swap(&mut hgl1, &mut hgl2);
        }

        let y_full = ctx.links.xsect_y_full[uj];
        let setting = ctx.links.setting[uj];
        let mut hcrest = ctx.nodes.invert_elev[un1] + ctx.link_subtypes.weirs.crest_height[wrow];
        let hcrown = hcrest + y_full;
        hcrest += (1.0 - setting) * y_full;
        let mut head = hgl1 - hcrest;
        if head <= FUDGE_W || hcrest >= hcrown {
            ctx.links.flow[uj] = 0.0;
            ctx.links.depth[uj] = 0.0;
            ctx.links.dqdh[uj] = 0.0;
            self.weirs_.surf_area[uk] = 0.0;
            return;
        }

        let cd_written = self.weirs_.c_disch1[uk];
        let cd2 = self.weirs_.c_disch2[uk];
        let length = ctx.links.xsect_w_max[uj];
        let wt = self.weirs_.weir_type[uk];
        let flap = ctx.links.has_flap_gate[uj] != 0;
        let mut q;

        const M3_PER_FT3: f64 = 0.028317;
        let us = ucf::get_unit_system(ctx.options.flow_units as i32);
        let ucf_l = ucf::UCF[ucf::LENGTH as usize][us as usize];
        let end_con = self.weirs_.end_con[uk];
        let slope = self.weirs_.slope[uk];
        let cd_curve = self.weirs_.cd_curve[uk];
        let y_bot = ctx.links.xsect_y_bot[uj];
        let weir_flow = |ctx: &mut SimulationContext, head_ft: f64, dir_f: f64, q1_out: &mut f64, q2_out: &mut f64| {
            *q1_out = 0.0;
            *q2_out = 0.0;
            if head_ft <= 0.0 {
                return;
            }
            let mut length_d = length * ucf_l;
            let h = head_ft * ucf_l;
            let mut cd = cd_written;
            if cd_curve >= 0 {
                let uci = cd_curve as usize;
                if uci < ctx.tables.tables.len() {
                    cd = table_lookup_cursor(&mut ctx.tables.tables[uci], h);
                }
            }
            let mut w_type = wt;
            if w_type == 2 && setting < 1.0 {
                w_type = 3;
            }
            match w_type {
                0 => {
                    length_d -= 0.1 * end_con * h;
                    length_d = cmax(length_d, 0.0);
                    *q1_out = cd * length_d * h.powf(1.5); // LIBM pow
                }
                1 => {
                    length_d -= 0.1 * end_con * h;
                    length_d = cmax(length_d, 0.0);
                    if dir_f < 0.0 {
                        *q1_out = cd * length_d * h.powf(1.5); // LIBM pow
                    } else {
                        *q1_out = cd * length_d.powf(0.83) * h.powf(1.67); // LIBM pow
                    }
                }
                2 => {
                    *q1_out = cd * slope * h.powf(2.5); // LIBM pow
                }
                3 => {
                    let y = (1.0 - setting) * y_full;
                    let w = if wt == 2 { 2.0 * slope * y } else { y_bot + 2.0 * slope * y };
                    *q1_out = cd * (w * ucf_l) * h.powf(1.5); // LIBM pow
                    *q2_out = cd2 * slope * h.powf(2.5); // LIBM pow
                }
                _ => {}
            }
            if us == 1 {
                *q1_out /= M3_PER_FT3;
                *q2_out /= M3_PER_FT3;
            }
        };

        let weir_open_area = |y: f64| -> f64 {
            let z = (1.0 - setting) * y_full;
            let zy = cmin(z + y, y_full);
            let area_of = |d: f64| -> f64 {
                match wt {
                    2 => slope * d * d,
                    3 => (y_bot + slope * d) * d,
                    _ => length * d,
                }
            };
            area_of(zy) - area_of(z)
        };

        {
            let h1_min = cmin(hgl1, hcrown);
            let y_sa = y_full - (hcrown - h1_min);
            let mut width_at_y = 0.0;
            if y_sa > 0.0 {
                width_at_y = match wt {
                    0 | 1 => length,
                    2 => 2.0 * slope * y_sa,
                    3 => y_bot + 2.0 * slope * y_sa,
                    _ => length,
                };
            }
            self.weirs_.surf_area[uk] = width_at_y * self.weirs_.length_eff[uk];
        }

        if hgl1 >= hcrown && self.weirs_.can_surcharge[uk] != 0 {
            let h_full = setting * y_full;
            let mut q1f = 0.0;
            let mut q2f = 0.0;
            weir_flow(ctx, h_full, 1.0, &mut q1f, &mut q2f);
            let h_half = h_full / 2.0;
            let c_surcharge = if h_half > 0.0 { (q1f + q2f) / h_half.sqrt() } else { 0.0 };

            let y_mid = (hcrest + hcrown) / 2.0;
            let mut h_orif = if hgl2 < y_mid { hgl1 - y_mid } else { hgl1 - hgl2 };
            let y_open = hcrown - hcrest;
            q = c_surcharge * h_orif.sqrt();
            if flap {
                let a = weir_open_area(y_open);
                if a > 0.0 {
                    let v = q / a;
                    let hloss = (4.0 / constants::GRAVITY) * v * v * (-1.15 * v / y_open.sqrt()).exp(); // LIBM exp
                    h_orif -= hloss;
                    h_orif = cmax(h_orif, 0.0);
                    q = c_surcharge * h_orif.sqrt();
                }
            }
            ctx.links.depth[uj] = y_open;
            ctx.links.dqdh[uj] = if h_orif > 0.0 { q / (2.0 * h_orif) } else { 0.0 };
        } else {
            if hgl1 >= hcrown {
                head = hcrown - hcrest;
            }
            let mut q1 = 0.0;
            let mut q2 = 0.0;
            weir_flow(ctx, head, dir, &mut q1, &mut q2);

            if flap {
                let area = weir_open_area(head);
                if area > 1.0e-6 {
                    let veloc = (q1 + q2) / area;
                    let h_loss = (4.0 / constants::GRAVITY) * veloc * veloc * (-1.15 * veloc / head.sqrt()).exp(); // LIBM exp
                    head -= h_loss;
                    if head < 0.0 {
                        head = 0.0;
                    }
                    weir_flow(ctx, head, dir, &mut q1, &mut q2);
                }
            }

            if head.abs() >= FUDGE_W {
                let q1h = (q1 / head).abs();
                let q2h = (q2 / head).abs();
                ctx.links.dqdh[uj] = match wt {
                    0 => 1.5 * q1h,
                    1 => (if dir < 0.0 { 1.5 } else { 1.67 }) * q1h,
                    2 => {
                        if q2h == 0.0 {
                            2.5 * q1h
                        } else {
                            1.5 * q1h + 2.5 * q2h
                        }
                    }
                    _ => 1.5 * q1h + 2.5 * q2h,
                };
            } else {
                ctx.links.dqdh[uj] = 0.0;
            }

            if hgl2 > hcrest {
                const WEIR_POWER: [f64; 4] = [1.5, 5.0 / 3.0, 2.5, 1.5];
                let ratio = (hgl2 - hcrest) / (hgl1 - hcrest);
                let wp = if (0..4).contains(&wt) { WEIR_POWER[wt as usize] } else { 1.5 };
                q1 *= (1.0 - ratio.powf(wp)).powf(0.385); // LIBM pow
                if q2 > 0.0 {
                    q2 *= (1.0 - ratio.powf(2.5)).powf(0.385); // LIBM pow
                }
            }

            q = q1 + q2;
            ctx.links.depth[uj] = cmin(hgl1 - hcrest, y_full);
        }

        if q < 0.0 {
            q = 0.0;
        }
        ctx.links.flow[uj] = q * dir;
    }

    /// computeOutletFlowK
    pub fn compute_outlet_flow_k(&mut self, ctx: &mut SimulationContext, k: i32) {
        let fu = ctx.options.flow_units as i32;
        let unit_sys = ucf::get_unit_system(fu);
        let ucf_len = ucf::UCF[ucf::LENGTH as usize][unit_sys as usize];
        let ucf_flow = ucf::QCF[fu as usize];
        const FUDGE_OUT: f64 = 0.0001;

        let uk = k as usize;
        let j = self.outlets_.link_idx[uk];
        let uj = j as usize;
        let n1 = ctx.links.node1[uj];
        let n2 = ctx.links.node2[uj];
        if n1 < 0 || n2 < 0 {
            ctx.links.flow[uj] = 0.0;
            ctx.links.depth[uj] = 0.0;
            return;
        }
        let un1 = n1 as usize;
        let un2 = n2 as usize;

        if ctx.links.setting[uj] == 0.0 {
            ctx.links.flow[uj] = 0.0;
            return;
        }

        let mut h1 = ctx.nodes.depth[un1] + ctx.nodes.invert_elev[un1];
        let mut h2 = tailwater_head(ctx, un1, un2);
        let dir: i32 = if h1 >= h2 { 1 } else { -1 };

        let mut y1 = ctx.nodes.depth[un1];
        if dir < 0 {
            std::mem::swap(&mut h1, &mut h2);
            y1 = ctx.nodes.depth[un2];
        }

        let olr_c = ctx.link_subtypes.outlet_row(j);
        let hcrest = ctx.nodes.invert_elev[un1] + if olr_c >= 0 { ctx.link_subtypes.outlets.crest_height[olr_c as usize] } else { 0.0 };

        let outlet_type = ctx.link_subtypes.outlets.outlet_type[ctx.link_subtypes.outlet_row(j) as usize] as i32;
        let depth_based = outlet_type == 1 || outlet_type == 3;

        let head = if depth_based { h1 - hcrest } else { h1 - cmax(h2, hcrest) };
        let blocked_by_flap = flap_gate_blocks(ctx, j, n1, n2, dir as f64);

        if head <= FUDGE_OUT || y1 <= FUDGE_OUT || blocked_by_flap {
            ctx.links.flow[uj] = 0.0;
            ctx.links.depth[uj] = 0.0;
            ctx.links.flow_class[uj] = FlowClass::Dry;
            ctx.links.dqdh[uj] = 0.0;
            return;
        }

        let lookup_val = head * ucf_len;

        let mut q = 0.0;
        let ci = self.outlets_.curve_idx[uk];
        if ci >= 0 && (ci as usize) < ctx.tables.tables.len() {
            q = table_lookup_cursor(&mut ctx.tables.tables[ci as usize], lookup_val);
            q /= ucf_flow;
        } else if lookup_val > 0.0 {
            q = self.outlets_.q_coeff[uk] * lookup_val.powf(self.outlets_.q_expon[uk]) / ucf_flow; // LIBM pow
        }

        q *= ctx.links.setting[uj];

        if q < 0.0 {
            q = 0.0;
        }
        ctx.links.flow[uj] = q * (dir as f64);
        ctx.links.depth[uj] = head;
        ctx.links.flow_class[uj] = FlowClass::Subcritical;
        ctx.links.dqdh[uj] = 0.0;
    }

    /// computeNonConduitFlowOne: one link's flow, in legacy findLinkFlows' per-link order.
    pub fn compute_non_conduit_flow_one(&mut self, ctx: &mut SimulationContext, dt: f64, node_new_surf_area: Option<&mut [f64]>, link_idx: i32) {
        let uj = link_idx as usize;
        if ctx.links.type_[uj] == LinkType::Conduit {
            ctx.links.dqdh[uj] = 0.0;
            let mut q = 0.0;
            if ctx.links.setting[uj] != 0.0 {
                let n1 = ctx.links.node1[uj];
                if n1 >= 0 {
                    let un1 = n1 as usize;
                    if ctx.nodes.type_[un1] == NodeType::Divider {
                        q = divider::get_outflow(ctx, n1, link_idx);
                    } else {
                        q = ctx.nodes.inflow[un1] + ctx.nodes.overflow[un1];
                    }
                }
                let qlim = ctx.links.q_limit[uj];
                if qlim > 0.0 && q > qlim {
                    q = qlim;
                }
            }
            ctx.links.flow[uj] = q;
            return;
        }

        let k = if uj < self.nc_group_k_.len() { self.nc_group_k_[uj] } else { -1 };
        if k < 0 {
            return;
        }

        ctx.links.dqdh[uj] = 0.0;

        match ctx.links.type_[uj] {
            LinkType::Pump => self.compute_pump_flow_k(ctx, dt, node_new_surf_area.as_deref(), k),
            LinkType::Orifice => self.compute_orifice_flow_k(ctx, node_new_surf_area, k),
            LinkType::Weir => self.compute_weir_flow_k(ctx, node_new_surf_area, k),
            LinkType::Outlet => self.compute_outlet_flow_k(ctx, k),
            _ => {}
        }
    }

    /// scatterHeldSurfArea: a bypassed orifice still scatters its held surface area.
    pub fn scatter_held_surf_area(&self, ctx: &SimulationContext, node_new_surf_area: Option<&mut [f64]>, link_idx: i32) {
        let uj = link_idx as usize;
        if ctx.links.type_[uj] != LinkType::Orifice {
            return;
        }
        let k = if uj < self.nc_group_k_.len() { self.nc_group_k_[uj] } else { -1 };
        if k < 0 {
            return;
        }
        let n1 = ctx.links.node1[uj];
        let n2 = ctx.links.node2[uj];
        if n1 < 0 || n2 < 0 {
            return;
        }
        self.scatter_orifice_surf_area(ctx, node_new_surf_area, k as usize, uj, n1 as usize, n2 as usize);
    }
}
