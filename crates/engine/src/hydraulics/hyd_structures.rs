// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/HydStructures.hpp and
// HydStructures.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Non-conduit hydraulic structures (pumps, orifices, weirs, outlets): the per-type SoA
// groups and StructureSolver::init.
//
// NOT YET TRANSLATED (run-time step functions, with the routing step): flapGateBlocks,
// tailwaterHead, computeAllFlows, updatePumpTargetSettings, computeNonConduitFlowOne,
// scatterHeldSurfArea, computePumpFlows/computePumpFlowK, computeOrificeFlowK,
// computeWeirFlowK, computeOutletFlowK, scatterOrificeSurfArea, computeOrificeFlows,
// computeWeirFlows, computeOutletFlows, and the file-local buildXSP.
//
// Differences in form, never in result:
//   * Members keep upstream's names and are `pub` (upstream's private ones too) so the
//     gate's state dumper can read them.

use crate::core::cmath::cmax;
use crate::core::constants;
use crate::core::context::SimulationContext;
use crate::data::link::{LinkType, XsectShape};

pub const GRAVITY: f64 = 32.2;
pub const FUDGE: f64 = 0.0001;

#[derive(Debug, Clone, Default)]
pub struct PumpGroup {
    pub count: i32,
    pub link_idx: Vec<i32>,
    pub curve_idx: Vec<i32>,
    pub curve_type: Vec<i32>,
    pub speed: Vec<f64>,
    pub y_on: Vec<f64>,
    pub y_off: Vec<f64>,
}

impl PumpGroup {
    pub fn resize(&mut self, n: i32) {
        self.count = n;
        let u = n as usize;
        self.link_idx.resize(u, 0);
        self.curve_idx.resize(u, -1);
        self.curve_type.resize(u, 0);
        self.speed.resize(u, 1.0);
        self.y_on.resize(u, 0.0);
        self.y_off.resize(u, 0.0);
    }
}

#[derive(Debug, Clone, Default)]
pub struct OrificeGroup {
    pub count: i32,
    pub link_idx: Vec<i32>,
    pub shape: Vec<i32>,
    pub surf_area: Vec<f64>,
    pub length_eff: Vec<f64>,
}

impl OrificeGroup {
    pub fn resize(&mut self, n: i32) {
        self.count = n;
        let u = n as usize;
        self.link_idx.resize(u, 0);
        self.shape.resize(u, 0);
        self.surf_area.resize(u, 0.0);
        self.length_eff.resize(u, 0.0);
    }
}

#[derive(Debug, Clone, Default)]
pub struct WeirGroup {
    pub count: i32,
    pub link_idx: Vec<i32>,
    pub weir_type: Vec<i32>,
    pub c_disch1: Vec<f64>,
    pub c_disch2: Vec<f64>,
    pub end_con: Vec<f64>,
    pub slope: Vec<f64>,
    pub cd_curve: Vec<i32>,
    pub road_width: Vec<f64>,
    pub road_surface: Vec<i32>,
    pub surf_area: Vec<f64>,
    pub length_eff: Vec<f64>,
    pub can_surcharge: Vec<u8>,
}

impl WeirGroup {
    pub fn resize(&mut self, n: i32) {
        self.count = n;
        let u = n as usize;
        self.link_idx.resize(u, 0);
        self.weir_type.resize(u, 0);
        self.c_disch1.resize(u, 0.0);
        self.c_disch2.resize(u, 0.0);
        self.end_con.resize(u, 0.0);
        self.slope.resize(u, 0.0);
        self.cd_curve.resize(u, -1);
        self.road_width.resize(u, 0.0);
        self.road_surface.resize(u, 0);
        self.surf_area.resize(u, 0.0);
        self.length_eff.resize(u, 0.0);
        self.can_surcharge.resize(u, 1);
    }
}

#[derive(Debug, Clone, Default)]
pub struct OutletGroup {
    pub count: i32,
    pub link_idx: Vec<i32>,
    pub curve_idx: Vec<i32>,
    pub q_coeff: Vec<f64>,
    pub q_expon: Vec<f64>,
}

impl OutletGroup {
    pub fn resize(&mut self, n: i32) {
        self.count = n;
        let u = n as usize;
        self.link_idx.resize(u, 0);
        self.curve_idx.resize(u, -1);
        self.q_coeff.resize(u, 0.0);
        self.q_expon.resize(u, 1.0);
    }
}

/// Batch solver for the non-conduit links.
#[derive(Debug, Clone, Default)]
pub struct StructureSolver {
    pub pumps_: PumpGroup,
    pub orifices_: OrificeGroup,
    pub weirs_: WeirGroup,
    pub outlets_: OutletGroup,
    /// All non-conduit link indices (DUMMY conduits included)
    pub nc_indices_: Vec<i32>,
    /// link index -> row in its type group (-1 = none)
    pub nc_group_k_: Vec<i32>,
}

impl StructureSolver {
    /// Group the links by type and populate the per-type SoA blocks.
    pub fn init(&mut self, ctx: &mut SimulationContext) {
        let mut n_pumps = 0i32;
        let mut n_orifices = 0i32;
        let mut n_weirs = 0i32;
        let mut n_outlets = 0i32;

        for j in 0..ctx.n_links() {
            match ctx.links.type_[j as usize] {
                LinkType::Pump => n_pumps += 1,
                LinkType::Orifice => n_orifices += 1,
                LinkType::Weir => n_weirs += 1,
                LinkType::Outlet => n_outlets += 1,
                _ => {}
            }
        }

        self.pumps_.resize(n_pumps);
        self.orifices_.resize(n_orifices);
        self.weirs_.resize(n_weirs);
        self.outlets_.resize(n_outlets);

        let mut ip = 0usize;
        let mut io = 0usize;
        let mut iw = 0usize;
        let mut ix = 0usize;
        for j in 0..ctx.n_links() {
            let uj = j as usize;
            match ctx.links.type_[uj] {
                LinkType::Pump => {
                    let uk = ip;
                    let pr = ctx.link_subtypes.pump_row(j) as usize;
                    let p = &ctx.link_subtypes.pumps;
                    self.pumps_.link_idx[uk] = j;
                    self.pumps_.curve_idx[uk] = p.curve[pr];
                    self.pumps_.speed[uk] = ctx.links.setting[uj];
                    self.pumps_.y_on[uk] = p.startup[pr];
                    self.pumps_.y_off[uk] = p.shutoff[pr];
                    let ci = p.curve[pr];
                    if ci >= 0 && ci < ctx.tables.tables.len() as i32 {
                        let tt = ctx.tables.tables[ci as usize].type_ as i32;
                        // CURVE_PUMP1..PUMP5 = 7..11 map to curve_type 1..5
                        if (7..=11).contains(&tt) {
                            self.pumps_.curve_type[uk] = tt - 6;
                        } else {
                            self.pumps_.curve_type[uk] = 6; // curve present but not a pump type
                        }
                    } else {
                        // No curve ('*' in [PUMPS]) => IDEAL_PUMP.
                        self.pumps_.curve_type[uk] = 6;
                    }
                    ctx.link_subtypes.pumps.curve_type[pr] = self.pumps_.curve_type[uk];
                    ip += 1;
                }
                LinkType::Orifice => {
                    let uk = io;
                    self.orifices_.link_idx[uk] = j;
                    let y_full = ctx.links.xsect_y_full[uj];
                    {
                        let route_step = ctx.options.link_validate_routing_step();
                        let l = 2.0 * route_step * (constants::GRAVITY * y_full).sqrt();
                        self.orifices_.length_eff[uk] = cmax(200.0, l);
                    }
                    io += 1;
                }
                LinkType::Weir => {
                    let uk = iw;
                    let w = &ctx.link_subtypes.weirs;
                    let wr = ctx.link_subtypes.weir_row(j) as usize;
                    self.weirs_.link_idx[uk] = j;
                    self.weirs_.c_disch1[uk] = w.cd[wr];
                    self.weirs_.c_disch2[uk] = w.cd2[wr];
                    self.weirs_.can_surcharge[uk] = w.can_surcharge[wr];
                    self.weirs_.weir_type[uk] = w.weir_type[wr] as i32;
                    self.weirs_.end_con[uk] = w.end_contractions[wr];
                    self.weirs_.cd_curve[uk] = w.cd_curve[wr];
                    self.weirs_.road_width[uk] = w.road_width[wr];
                    self.weirs_.road_surface[uk] = w.road_surface[wr] as i32;
                    self.weirs_.slope[uk] = ctx.links.xsect_s_bot[uj];
                    {
                        let y_full = ctx.links.xsect_y_full[uj];
                        let route_step = ctx.options.link_validate_routing_step();
                        let l = 2.0 * route_step * (constants::GRAVITY * y_full).sqrt();
                        self.weirs_.length_eff[uk] = cmax(200.0, l);
                    }
                    iw += 1;
                }
                LinkType::Outlet => {
                    let uk = ix;
                    let o = &ctx.link_subtypes.outlets;
                    let xr = ctx.link_subtypes.outlet_row(j) as usize;
                    self.outlets_.link_idx[uk] = j;
                    self.outlets_.q_coeff[uk] = o.coeff[xr];
                    self.outlets_.q_expon[uk] = o.expon[xr];
                    self.outlets_.curve_idx[uk] = o.curve[xr];
                    ix += 1;
                }
                _ => {}
            }
        }

        // Flat index of all non-conduit links (DUMMY conduits included).
        self.nc_indices_.clear();
        self.nc_indices_.reserve((n_pumps + n_orifices + n_weirs + n_outlets) as usize);
        for j in 0..ctx.n_links() {
            let uj = j as usize;
            if ctx.links.type_[uj] != LinkType::Conduit || ctx.links.xsect_shape[uj] == XsectShape::Dummy {
                self.nc_indices_.push(j);
            }
        }

        self.nc_group_k_.clear();
        self.nc_group_k_.resize(ctx.n_links() as usize, -1);
        for k in 0..self.pumps_.count {
            self.nc_group_k_[self.pumps_.link_idx[k as usize] as usize] = k;
        }
        for k in 0..self.orifices_.count {
            self.nc_group_k_[self.orifices_.link_idx[k as usize] as usize] = k;
        }
        for k in 0..self.weirs_.count {
            self.nc_group_k_[self.weirs_.link_idx[k as usize] as usize] = k;
        }
        for k in 0..self.outlets_.count {
            self.nc_group_k_[self.outlets_.link_idx[k as usize] as usize] = k;
        }
    }

    pub fn non_conduit_indices(&self) -> &[i32] {
        &self.nc_indices_
    }
}
