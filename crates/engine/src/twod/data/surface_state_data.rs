// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/data/SurfaceStateData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The 2D surface state: volumes, reconstructed depths and heads, gradients, forcing arrays,
// statistics. Members in upstream's order with upstream's names.
//
// Differences in form, never in result: upstream's non-owning pointers are not pointers here.
// Each names a member of the router or of the context, which SurfaceRouter2D::initialize
// attaches; the translation passes the targets explicitly where they are read, and keeps the
// attachment as a flag (true = attached to its one possible target):
//   boundary            -> SurfaceRouter2D::boundary_
//   nodes_1d            -> SimulationContext::nodes
//   node_row_conc       -> SurfaceRouter2D::node_row_conc_
//   node_coupling       -> SurfaceRouter2D::node_coupling_points_
//   coupling_head_cells -> SurfaceRouter2D::coupling_head_cells_
//   node_pond_area_m2   -> SurfaceRouter2D::node_pond_area_m2_
//   snow_context        -> the engine's SimulationContext
// node_pond_base has two possible targets, so it is a PondBase.

#![allow(clippy::all)]

use crate::twod::coupling::runtime_coupling::{RuntimeForcings, RuntimeSources};
use crate::twod::data::mesh_data::K_MAX_CELL_VERTS;
use crate::twod::data::surface_transport_state::SurfaceTransportState;
use crate::twod::heat::mesh_heat::MeshHeatState;
use crate::twod::heat::mesh_snow::MeshSnow;

/// The target of SurfaceStateData::node_pond_base.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PondBase {
    #[default]
    None = 0,
    /// &ctx.nodes.rpt_full_volume (dynamic wave)
    RptFullVolume = 1,
    /// &ctx.nodes.full_volume (other routing)
    FullVolume = 2,
}

#[derive(Debug, Clone, Default)]
pub struct SurfaceStateData {
    pub boundary: bool,
    pub nodes_1d: bool,
    pub node_row_conc: bool,
    pub node_coupling: bool,
    pub node_pond_base: PondBase,
    pub coupling_head_cells: bool,
    pub node_pond_area_m2: bool,
    pub coupling_window_s: f64,
    pub transport: SurfaceTransportState,
    pub heat: MeshHeatState,
    pub depth: Vec<f64>,
    pub head: Vec<f64>,
    pub volume: Vec<f64>,
    pub grad_hx: Vec<f64>,
    pub grad_hy: Vec<f64>,
    pub grad_hx_lim: Vec<f64>,
    pub grad_hy_lim: Vec<f64>,
    pub vert_head: Vec<f64>,
    pub vert_depth_signed: Vec<f64>,
    pub face_vx: Vec<f64>,
    pub face_vy: Vec<f64>,
    pub cell_continuity_err: Vec<f64>,
    pub edge_flux: Vec<f64>,
    pub snow: MeshSnow,
    pub snow_context: bool,
    pub snow_catch: Vec<f64>,
    pub rainfall: Vec<f64>,
    pub evap_rate: Vec<f64>,
    pub infil_rate: Vec<f64>,
    pub infil_applied: Vec<f64>,
    pub coupling_flux: Vec<f64>,
    pub coupling_applied: Vec<f64>,
    pub boundary_receipt_slots: Vec<i32>,
    pub boundary_in_applied: Vec<f64>,
    pub boundary_out_applied: Vec<f64>,
    pub boundary_requested: Vec<f64>,
    pub boundary_last_flow: Vec<f64>,
    pub boundary_last_dt: Vec<f64>,
    pub net_source: Vec<f64>,
    pub rainfall_forced: Vec<i8>,
    pub rainfall_persist: Vec<i8>,
    pub rainfall_force_val: Vec<f64>,
    pub evap_forced: Vec<i8>,
    pub evap_persist: Vec<i8>,
    pub evap_force_val: Vec<f64>,
    pub coupling_forced: Vec<i8>,
    pub coupling_persist: Vec<i8>,
    pub coupling_force_val: Vec<f64>,
    pub runtime_sources: RuntimeSources,
    pub runtime_forcings: RuntimeForcings,
    pub coupling_native: Vec<f64>,
    pub forcing_dirty: bool,
    pub forcing_ever_set: bool,
    pub old_depth: Vec<f64>,
    pub old_volume: Vec<f64>,
    pub stat_max_depth: Vec<f64>,
    pub stat_max_velocity: Vec<f64>,
    pub stat_max_cont_err: Vec<f64>,
    pub stat_cum_volume: Vec<f64>,
    pub evap_loss_total: f64,
    pub rk2_rain_pending: Vec<f64>,
}

fn assign<T: Clone>(v: &mut Vec<T>, n: usize, x: T) {
    v.clear();
    v.resize(n, x);
}

impl SurfaceStateData {
    /// liquidRainfall(i): the liquid part of the cell's rain (the snow owner's, or the runtime
    /// forcing's over the native rate).
    pub fn liquid_rainfall(&self, i: i32) -> f64 {
        if self.snow.owns(i) {
            self.snow.liquid_rates()[i as usize]
        } else {
            self.runtime_forcings.rate(1, i, self.rainfall[i as usize])
        }
    }

    pub fn resize(&mut self, n_triangles: i32, n_vertices: i32) {
        self.heat.clear();
        self.snow = MeshSnow::default();
        self.snow_context = false;
        self.snow_catch.clear();
        let nt = n_triangles as usize;
        let nv = n_vertices as usize;
        let n3 = nt * K_MAX_CELL_VERTS as usize;
        assign(&mut self.depth, nt, 0.0);
        assign(&mut self.head, nt, 0.0);
        assign(&mut self.volume, nt, 0.0);
        assign(&mut self.grad_hx, nt, 0.0);
        assign(&mut self.grad_hy, nt, 0.0);
        assign(&mut self.grad_hx_lim, nt, 0.0);
        assign(&mut self.grad_hy_lim, nt, 0.0);
        assign(&mut self.vert_head, nv, 0.0);
        assign(&mut self.vert_depth_signed, nv, 0.0);
        assign(&mut self.face_vx, nt, 0.0);
        assign(&mut self.face_vy, nt, 0.0);
        assign(&mut self.cell_continuity_err, nt, 0.0);
        assign(&mut self.edge_flux, n3, 0.0);
        assign(&mut self.rainfall, nt, 0.0);
        assign(&mut self.evap_rate, nt, 0.0);
        assign(&mut self.infil_rate, nt, 0.0);
        assign(&mut self.infil_applied, nt, 0.0);
        assign(&mut self.coupling_applied, nt, 0.0);
        self.boundary_receipt_slots.clear();
        self.boundary_in_applied.clear();
        self.boundary_out_applied.clear();
        self.boundary_requested.clear();
        self.boundary_last_flow.clear();
        self.boundary_last_dt.clear();
        assign(&mut self.coupling_flux, nt, 0.0);
        self.coupling_native.clear();
        self.runtime_sources = RuntimeSources::default();
        self.runtime_forcings = RuntimeForcings::default();
        self.rk2_rain_pending.clear();
        assign(&mut self.net_source, nt, 0.0);
        assign(&mut self.rainfall_forced, nt, 0i8);
        assign(&mut self.rainfall_persist, nt, 0i8);
        assign(&mut self.rainfall_force_val, nt, 0.0);
        assign(&mut self.evap_forced, nt, 0i8);
        assign(&mut self.evap_persist, nt, 0i8);
        assign(&mut self.evap_force_val, nt, 0.0);
        assign(&mut self.coupling_forced, nt, 0i8);
        assign(&mut self.coupling_persist, nt, 0i8);
        assign(&mut self.coupling_force_val, nt, 0.0);
        assign(&mut self.old_depth, nt, 0.0);
        assign(&mut self.old_volume, nt, 0.0);
        assign(&mut self.stat_max_depth, nt, 0.0);
        assign(&mut self.stat_max_velocity, nt, f64::NAN); // std::numeric_limits<double>::quiet_NaN()
        assign(&mut self.stat_max_cont_err, nt, 0.0);
        assign(&mut self.stat_cum_volume, nt, 0.0);
        self.evap_loss_total = 0.0;
    }
}
