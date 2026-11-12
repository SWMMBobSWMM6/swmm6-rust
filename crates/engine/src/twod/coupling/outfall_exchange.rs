// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/coupling/OutfallExchange.hpp and
// OutfallExchange.cpp (initialize, resetWindow) (HydroCouple/openswmm.engine, branch swmm6_rel,
// commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The outfall exchange's donor cells, SPREAD receivers and per-node ledgers, built at
// initialize. prepare / apply / recordTailwater (the step) come with wave 4. The two
// std::runtime_error throws are Err(what), caught by SWMMEngine::initHydraulics' handler.
// `users.at(c)` cannot go out of range: every cell comes from the mesh's own stencils.

#![allow(clippy::all)]

use crate::core::cmath::cmax;
use crate::core::context::SimulationContext;
use crate::core::options::Bytes;
use crate::twod::coupling::node_coupling::CouplingPoint;
use crate::twod::data::mesh_data::MeshData;

/// OutfallExchange::Point.
#[derive(Debug, Clone, Default)]
pub struct Point {
    pub index: i32,
    pub node: i32,
    pub cells: Vec<i32>,
    pub fraction: Vec<f64>,
    pub reserved: Vec<f64>,
    pub total: f64,
    pub area: f64,
    pub spread: Vec<i32>,
    pub bed_ref: f64,
}

#[derive(Debug, Clone, Default)]
pub struct OutfallExchange {
    pub points_: Vec<Point>,
    pub recv_: Vec<i32>,
    pub nodes_: Vec<i32>,
    pub degree_: Vec<i32>,
    pub node_area_: Vec<f64>,
    pub node_reserved_: Vec<f64>,
    pub tw_engaged_s_: Vec<f64>,
    pub tw_submerged_s_: Vec<f64>,
    pub tw_max_: Vec<f64>,
    pub pending_: Vec<f64>,
    pub prepared_dt_: f64,
    pub prepared_: bool,
}

fn assign<T: Clone>(v: &mut Vec<T>, n: usize, x: T) {
    v.clear();
    v.resize(n, x);
}

impl OutfallExchange {
    pub fn reset_window(&mut self) {
        for x in self.pending_.iter_mut() {
            *x = 0.0;
        }
    }

    /// initialize(mesh, cps, ctx).
    pub fn initialize(&mut self, mesh: &MeshData, cps: &[CouplingPoint], ctx: &mut SimulationContext) -> Result<(), Bytes> {
        self.points_.clear();
        self.recv_.clear();
        self.nodes_.clear();
        self.pending_.clear();
        self.prepared_dt_ = 0.0;
        self.prepared_ = false;
        ctx.surface_outfall_link_limit.clear();
        let nn = ctx.n_nodes() as usize;
        assign(&mut self.degree_, nn, 0);
        assign(&mut self.node_area_, nn, 0.0);
        assign(&mut self.node_reserved_, nn, 0.0);
        assign(&mut self.tw_engaged_s_, nn, 0.0);
        assign(&mut self.tw_submerged_s_, nn, 0.0);
        assign(&mut self.tw_max_, nn, -1.0e30);
        let mut users = vec![0i32; mesh.n_cells() as usize];
        let mut seen: Vec<u8> = Vec::new();
        for (k, cp) in cps.iter().enumerate() {
            if !cp.is_outfall {
                continue;
            }
            let mut p = Point::default();
            p.index = k as i32;
            p.node = cp.node_idx;
            if cp.vertex_idx < 0 {
                p.cells.push(cp.cell_idx);
            } else {
                let v = cp.vertex_idx as usize;
                for j in mesh.vert_stencil_ptr[v]..mesh.vert_stencil_ptr[v + 1] {
                    p.cells.push(mesh.vert_stencil_idx[j as usize]);
                }
            }
            if p.cells.is_empty() {
                return Err(b"2D outfall has no donor cells.".to_vec());
            }
            for &c in &p.cells {
                users[c as usize] += 1;
                p.area += mesh.tri_area[c as usize];
            }
            self.node_area_[p.node as usize] += p.area;
            p.reserved.resize(p.cells.len(), 0.0);
            // SPREAD r: the cells whose centroid lies within r of the point, gathered outward
            // from its own cells, receive its discharge.
            if cp.spread > 0.0 {
                if seen.is_empty() {
                    assign(&mut seen, mesh.n_cells() as usize, 0u8);
                }
                let cx = if cp.vertex_idx >= 0 { mesh.vx[cp.vertex_idx as usize] } else { mesh.tri_cx[cp.cell_idx as usize] };
                let cy = if cp.vertex_idx >= 0 { mesh.vy[cp.vertex_idx as usize] } else { mesh.tri_cy[cp.cell_idx as usize] };
                p.spread = p.cells.clone();
                p.bed_ref = f64::NEG_INFINITY;
                let mut marked: Vec<i32> = p.cells.clone();
                for &c in &p.cells {
                    seen[c as usize] = 1;
                    p.bed_ref = cmax(p.bed_ref, mesh.tri_cz[c as usize]);
                }
                let mut q = 0usize;
                while q < p.spread.len() {
                    let c = p.spread[q];
                    for e in 0..mesh.cell_vertex_count(c) {
                        let nb = mesh.cell_neighbour(c, e);
                        if nb < 0 || seen[nb as usize] != 0 {
                            continue;
                        }
                        seen[nb as usize] = 1;
                        marked.push(nb);
                        let dx = mesh.tri_cx[nb as usize] - cx;
                        let dy = mesh.tri_cy[nb as usize] - cy;
                        if dx * dx + dy * dy <= cp.spread * cp.spread {
                            p.spread.push(nb);
                        }
                    }
                    q += 1;
                }
                for &c in &marked {
                    seen[c as usize] = 0;
                }
            }
            if !self.nodes_.contains(&p.node) {
                self.nodes_.push(p.node);
            }
            self.points_.push(p);
        }
        if self.points_.is_empty() {
            return Ok(());
        }
        assign(&mut self.pending_, mesh.n_cells() as usize, 0.0);
        for j in 0..ctx.n_links() as usize {
            if ctx.links.node1[j] >= 0 {
                self.degree_[ctx.links.node1[j] as usize] += 1;
            }
            if ctx.links.node2[j] >= 0 {
                self.degree_[ctx.links.node2[j] as usize] += 1;
            }
        }
        for p in self.points_.iter_mut() {
            for &c in &p.cells {
                p.fraction.push(1.0 / users[c as usize] as f64);
            }
            if self.degree_[p.node as usize] == 0 {
                return Err(b"2D outfall has no incident hydraulic link.".to_vec());
            }
        }
        assign(&mut ctx.surface_outfall_link_limit, nn, f64::INFINITY);
        Ok(())
    }
}
