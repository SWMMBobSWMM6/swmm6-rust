// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/data/MeshData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The unified triangle/quad mesh, structure of arrays. Every cell owns kMaxCellVerts = 4
// vertex and neighbour slots (a triangle's fourth holds vertex -1, neighbour -2).
//
// In form only: cell_centroid_xy reads a vertex coordinate through `vcoord`, which gives NaN
// for an index outside the vertex arrays, where upstream reads out of bounds (UPSTREAM_NOTES
// U300). Every in-range index gives upstream's arithmetic, in upstream's order.

#![allow(clippy::all)]

use crate::core::options::Bytes;

pub const K_MAX_CELL_VERTS: i32 = 4;

/// MeshData::TriCouplingRow: one [2D_TRIANGLE_NODE_MAP] row.
#[derive(Debug, Clone)]
pub struct TriCouplingRow {
    pub tri: i32,
    pub node: i32,
    pub node_name: Bytes,
    pub cd: f64,
    pub cd_set: bool,
    pub area: f64,
    pub area_set: bool,
}

impl Default for TriCouplingRow {
    fn default() -> Self {
        TriCouplingRow { tri: -1, node: -1, node_name: Bytes::new(), cd: 0.65, cd_set: false, area: 1.0, area_set: false }
    }
}

#[derive(Debug, Clone, Default)]
pub struct MeshData {
    pub vx: Vec<f64>,
    pub vy: Vec<f64>,
    pub vz: Vec<f64>,
    pub vtag: Vec<Bytes>,
    pub cell_nv: Vec<u8>,
    pub cell_v: Vec<i32>,
    pub cell_nbr: Vec<i32>,
    pub tri_area: Vec<f64>,
    pub tri_cx: Vec<f64>,
    pub tri_cy: Vec<f64>,
    pub tri_cz: Vec<f64>,
    pub edge_length: Vec<f64>,
    pub edge_nx: Vec<f64>,
    pub edge_ny: Vec<f64>,
    pub edge_mx: Vec<f64>,
    pub edge_my: Vec<f64>,
    pub edge_mz: Vec<f64>,
    pub edge_dist_c: Vec<f64>,
    pub quad_vfr_z: Vec<f64>,
    pub quad_vfr_a: Vec<f64>,
    pub edge_conveyance: Vec<f64>,
    pub mannings_n: Vec<f64>,
    pub tri_init_depth: Vec<f64>,
    pub tri_init_u: Vec<f64>,
    pub tri_init_v: Vec<f64>,
    pub tri_tag: Vec<Bytes>,
    pub vert_stencil_ptr: Vec<i32>,
    pub vert_stencil_idx: Vec<i32>,
    pub vert_stencil_wt: Vec<f64>,
    pub vert_coupled_node: Vec<i32>,
    pub tri_coupled_node: Vec<i32>,
    pub vert_coupling_cd: Vec<f64>,
    pub vert_coupling_cd_set: Vec<u8>,
    pub vert_coupling_area: Vec<f64>,
    pub vert_coupling_area_set: Vec<u8>,
    pub tri_coupling_cd: Vec<f64>,
    pub tri_coupling_area: Vec<f64>,
    pub vert_coupled_node_name: Vec<Bytes>,
    pub tri_coupled_node_name: Vec<Bytes>,
    pub tri_couplings: Vec<TriCouplingRow>,
}

impl MeshData {
    pub fn n_vertices(&self) -> i32 {
        self.vx.len() as i32
    }
    pub fn n_triangles(&self) -> i32 {
        self.cell_nv.len() as i32
    }
    pub fn n_cells(&self) -> i32 {
        self.cell_nv.len() as i32
    }
    pub fn n_quads(&self) -> i32 {
        let mut q = 0;
        for &nv in &self.cell_nv {
            q += (nv == 4) as i32;
        }
        q
    }
    pub fn edge_stride(&self) -> i32 {
        if self.n_quads() > 0 {
            K_MAX_CELL_VERTS
        } else {
            3
        }
    }
    pub fn n_edge_slots(&self) -> i32 {
        self.n_cells() * K_MAX_CELL_VERTS
    }
    pub const fn slot(c: i32, k: i32) -> i32 {
        c * K_MAX_CELL_VERTS + k
    }
    pub const fn slot_local(s: i32) -> i32 {
        s % K_MAX_CELL_VERTS
    }
    pub const fn slot_cell(s: i32) -> i32 {
        s / K_MAX_CELL_VERTS
    }
    pub fn cell_vertex_count(&self, c: i32) -> i32 {
        self.cell_nv[c as usize] as i32
    }
    pub fn cell_vertex(&self, c: i32, k: i32) -> i32 {
        self.cell_v[Self::slot(c, k) as usize]
    }

    /// A vertex coordinate; NaN for an index outside the arrays (see the file notes).
    fn vcoord(v: &[f64], i: i32) -> f64 {
        if i >= 0 && (i as usize) < v.len() {
            v[i as usize]
        } else {
            f64::NAN
        }
    }

    /// cell_centroid_xy(): true-area centroid of a cell from its vertices.
    pub fn cell_centroid_xy(&self, c: i32, x: &mut f64, y: &mut f64) {
        if self.cell_vertex_count(c) == 3 {
            let v0 = self.cell_vertex(c, 0);
            let v1 = self.cell_vertex(c, 1);
            let v2 = self.cell_vertex(c, 2);
            *x = (Self::vcoord(&self.vx, v0) + Self::vcoord(&self.vx, v1) + Self::vcoord(&self.vx, v2)) / 3.0;
            *y = (Self::vcoord(&self.vy, v0) + Self::vcoord(&self.vy, v1) + Self::vcoord(&self.vy, v2)) / 3.0;
            return;
        }
        let mut px = [0.0f64; 4];
        let mut py = [0.0f64; 4];
        for k in 0..4 {
            let v = self.cell_vertex(c, k);
            px[k as usize] = Self::vcoord(&self.vx, v);
            py[k as usize] = Self::vcoord(&self.vy, v);
        }
        let a1 = (0.5 * ((px[1] - px[0]) * (py[2] - py[0]) - (px[2] - px[0]) * (py[1] - py[0]))).abs();
        let a2 = (0.5 * ((px[2] - px[0]) * (py[3] - py[0]) - (px[3] - px[0]) * (py[2] - py[0]))).abs();
        let area = a1 + a2;
        if area > 0.0 {
            *x = (a1 * ((px[0] + px[1] + px[2]) / 3.0) + a2 * ((px[0] + px[2] + px[3]) / 3.0)) / area;
            *y = (a1 * ((py[0] + py[1] + py[2]) / 3.0) + a2 * ((py[0] + py[2] + py[3]) / 3.0)) / area;
        } else {
            *x = 0.25 * (px[0] + px[1] + px[2] + px[3]);
            *y = 0.25 * (py[0] + py[1] + py[2] + py[3]);
        }
    }

    pub fn cell_neighbour(&self, c: i32, k: i32) -> i32 {
        self.cell_nbr[Self::slot(c, k) as usize]
    }
    pub fn is_quad(&self, c: i32) -> bool {
        self.cell_nv[c as usize] == 4
    }
    pub fn cell_edge_vertices(&self, c: i32, k: i32, a: &mut i32, b: &mut i32) {
        let nv = self.cell_nv[c as usize] as i32;
        *a = self.cell_v[Self::slot(c, (k + 1) % nv) as usize];
        *b = self.cell_v[Self::slot(c, (k + 2) % nv) as usize];
    }
    pub fn set_triangle(&mut self, c: i32, v0: i32, v1: i32, v2: i32) {
        self.cell_nv[c as usize] = 3;
        self.cell_v[Self::slot(c, 0) as usize] = v0;
        self.cell_v[Self::slot(c, 1) as usize] = v1;
        self.cell_v[Self::slot(c, 2) as usize] = v2;
        self.cell_v[Self::slot(c, 3) as usize] = -1;
    }
    pub fn set_quad(&mut self, c: i32, v0: i32, v1: i32, v2: i32, v3: i32) {
        self.cell_nv[c as usize] = 4;
        self.cell_v[Self::slot(c, 0) as usize] = v0;
        self.cell_v[Self::slot(c, 1) as usize] = v1;
        self.cell_v[Self::slot(c, 2) as usize] = v2;
        self.cell_v[Self::slot(c, 3) as usize] = v3;
    }

    /// resize_vertices(): every per-vertex array to `nv` (old entries kept).
    pub fn resize_vertices(&mut self, nv: i32) {
        let n = nv as usize;
        self.vx.resize(n, 0.0);
        self.vy.resize(n, 0.0);
        self.vz.resize(n, 0.0);
        self.vtag.resize(n, Bytes::new());
        self.vert_coupled_node.resize(n, -1);
        self.vert_coupled_node_name.resize(n, Bytes::new());
        self.vert_coupling_cd.resize(n, 0.65);
        self.vert_coupling_cd_set.resize(n, 0);
        self.vert_coupling_area.resize(n, 1.0);
        self.vert_coupling_area_set.resize(n, 0);
    }

    /// resize_triangles(): every per-cell and per-slot array to `nt` cells (old entries kept;
    /// the padding slot of each new cell is vertex -1, neighbour -2).
    pub fn resize_triangles(&mut self, nt: i32) {
        let n = nt as usize;
        let old_n = self.cell_nv.len();
        self.cell_nv.resize(n, 3);
        let n4 = n * K_MAX_CELL_VERTS as usize;
        self.cell_v.resize(n4, 0);
        self.cell_nbr.resize(n4, -1);
        for c in old_n..n {
            self.cell_v[c * K_MAX_CELL_VERTS as usize + 3] = -1;
            self.cell_nbr[c * K_MAX_CELL_VERTS as usize + 3] = -2;
        }
        self.tri_area.resize(n, 0.0);
        self.tri_cx.resize(n, 0.0);
        self.tri_cy.resize(n, 0.0);
        self.tri_cz.resize(n, 0.0);
        self.edge_length.resize(n4, 0.0);
        self.edge_nx.resize(n4, 0.0);
        self.edge_ny.resize(n4, 0.0);
        self.edge_mx.resize(n4, 0.0);
        self.edge_my.resize(n4, 0.0);
        self.edge_mz.resize(n4, 0.0);
        self.edge_dist_c.resize(n4, 0.0);
        self.edge_conveyance.resize(n4, 1.0);
        self.mannings_n.resize(n, 0.035);
        self.tri_init_depth.resize(n, 0.0);
        self.tri_init_u.resize(n, 0.0);
        self.tri_init_v.resize(n, 0.0);
        self.tri_tag.resize(n, Bytes::new());
        self.tri_coupled_node.resize(n, -1);
        self.tri_coupled_node_name.resize(n, Bytes::new());
        self.tri_coupling_cd.resize(n, 0.65);
        self.tri_coupling_area.resize(n, 1.0);
    }
}
