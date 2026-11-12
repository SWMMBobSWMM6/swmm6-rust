// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/solver/InertialEdges.hpp and InertialEdges.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The marcher's unique-interior-edge table (structure of arrays) and the per-cell CSR
// incidence, built once from the mesh. std::hypot is f64::hypot (glibc's hypot on the gate
// machine; see PLAN.md, "The C math library").

#![allow(clippy::all, non_snake_case)]

use crate::core::cmath::{cmax, cmin};
use crate::twod::data::mesh_data::{MeshData, K_MAX_CELL_VERTS};

#[derive(Debug, Clone, Default)]
pub struct InertialEdges {
    pub ne: i32,
    pub cL: Vec<i32>,
    pub cR: Vec<i32>,
    pub xi: Vec<f64>,
    pub conveyance: Vec<f64>,
    pub inv_dx: Vec<f64>,
    pub zface: Vec<f64>,
    pub ze_lo: Vec<f64>,
    pub ze_hi: Vec<f64>,
    pub slotL: Vec<i32>,
    pub slotR: Vec<i32>,
    pub nx: Vec<f64>,
    pub ny: Vec<f64>,
    pub mx: Vec<f64>,
    pub my: Vec<f64>,
    pub inv_dx_normal: Vec<f64>,
    pub n2_face: Vec<f64>,
    pub n_face: Vec<f64>,
    pub cell_lchar: Vec<f64>,
    pub cell_lpos: Vec<f64>,
    pub cell_ptr: Vec<i32>,
    pub cell_edge: Vec<i32>,
    pub cell_sign: Vec<i8>,
    pub cell_arm_x: Vec<f64>,
    pub cell_arm_y: Vec<f64>,
}

impl InertialEdges {
    pub fn empty(&self) -> bool {
        self.ne == 0
    }

    /// build(mesh).
    pub fn build(&mut self, mesh: &MeshData) {
        let nt = mesh.n_triangles();
        self.cL.clear();
        self.cR.clear();
        self.xi.clear();
        self.inv_dx.clear();
        self.zface.clear();
        self.ze_lo.clear();
        self.ze_hi.clear();
        self.slotL.clear();
        self.slotR.clear();
        self.nx.clear();
        self.ny.clear();
        self.mx.clear();
        self.my.clear();
        self.inv_dx_normal.clear();
        self.n2_face.clear();
        self.n_face.clear();
        self.cell_lchar.clear();
        self.cell_lpos.clear();

        // slot_edge[t][e] = unique-edge id incident to (t, local edge e), or -1.
        let mut slot_edge: Vec<[i32; K_MAX_CELL_VERTS as usize]> = vec![[-1, -1, -1, -1]; nt as usize];

        // 1. Enumerate unique interior edges (each once, nbr > t).
        for t in 0..nt {
            let ut = t as usize;
            let nv_t = mesh.cell_vertex_count(t);
            for e in 0..nv_t {
                let nb = mesh.cell_neighbour(t, e);
                if nb < 0 || nb < t {
                    continue; // boundary, or already counted
                }
                let un = nb as usize;
                let eid = self.cL.len() as i32;
                self.cL.push(t);
                self.cR.push(nb);
                self.xi.push(mesh.edge_length[MeshData::slot(t, e) as usize]);
                let ddx = (mesh.tri_cx[ut] - mesh.tri_cx[un]).hypot(mesh.tri_cy[ut] - mesh.tri_cy[un]); // LIBM (hypot)
                self.inv_dx.push(if ddx > 1.0e-12 { 1.0 / ddx } else { 0.0 });
                self.zface.push(cmax(mesh.tri_cz[ut], mesh.tri_cz[un]));
                {
                    let (mut va, mut vb) = (0, 0);
                    mesh.cell_edge_vertices(t, e, &mut va, &mut vb);
                    let za = mesh.vz[va as usize];
                    let zb = mesh.vz[vb as usize];
                    self.ze_lo.push(cmin(za, zb));
                    self.ze_hi.push(cmax(za, zb));
                }
                self.slotL.push(MeshData::slot(t, e));
                // Find nb's local edge facing t (the mirror slot) and record it.
                let nv_nb = mesh.cell_vertex_count(nb);
                let mut e2 = 0;
                while e2 < nv_nb {
                    if mesh.cell_neighbour(nb, e2) == t {
                        break;
                    }
                    e2 += 1;
                }
                self.slotR.push(MeshData::slot(nb, e2));
                slot_edge[ut][e as usize] = eid;
                if e2 < nv_nb {
                    slot_edge[un][e2 as usize] = eid;
                }
                let sl = MeshData::slot(t, e) as usize;
                self.nx.push(mesh.edge_nx[sl]);
                self.ny.push(mesh.edge_ny[sl]);
                self.mx.push(mesh.edge_mx[sl]);
                self.my.push(mesh.edge_my[sl]);
                {
                    let dxc = mesh.tri_cx[un] - mesh.tri_cx[ut];
                    let dyc = mesh.tri_cy[un] - mesh.tri_cy[ut];
                    let chord = dxc.hypot(dyc); // LIBM (hypot)
                    let mut dn = (dxc * mesh.edge_nx[sl] + dyc * mesh.edge_ny[sl]).abs();
                    dn = cmax(dn, 0.3 * chord); // near-degenerate floor
                    self.inv_dx_normal.push(if dn > 1.0e-12 { 1.0 / dn } else { 0.0 });
                }
                {
                    let nf = 0.5 * (mesh.mannings_n[ut] + mesh.mannings_n[un]);
                    self.n2_face.push(nf * nf);
                    self.n_face.push(nf);
                }
            }
        }
        self.ne = self.cL.len() as i32;
        self.conveyance.clear();
        self.conveyance.resize(self.ne as usize, 1.0);
        self.refresh_conveyance(mesh);

        // Per-cell characteristic length for the CFL bound.
        self.cell_lchar.clear();
        self.cell_lchar.resize(nt as usize, 0.0);
        {
            let mut s = vec![0.0f64; nt as usize];
            for e in 0..self.ne as usize {
                let sv = self.xi[e] * self.inv_dx_normal[e];
                s[self.cL[e] as usize] += sv;
                s[self.cR[e] as usize] += sv;
            }
            for t in 0..nt {
                let ut = t as usize;
                if s[ut] > 1.0e-30 {
                    self.cell_lchar[ut] = (2.0 * mesh.tri_area[ut] / s[ut]).sqrt();
                } else if mesh.cell_vertex_count(t) == 3 {
                    let mut xi_max = 0.0;
                    for e in 0..3 {
                        xi_max = cmax(xi_max, mesh.edge_length[MeshData::slot(t, e) as usize]);
                    }
                    self.cell_lchar[ut] = if xi_max > 0.0 { 2.0 * mesh.tri_area[ut] / xi_max } else { 0.0 };
                } else {
                    // Quad with no interior face: twice the smallest centroid -> edge distance.
                    let mut dmin = 1.0e300;
                    for e in 0..mesh.cell_vertex_count(t) {
                        dmin = cmin(dmin, mesh.edge_dist_c[MeshData::slot(t, e) as usize]);
                    }
                    self.cell_lchar[ut] = if dmin > 0.0 && dmin < 1.0e300 { 2.0 * dmin } else { 0.0 };
                }
            }
        }
        self.cell_lpos.clear();
        self.cell_lpos.resize(nt as usize, 0.0);
        for t in 0..nt {
            let ut = t as usize;
            let mut p = 0.0;
            for e in 0..mesh.cell_vertex_count(t) {
                p += mesh.edge_length[MeshData::slot(t, e) as usize];
            }
            self.cell_lpos[ut] = if p > 0.0 { 2.0 * mesh.tri_area[ut] / p } else { 0.0 };
        }

        // 2. Per-cell CSR incidence with orientation signs.
        self.cell_ptr.clear();
        self.cell_ptr.resize(nt as usize + 1, 0);
        for t in 0..nt as usize {
            let mut c = 0;
            for e in 0..K_MAX_CELL_VERTS as usize {
                if slot_edge[t][e] >= 0 {
                    c += 1;
                }
            }
            self.cell_ptr[t + 1] = self.cell_ptr[t] + c;
        }
        let total = self.cell_ptr[nt as usize] as usize;
        self.cell_edge.clear();
        self.cell_edge.resize(total, 0);
        self.cell_sign.clear();
        self.cell_sign.resize(total, 0);
        self.cell_arm_x.clear();
        self.cell_arm_x.resize(total, 0.0);
        self.cell_arm_y.clear();
        self.cell_arm_y.resize(total, 0.0);
        let mut fill = vec![0i32; nt as usize];
        for t in 0..nt as usize {
            for e in 0..K_MAX_CELL_VERTS as usize {
                let eid = slot_edge[t][e];
                if eid < 0 {
                    continue;
                }
                let pos = (self.cell_ptr[t] + fill[t]) as usize;
                fill[t] += 1;
                let ue = eid as usize;
                self.cell_edge[pos] = eid;
                self.cell_sign[pos] = if self.cL[ue] == t as i32 { 1 } else { -1 };
                self.cell_arm_x[pos] = self.mx[ue] - mesh.tri_cx[t];
                self.cell_arm_y[pos] = self.my[ue] - mesh.tri_cy[t];
            }
        }
    }

    /// refreshConveyance(mesh): the face conveyance factors from the mesh edge slots.
    pub fn refresh_conveyance(&mut self, mesh: &MeshData) -> bool {
        let mut changed = false;
        for e in 0..self.ne as usize {
            let slot = self.slotL[e] as usize;
            let factor = if slot < mesh.edge_conveyance.len() { mesh.edge_conveyance[slot] } else { 1.0 };
            changed |= self.conveyance[e] != factor;
            self.conveyance[e] = factor;
        }
        changed
    }
}
