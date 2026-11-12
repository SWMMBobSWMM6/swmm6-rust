// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/mesh/MeshBuilder.hpp and MeshBuilder.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// buildMeshTopology (cell neighbours, cell and edge geometry, the quad VFR side data) and
// validateMesh, as SurfaceRouter2D::initialize calls them (no progress callback).
//
// Differences in form, never in result:
//   * The edge map is a HashMap keyed by the sorted vertex pair; upstream's unordered_map with
//     the Cantor-pairing hash is only searched and inserted into, never iterated, so its order
//     is not observable.
//   * buildMeshTopology runs before validateMesh, so a cell naming a vertex outside the vertex
//     arrays reads out of bounds upstream (undefined). Here such a coordinate reads as NaN (as
//     MeshData::cell_centroid_xy does, U301). validateMesh then refuses the mesh on both sides;
//     what such a cell's garbage geometry would have been cannot reach a result, since the
//     only checks validateMesh makes before reaching that cell read other cells.

#![allow(clippy::all)]

use std::collections::HashMap;

use crate::core::options::Bytes;
use crate::twod::data::mesh_data::{MeshData, K_MAX_CELL_VERTS};
use crate::twod::mesh::quad_vfr::{quad_vfr_precompute, K_QUAD_VFR_Z};

/// A coordinate of vertex `i`, NaN outside the array (see the header).
fn at(v: &[f64], i: i32) -> f64 {
    if i >= 0 && (i as usize) < v.len() {
        v[i as usize]
    } else {
        f64::NAN
    }
}

fn make_edge_key(va: i32, vb: i32) -> (i32, i32) {
    (va.min(vb), va.max(vb))
}

/// buildMeshTopology(mesh).
pub fn build_mesh_topology(mesh: &mut MeshData) {
    let nt = mesh.n_triangles();

    // --- 1. Build edge-neighbour adjacency ---
    let mut edge_map: HashMap<(i32, i32), (i32, i32)> = HashMap::with_capacity(nt as usize * K_MAX_CELL_VERTS as usize);
    // Initialize neighbours to -1 (boundary) on real slots, -2 on padding.
    for t in 0..nt {
        let nv = mesh.cell_vertex_count(t);
        for k in 0..K_MAX_CELL_VERTS {
            mesh.cell_nbr[MeshData::slot(t, k) as usize] = if k < nv { -1 } else { -2 };
        }
    }
    for t in 0..nt {
        let nv = mesh.cell_vertex_count(t);
        for e in 0..nv {
            let (mut va, mut vb) = (0, 0);
            mesh.cell_edge_vertices(t, e, &mut va, &mut vb);
            let key = make_edge_key(va, vb);
            match edge_map.get(&key) {
                None => {
                    // First cell to claim this edge
                    edge_map.insert(key, (t, e));
                }
                Some(&(first_cell, first_edge)) => {
                    // Second cell: complete the pair
                    mesh.cell_nbr[MeshData::slot(t, e) as usize] = first_cell;
                    mesh.cell_nbr[MeshData::slot(first_cell, first_edge) as usize] = t;
                }
            }
        }
    }

    // --- 2. Compute cell geometry ---
    for t in 0..nt {
        let nv = mesh.cell_vertex_count(t);
        if nv == 3 {
            let v0 = mesh.cell_vertex(t, 0);
            let v1 = mesh.cell_vertex(t, 1);
            let v2 = mesh.cell_vertex(t, 2);
            let (x0, y0, z0) = (at(&mesh.vx, v0), at(&mesh.vy, v0), at(&mesh.vz, v0));
            let (x1, y1, z1) = (at(&mesh.vx, v1), at(&mesh.vy, v1), at(&mesh.vz, v1));
            let (x2, y2, z2) = (at(&mesh.vx, v2), at(&mesh.vy, v2), at(&mesh.vz, v2));
            // Centroid
            let (mut cx, mut cy) = (0.0, 0.0);
            mesh.cell_centroid_xy(t, &mut cx, &mut cy);
            mesh.tri_cx[t as usize] = cx;
            mesh.tri_cy[t as usize] = cy;
            mesh.tri_cz[t as usize] = (z0 + z1 + z2) / 3.0;
            // Area via cross product: 0.5 * ||(v1-v0) x (v2-v0)||
            let dx1 = x1 - x0;
            let dy1 = y1 - y0;
            let dx2 = x2 - x0;
            let dy2 = y2 - y0;
            mesh.tri_area[t as usize] = 0.5 * (dx1 * dy2 - dx2 * dy1).abs();
        } else {
            // Quadrilateral (v0 v1 v2 v3): shoelace area; centroid as the area-weighted
            // centroid of the two sub-triangles. Mean bed elevation over the four vertices.
            let mut x = [0.0f64; 4];
            let mut y = [0.0f64; 4];
            let mut z = [0.0f64; 4];
            for k in 0..4 {
                let v = mesh.cell_vertex(t, k);
                x[k as usize] = at(&mesh.vx, v);
                y[k as usize] = at(&mesh.vy, v);
                z[k as usize] = at(&mesh.vz, v);
            }
            let a1s = 0.5 * ((x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]));
            let a2s = 0.5 * ((x[2] - x[0]) * (y[3] - y[0]) - (x[3] - x[0]) * (y[2] - y[0]));
            let a1 = a1s.abs();
            let a2 = a2s.abs();
            let area = a1 + a2;
            let (mut cx, mut cy) = (0.0, 0.0);
            mesh.cell_centroid_xy(t, &mut cx, &mut cy);
            mesh.tri_cx[t as usize] = cx;
            mesh.tri_cy[t as usize] = cy;
            mesh.tri_cz[t as usize] = 0.25 * (z[0] + z[1] + z[2] + z[3]);
            mesh.tri_area[t as usize] = area;
        }
    }

    // Quad VFR side data (needs the connectivity + z only; sized on first quad).
    mesh.quad_vfr_z.clear();
    mesh.quad_vfr_a.clear();
    for t in 0..nt {
        compute_quad_vfr(mesh, t);
    }

    // --- 3. Compute edge geometry ---
    for t in 0..nt {
        let cx = mesh.tri_cx[t as usize];
        let cy = mesh.tri_cy[t as usize];
        let nv = mesh.cell_vertex_count(t);
        for e in 0..K_MAX_CELL_VERTS {
            let idx = MeshData::slot(t, e) as usize;
            if e >= nv {
                // Padding slot: zero geometry, so a stray read is inert.
                mesh.edge_mz[idx] = 0.0;
                mesh.edge_my[idx] = 0.0;
                mesh.edge_mx[idx] = 0.0;
                mesh.edge_length[idx] = 0.0;
                mesh.edge_ny[idx] = 0.0;
                mesh.edge_nx[idx] = 0.0;
                mesh.edge_dist_c[idx] = 0.0;
                continue;
            }
            let (mut va, mut vb) = (0, 0);
            mesh.cell_edge_vertices(t, e, &mut va, &mut vb);
            let (ax, ay, az) = (at(&mesh.vx, va), at(&mesh.vy, va), at(&mesh.vz, va));
            let (bx, by, bz) = (at(&mesh.vx, vb), at(&mesh.vy, vb), at(&mesh.vz, vb));
            // Edge midpoint
            mesh.edge_mx[idx] = 0.5 * (ax + bx);
            mesh.edge_my[idx] = 0.5 * (ay + by);
            mesh.edge_mz[idx] = 0.5 * (az + bz);
            // Edge length (planimetric)
            let dx = bx - ax;
            let dy = by - ay;
            mesh.edge_length[idx] = (dx * dx + dy * dy).sqrt();
            // Outward unit normal
            let mut nx = dy;
            let mut ny = -dx;
            let to_mid_x = mesh.edge_mx[idx] - cx;
            let to_mid_y = mesh.edge_my[idx] - cy;
            if nx * to_mid_x + ny * to_mid_y < 0.0 {
                nx = -nx;
                ny = -ny;
            }
            // Normalize
            let len = (nx * nx + ny * ny).sqrt();
            if len > 1.0e-15 {
                mesh.edge_nx[idx] = nx / len;
                mesh.edge_ny[idx] = ny / len;
            } else {
                mesh.edge_nx[idx] = 0.0;
                mesh.edge_ny[idx] = 0.0;
            }
            // Centroid -> edge normal distance (the ghost-cell arm).
            mesh.edge_dist_c[idx] = to_mid_x * mesh.edge_nx[idx] + to_mid_y * mesh.edge_ny[idx];
        }
    }
}

/// validateMesh(mesh): empty when the mesh is valid, else the first problem.
pub fn validate_mesh(mesh: &MeshData) -> Bytes {
    let nv_total = mesh.n_vertices();
    let nt = mesh.n_triangles();
    if nv_total < 3 {
        return b"Mesh must have at least 3 vertices".to_vec();
    }
    if nt < 1 {
        return b"Mesh must have at least 1 cell".to_vec();
    }
    for t in 0..nt {
        let nv = mesh.cell_vertex_count(t);
        let kind: &str = if nv == 4 { "Quad " } else { "Triangle " };
        if nv != 3 && nv != 4 {
            return format!("Cell {} has unsupported vertex count {}", t, nv).into_bytes();
        }
        let mut v = [-1i32; 4];
        for k in 0..nv {
            v[k as usize] = mesh.cell_vertex(t, k);
            if v[k as usize] < 0 || v[k as usize] >= nv_total {
                return format!("{}{} has out-of-range vertex index", kind, t).into_bytes();
            }
        }
        for a in 0..nv as usize {
            for b in a + 1..nv as usize {
                if v[a] == v[b] {
                    return format!("{}{} has duplicate vertex indices", kind, t).into_bytes();
                }
            }
        }
        if mesh.tri_area[t as usize] <= 0.0 {
            return format!("{}{} has zero or negative area", kind, t).into_bytes();
        }
        if nv == 4 {
            // Convex and simple: every consecutive edge-pair cross product has the same sign.
            let mut sign = 0;
            for k in 0..4usize {
                let p = v[k] as usize;
                let q = v[(k + 1) % 4] as usize;
                let r = v[(k + 2) % 4] as usize;
                let cr = (mesh.vx[q] - mesh.vx[p]) * (mesh.vy[r] - mesh.vy[q]) - (mesh.vy[q] - mesh.vy[p]) * (mesh.vx[r] - mesh.vx[q]);
                let s = if cr > 0.0 {
                    1
                } else if cr < 0.0 {
                    -1
                } else {
                    0
                };
                if s == 0 || (sign != 0 && s != sign) {
                    return format!("Quad {} (vertices {} {} {} {}) is not a convex quadrilateral", t, v[0], v[1], v[2], v[3]).into_bytes();
                }
                sign = s;
            }
        }
        if mesh.mannings_n[t as usize] <= 0.0 {
            return format!("{}{} has non-positive Manning's n", kind, t).into_bytes();
        }
    }
    Bytes::new()
}

/// computeQuadVfr(mesh, t): the quad VFR side data of cell t (no-op for triangles). Allocates
/// the arrays on first use so all-triangle meshes carry nothing.
fn compute_quad_vfr(mesh: &mut MeshData, t: i32) {
    if mesh.cell_vertex_count(t) != 4 {
        return;
    }
    let nt = mesh.n_triangles() as usize;
    let kz = K_QUAD_VFR_Z as usize;
    if mesh.quad_vfr_z.len() != nt * kz {
        mesh.quad_vfr_z.clear();
        mesh.quad_vfr_z.resize(nt * kz, 0.0);
    }
    if mesh.quad_vfr_a.len() != nt * 2 {
        mesh.quad_vfr_a.clear();
        mesh.quad_vfr_a.resize(nt * 2, 0.0);
    }
    let mut x = [0.0f64; 4];
    let mut y = [0.0f64; 4];
    let mut z = [0.0f64; 4];
    for k in 0..4 {
        let v = mesh.cell_vertex(t, k);
        x[k as usize] = at(&mesh.vx, v);
        y[k as usize] = at(&mesh.vy, v);
        z[k as usize] = at(&mesh.vz, v);
    }
    let mut a1 = 0.0;
    let mut a2 = 0.0;
    let base = t as usize * kz;
    quad_vfr_precompute(&x, &y, &z, &mut mesh.quad_vfr_z[base..base + kz], &mut a1, &mut a2);
    mesh.quad_vfr_a[t as usize * 2] = a1;
    mesh.quad_vfr_a[t as usize * 2 + 1] = a2;
}
