// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/mesh/VertexReconstruction.hpp and
// VertexReconstruction.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The pseudo-Laplacian vertex stencils and the vertex reconstructions (heads, render depths).
// Upstream's OpenMP loops are `schedule(static)` with each iteration writing only its own
// entries, so a serial loop gives the same bits at any thread count.

#![allow(clippy::all)]

use crate::twod::data::mesh_data::MeshData;
use crate::twod::data::surface_state_data::SurfaceStateData;
use crate::twod::mesh::quad_vfr::{quad_eta_from_mean_depth, K_QUAD_VFR_Z};
use crate::twod::mesh::vfr_closure::{vfr_eta_from_mean_depth, vfr_sort3};

/// buildVertexStencils(mesh).
pub fn build_vertex_stencils(mesh: &mut MeshData) {
    let nv = mesh.n_vertices();
    let nt = mesh.n_triangles();

    // Step 1: For each vertex, collect all triangles that share it
    let mut vert_triangles: Vec<Vec<i32>> = vec![Vec::new(); nv as usize];
    for t in 0..nt {
        let nvc = mesh.cell_vertex_count(t);
        for k in 0..nvc {
            vert_triangles[mesh.cell_vertex(t, k) as usize].push(t);
        }
    }

    // Step 2: Build CSR stencil
    mesh.vert_stencil_ptr.resize(nv as usize + 1, 0);
    mesh.vert_stencil_idx.clear();
    mesh.vert_stencil_wt.clear();

    let mut ptr = 0i32;
    for b in 0..nv as usize {
        mesh.vert_stencil_ptr[b] = ptr;
        let tris = &vert_triangles[b];
        let m = tris.len() as i32;
        if m == 0 {
            // Isolated vertex: no stencil
            continue;
        }

        let xb = mesh.vx[b];
        let yb = mesh.vy[b];

        // Compute moments (Eq. [20]-[21])
        let mut i_xx = 0.0;
        let mut i_yy = 0.0;
        let mut i_xy = 0.0;
        for &t in tris.iter() {
            let dx = mesh.tri_cx[t as usize] - xb;
            let dy = mesh.tri_cy[t as usize] - yb;
            i_xx += dx * dx;
            i_yy += dy * dy;
            i_xy += dx * dy;
        }

        // Determinant of the moment matrix
        let det = i_xx * i_yy - i_xy * i_xy;

        if det.abs() < 1.0e-30 {
            // Degenerate (collinear stencil): use uniform weights
            let w = 1.0 / m as f64;
            for &t in tris.iter() {
                mesh.vert_stencil_idx.push(t);
                mesh.vert_stencil_wt.push(w);
            }
            ptr += m;
            continue;
        }

        let inv_det = 1.0 / det;

        let mut sx = 0.0;
        let mut sy = 0.0;
        for &t in tris.iter() {
            sx += mesh.tri_cx[t as usize] - xb;
            sy += mesh.tri_cy[t as usize] - yb;
        }
        let r_x = -sx / m as f64;
        let r_y = -sy / m as f64;

        let lambda_x = (i_yy * r_x - i_xy * r_y) * inv_det;
        let lambda_y = (i_xx * r_y - i_xy * r_x) * inv_det;

        // Compute weights and clip negatives (Jawahar & Kamath boundary treatment)
        let mut weights = vec![0.0f64; m as usize];
        let mut w_sum = 0.0;
        for (i, &t) in tris.iter().enumerate() {
            let dx = mesh.tri_cx[t as usize] - xb;
            let dy = mesh.tri_cy[t as usize] - yb;
            weights[i] = (1.0 / m as f64) + lambda_x * dx + lambda_y * dy;
            if weights[i] < 0.0 {
                weights[i] = 0.0; // Clip extraneous weights
            }
            w_sum += weights[i];
        }

        // Renormalize to ensure partition of unity; all clipped -> uniform weights.
        if w_sum > 1.0e-30 {
            for w in weights.iter_mut() {
                *w /= w_sum;
            }
        } else {
            let w_uniform = 1.0 / m as f64;
            for w in weights.iter_mut() {
                *w = w_uniform;
            }
        }

        for (i, &t) in tris.iter().enumerate() {
            mesh.vert_stencil_idx.push(t);
            mesh.vert_stencil_wt.push(weights[i]);
        }
        ptr += m;
    }
    mesh.vert_stencil_ptr[nv as usize] = ptr;
}

/// reconstructVertexHeads(mesh, state, nthreads): each vertex's head from its stencil.
pub fn reconstruct_vertex_heads(mesh: &MeshData, state: &mut SurfaceStateData, _nthreads: i32) {
    let nv = mesh.n_vertices();
    for b in 0..nv as usize {
        let start = mesh.vert_stencil_ptr[b] as usize;
        let end = mesh.vert_stencil_ptr[b + 1] as usize;
        let mut h = 0.0;
        for k in start..end {
            h += mesh.vert_stencil_wt[k] * state.head[mesh.vert_stencil_idx[k] as usize];
        }
        state.vert_head[b] = h;
    }
}

/// cellFreeSurfaceElevationOf(mesh, t, mean_depth).
pub fn cell_free_surface_elevation_of(mesh: &MeshData, t: i32, mean_depth: f64) -> f64 {
    if mesh.cell_vertex_count(t) == 4 {
        // Quad: exact (eps = 0) B&S 2007 two-plane relation.
        let u = t as usize;
        let kz = K_QUAD_VFR_Z as usize;
        let zs = &mesh.quad_vfr_z[u * kz..u * kz + kz];
        if !(mean_depth > 0.0) {
            return if zs[0] < zs[3] { zs[0] } else { zs[3] };
        }
        return quad_eta_from_mean_depth(zs, mesh.quad_vfr_a[u * 2], mesh.quad_vfr_a[u * 2 + 1], mean_depth, 0.0);
    }
    cell_free_surface_elevation(
        mean_depth,
        mesh.vz[mesh.cell_vertex(t, 0) as usize],
        mesh.vz[mesh.cell_vertex(t, 1) as usize],
        mesh.vz[mesh.cell_vertex(t, 2) as usize],
    )
}

/// cellFreeSurfaceElevation(mean_depth, za, zb, zc): the exact (eps = 0) VFR closure.
pub fn cell_free_surface_elevation(mean_depth: f64, za: f64, zb: f64, zc: f64) -> f64 {
    let (mut z1, mut z2, mut z3) = (za, zb, zc);
    vfr_sort3(&mut z1, &mut z2, &mut z3);
    if !(mean_depth > 0.0) {
        return z1;
    }
    vfr_eta_from_mean_depth(z1, z2, z3, mean_depth, 0.0)
}

/// reconstructVertexRenderDepths(mesh, state, dry_depth, nthreads, eta_scratch).
pub fn reconstruct_vertex_render_depths(mesh: &MeshData, state: &mut SurfaceStateData, dry_depth: f64, _nthreads: i32, eta_scratch: &mut Vec<f64>) {
    let nv = mesh.n_vertices();
    let nt = mesh.n_triangles();
    eta_scratch.resize(nt as usize, 0.0);
    for t in 0..nt {
        let h = state.depth[t as usize];
        if !(h >= dry_depth) {
            continue;
        }
        eta_scratch[t as usize] = cell_free_surface_elevation_of(mesh, t, h);
    }
    for b in 0..nv as usize {
        let start = mesh.vert_stencil_ptr[b] as usize;
        let end = mesh.vert_stencil_ptr[b + 1] as usize;
        let mut vsum = 0.0;
        let mut wsum = 0.0;
        let mut eta_max = 0.0;
        let mut wet = false;
        for k in start..end {
            let t = mesh.vert_stencil_idx[k] as usize;
            let h = state.depth[t];
            if !(h >= dry_depth) {
                continue; // dry skip (NaN-robust)
            }
            let eta = eta_scratch[t];
            if !eta.is_finite() {
                continue; // nodata z must not spread
            }
            // Wetted-contact gate: a cell votes at this vertex only if its water surface
            // reaches the vertex's corner.
            if !(eta > mesh.vz[b]) {
                continue;
            }
            vsum += h * eta;
            wsum += h;
            if !wet || eta > eta_max {
                eta_max = eta;
                wet = true;
            }
        }
        if wsum > 0.0 {
            let mut eta_v = vsum / wsum;
            if eta_v > eta_max {
                eta_v = eta_max;
            }
            let d = eta_v - mesh.vz[b]; // SIGNED depth (VFR)
            state.vert_depth_signed[b] = if d.is_finite() { d } else { 0.0 };
        } else {
            state.vert_depth_signed[b] = 0.0; // no wet incident cell
        }
    }
}
