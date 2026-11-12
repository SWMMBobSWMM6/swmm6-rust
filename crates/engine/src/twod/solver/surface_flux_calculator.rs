// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/solver/SurfaceFluxCalculator.hpp and
// SurfaceFluxCalculator.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// So far: the header's two inline scalar sink ramps (evapSink and infilSink, the depth-limited
// evaporation and infiltration rates with the cubic Hermite wet/dry ramp; wave 2, checked bit
// for bit against upstream's header by gate/diff/twod_sinks.cpp and examples/twod_sinks.rs),
// and the diagnostics SurfaceRouter2D::initialize runs (wave 3): computeUnlimitedGradients,
// computeLimitedGradients and computeFaceVelocity, with the file's helpers tri_nbr, sq and
// edgeEndpointZ. The boundary-edge fluxes and computeCellContinuity come with the step (wave 4).
// Upstream's OpenMP loops are `schedule(static)` and each iteration writes only its own cell,
// so the serial loops give the same bits at any thread count.

use crate::core::cmath::cmax;
use crate::twod::data::boundary_data::{BoundaryData, BoundaryType};
use crate::twod::data::mesh_data::{MeshData, K_MAX_CELL_VERTS};
use crate::twod::data::solver_options_2d::{FaceDepth2D, Momentum2D, SolverOptions2D};
use crate::twod::data::surface_state_data::SurfaceStateData;
use crate::twod::solver::inertial_kernels::{face_depth_from_eta, face_flow_depth};

/// Depth-limited evaporation sink rate, m/s (evapSink).
#[inline]
pub fn evap_sink(rate: f64, depth: f64, dry_depth: f64) -> f64 {
    if rate <= 0.0 || depth <= 0.0 {
        return 0.0;
    }
    if depth >= dry_depth {
        return rate;
    }
    let t = depth / dry_depth;
    rate * t * t * (3.0 - 2.0 * t)
}

/// Depth-limited infiltration sink rate, m/s (infilSink); the same ramp as evap_sink.
#[inline]
pub fn infil_sink(rate: f64, depth: f64, dry_depth: f64) -> f64 {
    if rate <= 0.0 || depth <= 0.0 {
        return 0.0;
    }
    if depth >= dry_depth {
        return rate;
    }
    let t = depth / dry_depth;
    rate * t * t * (3.0 - 2.0 * t)
}

#[inline]
fn tri_nbr(mesh: &MeshData, t: i32, e: i32) -> i32 {
    if e >= 0 && e < mesh.cell_vertex_count(t) {
        mesh.cell_neighbour(t, e)
    } else {
        -1
    }
}

#[inline]
fn sq(x: f64) -> f64 {
    x * x
}

/// edgeEndpointZ(mesh, t, e, z_lo, z_hi): the endpoint bed elevations of local edge e of cell
/// t, sorted, returned as (z_lo, z_hi).
#[inline]
pub fn edge_endpoint_z(mesh: &MeshData, t: i32, e: i32) -> (f64, f64) {
    let (mut va, mut vb) = (0, 0);
    mesh.cell_edge_vertices(t, e, &mut va, &mut vb);
    let za = mesh.vz[va as usize];
    let zb = mesh.vz[vb as usize];
    let z_lo = if za < zb { za } else { zb };
    let z_hi = if za < zb { zb } else { za };
    (z_lo, z_hi)
}

/// computeUnlimitedGradients(mesh, state, nthreads): Green-Gauss head gradients.
pub fn compute_unlimited_gradients(mesh: &MeshData, state: &mut SurfaceStateData, _nthreads: i32) {
    let nt = mesh.n_triangles();
    for i in 0..nt {
        let ui = i as usize;
        let inv_area = if mesh.tri_area[ui] > 1.0e-30 { 1.0 / mesh.tri_area[ui] } else { 0.0 };
        let mut gx = 0.0;
        let mut gy = 0.0;
        let nv_i = mesh.cell_vertex_count(i);
        for e in 0..nv_i {
            let idx = MeshData::slot(i, e) as usize;
            let nbr = tri_nbr(mesh, i, e);
            // Head at edge midpoint: average of this cell and neighbour
            let h_edge = if nbr >= 0 {
                0.5 * (state.head[ui] + state.head[nbr as usize])
            } else {
                // Boundary: this cell's head (zero-gradient extrapolation)
                state.head[ui]
            };
            gx += h_edge * mesh.edge_nx[idx] * mesh.edge_length[idx];
            gy += h_edge * mesh.edge_ny[idx] * mesh.edge_length[idx];
        }
        state.grad_hx[ui] = gx * inv_area;
        state.grad_hy[ui] = gy * inv_area;
    }
}

/// computeLimitedGradients(mesh, state, epsilon, nthreads): Jawahar-Kamath weighted gradients.
pub fn compute_limited_gradients(mesh: &MeshData, state: &mut SurfaceStateData, epsilon: f64, _nthreads: i32) {
    let nt = mesh.n_triangles();
    let eps2 = epsilon * epsilon;
    for i in 0..nt {
        let ui = i as usize;
        let q0 = sq(state.grad_hx[ui]) + sq(state.grad_hy[ui]) + eps2;
        let mut q = [0.0f64; K_MAX_CELL_VERTS as usize];
        let mut gx_nbr = [0.0f64; K_MAX_CELL_VERTS as usize];
        let mut gy_nbr = [0.0f64; K_MAX_CELL_VERTS as usize];
        let nv_i = mesh.cell_vertex_count(i);
        for e in 0..nv_i as usize {
            let nbr = tri_nbr(mesh, i, e as i32);
            if nbr >= 0 {
                let un = nbr as usize;
                q[e] = sq(state.grad_hx[un]) + sq(state.grad_hy[un]) + eps2;
                gx_nbr[e] = state.grad_hx[un];
                gy_nbr[e] = state.grad_hy[un];
            } else {
                // Boundary: mirror the cell's own gradient.
                q[e] = q0;
                gx_nbr[e] = state.grad_hx[ui];
                gy_nbr[e] = state.grad_hy[ui];
            }
        }
        if nv_i == 4 {
            // Quad: the same JK product weights over 5 contributing gradients.
            let qq = [q0, q[0], q[1], q[2], q[3]];
            let gx = [state.grad_hx[ui], gx_nbr[0], gx_nbr[1], gx_nbr[2], gx_nbr[3]];
            let gy = [state.grad_hy[ui], gy_nbr[0], gy_nbr[1], gy_nbr[2], gy_nbr[3]];
            let mut num = [0.0f64; 5];
            let mut den = 0.0;
            for k in 0..5 {
                let mut pr = 1.0;
                for j in 0..5 {
                    if j != k {
                        pr *= qq[j];
                    }
                }
                num[k] = pr;
                den += pr;
            }
            let mut sx = 0.0;
            let mut sy = 0.0;
            for k in 0..5 {
                sx += num[k] / den * gx[k];
                sy += num[k] / den * gy[k];
            }
            state.grad_hx_lim[ui] = sx;
            state.grad_hy_lim[ui] = sy;
            continue;
        }
        let n0 = q[0] * q[1] * q[2]; // skip self
        let n1 = q0 * q[1] * q[2]; // skip nbr 0
        let n2 = q0 * q[0] * q[2]; // skip nbr 1
        let n3 = q0 * q[0] * q[1]; // skip nbr 2
        let denom = n0 + n1 + n2 + n3;
        let w0 = n0 / denom;
        let w1 = n1 / denom;
        let w2 = n2 / denom;
        let w3 = n3 / denom;
        state.grad_hx_lim[ui] = w0 * state.grad_hx[ui] + w1 * gx_nbr[0] + w2 * gx_nbr[1] + w3 * gx_nbr[2];
        state.grad_hy_lim[ui] = w0 * state.grad_hy[ui] + w1 * gy_nbr[0] + w2 * gy_nbr[1] + w3 * gy_nbr[2];
    }
}

/// computeFaceVelocity(mesh, state, opts): the output cell velocity, a least-squares fit of the
/// edge flows. `boundary` is the router's BoundaryData when state.boundary is attached.
pub fn compute_face_velocity(mesh: &MeshData, state: &mut SurfaceStateData, boundary: Option<&BoundaryData>, opts: &SolverOptions2D) {
    let nc = mesh.n_cells();
    let missing = f64::NAN; // std::numeric_limits<double>::quiet_NaN()
    for i in 0..nc {
        let ui = i as usize;
        let depth = state.depth[ui];
        state.face_vy[ui] = missing;
        state.face_vx[ui] = missing;
        if !depth.is_finite() {
            continue;
        }
        if depth <= 0.0 || depth < opts.dry_depth {
            state.face_vy[ui] = 0.0;
            state.face_vx[ui] = 0.0;
            continue;
        }
        let mut a00 = 0.0;
        let mut a01 = 0.0;
        let mut a11 = 0.0;
        let mut b0 = 0.0;
        let mut b1 = 0.0;
        let mut known = true;
        for e in 0..mesh.cell_vertex_count(i) {
            let idx = MeshData::slot(i, e) as usize;
            let j = mesh.cell_neighbour(i, e);
            let flow = state.edge_flux[idx];
            let length = mesh.edge_length[idx];
            if !flow.is_finite() {
                known = false;
                break;
            }
            let mut eta = state.head[ui];
            let mut hf = depth;
            if j >= 0 {
                eta = cmax(eta, state.head[j as usize]);
                hf = face_flow_depth(eta, eta, mesh.edge_mz[idx]);
            } else if let Some(b) = boundary.filter(|b| b.edge_bc_type[idx] == BoundaryType::SpecifiedStage as i8) {
                let stage = b.edge_bc_head[idx];
                eta = cmax(eta, stage);
                hf = if opts.momentum == Momentum2D::DiffusiveWave {
                    if state.head[ui] > stage {
                        depth
                    } else {
                        cmax(stage - mesh.tri_cz[ui], 0.0)
                    }
                } else {
                    face_flow_depth(eta, eta, mesh.tri_cz[ui])
                };
            }
            if opts.face_reconstruction == FaceDepth2D::VfrFace {
                let (zlo, zhi) = edge_endpoint_z(mesh, i, e);
                hf = face_depth_from_eta(eta, zlo, zhi);
            }
            let porosity = if j >= 0 { mesh.edge_conveyance[idx] } else { 1.0 };
            let area = length * hf * porosity;
            if !area.is_finite() {
                known = false;
                break;
            }
            if length <= 1e-12 || hf <= opts.dry_depth || area <= 0.0 {
                if flow != 0.0 {
                    known = false;
                    break;
                }
                continue;
            }
            let nx = mesh.edge_nx[idx];
            let ny = mesh.edge_ny[idx];
            a00 += area * nx * nx;
            a01 += area * nx * ny;
            a11 += area * ny * ny;
            b0 += nx * flow;
            b1 += ny * flow;
        }
        let det = a00 * a11 - a01 * a01;
        let trace = a00 + a11;
        if !known || !(det > 1e-10 * trace * trace) {
            continue;
        }
        let vx = (a11 * b0 - a01 * b1) / det;
        let vy = (a00 * b1 - a01 * b0) / det;
        if vx.is_finite() && vy.is_finite() {
            state.face_vx[ui] = vx;
            state.face_vy[ui] = vy;
        }
    }
}
