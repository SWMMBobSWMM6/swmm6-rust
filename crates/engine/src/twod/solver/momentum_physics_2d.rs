// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/solver/MomentumPhysics2D.hpp and
// MomentumPhysics2D.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The 2D momentum physics beyond the local-inertial face law: conservative Smagorinsky /
// mixing-length eddy-viscosity stresses, wall drag and the Coriolis rotation (upstream class
// openswmm::twoD::MomentumPhysics2D).
//
// So far: the class's members (with their initializers) and the static `rotate`, the only
// member that reads no mesh. NOT TRANSLATED yet, because each reads the mesh topology and
// geometry (MeshData, a wave-1 data struct), the solver's face list (InertialEdges, wave 3) or
// the boundary data (BoundaryData, wave 1): evaluate, advance, computeWallDistance, energy and
// wallDrag. They come with the step (wave 4); the reference library exports all five, so their
// differential can call upstream's own compiled code on a mesh built by upstream's MeshBuilder.
//
// Differences in form, never in result: rotate takes the two discharge vectors as `&mut [f64]`
// (upstream `std::vector<double>&`), and its std::cos / std::sin are f64::cos / f64::sin, which
// call glibc (checked by the differential against the library's own compiled rotate).
//
// Checked bit for bit against the reference library by gate/diff/twod_momentum.cpp and
// examples/twod_momentum.rs.

/// The private least-squares gradient stencil of one cell (MomentumPhysics2D::Stencil).
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Stencil {
    pub ids: [i32; 5],
    pub gx: [f64; 5],
    pub gy: [f64; 5],
    pub count: i32,
}

/// Conservative depth-averaged Smagorinsky face stresses, wall drag and Coriolis rotation
/// (MomentumPhysics2D). Field order and initializers are upstream's.
#[derive(Debug, Clone)]
pub struct MomentumPhysics2D {
    /// Per cell: m2/s.
    pub viscosity: Vec<f64>,
    /// Per cell: |S|, 1/s.
    pub strain: Vec<f64>,
    /// Per cell: |omega|, 1/s.
    pub vorticity: Vec<f64>,
    /// Distance from each centroid to the nearest mesh boundary edge (m); built on first use.
    pub wall_distance: Vec<f64>,
    pub last_limit: f64,
    pub substeps: u64,
    pub retries: u64,
    // Private upstream (scratch, not checkpoint state); public here for the gate dumpers.
    pub stencils_: Vec<Stencil>,
    pub connected_: Vec<u8>,
    pub bound_x_: Vec<f64>,
    pub bound_y_: Vec<f64>,
    pub rx_: Vec<f64>,
    pub ry_: Vec<f64>,
    pub sx_: Vec<f64>,
    pub sy_: Vec<f64>,
    pub tx_: Vec<f64>,
    pub ty_: Vec<f64>,
}

impl Default for MomentumPhysics2D {
    fn default() -> Self {
        Self {
            viscosity: Vec::new(),
            strain: Vec::new(),
            vorticity: Vec::new(),
            wall_distance: Vec::new(),
            last_limit: f64::INFINITY, // std::numeric_limits<double>::infinity()
            substeps: 0,
            retries: 0,
            stencils_: Vec::new(),
            connected_: Vec::new(),
            bound_x_: Vec::new(),
            bound_y_: Vec::new(),
            rx_: Vec::new(),
            ry_: Vec::new(),
            sx_: Vec::new(),
            sy_: Vec::new(),
            tx_: Vec::new(),
            ty_: Vec::new(),
        }
    }
}

impl MomentumPhysics2D {
    /// von Karman's constant, the wall cap of the depth mixing length (kappa).
    pub const KAPPA: f64 = 0.41;

    /// Rotate every cell discharge (qx, qy) by `angle` radians, clockwise in (x, y)
    /// (MomentumPhysics2D::rotate): q' = (c qx + s qy, c qy - s qx). No-op for a zero angle.
    pub fn rotate(qx: &mut [f64], qy: &mut [f64], angle: f64) {
        if angle == 0.0 {
            return;
        }
        let c = angle.cos(); // LIBM
        let s = angle.sin(); // LIBM
        for i in 0..qx.len() {
            let x = qx[i];
            let y = qy[i];
            qx[i] = c * x + s * y;
            qy[i] = c * y - s * x;
        }
    }
}
