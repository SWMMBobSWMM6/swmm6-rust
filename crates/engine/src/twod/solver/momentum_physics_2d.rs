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
// The class's members (with their initializers) and the static `rotate` came with the leaf
// numerics (wave 2); evaluate, advance, computeWallDistance, energy and wallDrag with the step
// (wave 4a, 12 November 2026), checked inside the step2d state gate.
//
// Differences in form, never in result: rotate takes the two discharge vectors as `&mut [f64]`
// (upstream `std::vector<double>&`), and its std::cos / std::sin are f64::cos / f64::sin, which
// call glibc (checked by the differential against the library's own compiled rotate).
// std::hypot is f64::hypot (glibc). `std::pow(x, 2)` is `x * x`, as g++ expands it (it folds a
// pow with exponent 2 into the one correctly rounded multiplication). std::max({...}) keeps the
// first largest under `<` (max3; U351). The exceptions evaluate and advance throw
// (std::runtime_error, propagating to SWMMEngine::step's catch) are Err(Terminated) with
// exception "throw std::runtime_error"; advance's Rollback guard is the explicit restore on Err.
// evaluate's output vectors are `&mut Vec<f64>` arguments; advance passes its own scratch
// members (rx_, ry_, sx_, sy_) by taking them out of the struct for the call.
//
// Checked bit for bit against the reference library by gate/diff/twod_momentum.cpp and
// examples/twod_momentum.rs.

use crate::core::cmath::{cclamp, cmax, cmin};
use crate::core::engine_init_hydraulics::Terminated;
use crate::twod::data::boundary_data::BoundaryData;
use crate::twod::data::mesh_data::MeshData;
use crate::twod::solver::inertial_edges::InertialEdges;
use crate::twod::solver::inertial_kernels::face_depth_from_eta;

fn runtime_error(what: &'static str) -> Terminated {
    Terminated { exception: "throw std::runtime_error", what }
}

/// std::max({a, b, c}): the first largest under `<` (U351).
#[inline]
fn max3(a: f64, b: f64, c: f64) -> f64 {
    let mut m = a;
    if m < b {
        m = b;
    }
    if m < c {
        m = c;
    }
    m
}

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

    /// Exact tangential quadratic decay on each selected exterior wall face
    /// (MomentumPhysics2D::wallDrag): the slots in order, then in reverse; returns the kinetic
    /// energy removed (per unit density).
    pub fn wall_drag(m: &MeshData, b: &BoundaryData, slots: &[i32], h: &[f64], eta: &[f64], qx: &mut [f64], qy: &mut [f64], dry: f64, seconds: f64, vfr: bool) -> f64 {
        let mut removed = 0.0;
        let mut apply = |slot: i32, qx: &mut [f64], qy: &mut [f64]| {
            let i = MeshData::slot_cell(slot) as usize;
            let us = slot as usize;
            if h[i] <= dry {
                return;
            }
            let mut hw = h[i];
            if vfr {
                let (mut a, mut z) = (0, 0);
                m.cell_edge_vertices(i as i32, MeshData::slot_local(slot), &mut a, &mut z);
                hw = face_depth_from_eta(eta[i], cmin(m.vz[a as usize], m.vz[z as usize]), cmax(m.vz[a as usize], m.vz[z as usize]));
            }
            let k = b.edge_wall_drag[us] * m.edge_length[us] * hw / (m.tri_area[i] * h[i]);
            if k == 0.0 {
                return;
            }
            let tx = -m.edge_ny[us];
            let ty = m.edge_nx[us];
            let ut = (qx[i] * tx + qy[i] * ty) / h[i];
            if ut == 0.0 {
                return;
            }
            let next = (1.0 / (1.0 / ut.abs() + 0.5 * seconds * k)).copysign(ut);
            let dq = h[i] * (next - ut);
            qx[i] += dq * tx;
            qy[i] += dq * ty;
            removed += 0.5 * m.tri_area[i] * h[i] * (ut - next) * (ut + next);
        };
        for &slot in slots {
            apply(slot, qx, qy);
        }
        for &slot in slots.iter().rev() {
            apply(slot, qx, qy);
        }
        removed
    }

    /// The viscosity of a cell from its velocity gradient (evaluate's shared tail).
    #[inline]
    fn cell_viscosity(&self, m: &MeshData, i: usize, h: f64, fixed_viscosity: f64, background: f64, depth_coefficient: f64, cs: f64) -> f64 {
        if fixed_viscosity >= 0.0 {
            fixed_viscosity
        } else if depth_coefficient > 0.0 {
            let l = cmin(depth_coefficient * h, Self::KAPPA * self.wall_distance[i]);
            background + l * l * self.vorticity[i]
        } else {
            background + cs * cs * m.tri_area[i] * self.strain[i]
        }
    }

    /// Evaluate the conservative shared-face dq/dt into rx, ry (MomentumPhysics2D::evaluate);
    /// returns the explicit step bound 0.4 / max rate (infinity when nothing mixes).
    pub fn evaluate(
        &mut self,
        m: &MeshData,
        edges: &InertialEdges,
        h: &[f64],
        eta: &[f64],
        qx: &[f64],
        qy: &[f64],
        dry: f64,
        cs: f64,
        rx: &mut Vec<f64>,
        ry: &mut Vec<f64>,
        fixed_viscosity: f64,
        background: f64,
        depth_coefficient: f64,
    ) -> Result<f64, Terminated> {
        let n = m.n_cells() as usize;
        rx.clear();
        rx.resize(n, 0.0);
        ry.clear();
        ry.resize(n, 0.0);
        self.viscosity.clear();
        self.viscosity.resize(n, 0.0);
        self.strain.clear();
        self.strain.resize(n, 0.0);
        self.vorticity.clear();
        self.vorticity.resize(n, 0.0);
        if depth_coefficient > 0.0 && self.wall_distance.len() != n {
            self.compute_wall_distance(m);
        }
        self.bound_x_.clear();
        self.bound_x_.resize(n, 0.0);
        self.bound_y_.clear();
        self.bound_y_.resize(n, 0.0);
        self.stencils_.clear();
        self.stencils_.resize(n, Stencil::default());
        self.connected_.clear();
        self.connected_.resize(edges.ne as usize, 0);
        for e in 0..edges.ne as usize {
            let i = edges.cL[e] as usize;
            let j = edges.cR[e] as usize;
            let sill = cmax(eta[i] - h[i], eta[j] - h[j]);
            self.connected_[e] = (edges.conveyance[e] > 0.0 && h[i] > dry && h[j] > dry && cmin(eta[i], eta[j]) - sill > dry) as u8;
        }
        for i in 0..n {
            if h[i] <= dry {
                continue;
            }
            let mut st = Stencil::default();
            st.ids[0] = i as i32;
            st.count = 1;
            let (mut xx, mut xy, mut yy) = (0.0, 0.0, 0.0);
            for p in edges.cell_ptr[i] as usize..edges.cell_ptr[i + 1] as usize {
                let e = edges.cell_edge[p] as usize;
                if self.connected_[e] == 0 {
                    continue;
                }
                let j = if edges.cL[e] == i as i32 { edges.cR[e] } else { edges.cL[e] } as usize;
                let dx = m.tri_cx[j] - m.tri_cx[i];
                let dy = m.tri_cy[j] - m.tri_cy[i];
                let w = 1.0 / (dx * dx + dy * dy);
                let c = st.count as usize;
                st.ids[c] = j as i32;
                st.gx[c] = w * dx;
                st.gy[c] = w * dy;
                st.count += 1;
                xx += w * dx * dx;
                xy += w * dx * dy;
                yy += w * dy * dy;
            }
            let det = xx * yy - xy * xy;
            if st.count < 3 || det <= 1e-10 * (xx + yy) * (xx + yy) {
                st.count = 1;
                self.stencils_[i] = st;
                continue;
            }
            let (mut ux, mut uy, mut vx, mut vy) = (0.0, 0.0, 0.0, 0.0);
            for k in 1..st.count as usize {
                let x = st.gx[k];
                let y = st.gy[k];
                st.gx[k] = (yy * x - xy * y) / det;
                st.gy[k] = (xx * y - xy * x) / det;
                st.gx[0] -= st.gx[k];
                st.gy[0] -= st.gy[k];
            }
            for k in 0..st.count as usize {
                let j = st.ids[k] as usize;
                let u = qx[j] / h[j];
                let v = qy[j] / h[j];
                ux += st.gx[k] * u;
                uy += st.gy[k] * u;
                vx += st.gx[k] * v;
                vy += st.gy[k] * v;
            }
            self.stencils_[i] = st;
            self.strain[i] = (2.0 * ux * ux + 2.0 * vy * vy + (uy + vx) * (uy + vx)).sqrt();
            self.vorticity[i] = (vx - uy).abs();
            self.viscosity[i] = self.cell_viscosity(m, i, h[i], fixed_viscosity, background, depth_coefficient, cs);
            if !self.viscosity[i].is_finite() {
                return Err(runtime_error("Nonfinite Smagorinsky viscosity."));
            }
        }
        // A corner with one wet neighbor: extend a connected neighbor's complete stencil.
        for i in 0..n {
            if h[i] <= dry || self.stencils_[i].count != 1 {
                continue;
            }
            for p in edges.cell_ptr[i] as usize..edges.cell_ptr[i + 1] as usize {
                let e = edges.cell_edge[p] as usize;
                if self.connected_[e] == 0 {
                    continue;
                }
                let j = if edges.cL[e] == i as i32 { edges.cR[e] } else { edges.cL[e] } as usize;
                if self.stencils_[j].count > 1 {
                    self.stencils_[i] = self.stencils_[j];
                    break;
                }
            }
            let st = self.stencils_[i];
            let (mut ux, mut uy, mut vx, mut vy) = (0.0, 0.0, 0.0, 0.0);
            for k in 0..st.count as usize {
                let j = st.ids[k] as usize;
                let u = qx[j] / h[j];
                let v = qy[j] / h[j];
                ux += st.gx[k] * u;
                uy += st.gy[k] * u;
                vx += st.gx[k] * v;
                vy += st.gy[k] * v;
            }
            self.strain[i] = (2.0 * ux * ux + 2.0 * vy * vy + (uy + vx) * (uy + vx)).sqrt();
            self.vorticity[i] = (vx - uy).abs();
            self.viscosity[i] = self.cell_viscosity(m, i, h[i], fixed_viscosity, background, depth_coefficient, cs);
        }
        for e in 0..edges.ne as usize {
            if self.connected_[e] == 0 {
                continue;
            }
            let i = edges.cL[e] as usize;
            let j = edges.cR[e] as usize;
            let nx = edges.nx[e];
            let ny = edges.ny[e];
            let dx = m.tri_cx[j] - m.tri_cx[i];
            let dy = m.tri_cy[j] - m.tri_cy[i];
            let dn = dx * nx + dy * ny;
            if !(dn > 0.0) {
                return Err(runtime_error("Invalid mixing face geometry."));
            }
            let mut ids = [0i32; 12];
            let mut gx = [0.0f64; 12];
            let mut gy = [0.0f64; 12];
            let mut count = 0usize;
            let mut add = |id: i32, x: f64, y: f64, ids: &mut [i32; 12], gx: &mut [f64; 12], gy: &mut [f64; 12]| {
                let mut k = 0usize;
                while k < count && ids[k] != id {
                    k += 1;
                }
                if k == count {
                    ids[k] = id;
                    count += 1;
                }
                gx[k] += x;
                gy[k] += y;
            };
            let fx = edges.mx[e];
            let fy = edges.my[e];
            let di = (fx - m.tri_cx[i]).hypot(fy - m.tri_cy[i]); // LIBM (glibc hypot)
            let dj = (fx - m.tri_cx[j]).hypot(fy - m.tri_cy[j]); // LIBM (glibc hypot)
            let wi = dj / (di + dj);
            let wj = 1.0 - wi;
            for cell in [i, j] {
                let st = self.stencils_[cell];
                let w = if cell == i { wi } else { wj };
                for k in 0..st.count as usize {
                    add(st.ids[k], w * st.gx[k], w * st.gy[k], &mut ids, &mut gx, &mut gy);
                }
            }
            add(i as i32, 0.0, 0.0, &mut ids, &mut gx, &mut gy);
            add(j as i32, 0.0, 0.0, &mut ids, &mut gx, &mut gy);
            let (mut ux, mut uy, mut vx, mut vy, mut bx, mut by) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            let nu = if fixed_viscosity >= 0.0 { fixed_viscosity } else { wi * self.viscosity[i] + wj * self.viscosity[j] };
            let wet = cmax(0.0, cmin(eta[i], eta[j]) - cmax(eta[i] - h[i], eta[j] - h[j]));
            let weight = edges.xi[e] * wet * nu * edges.conveyance[e];
            for k in 0..count {
                let id = ids[k];
                let uid = id as usize;
                let delta = (if id == j as i32 { 1.0 } else { 0.0 }) - (if id == i as i32 { 1.0 } else { 0.0 });
                let correction = (delta - gx[k] * dx - gy[k] * dy) / dn;
                gx[k] += nx * correction;
                gy[k] += ny * correction;
                ux += gx[k] * qx[uid] / h[uid];
                uy += gy[k] * qx[uid] / h[uid];
                vx += gx[k] * qy[uid] / h[uid];
                vy += gy[k] * qy[uid] / h[uid];
                bx += (2.0 * nx * gx[k] + ny * gy[k]).abs() + (ny * gx[k]).abs();
                by += (nx * gy[k]).abs() + (nx * gx[k] + 2.0 * ny * gy[k]).abs();
            }
            let tx = weight * (2.0 * ux * nx + (uy + vx) * ny);
            let ty = weight * (2.0 * vy * ny + (uy + vx) * nx);
            rx[i] += tx / m.tri_area[i];
            ry[i] += ty / m.tri_area[i];
            rx[j] -= tx / m.tri_area[j];
            ry[j] -= ty / m.tri_area[j];
            for cell in [i, j] {
                let mass = m.tri_area[cell] * h[cell];
                self.bound_x_[cell] += weight * bx / mass;
                self.bound_y_[cell] += weight * by / mass;
            }
        }
        let mut rate = 0.0;
        for i in 0..n {
            rate = max3(rate, self.bound_x_[i], self.bound_y_[i]);
        }
        if !rate.is_finite() {
            return Err(runtime_error("Nonfinite mixing timestep bound."));
        }
        self.last_limit = if rate > 0.0 { 0.4 / rate } else { f64::INFINITY };
        Ok(self.last_limit)
    }

    /// SSP-RK2 mixing at fixed water storage (MomentumPhysics2D::advance): substeps bounded by
    /// evaluate's limit, halved until the kinetic energy does not grow; on failure qx, qy are
    /// restored (the Rollback guard).
    pub fn advance(
        &mut self,
        m: &MeshData,
        edges: &InertialEdges,
        h: &[f64],
        eta: &[f64],
        qx: &mut Vec<f64>,
        qy: &mut Vec<f64>,
        dry: f64,
        cs: f64,
        seconds: f64,
        background: f64,
        depth_coefficient: f64,
    ) -> Result<(), Terminated> {
        if seconds <= 0.0 || (cs == 0.0 && background == 0.0 && depth_coefficient == 0.0) {
            return Ok(());
        }
        let old_x = qx.clone();
        let old_y = qy.clone();
        let r = self.advance_inner(m, edges, h, eta, qx, qy, dry, cs, seconds, background, depth_coefficient);
        if r.is_err() {
            *qx = old_x;
            *qy = old_y;
        }
        r
    }

    fn advance_inner(
        &mut self,
        m: &MeshData,
        edges: &InertialEdges,
        h: &[f64],
        eta: &[f64],
        qx: &mut Vec<f64>,
        qy: &mut Vec<f64>,
        dry: f64,
        cs: f64,
        seconds: f64,
        background: f64,
        depth_coefficient: f64,
    ) -> Result<(), Terminated> {
        let mut elapsed = 0.0;
        while elapsed < seconds {
            let mut rx = std::mem::take(&mut self.rx_);
            let mut ry = std::mem::take(&mut self.ry_);
            let bound = self.evaluate(m, edges, h, eta, qx, qy, dry, cs, &mut rx, &mut ry, -1.0, background, depth_coefficient);
            self.rx_ = rx;
            self.ry_ = ry;
            let bound = bound?;
            let mut dt = cmin(seconds - elapsed, bound);
            let e0 = Self::energy(m, h, qx, qy, dry);
            let mut accepted = false;
            for _retry in 0..40 {
                self.tx_.clone_from(qx);
                self.ty_.clone_from(qy);
                for j in 0..self.tx_.len() {
                    self.tx_[j] += dt * self.rx_[j];
                    self.ty_[j] += dt * self.ry_[j];
                }
                let tx = std::mem::take(&mut self.tx_);
                let ty = std::mem::take(&mut self.ty_);
                let mut sx = std::mem::take(&mut self.sx_);
                let mut sy = std::mem::take(&mut self.sy_);
                let stage_bound = self.evaluate(m, edges, h, eta, &tx, &ty, dry, cs, &mut sx, &mut sy, -1.0, background, depth_coefficient);
                self.tx_ = tx;
                self.ty_ = ty;
                self.sx_ = sx;
                self.sy_ = sy;
                let stage_bound = stage_bound?;
                if dt <= stage_bound {
                    for j in 0..self.tx_.len() {
                        self.tx_[j] = 0.5 * (qx[j] + self.tx_[j] + dt * self.sx_[j]);
                        self.ty_[j] = 0.5 * (qy[j] + self.ty_[j] + dt * self.sy_[j]);
                    }
                    let e1 = Self::energy(m, h, &self.tx_, &self.ty_, dry);
                    if e1.is_finite() && e1 <= e0 + 1e-13 * cmax(e0, 1e-30) {
                        accepted = true;
                        break;
                    }
                }
                dt *= 0.5;
                self.retries += 1;
            }
            if !accepted || elapsed + dt == elapsed {
                return Err(runtime_error("Smagorinsky mixing cannot find an admissible timestep."));
            }
            std::mem::swap(qx, &mut self.tx_);
            std::mem::swap(qy, &mut self.ty_);
            elapsed = if dt == seconds - elapsed { seconds } else { elapsed + dt };
            self.substeps += 1;
        }
        let mut rx = std::mem::take(&mut self.rx_);
        let mut ry = std::mem::take(&mut self.ry_);
        let r = self.evaluate(m, edges, h, eta, qx, qy, dry, cs, &mut rx, &mut ry, -1.0, background, depth_coefficient);
        self.rx_ = rx;
        self.ry_ = ry;
        r.map(|_| ())
    }

    /// Nearest-boundary-edge distance for every cell centroid (computeWallDistance).
    pub fn compute_wall_distance(&mut self, m: &MeshData) {
        let mut segments: Vec<[f64; 4]> = Vec::new();
        for c in 0..m.n_cells() {
            for k in 0..m.cell_vertex_count(c) {
                if m.cell_neighbour(c, k) < 0 {
                    let (mut a, mut b) = (0, 0);
                    m.cell_edge_vertices(c, k, &mut a, &mut b);
                    segments.push([m.vx[a as usize], m.vy[a as usize], m.vx[b as usize], m.vy[b as usize]]);
                }
            }
        }
        self.wall_distance.clear();
        self.wall_distance.resize(m.n_cells() as usize, f64::INFINITY);
        for c in 0..m.n_cells() as usize {
            let px = m.tri_cx[c];
            let py = m.tri_cy[c];
            let mut best = f64::INFINITY;
            for s in &segments {
                let dx = s[2] - s[0];
                let dy = s[3] - s[1];
                let len2 = dx * dx + dy * dy;
                let t = if len2 > 0.0 { cclamp(((px - s[0]) * dx + (py - s[1]) * dy) / len2, 0.0, 1.0) } else { 0.0 };
                best = cmin(best, (px - s[0] - t * dx).hypot(py - s[1] - t * dy)); // LIBM (glibc hypot)
            }
            self.wall_distance[c] = best;
        }
    }

    /// The wet cells' kinetic energy per unit density (MomentumPhysics2D::energy).
    pub fn energy(m: &MeshData, h: &[f64], qx: &[f64], qy: &[f64], dry: f64) -> f64 {
        let mut sum = 0.0;
        for i in 0..m.n_cells() as usize {
            if h[i] > dry {
                sum += 0.5 * m.tri_area[i] * (qx[i] * qx[i] + qy[i] * qy[i]) / h[i];
            }
        }
        sum
    }
}
