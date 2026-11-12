// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/mesh/RainfallInterpolator.hpp and
// RainfallInterpolator.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Static per-cell rain-gage weights (natural-neighbour / Laplace inside the gage hull,
// inverse distance outside it, nearest gage on request), built once at initialize; apply()
// (the per-step sparse weighted sum) came with the step (wave 4a, 12 November 2026).
//
// Differences in form, never in result:
//   * build's two exceptions (std::invalid_argument, std::runtime_error, caught by
//     SWMMEngine::initHydraulics) are Err(what).
//   * laplaceWeights' two unordered_map<int, Pt> are HashMaps: they are only inserted into
//     (the last insert for a key wins, as operator[] assignment does) and searched.

#![allow(clippy::all, non_snake_case)]

use std::collections::HashMap;

use crate::core::cmath::{cmax, cmin};
use crate::core::options::Bytes;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Method {
    NaturalNeighbour,
    NearestNeighbour,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i8)]
pub enum CellMethod {
    #[default]
    NaturalNeighbour = 0,
    InverseDistance = 1,
    Nearest = 2,
}

#[derive(Debug, Clone, Default)]
pub struct Diagnostics {
    pub unlocated: i32,
    pub invalid: i32,
    pub duplicates: i32,
    pub idwCells: i32,
}

#[derive(Debug, Clone, Default)]
pub struct RainfallInterpolator {
    pub diagnostics_: Diagnostics,
    pub cell_method_: Vec<CellMethod>,
    pub ready_: bool,
    pub nt_: i32,
    pub w_ptr_: Vec<i32>,
    pub w_gage_: Vec<i32>,
    pub w_val_: Vec<f64>,
}

#[derive(Debug, Clone, Copy)]
struct Pt {
    x: f64,
    y: f64,
}

#[derive(Debug, Clone, Copy)]
struct Tri {
    a: usize,
    b: usize,
    c: usize,
}

#[inline]
fn orient2d(a: &Pt, b: &Pt, c: &Pt) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
}

#[inline]
fn dist2(p: &Pt, q: &Pt) -> f64 {
    let dx = p.x - q.x;
    let dy = p.y - q.y;
    dx * dx + dy * dy
}

#[inline]
fn circumcenter(a: &Pt, b: &Pt, c: &Pt) -> Pt {
    let bx = b.x - a.x;
    let by = b.y - a.y;
    let cx = c.x - a.x;
    let cy = c.y - a.y;
    let d = 2.0 * (bx * cy - by * cx);
    let b2 = bx * bx + by * by;
    let c2 = cx * cx + cy * cy;
    Pt { x: a.x + (b2 * cy - c2 * by) / d, y: a.y + (bx * c2 - cx * b2) / d }
}

#[inline]
fn in_circle(a: &Pt, b: &Pt, c: &Pt, p: &Pt) -> bool {
    let ax = a.x - p.x;
    let ay = a.y - p.y;
    let bx = b.x - p.x;
    let by = b.y - p.y;
    let cx = c.x - p.x;
    let cy = c.y - p.y;
    let det = (ax * ax + ay * ay) * (bx * cy - cx * by) - (bx * bx + by * by) * (ax * cy - cx * ay) + (cx * cx + cy * cy) * (ax * by - bx * ay);
    det > 0.0
}

fn cavity_boundary(bad_edges: &[[usize; 2]]) -> Vec<[usize; 2]> {
    let mut boundary = Vec::new();
    for e in bad_edges {
        let mut shared = false;
        for f in bad_edges {
            if f[0] == e[1] && f[1] == e[0] {
                shared = true;
                break;
            }
        }
        if !shared {
            boundary.push(*e);
        }
    }
    boundary
}

/// Bowyer-Watson Delaunay triangulation of the sites (with a super-triangle).
fn delaunay(site: &[Pt]) -> Vec<Tri> {
    let m = site.len();
    if m < 3 {
        return Vec::new();
    }
    let mut minx = site[0].x;
    let mut maxx = site[0].x;
    let mut miny = site[0].y;
    let mut maxy = site[0].y;
    for s in site {
        minx = cmin(minx, s.x);
        maxx = cmax(maxx, s.x);
        miny = cmin(miny, s.y);
        maxy = cmax(maxy, s.y);
    }
    let dmax = cmax(maxx - minx, maxy - miny);
    if dmax <= 0.0 {
        return Vec::new(); // all coincident
    }
    let mut p: Vec<Pt> = site.to_vec();
    let midx = 0.5 * (minx + maxx);
    let midy = 0.5 * (miny + maxy);
    let r = 20.0 * dmax;
    p.push(Pt { x: midx - 2.0 * r, y: midy - r });
    p.push(Pt { x: midx + 2.0 * r, y: midy - r });
    p.push(Pt { x: midx, y: midy + 2.0 * r });
    let (s0, s1, s2) = (m, m + 1, m + 2);
    let ccw = |p: &Vec<Pt>, a: usize, b: usize, c: usize| -> Tri {
        let (mut b, mut c) = (b, c);
        if orient2d(&p[a], &p[b], &p[c]) < 0.0 {
            std::mem::swap(&mut b, &mut c);
        }
        Tri { a, b, c }
    };
    let mut tris: Vec<Tri> = vec![ccw(&p, s0, s1, s2)];
    for i in 0..m {
        let pt = p[i];
        let mut bad_edges: Vec<[usize; 2]> = Vec::new();
        let mut good: Vec<Tri> = Vec::with_capacity(tris.len());
        for t in &tris {
            if in_circle(&p[t.a], &p[t.b], &p[t.c], &pt) {
                bad_edges.push([t.a, t.b]);
                bad_edges.push([t.b, t.c]);
                bad_edges.push([t.c, t.a]);
            } else {
                good.push(*t);
            }
        }
        tris = good;
        for e in cavity_boundary(&bad_edges) {
            tris.push(ccw(&p, e[0], e[1], i));
        }
    }
    // Drop triangles touching a super vertex.
    let mut result = Vec::with_capacity(tris.len());
    for t in tris {
        if t.a < m && t.b < m && t.c < m {
            result.push(t);
        }
    }
    result
}

type Weights = Vec<(i32, f64)>; // (global gage idx, weight)

fn idw_all(p: &Pt, site: &[Pt], gid: &[i32]) -> Weights {
    let m = site.len();
    for j in 0..m {
        if dist2(p, &site[j]) < 1.0e-18 {
            return vec![(gid[j], 1.0)];
        }
    }
    let mut min_d2 = f64::INFINITY;
    for s in site {
        min_d2 = cmin(min_d2, dist2(p, s));
    }
    let mut w: Weights = Vec::with_capacity(m);
    let mut sum = 0.0;
    for j in 0..m {
        let wj = min_d2 / dist2(p, &site[j]); // rescaled 1/d^2; avoids overflow
        w.push((gid[j], wj));
        sum += wj;
    }
    for e in w.iter_mut() {
        e.1 /= sum;
    }
    w
}

fn laplace_weights(p: &Pt, site: &[Pt], gid: &[i32], tris: &[Tri]) -> Weights {
    let m = site.len();
    for j in 0..m {
        if dist2(p, &site[j]) < 1.0e-18 {
            return vec![(gid[j], 1.0)];
        }
    }
    // Inside-hull test: p must lie in some Delaunay triangle.
    let mut inside = false;
    for t in tris {
        let d0 = orient2d(&site[t.a], &site[t.b], p);
        let d1 = orient2d(&site[t.b], &site[t.c], p);
        let d2 = orient2d(&site[t.c], &site[t.a], p);
        if !((d0 < 0.0 || d1 < 0.0 || d2 < 0.0) && (d0 > 0.0 || d1 > 0.0 || d2 > 0.0)) {
            inside = true;
            break;
        }
    }
    if !inside {
        return Vec::new();
    }
    // Insertion cavity: triangles whose circumcircle contains p.
    let mut bad_edges: Vec<[usize; 2]> = Vec::new();
    for t in tris {
        if in_circle(&site[t.a], &site[t.b], &site[t.c], p) {
            bad_edges.push([t.a, t.b]);
            bad_edges.push([t.b, t.c]);
            bad_edges.push([t.c, t.a]);
        }
    }
    let boundary = cavity_boundary(&bad_edges);
    if boundary.is_empty() {
        return Vec::new();
    }
    // Circumcenter of each fan triangle (p, u, v); index by tail u and head v.
    let mut c_by_tail: HashMap<usize, Pt> = HashMap::new();
    let mut c_by_head: HashMap<usize, Pt> = HashMap::new();
    for e in &boundary {
        let u = &site[e[0]];
        let v = &site[e[1]];
        if orient2d(p, u, v).abs() < 1.0e-20 {
            return Vec::new(); // collinear -> IDW
        }
        let cc = circumcenter(p, u, v);
        c_by_tail.insert(e[0], cc);
        c_by_head.insert(e[1], cc);
    }
    let mut w: Weights = Vec::new();
    let mut sum = 0.0;
    for e in &boundary {
        let s = e[1]; // each boundary site is the head of exactly one edge
        let in_edge = c_by_head.get(&s); // fan triangle (p, prev, s)
        let out_edge = c_by_tail.get(&s); // fan triangle (p, s, next)
        let (Some(ie), Some(oe)) = (in_edge, out_edge) else {
            return Vec::new();
        };
        let facet = dist2(ie, oe).sqrt();
        let wj = facet / dist2(p, &site[s]).sqrt();
        w.push((gid[s], wj));
        sum += wj;
    }
    if !(sum > 0.0) || !sum.is_finite() {
        return Vec::new();
    }
    for e in w.iter_mut() {
        e.1 /= sum;
    }
    w
}

impl RainfallInterpolator {
    pub fn diagnostics(&self) -> &Diagnostics {
        &self.diagnostics_
    }
    pub fn ready(&self) -> bool {
        self.ready_
    }

    /// build(cx, cy, gage_x, gage_y, n_gages, gage_scale, method). Err(what) where upstream
    /// throws.
    pub fn build(&mut self, cx: &[f64], cy: &[f64], gage_x: &[f64], gage_y: &[f64], n_gages: i32, gage_scale: f64, method: Method) -> Result<(), Bytes> {
        self.ready_ = false;
        self.diagnostics_ = Diagnostics::default();
        if cx.len() != cy.len() || !gage_scale.is_finite() || gage_scale <= 0.0 {
            return Err(b"2D rainfall: invalid coordinate dimensions or scale".to_vec());
        }
        for i in 0..cx.len() {
            if !cx[i].is_finite() || !cy[i].is_finite() {
                return Err(format!("2D rainfall: nonfinite cell centroid {}", i).into_bytes());
            }
        }
        self.nt_ = cx.len() as i32;
        self.w_ptr_.clear();
        self.w_ptr_.resize(self.nt_ as usize + 1, 0);
        self.w_gage_.clear();
        self.w_val_.clear();
        self.cell_method_.clear();
        self.cell_method_.resize(self.nt_ as usize, CellMethod::Nearest);

        // Collect located gages (skip the un-located (0,0) sentinel and duplicates).
        let mut site: Vec<Pt> = Vec::new();
        let mut gid: Vec<i32> = Vec::new();
        let gx_n = gage_x.len() as i32;
        let gy_n = gage_y.len() as i32;
        for g in 0..n_gages {
            let gx = if g < gx_n { gage_x[g as usize] } else { 0.0 };
            let gy = if g < gy_n { gage_y[g as usize] } else { 0.0 };
            if gx == 0.0 && gy == 0.0 {
                self.diagnostics_.unlocated += 1;
                continue;
            }
            let s = Pt { x: gx * gage_scale, y: gy * gage_scale };
            if !s.x.is_finite() || !s.y.is_finite() {
                self.diagnostics_.invalid += 1;
                continue;
            }
            let mut dup = false;
            for e in &site {
                if dist2(e, &s) < 1.0e-12 {
                    dup = true;
                    break;
                }
            }
            if dup {
                self.diagnostics_.duplicates += 1;
                continue;
            }
            site.push(s);
            gid.push(g);
        }
        let m = site.len();
        if m == 0 || self.nt_ == 0 {
            return Ok(()); // ready_ stays false: the caller uses the SYSTEM mean
        }

        // One dimensionless local frame for both sites and cells.
        let origin = site[0];
        let mut scale = 0.0;
        for s in &site {
            scale = cmax(scale, cmax((s.x - origin.x).abs(), (s.y - origin.y).abs()));
        }
        if !(scale > 0.0) {
            scale = 1.0;
        }
        for s in site.iter_mut() {
            s.x = (s.x - origin.x) / scale;
            s.y = (s.y - origin.y) / scale;
        }

        let mut rows: Vec<Weights> = vec![Vec::new(); self.nt_ as usize];
        let tris = if method == Method::NaturalNeighbour { delaunay(&site) } else { Vec::new() };
        for i in 0..self.nt_ as usize {
            let p = Pt { x: (cx[i] - origin.x) / scale, y: (cy[i] - origin.y) / scale };
            let mut w: Weights;
            if m == 1 {
                w = vec![(gid[0], 1.0)];
            } else if method == Method::NearestNeighbour {
                let mut nearest = 0usize;
                let mut best = dist2(&p, &site[0]);
                for g in 1..m {
                    let d = dist2(&p, &site[g]);
                    if d < best {
                        best = d;
                        nearest = g;
                    }
                }
                w = vec![(gid[nearest], 1.0)];
            } else {
                self.cell_method_[i] = CellMethod::NaturalNeighbour;
                w = Vec::new();
                if !tris.is_empty() {
                    w = laplace_weights(&p, &site, &gid, &tris);
                }
                if w.is_empty() {
                    w = idw_all(&p, &site, &gid);
                    self.diagnostics_.idwCells += 1;
                    self.cell_method_[i] = CellMethod::InverseDistance;
                }
            }
            let mut sum = 0.0;
            for e in &w {
                if !e.1.is_finite() || e.1 < 0.0 {
                    return Err(format!("2D rainfall: invalid weight at cell {}", i).into_bytes());
                }
                sum += e.1;
            }
            if !(sum > 0.0) || !sum.is_finite() {
                return Err(format!("2D rainfall: no valid weights at cell {}", i).into_bytes());
            }
            for e in w.iter_mut() {
                e.1 /= sum;
            }
            rows[i] = w;
        }

        let mut total = 0usize;
        for r in &rows {
            total += r.len();
        }
        self.w_gage_.reserve(total);
        self.w_val_.reserve(total);
        for i in 0..self.nt_ as usize {
            self.w_ptr_[i] = self.w_gage_.len() as i32;
            for e in &rows[i] {
                self.w_gage_.push(e.0);
                self.w_val_.push(e.1);
            }
        }
        self.w_ptr_[self.nt_ as usize] = self.w_gage_.len() as i32;
        self.ready_ = true;
        Ok(())
    }

    /// apply(rain, out): every cell's rainfall as its precomputed weighted sum of the gages'.
    pub fn apply(&self, rain: &[f64], out: &mut Vec<f64>) {
        if out.len() as i32 != self.nt_ {
            out.clear();
            out.resize(self.nt_ as usize, 0.0);
        }
        for i in 0..self.nt_ as usize {
            let mut acc = 0.0;
            let beg = self.w_ptr_[i] as usize;
            let end = self.w_ptr_[i + 1] as usize;
            for k in beg..end {
                acc += self.w_val_[k] * rain[self.w_gage_[k] as usize];
            }
            out[i] = acc;
        }
    }
}
