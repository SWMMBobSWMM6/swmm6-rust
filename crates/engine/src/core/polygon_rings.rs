// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/PolygonRings.hpp (HydroCouple/openswmm.engine,
// branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache
// License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Rings of a flat spatial polygon and the even-odd point-in-polygon test (upstream namespace
// openswmm::polygon). Used by SurfaceRouter2D::initialize (the SUBCATCH_AQUIFER destination).

#![allow(dead_code)]

/// A ring: [first, first + count), excluding a repeated closing vertex.
#[derive(Debug, Clone, Copy)]
pub struct Ring {
    pub first: usize,
    pub count: usize,
}

/// rings(x, y): successive explicitly closed rings, then a final open ring of 3 or more.
pub fn rings(x: &[f64], y: &[f64]) -> Vec<Ring> {
    let mut out = Vec::new();
    if x.len() != y.len() {
        return out;
    }
    let mut first = 0usize;
    for i in 1..x.len() {
        if i - first >= 3 && x[i] == x[first] && y[i] == y[first] {
            out.push(Ring { first, count: i - first });
            first = i + 1;
        }
    }
    if x.len() >= first + 3 {
        out.push(Ring { first, count: x.len() - first });
    }
    out
}

/// contains(x, y, ring, px, py): the ray-crossing test on one ring.
pub fn contains_ring(x: &[f64], y: &[f64], ring: Ring, px: f64, py: f64) -> bool {
    let mut inside = false;
    let mut j = ring.count.wrapping_sub(1);
    for i in 0..ring.count {
        let a = ring.first + i;
        let b = ring.first + j;
        if (y[a] > py) != (y[b] > py) && px < (x[b] - x[a]) * (py - y[a]) / (y[b] - y[a]) + x[a] {
            inside = !inside;
        }
        j = i;
    }
    inside
}

/// contains(x, y, px, py): even-odd over every ring.
pub fn contains(x: &[f64], y: &[f64], px: f64, py: f64) -> bool {
    let mut inside = false;
    for ring in rings(x, y) {
        inside ^= contains_ring(x, y, ring, px, py);
    }
    inside
}
