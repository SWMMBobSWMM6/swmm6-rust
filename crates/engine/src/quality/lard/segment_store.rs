// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/quality/lard/SegmentStore.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin,
// licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The LAGRANGIAN engine's per-link segment ring (a fixed-capacity slab per link, front = the
// upstream end), species-major concentrations. The checkpoint serializer is not translated
// (engine checkpoints are not on the 1D run path).

use crate::core::cmath::cmin;

pub const K_MAX_SEGMENTS_PER_LINK: i32 = 100;
pub const K_MERGE_ATOL: f64 = 1.0e-6;
pub const K_MERGE_RTOL: f64 = 1.0e-4;

/// lard::SegmentStore
#[derive(Debug, Clone)]
pub struct SegmentStore {
    pub nl_: i32,
    pub ns_: i32,
    pub cap_: i32,
    pub head_: Vec<i32>,
    pub count_: Vec<i32>,
    pub vol_: Vec<f64>,
    pub conc_: Vec<f64>,
}

impl Default for SegmentStore {
    fn default() -> Self {
        SegmentStore { nl_: 0, ns_: 0, cap_: K_MAX_SEGMENTS_PER_LINK, head_: Vec::new(), count_: Vec::new(), vol_: Vec::new(), conc_: Vec::new() }
    }
}

fn assign<T: Clone>(v: &mut Vec<T>, n: usize, x: T) {
    v.clear();
    v.resize(n, x);
}

impl SegmentStore {
    pub fn resize(&mut self, n_links: i32, n_species: i32, cap: i32) {
        self.nl_ = n_links;
        self.ns_ = n_species;
        self.cap_ = cap;
        let slab = self.nl_ as usize * self.cap_ as usize;
        assign(&mut self.head_, self.nl_ as usize, 0);
        assign(&mut self.count_, self.nl_ as usize, 0);
        assign(&mut self.vol_, slab, 0.0);
        assign(&mut self.conc_, slab * (if self.ns_ > 0 { self.ns_ } else { 1 }) as usize, 0.0);
    }
    pub fn count(&self, l: i32) -> i32 {
        self.count_[l as usize]
    }
    pub fn species(&self) -> i32 {
        self.ns_
    }
    pub fn capacity(&self) -> i32 {
        self.cap_
    }
    pub fn total_volume(&self, l: i32) -> f64 {
        let mut v = 0.0;
        for i in 0..self.count(l) {
            v += self.vol_[self.phys(l, i)];
        }
        v
    }
    pub fn seg_volume(&self, l: i32, i: i32) -> f64 {
        self.vol_[self.phys(l, i)]
    }
    pub fn seg_conc(&self, l: i32, i: i32, s: i32) -> f64 {
        self.conc_[self.cidx(l, i, s)]
    }
    pub fn set_seg_conc(&mut self, l: i32, i: i32, s: i32, c: f64) {
        let k = self.cidx(l, i, s);
        self.conc_[k] = c;
    }
    pub fn mean_conc(&self, l: i32, out: &mut [f64]) {
        for s in 0..self.ns_ as usize {
            out[s] = 0.0;
        }
        let mut v = 0.0;
        for i in 0..self.count(l) {
            let sv = self.vol_[self.phys(l, i)];
            v += sv;
            for s in 0..self.ns_ {
                out[s as usize] += sv * self.conc_[self.cidx(l, i, s)];
            }
        }
        if v > 0.0 {
            for s in 0..self.ns_ as usize {
                out[s] /= v;
            }
        }
    }
    pub fn clear_link(&mut self, l: i32) {
        self.head_[l as usize] = 0;
        self.count_[l as usize] = 0;
    }
    pub fn drain_back(&mut self, l: i32, vol_target: f64, mass_out: &mut [f64]) -> f64 {
        let mut drained = 0.0;
        while vol_target - drained > 0.0 && self.count(l) > 0 {
            let i = self.count(l) - 1;
            let p = self.phys(l, i);
            let take = cmin(self.vol_[p], vol_target - drained);
            for s in 0..self.ns_ {
                mass_out[s as usize] += take * self.conc_[self.cidx(l, i, s)];
            }
            self.vol_[p] -= take;
            drained += take;
            if self.vol_[p] <= 0.0 {
                self.count_[l as usize] -= 1;
            }
        }
        drained
    }
    pub fn drain_front(&mut self, l: i32, vol_target: f64, mass_out: &mut [f64]) -> f64 {
        let mut drained = 0.0;
        while vol_target - drained > 0.0 && self.count(l) > 0 {
            let p = self.phys(l, 0);
            let take = cmin(self.vol_[p], vol_target - drained);
            for s in 0..self.ns_ {
                mass_out[s as usize] += take * self.conc_[self.cidx(l, 0, s)];
            }
            self.vol_[p] -= take;
            drained += take;
            if self.vol_[p] <= 0.0 {
                self.head_[l as usize] = (self.head_[l as usize] + 1) % self.cap_;
                self.count_[l as usize] -= 1;
            }
        }
        drained
    }
    pub fn push_front(&mut self, l: i32, volume: f64, c: &[f64]) {
        if volume <= 0.0 {
            return;
        }
        if self.count(l) > 0 && self.mergeable(l, 0, c) {
            self.merge_into(l, 0, volume, c);
            return;
        }
        if self.count(l) >= self.cap_ {
            self.merge_front_pair(l);
        }
        self.head_[l as usize] = (self.head_[l as usize] - 1 + self.cap_) % self.cap_;
        self.count_[l as usize] += 1;
        let p = self.phys(l, 0);
        self.vol_[p] = volume;
        for s in 0..self.ns_ {
            let k = self.cidx(l, 0, s);
            self.conc_[k] = c[s as usize];
        }
    }
    pub fn reverse(&mut self, l: i32) {
        let n = self.count(l);
        let mut i = 0;
        let mut j = n - 1;
        while i < j {
            let (pi, pj) = (self.phys(l, i), self.phys(l, j));
            self.vol_.swap(pi, pj);
            for s in 0..self.ns_ {
                let (ci, cj) = (self.cidx(l, i, s), self.cidx(l, j, s));
                self.conc_.swap(ci, cj);
            }
            i += 1;
            j -= 1;
        }
    }
    pub fn add_species(&mut self, s: i32, delta: f64) {
        for l in 0..self.nl_ {
            for i in 0..self.count(l) {
                let k = self.cidx(l, i, s);
                self.conc_[k] += delta;
            }
        }
    }
    pub fn decay_species(&mut self, s: i32, factor: f64) -> f64 {
        let mut removed = 0.0;
        for l in 0..self.nl_ {
            for i in 0..self.count(l) {
                let ci = self.cidx(l, i, s);
                removed += self.vol_[self.phys(l, i)] * self.conc_[ci] * (1.0 - factor);
                self.conc_[ci] *= factor;
            }
        }
        removed
    }

    fn phys(&self, l: i32, i: i32) -> usize {
        let slot = (self.head_[l as usize] + i) % self.cap_;
        l as usize * self.cap_ as usize + slot as usize
    }
    fn cidx(&self, l: i32, i: i32, s: i32) -> usize {
        s as usize * self.nl_ as usize * self.cap_ as usize + self.phys(l, i)
    }
    fn mergeable(&self, l: i32, i: i32, c: &[f64]) -> bool {
        for s in 0..self.ns_ {
            let cf = self.conc_[self.cidx(l, i, s)];
            if (c[s as usize] - cf).abs() > K_MERGE_ATOL + K_MERGE_RTOL * cf.abs() {
                return false;
            }
        }
        true
    }
    fn merge_into(&mut self, l: i32, i: i32, volume: f64, c: &[f64]) {
        let p = self.phys(l, i);
        let v0 = self.vol_[p];
        let vt = v0 + volume;
        if vt <= 0.0 {
            return;
        }
        for s in 0..self.ns_ {
            let ci = self.cidx(l, i, s);
            self.conc_[ci] = (self.conc_[ci] * v0 + c[s as usize] * volume) / vt;
        }
        self.vol_[p] = vt;
    }
    fn merge_front_pair(&mut self, l: i32) {
        if self.count(l) < 2 {
            return;
        }
        let p0 = self.phys(l, 0);
        let p1 = self.phys(l, 1);
        let v0 = self.vol_[p0];
        let v1 = self.vol_[p1];
        let vt = v0 + v1;
        if vt > 0.0 {
            for s in 0..self.ns_ {
                let c0 = self.cidx(l, 0, s);
                let c1 = self.cidx(l, 1, s);
                self.conc_[c1] = (self.conc_[c0] * v0 + self.conc_[c1] * v1) / vt;
            }
        }
        self.vol_[p1] = vt;
        self.head_[l as usize] = (self.head_[l as usize] + 1) % self.cap_;
        self.count_[l as usize] -= 1;
    }
}
