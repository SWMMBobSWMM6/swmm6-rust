// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/quality/lard/RwptDispersion.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin,
// licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The LAGRANGIAN engine's random-walk particle tracking dispersion ([OPTIONS] DISPERSION RWPT):
// a counter-based random stream (splitmix64 hashes of seed, link, substep, particle and draw),
// vertical shear and an Ito random walk per particle, and the segment exchanges the particles'
// crossings carry. The checkpoint serializer is not translated.
//
// Differences in form, never in result:
//   * The unsigned 64-bit arithmetic of rwpt_hash wraps, as C++'s does (wrapping_add,
//     wrapping_mul).
//   * std::min({a, b, c}) is cmin(cmin(a, b), c).

use crate::core::cmath::{cmax, cmin};
use crate::quality::lard::segment_store::SegmentStore;

pub const K_RWPT_PARTICLES_PER_LINK: i32 = 2000;
pub const K_KAPPA: f64 = 0.41;
pub const K_NU_WATER: f64 = 1.05e-5;
pub const K_DM: f64 = 1.3454e-8;
pub const K_ETA_MIN: f64 = 1.0e-3;
pub const K_RE_TURB: f64 = 2000.0;

pub fn rwpt_hash(z: u64) -> u64 {
    let mut z = z.wrapping_add(0x9e3779b97f4a7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}

pub fn rwpt_uniform(seed: u64, link: u64, stepc: u64, particle: u64, draw: u64) -> f64 {
    let mut k = rwpt_hash(seed ^ rwpt_hash(link));
    k = rwpt_hash(k ^ rwpt_hash(stepc));
    k = rwpt_hash(k ^ rwpt_hash(particle.wrapping_mul(8).wrapping_add(draw)));
    ((k >> 11) as f64 + 0.5) * (1.0 / 9007199254740992.0)
}

pub fn rwpt_normal(seed: u64, link: u64, stepc: u64, particle: u64, draw: u64) -> f64 {
    let u1 = rwpt_uniform(seed, link, stepc, particle, draw);
    let u2 = rwpt_uniform(seed, link, stepc, particle, draw.wrapping_add(1));
    (-2.0 * u1.ln()).sqrt() * (6.283185307179586 * u2).cos() // LIBM (log, cos)
}

pub fn rwpt_u_dev(eta: f64, ubar: f64, ustar: f64, turbulent: bool) -> f64 {
    if turbulent {
        return (ustar / K_KAPPA) * (1.0 + cmax(eta, K_ETA_MIN).ln()); // LIBM
    }
    1.5 * ubar * eta * (2.0 - eta) - ubar
}

pub fn rwpt_d_eta(eta: f64, h: f64, ustar: f64, turbulent: bool) -> f64 {
    if !turbulent {
        return K_DM / (h * h);
    }
    let e = cmin(cmax(eta, 0.0), 1.0);
    K_KAPPA * ustar * e * (1.0 - e) / h
}

pub fn rwpt_d_eta_grad(eta: f64, h: f64, ustar: f64, turbulent: bool) -> f64 {
    if !turbulent {
        return 0.0;
    }
    K_KAPPA * ustar * (1.0 - 2.0 * eta) / h
}

pub fn rwpt_hyd_radius(_area: f64, depth: f64, diam: f64, circular: bool) -> f64 {
    if !circular || diam <= 0.0 || depth <= 0.0 {
        return cmax(depth, 0.0);
    }
    let ratio = cmin(cmax(1.0 - 2.0 * depth / diam, -1.0), 1.0);
    let theta = 2.0 * ratio.acos(); // LIBM
    if theta <= 1.0e-6 {
        return cmax(depth, 0.0);
    }
    let a = diam * diam / 8.0 * (theta - theta.sin()); // LIBM
    let p = diam * theta / 2.0;
    if p > 0.0 {
        a / p
    } else {
        cmax(depth, 0.0)
    }
}

/// lard::RwptDispersion
#[derive(Debug, Clone, Default)]
pub struct RwptDispersion {
    pub nl_: i32,
    pub npart_: i32,
    pub zeta_: Vec<f64>,
    pub eta_: Vec<f64>,
    pub seg_idx_: Vec<i32>,
    pub bounds_: Vec<f64>,
}

impl RwptDispersion {
    pub fn resize(&mut self, n_links: i32, n_per_link: i32) {
        self.nl_ = n_links;
        self.npart_ = n_per_link;
        let total = self.nl_ as usize * self.npart_ as usize;
        self.zeta_.clear();
        self.zeta_.resize(total, -1.0);
        self.eta_.clear();
        self.eta_.resize(total, 0.5);
        self.seg_idx_.clear();
        self.seg_idx_.resize(total, 0);
    }

    /// disperse: one substep of the particle field of link l.
    #[allow(clippy::too_many_arguments)]
    pub fn disperse(&mut self, segs: &mut SegmentStore, l: i32, ubar: f64, h: f64, rh: f64, n_manning: f64, v_in: f64, q_dt: f64, dt: f64, stepc: u64, seed: u64, delta: &mut [f64]) {
        let nseg = segs.count(l);
        if nseg < 2 || h <= 1.0e-3 || ubar <= 0.0 {
            return;
        }
        let vol = segs.total_volume(l);
        if vol <= 0.0 {
            return;
        }
        let ns = segs.species();
        let sf = (n_manning * ubar) * (n_manning * ubar) / (2.208 * cmax(rh, 1.0e-3).powf(4.0 / 3.0)); // LIBM
        let ustar = (32.174 * cmax(rh, 1.0e-3) * sf).sqrt();
        let turbulent = (ubar * 4.0 * rh / K_NU_WATER) >= K_RE_TURB;
        self.bounds_.clear();
        self.bounds_.resize(nseg as usize + 1, 0.0);
        for j in 0..nseg as usize {
            self.bounds_[j + 1] = self.bounds_[j] + segs.seg_volume(l, j as i32);
        }
        let base = l as usize * self.npart_ as usize;
        let ul = l as u64;
        for p in 0..self.npart_ {
            let ip = base + p as usize;
            let mut z = self.zeta_[ip];
            let mut e = self.eta_[ip];
            if z < 0.0 {
                z = vol * (p as f64 + 0.5) / self.npart_ as f64;
                e = rwpt_uniform(seed, ul, stepc, p as u64, 7);
            } else {
                z += v_in;
                if z >= vol {
                    z = cmin(v_in, vol * rwpt_uniform(seed, ul, stepc, p as u64, 8));
                    e = rwpt_uniform(seed, ul, stepc, p as u64, 9);
                }
            }
            let j_from = self.locate(z, nseg);
            let dz = (rwpt_u_dev(e, ubar, ustar, turbulent) / ubar) * q_dt;
            let d_eta = rwpt_d_eta(e, h, ustar, turbulent);
            let drift = rwpt_d_eta_grad(e, h, ustar, turbulent);
            e += drift * dt + rwpt_normal(seed, ul, stepc, p as u64, 3) * cmax(2.0 * d_eta * dt, 0.0).sqrt();
            while e < 0.0 || e > 1.0 {
                if e < 0.0 {
                    e = -e;
                }
                if e > 1.0 {
                    e = 2.0 - e;
                }
            }
            z += dz;
            while z < 0.0 || z > vol {
                if z < 0.0 {
                    z = -z;
                }
                if z > vol {
                    z = 2.0 * vol - z;
                }
            }
            let j_to = self.locate(z, nseg);
            if j_to != j_from {
                let vshare = vol / self.npart_ as f64;
                let mut j = j_from;
                let stepdir = if j_to > j_from { 1 } else { -1 };
                while j != j_to {
                    let jn = j + stepdir;
                    let bnd = self.bounds_[std::cmp::max(j, jn) as usize];
                    let pen = cmin((z - bnd).abs(), vshare);
                    for s in 0..ns {
                        delta[s as usize] = segs.seg_conc(l, j, s) * pen;
                    }
                    Self::exchange(segs, l, j, jn, delta, ns);
                    j = jn;
                }
            }
            self.zeta_[ip] = z;
            self.eta_[ip] = e;
            self.seg_idx_[ip] = j_to;
        }
    }

    fn locate(&self, z: f64, nseg: i32) -> i32 {
        let mut j = 0;
        while j + 1 < nseg && z >= self.bounds_[j as usize + 1] {
            j += 1;
        }
        j
    }

    fn exchange(segs: &mut SegmentStore, l: i32, j: i32, jn: i32, delta: &[f64], ns: i32) {
        for s in 0..ns {
            let vj = segs.seg_volume(l, j);
            let vn = segs.seg_volume(l, jn);
            if vj <= 0.0 || vn <= 0.0 {
                continue;
            }
            let cj = segs.seg_conc(l, j, s);
            let cn = segs.seg_conc(l, jn, s);
            let take = cmin(cmin(delta[s as usize], cj * vj), cmax(0.0, cj - cn) * vn);
            if take <= 0.0 {
                continue;
            }
            segs.set_seg_conc(l, j, s, (cj * vj - take) / vj);
            segs.set_seg_conc(l, jn, s, (cn * vn + take) / vn);
        }
    }
}
