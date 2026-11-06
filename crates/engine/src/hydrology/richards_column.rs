// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/RichardsColumn.hpp and RichardsColumn.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 6 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The Richards-equation soil column (namespace richards): the van Genuchten-Mualem
// conductivity, the column types Cell, Column and Bottom, and advance(), the adaptive BDF1
// integrator with step doubling (Batch: rhs, implicit, run). theta, storage and pressure are
// translated in hydrology/lid_node_init.rs (richards_theta, richards_storage,
// richards_pressure), where initialize first needs them; the two valid() overloads are in
// hydrology/lid_node.rs. Material, Options, Transfer and Report are in data/lid_node.rs.
//
// Differences in form, never in result:
//   * perf::ProcessTimer (timing only) is left out.
//   * Batch<Packed> is the unpacked (Batch<false>) form only. The packed form copies the
//     same cell coefficients into private arrays (cellAt returns the same values), and
//     upstream builds it only on the OpenMP path, which needs at least 256 columns and more
//     than one thread; the reference is held to one thread. On success the OpenMP path
//     publishes what the serial walk publishes, group for group, and on failure it resets
//     the reports of later jobs to Report{}, which the serial walk never runs either.
//   * Batch borrows its columns and reports by pointer upstream; here it holds their indices
//     (`ids`) into the trial columns and the reports, which every method takes as arguments.
//     The Batch scratch members head, K and the lane mask `done` are passed to rhs
//     explicitly; f and fp are members as upstream (run's local `f` shadows the member, as
//     upstream's does).
//   * advance's duplicate-pointer test (std::set<const Column*> seen) has no counterpart: the
//     columns are a slice, so two entries are never the same column.
//   * std::nextafter(rate, 0.0) is `next_toward_zero`, exact for every finite rate.

use crate::core::cmath::{cclamp, cmax, cmin};
use crate::data::lid_node::richards::{Material, Options, Report, Transfer};
use crate::hydrology::lid_node::{richards_valid_material, richards_valid_options};
pub use crate::hydrology::lid_node_init::{richards_pressure as pressure, richards_storage as storage, richards_theta as theta};

/// richards::Cell
#[derive(Debug, Clone, Copy)]
pub struct Cell {
    pub volume: f64,
    pub area: f64,
    pub z: f64,
    pub dz: f64,
    /// Surface theta_s is the open storage fraction; porous theta_s is saturated pore water.
    pub theta_s: f64,
    pub ks: f64,
    pub wilting: f64,
    pub potential_et: f64,
    pub material: Material,
}

impl Default for Cell {
    fn default() -> Self {
        Cell {
            volume: 0.0,
            area: 0.0,
            z: 0.0,
            dz: 0.0,
            theta_s: 1.0,
            ks: 0.0,
            wilting: 0.0,
            potential_et: 0.0,
            material: Material::default(),
        }
    }
}

/// richards::Bottom
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Bottom {
    Sealed,
    FreeDrainage,
    Head,
}

/// richards::Column: a surface reservoir first, then porous cells top to bottom.
#[derive(Debug, Clone)]
pub struct Column {
    pub cells: Vec<Cell>,
    pub water: Vec<f64>,
    pub options: Options,
    pub bottom: Bottom,
    pub bottom_k: f64,
    pub bottom_head: f64,
    pub bottom_distance: f64,
    pub bottom_out_limit: f64,
    pub bottom_in_limit: f64,
    pub bottom_area: f64,
    pub bottom_factor: f64,
    /// accepted, held surface supply, m3/s
    pub inflow: f64,
}

impl Default for Column {
    fn default() -> Self {
        Column {
            cells: Vec::new(),
            water: Vec::new(),
            options: Options::default(),
            bottom: Bottom::Sealed,
            bottom_k: 0.0,
            bottom_head: 0.0,
            bottom_distance: 0.0,
            bottom_out_limit: f64::INFINITY,
            bottom_in_limit: f64::INFINITY,
            bottom_area: 0.0,
            bottom_factor: 1.0,
            inflow: 0.0,
        }
    }
}

/// richards::conductivity: van Genuchten-Mualem unsaturated conductivity (m/s).
pub fn conductivity(p: &Material, ks: f64, h: f64) -> f64 {
    if h >= 0.0 {
        return ks;
    }
    let m = 1.0 - 1.0 / p.n;
    let se = (1.0 + (-p.alpha * h).powf(p.n)).powf(-m); // LIBM pow
    let term = -(m * (-se.powf(1.0 / m)).ln_1p()).exp_m1(); // LIBM pow, log1p, expm1
    ks * se.powf(p.l) * term * term // LIBM pow
}

/// std::nextafter(x, 0.0)
fn next_toward_zero(x: f64) -> f64 {
    if x == 0.0 || x.is_nan() {
        return if x.is_nan() { x } else { 0.0 };
    }
    f64::from_bits(x.to_bits() - 1)
}

struct BoundaryRates {
    out: Vec<f64>,
    in_: Vec<f64>,
}

/// reservedRate: a remaining directional volume reservation as a safe step rate.
fn reserved_rate(limit: f64, spent: f64, step: f64) -> f64 {
    if !limit.is_finite() {
        return f64::INFINITY;
    }
    if step <= 0.0 || spent >= limit {
        return 0.0;
    }
    let mut rate = cmin((limit - spent) / step, f64::MAX);
    // Round toward the budget, using the same multiply/add as the receipt.
    while spent + step * rate > limit {
        rate = next_toward_zero(rate);
    }
    rate
}

/// Batch<false>: one group of columns with the same cell count, cell-major, lane-minor.
struct Batch {
    n: usize,
    lanes: usize,
    /// the trial column (and report) of each lane
    ids: Vec<usize>,
    y: Vec<f64>,
    volume: Vec<f64>,
    lower_bound: Vec<f64>,
    head: Vec<f64>,
    k: Vec<f64>,
    f: Vec<f64>,
    fp: Vec<f64>,
    residual: Vec<f64>,
    lo: Vec<f64>,
    di: Vec<f64>,
    up: Vec<f64>,
    delta: Vec<f64>,
    perturbed: Vec<f64>,
    epsilon: Vec<f64>,
    done: Vec<u8>,
}

/// Batch::rhs: conservative intercell and boundary water rates (m3/s), flux positive
/// downward. The Batch members head and K are `head` and `kv`.
#[allow(clippy::too_many_arguments)]
fn rhs(
    n: usize,
    lanes: usize,
    ids: &[usize],
    cols: &[Column],
    reps: &mut [Report],
    head: &mut Vec<f64>,
    kv: &mut Vec<f64>,
    w: &[f64],
    f: &mut Vec<f64>,
    mut faces: Option<&mut Vec<f64>>,
    mut et: Option<&mut Vec<f64>>,
    inactive: Option<&[u8]>,
    rates: Option<&BoundaryRates>,
    mut limited: Option<&mut Vec<i32>>,
) {
    let off = |a: usize| -> bool { inactive.is_some_and(|m| m[a] != 0) };
    f.clear();
    f.resize(w.len(), 0.0);
    head.clear();
    head.resize(w.len(), 0.0);
    kv.clear();
    kv.resize(w.len(), 0.0);
    if let Some(fc) = faces.as_deref_mut() {
        fc.clear();
        fc.resize(n * lanes, 0.0);
    }
    if let Some(e) = et.as_deref_mut() {
        e.clear();
        e.resize(w.len(), 0.0);
    }
    if let Some(l) = limited.as_deref_mut() {
        l.clear();
        l.resize(lanes, 0);
    }
    for a in 0..lanes {
        if !off(a) {
            reps[ids[a]].rhs += 1;
        }
    }
    for a in 0..lanes {
        if !off(a) && cols[ids[a]].inflow > 0.0 {
            f[a] += cols[ids[a]].inflow;
        }
    }
    for k in 0..n {
        for a in 0..lanes {
            if off(a) {
                continue;
            }
            let i = k * lanes + a;
            let c = &cols[ids[a]].cells[k];
            head[i] = c.z
                + if k != 0 {
                    pressure(&c.material, c.theta_s, w[i] / c.volume)
                } else {
                    cmax(0.0, w[i]) / (c.area * c.theta_s)
                };
            kv[i] = if k != 0 { conductivity(&c.material, c.ks, head[i] - c.z) } else { 0.0 };
            let stress = if k != 0 {
                cclamp((w[i] / c.volume - c.wilting) / cmax(c.theta_s - c.wilting, 1.0e-12), 0.0, 1.0)
            } else {
                cclamp(w[i] / (c.area * c.theta_s * 1.0e-6), 0.0, 1.0)
            };
            let e = c.potential_et * stress;
            f[i] -= e;
            if let Some(et) = et.as_deref_mut() {
                et[i] = e;
            }
        }
    }
    for k in 0..n.saturating_sub(1) {
        for a in 0..lanes {
            if off(a) {
                continue;
            }
            let i = k * lanes + a;
            let j = i + lanes;
            let c = &cols[ids[a]].cells[k];
            let d = &cols[ids[a]].cells[k + 1];
            let area = cmin(c.area, d.area);
            let mut conductance = 0.0;
            // Arithmetic conductivity within a material follows openRE (Ireson et al. 2023,
            // Eq. 11/14).
            if k == 0 {
                conductance = area * (d.ks + kv[j]) / d.dz;
            } else if c.material.alpha == d.material.alpha
                && c.material.n == d.material.n
                && c.material.theta_r == d.material.theta_r
                && c.material.l == d.material.l
                && c.ks == d.ks
            {
                conductance = area * (c.dz * kv[j] + d.dz * kv[i]) / (c.dz + d.dz) / (0.5 * (c.dz + d.dz));
            } else if kv[i] > 0.0 && kv[j] > 0.0 {
                conductance = area / (0.5 * c.dz / kv[i] + 0.5 * d.dz / kv[j]);
            }
            let mut q = conductance * (head[i] - head[j]);
            if k == 0 && q > 0.0 {
                q *= cclamp(w[i] / (c.area * c.theta_s * 1.0e-6), 0.0, 1.0);
            }
            f[i] -= q;
            f[j] += q;
            if let Some(fc) = faces.as_deref_mut() {
                fc[i] = q;
            }
        }
    }
    for a in 0..lanes {
        if off(a) {
            continue;
        }
        let i = (n - 1) * lanes + a;
        let col = &cols[ids[a]];
        let c = &col.cells[n - 1];
        let area = if col.bottom_area > 0.0 { col.bottom_area } else { c.area };
        let mut q = 0.0;
        if col.bottom == Bottom::FreeDrainage {
            q = cmin(kv[i], col.bottom_k) * area;
        } else if col.bottom == Bottom::Head && kv[i] > 0.0 && col.bottom_k > 0.0 {
            let resistance = 0.5 * c.dz / kv[i] + col.bottom_distance / col.bottom_k;
            q = area * (head[i] - col.bottom_head) / resistance;
        }
        q *= col.bottom_factor; // once, on the whole series conductance
        if let Some(r) = rates {
            if q > r.out[a] {
                q = r.out[a];
                if let Some(l) = limited.as_deref_mut() {
                    l[a] = 1;
                }
            } else if q < -r.in_[a] {
                q = -r.in_[a];
                if let Some(l) = limited.as_deref_mut() {
                    l[a] = -1;
                }
            }
        }
        f[i] -= q;
        if let Some(fc) = faces.as_deref_mut() {
            fc[i] = q;
        }
    }
}

impl Batch {
    /// Batch(c, r): copy the lanes' water and volumes; residual-volume lower bounds.
    fn new(cols: &[Column], ids: Vec<usize>) -> Batch {
        let n = cols[ids[0]].cells.len();
        let lanes = ids.len();
        let size = n * lanes;
        let mut b = Batch {
            n,
            lanes,
            ids,
            y: vec![0.0; size],
            volume: vec![0.0; size],
            lower_bound: vec![0.0; size],
            head: vec![0.0; size],
            k: vec![0.0; size],
            f: vec![0.0; size],
            fp: vec![0.0; size],
            residual: vec![0.0; size],
            lo: vec![0.0; size],
            di: vec![0.0; size],
            up: vec![0.0; size],
            delta: vec![0.0; size],
            perturbed: vec![0.0; size],
            epsilon: vec![0.0; size],
            done: vec![0; lanes],
        };
        for k in 0..n {
            for a in 0..lanes {
                let col = &cols[b.ids[a]];
                let cell = &col.cells[k];
                let i = k * lanes + a;
                b.y[i] = col.water[k];
                b.volume[i] = cell.volume;
                b.lower_bound[i] = if k != 0 { cell.material.theta_r * cell.volume } else { 0.0 };
            }
        }
        b
    }

    /// Batch::implicit: one lane-wise backward-Euler stage (up to 35 damped Newton
    /// iterations with three-colour finite-difference tridiagonal Jacobians).
    fn implicit(&mut self, cols: &[Column], reps: &mut [Report], start: &[f64], steps: &[f64], out: &mut Vec<f64>, success: &mut Vec<u8>, rates: Option<&BoundaryRates>) {
        let (n, lanes) = (self.n, self.lanes);
        out.clear();
        out.extend_from_slice(start);
        self.done.clear();
        self.done.resize(lanes, 0);
        success.clear();
        success.resize(lanes, 0);
        self.residual.iter_mut().for_each(|x| *x = 0.0);
        for a in 0..lanes {
            if steps[a] == 0.0 {
                self.done[a] = 1;
                success[a] = 1;
            }
        }
        for _it in 0..35 {
            {
                let mut f = std::mem::take(&mut self.f);
                rhs(n, lanes, &self.ids, cols, reps, &mut self.head, &mut self.k, out, &mut f, None, None, Some(&self.done), rates, None);
                self.f = f;
            }
            let mut all = true;
            for a in 0..lanes {
                if self.done[a] == 0 {
                    reps[self.ids[a]].newton += 1;
                    let mut norm = 0.0;
                    for k in 0..n {
                        let i = k * lanes + a;
                        self.residual[i] = out[i] - start[i] - steps[a] * self.f[i];
                        if !self.residual[i].is_finite() {
                            norm = f64::INFINITY;
                            break;
                        }
                        norm = cmax(norm, self.residual[i].abs() / self.volume[i]);
                    }
                    if norm.is_finite() && norm <= cmin(1.0e-11, cols[self.ids[a]].options.atol * 1.0e-3) {
                        self.done[a] = 1;
                        success[a] = 1;
                    } else {
                        all = false;
                    }
                }
            }
            if all {
                break;
            }
            self.lo.iter_mut().for_each(|x| *x = 0.0);
            self.up.iter_mut().for_each(|x| *x = 0.0);
            self.di.iter_mut().for_each(|x| *x = 1.0);
            // Three-color differentiation of a nearest-neighbor stencil.
            for color in 0..3usize {
                self.perturbed.clear();
                self.perturbed.extend_from_slice(out);
                let mut k = color;
                while k < n {
                    for a in 0..lanes {
                        if self.done[a] == 0 {
                            let i = k * lanes + a;
                            // Resolve the narrow saturated storage branch rather than
                            // perturbing across it with a fixed moisture increment.
                            let wet_distance = (out[i] - cols[self.ids[a]].cells[k].theta_s * self.volume[i]).abs();
                            self.epsilon[i] = 1.0e-8 * if k != 0 { cmax(wet_distance, 1.0e-4 * self.volume[i]) } else { self.volume[i] };
                            self.perturbed[i] += self.epsilon[i];
                        }
                    }
                    k += 3;
                }
                {
                    let mut fp = std::mem::take(&mut self.fp);
                    let perturbed = std::mem::take(&mut self.perturbed);
                    rhs(n, lanes, &self.ids, cols, reps, &mut self.head, &mut self.k, &perturbed, &mut fp, None, None, Some(&self.done), rates, None);
                    self.fp = fp;
                    self.perturbed = perturbed;
                }
                let mut k = color;
                while k < n {
                    for a in 0..lanes {
                        if self.done[a] == 0 {
                            let j = k * lanes + a;
                            let lo_row = if k >= 1 { k - 1 } else { 0 };
                            let hi_row = cmin(n - 1, k + 1);
                            for row in lo_row..=hi_row {
                                let i = row * lanes + a;
                                let value = -steps[a] * (self.fp[i] - self.f[i]) / self.epsilon[j];
                                if row == k {
                                    self.di[i] += value;
                                } else if row < k {
                                    self.up[i] = value;
                                } else {
                                    self.lo[i] = value;
                                }
                            }
                        }
                    }
                    k += 3;
                }
            }
            for k in 0..n {
                for a in 0..lanes {
                    let i = k * lanes + a;
                    self.delta[i] = if self.done[a] != 0 { 0.0 } else { -self.residual[i] };
                    if k != 0 {
                        let factor = self.lo[i] / self.di[i - lanes];
                        self.di[i] -= factor * self.up[i - lanes];
                        self.delta[i] -= factor * self.delta[i - lanes];
                    }
                }
            }
            for k in (0..n).rev() {
                for a in 0..lanes {
                    let i = k * lanes + a;
                    if k + 1 < n {
                        self.delta[i] -= self.up[i] * self.delta[i + lanes];
                    }
                    self.delta[i] /= self.di[i];
                }
            }
            for a in 0..lanes {
                if self.done[a] == 0 {
                    let mut damping = 1.0;
                    for k in 0..n {
                        let i = k * lanes + a;
                        if !self.delta[i].is_finite() {
                            self.done[a] = 1;
                            damping = 0.0;
                            break;
                        }
                        if self.delta[i] < 0.0 {
                            damping = cmin(damping, 0.9 * (out[i] - self.lower_bound[i]) / -self.delta[i]);
                        }
                    }
                    if damping <= 1.0e-10 {
                        self.done[a] = 1;
                        continue;
                    }
                    for k in 0..n {
                        out[k * lanes + a] += damping * self.delta[k * lanes + a];
                    }
                }
            }
        }
    }

    /// Batch::rates (run's lambda): the lanes' directional bottom caps for one stage.
    fn rates(&self, cols: &[Column], reps: &[Report], b: &mut BoundaryRates, steps: &[f64], first_faces: Option<&[f64]>) {
        let (n, lanes) = (self.n, self.lanes);
        b.out.resize(lanes, 0.0);
        b.in_.resize(lanes, 0.0);
        for a in 0..lanes {
            let first = match first_faces {
                Some(ff) => steps[a] * ff[(n - 1) * lanes + a],
                None => 0.0,
            };
            let col = &cols[self.ids[a]];
            let rep = &reps[self.ids[a]];
            b.out[a] = reserved_rate(col.bottom_out_limit, rep.bottom_outflow + cmax(0., first), steps[a]);
            b.in_[a] = reserved_rate(col.bottom_in_limit, rep.bottom_inflow + cmax(0., -first), steps[a]);
        }
    }

    /// Batch::run: integrate this group with independent adaptive lane clocks.
    fn run(&mut self, cols: &mut [Column], reps: &mut [Report], dt: f64) {
        let (n, lanes) = (self.n, self.lanes);
        let mut time = vec![0.0f64; lanes];
        let mut step = vec![0.0f64; lanes];
        let mut full: Vec<f64> = Vec::new();
        let mut half: Vec<f64> = Vec::new();
        let mut two: Vec<f64> = Vec::new();
        let mut f: Vec<f64> = Vec::new();
        let mut faces1: Vec<f64> = Vec::new();
        let mut faces2: Vec<f64> = Vec::new();
        let mut et1: Vec<f64> = Vec::new();
        let mut et2: Vec<f64> = Vec::new();
        let mut h: Vec<f64>;
        let mut good = vec![0u8; lanes];
        let mut good1 = vec![0u8; lanes];
        let mut good2 = vec![0u8; lanes];
        let mut inactive = vec![0u8; lanes];
        let mut limited1: Vec<i32> = Vec::new();
        let mut limited2: Vec<i32> = Vec::new();
        let mut full_rates = BoundaryRates { out: Vec::new(), in_: Vec::new() };
        let mut half_rates = BoundaryRates { out: Vec::new(), in_: Vec::new() };
        let mut second_rates = BoundaryRates { out: Vec::new(), in_: Vec::new() };
        for a in 0..lanes {
            step[a] = cmin(dt, cols[self.ids[a]].options.max_step);
        }
        let bounded = self.ids.iter().any(|&id| cols[id].bottom_out_limit.is_finite() || cols[id].bottom_in_limit.is_finite());
        let mut attempts = 0i32;
        loop {
            // *std::min_element(time.begin(), time.end()) < dt
            let mut smallest = time[0];
            for &t in time.iter().skip(1) {
                if t < smallest {
                    smallest = t;
                }
            }
            if !(smallest < dt) {
                break;
            }
            attempts += 1;
            if attempts > 200000 {
                for &id in self.ids.iter() {
                    reps[id].ok = false;
                    reps[id].error = b"Richards adaptive-step limit exceeded".to_vec();
                }
                return;
            }
            for a in 0..lanes {
                step[a] = cmin(step[a], cmax(0.0, dt - time[a]));
            }
            let c: &[Column] = cols;
            if bounded {
                self.rates(c, reps, &mut full_rates, &step, None);
            }
            self.implicit(c, reps, &self.y.clone(), &step, &mut full, &mut good, if bounded { Some(&full_rates) } else { None });
            h = step.clone();
            for v in h.iter_mut() {
                *v *= 0.5;
            }
            if bounded {
                self.rates(c, reps, &mut half_rates, &h, None);
            }
            self.implicit(c, reps, &self.y.clone(), &h, &mut half, &mut good1, if bounded { Some(&half_rates) } else { None });
            for a in 0..lanes {
                inactive[a] = (step[a] == 0.0) as u8;
            }
            if bounded {
                rhs(n, lanes, &self.ids, c, reps, &mut self.head, &mut self.k, &half, &mut f, Some(&mut faces1), Some(&mut et1), Some(&inactive), Some(&half_rates), Some(&mut limited1));
                self.rates(c, reps, &mut second_rates, &h, Some(&faces1));
            }
            self.implicit(c, reps, &half, &h, &mut two, &mut good2, if bounded { Some(&second_rates) } else { None });
            if !bounded {
                rhs(n, lanes, &self.ids, c, reps, &mut self.head, &mut self.k, &half, &mut f, Some(&mut faces1), Some(&mut et1), Some(&inactive), None, None);
            }
            rhs(
                n,
                lanes,
                &self.ids,
                c,
                reps,
                &mut self.head,
                &mut self.k,
                &two,
                &mut f,
                Some(&mut faces2),
                Some(&mut et2),
                Some(&inactive),
                if bounded { Some(&second_rates) } else { None },
                if bounded { Some(&mut limited2) } else { None },
            );
            for a in 0..lanes {
                if step[a] == 0.0 {
                    continue;
                }
                let id = self.ids[a];
                let opt = cols[id].options;
                let inflow = cols[id].inflow;
                let report = &mut reps[id];
                let mut error = 0.0;
                for k in 0..n {
                    let i = k * lanes + a;
                    let scale = opt.atol * self.volume[i] + opt.rtol * cmax(two[i].abs(), self.y[i].abs());
                    error = cmax(error, (two[i] - full[i]).abs() / scale);
                }
                if good[a] != 0 && good1[a] != 0 && good2[a] != 0 && error.is_finite() && error <= 1.0 {
                    report.accepted += 1;
                    report.min_step = if report.min_step > 0.0 { cmin(report.min_step, step[a]) } else { step[a] };
                    // Commit two accepted half steps and their endpoint flux integrals in
                    // chronological order, retaining directions.
                    for stage in 0..2 {
                        for k in 0..n {
                            if k == 0 && inflow > 0.0 {
                                let received = inflow * h[a];
                                report.inflow += received;
                                report.transfers.push(Transfer { from: -2, to: 0, volume: received });
                            }
                            let i = k * lanes + a;
                            let v = h[a] * if stage != 0 { faces2[i] } else { faces1[i] };
                            if k + 1 < n {
                                if v > 0.0 {
                                    report.transfers.push(Transfer { from: k as i32, to: (k + 1) as i32, volume: v });
                                } else if v < 0.0 {
                                    report.transfers.push(Transfer { from: (k + 1) as i32, to: k as i32, volume: -v });
                                }
                            } else {
                                report.bottom_volume += v;
                                if v > 0.0 {
                                    report.bottom_outflow += v;
                                } else if v < 0.0 {
                                    report.bottom_inflow -= v;
                                }
                                if v > 0.0 {
                                    report.transfers.push(Transfer { from: k as i32, to: -4, volume: v });
                                } else if v < 0.0 {
                                    report.transfers.push(Transfer { from: -5, to: k as i32, volume: -v });
                                }
                            }
                            let e = h[a] * if stage != 0 { et2[i] } else { et1[i] };
                            report.evaporation += e;
                            if e > 0.0 {
                                report.transfers.push(Transfer { from: k as i32, to: -3, volume: e });
                            }
                        }
                    }
                    if bounded {
                        for limited in [limited1[a], limited2[a]] {
                            report.bottom_out_limited += (limited == 1) as i32;
                            report.bottom_in_limited += (limited == -1) as i32;
                        }
                    }
                    for k in 0..n {
                        self.y[k * lanes + a] = two[k * lanes + a];
                    }
                    time[a] = cmin(dt, time[a] + step[a]);
                    step[a] = cmin(opt.max_step, step[a] * cclamp(0.8 / cmax(error, 1.0e-8).sqrt(), 0.5, 2.0));
                } else {
                    report.rejected += 1;
                    step[a] *= if error.is_finite() && error > 1.0 { cclamp(0.8 / error.sqrt(), 0.1, 0.5) } else { 0.5 };
                    if step[a] < 1.0e-9 * cmax(1.0, dt) {
                        report.ok = false;
                        report.error = b"Richards step underflow / Newton failure".to_vec();
                        return;
                    }
                }
            }
        }
        for a in 0..lanes {
            let id = self.ids[a];
            let col = &mut cols[id];
            let report = &mut reps[id];
            let mut initial = 0.0;
            for &w in col.water.iter() {
                initial += w;
            }
            let mut final_ = 0.0;
            for k in 0..n {
                col.water[k] = self.y[k * lanes + a];
                final_ += col.water[k];
            }
            report.balance = final_ - initial + report.bottom_volume + report.evaporation;
            if col.inflow > 0.0 {
                report.balance -= report.inflow;
            }
        }
    }
}

fn fail_all(reports: &mut [Report], what: &[u8]) {
    for r in reports.iter_mut() {
        r.ok = false;
        r.error = what.to_vec();
    }
}

/// richards::advance: advance independently adaptive Richards columns with all-column
/// rollback. On success every column's water is replaced; on any failure none is.
pub fn advance(columns: &mut [Column], dt: f64) -> Vec<Report> {
    let mut reports: Vec<Report> = vec![Report::default(); columns.len()];
    if !dt.is_finite() || dt <= 0.0 {
        fail_all(&mut reports, b"Invalid Richards interval");
        return reports;
    }
    // Work on copies: any failure is transactional across the complete batch.
    let mut trial: Vec<Column> = Vec::with_capacity(columns.len());
    for c in columns.iter() {
        if c.cells.len() < 2
            || c.water.len() != c.cells.len()
            || !richards_valid_options(&c.options)
            || !c.bottom_k.is_finite()
            || c.bottom_k < 0.0
            || !c.bottom_head.is_finite()
            || !c.bottom_distance.is_finite()
            || c.bottom_distance < 0.0
            || c.bottom_out_limit.is_nan()
            || c.bottom_out_limit < 0.0
            || c.bottom_in_limit.is_nan()
            || c.bottom_in_limit < 0.0
            || !c.bottom_area.is_finite()
            || c.bottom_area < 0.0
            || !c.bottom_factor.is_finite()
            || c.bottom_factor < 0.0
            || c.bottom_factor > 1.0
            || !c.inflow.is_finite()
            || c.inflow < 0.0
        {
            fail_all(&mut reports, b"Invalid Richards column");
            return reports;
        }
        for k in 0..c.cells.len() {
            let cell = &c.cells[k];
            let bad_state = if k != 0 { !richards_valid_material(&cell.material, cell.theta_s) || c.water[k] <= cell.material.theta_r * cell.volume } else { c.water[k] < 0.0 };
            if !cell.volume.is_finite()
                || cell.volume <= 0.0
                || !cell.area.is_finite()
                || cell.area <= 0.0
                || !cell.z.is_finite()
                || !cell.dz.is_finite()
                || cell.dz <= 0.0
                || !cell.ks.is_finite()
                || cell.ks < 0.0
                || !cell.theta_s.is_finite()
                || cell.theta_s <= 0.0
                || cell.theta_s > 1.0
                || !cell.wilting.is_finite()
                || cell.wilting < 0.0
                || cell.wilting >= cell.theta_s
                || !cell.potential_et.is_finite()
                || cell.potential_et < 0.0
                || !c.water[k].is_finite()
                || bad_state
            {
                fail_all(&mut reports, b"Invalid Richards state/material");
                return reports;
            }
        }
        trial.push(c.clone());
    }
    // std::map<std::size_t, std::vector<std::size_t>> groups: ascending cell count
    let mut groups: std::collections::BTreeMap<usize, Vec<usize>> = std::collections::BTreeMap::new();
    for a in 0..trial.len() {
        groups.entry(trial[a].cells.len()).or_default().push(a);
    }
    for ids in groups.values() {
        let mut first = 0usize;
        while first < ids.len() {
            let last = cmin(ids.len(), first + 64);
            let lane_ids: Vec<usize> = ids[first..last].to_vec();
            let mut batch = Batch::new(&trial, lane_ids);
            batch.run(&mut trial, &mut reports, dt);
            for j in first..last {
                if !reports[ids[j]].ok {
                    for report in reports.iter_mut() {
                        report.ok = false;
                        if report.error.is_empty() {
                            report.error = b"Richards batch rolled back".to_vec();
                        }
                    }
                    return reports;
                }
            }
            first += 64;
        }
    }
    for a in 0..trial.len() {
        columns[a].water = std::mem::take(&mut trial[a].water);
    }
    reports
}
