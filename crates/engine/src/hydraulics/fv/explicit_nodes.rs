// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/fv/ExplicitFvSolver.cpp and
// ExplicitFvSolver.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The explicit FV step's cells and nodes: relaxNodeFluxes, solveNodes, relaxOneNode,
// solveAlgebraicNode (the bracketed junction solve with its pass-through and degree-1
// paths), updateCells, refreshStructFlows, refreshDummyFlows, capWithdrawal, the accepted-
// water bookings (bookFaceWater, bookSpliceReceipt, bookStructureWater), updateNodes,
// settleAlgebraicNode, injectDivertedLateral, bookLateralSpills, frictionFor,
// applyNodeCapacity, and the header's inline helpers (nodeLateral, bookDivertedLateral,
// nodeQDummy, incidentPressurized).
//
// Differences in form, never in result:
//   * solveNodes' dt_of callable is (dt, tiered): tiered gives (1 << node_tier_[n]) * dt.
//   * solveAlgebraicNode's lambdas (residual, finalize, resistAt, prescribeDegree1) are
//     methods; the state they capture (h_best, r_best, h_last) is a struct passed along.
//   * refreshDummyFlows' arrived(node) callable is the Arrived enum (explicit_step.rs).

use crate::core::cmath::{cmax, cmin};
use crate::hydraulics::fv::explicit_step::{Arrived, FvStep};
use crate::hydraulics::fv::kernels as k;
use crate::hydraulics::fv::{FvGeometry, K_NODE_VIRTUAL};
use crate::hydraulics::hyd_closure_kernels::{fric_slope_dw, fric_slope_hw};
use crate::hydraulics::xsect_batch::XSectShape;

const K_ALG_Q_EPS: f64 = 1.0e-12;
const K_ALG_MAX_ITERS: i32 = 16;
const K_ALG_HEAD_TOL: f64 = 1.0e-6;
const K_ALG_MAX_EXPAND: i32 = 64;

/// solveAlgebraicNode's captured search state.
struct AlgSearch {
    h_best: f64,
    r_best: f64,
    h_last: f64,
}

impl<'a> FvStep<'a> {
    /// nodeLateral(forcing, un)
    #[inline]
    pub(crate) fn node_lateral(&self, un: usize) -> f64 {
        if !self.s.node_lat_div_.is_empty() && self.s.node_lat_div_[un] != 0 {
            return 0.0;
        }
        self.fc.node_lateral[un]
    }

    /// nodeQDummy(un)
    #[inline]
    pub(crate) fn node_q_dummy(&self, un: usize) -> f64 {
        if self.s.node_qdummy_.is_empty() { 0.0 } else { self.s.node_qdummy_[un] }
    }

    /// bookDivertedLateral(un, dt, forcing)
    pub(crate) fn book_diverted_lateral(&mut self, un: usize, dt: f64) {
        if self.s.node_lat_div_.is_empty() || self.s.node_lat_div_[un] == 0 {
            return;
        }
        let lat = self.fc.node_lateral[un];
        if lat > 0.0 {
            self.s.node_out_[un] += lat * dt;
        }
    }

    /// incidentPressurized(un)
    pub(crate) fn incident_pressurized(&self, un: usize) -> bool {
        let m = self.m;
        for p in m.node_face_ptr[un]..m.node_face_ptr[un + 1] {
            let uf = m.node_face_idx[p as usize] as usize;
            let c = if m.face_cl[uf] >= 0 { m.face_cl[uf] } else { m.face_cr[uf] };
            if c < 0 {
                continue;
            }
            let uc = c as usize;
            if self.tpa_cell(uc) {
                return true;
            }
            let g = &m.geom[m.cell_geom[uc] as usize];
            if g.is_open == 0 && self.st.cell_h[uc] >= g.y_crown {
                return true;
            }
        }
        false
    }

    /// relaxNodeFluxes(dt, forcing) (no pressurized fold: press_step_ is always false here)
    pub(crate) fn relax_node_fluxes(&mut self, dt: f64) {
        self.s.all_faces_live_ = true;
        let nn = self.m.n_nodes();
        let mut list = std::mem::take(&mut self.s.node_solve_list_);
        list.clear();
        for n in 0..nn {
            list.push(n);
        }
        self.solve_nodes(&list, dt, false);
        self.s.node_solve_list_ = list;
    }

    /// solveNodes(nodes, dt_of, forcing)
    pub(crate) fn solve_nodes(&mut self, nodes: &[i32], dt: f64, tiered: bool) {
        let mut heavy = std::mem::take(&mut self.s.node_heavy_);
        let mut light = std::mem::take(&mut self.s.node_light_);
        heavy.clear();
        light.clear();
        let have_pass = !self.s.node_pass_.is_empty();
        for &n in nodes.iter() {
            let un = n as usize;
            let fixed = self.fc.fixed(un);
            if !fixed && self.algebraic_active(n) && !(have_pass && self.s.node_pass_[un] != 0) {
                heavy.push(n);
            } else {
                light.push(n);
            }
        }
        let dt_of = |s: &crate::hydraulics::fv::explicit_solver::ExplicitFvSolver, n: i32| -> f64 {
            if tiered {
                (1i32 << s.node_tier_[n as usize]) as f64 * dt
            } else {
                dt
            }
        };
        for &n in light.iter() {
            let d = dt_of(self.s, n);
            if self.algebraic_active(n) {
                self.solve_algebraic_node(n, d);
            } else {
                self.relax_one_node(n, d);
            }
        }
        for &n in heavy.iter() {
            let d = dt_of(self.s, n);
            self.solve_algebraic_node(n, d);
        }
        self.s.node_heavy_ = heavy;
        self.s.node_light_ = light;
    }

    /// relaxOneNode(n, dt, forcing)
    pub(crate) fn relax_one_node(&mut self, n: i32, dt: f64) {
        let m = self.m;
        let un = n as usize;
        if m.node_kind[un] == K_NODE_VIRTUAL {
            return;
        }
        if self.fc.fixed(un) {
            return;
        }
        if !(self.st.node_surf_area[un] > 0.0) {
            return;
        }
        if !m.node_dummy_drain.is_empty() && m.node_dummy_drain[un] != 0 {
            return;
        }
        let b = m.node_face_ptr[un];
        let e = m.node_face_ptr[un + 1];
        let q_lat = self.node_lateral(un) + self.s.node_qstruct_[un] + self.node_q_dummy(un);
        let h0 = self.st.node_head[un];
        let mut sum_f = 0.0;
        let mut resist = 0.0;
        for p in b..e {
            let up = p as usize;
            let f = m.node_face_idx[up];
            let uf = f as usize;
            sum_f += m.node_face_sign[up] * self.s.f_mass_[uf];
            self.s.relax_beta_[uf] = 0.0;
            let cell = if m.face_cl[uf] >= 0 { m.face_cl[uf] } else { m.face_cr[uf] };
            if cell < 0 {
                continue;
            }
            let uc = cell as usize;
            let g = &m.geom[m.cell_geom[uc] as usize];
            let hg = h0 - m.face_zb[uf];
            if hg <= k::K_DRY_DEPTH {
                continue;
            }
            let ag = k::area_of_depth(g, hg);
            let tg = k::width_of_depth(g, hg);
            if ag <= k::K_DRY_AREA || tg <= 0.0 {
                continue;
            }
            let beta = (k::K_GRAVITY * ag * tg).sqrt();
            self.s.relax_beta_[uf] = beta;
            resist += beta;
        }
        if resist <= 0.0 {
            return;
        }
        let as_ = self.st.node_surf_area[un];
        let dh = dt * (sum_f + q_lat) / (as_ + dt * resist);
        if dh == 0.0 {
            return;
        }
        for p in b..e {
            let up = p as usize;
            let f = m.node_face_idx[up];
            let uf = f as usize;
            let a = self.s.relax_beta_[uf];
            if a <= 0.0 {
                continue;
            }
            self.s.f_mass_[uf] += -m.node_face_sign[up] * a * dh;
        }
    }

    /// solveAlgebraicNode's residual lambda: the node's net flux imbalance at head h.
    fn alg_residual(&mut self, un: usize, b: i32, e: i32, q_ext: f64, h: f64, sr: &mut AlgSearch) -> f64 {
        let m = self.m;
        self.st.node_head[un] = h;
        sr.h_last = h;
        let mut r = q_ext;
        for p in b..e {
            let up = p as usize;
            let f = m.node_face_idx[up];
            let uf = f as usize;
            if self.face_is_live(f) {
                let cl = m.face_cl[uf];
                let cr = m.face_cr[uf];
                let la = cl >= 0 && self.s.cell_active_[cl as usize] != 0;
                let ra = cr >= 0 && self.s.cell_active_[cr as usize] != 0;
                if la || ra {
                    self.compute_face_flux(f, true);
                }
            }
            r += m.node_face_sign[up] * self.s.f_mass_[uf];
        }
        if r.abs() < sr.r_best.abs() {
            sr.h_best = h;
            sr.r_best = r;
        }
        r
    }

    /// solveAlgebraicNode's finalize lambda: full fluxes at the last evaluated head.
    fn alg_finalize(&mut self, un: usize, b: i32, e: i32, sr: &AlgSearch) {
        let m = self.m;
        if sr.h_last != sr.h_last {
            return;
        }
        let keep = self.st.node_head[un];
        self.st.node_head[un] = sr.h_last;
        for p in b..e {
            let f = m.node_face_idx[p as usize];
            let uf = f as usize;
            if !self.face_is_live(f) {
                continue;
            }
            let cl = m.face_cl[uf];
            let cr = m.face_cr[uf];
            let la = cl >= 0 && self.s.cell_active_[cl as usize] != 0;
            let ra = cr >= 0 && self.s.cell_active_[cr as usize] != 0;
            if la || ra {
                self.compute_face_flux(f, false);
            }
        }
        self.st.node_head[un] = keep;
    }

    /// solveAlgebraicNode's resistAt lambda.
    fn alg_resist_at(&self, b: i32, e: i32, h: f64) -> f64 {
        let m = self.m;
        let mut resist = 0.0;
        for p in b..e {
            let uf = m.node_face_idx[p as usize] as usize;
            let cell = if m.face_cl[uf] >= 0 { m.face_cl[uf] } else { m.face_cr[uf] };
            if cell < 0 {
                continue;
            }
            let g = &m.geom[m.cell_geom[cell as usize] as usize];
            let hg = h - m.face_zb[uf];
            if hg <= k::K_DRY_DEPTH {
                continue;
            }
            let ag = k::area_of_depth(g, hg);
            let tg = k::width_of_depth(g, hg);
            if ag <= k::K_DRY_AREA || tg <= 0.0 {
                continue;
            }
            resist += (k::K_GRAVITY * ag * tg).sqrt();
        }
        resist
    }

    /// solveAlgebraicNode's prescribeDegree1 lambda.
    fn alg_prescribe_degree1(&mut self, un: usize, b: i32, e: i32, q_ext: f64, invert: f64) -> bool {
        let m = self.m;
        if e - b != 1 || q_ext == 0.0 {
            return false;
        }
        let up0 = b as usize;
        let f0 = m.node_face_idx[up0];
        let uf0 = f0 as usize;
        let cell = if m.face_cl[uf0] >= 0 { m.face_cl[uf0] } else { m.face_cr[uf0] };
        if cell < 0 || !self.face_is_live(f0) || self.s.cell_active_[cell as usize] == 0 || m.face_culvert[uf0] >= 0 || m.face_gate[uf0] != 0 {
            return false;
        }
        let eta_c = self.s.cell_eta_[cell as usize];
        let rim = invert + m.node_full_depth[un] + if m.node_can_pond[un] != 0 { 0.0 } else { m.node_sur_depth[un] };
        if m.node_full_depth[un] > 0.0 && eta_c >= rim {
            return false;
        }
        if self.s.deg1_guard_slot_ {
            let ucg = cell as usize;
            let gg = &m.geom[m.cell_geom[ucg] as usize];
            if gg.is_open == 0 && self.st.cell_h[ucg] >= gg.y_full {
                return false;
            }
        }
        let sgn = m.node_face_sign[up0];
        let inp = if m.face_cl[uf0] >= 0 { self.s.f_state_l_[uf0] } else { self.s.f_state_r_[uf0] };
        self.s.f_mass_[uf0] = -sgn * q_ext;
        self.s.f_mom_[uf0] = if inp.a > k::K_DRY_AREA { k::K_GRAVITY * inp.i1 } else { 0.0 };
        self.st.node_head[un] = cmax(invert, eta_c);
        true
    }

    /// solveAlgebraicNode(n, dt, forcing)
    pub(crate) fn solve_algebraic_node(&mut self, n: i32, dt: f64) {
        let m = self.m;
        let un = n as usize;
        if self.fc.fixed(un) {
            return;
        }
        let mut dummy_capacity = 0.0;
        if !m.node_dummy_drain.is_empty() && m.node_dummy_drain[un] != 0 {
            dummy_capacity = f64::INFINITY;
            for s in 0..m.struct_link.len() {
                if m.struct_is_dummy[s] == 0 || m.struct_n1[s] != n {
                    continue;
                }
                let cap = self.fc.link_q_cap[m.struct_link[s] as usize];
                if cap >= 0.0 {
                    dummy_capacity = cmin(dummy_capacity, cap);
                }
            }
            if n < self.fc.n_nodes {
                if let Some(sm) = self.fc.node_source_max.as_ref() {
                    dummy_capacity = cmin(dummy_capacity, sm[un]);
                }
            }
            if !dummy_capacity.is_finite() {
                self.st.node_head[un] = m.node_invert[un];
                return;
            }
        }
        if !self.s.node_pass_.is_empty() && self.s.node_pass_[un] != 0 {
            let b0 = m.node_face_ptr[un];
            let e0 = m.node_face_ptr[un + 1];
            let mut s = 0.0f64;
            let mut lowest = 0.0f64;
            let mut mm: i32 = 0;
            let mut any = false;
            for p in b0..e0 {
                let uf0 = m.node_face_idx[p as usize] as usize;
                let c = if m.face_cl[uf0] >= 0 { m.face_cl[uf0] } else { m.face_cr[uf0] };
                if c < 0 {
                    continue;
                }
                let uc = c as usize;
                let eta = self.s.cell_eta_[uc];
                if !any || eta < lowest {
                    lowest = eta;
                    any = true;
                }
                if !self.tpa_cell(uc) && self.st.cell_h[uc] <= k::K_DRY_DEPTH {
                    continue;
                }
                s += eta;
                mm += 1;
            }
            if mm > 0 {
                self.st.node_head[un] = s / mm as f64;
            } else if any {
                self.st.node_head[un] = lowest;
            }
            return;
        }
        let b = m.node_face_ptr[un];
        let e = m.node_face_ptr[un + 1];
        let invert = m.node_invert[un];
        let q_ext = self.node_lateral(un) + self.s.node_qstruct_[un] + self.s.node_carry_[un] / dt - dummy_capacity;

        let mut sr = AlgSearch { h_best: 0.0, r_best: f64::INFINITY, h_last: f64::NAN };

        let full = m.node_full_depth[un];
        let top;
        if full > 0.0 {
            top = if m.node_can_pond[un] != 0 { full } else { full + m.node_sur_depth[un] };
        } else {
            let mut t = 1.0f64;
            for p in b..e {
                let uf = m.node_face_idx[p as usize] as usize;
                let cell = if m.face_cl[uf] >= 0 { m.face_cl[uf] } else { m.face_cr[uf] };
                if cell < 0 {
                    continue;
                }
                let g: &FvGeometry = &m.geom[m.cell_geom[cell as usize] as usize];
                t = cmax(t, (m.face_zb[uf] - invert) + 2.0 * g.y_full);
            }
            top = t;
        }
        let mut lo = invert;
        let mut hi = invert + top;
        if self.s.tpa_ && m.node_sur_depth[un] > 0.0 {
            let mut any_flagged = false;
            for p in b..e {
                let uf = m.node_face_idx[p as usize] as usize;
                let cell = if m.face_cl[uf] >= 0 { m.face_cl[uf] } else { m.face_cr[uf] };
                if cell >= 0 && self.tpa_cell(cell as usize) {
                    any_flagged = true;
                    break;
                }
            }
            if any_flagged {
                lo = invert - 35.0;
            }
        }

        let h0 = self.st.node_head[un];
        let mut h = cmin(cmax(h0, lo), hi);
        let mut r = if h == h0 { self.alg_residual(un, b, e, q_ext, h0, &mut sr) } else { self.alg_residual(un, b, e, q_ext, h, &mut sr) };
        if r.abs() <= K_ALG_Q_EPS {
            self.st.node_head[un] = h0;
            self.alg_finalize(un, b, e, &sr);
            return;
        }
        let dir = if r > 0.0 { 1.0 } else { -1.0 };
        let res0 = self.alg_resist_at(b, e, h);
        let mut step = cmax(K_ALG_HEAD_TOL, if res0 > 0.0 { r.abs() / res0 } else { K_ALG_HEAD_TOL });
        let mut h_a = h;
        let mut r_a = r;
        let mut bracketed = false;
        for _it in 0..K_ALG_MAX_EXPAND {
            let mut h_b = h_a + dir * step;
            if h_b <= lo {
                h_b = lo;
            }
            if h_b >= hi {
                h_b = hi;
            }
            let r_b = self.alg_residual(un, b, e, q_ext, h_b, &mut sr);
            if r_b.abs() <= K_ALG_Q_EPS {
                self.alg_finalize(un, b, e, &sr);
                return;
            }
            if (r_a > 0.0) != (r_b > 0.0) {
                if dir > 0.0 {
                    lo = h_a;
                    hi = h_b;
                } else {
                    lo = h_b;
                    hi = h_a;
                }
                h = h_b;
                r = r_b;
                bracketed = true;
                break;
            }
            if h_b == lo || h_b == hi {
                self.alg_finalize(un, b, e, &sr);
                if h_b == hi && r_b > 0.0 && self.alg_prescribe_degree1(un, b, e, q_ext, invert) {
                    return;
                }
                return;
            }
            h_a = h_b;
            r_a = r_b;
            step *= 2.0;
        }
        if !bracketed {
            if sr.h_best != h_a {
                let hb = sr.h_best;
                self.alg_residual(un, b, e, q_ext, hb, &mut sr);
            }
            self.alg_finalize(un, b, e, &sr);
            return;
        }
        let root_budget = if self.s.conservative_node_balance_ { 64 } else { K_ALG_MAX_ITERS };
        let head_stop = if self.s.conservative_node_balance_ { 16.0 * f64::EPSILON * cmax(1.0, h.abs()) } else { K_ALG_HEAD_TOL };
        for _it in 0..root_budget {
            let resist = self.alg_resist_at(b, e, h);
            let mut h_new = if resist > 0.0 { h + r / resist } else { 0.5 * (lo + hi) };
            if !(h_new > lo && h_new < hi) {
                h_new = 0.5 * (lo + hi);
            }
            let dh = h_new - h;
            h = h_new;
            r = self.alg_residual(un, b, e, q_ext, h, &mut sr);
            if r > 0.0 {
                lo = h;
            } else {
                hi = h;
            }
            if dh.abs() <= head_stop || r.abs() <= K_ALG_Q_EPS {
                break;
            }
        }
        if sr.r_best.abs() < r.abs() {
            let hb = sr.h_best;
            self.alg_residual(un, b, e, q_ext, hb, &mut sr);
        }
        self.alg_finalize(un, b, e, &sr);
    }

    /// updateCells(dt, forcing) (no species, no unsteady friction: refused before the step)
    pub(crate) fn update_cells(&mut self, dt: f64) {
        let m = self.m;
        let nc = m.n_cells();
        for c in 0..nc {
            let uc = c as usize;
            if self.s.cell_active_[uc] == 0 {
                continue;
            }
            let g = &m.geom[m.cell_geom[uc] as usize];
            let dx = m.cell_dx[uc];
            let inv_dx = 1.0 / dx;
            let faces = [m.cell_face0[uc], m.cell_face1[uc]];
            let sides = [m.cell_side0[uc], m.cell_side1[uc]];
            let mut d_a = 0.0;
            let mut d_q = 0.0;
            let mut k_loss = 0.0;
            let mut q_thru = 0.0;
            for e in 0..2 {
                let uf = faces[e] as usize;
                let fa = self.s.f_mass_[uf];
                let fq = self.s.f_mom_[uf];
                if sides[e] == 0 {
                    d_a -= fa;
                    d_q -= m.face_dir_l[uf] as f64 * (fq + self.s.f_corr_l_[uf]);
                    q_thru += m.face_dir_l[uf] as f64 * fa;
                } else {
                    d_a += fa;
                    d_q += m.face_dir_r[uf] as f64 * (fq + self.s.f_corr_r_[uf]);
                    q_thru += m.face_dir_r[uf] as f64 * fa;
                }
                if m.face_node[uf] >= 0 {
                    k_loss += if m.face_cl[uf] < 0 { m.conduit_loss_inlet[m.cell_conduit[uc] as usize] } else { m.conduit_loss_outlet[m.cell_conduit[uc] as usize] };
                }
            }
            self.s.cell_q_int_[uc] += 0.5 * q_thru * dt;
            if self.s.opts_.order >= 2 && self.s.cell_ho2_[uc] != 0 {
                let half = 0.5 * m.cell_dx[uc];
                let zp = m.cell_zb[uc] + m.cell_dzdx[uc] * half;
                let zm = m.cell_zb[uc] - m.cell_dzdx[uc] * half;
                let eta_c = self.s.cell_eta_[uc];
                let hp = cmax(0.0, eta_c - zp);
                let hm = cmax(0.0, eta_c - zm);
                let i1p = if hp > 0.0 { k::i1_of_depth(g, hp, k::area_of_depth(g, hp)) } else { 0.0 };
                let i1m = if hm > 0.0 { k::i1_of_depth(g, hm, k::area_of_depth(g, hm)) } else { 0.0 };
                d_q += k::K_GRAVITY * (i1p - i1m);
            }
            let a_old = self.st.cell_a[uc];
            let mut a_new = a_old + dt * d_a * inv_dx;
            let mut q_new = self.st.cell_q[uc] + dt * d_q * inv_dx;
            if !self.s.cell_qlat_.is_empty() && self.s.cell_qlat_[uc] != 0.0 {
                let before_lateral = a_new;
                a_new = self.inject_diverted_lateral(uc, a_new, dt, inv_dx);
                if self.s.water_tracking_ {
                    self.s.water_.cell_lateral[uc] += (cmax(0.0, a_new) - cmax(0.0, before_lateral)) / inv_dx;
                }
            }
            // (species: n_species == 0)
            {
                let cr = m.cell_conduit[uc];
                if cr >= 0 {
                    let before_loss = a_new;
                    a_new -= dt * self.fc.conduit_loss[cr as usize];
                    if self.s.water_tracking_ {
                        self.s.water_.cell_loss[uc] += (cmax(0.0, before_loss) - cmax(0.0, a_new)) / inv_dx;
                    }
                }
            }
            if a_new < 0.0 {
                a_new = 0.0;
            }
            let tpa_cell = self.tpa_cell(uc);
            let h_new = if tpa_cell { k::tpa_depth_of_area(g, a_new) } else { k::depth_of_area(g, a_new) };
            self.st.cell_a[uc] = a_new;
            self.st.cell_h[uc] = h_new;
            self.s.cell_eta_[uc] = m.cell_zb[uc] + h_new;
            let radius = self.cache_closure_r(uc, g, h_new);
            if (!tpa_cell && h_new <= k::K_DRY_DEPTH) || (tpa_cell && a_new <= k::K_DRY_AREA) {
                self.st.cell_q[uc] = 0.0;
                self.s.cell_u_[uc] = 0.0;
                continue;
            }
            let u = q_new / a_new;
            let cdi = m.cell_conduit[uc] as usize;
            q_new = self.friction_for(g, m.conduit_roughness[cdi], m.conduit_rough_factor[cdi], q_new, u, if tpa_cell { g.y_full } else { h_new }, dt, radius);
            if k_loss > 0.0 {
                q_new = k::local_loss_update(q_new, u, k_loss, dx, dt);
            }
            self.st.cell_q[uc] = q_new;
            self.s.cell_u_[uc] = q_new / a_new;
        }
    }

    /// refreshStructFlows(forcing)
    pub(crate) fn refresh_struct_flows(&mut self) {
        let m = self.m;
        let nn = m.n_nodes() as usize;
        self.s.node_qstruct_.clear();
        self.s.node_qstruct_.resize(nn, 0.0);
        for us in 0..m.struct_link.len() {
            if !m.struct_is_dummy.is_empty() && m.struct_is_dummy[us] != 0 {
                continue;
            }
            let q = self.fc.structure_flow[m.struct_link[us] as usize];
            if q == 0.0 {
                continue;
            }
            self.s.node_qstruct_[m.struct_n1[us] as usize] -= q;
            self.s.node_qstruct_[m.struct_n2[us] as usize] += q;
        }
        if self.s.cell_qlat_.len() != m.n_cells() as usize {
            self.s.cell_qlat_.clear();
            self.s.cell_qlat_.resize(m.n_cells() as usize, 0.0);
        } else {
            self.s.cell_qlat_.iter_mut().for_each(|x| *x = 0.0);
        }
        if self.s.node_lat_div_.len() != nn {
            self.s.node_lat_div_.clear();
            self.s.node_lat_div_.resize(nn, 0);
        } else {
            self.s.node_lat_div_.iter_mut().for_each(|x| *x = 0);
        }
        for un in 0..nn {
            let drains = !m.node_dummy_drain.is_empty() && m.node_dummy_drain[un] != 0;
            let extracting = self.fc.node_withdrawal && self.s.withdrawal_rates_[un] > 0.0 && self.fc.node_lateral[un] < 0.0;
            let clean = self.s.node_pass_static_[un] != 0 && !drains && self.s.node_qstruct_[un] == 0.0 && self.s.node_carry_[un] == 0.0 && !self.incident_pressurized(un) && !extracting;
            let lat = self.fc.node_lateral[un];
            if clean && lat != 0.0 {
                let mut credited = false;
                for p in m.node_face_ptr[un]..m.node_face_ptr[un + 1] {
                    let uf = m.node_face_idx[p as usize] as usize;
                    let c = if m.face_cl[uf] >= 0 { m.face_cl[uf] } else { m.face_cr[uf] };
                    if c >= 0 {
                        self.s.cell_qlat_[c as usize] += 0.5 * lat;
                        credited = true;
                        if !self.s.cell_active_.is_empty() && self.s.cell_active_[c as usize] == 0 {
                            self.s.lists_valid_ = false;
                        }
                    }
                }
                if credited {
                    self.s.node_lat_div_[un] = 1;
                }
            } else if lat != 0.0 && m.node_kind[un] == K_NODE_VIRTUAL && !m.node_vj_face.is_empty() && m.node_vj_face[un] >= 0 {
                let uf = m.node_vj_face[un] as usize;
                let cl = m.face_cl[uf];
                let cr = m.face_cr[uf];
                if cl >= 0 && cr >= 0 {
                    self.s.cell_qlat_[cl as usize] += 0.5 * lat;
                    self.s.cell_qlat_[cr as usize] += 0.5 * lat;
                    self.s.node_lat_div_[un] = 1;
                    if !self.s.cell_active_.is_empty() && (self.s.cell_active_[cl as usize] == 0 || self.s.cell_active_[cr as usize] == 0) {
                        self.s.lists_valid_ = false;
                    }
                }
            }
            self.s.node_pass_[un] = if clean && (lat == 0.0 || self.s.node_lat_div_[un] != 0) { 1 } else { 0 };
            self.s.node_pub_[un] = if self.s.node_pub_static_[un] != 0 && self.s.node_qstruct_[un] == 0.0 && !drains { 1 } else { 0 };
        }
    }

    /// refreshDummyFlows(dt, forcing, arrived, book_receipts)
    pub(crate) fn refresh_dummy_flows(&mut self, dt: f64, arrived: Arrived, book_receipts: bool) {
        let m = self.m;
        if m.struct_is_dummy.is_empty() || dt <= 0.0 {
            return;
        }
        self.s.node_qdummy_.iter_mut().for_each(|x| *x = 0.0);
        for us in 0..m.struct_link.len() {
            if m.struct_is_dummy[us] == 0 {
                continue;
            }
            let u1 = m.struct_n1[us] as usize;
            let u2 = m.struct_n2[us] as usize;
            let arr = match arrived {
                Arrived::CellScaled(out_cell, dtc) => {
                    let mut volume = 0.0;
                    for p in m.node_face_ptr[u1]..m.node_face_ptr[u1 + 1] {
                        let f = m.node_face_idx[p as usize] as usize;
                        let c = if m.face_cl[f] >= 0 { m.face_cl[f] } else { m.face_cr[f] };
                        let q = m.node_face_sign[p as usize] * self.s.f_mass_[f];
                        volume += q * dtc * if q > 0.0 && c >= 0 { out_cell[c as usize] } else { 1.0 };
                    }
                    volume
                }
                Arrived::AccNvol => self.s.acc_nvol_[u1],
                Arrived::FaceSum(dtf) => {
                    let mut sum = 0.0;
                    for p in m.node_face_ptr[u1]..m.node_face_ptr[u1 + 1] {
                        let up = p as usize;
                        sum += m.node_face_sign[up] * self.s.f_mass_[m.node_face_idx[up] as usize];
                    }
                    sum * dtf
                }
            };
            let mut q = arr / dt + self.node_lateral(u1) + self.s.node_qstruct_[u1];
            if q < 0.0 {
                q = 0.0;
            }
            {
                let cap = self.fc.link_q_cap[m.struct_link[us] as usize];
                if cap >= 0.0 && q > cap {
                    q = cap;
                }
            }
            if (u1 as i32) < self.fc.n_nodes {
                if let Some(sm) = self.fc.node_source_max.as_ref() {
                    q = cmin(q, sm[u1]);
                }
            }
            self.s.node_qdummy_[u1] -= q;
            self.s.node_qdummy_[u2] += q;
            if book_receipts {
                self.s.dummy_vol_[us] += q * dt;
                if self.s.water_tracking_ {
                    self.s.water_.structure(us, q * dt);
                }
            }
        }
    }

    /// capWithdrawal(node, balance, dt, forcing)
    pub(crate) fn cap_withdrawal(&mut self, node: usize, mut balance: f64, dt: f64) -> f64 {
        if balance < 0.0 && self.fc.node_withdrawal {
            let refused = cmin(-balance, cmax(0.0, self.s.withdrawal_rates_[node]) * dt);
            self.s.unmet_withdrawal_[node] += refused;
            balance += refused;
        }
        balance
    }

    /// bookFaceWater(f, dt)
    pub(crate) fn book_face_water(&mut self, f: i32, dt: f64) {
        let m = self.m;
        let uf = f as usize;
        let v = self.s.f_mass_[uf] * dt;
        self.s.water_.face(uf, v);
        let mut node = m.face_node[uf];
        if node < 0 {
            node = m.face_vj_node[uf];
        }
        if node < 0 || self.s.node_lat_div_[node as usize] == 0 {
            return;
        }
        let volume = 0.5 * self.fc.node_lateral[node as usize] * dt;
        if m.face_cl[uf] >= 0 {
            self.s.water_.face_lateral_left[uf] += volume;
        }
        if m.face_cr[uf] >= 0 {
            self.s.water_.face_lateral_right[uf] += volume;
        }
    }

    /// bookSpliceReceipt(vj, flux, dt)
    pub(crate) fn book_splice_receipt(&mut self, vj: usize, flux: f64, dt: f64) {
        let mut h = 0.0;
        if self.s.node_lat_div_[vj] != 0 {
            h = 0.5 * self.fc.node_lateral[vj];
        }
        let left = flux - h;
        let right = -(flux + h);
        self.s.node_in_[vj] += (cmax(left, 0.0) + cmax(right, 0.0)) * dt;
        self.s.node_out_[vj] += (cmax(-left, 0.0) + cmax(-right, 0.0)) * dt;
    }

    /// bookStructureWater(dt, forcing)
    pub(crate) fn book_structure_water(&mut self, dt: f64) {
        let m = self.m;
        for s in 0..m.struct_link.len() {
            if !m.struct_is_dummy.is_empty() && m.struct_is_dummy[s] != 0 {
                continue;
            }
            let v = self.fc.structure_flow[m.struct_link[s] as usize] * dt;
            self.s.water_.structure(s, v);
        }
    }

    /// updateNodes(dt, forcing)
    pub(crate) fn update_nodes(&mut self, dt: f64) {
        let m = self.m;
        if self.s.water_tracking_ {
            for ai in 0..self.s.active_faces_.len() {
                let f = self.s.active_faces_[ai];
                self.book_face_water(f, dt);
            }
            self.book_structure_water(dt);
        }
        let nn = m.n_nodes();
        self.refresh_dummy_flows(dt, Arrived::FaceSum(dt), true);
        self.book_lateral_spills(None);

        for n in 0..nn {
            let un = n as usize;
            if m.node_kind[un] == K_NODE_VIRTUAL {
                let face = m.node_vj_face[un];
                if face >= 0 {
                    let fm = self.s.f_mass_[face as usize];
                    self.book_splice_receipt(un, fm, dt);
                }
                continue;
            }
            let mut sum_faces = 0.0;
            let b = m.node_face_ptr[un];
            let e = m.node_face_ptr[un + 1];
            for p in b..e {
                let up = p as usize;
                sum_faces += m.node_face_sign[up] * self.s.f_mass_[m.node_face_idx[up] as usize];
            }
            self.s.node_exch_[un] += sum_faces * dt;
            for p2 in b..e {
                let up2 = p2 as usize;
                let q = m.node_face_sign[up2] * self.s.f_mass_[m.node_face_idx[up2] as usize];
                if q > 0.0 {
                    self.s.node_in_[un] += q * dt;
                } else {
                    self.s.node_out_[un] -= q * dt;
                }
            }
            self.book_diverted_lateral(un, dt);

            if self.fc.fixed(un) {
                self.st.node_head[un] = self.fc.node_fixed_head[un];
                continue;
            }
            let q_lat = self.node_lateral(un) + self.s.node_qstruct_[un] + self.node_q_dummy(un);
            if self.algebraic_active(n) {
                let bal = self.cap_withdrawal(un, self.s.node_carry_[un] + dt * (sum_faces + q_lat), dt);
                self.settle_algebraic_node(n, bal);
                continue;
            }
            let v_prev = self.st.node_volume[un];
            let mut vol = self.cap_withdrawal(un, v_prev + dt * (sum_faces + q_lat), dt);
            if vol < 0.0 {
                if self.s.water_tracking_ {
                    self.s.water_.node_deficit[un] -= vol;
                }
                vol = 0.0;
            }
            let mut depth = self.node_depth_from_volume(n, vol);
            self.apply_node_capacity(n, v_prev, &mut vol, &mut depth);
            self.st.node_volume[un] = vol;
            self.st.node_head[un] = m.node_invert[un] + depth;
        }
    }

    /// settleAlgebraicNode(n, carry)
    pub(crate) fn settle_algebraic_node(&mut self, n: i32, mut carry: f64) {
        let m = self.m;
        let un = n as usize;
        let invert = m.node_invert[un];
        let depth = cmax(0.0, self.st.node_head[un] - invert);
        let full = m.node_full_depth[un];
        let published = self.st.node_volume[un];
        let ledger = cmax(0.0, self.s.node_carry_[un]);
        if published > ledger {
            carry += published - ledger;
        }
        if carry > 0.0 && full > 0.0 {
            let eps = 1.0e-9;
            if m.node_can_pond[un] != 0 {
                if depth >= full - eps {
                    self.st.node_volume[un] = self.s.node_vfull_[un] + carry;
                    let d = self.node_depth_from_volume(n, self.st.node_volume[un]);
                    self.st.node_head[un] = invert + d;
                    self.s.node_carry_[un] = 0.0;
                    return;
                }
            } else if depth >= full + m.node_sur_depth[un] - eps {
                self.s.flood_vol_[un] += carry;
                self.s.node_carry_[un] = 0.0;
                self.st.node_volume[un] = 0.0;
                return;
            }
        }
        if !self.s.node_pass_.is_empty() && self.s.node_pass_[un] != 0 && carry != 0.0 && !self.incident_pressurized(un) {
            let mut left = carry;
            let b = m.node_face_ptr[un];
            let e = m.node_face_ptr[un + 1];
            let deg = if e > b { e - b } else { 1 };
            for p in b..e {
                let uf = m.node_face_idx[p as usize] as usize;
                let c = if m.face_cl[uf] >= 0 { m.face_cl[uf] } else { m.face_cr[uf] };
                if c < 0 {
                    continue;
                }
                let uc = c as usize;
                let dx = m.cell_dx[uc];
                let mut a_new = self.st.cell_a[uc] + (carry / deg as f64) / dx;
                if a_new < 0.0 {
                    a_new = 0.0;
                }
                if self.s.water_tracking_ {
                    let add = (a_new - self.st.cell_a[uc]) * dx;
                    if m.face_cl[uf] == c {
                        self.s.water_.face_lateral_left[uf] += add;
                    } else {
                        self.s.water_.face_lateral_right[uf] += add;
                    }
                }
                left -= (a_new - self.st.cell_a[uc]) * dx;
                let g = &m.geom[m.cell_geom[uc] as usize];
                let h_new = k::depth_of_area(g, a_new);
                self.st.cell_a[uc] = a_new;
                self.st.cell_h[uc] = h_new;
                self.s.cell_eta_[uc] = m.cell_zb[uc] + h_new;
                self.s.cache_closure(self.st, uc, g, h_new);
                if h_new <= k::K_DRY_DEPTH {
                    self.st.cell_q[uc] = 0.0;
                    self.s.cell_u_[uc] = 0.0;
                } else {
                    self.s.cell_u_[uc] = self.st.cell_q[uc] / a_new;
                }
            }
            self.s.node_carry_[un] = if left.abs() > 1.0e-12 { left } else { 0.0 };
            self.st.node_volume[un] = cmax(0.0, self.s.node_carry_[un]);
            return;
        }
        self.s.node_carry_[un] = carry;
        self.st.node_volume[un] = cmax(0.0, carry);
    }

    /// injectDivertedLateral(uc, a_new, dt, inv_dx)
    pub(crate) fn inject_diverted_lateral(&mut self, uc: usize, a_new: f64, dt: f64, inv_dx: f64) -> f64 {
        let m = self.m;
        let add = dt * self.s.cell_qlat_[uc] * inv_dx;
        if add <= 0.0 || self.s.face_a_rim_.is_empty() {
            if self.s.water_tracking_ && a_new + add < 0.0 && add < 0.0 {
                let mut requested = 0.0;
                let faces = [m.cell_face0[uc], m.cell_face1[uc]];
                for &f in faces.iter() {
                    let mut n = m.face_node[f as usize];
                    if n < 0 {
                        n = m.face_vj_node[f as usize];
                    }
                    if n >= 0 && self.s.node_lat_div_[n as usize] != 0 {
                        requested += cmax(0.0, -0.5 * self.fc.node_lateral[n as usize] * dt);
                    }
                }
                let refused = cmin(requested, -(a_new + add) / inv_dx);
                if requested > 0.0 {
                    for &f in faces.iter() {
                        let mut n = m.face_node[f as usize];
                        if n < 0 {
                            n = m.face_vj_node[f as usize];
                        }
                        if n < 0 || self.s.node_lat_div_[n as usize] == 0 {
                            continue;
                        }
                        let share = cmax(0.0, -0.5 * self.fc.node_lateral[n as usize] * dt) / requested;
                        if m.face_cl[f as usize] == uc as i32 {
                            self.s.water_.face_lateral_left[f as usize] += refused * share;
                        } else {
                            self.s.water_.face_lateral_right[f as usize] += refused * share;
                        }
                    }
                }
            }
            return a_new + add;
        }
        let faces = [m.cell_face0[uc], m.cell_face1[uc]];
        let mut f_rim: i32 = -1;
        let mut a_rim = f64::INFINITY;
        for &f in faces.iter() {
            let uf = f as usize;
            if self.s.face_a_rim_[uf] < a_rim {
                a_rim = self.s.face_a_rim_[uf];
                f_rim = f;
            }
        }
        if f_rim < 0 || a_new + add <= a_rim {
            return a_new + add;
        }
        let room = cmax(0.0, a_rim - a_new);
        self.s.face_lat_spill_[f_rim as usize] += (add - room) / inv_dx;
        if self.s.water_tracking_ {
            self.s.water_.face_spill[f_rim as usize] += (add - room) / inv_dx;
        }
        a_new + room
    }

    /// bookLateralSpills(nodes)
    pub(crate) fn book_lateral_spills(&mut self, nodes: Option<&[i32]>) {
        if self.s.face_lat_spill_.is_empty() {
            return;
        }
        let m = self.m;
        let book = |this: &mut Self, n: i32| {
            let un = n as usize;
            if this.s.node_pass_static_[un] == 0 {
                return;
            }
            let mut spill = 0.0;
            for p in m.node_face_ptr[un]..m.node_face_ptr[un + 1] {
                let uf = m.node_face_idx[p as usize] as usize;
                spill += this.s.face_lat_spill_[uf];
                this.s.face_lat_spill_[uf] = 0.0;
            }
            if spill <= 0.0 {
                return;
            }
            if m.node_can_pond[un] != 0 {
                this.st.node_volume[un] = cmax(this.st.node_volume[un], this.s.node_vfull_[un]) + spill;
            } else {
                this.s.flood_vol_[un] += spill;
            }
        };
        match nodes {
            Some(list) => {
                for &n in list.iter() {
                    book(self, n);
                }
            }
            None => {
                for n in 0..m.n_nodes() {
                    book(self, n);
                }
            }
        }
    }

    /// frictionFor(g, roughness, rough_factor, q, u, h, dt, r)
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn friction_for(&self, g: &FvGeometry, roughness: f64, rough_factor: f64, q: f64, u: f64, h: f64, dt: f64, r: f64) -> f64 {
        if g.xs.type_ == XSectShape::ForceMain as i32 && h >= g.y_full {
            let absu = u.abs();
            if absu <= 0.0 {
                return q;
            }
            let sf = if roughness < 1.0 { fric_slope_dw(u, r, roughness) } else { fric_slope_hw(u, r, roughness) };
            return q / (1.0 + dt * k::K_GRAVITY * sf / absu);
        }
        k::friction_update(q, u, r, dt, rough_factor)
    }

    /// applyNodeCapacity(node, v_prev, vol, depth)
    pub(crate) fn apply_node_capacity(&mut self, node: i32, v_prev: f64, vol: &mut f64, depth: &mut f64) {
        let m = self.m;
        let un = node as usize;
        let full = m.node_full_depth[un];
        if full <= 0.0 || *depth <= full {
            return;
        }
        let v_full = self.node_volume_from_depth(node, full);
        if m.node_can_pond[un] != 0 {
            self.s.flood_vol_[un] += cmax(*vol - cmax(v_prev, v_full), 0.0);
            return;
        }
        let y_max = full + m.node_sur_depth[un];
        if *depth <= y_max {
            return;
        }
        let v_max = self.node_volume_from_depth(node, y_max);
        self.s.flood_vol_[un] += *vol - v_max;
        *vol = v_max;
        *depth = y_max;
    }
}
