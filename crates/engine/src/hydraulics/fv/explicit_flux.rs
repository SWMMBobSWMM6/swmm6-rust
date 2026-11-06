// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/fv/ExplicitFvSolver.cpp and
// transport/fvkernels/SpeciesTransportKernels.hpp (limitSlope) (HydroCouple/openswmm.engine,
// branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the
// Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The explicit FV step's face side: reconstructState (MUSCL-Hancock slopes), faceSide,
// computeFluxes, computeFaceFlux (HLL/HLLC with gates, culvert inlet control and the
// surface-outfall source cap), algebraicAvailableWater, limitPositivity, and cacheClosure
// with its radius output.
//
// Differences in form, never in result:
//   * faceSide's output references (out, i1_raw, z_side) are a returned tuple.
//   * The function-static out_cell / out_node of limitPositivity are local vectors: they
//     are assigned in full on every call, so nothing carries over.

use crate::core::cmath::{cmax, cmin};
use crate::core::options::fv::Limiter;
use crate::hydraulics::fv::explicit_solver::{FaceFlux, FaceState};
use crate::hydraulics::fv::explicit_step::{Arrived, FvStep};
use crate::hydraulics::fv::kernels as k;
use crate::hydraulics::fv::FvGeometry;
use crate::hydraulics::hyd_closure_kernels::culvert_inflow;

/// transport::fvkernels::limitSlope(a, b, lim)
pub fn limit_slope(a: f64, b: f64, lim: Limiter) -> f64 {
    if a * b <= 0.0 {
        return 0.0;
    }
    match lim {
        Limiter::Vanleer => 2.0 * a * b / (a + b),
        Limiter::Superbee => {
            let s = if a > 0.0 { 1.0 } else { -1.0 };
            s * cmax(cmin(2.0 * a.abs(), b.abs()), cmin(a.abs(), 2.0 * b.abs()))
        }
        Limiter::Minmod => {
            if a.abs() < b.abs() {
                a
            } else {
                b
            }
        }
    }
}

impl<'a> FvStep<'a> {
    /// cacheClosure(uc, g, h, &radius): the cell's closure cache, returning the radius.
    pub(crate) fn cache_closure_r(&mut self, uc: usize, g: &FvGeometry, h: f64) -> f64 {
        if h > k::K_DRY_DEPTH && !self.tpa_cell(uc) {
            let (e, radius) = k::closure_eval_r(&g.closure_tbl, h);
            self.s.cell_ah_[uc] = e.a;
            self.s.cell_t_[uc] = e.t;
            self.s.cell_i1_[uc] = e.i1;
            radius
        } else {
            self.s.cell_ah_[uc] = 0.0;
            self.s.cell_t_[uc] = 0.0;
            self.s.cell_i1_[uc] = 0.0;
            k::hyd_rad_of_depth(g, if self.tpa_cell(uc) { g.y_full } else { h })
        }
    }

    /// reconstructState()
    pub(crate) fn reconstruct_state(&mut self) {
        let m = self.m;
        let nc = m.n_cells();
        if self.s.opts_.order < 2 {
            self.s.cell_eta_slope_.iter_mut().for_each(|x| *x = 0.0);
            self.s.cell_u_slope_.iter_mut().for_each(|x| *x = 0.0);
            self.s.cell_ho2_.iter_mut().for_each(|x| *x = 0);
            return;
        }
        for c in 0..nc {
            let uc = c as usize;
            let h = self.st.cell_h[uc];
            let fall = m.cell_dzdx[uc].abs() * m.cell_dx[uc];
            self.s.cell_ho2_[uc] = if h > k::K_DRY_DEPTH && fall < 0.5 * h && !self.tpa_cell(uc) { 1 } else { 0 };
        }
        let lim = self.s.opts_.limiter;
        for ch in 0..m.n_chains() {
            let b = m.chain_ptr[ch as usize];
            let e = m.chain_ptr[ch as usize + 1];
            for i in b..e {
                let c = m.chain_cells[i as usize];
                let uc = c as usize;
                if i == b || i == e - 1 || self.s.cell_active_[uc] == 0 || self.s.cell_ho2_[uc] == 0 {
                    self.s.cell_eta_slope_[uc] = 0.0;
                    self.s.cell_u_slope_[uc] = 0.0;
                    continue;
                }
                let cm = m.chain_cells[(i - 1) as usize];
                let cp = m.chain_cells[(i + 1) as usize];
                let um = cm as usize;
                let up = cp as usize;
                if self.st.cell_h[um] <= k::K_DRY_DEPTH || self.st.cell_h[up] <= k::K_DRY_DEPTH || self.st.cell_h[uc] <= k::K_DRY_DEPTH {
                    self.s.cell_eta_slope_[uc] = 0.0;
                    self.s.cell_u_slope_[uc] = 0.0;
                    continue;
                }
                let dm = 0.5 * (m.cell_dx[um] + m.cell_dx[uc]);
                let dp = 0.5 * (m.cell_dx[uc] + m.cell_dx[up]);
                let dirc = m.chain_dir[i as usize] as f64;
                let dirm = m.chain_dir[(i - 1) as usize] as f64;
                let dirp = m.chain_dir[(i + 1) as usize] as f64;
                let d_eta_m = self.s.cell_eta_[uc] - self.s.cell_eta_[um];
                let d_eta_p = self.s.cell_eta_[up] - self.s.cell_eta_[uc];
                if d_eta_m.abs() < k::K_ETA_DEADBAND && d_eta_p.abs() < k::K_ETA_DEADBAND {
                    self.s.cell_eta_slope_[uc] = 0.0;
                } else {
                    self.s.cell_eta_slope_[uc] = limit_slope(d_eta_m / dm, d_eta_p / dp, lim) * dirc;
                }
                let um_axis = self.s.cell_u_[um] * dirm * dirc;
                let up_axis = self.s.cell_u_[up] * dirp * dirc;
                let gm_u = (self.s.cell_u_[uc] - um_axis) / dm;
                let gp_u = (up_axis - self.s.cell_u_[uc]) / dp;
                self.s.cell_u_slope_[uc] = limit_slope(gm_u, gp_u, lim);
            }
        }
    }

    /// faceSide(face, cell, node, zstar, dir, u_interior, out, i1_raw, z_side, measure_only),
    /// returned as (out, i1_raw, z_side).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn face_side(&self, face: i32, cell: i32, node: i32, zstar: f64, dir: i32, u_interior: f64, measure_only: bool) -> (FaceState, f64, f64) {
        let m = self.m;
        let uf = face as usize;
        let g: &FvGeometry;
        let mut eta;
        let mut h_raw;
        let u;
        let mut z_side;
        let mut centred = false;
        let mut uc_cached = 0usize;
        let mut node_head_ghost = false;
        let gf = &m.geom[m.face_geom[uf] as usize];

        if cell >= 0 {
            let uc = cell as usize;
            g = &m.geom[m.cell_geom[uc] as usize];
            if self.s.opts_.order < 2 || self.s.cell_ho2_[uc] == 0 {
                eta = self.s.cell_eta_[uc];
                h_raw = self.st.cell_h[uc];
                u = dir as f64 * self.s.cell_u_[uc];
                z_side = m.cell_zb[uc];
                centred = true;
                uc_cached = uc;
            } else {
                let is_left = m.face_cl[uf] == cell;
                let sgn = if is_left { m.face_dir_l[uf] as f64 } else { -(m.face_dir_r[uf] as f64) };
                let half = sgn * 0.5 * m.cell_dx[uc];
                z_side = m.cell_zb[uc] + m.cell_dzdx[uc] * half;
                eta = self.s.cell_eta_[uc] + self.s.cell_eta_slope_[uc] * half;
                h_raw = cmax(0.0, eta - z_side);
                u = dir as f64 * (self.s.cell_u_[uc] + self.s.cell_u_slope_[uc] * half);
            }
        } else {
            let other = if m.face_cl[uf] >= 0 { m.face_cl[uf] } else { m.face_cr[uf] };
            let uo = other as usize;
            g = &m.geom[m.cell_geom[uo] as usize];
            if node >= 0 {
                let und = node as usize;
                eta = self.st.node_head[und];
                h_raw = cmax(0.0, eta - m.face_zb[uf]);
                u = u_interior;
                z_side = m.face_zb[uf];
                node_head_ghost = true;
                if !self.s.node_pass_.is_empty() && self.s.node_pass_[und] != 0 && self.algebraic_active(node) {
                    let nb2 = m.node_face_ptr[und];
                    let ne2 = m.node_face_ptr[und + 1];
                    for p2 in nb2..ne2 {
                        let of = m.node_face_idx[p2 as usize];
                        if of == face {
                            continue;
                        }
                        let uof = of as usize;
                        let oc = if m.face_cl[uof] >= 0 { m.face_cl[uof] } else { m.face_cr[uof] };
                        if oc >= 0 {
                            let uoc = oc as usize;
                            eta = self.s.cell_eta_[uoc];
                            h_raw = self.st.cell_h[uoc];
                            z_side = m.cell_zb[uoc];
                            node_head_ghost = false;
                        }
                    }
                }
            } else {
                eta = self.s.cell_eta_[uo];
                h_raw = self.st.cell_h[uo];
                u = -u_interior;
                z_side = m.face_zb[uf];
            }
        }

        if measure_only {
            return (FaceState::default(), 0.0, z_side);
        }

        let mut tpa_side = false;
        if self.s.tpa_ {
            if cell >= 0 {
                tpa_side = self.tpa_cell(cell as usize);
            } else {
                let other2 = if m.face_cl[uf] >= 0 { m.face_cl[uf] } else { m.face_cr[uf] };
                if other2 >= 0 {
                    let interior_tpa = self.tpa_cell(other2 as usize);
                    if node >= 0 {
                        tpa_side = interior_tpa && m.node_sur_depth[node as usize] > 0.0;
                    } else {
                        tpa_side = interior_tpa;
                    }
                }
            }
        }
        if tpa_side {
            let mut out = FaceState::default();
            let i1_raw = k::tpa_i1_of_depth(g, eta - z_side);
            let h_star_t = eta - zstar;
            let mut a_t = k::tpa_area_of_depth(gf, h_star_t);
            if a_t < k::K_DRY_AREA {
                a_t = k::K_DRY_AREA;
            }
            out.a = a_t;
            out.u = u;
            out.q = a_t * u;
            out.c = k::celerity(a_t, gf.t_slot);
            out.i1 = k::tpa_i1_of_depth(gf, h_star_t);
            out.press = if self.s.mixed_wave_ && self.s.mw_rh_ { 1 } else { 0 };
            return (out, i1_raw, z_side);
        }

        let i1_raw = if cell >= 0 && h_raw > k::K_DRY_DEPTH {
            if centred {
                self.s.cell_i1_[uc_cached]
            } else {
                k::closure_i1(&g.closure_tbl, h_raw)
            }
        } else {
            0.0
        };

        let h_star = cmax(0.0, eta - zstar);
        if h_star <= k::K_DRY_DEPTH {
            return (FaceState::default(), i1_raw, z_side);
        }
        let finish = |o: &mut FaceState| {
            if !self.s.mixed_wave_ {
                return;
            }
            o.press = if self.s.mw_rh_ && gf.is_open == 0 && h_star >= gf.y_full { 1 } else { 0 };
            if node_head_ghost && node >= 0 && m.node_sur_depth[node as usize] <= 0.0 {
                o.press = if self.s.mw_rh_ { 2 } else { 0 };
            }
        };
        let mut out = FaceState::default();
        if centred && std::ptr::eq(gf, g) && h_star == h_raw {
            out.a = self.s.cell_ah_[uc_cached];
            out.u = u;
            out.q = out.a * u;
            out.c = k::celerity(out.a, self.s.cell_t_[uc_cached]);
            out.i1 = self.s.cell_i1_[uc_cached];
            finish(&mut out);
            return (out, i1_raw, z_side);
        }
        let e = k::closure_eval(&gf.closure_tbl, h_star);
        out.a = e.a;
        out.u = u;
        out.q = e.a * u;
        out.c = k::celerity(e.a, e.t);
        out.i1 = e.i1;
        finish(&mut out);
        (out, i1_raw, z_side)
    }

    /// computeFluxes()
    pub(crate) fn compute_fluxes(&mut self) {
        let n_act = self.s.active_faces_.len();
        for a in 0..n_act {
            let f = self.s.active_faces_[a];
            self.compute_face_flux(f, false);
        }
        self.s.total_flux_ += n_act as i64;
    }

    /// computeFaceFlux(f, mass_only)
    pub(crate) fn compute_face_flux(&mut self, f: i32, mass_only: bool) {
        let m = self.m;
        let uf = f as usize;
        let cl = m.face_cl[uf];
        let cr = m.face_cr[uf];
        let nd = m.face_node[uf];
        let (zl, zr);
        if cl >= 0 && cr >= 0 && (self.s.opts_.order < 2 || (self.s.cell_ho2_[cl as usize] == 0 && self.s.cell_ho2_[cr as usize] == 0)) {
            zl = m.cell_zb[cl as usize];
            zr = m.cell_zb[cr as usize];
        } else {
            zl = self.face_side(f, cl, nd, 0.0, m.face_dir_l[uf] as i32, 0.0, true).2;
            zr = self.face_side(f, cr, nd, 0.0, m.face_dir_r[uf] as i32, 0.0, true).2;
        }
        let zstar = cmax(zl, zr);
        let interior = if cl >= 0 { cl } else { cr };
        let mut u_int = if interior >= 0 { (if cl >= 0 { m.face_dir_l[uf] } else { m.face_dir_r[uf] }) as f64 * self.s.cell_u_[interior as usize] } else { 0.0 };
        if nd >= 0 && !self.s.node_pass_.is_empty() && self.s.node_pass_[nd as usize] != 0 && self.algebraic_active(nd) {
            let und2 = nd as usize;
            let nb2 = m.node_face_ptr[und2];
            let ne2 = m.node_face_ptr[und2 + 1];
            let mut s_this = 0.0;
            let mut s_far = 0.0;
            let mut oc: i32 = -1;
            let mut uof = 0usize;
            for p2 in nb2..ne2 {
                let up2 = p2 as usize;
                let of = m.node_face_idx[up2];
                if of == f {
                    s_this = m.node_face_sign[up2];
                    continue;
                }
                uof = of as usize;
                s_far = m.node_face_sign[up2];
                oc = if m.face_cl[uof] >= 0 { m.face_cl[uof] } else { m.face_cr[uof] };
            }
            if oc >= 0 && s_this != 0.0 && s_far != 0.0 {
                let odir = if m.face_cl[uof] >= 0 { m.face_dir_l[uof] as f64 } else { m.face_dir_r[uof] as f64 };
                u_int = -s_this * s_far * odir * self.s.cell_u_[oc as usize];
            }
        }

        let (mut l, i1l, _) = self.face_side(f, cl, nd, zstar, m.face_dir_l[uf] as i32, u_int, false);
        let (mut r, i1r, _) = self.face_side(f, cr, nd, zstar, m.face_dir_r[uf] as i32, u_int, false);

        if mass_only {
            let mut mass = k::riemann_mass_flux(&l, &r);
            let gate_m = m.face_gate[uf];
            if gate_m != 0 && ((mass > 0.0 && (gate_m & 1) != 0) || (mass < 0.0 && (gate_m & 2) != 0)) {
                mass = 0.0;
            }
            let cvr_m = m.face_culvert[uf];
            if cvr_m >= 0 && mass > 0.0 && nd >= 0 {
                let ucv = cvr_m as usize;
                let gc = &m.geom[m.conduit_section[ucv] as usize];
                let head = self.st.node_head[nd as usize] - m.face_zb[uf];
                let mut dqdh = 0.0;
                let q_cap = culvert_inflow(mass, head, gc.y_full, gc.xs.a_full, m.conduit_slope[ucv], &m.conduit_culvert_curve[ucv], m.conduit_culvert_mitered[ucv] != 0, &mut dqdh);
                if q_cap < mass {
                    mass = q_cap;
                    self.s.inlet_control_[ucv] = 1;
                }
            }
            if nd >= 0 && nd < self.fc.n_nodes {
                if let Some(sm) = self.fc.node_source_max.as_ref() {
                    let cap = sm[nd as usize];
                    if cl < 0 && mass > cap {
                        mass = cap;
                    }
                    if cr < 0 && mass < -cap {
                        mass = -cap;
                    }
                }
            }
            self.s.f_mass_[uf] = mass;
            return;
        }

        let mut fl: FaceFlux = k::riemann_flux(&l, &r);
        let gate = m.face_gate[uf];
        if gate != 0 && ((fl.mass > 0.0 && (gate & 1) != 0) || (fl.mass < 0.0 && (gate & 2) != 0)) {
            if cl < 0 {
                l = r;
                l.u = -r.u;
                l.q = -r.q;
            } else {
                r = l;
                r.u = -l.u;
                r.q = -l.q;
            }
            fl = k::riemann_flux(&l, &r);
            fl.mass = 0.0;
        }
        let cvr = m.face_culvert[uf];
        if cvr >= 0 && fl.mass > 0.0 && nd >= 0 {
            let ucv = cvr as usize;
            let gc = &m.geom[m.conduit_section[ucv] as usize];
            let head = self.st.node_head[nd as usize] - m.face_zb[uf];
            let mut dqdh = 0.0;
            let q_cap = culvert_inflow(fl.mass, head, gc.y_full, gc.xs.a_full, m.conduit_slope[ucv], &m.conduit_culvert_curve[ucv], m.conduit_culvert_mitered[ucv] != 0, &mut dqdh);
            if q_cap < fl.mass {
                let up = if cl < 0 { &l } else { &r };
                let a_up = if up.a > k::K_DRY_AREA { up.a } else { k::K_DRY_AREA };
                fl.mass = q_cap;
                fl.mom = q_cap * q_cap / a_up + k::K_GRAVITY * up.i1;
                self.s.inlet_control_[cvr as usize] = 1;
            }
        }
        if nd >= 0 && nd < self.fc.n_nodes {
            if let Some(sm) = self.fc.node_source_max.as_ref() {
                let cap = sm[nd as usize];
                let accepted = if cl < 0 { cmin(fl.mass, cap) } else { cmax(fl.mass, -cap) };
                if accepted != fl.mass {
                    let up = if cl < 0 { &l } else { &r };
                    fl.mass = accepted;
                    fl.mom = accepted * accepted / cmax(up.a, k::K_DRY_AREA) + k::K_GRAVITY * up.i1;
                }
            }
        }
        self.s.f_mass_[uf] = fl.mass;
        self.s.f_mom_[uf] = fl.mom;
        self.s.f_state_l_[uf] = l;
        self.s.f_state_r_[uf] = r;
        if self.st.n_species > 0 {
            self.s.f_sstar_[uf] = fl.sstar;
            self.s.f_flux_[uf] = fl;
        }
        self.s.f_corr_l_[uf] = if cl >= 0 { k::K_GRAVITY * (i1l - l.i1) } else { 0.0 };
        self.s.f_corr_r_[uf] = if cr >= 0 { k::K_GRAVITY * (i1r - r.i1) } else { 0.0 };
    }

    /// algebraicAvailableWater(node, dt, pending, cell_scales)
    pub(crate) fn algebraic_available_water(&self, node: usize, dt: f64, pending: f64, cell_scales: &[f64]) -> f64 {
        let m = self.m;
        let mut volume = self.st.node_volume[node] + cmin(0.0, self.s.node_carry_[node]) + pending;
        // (forcing_ is set)
        volume += dt * (self.node_lateral(node) + self.s.node_qstruct_[node] + self.node_q_dummy(node));
        for p in m.node_face_ptr[node]..m.node_face_ptr[node + 1] {
            let f = m.node_face_idx[p as usize] as usize;
            let incoming = m.node_face_sign[p as usize] * self.s.f_mass_[f];
            if incoming <= 0.0 {
                continue;
            }
            let c = if m.face_cl[f] >= 0 { m.face_cl[f] } else { m.face_cr[f] };
            volume += dt * incoming * if c >= 0 { cell_scales[c as usize] } else { 1.0 };
        }
        cmax(0.0, volume)
    }

    /// limitPositivity(dt)
    pub(crate) fn limit_positivity(&mut self, dt: f64) {
        let m = self.m;
        let nc = m.n_cells() as usize;
        let nn = m.n_nodes() as usize;
        let mut out_cell = vec![0.0f64; nc];
        let mut out_node = vec![0.0f64; nn];

        for ai in 0..self.s.active_faces_.len() {
            let f = self.s.active_faces_[ai];
            let uf = f as usize;
            let nd = m.face_node[uf];
            if nd >= 0 && nd < self.fc.n_nodes {
                if let Some(sm) = self.fc.node_source_max.as_ref() {
                    let cap = sm[nd as usize];
                    let left = m.face_cl[uf] < 0;
                    let accepted = if left { cmin(self.s.f_mass_[uf], cap) } else { cmax(self.s.f_mass_[uf], -cap) };
                    if accepted != self.s.f_mass_[uf] {
                        let up = if left { self.s.f_state_l_[uf] } else { self.s.f_state_r_[uf] };
                        self.s.f_mass_[uf] = accepted;
                        self.s.f_mom_[uf] = accepted * accepted / cmax(up.a, k::K_DRY_AREA) + k::K_GRAVITY * up.i1;
                        if self.st.n_species > 0 {
                            self.s.f_flux_[uf].mass = self.s.f_mass_[uf];
                            self.s.f_flux_[uf].mom = self.s.f_mom_[uf];
                        }
                    }
                }
            }
            let fa = self.s.f_mass_[uf];
            if fa == 0.0 {
                continue;
            }
            if fa > 0.0 {
                if m.face_cl[uf] >= 0 {
                    out_cell[m.face_cl[uf] as usize] += fa;
                } else if m.face_node[uf] >= 0 {
                    out_node[m.face_node[uf] as usize] += fa;
                }
            } else if m.face_cr[uf] >= 0 {
                out_cell[m.face_cr[uf] as usize] -= fa;
            } else if m.face_node[uf] >= 0 {
                out_node[m.face_node[uf] as usize] -= fa;
            }
        }

        for uc in 0..nc {
            if out_cell[uc] <= 0.0 {
                out_cell[uc] = 1.0;
                continue;
            }
            out_cell[uc] = k::positivity_scale(self.st.cell_a[uc] * m.cell_dx[uc], out_cell[uc], dt);
        }
        if m.struct_is_dummy.iter().any(|&v| v != 0) {
            self.refresh_dummy_flows(dt, Arrived::CellScaled(&out_cell, dt), false);
        }
        for un in 0..nn {
            if out_node[un] <= 0.0 {
                out_node[un] = 1.0;
                continue;
            }
            if self.fc.fixed(un) {
                out_node[un] = 1.0;
                continue;
            }
            if self.algebraic_active(un as i32) {
                out_node[un] = if self.s.conservative_node_balance_ { k::positivity_scale(self.algebraic_available_water(un, dt, 0.0, &out_cell), out_node[un], dt) } else { 1.0 };
                continue;
            }
            out_node[un] = k::positivity_scale(self.st.node_volume[un], out_node[un], dt);
        }

        for ai in 0..self.s.active_faces_.len() {
            let f = self.s.active_faces_[ai];
            let uf = f as usize;
            let fa = self.s.f_mass_[uf];
            if fa == 0.0 {
                continue;
            }
            let s = if fa > 0.0 {
                if m.face_cl[uf] >= 0 {
                    out_cell[m.face_cl[uf] as usize]
                } else if m.face_node[uf] >= 0 {
                    out_node[m.face_node[uf] as usize]
                } else {
                    1.0
                }
            } else if m.face_cr[uf] >= 0 {
                out_cell[m.face_cr[uf] as usize]
            } else if m.face_node[uf] >= 0 {
                out_node[m.face_node[uf] as usize]
            } else {
                1.0
            };
            if s < 1.0 {
                self.s.f_mass_[uf] *= s;
                self.s.f_mom_[uf] *= s;
            }
        }
    }
}
