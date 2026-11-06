// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/fvkernels/SpeciesTransportKernels.hpp and
// SpeciesTransportKernels.cpp, and speciesFlux from src/engine/hydraulics/fv/FvKernels.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The shared 1D FV species transport kernels: the advection reconstruction (MUSCL, UPWIND,
// QUICKEST-ULTIMATE), the Zalesak FCT flux limiter and the implicit per-chain dispersion solve.
//
// Differences in form, never in result:
//   * SpeciesKernelView holds Rust borrows instead of pointers; a caller builds one for the
//     reconstruction (which needs every scratch vector) and a smaller one for the dispersion
//     solve, which reads only the mesh, the state and the dispersion fields.
//   * dispersionSolve's `static thread_local` Thomas work arrays (aa, bb, ccv, rr, xx) are
//     per-call locals here: each chain assigns all five before any element is read.
//   * std::max({a, b, c}) / std::min({a, b, c}) (max_element / min_element: the first of the
//     largest) is cmax(cmax(a, b), c) / cmin(cmin(a, b), c), the same comparisons in the same
//     order.

use crate::core::cmath::{cclamp, cmax, cmin};
use crate::core::options::fv::{Limiter, ScalarScheme};
use crate::hydraulics::fv::explicit_solver::{FaceFlux, FaceState};
use crate::hydraulics::fv::kernels::K_DRY_AREA;
use crate::hydraulics::fv::{NetworkMeshData, NetworkStateData};

/// fv::kernels::speciesFlux (FvKernels.hpp): the scalar flux of one face.
pub fn species_flux(l: &FaceState, r: &FaceState, f: &FaceFlux, phi_l: f64, phi_r: f64, hllc: bool) -> f64 {
    if hllc {
        return f.mass * (if f.sstar >= 0.0 { phi_l } else { phi_r });
    }
    if f.sl >= 0.0 {
        return l.q * phi_l;
    }
    if f.sr <= 0.0 {
        return r.q * phi_r;
    }
    let dsr = f.sr - f.sl;
    (f.sr * l.q * phi_l - f.sl * r.q * phi_r + f.sl * f.sr * (r.a * phi_r - l.a * phi_l)) / dsr
}

/// limitSlope: van Leer, superbee or minmod.
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

/// SpeciesKernelView: everything the reconstruction and the FCT limiter read and write.
pub struct SpeciesKernelView<'a> {
    pub mesh: &'a NetworkMeshData,
    pub state: &'a mut NetworkStateData,
    pub scalar_scheme: ScalarScheme,
    pub limiter: Limiter,
    pub hllc: bool,
    pub f_mass: &'a [f64],
    pub f_sstar: &'a [f64],
    pub f_state_l: &'a [FaceState],
    pub f_state_r: &'a [FaceState],
    pub f_flux: &'a [FaceFlux],
    pub cell_u: &'a [f64],
    pub cell_active: &'a [u8],
    pub active_faces: &'a [i32],
    pub f_phi_l: &'a mut Vec<f64>,
    pub f_phi_r: &'a mut Vec<f64>,
    pub f_phi_flux: &'a mut Vec<f64>,
    pub cell_slope: &'a mut Vec<f64>,
    pub lo_flux: &'a mut Vec<f64>,
    pub anti_flux: &'a mut Vec<f64>,
    pub td: &'a mut Vec<f64>,
    pub anew: &'a mut Vec<f64>,
    pub rplus: &'a mut Vec<f64>,
    pub rminus: &'a mut Vec<f64>,
}

/// The part of the view dispersionSolve reads.
pub struct DispersionView<'a> {
    pub mesh: &'a NetworkMeshData,
    pub state: &'a mut NetworkStateData,
    pub dispersion: f64,
    pub cell_dispersion: Option<&'a [f64]>,
}

fn assign(v: &mut Vec<f64>, n: usize, x: f64) {
    v.clear();
    v.resize(n, x);
}

/// reconstructScalars: bounded face values and limited fluxes for every species row.
pub fn reconstruct_scalars(v: &mut SpeciesKernelView, dt: f64) {
    let ns = v.state.n_species;
    if ns <= 0 {
        return;
    }
    let mesh = v.mesh;
    let nc = mesh.n_cells();
    let nf = mesh.n_faces();
    for s in 0..ns {
        let sbase = s as usize * nc as usize;
        let fbase = s as usize * nf as usize;
        if v.scalar_scheme == ScalarScheme::Upwind {
            v.cell_slope.iter_mut().for_each(|x| *x = 0.0);
        } else {
            for ch in 0..mesh.n_chains() {
                let b = mesh.chain_ptr[ch as usize];
                let e = mesh.chain_ptr[ch as usize + 1];
                for i in b..e {
                    let c = mesh.chain_cells[i as usize];
                    let uc = c as usize;
                    if i == b || i == e - 1 {
                        v.cell_slope[uc] = 0.0;
                        continue;
                    }
                    let um = mesh.chain_cells[(i - 1) as usize] as usize;
                    let up = mesh.chain_cells[(i + 1) as usize] as usize;
                    let dm = 0.5 * (mesh.cell_dx[um] + mesh.cell_dx[uc]);
                    let dp = 0.5 * (mesh.cell_dx[uc] + mesh.cell_dx[up]);
                    let gm = (v.state.cell_phi[sbase + uc] - v.state.cell_phi[sbase + um]) / dm;
                    let gp = (v.state.cell_phi[sbase + up] - v.state.cell_phi[sbase + uc]) / dp;
                    v.cell_slope[uc] = limit_slope(gm, gp, v.limiter) * mesh.chain_dir[i as usize] as f64;
                }
            }
        }
        for &f in v.active_faces.iter() {
            let uf = f as usize;
            let cl = mesh.face_cl[uf];
            let cr = mesh.face_cr[uf];
            let el = if cl >= 0 { cl } else { cr };
            let er = if cr >= 0 { cr } else { cl };
            let ul = el as usize;
            let ur = er as usize;
            let mut phil = v.state.cell_phi[sbase + ul];
            let mut phir = v.state.cell_phi[sbase + ur];
            if cl >= 0 {
                let sign = mesh.face_dir_l[uf] as f64;
                phil += sign * v.cell_slope[ul] * 0.5 * mesh.cell_dx[ul];
            }
            if cr >= 0 {
                let sign = -(mesh.face_dir_r[uf] as f64);
                phir += sign * v.cell_slope[ur] * 0.5 * mesh.cell_dx[ur];
            }
            if v.scalar_scheme == ScalarScheme::QuickestUltimate {
                let fwd = v.f_sstar[uf] >= 0.0;
                let cu = if fwd { cl } else { cr };
                let cd = if fwd { cr } else { cl };
                if cu >= 0 && cd >= 0 {
                    let uu = cu as usize;
                    let chn = mesh.cell_chain[uu];
                    let pos = mesh.cell_chain_pos[uu];
                    let cb = mesh.chain_ptr[chn as usize];
                    let ce = mesh.chain_ptr[chn as usize + 1];
                    let dpos = mesh.cell_chain_pos[cd as usize];
                    let upos = pos + (pos - dpos);
                    if upos >= 0 && upos < ce - cb {
                        let cuu = mesh.chain_cells[(cb + upos) as usize];
                        let uuu = cuu as usize;
                        let pu = v.state.cell_phi[sbase + uu];
                        let pd = v.state.cell_phi[sbase + cd as usize];
                        let puu = v.state.cell_phi[sbase + uuu];
                        let dx = mesh.cell_dx[uu];
                        let cr_no = cmin(1.0, v.cell_u[uu].abs() * dt / cmax(dx, 1.0e-12));
                        let q = pd - 2.0 * pu + puu;
                        let mut pf = 0.5 * (pd + pu) - 0.5 * cr_no * (pd - pu) - (1.0 - cr_no * cr_no) / 6.0 * q;
                        let den = pd - puu;
                        if den.abs() > 1.0e-30 {
                            let un_ = (pu - puu) / den;
                            if un_ <= 0.0 || un_ >= 1.0 {
                                pf = pu;
                            } else {
                                let mut fn_ = (pf - puu) / den;
                                let hi = cmin(1.0, un_ / cmax(cr_no, 1.0e-12));
                                fn_ = cclamp(fn_, un_, hi);
                                pf = puu + fn_ * den;
                            }
                        } else {
                            pf = pu;
                        }
                        if fwd {
                            phil = pf;
                        } else {
                            phir = pf;
                        }
                    }
                }
            }
            if v.scalar_scheme != ScalarScheme::Upwind {
                let a = v.state.cell_phi[sbase + ul];
                let b = v.state.cell_phi[sbase + ur];
                let lo = cmin(a, b);
                let hi = cmax(a, b);
                phil = cclamp(phil, lo, hi);
                phir = cclamp(phir, lo, hi);
            }
            v.f_phi_l[fbase + uf] = phil;
            v.f_phi_r[fbase + uf] = phir;
        }
        limit_species_fluxes(v, s, dt);
    }
}

/// limitSpeciesFluxes: the face fluxes of one species row, Zalesak FCT.
pub fn limit_species_fluxes(v: &mut SpeciesKernelView, species: i32, dt: f64) {
    let mesh = v.mesh;
    let nc = mesh.n_cells();
    let nf = mesh.n_faces();
    let cbase = species as usize * nc as usize;
    let fbase = species as usize * nf as usize;
    assign(v.lo_flux, nf as usize, 0.0);
    assign(v.anti_flux, nf as usize, 0.0);
    for &f in v.active_faces.iter() {
        let uf = f as usize;
        let cl = mesh.face_cl[uf];
        let cr = mesh.face_cr[uf];
        let ul = (if cl >= 0 { cl } else { cr }) as usize;
        let ur = (if cr >= 0 { cr } else { cl }) as usize;
        let fm = v.f_mass[uf];
        let lo = fm * (if fm >= 0.0 { v.state.cell_phi[cbase + ul] } else { v.state.cell_phi[cbase + ur] });
        // adjustedFlux(f): the face flux record with the positivity-scaled mass flux
        let mut adj = v.f_flux[uf];
        adj.mass = v.f_mass[uf];
        let hi = species_flux(&v.f_state_l[uf], &v.f_state_r[uf], &adj, v.f_phi_l[fbase + uf], v.f_phi_r[fbase + uf], v.hllc);
        v.lo_flux[uf] = lo;
        v.anti_flux[uf] = hi - lo;
    }
    assign(v.td, nc as usize, 0.0);
    assign(v.anew, nc as usize, 0.0);
    for c in 0..nc {
        let uc = c as usize;
        if v.cell_active[uc] == 0 {
            v.td[uc] = v.state.cell_phi[cbase + uc];
            v.anew[uc] = v.state.cell_a[uc];
            continue;
        }
        let faces = [mesh.cell_face0[uc], mesh.cell_face1[uc]];
        let sides = [mesh.cell_side0[uc], mesh.cell_side1[uc]];
        let mut d_a = 0.0;
        let mut dm = 0.0;
        for e in 0..2 {
            let uf = faces[e] as usize;
            let sg = if sides[e] == 0 { -1.0 } else { 1.0 };
            d_a += sg * v.f_mass[uf];
            dm += sg * v.lo_flux[uf];
        }
        let inv_dx = 1.0 / mesh.cell_dx[uc];
        let a_new = cmax(0.0, v.state.cell_a[uc] + dt * d_a * inv_dx);
        v.anew[uc] = a_new;
        let m = v.state.cell_a[uc] * v.state.cell_phi[cbase + uc] + dt * dm * inv_dx;
        v.td[uc] = if a_new > 0.0 { m / a_new } else { v.state.cell_phi[cbase + uc] };
    }
    assign(v.rplus, nc as usize, 1.0);
    assign(v.rminus, nc as usize, 1.0);
    for c in 0..nc {
        let uc = c as usize;
        if v.cell_active[uc] == 0 {
            continue;
        }
        let faces = [mesh.cell_face0[uc], mesh.cell_face1[uc]];
        let sides = [mesh.cell_side0[uc], mesh.cell_side1[uc]];
        let mut pmax = cmax(v.state.cell_phi[cbase + uc], v.td[uc]);
        let mut pmin = cmin(v.state.cell_phi[cbase + uc], v.td[uc]);
        let mut pplus = 0.0;
        let mut pminus = 0.0;
        for e in 0..2 {
            let uf = faces[e] as usize;
            let cl = mesh.face_cl[uf];
            let cr = mesh.face_cr[uf];
            let nb = if cl == c { cr } else { cl };
            if nb >= 0 {
                let un = nb as usize;
                pmax = cmax(cmax(pmax, v.state.cell_phi[cbase + un]), v.td[un]);
                pmin = cmin(cmin(pmin, v.state.cell_phi[cbase + un]), v.td[un]);
            }
            let sg = if sides[e] == 0 { -1.0 } else { 1.0 };
            let a = sg * v.anti_flux[uf];
            if a > 0.0 {
                pplus += a;
            } else {
                pminus -= a;
            }
        }
        let cap = v.anew[uc] * mesh.cell_dx[uc] / dt;
        if pplus > 0.0 {
            v.rplus[uc] = cmin(1.0, (pmax - v.td[uc]) * cap / pplus);
        }
        if pminus > 0.0 {
            v.rminus[uc] = cmin(1.0, (v.td[uc] - pmin) * cap / pminus);
        }
        v.rplus[uc] = cmax(0.0, v.rplus[uc]);
        v.rminus[uc] = cmax(0.0, v.rminus[uc]);
    }
    for &f in v.active_faces.iter() {
        let uf = f as usize;
        let a = v.anti_flux[uf];
        let mut coef = 1.0;
        if a != 0.0 {
            let cl = mesh.face_cl[uf];
            let cr = mesh.face_cr[uf];
            let from_r = if cr >= 0 {
                if a > 0.0 {
                    v.rplus[cr as usize]
                } else {
                    v.rminus[cr as usize]
                }
            } else {
                1.0
            };
            let from_l = if cl >= 0 {
                if a > 0.0 {
                    v.rminus[cl as usize]
                } else {
                    v.rplus[cl as usize]
                }
            } else {
                1.0
            };
            coef = cmin(from_r, from_l);
        }
        v.f_phi_flux[fbase + uf] = v.lo_flux[uf] + coef * a;
    }
}

/// dispersionSolve: implicit per-chain longitudinal dispersion (Thomas solve per species).
pub fn dispersion_solve(v: &mut DispersionView, dt: f64) {
    let ns = v.state.n_species;
    if ns <= 0 {
        return;
    }
    if v.cell_dispersion.is_none() && v.dispersion <= 0.0 {
        return;
    }
    let mesh = v.mesh;
    let (cdisp, disp) = (v.cell_dispersion, v.dispersion);
    let cell_d = |c: i32| -> f64 {
        match cdisp {
            Some(cd) => cd[c as usize],
            None => disp,
        }
    };
    let nc = mesh.n_cells();
    let mut aa: Vec<f64> = Vec::new();
    let mut bb: Vec<f64> = Vec::new();
    let mut ccv: Vec<f64> = Vec::new();
    let mut rr: Vec<f64> = Vec::new();
    let mut xx: Vec<f64> = Vec::new();
    for ch in 0..mesh.n_chains() {
        let b = mesh.chain_ptr[ch as usize];
        let e = mesh.chain_ptr[ch as usize + 1];
        let m = e - b;
        if m < 2 {
            continue;
        }
        let um = m as usize;
        assign(&mut aa, um, 0.0);
        assign(&mut bb, um, 0.0);
        assign(&mut ccv, um, 0.0);
        assign(&mut rr, um, 0.0);
        assign(&mut xx, um, 0.0);
        for s in 0..ns {
            let base = s as usize * nc as usize;
            for i in 0..m {
                let c = mesh.chain_cells[(b + i) as usize];
                let uc = c as usize;
                let dx = mesh.cell_dx[uc];
                let a = cmax(v.state.cell_a[uc], K_DRY_AREA);
                let dc = cell_d(c);
                let mut lo = 0.0;
                let mut hi = 0.0;
                if i > 0 {
                    let cm = mesh.chain_cells[(b + i - 1) as usize];
                    let am = cmax(v.state.cell_a[cm as usize], K_DRY_AREA);
                    lo = 0.5 * (dc + cell_d(cm)) * 0.5 * (a + am) * dt / (dx * dx * a);
                }
                if i < m - 1 {
                    let cp = mesh.chain_cells[(b + i + 1) as usize];
                    let ap = cmax(v.state.cell_a[cp as usize], K_DRY_AREA);
                    hi = 0.5 * (dc + cell_d(cp)) * 0.5 * (a + ap) * dt / (dx * dx * a);
                }
                aa[i as usize] = -lo;
                ccv[i as usize] = -hi;
                bb[i as usize] = 1.0 + lo + hi;
                rr[i as usize] = v.state.cell_phi[base + uc];
            }
            for i in 1..um {
                let w = aa[i] / bb[i - 1];
                bb[i] -= w * ccv[i - 1];
                rr[i] -= w * rr[i - 1];
            }
            let ulast = um - 1;
            xx[ulast] = rr[ulast] / bb[ulast];
            let mut i = m - 2;
            while i >= 0 {
                let ui = i as usize;
                xx[ui] = (rr[ui] - ccv[ui] * xx[ui + 1]) / bb[ui];
                i -= 1;
            }
            for i in 0..m {
                let c = mesh.chain_cells[(b + i) as usize];
                v.state.cell_phi[base + c as usize] = xx[i as usize];
            }
        }
    }
}
