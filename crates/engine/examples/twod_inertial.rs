// Rust side of gate/diff/twod_inertial.cpp: the same inputs through
// twod::solver::inertial_kernels, printed the same way.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
#![allow(clippy::too_many_arguments)]
#[path = "common/twodgen.rs"]
mod twodgen;
use swmm6_rust::twod::mesh::quad_vfr::quad_vfr_precompute;
use swmm6_rust::twod::solver::inertial_kernels::*;
use twodgen::*;

fn args(o: &mut Out, a: &[f64]) {
    for &v in a {
        o.hx(v);
    }
    o.sep();
}

fn b2f(b: bool) -> f64 {
    if b {
        1.0
    } else {
        0.0
    }
}

fn main() {
    let mut o = Out::new();
    let mut g = Gen(0x2d0f1a7e5eed0003);
    let n = 40000;

    // ---- momentumBlendFraction
    let th = [1.0, 0.5, 0.9, 0.99, 0.999999, 1e-3, 0.0, -0.0, -0.5, 1.0f64.next_down(), 2.0, DMIN];
    let dt_ = [0.0, -0.0, 1e-12, 1e-9, 0.01, 1.0, 60.0, 3600.0, -1.0, DMIN];
    let rs = [0.0, -0.0, 1.0, 0.5, 1e-6, 1e6, -1.0, DMIN];
    grid(&[th.len(), dt_.len(), rs.len()], |i| {
        let (a, b, c) = (opq(th[i[0]]), opq(dt_[i[1]]), opq(rs[i[2]]));
        o.tag("mbf");
        args(&mut o, &[a, b, c]);
        o.hx(momentum_blend_fraction(a, b, c));
        o.nl();
    });
    for _ in 0..n {
        let a = g.draw(&th, 0.0, 1.0);
        let b = g.draw(&dt_, 0.0, 120.0);
        let c = g.draw(&rs, 0.0, 120.0);
        o.tag("mbf");
        args(&mut o, &[a, b, c]);
        o.hx(momentum_blend_fraction(opq(a), opq(b), opq(c)));
        o.nl();
    }

    // ---- qMagnitude
    let qs = [0.0, -0.0, 1.0, -1.0, 3.0, 4.0, 1e-160, 1e160, 1e-300, DMIN, 0.1, -2.5];
    grid(&[qs.len(), qs.len()], |i| {
        let (a, b) = (opq(qs[i[0]]), opq(qs[i[1]]));
        o.tag("qm");
        args(&mut o, &[a, b]);
        o.hx(q_magnitude(a, b));
        o.nl();
    });
    for _ in 0..n {
        let a = g.draw(&qs, -10.0, 10.0);
        let b = g.draw(&qs, -10.0, 10.0);
        o.tag("qm");
        args(&mut o, &[a, b]);
        o.hx(q_magnitude(opq(a), opq(b)));
        o.nl();
    }

    // ---- etaDepthScalar / volumeFromEtaScalar
    let ar = [1.0, 0.0, -0.0, 1e-30, 1e-30f64.next_up(), 1e-30f64.next_down(), 25.0, 1e4];
    let zz = [0.0, -0.0, 1.0, 1e-9, 100.0, -3.0, 0.5];
    let ep = [0.0, -0.0, 1e-6, 0.01, 0.1, 0.5, 1.0];
    let vv = [0.0, -0.0, -1.0, 1e-12, 1.0, 50.0, DMIN];
    for _ in 0..n {
        let area = g.draw(&ar, 0.0, 500.0);
        let cz = g.draw(&zz, -5.0, 50.0);
        let base = g.uni(-5.0, 50.0);
        let (z1, z2, z3);
        if g.below(3) == 0 {
            z1 = g.draw(&zz, -5.0, 5.0);
            z2 = g.draw(&zz, -5.0, 5.0);
            z3 = g.draw(&zz, -5.0, 5.0);
        } else {
            z1 = base + g.uni(0.0, 2.0);
            z2 = base + g.uni(0.0, 2.0);
            z3 = base + g.uni(0.0, 2.0);
        }
        let vfr = g.below(2) == 0;
        let eps = g.draw(&ep, 0.0, 0.3);
        let vol = g.draw(&vv, 0.0, 2.0 * area);
        let eta = if g.below(2) == 0 { g.uni(z1 - 0.5, z3 + 2.0) } else { g.draw(&zz, -5.0, 50.0) };
        let (e, d) = eta_depth_scalar(opq(area), opq(cz), opq(z1), opq(z2), opq(z3), vfr, opq(eps), opq(vol));
        o.tag("eds");
        args(&mut o, &[area, cz, z1, z2, z3, b2f(vfr), eps, vol]);
        o.hx(e);
        o.hx(d);
        o.nl();
        o.tag("vfe");
        args(&mut o, &[area, cz, z1, z2, z3, b2f(vfr), eps, eta]);
        o.hx(volume_from_eta_scalar(opq(area), opq(cz), opq(z1), opq(z2), opq(z3), vfr, opq(eps), opq(eta)));
        o.nl();
    }

    // ---- etaDepthQuadScalar / volumeFromEtaQuadScalar
    for _ in 0..n / 2 {
        let mut x = [0.0, 1.0, 1.0, 0.0];
        let mut y = [0.0, 0.0, 1.0, 1.0];
        let mut zq = [0.0f64; 4];
        let mut zs = [0.0f64; 6];
        let mut a1 = 0.0;
        let mut a2 = 0.0;
        let base = g.uni(-5.0, 50.0);
        let sc = g.uni(0.1, 30.0);
        for v in 0..4 {
            x[v] *= sc;
            y[v] *= sc;
            zq[v] = base + g.uni(0.0, 2.0);
        }
        if g.below(5) == 0 {
            let src = g.below(4) as usize; // C++17: the right operand of = is sequenced first
            let dst = g.below(4) as usize;
            zq[dst] = zq[src];
        }
        quad_vfr_precompute(&x, &y, &zq, &mut zs, &mut a1, &mut a2);
        let area = if g.below(8) == 0 { g.pick(&ar) } else { a1 + a2 };
        let cz = (zq[0] + zq[1] + zq[2] + zq[3]) / 4.0;
        let vfr = g.below(2) == 0;
        let eps = g.draw(&ep, 0.0, 0.3);
        let vol = g.draw(&vv, 0.0, 2.0 * area);
        let eta = if g.below(2) == 0 { g.uni(base - 0.5, base + 3.0) } else { g.draw(&zz, -5.0, 50.0) };
        for v in 0..6 {
            zs[v] = opq(zs[v]);
        }
        let (e, d) = eta_depth_quad_scalar(opq(area), opq(cz), &zs, opq(a1), opq(a2), vfr, opq(eps), opq(vol));
        o.tag("edq");
        args(&mut o, &[area, cz, zs[0], zs[1], zs[2], zs[3], zs[4], zs[5], a1, a2, b2f(vfr), eps, vol]);
        o.hx(e);
        o.hx(d);
        o.nl();
        o.tag("vfq");
        args(&mut o, &[area, cz, zs[0], zs[1], zs[2], zs[3], zs[4], zs[5], a1, a2, b2f(vfr), eps, eta]);
        o.hx(volume_from_eta_quad_scalar(opq(area), opq(cz), &zs, opq(a1), opq(a2), vfr, opq(eps), opq(eta)));
        o.nl();
    }

    // ---- faceFlowDepth, faceDepthFromEta, faceFlowDepthVfr
    let et = [0.0, -0.0, 1.0, 1.0 + 1e-9, 1.0 + 5e-10, 2.0, 0.5, -1.0, 1e-12, 3.0];
    grid(&[et.len(), et.len(), et.len()], |i| {
        let (a, b, c) = (opq(et[i[0]]), opq(et[i[1]]), opq(et[i[2]]));
        o.tag("ffd");
        args(&mut o, &[a, b, c]);
        o.hx(face_flow_depth(a, b, c));
        o.nl();
        o.tag("fde");
        args(&mut o, &[a, b, c]);
        o.hx(face_depth_from_eta(a, b, c));
        o.nl();
    });
    for _ in 0..n {
        let zlo = g.draw(&et, -2.0, 10.0);
        let m = g.below(4);
        let zhi = if m == 0 {
            zlo + g.uni(0.0, 2e-9)
        } else if m == 1 {
            g.draw(&et, -2.0, 10.0)
        } else {
            zlo + g.uni(0.0, 3.0)
        };
        let sel = [zlo, zhi, zlo.next_up(), zhi.next_up(), zhi.next_down()];
        let eta = if g.below(3) == 0 { sel[g.below(5) as usize] } else { g.uni(zlo - 1.0, zhi + 1.0) };
        let eta2 = g.draw(&et, zlo - 1.0, zhi + 1.0);
        o.tag("fde");
        args(&mut o, &[eta, zlo, zhi]);
        o.hx(face_depth_from_eta(opq(eta), opq(zlo), opq(zhi)));
        o.nl();
        o.tag("ffv");
        args(&mut o, &[eta, eta2, zlo, zhi]);
        o.hx(face_flow_depth_vfr(opq(eta), opq(eta2), opq(zlo), opq(zhi)));
        o.nl();
        o.tag("ffd");
        args(&mut o, &[eta, eta2, zlo]);
        o.hx(face_flow_depth(opq(eta), opq(eta2), opq(zlo)));
        o.nl();
    }

    // ---- inertialFaceUpdate
    let hf_ = [1e-3, 1.0, 0.0, -0.0, 1e-12, 1e-6, 5.0, DMIN, 1e-300, -1.0];
    let sl = [0.0, -0.0, 1e-3, -1e-3, 0.05, -0.05, 1e-12, 1.0];
    let qq = [0.0, -0.0, 0.5, -0.5, 2.0, 1e-9, -10.0];
    let nn = [0.0, 0.0004, 0.0025, 0.01, 1.0];
    for _ in 0..2 * n {
        let q = g.draw(&qq, -5.0, 5.0);
        let qhat = g.draw(&qq, -5.0, 5.0);
        let hf = g.draw(&hf_, 0.0, 3.0);
        let dt = g.draw(&dt_, 0.0, 10.0);
        let slope = g.draw(&sl, -0.1, 0.1);
        let n2 = g.draw(&nn, 0.0, 0.01);
        let qmag = g.draw(&qq, 0.0, 5.0);
        let adv = if g.below(2) == 0 { 0.0 } else { g.draw(&qq, -1.0, 1.0) };
        o.tag("ifu");
        args(&mut o, &[q, qhat, hf, dt, slope, n2, qmag, adv]);
        o.hx(inertial_face_update(opq(q), opq(qhat), opq(hf), opq(dt), opq(slope), opq(n2), opq(qmag), opq(adv)));
        o.nl();
    }
    for _ in 0..n / 4 {
        let q = g.draw(&qq, -5.0, 5.0);
        let qhat = g.draw(&qq, -5.0, 5.0);
        let hf = g.draw(&hf_, 0.0, 3.0);
        let dt = g.draw(&dt_, 0.0, 10.0);
        let slope = g.draw(&sl, -0.1, 0.1);
        let n2 = g.draw(&nn, 0.0, 0.01);
        let qmag = g.draw(&qq, 0.0, 5.0);
        o.tag("ifu0");
        args(&mut o, &[q, qhat, hf, dt, slope, n2, qmag]);
        o.hx(inertial_face_update(opq(q), opq(qhat), opq(hf), opq(dt), opq(slope), opq(n2), opq(qmag), 0.0));
        o.nl();
    }

    // ---- inertialAdvection
    let un = [0.0, -0.0, 1.0, -1.0, 0.3, -0.3, 1e-12, -1e-12];
    for _ in 0..n {
        let qf = g.draw(&qq, -5.0, 5.0);
        let un_l = g.draw(&un, -3.0, 3.0);
        let h_l = g.draw(&hf_, 0.0, 3.0);
        let un_r = g.draw(&un, -3.0, 3.0);
        let h_r = g.draw(&hf_, 0.0, 3.0);
        let idx = g.draw(&ar, 0.01, 2.0);
        o.tag("adv");
        args(&mut o, &[qf, un_l, h_l, un_r, h_r, idx]);
        o.hx(inertial_advection(opq(qf), opq(un_l), opq(h_l), opq(un_r), opq(h_r), opq(idx)));
        o.nl();
    }

    // ---- froudeCap
    let fr_ = [1.0, 0.0, -0.0, 0.9, 2.0, -1.0, 1e-9];
    grid(&[qq.len(), hf_.len(), fr_.len()], |i| {
        let (a, b, c) = (opq(qq[i[0]]), opq(hf_[i[1]]), opq(fr_[i[2]]));
        o.tag("fc");
        args(&mut o, &[a, b, c]);
        o.hx(froude_cap(a, b, c));
        o.nl();
    });
    for _ in 0..n {
        let q = g.draw(&qq, -20.0, 20.0);
        let hf = g.draw(&hf_, 0.0, 3.0);
        let fr = g.draw(&fr_, 0.0, 3.0);
        o.tag("fc");
        args(&mut o, &[q, hf, fr]);
        o.hx(froude_cap(opq(q), opq(hf), opq(fr)));
        o.nl();
    }

    // ---- cellCflDt
    let al = [0.5, 0.9, 1.0, 0.0, -0.5];
    let lc = [1.0, 10.0, 0.0, 1e-3, 1e4];
    let hh = [0.0, -0.0, 1e-30, 1.0e-25, 1.0198e-25, 1.02e-25, 1.0, 1e-6, -1.0];
    let sp = [0.0, -0.0, 1e-12, 1e-12f64.next_up(), 1.0, 5.0, -1.0];
    grid(&[al.len(), lc.len(), hh.len(), sp.len()], |i| {
        let (a, b, c, d) = (opq(al[i[0]]), opq(lc[i[1]]), opq(hh[i[2]]), opq(sp[i[3]]));
        o.tag("cfl");
        args(&mut o, &[a, b, c, d]);
        o.hx(cell_cfl_dt(a, b, c, d));
        o.nl();
    });
    for _ in 0..n {
        let a = g.draw(&al, 0.1, 1.0);
        let b = g.draw(&lc, 0.1, 100.0);
        let c = g.draw(&hh, 0.0, 3.0);
        let d = g.draw(&sp, 0.0, 5.0);
        o.tag("cfl");
        args(&mut o, &[a, b, c, d]);
        o.hx(cell_cfl_dt(opq(a), opq(b), opq(c), opq(d)));
        o.nl();
    }

    // ---- positivityScale
    let pv = [0.0, -0.0, -1.0, 1.0, 10.0, 1e-12, DMIN];
    let po = [0.0, -0.0, -1.0, 1.0, 0.1, 100.0, 1e-12];
    let pd = [0.0, -0.0, 1.0, 0.1, 10.0, -1.0];
    let pb = [1.0, 0.0, 0.5, 0.95];
    grid(&[pv.len(), po.len(), pd.len(), pb.len()], |i| {
        let (a, b, c, d) = (opq(pv[i[0]]), opq(po[i[1]]), opq(pd[i[2]]), opq(pb[i[3]]));
        o.tag("pos");
        args(&mut o, &[a, b, c, d]);
        o.hx(positivity_scale(a, b, c, d));
        o.nl();
    });
    for _ in 0..n {
        let a = g.draw(&pv, -1.0, 20.0);
        let b = g.draw(&po, 0.0, 20.0);
        let c = g.draw(&pd, 0.0, 5.0);
        let d = g.draw(&pb, 0.0, 1.0);
        o.tag("pos");
        args(&mut o, &[a, b, c, d]);
        o.hx(positivity_scale(opq(a), opq(b), opq(c), opq(d)));
        o.nl();
    }

    o.flush();
    eprintln!("twod_inertial: {} lines", o.lines);
}
