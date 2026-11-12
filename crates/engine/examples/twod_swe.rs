// Rust side of gate/diff/twod_swe.cpp: the same inputs through twod::solver::swe_kernels,
// printed the same way.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
#![allow(clippy::too_many_arguments)]
#[path = "common/twodgen.rs"]
mod twodgen;
use swmm6_rust::twod::solver::swe_kernels::*;
use twodgen::*;

fn args(o: &mut Out, a: &[f64]) {
    for &v in a {
        o.hx(v);
    }
    o.sep();
}

const HS: [f64; 9] = [0.0, -0.0, 1e-12, 1e-6, 0.5, 1.0, 2.0, -1.0, DMIN];
const US: [f64; 9] = [0.0, -0.0, 1.0, -1.0, 3.0, -3.0, 10.0, -10.0, 1e-12];
const NXS: [f64; 8] = [1.0, 0.0, -1.0, -0.0, 0.6, -0.8, 0.7071067811865476, -0.7071067811865476];
const NYS: [f64; 8] = [0.0, 1.0, -0.0, -1.0, 0.8, 0.6, 0.7071067811865476, 0.7071067811865476];

fn side(g: &mut Gen) -> FaceSide {
    let h = g.draw(&HS, 0.0, 3.0);
    let un = g.draw(&US, -8.0, 8.0);
    let ut = g.draw(&US, -8.0, 8.0);
    FaceSide { h, un, ut }
}

fn normal(g: &mut Gen) -> (f64, f64) {
    if g.below(2) == 0 {
        let k = g.below(NXS.len() as u64) as usize;
        (NXS[k], NYS[k])
    } else {
        let nx = g.uni(-1.0, 1.0);
        let mut ny = (1.0 - nx * nx).sqrt();
        if g.below(2) == 0 {
            ny = -ny;
        }
        (nx, ny)
    }
}

fn face_flux_line(o: &mut Out, tag: &str, inp: &[f64], f: &FaceFlux, a: f64, b: f64, c: f64, d: f64, r: bool) {
    o.tag(tag);
    args(o, inp);
    o.hx(f.mass);
    o.hx(f.mx);
    o.hx(f.my);
    o.hx(f.sstar);
    o.hx(a);
    o.hx(b);
    o.hx(c);
    o.hx(d);
    o.hi(r as i64);
    o.nl();
}

fn neg7() -> FaceFlux {
    FaceFlux { mass: -7.0, mx: -7.0, my: -7.0, sstar: -7.0 }
}

fn main() {
    let mut o = Out::new();
    let mut g = Gen(0x2d0f1a7e5eed0004);
    let n = 40000;

    // ---- physicalFlux, waveSpeeds, hllcFlux
    let hg = [0.0, 1e-12, 0.5, 2.0, -1.0];
    let ug = [0.0, -0.0, 1.0, -1.0, 10.0, -10.0];
    let tg = [0.0, 2.0];
    grid(&[hg.len(), ug.len(), tg.len(), hg.len(), ug.len(), tg.len()], |i| {
        let l = FaceSide { h: opq(hg[i[0]]), un: opq(ug[i[1]]), ut: opq(tg[i[2]]) };
        let r = FaceSide { h: opq(hg[i[3]]), un: opq(ug[i[4]]), ut: opq(tg[i[5]]) };
        let (sl, sr) = wave_speeds(&l, &r);
        let (fh, fn_, ft, ss) = hllc_flux(&l, &r);
        o.tag("hllc");
        args(&mut o, &[l.h, l.un, l.ut, r.h, r.un, r.ut]);
        o.hx(sl);
        o.hx(sr);
        o.hx(fh);
        o.hx(fn_);
        o.hx(ft);
        o.hx(ss);
        o.nl();
    });
    for _ in 0..2 * n {
        let mut l = side(&mut g);
        let mut r = side(&mut g);
        if g.below(8) == 0 {
            r.h = l.h;
        }
        l = FaceSide { h: opq(l.h), un: opq(l.un), ut: opq(l.ut) };
        r = FaceSide { h: opq(r.h), un: opq(r.un), ut: opq(r.ut) };
        let (pf, pn, pt) = physical_flux(&l);
        let (sl, sr) = wave_speeds(&l, &r);
        let (fh, fn_, ft, ss) = hllc_flux(&l, &r);
        o.tag("hllc");
        args(&mut o, &[l.h, l.un, l.ut, r.h, r.un, r.ut]);
        o.hx(sl);
        o.hx(sr);
        o.hx(fh);
        o.hx(fn_);
        o.hx(ft);
        o.hx(ss);
        o.hx(pf);
        o.hx(pn);
        o.hx(pt);
        o.nl();
    }

    // ---- faceFlux
    let eta_ = [0.0, 1.0, 1.5, 2.0, -0.0, 1e-12];
    let qx_ = [0.0, -0.0, 0.5, -0.5, 3.0, -3.0];
    let hd_ = [1e-6, 0.0, 1e-3, 1e-12];
    grid(&[eta_.len(), HS.len(), eta_.len(), HS.len(), 4], |i| {
        let (eta_l, h_l, eta_r, h_r) = (opq(eta_[i[0]]), opq(HS[i[1]]), opq(eta_[i[2]]), opq(HS[i[3]]));
        let (nx, ny, hd) = (opq(NXS[i[4]]), opq(NYS[i[4]]), opq(1e-6));
        let mut f = neg7();
        let (mut a, mut b, mut c, mut d) = (-7.0, -7.0, -7.0, -7.0);
        let r = face_flux(
            eta_l,
            h_l,
            opq(0.0),
            opq(0.0),
            eta_r,
            h_r,
            opq(0.0),
            opq(0.0),
            nx,
            ny,
            hd,
            &mut f,
            &mut a,
            &mut b,
            &mut c,
            &mut d,
        );
        face_flux_line(&mut o, "ff", &[eta_l, h_l, 0.0, 0.0, eta_r, h_r, 0.0, 0.0, nx, ny, hd], &f, a, b, c, d, r);
    });
    for _ in 0..2 * n {
        let z_l = g.draw(&eta_, -2.0, 5.0);
        let z_r = if g.below(4) == 0 { z_l } else { g.draw(&eta_, -2.0, 5.0) };
        let h_l = g.draw(&HS, 0.0, 3.0);
        let h_r = g.draw(&HS, 0.0, 3.0);
        let eta_l = if g.below(16) == 0 { g.draw(&eta_, -2.0, 5.0) } else { canon(z_l + h_l) };
        let eta_r = if g.below(16) == 0 { g.draw(&eta_, -2.0, 5.0) } else { canon(z_r + h_r) };
        let qx_l = g.draw(&qx_, -5.0, 5.0);
        let qy_l = g.draw(&qx_, -5.0, 5.0);
        let qx_r = g.draw(&qx_, -5.0, 5.0);
        let qy_r = g.draw(&qx_, -5.0, 5.0);
        let (nx, ny) = normal(&mut g);
        let hd = g.draw(&hd_, 0.0, 1e-2);
        let mut f = neg7();
        let (mut a, mut b, mut c, mut d) = (-7.0, -7.0, -7.0, -7.0);
        let r = face_flux(
            opq(eta_l),
            opq(h_l),
            opq(qx_l),
            opq(qy_l),
            opq(eta_r),
            opq(h_r),
            opq(qx_r),
            opq(qy_r),
            opq(nx),
            opq(ny),
            opq(hd),
            &mut f,
            &mut a,
            &mut b,
            &mut c,
            &mut d,
        );
        face_flux_line(&mut o, "ff", &[eta_l, h_l, qx_l, qy_l, eta_r, h_r, qx_r, qy_r, nx, ny, hd], &f, a, b, c, d, r);
    }

    // ---- faceFluxReconBed and faceFluxRecon
    for _ in 0..2 * n {
        let z_l = g.draw(&eta_, -2.0, 5.0);
        let z_r = if g.below(4) == 0 { z_l } else { g.draw(&eta_, -2.0, 5.0) };
        let h_lc = g.draw(&HS, 0.0, 3.0);
        let h_rc = g.draw(&HS, 0.0, 3.0);
        let z_lf = if g.below(3) == 0 { z_l } else { z_l + g.uni(-0.2, 0.2) };
        let z_rf = if g.below(3) == 0 { z_r } else { z_r + g.uni(-0.2, 0.2) };
        let eta_lf = if g.below(16) == 0 { g.draw(&eta_, -2.0, 5.0) } else { canon(z_l + h_lc + g.uni(-0.1, 0.1)) };
        let eta_rf = if g.below(16) == 0 { g.draw(&eta_, -2.0, 5.0) } else { canon(z_r + h_rc + g.uni(-0.1, 0.1)) };
        let ux_l = g.draw(&US, -5.0, 5.0);
        let uy_l = g.draw(&US, -5.0, 5.0);
        let ux_r = g.draw(&US, -5.0, 5.0);
        let uy_r = g.draw(&US, -5.0, 5.0);
        let (nx, ny) = normal(&mut g);
        let hd = g.draw(&hd_, 0.0, 1e-2);
        let mut f = neg7();
        let (mut a, mut b, mut c, mut d) = (-7.0, -7.0, -7.0, -7.0);
        let r = face_flux_recon_bed(
            opq(eta_lf),
            opq(ux_l),
            opq(uy_l),
            opq(z_l),
            opq(h_lc),
            opq(eta_rf),
            opq(ux_r),
            opq(uy_r),
            opq(z_r),
            opq(h_rc),
            opq(nx),
            opq(ny),
            opq(hd),
            &mut f,
            &mut a,
            &mut b,
            &mut c,
            &mut d,
            opq(z_lf),
            opq(z_rf),
        );
        face_flux_line(
            &mut o,
            "frb",
            &[eta_lf, ux_l, uy_l, z_l, h_lc, eta_rf, ux_r, uy_r, z_r, h_rc, nx, ny, hd, z_lf, z_rf],
            &f,
            a,
            b,
            c,
            d,
            r,
        );
        let mut f2 = neg7();
        let (mut a, mut b, mut c, mut d) = (-7.0, -7.0, -7.0, -7.0);
        let r = face_flux_recon(
            opq(eta_lf),
            opq(ux_l),
            opq(uy_l),
            opq(z_l),
            opq(h_lc),
            opq(eta_rf),
            opq(ux_r),
            opq(uy_r),
            opq(z_r),
            opq(h_rc),
            opq(nx),
            opq(ny),
            opq(hd),
            &mut f2,
            &mut a,
            &mut b,
            &mut c,
            &mut d,
        );
        face_flux_line(
            &mut o,
            "frc",
            &[eta_lf, ux_l, uy_l, z_l, h_lc, eta_rf, ux_r, uy_r, z_r, h_rc, nx, ny, hd],
            &f2,
            a,
            b,
            c,
            d,
            r,
        );
    }

    // ---- bjLimiter
    let bw = [0.0, -0.0, 1.0, 1e-14, -1e-14, 2e-14, -2e-14, 0.5, -1.0, 2.0];
    grid(&[bw.len(), bw.len(), bw.len(), bw.len()], |i| {
        let (a, b, c, d) = (opq(bw[i[0]]), opq(bw[i[1]]), opq(bw[i[2]]), opq(bw[i[3]]));
        o.tag("bj");
        args(&mut o, &[a, b, c, d]);
        o.hx(bj_limiter(a, b, c, d));
        o.nl();
    });
    for _ in 0..n {
        let w = g.draw(&bw, -2.0, 2.0);
        let wf = if g.below(4) == 0 { w + g.uni(-3e-14, 3e-14) } else { g.draw(&bw, -3.0, 3.0) };
        let wmin = if g.below(4) == 0 { g.draw(&bw, -3.0, 3.0) } else { w - g.uni(0.0, 1.0) };
        let wmax = if g.below(4) == 0 { g.draw(&bw, -3.0, 3.0) } else { w + g.uni(0.0, 1.0) };
        o.tag("bj");
        args(&mut o, &[w, wf, wmin, wmax]);
        o.hx(bj_limiter(opq(w), opq(wf), opq(wmin), opq(wmax)));
        o.nl();
    }

    // ---- frictionUpdate
    let nm = [0.0, 0.02, 0.05, 0.1, -0.03];
    let dts = [0.0, -0.0, 0.1, 1.0, 10.0, -1.0];
    for _ in 0..2 * n {
        let mut qx = g.draw(&qx_, -5.0, 5.0);
        let mut qy = g.draw(&qx_, -5.0, 5.0);
        let h = g.draw(&HS, 0.0, 3.0);
        let nn = g.draw(&nm, 0.0, 0.1);
        let dt = g.draw(&dts, 0.0, 30.0);
        o.tag("fr");
        args(&mut o, &[qx, qy, h, nn, dt]);
        qx = opq(qx);
        qy = opq(qy);
        friction_update(&mut qx, &mut qy, opq(h), opq(nn), opq(dt));
        o.hx(qx);
        o.hx(qy);
        o.nl();
    }

    // ---- criticalDepth
    for _ in 0..n {
        let q = g.draw(&qx_, -10.0, 10.0);
        o.tag("cd");
        args(&mut o, &[q]);
        o.hx(critical_depth(opq(q)));
        o.nl();
    }

    // ---- depthFromInvariantAndDischarge
    let rv = [0.0, 1.0, 5.0, 10.0, -1.0, 2.0 * 9.80665f64.sqrt(), 6.264, 20.0];
    let qv = [0.0, -0.0, -1.0, 1.0, -0.5, 0.5, -5.0, 5.0, 1e-9];
    let hv = [0.0, -0.0, 1.0, 0.5, 2.0, 1e-6, -1.0, 1e-12];
    grid(&[rv.len(), qv.len(), hv.len()], |i| {
        let (a, b, c) = (opq(rv[i[0]]), opq(qv[i[1]]), opq(hv[i[2]]));
        o.tag("dri");
        args(&mut o, &[a, b, c]);
        o.hx(depth_from_invariant_and_discharge(a, b, c));
        o.nl();
    });
    for _ in 0..2 * n {
        let h = g.draw(&hv, 0.0, 3.0);
        let u = g.draw(&US, -3.0, 3.0);
        let r = if g.below(4) == 0 { g.draw(&rv, -5.0, 20.0) } else { canon(u + 2.0 * (9.80665 * h).sqrt()) };
        let q = if g.below(4) == 0 { g.draw(&qv, -5.0, 5.0) } else { canon(u * h * g.uni(0.5, 1.5)) };
        o.tag("dri");
        args(&mut o, &[r, q, h]);
        o.hx(depth_from_invariant_and_discharge(opq(r), opq(q), opq(h)));
        o.nl();
    }

    o.flush();
    eprintln!("twod_swe: {} lines", o.lines);
}
