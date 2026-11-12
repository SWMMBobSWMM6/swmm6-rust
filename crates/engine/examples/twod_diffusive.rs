// Rust side of gate/diff/twod_diffusive.cpp: the same inputs through
// twod::solver::diffusive_kernels, printed the same way.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
#[path = "common/twodgen.rs"]
mod twodgen;
use swmm6_rust::twod::solver::diffusive_kernels::*;
use twodgen::*;

fn args(o: &mut Out, a: &[f64]) {
    for &v in a {
        o.hx(v);
    }
    o.sep();
}

fn main() {
    let mut o = Out::new();
    let mut g = Gen(0x2d0f1a7e5eed0005);
    let n = 60000;
    let hf = [0.0, -0.0, 1e-12, 1e-3, 0.3, 1.0, 4.0, -1.0, DMIN, 1e-200];
    let sl = [0.0, -0.0, 1e-6, -1e-6, 1e-6f64.next_down(), 1e-6f64.next_up(), 1e-3, -1e-3, 0.2, -0.2];
    let nf = [0.0, -0.0, 0.03, 0.1, -0.02, 1e-9];
    let se = [1e-6, 0.0, 1e-4, 1e-9];
    grid(&[hf.len(), sl.len(), nf.len(), se.len()], |i| {
        let (a, b, c, d) = (opq(hf[i[0]]), opq(sl[i[1]]), opq(nf[i[2]]), opq(se[i[3]]));
        o.tag("fq");
        args(&mut o, &[a, b, c, d]);
        o.hx(face_discharge(a, b, c, d));
        o.nl();
    });
    for _ in 0..n {
        let a = g.draw(&hf, 0.0, 3.0);
        let b = g.draw(&sl, -0.05, 0.05);
        let c = g.draw(&nf, 0.01, 0.1);
        let d = g.draw(&se, 0.0, 1e-4);
        o.tag("fq");
        args(&mut o, &[a, b, c, d]);
        o.hx(face_discharge(opq(a), opq(b), opq(c), opq(d)));
        o.nl();
    }
    let al = [0.25, 0.5, 1.0, 0.0];
    let lc = [1.0, 10.0, 0.0, 1e-3];
    let sa = [0.0, 1e-6, 1e-6f64.next_down(), 1e-6f64.next_up(), 1e-3, 0.2, -0.1];
    grid(&[al.len(), lc.len(), hf.len(), nf.len(), sa.len(), se.len()], |i| {
        let (a, b, c) = (opq(al[i[0]]), opq(lc[i[1]]), opq(hf[i[2]]));
        let (d, e, f) = (opq(nf[i[3]]), opq(sa[i[4]]), opq(se[i[5]]));
        o.tag("dt");
        args(&mut o, &[a, b, c, d, e, f]);
        o.hx(cell_diffusive_dt(a, b, c, d, e, f));
        o.nl();
    });
    for _ in 0..n {
        let a = g.draw(&al, 0.1, 1.0);
        let b = g.draw(&lc, 0.1, 100.0);
        let c = g.draw(&hf, 0.0, 3.0);
        let d = g.draw(&nf, 0.01, 0.1);
        let e = g.draw(&sa, 0.0, 0.1);
        let f = g.draw(&se, 0.0, 1e-4);
        o.tag("dt");
        args(&mut o, &[a, b, c, d, e, f]);
        o.hx(cell_diffusive_dt(opq(a), opq(b), opq(c), opq(d), opq(e), opq(f)));
        o.nl();
    }
    o.flush();
    eprintln!("twod_diffusive: {} lines", o.lines);
}
