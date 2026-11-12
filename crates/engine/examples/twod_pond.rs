// Rust side of gate/diff/twod_pond.cpp: the same inputs through twod::coupling::pond_model,
// printed the same way.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
#[path = "common/twodgen.rs"]
mod twodgen;
use swmm6_rust::twod::coupling::pond_model::*;
use twodgen::*;

fn args(o: &mut Out, a: &[f64]) {
    for &v in a {
        o.hx(v);
    }
    o.sep();
}

fn main() {
    let mut o = Out::new();
    let mut g = Gen(0x2d0f1a7e5eed0009);
    let n = 60000;
    let xx = [0.0, -0.0, -1.0, 0.05, 0.05f64.next_down(), 0.05f64.next_up(), 0.025, 1e-12, 1.0, DMIN, 0.1];
    let aq = [10.0, 0.0, -0.0, -5.0, 1.0, 1e6, 1e-9];
    let bd = [0.05, 0.0, -0.0, -0.05, 0.1, 1e-9, 1.0];
    grid(&[xx.len(), aq.len(), bd.len()], |i| {
        let (x, a, b) = (opq(xx[i[0]]), opq(aq[i[1]]), opq(bd[i[2]]));
        o.tag("pa");
        args(&mut o, &[x, a, b]);
        o.hx(area(x, a, b));
        o.hx(volume(x, a, b));
        o.nl();
        o.tag("pd");
        args(&mut o, &[x, a, b]);
        o.hx(depth_for_volume(x, a, b));
        o.nl();
    });
    for _ in 0..n {
        let a = g.draw(&aq, 0.0, 100.0);
        let b = g.draw(&bd, 0.0, 0.2);
        let x = if g.below(2) == 0 { g.draw(&xx, -0.01, 0.3) } else { canon(b * g.uni(0.0, 1.2)) };
        o.tag("pa");
        args(&mut o, &[x, a, b]);
        let ar = area(opq(x), opq(a), opq(b));
        let vo = volume(opq(x), opq(a), opq(b));
        o.hx(ar);
        o.hx(vo);
        o.hx(depth_for_volume(opq(vo), opq(a), opq(b)));
        o.nl();
        let v = if g.below(2) == 0 { g.draw(&xx, -0.01, 3.0) } else { canon(0.5 * a * b * g.uni(0.0, 1.1)) };
        o.tag("pd");
        args(&mut o, &[v, a, b]);
        o.hx(depth_for_volume(opq(v), opq(a), opq(b)));
        o.nl();
    }
    let dh = [0.0, -0.0, 1.0, -1.0, 1e-12, -0.3];
    let wn = [0.0, -0.0, 1.0, 30.0, -1.0, 1e-9];
    let a1 = [0.0, -0.0, 1.0, 12.5, -2.0, 1e-9];
    grid(&[dh.len(), wn.len(), a1.len()], |i| {
        let (a, b, c) = (opq(dh[i[0]]), opq(wn[i[1]]), opq(a1[i[2]]));
        o.tag("eq");
        args(&mut o, &[a, b, c]);
        o.hx(equalize_cap(a, b, c));
        o.nl();
    });
    for _ in 0..n {
        let a = g.draw(&dh, -2.0, 2.0);
        let b = g.draw(&wn, 0.0, 60.0);
        let c = g.draw(&a1, 0.0, 50.0);
        o.tag("eq");
        args(&mut o, &[a, b, c]);
        o.hx(equalize_cap(opq(a), opq(b), opq(c)));
        o.nl();
    }
    o.flush();
    eprintln!("twod_pond: {} lines", o.lines);
}
