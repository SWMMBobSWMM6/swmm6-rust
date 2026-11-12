// Rust side of gate/diff/twod_sinks.cpp: the same inputs through
// twod::solver::surface_flux_calculator::{evap_sink, infil_sink}, printed the same way.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
#[path = "common/twodgen.rs"]
mod twodgen;
use swmm6_rust::twod::solver::surface_flux_calculator::{evap_sink, infil_sink};
use twodgen::*;

fn args(o: &mut Out, a: &[f64]) {
    for &v in a {
        o.hx(v);
    }
    o.sep();
}

fn main() {
    let mut o = Out::new();
    let mut g = Gen(0x2d0f1a7e5eed0006);
    let n = 50000;
    let rt = [0.0, -0.0, 1e-7, 2.7e-8, -1e-7, 1.0, DMIN];
    let dp = [0.0, -0.0, 1e-4, 1e-4f64.next_down(), 1e-4f64.next_up(), 5e-5, 1e-12, 1.0, -1.0, DMIN];
    let dd = [1e-4, 0.0, -0.0, 1e-3, 1e-12, -1e-4];
    grid(&[rt.len(), dp.len(), dd.len()], |i| {
        let (a, b, c) = (opq(rt[i[0]]), opq(dp[i[1]]), opq(dd[i[2]]));
        o.tag("ev");
        args(&mut o, &[a, b, c]);
        o.hx(evap_sink(a, b, c));
        o.nl();
        o.tag("in");
        args(&mut o, &[a, b, c]);
        o.hx(infil_sink(a, b, c));
        o.nl();
    });
    for _ in 0..n {
        let a = g.draw(&rt, 0.0, 1e-5);
        let b = g.draw(&dp, 0.0, 2e-4);
        let c = g.draw(&dd, 0.0, 1e-3);
        o.tag("ev");
        args(&mut o, &[a, b, c]);
        o.hx(evap_sink(opq(a), opq(b), opq(c)));
        o.nl();
        o.tag("in");
        args(&mut o, &[a, b, c]);
        o.hx(infil_sink(opq(a), opq(b), opq(c)));
        o.nl();
    }
    o.flush();
    eprintln!("twod_sinks: {} lines", o.lines);
}
