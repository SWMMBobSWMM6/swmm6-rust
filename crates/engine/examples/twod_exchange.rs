// Rust side of gate/diff/twod_exchange.cpp: the same inputs through
// twod::coupling::exchange_law, printed the same way.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
#[path = "common/twodgen.rs"]
mod twodgen;
use swmm6_rust::twod::coupling::exchange_law::*;
use twodgen::*;

fn args(o: &mut Out, a: &[f64]) {
    for &v in a {
        o.hx(v);
    }
    o.sep();
}

fn main() {
    let mut o = Out::new();
    let mut g = Gen(0x2d0f1a7e5eed0008);
    let n = 60000;
    let xs = [0.0, -0.0, -1.0, 0.02, 0.02f64.next_down(), 0.02f64.next_up(), 0.01, 1e-12, 1.0, 4.0, DMIN, 1e-3];
    let es = [0.02, 0.0, 1e-3, 0.1, 1e-9, -0.02, 1.0];
    grid(&[xs.len(), es.len()], |i| {
        let (x, e) = (opq(xs[i[0]]), opq(es[i[1]]));
        o.tag("phi");
        args(&mut o, &[x, e]);
        o.hx(phi(x, e));
        o.hx(dphi(x, e));
        o.nl();
    });
    for _ in 0..n {
        let x = g.draw(&xs, -0.05, 0.5);
        let e = g.draw(&es, 1e-4, 0.1);
        o.tag("phi");
        args(&mut o, &[x, e]);
        o.hx(phi(opq(x), opq(e)));
        o.hx(dphi(opq(x), opq(e)));
        o.nl();
    }
    let ts = [0.0, -0.0, 1.0, 0.5, -0.5, 1.5, 1e-12, 1.0f64.next_down(), 1.0f64.next_up(), DMIN];
    for t in ts {
        o.tag("ss");
        args(&mut o, &[t]);
        o.hx(smoothstep(opq(t)));
        o.nl();
    }
    for _ in 0..n {
        let t = g.draw(&ts, -0.2, 1.2);
        o.tag("ss");
        args(&mut o, &[t]);
        o.hx(smoothstep(opq(t)));
        o.nl();
    }
    let ab = [0.0, -0.0, -1.0, 1.0, 2.0, 0.5, 1e-3, 1e-20, 1e20, 1e-300, 1e300, DMIN, 0.9999999];
    grid(&[ab.len(), ab.len()], |i| {
        let (a, b) = (opq(ab[i[0]]), opq(ab[i[1]]));
        let (q, wa, wb) = softmin(a, b);
        o.tag("sm");
        args(&mut o, &[a, b]);
        o.hx(q);
        o.hx(wa);
        o.hx(wb);
        o.nl();
    });
    for _ in 0..n {
        let a = g.draw(&ab, 0.0, 2.0);
        let b = if g.below(3) == 0 { canon(a * g.uni(0.5, 2.0)) } else { g.draw(&ab, 0.0, 2.0) };
        let (q, wa, wb) = softmin(opq(a), opq(b));
        o.tag("sm");
        args(&mut o, &[a, b]);
        o.hx(q);
        o.hx(wa);
        o.hx(wb);
        o.nl();
    }
    let zc = [0.0, 10.0, -0.0, 100.25];
    let lw = [0.0, 1.0, 3.14, 0.5];
    let ao = [0.0, 0.1, 0.785, 1e-4];
    let co = [0.54, 0.36, 0.167, 0.16, 0.0, 1.0];
    let he = [0.02, 1e-3, 0.1, 1e-9];
    let rh = [0.0, -0.0, 0.01, 0.02, 0.05, 0.1, 0.5, 1.0, -0.5, 1e-12, -1e-12, 2.0];
    for _ in 0..3 * n {
        let mut geo = Geometry::default();
        geo.z_crest = g.draw(&zc, -5.0, 50.0);
        geo.l_weir = g.draw(&lw, 0.0, 5.0);
        geo.a_orif = g.draw(&ao, 0.0, 2.0);
        let mut c = Coefficients::default();
        if g.below(3) != 0 {
            c.c_w = g.draw(&co, 0.0, 1.0);
            c.c_sw = g.draw(&co, 0.0, 1.0);
            c.c_oi = g.draw(&co, 0.0, 1.0);
            c.c_oo = g.draw(&co, 0.0, 1.0);
            c.h_eps = g.draw(&he, 1e-4, 0.1);
        }
        let z = geo.z_crest;
        let (h2, h1);
        match g.below(6) {
            0 => {
                h2 = canon(z + g.pick(&rh));
                h1 = canon(z + g.pick(&rh));
            }
            1 => {
                h2 = canon(z + g.uni(0.0, 2.0));
                h1 = canon(z + g.uni(-1.0, 2.0));
            }
            2 => {
                let hs = g.uni(1e-3, 1.0);
                h2 = canon(z + hs);
                h1 = canon(z + hs * g.uni(0.0, 0.12));
            }
            3 => {
                h2 = canon(z + g.uni(-1.0, 2.0));
                h1 = h2;
            }
            4 => {
                h2 = g.draw(&rh, -1.0, 60.0);
                h1 = g.draw(&rh, -1.0, 60.0);
            }
            _ => {
                let hd = g.uni(1e-3, 1.0);
                h1 = canon(z + hd);
                h2 = canon(z + hd * g.uni(-0.5, 0.15));
            }
        }
        let h2 = opq(h2);
        let h1 = opq(h1);
        let gq = Geometry { z_crest: opq(geo.z_crest), l_weir: opq(geo.l_weir), a_orif: opq(geo.a_orif) };
        let cq = Coefficients { c_w: opq(c.c_w), c_sw: opq(c.c_sw), c_oi: opq(c.c_oi), c_oo: opq(c.c_oo), h_eps: opq(c.h_eps) };
        let r = exchange(&gq, &cq, h2, h1);
        o.tag("ex");
        args(&mut o, &[geo.z_crest, geo.l_weir, geo.a_orif, c.c_w, c.c_sw, c.c_oi, c.c_oo, c.h_eps, h2, h1]);
        o.hx(r.q);
        o.hx(r.g);
        o.nl();
    }
    o.flush();
    eprintln!("twod_exchange: {} lines", o.lines);
}
