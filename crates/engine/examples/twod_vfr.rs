// Rust side of gate/diff/twod_vfr.cpp: the same inputs through twod::mesh::vfr_closure, printed
// the same way.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
#[path = "common/twodgen.rs"]
mod twodgen;
use swmm6_rust::twod::mesh::vfr_closure::*;
use twodgen::*;

#[derive(Clone, Copy)]
struct Tri {
    a: f64,
    b: f64,
    c: f64,
}
fn tri(a: f64, b: f64, c: f64) -> Tri {
    Tri { a, b, c }
}

fn structured() -> Vec<Tri> {
    vec![
        tri(0.0, 0.0, 0.0),
        tri(1.0, 1.0, 1.0),
        tri(-0.0, 0.0, 0.0),
        tri(0.0, 0.0, 1.0),
        tri(0.0, 1.0, 1.0),
        tri(0.0, 0.5, 1.0),
        tri(0.0, 0.25, 1.0),
        tri(0.0, 1e-10, 2e-10),
        tri(0.0, 5e-10, 1e-9),
        tri(0.0, 1e-9, 1.0),
        tri(0.0, 9.99e-10, 1.0),
        tri(0.0, 1.0 - 1e-9, 1.0),
        tri(0.0, 1.0 - 5e-10, 1.0),
        tri(100.0, 100.3, 101.7),
        tri(-5.0, -2.0, 3.0),
        tri(1000.0, 1000.0 + 1e-6, 1000.0 + 3e-6),
        tri(12345.678, 12345.679, 12345.68),
        tri(0.0, 2.0, 3.0),
        tri(-1e-3, 0.0, 1e-3),
        tri(0.0, 1e-6, 1e-3),
        tri(0.1, 0.2, 0.3),
        tri(1.0, 2.0, 1000.0),
        tri(-0.0, 1e-300, 2e-300),
        tri(0.0, 1e-9, 2e-9),
        tri(5.0, 5.0, 5.0 + 1e-9),
        tri(5.0, 5.0 + 1e-9, 5.0 + 1e-9),
        tri(DMIN, 2.0 * DMIN, 3.0 * DMIN),
        tri(3.0, 1.0, 2.0),
        tri(1.0, 0.0, 0.5),
    ]
}

fn stages(z: Tri) -> Vec<f64> {
    let z = tri(opq(z.a), opq(z.b), opq(z.c));
    let zbar = canon((z.a + z.b + z.c) / 3.0);
    let mut e = vec![
        z.a,
        z.b,
        z.c,
        zbar,
        canon(0.5 * (z.a + z.b)),
        canon(0.5 * (z.b + z.c)),
        z.a - 1.0,
        z.c + 1.0,
        z.a - 1e-12,
        z.c + 1e-12,
        -INF,
        INF,
    ];
    for v in [z.a, z.b, z.c] {
        e.push(v.next_up());
        e.push(v.next_down());
        e.push(v.next_up().next_up());
    }
    e
}

const EPS: [f64; 12] = [0.0, -0.0, 1e-6, 1e-4, 0.01, 0.05, 0.1, 0.3, 0.9, 1.0, -0.01, 0.999999];

fn line_tri(o: &mut Out, tag: &str, z: Tri) {
    o.tag(tag);
    o.hx(z.a);
    o.hx(z.b);
    o.hx(z.c);
}

fn one(o: &mut Out, z: Tri, eta: f64, eps: f64) {
    let (z1, z2, z3) = (opq(z.a), opq(z.b), opq(z.c));
    let eta = opq(eta);
    let eps = opq(eps);
    line_tri(o, "wf", z);
    o.hx(eta);
    o.sep();
    o.hx(vfr_wet_fraction(z1, z2, z3, eta));
    o.nl();
    line_tri(o, "mdx", z);
    o.hx(eta);
    o.sep();
    o.hx(vfr_mean_depth_from_eta_exact(z1, z2, z3, eta));
    o.nl();
    line_tri(o, "md", z);
    o.hx(eta);
    o.hx(eps);
    let d = vfr_mean_depth_from_eta(z1, z2, z3, eta, eps);
    o.sep();
    o.hx(d);
    o.hx(vfr_eta_from_mean_depth(z1, z2, z3, opq(d), eps));
    o.nl();
    line_tri(o, "de", z);
    o.hx(eta);
    o.hx(eps);
    o.sep();
    o.hx(vfr_d_eta_d_mean_depth(z1, z2, z3, eta, eps));
    o.nl();
}

fn per_eps(o: &mut Out, z: Tri, eps: f64) {
    let (z1, z2, z3) = (opq(z.a), opq(z.b), opq(z.c));
    let eps = opq(eps);
    line_tri(o, "st", z);
    o.hx(eps);
    o.sep();
    o.hx(vfr_stage_at_wet_fraction(z1, z2, z3, eps));
    o.nl();
    line_tri(o, "dry", z);
    o.hx(eps);
    o.sep();
    o.hx(vfr_dry_eta(z1, z2, z3, eps));
    o.nl();
    let z = tri(opq(z.a), opq(z.b), opq(z.c));
    let zbar = (z.a + z.b + z.c) / 3.0;
    let relief = z.c - z.a;
    let d21 = z.b - z.a;
    let mut ds = vec![0.0, -0.0, -1.0, DMIN, 1e-12, 1e-6, 0.5, z.c - zbar, d21 * d21 / (3.0 * relief), 1e3, INF, -INF];
    if eps > 0.0 && relief >= K_VFR_FLAT_RELIEF {
        let es = vfr_stage_at_wet_fraction(z.a, z.b, z.c, eps);
        ds.push(vfr_mean_depth_from_eta_exact(z.a, z.b, z.c, es));
    }
    for d in ds.iter_mut() {
        *d = canon(*d);
    }
    let n0 = ds.len();
    for i in 4..n0 {
        let v = ds[i];
        ds.push(v.next_up());
        ds.push(v.next_down());
    }
    for d in ds {
        let d = opq(d);
        line_tri(o, "ef", z);
        o.hx(d);
        o.hx(eps);
        o.sep();
        o.hx(vfr_eta_from_mean_depth(z1, z2, z3, d, eps));
        o.nl();
    }
}

fn main() {
    let mut o = Out::new();
    let mut g = Gen(0x2d0f1a7e5eed0001);
    let sz: [f64; 9] = [-1.0, 0.0, -0.0, 1e-10, 1.0, 2.0, 1.0 + 1e-12, 5.0, 100.0];
    grid(&[sz.len(), sz.len(), sz.len()], |i| {
        let (mut a, mut b, mut c) = (opq(sz[i[0]]), opq(sz[i[1]]), opq(sz[i[2]]));
        o.tag("sort");
        o.hx(a);
        o.hx(b);
        o.hx(c);
        o.sep();
        vfr_sort3(&mut a, &mut b, &mut c);
        o.hx(a);
        o.hx(b);
        o.hx(c);
        o.nl();
    });
    for z in structured() {
        for eps in EPS {
            per_eps(&mut o, z, eps);
            for eta in stages(z) {
                one(&mut o, z, eta, eps);
            }
        }
    }
    let spz = [0.0, -0.0, 1.0, 1e-9, 1e-10, 100.0, -3.0];
    for _k in 0..6000 {
        let mut z = match g.below(4) {
            0 => {
                let a = g.uni(-10.0, 100.0);
                let b = g.uni(-10.0, 100.0);
                let c = g.uni(-10.0, 100.0);
                tri(a, b, c)
            }
            1 => {
                let b = g.uni(-50.0, 500.0);
                let p = g.pos_ebits(1023 - 35, 1023 + 2);
                let q = g.pos_ebits(1023 - 35, 1023 + 2);
                tri(b, b + p, b + q)
            }
            2 => {
                let a = g.draw(&spz, -10.0, 10.0);
                let b = g.draw(&spz, -10.0, 10.0);
                let c = g.draw(&spz, -10.0, 10.0);
                tri(a, b, c)
            }
            _ => {
                let b = g.uni(0.0, 2000.0);
                let p = g.uni(0.0, 1e-8);
                let q = g.uni(0.0, 3.0);
                tri(b, b + p, b + q)
            }
        };
        if g.below(16) != 0 {
            vfr_sort3(&mut z.a, &mut z.b, &mut z.c);
        }
        let st = stages(z);
        for r in 0..4 {
            let eps = if g.below(3) == 0 { g.pick(&EPS) } else { g.uni(0.0, 0.5) };
            let eta = match g.below(3) {
                0 => g.pick(&st),
                1 => g.uni(z.a - 1.0, z.c + 1.0),
                _ => {
                    let v = g.pick(&st);
                    g.nudge(v)
                }
            };
            if r == 0 {
                per_eps(&mut o, z, eps);
            }
            one(&mut o, z, eta, eps);
        }
    }
    o.flush();
    eprintln!("twod_vfr: {} lines", o.lines);
}
