// Rust side of gate/diff/twod_momentum.cpp: the same inputs through
// twod::solver::momentum_physics_2d::MomentumPhysics2D (its Default and rotate), printed the same
// way.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
#[path = "common/twodgen.rs"]
mod twodgen;
use swmm6_rust::twod::solver::momentum_physics_2d::MomentumPhysics2D;
use twodgen::*;

fn main() {
    let mut o = Out::new();
    let mut g = Gen(0x2d0f1a7e5eed0007);
    {
        let mp = MomentumPhysics2D::default();
        o.tag("def");
        o.sep();
        o.hx(mp.last_limit);
        o.hi(mp.substeps as i64);
        o.hi(mp.retries as i64);
        o.hx(MomentumPhysics2D::KAPPA);
        o.hi((mp.viscosity.len() + mp.strain.len() + mp.vorticity.len() + mp.wall_distance.len()) as i64);
        o.nl();
    }
    let pi = 3.141592653589793f64;
    let an = [
        0.0,
        -0.0,
        pi,
        -pi,
        pi / 2.0,
        pi / 4.0,
        2.0 * pi,
        1e-8,
        -1e-8,
        1e-300,
        DMIN,
        1.0,
        -1.0,
        7.292e-5,
        1e6,
        1e22,
        -1e22,
        1.7e308,
        INF,
        -INF,
        QNAN,
        0.5,
        100.0,
        3.0e-3,
    ];
    let qv = [0.0, -0.0, 1.0, -1.0, 2.5, 1e-12, -3.0];
    let call = |g: &mut Gen, o: &mut Out, angle: f64, n: i32| {
        let mut qx = vec![1.0, 0.0];
        let mut qy = vec![0.0, 1.0];
        for _ in 0..n {
            qx.push(g.draw(&qv, -5.0, 5.0));
            qy.push(g.draw(&qv, -5.0, 5.0));
        }
        o.tag("rot");
        o.hx(angle);
        for k in 0..qx.len() {
            o.hx(qx[k]);
            o.hx(qy[k]);
        }
        o.sep();
        for k in 0..qx.len() {
            qx[k] = opq(qx[k]);
            qy[k] = opq(qy[k]);
        }
        MomentumPhysics2D::rotate(&mut qx, &mut qy, opq(angle));
        for k in 0..qx.len() {
            o.hx(qx[k]);
            o.hx(qy[k]);
        }
        o.nl();
    };
    for a in an {
        call(&mut g, &mut o, a, 6);
    }
    for _ in 0..60000 {
        let a = match g.below(5) {
            0 => g.uni(-2.0 * pi, 2.0 * pi),
            1 => g.ebits(1023 - 30, 1023 + 3),
            2 => g.ebits(1, 2046),
            3 => {
                let v = g.pick(&an);
                g.nudge(v)
            }
            _ => g.uni(-1e-3, 1e-3),
        };
        let n = g.below(9) as i32;
        call(&mut g, &mut o, a, n);
    }
    o.flush();
    eprintln!("twod_momentum: {} lines", o.lines);
}
