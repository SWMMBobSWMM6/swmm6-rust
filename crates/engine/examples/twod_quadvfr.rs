// Rust side of gate/diff/twod_quadvfr.cpp: the same inputs through twod::mesh::quad_vfr, printed
// the same way.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
#[path = "common/twodgen.rs"]
mod twodgen;
use swmm6_rust::twod::mesh::quad_vfr::*;
use swmm6_rust::twod::mesh::vfr_closure::vfr_sort3;
use twodgen::*;

const EPS: [f64; 11] = [0.0, -0.0, 1e-6, 1e-3, 0.01, 0.05, 0.2, 0.5, 0.99, 1.0, -0.1];

fn quad_funcs(o: &mut Out, g: &mut Gen, zs_in: &[f64], a1: f64, a2: f64, nstage: i32) {
    let mut zs = [0.0f64; 6];
    for k in 0..6 {
        zs[k] = opq(zs_in[k]);
    }
    let a1 = opq(a1);
    let a2 = opq(a2);
    let mut st = vec![-INF, INF];
    for k in 0..6 {
        st.push(zs[k]);
        st.push(zs[k].next_up());
        st.push(zs[k].next_down());
        st.push(zs[k] + 0.25);
        st.push(zs[k] - 0.25);
    }
    let zw = canon((a1 * (zs[0] + zs[1] + zs[2]) / 3.0 + a2 * (zs[3] + zs[4] + zs[5]) / 3.0) / (a1 + a2));
    st.push(zw);
    for r in 0..nstage {
        let eta = opq(if (r as usize) < st.len() { st[r as usize] } else { g.uni(zs[0] - 1.0, zs[5] + 1.0) });
        let eps = opq(g.pick(&EPS));
        let pre = |o: &mut Out, tag: &str| {
            o.tag(tag);
            for k in 0..6 {
                o.hx(zs[k]);
            }
            o.hx(a1);
            o.hx(a2);
        };
        pre(o, "qwf");
        o.hx(eta);
        o.sep();
        o.hx(quad_wet_fraction(&zs, a1, a2, eta));
        o.nl();
        pre(o, "qmd");
        o.hx(eta);
        o.hx(eps);
        let d = quad_mean_depth_from_eta(&zs, a1, a2, eta, eps);
        o.sep();
        o.hx(d);
        o.hx(quad_eta_from_mean_depth(&zs, a1, a2, opq(d), eps));
        o.nl();
        pre(o, "qsl");
        o.hx(eta);
        o.hx(eps);
        o.sep();
        o.hx(quad_mean_depth_slope(&zs, a1, a2, eta, eps));
        o.nl();
        pre(o, "qde");
        o.hx(eta);
        o.hx(eps);
        o.sep();
        o.hx(quad_d_eta_d_mean_depth(&zs, a1, a2, eta, eps));
        o.nl();
        let ztop = if zs[2] > zs[5] { zs[2] } else { zs[5] };
        let u = g.uni(0.0, 3.0);
        let tw = canon(ztop - zw);
        let dsel = [0.0, -0.0, -1.0, DMIN, 1e-9, 1e-3, tw, tw.next_down(), tw.next_up(), u, INF];
        let md = opq(dsel[g.below(11) as usize]);
        pre(o, "qef");
        o.hx(md);
        o.hx(eps);
        o.sep();
        o.hx(quad_eta_from_mean_depth(&zs, a1, a2, md, eps));
        o.nl();
    }
}

fn precompute(o: &mut Out, g: &mut Gen, x: &[f64], y: &[f64], z: &[f64], nstage: i32) {
    let mut xs = [0.0f64; 4];
    let mut ys = [0.0f64; 4];
    let mut zz = [0.0f64; 4];
    let mut zs = [0.0f64; 6];
    let mut a1 = -7.0;
    let mut a2 = -7.0;
    for k in 0..4 {
        xs[k] = opq(x[k]);
        ys[k] = opq(y[k]);
        zz[k] = opq(z[k]);
    }
    let kase = quad_vfr_precompute(&xs, &ys, &zz, &mut zs, &mut a1, &mut a2);
    o.tag("pre");
    for k in 0..4 {
        o.hx(xs[k]);
        o.hx(ys[k]);
        o.hx(zz[k]);
    }
    o.sep();
    for k in 0..6 {
        o.hx(zs[k]);
    }
    o.hx(a1);
    o.hx(a2);
    o.hi(kase as i64);
    o.nl();
    quad_funcs(o, g, &zs, a1, a2, nstage);
}

/// std::next_permutation over a small int array.
fn next_permutation(p: &mut [i32]) -> bool {
    let n = p.len();
    if n < 2 {
        return false;
    }
    let mut i = n - 1;
    while i > 0 && p[i - 1] >= p[i] {
        i -= 1;
    }
    if i == 0 {
        p.reverse();
        return false;
    }
    let mut j = n - 1;
    while p[j] <= p[i - 1] {
        j -= 1;
    }
    p.swap(i - 1, j);
    p[i..].reverse();
    true
}

fn main() {
    let mut o = Out::new();
    let mut g = Gen(0x2d0f1a7e5eed0002);
    let sqx = [0.0, 1.0, 1.0, 0.0];
    let sqy = [0.0, 0.0, 1.0, 1.0];
    let rvx = [0.0, 0.0, 1.0, 1.0];
    let rvy = [0.0, 1.0, 1.0, 0.0];
    let skx = [0.0, 2.5, 3.1, 0.4];
    let sky = [0.0, 0.2, 1.7, 1.3];
    let lnx = [0.0, 1.0, 2.0, 3.0];
    let lny = [0.0, 1.0, 2.0, 3.0000001];
    let xs: [&[f64]; 4] = [&sqx, &rvx, &skx, &lnx];
    let ys: [&[f64]; 4] = [&sqy, &rvy, &sky, &lny];
    let zv = [0.0, 1.0, 2.0];
    for q in 0..4 {
        grid(&[3, 3, 3, 3], |i| {
            let z = [zv[i[0]], zv[i[1]], zv[i[2]], zv[i[3]]];
            precompute(&mut o, &mut g, xs[q], ys[q], &z, 40);
        });
    }
    let l4 = [10.0, 10.5, 11.75, 13.0];
    let mut perm = [0i32, 1, 2, 3];
    loop {
        let z = [l4[perm[0] as usize], l4[perm[1] as usize], l4[perm[2] as usize], l4[perm[3] as usize]];
        for q in 0..3 {
            precompute(&mut o, &mut g, xs[q], ys[q], &z, 40);
        }
        if !next_permutation(&mut perm) {
            break;
        }
    }
    let spz = [0.0, 1.0, 1e-9, 1e-10, 100.0];
    for _k in 0..6000 {
        let w = g.uni(0.1, 50.0);
        let h = g.uni(0.1, 50.0);
        let x0 = g.uni(-1e3, 1e3);
        let y0 = g.uni(-1e3, 1e3);
        let mut x = [x0, x0 + w, x0 + w, x0];
        let mut y = [y0, y0, y0 + h, y0 + h];
        for v in 0..4 {
            x[v] += g.uni(-0.2, 0.2) * w;
            y[v] += g.uni(-0.2, 0.2) * h;
        }
        let base = g.uni(-10.0, 500.0);
        let mut z = [0.0f64; 4];
        let mode = g.below(4);
        for v in 0..4 {
            z[v] = match mode {
                0 => base + g.uni(0.0, 3.0),
                1 => base + g.pick(&spz),
                2 => g.draw(&spz, -5.0, 5.0),
                _ => base + g.uni(0.0, 1e-8),
            };
        }
        if g.below(2) == 0 {
            x.swap(1, 3);
            y.swap(1, 3);
            z.swap(1, 3);
        }
        precompute(&mut o, &mut g, &x, &y, &z, 8);
    }
    for _k in 0..3000 {
        let mut zs = [0.0f64; 6];
        for v in 0..6 {
            zs[v] = g.draw(&spz, -5.0, 5.0);
        }
        if g.below(4) != 0 {
            let (mut a, mut b, mut c) = (zs[0], zs[1], zs[2]);
            vfr_sort3(&mut a, &mut b, &mut c);
            zs[0] = a;
            zs[1] = b;
            zs[2] = c;
            let (mut a, mut b, mut c) = (zs[3], zs[4], zs[5]);
            vfr_sort3(&mut a, &mut b, &mut c);
            zs[3] = a;
            zs[4] = b;
            zs[5] = c;
        }
        let sa = [0.0, -0.0, 1.0, 0.5, 1e-12, 1e6];
        let a1 = g.draw(&sa, 0.0, 100.0);
        let a2 = g.draw(&sa, 0.0, 100.0);
        quad_funcs(&mut o, &mut g, &zs, a1, a2, 4);
    }
    o.flush();
    eprintln!("twod_quadvfr: {} lines", o.lines);
}
