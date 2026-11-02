// Rust side of gate/diff/xslookup.cpp: reads the file `xslookup gen` wrote, calls every
// translated function of hydraulics::xsect_lookup with the same arguments, and prints the
// same line format. Inputs on which upstream is undefined print "-" on both sides.
#![allow(clippy::int_plus_one)] // the harness's `v <= n - 1`, kept as written there
use std::io::{Read, Write};
use swmm6_rust::hydraulics::transect_types::{TransectData, N_TRANSECT_TBL};
use swmm6_rust::hydraulics::xsect_lookup::*;

fn hx(s: &mut String, v: f64) {
    s.push_str(&format!("{:016x} ", v.to_bits()));
}

// Same predicates as the harness: upstream's lookup_exact is defined (the cast fits an
// int and every read stays in the table), and lut_bucket's cast fits an int.
fn ok_lookup(x: f64, n: i32) -> bool {
    let delta = 1.0 / ((n - 1) as f64);
    let q = x / delta;
    if !(q > -1.0 && q < 2147483647.0) {
        return false;
    }
    let i = q as i32;
    if i >= n - 1 {
        return true;
    }
    if i < 2 && i + 2 > n - 1 {
        return false;
    }
    true
}
fn ok_bucket(l: &LocateLut, y: f64) -> bool {
    let v = (y - l.t0) * l.scale;
    v > -2147483648.0 && v < 2147483647.0
}
fn dump_lut(s: &mut String, l: &LocateLut) {
    hx(s, l.t0);
    hx(s, l.scale);
    s.push_str(&format!("{} ", l.j_last));
    for b in 0..=LocateLut::K_BUCKETS as usize {
        s.push_str(&format!("{:02x}", l.lo[b]));
    }
    s.push(' ');
}

fn main() {
    let mut out = String::new();
    {
        let td = TransectData::default();
        out.push_str(&format!("TD {} {} ", N_TRANSECT_TBL, td.name.len()));
        for v in [td.n_left, td.n_right, td.n_channel, td.length_factor] {
            hx(&mut out, v);
        }
        out.push_str(&format!("{} {} ", td.stations.len(), td.elevations.len()));
        for v in [td.x_left_bank, td.x_right_bank, td.y_full, td.w_max, td.a_full, td.r_full, td.s_max, td.a_max] {
            hx(&mut out, v);
        }
        for j in 0..N_TRANSECT_TBL as usize {
            hx(&mut out, td.area_tbl[j]);
            hx(&mut out, td.hrad_tbl[j]);
            hx(&mut out, td.width_tbl[j]);
        }
        dump_lut(&mut out, &td.area_lut);
        out.push_str(&format!("{}\n", LocateLut::K_BUCKETS));
    }
    let params = [0.0, -1.0, f64::NAN, 0.37, 1.0, 2.5, f64::INFINITY];

    let mut b = Vec::new();
    std::fs::File::open(std::env::args().nth(1).expect("usage: xslookup FILE"))
        .unwrap()
        .read_to_end(&mut b)
        .unwrap();
    let mut o = 0usize;
    let rd_i32 = |o: &mut usize| {
        let v = i32::from_le_bytes(b[*o..*o + 4].try_into().unwrap());
        *o += 4;
        v
    };
    let rd_f64 = |o: &mut usize, k: usize| -> Vec<f64> {
        let v = (0..k).map(|j| f64::from_le_bytes(b[*o + 8 * j..*o + 8 * j + 8].try_into().unwrap())).collect();
        *o += 8 * k;
        v
    };
    let nt = rd_i32(&mut o);
    let mut total = 0usize;
    let stdout = std::io::stdout();
    let mut w = std::io::BufWriter::new(stdout.lock());
    for tix in 0..nt {
        let n = rd_i32(&mut o);
        let t = rd_f64(&mut o, n as usize);
        let nq = rd_i32(&mut o);
        let q = rd_f64(&mut o, nq as usize);

        let mut jl: Vec<i32> = Vec::new();
        for v in [0, 1, 2, n / 2, n - 1] {
            if v >= 0 && v <= n - 1 && !jl.contains(&v) {
                jl.push(v);
            }
        }
        let mut ls = vec![LocateLut::default(); jl.len()];
        out.push_str(&format!("T {} {} ", tix, n));
        for k in 0..jl.len() {
            build_locate_lut(&mut ls[k], Some(&t), jl[k]);
            dump_lut(&mut out, &ls[k]);
        }
        let mut lnull = LocateLut { scale: 9.0, ..Default::default() }; // build must reset it
        build_locate_lut(&mut lnull, None, n - 1);
        dump_lut(&mut out, &lnull);
        let mut linv = LocateLut::default();
        let mut linv_null = LocateLut { t0: 3.0, ..Default::default() };
        build_invlookup_lut(&mut linv, Some(&t), n);
        build_invlookup_lut(&mut linv_null, None, n);
        dump_lut(&mut out, &linv);
        dump_lut(&mut out, &linv_null);
        let mut neff = n;
        if n >= 4 && t[(n - 3) as usize] > t[(n - 1) as usize] {
            neff = n - 2;
        }
        out.push_str(&format!("{}\n", neff));
        let lfull = *ls.last().unwrap(); // built for jLast == n-1

        for &x in &q {
            hx(&mut out, x);
            if ok_lookup(x, n) {
                hx(&mut out, lookup_exact(x, &t, n));
            } else {
                out.push_str("- ");
            }
            for &p in &params {
                let xn = norm_x(x, p);
                hx(&mut out, xn);
                if ok_lookup(xn, n) {
                    hx(&mut out, norm_lookup(x, p, &t, n));
                } else {
                    out.push_str("- ");
                }
            }
            out.push_str("| ");
            for k in 0..jl.len() {
                let j = jl[k];
                out.push_str(&format!(
                    "{} {} {} {} ",
                    locate_bisect(x, &t, j),
                    locate_maybe_lut(x, &t, j, Some(&ls[k])),
                    locate_maybe_lut(x, &t, j, None),
                    locate_maybe_lut(x, &t, j, Some(&lfull))
                ));
                if ls[k].scale > 0.0 {
                    out.push_str(&format!("{} ", locate_lut(x, &t, j, &ls[k])));
                } else {
                    out.push_str("- ");
                }
                if ok_bucket(&ls[k], x) {
                    out.push_str(&format!("{} ", lut_bucket(&ls[k], x)));
                } else {
                    out.push_str("- ");
                }
            }
            if n >= 4 {
                out.push_str(&format!("{} ", locate_maybe_lut(x, &t, neff - 1, Some(&linv))));
            }
            out.push('\n');
            total += 1;
        }
        w.write_all(out.as_bytes()).unwrap();
        out.clear();
    }
    w.flush().unwrap();
    eprintln!("{} tables, {} queries", nt, total);
}
