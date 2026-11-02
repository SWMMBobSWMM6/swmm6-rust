// Rust side of gate/diff/hydfuncs.cpp: reads the file `hydfuncs gen` wrote, runs every
// command with the translated hydraulics::{transect, street, force_main, hyd_closure_kernels,
// link, culvert} and math::find_root, and prints the same lines.
#![allow(clippy::too_many_arguments, clippy::needless_range_loop, clippy::field_reassign_with_default)]
use std::collections::HashMap;
use std::fmt::Write as _;
use std::io::Write;

use swmm6_rust::data::infra::TransectStore;
use swmm6_rust::data::link::{LinkData, LinkType, XsectShape};
use swmm6_rust::hydraulics::culvert;
use swmm6_rust::hydraulics::force_main::{self as forcemain, FrictionModel};
use swmm6_rust::hydraulics::hyd_closure_kernels::{self as hydkernels, CulvertCurve};
use swmm6_rust::hydraulics::link;
use swmm6_rust::hydraulics::street::{self, StreetParams};
use swmm6_rust::hydraulics::transect;
use swmm6_rust::hydraulics::transect_types::TransectData;
use swmm6_rust::hydraulics::xsect_batch::XSectParams;
use swmm6_rust::hydraulics::xsect_lookup::LocateLut;
use swmm6_rust::hydraulics::xsection as xsect;
use swmm6_rust::math::find_root as findroot;


// ---------------------------------------------------------------- output
struct Out(String);
impl Out {
    fn hx(&mut self, v: f64) {
        let _ = write!(self.0, "{:016x} ", v.to_bits());
    }
    fn hi(&mut self, v: i64) {
        let _ = write!(self.0, "{} ", v);
    }
    fn hs(&mut self, s: &str) {
        self.0.push_str(s);
        self.0.push(' ');
    }
    fn nl(&mut self) {
        self.0.push('\n');
        if self.0.len() > (1 << 22) {
            let _ = std::io::stdout().lock().write_all(self.0.as_bytes());
            self.0.clear();
        }
    }
    fn harr(&mut self, a: &[f64]) {
        for &v in a {
            self.hx(v);
        }
    }
    fn hname(&mut self, s: &[u8]) {
        if s.is_empty() {
            self.0.push_str("~ ");
            return;
        }
        for &c in s {
            let _ = write!(self.0, "{:02x}", c);
        }
        self.0.push(' ');
    }
}

// ---------------------------------------------------------------- input
fn rdh<'a, I: Iterator<Item = &'a str>>(it: &mut I) -> f64 {
    f64::from_bits(u64::from_str_radix(it.next().unwrap(), 16).unwrap())
}
fn rdi<'a, I: Iterator<Item = &'a str>>(it: &mut I) -> i64 {
    it.next().unwrap().parse().unwrap()
}
fn rdname<'a, I: Iterator<Item = &'a str>>(it: &mut I) -> Vec<u8> {
    let s = it.next().unwrap();
    if s == "~" {
        return Vec::new();
    }
    let b = s.as_bytes();
    let mut r = Vec::new();
    let mut i = 0;
    while i + 1 < b.len() {
        r.push(u8::from_str_radix(&s[i..i + 2], 16).unwrap());
        i += 2;
    }
    r
}
fn rdvec<'a, I: Iterator<Item = &'a str>>(it: &mut I) -> Vec<f64> {
    let n = rdi(it) as usize;
    (0..n).map(|_| rdh(it)).collect()
}

// ---------------------------------------------------------------- dumps
fn dump_lut(o: &mut Out, l: &LocateLut) {
    o.hx(l.t0);
    o.hx(l.scale);
    o.hi(l.j_last as i64);
    for b in 0..=LocateLut::K_BUCKETS as usize {
        o.hi(l.lo[b] as i64);
    }
}
fn dump_td(o: &mut Out, td: &TransectData) {
    o.hname(&td.name);
    o.hx(td.n_left);
    o.hx(td.n_right);
    o.hx(td.n_channel);
    o.hx(td.length_factor);
    o.hi(td.stations.len() as i64);
    o.harr(&td.stations);
    o.hi(td.elevations.len() as i64);
    o.harr(&td.elevations);
    o.hx(td.x_left_bank);
    o.hx(td.x_right_bank);
    o.hx(td.y_full);
    o.hx(td.w_max);
    o.hx(td.a_full);
    o.hx(td.r_full);
    o.hx(td.s_max);
    o.hx(td.a_max);
    o.harr(&td.area_tbl);
    o.harr(&td.hrad_tbl);
    o.harr(&td.width_tbl);
    dump_lut(o, &td.area_lut);
}
fn ptr_code<T: ?Sized>(p: Option<&T>, mine: Option<*const u8>) -> i64 {
    match p {
        None => 0,
        Some(r) => {
            if Some(r as *const T as *const u8) == mine {
                1
            } else {
                2
            }
        }
    }
}
fn dump_xs(o: &mut Out, xs: &XSectParams, td: Option<&TransectData>) {
    o.hi(xs.type_ as i64);
    o.hi(xs.culvert_code as i64);
    o.hi(xs.transect as i64);
    for v in [
        xs.y_full, xs.w_max, xs.yw_max, xs.a_full, xs.r_full, xs.s_full, xs.s_max, xs.y_bot, xs.a_bot,
        xs.s_bot, xs.r_bot,
    ] {
        o.hx(v);
    }
    o.hi(xs.transect_tbl_size as i64);
    let tp = |a: &[f64]| a.as_ptr() as *const u8;
    let c = ptr_code(xs.area_tbl.map(|s| &s[0]), td.map(|t| tp(&t.area_tbl)));
    o.hi(c);
    let c = ptr_code(xs.hrad_tbl.map(|s| &s[0]), td.map(|t| tp(&t.hrad_tbl)));
    o.hi(c);
    let c = ptr_code(xs.width_tbl.map(|s| &s[0]), td.map(|t| tp(&t.width_tbl)));
    o.hi(c);
    let c = ptr_code(xs.area_lut, td.map(|t| &t.area_lut as *const LocateLut as *const u8));
    o.hi(c);
}

// ---------------------------------------------------------------- sweeps
fn base_of(v: f64) -> f64 {
    if v.is_finite() && v > 0.0 {
        v
    } else {
        1.0
    }
}

struct Ctx {
    lists: HashMap<String, Vec<f64>>,
}
impl Ctx {
    fn l(&self, n: &str) -> &Vec<f64> {
        &self.lists[n]
    }
}

fn cullines(o: &mut Out, c: &Ctx, xs: &XSectParams, codes: &[i32], slopes: &[f64]) {
    let yb = base_of(xs.y_full);
    let qb = base_of(xs.a_full * (32.2 * yb).sqrt());
    for &code in codes {
        for &slope in slopes {
            for q0 in [1e30, 0.5 * qb, 0.0] {
                o.hs("CUL");
                o.hi(code as i64);
                o.hx(slope);
                o.hx(q0);
                for &fr in c.l("CD") {
                    let y = fr * yb;
                    let mut dqdh = 7.0;
                    let mut controls = true;
                    let q = culvert::get_inflow(q0, y, xs, slope, code, &mut dqdh, &mut controls);
                    o.hx(q);
                    o.hx(dqdh);
                    o.hi(if controls { 1 } else { 0 });
                }
                o.nl();
            }
        }
    }
}

fn link_sweep(o: &mut Out, c: &Ctx, xs: &XSectParams, level: i32) {
    let yb = base_of(xs.y_full);
    let d: Vec<f64> = c.l("LD").iter().map(|&fr| fr * yb).collect();
    for &flow in c.l("LQV") {
        for b in [1, 0, 3, -2] {
            o.hs("VEL");
            o.hx(flow);
            o.hi(b as i64);
            for &dd in &d {
                o.hx(link::get_velocity(xs, flow, dd, b));
            }
            o.nl();
        }
    }
    for &v in c.l("LV") {
        o.hs("FRO");
        o.hx(v);
        for &dd in &d {
            o.hx(link::get_froude(xs, v, dd));
        }
        o.nl();
    }
    o.hs("CAP");
    for &dd in &d {
        o.hx(link::get_capacity(xs, dd));
    }
    o.nl();
    let ln = c.l("LN");
    let lsl = c.l("LS");
    let nn = if level != 0 { ln.len() } else { 2 };
    let ns = if level != 0 { lsl.len() } else { 2 };
    for i in 0..nn {
        for j in 0..ns {
            let mut beta = -9.0;
            let mut rf = -9.0;
            let mut qf = -9.0;
            link::compute_conveyance(ln[i], lsl[j], xs.s_full, &mut beta, &mut rf, &mut qf);
            o.hs("CONV");
            o.hx(beta);
            o.hx(rf);
            o.hx(qf);
            o.nl();
            let qb = base_of(qf);
            for qmax in [-1.0, 0.0, 0.5 * qb, qb, f64::NAN] {
                o.hs("DFF");
                o.hx(qmax);
                for &fr in c.l("LQF") {
                    o.hx(link::get_depth_from_flow(xs, beta, fr * qb, qmax));
                }
                o.nl();
            }
        }
    }
    if level != 0 {
        let codes: Vec<i32> = (-1..=58).collect();
        let cs = c.l("CS").clone();
        cullines(o, c, xs, &codes, &cs);
    } else {
        cullines(o, c, xs, &[0, 1, 5, 12, 37, 50, 58], &[0.01]);
    }
}

fn tabulated_probe(o: &mut Out, td: &TransectData, idx: i32) {
    for t in [22, 25, 23, 1] {
        for yf in [td.y_full, 2.5] {
            let mut xs = XSectParams { type_: t, a_bot: 0.123, ..XSectParams::default() };
            link::apply_tabulated_xsect_params(&mut xs, td, yf, idx);
            o.hs("ATP");
            o.hi(t as i64);
            o.hx(yf);
            dump_xs(o, &xs, Some(td));
            o.nl();
        }
    }
}

fn sane_sweep(xs: &XSectParams) -> bool {
    xs.y_full > 0.0 && xs.a_full > 0.0 && xs.y_full.is_finite() && xs.a_full.is_finite()
}

// ---------------------------------------------------------------- root-finder test functions
struct Ev {
    n: i64,
    h: u64,
}
impl Ev {
    fn reset(&mut self) {
        self.n = 0;
        self.h = 1469598103934665603u64;
    }
    fn ev(&mut self, x: f64) {
        self.h = self.h.wrapping_mul(1099511628211u64) ^ x.to_bits();
        self.n += 1;
    }
}
fn tf(e: &mut Ev, fid: i32, x: f64, df: &mut f64) -> f64 {
    e.ev(x);
    match fid {
        0 => {
            *df = 3.0 * x * x - 2.0;
            x * x * x - 2.0 * x - 5.0
        }
        1 => {
            *df = -x.sin() - 1.0;
            x.cos() - x
        }
        2 => {
            *df = 2.0 * x;
            x * x - 2.0
        }
        3 => {
            *df = x.exp();
            x.exp() - 3.0
        }
        4 => {
            let s = (1.0 + x * x).sqrt();
            *df = 1.0 / (s * s * s);
            x / s
        }
        5 => {
            *df = 3.0 * (x - 1.0) * (x - 1.0);
            (x - 1.0) * (x - 1.0) * (x - 1.0)
        }
        6 => {
            *df = 0.0;
            0.0
        }
        7 => {
            *df = 0.0;
            if x < 0.5 {
                -1.0
            } else {
                1.0
            }
        }
        8 => {
            *df = 0.5 / x.sqrt();
            x.sqrt() - 1.0
        }
        9 => {
            *df = -1.0 / ((x - 0.3) * (x - 0.3));
            1.0 / (x - 0.3)
        }
        _ => {
            *df = 2.0 * x;
            x * x + 1.0
        }
    }
}

// ---------------------------------------------------------------- run
fn main() {
    let path = std::env::args().nth(1).expect("usage: hydfuncs FILE");
    let text = std::fs::read_to_string(&path).unwrap();
    let mut o = Out(String::new());
    let mut c = Ctx { lists: HashMap::new() };
    let mut ts = TransectStore::default();
    let mut pool: Vec<TransectData> = Vec::new();
    let mut shxs: Vec<XSectParams<'static>> = Vec::new();
    let mut ev = Ev { n: 0, h: 0 };

    for line in text.lines() {
        let mut it = line.split_whitespace();
        let tag = match it.next() {
            Some(t) => t,
            None => continue,
        };
        match tag {
            "L" => {
                let name = it.next().unwrap().to_string();
                let v = rdvec(&mut it);
                c.lists.insert(name, v);
            }
            "TS" => {
                ts = TransectStore::default();
            }
            "TR" => {
                ts.names.push(rdname(&mut it));
                ts.comments.push(Vec::new());
                ts.n_left.push(rdh(&mut it));
                ts.n_right.push(rdh(&mut it));
                ts.n_channel.push(rdh(&mut it));
                ts.x_left_bank.push(rdh(&mut it));
                ts.x_right_bank.push(rdh(&mut it));
                ts.x_left_encroachment.push(0.0);
                ts.x_right_encroachment.push(0.0);
                ts.x_factor.push(rdh(&mut it));
                ts.y_factor.push(rdh(&mut it));
                ts.length_factor.push(rdh(&mut it));
                ts.stations.push(rdvec(&mut it));
                ts.elevations.push(rdvec(&mut it));
            }
            "TE" => {
                let ucfs = rdvec(&mut it);
                for &u in &ucfs {
                    for i in 0..ts.count() {
                        let mut td = TransectData::default();
                        let ret = transect::build_from_store(&ts, i, u, &mut td);
                        o.hs("BFS");
                        o.hx(u);
                        o.hi(i as i64);
                        o.hi(if ret { 1 } else { 0 });
                        dump_td(&mut o, &td);
                        o.nl();
                        if ret {
                            let mut t2 = TransectData {
                                name: td.name.clone(),
                                n_left: td.n_left,
                                n_right: td.n_right,
                                n_channel: td.n_channel,
                                length_factor: td.length_factor,
                                stations: td.stations.clone(),
                                elevations: td.elevations.clone(),
                                x_left_bank: td.x_left_bank,
                                x_right_bank: td.x_right_bank,
                                ..TransectData::default()
                            };
                            transect::build_tables(&mut t2, false);
                            o.hs("BT0");
                            dump_td(&mut o, &t2);
                            o.nl();
                        }
                        let idx = pool.len() as i32;
                        pool.push(td);
                        let pt = pool.last().unwrap();
                        tabulated_probe(&mut o, pt, idx);
                        if ret
                            && pt.y_full > 0.0
                            && pt.a_full > 0.0
                            && pt.y_full.is_finite()
                            && pt.a_full.is_finite()
                        {
                            let mut xs = XSectParams { type_: 22, ..XSectParams::default() };
                            link::apply_tabulated_xsect_params(&mut xs, pt, pt.y_full, idx);
                            link_sweep(&mut o, &c, &xs, 0);
                        }
                    }
                }
            }
            "ST" => {
                let mut sp = StreetParams::default();
                sp.width = rdh(&mut it);
                sp.curb_height = rdh(&mut it);
                sp.slope = rdh(&mut it);
                sp.roughness = rdh(&mut it);
                sp.gutter_depression = rdh(&mut it);
                sp.gutter_width = rdh(&mut it);
                sp.sides = rdi(&mut it) as i32;
                sp.back_width = rdh(&mut it);
                sp.back_slope = rdh(&mut it);
                sp.back_roughness = rdh(&mut it);
                let mut td = TransectData::default();
                street::build_transect(&sp, &mut td);
                o.hs("STR");
                dump_td(&mut o, &td);
                o.nl();
                let idx = pool.len() as i32;
                pool.push(td);
                let pt = pool.last().unwrap();
                tabulated_probe(&mut o, pt, idx);
                if pt.y_full > 0.0 && pt.a_full > 0.0 && pt.y_full.is_finite() && pt.a_full.is_finite() {
                    let mut xs = XSectParams { type_: 25, ..XSectParams::default() };
                    link::apply_tabulated_xsect_params(&mut xs, pt, pt.y_full, idx);
                    link_sweep(&mut o, &c, &xs, 0);
                }
            }
            "CU" => {
                let name = rdname(&mut it);
                let yf = rdh(&mut it);
                let wm = rdh(&mut it);
                let n = rdi(&mut it) as i32;
                let x = rdvec(&mut it);
                let y = rdvec(&mut it);
                let mut td = TransectData {
                    name,
                    y_full: yf,
                    a_full: 0.0,
                    r_full: 0.0,
                    w_max: wm,
                    ..TransectData::default()
                };
                transect::build_custom_tables(&mut td, yf, &x, &y, n);
                o.hs("CUS");
                dump_td(&mut o, &td);
                o.nl();
                let idx = pool.len() as i32;
                pool.push(td);
                let pt = pool.last().unwrap();
                tabulated_probe(&mut o, pt, idx);
                let mut xs = XSectParams { type_: 23, ..XSectParams::default() };
                link::apply_tabulated_xsect_params(&mut xs, pt, yf, idx);
                if sane_sweep(&xs) && pt.area_lut.scale != 0.0 {
                    link_sweep(&mut o, &c, &xs, 0);
                }
            }
            "SH" => {
                let t = rdi(&mut it) as i32;
                let p = [rdh(&mut it), rdh(&mut it), rdh(&mut it), rdh(&mut it), 0.0];
                let u = rdh(&mut it);
                let full = rdi(&mut it) as i32;
                let mut xs = XSectParams::default();
                let rc = xsect::set_params(&mut xs, t, &p, u);
                o.hs("SHP");
                o.hi(t as i64);
                o.hi(rc as i64);
                dump_xs(&mut o, &xs, None);
                o.nl();
                shxs.push(xs);
                link_sweep(&mut o, &c, &xs, full);
            }
            "XS" => {
                use XsectShape::*;
                let all = [
                    Circular,
                    FilledCircular,
                    RectClosed,
                    RectOpen,
                    Trapezoidal,
                    Triangular,
                    Parabolic,
                    Power,
                    Modbaskethandle,
                    Eggshaped,
                    Horseshoe,
                    Gothic,
                    Catenary,
                    Semielliptical,
                    Baskethandle,
                    Semicircular,
                    RectTriang,
                    RectRound,
                    HorizEllipse,
                    VertEllipse,
                    Arch,
                    Irregular,
                    Custom,
                    ForceMain,
                    StreetXsect,
                    Dummy,
                ];
                o.hs("XSH");
                for (k, &s) in all.iter().enumerate() {
                    assert_eq!(s as i32, k as i32);
                    o.hi(link::translate_shape(s) as i64);
                }
                o.nl();
            }
            "BX" => {
                let n_p = pool.len() as i32;
                let n_s = shxs.len() as i32;
                let n = n_p + n_s + 8;
                let mut links = LinkData::default();
                links.resize(n);
                let shapes = [22, 23, 25, 1, 0, 24, 3];
                for i in 0..n {
                    let u = i as usize;
                    if i < n_s {
                        let s = &shxs[u];
                        links.xsect_batch_shape[u] = s.type_;
                        links.xsect_y_full[u] = s.y_full;
                        links.xsect_a_full[u] = s.a_full;
                        links.xsect_w_max[u] = s.w_max;
                        links.xsect_r_full[u] = s.r_full;
                        links.xsect_s_full[u] = s.s_full;
                        links.xsect_s_max[u] = s.s_max;
                        links.xsect_y_bot[u] = s.y_bot;
                        links.xsect_a_bot[u] = s.a_bot;
                        links.xsect_s_bot[u] = s.s_bot;
                        links.xsect_r_bot[u] = s.r_bot;
                        links.xsect_curve[u] = i - 2;
                    } else {
                        let fi = i as f64;
                        links.xsect_batch_shape[u] = shapes[(i % 7) as usize];
                        links.xsect_y_full[u] = fi * 0.25;
                        links.xsect_a_full[u] = fi * 0.5;
                        links.xsect_w_max[u] = fi * 0.125;
                        links.xsect_r_full[u] = fi * 0.0625;
                        links.xsect_s_full[u] = fi * 1.5;
                        links.xsect_s_max[u] = fi * 1.75;
                        links.xsect_y_bot[u] = fi * 0.03125;
                        links.xsect_a_bot[u] = fi * 2.0;
                        links.xsect_s_bot[u] = fi * 3.0;
                        links.xsect_r_bot[u] = fi * 4.0;
                        links.xsect_curve[u] = (i - n_s) - 3;
                    }
                    for with in [0, 1] {
                        let xs = link::build_xsect_params(&links, u, if with != 0 { Some(&pool[..]) } else { None });
                        let td = if xs.transect >= 0 && xs.transect < n_p {
                            Some(&pool[xs.transect as usize])
                        } else {
                            None
                        };
                        o.hs("BXP");
                        o.hi(i as i64);
                        o.hi(with as i64);
                        dump_xs(&mut o, &xs, td);
                        o.nl();
                    }
                }
            }
            "FB" => {
                use XsectShape::*;
                let shapes = [
                    Circular,
                    FilledCircular,
                    RectClosed,
                    RectOpen,
                    Trapezoidal,
                    Triangular,
                    Parabolic,
                    Power,
                    Modbaskethandle,
                    Eggshaped,
                    Horseshoe,
                    Gothic,
                    Catenary,
                    Semielliptical,
                    Baskethandle,
                    Semicircular,
                    RectTriang,
                    RectRound,
                    HorizEllipse,
                    VertEllipse,
                    Arch,
                    Irregular,
                    Custom,
                    ForceMain,
                    StreetXsect,
                    Dummy,
                ];
                let types = [LinkType::Conduit, LinkType::Pump, LinkType::Orifice, LinkType::Weir, LinkType::Outlet];
                let mut links = LinkData::default();
                links.resize(5 * 26);
                for i in 0..5 * 26usize {
                    links.type_[i] = types[i / 26];
                    links.xsect_shape[i] = shapes[i % 26];
                    links.xsect_y_bot[i] = i as f64 * 0.125 + 0.5;
                }
                o.hs("FBP");
                for i in 0..5 * 26usize {
                    for r in [false, true] {
                        o.hx(link::filled_circular_offset_bump(&links, i, r));
                    }
                }
                o.nl();
            }
            "FM" => {
                let fv = c.l("FV").clone();
                let frr = c.l("FR").clone();
                let fp = c.l("FP").clone();
                for &v in &fv {
                    for &r in &frr {
                        o.hs("FHW");
                        for &cc in &fp {
                            o.hx(forcemain::get_fric_slope_hw(v, r, cc));
                            o.hx(hydkernels::fric_slope_hw(v, r, cc));
                        }
                        o.nl();
                        o.hs("FDW");
                        for &cc in &fp {
                            o.hx(forcemain::get_fric_slope_dw(v, r, cc));
                            o.hx(hydkernels::fric_slope_dw(v, r, cc));
                        }
                        o.nl();
                    }
                }
                for &e in c.l("FE") {
                    for &h in c.l("FH") {
                        o.hs("FFF");
                        for &re in c.l("FRE") {
                            o.hx(forcemain::get_fric_factor(e, h, re));
                        }
                        o.nl();
                    }
                }
                for m in [FrictionModel::HazenWilliams, FrictionModel::DarcyWeisbach] {
                    for &rb in c.l("FB") {
                        for &yf in c.l("FY") {
                            o.hs("FEN");
                            for &s in c.l("FS") {
                                o.hx(forcemain::get_equiv_n(m, rb, yf, s, 0.013));
                            }
                            o.nl();
                        }
                        o.hs("FRF");
                        for &lf in c.l("FL") {
                            o.hx(forcemain::get_rough_factor(m, rb, lf));
                        }
                        o.nl();
                    }
                }
                let nb = fv.len().min(frr.len()).min(fp.len());
                for m in [FrictionModel::HazenWilliams, FrictionModel::DarcyWeisbach] {
                    for cnt in [nb as i32, 0, -1, 3] {
                        let mut out = vec![-7.0; nb];
                        forcemain::batch_fric_slope(&fv, &frr, &fp, &mut out, m, cnt);
                        o.hs("FBS");
                        o.harr(&out);
                        o.nl();
                    }
                }
            }
            "KC" => {
                let mut curves: Vec<CulvertCurve> = Vec::new();
                for code in 0..=57 {
                    let k = culvert::get_coeffs(code);
                    curves.push(CulvertCurve { K: k.K, M: k.M, C: k.C, Y: k.Y });
                }
                let cv = |k: f64, m: f64, cc: f64, y: f64| CulvertCurve { K: k, M: m, C: cc, Y: y };
                curves.push(cv(0.01, 0.0, 0.03, 0.7));
                curves.push(cv(0.01, 2.0, 0.0, 0.7));
                curves.push(cv(0.01, 2.0, -0.03, 0.7));
                curves.push(cv(-0.01, 2.0, 0.03, 0.7));
                curves.push(cv(f64::NAN, 2.0, 0.03, 0.7));
                curves.push(cv(0.5, 0.667, 0.03, 0.3));
                for cc in &curves {
                    for mit in [false, true] {
                        for &s in c.l("KS") {
                            for yf in [3.0, 0.0, f64::NAN] {
                                for q in [1e30, 5.0] {
                                    o.hs("KCI");
                                    for &h in c.l("KH") {
                                        let mut dq = 7.0;
                                        o.hx(hydkernels::culvert_inflow(q, h, yf, 7.0685835, s, cc, mit, &mut dq));
                                        o.hx(dq);
                                    }
                                    o.nl();
                                }
                            }
                            o.hs("KCA");
                            for af in [0.0, -1.0, f64::NAN, 1e-12] {
                                let mut dq = 7.0;
                                o.hx(hydkernels::culvert_inflow(1e30, 1.5, 3.0, af, s, cc, mit, &mut dq));
                                o.hx(dq);
                            }
                            o.nl();
                        }
                    }
                }
            }
            "GC" => {
                o.hs("GCO");
                for code in -2..=60 {
                    let k = culvert::get_coeffs(code);
                    for v in [k.form, k.K, k.M, k.C, k.Y] {
                        o.hx(v);
                    }
                }
                o.nl();
            }
            "XP" => {
                for t in [3, 1] {
                    for yf in [3.0, 0.0, -3.0, f64::NAN] {
                        for af in [7.0, 0.0, -1.0, f64::NAN] {
                            let mut xs = XSectParams::default();
                            let p = [3.0, 2.0, 0.0, 0.0, 0.0];
                            xsect::set_params(&mut xs, t, &p, 1.0);
                            xs.y_full = yf;
                            xs.a_full = af;
                            o.hs("XPP");
                            o.hi(t as i64);
                            o.hx(yf);
                            o.hx(af);
                            for &fr in c.l("LD") {
                                let d = fr * 3.0;
                                o.hx(link::get_capacity(&xs, d));
                                o.hx(link::get_velocity(&xs, 2.0, d, 1));
                                o.hx(link::get_froude(&xs, 1.5, d));
                            }
                            o.nl();
                            for code in [1, 12, 30, 37] {
                                o.hs("XPC");
                                o.hi(code as i64);
                                for &fr in c.l("CD") {
                                    let mut dq = 7.0;
                                    let mut ctl = true;
                                    o.hx(culvert::get_inflow(1e30, fr * 3.0, &xs, 0.01, code, &mut dq, &mut ctl));
                                    o.hx(dq);
                                    o.hi(if ctl { 1 } else { 0 });
                                }
                                o.nl();
                            }
                        }
                    }
                }
            }
            "HP" => {
                let h = c.l("FV").clone();
                for &q in c.l("LQF") {
                    o.hs("HPW");
                    for &h1 in &h {
                        for &h2 in &h {
                            o.hx(link::get_hyd_power(q, h1, h2));
                        }
                    }
                    o.nl();
                }
            }
            "FR" => {
                let fid = rdi(&mut it) as i32;
                let x1 = rdh(&mut it);
                let x2 = rdh(&mut it);
                let rts0 = rdh(&mut it);
                let xacc = rdh(&mut it);
                ev.reset();
                let mut rts = rts0;
                let n = findroot::newton(x1, x2, &mut rts, xacc, |x, f, df| {
                    *f = tf(&mut ev, fid, x, df);
                });
                o.hs("FRN");
                o.hi(n as i64);
                o.hx(rts);
                o.hi(ev.n);
                o.hi(ev.h as i64);
                ev.reset();
                let r = findroot::ridder(x1, x2, xacc, |x| {
                    let mut d = 0.0;
                    tf(&mut ev, fid, x, &mut d)
                });
                o.hs("FRR");
                o.hx(r);
                o.hi(ev.n);
                o.hi(ev.h as i64);
                o.nl();
            }
            _ => panic!("unknown command {}", tag),
        }
    }
    let _ = std::io::stdout().lock().write_all(o.0.as_bytes());
}
