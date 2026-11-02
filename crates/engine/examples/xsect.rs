// Rust side of gate/diff/xsect.cpp: reads the file `xsect gen` wrote, calls every
// translated function of hydraulics::{xsect_tables, xsect_kernels, xsection, xsect_batch}
// with the same arguments, and prints the same lines. Calls on which upstream is undefined
// or does not return print "-" on both sides (same predicates as the harness).
#![allow(clippy::too_many_arguments, clippy::needless_range_loop)]
use std::fmt::Write as _;
use std::io::Write;

use swmm6_rust::core::context::SimulationContext;
use swmm6_rust::hydraulics::transect_types::{TransectData, N_TRANSECT_TBL};
use swmm6_rust::hydraulics::xsect_batch::{self as xb, XSectGroups, XSectParams, XSectShape};
use swmm6_rust::hydraulics::xsect_kernels::{
    LutId, GRAVITY, PI, RECT_ALFMAX, RECT_ROUND_ALFMAX, RECT_TRIANG_ALFMAX, TINY,
};
use swmm6_rust::hydraulics::xsect_lookup::{build_invlookup_lut, LocateLut};
use swmm6_rust::hydraulics::xsect_tables as xt;
use swmm6_rust::hydraulics::xsection as xs_;

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
    fn dash(&mut self) {
        self.0.push_str("- ");
    }
    fn nl(&mut self) {
        self.0.push('\n');
    }
    fn harr(&mut self, a: &[f64]) {
        for &v in a {
            self.hx(v);
        }
    }
    fn prline(&mut self, tag: &str, vi: i32, v: &[f64]) {
        self.hs(tag);
        self.hi(vi as i64);
        for &x in v {
            self.hx(x);
        }
        self.nl();
    }
}

// ---------------------------------------------------------------- input
struct TRec {
    yf: f64,
    wm: f64,
    af: f64,
    rf: f64,
    sf: f64,
    sm: f64,
    ab: f64,
    ar: Vec<f64>,
    hr: Vec<f64>,
    wt: Vec<f64>,
}
struct PRec {
    type_: i32,
    p: [f64; 4],
    ucf: f64,
    tid: i32,
}

fn rdh(it: &mut std::str::SplitWhitespace) -> f64 {
    f64::from_bits(u64::from_str_radix(it.next().unwrap(), 16).unwrap())
}

// ---------------------------------------------------------------- predicates
fn is_tab_s(t: i32) -> bool {
    (15..=21).contains(&t)
}
fn is_rec_type(t: i32) -> bool {
    t == 0 || !(0..=25).contains(&t)
}
fn is_tabulated(t: i32) -> bool {
    t == 22 || t == 23 || t == 25
}
fn pofy_ok(xs: &XSectParams, y: f64) -> bool {
    if !(y > 0.0) {
        return true;
    }
    let dy1 = 0.02 * xs.y_full;
    if dy1.is_nan() {
        return true;
    }
    if !(dy1 > 0.0) {
        return false;
    }
    if dy1.is_infinite() {
        return true;
    }
    y.is_finite() && y / dy1 <= 1.0e4
}
fn pf_yofa(xs: &XSectParams, a: f64) -> f64 {
    (a / xs.r_bot).powf(1.0 / (xs.s_bot + 1.0))
}
fn pf_sane(xs: &XSectParams) -> bool {
    xs.y_full > 0.0
        && xs.y_full.is_finite()
        && xs.a_full > 0.0
        && xs.a_full.is_finite()
        && xs.r_bot > 0.0
        && xs.r_bot.is_finite()
        && xs.s_bot > 0.0
        && xs.s_bot.is_finite()
}
fn ok_rofy(xs: &XSectParams, y: f64) -> bool {
    if xs.type_ == 8 {
        return pofy_ok(xs, y);
    }
    true
}
fn ok_rofa(xs: &XSectParams, a: f64) -> bool {
    if is_rec_type(xs.type_) {
        return a <= 0.0;
    }
    if xs.type_ == 8 {
        return a <= 0.0 || pofy_ok(xs, pf_yofa(xs, a));
    }
    true
}
fn ok_sofa(xs: &XSectParams, a: f64) -> bool {
    if is_rec_type(xs.type_) {
        return a <= 0.0;
    }
    if xs.type_ == 8 {
        return a == 0.0 || ok_rofa(xs, a);
    }
    true
}
fn tab_ok(alpha: f64, n: i32) -> bool {
    let v = alpha / (1.0 / ((n - 1) as f64));
    v > -1.0 && v < 2147483648.0
}
fn ok_dsda(xs: &XSectParams, a: f64) -> bool {
    let t = xs.type_;
    if t == 1 || t == 24 {
        let alpha = a / xs.a_full;
        if alpha <= 1.0e-30 || alpha < 0.04 {
            return true;
        }
        return tab_ok(alpha, 51);
    }
    if is_tab_s(t) {
        return tab_ok(a / xs.a_full, 51);
    }
    if t == 3 || t == 4 || t == 9 || t == 10 || t == 11 || t == 5 || t == 6 {
        return true;
    }
    let alpha = a / xs.a_full;
    let mut alpha1 = alpha - 0.001;
    let alpha2 = alpha + 0.001;
    if alpha1 < 0.0 {
        alpha1 = 0.0;
    }
    let a1 = alpha1 * xs.a_full;
    let a2 = alpha2 * xs.a_full;
    ok_sofa(xs, a2) && ok_sofa(xs, a1)
}
fn ok_aofs(xs: &XSectParams, s: f64) -> bool {
    if s <= 0.0 {
        return true;
    }
    let t = xs.type_;
    if t == 0 {
        return true;
    }
    if !(0..=25).contains(&t) {
        return false;
    }
    if t == 8 {
        return pf_sane(xs);
    }
    true
}
fn ok_ycrit(xs: &XSectParams, q: f64) -> bool {
    if q <= 0.0 {
        return true;
    }
    let q2g = q * q / GRAVITY;
    if q2g == 0.0 {
        return true;
    }
    let t = xs.type_;
    if t == 0 || t == 3 || t == 4 || t == 6 || t == 7 || t == 8 {
        return true;
    }
    let mut y0 = 1.01 * (q2g / xs.y_full).powf(0.25);
    if y0 >= xs.y_full {
        y0 = 0.97 * xs.y_full;
    }
    let r = xs.a_full / (PI / 4.0 * xs.y_full * xs.y_full);
    if r >= 0.5 && r <= 2.0 {
        let dy = xs.y_full / 25.0;
        let v = y0 / dy;
        return v >= -100.0 && v <= 100.0;
    }
    true
}
fn batch_ok(xs: &XSectParams) -> bool {
    xs.type_ != 8 || pf_sane(xs)
}
fn base_of(v: f64) -> f64 {
    if v.is_finite() && v > 0.0 {
        v
    } else {
        1.0
    }
}

const LUT_IDS: [LutId; 15] = [
    LutId::YGothic,
    LutId::YCatenary,
    LutId::YSemiEllip,
    LutId::YSemiCirc,
    LutId::AHorizEllipse,
    LutId::AVertEllipse,
    LutId::AArch,
    LutId::SCirc,
    LutId::SEgg,
    LutId::SHorseshoe,
    LutId::SGothic,
    LutId::SCatenary,
    LutId::SSemiEllip,
    LutId::SBasketHandle,
    LutId::SSemiCirc,
];

fn table_list() -> Vec<(&'static str, &'static [f64], i32, i32)> {
    let t = xs_::host_tables();
    vec![
        ("A_Arch", t.A_Arch, t.N_A_Arch, 6),
        ("A_Baskethandle", t.A_Baskethandle, t.N_A_Baskethandle, -1),
        ("A_Circ", t.A_Circ, t.N_A_Circ, -1),
        ("A_Egg", t.A_Egg, t.N_A_Egg, -1),
        ("A_HorizEllipse", t.A_HorizEllipse, t.N_A_HorizEllipse, 4),
        ("A_Horseshoe", t.A_Horseshoe, t.N_A_Horseshoe, -1),
        ("A_VertEllipse", t.A_VertEllipse, t.N_A_VertEllipse, 5),
        ("R_Arch", t.R_Arch, t.N_R_Arch, -1),
        ("R_Baskethandle", t.R_Baskethandle, t.N_R_Baskethandle, -1),
        ("R_Circ", t.R_Circ, t.N_R_Circ, -1),
        ("R_Egg", t.R_Egg, t.N_R_Egg, -1),
        ("R_HorizEllipse", t.R_HorizEllipse, t.N_R_HorizEllipse, -1),
        ("R_Horseshoe", t.R_Horseshoe, t.N_R_Horseshoe, -1),
        ("R_VertEllipse", t.R_VertEllipse, t.N_R_VertEllipse, -1),
        ("S_BasketHandle", t.S_BasketHandle, t.N_S_BasketHandle, 13),
        ("S_Catenary", t.S_Catenary, t.N_S_Catenary, 11),
        ("S_Circ", t.S_Circ, t.N_S_Circ, 7),
        ("S_Egg", t.S_Egg, t.N_S_Egg, 8),
        ("S_Gothic", t.S_Gothic, t.N_S_Gothic, 10),
        ("S_Horseshoe", t.S_Horseshoe, t.N_S_Horseshoe, 9),
        ("S_SemiCirc", t.S_SemiCirc, t.N_S_SemiCirc, 14),
        ("S_SemiEllip", t.S_SemiEllip, t.N_S_SemiEllip, 12),
        ("W_Arch", t.W_Arch, t.N_W_Arch, -1),
        ("W_BasketHandle", t.W_BasketHandle, t.N_W_BasketHandle, -1),
        ("W_Catenary", t.W_Catenary, t.N_W_Catenary, -1),
        ("W_Circ", t.W_Circ, t.N_W_Circ, -1),
        ("W_Egg", t.W_Egg, t.N_W_Egg, -1),
        ("W_Gothic", t.W_Gothic, t.N_W_Gothic, -1),
        ("W_HorizEllipse", t.W_HorizEllipse, t.N_W_HorizEllipse, -1),
        ("W_Horseshoe", t.W_Horseshoe, t.N_W_Horseshoe, -1),
        ("W_SemiCirc", t.W_SemiCirc, t.N_W_SemiCirc, -1),
        ("W_SemiEllip", t.W_SemiEllip, t.N_W_SemiEllip, -1),
        ("W_VertEllipse", t.W_VertEllipse, t.N_W_VertEllipse, -1),
        ("Y_BasketHandle", t.Y_BasketHandle, t.N_Y_BasketHandle, -1),
        ("Y_Catenary", t.Y_Catenary, t.N_Y_Catenary, 1),
        ("Y_Circ", t.Y_Circ, t.N_Y_Circ, -1),
        ("Y_Egg", t.Y_Egg, t.N_Y_Egg, -1),
        ("Y_Gothic", t.Y_Gothic, t.N_Y_Gothic, 0),
        ("Y_Horseshoe", t.Y_Horseshoe, t.N_Y_Horseshoe, -1),
        ("Y_SemiCirc", t.Y_SemiCirc, t.N_Y_SemiCirc, 3),
        ("Y_SemiEllip", t.Y_SemiEllip, t.N_Y_SemiEllip, 2),
    ]
}

fn dump_lut(o: &mut Out, l: &LocateLut) {
    o.hx(l.t0);
    o.hx(l.scale);
    o.hi(l.j_last as i64);
    for b in 0..=LocateLut::K_BUCKETS as usize {
        o.hi(l.lo[b] as i64);
    }
}

fn dump_xs(o: &mut Out, xs: &XSectParams) {
    o.hi(xs.type_ as i64);
    o.hi(xs.culvert_code as i64);
    o.hi(xs.transect as i64);
    for v in [
        xs.y_full, xs.w_max, xs.yw_max, xs.a_full, xs.r_full, xs.s_full, xs.s_max, xs.y_bot, xs.a_bot,
        xs.s_bot, xs.r_bot,
    ] {
        o.hx(v);
    }
}

fn dump_groups(o: &mut Out, tag: &str, g_: &XSectGroups) {
    o.hs("G");
    o.hs(tag);
    o.hi(g_.num_groups() as i64);
    o.nl();
    for gi in 0..g_.num_groups() {
        let g = g_.group(gi);
        o.hs("g");
        o.hi(gi as i64);
        o.hi(g.shape as i64);
        o.hi(g.count as i64);
        o.hi(g.link_idx.len() as i64);
        for &v in &g.link_idx {
            o.hi(v as i64);
        }
        for col in [
            &g.y_full,
            &g.a_full,
            &g.r_full,
            &g.s_full,
            &g.w_max,
            &g.inv_y_full,
            &g.y_bot,
            &g.a_bot,
            &g.s_bot,
            &g.r_bot,
        ] {
            o.hi(col.len() as i64);
            for &v in col.iter() {
                o.hx(v);
            }
        }
        for tv in [&g.area_tables, &g.hrad_tables, &g.width_tables] {
            o.hi(tv.len() as i64);
            for p in tv.iter() {
                match p {
                    None => o.hs("n"),
                    Some(t) => {
                        o.hx(t[0]);
                        o.hx(t[17]);
                        o.hx(t[50]);
                    }
                }
            }
        }
        o.hi(g.transect_tbl_size as i64);
        o.hi(g.buf_d.borrow().len() as i64);
        o.hi(g.buf_r.borrow().len() as i64);
        o.hi(g.buf_r2.borrow().len() as i64);
        o.nl();
    }
    o.hs("find");
    for s in 0..=25 {
        match g_.find_group(XSectShape::from_code(s).unwrap()) {
            Some(p) => o.hi(p.count as i64),
            None => o.hi(-1),
        }
    }
    o.nl();
}

const SENT: f64 = -12345.678;

fn main() {
    let path = std::env::args().nth(1).expect("usage: xsect FILE");
    let text = std::fs::read_to_string(&path).unwrap();
    let mut fv: Vec<f64> = Vec::new();
    let mut av: Vec<f64> = Vec::new();
    let mut qv: Vec<f64> = Vec::new();
    let mut tv: Vec<TRec> = Vec::new();
    let mut pv: Vec<PRec> = Vec::new();
    let nn = N_TRANSECT_TBL as usize;
    for line in text.lines() {
        let mut it = line.split_whitespace();
        let tag = match it.next() {
            Some(t) => t,
            None => continue,
        };
        match tag {
            "F" | "A" | "Q" => {
                let n: usize = it.next().unwrap().parse().unwrap();
                let v = match tag {
                    "F" => &mut fv,
                    "A" => &mut av,
                    _ => &mut qv,
                };
                for _ in 0..n {
                    v.push(rdh(&mut it));
                }
            }
            "T" => {
                let _id: i32 = it.next().unwrap().parse().unwrap();
                let yf = rdh(&mut it);
                let wm = rdh(&mut it);
                let af = rdh(&mut it);
                let rf = rdh(&mut it);
                let sf = rdh(&mut it);
                let sm = rdh(&mut it);
                let ab = rdh(&mut it);
                let ar: Vec<f64> = (0..nn).map(|_| rdh(&mut it)).collect();
                let hr: Vec<f64> = (0..nn).map(|_| rdh(&mut it)).collect();
                let wt: Vec<f64> = (0..nn).map(|_| rdh(&mut it)).collect();
                tv.push(TRec {
                    yf,
                    wm,
                    af,
                    rf,
                    sf,
                    sm,
                    ab,
                    ar,
                    hr,
                    wt,
                });
            }
            "P" => {
                let type_: i32 = it.next().unwrap().parse().unwrap();
                let p = [rdh(&mut it), rdh(&mut it), rdh(&mut it), rdh(&mut it)];
                let ucf = rdh(&mut it);
                let tid: i32 = it.next().unwrap().parse().unwrap();
                pv.push(PRec { type_, p, ucf, tid });
            }
            _ => {}
        }
    }
    let n_t = tv.len() as i32;
    let mut o = Out(String::with_capacity(64 << 20));

    // ---- tables and bucket maps
    let tl = table_list();
    for (name, t, n, _) in &tl {
        o.hs("TBL");
        o.hs(name);
        o.hi(*n as i64);
        o.harr(&t[..*n as usize]);
        o.nl();
    }
    o.hs("TBL Amax");
    o.harr(&xt::Amax);
    o.nl();
    o.hs("TBL Ellipse");
    o.hi(xt::NumCodesEllipse as i64);
    o.harr(&xt::MinorAxis_Ellipse);
    o.harr(&xt::MajorAxis_Ellipse);
    o.harr(&xt::Afull_Ellipse);
    o.harr(&xt::Rfull_Ellipse);
    o.nl();
    o.hs("TBL Arch");
    o.hi(xt::NumCodesArch as i64);
    o.harr(&xt::Yfull_Arch);
    o.harr(&xt::Wmax_Arch);
    o.harr(&xt::Afull_Arch);
    o.harr(&xt::Rfull_Arch);
    o.nl();
    {
        let t = xs_::host_tables();
        o.hs("TBL hostAmax");
        o.harr(&t.Amax[..26]);
        o.nl();
        for (id, lid) in LUT_IDS.iter().enumerate() {
            o.hs("LUT");
            o.hi(id as i64);
            dump_lut(&mut o, t.lut(*lid).unwrap());
            o.nl();
        }
    }
    o.hs("CONST");
    for v in [
        TINY,
        PI,
        GRAVITY,
        RECT_ALFMAX,
        RECT_TRIANG_ALFMAX,
        RECT_ROUND_ALFMAX,
    ] {
        o.hx(v);
    }
    o.nl();

    // ---- lookup helpers on every table
    let mut fa = fv.clone();
    fa.extend_from_slice(&av);
    for (name, t, n, lut) in &tl {
        o.hs("L");
        o.hs(name);
        for &x in &fa {
            o.hx(xs_::lookup(x, t, *n));
        }
        o.nl();
        o.hs("I");
        o.hs(name);
        for &y in &fa {
            o.hx(xs_::inv_lookup(y, t, *n, None));
        }
        o.nl();
        if *lut >= 0 {
            let l = xs_::host_tables().lut(LUT_IDS[*lut as usize]);
            o.hs("IL");
            o.hs(name);
            for &y in &fa {
                o.hx(xs_::inv_lookup(y, t, *n, l));
            }
            o.nl();
        }
        o.hs("LOC");
        o.hs(name);
        for &y in &fa {
            o.hi(xs_::locate(y, t, *n - 1) as i64);
        }
        o.nl();
    }
    {
        let e = xs_::host_eval();
        o.hs("YC");
        for &x in &fa {
            o.hx(xs_::get_ycircular(x));
        }
        o.nl();
        o.hs("SC");
        for &x in &fa {
            o.hx(xs_::get_scircular(x));
        }
        o.nl();
        o.hs("AC");
        for &x in &fa {
            o.hx(e.get_acircular(x));
        }
        o.nl();
        o.hs("TA");
        for &x in &fa {
            o.hx(e.get_theta_of_alpha(x));
        }
        o.nl();
        o.hs("TP");
        for &x in &fa {
            o.hx(e.get_theta_of_psi(x));
        }
        o.nl();
    }

    // ---- transect LUTs
    let mut tluts = vec![LocateLut::default(); tv.len()];
    for t in 0..tv.len() {
        build_invlookup_lut(&mut tluts[t], Some(&tv[t].ar), N_TRANSECT_TBL);
        o.hs("TLUT");
        o.hi(t as i64);
        dump_lut(&mut o, &tluts[t]);
        o.nl();
    }

    // ---- per-element cases
    let mut net: Vec<XSectParams> = Vec::new();
    let mut curve: Vec<i32> = Vec::new();
    for (c, pr) in pv.iter().enumerate() {
        let mut xs = XSectParams::default();
        let ret = xs_::set_params(&mut xs, pr.type_, &pr.p, pr.ucf);
        if is_tabulated(pr.type_) && pr.tid >= 0 && pr.tid < n_t {
            let t = &tv[pr.tid as usize];
            xs.y_full = t.yf;
            xs.w_max = t.wm;
            xs.a_full = t.af;
            xs.r_full = t.rf;
            xs.s_full = t.sf;
            xs.s_max = t.sm;
            xs.a_bot = t.ab;
            xs.area_tbl = Some(&t.ar);
            xs.hrad_tbl = Some(&t.hr);
            xs.width_tbl = Some(&t.wt);
            xs.transect_tbl_size = N_TRANSECT_TBL;
            xs.area_lut = Some(&tluts[pr.tid as usize]);
        }
        o.hs("C");
        o.hi(c as i64);
        o.hi(ret as i64);
        dump_xs(&mut o, &xs);
        o.nl();

        let yb = base_of(xs.y_full);
        let ab = base_of(xs.a_full);
        let sb = base_of(xs.s_full);
        let qb0 = xs.a_full * (GRAVITY * xs.y_full).sqrt();
        let qb = base_of(qb0);
        let mut d: Vec<f64> = fv.iter().map(|f| f * yb).collect();
        let mut aa: Vec<f64> = fv.iter().map(|f| f * ab).collect();
        let mut ss: Vec<f64> = fv.iter().map(|f| f * sb).collect();
        let mut qq: Vec<f64> = fv.iter().map(|f| f * qb).collect();
        d.extend_from_slice(&av);
        aa.extend_from_slice(&av);
        ss.extend_from_slice(&av);
        qq.extend_from_slice(&qv);

        o.hs("AofY");
        for &y in &d {
            o.hx(xs_::get_a_of_y(&xs, y));
        }
        o.nl();
        o.hs("WofY");
        for &y in &d {
            o.hx(xs_::get_w_of_y(&xs, y));
        }
        o.nl();
        o.hs("RofY");
        for &y in &d {
            if ok_rofy(&xs, y) {
                o.hx(xs_::get_r_of_y(&xs, y))
            } else {
                o.dash()
            }
        }
        o.nl();
        o.hs("YofA");
        for &a in &aa {
            o.hx(xs_::get_y_of_a(&xs, a));
        }
        o.nl();
        o.hs("SofA");
        for &a in &aa {
            if ok_sofa(&xs, a) {
                o.hx(xs_::get_s_of_a(&xs, a))
            } else {
                o.dash()
            }
        }
        o.nl();
        o.hs("RofA");
        for &a in &aa {
            if ok_rofa(&xs, a) {
                o.hx(xs_::get_r_of_a(&xs, a))
            } else {
                o.dash()
            }
        }
        o.nl();
        o.hs("dSdA");
        for &a in &aa {
            if ok_dsda(&xs, a) {
                o.hx(xs_::get_ds_da(&xs, a))
            } else {
                o.dash()
            }
        }
        o.nl();
        o.hs("AofS");
        for &s in &ss {
            if ok_aofs(&xs, s) {
                o.hx(xs_::get_a_of_s(&xs, s))
            } else {
                o.dash()
            }
        }
        o.nl();
        o.hs("Ycrit");
        for &q in &qq {
            if ok_ycrit(&xs, q) {
                o.hx(xs_::get_ycrit(&xs, q))
            } else {
                o.dash()
            }
        }
        o.nl();
        o.hs("Amax");
        o.hx(xs_::get_amax(&xs));
        o.hi(if xs_::is_open(pr.type_) { 1 } else { 0 });
        o.nl();

        if batch_ok(&xs) {
            // The network keeps only the scalar fields (build ignores table pointers).
            net.push(XSectParams {
                area_tbl: None,
                hrad_tbl: None,
                width_tbl: None,
                area_lut: None,
                ..xs
            });
            curve.push(if is_tabulated(pr.type_) { pr.tid } else { -1 });
        }
    }
    for t in [
        -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
        100,
    ] {
        o.hs("open");
        o.hi(t as i64);
        o.hi(if xs_::is_open(t) { 1 } else { 0 });
        o.nl();
    }

    // ---- batch: network of the batch-safe cases
    let n = net.len() as i32;
    let un = net.len();
    let mut ctx = SimulationContext::default();
    ctx.links.resize(n);
    for i in 0..un {
        ctx.links.xsect_curve[i] = curve[i];
    }
    for t in &tv {
        let mut td = TransectData::default();
        td.y_full = t.yf;
        td.w_max = t.wm;
        td.a_full = t.af;
        td.r_full = t.rf;
        for j in 0..nn {
            td.area_tbl[j] = t.ar[j];
            td.hrad_tbl[j] = t.hr[j];
            td.width_tbl[j] = t.wt[j];
        }
        ctx.transect_tables.push(td);
    }
    let by: Vec<f64> = net.iter().map(|x| base_of(x.y_full)).collect();
    let ba: Vec<f64> = net.iter().map(|x| base_of(x.a_full)).collect();
    let mut fb: Vec<f64> = Vec::new();
    for (i, &f) in fv.iter().enumerate() {
        let keep = f.is_nan() || (f.is_finite() && f.abs() <= 100.0);
        if keep && (i % 3 == 0 || i < 21) {
            fb.push(f);
        }
    }
    let nv = fb.len() as i32;
    let dvec = |vi: i32| -> Vec<f64> { (0..un).map(|i| fb[(vi % nv) as usize] * by[i]).collect() };
    let avec = |vi: i32| -> Vec<f64> { (0..un).map(|i| fb[(vi % nv) as usize] * ba[i]).collect() };
    let mut masks: Vec<Vec<u8>> = vec![vec![0u8; un]; 4];
    for i in 0..un {
        masks[0][i] = (i % 3 == 0) as u8;
        masks[1][i] = 1;
        masks[2][i] = 0;
        masks[3][i] = ((i * 7) % 5 < 2) as u8;
    }

    o.hs("NET");
    o.hi(n as i64);
    o.hi(nv as i64);
    o.nl();
    let mut g = XSectGroups::default();
    g.build(&net, n);
    dump_groups(&mut o, "build", &g);
    {
        let d = dvec(5);
        let mut o1 = vec![SENT; un];
        let mut o2 = vec![SENT; un];
        let mut o3 = vec![SENT; un];
        g.compute_areas(&d, &mut o1, n);
        g.compute_hyd_rad(&d, &mut o2, n);
        g.compute_widths(&d, &mut o3, n);
        o.prline("preA", 5, &o1);
        o.prline("preR", 5, &o2);
        o.prline("preW", 5, &o3);
    }
    g.attach_transect_tables(&ctx);
    dump_groups(&mut o, "attach", &g);
    g.attach_transect_tables(&ctx);
    dump_groups(&mut o, "attach2", &g);

    for vi in 0..nv {
        let d = dvec(vi);
        let d2 = dvec(vi + 1);
        let dm = dvec(vi + 2);
        let a = avec(vi);
        let mut o1 = vec![SENT; un];
        let mut o2 = vec![SENT; un];
        g.compute_areas(&d, &mut o1, n);
        o.prline("bA", vi, &o1);
        o1.fill(SENT);
        g.compute_hyd_rad(&d, &mut o1, n);
        o.prline("bR", vi, &o1);
        o1.fill(SENT);
        g.compute_widths(&d, &mut o1, n);
        o.prline("bW", vi, &o1);
        o1.fill(SENT);
        g.compute_area_and_hyd_rad(&d, &mut o1, &mut o2, n);
        o.prline("bAR_a", vi, &o1);
        o.prline("bAR_r", vi, &o2);
        o1.fill(SENT);
        g.compute_section_factors(&a, &mut o1, n);
        o.prline("bSF", vi, &o1);
        o1.fill(SENT);
        g.compute_depths_from_area(&a, &mut o1, n);
        o.prline("bDA", vi, &o1);

        let triple = |o: &mut Out, g: &XSectGroups, tag: &str| {
            let mut a1 = vec![SENT; un];
            let mut a2 = vec![SENT; un];
            let mut am = vec![SENT; un];
            let mut h1 = vec![SENT; un];
            let mut hm = vec![SENT; un];
            g.compute_area_hyd_rad_triple(&d, &d2, &dm, &mut a1, &mut a2, &mut am, &mut h1, &mut hm, n);
            let mut w1 = vec![SENT; un];
            let mut w2 = vec![SENT; un];
            let mut wm = vec![SENT; un];
            g.compute_widths_triple(&d, &d2, &dm, &mut w1, &mut w2, &mut wm, n);
            o.prline(&format!("{}_a1", tag), vi, &a1);
            o.prline(&format!("{}_a2", tag), vi, &a2);
            o.prline(&format!("{}_am", tag), vi, &am);
            o.prline(&format!("{}_h1", tag), vi, &h1);
            o.prline(&format!("{}_hm", tag), vi, &hm);
            o.prline(&format!("{}_w1", tag), vi, &w1);
            o.prline(&format!("{}_w2", tag), vi, &w2);
            o.prline(&format!("{}_wm", tag), vi, &wm);
        };
        triple(&mut o, &g, "T");
        if vi % 5 == 0 {
            for m in 0..4 {
                g.set_bypass_mask(Some(&masks[m][..]));
                let tag = format!("M{}", m);
                triple(&mut o, &g, &tag);
                o1.fill(SENT);
                g.compute_areas(&d, &mut o1, n); // unmasked
                o.prline(&format!("{}_A", tag), vi, &o1);
            }
            g.set_bypass_mask(None);
            triple(&mut o, &g, "M_off");
        }
        if vi % 7 == 0 {
            for nt in [2, 3, 5, 0, -3] {
                for with_mask in [0, 1] {
                    if with_mask == 1 {
                        g.set_bypass_mask(Some(&masks[3][..]));
                    }
                    let mut a1 = vec![SENT; un];
                    let mut a2 = vec![SENT; un];
                    let mut am = vec![SENT; un];
                    let mut h1 = vec![SENT; un];
                    let mut hm = vec![SENT; un];
                    let mut w1 = vec![SENT; un];
                    let mut w2 = vec![SENT; un];
                    let mut wm = vec![SENT; un];
                    let workers = if nt < 1 { 1 } else { nt };
                    for tid in 0..workers {
                        g.compute_area_hyd_rad_triple_team(
                            &d, &d2, &dm, &mut a1, &mut a2, &mut am, &mut h1, &mut hm, n, tid, nt,
                        );
                        g.compute_widths_triple_team(&d, &d2, &dm, &mut w1, &mut w2, &mut wm, n, tid, nt);
                    }
                    if with_mask == 1 {
                        g.set_bypass_mask(None);
                    }
                    let t = format!("TM{}_{}", nt, with_mask);
                    o.prline(&format!("{}_a1", t), vi, &a1);
                    o.prline(&format!("{}_a2", t), vi, &a2);
                    o.prline(&format!("{}_am", t), vi, &am);
                    o.prline(&format!("{}_h1", t), vi, &h1);
                    o.prline(&format!("{}_hm", t), vi, &hm);
                    o.prline(&format!("{}_w1", t), vi, &w1);
                    o.prline(&format!("{}_w2", t), vi, &w2);
                    o.prline(&format!("{}_wm", t), vi, &wm);
                }
            }
        }
    }
    {
        let mut g2 = XSectGroups::default();
        g2.build(&net, n / 2);
        dump_groups(&mut o, "half", &g2);
        g2.set_bypass_mask(Some(&masks[0][..]));
        g2.build(&net, 0);
        dump_groups(&mut o, "zero", &g2);
        g2.build(&[], -5);
        dump_groups(&mut o, "neg", &g2);
    }

    // ---- every batch kernel over the network's columns
    let cy: Vec<f64> = net.iter().map(|x| x.y_full).collect();
    let ca: Vec<f64> = net.iter().map(|x| x.a_full).collect();
    let cr: Vec<f64> = net.iter().map(|x| x.r_full).collect();
    let cw: Vec<f64> = net.iter().map(|x| x.w_max).collect();
    let cyb: Vec<f64> = net.iter().map(|x| x.y_bot).collect();
    let csb: Vec<f64> = net.iter().map(|x| x.s_bot).collect();
    let crb: Vec<f64> = net.iter().map(|x| x.r_bot).collect();
    let ta: Vec<Option<&[f64]>> = (0..un)
        .map(|i| {
            let k = (i % 5) as i32;
            if k < n_t {
                Some(&tv[k as usize].ar[..])
            } else {
                None
            }
        })
        .collect();
    let tr: Vec<Option<&[f64]>> = (0..un)
        .map(|i| {
            let k2 = ((i + 1) % 5) as i32;
            if k2 < n_t {
                Some(&tv[k2 as usize].hr[..])
            } else {
                None
            }
        })
        .collect();
    let ht = xs_::host_tables();
    let atab: [(&str, &[f64], i32); 7] = [
        ("A_Circ", ht.A_Circ, ht.N_A_Circ),
        ("A_Egg", ht.A_Egg, ht.N_A_Egg),
        ("A_Horseshoe", ht.A_Horseshoe, ht.N_A_Horseshoe),
        ("A_Baskethandle", ht.A_Baskethandle, ht.N_A_Baskethandle),
        ("A_HorizEllipse", ht.A_HorizEllipse, ht.N_A_HorizEllipse),
        ("A_VertEllipse", ht.A_VertEllipse, ht.N_A_VertEllipse),
        ("A_Arch", ht.A_Arch, ht.N_A_Arch),
    ];
    let ytab: [(&str, &[f64], i32); 4] = [
        ("Y_Gothic", ht.Y_Gothic, ht.N_Y_Gothic),
        ("Y_Catenary", ht.Y_Catenary, ht.N_Y_Catenary),
        ("Y_SemiEllip", ht.Y_SemiEllip, ht.N_Y_SemiEllip),
        ("Y_SemiCirc", ht.Y_SemiCirc, ht.N_Y_SemiCirc),
    ];
    let ylut = [
        LutId::YGothic,
        LutId::YCatenary,
        LutId::YSemiEllip,
        LutId::YSemiCirc,
    ];
    let rtab: [(&str, &[f64], i32); 6] = [
        ("R_Egg", ht.R_Egg, ht.N_R_Egg),
        ("R_Horseshoe", ht.R_Horseshoe, ht.N_R_Horseshoe),
        ("R_Baskethandle", ht.R_Baskethandle, ht.N_R_Baskethandle),
        ("R_HorizEllipse", ht.R_HorizEllipse, ht.N_R_HorizEllipse),
        ("R_VertEllipse", ht.R_VertEllipse, ht.N_R_VertEllipse),
        ("R_Arch", ht.R_Arch, ht.N_R_Arch),
    ];
    let wtab: [(&str, &[f64], i32); 10] = [
        ("W_Egg", ht.W_Egg, ht.N_W_Egg),
        ("W_Horseshoe", ht.W_Horseshoe, ht.N_W_Horseshoe),
        ("W_Gothic", ht.W_Gothic, ht.N_W_Gothic),
        ("W_Catenary", ht.W_Catenary, ht.N_W_Catenary),
        ("W_SemiEllip", ht.W_SemiEllip, ht.N_W_SemiEllip),
        ("W_BasketHandle", ht.W_BasketHandle, ht.N_W_BasketHandle),
        ("W_SemiCirc", ht.W_SemiCirc, ht.N_W_SemiCirc),
        ("W_HorizEllipse", ht.W_HorizEllipse, ht.N_W_HorizEllipse),
        ("W_VertEllipse", ht.W_VertEllipse, ht.N_W_VertEllipse),
        ("W_Arch", ht.W_Arch, ht.N_W_Arch),
    ];
    let nt_tbl = N_TRANSECT_TBL;
    let mut vi = 0;
    while vi < nv {
        let d = dvec(vi);
        let mut ov = vec![SENT; un];
        let mut o2 = vec![SENT; un];
        let k = |o: &mut Out, ov: &mut Vec<f64>, tag: &str| {
            o.prline(tag, vi, ov);
            ov.fill(SENT);
        };
        xb::area_circular(&d, &cy, &ca, &mut ov, n);
        k(&mut o, &mut ov, "k_area_circular");
        xb::area_rect(&d, &cw, &mut ov, n);
        k(&mut o, &mut ov, "k_area_rect");
        xb::area_trapezoidal(&d, &cyb, &csb, &mut ov, n);
        k(&mut o, &mut ov, "k_area_trapezoidal");
        xb::area_triangular(&d, &csb, &mut ov, n);
        k(&mut o, &mut ov, "k_area_triangular");
        xb::area_parabolic(&d, &crb, &mut ov, n);
        k(&mut o, &mut ov, "k_area_parabolic");
        xb::area_powerfunc(&d, &csb, &crb, &mut ov, n);
        k(&mut o, &mut ov, "k_area_powerfunc");
        for (name, t, tn) in atab.iter() {
            xb::area_tabulated(&d, &cy, &ca, t, *tn, &mut ov, n);
            k(&mut o, &mut ov, &format!("k_area_tabulated_{}", name));
        }
        for j in 0..4 {
            let (name, t, tn) = ytab[j];
            xb::area_inv_tabulated(&d, &cy, &ca, t, tn, None, &mut ov, n);
            k(&mut o, &mut ov, &format!("k_area_inv_{}", name));
            xb::area_inv_tabulated(&d, &cy, &ca, t, tn, ht.lut(ylut[j]), &mut ov, n);
            k(&mut o, &mut ov, &format!("k_area_inv_lut_{}", name));
        }
        xb::perlink_tabulated(&d, &cy, &ca, &ta, nt_tbl, &mut ov, n);
        k(&mut o, &mut ov, "k_perlink");
        xb::perlink_tabulated_pair(&d, &cy, &ca, &cr, &ta, &tr, nt_tbl, &mut ov, &mut o2, n);
        o.prline("k_pair_b", vi, &o2);
        o2.fill(SENT);
        k(&mut o, &mut ov, "k_pair_a");
        xb::hydrad_circular(&d, &cy, &cr, &mut ov, n);
        k(&mut o, &mut ov, "k_hydrad_circular");
        xb::area_hydrad_circular(&d, &cy, &ca, &cr, &mut ov, &mut o2, n);
        o.prline("k_ahc_r", vi, &o2);
        o2.fill(SENT);
        k(&mut o, &mut ov, "k_ahc_a");
        xb::hydrad_trapezoidal(&d, &cyb, &csb, &crb, &mut ov, n);
        k(&mut o, &mut ov, "k_hydrad_trapezoidal");
        xb::hydrad_triangular(&d, &csb, &crb, &mut ov, n);
        k(&mut o, &mut ov, "k_hydrad_triangular");
        xb::hydrad_rect(&d, &cw, &mut ov, n);
        k(&mut o, &mut ov, "k_hydrad_rect");
        xb::hydrad_rect_closed(&d, &cw, &ca, &mut ov, n);
        k(&mut o, &mut ov, "k_hydrad_rect_closed");
        xb::hydrad_rect_open(&d, &cw, &csb, &mut ov, n);
        k(&mut o, &mut ov, "k_hydrad_rect_open");
        for (name, t, tn) in rtab.iter() {
            xb::hydrad_tabulated(&d, &cy, &cr, t, *tn, &mut ov, n);
            k(&mut o, &mut ov, &format!("k_hydrad_tabulated_{}", name));
        }
        xb::width_circular(&d, &cy, &cw, &mut ov, n);
        k(&mut o, &mut ov, "k_width_circular");
        xb::width_trapezoidal(&d, &cyb, &csb, &mut ov, n);
        k(&mut o, &mut ov, "k_width_trapezoidal");
        xb::width_triangular(&d, &csb, &mut ov, n);
        k(&mut o, &mut ov, "k_width_triangular");
        xb::width_rect(&cw, &mut ov, n);
        k(&mut o, &mut ov, "k_width_rect");
        xb::width_rect_closed(&d, &cy, &cw, &mut ov, n);
        k(&mut o, &mut ov, "k_width_rect_closed");
        for (name, t, tn) in wtab.iter() {
            xb::width_tabulated(&d, &cy, &cw, t, *tn, &mut ov, n);
            k(&mut o, &mut ov, &format!("k_width_tabulated_{}", name));
        }
        // count 0 and negative: no writes
        xb::area_circular(&d, &cy, &ca, &mut ov, 0);
        xb::perlink_tabulated(&d, &cy, &ca, &ta, nt_tbl, &mut ov, -4);
        k(&mut o, &mut ov, "k_count0");
        vi += 3;
    }

    std::io::stdout().write_all(o.0.as_bytes()).unwrap();
}
