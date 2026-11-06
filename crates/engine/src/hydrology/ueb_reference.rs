// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/surface/UebReference.cpp and
// UebReference.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright
// 2026 Caleb Buahin, licensed under the Apache License, Version 2.0. UebReference.cpp is itself
// an f2c translation of the Utah Energy Balance snow model (Tarboton and Luce; Mahat and
// Tarboton's modified force-restore surface temperature).
// Modified: translated from C++ to Rust, 6 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The modified force-restore point model (namespace force_restore: snowueb2_ and what it
// calls), entered through forceRestoreStep.
//
// Differences in form, never in result:
//   * The f2c pointer arguments are values and &mut, the one-based padded arrays keep their
//     padding slot so every index reads as upstream's. The thread_local globals (diagnostics,
//     age_dust, mean_surface, mean_bulk, visible_age, infrared_age, snow_fraction) and the
//     ts_save__ common block are one struct, `G`, set by force_restore_step as upstream sets
//     the globals before each call. The function-local `static thread_local` constants are
//     consts; predicorr_'s ncall counter is never read and is left out.
//   * surfeb_ and qcest_ assign their parameters es, z__ and fstab from the parameter table
//     (`*es = param[3]` ...), the values the caller already passed; the writes are left out.
//   * The Fortran list-directed writes (s_wsle ... e_wsle) only count: s_wsle increments
//     `diagnostics`, as upstream's stub does. s_stop throws std::runtime_error("UEB reference
//     rejected a negative melt rate."); here an Err with that message.
//   * f2c's min/max templates (std::min / std::max on the common type) are cmin / cmax.
//
// NOT YET TRANSLATED: the equilibrium-temperature model (namespace equilibrium, snow_ and its
// helpers, entered through equilibriumStep). ueb::advanceSubcatchments refuses a deck whose
// UEB parameters select it.

use crate::core::cmath::{cmax, cmin, cvt_i32};

/// The model's per-call globals (upstream's thread_local variables and the ts_save__ block).
#[derive(Debug, Clone, Copy, Default)]
pub struct G {
    pub diagnostics: i32,
    pub age_dust: f64,
    pub mean_surface: f64,
    pub mean_bulk: f64,
    pub visible_age: f64,
    pub infrared_age: f64,
    pub snow_fraction: f64,
    pub ts_old: f64,
    pub tave_old: f64,
    pub ts_ave: f64,
    pub tave_ave: f64,
}

const NEG_MELT: &str = "UEB reference rejected a negative melt rate.";

fn s_wsle(g: &mut G) {
    g.diagnostics += 1;
}

/// The arguments srftmp_ hands surfeb_ and qfm_ hands qcest_, in upstream's order (ts first).
#[derive(Clone, Copy)]
struct Seb<'a> {
    rkn: f64,
    ws: f64,
    tak: f64,
    z: f64,
    g: f64,
    fstab: f64,
    qp: f64,
    densa: f64,
    cp: f64,
    hneu: f64,
    pr: f64,
    ea: f64,
    tk: f64,
    dens: f64,
    cs: f64,
    rs: f64,
    tavek: f64,
    qsn: f64,
    qli: f64,
    fc: f64,
    es: f64,
    sbc: f64,
    qnetob: f64,
    iradfl: i32,
    param: &'a [f64; 23],
    itsmethod: i32,
    w: f64,
    dt: f64,
    ub: f64,
    refdepth: f64,
}

fn svpi(t: f64) -> f64 {
    let r = t * (t * (t * (t * (t * (t * 1.838826904e-10 + 4.838803174e-8) + 5.82472028e-6) + 4.176223716e-4) + 0.01886013408) + 0.503469897) + 6.109177956;
    r * 100.0
}

fn svpw(t: f64) -> f64 {
    let r = t * (t * (t * (t * (t * (t * 6.136820929e-11 + 2.034080948e-8) + 3.031240936e-6) + 2.650648471e-4) + 0.01428945805) + 0.4436518521) + 6.107799961;
    r * 100.0
}

fn svp(t: f64) -> f64 {
    if t >= 0.0 {
        svpw(t)
    } else {
        svpi(t)
    }
}

fn rkinst(rkn: f64, ws: f64, ta: f64, ts: f64, z: f64, g: f64, fstab: f64) -> f64 {
    let mut ret;
    if ws <= 0.0 {
        ret = 0.0;
    } else {
        let rich = g * (ta - ts) * z / (ws * ws * ta);
        if rich >= 0.0 {
            ret = rkn / (rich * 10.0 + 1.0);
        } else {
            let d3 = 1.0 - rich * 16.0;
            let d2 = d3.powf(0.75); // LIBM
            ret = rkn * cmin(3.0, d2);
        }
    }
    ret = rkn + fstab * (ret - rkn);
    ret
}

fn partsnow(g: &G, p: f64, ta: f64, tr: f64, ts: f64) -> f64 {
    if g.snow_fraction >= 0.0 {
        return p * g.snow_fraction;
    }
    if ta < ts {
        p
    } else if ta > tr {
        0.0
    } else {
        p * (tr - ta) / (tr - ts)
    }
}

#[allow(clippy::too_many_arguments)]
fn qpf(pr: f64, ta: f64, to: f64, ps: f64, rhow: f64, hf: f64, cw: f64, cs: f64) -> f64 {
    let (train, tsnow) = if ta > to { (ta, to) } else { (to, ta) };
    pr * rhow * (hf + cw * (train - to)) + ps * rhow * cs * (tsnow - to)
}

/// turbflux_: returns (qh, qe, e).
#[allow(clippy::too_many_arguments)]
fn turbflux(pr: f64, ra: f64, ta: f64, tk: f64, ts: f64, z: f64, g: f64, cp: f64, rkn: f64, ws: f64, ea: f64, rhow: f64, hneu: f64, fstab: f64) -> (f64, f64, f64) {
    let tak = ta + tk;
    let tsk = ts + tk;
    let rkin = rkinst(rkn, ws, tak, tsk, z, g, fstab);
    let rhoa = pr / (ra * tak);
    let qh = rhoa * (ta - ts) * cp * rkin;
    let es = svpi(ts);
    let qe = hneu * 0.622 / (ra * tak) * rkin * (ea - es);
    let e = -qe / (rhow * hneu);
    (qh, qe, e)
}

#[allow(clippy::too_many_arguments)]
fn fmelt(g: &mut G, ub: f64, rhow: f64, w: f64, hf: f64, lc: f64, rid: f64, ks: f64, prain: f64) -> Result<f64, String> {
    let mut ret;
    if ub < 0.0 {
        ret = 0.0;
    } else if w <= 0.0 {
        ret = prain;
        if prain <= 0.0 {
            ret = 0.0;
        }
    } else {
        let mut uu = ub / (rhow * w * hf);
        if uu > 0.99 {
            uu = 0.99;
        }
        let ss = if uu / (1.0 - uu) <= lc { 0.0 } else { (uu / ((1.0 - uu) * rid) - lc / rid) / (1.0 - lc / rid) };
        ret = ks * (ss * (ss * ss));
    }
    if ret < 0.0 {
        s_wsle(g);
        return Err(NEG_MELT.to_string());
    }
    Ok(ret)
}

#[allow(clippy::too_many_arguments)]
fn tavg(ub: f64, w: f64, rhow: f64, cs: f64, to: f64, rhog: f64, de: f64, cg: f64, hf: f64) -> f64 {
    let snhc = rhow * w * cs;
    let shc = rhog * de * cg;
    let chc = snhc + shc;
    if ub <= 0.0 {
        ub / chc
    } else {
        let al = ub / (rhow * hf);
        if w > al {
            to
        } else {
            (ub - w * rhow * hf) / chc
        }
    }
}

/// prehelp_: rates adjusted so the provisional water is not negative. Returns nothing; writes
/// fm, mr, e, q, qm, qe.
#[allow(clippy::too_many_arguments)]
fn prehelp(w1: f64, w: f64, dt: f64, fm: &mut f64, fm1: f64, fac: f64, ps: f64, prain: f64, e: &mut f64, rhow: f64, hf: f64, q: &mut f64, qm: &mut f64, mr: &mut f64, qe: &mut f64, hneu: f64) {
    *fm = (w1 - w) / dt * fac - fm1;
    *mr = cmax(0.0, ps + prain - *fm - *e);
    *e = ps + prain - *fm - *mr;
    let qother = *q + *qm - *qe;
    *qm = *mr * rhow * hf;
    *qe = -(*e) * rhow * hneu;
    *q = qother - *qm + *qe;
}

fn albedo(g: &G, tausn: f64, coszen: f64, d: f64, aep: f64, abg: f64, avo: f64, airo: f64) -> f64 {
    let b = 2.0;
    let cs = g.visible_age;
    let cn = g.infrared_age;
    let fage = tausn / (tausn + 1.0);
    let fzen = if coszen < 0.5 { 1.0 / b * ((b + 1.0) / (b * 2.0 * coszen + 1.0) - 1.0) } else { 0.0 };
    let avd = (1.0 - cs * fage) * avo;
    let avis = avd + fzen * 0.4 * (1.0 - avd);
    let aird = (1.0 - cn * fage) * airo;
    let anir = aird + fzen * 0.4 * (1.0 - aird);
    let mut ret = (avis + anir) / 2.0;
    if d < aep {
        let rr = (1.0 - d / aep) * (-d * 0.5 / aep).exp(); // LIBM
        ret = rr * abg + (1.0 - rr) * ret;
    }
    ret
}

fn agesn(g: &G, tausn: &mut f64, dt: f64, ps: f64, tsurf: f64, tk: f64, dnews: f64) {
    let tsk = tsurf + tk;
    let r1 = ((1.0 / tk - 1.0 / tsk) * 5e3).exp(); // LIBM
    let mut d1 = r1;
    d1 *= d1;
    let d2 = d1;
    d1 *= d1;
    let mut r2 = d2 * (d1 * d1);
    if r2 > 1.0 {
        r2 = 1.0;
    }
    let r3 = g.age_dust;
    *tausn += (r1 + r2 + r3) * 0.0036 * dt;
    if ps > 0.0 {
        if dnews > 0.0 {
            let x = *tausn * (1.0 - ps * dt / dnews);
            *tausn = cmax(x, 0.0);
        } else {
            *tausn = 0.0;
        }
    }
}

/// lane_: returns (value, ze).
#[allow(clippy::too_many_arguments)]
fn lane(lans: f64, lang: f64, zs: f64, rho: f64, rhog: f64, cs: f64, cg: f64, r: f64, w1day: f64) -> (f64, f64) {
    let fkappas = lans / (rho * cs);
    let fkappag = lang / (rhog * cg);
    let d1 = (fkappas * 2.0 / w1day).sqrt();
    let d2 = (fkappag * 2.0 / w1day).sqrt();
    let mut ret = cmin(r, zs / d1) * d1 / lans + cmax(0.0, r - zs / d1) * d2 / lang;
    let ze = cmin(r, zs / d1) * d1 + cmax(0.0, r - zs / d1) * d2;
    if ret != 0.0 {
        ret = 1.0 / ret * ze;
    }
    (ret, ze)
}

#[allow(clippy::too_many_arguments)]
fn refdep(flans: f64, a: f64, b: f64, hf: f64, rhom: f64, dt: f64, x1: f64) -> f64 {
    let temp = flans * flans + b * 2.0 * (-a * flans * dt / (hf * rhom) + flans * x1 + b * 0.5 * x1 * x1);
    if temp < 0.0 || a > 0.0 || b == 0.0 {
        0.0
    } else {
        (-flans + temp.sqrt()) / b
    }
}

/// grad_: (a, b) unchanged when t2 == t1.
fn grad(qc1: f64, qc2: f64, t1: f64, t2: f64, a: &mut f64, b: &mut f64) {
    if t2 - t1 != 0.0 {
        *b = (qc2 - qc1) / (t1 - t2);
        *a = qc1 + *b * t1;
    }
}

const RHOW: f64 = 1e3;
const W1DAY: f64 = 0.261799;

/// surfeb_: the surface energy balance at surface temperature ts (K).
fn surfeb(g: &G, s: &Seb, ts: f64) -> f64 {
    let param = s.param;
    let cg = param[4];
    let rho = param[7];
    let rhog = param[8];
    let lans = param[15];
    let lang = param[16];
    let wlf = param[17];
    let rd1 = param[18];
    let zs = s.w * RHOW / rho;
    let rkin = rkinst(s.rkn, s.ws, s.tak, ts, s.z, s.g, s.fstab);
    let qcs;
    if s.itsmethod == 1 {
        qcs = s.dens * s.cs * s.rs * (ts - s.tavek);
    } else if s.itsmethod == 2 {
        let fkappas = lans / (rho * s.cs);
        let d1 = (fkappas * 2.0 / W1DAY).sqrt();
        let lane_ze;
        if zs >= rd1 * d1 {
            lane_ze = lans / (rd1 * d1);
        } else {
            let (v, ze) = lane(lans, lang, zs, rho, rhog, s.cs, cg, rd1, W1DAY);
            lane_ze = v / ze;
        }
        qcs = lane_ze * (ts - s.tavek);
    } else if s.itsmethod == 3 {
        let fkappas = lans / (rho * s.cs);
        let d1 = (fkappas * 2.0 / W1DAY).sqrt();
        let (lane_ze, ze);
        if zs >= rd1 * d1 {
            lane_ze = lans / (rd1 * d1);
            ze = rd1 * d1;
        } else {
            let (v, z) = lane(lans, lang, zs, rho, rhog, s.cs, cg, rd1, W1DAY);
            ze = z;
            lane_ze = v / ze;
        }
        let de = ze / rd1;
        let lane_de = lane_ze / de * ze;
        qcs = lane_de * (ts - g.ts_old) / (W1DAY * s.dt) + lane_ze * (ts - s.tavek);
    } else {
        let fkappas = lans / (rho * s.cs);
        let fkappag = lang / (rhog * cg);
        let d1 = (fkappas * 2.0 / W1DAY).sqrt();
        let dlf = (fkappag * 2.0 / wlf).sqrt();
        let (lane_ze, ze);
        if zs >= rd1 * d1 {
            lane_ze = lans / (rd1 * d1);
            ze = rd1 * d1;
        } else {
            let (v, z) = lane(lans, lang, zs, rho, rhog, s.cs, cg, rd1, W1DAY);
            ze = z;
            lane_ze = v / ze;
        }
        let (lane_ze2, ze2);
        if zs >= rd1 * dlf {
            lane_ze2 = lans / (rd1 * dlf);
            ze2 = rd1 * dlf;
        } else {
            let (v, z) = lane(lans, lang, zs, rho, rhog, s.cs, cg, rd1, wlf);
            ze2 = z;
            lane_ze2 = v / ze2;
        }
        let de = ze / rd1;
        let lane_de = lane_ze / de * ze;
        let de2 = ze2 / rd1;
        let lane_de2 = lane_ze2 / de2 * ze2;
        if s.ub <= 0.0 || s.refdepth <= 0.0 {
            qcs = lane_de * (ts - g.ts_old) / (W1DAY * s.dt) + lane_ze * (ts - g.ts_ave) + lane_de2 * (g.ts_ave - g.tave_ave);
        } else if s.refdepth > rd1 * d1 {
            qcs = lane_de * (ts - g.ts_old) / (W1DAY * s.dt) + lane_ze * (ts - g.ts_ave) + lane_de2 * (g.ts_ave - g.tave_ave);
        } else {
            qcs = lane_ze * ze * (ts - s.tk) / s.refdepth;
        }
    }
    let mut ret = s.qp + rkin * s.densa * s.cp * (s.tak - ts) + s.hneu * s.densa * 0.622 * rkin / s.pr * (s.ea - svp(ts - s.tk)) - qcs;
    if s.iradfl == 0 {
        let d = ts * ts;
        ret = ret + s.qsn + s.qli - (1.0 - s.fc) * s.es * s.sbc * (d * d);
    } else {
        ret += s.qnetob;
    }
    ret
}

/// qcest_: surfeb_ without the conduction term.
fn qcest(s: &Seb, ts: f64) -> f64 {
    let rkin = rkinst(s.rkn, s.ws, s.tak, ts, s.z, s.g, s.fstab);
    let mut ret = s.qp + rkin * s.densa * s.cp * (s.tak - ts) + s.hneu * s.densa * 0.622 * rkin / s.pr * (s.ea - svp(ts - s.tk));
    if s.iradfl == 0 {
        let d = ts * ts;
        ret = ret + s.qsn + s.qli - (1.0 - s.fc) * s.es * s.sbc * (d * d);
    } else {
        ret += s.qnetob;
    }
    ret
}

/// srftmp_: the surface temperature (C), Newton-type iteration then bisection; returns
/// (tsurf, smelt).
fn srftmp(g: &mut G, s: &Seb) -> (f64, f64) {
    const TOL: f64 = 0.05;
    const NITERMAX: i32 = 10;
    let fff = (273.0 - TOL) / 273.0;
    let mut ts = s.tak;
    let mut er = TOL * 2.0;
    let mut niter = 0;
    let mut go_l11 = false;
    while er > TOL && niter < NITERMAX {
        let tslast = ts;
        let f1 = surfeb(g, s, ts);
        let f2 = surfeb(g, s, fff * ts);
        ts -= (1.0 - fff) * ts * f1 / (f1 - f2);
        if ts < s.tak - 50.0 {
            go_l11 = true;
            break;
        }
        er = (ts - tslast).abs();
        niter += 1;
    }
    let mut to_l10 = !go_l11 && er <= TOL;
    if !to_l10 {
        // L11: bisection
        let mut tlb = s.tak - 20.0;
        let mut tub = s.tak + 10.0;
        let mut flb = surfeb(g, s, tlb);
        let mut fub = surfeb(g, s, tub);
        let mut ibtowrite = 0;
        let mut bisect = true;
        if flb * fub > 0.0 {
            tlb = s.tak - 150.0;
            tub = s.tak + 100.0;
            flb = surfeb(g, s, tlb);
            fub = surfeb(g, s, tub);
            ibtowrite = 1;
            if flb * fub > 0.0 {
                for _ in 0..5 {
                    s_wsle(g);
                }
                ts = s.tk;
                to_l10 = true;
                bisect = false;
            } else {
                for _ in 0..6 {
                    s_wsle(g);
                }
            }
        }
        if bisect {
            let n = cvt_i32(((tub - tlb) / TOL).ln() / 2.0f64.ln()); // LIBM
            for _ in 1..=n {
                ts = (tub + tlb) * 0.5;
                let f1 = surfeb(g, s, ts);
                if f1 > 0.0 {
                    tlb = ts;
                } else {
                    tub = ts;
                }
            }
            if ibtowrite == 1 {
                for _ in 0..3 {
                    s_wsle(g);
                }
            }
        }
        let _ = to_l10;
    }
    // L10
    ts -= s.tk;
    if s.w > 0.0 && ts > 0.0 {
        let ret = 0.0;
        let smelt = surfeb(g, s, ret + s.tk);
        (ret, smelt)
    } else {
        (ts, 0.0)
    }
}

/// The outputs qfm_ writes.
#[derive(Debug, Clone, Copy, Default)]
struct Rates {
    fm: f64,
    q: f64,
    qm: f64,
    mr: f64,
    qe: f64,
    e: f64,
    tsurf: f64,
    tave: f64,
    qh: f64,
    qnet: f64,
    smelt: f64,
}

/// The inputs qfm_ reads that predicorr_ passes through unchanged.
#[derive(Clone, Copy)]
struct QfmIn<'a> {
    a: f64,
    ta: f64,
    prain: f64,
    ps: f64,
    ws: f64,
    rh: f64,
    qsi: f64,
    qli: f64,
    rkn: f64,
    iradfl: i32,
    qnetob: f64,
    rid: f64,
    param: &'a [f64; 23],
    sitev: &'a [f64; 6],
    itsmethod: i32,
    dt: f64,
}

/// qfm_: the coupled water and energy rates for inventory (ub, w).
fn qfm(g: &mut G, ub: f64, w: f64, x: &QfmIn, r: &mut Rates, refdepth: &mut f64, totalrefdepth: f64) -> Result<(), String> {
    const TO: f64 = 0.0;
    const TK: f64 = 273.15;
    const SBC: f64 = 2.0747e-7;
    const HF: f64 = 333.5;
    const HNEU: f64 = 2834.0;
    const CW: f64 = 4.18;
    const CS: f64 = 2.09;
    const CP: f64 = 1.005;
    const RA: f64 = 287.0;
    const GRAV: f64 = 9.81;
    let param = x.param;
    let es = param[3];
    let cg = param[4];
    let z = param[5];
    let rho = param[7];
    let rhog = param[8];
    let lc = param[9];
    let ks = param[10];
    let de = param[11];
    let lans = param[15];
    let rd1 = param[18];
    let fstab = param[19];
    let gsurf = param[22];
    let fc = x.sitev[1];
    let pr = x.sitev[3];
    let qg = x.sitev[4];
    let qsn = x.qsi * (1.0 - x.a);
    let tak = x.ta + TK;
    let tavek = r.tave + TK;
    let densa = pr / (RA * tak);
    let dens = rho;
    let rhom = lc * rho;
    let fkappas = lans / (rho * CS);
    let ds = (fkappas * 2.0 / W1DAY).sqrt();
    let rs = lans / (ds * rd1 * rho * CS);
    g.tave_old = r.tave;
    let ea = svpw(x.ta) * x.rh;
    let qp = qpf(x.prain, x.ta, TO, x.ps, RHOW, HF, CW, CS);
    r.tave = tavg(ub, w, RHOW, CS, TO, rhog, de, cg, HF);
    let mut seb = Seb { rkn: x.rkn, ws: x.ws, tak, z, g: GRAV, fstab, qp, densa, cp: CP, hneu: HNEU, pr, ea, tk: TK, dens, cs: CS, rs, tavek, qsn, qli: x.qli, fc, es, sbc: SBC, qnetob: x.qnetob, iradfl: x.iradfl, param, itsmethod: x.itsmethod, w, dt: x.dt, ub, refdepth: *refdepth };
    if ub > 0.0 && x.ps <= 0.0 && x.prain <= 0.0 && totalrefdepth <= rd1 * ds && w > 0.0 {
        let qc1 = qcest(&seb, TK);
        let qc2 = qcest(&seb, TK - 0.01);
        let mut var_a = 0.0;
        let mut var_b = 0.0;
        grad(qc1, qc2, 0.0, -0.01, &mut var_a, &mut var_b);
        let x1 = refdep(lans, var_a, var_b, HF, rhom, x.dt, *refdepth);
        *refdepth = x1;
    } else {
        *refdepth = 0.0;
    }
    // srftmp_ recomputes its own tavek from the UPDATED tave, and drops a positive qp over a
    // snow-free point (qcest_ above read the entry tavek and the unadjusted qp)
    seb.refdepth = *refdepth;
    seb.tavek = r.tave + TK;
    if seb.w <= 0.0 && seb.qp > 0.0 {
        seb.qp = 0.0;
    }
    let (tsurf, smelt) = srftmp(g, &seb);
    r.tsurf = tsurf;
    r.smelt = smelt;
    let mut d1 = r.tsurf + TK;
    d1 *= d1;
    let qle = (1.0 - fc) * es * SBC * (d1 * d1);
    let qlnet = x.qli - qle;
    let (qh, qe, e) = turbflux(pr, RA, x.ta, TK, r.tsurf, z, GRAV, CP, x.rkn, x.ws, ea, RHOW, HNEU, fstab);
    r.qh = qh;
    r.qe = qe;
    r.e = e;
    r.mr = fmelt(g, ub, RHOW, w, HF, lc, x.rid, ks, x.prain)?;
    r.qm = r.mr * RHOW * (HF + (r.tave - TO) * CW);
    r.mr += r.smelt / (HF * RHOW) * gsurf;
    r.qm += r.smelt * gsurf;
    if x.iradfl == 0 {
        r.qnet = x.qsi * (1.0 - x.a) + qlnet;
    } else {
        r.qnet = x.qnetob;
    }
    r.q = r.qnet + qp + qg + r.qh + r.qe - r.qm;
    r.fm = x.prain + x.ps - r.mr - r.e;
    Ok(())
}

/// predicorr_: the predictor-corrector over one step; updates ub, w and the rates.
fn predicorr(g: &mut G, ub: &mut f64, w: &mut f64, x: &QfmIn, r: &mut Rates, refdepth: &mut f64, totalrefdepth: f64) -> Result<(), String> {
    const WTOL: f64 = 0.025;
    const UTOL: f64 = 2e3;
    const HF: f64 = 333.5;
    const HNEU: f64 = 2834.0;
    let dt = x.dt;
    qfm(g, *ub, *w, x, r, refdepth, totalrefdepth)?;
    let mut w1 = *w + dt * r.fm;
    if w1 < 0.0 {
        w1 = 0.0;
        prehelp(w1, *w, dt, &mut r.fm, 0.0, 1.0, x.ps, x.prain, &mut r.e, RHOW, HF, &mut r.q, &mut r.qm, &mut r.mr, &mut r.qe, HNEU);
    }
    let mut ub1 = *ub + dt * r.q;
    let mut q1 = r.q;
    let fm1 = r.fm;
    let qh1 = r.qh;
    let qe1 = r.qe;
    let e1 = r.e;
    let mut mr1 = r.mr;
    let smelt1 = r.smelt;
    let mut qm1 = r.qm;
    let tsurf1 = r.tsurf;
    let qnet1 = r.qnet;
    qfm(g, ub1, w1, x, r, refdepth, totalrefdepth)?;
    let mut w2 = *w + dt / 2.0 * (fm1 + r.fm);
    if w2 < 0.0 {
        w2 = 0.0;
        prehelp(w2, *w, dt, &mut r.fm, fm1, 2.0, x.ps, x.prain, &mut r.e, RHOW, HF, &mut r.q, &mut r.qm, &mut r.mr, &mut r.qe, HNEU);
    }
    let mut ub2 = *ub + dt / 2.0 * (q1 + r.q);
    let mut niter = -3;
    let imax = 1;
    while ((w2 - w1).abs() > WTOL || (ub2 - ub1).abs() > UTOL) && niter < imax {
        w1 = w2;
        ub1 = ub2;
        qfm(g, ub1, w1, x, r, refdepth, totalrefdepth)?;
        w2 = *w + dt / 2.0 * (fm1 + r.fm);
        if w2 < 0.0 {
            w2 = 0.0;
            prehelp(w2, *w, dt, &mut r.fm, fm1, 2.0, x.ps, x.prain, &mut r.e, RHOW, HF, &mut r.q, &mut r.qm, &mut r.mr, &mut r.qe, HNEU);
        }
        ub2 = *ub + dt / 2.0 * (q1 + r.q);
        niter += 1;
        if niter >= 1 {
            let ae = (q1 + r.q + qm1 + r.qm) * 0.5 * dt;
            if *ub > 0.0 && ae > 0.0 && *w > 0.0 {
                let e2 = (r.e + e1) * 0.5;
                let rlf = (*ub + ae) / (RHOW * *w * HF);
                if rlf >= 1.0 {
                    r.mr = *w / dt + (x.ps + x.prain - e2);
                    if r.mr < 0.0 {
                        r.mr = 0.0;
                    }
                    r.qm = r.mr * RHOW * HF;
                    w2 = 0.0;
                    ub2 = *ub + ae - r.qm * dt;
                } else {
                    let rr = *w / *ub;
                    ub2 = (RHOW * HF * (*w + (x.ps + x.prain - e2) * dt) - ae - *ub) / (RHOW * HF * rr - 1.0);
                    w2 = rr * ub2;
                    if w2 < 0.0 {
                        w2 = 0.0;
                    }
                    r.mr = (*w - w2) / dt - e2 + x.ps + x.prain;
                    if r.mr < 0.0 {
                        r.mr = 0.0;
                        w2 = *w + (x.ps + x.prain - e2) * dt;
                        if w2 < 0.0 {
                            w2 = 0.0;
                        }
                    }
                    r.qm = r.mr * RHOW * HF;
                    ub2 = *ub + ae - r.qm * dt;
                }
                if r.mr < 0.0 {
                    s_wsle(g);
                }
                if w2 < 0.0 {
                    s_wsle(g);
                }
                r.q = ae / dt - r.qm;
                qm1 = r.qm;
                mr1 = r.mr;
                q1 = r.q;
            }
        }
    }
    *w = w2;
    *ub = ub2;
    r.qh = (r.qh + qh1) * 0.5;
    r.qe = (r.qe + qe1) * 0.5;
    r.e = (r.e + e1) * 0.5;
    r.mr = (r.mr + mr1) * 0.5;
    r.qm = (r.qm + qm1) * 0.5;
    r.tsurf = (r.tsurf + tsurf1) * 0.5;
    r.qnet = (r.qnet + qnet1) * 0.5;
    r.q = (r.q + q1) * 0.5;
    r.smelt = (r.smelt + smelt1) * 0.5 / (HF * RHOW);
    Ok(())
}

/// snowueb2_ for one step (nt = 1, nstepday = 1). `input` is the 15-slot table (slots 8..14),
/// `sitev`, `statev`, `param` and `outv` keep their padding slot 0; `iflag` slots 1..5.
#[allow(clippy::too_many_arguments)]
fn snowueb2(g: &mut G, dt: f64, input: &[f64; 15], sitev: &[f64; 6], statev: &mut [f64; 6], tsprevday: &mut f64, taveprevday: &mut f64, param: &[f64; 23], iflag: &[i32; 6], outv: &mut [f64; 15]) -> Result<(), String> {
    const TO: f64 = 0.0;
    const K: f64 = 0.4;
    const HFF: f64 = 3600.0;
    const RHOI: f64 = 917.0;
    const TK: f64 = 273.15;
    const HF: f64 = 333.5;
    const CS: f64 = 2.09;
    let tr = param[1];
    let ts = param[2];
    let cg = param[4];
    let z = param[5];
    let zo = param[6];
    let rho = param[7];
    let rhog = param[8];
    let lc = param[9];
    let de = param[11];
    let abg = param[12];
    let avo = param[13];
    let anir0 = param[14];
    let dnews = param[21];
    let mut ub = statev[1];
    let mut w = statev[2];
    let mut refdepth;
    let mut totalrefdepth;
    if ub <= 0.0 {
        refdepth = 0.0;
        totalrefdepth = 0.0;
    } else {
        refdepth = statev[4];
        totalrefdepth = statev[5];
    }
    let ub_old = ub;
    let fc = sitev[1];
    let df = sitev[2];
    let aep = sitev[5];
    let iradfl = iflag[1];
    let itsmethod = iflag[5];

    let d1 = (z / zo).ln(); // LIBM
    let cd = K * K * HFF / (d1 * d1) * (1.0 - fc * 0.8);
    let rrhoi = RHOI / RHOW;
    let rrho = rho / RHOW;
    let rid = 1.0 / rrho - 1.0 / rrhoi;

    // i = 1 (nt = 1)
    let ta = input[8];
    let p = input[9];
    let ws = input[10];
    let rh = input[11];
    let mut qsi = 0.0;
    let mut qli = 0.0;
    let mut qnetob = 0.0;
    if iradfl == 0 {
        qsi = input[12];
        qli = input[13];
    } else {
        qnetob = input[12];
    }
    let coszen = input[14];
    g.ts_ave = g.mean_surface + TK;
    g.tave_ave = g.mean_bulk + TK;
    g.ts_old = *tsprevday + TK;
    g.tave_old = *taveprevday + TK;
    if g.ts_old < 0.0 {
        s_wsle(g);
        g.ts_old = TK;
    }
    if g.tave_old < 0.0 {
        s_wsle(g);
        g.tave_old = TK;
    }
    if g.ts_ave < 0.0 {
        s_wsle(g);
        g.ts_ave = TK;
    }
    if g.tave_ave < 0.0 {
        s_wsle(g);
        g.tave_ave = TK;
    }
    let mut ps = partsnow(g, p, ta, tr, ts);
    let prain = p - ps;
    ps *= df;
    let a = if iflag[4] == 1 { albedo(g, statev[3], coszen, w / rrho, aep, abg, avo, anir0) } else { statev[3] };
    let rkn = cd * ws;
    qsi *= 1.0 - fc;
    qli *= 1.0 - fc;
    qnetob *= 1.0 - fc;
    let x = QfmIn { a, ta, prain, ps, ws, rh, qsi, qli, rkn, iradfl, qnetob, rid, param, sitev, itsmethod, dt };
    let mut r = Rates::default();
    predicorr(g, &mut ub, &mut w, &x, &mut r, &mut refdepth, totalrefdepth)?;
    r.tave = tavg(ub, w, RHOW, CS, TO, rhog, de, cg, HF);
    if r.tave > 0.0 {
        r.mr += w / dt;
        r.qm += w / dt * RHOW * HF;
        r.q -= w / dt * RHOW * HF;
        ub -= w * RHOW * HF;
        w = 0.0;
    }
    if iflag[4] == 1 {
        let mut tausn = statev[3];
        agesn(g, &mut tausn, dt, ps, r.tsurf, TK, dnews);
        statev[3] = tausn;
    }
    r.tave = tavg(ub, w, RHOW, CS, TO, rhog, de, cg, HF);
    if lc > 0.0 {
        if refdepth > 0.0 {
            totalrefdepth = refdepth;
        } else if r.mr > 0.0 || (ub > ub_old && ub > 0.0) {
            totalrefdepth = 0.0;
        }
    } else if r.mr > 0.0 || (ub > ub_old && ub > 0.0) {
        totalrefdepth = 0.0;
    }
    if totalrefdepth < 0.0 {
        totalrefdepth = 0.0;
    }
    *tsprevday = r.tsurf;
    *taveprevday = r.tave;
    statev[1] = ub;
    statev[2] = w;
    statev[4] = refdepth;
    statev[5] = totalrefdepth;
    outv[1] = prain;
    outv[2] = ps;
    outv[3] = a;
    outv[4] = r.qh;
    outv[5] = r.qe;
    outv[6] = r.e;
    outv[7] = r.mr;
    outv[8] = r.qm;
    outv[9] = r.q;
    outv[10] = r.fm;
    outv[11] = r.tave;
    outv[12] = r.tsurf;
    outv[13] = r.qnet;
    outv[14] = r.smelt;
    Ok(())
}

/// forceRestoreStep: one modified force-restore step; returns the diagnostics count.
#[allow(clippy::too_many_arguments)]
pub fn force_restore_step(dt: f64, input: &[f64; 15], site: &[f64; 6], state: &mut [f64; 6], surface: &mut f64, average: &mut f64, parameters: &[f64; 23], output: &mut [f64; 15], dust: f64, surface_mean: f64, bulk_mean: f64, cv: f64, cir: f64, fraction: f64) -> Result<i32, String> {
    let mut g = G { diagnostics: 0, age_dust: dust, mean_surface: surface_mean, mean_bulk: bulk_mean, visible_age: cv, infrared_age: cir, snow_fraction: fraction, ..G::default() };
    let flags: [i32; 6] = [0, 0, 0, 0, 1, 4];
    snowueb2(&mut g, dt, input, site, state, surface, average, parameters, &flags, output)?;
    Ok(g.diagnostics)
}
