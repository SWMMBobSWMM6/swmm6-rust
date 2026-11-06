// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/components/ReactionModule/
// ReactionLegacyBinding.cpp, ReactionIntegrator.cpp and ReactionIntegrator.hpp,
// ReactionExpression.cpp (evalReactionExpression), ReactionGeometry.hpp, ReactionSources.hpp
// (reactionBoosterMass), MsxLegacyTransport.cpp (armMsxRoutingLedger, routeLegacyMsx), and
// src/engine/quality/MsxLidReactions.cpp (enabledNode, heldMass, mobileWater, outlet,
// receive) and FvQualityTransport.cpp (retainedSpeciesMass) (HydroCouple/openswmm.engine,
// branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache
// License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The reaction species under QUALITY_SOLVER LEGACY, as QualitySolver::execute runs them each
// routing step: routeLegacyMsx (species transport on the CSTR mirror of the pollutant
// kernels, with its routing ledger), and reactLegacyNodes / reactLegacyLinks (exact
// first-order pollutant decay, then each element's species through ReactionIntegrator::step).
//
// Differences in form, never in result:
//   * Upstream keeps the integrator workspace and the element scratch in process-wide statics
//     (BindingScratch, MsxScratch). Within a call every value read from them is written first
//     (RxWorkspace's Jacobian flags are reset at each RATE group, terms are recomputed where
//     needed, the old-state snapshots are copied in), so here they are locals of each call.
//   * evalReactionExpression's fixed stack is an array of K_RX_MAX_STACK_DEPTH doubles.
//   * The shared-context (RxRateContext) argument of ReactionIntegrator::step is not carried:
//     the legacy binding passes none.
//
// NOT YET TRANSLATED (StepError::NotTranslated; SWMMEngine::stepRouting refuses such a deck
// before any state changes): the ROS2 and BDF2 solvers, EQUIL species (the joint Newton),
// booster sources (reactionBoosterMass with a booster row, reactionSourceFactor), and the
// reaction-species LID controls (msxlid::prepareNodes / reactNodes).

use crate::core::cmath::{cclamp, cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::engine_step_routing::StepError;
use crate::core::options::{Bytes, RoutingModel};
use crate::data::node::NodeType;
use crate::data::reaction::{ReactionAreaUnits, ReactionCoupling, ReactionData, ReactionExprForm, ReactionRateUnits, ReactionSolverKind};
use crate::data::reaction_tokens::{RxExprSpan, RxHydVar, RxToken, K_RX_MAX_STACK_DEPTH};
use crate::hydraulics::link;
use crate::hydraulics::xsection as xsect;
use crate::hydrology::lid_node_step as lidnode;
use crate::quality::quality_routing::{conduit_mixing_inflow, link_is_dry, link_takes_upstream_value, node_is_dry, node_is_reactor, storage_evap_factor, LEGACY_ZERO, ZERO_DEPTH, ZERO_VOLUME};
use crate::transport::reaction_legacy_binding::{ensure_msx_state, legacy_reactions_active};

const K_MAX_SUBSTEPS: i32 = 100000;
const K_MIN_STEP_FRAC: f64 = 1.0e-8;

/// The deck-level reaction features this file does not translate yet (checked by
/// SWMMEngine::stepRouting before the step changes anything).
pub fn reactions_not_translated(ctx: &SimulationContext) -> Option<&'static str> {
    let rx = &ctx.reactions;
    if !rx.lid_controls.is_empty() {
        return Some("reaction species LID controls");
    }
    if rx.sources.iter().any(|s| s.booster != 0) {
        return Some("reaction booster sources");
    }
    if legacy_reactions_active(ctx) && rx.n_species() > 0 {
        if rx.solver == ReactionSolverKind::Ros2 || rx.solver == ReactionSolverKind::Bdf2 {
            return Some("reaction solvers ROS2 and BDF2");
        }
        if rx.tank_form.iter().chain(rx.pipe_form.iter()).any(|f| *f == ReactionExprForm::Equil) {
            return Some("reaction EQUIL species");
        }
    }
    None
}

// ---------------------------------------------------------------------------
// ReactionExpression.cpp: evalReactionExpression
// ---------------------------------------------------------------------------

/// RxEvalEnv: what a compiled expression reads.
struct Env<'a> {
    species: &'a [f64],
    coefs: &'a [f64],
    terms: &'a [f64],
    hydvar: &'a [f64],
    pollutants: Option<&'a [f64]>,
}

/// evalReactionExpression: the postfix span on a fixed stack.
fn eval_expr(pool: &[RxToken], span: &RxExprSpan, env: &Env) -> f64 {
    if span.len <= 0 {
        return 0.0;
    }
    let mut st = [0.0f64; K_RX_MAX_STACK_DEPTH as usize];
    let mut sp: usize = 0;
    for t in &pool[span.begin as usize..(span.begin + span.len) as usize] {
        match t.op {
            RxToken::PUSH_NUM => {
                st[sp] = t.value;
                sp += 1;
            }
            RxToken::PUSH_SPECIES => {
                st[sp] = env.species[t.idx as usize];
                sp += 1;
            }
            RxToken::PUSH_COEF => {
                st[sp] = env.coefs[t.idx as usize];
                sp += 1;
            }
            RxToken::PUSH_TERM => {
                st[sp] = env.terms[t.idx as usize];
                sp += 1;
            }
            RxToken::PUSH_HYDVAR => {
                st[sp] = env.hydvar[t.idx as usize];
                sp += 1;
            }
            RxToken::PUSH_POLLUT => {
                st[sp] = env.pollutants.unwrap()[t.idx as usize];
                sp += 1;
            }
            RxToken::ADD => {
                sp -= 1;
                st[sp - 1] += st[sp];
            }
            RxToken::SUB => {
                sp -= 1;
                st[sp - 1] -= st[sp];
            }
            RxToken::MUL => {
                sp -= 1;
                st[sp - 1] *= st[sp];
            }
            RxToken::DIV => {
                sp -= 1;
                st[sp - 1] /= st[sp];
            }
            RxToken::POW => {
                sp -= 1;
                st[sp - 1] = st[sp - 1].powf(st[sp]); // LIBM
            }
            RxToken::NEG => st[sp - 1] = -st[sp - 1],
            RxToken::F_EXP => st[sp - 1] = st[sp - 1].exp(),     // LIBM
            RxToken::F_LOG => st[sp - 1] = st[sp - 1].ln(),      // LIBM
            RxToken::F_LOG10 => st[sp - 1] = st[sp - 1].log10(), // LIBM
            RxToken::F_SQRT => st[sp - 1] = st[sp - 1].sqrt(),
            RxToken::F_ABS => st[sp - 1] = st[sp - 1].abs(),
            RxToken::F_SGN => {
                st[sp - 1] = if st[sp - 1] > 0.0 {
                    1.0
                } else if st[sp - 1] < 0.0 {
                    -1.0
                } else {
                    0.0
                }
            }
            RxToken::F_STEP => st[sp - 1] = if st[sp - 1] > 0.0 { 1.0 } else { 0.0 },
            RxToken::F_SIN => st[sp - 1] = st[sp - 1].sin(), // LIBM
            RxToken::F_COS => st[sp - 1] = st[sp - 1].cos(), // LIBM
            RxToken::F_TAN => st[sp - 1] = st[sp - 1].tan(), // LIBM
            RxToken::F_MIN => {
                sp -= 1;
                st[sp - 1] = if st[sp - 1] < st[sp] { st[sp - 1] } else { st[sp] };
            }
            RxToken::F_MAX => {
                sp -= 1;
                st[sp - 1] = if st[sp - 1] > st[sp] { st[sp - 1] } else { st[sp] };
            }
            _ => {}
        }
    }
    st[0]
}

// ---------------------------------------------------------------------------
// ReactionIntegrator.cpp
// ---------------------------------------------------------------------------

fn unit_factor(u: ReactionRateUnits) -> f64 {
    match u {
        ReactionRateUnits::Sec => 1.0,
        ReactionRateUnits::Min => 1.0 / 60.0,
        ReactionRateUnits::Hr => 1.0 / 3600.0,
        _ => 1.0 / 86400.0,
    }
}

/// evalTerms: the needed named terms, in order (a term reads earlier ones).
fn eval_terms(rx: &ReactionData, species: &[f64], hydvar: &[f64], pollut: Option<&[f64]>, terms: &mut [f64], needed: &[u8], coefficients: &[f64]) {
    for i in 0..rx.term_expr.len() {
        if needed[i] != 0 {
            let v = {
                let env = Env { species, coefs: coefficients, terms, hydvar, pollutants: pollut };
                eval_expr(&rx.token_pool, &rx.term_expr[i], &env)
            };
            terms[i] = v;
        }
    }
}

/// errNorm: RMS of the error scaled by atol + rtol*|y| (per-species tolerances when set).
fn err_norm(rx: &ReactionData, idx: &[i32], count: usize, y: &[f64], e: &[f64]) -> f64 {
    if count == 0 {
        return 0.0;
    }
    let mut s = 0.0;
    for k in 0..count {
        let us = idx[k] as usize;
        let atol = if rx.species_atol[us] > 0.0 { rx.species_atol[us] } else { rx.atol };
        let rtol = if rx.species_rtol[us] > 0.0 { rx.species_rtol[us] } else { rx.rtol };
        let w = atol + rtol * y[k].abs();
        let r = e[k] / w;
        s += r * r;
    }
    (s / count as f64).sqrt()
}

fn finite_all(v: &[f64]) -> bool {
    v.iter().all(|x| x.is_finite())
}

/// RxStepReport (the fields the legacy binding reads).
struct StepReport {
    ok: bool,
    error: Bytes,
}

/// The integrator's per-call workspace (RxWorkspace, the members the RATE path uses).
struct Ws {
    term_needed: Vec<u8>,
    terms: Vec<f64>,
    y: Vec<f64>,
    y0: Vec<f64>,
    ytmp: Vec<f64>,
    rates: Vec<f64>,
    err: Vec<f64>,
    k: [Vec<f64>; 6],
    grp_out: Vec<f64>,
}

impl Ws {
    /// RxWorkspace::init
    fn new(rx: &ReactionData) -> Ws {
        let un = rx.n_species() as usize;
        Ws {
            term_needed: vec![0; rx.term_expr.len()],
            terms: vec![0.0; rx.term_expr.len()],
            y: vec![0.0; un],
            y0: vec![0.0; un],
            ytmp: vec![0.0; un],
            rates: vec![0.0; un],
            err: vec![0.0; un],
            k: [vec![0.0; un], vec![0.0; un], vec![0.0; un], vec![0.0; un], vec![0.0; un], vec![0.0; un]],
            grp_out: vec![0.0; un],
        }
    }
}

/// What the RATE right-hand side (the `fg` lambda) reads.
struct Rhs<'a> {
    rx: &'a ReactionData,
    spans: &'a [RxExprSpan],
    pollutants: Option<&'a [f64]>,
    coefficients: &'a [f64],
    ridx: &'a [i32],
    factor: f64,
    kinetic_scale: Option<&'a [f64]>,
    coupled: bool,
    gb: usize,
    gn: usize,
}

impl Rhs<'_> {
    /// The `fg` lambda: install y into the group's species (others frozen at entry when not
    /// coupled), evaluate every RATE expression, return the group's rates in `out`.
    fn fg(&self, y: &[f64], out: &mut [f64], species: &mut [f64], hydvar: &[f64], ws_y0: &[f64], terms: &mut [f64], needed: &[u8], rates: &mut [f64]) {
        let ytmp: Vec<f64> = y[..self.gn].to_vec();
        for k in 0..self.gn {
            species[self.ridx[self.gb + k] as usize] = ytmp[k];
        }
        // freeze_others(coupled ? -1 : gb)
        if !self.coupled {
            let active = self.gb;
            for k in 0..self.ridx.len() {
                if k != active {
                    species[self.ridx[k] as usize] = ws_y0[self.ridx[k] as usize];
                }
            }
        }
        // evalRates (one context)
        eval_terms(self.rx, species, hydvar, self.pollutants, terms, needed, self.coefficients);
        {
            let env = Env { species, coefs: self.coefficients, terms, hydvar, pollutants: self.pollutants };
            for k in 0..self.ridx.len() {
                let ks = match self.kinetic_scale {
                    Some(sc) => sc[self.ridx[k] as usize],
                    None => 1.,
                };
                rates[k] = self.factor * ks * eval_expr(&self.rx.token_pool, &self.spans[self.ridx[k] as usize], &env);
            }
        }
        for k in 0..self.gn {
            out[k] = rates[self.gb + k];
        }
    }
}

/// ReactionIntegrator::step (one context): RATE species by the explicit solvers (EUL, RK5),
/// then FORMULA species in declaration order.
#[allow(clippy::too_many_arguments)]
fn integrator_step(rx: &ReactionData, forms: &[ReactionExprForm], spans: &[RxExprSpan], dt: f64, species: &mut [f64], hydvar: &mut [f64], ws: &mut Ws, pollutants: Option<&[f64]>, kinetic_scale: Option<&[f64]>, coefficients: &[f64]) -> Result<StepReport, StepError> {
    let mut rep = StepReport { ok: true, error: Vec::new() };
    let n = rx.n_species();
    if !rx.compiled || n == 0 || dt <= 0.0 {
        return Ok(rep);
    }
    if forms.len() != n as usize || spans.len() != n as usize {
        rep.ok = false;
        rep.error = b"Reaction context shape differs from species".to_vec();
        return Ok(rep);
    }
    hydvar[RxHydVar::Dt as usize] = dt;
    let factor = unit_factor(rx.rate_units);
    let mut ridx: Vec<i32> = Vec::new();
    let mut eidx: Vec<i32> = Vec::new();
    let mut fidx: Vec<i32> = Vec::new();
    for s in 0..n {
        match forms[s as usize] {
            ReactionExprForm::Rate => ridx.push(s),
            ReactionExprForm::Equil => eidx.push(s),
            ReactionExprForm::Formula => fidx.push(s),
            _ => {}
        }
    }
    ws.term_needed.iter_mut().for_each(|x| *x = 0);
    {
        let mark = |e: &RxExprSpan, needed: &mut Vec<u8>| {
            for i in e.begin..e.begin + e.len {
                if rx.token_pool[i as usize].op == RxToken::PUSH_TERM {
                    needed[rx.token_pool[i as usize].idx as usize] = 1;
                }
            }
        };
        for s in 0..n as usize {
            if forms[s] != ReactionExprForm::None {
                mark(&spans[s], &mut ws.term_needed);
            }
        }
        let mut term = rx.term_expr.len() as i32 - 1;
        while term >= 0 {
            if ws.term_needed[term as usize] != 0 {
                mark(&rx.term_expr[term as usize], &mut ws.term_needed);
            }
            term -= 1;
        }
    }

    let nr = ridx.len();
    if nr > 0 {
        let coupled = rx.coupling == ReactionCoupling::Full || nr == 1;
        for s2 in 0..n as usize {
            ws.y0[s2] = species[s2];
        }
        for k in 0..nr {
            ws.grp_out[ridx[k] as usize] = ws.y0[ridx[k] as usize];
        }
        let groups = if coupled { 1 } else { nr };
        let mut g = 0;
        while g < groups && rep.ok {
            let gb = if coupled { 0 } else { g };
            let gn = if coupled { nr } else { 1 };
            let rhs = Rhs { rx, spans, pollutants, coefficients, ridx: &ridx, factor, kinetic_scale, coupled, gb, gn };
            let mut t = 0.0;
            let mut h = dt;
            let h_min = dt * K_MIN_STEP_FRAC;
            let mut guard = 0;
            for k in 0..gn {
                ws.y[k] = species[ridx[gb + k] as usize];
            }
            while t < dt && rep.ok {
                guard += 1;
                if guard > K_MAX_SUBSTEPS {
                    rep.ok = false;
                    rep.error = if rx.solver == ReactionSolverKind::Rk5 || rx.solver == ReactionSolverKind::Eul {
                        b"reaction step exceeded the substep cap: the explicit solver is resolving kinetics too stiff for it. Set [REACTION_OPTIONS] SOLVER to ROS2 or BDF2, or declare the fast species EQUIL.".to_vec()
                    } else {
                        b"reaction step exceeded the substep cap (stiffness beyond the solver \xe2\x80\x94 see the D-R7 escape hatch)".to_vec()
                    };
                    break;
                }
                if h < h_min {
                    rep.ok = false;
                    rep.error = b"reaction step size collapsed below the floor (stiffness beyond the solver \xe2\x80\x94 see D-R7 escape hatch)".to_vec();
                    break;
                }
                if t + h > dt {
                    h = dt - t;
                }
                let y_now = ws.y.clone();
                let mut k1 = std::mem::take(&mut ws.k[0]);
                rhs.fg(&y_now, &mut k1, species, hydvar, &ws.y0, &mut ws.terms, &ws.term_needed, &mut ws.rates);
                ws.k[0] = k1;
                if !finite_all(&ws.k[0][..gn]) {
                    rep.ok = false;
                    rep.error = b"non-finite reaction rate at accepted state".to_vec();
                    break;
                }
                let mut accepted = false;
                let mut h_next = h;
                match rx.solver {
                    ReactionSolverKind::Eul => {
                        for k in 0..gn {
                            ws.y[k] += h * ws.k[0][k];
                        }
                        accepted = true;
                    }
                    ReactionSolverKind::Rk5 => {
                        const B21: f64 = 0.2;
                        const B31: f64 = 3.0 / 40.0;
                        const B32: f64 = 9.0 / 40.0;
                        const B41: f64 = 0.3;
                        const B42: f64 = -0.9;
                        const B43: f64 = 1.2;
                        const B51: f64 = -11.0 / 54.0;
                        const B52: f64 = 2.5;
                        const B53: f64 = -70.0 / 27.0;
                        const B54: f64 = 35.0 / 27.0;
                        const B61: f64 = 1631.0 / 55296.0;
                        const B62: f64 = 175.0 / 512.0;
                        const B63: f64 = 575.0 / 13824.0;
                        const B64: f64 = 44275.0 / 110592.0;
                        const B65: f64 = 253.0 / 4096.0;
                        const C1: f64 = 37.0 / 378.0;
                        const C3: f64 = 250.0 / 621.0;
                        const C4: f64 = 125.0 / 594.0;
                        const C6: f64 = 512.0 / 1771.0;
                        const D1: f64 = 2825.0 / 27648.0;
                        const D3: f64 = 18575.0 / 48384.0;
                        const D4: f64 = 13525.0 / 55296.0;
                        const D5: f64 = 277.0 / 14336.0;
                        const D6: f64 = 0.25;
                        let y = ws.y.clone();
                        let stage = |ws: &mut Ws, kk: usize, species: &mut [f64], hydvar: &mut [f64]| {
                            let yt = ws.ytmp.clone();
                            let mut out = std::mem::take(&mut ws.k[kk]);
                            rhs.fg(&yt, &mut out, species, hydvar, &ws.y0, &mut ws.terms, &ws.term_needed, &mut ws.rates);
                            ws.k[kk] = out;
                        };
                        for k in 0..gn {
                            ws.ytmp[k] = y[k] + h * B21 * ws.k[0][k];
                        }
                        stage(ws, 1, species, hydvar);
                        for k in 0..gn {
                            ws.ytmp[k] = y[k] + h * (B31 * ws.k[0][k] + B32 * ws.k[1][k]);
                        }
                        stage(ws, 2, species, hydvar);
                        for k in 0..gn {
                            ws.ytmp[k] = y[k] + h * (B41 * ws.k[0][k] + B42 * ws.k[1][k] + B43 * ws.k[2][k]);
                        }
                        stage(ws, 3, species, hydvar);
                        for k in 0..gn {
                            ws.ytmp[k] = y[k] + h * (B51 * ws.k[0][k] + B52 * ws.k[1][k] + B53 * ws.k[2][k] + B54 * ws.k[3][k]);
                        }
                        stage(ws, 4, species, hydvar);
                        for k in 0..gn {
                            ws.ytmp[k] = y[k] + h * (B61 * ws.k[0][k] + B62 * ws.k[1][k] + B63 * ws.k[2][k] + B64 * ws.k[3][k] + B65 * ws.k[4][k]);
                        }
                        stage(ws, 5, species, hydvar);
                        for k in 0..gn {
                            let y5 = y[k] + h * (C1 * ws.k[0][k] + C3 * ws.k[2][k] + C4 * ws.k[3][k] + C6 * ws.k[5][k]);
                            let y4 = y[k] + h * (D1 * ws.k[0][k] + D3 * ws.k[2][k] + D4 * ws.k[3][k] + D5 * ws.k[4][k] + D6 * ws.k[5][k]);
                            ws.err[k] = y5 - y4;
                            ws.ytmp[k] = y5;
                        }
                        let en = err_norm(rx, &ridx[gb..], gn, &ws.ytmp, &ws.err);
                        if en <= 1.0 {
                            for k in 0..gn {
                                ws.y[k] = ws.ytmp[k];
                            }
                            accepted = true;
                            if en > 1.0e-30 {
                                h_next = cmin(dt, 0.9 * h * en.powf(-0.2)); // LIBM
                            }
                        } else {
                            h *= cmax(0.1, 0.9 * en.powf(-0.25)); // LIBM
                        }
                    }
                    _ => return Err(StepError::NotTranslated("reaction solvers ROS2 and BDF2")),
                }
                if accepted {
                    t += h;
                    h = h_next;
                } else if !finite_all(&ws.y[..gn]) {
                    rep.ok = false;
                    rep.error = b"non-finite species state during reaction integration".to_vec();
                }
            }
            if rep.ok {
                for k in 0..gn {
                    ws.grp_out[ridx[gb + k] as usize] = ws.y[k];
                }
            }
            g += 1;
        }
        if rep.ok {
            for k in 0..nr {
                species[ridx[k] as usize] = ws.grp_out[ridx[k] as usize];
            }
        }
    }
    if !rep.ok {
        return Ok(rep);
    }

    if !eidx.is_empty() {
        return Err(StepError::NotTranslated("reaction EQUIL species"));
    }

    // ---- 3. FORMULA: declaration order, current values
    for &s in fidx.iter() {
        eval_terms(rx, species, hydvar, pollutants, &mut ws.terms, &ws.term_needed, coefficients);
        let v = {
            let env = Env { species, coefs: coefficients, terms: &ws.terms, hydvar, pollutants };
            eval_expr(&rx.token_pool, &spans[s as usize], &env)
        };
        species[s as usize] = v;
    }

    if !finite_all(&species[..n as usize]) {
        rep.ok = false;
        rep.error = b"non-finite species state after reaction step".to_vec();
    }
    Ok(rep)
}

// ---------------------------------------------------------------------------
// ReactionGeometry.hpp
// ---------------------------------------------------------------------------

fn reaction_area_factor(units: ReactionAreaUnits) -> f64 {
    if units == ReactionAreaUnits::M2 {
        0.09290304
    } else if units == ReactionAreaUnits::Cm2 {
        929.0304
    } else {
        1.0
    }
}

fn reaction_wall_area(c: &SimulationContext, link_element: bool, e: i32, reference: bool) -> f64 {
    for row in c.reactions.wall_areas.iter() {
        if row.link == link_element && row.element == e {
            return row.area;
        }
    }
    if !link_element || e < 0 || e >= c.n_links() {
        return 0.0;
    }
    let r = c.link_subtypes.conduit_row(e);
    if r < 0 {
        return 0.0;
    }
    let ue = e as usize;
    let y = if reference { c.links.xsect_y_full[ue] } else { c.links.depth[ue] };
    if !(y > 0.0) {
        return 0.0;
    }
    let xs = link::build_xsect_params(&c.links, ue, Some(&c.transect_tables));
    let area = xsect::get_a_of_y(&xs, y);
    let radius = xsect::get_r_of_y(&xs, y);
    if !(radius > 0.0) {
        return 0.0;
    }
    let conduits = &c.link_subtypes.conduits;
    let ur = r as usize;
    let length = if ur < conduits.true_length.len() && conduits.true_length[ur] > 0.0 { conduits.true_length[ur] } else { conduits.length[ur] };
    area / radius * length * (cmax(1, conduits.barrels[ur]) as f64) * reaction_area_factor(c.reactions.area_units)
}

fn reaction_area_volume(c: &SimulationContext, link_element: bool, e: i32, volume: f64) -> f64 {
    if volume > 0.0353147 {
        reaction_wall_area(c, link_element, e, false) / (volume * 28.317)
    } else {
        0.0
    }
}

fn reaction_wall_scales(c: &SimulationContext, link_element: bool, e: i32, volume: f64, scales: &mut Vec<f64>) {
    scales.clear();
    scales.resize(c.reactions.n_species() as usize, 1.0);
    let reference = reaction_wall_area(c, link_element, e, true);
    let active = if volume > 0.0353147 { reaction_wall_area(c, link_element, e, false) } else { 0.0 };
    for s in 0..c.reactions.n_species() as usize {
        if c.reactions.species_is_wall[s] != 0 {
            scales[s] = if reference > 0.0 { cclamp(active / reference, 0.0, 1.0) } else { 0.0 };
        }
    }
}

// ---------------------------------------------------------------------------
// ReactionLegacyBinding.cpp
// ---------------------------------------------------------------------------

fn note_negative_clamp(ctx: &mut SimulationContext, where_: &[u8]) {
    let rx = &mut ctx.reactions;
    if rx.warned_negative_clamp {
        return;
    }
    rx.warned_negative_clamp = true;
    let mut w = b"Reaction step produced a negative concentration at ".to_vec();
    w.extend_from_slice(where_);
    w.extend_from_slice(b"; it was set to zero (mass added). The explicit solver is too coarse for these kinetics: use [REACTION_OPTIONS] SOLVER ROS2, COUPLING FULL or a shorter step.");
    ctx.warnings.push(w);
}

/// clampRateSpecies: zero negative RATE species; returns (and lists) what was added.
fn clamp_rate_species(rx: &ReactionData, tank: bool, block: &mut [f64], added: Option<&mut [f64]>) -> f64 {
    let forms = if tank { &rx.tank_form } else { &rx.pipe_form };
    let mut total = 0.0;
    let ns = rx.n_species();
    let mut added = added;
    for s in 0..ns as usize {
        let mut a = 0.0;
        if s < forms.len() && forms[s] == ReactionExprForm::Rate && block[s] < 0.0 {
            a = -block[s];
            block[s] = 0.0;
        }
        if let Some(ad) = added.as_deref_mut() {
            ad[s] = a;
        }
        total += a;
    }
    total
}

/// decayPollutantsExact: exact first-order pollutant decay on the elements (nodes skip
/// non-storage pools under a litre).
fn decay_pollutants_exact(ctx: &mut SimulationContext, dt: f64, tank: bool) {
    let np = ctx.n_pollutants();
    let n_elems = if tank { ctx.n_nodes() } else { ctx.n_links() };
    for p in 0..np {
        let k = ctx.pollutants.k_decay[p as usize];
        if k == 0.0 {
            continue;
        }
        let f = (-k * dt).exp(); // LIBM
        let mut removed = 0.0;
        for i in 0..n_elems {
            let idx = i as usize * np as usize + p as usize;
            let (conc, volume) = if tank { (&mut ctx.nodes.conc, &ctx.nodes.volume) } else { (&mut ctx.links.conc, &ctx.links.volume) };
            if idx >= conc.len() {
                continue;
            }
            let ui = i as usize;
            let v = if ui < volume.len() { volume[ui] } else { 0.0 };
            if tank && ui < ctx.nodes.type_.len() && ctx.nodes.type_[ui] != NodeType::Storage && v <= 0.0353147 {
                continue;
            }
            removed += conc[idx] * (1.0 - f) * v;
            conc[idx] *= f;
            if conc[idx] < 0.0 {
                conc[idx] = 0.0;
            }
        }
        if (p as usize) < ctx.mass_balance.qual_routing_reacted.len() {
            ctx.mass_balance.qual_routing_reacted[p as usize] += removed;
        }
    }
}

/// reactElements: each element's species through the integrator, with the routing ledger.
fn react_elements(ctx: &mut SimulationContext, dt: f64, tank: bool) -> Result<(), StepError> {
    let ns = ctx.reactions.n_species();
    if ns == 0 {
        return Ok(());
    }
    let uns = ns as usize;
    let mut ws = Ws::new(&ctx.reactions);
    let mut block = vec![0.0f64; uns];
    let mut added = vec![0.0f64; uns];
    let mut coef_scratch: Vec<f64> = Vec::new();
    let mut wall_scales: Vec<f64> = Vec::new();

    let np = ctx.n_pollutants();
    let mut hydvar = [0.0f64; RxHydVar::Count as usize];
    let n_elems = if tank { ctx.n_nodes() } else { ctx.n_links() };
    let what: &[u8] = if tank { b"node" } else { b"link" };

    for e in 0..n_elems {
        let ue = e as usize;
        let base = ue * uns;
        let vol_len = if tank { ctx.nodes.volume.len() } else { ctx.links.volume.len() };
        let vol_e = if ue < vol_len { if tank { ctx.nodes.volume[ue] } else { ctx.links.volume[ue] } } else { 0.0 };
        if tank && (!node_is_reactor(ctx, e) || ue >= vol_len || vol_e <= ZERO_VOLUME) {
            continue;
        }
        if ctx.reactions.wall_geometry && (ue >= vol_len || vol_e <= ZERO_VOLUME) {
            continue;
        }
        if tank && lidnode::active(ctx, e) {
            let row = ctx.node_subtypes.storage_row(e);
            if row >= 0 && ctx.reactions.lid_enabled(ctx.node_subtypes.storages.lid[row as usize].control) {
                continue;
            }
        }

        hydvar[RxHydVar::Hrt as usize] = if tank && ue < ctx.nodes.hrt.len() { ctx.nodes.hrt[ue] } else { 0.0 };
        let elem_temp = if tank { &ctx.heat_state.node_temp } else { &ctx.heat_state.link_temp };
        hydvar[RxHydVar::Temp as usize] = if ctx.options.heat_transport && ue < elem_temp.len() { elem_temp[ue] } else { ctx.reactions.default_temp_c };
        hydvar[RxHydVar::Av as usize] = if ctx.reactions.wall_geometry { reaction_area_volume(ctx, !tank, e, vol_e) } else { 0.0 };
        if ctx.reactions.wall_geometry {
            reaction_wall_scales(ctx, !tank, e, vol_e, &mut wall_scales);
        }

        {
            let msx = if tank { &ctx.reactions.msx_node_conc } else { &ctx.reactions.msx_link_conc };
            for s in 0..uns {
                block[s] = msx[base + s];
            }
        }

        let rep = {
            let rx = &ctx.reactions;
            let pollut_conc = if tank { &ctx.nodes.conc } else { &ctx.links.conc };
            let pollut: Option<&[f64]> = if np > 0 { Some(&pollut_conc[ue * np as usize..]) } else { None };
            let coefs: Vec<f64> = rx.coefficients(!tank, e, &mut coef_scratch).to_vec();
            let ks: Option<&[f64]> = if rx.wall_geometry { Some(&wall_scales) } else { None };
            let (forms, spans) = if tank { (&rx.tank_form, &rx.tank_expr) } else { (&rx.pipe_form, &rx.pipe_expr) };
            integrator_step(rx, forms, spans, dt, &mut block, &mut hydvar, &mut ws, pollut, ks, &coefs)?
        };

        if rep.ok {
            let v = vol_e;
            let book = ctx.reactions.msx_ledger.armed && ctx.reactions.msx_ledger.reacted.len() == uns;
            if book {
                for s in 0..uns {
                    if ctx.reactions.species_is_wall[s] == 0 {
                        let old = if tank { ctx.reactions.msx_node_conc[base + s] } else { ctx.reactions.msx_link_conc[base + s] };
                        ctx.reactions.msx_ledger.reacted[s] += (old - block[s]) * v;
                    }
                }
            }
            if clamp_rate_species(&ctx.reactions, tank, &mut block, Some(&mut added)) > 0.0 {
                note_negative_clamp(ctx, what);
                if book {
                    for s in 0..uns {
                        if ctx.reactions.species_is_wall[s] == 0 {
                            ctx.reactions.msx_ledger.clamped[s] += added[s] * v;
                        }
                    }
                }
            }
            let msx = if tank { &mut ctx.reactions.msx_node_conc } else { &mut ctx.reactions.msx_link_conc };
            for s in 0..uns {
                msx[base + s] = block[s];
            }
        } else if !ctx.reactions.warned_react_failure {
            ctx.reactions.warned_react_failure = true;
            let mut w = b"Reaction step failed at ".to_vec();
            w.extend_from_slice(what);
            w.extend_from_slice(format!(" index {}", e).as_bytes());
            w.extend_from_slice(b" (element state left unchanged): ");
            w.extend_from_slice(&rep.error);
            ctx.warnings.push(w);
        }
    }
    Ok(())
}

/// transport::reactLegacyNodes
pub fn react_legacy_nodes(ctx: &mut SimulationContext, dt: f64) -> Result<(), StepError> {
    if !legacy_reactions_active(ctx) || dt <= 0.0 {
        return Ok(());
    }
    ensure_msx_state(ctx);
    decay_pollutants_exact(ctx, dt, true);
    react_elements(ctx, dt, true)
}

/// transport::reactLegacyLinks
pub fn react_legacy_links(ctx: &mut SimulationContext, dt: f64) -> Result<(), StepError> {
    if !legacy_reactions_active(ctx) || dt <= 0.0 {
        return Ok(());
    }
    ensure_msx_state(ctx);
    decay_pollutants_exact(ctx, dt, false);
    react_elements(ctx, dt, false)
}

// ---------------------------------------------------------------------------
// MsxLidReactions.cpp: the read helpers routeLegacyMsx calls
// ---------------------------------------------------------------------------

fn enabled_node(c: &SimulationContext, row: i32) -> bool {
    let s = &c.node_subtypes.storages;
    row >= 0 && row < s.count() && !s.lid_state[row as usize].cells.is_empty() && c.reactions.lid_enabled(s.lid[row as usize].control)
}

/// msxlid::heldMass
pub fn msxlid_held_mass(c: &SimulationContext, node: i32, species: i32) -> f64 {
    let row = c.node_subtypes.storage_row(node);
    let ns = c.reactions.n_species();
    if row < 0 || row >= c.reactions.lid.node_mass.len() as i32 || species < 0 || species >= ns {
        return 0.0;
    }
    let mut total = 0.0;
    let mass = &c.reactions.lid.node_mass[row as usize];
    let mut i = species as usize;
    while i < mass.len() {
        total += mass[i];
        i += ns as usize;
    }
    total
}

fn msxlid_mobile_water(c: &SimulationContext, node: i32, fallback: f64) -> f64 {
    let row = c.node_subtypes.storage_row(node);
    if enabled_node(c, row) && row < c.reactions.lid.node_mobile_old.len() as i32 {
        c.reactions.lid.node_mobile_old[row as usize]
    } else {
        fallback
    }
}

fn msxlid_outlet(c: &SimulationContext, node: i32, link: i32, species: i32, fallback: f64) -> f64 {
    let row = c.node_subtypes.storage_row(node);
    let ns = c.reactions.n_species();
    if row < 0 || row >= c.reactions.lid.node_outlet.len() as i32 {
        return fallback;
    }
    let out = &c.reactions.lid.node_outlet[row as usize];
    let k = link as usize * ns as usize + species as usize;
    if k < out.len() && out[k].is_finite() {
        out[k]
    } else {
        fallback
    }
}

fn msxlid_receive(c: &mut SimulationContext, node: i32, link: i32, water: f64, species: i32, rate: f64, dt: f64) -> f64 {
    let row = c.node_subtypes.storage_row(node);
    let ns = c.reactions.n_species();
    if !enabled_node(c, row) || !(water > 0.0) {
        return 0.0;
    }
    let ports = c.node_subtypes.storages.lid_state[row as usize].quality_ports.clone();
    for p in ports.iter() {
        if p.link == link && p.volume > 0.0 {
            let sent = cmax(0., rate) * cclamp(p.volume / water, 0., 1.);
            c.reactions.lid.node_mass[row as usize][(p.cell * ns + species) as usize] += sent * dt;
            return sent;
        }
    }
    0.0
}

// ---------------------------------------------------------------------------
// MsxLegacyTransport.cpp
// ---------------------------------------------------------------------------

/// reactionBoosterMass: the booster sources' added mass (none without a booster row).
fn reaction_booster_mass(ctx: &SimulationContext, node: i32, species: i32, mixed_volume: f64, arriving_volume: f64) -> Result<f64, StepError> {
    if !(mixed_volume > 1e-12) || !(arriving_volume > 1e-12) {
        return Ok(0.0);
    }
    for source in ctx.reactions.sources.iter() {
        if source.booster != 0 && source.node == node && source.species == species {
            return Err(StepError::NotTranslated("reaction booster sources"));
        }
    }
    Ok(0.0)
}

/// armMsxRoutingLedger: zero and arm the routing ledger once, booking the opening inventory.
/// `fv_retained`: the FV transport's retained species mass (state_accessors).
fn arm_msx_routing_ledger(ctx: &mut SimulationContext, fv_retained: &[f64]) {
    ensure_msx_state(ctx);
    let nn = ctx.n_nodes();
    let nl = ctx.n_links();
    let uns = ctx.reactions.n_species() as usize;
    if !ctx.reactions.msx_ledger.armed || ctx.reactions.msx_ledger.init.len() != uns {
        ctx.reactions.msx_ledger.arm(uns);
        {
            let held = fv_retained;
            let rx = &mut ctx.reactions;
            let mut s = 0;
            while s < uns && s < held.len() {
                if rx.species_is_wall[s] == 0 {
                    rx.msx_ledger.init[s] += held[s];
                    rx.msx_ledger.final_dry[s] = held[s];
                }
                s += 1;
            }
        }
        for i in 0..nn {
            for s in 0..uns {
                if ctx.reactions.species_is_wall[s] == 0 {
                    let add = ctx.reactions.msx_node_conc[i as usize * uns + s] * ctx.nodes.old_volume[i as usize] + msxlid_held_mass(ctx, i, s as i32);
                    ctx.reactions.msx_ledger.init[s] += add;
                }
            }
        }
        for j in 0..nl as usize {
            for s in 0..uns {
                if ctx.reactions.species_is_wall[s] == 0 {
                    ctx.reactions.msx_ledger.init[s] += ctx.reactions.msx_link_conc[j * uns + s] * ctx.links.old_volume[j];
                }
            }
        }
    }
}

/// transport::routeLegacyMsx: species transport on the CSTR mirror (no reaction here).
pub fn route_legacy_msx(ctx: &mut SimulationContext, dt: f64, fv_retained: &[f64]) -> Result<(), StepError> {
    let nsp = ctx.reactions.n_species();
    if nsp <= 0 || dt <= 0.0 {
        return Ok(());
    }
    let nn = ctx.n_nodes();
    let nl = ctx.n_links();
    let uns = nsp as usize;
    ensure_msx_state(ctx);
    if ctx.reactions.msx_node_conc.len() < nn as usize * uns || ctx.reactions.msx_link_conc.len() < nl as usize * uns {
        return Ok(());
    }
    arm_msx_routing_ledger(ctx, fv_retained);
    let node_old: Vec<f64> = ctx.reactions.msx_node_conc[..nn as usize * uns].to_vec();
    let link_old: Vec<f64> = ctx.reactions.msx_link_conc[..nl as usize * uns].to_vec();
    let mut mass_in = vec![0.0f64; nn as usize * uns];

    for j in 0..nl {
        let uj = j as usize;
        let q = ctx.links.flow[uj].abs();
        if q <= 0.0 {
            continue;
        }
        let dn = if ctx.links.flow[uj] >= 0.0 { ctx.links.node2[uj] } else { ctx.links.node1[uj] };
        if dn < 0 || dn >= nn {
            continue;
        }
        let udn = dn as usize;
        let source = if ctx.links.flow[uj] >= 0.0 { ctx.links.node1[uj] } else { ctx.links.node2[uj] };
        for s in 0..uns {
            if ctx.reactions.species_is_wall[s] == 0 {
                let mut concentration = link_old[uj * uns + s];
                if link_takes_upstream_value(ctx, j) && (lidnode::active(ctx, source) || lidnode::active(ctx, dn)) {
                    concentration = msxlid_outlet(ctx, source, j, s as i32, node_old[source as usize * uns + s]);
                }
                let rate = q * concentration;
                let rcv = msxlid_receive(ctx, dn, j, q * dt, s as i32, rate, dt);
                mass_in[udn * uns + s] += rate - rcv;
            }
        }
    }

    let has_ext = ctx.reactions.msx_ext_mass_in.len() >= nn as usize * uns;
    for i in 0..nn {
        let ui = i as usize;
        let v_old = msxlid_mobile_water(ctx, i, ctx.nodes.old_volume[ui]);
        let v_in = ctx.nodes.qual_vol_in[ui];
        let q_in = if dt > 0.0 { v_in / dt } else { 0.0 };
        let zero_bf = ctx.options.outfall_backflow_zero && ctx.nodes.type_[ui] == NodeType::Outfall && v_in <= 0.0;
        let reactor = node_is_reactor(ctx, i);
        let dry = node_is_dry(ctx, i, q_in) && !lidnode::active(ctx, i);
        let f_evap = if reactor { storage_evap_factor(ctx, i) } else { 1.0 };
        let mut v_exfil = 0.0;
        if reactor {
            let srow = ctx.node_subtypes.storage_row(i);
            let st = &ctx.node_subtypes.storages;
            if srow >= 0 && (srow as usize) < st.exfil_loss.len() {
                v_exfil = st.exfil_loss[srow as usize];
            }
            if srow >= 0 && ctx.reactions.lid_enabled(st.lid[srow as usize].control) {
                v_exfil = cmax(0., v_exfil - st.lid_state[srow as usize].richards_bottom_loss);
            }
        }
        for s in 0..uns {
            if ctx.reactions.species_is_wall[s] != 0 {
                continue;
            }
            let c_old = node_old[ui * uns + s];
            let ext_rate = if has_ext { ctx.reactions.msx_ext_mass_in[ui * uns + s] } else { 0.0 };
            let mass = (mass_in[ui * uns + s] + ext_rate) * dt;
            ctx.reactions.msx_ledger.inflow[s] += ext_rate * dt;
            if !reactor {
                let mut c;
                if q_in > LEGACY_ZERO {
                    c = mass / v_in;
                } else if zero_bf {
                    c = 0.0;
                } else if ctx.nodes.depth[ui] > ZERO_DEPTH {
                    c = c_old;
                } else {
                    c = 0.0;
                }
                if q_in > LEGACY_ZERO {
                    c += reaction_booster_mass(ctx, i, s as i32, v_in, v_in)? / v_in;
                }
                ctx.reactions.msx_node_conc[ui * uns + s] = cmax(c, 0.0);
                continue;
            }
            ctx.reactions.msx_ledger.seepage[s] += v_exfil * c_old;
            let c1 = c_old * f_evap;
            let mut c_new;
            if q_in <= LEGACY_ZERO {
                c_new = if zero_bf { 0.0 } else { c1 };
            } else {
                let c_in = mass / v_in;
                let c_max = cmax(c1, c_in);
                c_new = (c1 * v_old + mass) / (v_old + v_in);
                c_new = cmin(c_new, c_max);
            }
            let mixed = v_old + v_in;
            if !dry && mixed > 1e-12 {
                c_new += reaction_booster_mass(ctx, i, s as i32, mixed, v_in)? / mixed;
            }
            if dry {
                ctx.reactions.msx_ledger.final_dry[s] += c_new * ctx.nodes.volume[ui];
                c_new = 0.0;
            }
            ctx.reactions.msx_node_conc[ui * uns + s] = cmax(c_new, 0.0);
        }
    }

    let is_steady = ctx.options.routing_model == RoutingModel::Steady;
    for j in 0..nl {
        let uj = j as usize;
        let up = if is_steady {
            ctx.links.node1[uj]
        } else if ctx.links.flow[uj] >= 0.0 {
            ctx.links.node1[uj]
        } else {
            ctx.links.node2[uj]
        };
        if up < 0 || up >= nn {
            continue;
        }
        let uup = up as usize;
        let passthrough = link_takes_upstream_value(ctx, j);
        let v_old = ctx.links.old_volume[uj];
        let q_in = if passthrough { 0.0 } else { conduit_mixing_inflow(ctx, j, dt) };
        let dry = !is_steady && !passthrough && link_is_dry(ctx, j);
        let mut q_seep = 0.0;
        let mut f_evap = 1.0;
        if !passthrough {
            let cr = ctx.link_subtypes.conduit_row(j);
            if cr >= 0 {
                let ucr = cr as usize;
                let cd = &ctx.link_subtypes.conduits;
                let barrels = cd.barrels[ucr] as f64;
                q_seep = cd.seep_loss_rate[ucr] * barrels;
                let v_evap = cd.evap_loss_rate[ucr] * barrels * dt;
                if v_evap > 0.0 && v_old > ZERO_VOLUME {
                    f_evap += v_evap / v_old;
                }
            }
        }
        for s in 0..uns {
            if ctx.reactions.species_is_wall[s] != 0 {
                continue;
            }
            let c_old = link_old[uj * uns + s];
            let c_up = msxlid_outlet(ctx, up, j, s as i32, ctx.reactions.msx_node_conc[uup * uns + s]);
            let mut c_new;
            if passthrough {
                c_new = c_up;
            } else if is_steady {
                ctx.reactions.msx_ledger.seepage[s] += q_seep * c_up * dt;
                c_new = c_up * f_evap;
            } else {
                ctx.reactions.msx_ledger.seepage[s] += q_seep * c_old * dt;
                let c1 = c_old * f_evap;
                if q_in <= LEGACY_ZERO {
                    c_new = c1;
                } else {
                    let v_in = q_in * dt;
                    let c_max = cmax(c1, c_up);
                    c_new = (c1 * v_old + c_up * v_in) / (v_old + v_in);
                    c_new = cmin(c_new, c_max);
                }
            }
            if dry {
                ctx.reactions.msx_ledger.final_dry[s] += c_new * ctx.links.volume[uj];
                c_new = 0.0;
            }
            ctx.reactions.msx_link_conc[uj * uns + s] = cmax(c_new, 0.0);
        }
    }
    Ok(())
}

/// FvQualityTransport::retainedSpeciesMass: per species, the held dry inventory row np + s
/// (zero for rows not held).
pub fn fv_retained_species_mass(held: &[Vec<f64>], n_pollutants: i32, n_species: i32) -> Vec<f64> {
    let mut mass = vec![0.0f64; if n_species > 0 { n_species as usize } else { 0 }];
    let np = n_pollutants as usize;
    let mut s = 0;
    while s < mass.len() && np + s < held.len() {
        mass[s] = held[np + s].iter().fold(0.0, |a, &x| a + x);
        s += 1;
    }
    mass
}
