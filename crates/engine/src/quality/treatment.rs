// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/quality/Treatment.cpp and Treatment.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// treatment::evaluate, the co-treatment overload (the one the routing step's treatment calls:
// QualitySolver::applyTreatment through applyNodeTreatment, and the node-LID layer treatment).
// The parser (parse, tokenize, shuntingYard) is in hydrology/lid_node.rs (module treatment).
//
// Differences in form, never in result:
//   * The std::stack<double> is a Vec<f64>; an empty stack at the end gives 0.0, as upstream.
//   * The optional arrays (cin, removal, cpollut: const double*, possibly null) are
//     Option<&[f64]>.

use crate::core::cmath::{cmax, cmin};
use crate::data::treatment_types::{TokenType, TreatExpr, TreatVar};

/// treatment::evaluate(expr, c, dt, hrt, q, v, d, cin, removal, n_pollut, area, cpollut): the
/// postfix stack machine with pollutant concentration (C_POLLUT) and removal (R_POLLUT)
/// references. STEP(0) is 0 here (strictly positive), unlike the scalar overload.
#[allow(clippy::too_many_arguments)]
pub fn evaluate(
    expr: &TreatExpr,
    c: f64,
    dt: f64,
    hrt: f64,
    q: f64,
    v: f64,
    d: f64,
    cin: Option<&[f64]>,
    removal: Option<&[f64]>,
    n_pollut: i32,
    area: f64,
    cpollut: Option<&[f64]>,
) -> f64 {
    let mut stk: Vec<f64> = Vec::new();

    for tok in expr.tokens.iter() {
        match tok.type_ {
            TokenType::Number => stk.push(tok.value),

            TokenType::Variable => match tok.var {
                TreatVar::C => stk.push(c),
                TreatVar::R => stk.push(0.0),
                TreatVar::Dt => stk.push(dt),
                TreatVar::Hrt => stk.push(hrt),
                TreatVar::Q => stk.push(q),
                TreatVar::V => stk.push(v),
                TreatVar::D => stk.push(d),
                TreatVar::Area => stk.push(area), // Gap #16
                TreatVar::CPollut => {
                    if tok.pollut_ref >= 0 && tok.pollut_ref < n_pollut {
                        if let Some(cp) = cpollut {
                            stk.push(cp[tok.pollut_ref as usize]);
                        } else if let Some(ci) = cin {
                            stk.push(ci[tok.pollut_ref as usize]);
                        } else {
                            stk.push(0.0);
                        }
                    } else {
                        stk.push(0.0);
                    }
                }
                TreatVar::RPollut => {
                    if let (Some(r), true) = (removal, tok.pollut_ref >= 0 && tok.pollut_ref < n_pollut) {
                        stk.push(cmax(0.0, r[tok.pollut_ref as usize]));
                    } else {
                        stk.push(0.0);
                    }
                }
            },

            TokenType::Neg => {
                let a = match stk.pop() {
                    Some(a) => a,
                    None => return 0.0,
                };
                stk.push(-a);
            }

            TokenType::Add | TokenType::Sub | TokenType::Mul | TokenType::Div | TokenType::Pow => {
                if stk.len() < 2 {
                    return 0.0;
                }
                let b = stk.pop().unwrap();
                let a = stk.pop().unwrap();
                match tok.type_ {
                    TokenType::Add => stk.push(a + b),
                    TokenType::Sub => stk.push(a - b),
                    TokenType::Mul => stk.push(a * b),
                    TokenType::Div => stk.push(if b != 0.0 { a / b } else { 0.0 }),
                    TokenType::Pow => stk.push(a.powf(b)), // LIBM
                    _ => {}
                }
            }

            TokenType::FuncExp | TokenType::FuncLog | TokenType::FuncSqrt | TokenType::FuncAbs | TokenType::FuncSgn | TokenType::FuncStep => {
                let a = match stk.pop() {
                    Some(a) => a,
                    None => return 0.0,
                };
                match tok.type_ {
                    TokenType::FuncExp => stk.push(a.exp()), // LIBM
                    TokenType::FuncLog => stk.push(if a > 0.0 { a.ln() } else { 0.0 }), // LIBM
                    TokenType::FuncSqrt => stk.push(if a >= 0.0 { a.sqrt() } else { 0.0 }),
                    TokenType::FuncAbs => stk.push(a.abs()),
                    TokenType::FuncSgn => stk.push(if a > 0.0 {
                        1.0
                    } else if a < 0.0 {
                        -1.0
                    } else {
                        0.0
                    }),
                    TokenType::FuncStep => stk.push(if a > 0.0 { 1.0 } else { 0.0 }),
                    _ => {}
                }
            }

            TokenType::FuncMin | TokenType::FuncMax => {
                if stk.len() < 2 {
                    return 0.0;
                }
                let b = stk.pop().unwrap();
                let a = stk.pop().unwrap();
                stk.push(if tok.type_ == TokenType::FuncMin { cmin(a, b) } else { cmax(a, b) });
            }

            _ => {}
        }
    }

    match stk.last() {
        Some(&x) => x,
        None => 0.0,
    }
}
