// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/components/ReactionModule/
// ReactionExpression.hpp and ReactionExpression.cpp (HydroCouple/openswmm.engine, branch
// swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache
// License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The reaction-expression compiler (lexer + shunting-yard into the flat RxToken pool),
// the per-context availability diagnostic and the vocabulary tables. The Tier-1 evaluator
// (evalReactionExpression) runs during quality routing and is not translated here.
//
// Differences of form only:
//   * RxSymbols holds slices instead of nullable vector pointers (None = null).
//   * A number is read with glibc strtod from the lexeme's position to the end of the C
//     string (core::cstdlib::strtod), as `std::strtod(s.c_str() + i, &end)` reads it.
//   * isdigit/isalpha/isalnum/toupper are the C locale's: ASCII only.

use crate::core::cstdlib::strtod;
use crate::core::options::Bytes;
use crate::data::reaction::ReactionData;
use crate::data::reaction_tokens::{RxExprSpan, RxHydVar, RxToken, K_RX_MAX_STACK_DEPTH};

pub struct RxSymbols<'a> {
    /// local slots
    pub species: Option<&'a [Bytes]>,
    pub coefs: Option<&'a [Bytes]>,
    /// only [0, max_term) resolvable
    pub terms: Option<&'a [Bytes]>,
    pub pollutants: Option<&'a [Bytes]>,
    /// forward-only rule: terms below this index only
    pub max_term: i32,
}

// ---------------------------------------------------------------------------
// Lexer
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LexKind {
    Number,
    Ident,
    Op,
    LParen,
    RParen,
    Comma,
    End,
}

struct Lexeme {
    kind: LexKind,
    /// for OP: + - * / ^
    op: u8,
    num: f64,
    ident: Bytes,
    /// 1-based start column
    col: i32,
}

struct Lexer<'a> {
    s: &'a [u8],
    i: usize,
}

impl<'a> Lexer<'a> {
    fn new(src: &'a [u8]) -> Self {
        Lexer { s: src, i: 0 }
    }

    fn next(&mut self, err: &mut Bytes) -> Lexeme {
        let s = self.s;
        while self.i < s.len() && (s[self.i] == b' ' || s[self.i] == b'\t') {
            self.i += 1;
        }
        let mut lx = Lexeme { kind: LexKind::End, op: 0, num: 0.0, ident: Vec::new(), col: self.i as i32 + 1 };
        if self.i >= s.len() {
            lx.kind = LexKind::End;
            return lx;
        }
        let c = s[self.i];
        if c.is_ascii_digit() || c == b'.' {
            let (v, n) = strtod(&s[self.i..]);
            lx.num = v;
            lx.kind = LexKind::Number;
            self.i += n;
            return lx;
        }
        if c.is_ascii_alphabetic() || c == b'_' {
            let b = self.i;
            while self.i < s.len() && (s[self.i].is_ascii_alphanumeric() || s[self.i] == b'_') {
                self.i += 1;
            }
            lx.kind = LexKind::Ident;
            lx.ident = s[b..self.i].to_vec();
            return lx;
        }
        match c {
            b'+' | b'-' | b'*' | b'/' | b'^' => {
                lx.kind = LexKind::Op;
                lx.op = c;
                self.i += 1;
                lx
            }
            b'(' => {
                lx.kind = LexKind::LParen;
                self.i += 1;
                lx
            }
            b')' => {
                lx.kind = LexKind::RParen;
                self.i += 1;
                lx
            }
            b',' => {
                lx.kind = LexKind::Comma;
                self.i += 1;
                lx
            }
            _ => {
                let mut m = b"unexpected character '".to_vec();
                m.push(c);
                m.push(b'\'');
                *err = m;
                lx
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Function table
// ---------------------------------------------------------------------------

struct FuncDef {
    name: &'static str,
    op: u8,
    arity: i32,
}

const K_FUNCS: [FuncDef; 13] = [
    FuncDef { name: "EXP", op: RxToken::F_EXP, arity: 1 },
    FuncDef { name: "LOG", op: RxToken::F_LOG, arity: 1 },
    FuncDef { name: "LOG10", op: RxToken::F_LOG10, arity: 1 },
    FuncDef { name: "SQRT", op: RxToken::F_SQRT, arity: 1 },
    FuncDef { name: "ABS", op: RxToken::F_ABS, arity: 1 },
    FuncDef { name: "SGN", op: RxToken::F_SGN, arity: 1 },
    FuncDef { name: "STEP", op: RxToken::F_STEP, arity: 1 },
    FuncDef { name: "SIN", op: RxToken::F_SIN, arity: 1 },
    FuncDef { name: "COS", op: RxToken::F_COS, arity: 1 },
    FuncDef { name: "TAN", op: RxToken::F_TAN, arity: 1 },
    FuncDef { name: "MIN", op: RxToken::F_MIN, arity: 2 },
    FuncDef { name: "MAX", op: RxToken::F_MAX, arity: 2 },
    FuncDef { name: "POW", op: RxToken::POW, arity: 2 },
];

fn upper(t: &[u8]) -> Bytes {
    t.iter().map(|c| c.to_ascii_uppercase()).collect()
}

fn find_func(ident: &[u8]) -> Option<&'static FuncDef> {
    let u = upper(ident);
    K_FUNCS.iter().find(|f| u == f.name.as_bytes())
}

const K_HYD_VAR_NAMES: [&str; 14] = ["D", "Q", "U", "RE", "US", "FF", "AV", "HRT", "DT", "TEMP", "GW_THETA", "GW_POROSITY", "GW_SE", "GW_DEPTH_M"];
const K_HYD_VAR_DESCS: [&str; 14] = [
    "Depth (ft); not yet populated, refused",
    "Flow (cfs); not yet populated, refused",
    "Velocity (ft/s); not yet populated, refused",
    "Reynolds number; not yet populated, refused",
    "Shear velocity (ft/s); not yet populated, refused",
    "Darcy-Weisbach friction factor; not yet populated, refused",
    "Wetted contact area per litre (AREA_UNITS/L; WALL_GEOMETRY YES)",
    "Hydraulic residence time (s)",
    "Reaction step (s)",
    "Water temperature (degC; TEMPERATURE option when heat is off)",
    "Groundwater water content (dimensionless; groundwater only)",
    "Groundwater porosity (dimensionless; groundwater only)",
    "Groundwater effective saturation (dimensionless; groundwater only)",
    "Groundwater reactor centre depth (m; groundwater only)",
];

fn find_hyd_var(ident: &[u8]) -> i32 {
    let u = upper(ident);
    for i in 0..RxHydVar::Count as i32 {
        if u == K_HYD_VAR_NAMES[i as usize].as_bytes() {
            return i;
        }
    }
    -1
}

// ---------------------------------------------------------------------------
// Shunting-yard
// ---------------------------------------------------------------------------

/// OpEntry{} value-initialized: every member zero/false/null.
#[derive(Clone, Copy, Default)]
#[allow(dead_code)] // right_assoc is recorded, as upstream records it, and never read
struct OpEntry {
    op: u8,
    prec: i32,
    right_assoc: bool,
    col: i32,
    is_func: bool,
    /// funcs
    arity: i32,
    /// funcs — for the arity diagnostic (D-R9)
    name: &'static str,
    /// '(' marker
    is_paren: bool,
    /// parens — argument-separator count (D-R9)
    commas: i32,
}

fn prec_of(op: u8) -> i32 {
    match op {
        b'+' | b'-' => 1,
        b'*' | b'/' => 2,
        b'^' => 3,
        _ => 0,
    }
}

fn bin_of(op: u8) -> u8 {
    match op {
        b'+' => RxToken::ADD,
        b'-' => RxToken::SUB,
        b'*' => RxToken::MUL,
        b'/' => RxToken::DIV,
        _ => RxToken::POW,
    }
}

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

/// Compile `src` into `pool` (appended on success), setting `out` to its span. Returns an
/// empty string on success, else the diagnostic, with `error_col` its 1-based column.
pub fn compile_reaction_expression(src: &[u8], symbols: &RxSymbols, pool: &mut Vec<RxToken>, out: &mut RxExprSpan, error_col: &mut i32) -> Bytes {
    *out = RxExprSpan { begin: pool.len() as i32, len: 0 };
    *error_col = 0;

    let mut output: Vec<RxToken> = Vec::new();
    let mut ops: Vec<OpEntry> = Vec::new();
    let mut lex = Lexer::new(src);

    let fail = |error_col: &mut i32, msg: Bytes, col: i32| -> Bytes {
        *error_col = col;
        msg
    };
    let pop_op_to_output = |output: &mut Vec<RxToken>, e: &OpEntry| {
        let t = RxToken { op: e.op, ..Default::default() };
        output.push(t);
    };

    let mut expect_operand = true; // start of expression / after operator
    let mut depth_guard = 0;

    loop {
        let mut lex_err: Bytes = Vec::new();
        let lx = lex.next(&mut lex_err);
        if !lex_err.is_empty() {
            return fail(error_col, lex_err, lx.col);
        }
        if lx.kind == LexKind::End {
            break;
        }

        match lx.kind {
            LexKind::Number => {
                if !expect_operand {
                    return fail(error_col, b"unexpected number".to_vec(), lx.col);
                }
                let t = RxToken { op: RxToken::PUSH_NUM, value: lx.num, ..Default::default() };
                output.push(t);
                expect_operand = false;
            }
            LexKind::Ident => {
                if !expect_operand {
                    return fail(error_col, cat(&[b"unexpected identifier '", &lx.ident, b"'"]), lx.col);
                }
                if let Some(f) = find_func(&lx.ident) {
                    // Must be followed by '('.
                    let mut e2: Bytes = Vec::new();
                    let la = lex.next(&mut e2);
                    if !e2.is_empty() {
                        return fail(error_col, e2, la.col);
                    }
                    if la.kind != LexKind::LParen {
                        return fail(error_col, cat(&[b"function '", &lx.ident, b"' needs '('"]), lx.col);
                    }
                    let fe = OpEntry { op: f.op, is_func: true, arity: f.arity, name: f.name, col: lx.col, ..Default::default() };
                    ops.push(fe);
                    let pe = OpEntry { is_paren: true, col: la.col, ..Default::default() };
                    ops.push(pe);
                    depth_guard += 1;
                    expect_operand = true;
                    continue;
                }
                let mut t = RxToken::default();
                let mut idx: i32 = -1;
                // Resolution order: species, coefficient, term (forward-only), pollutant,
                // hydraulic variable.
                if let Some(sp) = symbols.species {
                    for k in 0..sp.len() {
                        if sp[k] == lx.ident {
                            idx = k as i32;
                            t.op = RxToken::PUSH_SPECIES;
                            break;
                        }
                    }
                }
                if idx < 0 {
                    if let Some(cf) = symbols.coefs {
                        for k in 0..cf.len() {
                            if cf[k] == lx.ident {
                                idx = k as i32;
                                t.op = RxToken::PUSH_COEF;
                                break;
                            }
                        }
                    }
                }
                if idx < 0 {
                    if let Some(tm) = symbols.terms {
                        for k in 0..tm.len() {
                            if tm[k] == lx.ident {
                                if k as i32 >= symbols.max_term {
                                    return fail(
                                        error_col,
                                        cat(&[b"term '", &lx.ident, "' references later term \u{2014} term references are forward-only".as_bytes()]),
                                        lx.col,
                                    );
                                }
                                idx = k as i32;
                                t.op = RxToken::PUSH_TERM;
                                break;
                            }
                        }
                    }
                }
                if idx < 0 {
                    if let Some(pl) = symbols.pollutants {
                        for k in 0..pl.len() {
                            if pl[k] == lx.ident {
                                idx = k as i32;
                                t.op = RxToken::PUSH_POLLUT;
                                break;
                            }
                        }
                    }
                }
                if idx < 0 {
                    let hv = find_hyd_var(&lx.ident);
                    if hv >= 0 {
                        idx = hv;
                        t.op = RxToken::PUSH_HYDVAR;
                    }
                }
                if idx < 0 {
                    return fail(error_col, cat(&[b"undefined identifier '", &lx.ident, b"'"]), lx.col);
                }
                t.idx = idx;
                output.push(t);
                expect_operand = false;
            }
            LexKind::Op => {
                if expect_operand {
                    if lx.op == b'-' {
                        // D-R8: unary minus binds BELOW '^' (precedence 2, right-assoc).
                        let e = OpEntry { op: RxToken::NEG, prec: 2, right_assoc: true, col: lx.col, ..Default::default() };
                        ops.push(e);
                        continue;
                    }
                    if lx.op == b'+' {
                        continue; // unary plus: no-op
                    }
                    return fail(error_col, cat(&[b"unexpected operator '", &[lx.op], b"'"]), lx.col);
                }
                let p = prec_of(lx.op);
                let ra = lx.op == b'^';
                while let Some(back) = ops.last() {
                    if !(!back.is_paren && !back.is_func && (back.prec > p || (back.prec == p && !ra))) {
                        break;
                    }
                    let e = *back;
                    pop_op_to_output(&mut output, &e);
                    ops.pop();
                }
                let e = OpEntry { op: bin_of(lx.op), prec: p, right_assoc: ra, col: lx.col, ..Default::default() };
                ops.push(e);
                expect_operand = true;
            }
            LexKind::LParen => {
                if !expect_operand {
                    return fail(error_col, b"unexpected '('".to_vec(), lx.col);
                }
                let e = OpEntry { is_paren: true, col: lx.col, ..Default::default() };
                ops.push(e);
                depth_guard += 1;
                if depth_guard > K_RX_MAX_STACK_DEPTH {
                    return fail(error_col, b"expression too deep".to_vec(), lx.col);
                }
            }
            LexKind::RParen => {
                if expect_operand {
                    return fail(error_col, b"unexpected ')'".to_vec(), lx.col);
                }
                while let Some(back) = ops.last() {
                    if back.is_paren {
                        break;
                    }
                    let e = *back;
                    pop_op_to_output(&mut output, &e);
                    ops.pop();
                }
                if ops.is_empty() {
                    return fail(error_col, b"unbalanced ')'".to_vec(), lx.col);
                }
                let commas = ops.last().unwrap().commas;
                ops.pop(); // the '('
                depth_guard -= 1;
                // Function application?
                if let Some(back) = ops.last() {
                    if back.is_func {
                        let argc = commas + 1;
                        if argc != back.arity {
                            let m = cat(&[
                                b"function '",
                                back.name.as_bytes(),
                                b"' expects ",
                                back.arity.to_string().as_bytes(),
                                if back.arity == 1 { b" argument" } else { b" arguments" },
                                b", got ",
                                argc.to_string().as_bytes(),
                            ]);
                            let col = back.col;
                            return fail(error_col, m, col);
                        }
                        let e = *back;
                        pop_op_to_output(&mut output, &e);
                        ops.pop();
                    }
                }
            }
            LexKind::Comma => {
                if expect_operand {
                    return fail(error_col, b"unexpected ','".to_vec(), lx.col);
                }
                while let Some(back) = ops.last() {
                    if back.is_paren {
                        break;
                    }
                    let e = *back;
                    pop_op_to_output(&mut output, &e);
                    ops.pop();
                }
                if ops.is_empty() || ops.len() < 2 || !ops[ops.len() - 2].is_func {
                    return fail(error_col, b"',' outside a function call".to_vec(), lx.col);
                }
                ops.last_mut().unwrap().commas += 1; // the paren frame (D-R9)
                expect_operand = true;
            }
            LexKind::End => {}
        }
    }

    if expect_operand {
        return fail(error_col, b"expression ends with an operator".to_vec(), src.len() as i32 + 1);
    }
    while let Some(back) = ops.last() {
        if back.is_paren {
            let col = back.col;
            return fail(error_col, b"unbalanced '('".to_vec(), col);
        }
        let e = *back;
        pop_op_to_output(&mut output, &e);
        ops.pop();
    }
    if output.is_empty() {
        return fail(error_col, b"empty expression".to_vec(), 1);
    }

    // Verify stack discipline + final depth 1, and bound the depth.
    let mut depth = 0;
    let mut max_depth = 0;
    for t in &output {
        match t.op {
            RxToken::PUSH_NUM | RxToken::PUSH_SPECIES | RxToken::PUSH_COEF | RxToken::PUSH_TERM | RxToken::PUSH_HYDVAR | RxToken::PUSH_POLLUT => {
                depth += 1;
            }
            RxToken::NEG
            | RxToken::F_EXP
            | RxToken::F_LOG
            | RxToken::F_LOG10
            | RxToken::F_SQRT
            | RxToken::F_ABS
            | RxToken::F_SGN
            | RxToken::F_STEP
            | RxToken::F_SIN
            | RxToken::F_COS
            | RxToken::F_TAN => {
                if depth < 1 {
                    return fail(error_col, b"malformed expression".to_vec(), 1);
                }
            }
            _ => {
                // binary
                if depth < 2 {
                    return fail(error_col, b"malformed expression".to_vec(), 1);
                }
                depth -= 1;
            }
        }
        if depth > max_depth {
            max_depth = depth;
        }
    }
    if depth != 1 {
        return fail(error_col, b"malformed expression".to_vec(), 1);
    }
    if max_depth > K_RX_MAX_STACK_DEPTH {
        return fail(error_col, b"expression too deep".to_vec(), 1);
    }

    out.begin = pool.len() as i32;
    out.len = output.len() as i32;
    pool.extend_from_slice(&output);
    Vec::new()
}

// ---------------------------------------------------------------------------
// Vocabulary discovery (E-C1)
// ---------------------------------------------------------------------------

pub fn reaction_function_count() -> i32 {
    K_FUNCS.len() as i32
}

pub fn reaction_function_info(idx: i32) -> Option<(&'static str, i32)> {
    if idx < 0 || idx >= reaction_function_count() {
        return None;
    }
    Some((K_FUNCS[idx as usize].name, K_FUNCS[idx as usize].arity))
}

pub fn reaction_hyd_var_count() -> i32 {
    RxHydVar::Count as i32
}

/// (name, description), or None out of range (upstream returns false and leaves the
/// outputs unwritten).
pub fn reaction_hyd_var_info(idx: i32) -> Option<(&'static str, &'static str)> {
    if idx < 0 || idx >= reaction_hyd_var_count() {
        return None;
    }
    Some((K_HYD_VAR_NAMES[idx as usize], K_HYD_VAR_DESCS[idx as usize]))
}

pub fn reaction_hyd_var_populated(idx: i32) -> bool {
    idx >= RxHydVar::Av as i32 && idx <= RxHydVar::Temp as i32
}

/// Physics review A1: refuse a compiled expression that would silently read zero in the
/// given context (0 pipe/tank, 1 term, 2 groundwater, 3 wet pool, 4 LID). Empty if valid.
pub fn reaction_context_diagnostic(rx: &ReactionData, pool: &[RxToken], span: &RxExprSpan, context: i32) -> Bytes {
    let mut i = span.begin;
    while i < span.begin + span.len {
        if i < 0 || i as usize >= pool.len() {
            return b"invalid compiled reaction context span".to_vec();
        }
        let token = pool[i as usize];
        if token.op == RxToken::PUSH_HYDVAR {
            if token.idx == RxHydVar::Av as i32 && !rx.wall_geometry {
                return b"hydraulic variable 'AV' requires [REACTION_OPTIONS] WALL_GEOMETRY YES".to_vec();
            }
            let ground = token.idx >= RxHydVar::GwTheta as i32;
            let valid = if context >= 3 {
                token.idx == RxHydVar::Dt as i32 || token.idx == RxHydVar::Temp as i32
            } else if context == 2 {
                token.idx >= RxHydVar::Dt as i32
            } else {
                reaction_hyd_var_populated(token.idx) || (context == 1 && ground)
            };
            if !valid {
                let name: &str = match reaction_hyd_var_info(token.idx) {
                    Some((n, _)) => n,
                    None => "unknown",
                };
                return cat(&[b"hydraulic variable '", name.as_bytes(), b"' is unavailable in this reaction context"]);
            }
        }
        if context == 2 && token.op == RxToken::PUSH_SPECIES && rx.species_is_wall[token.idx as usize] != 0 {
            return b"groundwater reactions cannot depend on WALL species".to_vec();
        }
        if context >= 3 && token.op == RxToken::PUSH_SPECIES && rx.species_is_wall[token.idx as usize] != 0 {
            return b"wet-pool/LID reactions cannot depend on WALL species".to_vec();
        }
        if context >= 3 && token.op == RxToken::PUSH_POLLUT {
            return b"wet-pool/LID reactions cannot read pollutant concentrations (no local pollutant context)".to_vec();
        }
        if token.op == RxToken::PUSH_TERM {
            if token.idx < 0 || token.idx as usize >= rx.term_expr.len() {
                return b"invalid reaction term reference".to_vec();
            }
            let error = reaction_context_diagnostic(rx, &rx.token_pool, &rx.term_expr[token.idx as usize], context);
            if !error.is_empty() {
                return error;
            }
        }
        i += 1;
    }
    Vec::new()
}
