// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/LidNode.cpp, LidNodeTreatment.cpp and
// LidNodeRichards.cpp (the [LID_NODES], [LID_NODE_OUTLETS], [LID_LAYER_TREATMENT] and
// [LID_RICHARDS] readers, lidnode::validate and every helper they call), with
// richards::valid from src/engine/hydrology/RichardsColumn.cpp and treatment::validate /
// treatment::parse (with tokenize and shuntingYard) from src/engine/quality/Treatment.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// (LidNodeRichards.cpp and RichardsColumn.cpp carry no copyright line of their own; their
// doc comments give "Copyright (c) 2026 Caleb Buahin" and "@license Apache-2.0".)
//
// Differences in form only:
//   * readTreatment and readRichards extract with std::istringstream (`>>` into std::string,
//     int and double, std::ws, std::getline). `IStringStream` below reproduces libstdc++ 13
//     in the classic "C" locale: the sentry (skip C-locale white space, eofbit|failbit at the
//     end, failbit when the stream is not good), operator>>(string&), num_get's
//     _M_extract_int for `long` narrowed by operator>>(int&) (out of int range: failbit),
//     num_get's _M_extract_float followed by __convert_to_v (glibc strtod on the collected
//     text, which must be consumed whole; +-inf becomes +-DBL_MAX with failbit; an underflow
//     is accepted), std::ws and getline. Only eofbit and failbit are kept: a stringbuf never
//     sets badbit.
//   * `std::strtod` is core::cstdlib::strtod on the token up to its first NUL, as c_str()
//     would give it; `end != c_str` is "consumed > 0" and `*end == '\0'` is "consumed it all".
//   * std::stod (Treatment.cpp's tokenizers) throws invalid_argument when strtod converts
//     nothing and out_of_range when strtod sets ERANGE: on overflow, and on underflow (glibc,
//     with x86's tininess after rounding: a result that is subnormal or zero and inexact).
//     Here that is a result that is infinite, or zero from a literal holding a nonzero digit,
//     or subnormal. A subnormal is exact only for a literal of about 750 significant digits;
//     that case is not distinguished. A throw is Err(()). validate catches it (`catch (...)`);
//     the parse that follows only sees numbers validate already converted, so its Err cannot
//     occur from these readers.
//   * The pollutant lookup of treatment::parse(expr, result, names) is a closure, called only
//     where upstream calls its std::function.
//
// Not translated here (run-time; the hydrology and quality phases translate them):
//   LidNode.cpp: sync, active, heldVolume, infiltrationDeficit, reconcileInfiltration,
//     initializeInfiltration, initialize, prepareStepImpl, finishRows, prepareStep,
//     prepareStepAfterMesh, finishStep, finishRichardsNodes, bottomMode, bottomConductivity,
//     bottomCloggingFactor, portCell, resetPorts, portDepth, exchangePorts.
//   LidNodeTreatment.cpp: layerVolume, treatedConcentration, prepareQuality,
//     prepareOutletQuality, receiveQuality, outletQuality, heldMass.
//   LidNodeRichards.cpp: everything but readRichards (richardsMode, exchangeRichardsPorts,
//     refresh, refreshRichardsState, prepareRichardsStep, validateRichardsRouting and the
//     detail:: helpers).
//   RichardsColumn.cpp: everything but the two valid() overloads.
//   Treatment.cpp: both evaluate overloads, applyTreatment and the function-pointer parse.

use crate::core::cmath::cmax;
use crate::core::context::SimulationContext;
use crate::core::cstdlib::{c_str, strtod};
use crate::core::options::{Bytes, QualitySolverKind, RoutingModel};
use crate::core::unit_conversion as ucf;
use crate::data::lid_node::{richards, LidLayerTreatment, LidNodeBottom, LidNodeConfig, LidNodeLayer, LidNodeLayerKind, LidNodeOutlet};
use crate::input::tokenizer;

fn cat(parts: &[&[u8]]) -> Bytes {
    parts.concat()
}

/// ctype<char>::is(space, c) in the classic table, and isspace in the C locale.
#[inline]
fn c_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

// ============================================================================
// std::istringstream, as libstdc++ 13 extracts in the "C" locale
// ============================================================================

/// A std::istringstream over one line, read with `>>` (std::string, int, double), std::ws
/// and std::getline as libstdc++ 13 does in the classic locale. `eof` and `fail` are the
/// stream's eofbit and failbit. Public so that other readers, and a differential against
/// libstdc++, can use it.
pub struct IStringStream<'a> {
    s: &'a [u8],
    pos: usize,
    pub eof: bool,
    pub fail: bool,
}

impl<'a> IStringStream<'a> {
    pub fn new(s: &'a [u8]) -> Self {
        IStringStream { s, pos: 0, eof: false, fail: false }
    }

    pub fn good(&self) -> bool {
        !self.eof && !self.fail
    }

    /// `bool(stream)`: !fail().
    pub fn ok(&self) -> bool {
        !self.fail
    }

    /// basic_istream::sentry(in, noskip).
    fn sentry(&mut self, noskip: bool) -> bool {
        let mut err_eof = false;
        if self.good() && !noskip {
            while self.pos < self.s.len() && c_space(self.s[self.pos]) {
                self.pos += 1;
            }
            if self.pos == self.s.len() {
                err_eof = true;
            }
        }
        if self.good() && !err_eof {
            true
        } else {
            if err_eof {
                self.eof = true;
            }
            self.fail = true;
            false
        }
    }

    /// operator>>(istream&, std::string&).
    pub fn read_string(&mut self, out: &mut Bytes) -> &mut Self {
        let mut extracted = 0usize;
        let mut err_eof = false;
        if self.sentry(false) {
            out.clear();
            while self.pos < self.s.len() && !c_space(self.s[self.pos]) {
                out.push(self.s[self.pos]);
                self.pos += 1;
                extracted += 1;
            }
            if self.pos == self.s.len() {
                err_eof = true;
            }
        }
        if extracted == 0 {
            self.fail = true;
        }
        if err_eof {
            self.eof = true;
        }
        self
    }

    /// num_get<char>::_M_extract_int for `long`, basefield dec, C locale (no grouping,
    /// decimal point '.'). Returns (value, failbit, eofbit).
    fn extract_long(&mut self) -> (i64, bool, bool) {
        let s = self.s;
        let end = s.len();
        let mut c: u8 = 0;
        let mut testeof = self.pos == end;
        let mut negative = false;
        if !testeof {
            c = s[self.pos];
            negative = c == b'-';
            if negative || c == b'+' {
                self.pos += 1;
                if self.pos != end {
                    c = s[self.pos];
                } else {
                    testeof = true;
                }
            }
        }
        // leading zeros (base 10: any number of them)
        let mut found_zero = false;
        let mut sep_pos = 0i32;
        while !testeof {
            if c == b'.' {
                break;
            } else if c == b'0' {
                found_zero = true;
                sep_pos += 1;
            } else {
                // 'x'/'X' after a zero breaks too: the base stays 10
                break;
            }
            self.pos += 1;
            if self.pos != end {
                c = s[self.pos];
                if !found_zero {
                    break;
                }
            } else {
                testeof = true;
            }
        }
        let max: u64 = if negative { (i64::MIN as u64).wrapping_neg() } else { i64::MAX as u64 };
        let smax = max / 10;
        let mut result: u64 = 0;
        let mut testoverflow = false;
        while !testeof {
            if !(c >= b'0' && c < b'0' + 10) {
                break;
            }
            let digit = (c - b'0') as u64;
            if result > smax {
                testoverflow = true;
            } else {
                result = result.wrapping_mul(10);
                testoverflow |= result > max - digit;
                result = result.wrapping_add(digit);
                sep_pos += 1;
            }
            self.pos += 1;
            if self.pos != end {
                c = s[self.pos];
            } else {
                testeof = true;
            }
        }
        let (v, fail) = if sep_pos == 0 && !found_zero {
            (0i64, true)
        } else if testoverflow {
            (if negative { i64::MIN } else { i64::MAX }, true)
        } else {
            (if negative { result.wrapping_neg() as i64 } else { result as i64 }, false)
        };
        (v, fail, testeof)
    }

    /// basic_istream::operator>>(int&).
    pub fn read_int(&mut self, n: &mut i32) -> &mut Self {
        if self.sentry(false) {
            let (l, mut fail, eof) = self.extract_long();
            if l < i32::MIN as i64 {
                fail = true;
                *n = i32::MIN;
            } else if l > i32::MAX as i64 {
                fail = true;
                *n = i32::MAX;
            } else {
                *n = l as i32;
            }
            if fail {
                self.fail = true;
            }
            if eof {
                self.eof = true;
            }
        }
        self
    }

    /// num_get<char>::_M_extract_float, C locale: the characters a double may be built from.
    fn extract_float(&mut self, xtrc: &mut Bytes) -> bool {
        let s = self.s;
        let end = s.len();
        let mut c: u8 = 0;
        let mut testeof = self.pos == end;
        if !testeof {
            c = s[self.pos];
            let plus = c == b'+';
            if plus || c == b'-' {
                xtrc.push(if plus { b'+' } else { b'-' });
                self.pos += 1;
                if self.pos != end {
                    c = s[self.pos];
                } else {
                    testeof = true;
                }
            }
        }
        let mut found_mantissa = false;
        while !testeof {
            if c == b'.' {
                break;
            } else if c == b'0' {
                if !found_mantissa {
                    xtrc.push(b'0');
                    found_mantissa = true;
                }
                self.pos += 1;
                if self.pos != end {
                    c = s[self.pos];
                } else {
                    testeof = true;
                }
            } else {
                break;
            }
        }
        let mut found_dec = false;
        let mut found_sci = false;
        while !testeof {
            if c >= b'0' && c < b'0' + 10 {
                xtrc.push(c);
                found_mantissa = true;
            } else if c == b'.' && !found_dec && !found_sci {
                xtrc.push(b'.');
                found_dec = true;
            } else if (c == b'e' || c == b'E') && !found_sci && found_mantissa {
                xtrc.push(b'e');
                found_sci = true;
                self.pos += 1;
                if self.pos != end {
                    c = s[self.pos];
                    let plus = c == b'+';
                    if plus || c == b'-' {
                        xtrc.push(if plus { b'+' } else { b'-' });
                    } else {
                        continue;
                    }
                } else {
                    testeof = true;
                    break;
                }
            } else {
                break;
            }
            self.pos += 1;
            if self.pos != end {
                c = s[self.pos];
            } else {
                testeof = true;
            }
        }
        testeof
    }

    /// basic_istream::operator>>(double&): num_get::do_get, then __convert_to_v.
    pub fn read_double(&mut self, v: &mut f64) -> &mut Self {
        if self.sentry(false) {
            let mut xtrc: Bytes = Vec::with_capacity(32);
            let mut eof = self.extract_float(&mut xtrc);
            let mut fail = false;
            let (x, used) = strtod(&xtrc);
            if used == 0 || used != xtrc.len() {
                *v = 0.0;
                fail = true;
            } else if x == f64::INFINITY {
                *v = f64::MAX;
                fail = true;
            } else if x == f64::NEG_INFINITY {
                *v = -f64::MAX;
                fail = true;
            } else {
                *v = x;
            }
            if self.pos == self.s.len() {
                eof = true;
            }
            if fail {
                self.fail = true;
            }
            if eof {
                self.eof = true;
            }
        }
        self
    }

    /// `stream >> std::ws`.
    pub fn ws(&mut self) -> &mut Self {
        if self.sentry(true) {
            loop {
                if self.pos == self.s.len() {
                    self.eof = true;
                    break;
                }
                if !c_space(self.s[self.pos]) {
                    break;
                }
                self.pos += 1;
            }
        }
        self
    }

    /// std::getline(stream, str) with delimiter '\n'.
    pub fn getline(&mut self, out: &mut Bytes) -> &mut Self {
        let mut extracted = 0usize;
        let mut err_eof = false;
        if self.sentry(true) {
            out.clear();
            while self.pos < self.s.len() && self.s[self.pos] != b'\n' {
                out.push(self.s[self.pos]);
                self.pos += 1;
                extracted += 1;
            }
            if self.pos == self.s.len() {
                err_eof = true;
            } else {
                extracted += 1;
                self.pos += 1;
            }
        }
        if extracted == 0 {
            self.fail = true;
        }
        if err_eof {
            self.eof = true;
        }
        self
    }
}

// ============================================================================
// RichardsColumn.cpp: richards::valid
// ============================================================================

/// richards::valid(const Material&, double theta_s).
pub fn richards_valid_material(p: &richards::Material, ts: f64) -> bool {
    ts.is_finite()
        && ts > 0.0
        && ts <= 1.0
        && p.theta_r.is_finite()
        && p.theta_r >= 0.0
        && p.theta_r < ts
        && p.alpha.is_finite()
        && p.alpha > 0.0
        && p.n.is_finite()
        && p.n > 1.0
        && p.l.is_finite()
        && p.l >= 0.0
        && p.specific_storage.is_finite()
        && p.specific_storage > 0.0
}

/// richards::valid(const Options&).
pub fn richards_valid_options(o: &richards::Options) -> bool {
    o.cells_per_layer >= 1
        && o.cells_per_layer <= 256
        && o.atol.is_finite()
        && o.atol > 0.0
        && o.rtol.is_finite()
        && o.rtol > 0.0
        && o.max_step.is_finite()
        && o.max_step > 0.0
}

// ============================================================================
// Treatment.cpp: validate and parse (with tokenize and shuntingYard)
// ============================================================================

pub mod treatment {
    use super::{c_space, cat};
    use crate::core::cstdlib::strtod;
    use crate::core::options::Bytes;
    use crate::data::treatment_types::{Token, TokenType, TreatExpr, TreatVar};

    fn precedence(t: TokenType) -> i32 {
        match t {
            TokenType::Add | TokenType::Sub => 1,
            TokenType::Mul | TokenType::Div => 2,
            TokenType::Pow => 3,
            TokenType::Neg => 4,
            _ => 0,
        }
    }

    fn is_function(t: TokenType) -> bool {
        let v = t as i32;
        v >= TokenType::FuncExp as i32 && v <= TokenType::FuncStep as i32
    }

    fn is_operator(t: TokenType) -> bool {
        matches!(t, TokenType::Add | TokenType::Sub | TokenType::Mul | TokenType::Div | TokenType::Pow | TokenType::Neg)
    }

    fn is_right_assoc(t: TokenType) -> bool {
        t == TokenType::Pow || t == TokenType::Neg
    }

    /// func_map: function names (lower case).
    fn func_map(lower: &[u8]) -> Option<TokenType> {
        Some(match lower {
            b"exp" => TokenType::FuncExp,
            b"log" | b"ln" => TokenType::FuncLog,
            b"sqrt" => TokenType::FuncSqrt,
            b"min" => TokenType::FuncMin,
            b"max" => TokenType::FuncMax,
            b"abs" => TokenType::FuncAbs,
            b"sgn" => TokenType::FuncSgn,
            b"step" => TokenType::FuncStep,
            _ => return None,
        })
    }

    /// var_map: variable names (upper case).
    fn var_map(upper: &[u8]) -> Option<TreatVar> {
        Some(match upper {
            b"C" => TreatVar::C,
            b"R" => TreatVar::R,
            b"DT" => TreatVar::Dt,
            b"HRT" => TreatVar::Hrt,
            b"FLOW" | b"Q" => TreatVar::Q,
            b"V" => TreatVar::V,
            b"DEPTH" | b"D" => TreatVar::D,
            b"AREA" => TreatVar::Area,
            _ => return None,
        })
    }

    fn token(t: TokenType) -> Token {
        Token { type_: t, ..Token::default() }
    }

    /// std::stod: Err where it throws (see the module notes).
    fn stod(word: &[u8]) -> Result<f64, ()> {
        let (v, used) = strtod(word);
        if used == 0 {
            return Err(()); // invalid_argument
        }
        if v.is_infinite() {
            return Err(()); // out_of_range (overflow)
        }
        if v == 0.0 {
            // the consumed mantissa (before any exponent) holds a nonzero digit: underflow
            let mant = word[..used].iter().take_while(|&&c| c != b'e' && c != b'E');
            if mant.into_iter().any(|&c| (b'1'..=b'9').contains(&c)) {
                return Err(());
            }
        } else if v.abs() < f64::MIN_POSITIVE {
            return Err(()); // subnormal: tiny and (see notes) inexact
        }
        Ok(v)
    }

    /// The number scan both tokenizers use: digits, '.', 'e', 'E', and a sign right after
    /// an exponent letter. Returns the end index.
    fn scan_number(s: &[u8], mut i: usize) -> usize {
        while i < s.len()
            && (s[i].is_ascii_digit()
                || s[i] == b'.'
                || s[i] == b'e'
                || s[i] == b'E'
                || ((s[i] == b'+' || s[i] == b'-') && i > 0 && (s[i - 1] == b'e' || s[i - 1] == b'E')))
        {
            i += 1;
        }
        i
    }

    fn unary_minus_here(tokens: &[Token]) -> bool {
        match tokens.last() {
            None => true,
            Some(b) => b.type_ == TokenType::Lparen || is_operator(b.type_) || b.type_ == TokenType::Comma,
        }
    }

    /// tokenize(s, pollut_lookup): Ok(empty) for an unknown identifier or pollutant; Err(())
    /// where std::stod throws.
    fn tokenize(s: &[u8], pollut_lookup: Option<&dyn Fn(&[u8]) -> i32>) -> Result<Vec<Token>, ()> {
        let mut tokens: Vec<Token> = Vec::new();
        let mut i = 0usize;
        while i < s.len() {
            let c = s[i];
            if c_space(c) {
                i += 1;
                continue;
            }
            if c.is_ascii_digit() || c == b'.' {
                let start = i;
                i = scan_number(s, i);
                let mut t = token(TokenType::Number);
                t.value = stod(&s[start..i])?;
                tokens.push(t);
                continue;
            }
            if c.is_ascii_alphabetic() || c == b'_' {
                let start = i;
                while i < s.len() && (s[i].is_ascii_alphanumeric() || s[i] == b'_') {
                    i += 1;
                }
                let word = &s[start..i];
                let lower = word.to_ascii_lowercase();
                if let Some(ft) = func_map(&lower) {
                    tokens.push(token(ft));
                    continue;
                }
                let upper = word.to_ascii_uppercase();
                if let Some(v) = var_map(&upper) {
                    let mut t = token(TokenType::Variable);
                    t.var = v;
                    tokens.push(t);
                } else if pollut_lookup.is_some() && upper.len() > 2 && upper[1] == b'_' {
                    let lookup = pollut_lookup.unwrap();
                    let prefix = upper[0];
                    let pollut_name = &word[2..];
                    let pi = lookup(pollut_name);
                    if pi >= 0 {
                        let mut t = token(TokenType::Variable);
                        t.var = if prefix == b'R' { TreatVar::RPollut } else { TreatVar::CPollut };
                        t.pollut_ref = pi;
                        tokens.push(t);
                    } else {
                        return Ok(Vec::new()); // unknown pollutant name
                    }
                } else if pollut_lookup.is_some() && (pollut_lookup.unwrap())(word) >= 0 {
                    let mut t = token(TokenType::Variable);
                    t.var = TreatVar::CPollut;
                    t.pollut_ref = (pollut_lookup.unwrap())(word);
                    tokens.push(t);
                } else {
                    return Ok(Vec::new()); // unknown identifier
                }
                continue;
            }
            let ty = match c {
                b'+' => TokenType::Add,
                b'-' => {
                    if unary_minus_here(&tokens) {
                        TokenType::Neg
                    } else {
                        TokenType::Sub
                    }
                }
                b'*' => TokenType::Mul,
                b'/' => TokenType::Div,
                b'^' => TokenType::Pow,
                b'(' => TokenType::Lparen,
                b')' => TokenType::Rparen,
                b',' => TokenType::Comma,
                _ => {
                    i += 1; // skip unknown characters
                    continue;
                }
            };
            tokens.push(token(ty));
            i += 1;
        }
        Ok(tokens)
    }

    /// shuntingYard(infix, postfix): 0, or -1 for unpaired parentheses / a stray comma.
    fn shunting_yard(infix: &[Token], postfix: &mut Vec<Token>) -> i32 {
        postfix.clear();
        let mut ops: Vec<Token> = Vec::new();
        for tok in infix {
            if tok.type_ == TokenType::Number || tok.type_ == TokenType::Variable {
                postfix.push(*tok);
            } else if is_function(tok.type_) {
                ops.push(*tok);
            } else if tok.type_ == TokenType::Comma {
                while !ops.is_empty() && ops.last().unwrap().type_ != TokenType::Lparen {
                    postfix.push(ops.pop().unwrap());
                }
                if ops.is_empty() {
                    return -1;
                }
            } else if is_operator(tok.type_) {
                while let Some(top) = ops.last() {
                    let pops = is_operator(top.type_)
                        && ((is_right_assoc(tok.type_) && precedence(tok.type_) < precedence(top.type_))
                            || (!is_right_assoc(tok.type_) && precedence(tok.type_) <= precedence(top.type_)));
                    if !pops {
                        break;
                    }
                    postfix.push(ops.pop().unwrap());
                }
                ops.push(*tok);
            } else if tok.type_ == TokenType::Lparen {
                ops.push(*tok);
            } else if tok.type_ == TokenType::Rparen {
                while !ops.is_empty() && ops.last().unwrap().type_ != TokenType::Lparen {
                    postfix.push(ops.pop().unwrap());
                }
                if ops.is_empty() {
                    return -1;
                }
                ops.pop();
                if !ops.is_empty() && is_function(ops.last().unwrap().type_) {
                    postfix.push(ops.pop().unwrap());
                }
            }
        }
        while let Some(top) = ops.pop() {
            if top.type_ == TokenType::Lparen {
                return -1;
            }
            postfix.push(top);
        }
        0
    }

    /// find_first_not_of / find_last_not_of over " \t\r\n".
    fn lhs_ws(c: u8) -> bool {
        matches!(c, b' ' | b'\t' | b'\r' | b'\n')
    }

    /// parse(expr_str, result): the built-in-only parser. Err(()) where std::stod throws.
    pub fn parse(expr_str: &[u8], result: &mut TreatExpr) -> Result<i32, ()> {
        result.tokens.clear();
        result.is_removal = false;
        if expr_str.is_empty() {
            return Ok(-1);
        }
        let eq_pos = match expr_str.iter().position(|&c| c == b'=') {
            Some(p) => p,
            None => return Ok(-1),
        };
        let lhs = &expr_str[..eq_pos];
        let lhs_start = match lhs.iter().position(|&c| !lhs_ws(c)) {
            Some(p) => p,
            None => return Ok(-1),
        };
        let lhs_end = lhs.iter().rposition(|&c| !lhs_ws(c)).unwrap();
        let lhs = &lhs[lhs_start..lhs_end + 1];
        let lhs_char = lhs[0].to_ascii_uppercase();
        if lhs.len() == 1 && lhs_char == b'R' {
            result.is_removal = true;
        } else if lhs.len() == 1 && lhs_char == b'C' {
            result.is_removal = false;
        } else {
            return Ok(-1);
        }
        let rhs = &expr_str[eq_pos + 1..];
        let infix = tokenize(rhs, None)?;
        if infix.is_empty() && !rhs.is_empty() {
            return Ok(-1);
        }
        if infix.is_empty() {
            return Ok(-1);
        }
        if shunting_yard(&infix, &mut result.tokens) != 0 {
            result.tokens.clear();
            return Ok(-1);
        }
        Ok(0)
    }

    /// validate(expr_str, msg, col).
    pub fn validate(expr_str: &[u8], msg: &mut Bytes, col: &mut i32) -> i32 {
        msg.clear();
        *col = -1;
        let mut fail = |m: Bytes, c: i32| -> i32 {
            *msg = m;
            *col = c;
            -1
        };
        if expr_str.is_empty() {
            return fail(b"expression is empty".to_vec(), 0);
        }
        let eq_pos = match expr_str.iter().position(|&c| c == b'=') {
            Some(p) => p,
            None => {
                return fail(
                    "missing '=' — expression must have the form \"R = ...\" or \"C = ...\"".as_bytes().to_vec(),
                    0,
                )
            }
        };
        let lhs = &expr_str[..eq_pos];
        let lhs_start = match lhs.iter().position(|&c| !lhs_ws(c)) {
            Some(p) => p,
            None => return fail(b"left-hand side must be 'R' (removal) or 'C' (concentration)".to_vec(), 0),
        };
        let lhs_end = lhs.iter().rposition(|&c| !lhs_ws(c)).unwrap();
        let lhs = &lhs[lhs_start..lhs_end + 1];
        let lc = lhs[0].to_ascii_uppercase();
        if lhs.len() != 1 || (lc != b'R' && lc != b'C') {
            return fail(cat(&[b"left-hand side must be exactly 'R' or 'C', got '", lhs, b"'"]), lhs_start as i32);
        }
        let rhs = &expr_str[eq_pos + 1..];
        let rhs_off = eq_pos as i32 + 1;

        let mut tokens: Vec<Token> = Vec::new();
        let mut i = 0usize;
        while i < rhs.len() {
            let c = rhs[i];
            if c_space(c) {
                i += 1;
                continue;
            }
            if c.is_ascii_digit() || c == b'.' {
                let start = i;
                i = scan_number(rhs, i);
                let mut t = token(TokenType::Number);
                match stod(&rhs[start..i]) {
                    Ok(v) => t.value = v,
                    Err(()) => {
                        return fail(cat(&[b"malformed number '", &rhs[start..i], b"'"]), rhs_off + start as i32)
                    }
                }
                tokens.push(t);
                continue;
            }
            if c.is_ascii_alphabetic() || c == b'_' {
                let start = i;
                while i < rhs.len() && (rhs[i].is_ascii_alphanumeric() || rhs[i] == b'_') {
                    i += 1;
                }
                let word = &rhs[start..i];
                let lower = word.to_ascii_lowercase();
                if let Some(ft) = func_map(&lower) {
                    tokens.push(token(ft));
                    continue;
                }
                let upper = word.to_ascii_uppercase();
                if let Some(v) = var_map(&upper) {
                    let mut t = token(TokenType::Variable);
                    t.var = v;
                    tokens.push(t);
                    continue;
                }
                if upper.len() > 2 && upper[1] == b'_' && (upper[0] == b'R' || upper[0] == b'C') {
                    return fail(
                        cat(&[b"co-pollutant reference '", word, b"' is not supported here"]),
                        rhs_off + start as i32,
                    );
                }
                return fail(
                    cat(&[
                        b"unknown identifier '",
                        word,
                        b"' (variables: C R DT HRT Q V D AREA; functions: exp log ln sqrt min max abs sgn step)",
                    ]),
                    rhs_off + start as i32,
                );
            }
            let ty = match c {
                b'+' => TokenType::Add,
                b'-' => {
                    if unary_minus_here(&tokens) {
                        TokenType::Neg
                    } else {
                        TokenType::Sub
                    }
                }
                b'*' => TokenType::Mul,
                b'/' => TokenType::Div,
                b'^' => TokenType::Pow,
                b'(' => TokenType::Lparen,
                b')' => TokenType::Rparen,
                b',' => TokenType::Comma,
                _ => return fail(cat(&[b"unexpected character '", &[c], b"'"]), rhs_off + i as i32),
            };
            tokens.push(token(ty));
            i += 1;
        }
        if tokens.is_empty() {
            return fail(b"missing right-hand side after '='".to_vec(), rhs_off);
        }

        let mut postfix: Vec<Token> = Vec::new();
        if shunting_yard(&tokens, &mut postfix) != 0 {
            return fail(b"mismatched parentheses or misplaced comma".to_vec(), -1);
        }
        let mut depth = 0i32;
        for t in &postfix {
            match t.type_ {
                TokenType::Number | TokenType::Variable => depth += 1,
                TokenType::Neg
                | TokenType::FuncExp
                | TokenType::FuncLog
                | TokenType::FuncSqrt
                | TokenType::FuncAbs
                | TokenType::FuncSgn
                | TokenType::FuncStep => {
                    if depth < 1 {
                        return fail(b"function or operator is missing its operand".to_vec(), -1);
                    }
                }
                TokenType::FuncMin
                | TokenType::FuncMax
                | TokenType::Add
                | TokenType::Sub
                | TokenType::Mul
                | TokenType::Div
                | TokenType::Pow => {
                    if depth < 2 {
                        return fail(b"operator or min/max is missing an operand".to_vec(), -1);
                    }
                    depth -= 1;
                }
                _ => return fail(b"internal: unexpected token in parsed expression".to_vec(), -1),
            }
        }
        if depth != 1 {
            return fail(b"incomplete expression".to_vec(), -1);
        }
        // Drift guard. Its std::stod sees only the numbers converted above, so it cannot throw.
        let mut probe = TreatExpr::default();
        if parse(expr_str, &mut probe) != Ok(0) {
            return fail(b"invalid expression".to_vec(), -1);
        }
        0
    }

    /// parse_impl(expr_str, result, pollut_lookup). Err(()) where std::stod throws.
    fn parse_impl(expr_str: &[u8], result: &mut TreatExpr, pollut_lookup: Option<&dyn Fn(&[u8]) -> i32>) -> Result<i32, ()> {
        result.tokens.clear();
        result.pollutant_idx = -1;
        result.is_removal = false;
        if expr_str.is_empty() {
            return Ok(-1);
        }
        let eq = match expr_str.iter().position(|&c| c == b'=') {
            Some(p) => p,
            None => return Ok(-1),
        };
        if eq == 0 {
            return Ok(-1);
        }
        let lhs = expr_str[0].to_ascii_uppercase();
        if lhs == b'R' {
            result.is_removal = true;
        } else if lhs == b'C' {
            result.is_removal = false;
        } else {
            return Ok(-1);
        }
        let rhs = &expr_str[eq + 1..];
        let tokens = tokenize(rhs, pollut_lookup)?;
        if tokens.is_empty() && !rhs.is_empty() {
            return Ok(-1);
        }
        if shunting_yard(&tokens, &mut result.tokens) != 0 {
            return Ok(-1);
        }
        Ok(0)
    }

    /// parse(expr_str, result, pollut_names): names matched case-insensitively (toupper on
    /// each byte, C locale), first match wins. Err(()) where std::stod throws.
    pub fn parse_names(expr_str: &[u8], result: &mut TreatExpr, pollut_names: &[Bytes]) -> Result<i32, ()> {
        result.tokens.clear();
        result.pollutant_idx = -1;
        result.is_removal = false;
        let lookup = |name: &[u8]| -> i32 {
            for (i, p) in pollut_names.iter().enumerate() {
                if p.len() != name.len() {
                    continue;
                }
                let mut same = true;
                let mut k = 0;
                while k < name.len() && same {
                    same = p[k].to_ascii_uppercase() == name[k].to_ascii_uppercase();
                    k += 1;
                }
                if same {
                    return i as i32;
                }
            }
            -1
        };
        parse_impl(expr_str, result, Some(&lookup))
    }
}

// ============================================================================
// LidNode.cpp: helpers
// ============================================================================

/// number(): one complete finite token through strtod. `v` always receives strtod's result.
fn number(s: &[u8], v: &mut f64) -> bool {
    let cs = c_str(s);
    let (x, used) = strtod(cs);
    *v = x;
    used != 0 && used == cs.len() && x.is_finite()
}

/// rainDepth(): display rain-depth units per native ft.
fn rain_depth(ctx: &SimulationContext) -> f64 {
    ucf::UCF[ucf::RAINDEPTH as usize][ucf::get_unit_system(ctx.options.flow_units as i32) as usize]
}

/// outletNode(): the one endpoint of `link` with an assigned LID control, or -1.
fn outlet_node(ctx: &SimulationContext, link: i32) -> i32 {
    if link < 0 || link >= ctx.n_links() {
        return -1;
    }
    let mut found = -1;
    for node in [ctx.links.node1[link as usize], ctx.links.node2[link as usize]] {
        let r = ctx.node_subtypes.storage_row(node);
        if r >= 0 && ctx.node_subtypes.storages.lid[r as usize].control >= 0 {
            if found >= 0 {
                return -1; // three-column grammar is ambiguous for two LID ends
            }
            found = node;
        }
    }
    found
}

/// anchorHeight(): a one-based physical layer interface above the bottom of the stack.
fn anchor_height(stack: &[LidNodeLayer], layer: i32, top: bool) -> f64 {
    let mut n = 0;
    let total = thickness(stack);
    let mut down = 0.0;
    for l in stack {
        if l.kind == LidNodeLayerKind::Bottom {
            continue;
        }
        n += 1;
        if n == layer {
            return total - down - (if top { 0.0 } else { l.params[0] });
        }
        down += l.params[0];
    }
    -1.0
}

fn make_layer(kind: LidNodeLayerKind, params: &[f64]) -> LidNodeLayer {
    let mut l = LidNodeLayer { kind, ..LidNodeLayer::default() };
    l.params[..params.len()].copy_from_slice(params);
    l
}

/// layers(): the explicit NODE stack, or a BC/RG/IT/PP control adapted to one.
pub fn layers(ctx: &SimulationContext, control: i32) -> Vec<LidNodeLayer> {
    let c = &ctx.lid_controls;
    if control < 0 || control >= c.count() {
        return Vec::new();
    }
    let uc = control as usize;
    let type_ = c.lid_type[uc].as_slice();
    if type_ == b"NODE" {
        return if uc < c.node_layers.len() { c.node_layers[uc].clone() } else { Vec::new() };
    }
    if type_ != b"BC" && type_ != b"RG" && type_ != b"IT" && type_ != b"PP" {
        return Vec::new();
    }
    let mut out: Vec<LidNodeLayer> = Vec::new();
    let s = &c.surface[uc];
    if s[0] > 0.0 {
        out.push(make_layer(LidNodeLayerKind::Surface, &[s[0], s[1]]));
    }
    // Pavement void ratio is converted once; impermeable paver fraction reduces
    // available pore space. Its hydraulic layer remains above the soil.
    let p = &c.pavement[uc];
    if type_ == b"PP" && p[0] > 0.0 {
        out.push(make_layer(LidNodeLayerKind::Aggregate, &[p[0], p[1] / (1.0 + p[1]) * (1.0 - p[2]), p[3]]));
    }
    let m = &c.soil[uc];
    if m[0] > 0.0 {
        out.push(make_layer(LidNodeLayerKind::Media, m));
    }
    let g = &c.storage[uc];
    if g[0] > 0.0 {
        // Standard storage K is native-soil seepage, not gravel conductivity.
        // Use effectively free drainage within the aggregate.
        out.push(make_layer(LidNodeLayerKind::Aggregate, &[g[0], g[1] / (1.0 + g[1]), 1.0e12]));
        out.push(make_layer(LidNodeLayerKind::Bottom, &[g[2], g[3]]));
    }
    out
}

/// validLayer(): one layer's parameters.
pub fn valid_layer(l: &LidNodeLayer) -> bool {
    let p = &l.params;
    for &v in p.iter() {
        if !v.is_finite() || v < 0.0 {
            return false;
        }
    }
    let used = match l.kind {
        LidNodeLayerKind::Media => 7,
        LidNodeLayerKind::Aggregate => 3,
        _ => 2,
    };
    for i in used..7 {
        if p[i] != 0.0 {
            return false;
        }
    }
    match l.kind {
        LidNodeLayerKind::Surface => p[0] > 0.0 && p[1] < 1.0,
        LidNodeLayerKind::Media => p[0] > 0.0 && p[1] > 0.0 && p[1] <= 1.0 && p[3] < p[2] && p[2] < p[1],
        LidNodeLayerKind::Aggregate => p[0] > 0.0 && p[1] > 0.0 && p[1] <= 1.0,
        LidNodeLayerKind::Bottom => true,
    }
}

/// validateStack(): "" or the first problem with the stack.
pub fn validate_stack(stack: &[LidNodeLayer]) -> Bytes {
    let mut porous = false;
    for (i, l) in stack.iter().enumerate() {
        if !valid_layer(l) {
            return cat(&[b"invalid layer parameters at row ", (i + 1).to_string().as_bytes()]);
        }
        if l.kind == LidNodeLayerKind::Surface && i != 0 {
            return b"SURFACE must occur once, first".to_vec();
        }
        if l.kind == LidNodeLayerKind::Bottom && i + 1 != stack.len() {
            return b"BOTTOM must occur once, last".to_vec();
        }
        porous |= l.kind == LidNodeLayerKind::Media || l.kind == LidNodeLayerKind::Aggregate;
    }
    if !porous {
        return b"at least one MEDIA or AGGREGATE layer is required".to_vec();
    }
    if stack[0].flow.enabled {
        if stack[0].kind != LidNodeLayerKind::Surface {
            return b"Richards requires a SURFACE ponding layer".to_vec();
        }
        if !richards_valid_options(&stack[0].flow) {
            return b"invalid Richards numerical options".to_vec();
        }
        for l in stack {
            if l.kind == LidNodeLayerKind::Media || l.kind == LidNodeLayerKind::Aggregate {
                if !richards_valid_material(&l.retention, l.params[1]) {
                    return b"Richards requires valid retention and specific storage on every porous layer".to_vec();
                }
                if l.kind == LidNodeLayerKind::Media && l.retention.theta_r >= l.params[3] {
                    return b"Richards residual water content must be below the media wilting point".to_vec();
                }
            }
        }
    }
    Bytes::new()
}

/// thickness(): the summed thickness of every layer but BOTTOM.
pub fn thickness(stack: &[LidNodeLayer]) -> f64 {
    let mut d = 0.0;
    for l in stack {
        if l.kind != LidNodeLayerKind::Bottom {
            d += l.params[0];
        }
    }
    d
}

/// portOffset(): the link's port elevation above `node`'s invert (weir and rating-outlet
/// crests come from their subtype tables).
pub fn port_offset(ctx: &SimulationContext, link: i32, node: i32) -> f64 {
    let ul = link as usize;
    let mut offset = ctx.links.offset1[ul];
    let wr = ctx.link_subtypes.weir_row(link);
    let out = ctx.link_subtypes.outlet_row(link);
    if wr >= 0 {
        offset = ctx.link_subtypes.weirs.crest_height[wr as usize];
    } else if out >= 0 {
        offset = ctx.link_subtypes.outlets.crest_height[out as usize];
    } else {
        return if ctx.links.node1[ul] == node { offset } else { ctx.links.offset2[ul] };
    }
    // Keep a crest authored relative to its own node exact. Adding its datum
    // before subtracting the same datum can move an interface by one ULP.
    offset + (ctx.nodes.invert_elev[ctx.links.node1[ul] as usize] - ctx.nodes.invert_elev[node as usize])
}

// ============================================================================
// Section readers
// ============================================================================

/// lidnode::readNodes: [LID_NODES] `node control init_sat [BOTTOM CURRENT|SEALED|MESH_AQUIFER]`.
pub fn read_nodes(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let t = tokenizer::tokenize(line);
        if t.is_empty() {
            continue;
        }
        let mut sat = 0.0;
        let mut bottom = LidNodeBottom::Current;
        let mut valid_bottom = t.len() == 3;
        if t.len() == 5 && tokenizer::to_upper(&t[3]) == b"BOTTOM" {
            let mode = tokenizer::to_upper(&t[4]);
            valid_bottom = mode == b"CURRENT" || mode == b"SEALED" || mode == b"MESH_AQUIFER";
            if mode == b"SEALED" {
                bottom = LidNodeBottom::Sealed;
            }
            if mode == b"MESH_AQUIFER" {
                bottom = LidNodeBottom::MeshAquifer;
            }
        }
        if !valid_bottom || !number(&t[2], &mut sat) || sat < 0.0 || sat > 100.0 {
            ctx.errors.push(cat(&[b"Invalid [LID_NODES] row: ", line]));
            continue;
        }
        let n = ctx.node_names.find(&t[0]);
        let c = ctx.lid_names.find(&t[1]);
        if n < 0 || c < 0 {
            ctx.deferred_section_rows.push((b"LID_NODES".to_vec(), line.clone()));
            continue;
        }
        let r = ctx.node_subtypes.storage_row(n);
        if r < 0 {
            ctx.errors.push(cat(&[b"[LID_NODES] requires a storage node: ", &t[0]]));
            continue;
        }
        let config = &mut ctx.node_subtypes.storages.lid[r as usize];
        if config.control >= 0 {
            ctx.errors.push(cat(&[b"Duplicate [LID_NODES] assignment: ", &t[0]]));
            continue;
        }
        *config = LidNodeConfig { control: c, initial_saturation: sat, bottom };
    }
}

/// lidnode::readOutlets: [LID_NODE_OUTLETS] `link layer TOP|BOTTOM`.
pub fn read_outlets(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let t = tokenizer::tokenize(line);
        if t.is_empty() {
            continue;
        }
        let mut layer = 0.0;
        let pos: Bytes = if t.len() == 3 { tokenizer::to_upper(&t[2]) } else { Bytes::new() };
        if t.len() != 3
            || !number(&t[1], &mut layer)
            || layer < 1.0
            || layer > 2147483647.0
            || layer != layer.floor()
            || (pos != b"TOP" && pos != b"BOTTOM")
        {
            ctx.errors.push(cat(&[b"Invalid [LID_NODE_OUTLETS] row: ", line]));
            continue;
        }
        let link = ctx.link_names.find(&t[0]);
        if link < 0 {
            ctx.deferred_section_rows.push((b"LID_NODE_OUTLETS".to_vec(), line.clone()));
            continue;
        }
        let out = &mut ctx.lid_node_outlets;
        if out.iter().any(|a| a.link == link) {
            ctx.errors.push(cat(&[b"Duplicate LID outlet anchor: ", &t[0]]));
            continue;
        }
        out.push(LidNodeOutlet { link, layer: layer as i32, top: pos == b"TOP" });
    }
}

/// lidnode::validTreatment: false with `error` set when the rule cannot be used.
pub fn valid_treatment(ctx: &SimulationContext, rule: &LidLayerTreatment, error: &mut Bytes) -> bool {
    if ctx.pollutant_names.find(&rule.pollutant) < 0 {
        *error = cat(&[b"Unknown pollutant: ", &rule.pollutant]);
        return false;
    }
    if !rule.removal.is_finite() || rule.removal < 0.0 || rule.removal > 1.0 || !rule.decay.is_finite() || rule.decay < 0.0 {
        *error = b"Removal must be between 0 and 100 percent; decay must be finite and nonnegative".to_vec();
        return false;
    }
    if !rule.expression.is_empty() {
        let mut expr = crate::data::treatment_types::TreatExpr::default();
        let mut col = -1;
        // A std::stod throw inside parse cannot happen once validate has passed (see notes).
        if treatment::validate(&rule.expression, error, &mut col) != 0
            || treatment::parse_names(&rule.expression, &mut expr, ctx.pollutant_names.names()) != Ok(0)
        {
            if error.is_empty() {
                *error = b"Invalid treatment expression".to_vec();
            }
            return false;
        }
    }
    true
}

/// lidnode::readTreatment: [LID_LAYER_TREATMENT]
/// `control layer pollutant removal% decay/day [R= or C= expression | -]`.
pub fn read_treatment(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let mut stream = IStringStream::new(line);
        let mut control = Bytes::new();
        let mut pollutant = Bytes::new();
        let mut expression = Bytes::new();
        let mut layer: i32 = 0;
        let mut removal = 0.0;
        let mut decay = 0.0;
        stream.read_string(&mut control).read_int(&mut layer).read_string(&mut pollutant);
        stream.read_double(&mut removal).read_double(&mut decay);
        if !stream.ok() {
            ctx.errors.push(cat(&[b"Invalid [LID_LAYER_TREATMENT] row: ", line]));
            continue;
        }
        stream.ws().getline(&mut expression);
        let c = ctx.lid_names.find(&control);
        if c < 0
            || c >= ctx.lid_controls.node_layers.len() as i32
            || ctx.lid_controls.node_layers[c as usize].is_empty()
            || ctx.pollutant_names.find(&pollutant) < 0
        {
            ctx.deferred_section_rows.push((b"LID_LAYER_TREATMENT".to_vec(), line.clone()));
            continue;
        }
        if c < 0
            || c >= ctx.lid_controls.node_layers.len() as i32
            || layer < 1
            || layer > ctx.lid_controls.node_layers[c as usize].len() as i32
            || ctx.lid_controls.node_layers[c as usize][(layer - 1) as usize].kind == LidNodeLayerKind::Bottom
        {
            ctx.errors.push(cat(&[b"Invalid LID treatment control/layer: ", line]));
            continue;
        }
        let rule = LidLayerTreatment {
            pollutant: pollutant.clone(),
            removal: removal / 100.0,
            decay,
            expression: if expression == b"-" { Bytes::new() } else { expression },
        };
        let mut error = Bytes::new();
        if !valid_treatment(ctx, &rule, &mut error) {
            ctx.errors.push(cat(&[b"LID treatment: ", &error]));
            continue;
        }
        let p = ctx.pollutant_names.find(&pollutant);
        let rules = &ctx.lid_controls.node_layers[c as usize][(layer - 1) as usize].treatment;
        if rules.iter().any(|r| ctx.pollutant_names.find(&r.pollutant) == p) {
            ctx.errors.push(cat(&[b"Duplicate LID layer treatment: ", line]));
            continue;
        }
        ctx.lid_controls.node_layers[c as usize][(layer - 1) as usize].treatment.push(rule);
    }
}

/// lidnode::readRichards: [LID_RICHARDS]
/// `control OPTIONS cells_per_layer atol rtol max_step` or
/// `control layer theta_r alpha n l specific_storage`.
pub fn read_richards(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let t = tokenizer::tokenize(line);
        if t.is_empty() {
            continue;
        }
        let control = ctx.lid_names.find(&t[0]);
        if control < 0
            || control >= ctx.lid_controls.node_layers.len() as i32
            || ctx.lid_controls.node_layers[control as usize].is_empty()
        {
            ctx.deferred_section_rows.push((b"LID_RICHARDS".to_vec(), line.clone()));
            continue;
        }
        let mut input = IStringStream::new(line);
        let mut name = Bytes::new();
        let mut tag = Bytes::new();
        let mut trailing = Bytes::new();
        input.read_string(&mut name).read_string(&mut tag);
        if tokenizer::to_upper(&tag) == b"OPTIONS" {
            let mut o = richards::Options::default();
            o.enabled = true;
            let bad = {
                input.read_int(&mut o.cells_per_layer).read_double(&mut o.atol);
                input.read_double(&mut o.rtol).read_double(&mut o.max_step);
                !input.ok()
            } || input.read_string(&mut trailing).ok()
                || !richards_valid_options(&o)
                || ctx.lid_controls.node_layers[control as usize][0].flow.enabled;
            if bad {
                ctx.errors.push(cat(&[b"Invalid or duplicate [LID_RICHARDS] OPTIONS: ", line]));
                continue;
            }
            ctx.lid_controls.node_layers[control as usize][0].flow = o;
        } else {
            let stack = &ctx.lid_controls.node_layers[control as usize];
            let mut layer_stream = IStringStream::new(&tag);
            let mut layer: i32 = 0;
            let mut p = richards::Material::default();
            let bad = !layer_stream.read_int(&mut layer).ok()
                || layer_stream.read_string(&mut trailing).ok()
                || {
                    input.read_double(&mut p.theta_r).read_double(&mut p.alpha).read_double(&mut p.n);
                    input.read_double(&mut p.l).read_double(&mut p.specific_storage);
                    !input.ok()
                }
                || input.read_string(&mut trailing).ok()
                || layer < 1
                || layer > stack.len() as i32
                || (stack[(layer - 1) as usize].kind != LidNodeLayerKind::Media
                    && stack[(layer - 1) as usize].kind != LidNodeLayerKind::Aggregate)
                || !richards_valid_material(&p, stack[(layer - 1) as usize].params[1])
                || stack[(layer - 1) as usize].retention.alpha > 0.0;
            if bad {
                ctx.errors.push(cat(&[b"Invalid or duplicate [LID_RICHARDS] material: ", line]));
                continue;
            }
            ctx.lid_controls.node_layers[control as usize][(layer - 1) as usize].retention = p;
        }
    }
}

// ============================================================================
// lidnode::validate (called by PostParseResolver)
// ============================================================================

pub fn validate(ctx: &mut SimulationContext) {
    for c in 0..ctx.lid_controls.count() {
        if ctx.lid_controls.lid_type[c as usize] != b"NODE" {
            continue;
        }
        for layer in layers(ctx, c) {
            for rule in &layer.treatment {
                let mut diagnostic = Bytes::new();
                if !valid_treatment(ctx, rule, &mut diagnostic) {
                    ctx.errors.push(cat(&[b"LID treatment: ", &diagnostic]));
                }
            }
        }
        let error = validate_stack(&layers(ctx, c));
        if !error.is_empty() {
            let e = cat(&[b"LID control ", &ctx.lid_controls.names[c as usize], b": ", &error]);
            ctx.errors.push(e);
        }
    }
    let has_lid = ctx.node_subtypes.storages.lid.iter().any(|c| c.control >= 0);
    if has_lid && ctx.options.routing_model != RoutingModel::Dynwave {
        ctx.errors.push(b"LID storage nodes require FLOW_ROUTING DYNWAVE".to_vec());
    }
    if has_lid && ctx.n_pollutants() > 0 && ctx.options.quality_solver != QualitySolverKind::Legacy {
        ctx.errors.push(b"LID layer pollutant routing requires QUALITY_SOLVER LEGACY".to_vec());
    }
    for r in 0..ctx.node_subtypes.storages.count() {
        let cfg = ctx.node_subtypes.storages.lid[r as usize];
        if cfg.control < 0 {
            continue;
        }
        let n = ctx.node_subtypes.storages.node_idx[r as usize];
        let stack = layers(ctx, cfg.control);
        let error = validate_stack(&stack);
        if !error.is_empty() {
            let e = cat(&[b"LID node ", ctx.node_names.name_of(n), b": ", &error]);
            ctx.errors.push(e);
            continue;
        }
        let mesh_heat = ctx.options.heat_transport && cfg.bottom == LidNodeBottom::MeshAquifer;
        if stack[0].flow.enabled && (ctx.options.water_age || (ctx.options.heat_transport && !mesh_heat)) {
            let e = cat(&[
                b"Richards LID node ",
                ctx.node_names.name_of(n),
                b": water age and MSX require porous-cell transport adapters; heat is supported only for a qualified Mesh Richards/groundwater route",
            ]);
            ctx.errors.push(e);
        }
        // Physics review E9/E12 (2026-10-15): reaction species have no layer
        // adapter. They mixed as a plain tank that bypassed retained cells,
        // ports and layer treatment, and capture removed their water but not
        // their mass, so a decaying species could end above its source.
        else if ctx.reactions.configured
            && (!ctx.reactions.lid_enabled(cfg.control)
                || ctx.options.quality_solver != QualitySolverKind::Legacy
                || cfg.bottom == LidNodeBottom::MeshAquifer
                || ctx.reactions.species_is_wall.iter().any(|&v| v != 0))
        {
            let e = cat(&[
                b"LID node ",
                ctx.node_names.name_of(n),
                b": reaction species require an explicit [REACTION_LIDS] porous-cell transport adapter, BULK species, LEGACY quality and a non-Mesh bottom",
            ]);
            ctx.errors.push(e);
        }
        let d = thickness(&stack) / rain_depth(ctx);
        if (d - ctx.nodes.full_depth[n as usize]).abs() > 1.0e-7 * cmax(1.0, d) {
            let e = cat(&[b"LID node ", ctx.node_names.name_of(n), b": MaxDepth must equal the sum of layer thicknesses"]);
            ctx.errors.push(e);
        }
        if ctx.lid_controls.lid_type[cfg.control as usize] != b"NODE" {
            let surf = &ctx.lid_controls.surface[cfg.control as usize];
            let drain = &ctx.lid_controls.drain[cfg.control as usize];
            if surf[2] != 0.0 || surf[3] != 0.0 || surf[4] != 0.0 || drain.iter().any(|&v| v != 0.0) {
                let w = cat(&[
                    b"LID node ",
                    ctx.node_names.name_of(n),
                    b": surface roughness/slopes and DRAIN parameters are ignored; use hydraulic links",
                ]);
                ctx.warnings.push(w);
            }
        }
    }
    for k in 0..ctx.lid_node_outlets.len() {
        let anchor = ctx.lid_node_outlets[k];
        let n = outlet_node(ctx, anchor.link);
        if n < 0 {
            ctx.errors.push(b"LID outlet anchor requires exactly one LID endpoint".to_vec());
            continue;
        }
        let r = ctx.node_subtypes.storage_row(n);
        let h = anchor_height(&layers(ctx, ctx.node_subtypes.storages.lid[r as usize].control), anchor.layer, anchor.top);
        if h < 0.0 {
            ctx.errors.push(b"LID outlet anchor layer is outside the stack".to_vec());
            continue;
        }
        let offset = port_offset(ctx, anchor.link, n);
        if (offset - h / rain_depth(ctx)).abs() > 1.0e-7 {
            let w = cat(&[
                b"LID outlet ",
                ctx.link_names.name_of(anchor.link),
                b": offset does not match the anchored layer interface",
            ]);
            ctx.warnings.push(w);
        }
    }
}
