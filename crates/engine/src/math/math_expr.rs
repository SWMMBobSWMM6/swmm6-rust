// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/math/MathExpr.hpp and MathExpr.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The expression types (TokenType, Token, Expression), the shunting-yard parser, the
// function-name list, bind_variables and compute_max_stack_depth. SWMMEngine::initHydrology
// compiles the [GWF] expressions with parse + bind_variables.
//
// Differences of form:
//   * tokenize throws std::invalid_argument / std::out_of_range from std::stod on a bad
//     number; here tokenize and parse_stod return Err(StodError) naming the exception, and
//     parse returns Err(()), where upstream throws (nothing on upstream's paths catches
//     it, see UPSTREAM_NOTES U1; ControlEngine::addExpression in initHydraulics reaches
//     std::terminate this way). std::stod is read as
//     hydrology/groundwater_parse.rs reads it (from_chars on the scanned word; ERANGE on
//     overflow and on a nonzero literal that rounds to zero or a subnormal).
//   * function_names() builds its list once from an unordered_map sorted by TokenType;
//     the order is the TokenType order, so the list is a constant here.
//   * The character classes are the C locale's (ASCII), as in groundwater_parse.rs.
//   * NOT YET TRANSLATED (run-time step functions, phase 1c): both evaluate overloads and
//     evaluate_fast.

use crate::core::charconv::{from_chars_double, Errc};
use crate::core::cstdlib::c_isspace;
use crate::core::options::Bytes;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum TokenType {
    Number = 0,
    Variable = 1,
    Add = 2,
    Sub = 3,
    Mul = 4,
    Div = 5,
    Pow = 6,
    Neg = 7,
    LParen = 8,
    RParen = 9,
    FuncAbs = 10,
    FuncSgn = 11,
    FuncSqrt = 12,
    FuncLog = 13,
    FuncExp = 14,
    FuncSin = 15,
    FuncCos = 16,
    FuncTan = 17,
    FuncAsin = 18,
    FuncAcos = 19,
    FuncAtan = 20,
    FuncStep = 21,
    FuncMin = 22,
    FuncMax = 23,
    Comma = 24,
    // Legacy mathexpr.c functions not in original refactor (Gap 49)
    FuncCot = 25,
    FuncSinh = 26,
    FuncCosh = 27,
    FuncTanh = 28,
    FuncCoth = 29,
    FuncLog10 = 30,
    FuncAcot = 31,
}

#[derive(Debug, Clone)]
pub struct Token {
    pub type_: TokenType,
    /// For NUMBER tokens
    pub value: f64,
    /// For VARIABLE tokens (index into variable table)
    pub var_idx: i32,
    /// For VARIABLE tokens (name for lookup)
    pub var_name: Bytes,
}

impl Default for Token {
    fn default() -> Self {
        Token { type_: TokenType::Number, value: 0.0, var_idx: -1, var_name: Bytes::new() }
    }
}

/// The exception std::stod would throw.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StodError {
    /// std::invalid_argument (nothing converts)
    InvalidArgument,
    /// std::out_of_range (strtod reports ERANGE)
    OutOfRange,
}

/// Compiled expression (postfix token list).
#[derive(Debug, Clone, Default)]
pub struct Expression {
    /// Postfix (RPN) token sequence
    pub postfix: Vec<Token>,
    pub valid: bool,
}

fn precedence(t: TokenType) -> i32 {
    match t {
        TokenType::Add | TokenType::Sub => 1,
        TokenType::Mul | TokenType::Div => 2,
        TokenType::Pow => 3,
        TokenType::Neg => 4,
        _ => 0,
    }
}

/// Ordinal >= FUNC_ABS (COMMA included, as upstream's numeric test includes it).
fn is_function(t: TokenType) -> bool {
    t as i32 >= TokenType::FuncAbs as i32
}

fn is_operator(t: TokenType) -> bool {
    matches!(t, TokenType::Add | TokenType::Sub | TokenType::Mul | TokenType::Div | TokenType::Pow | TokenType::Neg)
}

fn is_right_assoc(t: TokenType) -> bool {
    t == TokenType::Pow || t == TokenType::Neg
}

/// func_map, in TokenType order (which is function_names()'s order).
const FUNCS: [(&[u8], TokenType); 21] = [
    (b"abs", TokenType::FuncAbs),
    (b"sgn", TokenType::FuncSgn),
    (b"sqrt", TokenType::FuncSqrt),
    (b"log", TokenType::FuncLog),
    (b"exp", TokenType::FuncExp),
    (b"sin", TokenType::FuncSin),
    (b"cos", TokenType::FuncCos),
    (b"tan", TokenType::FuncTan),
    (b"asin", TokenType::FuncAsin),
    (b"acos", TokenType::FuncAcos),
    (b"atan", TokenType::FuncAtan),
    (b"step", TokenType::FuncStep),
    (b"min", TokenType::FuncMin),
    (b"max", TokenType::FuncMax),
    (b"cot", TokenType::FuncCot),
    (b"sinh", TokenType::FuncSinh),
    (b"cosh", TokenType::FuncCosh),
    (b"tanh", TokenType::FuncTanh),
    (b"coth", TokenType::FuncCoth),
    (b"log10", TokenType::FuncLog10),
    (b"acot", TokenType::FuncAcot),
];

/// function_names(): the built-in names ordered by TokenType.
pub fn function_names() -> Vec<Bytes> {
    FUNCS.iter().map(|(n, _)| n.to_vec()).collect()
}

/// std::stod on a scanned word: Err where it throws.
fn stod(word: &[u8]) -> Result<f64, StodError> {
    let (r, _used) = from_chars_double(word);
    match r {
        Err(Errc::InvalidArgument) => Err(StodError::InvalidArgument),
        Err(Errc::ResultOutOfRange) => Err(StodError::OutOfRange),
        Ok(v) => {
            if v != 0.0 && v.abs() < f64::MIN_POSITIVE {
                return Err(StodError::OutOfRange); // subnormal: tiny and inexact (ERANGE)
            }
            Ok(v)
        }
    }
}

fn tokenize(s: &[u8]) -> Result<Vec<Token>, StodError> {
    let mut tokens: Vec<Token> = Vec::new();
    let mut i = 0usize;
    while i < s.len() {
        let c = s[i];
        if c_isspace(c) {
            i += 1;
            continue;
        }
        if c.is_ascii_digit() || c == b'.' {
            let start = i;
            while i < s.len() && (s[i].is_ascii_digit() || s[i] == b'.' || s[i] == b'e' || s[i] == b'E' || ((s[i] == b'+' || s[i] == b'-') && i > 0 && (s[i - 1] == b'e' || s[i - 1] == b'E'))) {
                i += 1;
            }
            let t = Token { type_: TokenType::Number, value: stod(&s[start..i])?, ..Token::default() };
            tokens.push(t);
        } else if c.is_ascii_alphabetic() || c == b'_' {
            let start = i;
            while i < s.len() && (s[i].is_ascii_alphanumeric() || s[i] == b'_') {
                i += 1;
            }
            let word = &s[start..i];
            let lower = word.to_ascii_lowercase();
            match FUNCS.iter().find(|(n, _)| *n == lower.as_slice()) {
                Some((_, ty)) => tokens.push(Token { type_: *ty, ..Token::default() }),
                None => tokens.push(Token { type_: TokenType::Variable, var_name: word.to_vec(), ..Token::default() }),
            }
        } else {
            let type_ = match c {
                b'+' => TokenType::Add,
                b'-' => {
                    // Unary minus if at start or after operator/lparen
                    match tokens.last() {
                        None => TokenType::Neg,
                        Some(b) if b.type_ == TokenType::LParen || is_operator(b.type_) || b.type_ == TokenType::Comma => TokenType::Neg,
                        Some(_) => TokenType::Sub,
                    }
                }
                b'*' => TokenType::Mul,
                b'/' => TokenType::Div,
                b'^' => TokenType::Pow,
                b'(' => TokenType::LParen,
                b')' => TokenType::RParen,
                b',' => TokenType::Comma,
                _ => {
                    i += 1;
                    continue;
                }
            };
            tokens.push(Token { type_, ..Token::default() });
            i += 1;
        }
    }
    Ok(tokens)
}

/// mathexpr::parse: Ok(0) and `result.valid` on success, Ok(-1) on mismatched
/// parentheses (with the postfix built so far), Err(()) where upstream's tokenizer
/// throws (after `result` was cleared).
pub fn parse(expr_str: &[u8], result: &mut Expression) -> Result<i32, ()> {
    parse_stod(expr_str, result).map_err(|_| ())
}

/// parse, with the exception std::stod would throw as the error.
pub fn parse_stod(expr_str: &[u8], result: &mut Expression) -> Result<i32, StodError> {
    result.postfix.clear();
    result.valid = false;

    let tokens = tokenize(expr_str)?;
    let mut ops: Vec<Token> = Vec::new();

    for tok in &tokens {
        if tok.type_ == TokenType::Number || tok.type_ == TokenType::Variable {
            result.postfix.push(tok.clone());
        } else if is_function(tok.type_) {
            ops.push(tok.clone());
        } else if tok.type_ == TokenType::Comma {
            // unreachable, as upstream: COMMA's ordinal makes is_function true
            while let Some(top) = ops.last() {
                if top.type_ == TokenType::LParen {
                    break;
                }
                result.postfix.push(ops.pop().unwrap());
            }
        } else if is_operator(tok.type_) {
            while let Some(top) = ops.last() {
                let go = is_operator(top.type_) && ((is_right_assoc(tok.type_) && precedence(tok.type_) < precedence(top.type_)) || (!is_right_assoc(tok.type_) && precedence(tok.type_) <= precedence(top.type_)));
                if !go {
                    break;
                }
                result.postfix.push(ops.pop().unwrap());
            }
            ops.push(tok.clone());
        } else if tok.type_ == TokenType::LParen {
            ops.push(tok.clone());
        } else if tok.type_ == TokenType::RParen {
            while let Some(top) = ops.last() {
                if top.type_ == TokenType::LParen {
                    break;
                }
                result.postfix.push(ops.pop().unwrap());
            }
            if ops.is_empty() {
                return Ok(-1); // mismatched parens
            }
            ops.pop(); // pop LPAREN
            if let Some(top) = ops.last() {
                if is_function(top.type_) {
                    result.postfix.push(ops.pop().unwrap());
                }
            }
        }
    }

    while let Some(top) = ops.last() {
        if top.type_ == TokenType::LParen {
            return Ok(-1);
        }
        result.postfix.push(ops.pop().unwrap());
    }

    result.valid = true;
    Ok(0)
}

/// Bind every VARIABLE token to the first table entry that is a case-insensitive prefix
/// of its name (legacy findmatch()/match()); returns the number bound. Unmatched tokens
/// get var_idx -1. Table names are C strings: each ends at its first NUL.
pub fn bind_variables(expr: &mut Expression, name_table: &[&[u8]], n_vars: i32) -> i32 {
    let mut bound = 0;
    for tok in expr.postfix.iter_mut() {
        if tok.type_ != TokenType::Variable {
            continue;
        }
        tok.var_idx = -1;
        for v in 0..n_vars {
            let tbl = crate::core::cstdlib::c_str(name_table[v as usize]);
            let tn = &tok.var_name;
            let mut matched = true;
            let mut k = 0usize;
            while k < tn.len() && k < tbl.len() {
                if tn[k].to_ascii_uppercase() != tbl[k].to_ascii_uppercase() {
                    matched = false;
                    break;
                }
                k += 1;
            }
            if matched && k == tbl.len() {
                // whole table name is a prefix of the token
                tok.var_idx = v;
                bound += 1;
                break;
            }
        }
    }
    bound
}

pub fn compute_max_stack_depth(expr: &Expression) -> i32 {
    let mut depth = 0;
    let mut max_depth = 0;
    for tok in &expr.postfix {
        match tok.type_ {
            TokenType::Number | TokenType::Variable => depth += 1,
            TokenType::Neg => {} // pop 1 push 1 = net 0
            TokenType::Add | TokenType::Sub | TokenType::Mul | TokenType::Div | TokenType::Pow | TokenType::FuncMin | TokenType::FuncMax => depth -= 1,
            _ => {} // unary functions: pop 1 push 1 = net 0
        }
        if depth > max_depth {
            max_depth = depth;
        }
    }
    max_depth
}
