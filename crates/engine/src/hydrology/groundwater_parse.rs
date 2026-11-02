// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/Groundwater.hpp and Groundwater.cpp
// (GW_VAR_NAMES and the [GWF] expression check gwf_validate only), with the parts of
// src/engine/math/MathExpr.hpp and MathExpr.cpp that gwf_validate's drift guard runs
// (the function-name list, tokenize and the shunting-yard parse)
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// What HydrologyHandler's handle_gwf calls, and nothing else of Groundwater.cpp (the GW
// solver is not here). Upstream signals a syntax failure by throwing GwfFail and catches
// it in gwf_validate; here the parser returns Result<_, GwfFail>. Upstream's catch(...)
// also turns any exception from the drift guard (std::stod's invalid_argument /
// out_of_range inside mathexpr's tokenizer) into "internal: expression rejected by the
// engine parser"; `stod` below reports those two cases as an Err.
//
// The character classes are the C locale's: isspace is space, \t, \n, \v, \f, \r; isdigit,
// isalpha and isalnum are ASCII only (mathexpr calls them with a plain, signed char; in
// glibc's C locale a negative value is no class either, so the two agree).
//
// strtod/stod on a scanned number: the word holds only digits, '.', 'e', 'E' and a sign
// right after an exponent letter, and starts with a digit or '.', so strtod reads exactly
// the decimal prefix std::from_chars reads (core::charconv), and the full-consumption test
// is the same. std::stod throws out_of_range when strtod sets ERANGE: on overflow, and on
// underflow (glibc: a result that is tiny and inexact). Here underflow is a nonzero literal
// that rounds to zero or to a subnormal. A subnormal is exact only for a literal of at
// least 700 significant digits; that case, and the half-ulp window just below DBL_MIN
// that rounds up to DBL_MIN, are not distinguished.

use crate::core::charconv::{from_chars_double, Errc};
use crate::core::options::Bytes;

/// Groundwater.hpp GW_VAR_NAMES.
pub const GW_VAR_NAMES: [&[u8]; 11] = [b"HGW", b"HSW", b"HCB", b"HGS", b"KS", b"K", b"THETA", b"PHI", b"FI", b"FU", b"A"];

/// MathExpr.cpp function_names(): the func_map keys in TokenType order.
pub const MATHEXPR_FUNCTION_NAMES: [&[u8]; 21] = [
    b"abs", b"sgn", b"sqrt", b"log", b"exp", b"sin", b"cos", b"tan", b"asin", b"acos", b"atan", b"step", b"min",
    b"max", b"cot", b"sinh", b"cosh", b"tanh", b"coth", b"log10", b"acot",
];

#[inline]
fn c_isspace(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

fn cat(parts: &[&[u8]]) -> Bytes {
    parts.concat()
}

// ============================================================================
// [GWF] expression validation — diagnostic peer of mathexpr::parse
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GwfTok {
    Num,
    Var,
    Func,
    Op,
    LParen,
    RParen,
    Comma,
    End,
}

struct GwfToken {
    kind: GwfTok,
    /// as typed (identifiers keep the user's case)
    text: Bytes,
    col: i32,
}

struct GwfFail {
    msg: Bytes,
    col: i32,
}

fn fail(msg: Bytes, col: i32) -> GwfFail {
    GwfFail { msg, col }
}

/// gwf_is_variable: some whole GW_VAR_NAMES entry is a prefix of `upper`.
fn gwf_is_variable(upper: &[u8]) -> bool {
    for v in GW_VAR_NAMES.iter() {
        let mut k = 0usize;
        while k < v.len() && k < upper.len() && upper[k] == v[k] {
            k += 1;
        }
        if k == v.len() {
            return true;
        }
    }
    false
}

/// gwf_is_function: exact lowercase name in mathexpr::function_names().
fn gwf_is_function(lower: &[u8]) -> bool {
    MATHEXPR_FUNCTION_NAMES.iter().any(|f| lower == *f)
}

/// The scanned-number word's length as strtod consumes it.
fn strtod_used(word: &[u8]) -> usize {
    match from_chars_double(word) {
        (Err(Errc::InvalidArgument), _) => 0,
        (_, used) => used,
    }
}

/// Whether a scanned number ends (strict tokenizer and mathexpr scanner alike).
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

fn gwf_tokenize(s: &[u8]) -> Result<Vec<GwfToken>, GwfFail> {
    let mut out: Vec<GwfToken> = Vec::new();
    let mut i = 0usize;
    while i < s.len() {
        let c = s[i];
        if c_isspace(c) {
            i += 1;
            continue;
        }
        let col = i as i32;

        if c.is_ascii_digit() || c == b'.' {
            let start = i;
            i = scan_number(s, i);
            let word = s[start..i].to_vec();
            let used = strtod_used(&word);
            if used == 0 || used != word.len() {
                return Err(fail(cat(&[b"malformed number '", &word, b"'"]), col));
            }
            out.push(GwfToken { kind: GwfTok::Num, text: word, col });
            continue;
        }

        if c.is_ascii_alphabetic() || c == b'_' {
            let start = i;
            while i < s.len() && (s[i].is_ascii_alphanumeric() || s[i] == b'_') {
                i += 1;
            }
            let word = s[start..i].to_vec();
            let lower = word.to_ascii_lowercase();
            let upper = word.to_ascii_uppercase();
            if gwf_is_function(&lower) {
                out.push(GwfToken { kind: GwfTok::Func, text: word, col });
            } else if gwf_is_variable(&upper) {
                out.push(GwfToken { kind: GwfTok::Var, text: word, col });
            } else {
                return Err(fail(cat(&[b"unknown variable '", &word, b"'"]), col));
            }
            continue;
        }

        match c {
            b'+' | b'-' | b'*' | b'/' | b'^' => out.push(GwfToken { kind: GwfTok::Op, text: vec![c], col }),
            b'(' => out.push(GwfToken { kind: GwfTok::LParen, text: b"(".to_vec(), col }),
            b')' => out.push(GwfToken { kind: GwfTok::RParen, text: b")".to_vec(), col }),
            b',' => out.push(GwfToken { kind: GwfTok::Comma, text: b",".to_vec(), col }),
            _ => return Err(fail(cat(&[b"unexpected character '", &[c], b"'"]), col)),
        }
        i += 1;
    }
    out.push(GwfToken { kind: GwfTok::End, text: Vec::new(), col: s.len() as i32 });
    Ok(out)
}

// expr   := term (('+'|'-') term)*
// term   := unary (('*'|'/') unary)*
// unary  := '-' unary | power
// power  := primary ('^' unary)?
// primary:= NUM | VAR | FUNC '(' expr (',' expr)* ')' | '(' expr ')'
struct GwfParser<'a> {
    t: &'a [GwfToken],
    p: usize,
    paren_depth: i32,
    func_depth: i32,
}

impl<'a> GwfParser<'a> {
    fn cur(&self) -> &'a GwfToken {
        &self.t[self.p]
    }
    fn is_op(&self, s: &[u8]) -> bool {
        self.cur().kind == GwfTok::Op && self.cur().text == s
    }
    fn expr(&mut self) -> Result<(), GwfFail> {
        self.term()?;
        while self.is_op(b"+") || self.is_op(b"-") {
            self.p += 1;
            self.term()?;
        }
        Ok(())
    }
    fn term(&mut self) -> Result<(), GwfFail> {
        self.unary()?;
        while self.is_op(b"*") || self.is_op(b"/") {
            self.p += 1;
            self.unary()?;
        }
        Ok(())
    }
    fn unary(&mut self) -> Result<(), GwfFail> {
        if self.is_op(b"-") {
            self.p += 1;
            return self.unary();
        }
        self.power()
    }
    fn power(&mut self) -> Result<(), GwfFail> {
        self.primary()?;
        if self.is_op(b"^") {
            self.p += 1;
            self.unary()?;
        }
        Ok(())
    }
    fn primary(&mut self) -> Result<(), GwfFail> {
        let k = self.cur();
        match k.kind {
            GwfTok::Num | GwfTok::Var => {
                self.p += 1;
                Ok(())
            }
            GwfTok::Func => {
                self.p += 1;
                if self.cur().kind != GwfTok::LParen {
                    return Err(fail(cat(&[b"function '", &k.text, b"' must be followed by '('"]), self.cur().col));
                }
                let open = self.cur();
                self.p += 1;
                self.paren_depth += 1;
                self.func_depth += 1;
                let mut nargs = 1;
                self.expr()?;
                while self.cur().kind == GwfTok::Comma {
                    self.p += 1;
                    nargs += 1;
                    self.expr()?;
                }
                self.close_paren(open)?;
                self.paren_depth -= 1;
                self.func_depth -= 1;
                let lower = k.text.to_ascii_lowercase();
                let want = if lower == b"min" || lower == b"max" { 2 } else { 1 };
                if nargs != want {
                    let n = want.to_string();
                    let unit: &[u8] = if want == 1 { b" argument" } else { b" arguments" };
                    return Err(fail(cat(&[b"function '", &k.text, b"' takes ", n.as_bytes(), unit]), k.col));
                }
                Ok(())
            }
            GwfTok::LParen => {
                let open = self.cur();
                self.p += 1;
                self.paren_depth += 1;
                self.expr()?;
                if self.cur().kind == GwfTok::Comma {
                    return Err(fail(b"unexpected ','".to_vec(), self.cur().col));
                }
                self.close_paren(open)?;
                self.paren_depth -= 1;
                Ok(())
            }
            GwfTok::Op => Err(fail(cat(&[b"missing operand before '", &k.text, b"'"]), k.col)),
            GwfTok::RParen => {
                if self.paren_depth == 0 {
                    return Err(fail(b"unbalanced parenthesis".to_vec(), k.col));
                }
                Err(fail(b"missing operand before ')'".to_vec(), k.col))
            }
            GwfTok::Comma => {
                if self.func_depth == 0 {
                    return Err(fail(b"unexpected ','".to_vec(), k.col));
                }
                Err(fail(b"missing operand before ','".to_vec(), k.col))
            }
            GwfTok::End => {
                let prev = &self.t[self.p - 1];
                if prev.kind == GwfTok::LParen {
                    return Err(fail(b"unbalanced parenthesis".to_vec(), prev.col));
                }
                Err(fail(cat(&[b"missing operand after '", &prev.text, b"'"]), prev.col))
            }
        }
    }
    // After an expression inside '(': the next token must close it.
    fn close_paren(&mut self, open: &GwfToken) -> Result<(), GwfFail> {
        if self.cur().kind == GwfTok::RParen {
            self.p += 1;
            return Ok(());
        }
        if self.cur().kind == GwfTok::End {
            return Err(fail(b"unbalanced parenthesis".to_vec(), open.col));
        }
        Err(fail(cat(&[b"missing operator before '", &self.cur().text, b"'"]), self.cur().col))
    }
}

/// gwf_validate: 0 if `expr` is a well-formed [GWF] expression, else -1 with `msg` and the
/// zero-based byte column `col` (-1 when the drift guard rejected it).
pub fn gwf_validate(expr: &[u8], msg: &mut Bytes, col: &mut i32) -> i32 {
    msg.clear();
    *col = -1;
    let result = (|| -> Result<(), GwfFail> {
        let tokens = gwf_tokenize(expr)?;
        if tokens[0].kind == GwfTok::End {
            return Err(fail(b"expression is empty".to_vec(), 0));
        }
        let mut parser = GwfParser { t: &tokens, p: 0, paren_depth: 0, func_depth: 0 };
        parser.expr()?;
        let rest = parser.cur();
        if rest.kind != GwfTok::End {
            if rest.kind == GwfTok::RParen {
                return Err(fail(b"unbalanced parenthesis".to_vec(), rest.col));
            }
            if rest.kind == GwfTok::Comma {
                return Err(fail(b"unexpected ','".to_vec(), rest.col));
            }
            return Err(fail(cat(&[b"missing operator before '", &rest.text, b"'"]), rest.col));
        }
        // Drift guard: the production parser must accept what we accept.
        match mathexpr_probe::parse(expr) {
            Ok(probe) if probe.valid => Ok(()),
            // a nonzero return, an invalid probe, or (upstream's catch(...)) a stod exception
            _ => Err(fail(b"internal: expression rejected by the engine parser".to_vec(), -1)),
        }
    })();
    match result {
        Ok(()) => 0,
        Err(f) => {
            *msg = f.msg;
            *col = f.col;
            -1
        }
    }
}

/// MathExpr.cpp tokenize + parse, as far as gwf_validate's drift guard uses them: the
/// postfix is built (and discarded by the caller); evaluation is not here.
mod mathexpr_probe {
    use super::{c_isspace, scan_number, MATHEXPR_FUNCTION_NAMES};
    use crate::core::charconv::{from_chars_double, Errc};

    /// MathExpr.hpp TokenType (the values matter: isFunction is `>= FUNC_ABS`).
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    #[repr(i32)]
    enum TokenType {
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
        FuncCot = 25,
        FuncSinh = 26,
        FuncCosh = 27,
        FuncTanh = 28,
        FuncCoth = 29,
        FuncLog10 = 30,
        FuncAcot = 31,
    }

    const FUNC_TYPES: [TokenType; 21] = [
        TokenType::FuncAbs,
        TokenType::FuncSgn,
        TokenType::FuncSqrt,
        TokenType::FuncLog,
        TokenType::FuncExp,
        TokenType::FuncSin,
        TokenType::FuncCos,
        TokenType::FuncTan,
        TokenType::FuncAsin,
        TokenType::FuncAcos,
        TokenType::FuncAtan,
        TokenType::FuncStep,
        TokenType::FuncMin,
        TokenType::FuncMax,
        TokenType::FuncCot,
        TokenType::FuncSinh,
        TokenType::FuncCosh,
        TokenType::FuncTanh,
        TokenType::FuncCoth,
        TokenType::FuncLog10,
        TokenType::FuncAcot,
    ];

    #[derive(Debug, Clone)]
    struct Token {
        type_: TokenType,
        #[allow(dead_code)]
        value: f64,
        #[allow(dead_code)]
        var_name: Vec<u8>,
    }

    impl Token {
        fn of(type_: TokenType) -> Token {
            Token { type_, value: 0.0, var_name: Vec::new() }
        }
    }

    pub struct Expression {
        #[allow(dead_code)]
        postfix: Vec<Token>,
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
    fn is_function(t: TokenType) -> bool {
        t as i32 >= TokenType::FuncAbs as i32
    }
    fn is_operator(t: TokenType) -> bool {
        matches!(t, TokenType::Add | TokenType::Sub | TokenType::Mul | TokenType::Div | TokenType::Pow | TokenType::Neg)
    }
    fn is_right_assoc(t: TokenType) -> bool {
        t == TokenType::Pow || t == TokenType::Neg
    }

    /// std::stod on a scanned word: Err where it throws (no conversion: invalid_argument;
    /// ERANGE: out_of_range). See the module notes on underflow.
    fn stod(word: &[u8]) -> Result<f64, ()> {
        let (r, used) = from_chars_double(word);
        match r {
            Err(Errc::InvalidArgument) => Err(()),
            Err(Errc::ResultOutOfRange) => Err(()),
            Ok(v) => {
                let _ = used;
                if v != 0.0 && v.abs() < f64::MIN_POSITIVE {
                    return Err(()); // subnormal: tiny and (see notes) inexact
                }
                Ok(v)
            }
        }
    }

    fn tokenize(s: &[u8]) -> Result<Vec<Token>, ()> {
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
                i = scan_number(s, i);
                let mut t = Token::of(TokenType::Number);
                t.value = stod(&s[start..i])?;
                tokens.push(t);
            } else if c.is_ascii_alphabetic() || c == b'_' {
                let start = i;
                while i < s.len() && (s[i].is_ascii_alphanumeric() || s[i] == b'_') {
                    i += 1;
                }
                let word = &s[start..i];
                let lower = word.to_ascii_lowercase();
                match MATHEXPR_FUNCTION_NAMES.iter().position(|f| *f == lower.as_slice()) {
                    Some(k) => tokens.push(Token::of(FUNC_TYPES[k])),
                    None => {
                        let mut t = Token::of(TokenType::Variable);
                        t.var_name = word.to_vec();
                        tokens.push(t);
                    }
                }
            } else {
                let type_ = match c {
                    b'+' => TokenType::Add,
                    b'-' => {
                        // Unary minus if at start or after operator/lparen
                        match tokens.last() {
                            None => TokenType::Neg,
                            Some(b) if b.type_ == TokenType::LParen || is_operator(b.type_) || b.type_ == TokenType::Comma => {
                                TokenType::Neg
                            }
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
                tokens.push(Token::of(type_));
                i += 1;
            }
        }
        Ok(tokens)
    }

    /// mathexpr::parse. Ok(expression) for upstream's 0 or -1 return (valid tells them
    /// apart); Err where upstream's tokenizer throws.
    pub fn parse(expr_str: &[u8]) -> Result<Expression, ()> {
        let mut result = Expression { postfix: Vec::new(), valid: false };
        let tokens = tokenize(expr_str)?;
        let mut ops: Vec<Token> = Vec::new();
        for tok in &tokens {
            if tok.type_ == TokenType::Number || tok.type_ == TokenType::Variable {
                result.postfix.push(tok.clone());
            } else if is_function(tok.type_) {
                ops.push(tok.clone());
            } else if tok.type_ == TokenType::Comma {
                // unreachable in upstream as well: COMMA's ordinal makes is_function true
                while let Some(top) = ops.last() {
                    if top.type_ == TokenType::LParen {
                        break;
                    }
                    result.postfix.push(ops.pop().unwrap());
                }
            } else if is_operator(tok.type_) {
                while let Some(top) = ops.last() {
                    let go = is_operator(top.type_)
                        && ((is_right_assoc(tok.type_) && precedence(tok.type_) < precedence(top.type_))
                            || (!is_right_assoc(tok.type_) && precedence(tok.type_) <= precedence(top.type_)));
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
                    return Ok(result); // mismatched parens: -1
                }
                ops.pop(); // pop LPAREN
                if let Some(top) = ops.last() {
                    if is_function(top.type_) {
                        result.postfix.push(ops.pop().unwrap());
                    }
                }
            }
        }
        while let Some(top) = ops.pop() {
            if top.type_ == TokenType::LParen {
                return Ok(result);
            }
            result.postfix.push(top);
        }
        result.valid = true;
        Ok(result)
    }
}
