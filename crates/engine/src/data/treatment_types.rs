// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/quality/Treatment.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The types of upstream's namespace openswmm::treatment: the token and variable enums, the
// postfix Token and the parsed TreatExpr. The functions Treatment.hpp declares are defined
// in Treatment.cpp and are translated with it; each is named at the end of this file.

/// Expression token types (`enum class TokenType : int`).
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
    /// Unary negation.
    Neg = 7,
    FuncExp = 8,
    FuncLog = 9,
    FuncSqrt = 10,
    FuncMin = 11,
    FuncMax = 12,
    FuncAbs = 13,
    FuncSgn = 14,
    /// Heaviside step.
    FuncStep = 15,
    /// Left parenthesis (parser internal, never in output).
    Lparen = 16,
    /// Right parenthesis (parser internal, never in output).
    Rparen = 17,
    /// Comma separator for min/max (parser internal).
    Comma = 18,
}

/// Variables a treatment expression can reference (`enum class TreatVar : int`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum TreatVar {
    /// Current concentration (this pollutant).
    C = 0,
    /// Removal fraction (this pollutant, unused in expressions).
    R = 1,
    /// Timestep (sec).
    Dt = 2,
    /// Hydraulic residence time (hours, per legacy convention).
    Hrt = 3,
    /// Flow rate in user units (Q*UCF(FLOW), matching legacy pvFLOW).
    Q = 4,
    /// Volume (ft3 or m3).
    V = 5,
    /// Depth in user units (y*UCF(LENGTH), matching legacy pvDEPTH).
    D = 6,
    /// Concentration of another pollutant (uses pollut_ref).
    CPollut = 7,
    /// Removal fraction of another pollutant (uses pollut_ref, co-treatment).
    RPollut = 8,
    /// Node surface area in user units squared ((a1+a2)/2*UCF(LENGTH)^2, legacy pvAREA).
    Area = 9,
}

/// One token of a postfix expression.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Token {
    pub type_: TokenType,
    /// For NUMBER tokens.
    pub value: f64,
    /// For VARIABLE tokens.
    pub var: TreatVar,
    /// Pollutant index for C_POLLUT/R_POLLUT.
    pub pollut_ref: i32,
}

impl Default for Token {
    fn default() -> Self {
        Token { type_: TokenType::Number, value: 0.0, var: TreatVar::C, pollut_ref: -1 }
    }
}

/// A parsed treatment expression, ready to evaluate.
#[derive(Debug, Clone, PartialEq)]
pub struct TreatExpr {
    /// Postfix (RPN) token list.
    pub tokens: Vec<Token>,
    /// Which pollutant this expression applies to.
    pub pollutant_idx: i32,
    /// True if the expression computes R (removal), false if C.
    pub is_removal: bool,
}

impl Default for TreatExpr {
    fn default() -> Self {
        TreatExpr { tokens: Vec::new(), pollutant_idx: -1, is_removal: false }
    }
}

// Functions declared in Treatment.hpp, defined in Treatment.cpp and translated with it:
//
//   int parse(const std::string& expr_str, TreatExpr& result);
//   int parse(const std::string& expr_str, TreatExpr& result,
//             int (*pollut_lookup)(const std::string& name));
//   int validate(const std::string& expr_str, std::string& msg, int& col);
//   int parse(const std::string& expr_str, TreatExpr& result,
//             const std::vector<std::string>& pollut_names);
//   double evaluate(const TreatExpr& expr, double c, double dt, double hrt, double q,
//                   double v, double d, double area = 0.0);
//   double evaluate(const TreatExpr& expr, double c, double dt, double hrt, double q,
//                   double v, double d, const double* cin, const double* removal,
//                   int n_pollut, double area = 0.0, const double* cpollut = nullptr);
//   double applyTreatment(const TreatExpr& expr, double c_in, double dt, double hrt,
//                         double q, double v, double d, double area = 0.0);
