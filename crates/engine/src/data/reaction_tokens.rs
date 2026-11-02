// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/ReactionTokens.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// POD reaction-VM instruction and span types (phase R2, D-L3 flat pool). Upstream nests
// an unscoped `enum Op : uint8_t` in RxToken and stores the field as a plain uint8_t,
// comparing it as `token.op == RxToken::PUSH_TERM`. Here the enum is RxTokenOp (CamelCase
// variants, same u8 discriminants), the field stays u8, and RxToken also carries the
// opcodes as associated u8 consts under upstream's spelling so `token.op ==
// RxToken::PUSH_TERM` and `match token.op { RxToken::ADD => .. }` read as upstream does.
// kRxMaxStackDepth is K_RX_MAX_STACK_DEPTH.

/// Hydraulic variables available in reaction expressions. Order is the ABI of the
/// evaluator's hydvar array.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum RxHydVar {
    /// depth (ft)
    D = 0,
    /// flow (cfs)
    Q = 1,
    /// velocity (ft/s)
    U = 2,
    /// Reynolds number
    Re = 3,
    /// shear velocity (ft/s)
    Us = 4,
    /// Darcy-Weisbach friction factor
    Ff = 5,
    /// wetted surface area per volume (1/ft)
    Av = 6,
    /// hydraulic residence time (s)
    Hrt = 7,
    /// reaction step (s)
    Dt = 8,
    /// local water temperature (degC)
    Temp = 9,
    /// water content (dimensionless), groundwater only
    GwTheta = 10,
    /// porosity (dimensionless), groundwater only
    GwPorosity = 11,
    /// effective saturation (dimensionless), groundwater only
    GwSe = 12,
    /// reactor centre depth below the ground surface (m)
    GwDepthM = 13,
    Count = 14,
}

/// RxToken::Op
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum RxTokenOp {
    /// value
    PushNum = 0,
    /// idx = local species slot
    PushSpecies = 1,
    /// idx = coefficient slot
    PushCoef = 2,
    /// idx = term slot (evaluated earlier in-order)
    PushTerm = 3,
    /// idx = RxHydVar
    PushHydvar = 4,
    /// idx = pollutant index (READ-ONLY coupling, R4)
    PushPollut = 5,
    Add = 6,
    Sub = 7,
    Mul = 8,
    Div = 9,
    Pow = 10,
    Neg = 11,
    FExp = 12,
    FLog = 13,
    FLog10 = 14,
    FSqrt = 15,
    FAbs = 16,
    FSgn = 17,
    FStep = 18,
    FSin = 19,
    FCos = 20,
    FTan = 21,
    FMin = 22,
    /// binary
    FMax = 23,
}

/// One VM instruction. POD, contiguous in the shared pool.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct RxToken {
    pub op: u8,
    pub idx: i32,
    pub value: f64,
}

// op = PUSH_NUM (0), idx = 0, value = 0.0 are Rust's defaults, so Default is derived.

impl RxToken {
    pub const PUSH_NUM: u8 = RxTokenOp::PushNum as u8;
    pub const PUSH_SPECIES: u8 = RxTokenOp::PushSpecies as u8;
    pub const PUSH_COEF: u8 = RxTokenOp::PushCoef as u8;
    pub const PUSH_TERM: u8 = RxTokenOp::PushTerm as u8;
    pub const PUSH_HYDVAR: u8 = RxTokenOp::PushHydvar as u8;
    pub const PUSH_POLLUT: u8 = RxTokenOp::PushPollut as u8;
    pub const ADD: u8 = RxTokenOp::Add as u8;
    pub const SUB: u8 = RxTokenOp::Sub as u8;
    pub const MUL: u8 = RxTokenOp::Mul as u8;
    pub const DIV: u8 = RxTokenOp::Div as u8;
    pub const POW: u8 = RxTokenOp::Pow as u8;
    pub const NEG: u8 = RxTokenOp::Neg as u8;
    pub const F_EXP: u8 = RxTokenOp::FExp as u8;
    pub const F_LOG: u8 = RxTokenOp::FLog as u8;
    pub const F_LOG10: u8 = RxTokenOp::FLog10 as u8;
    pub const F_SQRT: u8 = RxTokenOp::FSqrt as u8;
    pub const F_ABS: u8 = RxTokenOp::FAbs as u8;
    pub const F_SGN: u8 = RxTokenOp::FSgn as u8;
    pub const F_STEP: u8 = RxTokenOp::FStep as u8;
    pub const F_SIN: u8 = RxTokenOp::FSin as u8;
    pub const F_COS: u8 = RxTokenOp::FCos as u8;
    pub const F_TAN: u8 = RxTokenOp::FTan as u8;
    pub const F_MIN: u8 = RxTokenOp::FMin as u8;
    pub const F_MAX: u8 = RxTokenOp::FMax as u8;
}

/// Compiled expression = span into the shared pool. len 0 ⇒ no expression.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct RxExprSpan {
    pub begin: i32,
    pub len: i32,
}

/// Maximum operand-stack depth; the compiler rejects deeper expressions.
pub const K_RX_MAX_STACK_DEPTH: i32 = 32;
