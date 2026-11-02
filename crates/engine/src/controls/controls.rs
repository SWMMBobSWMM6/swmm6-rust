// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/controls/Controls.hpp and Controls.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The control-rule engine's data and its build half: the rule records, the premise SoA,
// ControlEngine::init, clearRules, parseRuleText (the full [CONTROLS] text parser),
// buildPremiseSoA, addNamedVariable and addExpression, with the parser's file-local
// helpers.
//
// NOT YET TRANSLATED (run-time functions, with the controls step): evaluate,
// batchEvaluateGroup, applyPendingActions, updateActionValue, evaluatePremise,
// compareValues, compareTimes, getVariableValue, computePIDSetting,
// getNamedVariableValue.
//
// Differences in form, never in result:
//   * parse_rule_text returns Err(StodError) where upstream lets std::stod's exception
//     escape: addExpression -> mathexpr::parse is called with no `try` (finding U1), and
//     SWMMEngine::initHydraulics is noexcept, so upstream terminates. The caller reports
//     it instead of aborting.
//   * `expr_index_` (unordered_map with CiHash/CiEqual) is only looked up: a HashMap keyed
//     by the ASCII-uppercased name, where `map[name] = idx` overwrites the value.
//   * buildPremiseSoA ITERATES `std::unordered_map<int,int> type_counts` to create the
//     premise groups, so the group order is libstdc++'s: reproduced with
//     core::std_unordered_map::UnorderedMapOrder and std::hash<int>.
//   * std::istringstream >> token splits at the C locale's isspace bytes; std::getline
//     at '\n' (a final line without one is still read; nothing follows a trailing '\n').
//   * std::stod(s, &pos) is core::cstdlib::strtod on s up to its first NUL plus the
//     ERANGE rules glibc applies (see stod_pos); std::strtol is strtol_base10 below.
//   * `std::vector<bool>` members are Vec<bool>.
//   * Members keep upstream's names and are `pub` (upstream's private ones too) so the
//     gate's state dumper can read them.

use std::collections::HashMap;

use crate::core::context::SimulationContext;
use crate::core::cstdlib::{c_str, strtod};
use crate::core::datetime;
use crate::core::error_codes::{ERR_CONTROL_RULE, ERR_ITEMS, ERR_KEYWORD, ERR_NAME, ERR_NUMBER};
use crate::core::options::Bytes;
use crate::core::std_unordered_map::{hash_int, UnorderedMapOrder};
use crate::data::gage::GageData;
use crate::data::link::LinkType;
use crate::data::table::{sscanf, ScanArg};
use crate::math::math_expr::{self, Expression, StodError};

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CompareOp {
    #[default]
    Eq = 0,
    Ne = 1,
    Lt = 2,
    Le = 3,
    Gt = 4,
    Ge = 5,
}

impl CompareOp {
    fn from_i32(v: i32) -> CompareOp {
        match v {
            1 => CompareOp::Ne,
            2 => CompareOp::Lt,
            3 => CompareOp::Le,
            4 => CompareOp::Gt,
            5 => CompareOp::Ge,
            _ => CompareOp::Eq,
        }
    }
}

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LogicOp {
    #[default]
    And = 0,
    Or = 1,
}

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ConditionVar {
    #[default]
    NodeDepth = 0,
    NodeMaxdepth = 1,
    NodeHead = 2,
    NodeVolume = 3,
    NodeInflow = 4,
    LinkFlow = 10,
    LinkDepth = 11,
    LinkSetting = 12,
    LinkStatus = 13,
    LinkFullflow = 14,
    LinkFulldepth = 15,
    LinkVelocity = 16,
    LinkLength = 17,
    LinkSlope = 18,
    LinkTimeopen = 19,
    LinkTimeclosed = 20,
    GageRain = 30,
    GageRainPast = 31,
    SimTime = 40,
    SimDate = 41,
    ClockTime = 42,
    SimDay = 43,
    SimMonth = 44,
    SimDayofyear = 45,
}

impl ConditionVar {
    fn from_i32(v: i32) -> ConditionVar {
        use ConditionVar::*;
        match v {
            1 => NodeMaxdepth,
            2 => NodeHead,
            3 => NodeVolume,
            4 => NodeInflow,
            10 => LinkFlow,
            11 => LinkDepth,
            12 => LinkSetting,
            13 => LinkStatus,
            14 => LinkFullflow,
            15 => LinkFulldepth,
            16 => LinkVelocity,
            17 => LinkLength,
            18 => LinkSlope,
            19 => LinkTimeopen,
            20 => LinkTimeclosed,
            30 => GageRain,
            31 => GageRainPast,
            40 => SimTime,
            41 => SimDate,
            42 => ClockTime,
            43 => SimDay,
            44 => SimMonth,
            45 => SimDayofyear,
            _ => NodeDepth,
        }
    }
}

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ActionType {
    #[default]
    Numeric = 0,
    Curve = 1,
    Timeseries = 2,
    Pid = 3,
}

#[derive(Debug, Clone, Copy)]
pub struct Premise {
    pub logic: LogicOp,
    pub is_expression: bool,
    pub expr_idx: i32,
    pub lhs_var: ConditionVar,
    pub lhs_idx: i32,
    pub lhs_param: i32,
    pub op: CompareOp,
    pub rhs_is_variable: bool,
    pub rhs_var: ConditionVar,
    pub rhs_idx: i32,
    pub rhs_param: i32,
    pub rhs_value: f64,
}

impl Default for Premise {
    fn default() -> Self {
        Premise {
            logic: LogicOp::And,
            is_expression: false,
            expr_idx: -1,
            lhs_var: ConditionVar::NodeDepth,
            lhs_idx: -1,
            lhs_param: 0,
            op: CompareOp::Eq,
            rhs_is_variable: false,
            rhs_var: ConditionVar::NodeDepth,
            rhs_idx: -1,
            rhs_param: 0,
            rhs_value: 0.0,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Action {
    pub link_idx: i32,
    pub rule_idx: i32,
    pub type_: ActionType,
    pub value: f64,
    pub curve_idx: i32,
    pub tseries_idx: i32,
    pub pid_idx: i32,
}

impl Default for Action {
    fn default() -> Self {
        Action { link_idx: -1, rule_idx: -1, type_: ActionType::Numeric, value: 0.0, curve_idx: -1, tseries_idx: -1, pid_idx: -1 }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct PIDState {
    pub kp: f64,
    pub ki: f64,
    pub kd: f64,
    pub setpoint: f64,
    pub e1: f64,
    pub e2: f64,
    pub action_idx: i32,
}

impl Default for PIDState {
    fn default() -> Self {
        PIDState { kp: 0.0, ki: 0.0, kd: 0.0, setpoint: 0.0, e1: 0.0, e2: 0.0, action_idx: -1 }
    }
}

#[derive(Debug, Clone)]
pub struct NamedVariable {
    pub name: Bytes,
    pub var: ConditionVar,
    pub idx: i32,
    pub param: i32,
}

#[derive(Debug, Clone, Default)]
pub struct Rule {
    pub name: Bytes,
    pub premises: Vec<Premise>,
    pub then_actions: Vec<Action>,
    pub else_actions: Vec<Action>,
    pub priority: f64,
}

#[derive(Debug, Clone)]
pub struct ParseError {
    pub line: i32,
    pub message: Bytes,
    pub code: i32,
    pub token: Bytes,
}

impl Default for ParseError {
    fn default() -> Self {
        ParseError { line: -1, message: Vec::new(), code: 0, token: Vec::new() }
    }
}

/// All premises testing one variable type (batch evaluation).
#[derive(Debug, Clone, Default)]
pub struct PremiseSoA {
    pub var_type: ConditionVar,
    pub count: i32,
    pub rule_idx: Vec<i32>,
    pub premise_idx: Vec<i32>,
    pub obj_idx: Vec<i32>,
    pub param: Vec<i32>,
    pub flat_idx: Vec<i32>,
    pub op: Vec<i32>,
    pub rhs_value: Vec<f64>,
    pub rhs_is_variable: Vec<bool>,
    pub is_expression: Vec<bool>,
    pub lhs_values: Vec<f64>,
    pub results: Vec<bool>,
}

#[derive(Debug, Clone, Copy)]
pub struct PendingAction {
    pub link_idx: i32,
    pub value: f64,
    pub priority: f64,
    pub rule_idx: i32,
    pub type_: ActionType,
}

#[derive(Debug, Clone, Default)]
pub struct ControlEngine {
    pub last_action_count_: i32,
    pub rules_: Vec<Rule>,
    pub pid_states_: Vec<PIDState>,
    pub named_vars_: Vec<NamedVariable>,
    pub expressions_: Vec<Expression>,
    /// case-insensitive name -> index (see the module notes)
    pub expr_index_: HashMap<Bytes, i32>,
    pub last_parse_error_: ParseError,
    pub premise_groups_: Vec<PremiseSoA>,
    pub premise_results_: Vec<u8>,
    pub rule_premise_offset_: Vec<i32>,
    pub total_premises_: i32,
    pub rule_results_: Vec<bool>,
    pub link_time_last_set_: Vec<f64>,
    pub control_value_: f64,
    pub set_point_: f64,
    pub eval_date_: f64,
    pub cur_date_: f64,
    pub cur_time_: f64,
    pub elapsed_days_: f64,
    pub pending_actions_: Vec<PendingAction>,
    pub action_slots_: Vec<i32>,
}

// ============================================================================
// Parser helpers (Controls.cpp's anonymous namespace)
// ============================================================================

fn to_upper(s: &[u8]) -> Bytes {
    s.iter().map(|c| c.to_ascii_uppercase()).collect()
}

fn ci_key(s: &[u8]) -> Bytes {
    to_upper(s)
}

fn ieq(a: &[u8], b: &[u8]) -> bool {
    a.eq_ignore_ascii_case(b)
}

fn c_isspace(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// istringstream >> token, repeatedly.
fn tokenize(line: &[u8]) -> Vec<Bytes> {
    let mut toks = Vec::new();
    let mut i = 0;
    while i < line.len() {
        while i < line.len() && c_isspace(line[i]) {
            i += 1;
        }
        if i >= line.len() {
            break;
        }
        let start = i;
        while i < line.len() && !c_isspace(line[i]) {
            i += 1;
        }
        toks.push(line[start..i].to_vec());
    }
    toks
}

/// std::stod(s, &pos): Ok((value, pos)) or the exception it throws. ERANGE follows glibc:
/// overflow (an infinite result from a non-"inf" literal) and underflow (a zero result
/// from a mantissa holding a nonzero digit, or any subnormal result; an exact subnormal
/// literal is not distinguished).
fn stod_pos(s: &[u8]) -> Result<(f64, usize), StodError> {
    let (v, used) = strtod(s);
    if used == 0 {
        return Err(StodError::InvalidArgument);
    }
    let t = &c_str(s)[..used];
    let mut k = 0;
    while k < t.len() && c_isspace(t[k]) {
        k += 1;
    }
    if k < t.len() && (t[k] == b'+' || t[k] == b'-') {
        k += 1;
    }
    let body = &t[k..];
    if v.is_nan() {
        return Ok((v, used));
    }
    if v.is_infinite() {
        let lit = body.first().is_some_and(|c| c.to_ascii_lowercase() == b'i');
        return if lit { Ok((v, used)) } else { Err(StodError::OutOfRange) };
    }
    let is_hex = body.len() >= 2 && body[0] == b'0' && body[1].to_ascii_lowercase() == b'x';
    if v == 0.0 {
        let nonzero = if is_hex {
            body[2..].iter().take_while(|&&c| c.to_ascii_lowercase() != b'p').any(|&c| c.is_ascii_hexdigit() && c != b'0')
        } else {
            body.iter().take_while(|&&c| c != b'e' && c != b'E').any(|&c| (b'1'..=b'9').contains(&c))
        };
        if nonzero {
            return Err(StodError::OutOfRange);
        }
    } else if v.abs() < f64::MIN_POSITIVE {
        return Err(StodError::OutOfRange);
    }
    Ok((v, used))
}

/// tryParseDouble: a complete std::stod conversion. `val` changes whenever stod returns,
/// even when trailing bytes make the result false.
fn try_parse_double(s: &[u8], val: &mut f64) -> bool {
    match stod_pos(s) {
        Ok((v, pos)) => {
            *val = v;
            pos == s.len()
        }
        Err(_) => false,
    }
}

/// std::strtol(s, nullptr, 10): leading white space, a sign, decimal digits; saturates
/// at LONG_MIN/LONG_MAX. The C string ends at its first NUL.
fn strtol_base10(s: &[u8]) -> i64 {
    let s = c_str(s);
    let mut i = 0;
    while i < s.len() && c_isspace(s[i]) {
        i += 1;
    }
    let mut neg = false;
    if i < s.len() && (s[i] == b'+' || s[i] == b'-') {
        neg = s[i] == b'-';
        i += 1;
    }
    let mut v: i128 = 0;
    while i < s.len() && s[i].is_ascii_digit() {
        if v <= i64::MAX as i128 + 1 {
            v = v * 10 + (s[i] - b'0') as i128;
        }
        i += 1;
    }
    let v = if neg { -v } else { v };
    if v > i64::MAX as i128 {
        i64::MAX
    } else if v < i64::MIN as i128 {
        i64::MIN
    } else {
        v as i64
    }
}

fn match_rel_op(tok: &[u8]) -> i32 {
    match tok {
        b"=" => CompareOp::Eq as i32,
        b"<>" => CompareOp::Ne as i32,
        b"<" => CompareOp::Lt as i32,
        b"<=" => CompareOp::Le as i32,
        b">" => CompareOp::Gt as i32,
        b">=" => CompareOp::Ge as i32,
        _ => -1,
    }
}

fn node_attribute(attr: &[u8]) -> i32 {
    match attr {
        b"DEPTH" => ConditionVar::NodeDepth as i32,
        b"MAXDEPTH" => ConditionVar::NodeMaxdepth as i32,
        b"HEAD" => ConditionVar::NodeHead as i32,
        b"VOLUME" => ConditionVar::NodeVolume as i32,
        b"INFLOW" => ConditionVar::NodeInflow as i32,
        _ => -1,
    }
}

fn link_attribute(attr: &[u8]) -> i32 {
    match attr {
        b"FLOW" => ConditionVar::LinkFlow as i32,
        b"DEPTH" => ConditionVar::LinkDepth as i32,
        b"SETTING" => ConditionVar::LinkSetting as i32,
        b"STATUS" => ConditionVar::LinkStatus as i32,
        b"FULLFLOW" => ConditionVar::LinkFullflow as i32,
        b"FULLDEPTH" => ConditionVar::LinkFulldepth as i32,
        b"VELOCITY" => ConditionVar::LinkVelocity as i32,
        b"LENGTH" => ConditionVar::LinkLength as i32,
        b"SLOPE" => ConditionVar::LinkSlope as i32,
        b"TIMEOPEN" => ConditionVar::LinkTimeopen as i32,
        b"TIMECLOSED" => ConditionVar::LinkTimeclosed as i32,
        _ => -1,
    }
}

fn sim_attribute(attr: &[u8]) -> i32 {
    match attr {
        b"TIME" => ConditionVar::SimTime as i32,
        b"DATE" => ConditionVar::SimDate as i32,
        b"CLOCKTIME" => ConditionVar::ClockTime as i32,
        b"DAY" => ConditionVar::SimDay as i32,
        b"MONTH" => ConditionVar::SimMonth as i32,
        b"DAYOFYEAR" => ConditionVar::SimDayofyear as i32,
        _ => -1,
    }
}

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

/// The outputs of parsePremiseVariable (its reference parameters).
struct VarOut {
    cv: ConditionVar,
    obj_idx: i32,
    extra_param: i32,
    err: Bytes,
    err_code: i32,
    err_tok: Bytes,
}

/// parsePremiseVariable. `out` carries the reference parameters' prior values in and
/// their (possibly partly changed) values out, as upstream's references do.
fn parse_premise_variable(toks: &[Bytes], k: &mut i32, ctx: &SimulationContext, out: &mut VarOut) -> bool {
    out.extra_param = 0;
    out.err_code = ERR_ITEMS;
    out.err_tok.clear();
    let n = toks.len() as i32;
    if *k >= n {
        out.err = b"condition clause is empty".to_vec();
        return false;
    }
    let obj_type = to_upper(&toks[*k as usize]);
    let incomplete = |toks: &[Bytes]| -> Bytes { cat(&[b"condition clause is incomplete after '", toks.last().unwrap(), b"'"]) };

    if obj_type == b"NODE" {
        if *k + 2 >= n {
            out.err = incomplete(toks);
            return false;
        }
        let t1 = &toks[(*k + 1) as usize];
        out.obj_idx = ctx.node_names.find(t1);
        if out.obj_idx < 0 {
            out.err_code = ERR_NAME;
            out.err_tok = t1.clone();
            out.err = cat(&[b"'", t1, b"' is not the name of a defined NODE"]);
            return false;
        }
        let t2 = &toks[(*k + 2) as usize];
        let a = node_attribute(&to_upper(t2));
        if a < 0 {
            out.err_code = ERR_KEYWORD;
            out.err_tok = t2.clone();
            out.err = cat(&[b"'", t2, b"' is not a valid NODE attribute (expected DEPTH, MAXDEPTH, HEAD, VOLUME or INFLOW)"]);
            return false;
        }
        out.cv = ConditionVar::from_i32(a);
        *k += 3;
        return true;
    }
    if obj_type == b"LINK" || obj_type == b"CONDUIT" || obj_type == b"PUMP" || obj_type == b"ORIFICE" || obj_type == b"WEIR" || obj_type == b"OUTLET" {
        if *k + 2 >= n {
            out.err = incomplete(toks);
            return false;
        }
        let t1 = &toks[(*k + 1) as usize];
        out.obj_idx = ctx.link_names.find(t1);
        if out.obj_idx < 0 {
            out.err_code = ERR_NAME;
            out.err_tok = t1.clone();
            out.err = cat(&[b"'", t1, b"' is not the name of a defined ", &obj_type]);
            return false;
        }
        let t2 = &toks[(*k + 2) as usize];
        let a = link_attribute(&to_upper(t2));
        if a < 0 {
            out.err_code = ERR_KEYWORD;
            out.err_tok = t2.clone();
            out.err = cat(&[b"'", t2, b"' is not a valid ", &obj_type, b" attribute"]);
            return false;
        }
        out.cv = ConditionVar::from_i32(a);
        *k += 3;
        return true;
    }
    if obj_type == b"GAGE" {
        if *k + 2 >= n {
            out.err = incomplete(toks);
            return false;
        }
        let t1 = &toks[(*k + 1) as usize];
        out.obj_idx = ctx.gage_names.find(t1);
        if out.obj_idx < 0 {
            out.err_code = ERR_NAME;
            out.err_tok = t1.clone();
            out.err = cat(&[b"'", t1, b"' is not the name of a defined GAGE"]);
            return false;
        }
        let t2 = &toks[(*k + 2) as usize];
        let attr = to_upper(t2);
        if attr == b"INTENSITY" {
            out.cv = ConditionVar::GageRain;
        } else {
            let hours = strtol_base10(&attr);
            if hours < 1 || hours > GageData::MAXPASTRAIN as i64 {
                out.err_code = ERR_KEYWORD;
                out.err_tok = t2.clone();
                out.err = cat(&[b"'", t2, b"' is not a valid GAGE attribute (expected INTENSITY or n-HR_DEPTH with n from 1 to 48)"]);
                return false;
            }
            out.cv = ConditionVar::GageRainPast;
            out.extra_param = hours as i32;
        }
        *k += 3;
        return true;
    }
    if obj_type == b"SIMULATION" {
        if *k + 1 >= n {
            out.err = incomplete(toks);
            return false;
        }
        let t1 = &toks[(*k + 1) as usize];
        let a = sim_attribute(&to_upper(t1));
        if a < 0 {
            out.err_code = ERR_KEYWORD;
            out.err_tok = t1.clone();
            out.err = cat(&[b"'", t1, b"' is not a valid SIMULATION attribute (expected TIME, DATE, CLOCKTIME, DAY, MONTH or DAYOFYEAR)"]);
            return false;
        }
        out.cv = ConditionVar::from_i32(a);
        out.obj_idx = -1;
        *k += 2;
        return true;
    }
    let t0 = &toks[*k as usize];
    out.err_code = ERR_KEYWORD;
    out.err_tok = t0.clone();
    out.err = cat(&[b"'", t0, b"' is not a known object, variable or expression in this condition"]);
    false
}

fn parse_status_value(tok: &[u8], val: &mut f64) -> bool {
    let u = to_upper(tok);
    if u == b"OFF" || u == b"CLOSED" {
        *val = 0.0;
        return true;
    }
    if u == b"ON" || u == b"OPEN" {
        *val = 1.0;
        return true;
    }
    false
}

fn parse_time_token(s: &[u8], out_days: &mut f64) -> bool {
    // Try decimal hours first
    if let Ok((hours, pos)) = stod_pos(s) {
        if pos == s.len() {
            *out_days = hours / 24.0;
            return true;
        }
    }
    // Parse HH:MM:SS (MM and SS optional, legacy sscanf behaviour)
    let mut hr = 0i32;
    let mut mn = 0i32;
    let mut sc = 0i32;
    let n = sscanf(s, b"%d:%d:%d", &mut [ScanArg::Int(&mut hr), ScanArg::Int(&mut mn), ScanArg::Int(&mut sc)]);
    if n >= 1 && hr >= 0 && mn >= 0 && sc >= 0 {
        *out_days = datetime::encode_time(hr, mn, sc);
        return true;
    }
    false
}

fn parse_date_token(s: &[u8], out_date: &mut f64) -> bool {
    let mut m = 0i32;
    let mut d = 0i32;
    let mut y = 0i32;
    let mut sep1 = 0u8;
    let mut sep2 = 0u8;
    if sscanf(s, b"%d%c%d%c%d", &mut [ScanArg::Int(&mut m), ScanArg::Char(&mut sep1), ScanArg::Int(&mut d), ScanArg::Char(&mut sep2), ScanArg::Int(&mut y)]) == 5 {
        *out_date = datetime::encode_date(y, m, d);
        return *out_date != -(datetime::DATE_DELTA as f64);
    }
    false
}

fn parse_day_of_year_token(s: &[u8], out_doy: &mut f64) -> bool {
    let mut m = 0i32;
    let mut d = 0i32;
    let mut sep = 0u8;
    if sscanf(s, b"%d%c%d", &mut [ScanArg::Int(&mut m), ScanArg::Char(&mut sep), ScanArg::Int(&mut d)]) == 3 {
        let date = datetime::encode_date(1947, m, d);
        if date != -(datetime::DATE_DELTA as f64) {
            *out_doy = datetime::day_of_year(date) as f64;
            return true;
        }
        return false;
    }
    let mut v = 0.0;
    if !try_parse_double(s, &mut v) {
        return false;
    }
    if !(1.0..=365.0).contains(&v) {
        return false;
    }
    *out_doy = v;
    true
}

fn parse_premise_rhs(tok: &[u8], lhs_var: ConditionVar, val: &mut f64) -> bool {
    match lhs_var {
        ConditionVar::LinkStatus => parse_status_value(tok, val),
        ConditionVar::SimTime | ConditionVar::ClockTime | ConditionVar::LinkTimeopen | ConditionVar::LinkTimeclosed => parse_time_token(tok, val),
        ConditionVar::SimDate => parse_date_token(tok, val),
        ConditionVar::SimDayofyear => parse_day_of_year_token(tok, val),
        ConditionVar::SimDay => {
            if !try_parse_double(tok, val) {
                return false;
            }
            *val >= 1.0 && *val <= 7.0
        }
        ConditionVar::SimMonth => {
            if !try_parse_double(tok, val) {
                return false;
            }
            *val >= 1.0 && *val <= 12.0
        }
        _ => try_parse_double(tok, val),
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum ParseState {
    Idle,
    InPremises,
    InThen,
    InElse,
}

impl ControlEngine {
    /// init(rules): install pre-built rules.
    pub fn init(&mut self, rules: &[Rule]) {
        self.rules_ = rules.to_vec();
        self.pid_states_.clear();
        for r in 0..self.rules_.len() {
            for a in self.rules_[r].then_actions.iter_mut() {
                a.rule_idx = r as i32;
            }
            for a in self.rules_[r].else_actions.iter_mut() {
                a.rule_idx = r as i32;
            }
        }
        self.build_premise_soa();
        self.link_time_last_set_.clear();
    }

    pub fn last_parse_error(&self) -> &ParseError {
        &self.last_parse_error_
    }

    pub fn rules(&self) -> &[Rule] {
        &self.rules_
    }

    pub fn last_action_count(&self) -> i32 {
        self.last_action_count_
    }

    /// Group the premises by LHS variable type (group order: libstdc++'s unordered_map).
    pub fn build_premise_soa(&mut self) {
        self.premise_groups_.clear();

        self.rule_premise_offset_.clear();
        self.rule_premise_offset_.resize(self.rules_.len() + 1, 0);
        for r in 0..self.rules_.len() {
            self.rule_premise_offset_[r + 1] = self.rule_premise_offset_[r] + self.rules_[r].premises.len() as i32;
        }
        self.total_premises_ = *self.rule_premise_offset_.last().unwrap();
        self.premise_results_.clear();
        self.premise_results_.resize(self.total_premises_ as usize, 0);

        // Count premises per variable type (unordered_map<int,int>).
        let mut tc_keys: Vec<i32> = Vec::new();
        let mut tc_vals: Vec<i32> = Vec::new();
        let mut tc_order = UnorderedMapOrder::new();
        for r in 0..self.rules_.len() {
            for prem in &self.rules_[r].premises {
                if !prem.is_expression {
                    let key = prem.lhs_var as i32;
                    match tc_keys.iter().position(|&k| k == key) {
                        Some(p) => tc_vals[p] += 1,
                        None => {
                            tc_order.insert(hash_int(key));
                            tc_keys.push(key);
                            tc_vals.push(1);
                        }
                    }
                }
            }
        }

        // Create groups, in the map's iteration order.
        let mut type_to_group: HashMap<i32, usize> = HashMap::new();
        for p in tc_order.order() {
            let (key, cnt) = (tc_keys[p], tc_vals[p]);
            type_to_group.insert(key, self.premise_groups_.len());
            let cn = cnt as usize;
            let g = PremiseSoA {
                var_type: ConditionVar::from_i32(key),
                count: 0,
                rule_idx: Vec::with_capacity(cn),
                premise_idx: Vec::with_capacity(cn),
                obj_idx: Vec::with_capacity(cn),
                param: Vec::with_capacity(cn),
                flat_idx: Vec::with_capacity(cn),
                op: Vec::with_capacity(cn),
                rhs_value: Vec::with_capacity(cn),
                rhs_is_variable: Vec::with_capacity(cn),
                is_expression: Vec::with_capacity(cn),
                lhs_values: Vec::new(),
                results: Vec::new(),
            };
            self.premise_groups_.push(g);
        }

        // Fill groups
        for r in 0..self.rules_.len() {
            let rule = &self.rules_[r];
            for (p, prem) in rule.premises.iter().enumerate() {
                if prem.is_expression {
                    continue;
                }
                let gi = match type_to_group.get(&(prem.lhs_var as i32)) {
                    None => continue,
                    Some(&g) => g,
                };
                let g = &mut self.premise_groups_[gi];
                g.rule_idx.push(r as i32);
                g.premise_idx.push(p as i32);
                g.obj_idx.push(prem.lhs_idx);
                g.param.push(prem.lhs_param);
                g.flat_idx.push(self.rule_premise_offset_[r] + p as i32);
                g.op.push(prem.op as i32);
                g.rhs_value.push(prem.rhs_value);
                g.rhs_is_variable.push(prem.rhs_is_variable);
                g.is_expression.push(false);
                g.count += 1;
            }
        }

        for g in self.premise_groups_.iter_mut() {
            let un = g.count as usize;
            g.lhs_values.resize(un, 0.0);
            g.results.resize(un, false);
        }

        self.rule_results_.resize(self.rules_.len(), true);
    }

    pub fn add_named_variable(&mut self, name: &[u8], var: ConditionVar, idx: i32, param: i32) {
        self.named_vars_.push(NamedVariable { name: name.to_vec(), var, idx, param });
    }

    /// addExpression: Ok(index), Ok(-1) when the formula does not parse, Err where
    /// std::stod throws (see the module notes).
    pub fn add_expression(&mut self, name: &[u8], formula: &[u8]) -> Result<i32, StodError> {
        let mut expr = Expression::default();
        if math_expr::parse_stod(formula, &mut expr)? != 0 {
            return Ok(-1);
        }
        let idx = self.expressions_.len() as i32;
        self.expressions_.push(expr);
        self.expr_index_.insert(ci_key(name), idx);
        Ok(idx)
    }

    pub fn clear_rules(&mut self) {
        self.rules_.clear();
        self.pid_states_.clear();
        self.named_vars_.clear();
        self.expressions_.clear();
        self.expr_index_.clear();
        self.premise_groups_.clear();
        self.premise_results_.clear();
        self.rule_premise_offset_.clear();
        self.total_premises_ = 0;
        self.rule_results_.clear();
        self.pending_actions_.clear();
        self.action_slots_.clear();
        self.last_parse_error_ = ParseError::default();
    }

    /// parseRuleText: the number of rules parsed, -1 on a rejection (last_parse_error_
    /// says why), or Err where upstream's std::stod exception escapes.
    #[allow(unused_assignments)]
    pub fn parse_rule_text(&mut self, text: &[u8], ctx: &SimulationContext) -> Result<i32, StodError> {
        self.last_parse_error_ = ParseError::default();

        // Split into lines (std::getline at '\n'), trimming " \t\r\n", keeping 1-based
        // line numbers of the original text.
        let mut lines: Vec<(i32, Bytes)> = Vec::new();
        {
            let mut lineno = 0;
            let mut pos = 0usize;
            while pos < text.len() {
                let end = text[pos..].iter().position(|&c| c == b'\n').map(|e| pos + e);
                let line = match end {
                    Some(e) => &text[pos..e],
                    None => &text[pos..],
                };
                pos = match end {
                    Some(e) => e + 1,
                    None => text.len(),
                };
                lineno += 1;
                let ws = |c: &u8| matches!(*c, b' ' | b'\t' | b'\r' | b'\n');
                let start = match line.iter().position(|c| !ws(c)) {
                    None => continue,
                    Some(s) => s,
                };
                let last = line.iter().rposition(|c| !ws(c)).unwrap();
                let trimmed = line[start..=last].to_vec();
                if !trimmed.is_empty() {
                    lines.push((lineno, trimmed));
                }
            }
        }

        if lines.is_empty() {
            return Ok(0);
        }

        let mut current_line: i32 = -1;
        macro_rules! fail {
            ($why:expr) => {
                fail!($why, ERR_CONTROL_RULE, Vec::<u8>::new())
            };
            ($why:expr, $code:expr) => {
                fail!($why, $code, Vec::<u8>::new())
            };
            ($why:expr, $code:expr, $tok:expr) => {{
                self.last_parse_error_.line = current_line;
                self.last_parse_error_.message = $why;
                self.last_parse_error_.code = $code;
                self.last_parse_error_.token = $tok;
                return Ok(-1);
            }};
        }

        let mut state = ParseState::Idle;
        let mut current_rule = Rule::default();
        let mut rules_parsed = 0i32;

        for (number, text_line) in &lines {
            current_line = *number;
            let toks = tokenize(text_line);
            if toks.is_empty() {
                continue;
            }
            let nt = toks.len() as i32;

            let keyword = to_upper(&toks[0]);

            // ---- VARIABLE keyword ----
            if keyword == b"VARIABLE" {
                if toks.len() < 6 || toks[2] != b"=" {
                    fail!(b"expected: VARIABLE <name> = <object> <id> <attribute>".to_vec(), ERR_ITEMS);
                }
                let mut k = 3;
                let mut vo = VarOut { cv: ConditionVar::NodeDepth, obj_idx: -1, extra_param: 0, err: Vec::new(), err_code: ERR_ITEMS, err_tok: Vec::new() };
                if !parse_premise_variable(&toks, &mut k, ctx, &mut vo) {
                    fail!(cat(&[b"in VARIABLE '", &toks[1], b"': ", &vo.err]), vo.err_code, vo.err_tok.clone());
                }
                let (name, cv, idx, extra) = (toks[1].clone(), vo.cv, vo.obj_idx, vo.extra_param);
                self.add_named_variable(&name, cv, idx, extra);
                continue;
            }

            // ---- EXPRESSION keyword ----
            if keyword == b"EXPRESSION" {
                if toks.len() < 4 || toks[2] != b"=" {
                    fail!(b"expected: EXPRESSION <name> = <formula>".to_vec(), ERR_ITEMS);
                }
                let mut formula: Bytes = Vec::new();
                for (i, t) in toks.iter().enumerate().skip(3) {
                    if i > 3 {
                        formula.push(b' ');
                    }
                    formula.extend_from_slice(t);
                }
                let name = toks[1].clone();
                if self.add_expression(&name, &formula)? < 0 {
                    fail!(cat(&[b"could not parse the formula for EXPRESSION '", &toks[1], b"'"]));
                }
                continue;
            }

            // ---- RULE keyword ----
            if keyword == b"RULE" {
                if state != ParseState::Idle {
                    // finishRule()
                    if !current_rule.name.is_empty() {
                        self.rules_.push(std::mem::take(&mut current_rule));
                        rules_parsed += 1;
                    }
                    current_rule = Rule::default();
                    state = ParseState::Idle;
                }
                if toks.len() < 2 {
                    fail!(b"RULE is missing its name".to_vec(), ERR_ITEMS);
                }
                current_rule.name = toks[1].clone();
                state = ParseState::Idle; // wait for IF
                continue;
            }

            // ---- PRIORITY keyword ----
            if keyword == b"PRIORITY" {
                if toks.len() < 2 {
                    fail!(b"PRIORITY is missing its value".to_vec(), ERR_ITEMS);
                }
                let mut p = 0.0;
                if !try_parse_double(&toks[1], &mut p) {
                    fail!(cat(&[b"PRIORITY value '", &toks[1], b"' is not a number"]), ERR_NUMBER, toks[1].clone());
                }
                current_rule.priority = p;
                continue;
            }

            // ---- IF / AND / OR keywords (premises) ----
            if keyword == b"IF" || keyword == b"AND" || keyword == b"OR" {
                let mut is_premise_context = keyword == b"IF" || state == ParseState::InPremises;

                if keyword == b"IF" {
                    state = ParseState::InPremises;
                }

                if (keyword == b"AND" || keyword == b"OR") && (state == ParseState::InThen || state == ParseState::InElse) {
                    is_premise_context = false;
                }

                if is_premise_context {
                    let mut prem = Premise { logic: if keyword == b"OR" { LogicOp::Or } else { LogicOp::And }, ..Premise::default() };

                    let mut k: i32 = 1;

                    let mut lhs = VarOut { cv: ConditionVar::NodeDepth, obj_idx: -1, extra_param: 0, err: Vec::new(), err_code: ERR_ITEMS, err_tok: Vec::new() };
                    if k < nt {
                        let tk = &toks[k as usize];
                        let ei = self.expr_index_.get(&ci_key(tk)).copied().unwrap_or(-1);
                        if ei >= 0 {
                            prem.is_expression = true;
                            prem.expr_idx = ei;
                            k += 1;
                        } else if let Some(nv) = self.named_vars_.iter().find(|nv| ieq(&nv.name, tk)) {
                            lhs.cv = nv.var;
                            lhs.obj_idx = nv.idx;
                            lhs.extra_param = nv.param;
                            prem.lhs_var = lhs.cv;
                            prem.lhs_idx = lhs.obj_idx;
                            prem.lhs_param = lhs.extra_param;
                            k += 1;
                        } else if !parse_premise_variable(&toks, &mut k, ctx, &mut lhs) {
                            fail!(lhs.err.clone(), lhs.err_code, lhs.err_tok.clone());
                        } else {
                            prem.lhs_var = lhs.cv;
                            prem.lhs_idx = lhs.obj_idx;
                            prem.lhs_param = lhs.extra_param;
                        }
                    } else {
                        fail!(cat(&[b"condition clause is empty after ", &keyword]), ERR_ITEMS);
                    }
                    let lhs_cv = lhs.cv;

                    // Relational operator
                    if k >= nt {
                        fail!(b"condition is missing a relational operator".to_vec(), ERR_ITEMS);
                    }
                    let rel = match_rel_op(&toks[k as usize]);
                    if rel < 0 {
                        fail!(cat(&[b"'", &toks[k as usize], b"' is not a relational operator (expected =, <>, <, <=, > or >=)"]), ERR_KEYWORD, toks[k as usize].clone());
                    }
                    prem.op = CompareOp::from_i32(rel);
                    k += 1;

                    // RHS: a variable or a constant value
                    if k >= nt {
                        fail!(b"condition is missing its right-hand value".to_vec(), ERR_ITEMS);
                    }

                    let saved_k = k;
                    let mut rhs = VarOut { cv: ConditionVar::NodeDepth, obj_idx: -1, extra_param: 0, err: Vec::new(), err_code: 0, err_tok: Vec::new() };
                    if let Some(nv) = self.named_vars_.iter().find(|nv| ieq(&nv.name, &toks[k as usize])) {
                        prem.rhs_is_variable = true;
                        prem.rhs_var = nv.var;
                        prem.rhs_idx = nv.idx;
                        prem.rhs_param = nv.param;
                        k += 1;
                    } else if parse_premise_variable(&toks, &mut k, ctx, &mut rhs) {
                        prem.rhs_is_variable = true;
                        prem.rhs_var = rhs.cv;
                        prem.rhs_idx = rhs.obj_idx;
                        prem.rhs_param = rhs.extra_param;
                    } else {
                        k = saved_k;
                        let val_tok = toks[k as usize].clone();
                        let mut val = 0.0;
                        if !parse_premise_rhs(&val_tok, lhs_cv, &mut val) {
                            fail!(cat(&[b"'", &val_tok, b"' is not a valid value for this condition variable"]), ERR_NUMBER, val_tok);
                        }
                        prem.rhs_value = val;
                        k += 1;
                    }
                    let _ = k;

                    current_rule.premises.push(prem);
                    continue;
                }
                // Fall through: AND after THEN/ELSE is handled as an action below
            }

            if keyword == b"THEN" {
                state = ParseState::InThen;
            }
            if keyword == b"ELSE" {
                state = ParseState::InElse;
            }

            // ---- Parse action (THEN/ELSE/AND after THEN or ELSE) ----
            if state == ParseState::InThen || state == ParseState::InElse {
                let mut k: i32 = 1;

                if k >= nt {
                    continue;
                }

                let obj_type = to_upper(&toks[k as usize]);
                if obj_type != b"LINK" && obj_type != b"CONDUIT" && obj_type != b"PUMP" && obj_type != b"ORIFICE" && obj_type != b"WEIR" && obj_type != b"OUTLET" {
                    fail!(cat(&[b"'", &toks[k as usize], b"' cannot start an action (expected LINK, CONDUIT, PUMP, ORIFICE, WEIR or OUTLET)"]), ERR_KEYWORD, toks[k as usize].clone());
                }
                k += 1;

                // Link name
                if k >= nt {
                    fail!(cat(&[&obj_type, b" action is missing its link name"]), ERR_ITEMS);
                }
                let link_idx = ctx.link_names.find(&toks[k as usize]);
                if link_idx < 0 {
                    fail!(cat(&[b"no link named '", &toks[k as usize], b"' exists in the model"]), ERR_NAME, toks[k as usize].clone());
                }
                let link_type = ctx.links.type_[link_idx as usize];
                if (obj_type == b"CONDUIT" && link_type != LinkType::Conduit)
                    || (obj_type == b"PUMP" && link_type != LinkType::Pump)
                    || (obj_type == b"ORIFICE" && link_type != LinkType::Orifice)
                    || (obj_type == b"WEIR" && link_type != LinkType::Weir)
                    || (obj_type == b"OUTLET" && link_type != LinkType::Outlet)
                {
                    fail!(cat(&[b"link '", &toks[k as usize], b"' is not a ", &obj_type]), ERR_NAME, toks[k as usize].clone());
                }
                k += 1;

                // Attribute
                if k >= nt {
                    fail!(b"action is missing its attribute (STATUS or SETTING)".to_vec(), ERR_ITEMS);
                }
                let attr = to_upper(&toks[k as usize]);
                if attr != b"STATUS" && attr != b"SETTING" {
                    fail!(cat(&[b"'", &toks[k as usize], b"' is not a valid action attribute (expected STATUS or SETTING)"]), ERR_KEYWORD, toks[k as usize].clone());
                }
                if obj_type == b"LINK" {
                    fail!(b"a LINK action is not accepted; name the link's type (CONDUIT, PUMP, ORIFICE, WEIR or OUTLET)".to_vec(), ERR_KEYWORD, obj_type.clone());
                }
                if obj_type == b"CONDUIT" && attr != b"STATUS" {
                    fail!(b"CONDUIT actions take STATUS only".to_vec(), ERR_KEYWORD, attr.clone());
                }
                if (obj_type == b"ORIFICE" || obj_type == b"WEIR" || obj_type == b"OUTLET") && attr != b"SETTING" {
                    fail!(cat(&[&obj_type, b" actions take SETTING only"]), ERR_KEYWORD, attr.clone());
                }
                k += 1;

                // '='
                if k >= nt {
                    fail!(cat(&[b"action is missing '=' after ", &attr]), ERR_ITEMS);
                }
                if toks[k as usize] != b"=" {
                    fail!(cat(&[b"expected '=' after ", &attr, b", found '", &toks[k as usize], b"'"]));
                }
                k += 1;

                // Value / control type
                if k >= nt {
                    fail!(b"action is missing its value after '='".to_vec(), ERR_ITEMS);
                }

                let mut action = Action { link_idx, ..Action::default() };

                let val_tok = to_upper(&toks[k as usize]);

                if val_tok == b"CURVE" {
                    k += 1;
                    if k >= nt {
                        fail!(b"CURVE action is missing its curve name".to_vec(), ERR_ITEMS);
                    }
                    let ci = ctx.find_curve(&toks[k as usize]);
                    if ci < 0 {
                        fail!(cat(&[b"no curve named '", &toks[k as usize], b"' exists in the model"]), ERR_NAME, toks[k as usize].clone());
                    }
                    action.type_ = ActionType::Curve;
                    action.curve_idx = ci;
                } else if val_tok == b"TIMESERIES" {
                    k += 1;
                    if k >= nt {
                        fail!(b"TIMESERIES action is missing its series name".to_vec(), ERR_ITEMS);
                    }
                    let ti = ctx.find_timeseries(&toks[k as usize]);
                    if ti < 0 {
                        fail!(cat(&[b"no time series named '", &toks[k as usize], b"' exists in the model"]), ERR_NAME, toks[k as usize].clone());
                    }
                    action.type_ = ActionType::Timeseries;
                    action.tseries_idx = ti;
                } else if val_tok == b"PID" {
                    action.type_ = ActionType::Pid;
                    k += 1;
                    if k + 2 >= nt {
                        fail!(b"PID action requires three coefficients: Kp Ki Kd".to_vec(), ERR_ITEMS);
                    }
                    let mut kp = 0.0;
                    let mut ki = 0.0;
                    let mut kd = 0.0;
                    if !try_parse_double(&toks[k as usize], &mut kp) {
                        fail!(cat(&[b"PID coefficient Kp '", &toks[k as usize], b"' is not a number"]), ERR_NUMBER, toks[k as usize].clone());
                    }
                    if !try_parse_double(&toks[(k + 1) as usize], &mut ki) {
                        fail!(cat(&[b"PID coefficient Ki '", &toks[(k + 1) as usize], b"' is not a number"]), ERR_NUMBER, toks[(k + 1) as usize].clone());
                    }
                    if !try_parse_double(&toks[(k + 2) as usize], &mut kd) {
                        fail!(cat(&[b"PID coefficient Kd '", &toks[(k + 2) as usize], b"' is not a number"]), ERR_NUMBER, toks[(k + 2) as usize].clone());
                    }
                    let mut pid = PIDState { kp, ki, kd, ..PIDState::default() };
                    action.pid_idx = self.pid_states_.len() as i32;
                    pid.action_idx = action.pid_idx;
                    self.pid_states_.push(pid);
                    action.value = 0.0;
                } else {
                    action.type_ = ActionType::Numeric;
                    if attr == b"STATUS" {
                        if !parse_status_value(&val_tok, &mut action.value) {
                            fail!(cat(&[b"'", &toks[k as usize], b"' is not a valid STATUS (expected ON, OFF, OPEN or CLOSED)"]), ERR_KEYWORD, toks[k as usize].clone());
                        }
                    } else {
                        let mut v = 0.0;
                        if !try_parse_double(&toks[k as usize], &mut v) {
                            fail!(cat(&[b"action value '", &toks[k as usize], b"' is not a number"]), ERR_NUMBER, toks[k as usize].clone());
                        }
                        action.value = v;
                    }
                }

                if state == ParseState::InThen {
                    current_rule.then_actions.push(action);
                } else {
                    current_rule.else_actions.push(action);
                }
                continue;
            }
        }

        // Finish the last rule
        if !current_rule.name.is_empty() {
            self.rules_.push(std::mem::take(&mut current_rule));
            rules_parsed += 1;
        }

        // Re-index rule_idx on all actions
        for r in 0..self.rules_.len() {
            for a in self.rules_[r].then_actions.iter_mut() {
                a.rule_idx = r as i32;
            }
            for a in self.rules_[r].else_actions.iter_mut() {
                a.rule_idx = r as i32;
            }
        }

        self.build_premise_soa();

        Ok(rules_parsed)
    }
}
