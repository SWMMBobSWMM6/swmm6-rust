// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/handlers/OptionsHandler.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// [OPTIONS]: key-value rows into ctx.options. Unknown keys go to options.ext_options with
// a warning; a value outside a legacy option's keyword list refuses the deck with
// ERR_KEYWORD (205), as legacy project_readOption does.
//
// Upstream reads several values with the C library (strtod, atof, atoi, sscanf). Those
// calls are reproduced below as glibc 2.39 performs them in the C locale, on the value's
// c_str() (so an embedded NUL ends it):
//   * c_strtod: leading C-locale white space, an optional sign, then INF/INFINITY,
//     NAN/NAN(n-char-sequence) (the sequence becomes the NaN payload as glibc sets it),
//     a hexadecimal float, or a decimal number whose exponent is taken only when
//     complete. Decimal values are correctly rounded, as glibc's are. No conversion
//     consumes nothing (end == start, not past the white space).
//   * c_atoi: strtol base 10 saturating at the 64-bit long range, then truncated to int.
//   * Scanf: the "%d", "%c", "%3s", "%lf" and "%n" conversions and literal characters
//     of the formats used here, with glibc's return value (EOF on an input failure before
//     the first conversion) and its greedy float scan (which consumes "1e" whole and
//     converts the "1").
// `static_cast<int>` / `static_cast<long>` of a double that may be NaN or out of range
// use cvt_i32 / cvt_i64 (the reference build's cvttsd2si result).

#![allow(clippy::all)]

use crate::core::charconv::from_chars_double;
use crate::core::cmath::{cmax, cvt_i32};
use crate::core::context::SimulationContext;
use crate::core::datetime;
use crate::core::error_codes::{format_error, ERR_DATETIME, ERR_KEYWORD, ERR_NUMBER};
use crate::core::options::fv::{Backend, Limiter, RiemannSolver, ScalarScheme, StructureCoupling, TimeIntegration};
use crate::core::options::{Bytes, FlowUnits, InfiltrationModel, NodeContinuity, QualitySolverKind, RoutingModel};
use crate::input::parse_utils::{parse_time_seconds, to_double};
use crate::input::tokenizer;
use crate::core::cstdlib::{atof as c_atof, atoi as c_atoi, c_isspace, c_str, strtod as c_strtod};

// ============================================================================
// The C library, as the reference build runs it
// ============================================================================

/// C-locale `isspace`.
#[inline]



/// `static_cast<long>(x)` as the reference build executes it (64-bit `cvttsd2si`).
fn cvt_i64(x: f64) -> i64 {
    if x.is_nan() || x >= 9223372036854775808.0 || x < -9223372036854775808.0 {
        i64::MIN
    } else {
        x as i64
    }
}






/// Why a scanf conversion stopped: an input failure (end of input) or a matching failure.
enum ScanFail {
    Input,
    Match,
}

/// A glibc `sscanf` over one string, one conversion at a time. `pos` is `%n`.
struct Scanf<'a> {
    s: &'a [u8],
    pos: usize,
}

impl<'a> Scanf<'a> {
    fn new(s: &'a [u8]) -> Self {
        Scanf { s: c_str(s), pos: 0 }
    }
    fn peek(&self) -> Option<u8> {
        self.s.get(self.pos).copied()
    }
    /// The white space skip before %d, %s and %lf.
    fn skip_space(&mut self) -> Result<(), ScanFail> {
        while let Some(c) = self.peek() {
            if !c_isspace(c) {
                return Ok(());
            }
            self.pos += 1;
        }
        Err(ScanFail::Input)
    }
    /// A literal format character.
    fn lit(&mut self, want: u8) -> Result<(), ScanFail> {
        match self.peek() {
            None => Err(ScanFail::Input),
            Some(c) if c == want => {
                self.pos += 1;
                Ok(())
            }
            Some(_) => Err(ScanFail::Match),
        }
    }
    /// %d
    fn int(&mut self) -> Result<i32, ScanFail> {
        self.skip_space()?;
        let start = self.pos;
        if let Some(c) = self.peek() {
            if c == b'-' || c == b'+' {
                self.pos += 1;
            }
        }
        let ds = self.pos;
        while let Some(c) = self.peek() {
            if !c.is_ascii_digit() {
                break;
            }
            self.pos += 1;
        }
        if self.pos == ds {
            return Err(ScanFail::Match);
        }
        Ok(c_atoi(&self.s[start..self.pos]))
    }
    /// %c
    fn chr(&mut self) -> Result<u8, ScanFail> {
        match self.peek() {
            None => Err(ScanFail::Input),
            Some(c) => {
                self.pos += 1;
                Ok(c)
            }
        }
    }
    /// %3s (into a char[4], NUL-terminated)
    fn str3(&mut self) -> Result<Bytes, ScanFail> {
        self.skip_space()?;
        let mut out = Vec::new();
        while out.len() < 3 {
            match self.peek() {
                Some(c) if !c_isspace(c) => {
                    out.push(c);
                    self.pos += 1;
                }
                _ => break,
            }
        }
        Ok(out)
    }
    /// %lf: glibc's character collection, then strtod over what was collected.
    fn dbl(&mut self) -> Result<f64, ScanFail> {
        self.skip_space()?;
        let mut buf: Vec<u8> = Vec::new();
        let mut got_sign = false;
        let mut got_digit = false;
        let mut got_dot = false;
        let mut got_e = false;
        let mut hexa = false;
        // c = inchar()
        let mut c = self.next();
        if c == Some(b'-') || c == Some(b'+') {
            got_sign = true;
            buf.push(c.unwrap());
            c = self.next();
            if c.is_none() {
                return Err(ScanFail::Match);
            }
        }
        let lower = |c: Option<u8>| c.map(|x| x.to_ascii_lowercase());
        let mut special = false;
        if lower(c) == Some(b'n') {
            buf.push(c.unwrap());
            c = self.next();
            if lower(c) != Some(b'a') {
                return Err(ScanFail::Match);
            }
            buf.push(c.unwrap());
            c = self.next();
            if lower(c) != Some(b'n') {
                return Err(ScanFail::Match);
            }
            buf.push(c.unwrap());
            special = true;
        } else if lower(c) == Some(b'i') {
            buf.push(c.unwrap());
            for want in [b'n', b'f'] {
                c = self.next();
                if lower(c) != Some(want) {
                    return Err(ScanFail::Match);
                }
                buf.push(c.unwrap());
            }
            c = self.next();
            if c.is_some() {
                if lower(c) == Some(b'i') {
                    buf.push(c.unwrap());
                    for want in [b'n', b'i', b't', b'y'] {
                        c = self.next();
                        if lower(c) != Some(want) {
                            return Err(ScanFail::Match);
                        }
                        buf.push(c.unwrap());
                    }
                } else {
                    self.unget(c);
                }
            }
            special = true;
        }
        if !special {
            let mut exp_char = b'e';
            if c == Some(b'0') {
                buf.push(b'0');
                c = self.next();
                if lower(c) == Some(b'x') {
                    buf.push(c.unwrap());
                    hexa = true;
                    exp_char = b'p';
                    c = self.next();
                } else {
                    got_digit = true;
                }
            }
            loop {
                match c {
                    Some(ch) if ch.is_ascii_digit() => {
                        buf.push(ch);
                        got_digit = true;
                    }
                    Some(ch) if !got_e && hexa && ch.is_ascii_hexdigit() => {
                        buf.push(ch);
                        got_digit = true;
                    }
                    Some(ch) if got_e && buf.last() == Some(&exp_char) && (ch == b'-' || ch == b'+') => {
                        buf.push(ch);
                    }
                    Some(ch) if got_digit && !got_e && ch.to_ascii_lowercase() == exp_char => {
                        buf.push(exp_char);
                        got_e = true;
                        got_dot = true;
                    }
                    Some(b'.') if !got_dot => {
                        buf.push(b'.');
                        got_dot = true;
                    }
                    _ => {
                        self.unget(c);
                        break;
                    }
                }
                c = self.next();
                if c.is_none() {
                    break;
                }
            }
            if buf.len() == got_sign as usize || (hexa && buf.len() == 2 + got_sign as usize) {
                return Err(ScanFail::Match);
            }
        }
        let (d, used) = c_strtod(&buf);
        if used == 0 {
            return Err(ScanFail::Match);
        }
        Ok(d)
    }
    fn next(&mut self) -> Option<u8> {
        let c = self.peek();
        if c.is_some() {
            self.pos += 1;
        }
        c
    }
    fn unget(&mut self, c: Option<u8>) {
        if c.is_some() {
            self.pos -= 1;
        }
    }
}

/// The value sscanf returns after `done` assignments and the failure `f`.
fn scanf_result(done: i32, f: ScanFail) -> i32 {
    match f {
        ScanFail::Input if done == 0 => -1, // EOF
        _ => done,
    }
}

/// `sscanf(s, "%d%c%d%c%d", &a, &c1, &b, &c2, &y)`.
fn sscanf_d_c_d_c_d(s: &[u8], a: &mut i32, c1: &mut u8, b: &mut i32, c2: &mut u8, y: &mut i32) -> i32 {
    let mut sc = Scanf::new(s);
    let mut done = 0;
    macro_rules! step {
        ($e:expr, $dst:expr) => {
            match $e {
                Ok(v) => {
                    *$dst = v;
                    done += 1;
                }
                Err(f) => return scanf_result(done, f),
            }
        };
    }
    step!(sc.int(), a);
    step!(sc.chr(), c1);
    step!(sc.int(), b);
    step!(sc.chr(), c2);
    step!(sc.int(), y);
    done
}

/// `sscanf(s, "%3s%c%d%c%d", text, &c1, &b, &c2, &y)`.
fn sscanf_3s_c_d_c_d(s: &[u8], text: &mut Bytes, c1: &mut u8, b: &mut i32, c2: &mut u8, y: &mut i32) -> i32 {
    let mut sc = Scanf::new(s);
    let mut done = 0;
    macro_rules! step {
        ($e:expr, $dst:expr) => {
            match $e {
                Ok(v) => {
                    *$dst = v;
                    done += 1;
                }
                Err(f) => return scanf_result(done, f),
            }
        };
    }
    step!(sc.str3(), text);
    step!(sc.chr(), c1);
    step!(sc.int(), b);
    step!(sc.chr(), c2);
    step!(sc.int(), y);
    done
}

/// `sscanf(s, "%d:%d:%d", &h, &m, &sec)`.
fn sscanf_hms(s: &[u8], h: &mut i32, m: &mut i32, sec: &mut i32) -> i32 {
    let mut sc = Scanf::new(s);
    let mut done = 0;
    macro_rules! step {
        ($e:expr, $dst:expr) => {
            match $e {
                Ok(v) => {
                    *$dst = v;
                    done += 1;
                }
                Err(f) => return scanf_result(done, f),
            }
        };
    }
    macro_rules! lit {
        ($c:expr) => {
            if let Err(f) = sc.lit($c) {
                return scanf_result(done, f);
            }
        };
    }
    step!(sc.int(), h);
    lit!(b':');
    step!(sc.int(), m);
    lit!(b':');
    step!(sc.int(), sec);
    done
}

/// `sscanf(s, "%d:%d:%lf%n", &h, &m, &sec, &n)` (`with_seconds`) or
/// `sscanf(s, "%d:%d%n", &h, &m, &n)`.
fn sscanf_clock_n(s: &[u8], with_seconds: bool, h: &mut i32, m: &mut i32, sec: &mut f64, n: &mut i32) -> i32 {
    let mut sc = Scanf::new(s);
    let mut done = 0;
    macro_rules! step {
        ($e:expr, $dst:expr) => {
            match $e {
                Ok(v) => {
                    *$dst = v;
                    done += 1;
                }
                Err(f) => return scanf_result(done, f),
            }
        };
    }
    macro_rules! lit {
        ($c:expr) => {
            if let Err(f) = sc.lit($c) {
                return scanf_result(done, f);
            }
        };
    }
    step!(sc.int(), h);
    lit!(b':');
    step!(sc.int(), m);
    if with_seconds {
        lit!(b':');
        step!(sc.dbl(), sec);
    }
    *n = sc.pos as i32;
    done
}

/// `sscanf(s, "%d:%d:%lf%c", &h, &m, &sec, &tail)`.
fn sscanf_clock_tail(s: &[u8], h: &mut i32, m: &mut i32, sec: &mut f64, tail: &mut u8) -> i32 {
    let mut sc = Scanf::new(s);
    let mut done = 0;
    macro_rules! step {
        ($e:expr, $dst:expr) => {
            match $e {
                Ok(v) => {
                    *$dst = v;
                    done += 1;
                }
                Err(f) => return scanf_result(done, f),
            }
        };
    }
    macro_rules! lit {
        ($c:expr) => {
            if let Err(f) = sc.lit($c) {
                return scanf_result(done, f);
            }
        };
    }
    step!(sc.int(), h);
    lit!(b':');
    step!(sc.int(), m);
    lit!(b':');
    step!(sc.dbl(), sec);
    step!(sc.chr(), tail);
    done
}

// ============================================================================
// Helper: normalize uppercase key
// ============================================================================

/// Copy and upper-case an option token with the tokenizer's (ASCII) rules.
fn norm(sv: &[u8]) -> Bytes {
    tokenizer::to_upper(sv)
}

// legacy findmatch / match (input.c): a keyword matches when it is a
// case-insensitive PREFIX of the token — `STEADYFLOW` is STEADY, `DYNWAVE2`
// DYNWAVE. `tok` is already upper-cased by norm().
fn legacy_prefix(tok: &[u8], keyword: &[u8]) -> bool {
    let n = keyword.len();
    tok.len() >= n && &tok[..n] == keyword
}

// legacy OptionWords (keywords.c): project_readOption warns only for a key
// none of these prefixes.
fn legacy_knows_option(key: &[u8]) -> bool {
    const K_WORDS: &[&[u8]] = &[
        b"FLOW_UNITS",
        b"INFILTRATION",
        b"FLOW_ROUTING",
        b"START_DATE",
        b"START_TIME",
        b"END_DATE",
        b"END_TIME",
        b"REPORT_START_DATE",
        b"REPORT_START_TIME",
        b"SWEEP_START",
        b"SWEEP_END",
        b"DRY_DAYS",
        b"WET_STEP",
        b"DRY_STEP",
        b"ROUTING_STEP",
        b"RULE_STEP",
        b"REPORT_STEP",
        b"ALLOW_PONDING",
        b"INERTIAL_DAMPING",
        b"SLOPE_WEIGHTING",
        b"VARIABLE_STEP",
        b"NORMAL_FLOW_LIMITED",
        b"LENGTHENING_STEP",
        b"MIN_SURFAREA",
        b"COMPATIBILITY",
        b"SKIP_STEADY_STATE",
        b"TEMPDIR",
        b"IGNORE_RAINFALL",
        b"FORCE_MAIN_EQUATION",
        b"LINK_OFFSETS",
        b"MIN_SLOPE",
        b"IGNORE_SNOWMELT",
        b"IGNORE_GROUNDWATER",
        b"IGNORE_ROUTING",
        b"IGNORE_QUALITY",
        b"MAX_TRIALS",
        b"HEAD_TOLERANCE",
        b"SYS_FLOW_TOL",
        b"LAT_FLOW_TOL",
        b"IGNORE_RDII",
        b"MINIMUM_STEP",
        b"THREADS",
        b"SURCHARGE_METHOD",
        b"OUTFALL_BACKFLOW_QUALITY",
    ];
    for w in K_WORDS {
        if legacy_prefix(key, w) {
            return true;
        }
    }
    false
}

// legacy keywords.c word lists for the [OPTIONS] values project_readOption
// resolves with findmatch, refusing the deck with ERR_KEYWORD when none
// matches. Where no value table follows, the option stores the word's index.
const K_FLOW_UNIT_WORDS: &[&[u8]] = &[b"CFS", b"GPM", b"MGD", b"CMS", b"LPS", b"MLD"];
const K_FLOW_UNIT_VALUES: [FlowUnits; 6] =
    [FlowUnits::Cfs, FlowUnits::Gpm, FlowUnits::Mgd, FlowUnits::Cms, FlowUnits::Lps, FlowUnits::Mld];
const K_INFIL_MODEL_WORDS: &[&[u8]] =
    &[b"HORTON", b"MODIFIED_HORTON", b"GREEN_AMPT", b"MODIFIED_GREEN_AMPT", b"CURVE_NUMBER"];
const K_INFIL_MODEL_VALUES: [InfiltrationModel; 5] = [
    InfiltrationModel::Horton,
    InfiltrationModel::ModHorton,
    InfiltrationModel::GreenAmpt,
    InfiltrationModel::ModGreenAmpt,
    InfiltrationModel::CurveNumber,
];
const K_LINK_OFFSET_WORDS: &[&[u8]] = &[b"DEPTH", b"ELEVATION"];
const K_FORCE_MAIN_EQN_WORDS: &[&[u8]] = &[b"H-W", b"D-W"];
const K_INERT_DAMPING_WORDS: &[&[u8]] = &[b"NONE", b"PARTIAL", b"FULL"];
// Index 3 is legacy NONE, which this engine calls NEITHER.
const K_NORMAL_FLOW_WORDS: &[&[u8]] = &[b"SLOPE", b"FROUDE", b"BOTH", b"NONE"];
const K_SURCHARGE_WORDS: &[&[u8]] = &[b"EXTRAN", b"SLOT"];
const K_NO_YES_WORDS: &[&[u8]] = &[b"NO", b"YES"];

/// The index of the first legacy keyword that is a prefix of `tok` (legacy findmatch), or
/// -1 when none is.
fn legacy_findmatch(tok: &[u8], words: &[&[u8]]) -> i32 {
    for (i, w) in words.iter().enumerate() {
        if legacy_prefix(tok, w) {
            return i as i32;
        }
    }
    -1
}

/// A legacy YES/NO option value: this engine's TRUE/FALSE/1/0, or legacy NoYesWords by
/// prefix. `out` is set only when the value is recognised.
fn parse_legacy_yes_no(val: &[u8], out: &mut bool) -> bool {
    let v = norm(val);
    if v == b"TRUE" || v == b"1" {
        *out = true;
        return true;
    }
    if v == b"FALSE" || v == b"0" {
        *out = false;
        return true;
    }
    let m = legacy_findmatch(&v, K_NO_YES_WORDS);
    if m < 0 {
        return false;
    }
    *out = m == 1;
    true
}

/// legacy getDouble: strtod must consume the complete value.
fn legacy_number(val: &[u8], out: &mut f64) -> bool {
    let (v, used) = c_strtod(val);
    *out = v;
    used == c_str(val).len()
}

/// A legacy numeric option with its own bounds; on refusal `out` is kept and ERR_NUMBER
/// recorded. (Upstream takes the context; here its error list, so that `out` can be a
/// field of ctx.options.) The defaults upstream gives are -inf, +inf, false.
fn read_legacy_number(errors: &mut Vec<Bytes>, val: &[u8], out: &mut f64, lower: f64, upper: f64, exclusive_upper: bool) {
    let mut value = 0.0;
    if !legacy_number(val, &mut value) || value < lower || value > upper || (exclusive_upper && value == upper) {
        errors.push(format_error(ERR_NUMBER, val));
    } else {
        *out = value;
    }
}

/// datetime_strToDate in the legacy input reader's M/D/Y format.
fn legacy_date(val: &[u8], out: &mut f64) -> bool {
    if !val.iter().any(|&c| c == b'/' || c == b'-') {
        return false;
    }
    let (mut year, mut month, mut day) = (0i32, 0i32, 0i32);
    let (mut sep1, mut sep2) = (0u8, 0u8);
    let mut month_text: Bytes = Vec::new();
    if sscanf_d_c_d_c_d(val, &mut month, &mut sep1, &mut day, &mut sep2, &mut year) < 3 {
        month = 0;
        if sscanf_3s_c_d_c_d(val, &mut month_text, &mut sep1, &mut day, &mut sep2, &mut year) < 3 {
            return false;
        }
    }
    if month == 0 {
        const MONTHS: [&[u8]; 12] =
            [b"JAN", b"FEB", b"MAR", b"APR", b"MAY", b"JUN", b"JUL", b"AUG", b"SEP", b"OCT", b"NOV", b"DEC"];
        for i in 0..12 {
            if norm(&month_text) == MONTHS[i] {
                month = i as i32 + 1;
                break;
            }
        }
    }
    let date = datetime::encode_date(year, month, day);
    if date == -(datetime::DATE_DELTA as f64) {
        return false;
    }
    *out = date;
    true
}

/// legacy datetime_strToTime: decimal hours or permissive H[:M[:S]].
fn legacy_time(val: &[u8], days: &mut f64) -> bool {
    let mut hours = 0.0;
    if legacy_number(val, &mut hours) {
        *days = hours / 24.0;
        return true;
    }
    let (mut h, mut m, mut s) = (0i32, 0i32, 0i32);
    if sscanf_hms(val, &mut h, &mut m, &mut s) < 1 || h < 0 || m < 0 || s < 0 {
        return false;
    }
    *days = (3600.0 * h as f64 + 60.0 * m as f64 + s as f64) / datetime::SECS_PER_DAY;
    true
}

/// A legacy time step in whole seconds, multi-day durations included.
fn legacy_step_seconds(days: f64) -> f64 {
    let (h, m, s) = datetime::decode_time(days);
    s as f64 + 60.0 * m as f64 + 3600.0 * h as f64 + 86400.0 * days.trunc()
}

/// `std::numeric_limits<int>::max() / datetime::SecsPerDay`
const MAX_STEP_DAYS: f64 = i32::MAX as f64 / datetime::SECS_PER_DAY;

/// `opt.ext_options[key] = val`
fn ext_set(ctx: &mut SimulationContext, key: &[u8], val: &[u8]) {
    ctx.options.ext_options.insert(key.to_vec(), val.to_vec());
}

/// `openswmm::from_chars_double(val, out)`: `out` is written only on success.
fn from_chars_into(val: &[u8], out: &mut f64) {
    if let Ok(v) = from_chars_double(val).0 {
        *out = v;
    }
}

/// `ctx.warnings.push_back(w); if (ctx.warning_code == 0) ctx.warning_code = 101;`
fn warn_101(ctx: &mut SimulationContext, w: Bytes) {
    ctx.warnings.push(w);
    if ctx.warning_code == 0 {
        ctx.warning_code = 101;
    }
}

// ============================================================================
// handle_options() — registered as built-in handler for "OPTIONS"
// ============================================================================

pub fn handle_options(ctx: &mut SimulationContext, lines: &[Bytes]) {
    let inf = f64::INFINITY;

    // We need to accumulate START_DATE/START_TIME because they may appear
    // in any order relative to each other.
    let (mut start_date_part, mut start_time_part) = (0.0f64, 0.0f64);
    let (mut end_date_part, mut end_time_part) = (0.0f64, 0.0f64);
    let (mut rpt_date_part, mut rpt_time_part) = (0.0f64, 0.0f64);
    let (mut got_start_date, mut got_start_time) = (false, false);
    let (mut got_end_date, mut got_end_time) = (false, false);
    let (mut got_rpt_date, mut got_rpt_time) = (false, false);

    for raw in lines {
        let tokens = tokenizer::tokenize(raw);
        if tokens.len() < 2 {
            continue;
        }

        let key = norm(&tokens[0]);
        let key: &[u8] = &key;
        // Value may be a single token or a multi-token (e.g., quoted CRS string).
        let val: &[u8] = &tokens[1];
        let opt = &mut ctx.options;

        // -----------------------------------------------------------------
        // Flow / routing units
        // -----------------------------------------------------------------
        if key == b"FLOW_UNITS" {
            let uv = norm(val);
            if uv == b"CFS" {
                opt.flow_units = FlowUnits::Cfs;
            } else if uv == b"GPM" {
                opt.flow_units = FlowUnits::Gpm;
            } else if uv == b"MGD" {
                opt.flow_units = FlowUnits::Mgd;
            } else if uv == b"CMS" {
                opt.flow_units = FlowUnits::Cms;
            } else if uv == b"LPS" {
                opt.flow_units = FlowUnits::Lps;
            } else if uv == b"MLD" {
                opt.flow_units = FlowUnits::Mld;
            }
            // Legacy findmatch also takes a keyword followed by anything (LPSX
            // reads as LPS) and refuses the rest.
            else {
                let m = legacy_findmatch(&uv, K_FLOW_UNIT_WORDS);
                if m >= 0 {
                    opt.flow_units = K_FLOW_UNIT_VALUES[m as usize];
                } else {
                    ctx.errors.push(format_error(ERR_KEYWORD, val));
                }
            }
        } else if key == b"INFILTRATION" {
            let iv = norm(val);
            if iv == b"HORTON" {
                opt.infiltration = InfiltrationModel::Horton;
            } else if iv == b"MOD_HORTON" || iv == b"MODIFIED_HORTON" {
                opt.infiltration = InfiltrationModel::ModHorton;
            } else if iv == b"GREEN_AMPT" {
                opt.infiltration = InfiltrationModel::GreenAmpt;
            }
            // Legacy InfilModelWords spells it MODIFIED_GREEN_AMPT;
            // MOD_GREEN_AMPT is this engine's own alias.
            else if iv == b"MOD_GREEN_AMPT" || iv == b"MODIFIED_GREEN_AMPT" {
                opt.infiltration = InfiltrationModel::ModGreenAmpt;
            } else if iv == b"CURVE_NUMBER" {
                opt.infiltration = InfiltrationModel::CurveNumber;
            } else {
                let m = legacy_findmatch(&iv, K_INFIL_MODEL_WORDS);
                if m >= 0 {
                    opt.infiltration = K_INFIL_MODEL_VALUES[m as usize];
                } else {
                    ctx.errors.push(format_error(ERR_KEYWORD, val));
                }
            }
        } else if key == b"FLOW_ROUTING" {
            let rv = norm(val);
            // legacy project_readOption ROUTE_MODEL: findmatch over
            // RouteModelWords {NONE, STEADY, KINWAVE, XKINWAVE, DYNWAVE}, then
            // OldRouteModelWords {NONE, NF, KW, EKW, DW}, each a prefix test
            // in that order.
            if legacy_prefix(&rv, b"NONE") {
                opt.ignore_routing = true;
            } else if legacy_prefix(&rv, b"STEADY") {
                opt.routing_model = RoutingModel::Steady;
            } else if legacy_prefix(&rv, b"KINWAVE") {
                opt.routing_model = RoutingModel::Kinwave;
            } else if legacy_prefix(&rv, b"XKINWAVE") {
                opt.routing_model = RoutingModel::Kinwave;
            } else if legacy_prefix(&rv, b"DYNWAVE") {
                opt.routing_model = RoutingModel::Dynwave;
            } else if legacy_prefix(&rv, b"NF") {
                opt.routing_model = RoutingModel::Steady;
            } else if legacy_prefix(&rv, b"KW") {
                opt.routing_model = RoutingModel::Kinwave;
            } else if legacy_prefix(&rv, b"EKW") {
                opt.routing_model = RoutingModel::Kinwave;
            } else if legacy_prefix(&rv, b"DW") {
                opt.routing_model = RoutingModel::Dynwave;
            }
            // This engine's own spellings.
            else if rv == b"KINEMATIC_WAVE" {
                opt.routing_model = RoutingModel::Kinwave;
            } else if rv == b"DYNAMIC_WAVE" {
                opt.routing_model = RoutingModel::Dynwave;
            } else if rv == b"FV" || rv == b"FINITE_VOLUME" {
                opt.routing_model = RoutingModel::Fv;
            } else if rv == b"NO_ROUTING" {
                opt.ignore_routing = true;
            } else {
                ctx.errors.push(format_error(ERR_KEYWORD, val));
            }
        } else if key == b"QUALITY_SOLVER" {
            let qv = norm(val);
            if qv == b"LEGACY" {
                opt.quality_solver = QualitySolverKind::Legacy;
            } else if qv == b"EULERIAN_ARD" || qv == b"ARD" {
                opt.quality_solver = QualitySolverKind::EulerianArd;
            } else if qv == b"LAGRANGIAN" || qv == b"LARD" {
                opt.quality_solver = QualitySolverKind::Lagrangian;
            } else {
                ext_set(ctx, key, val);
            }
        } else if key == b"WATER_AGE" {
            opt.water_age = norm(val) == b"ON" || norm(val) == b"YES";
        } else if key == b"OUTFALL_BACKFLOW_QUALITY" {
            let bv = norm(val);
            if legacy_prefix(&bv, b"LAST") {
                opt.outfall_backflow_zero = false;
            } else if legacy_prefix(&bv, b"ZERO") {
                opt.outfall_backflow_zero = true;
            } else {
                ctx.errors.push(format_error(ERR_KEYWORD, val));
            }
        } else if key == b"HEAT_TRANSPORT" {
            opt.heat_transport = norm(val) == b"ON" || norm(val) == b"YES";
        } else if key == b"WATER_DENSITY" {
            from_chars_into(val, &mut opt.water_density);
        } else if key == b"WATER_SPECIFIC_HEAT_CAPACITY" {
            from_chars_into(val, &mut opt.water_specific_heat);
        } else if key == b"WIND_FUNC_COEFF_A" {
            from_chars_into(val, &mut opt.wind_func_coeff_a);
        } else if key == b"WIND_FUNC_COEFF_B" {
            from_chars_into(val, &mut opt.wind_func_coeff_b);
        } else if key == b"PRESSURE_RATIO" {
            from_chars_into(val, &mut opt.pressure_ratio);

            // -----------------------------------------------------------------
            // Timesteps
            // -----------------------------------------------------------------
        } else if key == b"ROUTING_STEP" {
            let (mut step, mut days) = (0.0f64, 0.0f64);
            if !legacy_number(val, &mut step) {
                if !legacy_time(val, &mut days) || !days.is_finite() || days.abs() > MAX_STEP_DAYS {
                    ctx.errors.push(format_error(ERR_NUMBER, val));
                    continue;
                }
                step = legacy_step_seconds(days);
            }
            if step <= 0.0 {
                ctx.errors.push(format_error(ERR_NUMBER, val));
            } else {
                opt.routing_step = step;
            }
        } else if key == b"MINIMUM_STEP" {
            if !val.contains(&b':') {
                read_legacy_number(&mut ctx.errors, val, &mut opt.min_routing_step, 0.0, inf, false);
            } else {
                // GUI/API clock spelling, including fractional seconds.
                let (mut hours, mut minutes, mut consumed) = (0i32, 0i32, 0i32);
                let mut seconds = 0.0f64;
                let mut fields = sscanf_clock_n(val, true, &mut hours, &mut minutes, &mut seconds, &mut consumed);
                if fields == 2 {
                    fields = sscanf_clock_n(val, false, &mut hours, &mut minutes, &mut seconds, &mut consumed);
                }
                if fields < 2
                    || consumed != val.len() as i32
                    || hours < 0
                    || minutes < 0
                    || minutes >= 60
                    || !seconds.is_finite()
                    || seconds < 0.0
                    || seconds >= 60.0
                {
                    ctx.errors.push(format_error(ERR_NUMBER, val));
                } else {
                    opt.min_routing_step = hours as f64 * 3600.0 + minutes as f64 * 60.0 + seconds;
                }
            }
        } else if key == b"DRY_DAYS" {
            // project_readOption deliberately uses atof, accepting numeric prefixes and junk as 0.
            let value = c_atof(val);
            if value < 0.0 {
                ctx.errors.push(format_error(ERR_NUMBER, val));
            } else {
                opt.dry_days = value;
            }
        } else if key == b"DRY_STEP" || key == b"WET_STEP" || key == b"REPORT_STEP" || key == b"RULE_STEP" {
            // (upstream: `double days = 0.0, step = 0.0;` -- every path assigns step)
            let mut days = 0.0f64;
            let step: f64;
            // Explicit fractional seconds extend the clock spelling without
            // changing legacy decimal-hour tokens or RULE_STEP truncation.
            let last_colon = val.iter().rposition(|&c| c == b':');
            let fractional_clock = key != b"RULE_STEP"
                && match last_colon {
                    Some(lc) => val[lc..].contains(&b'.'),
                    None => false,
                };
            if fractional_clock {
                let (mut h, mut m) = (0i32, 0i32);
                let mut seconds = 0.0f64;
                let mut tail = 0u8;
                if sscanf_clock_tail(val, &mut h, &mut m, &mut seconds, &mut tail) != 3
                    || h < 0
                    || m < 0
                    || seconds < 0.0
                    || !seconds.is_finite()
                {
                    ctx.errors.push(format_error(ERR_NUMBER, val));
                    continue;
                }
                step = 3600.0 * h as f64 + 60.0 * m as f64 + seconds;
                if !step.is_finite() || step > i32::MAX as f64 {
                    ctx.errors.push(format_error(ERR_NUMBER, val));
                    continue;
                }
            } else {
                if !legacy_time(val, &mut days) || !days.is_finite() || days.abs() > MAX_STEP_DAYS {
                    ctx.errors.push(format_error(ERR_DATETIME, val));
                    continue;
                }
                step = legacy_step_seconds(days);
            }
            if if key == b"RULE_STEP" { step < 0.0 } else { step <= 0.0 } {
                ctx.errors.push(format_error(ERR_NUMBER, val));
                continue;
            }
            if key == b"DRY_STEP" {
                opt.dry_step = step;
            } else if key == b"WET_STEP" {
                opt.wet_step = step;
            } else if key == b"REPORT_STEP" {
                opt.report_step = step;
            } else {
                opt.rule_step = step;
            }
        } else if key == b"QUALITY_STEP" {
            opt.quality_step = parse_time_seconds(val);
        } else if key == b"MAX_SEGMENTS_PER_LINK" {
            let mut d = 0.0;
            from_chars_into(val, &mut d);
            opt.max_segments_per_link = cvt_i32(d);
        } else if key == b"DISPERSION" {
            let dv = norm(val);
            if dv == b"RWPT" {
                opt.lard_rwpt = true;
            } else if dv == b"OFF" || dv == b"NONE" {
                opt.lard_rwpt = false;
            } else {
                ext_set(ctx, key, val);
            }
        } else if key == b"RWPT_SEED" {
            let mut d = 0.0;
            from_chars_into(val, &mut d);
            opt.rwpt_seed = cvt_i32(d);

            // -----------------------------------------------------------------
            // Simulation dates / times
            // -----------------------------------------------------------------
        } else if key == b"START_DATE" {
            if legacy_date(val, &mut start_date_part) {
                got_start_date = true;
            } else {
                ctx.errors.push(format_error(ERR_DATETIME, val));
            }
        } else if key == b"START_TIME" {
            if legacy_time(val, &mut start_time_part) {
                got_start_time = true;
            } else {
                ctx.errors.push(format_error(ERR_DATETIME, val));
            }
        } else if key == b"END_DATE" {
            if legacy_date(val, &mut end_date_part) {
                got_end_date = true;
            } else {
                ctx.errors.push(format_error(ERR_DATETIME, val));
            }
        } else if key == b"END_TIME" {
            if legacy_time(val, &mut end_time_part) {
                got_end_time = true;
            } else {
                ctx.errors.push(format_error(ERR_DATETIME, val));
            }
        } else if key == b"REPORT_START_DATE" {
            if legacy_date(val, &mut rpt_date_part) {
                got_rpt_date = true;
            } else {
                ctx.errors.push(format_error(ERR_DATETIME, val));
            }
        } else if key == b"REPORT_START_TIME" {
            if legacy_time(val, &mut rpt_time_part) {
                got_rpt_time = true;
            } else {
                ctx.errors.push(format_error(ERR_DATETIME, val));
            }

            // -----------------------------------------------------------------
            // Solver settings
            // -----------------------------------------------------------------
        } else if key == b"MAX_TRIALS" {
            // Legacy atoi accepts a numeric prefix (and maps nonnumeric text to zero).
            let value = c_atoi(val);
            if value < 0 {
                ctx.errors.push(format_error(ERR_NUMBER, val));
            } else {
                opt.max_trials = value;
            }
        } else if key == b"HEAD_TOLERANCE" {
            read_legacy_number(&mut ctx.errors, val, &mut opt.head_tol, -inf, inf, false);
        } else if key == b"SYS_FLOW_TOL" {
            let mut pct = 0.0;
            if !legacy_number(val, &mut pct) {
                ctx.errors.push(format_error(ERR_NUMBER, val));
            } else {
                opt.sys_flow_tol = pct / 100.0; // input is percent, store as fraction
            }
        } else if key == b"LAT_FLOW_TOL" {
            let mut pct = 0.0;
            if !legacy_number(val, &mut pct) {
                ctx.errors.push(format_error(ERR_NUMBER, val));
            } else {
                opt.lat_flow_tol = pct / 100.0; // input is percent, store as fraction
            }
        } else if key == b"VARIABLE_STEP" {
            read_legacy_number(&mut ctx.errors, val, &mut opt.variable_step, 0.0, 2.0, false);
        } else if key == b"THREADS" {
            // Legacy atoi accepts a numeric prefix (and maps nonnumeric text to zero).
            let value = c_atoi(val);
            if value < 0 {
                ctx.errors.push(format_error(ERR_NUMBER, val));
            } else {
                opt.num_threads = value;
            }

            // -----------------------------------------------------------------
            // Boolean flags
            // -----------------------------------------------------------------
        } else if key == b"ALLOW_PONDING" {
            if !parse_legacy_yes_no(val, &mut opt.allow_ponding) {
                ctx.errors.push(format_error(ERR_KEYWORD, val));
            }
        } else if key == b"IGNORE_RAINFALL" {
            if !parse_legacy_yes_no(val, &mut opt.ignore_rainfall) {
                ctx.errors.push(format_error(ERR_KEYWORD, val));
            }
        } else if key == b"IGNORE_SNOWMELT" {
            if !parse_legacy_yes_no(val, &mut opt.ignore_snow_melt) {
                ctx.errors.push(format_error(ERR_KEYWORD, val));
            }
        } else if key == b"IGNORE_GROUNDWATER" {
            if !parse_legacy_yes_no(val, &mut opt.ignore_groundwater) {
                ctx.errors.push(format_error(ERR_KEYWORD, val));
            }
        } else if key == b"IGNORE_RDII" {
            if !parse_legacy_yes_no(val, &mut opt.ignore_rdii) {
                ctx.errors.push(format_error(ERR_KEYWORD, val));
            }
        } else if key == b"IGNORE_ROUTING" {
            if !parse_legacy_yes_no(val, &mut opt.ignore_routing) {
                ctx.errors.push(format_error(ERR_KEYWORD, val));
            }
        } else if key == b"IGNORE_QUALITY" {
            if !parse_legacy_yes_no(val, &mut opt.ignore_quality) {
                ctx.errors.push(format_error(ERR_KEYWORD, val));
            }
        } else if key == b"SLOPE_WEIGHTING" {
            // Legacy validates this retired option even though it has no solver effect.
            let mut value = false;
            if !parse_legacy_yes_no(val, &mut value) {
                ctx.errors.push(format_error(ERR_KEYWORD, val));
            } else {
                ext_set(ctx, key, val);
            }
        } else if key == b"IGNORE_2D" {
            opt.ignore_2d = tokenizer::parse_boolean(val);

            // -----------------------------------------------------------------
            // CRS (R06)
            // -----------------------------------------------------------------
        } else if key == b"CRS" {
            opt.crs = val.to_vec();
            ctx.spatial.crs = val.to_vec();
            // Heuristic: geographic if EPSG:4326 or contains "latlong"/"geographic"
            let lv = norm(val);
            ctx.spatial.is_geographic =
                lv == b"EPSG:4326" || contains(&lv, b"LATLONG") || contains(&lv, b"GEOGRAPHIC");

            // -----------------------------------------------------------------
            // Legacy-compatibility keys that we silently accept but don't use
            // -----------------------------------------------------------------
        } else if key == b"SWEEP_START" || key == b"SWEEP_END" {
            // Same non-leap anchor and 24-character input limit as legacy project_readOption.
            let mut date = 0.0;
            let mut s24 = val[..val.len().min(24)].to_vec();
            s24.extend_from_slice(b"/1947");
            if !legacy_date(&s24, &mut date) {
                ctx.errors.push(format_error(ERR_DATETIME, val));
            } else if key == b"SWEEP_START" {
                opt.sweep_start = datetime::day_of_year(date);
            } else {
                opt.sweep_end = datetime::day_of_year(date);
            }
        } else if key == b"NORMAL_FLOW_LIMITED" {
            let nv = norm(val);
            if nv == b"SLOPE" {
                opt.normal_flow_ltd = 0;
            } else if nv == b"FROUDE" {
                opt.normal_flow_ltd = 1;
            } else if nv == b"BOTH" {
                opt.normal_flow_ltd = 2;
            } else if nv == b"NEITHER" {
                opt.normal_flow_ltd = 3;
            }
            // Legacy prefix forms, and legacy NONE (index 3, NEITHER here).
            else {
                let m = legacy_findmatch(&nv, K_NORMAL_FLOW_WORDS);
                if m >= 0 {
                    opt.normal_flow_ltd = m;
                } else {
                    ctx.errors.push(format_error(ERR_KEYWORD, val));
                }
            }
        } else if key == b"FORCE_MAIN_EQUATION" {
            let fv = norm(val);
            if fv == b"H-W" || fv == b"HW" || fv == b"HAZEN-WILLIAMS" {
                opt.force_main_eqn = 0;
            } else if fv == b"D-W" || fv == b"DW" || fv == b"DARCY-WEISBACH" {
                opt.force_main_eqn = 1;
            } else {
                let m = legacy_findmatch(&fv, K_FORCE_MAIN_EQN_WORDS);
                if m >= 0 {
                    opt.force_main_eqn = m;
                } else {
                    ctx.errors.push(format_error(ERR_KEYWORD, val));
                }
            }
        } else if key == b"INERTIAL_DAMPING" {
            let iv = norm(val);
            if iv == b"NONE" {
                opt.inertial_damping = 0;
            } else if iv == b"PARTIAL" {
                opt.inertial_damping = 1;
            } else if iv == b"FULL" {
                opt.inertial_damping = 2;
            } else {
                let m = legacy_findmatch(&iv, K_INERT_DAMPING_WORDS);
                if m >= 0 {
                    opt.inertial_damping = m;
                } else {
                    ctx.errors.push(format_error(ERR_KEYWORD, val));
                }
            }
        } else if key == b"SURCHARGE_METHOD" {
            let sv = norm(val);
            if sv == b"EXTRAN" {
                opt.surcharge_method = 0;
            } else if sv == b"SLOT" {
                opt.surcharge_method = 1;
            } else if sv == b"DYNAMIC_SLOT" {
                opt.surcharge_method = 2;
            } else if sv == b"TPA" {
                opt.surcharge_method = 3; // issue #156
            } else {
                let m = legacy_findmatch(&sv, K_SURCHARGE_WORDS);
                if m >= 0 {
                    opt.surcharge_method = m;
                } else {
                    ctx.errors.push(format_error(ERR_KEYWORD, val));
                }
            }
        } else if key == b"TPA_CELERITY" {
            // issue #156
            opt.tpa_celerity = to_double(val, 0.0);
        } else if key == b"DPS_CELERITY" {
            opt.dps_target_celerity = to_double(val, 0.0);
        } else if key == b"DPS_ALPHA" {
            opt.dps_alpha = to_double(val, 0.0);
        } else if key == b"DPS_DECAY_TIME" {
            opt.dps_decay_time = to_double(val, 0.0);

            // Unsteady friction (issue #156). Accepted-and-INERT.
        } else if key == b"UNSTEADY_FRICTION" {
            let uv = norm(val);
            if uv == b"NONE" {
                opt.unsteady_friction = 0;
            } else if uv == b"VITKOVSKY" {
                opt.unsteady_friction = 1;
            }
        } else if key == b"UF_K3" {
            opt.uf_k3 = to_double(val, 0.0);
        } else if key == b"REPORT_SIGNED_HEADS" {
            // issue #156 O-6
            let sv2 = norm(val);
            opt.report_signed_heads = if sv2 == b"YES" || sv2 == b"TRUE" || sv2 == b"ON" || sv2 == b"1" { 1 } else { 0 };

            // -----------------------------------------------------------------
            // Explicit finite-volume solver (FLOW_ROUTING FV). Accepted and inert
            // under any other routing model.
            // -----------------------------------------------------------------
        } else if key == b"FV_CELL_LENGTH" {
            opt.fv.cell_length = to_double(val, 0.0);
        } else if key == b"FV_MIN_CELLS" {
            opt.fv.min_cells = cmax(1, cvt_i32(to_double(val, 0.0)));
        } else if key == b"FV_CFL" {
            opt.fv.cfl = to_double(val, 0.0);
        } else if key == b"FV_RIEMANN" {
            let rv2 = norm(val);
            if rv2 == b"HLL" {
                opt.fv.riemann = RiemannSolver::Hll;
            } else if rv2 == b"HLLC" {
                opt.fv.riemann = RiemannSolver::Hllc;
            }
        } else if key == b"FV_ORDER" {
            opt.fv.order = cvt_i32(to_double(val, 0.0));
        } else if key == b"FV_LIMITER" {
            let lv = norm(val);
            if lv == b"MINMOD" {
                opt.fv.limiter = Limiter::Minmod;
            } else if lv == b"VANLEER" {
                opt.fv.limiter = Limiter::Vanleer;
            } else if lv == b"SUPERBEE" {
                opt.fv.limiter = Limiter::Superbee;
            }
        } else if key == b"FV_SCALAR_SCHEME" {
            let sv2 = norm(val);
            if sv2 == b"UPWIND" {
                opt.fv.scalar_scheme = ScalarScheme::Upwind;
            } else if sv2 == b"MUSCL" {
                opt.fv.scalar_scheme = ScalarScheme::Muscl;
            } else if sv2 == b"QUICKEST_ULTIMATE" || sv2 == b"QUICKEST" {
                opt.fv.scalar_scheme = ScalarScheme::QuickestUltimate;
            }
        } else if key == b"FV_TIME_INTEGRATION" {
            let tv = norm(val);
            if tv == b"EULER" {
                opt.fv.time_integration = TimeIntegration::Euler;
            } else if tv == b"RK2" {
                opt.fv.time_integration = TimeIntegration::Rk2;
            }
        } else if key == b"FV_SLOT_CELERITY" {
            opt.fv.slot_celerity = to_double(val, 0.0);
        } else if key == b"FV_PRESSURIZED_IMPLICIT" {
            let pv = norm(val);
            opt.fv.pressurized_implicit = pv == b"YES" || pv == b"TRUE" || pv == b"ON" || pv == b"1";
        } else if key == b"FV_PRESSURE_CLOSURE" {
            // issue #156 Phase 4
            let pc = norm(val);
            if pc == b"SLOT" {
                opt.fv.pressure_closure = 0;
            } else if pc == b"TPA" {
                opt.fv.pressure_closure = 1;
            }
        } else if key == b"FV_DISPERSION" {
            opt.fv.dispersion = to_double(val, 0.0);
        } else if key == b"FV_STRUCTURE_COUPLING" {
            let cv = norm(val);
            if cv == b"SUBSTEP" {
                opt.fv.structure_coupling = StructureCoupling::Substep;
            } else if cv == b"ROUTING_STEP" {
                opt.fv.structure_coupling = StructureCoupling::RoutingStep;
            }
        } else if key == b"FV_NODE_COUPLING" {
            // RETIRED 2026-08-29 with FV_NODE_DT and FV_NODE_PICARD below.
            if norm(val) == b"EXPLICIT" {
                warn_101(
                    ctx,
                    b"WARNING: FV_NODE_COUPLING EXPLICIT is retired and will be \
treated as SEMI_IMPLICIT - storage nodes are always coupled \
semi-implicitly (plain junctions are algebraic interfaces \
under either)."
                        .to_vec(),
                );
            }
        } else if key == b"FV_COMPACTION" {
            let bv = norm(val);
            opt.fv.compaction = !(bv == b"NO" || bv == b"FALSE" || bv == b"0" || bv == b"OFF");
        } else if key == b"FV_BACKEND" {
            let bv = norm(val);
            if bv == b"CPU" {
                opt.fv.backend = Backend::Cpu;
            } else if bv == b"AUTO" {
                opt.fv.backend = Backend::Auto;
            } else if bv == b"OMP" {
                opt.fv.backend = Backend::Omp;
            } else if bv == b"CUDA" {
                opt.fv.backend = Backend::Cuda;
            } else if bv == b"HIP" {
                opt.fv.backend = Backend::Hip;
            } else if bv == b"SYCL" {
                opt.fv.backend = Backend::Sycl;
            }
        } else if key == b"FV_MIN_PARALLEL_CELLS" {
            opt.fv.min_parallel_cells = cvt_i64(to_double(val, 0.0));
        } else if key == b"FV_LTS" {
            let bv = norm(val);
            opt.fv.lts = !(bv == b"NO" || bv == b"FALSE" || bv == b"0" || bv == b"OFF");
        } else if key == b"FV_NODE_FEEDBACK_DT" {
            let bv = norm(val);
            opt.fv.node_feedback_dt = !(bv == b"NO" || bv == b"FALSE" || bv == b"0" || bv == b"OFF");
        } else if key == b"FV_LTS_MAX_TIERS" {
            opt.fv.lts_max_tiers = cmax(1, cvt_i32(to_double(val, 0.0)));
        } else if key == b"FV_CFL_CENSUS_INTERVAL" {
            opt.fv.cfl_census_interval = cmax(1, cvt_i32(to_double(val, 0.0)));
        } else if key == b"FV_NODE_DT" {
            // RETIRED 2026-08-29 -- see FV_NODE_COUPLING above.
            if norm(val) == b"NONE" {
                warn_101(
                    ctx,
                    b"WARNING: FV_NODE_DT NONE is retired and will be treated as \
STABILITY - the node accuracy bound is always armed (it \
binds storage nodes only; algebraic junctions and outfalls \
are exempt)."
                        .to_vec(),
                );
            }
        } else if key == b"FV_NODE_PICARD" {
            // RETIRED 2026-08-29 -- see FV_NODE_COUPLING above.
            if cmax(1, cvt_i32(to_double(val, 0.0))) > 1 {
                let mut w = b"WARNING: FV_NODE_PICARD ".to_vec();
                w.extend_from_slice(val);
                w.extend_from_slice(
                    b" is retired and will be \
treated as 1 - the semi-implicit node correction is always \
a single sweep.",
                );
                warn_101(ctx, w);
            }
        } else if key == b"FV_NODE_CELL_COUPLING" || key == b"FV_JUNCTION_MODEL" {
            // Retired options, accepted and ignored so existing project files
            // still parse.
        } else if key == b"NODE_CONTINUITY" {
            let nc = norm(val);
            if nc == b"EXPLICIT" {
                opt.node_continuity = NodeContinuity::Explicit;
            } else if nc == b"SEMI_IMPLICIT" {
                opt.node_continuity = NodeContinuity::SemiImplicit;
            }
        } else if key == b"VIRTUAL_JUNCTION_MOMENTUM" {
            // RETIRED 2026-08-14. FULL is accepted and warned for one release.
            let vm = norm(val);
            opt.virtual_junction_momentum = 0;
            if vm == b"FULL" {
                warn_101(
                    ctx,
                    b"WARNING: VIRTUAL_JUNCTION_MOMENTUM FULL is retired and \
will be treated as BASIC - the cross-junction momentum \
term was not mass conserving."
                        .to_vec(),
                );
            }
        } else if key == b"ANDERSON_ACCEL" {
            let av = norm(val);
            opt.anderson_accel = av == b"YES" || av == b"TRUE" || av == b"1";
        } else if key == b"LINK_OFFSETS" {
            let lv = norm(val);
            if lv == b"DEPTH" {
                opt.link_offsets = 0;
            } else if lv == b"ELEVATION" {
                opt.link_offsets = 1;
            } else {
                let m = legacy_findmatch(&lv, K_LINK_OFFSET_WORDS);
                if m >= 0 {
                    opt.link_offsets = m;
                } else {
                    ctx.errors.push(format_error(ERR_KEYWORD, val));
                }
            }
        } else if key == b"MIN_SLOPE" {
            read_legacy_number(&mut ctx.errors, val, &mut opt.min_slope, 0.0, 100.0, true);
        } else if key == b"MIN_SURFAREA" {
            read_legacy_number(&mut ctx.errors, val, &mut opt.min_surf_area, 0.0, inf, false);
        } else if key == b"LENGTHENING_STEP" {
            let (mut step, mut days) = (0.0f64, 0.0f64);
            if !legacy_number(val, &mut step) {
                if !legacy_time(val, &mut days) || !days.is_finite() || days.abs() > MAX_STEP_DAYS {
                    ctx.errors.push(format_error(ERR_NUMBER, val));
                    continue;
                }
                step = legacy_step_seconds(days);
            }
            opt.lengthening_step = cmax(0.0, step);
        } else if key == b"SKIP_STEADY_STATE" {
            if !parse_legacy_yes_no(val, &mut opt.skip_steady_state) {
                ctx.errors.push(format_error(ERR_KEYWORD, val));
            }
        } else if key == b"WRITE_ABSOLUTE_PATHS" {
            // Slice IO-4 portability opt-out — see SimulationOptions.hpp.
            opt.write_absolute_paths = tokenizer::parse_boolean(val);
        } else if key.starts_with(b"IO_") {
            ext_set(ctx, key, val);
        } else if key == b"COMPATIBILITY" {
            // Not acted on by the v6 engine; kept verbatim for the round trip.
            if val == b"3" || val == b"4" || val == b"5" {
                ext_set(ctx, key, val);
            } else {
                ctx.errors.push(format_error(ERR_KEYWORD, val));
            }

            // -----------------------------------------------------------------
            // Unknown key → ext_options (R05)
            // -----------------------------------------------------------------
        } else if legacy_knows_option(key) {
            // A legacy keyword v6 does not act on (e.g. TEMPDIR): kept for the
            // round trip, and no warning, as legacy accepts it.
            ext_set(ctx, key, val);
        } else {
            ext_set(ctx, key, val);
            // Leading "\n  " as legacy writes it (report_writeLine).
            let mut w = b"\n  WARNING: Unknown option keyword '".to_vec();
            w.extend_from_slice(&tokens[0]);
            w.extend_from_slice(b"' in [OPTIONS] section - option will be ignored.");
            warn_101(ctx, w);
        }
    }

    let opt = &mut ctx.options;
    // Combine date + time parts
    if got_start_date || got_start_time {
        opt.start_date = start_date_part + start_time_part;
    }
    if got_end_date || got_end_time {
        opt.end_date = end_date_part + end_time_part;
    }
    // PARITY swmm5.c:3198-3200: legacy forms TotalDuration from the SEPARATE
    // date and time parts — floor((EndDate-StartDate)*86400 +
    // (EndTime-StartTime)*86400) × 1000. Only computable here when both dates
    // were parsed from their parts.
    if got_start_date && got_end_date {
        opt.total_duration_ms =
            ((end_date_part - start_date_part) * 86400.0 + (end_time_part - start_time_part) * 86400.0).floor() * 1000.0;
    }
    if got_rpt_date || got_rpt_time {
        opt.report_start = rpt_date_part + rpt_time_part;
    } else {
        opt.report_start = opt.start_date;
    }
}

/// `std::string::find(needle) != npos`
fn contains(hay: &[u8], needle: &[u8]) -> bool {
    hay.windows(needle.len()).any(|w| w == needle)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Recorded from glibc 2.39 (Ubuntu 24.04): strtod bits and bytes consumed, atoi, and
    // sscanf(s, "%d:%d:%lf%n", &h, &m, &sec, &n) as (return, h, sec bits, m, n).
    const CASES: &[(&str, u64, usize, i32, i32, i32, u64, i32, i32)] = &[
        ("0x1.CE4A9EEd9CEcA8p-1024", 0x0007392a7bb673b2, 24, 0, 1, 0, 0, 0, 0),
        ("-0x1.8p-1075", 0x8000000000000001, 12, 0, 1, 0, 0, 0, 0),
        ("0x1.00000000000008p0", 0x3ff0000000000000, 20, 0, 1, 0, 0, 0, 0),
        ("1e", 0x3ff0000000000000, 1, 1, 1, 1, 0, 0, 0),
        ("  +1.5e+x", 0x3ff8000000000000, 6, 1, 1, 1, 0, 0, 0),
        ("nan(123)", 0x7ff800000000007b, 8, 0, 0, 0, 0, 0, 0),
        ("-inf", 0xfff0000000000000, 4, 0, 0, 0, 0, 0, 0),
        ("0x", 0, 1, 0, 1, 0, 0, 0, 0),
        ("0xg", 0, 1, 0, 1, 0, 0, 0, 0),
        ("2.4703282292062328e-324", 1, 23, 2, 1, 2, 0, 0, 0),
        ("12:30:1.5e", 0x4028000000000000, 2, 12, 3, 12, 0x3ff8000000000000, 30, 10),
        ("1:2:x", 0x3ff0000000000000, 1, 1, 2, 1, 0, 2, 0),
        ("99999999999", 0x42374876e7ff0000, 11, 1215752191, 1, 1215752191, 0, 0, 0),
        (".e1", 0, 0, 0, 0, 0, 0, 0, 0),
        ("0x1p1024", 0x7ff0000000000000, 8, 0, 1, 0, 0, 0, 0),
    ];

    #[test]
    fn c_library_matches_glibc() {
        for &(s, bits, used, atoi, r, h, sec, m, n) in CASES {
            let (v, u) = c_strtod(s.as_bytes());
            assert_eq!((v.to_bits(), u), (bits, used), "strtod {s}");
            assert_eq!(c_atoi(s.as_bytes()), atoi, "atoi {s}");
            let (mut hh, mut mm, mut ss, mut nn) = (0, 0, 0.0f64, 0);
            let rr = sscanf_clock_n(s.as_bytes(), true, &mut hh, &mut mm, &mut ss, &mut nn);
            assert_eq!((rr, hh, ss.to_bits(), mm, nn), (r, h, sec, m, n), "sscanf {s}");
        }
    }
}
