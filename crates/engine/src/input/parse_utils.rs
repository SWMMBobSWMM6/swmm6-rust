// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/InputParseUtils.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Every number in a section row is read here, through std::from_chars semantics
// (core::charconv): a prefix is accepted where the C++ accepts one, and a failed parse
// leaves the default in place, as `double v = def; from_chars(..., v)` does.

use crate::core::charconv::{from_chars_double, from_chars_i32, from_chars_u32};
use crate::core::datetime;
use crate::core::error_codes::{format_error, format_error_detail, ERR_DUP_NAME};
use crate::core::options::Bytes;
use crate::data::name_index::NameIndex;

/// The leading number of `sv`, or `def` if there is none (or it is out of range).
pub fn to_double(sv: &[u8], def: f64) -> f64 {
    match from_chars_double(sv).0 {
        Ok(v) => v,
        Err(_) => def,
    }
}

/// The whole of `sv` as one number. Empty is 0 (and true); a single leading '+' is
/// allowed, as strtod allows it; anything left over is false and leaves `out` alone.
pub fn parse_double_strict(sv: &[u8], out: &mut f64) -> bool {
    let mut sv = sv;
    if sv.is_empty() {
        *out = 0.0;
        return true;
    }
    if sv[0] == b'+' {
        sv = &sv[1..];
        if sv.is_empty() || sv[0] == b'+' || sv[0] == b'-' {
            return false;
        }
    }
    let (r, used) = from_chars_double(sv);
    match r {
        Ok(v) if used == sv.len() => {
            *out = v;
            true
        }
        _ => false,
    }
}

/// The leading integer of `sv`, or `def`.
pub fn to_int(sv: &[u8], def: i32) -> i32 {
    match from_chars_i32(sv).0 {
        Ok(v) => v,
        Err(_) => def,
    }
}

/// A reader over a byte slice with the C++ lambdas' `read_uint`.
struct Cursor<'a> {
    s: &'a [u8],
    p: usize,
}

impl<'a> Cursor<'a> {
    fn read_uint(&mut self, out: &mut u32) -> bool {
        let (r, used) = from_chars_u32(&self.s[self.p..]);
        match r {
            Ok(v) => {
                *out = v;
                self.p += used;
                true
            }
            Err(_) => false,
        }
    }
    fn at(&self, c: u8) -> bool {
        self.p < self.s.len() && self.s[self.p] == c
    }
}

/// M/D/Y with '/' or '-' between, as a date; 0.0 if it does not read.
pub fn parse_date(sv: &[u8]) -> f64 {
    let (mut m, mut d, mut y) = (0u32, 0u32, 0u32);
    let mut c = Cursor { s: sv, p: 0 };
    let read_sep = |c: &mut Cursor| -> bool {
        if c.at(b'/') || c.at(b'-') {
            c.p += 1;
            true
        } else {
            false
        }
    };
    if !c.read_uint(&mut m) {
        return 0.0;
    }
    if !read_sep(&mut c) {
        return 0.0;
    }
    if !c.read_uint(&mut d) {
        return 0.0;
    }
    if !read_sep(&mut c) {
        return 0.0;
    }
    if !c.read_uint(&mut y) {
        return 0.0;
    }
    datetime::encode_date(y as i32, m as i32, d as i32)
}

/// Seconds from a plain number, or from H, H:M or H:M:S (S may be fractional).
pub fn parse_time_seconds(sv: &[u8]) -> f64 {
    let (r, used) = from_chars_double(sv);
    if let Ok(v) = r {
        if used == sv.len() {
            return v;
        }
    }
    let (mut h, mut m) = (0u32, 0u32);
    let mut s = 0.0f64;
    let mut c = Cursor { s: sv, p: 0 };
    if !c.read_uint(&mut h) {
        return 0.0;
    }
    if c.at(b':') {
        c.p += 1;
        if !c.read_uint(&mut m) {
            return h as f64;
        }
        if c.at(b':') {
            c.p += 1;
            if let Ok(v) = from_chars_double(&sv[c.p..]).0 {
                s = v;
            }
        }
    }
    h as f64 * 3600.0 + m as f64 * 60.0 + s
}

/// A time of day as a fraction of a day: decimal hours, or H[:M[:S]].
pub fn parse_time_day_fraction(sv: &[u8]) -> f64 {
    let (r, used) = from_chars_double(sv);
    if let Ok(dec_hours) = r {
        if used == sv.len() {
            return dec_hours / 24.0;
        }
    }
    let (mut h, mut m, mut s) = (0u32, 0u32, 0u32);
    let mut c = Cursor { s: sv, p: 0 };
    if !c.read_uint(&mut h) {
        return 0.0;
    }
    if c.at(b':') {
        c.p += 1;
    }
    c.read_uint(&mut m);
    if c.at(b':') {
        c.p += 1;
    }
    c.read_uint(&mut s);
    datetime::encode_time(h as i32, m as i32, s as i32)
}

pub fn parse_datetime(date_sv: &[u8], time_sv: &[u8]) -> f64 {
    parse_date(date_sv) + parse_time_day_fraction(time_sv)
}

/// Add a name to an index; on a duplicate, record ERROR 207 (naming the existing spelling
/// when it differs only in case) and return -1.
pub fn add_unique(reg: &mut NameIndex, name: &[u8], errors: &mut Vec<Bytes>) -> i32 {
    let idx = reg.try_add(name);
    if idx < 0 {
        match reg.canonical(name) {
            Some(canon) if canon != name => {
                let mut detail = b"(matches existing ID '".to_vec();
                detail.extend_from_slice(canon);
                detail.extend_from_slice(b"'; names are case-insensitive)");
                errors.push(format_error_detail(ERR_DUP_NAME, name, &detail));
            }
            _ => errors.push(format_error(ERR_DUP_NAME, name)),
        }
    }
    idx
}
