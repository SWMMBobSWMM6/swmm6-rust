// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/data/Report2DVars.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The REPORT_2D_VARIABLES / REPORT_2D_SPECIES / REPORT_2D_STEP grammar shared by
// [2D_OPTIONS], the C API and the writer (namespace report2d upstream).

#![allow(clippy::all)]

use crate::core::cprintf::{format, Arg};
use crate::core::options::Bytes;
use crate::twod::data::solver_options_2d::report2d::*;

/// tokens(): the group keywords in bit order.
pub const TOKENS: [&[u8]; 16] = [
    b"DEPTH",
    b"VELOCITY",
    b"EDGE_FLUX",
    b"NODE_HEAD",
    b"SPECIES",
    b"RAINFALL",
    b"INFILTRATION",
    b"COUPLING",
    b"GRADIENTS",
    b"CONTINUITY",
    b"ENVELOPES",
    b"BUILDUP",
    b"GROUNDWATER",
    b"GW_DETAILED",
    b"HEAT",
    b"TURBULENCE",
];

/// iequalsTok(): byte case-insensitive equality with a keyword.
pub fn iequals_tok(a: &[u8], b: &[u8]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_ascii_uppercase() == y.to_ascii_uppercase())
}

/// splitList(): split on ',', ' ', '\t', '\r', '\n', dropping empty pieces.
pub fn split_list(text: &[u8]) -> Vec<Bytes> {
    let mut toks: Vec<Bytes> = Vec::new();
    let mut cur: Bytes = Bytes::new();
    for &c in text {
        if c == b',' || c == b' ' || c == b'\t' || c == b'\r' || c == b'\n' {
            if !cur.is_empty() {
                toks.push(std::mem::take(&mut cur));
            }
        } else {
            cur.push(c);
        }
    }
    if !cur.is_empty() {
        toks.push(cur);
    }
    toks
}

/// parseMask(): "" on success (mask set, DEPTH always on), otherwise the diagnostic. `mask` is
/// zeroed first and keeps whatever groups preceded an unknown token.
pub fn parse_mask(text: &[u8], mask: &mut u32) -> Bytes {
    *mask = 0;
    let toks = split_list(text);
    if toks.is_empty() {
        return b"REPORT_2D_VARIABLES needs DEFAULT|MINIMAL|ALL or a token list".to_vec();
    }
    for t in &toks {
        if iequals_tok(t, b"DEFAULT") {
            *mask |= DEFAULT_MASK;
            continue;
        }
        if iequals_tok(t, b"MINIMAL") {
            *mask |= MINIMAL_MASK;
            continue;
        }
        if iequals_tok(t, b"ALL") {
            *mask |= ALL_MASK;
            continue;
        }
        if iequals_tok(t, b"NONE") {
            continue;
        }
        let mut found = false;
        for (i, name) in TOKENS.iter().enumerate() {
            if iequals_tok(t, name) {
                *mask |= 1u32 << i;
                found = true;
                break;
            }
        }
        if !found {
            let mut all: Bytes = Bytes::new();
            for n in TOKENS.iter() {
                if !all.is_empty() {
                    all.push(b'|');
                }
                all.extend_from_slice(n);
            }
            let mut m = b"Unknown REPORT_2D_VARIABLES token '".to_vec();
            m.extend_from_slice(t);
            m.extend_from_slice(b"' (expected DEFAULT|MINIMAL|ALL|NONE or ");
            m.extend_from_slice(&all);
            m.push(b')');
            return m;
        }
    }
    *mask |= DEPTH;
    Bytes::new()
}

/// formatMask()
pub fn format_mask(mut mask: u32) -> Bytes {
    mask &= ALL_MASK;
    if mask == DEFAULT_MASK {
        return b"DEFAULT".to_vec();
    }
    if mask == MINIMAL_MASK {
        return b"MINIMAL".to_vec();
    }
    if mask == ALL_MASK {
        return b"ALL".to_vec();
    }
    let mut out: Bytes = Bytes::new();
    for (i, n) in TOKENS.iter().enumerate() {
        if mask & (1u32 << i) != 0 {
            if !out.is_empty() {
                out.push(b' ');
            }
            out.extend_from_slice(n);
        }
    }
    if out.is_empty() {
        b"DEPTH".to_vec()
    } else {
        out
    }
}

/// formatSpecies()
pub fn format_species(species: &[Bytes]) -> Bytes {
    if species.is_empty() {
        return b"ALL".to_vec();
    }
    let mut out: Bytes = Bytes::new();
    for n in species {
        if !out.is_empty() {
            out.push(b' ');
        }
        out.extend_from_slice(n);
    }
    out
}

/// parseSpecies(): the species list; a lone ALL is the empty list.
pub fn parse_species(text: &[u8]) -> Vec<Bytes> {
    let mut v = split_list(text);
    if v.len() == 1 && iequals_tok(&v[0], b"ALL") {
        v.clear();
    }
    v
}

/// formatStep(): HH:MM:SS for whole seconds, %.17g otherwise, 00:00:00 for <= 0 or NaN.
pub fn format_step(seconds: f64) -> Bytes {
    if !(seconds > 0.0) {
        return b"00:00:00".to_vec();
    }
    if seconds.floor() != seconds || seconds >= i64::MAX as f64 {
        return format(b"%.17g", &[Arg::F(seconds)]);
    }
    let total: i64 = if seconds > 0.0 { seconds as i64 } else { 0 };
    format(b"%02ld:%02ld:%02ld", &[Arg::I(total / 3600), Arg::I((total % 3600) / 60), Arg::I(total % 60)])
}
