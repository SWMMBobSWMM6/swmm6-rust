// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/Tokenizer.cpp and Tokenizer.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The input is bytes, as in the C++ (a std::string_view of char): a name may hold any byte
// but the separators, Latin-1 included, and is never re-encoded. White space is a space or
// a tab only; std::isspace (C locale) is used only by trim(). Upper-casing is ASCII only,
// as std::toupper is in the C locale.

use crate::core::charconv::from_chars_double;

#[inline]
fn is_space(c: u8) -> bool {
    c == b' ' || c == b'\t'
}

/// C-locale `std::isspace`.
#[inline]
fn is_c_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// The line up to its first `;` that is not inside double quotes.
pub fn strip_comment(line: &[u8]) -> &[u8] {
    let mut in_quote = false;
    for (i, &c) in line.iter().enumerate() {
        if c == b'"' {
            in_quote = !in_quote;
        } else if c == b';' && !in_quote {
            return &line[..i];
        }
    }
    line
}

fn is_pure_csv_line(line: &[u8]) -> bool {
    let mut in_quotes = false;
    for &c in line {
        if c == b'"' {
            in_quotes = !in_quotes;
            continue;
        }
        if !in_quotes && is_space(c) {
            return false;
        }
    }
    true
}

/// Whether the comma at `pos` separates fields. On a line with no white space outside
/// quotes every comma does; otherwise only one that opens or closes the line or touches
/// white space. A comma between two other characters is part of a name, as in legacy SWMM.
fn comma_separates(line: &[u8], pos: usize, pure_csv: bool) -> bool {
    if pure_csv {
        return true;
    }
    if pos == 0 || pos + 1 >= line.len() {
        return true;
    }
    is_space(line[pos - 1]) || is_space(line[pos + 1])
}

/// Split one line into tokens (the C++ `Tokenizer::tokenize`): comment removed, a quoted
/// token taken whole without its quotes, empty CSV fields kept.
pub fn tokenize(line: &[u8]) -> Vec<Vec<u8>> {
    let stripped = strip_comment(line);
    let pure_csv = is_pure_csv_line(stripped);
    let mut tokens: Vec<Vec<u8>> = Vec::with_capacity(8);
    let n = stripped.len();
    let mut i = 0;
    while i < n {
        while i < n && (stripped[i] == b' ' || stripped[i] == b'\t') {
            i += 1;
        }
        if i >= n {
            break;
        }
        if stripped[i] == b',' && comma_separates(stripped, i, pure_csv) {
            // An empty field only where this comma follows another or opens the line.
            if i == 0 || stripped[i - 1] == b',' {
                tokens.push(Vec::new());
            }
            i += 1;
            continue;
        }
        if stripped[i] == b'"' {
            i += 1;
            let mut tok = Vec::new();
            while i < n && stripped[i] != b'"' {
                tok.push(stripped[i]);
                i += 1;
            }
            if i < n {
                i += 1;
            }
            tokens.push(tok);
            while i < n && (stripped[i] == b' ' || stripped[i] == b'\t') {
                i += 1;
            }
            if i < n && stripped[i] == b',' {
                i += 1;
            }
            continue;
        }
        // A regular token; embedded quotes group and stay in it.
        let start = i;
        let mut in_quote = false;
        while i < n {
            let c = stripped[i];
            if c == b'"' {
                in_quote = !in_quote;
            } else if !in_quote && (is_space(c) || (c == b',' && comma_separates(stripped, i, pure_csv))) {
                break;
            }
            i += 1;
        }
        tokens.push(stripped[start..i].to_vec());
        if i < n && stripped[i] == b',' {
            i += 1;
        }
    }
    tokens
}

/// The C++ `Tokenizer::tokenize_views_into`: no quote handling and no empty fields.
pub fn tokenize_views_into<'a>(line: &'a [u8], tokens: &mut Vec<&'a [u8]>) {
    let stripped = strip_comment(line);
    let pure_csv = is_pure_csv_line(stripped);
    tokens.clear();
    let n = stripped.len();
    let mut i = 0;
    while i < n {
        while i < n && (stripped[i] == b' ' || stripped[i] == b'\t') {
            i += 1;
        }
        if i >= n {
            break;
        }
        if stripped[i] == b',' && comma_separates(stripped, i, pure_csv) {
            i += 1;
            continue;
        }
        let start = i;
        while i < n && !is_space(stripped[i]) && !(stripped[i] == b',' && comma_separates(stripped, i, pure_csv)) {
            i += 1;
        }
        if i > start {
            tokens.push(&stripped[start..i]);
        }
        if i < n && stripped[i] == b',' {
            i += 1;
        }
    }
}

pub fn tokenize_views(line: &[u8]) -> Vec<&[u8]> {
    let mut v = Vec::new();
    tokenize_views_into(line, &mut v);
    v
}

pub fn to_upper_inplace(s: &mut [u8]) {
    s.make_ascii_uppercase();
}

pub fn to_upper(s: &[u8]) -> Vec<u8> {
    s.to_ascii_uppercase()
}

pub fn trim(s: &[u8]) -> &[u8] {
    let mut lo = 0;
    while lo < s.len() && is_c_space(s[lo]) {
        lo += 1;
    }
    let mut hi = s.len();
    while hi > lo && is_c_space(s[hi - 1]) {
        hi -= 1;
    }
    &s[lo..hi]
}

/// Whether the whole token (after trim) is one number as std::from_chars reads it.
pub fn is_numeric(sv: &[u8]) -> bool {
    let sv = trim(sv);
    if sv.is_empty() {
        return false;
    }
    let (r, used) = from_chars_double(sv);
    r.is_ok() && used == sv.len()
}

fn upper_short(sv: &[u8]) -> Option<Vec<u8>> {
    // The C++ copies into char buf[8] and gives up on 8 bytes or more.
    if sv.len() >= 8 {
        None
    } else {
        Some(sv.to_ascii_uppercase())
    }
}

pub fn is_boolean(sv: &[u8]) -> bool {
    match upper_short(sv) {
        Some(u) => matches!(u.as_slice(), b"YES" | b"NO" | b"TRUE" | b"FALSE" | b"1" | b"0"),
        None => false,
    }
}

pub fn parse_boolean(sv: &[u8]) -> bool {
    match upper_short(sv) {
        Some(u) => matches!(u.as_slice(), b"YES" | b"TRUE" | b"1"),
        None => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn t(s: &str) -> Vec<String> {
        tokenize(s.as_bytes()).into_iter().map(|v| String::from_utf8_lossy(&v).into_owned()).collect()
    }
    #[test]
    fn commas_and_quotes() {
        assert_eq!(t("J1  100  5 ; comment"), ["J1", "100", "5"]);
        assert_eq!(t("a,,b"), ["a", "", "b"]);
        assert_eq!(t("A ,B"), ["A", "B"]);
        assert_eq!(t("N1,2 3"), ["N1,2", "3"]);
        assert_eq!(t("\"My Node\" 4"), ["My Node", "4"]);
        assert_eq!(t("config=\"a b/c.rxn\" x"), ["config=\"a b/c.rxn\"", "x"]);
        assert_eq!(t("x \";\" y"), ["x", ";", "y"]);
    }
}
