// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/charconv_compat.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Upstream calls std::from_chars(first, last, double&) (chars_format::general). On the
// reference build (Linux, GCC 13 libstdc++) that is the path taken, so this module
// reproduces libstdc++'s observable behavior, measured on that build:
//
//   * an optional '-', never '+', never leading white space;
//   * "inf", "infinity" and "nan" in any case, "nan(...)" with its payload;
//   * otherwise digits with at most one '.', at least one digit, then an exponent only if
//     it is complete ("1e" and "1e+" stop before the 'e');
//   * the longest such prefix is used; nothing after it is examined;
//   * a finite literal that rounds to infinity, or a nonzero one that rounds to zero, is
//     out of range: the value is not stored (2e-324 is out of range, 2.5e-324 gives the
//     smallest subnormal, 0e999999 is zero).
//
// Rounding is correct rounding in both libraries, so equal prefixes give equal bits.

/// The `std::errc` outcomes `from_chars` can report.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Errc {
    InvalidArgument,
    ResultOutOfRange,
}

/// `std::from_chars` for `double`: the value, or why not, and how many bytes were used
/// (0 on `InvalidArgument`; the matched length on `ResultOutOfRange`, as libstdc++ reports).
pub fn from_chars_double(s: &[u8]) -> (Result<f64, Errc>, usize) {
    let n = s.len();
    let mut i = 0;
    let neg = i < n && s[i] == b'-';
    if neg {
        i += 1;
    }
    // inf / infinity / nan
    let rest = &s[i..];
    let lower = |k: usize| -> Vec<u8> { rest.iter().take(k).map(|c| c.to_ascii_lowercase()).collect() };
    if lower(3) == b"inf" {
        let used = if lower(8) == b"infinity" { i + 8 } else { i + 3 };
        return (Ok(if neg { f64::NEG_INFINITY } else { f64::INFINITY }), used);
    }
    if lower(3) == b"nan" {
        let mut used = i + 3;
        if used < n && s[used] == b'(' {
            let mut j = used + 1;
            while j < n && (s[j].is_ascii_alphanumeric() || s[j] == b'_') {
                j += 1;
            }
            if j < n && s[j] == b')' {
                used = j + 1;
            }
        }
        let v = if neg { -f64::NAN } else { f64::NAN };
        return (Ok(v), used);
    }
    // decimal significand
    let start = i;
    let mut digits = 0usize;
    let mut nonzero = false;
    while i < n && s[i].is_ascii_digit() {
        nonzero |= s[i] != b'0';
        digits += 1;
        i += 1;
    }
    if i < n && s[i] == b'.' {
        i += 1;
        while i < n && s[i].is_ascii_digit() {
            nonzero |= s[i] != b'0';
            digits += 1;
            i += 1;
        }
    }
    if digits == 0 {
        return (Err(Errc::InvalidArgument), 0);
    }
    // exponent, only if complete
    if i < n && (s[i] == b'e' || s[i] == b'E') {
        let mut j = i + 1;
        if j < n && (s[j] == b'+' || s[j] == b'-') {
            j += 1;
        }
        let k = j;
        while j < n && s[j].is_ascii_digit() {
            j += 1;
        }
        if j > k {
            i = j;
        }
    }
    let text = std::str::from_utf8(&s[start..i]).expect("ASCII digits");
    let mag: f64 = text.parse().expect("validated decimal literal");
    if mag.is_infinite() || (mag == 0.0 && nonzero) {
        return (Err(Errc::ResultOutOfRange), i);
    }
    (Ok(if neg { -mag } else { mag }), i)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Every case below was run through std::from_chars on the reference build
    // (g++ 13.3.0, libstdc++) and the result recorded here.
    #[test]
    fn matches_libstdcxx() {
        let ok = |s: &str, used: usize, v: f64| {
            let (r, u) = from_chars_double(s.as_bytes());
            assert_eq!(u, used, "{s}");
            let got = r.unwrap_or_else(|e| panic!("{s}: {e:?}"));
            if v.is_nan() {
                assert!(got.is_nan(), "{s}");
            } else {
                assert_eq!(got.to_bits(), v.to_bits(), "{s}");
            }
        };
        let err = |s: &str, used: usize, e: Errc| {
            assert_eq!(from_chars_double(s.as_bytes()), (Err(e), used), "{s}");
        };
        ok("1", 1, 1.0);
        ok("-1", 2, -1.0);
        err("+1", 0, Errc::InvalidArgument);
        err(" 1", 0, Errc::InvalidArgument);
        err("1e400", 5, Errc::ResultOutOfRange);
        err("-1e400", 6, Errc::ResultOutOfRange);
        err("1e-400", 6, Errc::ResultOutOfRange);
        ok("4.9e-324", 8, 4.9406564584124654e-324);
        err("2e-324", 6, Errc::ResultOutOfRange);
        ok("2.5e-324", 8, 4.9406564584124654e-324);
        ok("1e-320", 6, 9.9998886718268301e-321);
        ok("inf", 3, f64::INFINITY);
        ok("INF", 3, f64::INFINITY);
        ok("Infinity", 8, f64::INFINITY);
        ok("infin", 3, f64::INFINITY);
        ok("nan", 3, f64::NAN);
        ok("NaN(abc)", 8, f64::NAN);
        ok("nan(", 3, f64::NAN);
        ok("1e", 1, 1.0);
        ok("1e+", 1, 1.0);
        ok("1.", 2, 1.0);
        err(".", 0, Errc::InvalidArgument);
        ok(".5", 2, 0.5);
        ok("-.5e2x", 5, -50.0);
        ok("0x10", 1, 0.0);
        ok("1,5", 1, 1.0);
        ok("1d5", 1, 1.0);
        ok("0e999999", 8, 0.0);
        ok("00012", 5, 12.0);
        ok("1e-310", 6, 9.9999999999999694e-311);
    }
}
