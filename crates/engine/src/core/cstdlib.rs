// SPDX-License-Identifier: Apache-2.0
//
// C library number parsing that OpenSWMM Engine calls directly (std::strtod), written for
// the translation of src/engine/input/handlers/NodesHandler.cpp (HydroCouple/openswmm.engine,
// branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the
// Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// This is not upstream code but the behavior of glibc's strtod in the "C" locale with the
// default rounding mode, which upstream relies on:
//
//   * leading isspace() bytes are skipped; one optional '+' or '-';
//   * "inf"/"infinity" and "nan"/"nan(chars)" in any case;
//   * "0x"/"0X" hexadecimal significands with an optional binary 'p' exponent; a bare "0x"
//     with no hex digit converts the "0" and stops at the 'x';
//   * otherwise decimal digits with at most one '.', at least one digit, and an exponent
//     only when it is complete ("1e" and "1e+" stop before the 'e');
//   * correctly rounded; overflow gives infinity and underflow gives zero or a subnormal
//     (errno is not modelled; no caller reads it).
//
// Decimal values are rounded by Rust's own parser, which is also correctly rounded, so
// equal prefixes give equal bits.

/// C `isspace` in the "C" locale.
fn is_c_space(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

fn starts_with_ci(s: &[u8], word: &[u8]) -> bool {
    s.len() >= word.len() && s[..word.len()].eq_ignore_ascii_case(word)
}

fn hex_val(c: u8) -> Option<u64> {
    match c {
        b'0'..=b'9' => Some((c - b'0') as u64),
        b'a'..=b'f' => Some((c - b'a' + 10) as u64),
        b'A'..=b'F' => Some((c - b'A' + 10) as u64),
        _ => None,
    }
}

/// 2^k as a double, for -1022 <= k <= 1023.
fn pow2(k: i32) -> f64 {
    f64::from_bits(((k + 1023) as u64) << 52)
}

/// m * 2^e for an m that fits in 53 bits, where the product is exact or overflows.
fn scale2(m: f64, mut e: i64) -> f64 {
    let mut v = m;
    while e > 1000 {
        v *= pow2(1000);
        e -= 1000;
        if v.is_infinite() {
            return v;
        }
    }
    while e < -1000 {
        v *= pow2(-1000);
        e += 1000;
    }
    v * pow2(e as i32)
}

/// Parse the hexadecimal significand starting at `s[i]` (just past "0x"); returns the
/// magnitude and the end index, or None if no hex digit follows.
fn parse_hex(s: &[u8], mut i: usize) -> Option<(f64, usize)> {
    let n = s.len();
    let mut mant: u64 = 0;
    let mut sticky = false;
    let mut exp2: i64 = 0;
    let mut any = false;
    while i < n {
        match hex_val(s[i]) {
            Some(d) => {
                any = true;
                if mant >> 60 == 0 {
                    mant = (mant << 4) | d;
                } else {
                    sticky |= d != 0;
                    exp2 += 4;
                }
                i += 1;
            }
            None => break,
        }
    }
    if i < n && s[i] == b'.' {
        let mut j = i + 1;
        let mut frac_any = false;
        while j < n {
            match hex_val(s[j]) {
                Some(d) => {
                    frac_any = true;
                    if mant >> 60 == 0 {
                        mant = (mant << 4) | d;
                        exp2 -= 4;
                    } else {
                        sticky |= d != 0;
                    }
                    j += 1;
                }
                None => break,
            }
        }
        if any || frac_any {
            any = true;
            i = j;
        }
    }
    if !any {
        return None;
    }
    // binary exponent: only when complete
    if i < n && (s[i] == b'p' || s[i] == b'P') {
        let mut j = i + 1;
        let mut eneg = false;
        if j < n && (s[j] == b'+' || s[j] == b'-') {
            eneg = s[j] == b'-';
            j += 1;
        }
        if j < n && s[j].is_ascii_digit() {
            let mut e: i64 = 0;
            while j < n && s[j].is_ascii_digit() {
                if e < 1_000_000_000 {
                    e = e * 10 + (s[j] - b'0') as i64;
                }
                j += 1;
            }
            exp2 += if eneg { -e } else { e };
            i = j;
        }
    }
    if mant == 0 {
        return Some((0.0, i));
    }
    // Round mant * 2^exp2 to 53 bits, or to the subnormal grid (lsb 2^-1074).
    let bits = 64 - mant.leading_zeros() as i64;
    let lsb = std::cmp::max(exp2 + bits - 53, -1074);
    let shift = lsb - exp2;
    // Round half to even on the dropped bits, with `sticky` standing for nonzero hex
    // digits beyond the 64 bits kept.
    let m = if shift > 64 {
        0 // below half the smallest step
    } else if shift == 64 {
        let above_half = (mant >> 63) == 1 && (mant & !(1u64 << 63) != 0 || sticky);
        if above_half {
            1
        } else {
            0
        }
    } else if shift > 0 {
        let kept = mant >> shift;
        let dropped = mant & ((1u64 << shift) - 1);
        let half = 1u64 << (shift - 1);
        let up = dropped > half || (dropped == half && (sticky || kept & 1 == 1));
        if up {
            kept + 1
        } else {
            kept
        }
    } else {
        mant << (-shift)
    };
    Some((scale2(m as f64, lsb), i))
}

/// `std::strtod` on the bytes of `s` (the caller passes the C string, i.e. stops at the
/// first NUL). Returns the value and the number of bytes consumed; 0 consumed means no
/// conversion was performed and the value is 0.0.
pub fn strtod(s: &[u8]) -> (f64, usize) {
    let n = s.len();
    let mut i = 0;
    while i < n && is_c_space(s[i]) {
        i += 1;
    }
    let mut neg = false;
    if i < n && (s[i] == b'+' || s[i] == b'-') {
        neg = s[i] == b'-';
        i += 1;
    }
    let sign = |v: f64| if neg { -v } else { v };
    let rest = &s[i..];
    if starts_with_ci(rest, b"inf") {
        let used = if starts_with_ci(rest, b"infinity") { 8 } else { 3 };
        return (sign(f64::INFINITY), i + used);
    }
    if starts_with_ci(rest, b"nan") {
        let mut j = i + 3;
        if j < n && s[j] == b'(' {
            let mut k = j + 1;
            while k < n && (s[k].is_ascii_alphanumeric() || s[k] == b'_') {
                k += 1;
            }
            if k < n && s[k] == b')' {
                j = k + 1;
            }
        }
        return (sign(f64::NAN), j);
    }
    if rest.len() >= 2 && rest[0] == b'0' && (rest[1] == b'x' || rest[1] == b'X') {
        return match parse_hex(s, i + 2) {
            Some((v, end)) => (sign(v), end),
            None => (sign(0.0), i + 1),
        };
    }
    // decimal
    let start = i;
    let mut digits = 0usize;
    while i < n && s[i].is_ascii_digit() {
        i += 1;
        digits += 1;
    }
    if i < n && s[i] == b'.' {
        let mut j = i + 1;
        let mut fd = 0usize;
        while j < n && s[j].is_ascii_digit() {
            j += 1;
            fd += 1;
        }
        if digits + fd > 0 {
            i = j;
            digits += fd;
        }
    }
    if digits == 0 {
        return (0.0, 0);
    }
    if i < n && (s[i] == b'e' || s[i] == b'E') {
        let mut j = i + 1;
        if j < n && (s[j] == b'+' || s[j] == b'-') {
            j += 1;
        }
        if j < n && s[j].is_ascii_digit() {
            while j < n && s[j].is_ascii_digit() {
                j += 1;
            }
            i = j;
        }
    }
    // The matched text is ASCII digits, '.', 'e', sign: valid UTF-8 and valid for Rust's
    // parser, which rounds correctly.
    let text = std::str::from_utf8(&s[start..i]).unwrap_or("0");
    let v: f64 = text.parse().unwrap_or(0.0);
    (sign(v), i)
}

#[cfg(test)]
mod tests {
    use super::strtod;
    #[test]
    fn cases() {
        assert_eq!(strtod(b"1.5"), (1.5, 3));
        assert_eq!(strtod(b"  +2"), (2.0, 4));
        assert_eq!(strtod(b"1e"), (1.0, 1));
        assert_eq!(strtod(b"1e+"), (1.0, 1));
        assert_eq!(strtod(b"abc"), (0.0, 0));
        assert_eq!(strtod(b"0x"), (0.0, 1));
        assert_eq!(strtod(b"0x10"), (16.0, 4));
        assert_eq!(strtod(b"0x1.8p1"), (3.0, 7));
        assert_eq!(strtod(b"-0x.8"), (-0.5, 5));
        assert_eq!(strtod(b"0x1p-1074").0.to_bits(), 1);
        assert_eq!(strtod(b"0x1p-1075").0.to_bits(), 0);
        assert_eq!(strtod(b"0x3p-1076").0.to_bits(), 1);
        assert_eq!(strtod(b"0x1p1024").0, f64::INFINITY);
        assert_eq!(strtod(b"0x1.fffffffffffff8p0").0, 2.0);
        assert_eq!(strtod(b"0x1.fffffffffffff7p0").0, f64::from_bits(0x3FFF_FFFF_FFFF_FFFF));
        assert!(strtod(b"inFinity").0.is_infinite() && strtod(b"inFinity").1 == 8);
        assert_eq!(strtod(b"nan(abc)x").1, 8);
        assert_eq!(strtod(b"2e-324").0, 0.0);
        assert_eq!(strtod(b".e1"), (0.0, 0));
    }
}
