// SPDX-License-Identifier: Apache-2.0
//
// The C library calls upstream makes on input text (std::strtod, atof, atoi), as glibc 2.39
// performs them in the C locale on the token's c_str(), so an embedded NUL ends it.
//   * strtod: leading C-locale white space, an optional sign, then INF/INFINITY,
//     NAN/NAN(n-char-sequence) (the sequence becomes the NaN payload as glibc sets it), a
//     hexadecimal float (including glibc's own rounding of some subnormal hex values), or a
//     decimal number whose exponent is taken only when complete; decimal values are
//     correctly rounded. No conversion consumes nothing.
//   * atoi: strtol base 10 saturating at the 64-bit long range, then truncated to int.
// Checked against glibc on 180,000 strings (by the agent that translated OptionsHandler).
// Written for swmm6-rust, 2 November 2026, by Claude (Anthropic) for Robert Dickinson.

pub fn c_isspace(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// `s.c_str()` as a C function sees it: the bytes before the first NUL.
pub fn c_str(s: &[u8]) -> &[u8] {
    match s.iter().position(|&c| c == 0) {
        Some(p) => &s[..p],
        None => s,
    }
}

/// glibc strtoull(s, &end, 0) over a NaN's n-char-sequence: the value when the whole
/// sequence is consumed, None otherwise. (The sequence holds only [0-9A-Za-z_].)
fn nan_payload(seq: &[u8]) -> Option<u64> {
    let (base, digits): (u32, &[u8]) = if seq.len() >= 2 && seq[0] == b'0' && (seq[1] == b'x' || seq[1] == b'X') {
        // "0x" with no hex digit after it converts the "0" and stops at the 'x'
        if seq.len() == 2 || !seq[2].is_ascii_hexdigit() {
            return None;
        }
        (16, &seq[2..])
    } else if !seq.is_empty() && seq[0] == b'0' {
        (8, seq)
    } else {
        (10, seq)
    };
    if digits.is_empty() {
        // strtoull("") converts nothing: end == start == the sequence's end
        return Some(0);
    }
    let mut v: u64 = 0;
    let mut overflow = false;
    for &c in digits {
        let d = match (c as char).to_digit(base) {
            Some(d) => d as u64,
            None => return None,
        };
        match v.checked_mul(base as u64).and_then(|x| x.checked_add(d)) {
            Some(x) => v = x,
            None => overflow = true,
        }
    }
    Some(if overflow { u64::MAX } else { v })
}

/// A hexadecimal significand (`int_digits`.`frac_digits`) times 2^`bexp`, as glibc 2.39's
/// strtod builds and rounds it (stdlib/strtod_l.c, the base-16 path and round_and_return):
/// the first 53 significant bits form the mantissa, the 54th is the round bit, any later
/// nonzero bit is sticky, and the result is rounded to nearest-even. For a subnormal result
/// glibc shifts the mantissa right and takes its new round bit from the mantissa, keeping
/// the old sticky bits but dropping the old round bit (round_and_return folds in only the
/// bits *below* it). So a hex subnormal can round differently from correct rounding, e.g.
/// 0x1.CE4A9EEd9CEcA8p-1024 gives 0x0007392a7bb673b2 where the exact value rounds to ...b3.
/// That is reproduced here; decimal input is not affected (glibc's decimal path is correctly
/// rounded, subnormals included).
fn hex_value(int_digits: &[u8], frac_digits: &[u8], bexp: i64) -> f64 {
    const MANT_DIG: i64 = 53;
    const MIN_EXP: i64 = -1021;
    const MAX_EXP: i64 = 1024;
    const DENORM_EXP: i64 = MIN_EXP - 2;
    let int_len = int_digits.len() as i64;
    let all: Vec<u32> = int_digits.iter().chain(frac_digits.iter()).map(|&c| (c as char).to_digit(16).unwrap()).collect();
    let k = match all.iter().position(|&d| d != 0) {
        None => return 0.0,
        Some(k) => k,
    };
    let d0 = all[k];
    let bits0 = 32 - d0.leading_zeros() as i64; // nbits[val]
    // exponent of the leading bit
    let mut exponent: i64 = bexp + (bits0 - 1) + 4 * (int_len - 1 - k as i64);
    // the mantissa (53 bits, leading bit at 52), the round bit, the sticky bits
    let mut m: u64 = 0;
    let mut nb: i64 = 0;
    let mut round_bit_set = false;
    let mut sticky = false;
    for (j, &d) in all[k..].iter().enumerate() {
        let width = if j == 0 { bits0 } else { 4 };
        for q in (0..width).rev() {
            let bit = (d >> q) & 1 == 1;
            if nb < MANT_DIG {
                m = (m << 1) | bit as u64;
            } else if nb == MANT_DIG {
                round_bit_set = bit;
            } else if bit {
                sticky = true;
            }
            nb += 1;
        }
        if nb > MANT_DIG + 1 && sticky {
            break;
        }
    }
    if nb < MANT_DIG {
        m <<= (MANT_DIG - nb) as u32;
    }
    // round_and_return
    let mut retval = m;
    let half_bit;
    let more_bits_nonzero;
    if exponent < MIN_EXP - 1 {
        if exponent < MIN_EXP - 1 - MANT_DIG {
            return 0.0; // underflow_value
        }
        let shift = MIN_EXP - 1 - exponent; // 1..=53
        let more_bits = sticky; // the old round bit is not carried over
        if shift == MANT_DIG {
            half_bit = (m >> 52) & 1 == 1;
            more_bits_nonzero = more_bits || (m & ((1u64 << 52) - 1)) != 0;
            retval = 0;
        } else {
            let rb = (shift - 1) as u32;
            half_bit = (m >> rb) & 1 == 1;
            more_bits_nonzero = more_bits || (m & ((1u64 << rb) - 1)) != 0;
            retval = m >> shift as u32;
        }
        exponent = DENORM_EXP;
    } else {
        half_bit = round_bit_set;
        more_bits_nonzero = sticky;
    }
    if exponent >= MAX_EXP {
        return f64::INFINITY;
    }
    // round_away, FE_TONEAREST
    if half_bit && ((retval & 1) != 0 || more_bits_nonzero) {
        retval += 1;
        if retval & (1u64 << 53) != 0 {
            exponent += 1;
            retval >>= 1;
            retval |= 1u64 << 52;
        } else if exponent == DENORM_EXP && retval & (1u64 << 52) != 0 {
            // The number was denormalized but now normalized.
            exponent = MIN_EXP - 1;
        }
    }
    if exponent >= MAX_EXP {
        return f64::INFINITY;
    }
    // MPN2FLOAT
    let biased = (exponent + 1023) as u64;
    f64::from_bits((biased << 52) | (retval & ((1u64 << 52) - 1)))
}

/// glibc `strtod(s, &end)` in the C locale: (value, bytes consumed). Nothing converted
/// gives (0.0, 0).
pub fn strtod(s: &[u8]) -> (f64, usize) {
    let s = c_str(s);
    let n = s.len();
    let at = |k: usize| -> u8 { if k < n { s[k] } else { 0 } };
    let mut i = 0;
    while i < n && c_isspace(s[i]) {
        i += 1;
    }
    let mut negative = false;
    if at(i) == b'-' {
        negative = true;
        i += 1;
    } else if at(i) == b'+' {
        i += 1;
    }
    let sign = |v: f64| if negative { -v } else { v };
    let c = at(i);
    if !(c == b'.' && at(i + 1).is_ascii_digit()) && !c.is_ascii_digit() {
        if starts_with_ci(&s[i..], b"inf") {
            let used = if starts_with_ci(&s[i + 3..], b"inity") { i + 8 } else { i + 3 };
            return (sign(f64::INFINITY), used);
        }
        if starts_with_ci(&s[i..], b"nan") {
            let mut retval = f64::NAN;
            let mut cp = i + 3;
            if at(cp) == b'(' {
                let start = cp + 1;
                let mut e = start;
                while e < n && (s[e].is_ascii_alphanumeric() || s[e] == b'_') {
                    e += 1;
                }
                if at(e) == b')' {
                    if let Some(mant) = nan_payload(&s[start..e]) {
                        retval = f64::from_bits(f64::NAN.to_bits() | (mant & 0x000F_FFFF_FFFF_FFFF & !(1u64 << 51)));
                    }
                    cp = e + 1;
                }
            }
            return (sign(retval), cp);
        }
        return (0.0, 0);
    }

    // Hexadecimal?
    let mut base16 = false;
    if c == b'0' && (at(i + 1) == b'x' || at(i + 1) == b'X') {
        base16 = true;
        i += 2;
    }
    let is_dig = |ch: u8| ch.is_ascii_digit() || (base16 && ch.is_ascii_hexdigit());
    let start_of_digits = i;
    let mut cp = i;
    while at(cp) == b'0' {
        cp += 1;
    }
    // "If no other digit but a '0' is found the result is 0.0."
    let ch = at(cp);
    let lo = ch.to_ascii_lowercase();
    let ok = is_dig(ch)
        || (ch == b'.'
            && (!base16 || cp != start_of_digits || is_dig(at(cp + 1))))
        || (base16 && cp != start_of_digits && lo == b'p')
        || (!base16 && lo == b'e');
    if !ok {
        let end = if cp == start_of_digits {
            if base16 {
                cp - 1
            } else {
                0
            }
        } else {
            cp
        };
        return (sign(0.0), end);
    }

    // Integer digits (from the start, leading zeros included), fraction, exponent.
    let int_start = start_of_digits;
    cp = start_of_digits;
    while is_dig(at(cp)) {
        cp += 1;
    }
    let int_digits = &s[int_start..cp];
    let mut frac_digits: &[u8] = &[];
    if at(cp) == b'.' {
        cp += 1;
        let fs = cp;
        while is_dig(at(cp)) {
            cp += 1;
        }
        frac_digits = &s[fs..cp];
    }
    let expp = cp;
    let mut exponent: i64 = 0;
    let lowc = at(cp).to_ascii_lowercase();
    if (base16 && lowc == b'p') || (!base16 && lowc == b'e') {
        let mut k = cp + 1;
        let mut exp_negative = false;
        if at(k) == b'-' {
            exp_negative = true;
            k += 1;
        } else if at(k) == b'+' {
            k += 1;
        }
        if at(k).is_ascii_digit() {
            while at(k).is_ascii_digit() {
                // saturate far beyond any finite result
                if exponent < 1_000_000_000 {
                    exponent = exponent * 10 + (at(k) - b'0') as i64;
                }
                k += 1;
            }
            if exp_negative {
                exponent = -exponent;
            }
            cp = k;
        } else {
            cp = expp;
        }
    }

    let v = if base16 {
        hex_value(int_digits, frac_digits, exponent)
    } else {
        let all_zero = int_digits.iter().chain(frac_digits.iter()).all(|&d| d == b'0');
        if all_zero {
            0.0
        } else {
            // Correctly rounded, as glibc's is. Normalize the exponent so the literal is
            // short: drop leading zeros and shift the decimal point.
            let mut digits: Vec<u8> = Vec::with_capacity(int_digits.len() + frac_digits.len());
            digits.extend_from_slice(int_digits);
            digits.extend_from_slice(frac_digits);
            let first = digits.iter().position(|&d| d != b'0').unwrap();
            let sig = &digits[first..];
            // value = 0.sig * 10^(int_len - first + exponent)
            let e10 = int_digits.len() as i64 - first as i64 + exponent;
            if e10 > 400 {
                f64::INFINITY
            } else if e10 < -400 {
                0.0
            } else {
                let mut lit = b"0.".to_vec();
                lit.extend_from_slice(sig);
                lit.extend_from_slice(format!("e{}", e10).as_bytes());
                std::str::from_utf8(&lit).unwrap().parse::<f64>().unwrap()
            }
        }
    };
    (sign(v), cp)
}

/// glibc `atof(s)`.
pub fn atof(s: &[u8]) -> f64 {
    strtod(s).0
}

/// glibc `atoi(s)`: `(int) strtol(s, NULL, 10)`.
pub fn atoi(s: &[u8]) -> i32 {
    let s = c_str(s);
    let n = s.len();
    let mut i = 0;
    while i < n && c_isspace(s[i]) {
        i += 1;
    }
    let mut neg = false;
    if i < n && (s[i] == b'-' || s[i] == b'+') {
        neg = s[i] == b'-';
        i += 1;
    }
    let mut v: i128 = 0;
    while i < n && s[i].is_ascii_digit() {
        if v <= i64::MAX as i128 + 1 {
            v = v * 10 + (s[i] - b'0') as i128;
        }
        i += 1;
    }
    let v = if neg { -v } else { v };
    let l: i64 = if v > i64::MAX as i128 {
        i64::MAX
    } else if v < i64::MIN as i128 {
        i64::MIN
    } else {
        v as i64
    };
    l as i32
}

pub fn starts_with_ci(s: &[u8], word: &[u8]) -> bool {
    s.len() >= word.len() && s[..word.len()].eq_ignore_ascii_case(word)
}
