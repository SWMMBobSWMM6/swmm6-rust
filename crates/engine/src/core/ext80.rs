// SPDX-License-Identifier: Apache-2.0
//
// Written for swmm6-rust, 6 November 2026, by Claude (Anthropic) at the request of Robert
// Dickinson. Not an OpenSWMM release.
//
// The x87 80-bit extended format (`long double` on x86-64 Linux, 64-bit significand, round
// to nearest even), in software. Safe Rust has no long double; upstream accumulates a dot
// product in one (hydraulics/fv/NativeGmres.cpp `dot`: `long double sum = 0.0; sum +=
// static_cast<long double>(a[i]) * b[i];`), so the translation reproduces each extended
// operation exactly: the product of two doubles rounded to 64 significand bits, the sum
// rounded to 64 bits, and the final conversion to double rounded to 53 bits (or fewer, in
// the subnormal range). Checked against g++'s long double by gate/diff/ext80dot.
//
// Only finite values are handled (the callers pass finite vectors); the exponent range of
// the format is far wider than any product of two doubles, so no extended overflow or
// underflow can occur.

/// A finite extended value: (-1)^neg * sig * 2^(exp - 63), sig normalized (bit 63 set) or 0.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Ext80 {
    pub neg: bool,
    pub exp: i32,
    pub sig: u64,
}

impl Ext80 {
    pub const ZERO: Ext80 = Ext80 { neg: false, exp: 0, sig: 0 };

    /// (long double)x: exact.
    pub fn from_f64(x: f64) -> Ext80 {
        let bits = x.to_bits();
        let neg = (bits >> 63) != 0;
        let e = ((bits >> 52) & 0x7ff) as i32;
        let frac = bits & ((1u64 << 52) - 1);
        if e == 0 {
            if frac == 0 {
                return Ext80 { neg, exp: 0, sig: 0 };
            }
            let lz = frac.leading_zeros() as i32;
            // frac * 2^-1074 = (frac << lz) * 2^(-1074 - lz), bit 63 set after the shift
            return Ext80 { neg, exp: -1011 - lz, sig: frac << lz };
        }
        Ext80 { neg, exp: e - 1023, sig: ((1u64 << 52) | frac) << 11 }
    }

    pub fn is_zero(&self) -> bool {
        self.sig == 0
    }

    /// v * 2^e0 (v != 0) rounded to a 64-bit significand, nearest even. The lowest bit of v may
    /// be a jammed sticky bit; v has at least two bits below the rounding point whenever bits
    /// are dropped, so the rounding is exact.
    fn round_u128(neg: bool, v: u128, e0: i32) -> Ext80 {
        let p = 127 - v.leading_zeros() as i32; // msb position
        if p <= 63 {
            let sig = (v << (63 - p)) as u64;
            return Ext80 { neg, exp: e0 + p, sig };
        }
        let s = (p - 63) as u32;
        let mut sig = (v >> s) as u64;
        let rem = v & ((1u128 << s) - 1);
        let half = 1u128 << (s - 1);
        let mut exp = e0 + p;
        if rem > half || (rem == half && (sig & 1) == 1) {
            let (n, carry) = sig.overflowing_add(1);
            if carry {
                sig = 1u64 << 63;
                exp += 1;
            } else {
                sig = n;
            }
        }
        Ext80 { neg, exp, sig }
    }

    /// (long double)a * (long double)b
    pub fn mul_f64(a: f64, b: f64) -> Ext80 {
        let x = Ext80::from_f64(a);
        let y = Ext80::from_f64(b);
        let neg = x.neg != y.neg;
        if x.is_zero() || y.is_zero() {
            return Ext80 { neg, exp: 0, sig: 0 };
        }
        let p = (x.sig as u128) * (y.sig as u128);
        // value = p * 2^(x.exp - 63 + y.exp - 63)
        Ext80::round_u128(neg, p, x.exp + y.exp - 126)
    }

    /// self + other, rounded to nearest even (an exact zero sum is +0).
    pub fn add(self, other: Ext80) -> Ext80 {
        if self.is_zero() && other.is_zero() {
            return Ext80 { neg: self.neg && other.neg, exp: 0, sig: 0 };
        }
        if self.is_zero() {
            return other;
        }
        if other.is_zero() {
            return self;
        }
        // a: the larger magnitude
        let (a, b) = if (self.exp, self.sig) >= (other.exp, other.sig) { (self, other) } else { (other, self) };
        let d = (a.exp - b.exp) as u32;
        let big = (a.sig as u128) << 63; // value big * 2^(a.exp - 126)
        let small_full = (b.sig as u128) << 63;
        let small = if d == 0 {
            small_full
        } else if d >= 127 {
            1 // only the sticky bit survives
        } else {
            let kept = small_full >> d;
            let lost = small_full & ((1u128 << d) - 1);
            kept | (lost != 0) as u128
        };
        let e0 = a.exp - 126;
        if a.neg == b.neg {
            Ext80::round_u128(a.neg, big + small, e0)
        } else {
            let r = big - small;
            if r == 0 {
                return Ext80::ZERO;
            }
            Ext80::round_u128(a.neg, r, e0)
        }
    }

    /// (double)self, rounded to nearest even (finite, in range).
    pub fn to_f64(self) -> f64 {
        if self.is_zero() {
            return if self.neg { -0.0 } else { 0.0 };
        }
        // value = sig * 2^(exp - 63); the double's lowest kept bit is 2^max(exp - 52, -1074)
        let lsb = (self.exp - 52).max(-1074);
        let drop = (lsb - (self.exp - 63)) as u32; // >= 11
        let mut m: u64;
        let mut l = lsb;
        if drop >= 64 {
            // everything is below the lowest kept bit: 0 or the smallest unit
            let half_pos = drop - 1;
            m = if half_pos < 64 {
                let half = 1u64 << half_pos;
                if self.sig > half { 1 } else { 0 } // == half rounds to even (0)
            } else {
                0
            };
        } else {
            m = self.sig >> drop;
            let rem = self.sig & ((1u64 << drop) - 1);
            let half = 1u64 << (drop - 1);
            if rem > half || (rem == half && (m & 1) == 1) {
                m += 1;
            }
        }
        if m >= (1u64 << 53) {
            m >>= 1;
            l += 1;
        }
        let sign = (self.neg as u64) << 63;
        if m == 0 {
            return f64::from_bits(sign);
        }
        if m >= (1u64 << 52) {
            let biased = l + 52 + 1023;
            if biased >= 2047 {
                return f64::from_bits(sign | (0x7ffu64 << 52));
            }
            return f64::from_bits(sign | ((biased as u64) << 52) | (m - (1u64 << 52)));
        }
        // subnormal (l == -1074)
        f64::from_bits(sign | m)
    }
}

/// NativeGmres.cpp `dot`: `long double sum = 0.0; for (...) sum += (long double)a[i] * b[i];
/// return (double)sum;`
pub fn dot_ext80(a: &[f64], b: &[f64]) -> f64 {
    let mut sum = Ext80::ZERO;
    for i in 0..a.len() {
        sum = sum.add(Ext80::mul_f64(a[i], b[i]));
    }
    sum.to_f64()
}
