// SPDX-License-Identifier: Apache-2.0
//
// The C++ standard library's comparisons, as upstream calls them. Rust's f64::min/max
// and clamp treat NaN differently (they return the other operand), so a translation that
// used them would differ from upstream whenever a NaN reaches one. These follow the
// definitions in the C++ standard exactly:
//   std::min(a, b)        = (b < a) ? b : a
//   std::max(a, b)        = (a < b) ? b : a
//   std::clamp(v, lo, hi) = (v < lo) ? lo : (hi < v) ? hi : v
// Written for swmm6-rust, 2 November 2026, by Claude (Anthropic) for Robert Dickinson.

#[inline(always)]
pub fn cmin<T: PartialOrd>(a: T, b: T) -> T {
    if b < a {
        b
    } else {
        a
    }
}

#[inline(always)]
pub fn cmax<T: PartialOrd>(a: T, b: T) -> T {
    if a < b {
        b
    } else {
        a
    }
}

#[inline(always)]
pub fn cclamp<T: PartialOrd>(v: T, lo: T, hi: T) -> T {
    if v < lo {
        lo
    } else if hi < v {
        hi
    } else {
        v
    }
}

/// `static_cast<int>(x)` as the reference build executes it (x86-64 `cvttsd2si`): truncation
/// toward zero for in-range values, and INT_MIN (0x80000000) for NaN, infinities and
/// anything out of range. Rust's `as i32` saturates instead (NaN -> 0, +inf -> i32::MAX), so
/// use this wherever the operand can be non-finite or out of range. (In C++ such a cast is
/// undefined; this is what g++ -O2 emits on the reference machine.)
#[inline(always)]
pub fn cvt_i32(x: f64) -> i32 {
    if x.is_nan() || x >= 2147483648.0 || x < -2147483648.0 {
        i32::MIN
    } else {
        x as i32
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn nan_follows_cxx() {
        assert!(cmin(f64::NAN, 1.0).is_nan()); // (1 < NaN) false -> a
        assert_eq!(cmin(1.0, f64::NAN), 1.0); // (NaN < 1) false -> a
        assert!(cmax(f64::NAN, 1.0).is_nan());
        assert_eq!(cmax(1.0, f64::NAN), 1.0);
        assert!(cclamp(f64::NAN, 0.0, 1.0).is_nan());
        // signed zeros: std::min(0.0, -0.0) is 0.0 (first), Rust's min may give either
        assert!(cmin(0.0f64, -0.0).is_sign_positive());
        assert!(cmin(-0.0f64, 0.0).is_sign_negative());
    }
}

/// glibc's cbrt (sysdeps/ieee754/dbl-64/s_cbrt.c, glibc 2.39), the one the reference build
/// links. Rust's f64::cbrt is not glibc's (it differs in the last bit on about half of all
/// inputs, measured on 200,000 values), so std::cbrt translates to this. The polynomial,
/// the Newton-Halley step and the order of every operation are glibc's; frexp and ldexp are
/// exact bit operations (a cube root of a finite nonzero double is always a normal double).
pub fn glibc_cbrt(x: f64) -> f64 {
    const CBRT2: f64 = 1.2599210498948731648; // 2^(1/3)
    const SQR_CBRT2: f64 = 1.5874010519681994748; // 2^(2/3)
    const FACTOR: [f64; 5] = [1.0 / SQR_CBRT2, 1.0 / CBRT2, 1.0, CBRT2, SQR_CBRT2];
    // xm = frexp(fabs(x), &xe); glibc's frexp gives xe == 0 for Inf and NaN
    let ax = x.abs();
    if ax == 0.0 || !ax.is_finite() {
        return x + x;
    }
    let bits = ax.to_bits();
    let mut e = ((bits >> 52) & 0x7ff) as i32;
    let mut frac = bits & ((1u64 << 52) - 1);
    if e == 0 {
        // subnormal: normalize
        let lz = frac.leading_zeros() as i32 - 11; // shift that puts the top bit at bit 52
        frac = (frac << lz) & ((1u64 << 52) - 1);
        e = 1 - lz;
    }
    let xe = e - 1022; // ax = xm * 2^xe, xm in [0.5, 1)
    let xm = f64::from_bits((1022u64 << 52) | frac);
    let u = 0.354895765043919860 + (1.50819193781584896 + (-2.11499494167371287 + (2.44693122563534430 + (-1.83469277483613086 + (0.784932344976639262 - 0.145263899385486377 * xm) * xm) * xm) * xm) * xm) * xm;
    let t2 = u * u * u;
    let ym = u * (t2 + 2.0 * xm) / (2.0 * t2 + xm) * FACTOR[(2 + xe % 3) as usize];
    let y = if x > 0.0 { ym } else { -ym };
    // ldexp(y, xe / 3): xe / 3 is within [-358, 342], so 2^n is a normal double and exact
    let n = xe / 3;
    y * f64::from_bits(((n + 1023) as u64) << 52)
}
