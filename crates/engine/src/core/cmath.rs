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
