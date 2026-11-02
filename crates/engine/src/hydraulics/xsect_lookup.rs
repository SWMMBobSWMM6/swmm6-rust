// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/XSectLookup.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Bit-exact geometry-table interpolation (legacy xsect.c lookup()) and the bucket LUT
// that narrows locate()'s bisection. Every product and quotient keeps upstream's
// grouping: IEEE division by `delta`, never a reciprocal multiply, and no fused
// multiply-add.
//
// Differences in form, never in result on upstream's defined domain:
//   * A `const double*` table plus an `int` count becomes a `&[f64]` plus the same `i32`
//     count; the count is still the logical length, so a longer slice behaves as the
//     pointer did. Indexing is bounds-checked, so an out-of-range read that is undefined
//     behaviour in C++ (lookup_exact with n < 4 and a low x, or a negative x below
//     -delta) panics here instead.
//   * `static_cast<int>` of a double becomes `as i32`. The two agree whenever the value
//     fits an int. Outside that (NaN, infinities, overflow) C++ is undefined; x86 yields
//     INT_MIN, while Rust saturates (NaN -> 0), as ARM does. Upstream's callers guard
//     against that case (see xsect::lookup in XSectKernels.hpp).
//   * A nullable `const double*` / `const LocateLut*` argument becomes an `Option`.
//
// Compile-time switch: upstream's CMake option OPENSWMM_FAST_XSECT_LOOKUP defines
// SWMM_XSECT_FAST_LOOKUP, which makes norm_x() multiply by a precomputed 1/y_full and
// routes norm_lookup() through lookup_fast() (index and interpolate with `* inv_delta`,
// pre-clamping x to [0, 1]). The reference build has it OFF, so only the default IEEE
// divide branch is translated. lookup_fast() itself is always compiled upstream but is
// wired into the kernels only under that switch (otherwise only the unit test
// tests/unit/engine/test_xsect_parity.cpp calls it); it is not translated here.

// The negated comparisons (`!(a >= b)`, `!(y > t0)`) are upstream's NaN handling and
// stay as written; so do upstream's range tests.
#![allow(clippy::neg_cmp_op_on_partial_ord, clippy::manual_range_contains)]

/// Interpolate a normalized geometry table with the legacy arithmetic grouping
/// (upstream `lookup_exact`). `t` holds `n` equally spaced entries on 0..1; `n >= 4`
/// for the low-end quadratic lookahead. Values at or above the upper extent return
/// the final entry; the result is floored at 0.
#[inline]
pub fn lookup_exact(x: f64, t: &[f64], n: i32) -> f64 {
    let delta = 1.0 / ((n - 1) as f64); // same rounding as legacy
    let i = (x / delta) as i32; // DIVIDE (legacy), not x*inv_delta
    if i >= n - 1 {
        return t[(n - 1) as usize];
    }
    let iu = i as usize;

    let x0 = (i as f64) * delta;
    let mut y = t[iu] + (x - x0) * (t[iu + 1] - t[iu]) / delta; // grouping matches legacy

    if i < 2 {
        // quadratic refinement for low x
        let x1 = ((i as f64) + 1.0) * delta;
        let y2 = y
            + (x - x0) * (x - x1) / (delta * delta)
                * (t[iu] / 2.0 - t[iu + 1] + t[iu + 2] / 2.0);
        if y2 > 0.0 {
            y = y2;
        }
    }
    if y < 0.0 {
        y = 0.0;
    }
    y
}

// lookup_fast(x, t, n) — upstream's opt-in reciprocal-multiply variant (NOT bit-exact,
// SWMM_XSECT_FAST_LOOKUP); not translated, see the header note.

/// Normalize a depth with the default (bit-exact) operand convention: `param` is the
/// full depth and x = depth/param by IEEE divide, or 0 unless both are positive.
///
/// Upstream under SWMM_XSECT_FAST_LOOKUP instead takes `param` as 1/y_full and returns
/// depth*param clamped to [0, 1]; that branch is not translated.
#[inline]
pub fn norm_x(depth: f64, param: f64) -> f64 {
    if param > 0.0 && depth > 0.0 {
        depth / param // param == y_full
    } else {
        0.0
    }
}

/// Normalize a depth and interpolate (default build: `lookup_exact`; upstream's
/// SWMM_XSECT_FAST_LOOKUP build would call `lookup_fast`).
#[inline]
pub fn norm_lookup(depth: f64, param: f64, t: &[f64], n: i32) -> f64 {
    lookup_exact(norm_x(depth, param), t, n)
}

// ============================================================================
// Bucket LUT for locate() — plan XSECT_LOOKUP_ACCEL §4 item A1.
//
// The table's VALUE range is split into K_BUCKETS equal spans; for each span `lo`
// holds an index bracket that provably contains locate()'s answer, so a short
// bisection inside it returns the same index plain bisection would. A table that
// fails the monotonicity check (or is flat) leaves `scale` at 0, which callers read
// as "no LUT" and fall back to plain bisection.
// ============================================================================

/// Precomputed value→index-bracket map for one geometry table (upstream `LocateLut`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LocateLut {
    /// table[0] — the low end of the value range
    pub t0: f64,
    /// K_BUCKETS / (table[jLast] - t0); 0 == disabled
    pub scale: f64,
    /// the `jLast` this map was built for
    pub j_last: i32,
    /// lo[b] = max{ j : bucket(table[j]) < b }
    pub lo: [u8; LocateLut::K_BUCKETS as usize + 1],
}

impl LocateLut {
    pub const K_BUCKETS: i32 = 64;
}

// By hand: Rust derives Default for arrays of at most 32 elements only.
impl Default for LocateLut {
    fn default() -> Self {
        LocateLut { t0: 0.0, scale: 0.0, j_last: 0, lo: [0u8; LocateLut::K_BUCKETS as usize + 1] }
    }
}

/// Map a value to one of 64 locate-acceleration buckets: truncate, then clamp.
#[inline]
pub fn lut_bucket(l: &LocateLut, y: f64) -> i32 {
    let mut b = ((y - l.t0) * l.scale) as i32;
    if b < 0 {
        b = 0;
    } else if b >= LocateLut::K_BUCKETS {
        b = LocateLut::K_BUCKETS - 1;
    }
    b
}

/// Find a lower bracketing index in a nondecreasing table by bisection
/// (0 at/below the first entry, `j_last` at/above the last; NaN resolves to 0).
#[inline]
pub fn locate_bisect(y: f64, table: &[f64], j_last: i32) -> i32 {
    let mut j1: i32 = 0;
    let mut j2: i32 = j_last;
    if y <= table[0] {
        return 0;
    }
    if y >= table[j_last as usize] {
        return j_last;
    }
    while j2 - j1 > 1 {
        let j = (j1 + j2) >> 1;
        if y >= table[j as usize] {
            j1 = j;
        } else {
            j2 = j;
        }
    }
    j1
}

/// Build the byte-indexed acceleration map for a monotone table. `None` (upstream's
/// NULL), `j_last` outside 2..=254, a non-monotone table (including NaN) or a flat one
/// leaves the map disabled (`scale == 0`).
#[inline]
pub fn build_locate_lut(l: &mut LocateLut, table: Option<&[f64]>, j_last: i32) {
    *l = LocateLut::default();
    let table = match table {
        Some(t) if !(j_last < 2 || j_last > 254) => t,
        _ => return,
    };
    for j in 1..=j_last {
        let ju = j as usize;
        if !(table[ju] >= table[ju - 1]) {
            return; // non-monotone (or NaN) — no LUT
        }
    }

    let span = table[j_last as usize] - table[0];
    if !(span > 0.0) {
        return; // flat table — nothing to bracket
    }

    l.t0 = table[0];
    l.scale = (LocateLut::K_BUCKETS as f64) / span;
    l.j_last = j_last;
    for b in 0..=LocateLut::K_BUCKETS {
        let mut best: i32 = 0;
        for j in 0..=j_last {
            if lut_bucket(l, table[j as usize]) >= b {
                break; // bucket(table[.]) is sorted
            }
            best = j;
        }
        l.lo[b as usize] = best as u8;
    }
}

/// Find the same lower bracket as `locate_bisect`, narrowed first by a prebuilt map.
/// `!(y > table[0])` also catches NaN, which plain bisection resolves to 0.
#[inline]
pub fn locate_lut(y: f64, table: &[f64], j_last: i32, l: &LocateLut) -> i32 {
    if !(y > table[0]) {
        return 0;
    }
    if y >= table[j_last as usize] {
        return j_last;
    }

    let b = lut_bucket(l, y);
    let mut j1: i32 = l.lo[b as usize] as i32;
    let mut j2: i32 = (l.lo[(b + 1) as usize] as i32) + 1;
    if j2 > j_last {
        j2 = j_last;
    }
    while j2 - j1 > 1 {
        let jm = (j1 + j2) >> 1;
        if y >= table[jm as usize] {
            j1 = jm;
        } else {
            j2 = jm;
        }
    }
    j1
}

/// Use a matching enabled bucket map (`scale > 0`, same `j_last`) or fall back to plain
/// bisection. `None` is upstream's null pointer.
#[inline]
pub fn locate_maybe_lut(y: f64, table: &[f64], j_last: i32, l: Option<&LocateLut>) -> i32 {
    if let Some(m) = l {
        if m.scale > 0.0 && m.j_last == j_last {
            return locate_lut(y, table, j_last, m);
        }
    }
    locate_bisect(y, table, j_last)
}

/// Build a map for the truncated extent the inverse geometry lookup searches: when
/// table[n-3] > table[n-1] the two declining top rows are excluded first.
#[inline]
pub fn build_invlookup_lut(l: &mut LocateLut, table: Option<&[f64]>, n_items: i32) {
    let table = match table {
        Some(t) if n_items >= 4 => t,
        _ => {
            *l = LocateLut::default();
            return;
        }
    };
    let mut n = n_items;
    if table[(n - 3) as usize] > table[(n - 1) as usize] {
        n -= 2;
    }
    build_locate_lut(l, Some(table), n - 1);
}
