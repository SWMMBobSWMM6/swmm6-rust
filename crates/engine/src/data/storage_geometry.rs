// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/StorageGeometry.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Storage-shape raw parameters (p1, p2, p3) to the area-relation coefficients of
// A(d) = c + a*d + b*d^2, and the INP/GeoPackage keyword for each shape. Every product
// and quotient is evaluated left to right as upstream writes it; no transcendental call
// is made. storage_shape_coeffs keeps upstream's out-parameters (a, b, c are left
// untouched on refusal) and its local names A, B, Z. Keywords are byte strings. Upstream's
// fallback `return "FUNCTIONAL"` after the switch serves enum values outside the six
// enumerators, which a Rust StorageShape cannot hold, so the match is exhaustive instead.

use crate::data::node::StorageShape;

/// Pi to the precision legacy uses (consts.h `#define PI 3.141592654`), not full double pi.
#[allow(clippy::approx_constant)]
pub const K_STORAGE_PI: f64 = 3.141592654;

/// True for CYLINDRICAL, CONICAL, PARABOLOID or PYRAMIDAL (a/b/c are quadratic coefficients).
pub const fn storage_shape_is_geometric(s: StorageShape) -> bool {
    matches!(
        s,
        StorageShape::Cylindrical
            | StorageShape::Conical
            | StorageShape::Paraboloid
            | StorageShape::Pyramidal
    )
}

/// Legacy node.c:700-708: length and width > 0, slope/height >= 0, and a PARABOLOID's
/// height not 0. TABULAR/FUNCTIONAL always validate.
pub fn storage_shape_params_valid(s: StorageShape, p1: f64, p2: f64, p3: f64) -> bool {
    if !storage_shape_is_geometric(s) {
        return true;
    }
    if p1 <= 0.0 || p2 <= 0.0 || p3 < 0.0 {
        return false;
    }
    if s == StorageShape::Paraboloid && p3 == 0.0 {
        return false;
    }
    true
}

/// Derive a (=a1), b (=a2), c (=a0) for a geometric shape (legacy node.c:713-752).
/// Returns false, leaving a/b/c untouched, for TABULAR/FUNCTIONAL or invalid params.
#[allow(non_snake_case)]
pub fn storage_shape_coeffs(
    s: StorageShape,
    p1: f64,
    p2: f64,
    p3: f64,
    a: &mut f64,
    b: &mut f64,
    c: &mut f64,
) -> bool {
    if !storage_shape_is_geometric(s) {
        return false;
    }
    if !storage_shape_params_valid(s, p1, p2, p3) {
        return false;
    }

    let A = p1 / 2.0; // base/top semi-axis length
    let B = p2 / 2.0; // base/top semi-axis width
    let Z = p3; // side slope (run/rise), or height for PARABOLOID

    match s {
        StorageShape::Cylindrical => {
            *a = 0.0;
            *b = 0.0;
            *c = K_STORAGE_PI * A * B;
            true
        }
        StorageShape::Conical => {
            *a = 2.0 * K_STORAGE_PI * B * Z;
            *b = K_STORAGE_PI * B / A * Z * Z;
            *c = K_STORAGE_PI * A * B;
            true
        }
        StorageShape::Paraboloid => {
            *a = K_STORAGE_PI * A * B / Z; // Z != 0 guaranteed by params_valid
            *b = 0.0;
            *c = 0.0;
            true
        }
        StorageShape::Pyramidal => {
            *a = 2.0 * (p1 + p2) * Z; // 2(L + W)Z: full base dims, not semi-axes
            *b = 4.0 * Z * Z;
            *c = p1 * p2; // L*W
            true
        }
        _ => false,
    }
}

/// Canonical INP/GeoPackage keyword for a shape. PARABOLOID's keyword is PARABOLIC, as
/// legacy writes and reads it (text.h:195).
pub fn storage_shape_keyword(s: StorageShape) -> &'static [u8] {
    match s {
        StorageShape::Tabular => b"TABULAR",
        StorageShape::Functional => b"FUNCTIONAL",
        StorageShape::Cylindrical => b"CYLINDRICAL",
        StorageShape::Conical => b"CONICAL",
        StorageShape::Paraboloid => b"PARABOLIC",
        StorageShape::Pyramidal => b"PYRAMIDAL",
    }
}

/// Parse a complete, already upper-cased shape keyword; PARABOLOID is accepted as an alias
/// for PARABOLIC. `out` is assigned on success and left as it was on refusal.
pub fn storage_shape_from_keyword(kw: &[u8], out: &mut StorageShape) -> bool {
    if kw == b"TABULAR" {
        *out = StorageShape::Tabular;
        return true;
    }
    if kw == b"FUNCTIONAL" {
        *out = StorageShape::Functional;
        return true;
    }
    if kw == b"CYLINDRICAL" {
        *out = StorageShape::Cylindrical;
        return true;
    }
    if kw == b"CONICAL" {
        *out = StorageShape::Conical;
        return true;
    }
    if kw == b"PARABOLIC" || kw == b"PARABOLOID" {
        *out = StorageShape::Paraboloid;
        return true;
    }
    if kw == b"PYRAMIDAL" {
        *out = StorageShape::Pyramidal;
        return true;
    }
    false
}

/// True for an integer in the inclusive TABULAR..PYRAMIDAL range.
pub const fn storage_shape_is_valid_code(v: i32) -> bool {
    v >= StorageShape::Tabular as i32 && v <= StorageShape::Pyramidal as i32
}
