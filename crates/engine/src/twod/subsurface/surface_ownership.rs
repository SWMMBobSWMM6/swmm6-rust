// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/subsurface/SurfaceOwnership.hpp and the
// parseSurfaceOwnerLine function of SubsurfaceSections.cpp (HydroCouple/openswmm.engine,
// branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache
// License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// [2D_SURFACE_OWNERSHIP] records. resolveSurfaceOwnership and the preview structs are run-time
// (phase 2b) and not translated here.

#![allow(clippy::all)]

use crate::core::cstdlib::c_str;
use crate::core::options::Bytes;

#[derive(Debug, Clone, Default, PartialEq)]
pub struct SurfaceOwnerRecord {
    pub subcatch: Bytes,
}

/// SubsurfaceSections.cpp iequals(std::string, const char*) with the C string `b`.
fn iequals_cstr(a: &[u8], b: &[u8]) -> bool {
    let b = c_str(b);
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_ascii_uppercase() == y.to_ascii_uppercase())
}

/// parseSurfaceOwnerLine(): `Subcatchment SUBCATCH UNIFORM`, no case-insensitive duplicate.
pub fn parse_surface_owner_line(t: &[Bytes], rows: &mut Vec<SurfaceOwnerRecord>) -> Bytes {
    if t.len() != 3 || !iequals_cstr(&t[1], b"SUBCATCH") || !iequals_cstr(&t[2], b"UNIFORM") {
        return b"expected Subcatchment SUBCATCH UNIFORM; mesh conversion and unreviewed distribution are unavailable.".to_vec();
    }
    for r in rows.iter() {
        if iequals_cstr(&r.subcatch, &t[0]) {
            let mut m = b"duplicate ownership record for ".to_vec();
            m.extend_from_slice(&t[0]);
            return m;
        }
    }
    rows.push(SurfaceOwnerRecord { subcatch: t[0].clone() });
    Bytes::new()
}
