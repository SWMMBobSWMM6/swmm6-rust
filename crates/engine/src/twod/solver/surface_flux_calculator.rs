// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/solver/SurfaceFluxCalculator.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// So far only the header's two inline scalar sink ramps (upstream namespace openswmm::twoD):
// evapSink and infilSink, the depth-limited evaporation and infiltration rates with the cubic
// Hermite wet/dry ramp. The flux calculator itself (computeLimitedGradients, the boundary-edge
// fluxes and SurfaceFluxCalculator.cpp) needs the mesh and the surface state and comes with the
// step (phase 2a wave 4).
//
// Checked bit for bit against upstream's header by gate/diff/twod_sinks.cpp and
// examples/twod_sinks.rs.

/// Depth-limited evaporation sink rate, m/s (evapSink).
#[inline]
pub fn evap_sink(rate: f64, depth: f64, dry_depth: f64) -> f64 {
    if rate <= 0.0 || depth <= 0.0 {
        return 0.0;
    }
    if depth >= dry_depth {
        return rate;
    }
    let t = depth / dry_depth;
    rate * t * t * (3.0 - 2.0 * t)
}

/// Depth-limited infiltration sink rate, m/s (infilSink); the same ramp as evap_sink.
#[inline]
pub fn infil_sink(rate: f64, depth: f64, dry_depth: f64) -> f64 {
    if rate <= 0.0 || depth <= 0.0 {
        return 0.0;
    }
    if depth >= dry_depth {
        return rate;
    }
    let t = depth / dry_depth;
    rate * t * t * (3.0 - 2.0 * t)
}
