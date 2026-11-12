// SPDX-License-Identifier: Apache-2.0
//
// The 2D surface routing module: upstream src/engine/2d (HydroCouple/openswmm.engine, branch
// swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache License,
// Version 2.0. Compiled only with the `twod` feature, as upstream compiles it only with
// OPENSWMM_HAS_2D. Each submodule names the upstream files it derives from.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) at the request of Robert
// Dickinson. Not an OpenSWMM release.

use std::cell::RefCell;
use std::rc::Rc;

// Phase 2a wave 2: the leaf numerics (each gated against upstream by a gate/diff/twod_*
// differential; gate/twodkern.sh builds and runs them all).
pub mod coupling;
pub mod mesh;
pub mod solver;
// Phase 2a wave 1: the 2D sections, their data and the router members open() fills.
pub mod data;
pub mod gw;
pub mod heat;
pub mod infil;
pub mod input;
pub mod quality;
pub mod subsurface;
pub mod surface_router_2d;
// Phase 2a wave 3: SurfaceRouter2D::initialize and what it calls.
pub mod surface_router_2d_init;
// Phase 2a wave 4a: the 2D step on an uncoupled mesh.
pub mod surface_router_2d_step;

/// A member of the 2D router that upstream reaches through raw pointers or captured references
/// as well as through its owner (SimulationContext::twod_io, the section handlers): one
/// shared, interior-mutable value, held by the router and cloned into each user.
pub type Shared<T> = Rc<RefCell<T>>;

/// A new shared value.
pub fn shared<T>(v: T) -> Shared<T> {
    Rc::new(RefCell::new(v))
}
