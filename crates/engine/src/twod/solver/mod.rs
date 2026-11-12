// SPDX-License-Identifier: Apache-2.0
//
// The 2D surface solver (upstream src/engine/2d/solver; HydroCouple/openswmm.engine, branch
// swmm6_rel, commit 88329d6, Copyright 2026 Caleb Buahin, Apache License 2.0). Each submodule
// names the upstream file it derives from.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) at the request of Robert
// Dickinson. Not an OpenSWMM release.

pub mod diffusive_kernels;
pub mod explicit_inertial_solver;
pub mod inertial_edges;
pub mod inertial_kernels;
pub mod momentum_physics_2d;
pub mod surface_flux_calculator;
pub mod surface_solver_factory;
pub mod swe_kernels;
