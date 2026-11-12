// SPDX-License-Identifier: Apache-2.0
//
// 2D heat (upstream src/engine/2d/heat; HydroCouple/openswmm.engine, branch swmm6_rel, commit
// 88329d6, Copyright 2026 Caleb Buahin, Apache License 2.0). Heat on the mesh is phase 2b: what
// is here is the state SurfaceRouter2D::initialize resets and the inactive paths it takes.
// Each submodule names the upstream file it derives from.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) at the request of Robert
// Dickinson. Not an OpenSWMM release.

pub mod channel_ground;
pub mod ground_heat;
pub mod mesh_heat;
pub mod mesh_snow;
