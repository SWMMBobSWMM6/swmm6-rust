// SPDX-License-Identifier: Apache-2.0
//
// The 2D subsurface authoring (upstream src/engine/2d/subsurface): the [2D_AQUIFER*] and
// [2D_SURFACE_OWNERSHIP] rows open() parses. Each submodule names the upstream files it derives
// from. The kernel itself is phase 2b.

pub mod subsurface_data;
pub mod subsurface_sections;
pub mod subsurface_solver;
pub mod surface_ownership;
