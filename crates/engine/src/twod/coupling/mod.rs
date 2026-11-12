// SPDX-License-Identifier: Apache-2.0
//
// 1D-2D coupling (upstream src/engine/2d/coupling; HydroCouple/openswmm.engine, branch
// swmm6_rel, commit 88329d6, Copyright 2026 Caleb Buahin, Apache License 2.0). Each submodule
// names the upstream file it derives from.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) at the request of Robert
// Dickinson. Not an OpenSWMM release.

pub mod exchange_law;
pub mod node_coupling;
pub mod outfall_exchange;
pub mod pond_model;
pub mod runtime_coupling;
