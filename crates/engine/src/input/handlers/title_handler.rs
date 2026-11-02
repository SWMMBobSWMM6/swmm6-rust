// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/handlers/TitleHandler.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.

#![allow(clippy::all)]

use crate::core::context::SimulationContext;
use crate::core::options::Bytes;

/// [TITLE]: replace the title notes with an exact copy of the section lines (blank and
/// comment text preserved, no tokenizing, no three-line cap).
pub fn handle_title(ctx: &mut SimulationContext, lines: &[Bytes]) {
    ctx.title_notes.clear();
    ctx.title_notes.reserve(lines.len());
    for line in lines {
        ctx.title_notes.push(line.clone());
    }
}
