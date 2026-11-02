// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/handlers/PluginsHandler.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.

#![allow(clippy::all)]

use crate::core::context::{PluginSpec, SimulationContext};
use crate::core::options::Bytes;
use crate::input::tokenizer;

/// [PLUGINS]: each nonempty row appends one PluginSpec; the first token is the library
/// path, the rest are its init arguments. Nothing is loaded or resolved here.
pub fn handle_plugins(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.is_empty() {
            continue;
        }

        let mut spec = PluginSpec::default();
        spec.path = tok[0].clone(); // First token is the library path

        // Remaining tokens are init_args
        for i in 1..tok.len() {
            spec.init_args.push(tok[i].clone());
        }

        ctx.plugin_specs.push(spec);
    }
}
