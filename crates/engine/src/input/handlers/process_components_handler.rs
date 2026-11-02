// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/handlers/ProcessComponentsHandler.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.

#![allow(clippy::all)]

use crate::core::context::{ProcessComponentSpec, SimulationContext};
use crate::core::options::Bytes;
use crate::input::tokenizer;

/// Split one component argument at its first '='. A nonempty key is required; one
/// matching pair of surrounding double quotes is removed from the value. On refusal
/// `key` and `val` are left as they were.
fn split_kv(tok: &[u8], key: &mut Bytes, val: &mut Bytes) -> bool {
    let eq = match tok.iter().position(|&c| c == b'=') {
        None => return false,
        Some(0) => return false,
        Some(p) => p,
    };
    *key = tok[..eq].to_vec();
    *val = tok[eq + 1..].to_vec();
    if val.len() >= 2 && val[0] == b'"' && val[val.len() - 1] == b'"' {
        *val = val[1..val.len() - 1].to_vec();
    }
    true
}

/// [PROCESS_COMPONENTS]: append one registration per nonempty row. CONFIG (any case) sets
/// the config path, last one winning; other key=value pairs keep spelling and order. A
/// malformed argument appends an error and the row still registers.
pub fn handle_process_components(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.is_empty() {
            continue;
        }

        let mut spec = ProcessComponentSpec::default();
        spec.id = tok[0].clone();

        for i in 1..tok.len() {
            let mut key: Bytes = Vec::new();
            let mut val: Bytes = Vec::new();
            if !split_kv(&tok[i], &mut key, &mut val) {
                let mut e = b"[PROCESS_COMPONENTS] '".to_vec();
                e.extend_from_slice(&spec.id);
                e.extend_from_slice(b"': argument '");
                e.extend_from_slice(&tok[i]);
                e.extend_from_slice(b"' is not a key=\"value\" pair.");
                ctx.errors.push(e);
                continue;
            }
            if tokenizer::to_upper(&key) == b"CONFIG" {
                spec.config_path = val;
            } else {
                spec.args.push((key, val));
            }
        }

        ctx.process_component_specs.push(spec);
    }
}
