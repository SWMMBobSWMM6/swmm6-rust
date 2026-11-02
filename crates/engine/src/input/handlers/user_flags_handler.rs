// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/handlers/UserFlagsHandler.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.

#![allow(clippy::all)]

use crate::core::context::SimulationContext;
use crate::core::options::Bytes;
use crate::core::user_flags::{UserFlagDef, UserFlagType};
use crate::input::tokenizer;

/// [USER_FLAGS]: `Name Type [Description]`. Names and types are upper-cased; a duplicate
/// name replaces the earlier definition. An unknown type becomes STRING and sets
/// warning_code 102 when no warning is recorded yet.
pub fn handle_user_flags(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 2 {
            continue; // Need at least Name + Type
        }

        let mut def = UserFlagDef::default();
        def.name = tokenizer::to_upper(&tok[0]);
        let type_str = tokenizer::to_upper(&tok[1]);

        if type_str == b"BOOLEAN" {
            def.type_ = UserFlagType::Boolean;
        } else if type_str == b"INTEGER" {
            def.type_ = UserFlagType::Integer;
        } else if type_str == b"REAL" {
            def.type_ = UserFlagType::Real;
        } else if type_str == b"STRING" {
            def.type_ = UserFlagType::String;
        } else {
            // Unknown type — default to STRING and emit a warning
            def.type_ = UserFlagType::String;
            if ctx.warning_code == 0 {
                ctx.warning_code = 102;
            }
        }

        // Optional description in tok[2] (already unquoted by tokenizer)
        if tok.len() > 2 {
            def.description = tok[2].clone();
        }

        ctx.user_flags.define(def);
    }
}
