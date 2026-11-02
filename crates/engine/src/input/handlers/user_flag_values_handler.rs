// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/handlers/UserFlagValuesHandler.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// std::from_chars into an int or double that starts at zero leaves the zero in place on
// any failure, so a failed parse is written here as the zero default.

#![allow(clippy::all)]

use crate::core::charconv::{from_chars_double, from_chars_i32};
use crate::core::context::SimulationContext;
use crate::core::options::Bytes;
use crate::core::user_flags::{UserFlagType, UserFlagValue};
use crate::input::tokenizer;

/// [USER_FLAG_VALUES]: `ObjectType ObjectName FlagName Value`, typed by the schema
/// registered so far. An undefined flag reads as STRING and sets warning_code 103 when no
/// warning is recorded yet. No check that the object exists.
pub fn handle_user_flag_values(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 4 {
            continue; // ObjectType ObjectName FlagName Value
        }

        let obj_type = tokenizer::to_upper(&tok[0]);
        let obj_name: &[u8] = &tok[1]; // preserve case
        let flag_name = tokenizer::to_upper(&tok[2]);
        let raw_val: &[u8] = &tok[3];

        // Determine the type from the schema definition (if registered).
        // If not defined yet, default to STRING and warn.
        let mut type_ = UserFlagType::String;
        let defined = ctx.user_flags.is_defined(&flag_name);
        if defined {
            type_ = ctx.user_flags.get_def(&flag_name).expect("defined").type_;
        } else if ctx.warning_code == 0 {
            ctx.warning_code = 103;
        }

        let value = match type_ {
            UserFlagType::Boolean => UserFlagValue::Bool(tokenizer::parse_boolean(raw_val)),

            UserFlagType::Integer => {
                let mut v: i32 = 0;
                if let Ok(x) = from_chars_i32(raw_val).0 {
                    v = x;
                }
                UserFlagValue::Int(v)
            }

            UserFlagType::Real => {
                let mut v: f64 = 0.0;
                if let Ok(x) = from_chars_double(raw_val).0 {
                    v = x;
                }
                UserFlagValue::Double(v)
            }

            UserFlagType::String => UserFlagValue::String(raw_val.to_vec()),
        };

        ctx.user_flags.set(&obj_type, obj_name, &flag_name, value);
    }
}
