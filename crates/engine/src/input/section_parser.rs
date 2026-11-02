// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/SectionParser.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Pairs each data line with the comment lines just above it. A ";;" line (the column
// header the GUI writes) discards the pending comment. Comment lines join with the two
// characters backslash and 'n' (upstream writes "\\n" in its source), not a newline.

use crate::core::options::Bytes;

#[derive(Debug, Clone, Default)]
pub struct ParsedLine {
    pub data: Bytes,
    pub comment: Bytes,
}

pub fn parse_section(raw_lines: &[Bytes]) -> Vec<ParsedLine> {
    let mut result = Vec::new();
    let mut pending: Bytes = Vec::new();
    for raw in raw_lines {
        if raw.is_empty() {
            continue;
        }
        if raw[0] == b';' {
            if raw.len() >= 2 && raw[1] == b';' {
                pending.clear();
            } else {
                if !pending.is_empty() {
                    pending.extend_from_slice(b"\\n");
                }
                pending.extend_from_slice(&raw[1..]);
            }
            continue;
        }
        result.push(ParsedLine { data: raw.clone(), comment: std::mem::take(&mut pending) });
        pending.clear();
    }
    result
}
