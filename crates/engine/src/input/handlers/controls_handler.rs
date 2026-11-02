// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/handlers/ControlsHandler.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// [CONTROLS] and [REPORT] (legacy reference: input.c readControl(), readReport()).

#![allow(clippy::all)]

use crate::core::context::SimulationContext;
use crate::core::options::Bytes;
use crate::input::tokenizer;

/// [CONTROLS]: split the section into raw text blocks, one per RULE, joined with '\n'.
/// Text before the first RULE forms a block of its own. No grammar check here.
pub fn handle_controls(ctx: &mut SimulationContext, lines: &[Bytes]) {
    // Each rule starts with "RULE name" and ends before the next "RULE" or
    // end of section. We accumulate lines into a single text block per rule.

    let mut current_block: Bytes = Vec::new();

    fn flush_block(ctx: &mut SimulationContext, current_block: &mut Bytes) {
        if !current_block.is_empty() {
            ctx.control_rules.rule_text.push(current_block.clone());
            current_block.clear();
        }
    }

    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.is_empty() {
            continue;
        }

        let keyword = tokenizer::to_upper(&tok[0]);

        if keyword == b"RULE" {
            // Flush the previous rule block (if any) and start a new one
            flush_block(ctx, &mut current_block);
        }

        // Append this line to the current block
        if !current_block.is_empty() {
            current_block.push(b'\n');
        }
        current_block.extend_from_slice(line);
    }

    // Flush the last rule block
    flush_block(ctx, &mut current_block);
}

/// A report selector: 0 for missing/NONE, 1 for ALL, 2 for an explicit name list (the
/// names are appended to `out`, which NONE and ALL never clear).
fn parse_report_list(tok: &[Bytes], start: usize, out: &mut Vec<Bytes>) -> i32 {
    if tok.len() <= start {
        return 0;
    }

    let val = tokenizer::to_upper(&tok[start]);
    if val == b"NONE" {
        return 0;
    }
    if val == b"ALL" {
        return 1;
    }

    // Individual names
    for i in start..tok.len() {
        out.push(tok[i].clone());
    }
    2
}

/// [REPORT]: object selectors and the INPUT/CONTINUITY/FLOWSTATS/CONTROLS/AVERAGES/
/// LINK_STEPS/DISABLED flags. Rows under two tokens and unknown keys are skipped.
pub fn handle_report(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 2 {
            continue;
        }

        let keyword = tokenizer::to_upper(&tok[0]);
        let opt = &mut ctx.options;

        if keyword == b"SUBCATCHMENTS" {
            opt.rpt_subcatchments = parse_report_list(&tok, 1, &mut opt.rpt_subcatch_names);
        } else if keyword == b"NODES" {
            opt.rpt_nodes = parse_report_list(&tok, 1, &mut opt.rpt_node_names);
        } else if keyword == b"LINKS" {
            opt.rpt_links = parse_report_list(&tok, 1, &mut opt.rpt_link_names);
        } else if keyword == b"INPUT" {
            opt.rpt_input = tokenizer::parse_boolean(&tok[1]);
        } else if keyword == b"CONTINUITY" {
            opt.rpt_continuity = tokenizer::parse_boolean(&tok[1]);
        } else if keyword == b"FLOWSTATS" {
            opt.rpt_flowstats = tokenizer::parse_boolean(&tok[1]);
        } else if keyword == b"CONTROLS" {
            opt.rpt_controls = tokenizer::parse_boolean(&tok[1]);
        } else if keyword == b"AVERAGES" {
            opt.rpt_averages = tokenizer::parse_boolean(&tok[1]);
        } else if keyword == b"LINK_STEPS" {
            opt.rpt_link_steps = tokenizer::parse_boolean(&tok[1]);
        } else if keyword == b"DISABLED" {
            opt.rpt_disabled = tokenizer::parse_boolean(&tok[1]);
        }
    }
}
