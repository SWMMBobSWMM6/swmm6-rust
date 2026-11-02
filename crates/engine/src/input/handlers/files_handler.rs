// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/handlers/FilesHandler.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// [FILES] — secondary file references (legacy reference: iface.c).

#![allow(clippy::all)]

use crate::core::context::{FileMode, HotstartFormat, HotstartSaveEntry, SimulationContext};
use crate::core::error_codes::{format_error, format_error_detail, ERR_DATETIME, ERR_ITEMS};
use crate::core::options::Bytes;
use crate::input::parse_utils::{parse_date, parse_datetime, parse_double_strict};
use crate::input::tokenizer;

/// Remove one matching pair of leading/trailing single or double quotes. The tokenizer
/// already removes the quotes of a token that starts with one; this catches the rest, as
/// the legacy parser does.
fn unquote(s: &[u8]) -> Bytes {
    let mut s = s.to_vec();
    if s.len() >= 2 && (s[0] == b'"' || s[0] == b'\'') && s[s.len() - 1] == s[0] {
        s = s[1..s.len() - 1].to_vec();
    }
    s
}

/// The optional date+time pair on a SAVE HOTSTART row; 0.0 (write at end of run) when
/// fewer than two tokens remain.
fn parse_optional_datetime(tok: &[Bytes], idx: usize) -> f64 {
    if idx + 1 >= tok.len() {
        return 0.0;
    }
    parse_datetime(&tok[idx], &tok[idx + 1])
}

/// A time token legacy datetime_strToTime accepts: decimal hours, or H[:M[:S]] in digits.
fn is_time_token(s: &[u8]) -> bool {
    let mut hours = 0.0;
    if parse_double_strict(s, &mut hours) {
        return true;
    }
    let mut fields = 0;
    let mut digit = false;
    for &c in s {
        if c >= b'0' && c <= b'9' {
            digit = true;
            continue;
        }
        // `c != ':' || !digit || ++fields > 2`, short-circuit as in C++
        if c != b':' || !digit {
            return false;
        }
        fields += 1;
        if fields > 2 {
            return false;
        }
        digit = false;
    }
    digit
}

/// [FILES]: `SAVE|USE  RAINFALL|RUNOFF|RDII|INFLOWS|OUTFLOWS|HOTSTART  path [...]`.
pub fn handle_files(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 3 {
            continue; // need at least mode + kind + path
        }

        let mode_word = tokenizer::to_upper(&tok[0]);
        let kind_word = tokenizer::to_upper(&tok[1]);
        let path = unquote(&tok[2]);

        let mode;
        if mode_word == b"SAVE" {
            mode = FileMode::Save;
        } else if mode_word == b"USE" {
            mode = FileMode::Use;
        } else {
            continue; // unrecognised mode — skip silently (legacy behaviour)
        }

        let files = &mut ctx.files;
        if kind_word == b"RAINFALL" {
            files.rainfall_mode = mode;
            files.rainfall_path.set(&path);
        } else if kind_word == b"RUNOFF" {
            files.runoff_mode = mode;
            files.runoff_path.set(&path);
        } else if kind_word == b"RDII" {
            files.rdii_mode = mode;
            files.rdii_path.set(&path);
        } else if kind_word == b"INFLOWS" {
            // Legacy: USE only — SAVE INFLOWS is a fatal input error
            // (iface.c iface_readFileParams returns ERR_ITEMS).
            if mode != FileMode::Use {
                ctx.errors.push(format_error_detail(ERR_ITEMS, b"", b"[FILES] SAVE INFLOWS"));
                continue;
            }
            files.inflows_path.set(&path);
        } else if kind_word == b"OUTFLOWS" {
            // Legacy: SAVE only — USE OUTFLOWS is a fatal input error.
            if mode != FileMode::Save {
                ctx.errors.push(format_error_detail(ERR_ITEMS, b"", b"[FILES] USE OUTFLOWS"));
                continue;
            }
            files.outflows_path.set(&path);
        } else if kind_word == b"HOTSTART" {
            if mode == FileMode::Save {
                // Legacy iface.c reads a datetime only from a five-token row
                // and rejects one whose date or time does not parse.
                if tok.len() == 5 {
                    if parse_date(&tok[3]) <= 0.0 {
                        ctx.errors.push(format_error(ERR_DATETIME, &tok[3]));
                        continue;
                    }
                    if !is_time_token(&tok[4]) {
                        ctx.errors.push(format_error(ERR_DATETIME, &tok[4]));
                        continue;
                    }
                }
                let mut entry = HotstartSaveEntry::default();
                entry.path.set(&path);
                entry.datetime = parse_optional_datetime(&tok, 3);
                // Optional trailing format keyword, after the datetime if
                // there is one: `SAVE HOTSTART f [date time] LEGACY|NATIVE`.
                // It must trail: legacy SWMM 5 reads a datetime only from a
                // five-token row and ignores later tokens, so a four- or
                // six-token row still opens there. Any other token in that
                // position is ignored, with a warning.
                let fmt_at: usize = if tok.len() >= 6 { 5 } else { 3 };
                if tok.len() == fmt_at + 1 {
                    let fmt = tokenizer::to_upper(&tok[fmt_at]);
                    if fmt == b"LEGACY" {
                        entry.format = HotstartFormat::Legacy;
                    } else if fmt == b"NATIVE" {
                        entry.format = HotstartFormat::Native;
                    } else {
                        let mut w = b"\n  WARNING: Unknown SAVE HOTSTART format '".to_vec();
                        w.extend_from_slice(&tok[fmt_at]);
                        w.extend_from_slice(
                            b"' in [FILES] section - LEGACY or NATIVE expected; \
the default format will be used.",
                        );
                        ctx.warnings.push(w);
                        if ctx.warning_code == 0 {
                            ctx.warning_code = 101;
                        }
                    }
                }
                files.hotstart_saves.push(entry);
            } else {
                // USE
                files.hotstart_use_path.set(&path);
            }
        }
        // Unknown kinds are silently ignored — preserves the legacy
        // behaviour of skipping rows that aren't part of the spec.
    }
}
