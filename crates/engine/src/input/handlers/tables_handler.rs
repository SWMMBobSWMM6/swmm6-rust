// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/handlers/TablesHandler.cpp and
// TablesHandler.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Section handlers for [TIMESERIES] and [CURVES]. A `Name FILE "path"` row only records
// the path in Table::file_path (FilePathPair assignment from a string: original = path,
// absolute cleared); nothing here opens the file — that is PostParseResolver's and
// TableData.cpp's job. The per-series anchor map (an unordered_map upstream) is only
// looked up and assigned, never iterated, so a HashMap gives the same result. std::fmod
// is `%` on f64 and std::floor is `floor`; both are exact.

use std::collections::HashMap;

use crate::core::charconv::{from_chars_double, from_chars_u32};
use crate::core::cmath::cvt_i32;
use crate::core::context::SimulationContext;
use crate::core::datetime;
use crate::core::error_codes::{format_error, ERR_KEYWORD};
use crate::core::options::Bytes;
use crate::data::name_index::ieq;
use crate::data::table::TableType;
use crate::input::parse_utils::{parse_datetime, to_double};
use crate::input::section_parser::parse_section;
use crate::input::tokenizer;

// ============================================================================
// handle_timeseries()
// ============================================================================

/// Per-series date anchor (legacy Tseries[j].lastDate).
#[derive(Clone, Copy, Default)]
struct SeriesParseState {
    last_date: f64,
    seen_date: bool,
}

pub fn handle_timeseries(ctx: &mut SimulationContext, lines: &[Bytes]) {
    // Rows are: Name Date Time Value | Name Time Value | Name FILE "path".
    let mut current_name: Bytes = Vec::new();
    let mut last_date = 0.0f64;
    let mut current_idx: i32 = -1;

    // The anchor is kept per series; locals on the hot path, parked in the map when the
    // row's series changes. Date-less rows before a series' first explicit date are
    // elapsed times stored relative to 0 and counted in Table::n_relative.
    let mut series_state: HashMap<i32, SeriesParseState> = HashMap::new();
    let mut seen_date = false;

    for pl in parse_section(lines) {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.is_empty() {
            continue;
        }

        // First column: name (non-empty) or continuation (empty — same name)
        let maybe_name = &tok[0];
        let mut is_new_table = false;
        if !maybe_name.is_empty() && !(current_idx >= 0 && ieq(maybe_name, &current_name)) {
            // Park the outgoing series' anchor state before switching.
            if current_idx >= 0 {
                series_state.insert(current_idx, SeriesParseState { last_date, seen_date });
            }
            current_name = maybe_name.clone();
            // Kind-scoped: a curve with the same name is a different object.
            current_idx = ctx.find_timeseries(&current_name);
            if current_idx < 0 {
                current_idx = ctx.tables.add(&current_name, TableType::Timeseries);
                is_new_table = true;
            }
            if let Some(st) = series_state.get(&current_idx) {
                last_date = st.last_date;
                seen_date = st.seen_date;
            } else if current_idx < ctx.n_tables() && !ctx.tables[current_idx].x.is_empty() {
                // Series resumed from an earlier [TIMESERIES] section: reconstruct the
                // anchor from its stored rows.
                let prev = &ctx.tables[current_idx];
                seen_date = prev.n_relative < prev.x.len() as i32;
                last_date = prev.x[prev.x.len() - 1].floor();
            } else {
                last_date = 0.0;
                seen_date = false;
            }
        }
        // Attach comment to the table on its first (name-introducing) row
        if is_new_table && !pl.comment.is_empty() && current_idx >= 0 {
            ctx.tables[current_idx].comment = pl.comment.clone();
        }

        if current_idx < 0 || current_name.is_empty() {
            continue;
        }
        let tbl = &mut ctx.tables[current_idx];

        if tok.len() < 2 {
            continue;
        }

        // Detect FILE reference
        let tok1_upper = tokenizer::to_upper(&tok[1]);
        if tok1_upper == b"FILE" {
            // External file: tok[2] is the quoted path, optionally with a ':column'
            // suffix; stored verbatim (outer quotes removed).
            if tok.len() > 2 {
                let mut path: &[u8] = &tok[2];
                if path.len() >= 2 && path[0] == b'"' && path[path.len() - 1] == b'"' {
                    path = &path[1..path.len() - 1];
                }
                tbl.file_path.set(path);
            }
            continue;
        }

        // Detect whether tok[1] looks like a date (contains '/')
        let has_date = tok[1].contains(&b'/');

        let x: f64;
        let y: f64;

        if has_date && tok.len() >= 4 {
            // Name  Date  Time  Value
            x = parse_datetime(&tok[1], &tok[2]);
            // Carry only the DATE (midnight) to subsequent date-less rows.
            last_date = x.floor();
            seen_date = true;
            y = to_double(&tok[3], 0.0);
        } else if !has_date && tok.len() >= 3 {
            // Name  Time  Value: decimal hours if the whole token is a number, otherwise
            // integer HH[:MM[:SS]].
            let time_frac: f64;
            let t = &tok[1];
            let (dec, used) = from_chars_double(t);
            match dec {
                Ok(dec_hours) if used == t.len() => {
                    // Whole token is a number -> decimal hours -> fraction of day.
                    time_frac = dec_hours / 24.0;
                }
                _ => {
                    let (mut th, mut tm) = (0u32, 0u32);
                    let mut ts = 0.0f64;
                    let mut tp = 0usize;
                    let rut = |tp: &mut usize, out: &mut u32| {
                        let (r, n) = from_chars_u32(&t[*tp..]);
                        if let Ok(v) = r {
                            *out = v;
                            *tp += n;
                        }
                    };
                    rut(&mut tp, &mut th);
                    if tp < t.len() && t[tp] == b':' {
                        tp += 1;
                        rut(&mut tp, &mut tm);
                    }
                    if tp < t.len() && t[tp] == b':' {
                        tp += 1;
                        let (r, n) = from_chars_double(&t[tp..]);
                        if let Ok(v) = r {
                            ts = v;
                            tp += n;
                        }
                    }
                    let _ = tp;
                    time_frac = datetime::encode_time(th as i32, tm as i32, cvt_i32(ts));
                }
            }

            // If the time wraps back to (or before) the previous time of day, increment
            // the day.
            if !tbl.x.is_empty() && time_frac <= tbl.x[tbl.x.len() - 1] % 1.0 {
                last_date += 1.0;
            }
            x = last_date + time_frac;
            y = to_double(&tok[2], 0.0);
            // Date-less rows before the series' first explicit date are counted so the
            // resolver adds start_date to exactly these rows.
            if !seen_date {
                tbl.n_relative = tbl.x.len() as i32 + 1;
            }
        } else {
            continue;
        }

        tbl.x.push(x);
        tbl.y.push(y);
    }
}

// ============================================================================
// handle_curves()
// ============================================================================

/// CURVE_TYPE_MAP: exact (already upper-cased) curve type keyword.
fn curve_type(word: &[u8]) -> Option<TableType> {
    match word {
        b"STORAGE" => Some(TableType::CurveStorage),
        b"DIVERSION" => Some(TableType::CurveDiversion),
        b"RATING" => Some(TableType::CurveRating),
        b"SHAPE" => Some(TableType::CurveShape),
        b"CONTROL" => Some(TableType::CurveControl),
        b"TIDAL" => Some(TableType::CurveTidal),
        b"WEIR" => Some(TableType::CurveWeir),
        b"PUMP1" => Some(TableType::CurvePump1),
        b"PUMP2" => Some(TableType::CurvePump2),
        b"PUMP3" => Some(TableType::CurvePump3),
        b"PUMP4" => Some(TableType::CurvePump4),
        b"PUMP5" => Some(TableType::CurvePump5),
        _ => None,
    }
}

pub fn handle_curves(ctx: &mut SimulationContext, lines: &[Bytes]) {
    let mut current_name: Bytes = Vec::new();
    let mut current_idx: i32 = -1;
    let mut current_type = TableType::CurveRating;

    for pl in parse_section(lines) {
        let tok = tokenizer::tokenize(&pl.data);
        if tok.is_empty() {
            continue;
        }

        let maybe_name = &tok[0];
        let mut is_new_table = false;
        // Skip the lookup when this row names the table we are already on.
        if !maybe_name.is_empty() && !(current_idx >= 0 && ieq(maybe_name, &current_name)) {
            current_name = maybe_name.clone();
            // Kind-scoped: a timeseries with the same name is a different object.
            current_idx = ctx.find_curve(&current_name);
            if current_idx < 0 {
                // Type appears in tok[1] on a new curve's first row; an unknown type is
                // ERR_KEYWORD but the table is still created with current_type.
                if tok.len() > 1 {
                    match curve_type(&tokenizer::to_upper(&tok[1])) {
                        Some(t) => current_type = t,
                        None => ctx.errors.push(format_error(ERR_KEYWORD, &tok[1])),
                    }
                }
                current_idx = ctx.tables.add(&current_name, current_type);
                is_new_table = true;
            }
        }

        if current_idx < 0 {
            continue;
        }
        let tbl = &mut ctx.tables[current_idx];
        if is_new_table && !pl.comment.is_empty() {
            tbl.comment = pl.comment.clone();
        }

        // After a recognised type token, the remaining tokens are x-y pairs.
        let mut data_start = 1usize;
        if tok.len() >= 2 && curve_type(&tokenizer::to_upper(&tok[1])).is_some() {
            data_start = 2; // type token consumed
        }

        // Pairs: x y [x y ...]
        let mut i = data_start;
        while i + 1 < tok.len() {
            tbl.x.push(to_double(&tok[i], 0.0));
            tbl.y.push(to_double(&tok[i + 1], 0.0));
            i += 2;
        }
    }
}
