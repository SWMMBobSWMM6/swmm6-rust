// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/handlers/InflowsHandler.cpp and
// InflowsHandler.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Section handlers for [PATTERNS], [INFLOWS], [DWF], [RDII], [HYDROGRAPHS] and
// [RDII_DECAY]. Upstream's default-initialized UnitHydEntry has every member assigned
// before use except gage_name (a std::string, so empty); Default here matches.

use crate::core::charconv::from_chars_i32;
use crate::core::context::SimulationContext;
use crate::core::error_codes::{format_error, ERR_ITEMS, ERR_KEYWORD};
use crate::core::options::Bytes;
use crate::data::inflow::{PatternData, RDIIDecayEntry, UnitHydEntry};
use crate::input::parse_utils::to_double;
use crate::input::tokenizer;

/// One complete case-insensitive pattern-type token: 0 MONTHLY, 1 DAILY, 2 HOURLY,
/// 3 WEEKEND, or -1.
fn parse_pattern_type(sv: &[u8]) -> i32 {
    let upper = sv.to_ascii_uppercase();
    if upper == b"MONTHLY" {
        return 0;
    }
    if upper == b"DAILY" {
        return 1;
    }
    if upper == b"HOURLY" {
        return 2;
    }
    if upper == b"WEEKEND" {
        return 3;
    }
    -1 // not a type keyword
}

/// An authored pattern by case-insensitive identifier, or -1.
fn find_pattern(pat: &PatternData, name: &[u8]) -> i32 {
    pat.find(name)
}

// ============================================================================
// handle_patterns()
// ============================================================================

pub fn handle_patterns(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 2 {
            continue;
        }
        // Name  [Type]  Multipliers...

        let name = &tok[0];
        let existing = find_pattern(&ctx.patterns, name);

        // Determine if the second token is a type keyword or a numeric value.
        let type_code = parse_pattern_type(&tok[1]);

        if type_code >= 0 && existing < 0 {
            // New pattern with explicit type keyword
            let mut facs: Vec<f64> = Vec::new();
            for t in &tok[2..] {
                facs.push(to_double(t, 1.0));
            }
            ctx.patterns.add(name, type_code, &facs);
        } else if type_code >= 0 && existing >= 0 {
            // Same name with type keyword again — append multipliers
            let facs = &mut ctx.patterns.factors[existing as usize];
            for t in &tok[2..] {
                facs.push(to_double(t, 1.0));
            }
        } else if existing >= 0 {
            // Continuation line — second token is numeric, append all values
            let facs = &mut ctx.patterns.factors[existing as usize];
            for t in &tok[1..] {
                facs.push(to_double(t, 1.0));
            }
        } else {
            // New pattern without explicit type — default to MONTHLY
            let mut facs: Vec<f64> = Vec::new();
            for t in &tok[1..] {
                facs.push(to_double(t, 1.0));
            }
            ctx.patterns.add(name, 0, &facs);
        }
    }
}

// ============================================================================
// handle_inflows()
// ============================================================================

pub fn handle_inflows(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 3 {
            continue;
        }
        // Node  Constituent  TimeSeries  [Type]  [Mfactor]  [Sfactor]  [Baseline]  [Pattern]

        // The node may be defined in a later section: keep the raw name.
        let node_idx = ctx.node_names.find(&tok[0]);

        let constituent = &tok[1];
        let ts_name = &tok[2];

        let inflow_type: Bytes = if tok.len() > 3 { tokenizer::to_upper(&tok[3]) } else { b"FLOW".to_vec() };
        let m_factor = if tok.len() > 4 { to_double(&tok[4], 1.0) } else { 1.0 };
        let s_factor = if tok.len() > 5 { to_double(&tok[5], 1.0) } else { 1.0 };
        let baseline = if tok.len() > 6 { to_double(&tok[6], 0.0) } else { 0.0 };
        let pat: &[u8] = if tok.len() > 7 { &tok[7] } else { b"" };

        ctx.ext_inflows.add(node_idx, constituent, ts_name, &inflow_type, m_factor, s_factor, baseline, pat, &tok[0]);
    }
}

// ============================================================================
// handle_dwf()
// ============================================================================

pub fn handle_dwf(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 3 {
            continue;
        }
        // Node  Constituent  AvgValue  [Pat1]  [Pat2]  [Pat3]  [Pat4]

        // Deferred node resolution — see handle_inflows() note.
        let node_idx = ctx.node_names.find(&tok[0]);

        let constituent = &tok[1];
        let avg_value = to_double(&tok[2], 0.0);

        let p1: &[u8] = if tok.len() > 3 { &tok[3] } else { b"" };
        let p2: &[u8] = if tok.len() > 4 { &tok[4] } else { b"" };
        let p3: &[u8] = if tok.len() > 5 { &tok[5] } else { b"" };
        let p4: &[u8] = if tok.len() > 6 { &tok[6] } else { b"" };

        ctx.dwf_inflows.add(node_idx, constituent, avg_value, p1, p2, p3, p4, &tok[0]);
    }
}

// ============================================================================
// handle_rdii()
// ============================================================================

pub fn handle_rdii(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 3 {
            continue;
        }
        // Node  UHgroup  SewerArea

        // Deferred node resolution — see handle_inflows() note.
        let node_idx = ctx.node_names.find(&tok[0]);

        let uh_name = &tok[1];
        let sewer_area = to_double(&tok[2], 0.0);

        ctx.rdii_assigns.add(node_idx, uh_name, sewer_area, &tok[0]);
    }
}

// ============================================================================
// handle_hydrographs()
// ============================================================================
//   UHgroup  RainGage                                   (gage assignment line)
//   UHgroup  Month  Response  R  T  K  [Dmax  Drecov  Dinit]

const MONTHS: [&[u8]; 12] = [
    b"JAN", b"FEB", b"MAR", b"APR", b"MAY", b"JUN", b"JUL", b"AUG", b"SEP", b"OCT", b"NOV", b"DEC",
];

pub fn handle_hydrographs(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.is_empty() {
            continue;
        }
        // Legacy rdii_readUnitHydParams: a gage row has 2 items, a response row at
        // least 6; anything else is ERR_ITEMS.
        if tok.len() == 1 || (tok.len() > 2 && tok.len() < 6) {
            ctx.errors.push(format_error(ERR_ITEMS, &tok[0]));
            continue;
        }

        let uh_name = &tok[0];

        // If only 2 tokens: UHgroup  RainGage (gage assignment line)
        if tok.len() == 2 {
            ctx.unit_hyds.add_gage(uh_name, &tok[1]);
            continue;
        }

        // Month: case-insensitive three-letter prefix, ALL (-1), or a whole integer 1..12.
        let mut month: i32 = -1;
        let month_upper = tokenizer::to_upper(&tok[1]);
        if month_upper != b"ALL" {
            for (m, mname) in MONTHS.iter().enumerate() {
                // std::string::compare(0, 3, s) == 0: the first 3 bytes (or fewer, if
                // shorter) equal s.
                let n = std::cmp::min(3, month_upper.len());
                if &month_upper[..n] == *mname {
                    month = m as i32;
                    break;
                }
            }
            // Preserve the engine's numeric-month extension.
            if month < 0 {
                let (r, used) = from_chars_i32(&tok[1]);
                if let Ok(numeric_month) = r {
                    if used == tok[1].len() && (1..=12).contains(&numeric_month) {
                        month = numeric_month - 1;
                    }
                }
            }
            if month < 0 {
                ctx.errors.push(format_error(ERR_KEYWORD, &tok[1]));
                continue;
            }
        }

        // Parse response type: Short=0, Medium=1, Long=2
        let mut response: i32 = -1;
        let resp_upper = tokenizer::to_upper(&tok[2]);
        if resp_upper == b"SHORT" {
            response = 0;
        } else if resp_upper == b"MEDIUM" {
            response = 1;
        } else if resp_upper == b"LONG" {
            response = 2;
        }
        if response < 0 {
            continue;
        }

        let entry = UnitHydEntry {
            name: uh_name.clone(),
            month,
            response,
            r: to_double(&tok[3], 0.0),
            t: to_double(&tok[4], 0.0), // hours
            k: to_double(&tok[5], 0.0), // tBase/tPeak ratio
            dmax: if tok.len() > 6 { to_double(&tok[6], 0.0) } else { 0.0 },
            drecov: if tok.len() > 7 { to_double(&tok[7], 0.0) } else { 0.0 },
            dinit: if tok.len() > 8 { to_double(&tok[8], 0.0) } else { 0.0 },
            ..UnitHydEntry::default()
        };

        ctx.unit_hyds.add(&entry);
    }
}

// ============================================================================
// handle_rdii_decay()
// ============================================================================
//   UHGroup  Response  k_dep  k_0  k_T  T_ref  theta_rec  T_freeze  [SNOW snow_T snow_ddf]

pub fn handle_rdii_decay(ctx: &mut SimulationContext, lines: &[Bytes]) {
    for line in lines {
        let tok = tokenizer::tokenize(line);
        if tok.len() < 8 {
            continue;
        }

        let mut e = RDIIDecayEntry { uh_name: tok[0].clone(), ..RDIIDecayEntry::default() };

        let resp_upper = tokenizer::to_upper(&tok[1]);
        if resp_upper == b"SHORT" {
            e.response = 0;
        } else if resp_upper == b"MEDIUM" {
            e.response = 1;
        } else if resp_upper == b"LONG" {
            e.response = 2;
        } else {
            continue;
        }

        e.k_dep = to_double(&tok[2], 0.0);
        e.k_0 = to_double(&tok[3], 0.0);
        e.k_T = to_double(&tok[4], 0.0);
        e.T_ref = to_double(&tok[5], 0.0);
        e.theta_rec = to_double(&tok[6], 0.0);
        e.T_freeze = to_double(&tok[7], 0.0);

        if e.k_dep < 0.0 || e.k_0 < 0.0 || e.k_T < 0.0 {
            continue;
        }

        // Optional degree-day snow model: SNOW snow_T snow_ddf
        if tok.len() > 8 && tokenizer::to_upper(&tok[8]) == b"SNOW" {
            if tok.len() < 11 {
                continue; // malformed snow clause
            }
            e.snow_on = true;
            e.snow_T = to_double(&tok[9], 0.0);
            e.snow_ddf = to_double(&tok[10], 0.0);
            if e.snow_ddf < 0.0 {
                continue;
            }
        }

        ctx.rdii_decay.add(&e);
    }
}
