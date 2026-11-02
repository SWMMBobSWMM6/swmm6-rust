// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/InputReader.hpp and InputReader.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Reads an .inp in one pass: a section's lines are gathered and handed to its handler
// when the next header (or the end of the file) arrives, in file order. Rows a handler
// could not resolve (an object declared further down) are replayed once at the end.
//
// Lines are split as std::getline splits them: on '\n'; a final line without one still
// counts; one trailing '\r' is removed.

use crate::core::context::{EngineState, SimulationContext};
use crate::core::error_codes::{format_error, ERR_NAME};
use crate::core::options::Bytes;
use crate::input::section_registry::SectionRegistry;
use crate::input::tokenizer;

/// Legacy SectWords (keywords.c). legacy findmatch() accepts a header that merely BEGINS
/// with one of them ("[PROFILES]" is "[PROFILE"), so such a header draws no
/// unknown-section warning even when there is no handler for it.
const K_LEGACY_SECTION_WORDS: &[&[u8]] = &[
    b"TITLE", b"OPTION", b"FILE", b"RAINGAGE", b"TEMPERATURE", b"EVAP", b"SUBCATCHMENT",
    b"SUBAREA", b"INFIL", b"AQUIFER", b"GROUNDWATER", b"SNOWPACK", b"JUNC", b"OUTFALL",
    b"STORAGE", b"DIVIDER", b"CONDUIT", b"PUMP", b"ORIFICE", b"WEIR", b"OUTLET", b"XSECT",
    b"TRANSECT", b"LOSS", b"CONTROL", b"POLLUT", b"LANDUSE", b"BUILDUP", b"WASHOFF",
    b"COVERAGE", b"INFLOW", b"DWF", b"PATTERN", b"RDII", b"HYDROGRAPH", b"LOADING",
    b"TREATMENT", b"CURVE", b"TIMESERIES", b"REPORT", b"COORDINATE", b"VERTICES",
    b"POLYGON", b"LABEL", b"SYMBOL", b"BACKDROP", b"TAG", b"PROFILE", b"MAP",
    b"LID_CONTROL", b"LID_USAGE", b"GWF", b"ADJUSTMENT", b"EVENT", b"STREET",
    b"INLET_USAGE", b"INLET",
];

fn legacy_knows_section(tag: &[u8]) -> bool {
    K_LEGACY_SECTION_WORDS.iter().any(|w| tag.starts_with(w))
}

pub struct InputReader<'r> {
    registry: &'r SectionRegistry,
    lines_read: i32,
    skipped_sections: Vec<Bytes>,
}

impl<'r> InputReader<'r> {
    pub fn new(registry: &'r SectionRegistry) -> Self {
        InputReader { registry, lines_read: 0, skipped_sections: Vec::new() }
    }

    pub fn lines_read(&self) -> i32 {
        self.lines_read
    }

    pub fn skipped_sections(&self) -> &[Bytes] {
        &self.skipped_sections
    }

    pub fn read(&mut self, path: &[u8], ctx: &mut SimulationContext) -> bool {
        let p = crate::core::file_io::path_from_bytes(path);
        let data = match std::fs::read(&p) {
            Ok(d) => d,
            Err(_) => {
                ctx.error_code = 2; // public SWMM_ERR_INPFILE
                let mut m = b"InputReader: cannot open file '".to_vec();
                m.extend_from_slice(path);
                m.push(b'\'');
                ctx.error_message = m;
                ctx.state = EngineState::ErrorState;
                return false;
            }
        };
        ctx.inp_file_path = path.to_vec();
        self.read_stream(&data, ctx)
    }

    pub fn read_stream(&mut self, data: &[u8], ctx: &mut SimulationContext) -> bool {
        self.lines_read = 0;
        self.skipped_sections.clear();
        ctx.unknown_section_headers.clear();
        // An .inp reports no subcatchments, nodes or links unless [REPORT] names them.
        ctx.options.rpt_subcatchments = 0;
        ctx.options.rpt_nodes = 0;
        ctx.options.rpt_links = 0;

        let mut current_tag: Bytes = Vec::new();
        let mut section_lines: Vec<Bytes> = Vec::new();

        let mut lines: Vec<&[u8]> = data.split(|&c| c == b'\n').collect();
        if data.is_empty() || data.last() == Some(&b'\n') {
            lines.pop();
        }
        for raw in lines {
            self.lines_read += 1;
            let mut raw_line = raw;
            if raw_line.last() == Some(&b'\r') {
                raw_line = &raw_line[..raw_line.len() - 1];
            }
            let trimmed_raw = tokenizer::trim(raw_line);
            if trimmed_raw.is_empty() {
                continue;
            }
            // A clean "[NAME]" (optionally followed by a comment) is a header; a line that
            // starts with '[' but carries more after the ']' is data.
            if trimmed_raw[0] == b'[' {
                let header = tokenizer::trim(tokenizer::strip_comment(trimmed_raw));
                if !header.is_empty() && *header.last().unwrap() == b']' {
                    if !current_tag.is_empty() {
                        self.dispatch_section(&current_tag, &section_lines, ctx);
                        section_lines.clear();
                    }
                    current_tag = Self::parse_section_header(header);
                    if !current_tag.is_empty() && !self.registry.has(&current_tag) && !legacy_knows_section(&current_tag) {
                        // legacy names the first token of the line, as written
                        let end = trimmed_raw.iter().position(|&c| c == b' ' || c == b'\t').unwrap_or(trimmed_raw.len());
                        ctx.unknown_section_headers.push((trimmed_raw[..end].to_vec(), self.lines_read as i64));
                    }
                    continue;
                }
            }
            // A comment line is kept inside a section, verbatim (trimmed).
            if trimmed_raw[0] == b';' {
                if !current_tag.is_empty() {
                    section_lines.push(trimmed_raw.to_vec());
                }
                continue;
            }
            // [TITLE] keeps the raw line, leading blanks and ';' included.
            if current_tag == b"TITLE" {
                section_lines.push(raw_line.to_vec());
                continue;
            }
            let stripped = tokenizer::strip_comment(trimmed_raw);
            let trimmed = tokenizer::trim(stripped);
            if !trimmed.is_empty() && !current_tag.is_empty() {
                section_lines.push(trimmed.to_vec());
            }
            // lines before the first header are ignored
        }
        if !current_tag.is_empty() {
            self.dispatch_section(&current_tag, &section_lines, ctx);
            section_lines.clear();
        }
        self.replay_deferred_rows(ctx);
        ctx.error_code == 0
    }

    fn replay_deferred_rows(&mut self, ctx: &mut SimulationContext) {
        if ctx.deferred_section_rows.is_empty() {
            return;
        }
        let pending = std::mem::take(&mut ctx.deferred_section_rows);
        // group by tag, file order within each tag
        let mut tags: Vec<Bytes> = Vec::new();
        for (t, _) in &pending {
            if !tags.contains(t) {
                tags.push(t.clone());
            }
        }
        for tag in &tags {
            let lines: Vec<Bytes> = pending.iter().filter(|(t, _)| t == tag).map(|(_, l)| l.clone()).collect();
            self.dispatch_section(tag, &lines, ctx);
        }
        // anything stashed again names an object that does not exist: ERROR 209
        let rows = std::mem::take(&mut ctx.deferred_section_rows);
        for (_, line) in &rows {
            let tok = tokenizer::tokenize(line);
            let name: &[u8] = if tok.is_empty() { line } else { &tok[0] };
            ctx.errors.push(format_error(ERR_NAME, name));
        }
        ctx.deferred_section_rows.clear();
    }

    fn parse_section_header(line: &[u8]) -> Bytes {
        if line.len() < 3 || line[0] != b'[' {
            return Vec::new();
        }
        let close = match line.iter().position(|&c| c == b']') {
            Some(c) => c,
            None => return Vec::new(),
        };
        let name = tokenizer::trim(&line[1..close]);
        if name.is_empty() {
            return Vec::new();
        }
        tokenizer::to_upper(name)
    }

    fn dispatch_section(&mut self, tag: &[u8], lines: &[Bytes], ctx: &mut SimulationContext) {
        if self.registry.has(tag) {
            self.registry.dispatch(tag, ctx, lines);
        } else if !legacy_knows_section(tag) {
            if !self.skipped_sections.iter().any(|s| s.as_slice() == tag) {
                self.skipped_sections.push(tag.to_vec());
            }
            if ctx.warning_code == 0 {
                ctx.warning_code = 100; // SWMM_WARN_UNKNOWN_SECTION
            }
        }
    }
}
