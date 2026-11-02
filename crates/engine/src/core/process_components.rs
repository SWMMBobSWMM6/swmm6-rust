// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/plugins/ProcessComponentRegistry.hpp and
// ProcessComponentRegistry.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit
// 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// [PROCESS_COMPONENTS]: the component registry, the component config-file reader and
// resolve_process_components (called from SWMMEngine::open after the .inp read).
//
// Differences of form only:
//   * Upstream's registry is a process singleton that SWMMEngine's constructor fills once
//     (std::call_once) with the four implemented components, over the six planned
//     placeholders. Nothing else ever registers, so here instance() returns that filled
//     registry, built once (OnceLock). entries_ is a std::map; a BTreeMap keeps the same
//     byte order for known_ids().
//   * Hooks are plain function pointers. Save hooks are not part of open() and are left out.
//   * resolve_process_components hands the hook a copy of the spec (with the resolved path
//     already set), since the hook also receives the context that owns the spec.
//   * std::ifstream on a directory opens and then reads nothing; File::open does the same
//     on Linux, and a failed read is treated as end of file, as getline does.

use std::collections::BTreeMap;
use std::io::Read;
use std::sync::OnceLock;

use crate::core::context::{ProcessComponentSpec, SimulationContext};
use crate::core::fs_path::FsPath;
use crate::core::options::Bytes;

/// One parsed component config file (or the embedded [REACTION_*] sections).
#[derive(Debug, Clone, Default)]
pub struct ComponentConfigSections {
    /// resolved absolute/effective path
    pub source_path: Bytes,
    pub sections: Vec<(Bytes, Vec<Bytes>)>,
}

impl ComponentConfigSections {
    pub fn find(&self, upper_tag: &[u8]) -> Option<&Vec<Bytes>> {
        for s in &self.sections {
            if s.0 == upper_tag {
                return Some(&s.1);
            }
        }
        None
    }
}

/// ComponentConfigApply
pub type ComponentConfigApply = fn(&mut SimulationContext, &ProcessComponentSpec, &ComponentConfigSections, &mut Vec<Bytes>);

#[derive(Clone)]
pub struct ProcessComponentEntry {
    pub description: Bytes,
    /// non-empty => planned, not yet implemented
    pub pending_phase: Bytes,
    /// None for planned entries
    pub apply: Option<ComponentConfigApply>,
}

// ===========================================================================
// Config-file reader (IO2)
// ===========================================================================

/// toupper over each byte in the C locale.
fn upper(s: &[u8]) -> Bytes {
    s.iter().map(|c| c.to_ascii_uppercase()).collect()
}

/// The first ';' begins a comment even inside quotes; trailing space/tab/CR removed.
fn strip_comment(line: &[u8]) -> Bytes {
    let mut s: Bytes = match line.iter().position(|&c| c == b';') {
        Some(sc) => line[..sc].to_vec(),
        None => line.to_vec(),
    };
    while let Some(&b) = s.last() {
        if b == b' ' || b == b'\t' || b == b'\r' {
            s.pop();
        } else {
            break;
        }
    }
    s
}

fn find_from(s: &[u8], c: u8, from: usize) -> Option<usize> {
    if from > s.len() {
        return None;
    }
    s[from..].iter().position(|&x| x == c).map(|p| p + from)
}

/// A separator, or an exact .so/.dylib/.dll suffix after a nonempty name.
fn looks_like_library(id: &[u8]) -> bool {
    if id.contains(&b'/') || id.contains(&b'\\') {
        return true;
    }
    for e in [&b".so"[..], &b".dylib"[..], &b".dll"[..]] {
        if id.len() > e.len() && &id[id.len() - e.len()..] == e {
            return true;
        }
    }
    false
}

/// std::getline over a byte buffer: split on '\n'; a final line without one still counts.
fn getline_split(data: &[u8]) -> Vec<&[u8]> {
    let mut lines: Vec<&[u8]> = data.split(|&c| c == b'\n').collect();
    if data.is_empty() || data.last() == Some(&b'\n') {
        lines.pop();
    }
    lines
}

/// Append parsed bracketed sections from `data`. Empty on accepted input.
pub fn parse_component_config_stream(data: &[u8], label: &[u8], out: &mut ComponentConfigSections) -> Bytes {
    let mut have_current = false;
    for raw in getline_split(data) {
        let line = strip_comment(raw);
        if line.is_empty() {
            continue;
        }

        // Section header?
        let lb = line.iter().position(|&c| c == b'[');
        if let Some(lb) = lb {
            let first_not_blank = line.iter().position(|&c| c != b' ' && c != b'\t');
            if find_from(&line, b']', lb).is_some() && first_not_blank == Some(lb) {
                let rb = find_from(&line, b']', lb).unwrap();
                let tag = upper(&line[lb + 1..rb]);
                if tag == b"PROCESS_COMPONENTS" {
                    let mut m = b"Component config '".to_vec();
                    m.extend_from_slice(label);
                    m.extend_from_slice(
                        "' contains a nested [PROCESS_COMPONENTS] section \u{2014} component files cannot register components (no recursion; TRANSPORT_IO_PLUGIN_CONFIG_PLAN \u{a7}3.3).".as_bytes(),
                    );
                    return m;
                }
                out.sections.push((tag, Vec::new()));
                have_current = true;
                continue;
            }
        }
        if !have_current {
            let mut m = b"Component config '".to_vec();
            m.extend_from_slice(label);
            m.extend_from_slice(b"' has content before its first [SECTION] header.");
            return m;
        }
        out.sections.last_mut().unwrap().1.push(line);
    }
    Vec::new()
}

pub fn parse_component_config_text(text: &[u8], label: &[u8], out: &mut ComponentConfigSections) -> Bytes {
    out.sections.clear();
    out.source_path = label.to_vec();
    parse_component_config_stream(text, label, out)
}

pub fn read_component_config(path: &[u8], base_dir: &[u8], out: &mut ComponentConfigSections) -> Bytes {
    out.sections.clear();

    let mut p = FsPath::new(path);
    if !p.is_absolute() && !base_dir.is_empty() {
        p = FsPath::new(base_dir).join(&p);
    }
    out.source_path = p.native().to_vec();

    let file = std::fs::File::open(crate::core::file_io::path_from_bytes(p.native()));
    let mut file = match file {
        Ok(f) => f,
        Err(_) => {
            let mut m = b"Process component config file not found or unreadable: '".to_vec();
            m.extend_from_slice(&out.source_path);
            m.extend_from_slice(b"'.");
            return m;
        }
    };
    let mut data = Vec::new();
    if file.read_to_end(&mut data).is_err() {
        data.clear(); // a directory: the first getline fails
    }
    let label = out.source_path.clone();
    parse_component_config_stream(&data, &label, out)
}

// ===========================================================================
// Registry (IO1)
// ===========================================================================

pub struct ProcessComponentRegistry {
    entries_: BTreeMap<Bytes, ProcessComponentEntry>,
}

impl ProcessComponentRegistry {
    /// The planned built-in ids (TRANSPORT_IO_PLUGIN_CONFIG_PLAN §2).
    fn new() -> Self {
        let planned: [(&str, &str, &str); 6] = [
            ("org.hydrocouple.openswmm.reactions", "Multispecies reaction system (EPANET-MSX conventions)", "R1 (reactions plan)"),
            ("org.hydrocouple.openswmm.transport.ard", "Eulerian ARD transport engine configuration", "E5 (Eulerian ARD plan)"),
            ("org.hydrocouple.openswmm.transport.lard", "Lagrangian (LARD) quality engine", "T5 (LARD plan)"),
            ("org.hydrocouple.openswmm.heat", "1D heat transport flux modules", "H1 (heat plan)"),
            ("org.hydrocouple.openswmm.waterage", "Water age tracking coordinator", "A1 (water age plan)"),
            ("org.hydrocouple.openswmm.integrated2d", "Unified 2D surface + two-zone groundwater component", "S1/G1 (2D + GW plans)"),
        ];
        let mut entries_ = BTreeMap::new();
        for (id, desc, phase) in planned {
            entries_.insert(
                id.as_bytes().to_vec(),
                ProcessComponentEntry { description: desc.as_bytes().to_vec(), pending_phase: phase.as_bytes().to_vec(), apply: None },
            );
        }
        ProcessComponentRegistry { entries_ }
    }

    /// The process registry as SWMMEngine's constructor leaves it: the planned ids, with
    /// the four implemented components registered over them in the constructor's order
    /// (reactions, ARD, water age, heat).
    pub fn instance() -> &'static ProcessComponentRegistry {
        static REG: OnceLock<ProcessComponentRegistry> = OnceLock::new();
        REG.get_or_init(|| {
            let mut reg = ProcessComponentRegistry::new();
            crate::transport::reactions_component::register_reactions_component(&mut reg);
            crate::transport::ard_config::register_ard_component(&mut reg);
            crate::transport::water_age_component::register_water_age_component(&mut reg);
            crate::transport::heat_component::register_heat_component(&mut reg);
            reg
        })
    }

    pub fn register_component(&mut self, id: &[u8], description: &[u8], apply: ComponentConfigApply) {
        self.entries_.insert(
            id.to_vec(),
            ProcessComponentEntry { description: description.to_vec(), pending_phase: Vec::new(), apply: Some(apply) },
        );
    }

    pub fn find(&self, id: &[u8]) -> Option<&ProcessComponentEntry> {
        self.entries_.get(id)
    }

    pub fn known_ids(&self) -> Vec<Bytes> {
        self.entries_.keys().cloned().collect()
    }
}

// ===========================================================================
// Resolution (IO1) — called from SWMMEngine::open after the input read
// ===========================================================================

pub fn resolve_process_components(ctx: &mut SimulationContext, base_dir: &[u8]) -> Vec<Bytes> {
    let mut errors: Vec<Bytes> = Vec::new();
    let reg = ProcessComponentRegistry::instance();

    // Duplicate registrations are refused.
    {
        let mut seen: Vec<Bytes> = Vec::new();
        for spec in &ctx.process_component_specs {
            for s0 in &seen {
                if *s0 == spec.id {
                    let mut m = b"[PROCESS_COMPONENTS] duplicate registration of '".to_vec();
                    m.extend_from_slice(&spec.id);
                    m.extend_from_slice("' \u{2014} each component id may appear once.".as_bytes());
                    errors.push(m);
                    break;
                }
            }
            seen.push(spec.id.clone());
        }
        if !errors.is_empty() {
            return errors;
        }
    }

    for k in 0..ctx.process_component_specs.len() {
        let id = ctx.process_component_specs[k].id.clone();
        if looks_like_library(&id) {
            let mut m = b"[PROCESS_COMPONENTS] '".to_vec();
            m.extend_from_slice(&id);
            m.extend_from_slice(
                "': library-loaded HydroCouple components are not available yet (arrives with plan phase HC2 \u{2014} IMPLEMENTATION_ROADMAP.md).".as_bytes(),
            );
            errors.push(m);
            continue;
        }
        let entry = match reg.find(&id) {
            Some(e) => e,
            None => {
                let mut msg = b"[PROCESS_COMPONENTS] unknown component id '".to_vec();
                msg.extend_from_slice(&id);
                msg.extend_from_slice(b"'. Known ids:");
                for kid in reg.known_ids() {
                    msg.push(b' ');
                    msg.extend_from_slice(&kid);
                }
                errors.push(msg);
                continue;
            }
        };
        if !entry.pending_phase.is_empty() {
            let mut m = b"[PROCESS_COMPONENTS] '".to_vec();
            m.extend_from_slice(&id);
            m.extend_from_slice(b"' (");
            m.extend_from_slice(&entry.description);
            m.extend_from_slice(") is recognized but not yet implemented \u{2014} arrives with plan phase ".as_bytes());
            m.extend_from_slice(&entry.pending_phase);
            m.push(b'.');
            errors.push(m);
            continue;
        }
        if ctx.process_component_specs[k].config_path.is_empty() {
            let mut m = b"[PROCESS_COMPONENTS] '".to_vec();
            m.extend_from_slice(&id);
            m.extend_from_slice("': missing required config=\"\u{2026}\" argument.".as_bytes());
            errors.push(m);
            continue;
        }
        let mut config = ComponentConfigSections::default();
        let cp = ctx.process_component_specs[k].config_path.clone();
        let err = read_component_config(&cp, base_dir, &mut config);
        if !err.is_empty() {
            errors.push(err);
            continue;
        }
        // IO3: remember where the config was actually read from.
        ctx.process_component_specs[k].resolved_config_path = config.source_path.clone();
        let spec = ctx.process_component_specs[k].clone();
        (entry.apply.unwrap())(ctx, &spec, &config, &mut errors);
    }
    errors
}
