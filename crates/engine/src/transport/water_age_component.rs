// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/components/WaterAgeModule/
// WaterAgeComponent.hpp and WaterAgeComponent.cpp (HydroCouple/openswmm.engine, branch
// swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache
// License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The `org.hydrocouple.openswmm.waterage` component's config apply hook (model.age,
// [WATER_AGE_SOURCES]). The serializer is not part of open() and is not translated here.

use crate::core::context::{ProcessComponentSpec, SimulationContext};
use crate::core::options::Bytes;
use crate::core::process_components::{ComponentConfigSections, ProcessComponentRegistry};
use crate::data::water_age::{WaterAgeConfigData, WaterAgeSource};
use crate::transport::ard_config::{parse_value, tokenize_ws, upper};

const K_AGE_ID: &str = "org.hydrocouple.openswmm.waterage";

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

fn source_of(u: &[u8]) -> i32 {
    match u {
        b"RAINFALL" => WaterAgeSource::Rainfall as i32,
        b"DWF" => WaterAgeSource::Dwf as i32,
        b"GW" => WaterAgeSource::Gw as i32,
        b"RDII" => WaterAgeSource::Rdii as i32,
        b"EXTERNAL_INFLOW" => WaterAgeSource::ExternalInflow as i32,
        b"IFACE" => WaterAgeSource::Iface as i32,
        b"INITIAL_STATE" => WaterAgeSource::InitialState as i32,
        _ => -1,
    }
}

fn apply_water_age_sections(ctx: &mut SimulationContext, _spec: &ProcessComponentSpec, config: &ComponentConfigSections, errors: &mut Vec<Bytes>) {
    ctx.water_age_config = WaterAgeConfigData::default(); // reopen hygiene
    let errors_before = errors.len();

    for sec in &config.sections {
        if sec.0 != b"WATER_AGE_SOURCES" {
            errors.push(cat(&[b"model.age: unknown section [", &sec.0, b"] (A1a recognizes: [WATER_AGE_SOURCES])."]));
            continue;
        }
        for line in &sec.1 {
            let toks = tokenize_ws(line);
            if toks.len() < 3 {
                errors.push(cat(&[b"[WATER_AGE_SOURCES] expects '<source> <scope> [name] <hours>': '", line, b"'."]));
                continue;
            }
            let src = source_of(&upper(&toks[0]));
            if src < 0 {
                errors.push(cat(&[
                    b"[WATER_AGE_SOURCES] unknown source '",
                    &toks[0],
                    b"' (RAINFALL, DWF, GW, RDII, EXTERNAL_INFLOW, IFACE, INITIAL_STATE).",
                ]));
                continue;
            }
            let scope = upper(&toks[1]);

            if scope == b"SUBCATCH" {
                errors.push(b"[WATER_AGE_SOURCES] SUBCATCH scope arrives with plan phase A3 (watershed age states).".to_vec());
                continue;
            }
            if scope == b"EDGE_BC" {
                errors.push(b"[WATER_AGE_SOURCES] EDGE_BC scope arrives with phase T6 (2D transport).".to_vec());
                continue;
            }
            let has_name = scope == b"NODE";
            let vpos: usize = if has_name { 3 } else { 2 };
            if toks.len() <= vpos {
                errors.push(cat(&[b"[WATER_AGE_SOURCES] malformed row: '", line, b"'."]));
                continue;
            }
            let vtok = &toks[vpos];
            if upper(vtok).starts_with(b"TIMESERIES") {
                errors.push(
                    "[WATER_AGE_SOURCES] TIMESERIES ages are prescribed as [INFLOWS] rows naming __WATER_AGE__ (amendment D-Y4 / Z1) \u{2014} this table takes constant hours."
                        .as_bytes()
                        .to_vec(),
                );
                continue;
            }
            if toks.len() != vpos + 1 {
                errors.push(cat(&[b"[WATER_AGE_SOURCES] malformed row: '", line, b"'."]));
                continue;
            }
            let mut hours = 0.0;
            if !parse_value(vtok, &mut hours) {
                errors.push(cat(&[b"[WATER_AGE_SOURCES] '", &toks[0], b"': '", vtok, b"' is not an age in hours."]));
                continue;
            }
            if hours < 0.0 {
                ctx.warnings.push(cat(&[
                    b"[WATER_AGE_SOURCES] '",
                    &toks[0],
                    b"' is negative (",
                    vtok,
                    b" h): the pathway EXTRACTS age-volume (water reads younger). Extraction clamps so age never goes below zero (D-NS1).",
                ]));
            }
            let seconds = hours * 3600.0;

            if scope == b"GLOBAL" {
                ctx.water_age_config.global_age[src as usize] = seconds;
            } else if scope == b"NODE" {
                if src != WaterAgeSource::Dwf as i32 && src != WaterAgeSource::ExternalInflow as i32 {
                    errors.push(cat(&[b"[WATER_AGE_SOURCES] NODE scope applies to DWF and EXTERNAL_INFLOW in A1a; '", &toks[0], b"' takes GLOBAL."]));
                    continue;
                }
                let nd = ctx.node_names.find(&toks[2]);
                if nd < 0 {
                    errors.push(cat(&[b"[WATER_AGE_SOURCES] unknown node '", &toks[2], b"'."]));
                    continue;
                }
                let mut dup = false;
                for i in 0..ctx.water_age_config.node_over_source.len() {
                    if ctx.water_age_config.node_over_source[i] == src && ctx.water_age_config.node_over_node[i] == nd {
                        errors.push(cat(&[b"[WATER_AGE_SOURCES] duplicate NODE row for '", &toks[0], b"' at '", &toks[2], b"'."]));
                        dup = true;
                        break;
                    }
                }
                if dup {
                    continue;
                }
                ctx.water_age_config.node_over_source.push(src);
                ctx.water_age_config.node_over_node.push(nd);
                ctx.water_age_config.node_over_age.push(seconds);
            } else {
                errors.push(cat(&[b"[WATER_AGE_SOURCES] unknown scope '", &toks[1], b"' (GLOBAL or NODE in A1a)."]));
            }
        }
    }

    if errors.len() != errors_before {
        ctx.water_age_config = WaterAgeConfigData::default(); // never half-apply
        return;
    }
    ctx.water_age_config.configured = true;

    if !ctx.options.water_age {
        ctx.warnings.push(
            "A waterage component is configured but [OPTIONS] WATER_AGE is OFF \u{2014} no age is tracked this simulation. Set WATER_AGE ON."
                .as_bytes()
                .to_vec(),
        );
    }
}

pub fn register_water_age_component(reg: &mut ProcessComponentRegistry) {
    reg.register_component(
        K_AGE_ID.as_bytes(),
        b"Water age tracking coordinator (A1a: [WATER_AGE_SOURCES] GLOBAL + NODE constant ages)",
        apply_water_age_sections,
    );
}
