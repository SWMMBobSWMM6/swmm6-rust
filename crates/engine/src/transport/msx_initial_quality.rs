// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/MsxInitialQuality.cpp and
// MsxInitialQuality.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// [INITIAL_QUALITY] rows naming a reactions-component (MSX) species, mirrored into the
// component's per-element seed table (namespace openswmm::transport).

use crate::core::context::SimulationContext;
use crate::core::options::Bytes;
use crate::data::initial_quality::InitialQualityData;

/// True when [INITIAL_QUALITY] holds a resolved MSX row for (scope, element, species).
pub fn initial_quality_has_msx_row(ctx: &SimulationContext, is_link: bool, elem_idx: i32, species: i32) -> bool {
    let iq = &ctx.initial_quality;
    let want = InitialQualityData::msx_kind(species);
    for r in 0..iq.count() {
        let ur = r as usize;
        if (iq.is_link[ur] != 0) == is_link && iq.elem_idx[ur] == elem_idx && iq.kind[ur] == want {
            return true;
        }
    }
    false
}

/// Append every resolved MSX row to ReactionData::init_elem_*, unless the table already
/// holds that (scope, element, species); a held row with a different value is an error.
pub fn mirror_initial_quality_msx_rows(ctx: &mut SimulationContext) -> Vec<Bytes> {
    let mut errs: Vec<Bytes> = Vec::new();
    let rx = &mut ctx.reactions;
    let iq = &ctx.initial_quality;
    if !rx.configured {
        return errs;
    }

    for r in 0..iq.count() {
        let ur = r as usize;
        let m = InitialQualityData::msx_species(iq.kind[ur]);
        if m < 0 || m >= rx.n_species() {
            continue;
        }
        let ei = iq.elem_idx[ur];
        let link = iq.is_link[ur] != 0;
        if ei < 0 {
            continue;
        }

        let mut found = false;
        for k in 0..rx.init_elem_idx.len() {
            if (rx.init_elem_is_link[k] != 0) == link && rx.init_elem_idx[k] == ei && rx.init_elem_species[k] == m {
                if rx.init_elem_value[k] != iq.value[ur] {
                    let mut e = b"[INITIAL_QUALITY] '".to_vec();
                    e.extend_from_slice(&iq.constituent[ur]);
                    e.extend_from_slice(b"' at ");
                    e.extend_from_slice(if link { b"link '" } else { b"node '" });
                    e.extend_from_slice(&iq.elem_name[ur]);
                    e.extend_from_slice(b"' conflicts with a [REACTION_QUALITY] ");
                    e.extend_from_slice(if link { b"LINK" } else { b"NODE" });
                    e.extend_from_slice(" row in the reactions config \u{2014} keep one.".as_bytes());
                    errs.push(e);
                }
                found = true;
                break;
            }
        }
        if !found {
            rx.init_elem_is_link.push(if link { 1 } else { 0 });
            rx.init_elem_idx.push(ei);
            rx.init_elem_species.push(m);
            rx.init_elem_value.push(iq.value[ur]);
        }
    }
    errs
}
