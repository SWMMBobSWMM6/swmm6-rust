// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/InitialQualityData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// SoA store for [INITIAL_QUALITY] per-element initial concentrations. Upstream's add()
// has defaulted trailing parameters (ei = -1, k = kKindUnresolved, fromFile = false);
// Rust has none, so add() takes the four leading ones and applies those defaults, and
// add_with() takes all seven. The kKind* consts are K_KIND_*; msxKind/msxSpecies are
// msx_kind/msx_species.

use crate::core::options::Bytes;

#[derive(Debug, Clone, Default)]
pub struct InitialQualityData {
    /// 0 = NODE row, 1 = LINK row
    pub is_link: Vec<u8>,
    /// Raw element name (post-parse re-resolution)
    pub elem_name: Vec<Bytes>,
    /// Resolved node/link index (-1 until resolved)
    pub elem_idx: Vec<i32>,
    /// Raw constituent name (writer round-trip)
    pub constituent: Vec<Bytes>,
    /// >=0 pollutant idx; K_KIND_* otherwise
    pub kind: Vec<i32>,
    /// Raw user-units value
    pub value: Vec<f64>,
    /// U2: 1 when the row came from the `[INITIAL_QUALITY] FILE <csv>` sidecar.
    pub from_file: Vec<u8>,
    /// U2: `[INITIAL_QUALITY] FILE <path>` as written (relative to the .inp); empty when none.
    pub file: Bytes,
}

impl InitialQualityData {
    /// __WATER_AGE__ (hours)
    pub const K_KIND_WATER_AGE: i32 = -1;
    /// __TEMPERATURE__ (degC)
    pub const K_KIND_TEMPERATURE: i32 = -2;
    /// Not yet classified
    pub const K_KIND_UNRESOLVED: i32 = -3;
    /// U2: MSX species m is encoded as kind = K_KIND_MSX_FIRST - m.
    pub const K_KIND_MSX_FIRST: i32 = -100;

    pub const fn msx_kind(species: i32) -> i32 {
        Self::K_KIND_MSX_FIRST - species
    }

    pub const fn msx_species(k: i32) -> i32 {
        if k <= Self::K_KIND_MSX_FIRST {
            Self::K_KIND_MSX_FIRST - k
        } else {
            -1
        }
    }

    pub fn count(&self) -> i32 {
        self.is_link.len() as i32
    }

    /// add(link, name, cons, v) with upstream's defaults ei = -1, k = K_KIND_UNRESOLVED,
    /// fromFile = false.
    pub fn add(&mut self, link: bool, name: &[u8], cons: &[u8], v: f64) {
        self.add_with(link, name, cons, v, -1, Self::K_KIND_UNRESOLVED, false);
    }

    /// Append a row without resolving its names.
    #[allow(clippy::too_many_arguments)]
    pub fn add_with(
        &mut self,
        link: bool,
        name: &[u8],
        cons: &[u8],
        v: f64,
        ei: i32,
        k: i32,
        from_file: bool,
    ) {
        self.is_link.push(if link { 1 } else { 0 });
        self.elem_name.push(name.to_vec());
        self.elem_idx.push(ei);
        self.constituent.push(cons.to_vec());
        self.kind.push(k);
        self.value.push(v);
        self.from_file.push(if from_file { 1 } else { 0 });
    }

    /// Remove row idx; negative or >= count() is a no-op. from_file is erased only when
    /// it reaches that row.
    pub fn erase(&mut self, idx: i32) {
        if idx < 0 || idx >= self.count() {
            return;
        }
        let u = idx as usize;
        self.is_link.remove(u);
        self.elem_name.remove(u);
        self.elem_idx.remove(u);
        self.constituent.remove(u);
        self.kind.remove(u);
        self.value.remove(u);
        if u < self.from_file.len() {
            self.from_file.remove(u);
        }
    }
}
