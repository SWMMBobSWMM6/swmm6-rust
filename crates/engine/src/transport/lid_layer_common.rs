// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/components/LidLayerCommon.hpp and
// LidLayerCommon.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The LID layer adjacency shared by the LID species, age and heat modules (namespace
// transport): the layer selectors, donorsFor, layerVolumes and buildOffsets. Complete.
//
// Differences of form: `double (&v)[kNL]` is `&mut [f64; K_NL]`.

use crate::data::lid_layer_species::{LidLayer, LidLayerSpeciesState};
use crate::hydrology::lid::{LIDGroupSoA, LIDSolver, LIDType};

pub const K_NL: usize = LidLayerSpeciesState::K_LAYER_COUNT as usize;

/// -1 means "external to the unit"; -2 means "this type has no such layer".
pub const K_EXTERNAL: i32 = -1;
pub const K_ABSENT: i32 = -2;

pub const K_SURF: i32 = LidLayer::Surface as i32;
pub const K_PAVE: i32 = LidLayer::Pavement as i32;
pub const K_SOIL: i32 = LidLayer::Soil as i32;
pub const K_STOR: i32 = LidLayer::Storage as i32;

/// Where each layer's water comes from (layer index, K_EXTERNAL or K_ABSENT).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Donors {
    pub surface: i32,
    pub pavement: i32,
    pub soil: i32,
    pub storage: i32,
}

fn d(surface: i32, pavement: i32, soil: i32, storage: i32) -> Donors {
    Donors { surface, pavement, soil, storage }
}

pub fn donors_for(t: LIDType, has_soil: bool) -> Donors {
    match t {
        LIDType::BioCell | LIDType::RainGarden => d(K_EXTERNAL, K_ABSENT, K_SURF, K_SOIL),
        LIDType::GreenRoof => d(K_EXTERNAL, K_ABSENT, K_SURF, K_SOIL),
        LIDType::InfilTrench => d(K_EXTERNAL, K_ABSENT, K_ABSENT, K_SURF),
        LIDType::PermPavement => d(K_EXTERNAL, K_SURF, K_PAVE, if has_soil { K_SOIL } else { K_PAVE }),
        LIDType::RainBarrel => d(K_ABSENT, K_ABSENT, K_ABSENT, K_EXTERNAL),
        LIDType::VegSwale | LIDType::RoofDiscon => d(K_EXTERNAL, K_ABSENT, K_ABSENT, K_ABSENT),
    }
}

/// Void-weighted water held in each layer, as a depth per unit area (ft).
pub fn layer_volumes(g: &LIDGroupSoA, ui: usize, v: &mut [f64; K_NL]) {
    v[K_SURF as usize] = g.surf_depth[ui] * g.surf_void_frac[ui];
    v[K_PAVE as usize] = g.pave_depth[ui] * g.pave_void[ui] * (1.0 - g.pave_imperv_frac[ui]);
    v[K_SOIL as usize] = g.soil_moist[ui] * g.soil_thick[ui];
    v[K_STOR as usize] = g.stor_depth[ui] * g.stor_void[ui];
}

/// Flat index of each per-type group's first unit; `last() == n_units`.
pub fn build_offsets(solver: &LIDSolver) -> Vec<i32> {
    let mut off: Vec<i32> = Vec::new();
    off.reserve(solver.num_groups() as usize + 1);
    let mut running = 0i32;
    for t in 0..solver.num_groups() {
        off.push(running);
        running += solver.group(t).count;
    }
    off.push(running);
    off
}
