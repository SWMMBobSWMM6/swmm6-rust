// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/LidLayerSpeciesData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Phase A4: per-(LID unit, layer, species) transported state. kLayerCount is
// K_LAYER_COUNT; ensureSized is ensure_sized. The borrowed `const std::vector<int>&`
// offsets are a slice, copied into owned storage as upstream's copy-assignment does.
// Index arithmetic is size_t in upstream and usize here, in the same order.

/// A4: layer index within a LID unit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum LidLayer {
    Surface = 0,
    Pavement = 1,
    Soil = 2,
    Storage = 3,
    Count = 4,
}

/// Species rows carried on the LID layers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum LidSpecies {
    Age = 0,
    /// H5b.
    Temperature = 1,
    Count = 2,
}

/// Per-(unit, layer, species) state for every LID unit in the model.
#[derive(Debug, Clone, Default)]
pub struct LidLayerSpeciesState {
    pub n_units: i32,
    pub n_species: i32,
    /// [type] flat index of that group's first unit; size = numGroups() + 1.
    pub group_offset: Vec<i32>,
    /// [(unit * K_LAYER_COUNT + layer) * n_species + species]
    pub value: Vec<f64>,
    /// [unit * K_LAYER_COUNT + layer] — stored water at the END of the previous step (ft).
    pub vol_prev: Vec<f64>,
    /// [unit * n_species + species] — what leaves through the underdrain.
    pub drain_value: Vec<f64>,
    /// [unit * n_species + species] — the value of the water arriving at the unit.
    pub inflow_value: Vec<f64>,
}

// n_units = n_species = 0 are Rust's defaults, so Default is derived.

impl LidLayerSpeciesState {
    pub const K_LAYER_COUNT: i32 = LidLayer::Count as i32;

    /// n_units > 0 && n_species > 0; does not validate arrays or offsets.
    pub fn active(&self) -> bool {
        self.n_units > 0 && self.n_species > 0
    }

    /// (unit*K_LAYER_COUNT+layer)*n_species+species, unchecked.
    pub fn layer_index(&self, unit: i32, layer: LidLayer, species: i32) -> usize {
        ((unit as usize) * (Self::K_LAYER_COUNT as usize) + (layer as i32 as usize))
            * (self.n_species as usize)
            + (species as usize)
    }

    /// Replace LID layer state with a zeroed shape and copied group offsets.
    pub fn resize(&mut self, units: i32, species: i32, offsets: &[i32]) {
        self.n_units = units;
        self.n_species = species;
        self.group_offset = offsets.to_vec();
        let n_layers = (units as usize) * (Self::K_LAYER_COUNT as usize);
        self.value.clear();
        self.value.resize(n_layers * (species as usize), 0.0);
        self.vol_prev.clear();
        self.vol_prev.resize(n_layers, 0.0);
        self.drain_value.clear();
        self.drain_value.resize((units as usize) * (species as usize), 0.0);
        self.inflow_value.clear();
        self.inflow_value.resize((units as usize) * (species as usize), 0.0);
    }

    /// Resize only when the counts or the value length differ.
    pub fn ensure_sized(&mut self, units: i32, species: i32, offsets: &[i32]) {
        let want = (units as usize) * (Self::K_LAYER_COUNT as usize) * (species as usize);
        if self.n_units == units && self.n_species == species && self.value.len() == want {
            return;
        }
        self.resize(units, species, offsets);
    }

    pub fn clear(&mut self) {
        *self = LidLayerSpeciesState::default();
    }
}
