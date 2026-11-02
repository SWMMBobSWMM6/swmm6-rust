// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/HeatOverrideData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The per-element heat override tokens (plan PE): the element kind and index a flux is
// evaluated for, the scope and attribute of one override row, and the row itself. The
// dense storage the rows resolve into is HeatOverrideData in heat.rs, as upstream keeps
// it in HeatData.hpp. Method names are snake_case (isSedimentAttr is is_sediment_attr).
// A Rust enum cannot hold an out-of-range ordinal, so is_sediment_attr only ever sees
// declared attributes; upstream's predicate is a plain threshold and agrees on all of them.

use crate::core::options::Bytes;

/// What kind of element a flux is being evaluated for (D-PE1). ARD cells and LARD
/// parcels resolve to their parent LINK.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum HeatElemKind {
    Node = 0,
    Link = 1,
    Subcatch = 2,
    Lid = 3,
    /// Surface mesh cell; CELL2D and MESH_TAG overrides are independent of the network
    /// LINK/NODE/TAG namespaces.
    Cell2d = 4,
}

/// Which element a flux evaluator is being called for. Defaults to an INVALID index
/// (LINK, -1) deliberately, so a binding that forgets the element gets the global config.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct HeatElement {
    pub kind: HeatElemKind,
    pub index: i32,
}

impl Default for HeatElement {
    fn default() -> Self {
        HeatElement { kind: HeatElemKind::Link, index: -1 }
    }
}

impl HeatElement {
    pub fn node(i: i32) -> HeatElement {
        HeatElement { kind: HeatElemKind::Node, index: i }
    }
    pub fn link(i: i32) -> HeatElement {
        HeatElement { kind: HeatElemKind::Link, index: i }
    }
    pub fn subcatch(i: i32) -> HeatElement {
        HeatElement { kind: HeatElemKind::Subcatch, index: i }
    }
    pub fn cell2d(i: i32) -> HeatElement {
        HeatElement { kind: HeatElemKind::Cell2d, index: i }
    }
}

/// Scope of one override row (D-PE3). Ordered by specificity; 4 and 5 are reserved.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum HeatScope {
    Global = 0,
    Tag = 1,
    Node = 2,
    Link = 3,
    Cell2d = 6,
    MeshTag = 7,
}

/// Which attribute an override row sets (RadiativeConfig, then SedimentConfig fields).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum HeatAttr {
    // RadiativeConfig
    Albedo = 0,
    ShadeFactor = 1,
    SkyView = 2,
    EmissWater = 3,
    EmissLandcover = 4,
    /// NEW (PE2) — NaN means "use air temperature"
    LandcoverTemp = 5,
    // SedimentConfig
    SedThermalDiffusivity = 6,
    SedSoluteDiffusivity = 7,
    SedBedThickness = 8,
    SedGroundDepth = 9,
    SedGroundTemp = 10,
    SedHyporheicVelocity = 11,
    SedDensity = 12,
    SedSpecificHeat = 13,
    SedMeshPorosity = 14,
    SedGroundSource = 15,
    Count = 16,
}

/// True for every ordinal >= SED_THERMAL_DIFFUSIVITY, including COUNT_.
#[inline]
pub fn is_sediment_attr(a: HeatAttr) -> bool {
    (a as i32) >= (HeatAttr::SedThermalDiffusivity as i32)
}

/// One parsed override row, kept verbatim so the serializer can render what the user
/// wrote rather than a resolved expansion of it.
#[derive(Debug, Clone, PartialEq)]
pub struct HeatOverrideRow {
    pub attr: HeatAttr,
    pub scope: HeatScope,
    /// tag or element name; empty at GLOBAL scope
    pub name: Bytes,
    pub value: f64,
}

impl Default for HeatOverrideRow {
    fn default() -> Self {
        HeatOverrideRow {
            attr: HeatAttr::Albedo,
            scope: HeatScope::Global,
            name: Vec::new(),
            value: 0.0,
        }
    }
}
