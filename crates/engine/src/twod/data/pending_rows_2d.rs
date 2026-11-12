// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/data/PendingRows2D.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The authored 2D rows SurfaceRouter2D::initialize drains once the mesh exists.

#![allow(clippy::all)]

use crate::core::options::Bytes;

/// [2D_BOUNDARY_CONDITIONS] row.
#[derive(Debug, Clone, Default)]
pub struct PendingBoundaryRow {
    pub tri: i32,
    pub edge: i32,
    pub bc_type: i32,
    pub param1: f64,
    pub name: Bytes,
    pub group: Bytes,
}

/// [2D_EDGE_CONVEYANCE] row.
#[derive(Debug, Clone)]
pub struct PendingEdgeConveyanceRow {
    pub v_from: i32,
    pub v_to: i32,
    pub conveyance: f64,
}

impl Default for PendingEdgeConveyanceRow {
    fn default() -> Self {
        PendingEdgeConveyanceRow { v_from: 0, v_to: 0, conveyance: 1.0 }
    }
}

/// [2D_INITIAL_QUALITY] row.
#[derive(Debug, Clone)]
pub struct PendingInitialQualityRow {
    pub tri: i32,
    pub tag: Bytes,
    pub all: bool,
    pub species: Bytes,
    pub conc: f64,
}

impl Default for PendingInitialQualityRow {
    fn default() -> Self {
        PendingInitialQualityRow { tri: -1, tag: Bytes::new(), all: false, species: Bytes::new(), conc: 0.0 }
    }
}

/// [2D_INITIAL_BED_QUALITY] row.
#[derive(Debug, Clone)]
pub struct PendingInitialBedQualityRow {
    pub tri: i32,
    pub tag: Bytes,
    pub all: bool,
    pub species: Bytes,
    pub conc: f64,
}

impl Default for PendingInitialBedQualityRow {
    fn default() -> Self {
        PendingInitialBedQualityRow { tri: -1, tag: Bytes::new(), all: false, species: Bytes::new(), conc: 0.0 }
    }
}

/// [2D_BOUNDARY_QUALITY] row.
#[derive(Debug, Clone, Default)]
pub struct PendingBoundaryQualityRow {
    pub tri: i32,
    pub edge: i32,
    pub species: Bytes,
    pub conc: f64,
}
