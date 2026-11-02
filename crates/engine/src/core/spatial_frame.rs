// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/SpatialFrame.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The CRS string, map extents, and coordinates of nodes, links, subcatchments and gages,
// each array indexed like its data store. Strings are byte strings.
//
// Upstream's abstract class ICouplingPoint is a trait here. Upstream stores borrowed raw
// pointers (`std::vector<ICouplingPoint*>`, "not owned"); safe Rust has no borrowed
// pointer that outlives its borrow, so each entry is a std::sync::Weak to a Mutex, which
// owns nothing and keeps the frame Send/Sync. A null pointer is `None` and is ignored, as
// upstream ignores NULL. No upstream engine code registers or reads a coupling point
// (SurfaceRouter2D keeps its own), so this choice cannot reach a result.

use crate::core::options::Bytes;
use std::sync::{Mutex, Weak};

/// A coupling point between a 1D object and a 2D mesh cell, implemented by a 2D component.
pub trait ICouplingPoint {
    /// Zero-based index of the 1D object, interpreted according to object_type().
    fn object_index(&self) -> i32;
    /// Category token such as NODE, LINK or SUBCATCH.
    fn object_type(&self) -> &[u8];
    /// Deliver a 1D-side sample to the 2D implementation.
    fn exchange_to_2d(&mut self, value: f64);
    /// Read the sample provided by the 2D implementation.
    fn receive_from_2d(&self) -> f64;
}

/// A non-owning reference to a registered coupling point.
pub type CouplingPointRef = Weak<Mutex<dyn ICouplingPoint + Send>>;

/// Spatial frame: CRS and georeferenced coordinates (model units when `crs` is empty).
#[derive(Debug, Clone, Default)]
pub struct SpatialFrame {
    // ---- coordinate reference system
    /// EPSG code or PROJ string; empty means no CRS (legacy behavior).
    pub crs: Bytes,
    /// True if the CRS represents geographic coordinates (lon/lat).
    pub is_geographic: bool,

    // ---- map extents (from [MAP] DIMENSIONS)
    pub map_x1: f64,
    pub map_y1: f64,
    pub map_x2: f64,
    pub map_y2: f64,
    /// [MAP] UNITS keyword: typically NONE, FEET, METERS or DEGREES.
    pub map_units: Bytes,

    // ---- node coordinates (indexed parallel to NodeData)
    pub node_x: Vec<f64>,
    pub node_y: Vec<f64>,
    /// 1 when the node's coordinates were authored; 0 means the (0, 0) fill. Sized with node_x.
    pub node_has_xy: Vec<u8>,

    // ---- link coordinates (one point per link)
    pub link_x: Vec<f64>,
    pub link_y: Vec<f64>,

    // ---- subcatchment centroid coordinates
    pub subcatch_x: Vec<f64>,
    pub subcatch_y: Vec<f64>,

    // ---- link vertices (polyline points between from-node and to-node)
    pub link_vertices_x: Vec<Vec<f64>>,
    pub link_vertices_y: Vec<Vec<f64>>,

    // ---- subcatchment polygon vertices (rings concatenated with explicit closures; the
    // final ring may be open)
    pub subcatch_polygon_x: Vec<Vec<f64>>,
    pub subcatch_polygon_y: Vec<Vec<f64>>,

    // ---- gage coordinates
    pub gage_x: Vec<f64>,
    pub gage_y: Vec<f64>,

    // ---- 2D coupling interface (R30): registered, not owned
    pub coupling_points: Vec<CouplingPointRef>,
}

// The member initializers (is_geographic = false, map extents = 0.0) are Rust's defaults,
// so Default is derived.

impl SpatialFrame {
    /// Append a coupling point without deduplication; `None` (NULL) is ignored.
    pub fn register_coupling_point(&mut self, point: Option<CouplingPointRef>) {
        if let Some(p) = point {
            self.coupling_points.push(p);
        }
    }

    pub fn has_2d_coupling(&self) -> bool {
        !self.coupling_points.is_empty()
    }

    /// Visits outer and inner geometry vectors; node_has_xy and coupling_points are not
    /// compacted (as upstream).
    pub fn shrink_to_fit(&mut self) {
        self.node_x.shrink_to_fit();
        self.node_y.shrink_to_fit();
        self.link_x.shrink_to_fit();
        self.link_y.shrink_to_fit();
        self.subcatch_x.shrink_to_fit();
        self.subcatch_y.shrink_to_fit();
        self.gage_x.shrink_to_fit();
        self.gage_y.shrink_to_fit();

        self.link_vertices_x.shrink_to_fit();
        for v in &mut self.link_vertices_x {
            v.shrink_to_fit();
        }
        self.link_vertices_y.shrink_to_fit();
        for v in &mut self.link_vertices_y {
            v.shrink_to_fit();
        }

        self.subcatch_polygon_x.shrink_to_fit();
        for v in &mut self.subcatch_polygon_x {
            v.shrink_to_fit();
        }
        self.subcatch_polygon_y.shrink_to_fit();
        for v in &mut self.subcatch_polygon_y {
            v.shrink_to_fit();
        }
    }
}
