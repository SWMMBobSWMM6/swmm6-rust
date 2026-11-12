// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/data/BoundaryData.hpp and the RuntimeReceipt
// struct of src/engine/2d/coupling/RuntimeCoupling.hpp (HydroCouple/openswmm.engine, branch
// swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache License,
// Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Per-edge-slot boundary conditions, filled by SurfaceRouter2D::initialize from the pending
// [2D_BOUNDARY_CONDITIONS] rows (open() leaves the store empty). The std::map members are
// BTreeMap (the same key order). RuntimeReceipt lives here until the wave that translates
// RuntimeCoupling.hpp moves it next to RuntimeSource.

#![allow(clippy::all)]

use std::cell::RefCell;
use std::collections::BTreeMap;

use crate::core::options::Bytes;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i8)]
pub enum BoundaryType {
    #[default]
    Wall = 0,
    NormalFlow = 1,
    SpecifiedStage = 2,
    SpecifiedFlow = 3,
    RatingCurve = 4,
}

/// coupling/RuntimeCoupling.hpp RuntimeReceipt.
#[derive(Debug, Clone, Default)]
pub struct RuntimeReceipt {
    pub requested: f64,
    pub applied: f64,
    pub water_in: f64,
    pub water_out: f64,
    pub heat_requested: f64,
    pub heat_applied: f64,
    pub last_flow: f64,
    pub last_dt: f64,
    pub requested_mass: Vec<f64>,
    pub applied_mass: Vec<f64>,
}

/// BoundaryData::RuntimeOriginal (no member initializers upstream; always built whole).
#[derive(Debug, Clone, Default)]
pub struct RuntimeOriginal {
    pub type_: i8,
    pub head: f64,
    pub flow: f64,
    pub slope: f64,
    pub stage_ts: i32,
    pub flow_ts: i32,
    pub rating: i32,
    pub stage_name: Bytes,
    pub flow_name: Bytes,
    pub rating_name: Bytes,
}

#[derive(Debug, Clone, Default)]
pub struct BoundaryData {
    pub edge_bc_type: Vec<i8>,
    pub edge_wall_drag: Vec<f64>,
    pub edge_is_exterior: Vec<u8>,
    pub edge_bed_slope: Vec<f64>,
    pub edge_bc_head: Vec<f64>,
    pub edge_bc_tseries: Vec<i32>,
    pub edge_bc_tseries_name: Vec<Bytes>,
    pub edge_bc_cum_flux: Vec<f64>,
    pub edge_bc_flow: Vec<f64>,
    pub edge_bc_flow_tseries: Vec<i32>,
    pub edge_bc_flow_tseries_name: Vec<Bytes>,
    pub edge_bc_rating_curve: Vec<i32>,
    pub edge_bc_rating_curve_name: Vec<Bytes>,
    pub runtime_original: BTreeMap<i32, RuntimeOriginal>,
    pub runtime_owner: BTreeMap<i32, Bytes>,
    /// `mutable` upstream: ensureReceipt fills it through a const method
    pub receipts: RefCell<BTreeMap<i32, RuntimeReceipt>>,
}

impl BoundaryData {
    pub fn ensure_receipt(&self, slot: i32, species: i32) {
        let mut r = self.receipts.borrow_mut();
        let receipt = r.entry(slot).or_default();
        if receipt.applied_mass.len() != species as usize {
            receipt.requested_mass.clear();
            receipt.requested_mass.resize(species as usize, 0.0);
            receipt.applied_mass.clear();
            receipt.applied_mass.resize(species as usize, 0.0);
        }
    }

    pub fn save_runtime(&mut self, slot: i32) {
        if self.runtime_original.contains_key(&slot) {
            return;
        }
        let s = slot as usize;
        let b = RuntimeOriginal {
            type_: self.edge_bc_type[s],
            head: self.edge_bc_head[s],
            flow: self.edge_bc_flow[s],
            slope: self.edge_bed_slope[s],
            stage_ts: self.edge_bc_tseries[s],
            flow_ts: self.edge_bc_flow_tseries[s],
            rating: self.edge_bc_rating_curve[s],
            stage_name: self.edge_bc_tseries_name[s].clone(),
            flow_name: self.edge_bc_flow_tseries_name[s].clone(),
            rating_name: self.edge_bc_rating_curve_name[s].clone(),
        };
        self.runtime_original.insert(slot, b);
    }

    pub fn clear_runtime(&mut self, slot: i32) {
        let b = match self.runtime_original.remove(&slot) {
            Some(b) => b,
            None => return,
        };
        let s = slot as usize;
        self.edge_bc_type[s] = b.type_;
        self.edge_bc_head[s] = b.head;
        self.edge_bc_flow[s] = b.flow;
        self.edge_bed_slope[s] = b.slope;
        self.edge_bc_tseries[s] = b.stage_ts;
        self.edge_bc_flow_tseries[s] = b.flow_ts;
        self.edge_bc_rating_curve[s] = b.rating;
        self.edge_bc_tseries_name[s] = b.stage_name;
        self.edge_bc_flow_tseries_name[s] = b.flow_name;
        self.edge_bc_rating_curve_name[s] = b.rating_name;
        self.runtime_owner.remove(&slot);
    }

    /// resize(): every per-slot array assigned to `n_edges` (the name arrays resized, old
    /// names kept), the runtime maps cleared.
    pub fn resize(&mut self, n_edges: i32) {
        self.runtime_original.clear();
        self.runtime_owner.clear();
        self.receipts.borrow_mut().clear();
        let n = n_edges as usize;
        self.edge_bc_type.clear();
        self.edge_bc_type.resize(n, BoundaryType::Wall as i8);
        self.edge_bed_slope.clear();
        self.edge_bed_slope.resize(n, 0.0);
        self.edge_wall_drag.clear();
        self.edge_wall_drag.resize(n, 0.0);
        self.edge_is_exterior.clear();
        self.edge_is_exterior.resize(n, 0);
        self.edge_bc_head.clear();
        self.edge_bc_head.resize(n, 0.0);
        self.edge_bc_tseries.clear();
        self.edge_bc_tseries.resize(n, -1);
        self.edge_bc_tseries_name.resize(n, Bytes::new());
        self.edge_bc_cum_flux.clear();
        self.edge_bc_cum_flux.resize(n, 0.0);
        self.edge_bc_flow.clear();
        self.edge_bc_flow.resize(n, 0.0);
        self.edge_bc_flow_tseries.clear();
        self.edge_bc_flow_tseries.resize(n, -1);
        self.edge_bc_flow_tseries_name.resize(n, Bytes::new());
        self.edge_bc_rating_curve.clear();
        self.edge_bc_rating_curve.resize(n, -1);
        self.edge_bc_rating_curve_name.resize(n, Bytes::new());
    }

    pub fn size(&self) -> i32 {
        self.edge_bc_type.len() as i32
    }
}
