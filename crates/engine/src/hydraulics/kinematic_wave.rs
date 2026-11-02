// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/KinematicWave.hpp and
// KinematicWave.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The kinematic-wave solver's state and its initialization (KWSolver::init,
// setLinkOrder).
//
// NOT YET TRANSLATED (run-time step functions, with the routing step): KWSolver::execute,
// KWSolver::solveConduit, capConduitLoss, finishRouting, the storage-node update and the
// file-local buildXSP_KW / inflow helpers of KinematicWave.cpp.
//
// Differences in form, never in result:
//   * init's `const XSectGroups&` argument is unused upstream (/*groups*/), so it is
//     dropped.
//   * The members keep upstream's names, trailing underscores included, and are `pub`
//     (upstream's private n_conduits_ too) so the gate's state dumper can read them, as
//     upstream's harness does through `#define private public`.

/// Distance weighting factor
pub const WX: f64 = 0.6;
/// Time weighting factor
pub const WT: f64 = 0.6;
/// Newton convergence tolerance
pub const EPSIL: f64 = 0.001;

/// The kinematic-wave solver over the whole conduit set.
#[derive(Debug, Clone, Default)]
pub struct KWSolver {
    /// Topological link order (upstream to downstream)
    pub sorted_links_: Vec<i32>,
    /// Per-node "already converged this step" flag for storage units
    pub storage_updated_: Vec<i8>,

    /// Previous inlet flow (cfs)
    pub q1_: Vec<f64>,
    /// Previous inlet area (ft2)
    pub a1_: Vec<f64>,
    /// Previous outlet flow (cfs)
    pub q2_: Vec<f64>,
    /// Previous outlet area (ft2)
    pub a2_: Vec<f64>,

    /// Inflow to each conduit (cfs)
    pub q_in_: Vec<f64>,
    /// Inlet area from inflow (ft2)
    pub a_in_: Vec<f64>,
    /// Computed outflow (cfs)
    pub q_out_: Vec<f64>,
    /// Outlet area from Newton solve (ft2)
    pub a_out_: Vec<f64>,
    /// Section factor at inlet
    pub sf_in_: Vec<f64>,
    /// End-of-step upstream flow depth (ft)
    pub y1_: Vec<f64>,
    /// End-of-step downstream flow depth (ft)
    pub y2_: Vec<f64>,

    pub n_conduits_: i32,
}

impl KWSolver {
    pub fn init(&mut self, n_conduits: i32) {
        self.n_conduits_ = n_conduits;
        let un = n_conduits as usize;

        self.q1_.clear();
        self.q1_.resize(un, 0.0);
        self.a1_.clear();
        self.a1_.resize(un, 0.0);
        self.q2_.clear();
        self.q2_.resize(un, 0.0);
        self.a2_.clear();
        self.a2_.resize(un, 0.0);

        self.q_in_.resize(un, 0.0);
        self.a_in_.resize(un, 0.0);
        self.q_out_.resize(un, 0.0);
        self.a_out_.resize(un, 0.0);
        self.sf_in_.resize(un, 0.0);
        self.y1_.clear();
        self.y1_.resize(un, 0.0);
        self.y2_.clear();
        self.y2_.resize(un, 0.0);
    }

    pub fn set_link_order(&mut self, sorted_links: &[i32]) {
        self.sorted_links_ = sorted_links.to_vec();
    }
}
