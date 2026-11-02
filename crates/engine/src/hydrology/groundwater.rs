// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/Groundwater.hpp and Groundwater.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The groundwater solver's per-subcatchment store (GWSoA, resize) and GWSolver::init,
// which SWMMEngine::initHydrology calls before filling the store from the aquifers.
// GW_VAR_NAMES and gwf_validate were translated earlier in groundwater_parse.rs (the
// [GWF] handler needs them); they are re-exported here.
//
// Differences of form:
//   * The solver member is public so the gate dumper can read it; upstream's is private.
//   * NOT YET TRANSLATED (run-time step functions, phase 1c): GWSolver::execute,
//     batchUpperPerc, batchGWFlow and the file-static ODE helpers (GWContext,
//     getFluxes, getDxDt, ...).

use crate::math::math_expr::Expression;

pub use crate::hydrology::groundwater_parse::{gwf_validate, GW_VAR_NAMES};

/// GW expression variable indices (legacy GWvariables enum).
pub const GWV_HGW: i32 = 0;
pub const GWV_HSW: i32 = 1;
pub const GWV_HCB: i32 = 2;
pub const GWV_HGS: i32 = 3;
pub const GWV_KS: i32 = 4;
pub const GWV_K: i32 = 5;
pub const GWV_THETA: i32 = 6;
pub const GWV_PHI: i32 = 7;
pub const GWV_FI: i32 = 8;
pub const GWV_FU: i32 = 9;
pub const GWV_A: i32 = 10;
pub const GWV_MAX: i32 = 11;

/// Short English gloss per GW_VAR_NAMES entry (editor tooltips / completers).
pub const GW_VAR_DESCRIPTIONS: [&str; 11] = [
    "Water table height above aquifer bottom (ft or m)",
    "Surface-water head at the receiving node above aquifer bottom (ft or m)",
    "Channel bottom height above aquifer bottom (Hstar) (ft or m)",
    "Ground surface height above aquifer bottom (ft or m)",
    "Saturated hydraulic conductivity (in/hr or mm/hr)",
    "Unsaturated hydraulic conductivity (in/hr or mm/hr)",
    "Upper-zone moisture content (fraction)",
    "Soil porosity (fraction)",
    "Surface infiltration rate (in/hr or mm/hr)",
    "Upper-zone percolation rate (in/hr or mm/hr)",
    "Subcatchment area (ac or ha)",
];

/// Per-subcatchment GW state.
#[derive(Debug, Clone, Default)]
pub struct GWSoA {
    pub n_subcatch: i32,

    // Aquifer properties (set at init, constant during sim)
    pub porosity: Vec<f64>,
    pub field_cap: Vec<f64>,
    pub wilt_point: Vec<f64>,
    /// Saturated conductivity (ft/sec)
    pub k_sat: Vec<f64>,
    /// Exponential decay slope
    pub k_slope: Vec<f64>,
    pub tension_slope: Vec<f64>,
    pub upper_evap_frac: Vec<f64>,
    /// Pattern index for monthly evap adjustment (-1 = none)
    pub upper_evap_pat: Vec<i32>,
    pub lower_evap_depth: Vec<f64>,
    /// Deep percolation coeff
    pub lower_loss_coeff: Vec<f64>,
    /// Aquifer thickness (ft)
    pub total_depth: Vec<f64>,
    /// Effective aquifer bottom elevation (ft)
    pub bottom_elev: Vec<f64>,

    // Lateral flow coefficients
    pub a1: Vec<f64>,
    pub b1: Vec<f64>,
    pub a2: Vec<f64>,
    pub b2: Vec<f64>,
    pub a3: Vec<f64>,
    /// Threshold water table height (ft)
    pub h_star: Vec<f64>,

    /// Upper zone moisture content (0-phi)
    pub theta: Vec<f64>,
    /// Lower zone depth (ft)
    pub lower_depth: Vec<f64>,

    /// Previous step GW flow (for trapezoidal avg)
    pub old_flow: Vec<f64>,

    /// Lateral GW flux over full subcatchment area (ft/s)
    pub gw_flow: Vec<f64>,
    pub upper_evap: Vec<f64>,
    pub lower_evap: Vec<f64>,
    pub deep_loss: Vec<f64>,
    /// Gap #40: max infiltration volume upper zone can accept in next step (ft).
    pub max_infil_vol: Vec<f64>,

    /// Per-subcatch compiled lateral flow expression ([GWF] LATERAL).
    pub lateral_expr: Vec<Expression>,
    /// Per-subcatch compiled deep percolation expression ([GWF] DEEP).
    pub deep_expr: Vec<Expression>,
}

fn assign<T: Clone>(v: &mut Vec<T>, n: usize, d: T) {
    v.clear();
    v.resize(n, d);
}

impl GWSoA {
    pub fn resize(&mut self, n: i32) {
        self.n_subcatch = n;
        let un = n as usize;

        assign(&mut self.porosity, un, 0.4);
        assign(&mut self.field_cap, un, 0.2);
        assign(&mut self.wilt_point, un, 0.1);
        assign(&mut self.k_sat, un, 0.0);
        assign(&mut self.k_slope, un, 0.0);
        assign(&mut self.tension_slope, un, 0.0);
        assign(&mut self.upper_evap_frac, un, 0.0);
        assign(&mut self.upper_evap_pat, un, -1);
        assign(&mut self.lower_evap_depth, un, 0.0);
        assign(&mut self.lower_loss_coeff, un, 0.0);
        assign(&mut self.total_depth, un, 0.0);
        assign(&mut self.bottom_elev, un, 0.0);

        assign(&mut self.a1, un, 0.0);
        assign(&mut self.b1, un, 0.0);
        assign(&mut self.a2, un, 0.0);
        assign(&mut self.b2, un, 0.0);
        assign(&mut self.a3, un, 0.0);
        assign(&mut self.h_star, un, 0.0);

        assign(&mut self.theta, un, 0.2);
        assign(&mut self.lower_depth, un, 0.0);
        assign(&mut self.old_flow, un, 0.0);

        assign(&mut self.gw_flow, un, 0.0);
        assign(&mut self.upper_evap, un, 0.0);
        assign(&mut self.lower_evap, un, 0.0);
        assign(&mut self.deep_loss, un, 0.0);
        assign(&mut self.max_infil_vol, un, 0.0);

        // resize, not assign: existing expressions are kept.
        self.lateral_expr.resize(un, Expression::default());
        self.deep_expr.resize(un, Expression::default());
    }
}

/// The groundwater solver.
#[derive(Debug, Clone, Default)]
pub struct GWSolver {
    pub soa_: GWSoA,
}

impl GWSolver {
    pub fn init(&mut self, n_subcatch: i32) {
        self.soa_.resize(n_subcatch);
    }

    pub fn state(&self) -> &GWSoA {
        &self.soa_
    }

    pub fn state_mut(&mut self) -> &mut GWSoA {
        &mut self.soa_
    }
}
