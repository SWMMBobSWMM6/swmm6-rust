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
//   * GWSolver::execute with its file-static helpers (GWContext, getFluxes, getDxDt) was
//     translated on 5 November 2026 for the runoff step. Custom [GWF] LATERAL / DEEP
//     expressions need mathexpr::evaluate_fast, NOT YET TRANSLATED: execute returns
//     Err("mathexpr::evaluate_fast") where upstream would first evaluate one. The OpenMP
//     split (>= 512 subcatchments) writes disjoint rows; the serial loop gives the same
//     result.

use crate::core::cmath::{cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::datetime;
use crate::core::unit_conversion as ucf;
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

// ============================================================================
// GWSolver::execute (Groundwater.cpp), translated 5 November 2026
// ============================================================================

/// Per-subcatchment context for the ODE callback (legacy gwater.c's static variables).
#[derive(Debug, Clone, Copy, Default)]
struct GWContext {
    porosity: f64,
    field_cap: f64,
    wilt_point: f64,
    k_sat: f64,
    k_slope: f64,
    tension_slope: f64,
    upper_evap_frac: f64,
    lower_evap_depth: f64,
    lower_loss_coeff: f64,
    total_depth: f64,
    a1: f64,
    b1: f64,
    a2: f64,
    b2: f64,
    a3: f64,
    h_star: f64,
    infil: f64,
    avail_evap: f64,
    max_evap: f64,
    sw_head: f64,
    dt: f64,
    area: f64,
    ucf_length: f64,
    ucf_gwflow: f64,
    ucf_rainfall: f64,
    ucf_landarea: f64,
    max_upper_perc: f64,
    max_gw_flow_pos: f64,
    max_gw_flow_neg: f64,
    upper_evap: f64,
    lower_evap: f64,
    upper_perc: f64,
    deep_loss: f64,
    gw_flow: f64,
}

/// getFluxes (legacy getFluxes / getEvapRates / getUpperPerc / getGWFlow). The custom
/// expressions are not evaluated here (execute stops before an aquifer that has one).
fn get_fluxes(c: &mut GWContext, theta: f64, lower_depth: f64) {
    let mut lower_depth = cmax(lower_depth, 0.0);
    lower_depth = cmin(lower_depth, c.total_depth);
    let upper_depth = c.total_depth - lower_depth;

    c.upper_evap = 0.0;
    c.lower_evap = 0.0;
    if c.infil <= 0.0 {
        if theta > c.wilt_point {
            c.upper_evap = c.upper_evap_frac * c.max_evap;
            c.upper_evap = cmin(c.upper_evap, c.avail_evap);
        }
        if c.lower_evap_depth > 0.0 {
            let mut frac = (c.lower_evap_depth - upper_depth) / c.lower_evap_depth;
            frac = cmax(0.0, cmin(frac, 1.0));
            c.lower_evap = frac * (1.0 - c.upper_evap_frac) * c.max_evap;
            c.lower_evap = cmin(c.lower_evap, c.avail_evap - c.upper_evap);
        }
    }

    c.upper_perc = 0.0;
    if upper_depth > 0.0 && theta > c.field_cap {
        let delta = theta - c.porosity;
        let hydcon = c.k_sat * (delta * c.k_slope).exp(); // LIBM
        let delta2 = theta - c.field_cap;
        let dhdz = 1.0 + c.tension_slope * 2.0 * delta2 / upper_depth;
        c.upper_perc = hydcon * dhdz;
    }
    c.upper_perc = cmin(c.upper_perc, c.max_upper_perc);

    c.deep_loss = if c.total_depth > 0.0 { c.lower_loss_coeff * lower_depth / c.total_depth } else { 0.0 };
    c.deep_loss = cmin(c.deep_loss, lower_depth / c.dt);

    c.gw_flow = 0.0;
    if lower_depth > c.h_star {
        let ucf_len = c.ucf_length;
        let hdiff1 = (lower_depth - c.h_star) * ucf_len;
        let t1 = if c.b1 == 0.0 {
            c.a1
        } else if c.b1 == 1.0 {
            c.a1 * hdiff1
        } else {
            c.a1 * hdiff1.powf(c.b1) // LIBM
        };
        let t2 = if c.b2 == 0.0 {
            c.a2
        } else if c.sw_head > c.h_star {
            let hdiff2 = (c.sw_head - c.h_star) * ucf_len;
            if c.b2 == 1.0 {
                c.a2 * hdiff2
            } else {
                c.a2 * hdiff2.powf(c.b2) // LIBM
            }
        } else {
            0.0
        };
        let t3 = c.a3 * lower_depth * c.sw_head * ucf_len * ucf_len;
        c.gw_flow = (t1 - t2 + t3) / c.ucf_gwflow;
        if c.gw_flow < 0.0 && c.a3 != 0.0 {
            c.gw_flow = 0.0;
        }
    }
    if c.gw_flow >= 0.0 {
        c.gw_flow = cmin(c.gw_flow, c.max_gw_flow_pos);
    } else {
        c.gw_flow = cmax(c.gw_flow, c.max_gw_flow_neg);
    }
}

/// getDxDt: the two-variable derivative (moisture fraction, saturated depth).
fn get_dxdt(c: &mut GWContext, x: &[f64; 2], dxdt: &mut [f64; 2]) {
    let theta_val = x[0];
    let lower_d = x[1];
    get_fluxes(c, theta_val, lower_d);
    let upper_d = c.total_depth - lower_d;
    dxdt[0] = if upper_d > 0.0 { (c.infil - c.upper_evap - c.upper_perc) / upper_d } else { 0.0 };
    let denom = c.porosity - theta_val;
    dxdt[1] = if denom > 0.0 { (c.upper_perc - c.deep_loss - c.lower_evap - c.gw_flow) / denom } else { 0.0 };
}

impl GWSolver {
    /// GWSolver::execute: per-subcatchment RKF45 integration of the two-zone aquifer, then
    /// the groundwater mass balance (legacy gwater_getGroundwater / updateMassBal).
    #[allow(clippy::too_many_arguments)]
    pub fn execute(&mut self, ctx: &mut SimulationContext, dt: f64, max_evap: f64, infil_rate: &[f64], sw_head: &[f64], frac_perv: &[f64], perv_evap_rate: &[f64], evap_pattern_month: i32) -> Result<(), &'static str> {
        let n = self.soa_.n_subcatch;
        if n == 0 {
            return Ok(());
        }
        for i in 0..n {
            self.update(ctx, i, dt, max_evap, infil_rate, sw_head, frac_perv, perv_evap_rate, evap_pattern_month)?;
        }

        let mb = &mut ctx.mass_balance;
        let ucf_land = ucf::ucf(ucf::LANDAREA, &ctx.options);
        for i in 0..n as usize {
            if self.soa_.total_depth[i] <= 0.0 {
                continue;
            }
            let area = ctx.subcatches.area[i] / ucf_land;
            let ft2sec = area * dt;
            mb.gw_infil += infil_rate[i] * ft2sec;
            mb.gw_upper_evap += self.soa_.upper_evap[i] * ft2sec;
            mb.gw_lower_evap += self.soa_.lower_evap[i] * ft2sec;
            mb.gw_lower_perc += self.soa_.deep_loss[i] * ft2sec;
            mb.gw_lateral_flow += 0.5 * (self.soa_.old_flow[i] + self.soa_.gw_flow[i]) * ft2sec;
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn update(&mut self, ctx: &SimulationContext, i: i32, dt: f64, max_evap: f64, infil_rate: &[f64], sw_head: &[f64], frac_perv: &[f64], perv_evap_rate: &[f64], evap_pattern_month: i32) -> Result<(), &'static str> {
        let ui = i as usize;
        let soa = &mut self.soa_;
        if soa.total_depth[ui] <= 0.0 {
            soa.gw_flow[ui] = 0.0;
            soa.upper_evap[ui] = 0.0;
            soa.lower_evap[ui] = 0.0;
            soa.deep_loss[ui] = 0.0;
            return Ok(());
        }
        if frac_perv[ui] <= 0.0 {
            return Ok(());
        }

        let mut c = GWContext { porosity: soa.porosity[ui], field_cap: soa.field_cap[ui], wilt_point: soa.wilt_point[ui], k_sat: soa.k_sat[ui], k_slope: soa.k_slope[ui], tension_slope: soa.tension_slope[ui], ..Default::default() };
        let mut evap_frac = soa.upper_evap_frac[ui];
        let pat_idx = soa.upper_evap_pat[ui];
        if pat_idx >= 0 && pat_idx < ctx.patterns.count() {
            let month = if evap_pattern_month > 0 { evap_pattern_month } else { datetime::month_of_year(ctx.current_date) };
            let facs = &ctx.patterns.factors[pat_idx as usize];
            if month >= 1 && month <= facs.len() as i32 {
                evap_frac *= facs[(month - 1) as usize];
            }
        }
        c.upper_evap_frac = evap_frac;
        c.lower_evap_depth = soa.lower_evap_depth[ui];
        c.lower_loss_coeff = soa.lower_loss_coeff[ui];
        c.total_depth = soa.total_depth[ui];
        c.a1 = soa.a1[ui];
        c.b1 = soa.b1[ui];
        c.a2 = soa.a2[ui];
        c.b2 = soa.b2[ui];
        c.a3 = soa.a3[ui];
        c.h_star = soa.h_star[ui];
        c.infil = infil_rate[ui];
        let fp = frac_perv[ui];
        c.max_evap = ctx.forcing.effective_evap_rate(ui, max_evap) * fp;
        c.avail_evap = cmax(c.max_evap - perv_evap_rate[ui], 0.0);
        c.sw_head = sw_head[ui];
        c.dt = dt;
        c.area = ctx.subcatches.area[ui];
        let unit_sys = ucf::get_unit_system(ctx.options.flow_units as i32) as usize;
        c.ucf_length = ucf::UCF[ucf::LENGTH as usize][unit_sys];
        c.ucf_gwflow = ucf::UCF[ucf::GWFLOW as usize][unit_sys];
        c.ucf_rainfall = ucf::UCF[ucf::RAINFALL as usize][unit_sys];
        c.ucf_landarea = ucf::UCF[ucf::LANDAREA as usize][unit_sys];
        if soa.lateral_expr[ui].valid || soa.deep_expr[ui].valid {
            return Err("mathexpr::evaluate_fast ([GWF] expression)");
        }

        let theta_0 = soa.theta[ui];
        let lower_d_0 = soa.lower_depth[ui];
        let upper_d_0 = c.total_depth - lower_d_0;

        let mut v_upper = upper_d_0 * (theta_0 - c.field_cap);
        v_upper = cmax(0.0, v_upper);
        c.max_upper_perc = v_upper / dt;

        c.max_gw_flow_pos = lower_d_0 * c.porosity / dt;

        let gw_capacity_neg = (c.total_depth - lower_d_0) * (c.porosity - theta_0) / dt;
        let mut gw_node = ctx.subcatches.gw_node[ui];
        if gw_node < 0 {
            gw_node = ctx.subcatches.outlet_node[ui];
        }
        let mut node_flow = 0.0;
        if gw_node >= 0 && gw_node < ctx.n_nodes() {
            let un = gw_node as usize;
            let area_ft2 = ctx.subcatches.area[ui] / ucf::ucf(ucf::LANDAREA, &ctx.options);
            if area_ft2 > 0.0 {
                node_flow = (ctx.nodes.inflow[un] + ctx.nodes.volume[un] / dt) / area_ft2;
            }
        }
        c.max_gw_flow_neg = -cmin(gw_capacity_neg, node_flow);

        let mut x = [theta_0, lower_d_0];
        // Legacy ignores the ODE return code.
        let _ = crate::math::ode_solver::integrate(&mut x, 0.0, dt, crate::core::constants::ode::GWTOL, dt, |_t, y: &[f64; 2], dydx: &mut [f64; 2]| get_dxdt(&mut c, y, dydx));

        x[0] = cmax(x[0], c.wilt_point);
        if x[0] >= c.porosity {
            x[0] = c.porosity - crate::core::constants::XTOL;
            x[1] = c.total_depth - crate::core::constants::XTOL;
        }
        x[1] = cmax(x[1], 0.0);
        if x[1] >= c.total_depth {
            x[1] = c.total_depth - crate::core::constants::XTOL;
        }

        soa.theta[ui] = x[0];
        soa.lower_depth[ui] = x[1];

        let fp_gw = frac_perv[ui];
        if fp_gw > 0.0 && c.total_depth > 0.0 {
            soa.max_infil_vol[ui] = (c.total_depth - x[1]) * (c.porosity - x[0]) / fp_gw;
            soa.max_infil_vol[ui] = cmax(soa.max_infil_vol[ui], 0.0);
        } else {
            soa.max_infil_vol[ui] = 0.0;
        }

        get_fluxes(&mut c, x[0], x[1]);
        soa.old_flow[ui] = soa.gw_flow[ui];
        soa.gw_flow[ui] = c.gw_flow;
        soa.upper_evap[ui] = c.upper_evap;
        soa.lower_evap[ui] = c.lower_evap;
        soa.deep_loss[ui] = c.deep_loss;
        Ok(())
    }
}
