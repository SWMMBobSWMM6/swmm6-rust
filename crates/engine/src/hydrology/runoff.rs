// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/Runoff.hpp and Runoff.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Subcatchment runoff (three-subarea nonlinear reservoir): the RunoffSoA store with
// resize and computeAlpha, the RunoffSolver members, and RunoffSolver::init, which
// SWMMEngine::initHydrology calls.
//
// Differences of form:
//   * init's `throw std::invalid_argument` on a bad reviewed-area override becomes
//     Err(message), returned before the solver is touched, as upstream throws before
//     clearing it. initHydrology passes no overrides, so it never fails there.
//   * The solver's members are public so the gate dumper (an example crate) can read
//     them; upstream's are private.
//   * Run-time functions translated 5 November 2026 for the runoff step:
//     updatePondedDepth, infilGetInfil, nativeInfilFullLid and execute. execute is
//     translated for the call SWMMEngine::stepRunoff makes (no infiltration boundary, no
//     completed-source forcing, no spatial fractions: the 2D coupling's arguments, all
//     null in a 1D build); the source-owner map is taken as stepRunoff passes it
//     (ordinarySourceOwners(), None in a 1D build). The per-subarea lambda processSubarea
//     is the method process_subarea, with the subarea's ponded depth copied in and out
//     (nothing else reads it during the call) and the lambda's captured accumulators in
//     SubareaAcc. The OpenMP split (>= 512 subcatchments) touches disjoint rows and gives
//     the same result; the translation runs the serial loop.
//   * NOT YET TRANSLATED: pendingRoutingVolume, reviewedAreas (2D), checkpoint,
//     validateCheckpoint.

use crate::core::cmath::{cclamp, cmax, cmin};
use crate::data::forcing::ForcingMode;
use crate::hydrology::gage;
use crate::hydrology::infil_bank::Factors;
use crate::hydrology::ueb_types::SnowKind;
use crate::core::constants;
use crate::core::context::SimulationContext;
use crate::core::unit_conversion as ucf;
use crate::hydrology::infil_bank::InfilBank;

// ============================================================================
// Completed-source types (declared in Runoff.hpp; used by execute)
// ============================================================================

/// Completed-interval rates, already resolved to internal ft/s.
#[derive(Debug, Clone, Copy)]
pub struct RunoffSourceForcing {
    pub subcatch: i32,
    pub rain: f64,
    pub pet: f64,
    /// ft/s by IMPERV0/IMPERV1/PERV; never rain
    pub surface_return: [f64; 3],
    /// before subcatchment scales and phase partition
    pub raw_precipitation: f64,
    /// optional exposed atmospheric footprints
    pub pet_imperv: f64,
    pub pet_perv: f64,
}

impl Default for RunoffSourceForcing {
    fn default() -> Self {
        RunoffSourceForcing { subcatch: -1, rain: 0.0, pet: 0.0, surface_return: [0.0; 3], raw_precipitation: -1.0, pet_imperv: -1.0, pet_perv: -1.0 }
    }
}

/// Actual completed subarea water budget (ft3; area in ft2). Every initializer is zero.
#[derive(Debug, Clone, Copy, Default)]
pub struct CompletedSubareaWater {
    pub area: f64,
    pub opening: f64,
    pub rain: f64,
    pub runon: f64,
    pub evaporation: f64,
    pub infiltration: f64,
    pub spatial_infiltration: f64,
    pub outflow: f64,
    pub outlet: f64,
    pub closing: f64,
}

// ============================================================================
// Constants
// ============================================================================

/// Manning's exponent (legacy subcatch.c:55 literal, not 5.0/3.0)
pub const MEXP: f64 = 1.6666667;
/// ODE solver tolerance (matching legacy)
pub const ODETOL: f64 = 0.0001;
/// Manning's US customary constant
pub const PHI: f64 = 1.486;

// ============================================================================
// Per-subcatchment subarea state (SoA)
// ============================================================================

#[derive(Debug, Clone, Default)]
pub struct RunoffSoA {
    pub n_subcatch: i32,

    /// Subcatchment area (ft2)
    pub area: Vec<f64>,
    /// Subcatchment width (ft)
    pub width: Vec<f64>,
    /// Average slope (ft/ft)
    pub slope: Vec<f64>,
    /// Impervious fraction (0-1)
    pub imperv_pct: Vec<f64>,
    /// IMPERV0 (no depression storage)
    pub frac_imperv0: Vec<f64>,
    /// IMPERV1 (with depression storage)
    pub frac_imperv1: Vec<f64>,

    pub alpha_imperv: Vec<f64>,
    pub alpha_perv: Vec<f64>,

    pub ds_imperv: Vec<f64>,
    pub ds_perv: Vec<f64>,

    pub n_imperv: Vec<f64>,
    pub n_perv: Vec<f64>,

    pub depth_imperv0: Vec<f64>,
    pub depth_imperv1: Vec<f64>,
    pub depth_perv: Vec<f64>,

    pub old_runoff_imperv0: Vec<f64>,
    pub old_runoff_imperv1: Vec<f64>,
    pub old_runoff_perv: Vec<f64>,

    /// Total runoff rate (cfs)
    pub runoff: Vec<f64>,
    pub evap_loss: Vec<f64>,
    pub infil_loss: Vec<f64>,
    /// Pervious-subarea evaporation this step (ft3) — legacy Vpevap
    pub perv_evap_vol: Vec<f64>,
    pub actual_perv_evap_vol: Vec<f64>,
    /// Completed partition diagnostics: ft3 spatial intake and ft/s outside native capacity.
    pub spatial_infil_vol: Vec<f64>,
    pub native_infil_rate: Vec<f64>,
    /// Non-LID infiltration this step (ft3) — legacy Vinfil
    pub infil_vol: Vec<f64>,
    pub subarea_runoff_rate: Vec<f64>,

    pub imperv_runoff_cfs: Vec<f64>,
    pub perv_runoff_cfs: Vec<f64>,
    pub outflow_vol: Vec<f64>,
    pub completed: Vec<[CompletedSubareaWater; 3]>,
}

impl RunoffSoA {
    /// Install n rows: geometry/depth/rates/diagnostics zero, Manning n 0.01 impervious
    /// and 0.1 pervious.
    pub fn resize(&mut self, n: i32) {
        self.n_subcatch = n;
        let un = n as usize;
        fn assign<T: Clone>(v: &mut Vec<T>, n: usize, d: T) {
            v.clear();
            v.resize(n, d);
        }
        assign(&mut self.area, un, 0.0);
        assign(&mut self.width, un, 0.0);
        assign(&mut self.slope, un, 0.0);
        assign(&mut self.imperv_pct, un, 0.0);
        assign(&mut self.frac_imperv0, un, 0.0);
        assign(&mut self.frac_imperv1, un, 0.0);

        assign(&mut self.alpha_imperv, un, 0.0);
        assign(&mut self.alpha_perv, un, 0.0);
        assign(&mut self.ds_imperv, un, 0.0);
        assign(&mut self.ds_perv, un, 0.0);
        assign(&mut self.n_imperv, un, 0.01);
        assign(&mut self.n_perv, un, 0.1);

        assign(&mut self.depth_imperv0, un, 0.0);
        assign(&mut self.depth_imperv1, un, 0.0);
        assign(&mut self.depth_perv, un, 0.0);

        assign(&mut self.old_runoff_imperv0, un, 0.0);
        assign(&mut self.old_runoff_imperv1, un, 0.0);
        assign(&mut self.old_runoff_perv, un, 0.0);

        assign(&mut self.runoff, un, 0.0);
        assign(&mut self.evap_loss, un, 0.0);
        assign(&mut self.infil_loss, un, 0.0);
        assign(&mut self.perv_evap_vol, un, 0.0);
        assign(&mut self.actual_perv_evap_vol, un, 0.0);
        assign(&mut self.infil_vol, un, 0.0);
        assign(&mut self.spatial_infil_vol, un, 0.0);
        assign(&mut self.native_infil_rate, un, 0.0);
        assign(&mut self.subarea_runoff_rate, un, 0.0);
        assign(&mut self.imperv_runoff_cfs, un, 0.0);
        assign(&mut self.perv_runoff_cfs, un, 0.0);
        assign(&mut self.outflow_vol, un, 0.0);
        assign(&mut self.completed, un, [CompletedSubareaWater::default(); 3]);
    }

    /// Alpha = PHI * width / area * sqrt(slope) / N (legacy subcatch.c:398-399 operand
    /// order). Rows with total area <= 0 keep their previous coefficients.
    pub fn compute_alpha(&mut self) {
        for i in 0..self.n_subcatch {
            let ui = i as usize;
            let sq_slope = self.slope[ui].sqrt();
            let fi = self.imperv_pct[ui];
            let fp = 1.0 - fi;
            if self.area[ui] > 0.0 {
                let area_imperv = self.area[ui] * fi;
                let area_perv = self.area[ui] * fp;
                self.alpha_imperv[ui] = if self.n_imperv[ui] > 0.0 && area_imperv > 0.0 { PHI * self.width[ui] / area_imperv * sq_slope / self.n_imperv[ui] } else { 0.0 };
                self.alpha_perv[ui] = if self.n_perv[ui] > 0.0 && area_perv > 0.0 { PHI * self.width[ui] / area_perv * sq_slope / self.n_perv[ui] } else { 0.0 };
            }
        }
    }
}

// ============================================================================
// Runoff solver
// ============================================================================

#[derive(Debug, Clone, Default)]
pub struct RunoffSolver {
    pub soa_: RunoffSoA,
    /// Reviewed denominator, absent on the legacy path.
    pub spatial_full_area_ft2_: Vec<f64>,

    /// Infiltration state (one per subcatchment)
    pub infil_bank_: InfilBank,
    /// The InfilFactor applied to each subcatchment this step (pattern or global)
    pub infil_factor_used_: Vec<f64>,
    /// Volume history for selected interval routing.
    pub completed_step_seconds_: Vec<f64>,

    // Working buffers (reused each step, sized to n_subcatch)
    pub precip_: Vec<f64>,
    pub evap_rate_: Vec<f64>,
    pub infil_rate_: Vec<f64>,
}

impl RunoffSolver {
    /// Rebuild owned runoff/native infiltration state for one model. `spatial_areas`
    /// are optional (zero-based source, m2 non-LID area) overrides; initHydrology passes
    /// none.
    pub fn init(&mut self, ctx: &mut SimulationContext, spatial_areas: &[(i32, f64)]) -> Result<(), &'static str> {
        let n = ctx.n_subcatches();
        let mut external_area: Vec<f64> = vec![-1.0; n as usize];
        for &(i, area) in spatial_areas {
            if i < 0 || i >= n || !area.is_finite() || area < 0.0 || external_area[i as usize] >= 0.0 {
                return Err("Invalid or duplicate reviewed non-LID source area.");
            }
            external_area[i as usize] = area;
        }
        self.soa_.resize(n);
        self.completed_step_seconds_.clear();
        self.completed_step_seconds_.resize(n as usize, 0.0);
        self.spatial_full_area_ft2_.clear();
        if !spatial_areas.is_empty() {
            self.spatial_full_area_ft2_.clear();
            self.spatial_full_area_ft2_.resize(n as usize, -1.0);
        }

        let ucf_area = ucf::ucf(ucf::LANDAREA, &ctx.options);
        let ucf_depth = ucf::ucf(ucf::RAINDEPTH, &ctx.options);

        let sc = &ctx.subcatches;
        for i in 0..n {
            let ui = i as usize;
            // Gap #23: subareas use the non-LID area (legacy subcatch_validate).
            let full_area_ft2 = sc.area[ui] / ucf_area;
            let lid_area_ft2 = sc.total_lid_area_ft2[ui]; // already in ft2
            self.soa_.area[ui] = cmax(0.0, full_area_ft2 - lid_area_ft2);
            if external_area[ui] >= 0.0 {
                self.soa_.area[ui] = external_area[ui] / (0.3048 * 0.3048);
                self.spatial_full_area_ft2_[ui] = self.soa_.area[ui] + lid_area_ft2;
            }
            self.soa_.width[ui] = sc.width[ui];
            self.soa_.slope[ui] = sc.slope[ui];
            self.soa_.imperv_pct[ui] = sc.frac_imperv[ui];
            // PARITY subcatch.c:268-270 — multiply by the AUTHORED percent, divide by
            // 100 afterwards.
            let fi_ = sc.frac_imperv[ui];
            let pz_ = sc.pct_zero[ui];
            self.soa_.frac_imperv0[ui] = fi_ * pz_ / 100.0;
            self.soa_.frac_imperv1[ui] = fi_ * (1.0 - pz_ / 100.0);
            self.soa_.n_imperv[ui] = sc.n_imperv[ui];
            self.soa_.n_perv[ui] = sc.n_perv[ui];
            self.soa_.ds_imperv[ui] = sc.ds_imperv[ui] / ucf_depth;
            self.soa_.ds_perv[ui] = sc.ds_perv[ui] / ucf_depth;
        }
        self.soa_.compute_alpha();

        self.infil_bank_.init(n);
        self.infil_factor_used_.clear();
        self.infil_factor_used_.resize(n as usize, 1.0);
        for i in 0..n {
            let ui = i as usize;
            // Preserve legacy subcatchment zero-drying-time normalization.
            let dry_time = if sc.infil_p4[ui] == 0.0 { constants::TINY } else { sc.infil_p4[ui] };
            let p: [f64; 5] = [sc.infil_p1[ui], sc.infil_p2[ui], sc.infil_p3[ui], dry_time, sc.infil_p5[ui]];
            self.infil_bank_.set_method(i, sc.infil_model[ui], &p, &ctx.options);
        }
        Ok(())
    }

    /// The last recorded native infiltration multiplier for one source.
    pub fn infil_factor_used(&self, i: i32) -> f64 {
        self.infil_factor_used_[i as usize]
    }

    pub fn soa(&self) -> &RunoffSoA {
        &self.soa_
    }

    pub fn soa_mut(&mut self) -> &mut RunoffSoA {
        &mut self.soa_
    }
}

// RunoffSolver::infil_get_state / infil_set_state (translated 5 November 2026 for the legacy
// hot-start runoff section).
impl RunoffSolver {
    pub fn infil_get_state(&self, i: i32, model: &mut i32, state: &mut [f64; 6]) {
        self.infil_bank_.pack(i, model, state);
    }
    pub fn infil_set_state(&mut self, i: i32, model: i32, state: &[f64; 6]) {
        self.infil_bank_.unpack(i, model, state);
    }
}

/// processSubarea's captured accumulators (legacy Vevap, Vpevap, Vinfil, Voutflow, ...).
#[derive(Debug, Clone, Copy, Default)]
struct SubareaAcc {
    vevap: f64,
    vpevap: f64,
    actual_perv_evap: f64,
    vinfil: f64,
    voutflow: f64,
}

// RunoffSolver run-time functions (translated 5 November 2026 for the runoff step).
impl RunoffSolver {
    /// updatePondedDepth: legacy subcatch.c updatePondedDepth + odesolve_integrate.
    /// `t_runoff` receives legacy's `*dt = tx`.
    pub fn update_ponded_depth(depth: &mut f64, inflow: f64, alpha: f64, d_store: f64, dt: f64, t_runoff: &mut f64) {
        let mut tx = dt;

        if *depth + inflow * tx <= d_store {
            *depth += inflow * tx;
        } else {
            let dx = d_store - *depth;
            if dx > 0.0 && inflow > 0.0 {
                tx -= dx / inflow;
                *depth = d_store;
            }

            if alpha > 0.0 && tx > 0.0 {
                let captured_inflow = inflow;
                let captured_alpha = alpha;
                let captured_dstore = d_store;
                let mut y = [*depth];
                crate::math::ode_solver::integrate(&mut y, 0.0, tx, ODETOL, tx, |_t, d: &[f64; 1], dddt: &mut [f64; 1]| {
                    let rx = d[0] - captured_dstore;
                    let outflow = if rx > 0.0 { captured_alpha * rx.powf(MEXP) } else { 0.0 }; // LIBM
                    dddt[0] = captured_inflow - outflow;
                });
                *depth = y[0];
            } else {
                tx = cmax(tx, 0.0);
                *depth += inflow * tx;
            }
        }

        *depth = cmax(*depth, 0.0);
        *t_runoff = tx;
    }

    /// infilGetInfil (legacy infil_getInfil).
    #[allow(clippy::too_many_arguments)]
    pub fn infil_get_infil(&mut self, _ctx: &SimulationContext, i: i32, precip: f64, runon: f64, depth: f64, dt: f64, local_infil: f64, recovery_factor: f64) -> f64 {
        self.infil_bank_.rate_feet(i, precip, runon, depth, dt, Factors { infiltration: local_infil, recovery: recovery_factor })
    }

    /// nativeInfilFullLid (legacy findNativeInfil for a subcatchment without pervious
    /// non-LID area): advances the subcatchment's infiltration state.
    pub fn native_infil_full_lid(&mut self, ctx: &SimulationContext, i: i32, dt: f64, recovery_factor: f64) -> f64 {
        let ui = i as usize;
        let rain = ctx.subcatches.rainfall[ui];
        let runon = ctx.subcatches.runon_rate[ui];
        let depth = self.soa_.depth_perv[ui];
        let f = self.infil_factor_used_[ui];
        self.infil_get_infil(ctx, i, rain, runon, depth, dt, f, recovery_factor)
    }

    /// execute: one runoff step for every subcatchment, as SWMMEngine::stepRunoff calls it
    /// (no boundary, completed forcing or spatial fractions; `source_owners` as given).
    #[allow(clippy::too_many_arguments)]
    pub fn execute(&mut self, ctx: &mut SimulationContext, dt: f64, evap_rate_in: f64, infil_factor: f64, recovery_factor: f64, month: i32, source_owners: Option<&Vec<i32>>) -> Result<(), &'static str> {
        let n = self.soa_.n_subcatch;
        if let Some(o) = source_owners {
            if o.len() != n as usize {
                return Err("Ordinary runoff ownership requires a matching exclusive source map.");
            }
        }
        if n == 0 {
            return Ok(());
        }
        let count = n;
        let owned = |i: i32| -> bool { matches!(source_owners, Some(o) if o[i as usize] >= 0) };

        let un = n as usize;
        self.precip_.resize(un, 0.0);
        self.evap_rate_.resize(un, 0.0);
        self.infil_rate_.resize(un, 0.0);
        let ucf_landarea = ucf::ucf(ucf::LANDAREA, &ctx.options);

        // ----- Step 1: Rainfall -> net precip (ft/sec)
        for i in 0..count {
            if owned(i) {
                continue;
            }
            let ui = i as usize;
            let p = gage::split_precip(ctx, ui);
            let mut rain = p.rainfall;
            if ui < ctx.forcing.subcatch_rainfall_mode.len() && ctx.forcing.subcatch_rainfall_mode[ui] != ForcingMode::None {
                let mut rain_inhr = p.rainfall * ucf::ucf(ucf::RAINFALL, &ctx.options);
                rain_inhr = ctx.forcing.effective_rainfall(ui, rain_inhr);
                rain = rain_inhr / ucf::ucf(ucf::RAINFALL, &ctx.options);
            }
            let snow = ctx.forcing.effective_snowfall(ui, p.snowfall);
            self.precip_[ui] = rain + snow;
            ctx.subcatches.rainfall[ui] = self.precip_[ui];
        }

        // ----- Step 2: Evaporation rate
        for i in 0..count {
            if owned(i) {
                continue;
            }
            let ui = i as usize;
            let is_dry_only = ctx.options.evap_dry_only;
            let rain = self.precip_[ui] * ucf::ucf(ucf::RAINFALL, &ctx.options);
            let broadcast = if is_dry_only && rain > 0.0 { 0.0 } else { evap_rate_in };
            self.evap_rate_[ui] = ctx.forcing.effective_evap_rate(ui, broadcast);
        }

        // ----- Step 3: Per-subcatchment subarea processing
        for i in 0..count {
            if owned(i) {
                continue;
            }
            self.process_source(ctx, i, dt, infil_factor, recovery_factor, month, ucf_landarea);
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn process_source(&mut self, ctx: &mut SimulationContext, i: i32, dt: f64, infil_factor: f64, recovery_factor: f64, month: i32, ucf_landarea: f64) {
        let ui = i as usize;
        let fi = self.soa_.imperv_pct[ui];
        let fp = 1.0 - fi;
        let f0 = self.soa_.frac_imperv0[ui];
        let f1 = self.soa_.frac_imperv1[ui];
        let total_area = self.soa_.area[ui];

        let precip = self.precip_[ui];
        let mut snow_open_imperv = 1.0;
        let mut snow_open_perv = 1.0;
        let snowpack = ctx.subcatches.snowpack[ui];
        if !ctx.options.ignore_snow_melt && snowpack >= 0 && ctx.snowpacks.config.model_for(snowpack) == SnowKind::Ueb {
            let flux = &ctx.snowpacks.runtime.last;
            if flux.len() >= 3 * (ui + 1) {
                let plow = if snowpack < ctx.snowpacks.plowable.len() as i32 { ctx.snowpacks.plowable[snowpack as usize][6] } else { 0.0 };
                snow_open_imperv = cclamp(1.0 - plow * flux[ui * 3].active_fraction - (1.0 - plow) * flux[ui * 3 + 1].active_fraction, 0., 1.);
                snow_open_perv = cclamp(1.0 - flux[ui * 3 + 2].active_fraction, 0., 1.);
            }
        }

        let alpha_i = self.soa_.alpha_imperv[ui];
        let alpha_p = self.soa_.alpha_perv[ui];

        let mut acc = SubareaAcc::default();
        self.soa_.native_infil_rate[ui] = 0.0;

        // --- Subarea inflows (legacy subcatch_getRunon), in its order
        let route_mode = ctx.subcatches.subarea_routing[ui];
        let pct = ctx.subcatches.pct_routed[ui];
        let f_outlet = 1.0 - pct;
        let history_scale = 1.0;
        let runon_sub = if total_area > 0.0 { ctx.subcatches.runon_rate[ui] } else { 0.0 };
        let runon_imperv0 = runon_sub;
        let mut runon_imperv1 = runon_sub;
        let mut runon_perv = runon_sub;

        if route_mode == 2 && fp > 0.0 {
            let q1 = self.soa_.old_runoff_imperv0[ui] * f0;
            let q2 = self.soa_.old_runoff_imperv1[ui] * f1;
            let q = q1 + q2;
            runon_perv += q * (1.0 - f_outlet) / fp * history_scale;
        } else if route_mode == 1 && f1 > 0.0 {
            let q = self.soa_.old_runoff_perv[ui];
            runon_imperv1 += q * (1.0 - f_outlet) * fp / f1 * history_scale;
        }

        {
            let q_ret = ctx.subcatches.lid_return_to_perv_cfs[ui];
            if q_ret > 0.0 && total_area > 0.0 && fp > 0.0 {
                let perv_area = total_area * fp;
                runon_perv += q_ret / perv_area;
            }
            ctx.subcatches.lid_return_to_perv_cfs[ui] = 0.0;
        }

        let mut precip_imperv = precip;
        let mut precip_perv = precip;
        if ctx.subcatches.snowpack[ui] >= 0 && !ctx.options.ignore_snow_melt {
            let sni = ctx.subcatches.snow_net_imperv[ui];
            let snp = ctx.subcatches.snow_net_perv[ui];
            if sni >= 0.0 {
                precip_imperv = sni;
            }
            if snp >= 0.0 {
                precip_perv = snp;
            }
        }

        let mut evap_rate = self.evap_rate_[ui] * snow_open_imperv;
        let n_imperv = self.soa_.n_imperv[ui];
        let mut d = self.soa_.depth_imperv0[ui];
        let runoff0 = self.process_subarea(ctx, i, &mut d, alpha_i, 0.0, f0, false, runon_imperv0, n_imperv, precip_imperv, evap_rate, dt, infil_factor, recovery_factor, month, total_area, &mut acc);
        self.soa_.depth_imperv0[ui] = d;
        let ds_imperv = self.soa_.ds_imperv[ui];
        let mut d = self.soa_.depth_imperv1[ui];
        let runoff1 = self.process_subarea(ctx, i, &mut d, alpha_i, ds_imperv, f1, false, runon_imperv1, n_imperv, precip_imperv, evap_rate, dt, infil_factor, recovery_factor, month, total_area, &mut acc);
        self.soa_.depth_imperv1[ui] = d;

        let mut alpha_p_adj = alpha_p;
        let mut ds_perv_adj = self.soa_.ds_perv[ui];
        if month >= 0 {
            let umon = month as usize;
            let monthly = |pat: &Vec<i32>, out: &mut f64| -> bool {
                if ui >= pat.len() {
                    return false;
                }
                let pi = pat[ui];
                if pi < 0 {
                    return false;
                }
                let upi = pi as usize;
                if upi >= ctx.patterns.types.len() || ctx.patterns.types[upi] != 0 {
                    return false;
                }
                let facs = &ctx.patterns.factors[upi];
                if umon >= facs.len() {
                    return false;
                }
                *out = facs[umon];
                true
            };
            let mut f = 0.0;
            if monthly(&ctx.subcatch_d_store_pattern, &mut f) && f >= 0.0 {
                ds_perv_adj *= f;
            }
            if monthly(&ctx.subcatch_n_perv_pattern, &mut f) {
                alpha_p_adj = if f <= 0.0 { 0.0 } else { alpha_p_adj / f };
            }
        }

        evap_rate = self.evap_rate_[ui] * snow_open_perv;
        let n_perv = self.soa_.n_perv[ui];
        let mut d = self.soa_.depth_perv[ui];
        let runoff_p = self.process_subarea(ctx, i, &mut d, alpha_p_adj, ds_perv_adj, fp, true, runon_perv, n_perv, precip_perv, evap_rate, dt, infil_factor, recovery_factor, month, total_area, &mut acc);
        self.soa_.depth_perv[ui] = d;

        self.soa_.old_runoff_imperv0[ui] = runoff0;
        self.soa_.old_runoff_imperv1[ui] = runoff1;
        self.soa_.old_runoff_perv[ui] = runoff_p;
        self.completed_step_seconds_[ui] = dt;

        {
            let mut sum = 0.0;
            if total_area > 0.0 {
                sum += runoff0 * (total_area * f0);
                sum += runoff1 * (total_area * f1);
                sum += runoff_p * (total_area * fp);
            }
            let full_area = if !self.spatial_full_area_ft2_.is_empty() && self.spatial_full_area_ft2_[ui] >= 0.0 { self.spatial_full_area_ft2_[ui] } else { ctx.subcatches.area[ui] / ucf_landarea };
            self.soa_.subarea_runoff_rate[ui] = if full_area > 0.0 { sum / full_area } else { 0.0 };
        }

        {
            let mut q_i = runoff0 * f0 + runoff1 * f1;
            if route_mode == 2 && self.soa_.imperv_pct[ui] < 1.0 {
                q_i *= 1.0 - pct;
            }
            self.soa_.imperv_runoff_cfs[ui] = q_i * total_area;
            let mut q_p = runoff_p * fp;
            if route_mode == 1 && self.soa_.imperv_pct[ui] > 0.0 {
                q_p *= 1.0 - pct;
            }
            self.soa_.perv_runoff_cfs[ui] = q_p * total_area;
        }

        // ----- Step 4: loss rates and net runoff
        let evap_loss = if total_area > 0.0 { acc.vevap / dt / total_area } else { 0.0 };
        let infil_loss = if total_area > 0.0 { acc.vinfil / dt / total_area } else { 0.0 };
        let new_runoff = acc.voutflow / dt;

        self.soa_.runoff[ui] = new_runoff;
        self.soa_.outflow_vol[ui] = acc.voutflow;

        ctx.subcatches.runoff[ui] = new_runoff;
        ctx.subcatches.subarea_runoff_rate[ui] = self.soa_.subarea_runoff_rate[ui];
        {
            let fi = self.soa_.imperv_pct[ui];
            ctx.subcatches.ponded_depth[ui] = self.soa_.depth_imperv0[ui] * self.soa_.frac_imperv0[ui] + self.soa_.depth_imperv1[ui] * self.soa_.frac_imperv1[ui] + self.soa_.depth_perv[ui] * (1.0 - fi);
        }
        ctx.subcatches.evap_loss[ui] = evap_loss;
        ctx.subcatches.infil_loss[ui] = infil_loss;
        self.soa_.perv_evap_vol[ui] = acc.vpevap;
        self.soa_.actual_perv_evap_vol[ui] = acc.actual_perv_evap;
        self.soa_.infil_vol[ui] = acc.vinfil;
        self.soa_.spatial_infil_vol[ui] = 0.0; // Vspatial (no boundary on this path)

        ctx.subcatches.stat_evap_vol[ui] += acc.vevap;
        ctx.subcatches.stat_infil_vol[ui] += acc.vinfil;
        let area_i0 = total_area * f0;
        let area_i1 = total_area * f1;
        let area_pv = total_area * fp;
        ctx.subcatches.stat_imperv_vol[ui] += (runoff0 * area_i0 + runoff1 * area_i1) * dt;
        ctx.subcatches.stat_perv_vol[ui] += runoff_p * area_pv * dt;
    }

    /// processSubarea (legacy getSubareaRunoff) on the path without boundary or spatial
    /// fractions. `depth` is the subarea's ponded depth, copied in and out by the caller.
    #[allow(clippy::too_many_arguments)]
    fn process_subarea(&mut self, ctx: &mut SimulationContext, i: i32, depth: &mut f64, alpha: f64, d_store: f64, frac: f64, is_pervious: bool, runon_in: f64, subarea_n: f64, precip: f64, evap_rate: f64, dt: f64, infil_factor: f64, recovery_factor: f64, month: i32, total_area: f64, acc: &mut SubareaAcc) -> f64 {
        let ui = i as usize;
        if frac <= 0.0 {
            return 0.0;
        }
        let subarea_area = total_area * frac;

        let mut surf_moisture = *depth / dt;
        let surf_evap = cmin(surf_moisture, evap_rate);

        let mut infil = 0.0;
        if is_pervious {
            let runon = runon_in;
            let mut local_infil = infil_factor;
            if month >= 0 && ui < ctx.subcatch_infil_pattern.len() {
                let pi = ctx.subcatch_infil_pattern[ui];
                if pi >= 0 && (pi as usize) < ctx.patterns.factors.len() {
                    let facs = &ctx.patterns.factors[pi as usize];
                    let umon = month as usize;
                    if umon < facs.len() {
                        local_infil = facs[umon];
                    }
                }
            }

            self.infil_factor_used_[ui] = local_infil;
            let native = self.infil_get_infil(ctx, i, precip, runon, *depth, dt, local_infil, recovery_factor);
            self.soa_.native_infil_rate[ui] = native;
            infil = self.soa_.native_infil_rate[ui];
        }

        // Gap #40: limit pervious infiltration by GW upper zone capacity.
        if is_pervious && dt > 0.0 {
            let max_iv = ctx.subcatches.gw_max_infil_vol[ui];
            if max_iv < 1.0e30 {
                infil = cmin(infil, max_iv / dt);
            }
        }

        let inflow = precip + runon_in;
        surf_moisture += inflow;

        acc.vevap += surf_evap * subarea_area * dt;
        if is_pervious {
            acc.actual_perv_evap += surf_evap * subarea_area * dt;
        }
        if is_pervious {
            acc.vpevap += acc.vevap;
        }
        acc.vinfil += infil * subarea_area * dt;

        let mut runoff_rate = 0.0;
        let mut t_runoff = dt;
        if surf_evap + infil >= surf_moisture {
            *depth = 0.0;
        } else {
            let net_inflow = inflow - (surf_evap + infil);
            Self::update_ponded_depth(depth, net_inflow, alpha, d_store, dt, &mut t_runoff);
        }

        {
            let x_depth = *depth - d_store;
            if x_depth > 1.0e-10 {
                if subarea_n > 0.0 {
                    runoff_rate = alpha * x_depth.powf(MEXP); // LIBM
                } else {
                    runoff_rate = x_depth / t_runoff;
                    *depth = d_store;
                }
            }
        }

        let mut f_outlet = 1.0;
        let mut route_mode = ctx.subcatches.subarea_routing[ui];
        let pct = ctx.subcatches.pct_routed[ui];
        if self.soa_.imperv_pct[ui] == 0.0 || self.soa_.imperv_pct[ui] == 1.0 {
            route_mode = 0;
        }
        if route_mode == 2 && !is_pervious {
            f_outlet = 1.0 - pct;
        } else if route_mode == 1 && is_pervious {
            f_outlet = 1.0 - pct;
        }
        acc.voutflow += f_outlet * runoff_rate * subarea_area * dt;

        runoff_rate
    }
}
