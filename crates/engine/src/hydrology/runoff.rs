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
//   * NOT YET TRANSLATED (run-time step functions, phase 1c): updatePondedDepth,
//     pendingRoutingVolume, infilGetInfil, nativeInfilFullLid, reviewedAreas, execute,
//     checkpoint, validateCheckpoint, infil_get_state, infil_set_state.

use crate::core::cmath::cmax;
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
