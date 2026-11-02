// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/quality/QualityRouting.hpp and QualityRouting.cpp,
// FvQualityTransport.hpp and FvQualityTransport.cpp, FvTransportWater.hpp and
// ConservativeMixing.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The LEGACY quality router (namespace quality), init part: the QualitySolver and the data it
// owns (the FV quality-transport state and its conservative-mixing workspaces), the module
// constants, QualitySolver::init, and the one read SWMMEngine::initMassBalance makes through
// the state accessor open() installs (fvRetainedPollutantMass).
//
// Differences of form:
//   * ConservativeMixing's `fv::NativeGmres gmres_` workspace is not carried: it is a
//     run-time solver workspace (hydraulics/fv/NativeGmres.hpp) that nothing reads before the
//     first FV quality step, and it comes with the FV solver. NOT YET TRANSLATED.
//   * ConservativeMixing::Result is a nested struct upstream; here it is
//     ConservativeMixingResult.
//
// Run-time functions, translated with the quality phase (not here):
//   QualityRouting.cpp: linkTakesUpstreamValue, conduitMixingInflow, linkIsDry,
//     nodeIsReactor, nodeIsDry, storageEvapFactor, updateUntreatedStorageHrt,
//     applyNodeTreatment, QualitySolver::{addQualityOnlyInflows, assembleExternalLoads,
//     addCouplingLoads, addExtInflowLoads, execute, addWetWeatherLoads, addRdiiLoads,
//     addDwfLoads, addGwLoads, addIfaceLoads, accumulateLinkLoads, mixAtNodes, applyDecay,
//     updateLinkQuality, applyTreatment} and their file-local helpers.
//   FvQualityTransport.cpp: execute, applyPollutantDecay, dryRestartState,
//     restoreDryRestart, retainedSpeciesMass (and the inline hasRetainedDryMass,
//     beginWindow, accepted).
//   ConservativeMixing.cpp: ConservativeMixing::{prepare, solve}.

// ============================================================================
// Constants
// ============================================================================

/// 1 liter in ft3
pub const ZERO_VOLUME: f64 = 0.0353147;
/// 1 mm in ft
pub const ZERO_DEPTH: f64 = 0.003281;
/// Legacy's effective-zero flow (`ZERO` in src/legacy/engine/consts.h).
pub const LEGACY_ZERO: f64 = 1.0e-10;

// ============================================================================
// ConservativeMixing.hpp
// ============================================================================

/// One accepted internal water transfer, ft3. Donor and receiver are control volume
/// indices; the same transported amount must be used on both sides. (No member
/// initializers upstream: a value-initialized one is all zeros, as Default gives.)
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct WaterTransfer {
    pub donor: usize,
    pub receiver: usize,
    pub volume: f64,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct MixingWater {
    pub old_volume: Vec<f64>,
    pub new_volume: Vec<f64>,
    pub source_volume: Vec<f64>,
    pub exported: Vec<f64>,
    pub seeped: Vec<f64>,
    pub evaporated: Vec<f64>,
    pub transfers: Vec<WaterTransfer>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MixingScalar {
    pub old_concentration: Vec<f64>,
    pub source_mass: Vec<f64>,
    pub held_dry_mass: Vec<f64>,
    /// false for temperature-volume
    pub nonnegative: bool,
}

impl Default for MixingScalar {
    fn default() -> Self {
        MixingScalar { old_concentration: Vec::new(), source_mass: Vec::new(), held_dry_mass: Vec::new(), nonnegative: true }
    }
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct MixedScalar {
    pub concentration: Vec<f64>,
    pub held_dry_mass: Vec<f64>,
    pub transfer_mass: Vec<f64>,
    pub exported_mass: Vec<f64>,
    pub seeped_mass: Vec<f64>,
}

/// ConservativeMixing::Result
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ConservativeMixingResult {
    pub converged: bool,
    pub max_water_residual: f64,
    pub mass_residual: f64,
    pub iterations: i32,
}

/// Native conservative complete-mix transport on accepted water receipts (data only; see
/// the file notes for gmres_).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ConservativeMixing {
    pub water_: MixingWater,
    pub diagonal_: Vec<f64>,
    pub rhs_: Vec<f64>,
    pub image_: Vec<f64>,
    pub correction_: Vec<f64>,
    pub guess_: Vec<f64>,
    // gmres_: fv::NativeGmres — NOT YET TRANSLATED (see the file notes).
    pub candidate_: MixedScalar,
    pub incoming_: Vec<f64>,
    pub outgoing_: Vec<f64>,
    pub prepared_: bool,
    pub water_residual_: f64,
}

// ============================================================================
// FvTransportWater.hpp
// ============================================================================

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct BoundaryWaterImport {
    pub node: usize,
    pub receiver: usize,
    pub volume: f64,
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct FvTransportWater {
    pub water: MixingWater,
    pub boundary_imports: Vec<BoundaryWaterImport>,
    pub withdrawal: Vec<f64>,
    pub flooding: Vec<f64>,
    pub ready: bool,
}

// ============================================================================
// FvQualityTransport.hpp
// ============================================================================

#[derive(Debug, Clone, Default, PartialEq)]
pub struct FvQualityTransport {
    pub water: FvTransportWater,
    pub solute_mixer_: ConservativeMixing,
    pub thermal_mixer_: ConservativeMixing,
    pub held_: Vec<Vec<f64>>,
    pub results_: Vec<MixedScalar>,
    pub inputs_: Vec<MixingScalar>,
    pub thermal_water_: MixingWater,
    pub accepted_: bool,
}

impl FvQualityTransport {
    /// retainedPollutantMass: per pollutant, the sum (std::accumulate from 0.0, in order) of
    /// the held dry inventory row; zero for rows not held.
    pub fn retained_pollutant_mass(&self, n_pollutants: i32) -> Vec<f64> {
        let mut mass = vec![0.0f64; n_pollutants as usize];
        let mut p = 0usize;
        while p < mass.len() && p < self.held_.len() {
            mass[p] = self.held_[p].iter().fold(0.0, |a, &x| a + x);
            p += 1;
        }
        mass
    }
}

// ============================================================================
// QualitySolver
// ============================================================================

#[derive(Debug, Clone, Default, PartialEq)]
pub struct QualitySolver {
    pub n_pollutants_: i32,
    pub fv_transport_: FvQualityTransport,
    /// Per-node LATERAL carrier volume added into qual_vol_in this step.
    pub qual_vol_lat_: Vec<f64>,
    /// Negative lateral components seen this step (ft3).
    pub qual_vol_lat_neg_: Vec<f64>,
}

impl QualitySolver {
    /// QualitySolver::init: sets the pollutant row count and resets the FV transport state.
    /// The node and link counts are unused; qual_vol_lat_ and qual_vol_lat_neg_ are left as
    /// they are.
    pub fn init(&mut self, n_nodes: i32, n_links: i32, n_pollutants: i32) {
        self.n_pollutants_ = n_pollutants;
        self.fv_transport_ = FvQualityTransport::default();
        let _ = n_nodes;
        let _ = n_links;
    }

    /// fvRetainedPollutantMass(ctx): fv_transport_.retainedPollutantMass(ctx).
    pub fn fv_retained_pollutant_mass(&self, ctx: &crate::core::context::SimulationContext) -> Vec<f64> {
        self.fv_transport_.retained_pollutant_mass(ctx.n_pollutants())
    }
}
