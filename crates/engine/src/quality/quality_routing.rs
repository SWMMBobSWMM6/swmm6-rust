// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/quality/QualityRouting.hpp and QualityRouting.cpp,
// NegativeSources.hpp, FvQualityTransport.hpp and FvQualityTransport.cpp, FvTransportWater.hpp
// and ConservativeMixing.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026 (init part) and 5 November 2026 (the
// routing step's quality stage), by Claude (Anthropic) at the request of Robert Dickinson.
// Not an OpenSWMM release.
//
// The LEGACY quality router (namespace quality): the QualitySolver and the data it owns (the FV
// quality-transport state and its conservative-mixing workspaces), the module constants,
// QualitySolver::init, the read SWMMEngine::initMassBalance makes (fvRetainedPollutantMass),
// and the routing step's quality stage: the shared kernel rules (linkTakesUpstreamValue,
// conduitMixingInflow, linkIsDry, nodeIsReactor, storageEvapFactor,
// updateUntreatedStorageHrt, nodeIsDry), the file-local loadersNeeded, addAgeVolume,
// addTempVolume, nodeHasTreatmentFor and applyLinkQualityForcing, and
// QualitySolver::{execute, assembleExternalLoads, addQualityOnlyInflows, addCouplingLoads,
// addExtInflowLoads, addWetWeatherLoads, addRdiiLoads, addDwfLoads, addGwLoads, addIfaceLoads,
// accumulateLinkLoads, mixAtNodes, applyDecay, updateLinkQuality, applyTreatment} with
// applyNodeTreatment; NegativeSources.hpp's bookNegativeSourceClamp.
//
// Differences of form:
//   * ConservativeMixing's `fv::NativeGmres gmres_` workspace is not carried: it is a
//     run-time solver workspace (hydraulics/fv/NativeGmres.hpp) that nothing reads before the
//     first FV quality step, and it comes with the FV solver. NOT YET TRANSLATED.
//   * ConservativeMixing::Result is a nested struct upstream; here it is
//     ConservativeMixingResult.
//   * execute returns Err(StepError) where upstream would throw (std::stod inside a node-LID
//     layer treatment expression) or where it would enter a path not translated yet.
//   * applyNodeTreatment's `cin` (a const double* that may alias treatment.cin) is
//     Option<&[f64]>: None is the aliasing call the LEGACY pass makes.
//
// NOT YET TRANSLATED (each returns StepError::NotTranslated where upstream would first act;
// SWMMEngine::stepRouting refuses such a deck before any state changes, so none is reached):
//   the FV branch of execute (FvQualityTransport::execute, applyPollutantDecay, the bed
//   exchange), msxlid::prepareNodes / reactNodes (reaction-species LID controls),
//   heat::applyBedSoluteExchange. (routeLegacyAge, routeLegacyHeat, routeLegacyMsx and
//   reactLegacyNodes / reactLegacyLinks are in transport/water_age_legacy.rs,
//   heat_legacy.rs and reaction_legacy_step.rs.)
//   FvQualityTransport.cpp: execute, applyPollutantDecay, dryRestartState,
//     restoreDryRestart, retainedSpeciesMass (and the inline hasRetainedDryMass,
//     beginWindow, accepted).
//   ConservativeMixing.cpp: ConservativeMixing::{prepare, solve}.

use crate::core::cmath::{cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::engine_step_routing::StepError;
use crate::core::options::RoutingModel;
use crate::core::unit_conversion as ucf;
use crate::data::forcing::ForcingMode;
use crate::data::heat::HeatSource;
use crate::data::link::{LinkType, XsectShape};
use crate::data::node::NodeType;
use crate::data::treatment_types::TreatVar;
use crate::data::water_age::WaterAgeSource;
use crate::hydraulics::node;
use crate::hydrology::lid_node_quality as lidq;
use crate::hydrology::lid_node_step as lidnode;
use crate::quality::treatment;
use crate::transport::ard_config::ard_boundaries_need_external_volumes;
use crate::transport::reaction_legacy_binding::legacy_reactions_active;

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

// ============================================================================
// The shared shape of the legacy mixing kernels (QualityRouting.hpp)
// ============================================================================

/// linkTakesUpstreamValue: a non-conduit link or a DUMMY conduit holds no water.
pub fn link_takes_upstream_value(ctx: &SimulationContext, link: i32) -> bool {
    let uj = link as usize;
    ctx.links.type_[uj] != LinkType::Conduit || ctx.links.xsect_shape[uj] == XsectShape::Dummy
}

/// conduitMixingInflow: the mixing inflow rate a conduit sees over one step.
pub fn conduit_mixing_inflow(ctx: &SimulationContext, link: i32, dt: f64) -> f64 {
    let uj = link as usize;
    let model = ctx.options.routing_model;
    let mut q_in = (if model == RoutingModel::Kinwave { ctx.links.kw_inflow[uj] } else { ctx.links.flow[uj] }).abs();

    if model != RoutingModel::Dynwave && model != RoutingModel::Fv {
        return cmax(q_in, 0.0);
    }
    if dt <= 0.0 {
        return cmax(q_in, 0.0);
    }

    let cr = ctx.link_subtypes.conduit_row(link);
    let mut v_losses = 0.0;
    if cr >= 0 {
        let ucr = cr as usize;
        let cd = &ctx.link_subtypes.conduits;
        let barrels = cd.barrels[ucr] as f64;
        v_losses = cd.seep_loss_rate[ucr] * barrels * dt + cd.evap_loss_rate[ucr] * barrels * dt;
    }
    q_in += (ctx.links.volume[uj] + v_losses - ctx.links.old_volume[uj]) / dt;
    cmax(q_in, 0.0)
}

/// linkIsDry: under a litre of water or under a millimetre of depth.
pub fn link_is_dry(ctx: &SimulationContext, link: i32) -> bool {
    let uj = link as usize;
    ctx.links.volume[uj] < ZERO_VOLUME || ctx.links.depth[uj] <= ZERO_DEPTH
}

/// nodeIsReactor: a STORAGE node, or one that started the step holding more than a litre.
pub fn node_is_reactor(ctx: &SimulationContext, node: i32) -> bool {
    let ui = node as usize;
    ctx.nodes.type_[ui] == NodeType::Storage || ctx.nodes.old_volume[ui] > ZERO_VOLUME
}

/// storageEvapFactor: 1 + vEvap/v1 for a storage node, 1 for any other.
pub fn storage_evap_factor(ctx: &SimulationContext, node: i32) -> f64 {
    let mut f_evap = 1.0;
    let ui = node as usize;
    if ui >= ctx.nodes.type_.len() || ctx.nodes.type_[ui] != NodeType::Storage {
        return f_evap;
    }
    let srow = ctx.node_subtypes.storage_row(node);
    if srow < 0 {
        return f_evap;
    }
    let us = srow as usize;
    let st = &ctx.node_subtypes.storages;
    let v_old = ctx.nodes.old_volume[ui];
    let v_evap = if us < st.evap_loss.len() { cmax(0.0, st.evap_loss[us] - st.lid_state[us].evap_volume) } else { 0.0 };
    if v_evap > 0.0 && v_old > ZERO_VOLUME {
        f_evap += v_evap / v_old;
    }
    f_evap
}

/// updateUntreatedStorageHrt: residence time at untreated storage nodes (reactions only).
pub fn update_untreated_storage_hrt(ctx: &mut SimulationContext, dt: f64) {
    if dt <= 0.0 || !legacy_reactions_active(ctx) {
        return;
    }
    let any_treatment = ctx.n_pollutants() > 0 && ctx.treatment.has_any();
    for j in 0..ctx.n_nodes() {
        let uj = j as usize;
        if ctx.nodes.type_[uj] != NodeType::Storage || uj >= ctx.nodes.hrt.len() {
            continue;
        }
        if any_treatment && uj < ctx.treatment.has_treatment.len() && ctx.treatment.has_treatment[uj] {
            continue;
        }
        let v_old = ctx.nodes.old_volume[uj];
        if v_old < LEGACY_ZERO {
            ctx.nodes.hrt[uj] = 0.0;
        } else {
            let q = if uj < ctx.nodes.inflow.len() { ctx.nodes.inflow[uj].abs() } else { 0.0 };
            ctx.nodes.hrt[uj] = cmax((ctx.nodes.hrt[uj] + dt) * v_old / (v_old + q * dt), 0.0);
        }
    }
}

/// nodeIsDry: essentially empty AND taking (almost) nothing in.
pub fn node_is_dry(ctx: &SimulationContext, node: i32, q_in: f64) -> bool {
    let ui = node as usize;
    (ctx.nodes.volume[ui] <= ZERO_VOLUME || ctx.nodes.depth[ui] <= ZERO_DEPTH) && q_in <= LEGACY_ZERO
}

// ---------------------------------------------------------------------------
// File-local helpers (QualityRouting.cpp's anonymous namespace)
// ---------------------------------------------------------------------------

/// loadersNeeded: any consumer of the external loads or carrier volumes.
fn loaders_needed(np: i32, ctx: &SimulationContext) -> bool {
    np > 0 || ard_boundaries_need_external_volumes(ctx) || ctx.options.water_age || ctx.options.heat_transport || legacy_reactions_active(ctx)
}

/// addAgeVolume: q * source age into the node's age-volume rate (WATER_AGE only).
fn add_age_volume(ctx: &mut SimulationContext, node: i32, q: f64, src: WaterAgeSource) {
    if !ctx.options.water_age {
        return;
    }
    let un = node as usize;
    if un >= ctx.water_age_state.node_age_vol_in.len() {
        return;
    }
    let mut age = ctx.water_age_config.source_age(src, node);
    if src == WaterAgeSource::ExternalInflow {
        let ov = &ctx.water_age_state.node_ext_inflow_age;
        if un < ov.len() && !ov[un].is_nan() {
            age = ov[un];
        }
    }
    ctx.water_age_state.node_age_vol_in[un] += q * age;
}

/// addTempVolume: q * source temperature into the node's temperature-volume rate.
fn add_temp_volume(ctx: &mut SimulationContext, node: i32, q: f64, src: HeatSource) {
    if !ctx.options.heat_transport {
        return;
    }
    let un = node as usize;
    if un >= ctx.heat_state.node_temp_vol_in.len() {
        return;
    }
    let t = ctx.heat_config.source_temp(src, node);
    ctx.heat_state.node_temp_vol_in[un] += q * t;
}

/// nodeHasTreatmentFor: a non-empty treatment expression for (node, p).
fn node_has_treatment_for(ctx: &SimulationContext, node: i32, p: i32) -> bool {
    let tr = &ctx.treatment;
    if tr.n_pollutants <= 0 {
        return false;
    }
    let idx = node as usize * tr.n_pollutants as usize + p as usize;
    idx < tr.expressions.len() && !tr.expressions[idx].is_empty()
}

/// applyLinkQualityForcing: runtime link concentration OVERRIDE / ADD.
fn apply_link_quality_forcing(ctx: &mut SimulationContext, n_pollutants: i32, dt: f64) {
    if n_pollutants <= 0 {
        return;
    }
    if ctx.forcing.link_quality_mode.is_empty() {
        return;
    }
    for j in 0..ctx.n_links() {
        let uj = j as usize;
        for p in 0..n_pollutants {
            let flat = uj * n_pollutants as usize + p as usize;
            if flat >= ctx.forcing.link_quality_mode.len() || flat >= ctx.links.conc.len() {
                continue;
            }
            if ctx.forcing.link_quality_mode[flat] == ForcingMode::Override {
                ctx.links.conc[flat] = cmax(ctx.forcing.link_quality_value[flat], 0.0);
            } else if ctx.forcing.link_quality_mode[flat] == ForcingMode::Add {
                ctx.links.conc[flat] = cmax(ctx.links.conc[flat] + ctx.forcing.link_quality_value[flat], 0.0);
                ctx.mass_balance.routing_forcing_qual_inflow[p as usize] += ctx.forcing.link_quality_value[flat] * dt;
            }
        }
    }
}

/// NegativeSources.hpp bookNegativeSourceClamp: count an extraction clamp and return the unmet
/// mass to the external-inflow row it was booked into.
pub fn book_negative_source_clamp(ctx: &mut SimulationContext, node: i32, p: i32, shortfall: f64) {
    let st = &mut ctx.negsrc;
    st.clamp_events += 1;
    st.shortfall_mass += shortfall;
    if p >= 0 && (p as usize) < ctx.mass_balance.qual_routing_ex_in.len() {
        ctx.mass_balance.qual_routing_ex_in[p as usize] += shortfall;
    }
    let st = &mut ctx.negsrc;
    if !st.first_clamp_recorded {
        st.first_clamp_recorded = true;
        st.first_node = node;
    }
}

/// NegativeSources.hpp bookNegativeCellSourceClamp: count an ARD cell-source extraction clamp
/// (counted only; no ledger row).
pub fn book_negative_cell_source_clamp(ctx: &mut SimulationContext, conduit: i32, shortfall: f64) {
    let st = &mut ctx.negsrc;
    st.clamp_events += 1;
    st.shortfall_mass += shortfall;
    if !st.first_clamp_recorded {
        st.first_clamp_recorded = true;
        st.first_node = conduit;
    }
}

/// A run-time path not translated yet (see the file notes).
fn nt(what: &'static str) -> StepError {
    StepError::NotTranslated(what)
}

impl QualitySolver {
    /// addQualityOnlyInflows: surface-outfall and gaining-conduit receipts (water that already
    /// arrived hydraulically, with the mass its donor gave up).
    pub fn add_quality_only_inflows(&mut self, ctx: &mut SimulationContext, dt: f64) {
        if !(dt > 0.0) {
            return;
        }
        let mut q = std::mem::take(&mut ctx.surface_outfall_quality);
        apply_quality_only(ctx, &mut q, dt, true);
        ctx.surface_outfall_quality = q;
        let mut q = std::mem::take(&mut ctx.link_gain_quality);
        apply_quality_only(ctx, &mut q, dt, false);
        ctx.link_gain_quality = q;
    }

    /// assembleExternalLoads: reset the node assembly arrays and run every external loader.
    pub fn assemble_external_loads(&mut self, ctx: &mut SimulationContext, dt: f64) {
        if !loaders_needed(self.n_pollutants_, ctx) {
            return;
        }

        ctx.nodes.qual_mass_in.iter_mut().for_each(|x| *x = 0.0);

        if ctx.options.water_age {
            let (nn, nl, ns) = (ctx.n_nodes(), ctx.n_links(), ctx.n_subcatches());
            let ws = &mut ctx.water_age_state;
            if ws.node_age_vol_in.len() != nn as usize {
                ws.resize(nn, nl, ns);
            }
            ws.node_age_vol_in.iter_mut().for_each(|x| *x = 0.0);
        }
        if ctx.options.heat_transport {
            let (nn, nl) = (ctx.n_nodes(), ctx.n_links());
            let t0 = ctx.heat_config.global_temp[HeatSource::InitialState as usize];
            let hs = &mut ctx.heat_state;
            if hs.node_temp_vol_in.len() != nn as usize {
                hs.resize(nn, nl, t0);
            }
            hs.node_temp_vol_in.iter_mut().for_each(|x| *x = 0.0);
        }
        ctx.nodes.qual_vol_in.iter_mut().for_each(|x| *x = 0.0);
        let nn = ctx.n_nodes() as usize;
        self.qual_vol_lat_.clear();
        self.qual_vol_lat_.resize(nn, 0.0);
        self.qual_vol_lat_neg_.clear();
        self.qual_vol_lat_neg_.resize(nn, 0.0);

        self.add_wet_weather_loads(ctx, dt);
        self.add_rdii_loads(ctx, dt);
        self.add_dwf_loads(ctx, dt);
        self.add_gw_loads(ctx, dt);
        self.add_iface_loads(ctx, dt);
        self.add_ext_inflow_loads(ctx, dt);
        self.add_coupling_loads(ctx, dt);
        self.add_quality_only_inflows(ctx, dt);

        if ctx.options.routing_model == RoutingModel::Fv {
            for i in 0..ctx.n_nodes() {
                let ui = i as usize;
                let mut q = ctx.nodes.user_lat_flow[ui];
                if ui < ctx.forcing.node_lat_inflow_mode.len() {
                    if ctx.forcing.node_lat_inflow_mode[ui] == ForcingMode::Override {
                        q = ctx.forcing.node_lat_inflow_value[ui];
                    } else if ctx.forcing.node_lat_inflow_mode[ui] == ForcingMode::Add {
                        q += ctx.forcing.node_lat_inflow_value[ui];
                    }
                }
                if q > 0.0 {
                    ctx.nodes.qual_vol_in[ui] += q * dt;
                    add_age_volume(ctx, i, q, WaterAgeSource::ExternalInflow);
                    add_temp_volume(ctx, i, q, HeatSource::ExternalInflow);
                }
            }
            return;
        }

        // legacy's floored net lateral denominator (routing.c:480), only where a negative
        // lateral component arrived
        let nodes = &mut ctx.nodes;
        let nn = cmin(self.qual_vol_lat_.len(), nodes.qual_vol_in.len());
        for ui in 0..nn {
            if !(ui < self.qual_vol_lat_neg_.len() && self.qual_vol_lat_neg_[ui] < 0.0) {
                continue;
            }
            let coupling = if ui < nodes.coupling_inflow.len() { nodes.coupling_inflow[ui] } else { 0.0 };
            let lat_net_pos = if ui < nodes.lat_flow.len() { cmax(0.0, nodes.lat_flow[ui] - coupling) * dt } else { 0.0 };
            let corr = lat_net_pos - self.qual_vol_lat_[ui];
            if corr < 0.0 {
                nodes.qual_vol_in[ui] += corr;
            }
        }
    }

    /// addCouplingLoads: the 2D surface's junction drain (water and species mass).
    pub fn add_coupling_loads(&mut self, ctx: &mut SimulationContext, dt: f64) {
        let np = self.n_pollutants_;
        if !loaders_needed(np, ctx) {
            return;
        }
        if ctx.nodes.coupling_inflow.is_empty() {
            return;
        }
        for i in 0..ctx.n_nodes() {
            let ui = i as usize;
            let q = ctx.nodes.coupling_inflow[ui];
            if q > 0.0 {
                ctx.nodes.qual_vol_in[ui] += q * dt;
                if ctx.nodes.coupling_tuple_age && ctx.options.water_age && ui < ctx.water_age_state.node_age_vol_in.len() && ui < ctx.nodes.coupling_age_vol_inflow.len() {
                    ctx.water_age_state.node_age_vol_in[ui] += ctx.nodes.coupling_age_vol_inflow[ui];
                } else {
                    add_age_volume(ctx, i, q, WaterAgeSource::ExternalInflow);
                }
                if ctx.nodes.coupling_tuple_temp && ctx.options.heat_transport && ui < ctx.heat_state.node_temp_vol_in.len() && ui < ctx.nodes.coupling_temp_vol_inflow.len() {
                    ctx.heat_state.node_temp_vol_in[ui] += ctx.nodes.coupling_temp_vol_inflow[ui];
                } else {
                    add_temp_volume(ctx, i, q, HeatSource::ExternalInflow);
                }
            }
            if ctx.nodes.coupling_qual_inflow.is_empty() {
                continue;
            }
            for p in 0..np {
                let nd_idx = ui * np as usize + p as usize;
                if nd_idx >= ctx.nodes.coupling_qual_inflow.len() {
                    continue;
                }
                let mass_rate = ctx.nodes.coupling_qual_inflow[nd_idx];
                if !(mass_rate > 0.0) {
                    continue;
                }
                if nd_idx < ctx.nodes.qual_mass_in.len() {
                    ctx.nodes.qual_mass_in[nd_idx] += mass_rate;
                }
                let pi = p as usize;
                if pi < ctx.mass_balance.qual_routing_ex_in.len() {
                    ctx.mass_balance.qual_routing_ex_in[pi] += mass_rate * dt;
                }
            }
        }
    }

    /// addExtInflowLoads: direct [INFLOWS] water and CONCEN/MASS loads, and the runtime-API
    /// forced quality mass flux.
    pub fn add_ext_inflow_loads(&mut self, ctx: &mut SimulationContext, dt: f64) {
        let np = self.n_pollutants_;
        if !loaders_needed(np, ctx) {
            return;
        }
        for i in 0..ctx.n_nodes() {
            let ui = i as usize;
            let q = ctx.nodes.ext_inflow[ui];
            if q < 0.0 && ui < self.qual_vol_lat_neg_.len() {
                self.qual_vol_lat_neg_[ui] += q * dt;
            }
            if q > 0.0 {
                ctx.nodes.qual_vol_in[ui] += q * dt;
                self.qual_vol_lat_[ui] += q * dt;
                add_age_volume(ctx, i, q, WaterAgeSource::ExternalInflow);
                add_temp_volume(ctx, i, q, HeatSource::ExternalInflow);
            }

            if ctx.nodes.ext_qual_mass.is_empty() {
                continue;
            }
            for p in 0..np {
                let nd_idx = ui * np as usize + p as usize;
                if nd_idx >= ctx.nodes.ext_qual_mass.len() {
                    continue;
                }
                let mass_rate = ctx.nodes.ext_qual_mass[nd_idx];
                if mass_rate == 0.0 {
                    continue;
                }
                if nd_idx < ctx.nodes.qual_mass_in.len() {
                    ctx.nodes.qual_mass_in[nd_idx] += mass_rate;
                }
                let pi = p as usize;
                if pi < ctx.mass_balance.qual_routing_ex_in.len() {
                    ctx.mass_balance.qual_routing_ex_in[pi] += mass_rate * dt;
                }
            }
        }

        if !ctx.nodes.user_conc_mass_flux.is_empty() {
            for i in 0..ctx.n_nodes() {
                let ui = i as usize;
                for p in 0..np {
                    let nd_idx = ui * np as usize + p as usize;
                    if nd_idx >= ctx.nodes.user_conc_mass_flux.len() {
                        continue;
                    }
                    let w = ctx.nodes.user_conc_mass_flux[nd_idx];
                    if w == 0.0 {
                        continue;
                    }
                    if w < 0.0 && !ctx.negsrc.api_warned {
                        ctx.negsrc.api_warned = true;
                        ctx.warnings.push(b"D-NS1: a negative quality mass flux was applied via the runtime API (extraction). It is clamped per step to the mass the element holds.".to_vec());
                    }
                    if nd_idx < ctx.nodes.qual_mass_in.len() {
                        ctx.nodes.qual_mass_in[nd_idx] += w;
                    }
                    let pi = p as usize;
                    if pi < ctx.mass_balance.qual_routing_ex_in.len() {
                        ctx.mass_balance.qual_routing_ex_in[pi] += w * dt;
                    }
                    if pi < ctx.mass_balance.routing_forcing_qual_inflow.len() {
                        ctx.mass_balance.routing_forcing_qual_inflow[pi] += w * dt;
                    }
                }
            }
        }
    }

    /// QualitySolver::execute: one quality routing step (LEGACY).
    pub fn execute(&mut self, ctx: &mut SimulationContext, dt: f64) -> Result<(), StepError> {
        if self.n_pollutants_ <= 0 && !legacy_reactions_active(ctx) && !ctx.options.water_age && !ctx.options.heat_transport {
            return Ok(());
        }

        self.assemble_external_loads(ctx, dt);
        if ctx.options.routing_model == RoutingModel::Fv {
            // fv_transport_.execute(ctx, dt) and what follows it
            return Err(nt("FV quality transport (FvQualityTransport::execute)"));
        }
        lidq::prepare_quality(ctx, dt)?;
        // msxlid::prepareNodes
        if !ctx.reactions.lid_controls.is_empty() && dt > 0.0 {
            return Err(nt("reaction species LID controls"));
        }
        let coupled_lids = self.n_pollutants_ > 0 && ctx.node_subtypes.storages.lid_state.iter().any(|state| !state.cells.is_empty());
        if coupled_lids {
            let mass_in = ctx.nodes.qual_mass_in.clone();
            let vol_in = ctx.nodes.qual_vol_in.clone();
            let reacted = ctx.mass_balance.qual_routing_reacted.clone();
            let seeped = ctx.mass_balance.qual_routing_seep.clone();
            let dry = ctx.mass_balance.qual_routing_final_dry.clone();
            let external = ctx.mass_balance.qual_routing_ex_in.clone();
            let negative_sources = ctx.negsrc.clone();
            let mut held: Vec<Vec<f64>> = Vec::new();
            for state in ctx.node_subtypes.storages.lid_state.iter() {
                held.push(state.quality_mass.clone());
            }
            ctx.nodes.conc = ctx.nodes.conc_old.clone();
            let mut converged = false;
            for _iteration in 0..100 {
                let previous = ctx.nodes.conc.clone();
                ctx.nodes.qual_mass_in = mass_in.clone();
                ctx.nodes.qual_vol_in = vol_in.clone();
                ctx.mass_balance.qual_routing_reacted = reacted.clone();
                ctx.mass_balance.qual_routing_seep = seeped.clone();
                ctx.mass_balance.qual_routing_final_dry = dry.clone();
                ctx.mass_balance.qual_routing_ex_in = external.clone();
                ctx.negsrc = negative_sources.clone();
                for r in 0..held.len() {
                    ctx.node_subtypes.storages.lid_state[r].quality_mass = held[r].clone();
                }
                lidq::prepare_outlet_quality(ctx, dt, false)?;
                self.accumulate_link_loads(ctx, dt);
                self.mix_at_nodes(ctx, dt);
                converged = true;
                for i in 0..previous.len() {
                    let c = ctx.nodes.conc[i];
                    if (previous[i] - c).abs() > 1.0e-11 * cmax(cmax(1.0, previous[i].abs()), c.abs()) {
                        converged = false;
                        break;
                    }
                }
                if converged {
                    break;
                }
            }
            if !converged {
                let warning: &[u8] = b"LID quality mixtures did not converge; inspect quality continuity and reduce the routing step";
                if !ctx.warnings.iter().any(|w| w.as_slice() == warning) {
                    ctx.warnings.push(warning.to_vec());
                }
            }
            // Iteration is provisional: layer reactions are booked only once.
            lidq::prepare_outlet_quality(ctx, dt, true)?;
        } else {
            self.accumulate_link_loads(ctx, dt);
            self.mix_at_nodes(ctx, dt);
        }
        update_untreated_storage_hrt(ctx, dt); // review A6
        // R4b: MSX element state advects on the CSTR mirror family, before the react stages
        {
            let held = crate::transport::reaction_legacy_step::fv_retained_species_mass(&self.fv_transport_.held_, ctx.n_pollutants(), ctx.reactions.n_species());
            crate::transport::reaction_legacy_step::route_legacy_msx(ctx, dt, &held)?;
        }
        self.apply_treatment(ctx, dt, true);
        if legacy_reactions_active(ctx) {
            crate::transport::reaction_legacy_step::react_legacy_nodes(ctx, dt)?;
            // msxlid::reactNodes: a no-op without reaction-species LID controls (refused above)
            if !ctx.reactions.lid_controls.is_empty() {
                return Err(nt("reaction species LID controls"));
            }
        } else {
            self.apply_decay(ctx, dt);
        }
        self.update_link_quality(ctx, dt);
        if legacy_reactions_active(ctx) {
            crate::transport::reaction_legacy_step::react_legacy_links(ctx, dt)?;
        }
        apply_link_quality_forcing(ctx, self.n_pollutants_, dt);

        // A1b: the LEGACY age mirror runs last, on the fully accumulated qual_vol_in
        crate::transport::water_age_legacy::route_legacy_age(ctx, dt);
        // H1: the LEGACY temperature mirror beside it
        crate::transport::heat_legacy::route_legacy_heat(ctx, dt)?;
        // H6b: the bed's solute exchange
        if crate::transport::heat_legacy::bed_exchange_enabled(ctx) {
            return Err(nt("bed solute exchange (heat::applyBedSoluteExchange)"));
        }
        Ok(())
    }

    /// addWetWeatherLoads: subcatchment runoff and its washoff load at the outlet nodes, at
    /// legacy's interpolation point between the runoff times; then the LID drain loads.
    pub fn add_wet_weather_loads(&mut self, ctx: &mut SimulationContext, dt: f64) {
        let np = self.n_pollutants_;

        for i in 0..ctx.n_subcatches() {
            let ui = i as usize;
            let out_node = ctx.subcatches.outlet_node[ui];
            if out_node < 0 || out_node >= ctx.n_nodes() {
                continue;
            }
            let ud = out_node as usize;

            let f = ctx.runoff_interp_f;
            let f1 = 1.0 - f;
            let q_new = ctx.subcatches.runoff[ui];
            let q_old = ctx.subcatches.old_runoff[ui];
            let q = f1 * q_old + f * q_new;
            if q <= 0.0 {
                continue;
            }

            ctx.nodes.qual_vol_in[ud] += q * dt;
            if ud < self.qual_vol_lat_.len() {
                self.qual_vol_lat_[ud] += q * dt;
            }
            {
                let have_sc = ctx.options.water_age && ui < ctx.water_age_state.subcatch_runoff_age.len();
                if have_sc {
                    let a = ctx.water_age_state.subcatch_runoff_age[ui];
                    let acc = &mut ctx.water_age_state.node_age_vol_in;
                    if ud < acc.len() {
                        acc[ud] += q * a;
                    }
                } else {
                    add_age_volume(ctx, out_node, q, WaterAgeSource::Rainfall);
                }
            }
            {
                let have_sc = ctx.options.heat_transport && ui < ctx.heat_state.subcatch_runoff_temp.len();
                if have_sc {
                    let t = ctx.heat_state.subcatch_runoff_temp[ui];
                    let acc = &mut ctx.heat_state.node_temp_vol_in;
                    if ud < acc.len() {
                        acc[ud] += q * t;
                    }
                } else {
                    add_temp_volume(ctx, out_node, q, HeatSource::Rainfall);
                }
            }

            for p in 0..np {
                let sc_idx = ui * np as usize + p as usize;
                let nd_idx = ud * np as usize + p as usize;

                let c_new = if sc_idx < ctx.subcatches.conc.len() { ctx.subcatches.conc[sc_idx] } else { 0.0 };
                let c_old = if sc_idx < ctx.subcatches.conc_old.len() { ctx.subcatches.conc_old[sc_idx] } else { 0.0 };

                let mass_rate = f1 * (q_old * c_old) + f * (q_new * c_new);

                if mass_rate > 0.0 && nd_idx < ctx.nodes.qual_mass_in.len() {
                    ctx.nodes.qual_mass_in[nd_idx] += mass_rate;
                }
                let pi = p as usize;
                if mass_rate > 0.0 && pi < ctx.mass_balance.qual_routing_wet.len() {
                    ctx.mass_balance.qual_routing_wet[pi] += mass_rate * dt;
                }
            }
        }

        // Gap #26: LID drain quality
        let f_dr = ctx.runoff_interp_f;
        let f1_dr = 1.0 - f_dr;
        for j in 0..ctx.n_nodes() {
            let uj = j as usize;
            let v_new = if uj < ctx.nodes.lid_drain_qual_vol.len() { ctx.nodes.lid_drain_qual_vol[uj] } else { 0.0 };
            let v_old = if uj < ctx.nodes.lid_drain_qual_vol_old.len() { ctx.nodes.lid_drain_qual_vol_old[uj] } else { 0.0 };
            let drain_vol_rate = f1_dr * v_old + f_dr * v_new;
            if drain_vol_rate <= 0.0 {
                continue;
            }

            ctx.nodes.qual_vol_in[uj] += drain_vol_rate * dt;
            if uj < self.qual_vol_lat_.len() {
                self.qual_vol_lat_[uj] += drain_vol_rate * dt;
            }
            if ctx.options.water_age && uj < ctx.water_age_state.node_lid_drain_age_vol_in.len() && ctx.lid_layer_state.active() {
                let a = ctx.water_age_state.node_lid_drain_age_vol_in[uj];
                let acc = &mut ctx.water_age_state.node_age_vol_in;
                if uj < acc.len() {
                    acc[uj] += a;
                }
            } else {
                add_age_volume(ctx, j, drain_vol_rate, WaterAgeSource::Rainfall);
            }
            if ctx.options.heat_transport && uj < ctx.heat_state.node_lid_drain_temp_vol_in.len() && ctx.lid_layer_state.active() {
                let t = ctx.heat_state.node_lid_drain_temp_vol_in[uj];
                let tacc = &mut ctx.heat_state.node_temp_vol_in;
                if uj < tacc.len() {
                    tacc[uj] += t;
                }
            } else {
                add_temp_volume(ctx, j, drain_vol_rate, HeatSource::Rainfall);
            }

            for p in 0..np {
                let nd_idx = uj * np as usize + p as usize;
                let w_new = if nd_idx < ctx.nodes.lid_drain_qual_load.len() { ctx.nodes.lid_drain_qual_load[nd_idx] } else { 0.0 };
                let w_old = if nd_idx < ctx.nodes.lid_drain_qual_load_old.len() { ctx.nodes.lid_drain_qual_load_old[nd_idx] } else { 0.0 };
                let load = f1_dr * w_old + f_dr * w_new;
                if load > 0.0 && nd_idx < ctx.nodes.qual_mass_in.len() {
                    ctx.nodes.qual_mass_in[nd_idx] += load;
                }
                let pi = p as usize;
                if load > 0.0 && pi < ctx.mass_balance.qual_routing_wet.len() {
                    ctx.mass_balance.qual_routing_wet[pi] += load * dt;
                }
            }
        }
    }

    /// addRdiiLoads: RDII water and its pollutant loads (q * c_rdii).
    fn add_rdii_loads(&mut self, ctx: &mut SimulationContext, dt: f64) {
        let np = self.n_pollutants_;
        if !loaders_needed(np, ctx) {
            return;
        }
        for i in 0..ctx.n_nodes() {
            let ui = i as usize;
            let q = ctx.nodes.rdii_inflow[ui];
            if q <= 0.0 {
                continue;
            }
            ctx.nodes.qual_vol_in[ui] += q * dt;
            self.qual_vol_lat_[ui] += q * dt;
            add_age_volume(ctx, i, q, WaterAgeSource::Rdii);
            add_temp_volume(ctx, i, q, HeatSource::Rdii);

            for p in 0..np {
                let c_rdii = ctx.pollutants.c_rdii[p as usize];
                if c_rdii <= 0.0 {
                    continue;
                }
                let nd_idx = ui * np as usize + p as usize;
                if nd_idx < ctx.nodes.qual_mass_in.len() {
                    ctx.nodes.qual_mass_in[nd_idx] += q * c_rdii;
                }
            }
            for p in 0..np {
                let c_rdii = ctx.pollutants.c_rdii[p as usize];
                if c_rdii <= 0.0 {
                    continue;
                }
                let pi = p as usize;
                if pi < ctx.mass_balance.qual_routing_ii_in.len() {
                    ctx.mass_balance.qual_routing_ii_in[pi] += q * c_rdii * dt;
                }
            }
        }
    }

    /// addDwfLoads: dry weather water, the default concentration and the per-node [DWF]
    /// pollutant adjustment.
    fn add_dwf_loads(&mut self, ctx: &mut SimulationContext, dt: f64) {
        let np = self.n_pollutants_;
        if !loaders_needed(np, ctx) {
            return;
        }
        for i in 0..ctx.n_nodes() {
            let ui = i as usize;
            let q = ctx.nodes.dwf_inflow[ui];
            if q < 0.0 && ui < self.qual_vol_lat_neg_.len() {
                self.qual_vol_lat_neg_[ui] += q * dt;
            }
            if q <= 0.0 {
                continue;
            }
            ctx.nodes.qual_vol_in[ui] += q * dt;
            self.qual_vol_lat_[ui] += q * dt;
            add_age_volume(ctx, i, q, WaterAgeSource::Dwf);
            add_temp_volume(ctx, i, q, HeatSource::Dwf);

            for p in 0..np {
                let c_dwf = ctx.pollutants.c_dwf[p as usize];
                let nd_idx = ui * np as usize + p as usize;
                let w_adj = if nd_idx < ctx.nodes.dwf_qual_mass.len() { ctx.nodes.dwf_qual_mass[nd_idx] } else { 0.0 };
                let w = (if c_dwf > 0.0 { q * c_dwf } else { 0.0 }) + w_adj;
                if w == 0.0 {
                    continue;
                }
                if nd_idx < ctx.nodes.qual_mass_in.len() {
                    ctx.nodes.qual_mass_in[nd_idx] += w;
                }
                let pi = p as usize;
                if pi < ctx.mass_balance.qual_routing_dw_in.len() {
                    ctx.mass_balance.qual_routing_dw_in[pi] += w * dt;
                }
            }
        }
    }

    /// addGwLoads: groundwater inflow and its pollutant loads (q * c_gw).
    fn add_gw_loads(&mut self, ctx: &mut SimulationContext, dt: f64) {
        let np = self.n_pollutants_;
        if !loaders_needed(np, ctx) {
            return;
        }
        for i in 0..ctx.n_nodes() {
            let ui = i as usize;
            let q = ctx.nodes.gw_inflow[ui];
            if q < 0.0 && ui < self.qual_vol_lat_neg_.len() {
                self.qual_vol_lat_neg_[ui] += q * dt;
            }
            if q <= 0.0 {
                continue;
            }
            ctx.nodes.qual_vol_in[ui] += q * dt;
            self.qual_vol_lat_[ui] += q * dt;
            add_age_volume(ctx, i, q, WaterAgeSource::Gw);
            add_temp_volume(ctx, i, q, HeatSource::Gw);

            for p in 0..np {
                let c_gw = ctx.pollutants.c_gw[p as usize];
                if c_gw <= 0.0 {
                    continue;
                }
                let nd_idx = ui * np as usize + p as usize;
                if nd_idx < ctx.nodes.qual_mass_in.len() {
                    ctx.nodes.qual_mass_in[nd_idx] += q * c_gw;
                }
            }
            for p in 0..np {
                let c_gw = ctx.pollutants.c_gw[p as usize];
                if c_gw <= 0.0 {
                    continue;
                }
                let pi = p as usize;
                if pi < ctx.mass_balance.qual_routing_gw_in.len() {
                    ctx.mass_balance.qual_routing_gw_in[pi] += q * c_gw * dt;
                }
            }
        }
    }

    /// addIfaceLoads: routing interface file water and loads.
    fn add_iface_loads(&mut self, ctx: &mut SimulationContext, dt: f64) {
        let np = self.n_pollutants_;
        if !loaders_needed(np, ctx) {
            return;
        }
        if ctx.nodes.iface_qual_mass.is_empty() {
            return;
        }
        for i in 0..ctx.n_nodes() {
            let ui = i as usize;
            let q = ctx.nodes.iface_inflow[ui];
            if q < 0.0 && ui < self.qual_vol_lat_neg_.len() {
                self.qual_vol_lat_neg_[ui] += q * dt;
            }
            if q <= 0.0 {
                continue;
            }
            ctx.nodes.qual_vol_in[ui] += q * dt;
            self.qual_vol_lat_[ui] += q * dt;
            add_age_volume(ctx, i, q, WaterAgeSource::Iface);
            add_temp_volume(ctx, i, q, HeatSource::Iface);

            for p in 0..np {
                let nd_idx = ui * np as usize + p as usize;
                let mass_rate = if nd_idx < ctx.nodes.iface_qual_mass.len() { ctx.nodes.iface_qual_mass[nd_idx] } else { 0.0 };
                if mass_rate <= 0.0 {
                    continue;
                }
                if nd_idx < ctx.nodes.qual_mass_in.len() {
                    ctx.nodes.qual_mass_in[nd_idx] += mass_rate;
                }
            }
            for p in 0..np {
                let nd_idx = ui * np as usize + p as usize;
                let mass_rate = if nd_idx < ctx.nodes.iface_qual_mass.len() { ctx.nodes.iface_qual_mass[nd_idx] } else { 0.0 };
                if mass_rate <= 0.0 {
                    continue;
                }
                let pi = p as usize;
                if pi < ctx.mass_balance.qual_routing_ex_in.len() {
                    ctx.mass_balance.qual_routing_ex_in[pi] += mass_rate * dt;
                }
            }
        }
    }

    /// accumulateLinkLoads: each link's mass flow into its downstream node.
    fn accumulate_link_loads(&mut self, ctx: &mut SimulationContext, dt: f64) {
        let np = self.n_pollutants_;
        for j in 0..ctx.n_links() {
            let uj = j as usize;
            let q = ctx.links.flow[uj].abs();
            if q <= 0.0 {
                continue;
            }
            let downstream = if ctx.links.flow[uj] >= 0.0 { ctx.links.node2[uj] } else { ctx.links.node1[uj] };
            if downstream < 0 || downstream >= ctx.n_nodes() {
                continue;
            }
            let ud = downstream as usize;

            ctx.nodes.qual_vol_in[ud] += q * dt;

            for p in 0..np {
                let lp = uj * np as usize + p as usize;
                let np_idx = ud * np as usize + p as usize;

                let mut c_link = if lp < ctx.links.conc_old.len() { ctx.links.conc_old[lp] } else { 0.0 };
                if link_takes_upstream_value(ctx, j) {
                    let source = if ctx.links.flow[uj] >= 0.0 { ctx.links.node1[uj] } else { ctx.links.node2[uj] };
                    if lidnode::active(ctx, source) || lidnode::active(ctx, downstream) {
                        c_link = ctx.nodes.conc[(source * np + p) as usize];
                        if ctx.nodes.type_[source as usize] == NodeType::Outfall {
                            c_link = if ctx.options.outfall_backflow_zero { 0.0 } else { ctx.nodes.conc_old[(source * np + p) as usize] };
                        }
                    }
                    c_link = lidq::outlet_quality(ctx, source, j, p, c_link);
                }
                let mass = q * c_link;
                if link_takes_upstream_value(ctx, j) && lidnode::active(ctx, downstream) {
                    let source = if ctx.links.flow[uj] >= 0.0 { ctx.links.node1[uj] } else { ctx.links.node2[uj] };
                    if ctx.nodes.type_[source as usize] == NodeType::Outfall {
                        ctx.mass_balance.qual_routing_ex_in[p as usize] += mass * dt;
                    }
                }
                if np_idx < ctx.nodes.qual_mass_in.len() {
                    ctx.nodes.qual_mass_in[np_idx] += mass;
                    lidq::receive_quality(ctx, downstream, j, q * dt, p, mass, dt);
                }
            }
            let r = ctx.node_subtypes.storage_row(downstream);
            if r >= 0 {
                for port in ctx.node_subtypes.storages.lid_state[r as usize].quality_ports.iter() {
                    if port.link == j && port.volume > 0.0 {
                        ctx.nodes.qual_vol_in[ud] -= port.volume;
                    }
                }
            }
        }
        let _ = dt;
    }

    /// mixAtNodes: legacy findStorageQual (a mixed reactor) / findNodeQual (flow-through).
    fn mix_at_nodes(&mut self, ctx: &mut SimulationContext, dt: f64) {
        let np = self.n_pollutants_;
        let zero_bf_opt = ctx.options.outfall_backflow_zero;
        let reactions_active = legacy_reactions_active(ctx);
        let any_treatment = ctx.treatment.has_any();

        for i in 0..ctx.n_nodes() {
            let ui = i as usize;
            let v_old = ctx.nodes.old_volume[ui];
            let v_in = ctx.nodes.qual_vol_in[ui];
            let is_storage = ctx.nodes.type_[ui] == NodeType::Storage;
            let zero_bf = zero_bf_opt && ctx.nodes.type_[ui] == NodeType::Outfall;

            if !node_is_reactor(ctx, i) {
                // ---- findNodeQual
                let q_node = if dt > 0.0 { v_in / dt } else { 0.0 };
                for p in 0..np {
                    let idx = ui * np as usize + p as usize;
                    if idx >= ctx.nodes.conc.len() {
                        continue;
                    }
                    if q_node > LEGACY_ZERO {
                        let mut mass_in = if idx < ctx.nodes.qual_mass_in.len() { ctx.nodes.qual_mass_in[idx] * dt } else { 0.0 };
                        if mass_in < 0.0 {
                            let avail = ctx.nodes.conc_old[idx] * v_old;
                            if mass_in < -avail {
                                book_negative_source_clamp(ctx, i, p, -(avail + mass_in));
                                mass_in = -avail;
                            }
                        }
                        ctx.nodes.conc[idx] = mass_in / v_in;
                    } else if zero_bf {
                        ctx.nodes.conc[idx] = 0.0;
                    } else if ctx.nodes.depth[ui] > ZERO_DEPTH {
                        ctx.nodes.conc[idx] = ctx.nodes.conc_old[idx];
                    } else {
                        ctx.nodes.conc[idx] = if idx < ctx.nodes.qual_mass_in.len() { ctx.nodes.qual_mass_in[idx] } else { 0.0 };
                    }
                }
                continue;
            }

            // ---- findStorageQual: a mixed reactor.
            let f_evap = storage_evap_factor(ctx, i);
            let mut q_exfil = 0.0;
            if is_storage {
                let srow = ctx.node_subtypes.storage_row(i);
                if srow >= 0 {
                    let us = srow as usize;
                    let st = &ctx.node_subtypes.storages;
                    let v_exfil = if us < st.exfil_loss.len() { st.exfil_loss[us] } else { 0.0 };
                    q_exfil = if dt > 0.0 { v_exfil / dt } else { 0.0 };
                }
            }

            let q_in = if dt > 0.0 { v_in / dt } else { 0.0 };
            let dry = node_is_dry(ctx, i, q_in) && !lidnode::active(ctx, i);

            for p in 0..np {
                let idx = ui * np as usize + p as usize;
                if idx >= ctx.nodes.conc.len() {
                    continue;
                }
                let up = p as usize;

                let mut c1 = ctx.nodes.conc_old[idx];

                let seeped = q_exfil * c1 * dt;

                c1 *= f_evap;

                let mut k = if up < ctx.pollutants.k_decay.len() { ctx.pollutants.k_decay[up] } else { 0.0 };
                if reactions_active {
                    k = 0.0;
                }
                if k != 0.0 && any_treatment && node_has_treatment_for(ctx, i, p) {
                    k = 0.0;
                }
                let mut reacted = 0.0;
                if k != 0.0 {
                    let c2 = cmax(0.0, c1 * (1.0 - k * dt));
                    reacted = (c1 - c2) * v_old;
                    c1 = c2;
                }

                let mut mass_in = if idx < ctx.nodes.qual_mass_in.len() { ctx.nodes.qual_mass_in[idx] * dt } else { 0.0 };
                if mass_in < 0.0 {
                    let avail = ctx.nodes.conc_old[idx] * v_old;
                    if mass_in < -avail {
                        book_negative_source_clamp(ctx, i, p, -(avail + mass_in));
                        mass_in = -avail;
                    }
                }

                let mut c_new;
                if q_in <= LEGACY_ZERO {
                    c_new = c1;
                } else {
                    let c_in = mass_in / v_in;
                    let c_max = cmax(c1, c_in);
                    c_new = (c1 * v_old + mass_in) / (v_old + v_in);
                    c_new = cmin(c_new, c_max);
                    c_new = cmax(c_new, 0.0);
                }

                if dry {
                    if up < ctx.mass_balance.qual_routing_final_dry.len() {
                        ctx.mass_balance.qual_routing_final_dry[up] += c_new * ctx.nodes.volume[ui];
                    }
                    c_new = 0.0;
                }

                c_new = cmax(c_new, 0.0);
                ctx.nodes.conc[idx] = c_new;

                if reacted > 0.0 && up < ctx.mass_balance.qual_routing_reacted.len() {
                    ctx.mass_balance.qual_routing_reacted[up] += reacted;
                }
                if seeped != 0.0 && up < ctx.mass_balance.qual_routing_seep.len() {
                    ctx.mass_balance.qual_routing_seep[up] += seeped;
                }
            }
        }
    }

    /// applyDecay: a no-op stage (the first-order reaction lives inside the mixing kernels).
    fn apply_decay(&mut self, ctx: &mut SimulationContext, dt: f64) {
        let _ = (ctx, dt);
    }

    /// updateLinkQuality: legacy findLinkQual (DW/KW) and findSFLinkQual (STEADY).
    pub fn update_link_quality(&mut self, ctx: &mut SimulationContext, dt: f64) {
        let np = self.n_pollutants_;
        let is_steady = ctx.options.routing_model == RoutingModel::Steady;
        let reactions_active = legacy_reactions_active(ctx);

        for j in 0..ctx.n_links() {
            let uj = j as usize;

            let upstream = if ctx.links.flow[uj] >= 0.0 { ctx.links.node1[uj] } else { ctx.links.node2[uj] };
            if upstream < 0 || upstream >= ctx.n_nodes() {
                continue;
            }
            let un = upstream as usize;

            if link_takes_upstream_value(ctx, j) {
                for p in 0..np {
                    let li = uj * np as usize + p as usize;
                    let ni = un * np as usize + p as usize;
                    if li >= ctx.links.conc.len() || ni >= ctx.nodes.conc.len() {
                        break;
                    }
                    ctx.links.conc[li] = lidq::outlet_quality(ctx, upstream, j, p, ctx.nodes.conc[ni]);
                }
                continue;
            }

            let cr = ctx.link_subtypes.conduit_row(j);
            let ucr = (if cr < 0 { 0 } else { cr }) as usize;
            let cd = &ctx.link_subtypes.conduits;
            let barrels = if cr >= 0 { cd.barrels[ucr] as f64 } else { 1.0 };
            let q_seep = if cr >= 0 { cd.seep_loss_rate[ucr] * barrels } else { 0.0 };
            let v_evap = if cr >= 0 { cd.evap_loss_rate[ucr] * barrels * dt } else { 0.0 };

            let v_old = ctx.links.old_volume[uj];

            let mut f_evap = 1.0;
            if v_evap > 0.0 && v_old > ZERO_VOLUME {
                f_evap += v_evap / v_old;
            }

            let q_in = conduit_mixing_inflow(ctx, j, dt);

            let dry = !is_steady && link_is_dry(ctx, j);

            for p in 0..np {
                let li = uj * np as usize + p as usize;
                let ni = un * np as usize + p as usize;
                if li >= ctx.links.conc.len() || ni >= ctx.nodes.conc.len() {
                    continue;
                }
                let up = p as usize;

                let mut k = if up < ctx.pollutants.k_decay.len() { ctx.pollutants.k_decay[up] } else { 0.0 };
                if reactions_active {
                    k = 0.0;
                }

                let mut c_new;
                let mut reacted = 0.0;
                let seeped;

                if is_steady {
                    let n1 = ctx.links.node1[uj];
                    if n1 < 0 || n1 >= ctx.n_nodes() {
                        continue;
                    }
                    let n1i = n1 as usize * np as usize + up;
                    if n1i >= ctx.nodes.conc.len() {
                        continue;
                    }

                    let mut c1 = ctx.nodes.conc[n1i];
                    seeped = q_seep * c1 * dt;
                    c1 *= f_evap;
                    c_new = c1;
                    if k > 0.0 {
                        c_new = cmax(0.0, c1 * (-k * dt).exp()); // LIBM
                        reacted = (c1 - c_new) * ctx.links.flow[uj] * dt;
                    }
                } else {
                    let mut c1 = ctx.links.conc_old[li];
                    seeped = q_seep * c1 * dt;
                    c1 *= f_evap;

                    let c2 = c1 * cmax(1.0 - k * dt, 0.0);
                    reacted = (c1 - c2) * v_old;

                    if q_in <= LEGACY_ZERO {
                        c_new = c2;
                    } else {
                        let w_in = lidq::outlet_quality(ctx, upstream, j, p, ctx.nodes.conc[ni]) * q_in;
                        let v_in = q_in * dt;
                        let c_in = w_in * dt / v_in;
                        let c_max = cmax(c2, c_in);
                        c_new = (c2 * v_old + w_in * dt) / (v_old + v_in);
                        c_new = cmin(c_new, c_max);
                        c_new = cmax(c_new, 0.0);
                    }

                    if dry {
                        if up < ctx.mass_balance.qual_routing_final_dry.len() {
                            ctx.mass_balance.qual_routing_final_dry[up] += c_new * ctx.links.volume[uj];
                        }
                        c_new = 0.0;
                    }
                }

                c_new = cmax(c_new, 0.0);
                ctx.links.conc[li] = c_new;

                if reacted > 0.0 && up < ctx.mass_balance.qual_routing_reacted.len() {
                    ctx.mass_balance.qual_routing_reacted[up] += reacted;
                }
                if seeped != 0.0 && up < ctx.mass_balance.qual_routing_seep.len() {
                    ctx.mass_balance.qual_routing_seep[up] += seeped;
                }
            }
        }
    }

    /// applyTreatment: [TREATMENT] at every node that has one (legacy treatmnt_treat).
    /// full_inflow: the caller assembled the node's complete inflow (the LEGACY pass).
    pub fn apply_treatment(&mut self, ctx: &mut SimulationContext, dt: f64, full_inflow: bool) {
        let np = self.n_pollutants_;
        if np <= 0 {
            return;
        }
        if !ctx.treatment.has_any() {
            return;
        }
        let nn = ctx.n_nodes();
        for j in 0..nn {
            let uj = j as usize;
            if !ctx.treatment.has_treatment[uj] {
                continue;
            }
            let vol_in = ctx.nodes.qual_vol_in[uj];
            let q_raw = if dt > 0.0 { vol_in / dt } else { 0.0 };
            for p in 0..np {
                let mi = uj * np as usize + p as usize;
                let have = mi < ctx.nodes.qual_mass_in.len();
                if full_inflow {
                    ctx.treatment.cin[p as usize] = if q_raw > 0.0 && have { ctx.nodes.qual_mass_in[mi] / q_raw } else { 0.0 };
                } else {
                    ctx.treatment.cin[p as usize] = if vol_in > 0.0 && have { ctx.nodes.qual_mass_in[mi] / vol_in } else { 0.0 };
                }
            }
            apply_node_treatment(ctx, j, dt, q_raw, None);
        }
    }
}

/// addQualityOnlyInflows' `apply` lambda for one receipt container. `ex` selects the ledger
/// row: qual_routing_ex_in (surface outfall) or qual_routing_gw_in (gaining conduits).
fn apply_quality_only(ctx: &mut SimulationContext, q: &mut crate::core::context::QualityOnlyInflow, dt: f64, ex: bool) {
    if !q.any {
        return;
    }
    let np = ctx.n_pollutants();
    let ns = ctx.reactions.n_species();
    let nn = cmin(q.vol.len(), ctx.n_nodes() as usize);
    let want = ctx.n_nodes() as usize * if ns > 0 { ns as usize } else { 0 };
    if ns > 0 && ctx.reactions.msx_ext_mass_in.len() != ctx.n_nodes() as usize * ns as usize {
        ctx.reactions.msx_ext_mass_in.clear();
        ctx.reactions.msx_ext_mass_in.resize(want, 0.0);
    }
    for i in 0..nn {
        if q.vol[i] <= 0.0 {
            continue;
        }
        ctx.nodes.qual_vol_in[i] += q.vol[i];
        for p in 0..np {
            let k = i * np as usize + p as usize;
            if k >= q.mass.len() || k >= ctx.nodes.qual_mass_in.len() {
                continue;
            }
            ctx.nodes.qual_mass_in[k] += q.mass[k] / dt;
            let row = if ex { &mut ctx.mass_balance.qual_routing_ex_in } else { &mut ctx.mass_balance.qual_routing_gw_in };
            if (p as usize) < row.len() {
                row[p as usize] += q.mass[k];
            }
        }
        for m in 0..ns {
            let k = i * ns as usize + m as usize;
            if k < q.msx.len() && k < ctx.reactions.msx_ext_mass_in.len() {
                ctx.reactions.msx_ext_mass_in[k] += q.msx[k] / dt;
            }
        }
        if ctx.options.water_age && i < ctx.water_age_state.node_age_vol_in.len() {
            ctx.water_age_state.node_age_vol_in[i] += q.age_vol[i] / dt;
        }
        if ctx.options.heat_transport && i < ctx.heat_state.node_temp_vol_in.len() {
            ctx.heat_state.node_temp_vol_in[i] += q.temp_vol[i] / dt;
        }
    }
    q.clear();
}

/// applyTreatment's getRemoval lambda (legacy treatmnt.c getRemoval with the node-wide cycle
/// latch).
#[allow(clippy::too_many_arguments)]
fn get_removal(
    p: i32,
    uj: usize,
    np: i32,
    compiled: &[crate::data::treatment_types::TreatExpr],
    removal: &mut Vec<f64>,
    cin: &[f64],
    conc: &[f64],
    cpollut: &[f64],
    cyclic_err: &mut bool,
    args: (f64, f64, f64, f64, f64, f64),
) -> f64 {
    let (dt, hrt_hours, q, v, d, area) = args;
    let up = p as usize;
    if removal[up] > 1.0 || *cyclic_err {
        *cyclic_err = true;
        return 0.0;
    }
    if removal[up] >= 0.0 && removal[up] <= 1.0 {
        return removal[up];
    }

    let idx = uj * np as usize + up;
    if idx >= compiled.len() || compiled[idx].tokens.is_empty() {
        removal[up] = 0.0;
        return 0.0;
    }

    removal[up] = 10.0;

    let expr = &compiled[idx];
    let ci_idx = uj * np as usize + up;
    let c_node = if ci_idx < conc.len() { conc[ci_idx] } else { 0.0 };
    let c_in = if expr.is_removal { cin[up] } else { c_node };

    if c_node <= 0.0 {
        removal[up] = 0.0;
        return 0.0;
    }

    for tok in expr.tokens.iter() {
        if tok.var == TreatVar::RPollut && tok.pollut_ref >= 0 && tok.pollut_ref < np {
            get_removal(tok.pollut_ref, uj, np, compiled, removal, cin, conc, cpollut, cyclic_err, args);
        }
    }

    let mut r_view = removal.clone();
    for x in r_view.iter_mut() {
        if *x < 0.0 || *x > 1.0 {
            *x = 0.0;
        }
    }

    let mut result = treatment::evaluate(expr, c_in, dt, hrt_hours, q, v, d, Some(cin), Some(&r_view), np, area, Some(cpollut));
    result = cmax(result, 0.0);

    if expr.is_removal {
        removal[up] = cmin(result, 1.0);
    } else {
        result = cmin(result, c_node);
        removal[up] = if c_node > 0.0 { 1.0 - result / c_node } else { 0.0 };
    }
    removal[up]
}

/// quality::applyNodeTreatment: one node's [TREATMENT] against the caller's inflow figures.
/// `cin` None: the caller's concentrations are already in treatment.cin (the LEGACY pass).
pub fn apply_node_treatment(ctx: &mut SimulationContext, j: i32, dt: f64, q_raw: f64, cin: Option<&[f64]>) {
    let np = ctx.n_pollutants();
    if np <= 0 {
        return;
    }
    if !ctx.treatment.has_any() {
        return;
    }
    let uj = j as usize;
    if uj >= ctx.treatment.has_treatment.len() || !ctx.treatment.has_treatment[uj] {
        return;
    }
    if let Some(c) = cin {
        for p in 0..np as usize {
            ctx.treatment.cin[p] = c[p];
        }
    }

    let unit_sys = ctx.options.flow_units as i32;
    let ucf_flow = ucf::ucf(ucf::FLOW, &ctx.options);
    let ucf_length = ucf::ucf(ucf::LENGTH, &ctx.options);

    let mut hrt_hours = ctx.nodes.hrt[uj] / 3600.0;
    let q = q_raw * ucf_flow;
    let v = (ctx.nodes.old_volume[uj] + ctx.nodes.volume[uj]) * 0.5;
    let d_ft = (ctx.nodes.depth[uj] + ctx.nodes.old_depth[uj]) * 0.5;
    let d = d_ft * ucf_length;

    let us = ucf::get_unit_system(unit_sys);
    let a1 = node::get_surf_area(&ctx.nodes, j, ctx.nodes.old_depth[uj], Some(&ctx.tables), us, Some(&ctx.node_subtypes));
    let a2 = node::get_surf_area(&ctx.nodes, j, ctx.nodes.depth[uj], Some(&ctx.tables), us, Some(&ctx.node_subtypes));
    let area = (a1 + a2) * 0.5 * ucf_length * ucf_length;

    if ctx.nodes.type_[uj] == NodeType::Storage {
        let v_old = ctx.nodes.old_volume[uj];
        if v_old < LEGACY_ZERO {
            ctx.nodes.hrt[uj] = 0.0;
        } else {
            let qdt = q_raw.abs() * dt;
            ctx.nodes.hrt[uj] = cmax((ctx.nodes.hrt[uj] + dt) * v_old / (v_old + qdt), 0.0);
        }
        hrt_hours = ctx.nodes.hrt[uj] / 3600.0;
    }

    let unp = np as usize;
    for p in 0..unp {
        ctx.treatment.removal[p] = -1.0;
    }

    let mut cpollut = vec![0.0f64; unp];
    for p in 0..unp {
        let ci = uj * unp + p;
        let p_is_removal = ci < ctx.treatment.compiled.len() && ctx.treatment.compiled[ci].is_removal;
        cpollut[p] = if p_is_removal { ctx.treatment.cin[p] } else if ci < ctx.nodes.conc.len() { ctx.nodes.conc[ci] } else { 0.0 };
    }

    let mut cyclic_err = false;
    let args = (dt, hrt_hours, q, v, d, area);
    {
        let tr = &mut ctx.treatment;
        let conc = &ctx.nodes.conc;
        for p in 0..unp {
            let ci2 = uj * unp + p;
            let no_eqn = ci2 >= tr.compiled.len() || tr.compiled[ci2].tokens.is_empty();
            if no_eqn {
                tr.removal[p] = 0.0;
            } else if tr.compiled[ci2].is_removal && q_raw.abs() <= LEGACY_ZERO {
                tr.removal[p] = 0.0;
            } else {
                get_removal(p as i32, uj, np, &tr.compiled, &mut tr.removal, &tr.cin, conc, &cpollut, &mut cyclic_err, args);
            }
        }
    }

    let v_old = ctx.nodes.old_volume[uj];
    for p in 0..unp {
        let rr = ctx.treatment.removal[p];
        if rr <= 0.0 {
            continue;
        }
        let ci = uj * unp + p;
        if ci >= ctx.nodes.conc.len() {
            continue;
        }
        let is_removal = ctx.treatment.compiled[ci].is_removal;
        let c_node = ctx.nodes.conc[ci];
        let mut c_out;
        if is_removal {
            let c_in = ctx.treatment.cin[p];
            c_out = if c_in > 0.0 { (1.0 - rr) * c_in } else { c_node };
            c_out = cmin(c_out, c_node);
        } else {
            c_out = (1.0 - rr) * c_node;
        }
        c_out = cmax(c_out, 0.0);

        let c_in_p = ctx.treatment.cin[p];
        let c_old_p = if ci < ctx.nodes.conc_old.len() { ctx.nodes.conc_old[ci] } else { c_node };
        let mut mass_lost = 0.0;
        if dt > 0.0 {
            mass_lost = (c_in_p * q_raw * dt + c_old_p * v_old - c_out * (q_raw * dt + v_old)) / dt;
            mass_lost = cmax(0.0, mass_lost);
        }
        if mass_lost > 0.0 && p < ctx.mass_balance.qual_routing_reacted.len() {
            ctx.mass_balance.qual_routing_reacted[p] += mass_lost * dt;
        }
        ctx.nodes.conc[ci] = c_out;
    }
}
