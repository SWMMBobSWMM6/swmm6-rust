// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/SWMMEngine.cpp (SWMMEngine::initQuality,
// refreshLanduseParams, initGeometry, initMassBalance, both reportedNodeVolume overloads,
// ordinarySourceOwners, ownsSpatialSource and the file-local snowCoverVolumeFt3)
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Three legs of SWMMEngine::init_modules (1D build: OPENSWMM_HAS_2D off). Like the other
// init legs (core/engine_init_hydraulics.rs, core/engine_init_hydrology.rs) it is a sibling
// module of core::engine, whose members it reaches as pub(crate).
//
// Differences of form:
//   * initQuality is noexcept upstream, and treatment::parse can throw (std::stod on a
//     [TREATMENT] number such as `.` or `1e999`), which calls std::terminate (U57). Here
//     init_quality returns Err(Terminated) at that point, as init_hydraulics does for U65:
//     the translation stops where upstream stops, and never aborts the process.
//   * initQuality's cycle check is a recursive std::function lambda; here `dfs` is a nested
//     fn over the same colour and cycle state, with the same visiting order.
//   * initMassBalance reads `ctx_.state_accessors.get_fv_retained_pollutant_mass`, which
//     open() always installs (non-metadata path) as `[this]{ return
//     quality_.fvRetainedPollutantMass(ctx_); }`; initialize(), the only caller, refuses
//     metadata handles. ctx_.state_accessors is not wired in the translation yet (see
//     engine.rs), so the closure's body is called directly.
//   * ordinarySourceOwners returns nullptr without OPENSWMM_HAS_2D, so ownsSpatialSource is
//     always false here; both are kept for the call sites.
//   * runoff_, lid_, snow_ and groundwater_ are the hydrology solvers (initHydrology fills
//     them); initQuality and initMassBalance only read them.

use crate::core::engine::{SwmmEngine, SWMM_ERR_PARSE};
use crate::core::engine_init_hydraulics::Terminated;
use crate::math::math_expr::StodError;
use crate::core::constants::TINY;
use crate::core::context::SimulationContext;
use crate::core::error_codes::{format_error, ERR_CYCLIC_TREATMENT};
use crate::core::options::RoutingModel;
use crate::core::unit_conversion::{self as ucf, LANDAREA};
use crate::data::link::{LinkType, XsectShape};
use crate::data::node::NodeType;
use crate::data::treatment_types::{TokenType, TreatExpr, TreatVar};
use crate::hydrology::lid_node::treatment;
use crate::hydrology::lid_node_held as lidnode;
use crate::hydrology::snow::{SnowSoA, N_SUBAREAS};
use crate::quality::landuse::{BuildupType, WashoffType};

/// snowCoverVolumeFt3 (file-local): snow water equivalent plus free water held in every
/// enabled pack (ft3), over the snow-covered (non-LID) project area.
pub(crate) fn snow_cover_volume_ft3(ctx: &SimulationContext, soa: &SnowSoA) -> f64 {
    let mut total = 0.0;
    for i in 0..ctx.n_subcatches() {
        let ui = i as usize;
        if ctx.subcatches.snowpack[ui] < 0 {
            continue;
        }
        if soa.enabled[ui] == 0 {
            continue;
        }
        let mut depth = 0.0;
        for k in 0..N_SUBAREAS {
            let idx = (i * N_SUBAREAS + k) as usize;
            if idx >= soa.wsnow.len() {
                break;
            }
            depth += (soa.wsnow[idx] + soa.fw[idx]) * soa.fArea[idx];
        }
        let lid_ft2 = if ui < ctx.subcatches.total_lid_area_ft2.len() { ctx.subcatches.total_lid_area_ft2[ui] } else { 0.0 };
        let mut area_ft2 = ctx.subcatches.area[ui] / ucf::ucf(LANDAREA, &ctx.options) - lid_ft2;
        if area_ft2 < 0.0 {
            area_ft2 = 0.0;
        }
        total += depth * area_ft2;
    }
    total
}

/// initQuality's DFS (a recursive std::function upstream).
fn dfs(p: i32, dep: &[Vec<i32>], color: &mut [i32], cycle: &mut bool) {
    if *cycle {
        return;
    }
    color[p as usize] = 1;
    for &q in &dep[p as usize] {
        if color[q as usize] == 1 {
            *cycle = true;
            return;
        }
        if color[q as usize] == 0 {
            dfs(q, dep, color, cycle);
        }
    }
    color[p as usize] = 2;
}

impl SwmmEngine {
    // ========================================================================
    // initQuality()
    // ========================================================================

    pub(crate) fn init_quality(&mut self) -> Result<(), Terminated> {
        // 7. Quality solver
        self.quality_.init(self.ctx_.n_nodes(), self.ctx_.n_links(), self.ctx_.n_pollutants());

        // Seed initial concentrations from [POLLUTANTS] Cinit (legacy qualrout_init): a wet
        // node/link starts at the pollutant's initial concentration, a dry one at zero.
        {
            let np = self.ctx_.n_pollutants();
            if np > 0 && (self.ctx_.files.hotstart_use_path.is_empty() || self.hotstart_deferred_) {
                // Legacy qualrout.c: static const double ZeroDepth = 0.003281 (1 mm).
                const ZERO_DEPTH: f64 = 0.003281;
                let ctx = &mut self.ctx_;
                for i in 0..ctx.n_nodes() {
                    let ui = i as usize;
                    let wet = ctx.nodes.depth[ui] > ZERO_DEPTH;
                    for p in 0..np {
                        let idx = ui * (np as usize) + p as usize;
                        if idx >= ctx.nodes.conc.len() {
                            continue;
                        }
                        let c = if wet { ctx.pollutants.init_conc[p as usize] } else { 0.0 };
                        ctx.nodes.conc[idx] = c;
                        ctx.nodes.conc_old[idx] = c;
                    }
                }
                for j in 0..ctx.n_links() {
                    let uj = j as usize;
                    let wet = ctx.links.depth[uj] > ZERO_DEPTH;
                    for p in 0..np {
                        let idx = uj * (np as usize) + p as usize;
                        if idx >= ctx.links.conc.len() {
                            continue;
                        }
                        let c = if wet { ctx.pollutants.init_conc[p as usize] } else { 0.0 };
                        ctx.links.conc[idx] = c;
                        ctx.links.conc_old[idx] = c;
                    }
                }

                // [INITIAL_QUALITY] per-element pollutant overrides, after the global Cinit
                // seed, wet or dry. Reserved-species rows (kind < 0) are consumed elsewhere.
                {
                    let iq = &ctx.initial_quality;
                    for r in 0..iq.count() {
                        let ur = r as usize;
                        let k = iq.kind[ur];
                        if k < 0 || k >= np {
                            continue; // reserved/unresolved
                        }
                        let ei = iq.elem_idx[ur];
                        if ei < 0 {
                            continue;
                        }
                        let (conc, conc_old) = if iq.is_link[ur] != 0 {
                            (&mut ctx.links.conc, &mut ctx.links.conc_old)
                        } else {
                            (&mut ctx.nodes.conc, &mut ctx.nodes.conc_old)
                        };
                        let idx = (ei as usize) * (np as usize) + k as usize;
                        if idx >= conc.len() {
                            continue;
                        }
                        conc[idx] = iq.value[ur];
                        conc_old[idx] = iq.value[ur];
                    }
                }
            }
        }

        // 10. Treatment: resize for nodes x pollutants + compile expressions
        if self.ctx_.n_pollutants() > 0 && self.ctx_.n_nodes() > 0 {
            if self.ctx_.treatment.n_nodes == 0 {
                let (nn, np) = (self.ctx_.n_nodes(), self.ctx_.n_pollutants());
                self.ctx_.treatment.resize(nn, np);
            }
            // Compile treatment expressions
            let np = self.ctx_.n_pollutants();
            for n in 0..self.ctx_.n_nodes() {
                for p in 0..np {
                    let idx = (n * np + p) as usize;
                    let ctx = &mut self.ctx_;
                    if ctx.treatment.expressions[idx].is_empty() {
                        continue;
                    }
                    let mut te = TreatExpr::default();
                    // Err: std::stod threw inside a noexcept function upstream (terminate);
                    // see the file notes.
                    let r = treatment::parse_names(&ctx.treatment.expressions[idx], &mut te, ctx.pollutant_names.names());
                    if let Err(e) = r {
                        return Err(Terminated {
                            what: "stod",
                            exception: match e {
                                StodError::InvalidArgument => "std::invalid_argument",
                                StodError::OutOfRange => "std::out_of_range",
                            },
                        });
                    }
                    if r == Ok(0) {
                        te.pollutant_idx = p;
                        ctx.treatment.compiled[idx] = te;
                        ctx.treatment.has_treatment[n as usize] = true;
                    }
                }
            }
        }

        // Gap #85: cyclic treatment expression detection.
        if self.ctx_.n_pollutants() > 0 && self.ctx_.n_nodes() > 0 {
            let np = self.ctx_.n_pollutants();
            let nn = self.ctx_.n_nodes();
            for n in 0..nn {
                let un = n as usize;
                if !self.ctx_.treatment.has_treatment[un] {
                    continue;
                }
                // dep[p] = the pollutants p's expression depends on
                let mut dep: Vec<Vec<i32>> = vec![Vec::new(); np as usize];
                for p in 0..np {
                    let idx = (n * np + p) as usize;
                    let te = &self.ctx_.treatment.compiled[idx];
                    for tok in &te.tokens {
                        if tok.type_ == TokenType::Variable
                            && tok.var == TreatVar::RPollut
                            && tok.pollut_ref >= 0
                            && tok.pollut_ref < np
                            && tok.pollut_ref != p
                        {
                            dep[p as usize].push(tok.pollut_ref);
                        }
                    }
                }
                // DFS cycle detection (0=white, 1=gray, 2=black)
                let mut color: Vec<i32> = vec![0; np as usize];
                let mut cycle = false;
                let mut p = 0;
                while p < np && !cycle {
                    if color[p as usize] == 0 {
                        dfs(p, &dep, &mut color, &mut cycle);
                    }
                    p += 1;
                }
                if cycle {
                    let e = format_error(ERR_CYCLIC_TREATMENT, &self.ctx_.node_names.names()[un]);
                    self.ctx_.errors.push(e);
                    let last = self.ctx_.errors.last().unwrap().clone();
                    self.set_error(SWMM_ERR_PARSE, &last);
                }
            }
        }

        // 11. Buildup/washoff: resize for landuses x pollutants
        if self.ctx_.n_landuses() > 0 && self.ctx_.n_pollutants() > 0 {
            if self.ctx_.buildup.n_landuses == 0 {
                let (nlu, np) = (self.ctx_.n_landuses(), self.ctx_.n_pollutants());
                self.ctx_.buildup.resize(nlu, np);
                self.ctx_.washoff.resize(nlu, np);
            }
        }

        // 11b. Initialize landuse quality solver and surface quality state
        if self.ctx_.n_pollutants() > 0 {
            let nlu = crate::core::cmath::cmax(self.ctx_.n_landuses(), 1);
            let np = self.ctx_.n_pollutants();
            self.landuse_solver_.init(nlu, np);
            self.surface_quality_.resize(self.ctx_.n_subcatches(), nlu, np);
            let nsc = self.ctx_.n_subcatches();
            self.ctx_.subcatches.resize_total_load(nsc, np);

            // Transfer parsed BuildupData/WashoffData into LanduseSolver params
            self.refresh_landuse_params();

            // Initialize quality mass balance vectors
            self.ctx_.mass_balance.resize_quality(np);

            // Initial surface buildup (legacy landuse_getInitBuildup): a [LOADINGS] initial
            // loading seeds initBuildup*fArea; otherwise the land use's buildup function
            // runs over the antecedent dry period. Stored per land use, per normalizer unit.
            let dry_days = self.ctx_.options.dry_days;
            if self.ctx_.n_landuses() > 0 {
                for i in 0..self.ctx_.n_subcatches() {
                    let ui = i as usize;
                    for lu in 0..nlu {
                        let ctx = &mut self.ctx_;
                        // legacy landuse_initState: lastSwept = StartDateTime - sweepDays0
                        {
                            let sw_idx = ui * (nlu as usize) + lu as usize;
                            if sw_idx < ctx.subcatches.sweep_last_swept.len() {
                                ctx.subcatches.sweep_last_swept[sw_idx] = ctx.landuses.last_swept[lu as usize];
                            }
                        }
                        let cov_idx = ui * (nlu as usize) + lu as usize;
                        let frac = if cov_idx < ctx.subcatches.coverage.len() { ctx.subcatches.coverage[cov_idx] / 100.0 } else { 0.0 };
                        if frac <= 0.0 {
                            continue;
                        }

                        // legacy fArea = fraction * area (landuse.c:377)
                        let f_area = frac * ctx.subcatches.area[ui];

                        for p in 0..np {
                            let k = (lu * np + p) as usize;
                            let bp = self.landuse_solver_.buildup_params[k];

                            let norm = if bp.normalizer == 0 { f_area } else { frac * ctx.subcatches.curb_length[ui] };

                            let li = ui * (np as usize) + p as usize;
                            let init_load = if li < ctx.subcatches.init_loading.len() { ctx.subcatches.init_loading[li] } else { 0.0 };

                            let mut mass = 0.0; // per normalizer unit
                            if init_load > 0.0 {
                                if norm > 0.0 {
                                    mass = init_load * f_area / norm;
                                }
                            } else if bp.type_ != BuildupType::None && dry_days > 0.0 {
                                let (c0, c1, c2) = (bp.coeff[0], bp.coeff[1], bp.coeff[2]);
                                match bp.type_ {
                                    BuildupType::Power => {
                                        mass = crate::core::cmath::cmin(c1 * dry_days.powf(c2), c0); // LIBM
                                    }
                                    BuildupType::Expon => {
                                        mass = c0 * (1.0 - (-c1 * dry_days).exp()); // LIBM
                                    }
                                    BuildupType::Satur => {
                                        mass = if c2 + dry_days > 0.0 { c0 * dry_days / (c2 + dry_days) } else { 0.0 };
                                    }
                                    _ => {}
                                }
                            } else {
                                continue;
                            }

                            let bu = self.surface_quality_.bu_idx(i, lu, p);
                            self.surface_quality_.buildup[bu] = mass;
                            ctx.mass_balance.qual_init_buildup[p as usize] += mass * norm;
                        }
                    }
                }
            }
        }

        // MSX element state is initial state (seeded here, not on the first step).
        if self.ctx_.reactions.n_species() > 0 {
            crate::transport::reaction_legacy_binding::ensure_msx_state(&mut self.ctx_);
        }
        // BW-MSX: initial buildup of the reactions component's species.
        if !self.ctx_.options.ignore_quality {
            crate::transport::msx_surface::init_state(&mut self.ctx_);
        }
        if !self.ctx_.options.ignore_quality {
            crate::quality::msx_lid::initialize(&mut self.ctx_);
        }
        if !self.ctx_.options.ignore_quality {
            crate::quality::msx_lid::initialize_units(&mut self.ctx_, &self.lid_);
        }
        Ok(())
    }

    /// refreshLanduseParams: re-derive the land-use buildup/washoff parameter cache from the
    /// live context (BuildupData/WashoffData). The buildup pool is left untouched.
    pub(crate) fn refresh_landuse_params(&mut self) {
        let np = self.ctx_.n_pollutants();
        let nlu = self.ctx_.n_landuses();
        if np <= 0 || nlu <= 0 {
            return;
        }
        if self.landuse_solver_.buildup_params.len() as i32 != nlu * np {
            return;
        }
        for lu in 0..nlu {
            for p in 0..np {
                let k = (lu * np + p) as usize;
                let ctx = &self.ctx_;
                let bp = &mut self.landuse_solver_.buildup_params[k];
                bp.type_ = BuildupType::from_i32(ctx.buildup.func_type[k]);
                bp.coeff[0] = ctx.buildup.coeff1[k]; // max buildup
                bp.coeff[1] = ctx.buildup.coeff2[k]; // rate constant
                bp.coeff[2] = ctx.buildup.coeff3[k]; // exponent/half-sat
                bp.normalizer = ctx.buildup.normalizer[k];
                bp.max_days = 0.0;
                // max_days: time to reach 99.9% of max buildup
                if bp.type_ == BuildupType::Expon && bp.coeff[1] > 0.0 {
                    bp.max_days = -(0.001f64).ln() / bp.coeff[1]; // LIBM
                } else if bp.type_ == BuildupType::Power && bp.coeff[1] > 0.0 && bp.coeff[2] > 0.0 && bp.coeff[0] > 0.0 {
                    bp.max_days = (bp.coeff[0] / bp.coeff[1]).powf(1.0 / bp.coeff[2]); // LIBM
                } else if bp.type_ == BuildupType::Satur && bp.coeff[2] > 0.0 && bp.coeff[0] > 0.0 {
                    bp.max_days = 999.0 * bp.coeff[2]; // asymptotic
                }

                let wp = &mut self.landuse_solver_.washoff_params[k];
                wp.type_ = WashoffType::from_i32(ctx.washoff.func_type[k]);
                wp.coeff = ctx.washoff.coeff[k];
                wp.expon = ctx.washoff.expon[k];
                wp.sweep_effic = ctx.washoff.sweep_effic[k];
                wp.bmp_effic = ctx.washoff.bmp_effic[k];
            }
        }
    }

    // ========================================================================
    // initGeometry()
    // ========================================================================

    pub(crate) fn init_geometry(&mut self) {
        let ctx = &mut self.ctx_;
        // 12. Node crown elevations (legacy dynwave_init): highest connecting link crown.
        for i in 0..ctx.n_nodes() {
            ctx.nodes.crown_elev[i as usize] = ctx.nodes.invert_elev[i as usize];
        }
        for j in 0..ctx.n_links() {
            let uj = j as usize;
            let n1 = ctx.links.node1[uj];
            let n2 = ctx.links.node2[uj];
            let un1 = n1 as usize;
            let un2 = n2 as usize;
            // PARITY dynwave.c:143-152: crowns accumulate over ALL link types.
            let mut off1 = ctx.links.offset1[uj];
            let off2 = ctx.links.offset2[uj];
            match ctx.links.type_[uj] {
                LinkType::Orifice => {} // offset1 raised, offset2 as parsed
                LinkType::Weir => {
                    let wr = ctx.link_subtypes.weir_row(j);
                    off1 = if wr >= 0 { ctx.link_subtypes.weirs.crest_height[wr as usize] } else { 0.0 };
                }
                LinkType::Outlet => {
                    let olr = ctx.link_subtypes.outlet_row(j);
                    off1 = if olr >= 0 { ctx.link_subtypes.outlets.crest_height[olr as usize] } else { 0.0 };
                }
                _ => {} // CONDUIT: parsed offsets; PUMP: zeros
            }
            // legacy xsect.yFull of a DUMMY section is TINY (1e-6), not 0.
            let mut yf = ctx.links.xsect_y_full[uj];
            if ctx.links.type_[uj] == LinkType::Outlet
                || (ctx.links.type_[uj] == LinkType::Conduit && ctx.links.xsect_shape[uj] == XsectShape::Dummy)
            {
                yf = TINY;
            }
            let z1 = ctx.nodes.invert_elev[un1] + off1 + yf;
            let z2 = ctx.nodes.invert_elev[un2] + off2 + yf;
            if z1 > ctx.nodes.crown_elev[un1] {
                ctx.nodes.crown_elev[un1] = z1;
            }
            if z2 > ctx.nodes.crown_elev[un2] {
                ctx.nodes.crown_elev[un2] = z2;
            }
        }

        // 12b. legacy Node.degree, in full (toposort_sortLinks + flowrout.c
        // validateGeneralLayout): the OUTFLOW-link count, negated under the dynamic wave for
        // a node with no INFLOW links.
        {
            let nn = ctx.n_nodes();
            ctx.nodes.degree.fill(0);
            let mut inflow_links: Vec<i32> = vec![0; nn as usize];
            for j in 0..ctx.n_links() {
                let uj = j as usize;
                let n1 = ctx.links.node1[uj];
                let n2 = ctx.links.node2[uj];
                if n1 < 0 || n1 >= nn || n2 < 0 || n2 >= nn {
                    continue;
                }
                // outflow count (toposort.c:75-91)
                let mut n = if ctx.links.direction[uj] < 0 { n2 } else { n1 };
                if ctx.nodes.type_[n as usize] == NodeType::Outfall {
                    n = if ctx.links.direction[uj] < 0 { n1 } else { n2 };
                }
                ctx.nodes.degree[n as usize] += 1;
                // inflow count (flowrout.c validateGeneralLayout)
                let mut dn = n1;
                if ctx.nodes.type_[dn as usize] != NodeType::Outfall {
                    dn = n2;
                }
                inflow_links[dn as usize] += 1;
            }
            if ctx.options.routing_model == RoutingModel::Dynwave || ctx.options.routing_model == RoutingModel::Fv {
                for i in 0..nn {
                    let ui = i as usize;
                    if inflow_links[ui] == 0 {
                        ctx.nodes.degree[ui] = -ctx.nodes.degree[ui];
                    }
                }
            }
        }

        // 13. Node full volumes
        let unit_sys = ucf::get_unit_system(ctx.options.flow_units as i32);
        for j in 0..ctx.n_nodes() {
            let uj = j as usize;
            let v = crate::hydraulics::node::get_volume(&ctx.nodes, j, ctx.nodes.full_depth[uj], Some(&ctx.tables), unit_sys, Some(&ctx.node_subtypes));
            ctx.nodes.full_volume[uj] = v;
        }
    }

    // ========================================================================
    // initMassBalance()
    // ========================================================================

    /// reportedNodeVolume(i): the legacy-convention reported volume at the node's state.
    pub(crate) fn reported_node_volume(&self, i: i32) -> f64 {
        let ui = i as usize;
        self.reported_node_volume_at(i, self.ctx_.nodes.depth[ui], self.ctx_.nodes.volume[ui])
    }

    /// reportedNodeVolume(i, depth, volume).
    pub(crate) fn reported_node_volume_at(&self, i: i32, depth: f64, volume: f64) -> f64 {
        let ctx = &self.ctx_;
        let ui = i as usize;
        if ctx.nodes.storage_area[ui] > 0.0 {
            return volume;
        }
        if ctx.nodes.type_[ui] == NodeType::Storage {
            return volume + lidnode::held_volume(ctx, i); // mobile + retained LID water
        }
        // Kinematic wave / steady flow: the accumulated net-inflow volume, as it is.
        if ctx.options.routing_model == RoutingModel::Kinwave || ctx.options.routing_model == RoutingModel::Steady {
            return volume;
        }
        let fd = ctx.nodes.full_depth[ui];
        if !(fd > 0.0) {
            return 0.0;
        }
        // A PONDING node holds real water above its rim; report the volume.
        let can_pond = ctx.options.allow_ponding || (ui < ctx.coupled_node.len() && ctx.coupled_node[ui] != 0);
        if can_pond && ctx.nodes.ponded_area[ui] > 0.0 && depth > fd {
            return volume;
        }
        // Legacy convention: a plain junction contributes ZERO (rpt_full_volume is 0).
        ctx.nodes.rpt_full_volume[ui] * (depth / fd)
    }

    /// ordinarySourceOwners: nullptr without OPENSWMM_HAS_2D.
    pub(crate) fn ordinary_source_owners(&self) -> Option<&Vec<i32>> {
        None
    }

    pub(crate) fn owns_spatial_source(&self, source: i32) -> bool {
        match self.ordinary_source_owners() {
            Some(owners) => source >= 0 && (source as usize) < owners.len() && owners[source as usize] >= 0,
            None => false,
        }
    }

    pub(crate) fn init_mass_balance(&mut self) {
        // 14. Mass balance: record initial storage (nodes + links, matching legacy)
        self.ctx_.mass_balance.reset();
        let nn = self.ctx_.n_nodes();
        let np = self.ctx_.n_pollutants();
        {
            let oe = &mut self.ctx_.outfall_exports;
            oe.water_m3.clear();
            oe.water_m3.resize(nn as usize, 0.0);
            oe.reverse_m3.clear();
            oe.reverse_m3.resize(nn as usize, 0.0);
            oe.heat_j.clear();
            oe.heat_j.resize(nn as usize, 0.0);
            oe.pollutant_amounts.clear();
            oe.pollutant_amounts.resize((nn as usize).wrapping_mul(np as usize), 0.0);
        }
        self.ctx_.negsrc.reset(); // D-NS1 (X6): clamp counters share the lifecycle
        for j in 0..self.ctx_.n_nodes() {
            // Legacy-convention node volume (junctions => 0)
            let v = self.reported_node_volume(j);
            self.ctx_.mass_balance.routing_init_storage += v;
        }
        for j in 0..self.ctx_.n_links() {
            let v = self.ctx_.links.volume[j as usize];
            self.ctx_.mass_balance.routing_init_storage += v;
        }

        // Initial stored pollutant mass, from the concentrations initQuality() seeded
        // (legacy massbal_open), with the same legacy-convention node volume.
        {
            let np = self.ctx_.n_pollutants();
            for p in 0..np {
                let up = p as usize;
                if up >= self.ctx_.mass_balance.qual_routing_init.len() {
                    break;
                }
                let mut m = 0.0;
                for j in 0..self.ctx_.n_nodes() {
                    let idx = (j as usize) * (np as usize) + up;
                    if idx < self.ctx_.nodes.conc.len() {
                        m += self.ctx_.nodes.conc[idx] * self.reported_node_volume(j) + lidnode::held_mass(&self.ctx_, j, p);
                    }
                }
                for j in 0..self.ctx_.n_links() {
                    let uj = j as usize;
                    let idx = uj * (np as usize) + up;
                    if idx < self.ctx_.links.conc.len() {
                        m += self.ctx_.links.conc[idx] * self.ctx_.links.volume[uj];
                    }
                }
                self.ctx_.mass_balance.qual_routing_init[up] = m;
            }
        }

        // A classic FV restart can carry disconnected scalar inventory; it belongs to both
        // opening storage and the current dry-storage ledger. (The state accessor open()
        // installs; see the file notes.)
        {
            let held = self.quality_.fv_retained_pollutant_mass(&self.ctx_);
            let mb = &mut self.ctx_.mass_balance;
            let mut p = 0usize;
            while p < held.len() && p < mb.qual_routing_init.len() {
                mb.qual_routing_init[p] += held[p];
                mb.qual_routing_final_dry[p] = held[p];
                p += 1;
            }
        }

        // Initial runoff storage: ponded depth x subarea area in ft2. Only a hot start
        // begins with water ponded.
        for j in 0..self.ctx_.n_subcatches() {
            if self.owns_spatial_source(j) {
                continue;
            }
            let uj = j as usize;
            let v = self.ctx_.subcatches.ponded_depth[uj] * self.runoff_.soa().area[uj];
            self.ctx_.mass_balance.runoff_init_store += v;
        }
        // Initial water in LID units is part of runoff storage.
        self.ctx_.mass_balance.runoff_init_store += self.lid_.stored_volume(None);

        // F8: the pack the deck starts with (legacy massbal.c:124-129).
        self.ctx_.mass_balance.runoff_init_snow = snow_cover_volume_ft3(&self.ctx_, self.snow_.state());

        // Initial groundwater storage: theta * upperDepth + porosity * lowerDepth, per area.
        {
            let gw = self.groundwater_.state();
            let ctx = &mut self.ctx_;
            for i in 0..ctx.n_subcatches() {
                let ui = i as usize;
                if gw.total_depth[ui] <= 0.0 {
                    continue;
                }
                let upper_d = gw.total_depth[ui] - gw.lower_depth[ui];
                let vol = gw.theta[ui] * upper_d + gw.porosity[ui] * gw.lower_depth[ui];
                let area = ctx.subcatches.area[ui] / ucf::ucf(LANDAREA, &ctx.options);
                ctx.mass_balance.gw_init_storage += vol * area;
            }
        }
    }
}
