// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Routing.hpp and Routing.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The top-level routing orchestrator's state and its initialization: RouteModel,
// applyConduitLengthening, Router::init (with the file-local buildDivertLinks),
// Router::initFv, Router::enableLinkStepStats, and the small inline accessors
// (hasCycle, setDWNumThreads, setStructureSolver, lastStepConverged, fv* getters).
//
// NOT YET TRANSLATED (run-time step functions, with the routing step): step, stepImpl,
// getAdaptiveStep, fvStepSink, accumulateLinkStepStats, initNodeFlows,
// computeConduitLosses, executeSteadyFlow, updateLinkStates, stepFv,
// accumulateStructureFlows, refreshFvBoundaryFlows, publishFv, acceptedQualityWater,
// restoreFvClassicState, privateDynamicWaveWindow, publishDynamicWaveWindow, the
// checkpoint session (enableCheckpointing, beginCheckpointStep,
// acceptCheckpointBoundary, saveCheckpoint, loadCheckpoint and their private helpers in
// RoutingCheckpoint.cpp), and the file-local buildXSP of Routing.cpp (used by the step).
//
// Differences in form, never in result:
//   * Router::init's local translateShape lambda is link::translate_shape: the two maps
//     are identical, and the lambda's `default:` arm cannot be reached with Rust's
//     closed XsectShape.
//   * `SimulationContext* checkpoint_context_` (a non-owning identity pointer, reset to
//     null by init) is `checkpoint_context_attached_: bool`.
//   * `hydstruct::StructureSolver* structures_` (non-owning, set by SWMMEngine to its own
//     hydstruct_) is `structures_attached_: bool`: the structure solver stays an engine
//     member, and the KW/STEADY step will take it as an argument.
//   * `std::unique_ptr<fv::INetworkSolver> fv_solver_` is an Option of the one backend the
//     1D reference build can construct (see hydraulics::fv).
//   * Members keep upstream's names, trailing underscores included, and are `pub` so the
//     gate's state dumper can read them, as upstream's harness does through
//     `#define private public`.
//   * applyConduitLengthening's static_cast<forcemain::FrictionModel>(force_main_eqn)
//     with an out-of-range value makes upstream's getRoughFactor return 0.0; so here.

use crate::core::cmath::{cmax, cmin};
use crate::core::constants;
use crate::core::context::SimulationContext;
use crate::core::options::{fv::FvOptions, Bytes};
use crate::data::link::{LinkData, LinkType, XsectShape};
use crate::data::node::NodeType;
use crate::hydraulics::dynamic_wave::{DWSolver, SurchargeMethod};
use crate::hydraulics::fv::{ExplicitFvSolver, NetworkMeshData, NetworkStateData};
use crate::hydraulics::kinematic_wave::KWSolver;
use crate::hydraulics::xsect_batch::{XSectGroups, XSectParams};
use crate::hydraulics::{force_main, link, outfall, topo_sort, xsection as xsect};
use crate::input::post_parse_resolver::{conduit_manning_n, conduit_true_length};

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RouteModel {
    /// Steady-state (pass-through)
    Steady = 0,
    /// Kinematic wave
    Kinwave = 1,
    /// Dynamic wave (St. Venant, implicit Picard)
    #[default]
    Dynwave = 2,
    /// Explicit conservative finite volume (Godunov)
    Fv = 3,
}

/// forcemain::getRoughFactor with upstream's static_cast of an int model.
fn rough_factor_for(eqn: i32, r_bot: f64, factor: f64) -> f64 {
    match eqn {
        0 => force_main::get_rough_factor(force_main::FrictionModel::HazenWilliams, r_bot, factor),
        1 => force_main::get_rough_factor(force_main::FrictionModel::DarcyWeisbach, r_bot, factor),
        _ => 0.0,
    }
}

/// Refresh routing lengths and apply dynamic-wave/FV Courant lengthening, with the
/// conveyance re-derivation that goes with it (legacy conduit_validate).
pub fn apply_conduit_lengthening(ctx: &mut SimulationContext, model: RouteModel) {
    let n_links = ctx.n_links();

    // Refresh the routing length (legacy conduit_getLength).
    for j in 0..n_links {
        if ctx.links.type_[j as usize] != LinkType::Conduit {
            continue;
        }
        let cr = ctx.link_subtypes.conduit_row(j);
        if cr < 0 {
            continue;
        }
        let tl = conduit_true_length(ctx, j);
        ctx.link_subtypes.conduits.true_length[cr as usize] = tl;
    }

    // The Courant lengthening factor per conduit row (legacy's local `lengthFactor`).
    let mut courant_factor = vec![1.0f64; ctx.link_subtypes.conduits.count() as usize];

    {
        use constants::{GRAVITY, PHI};
        let route_step = ctx.options.link_validate_routing_step();
        let lengthening_step = ctx.options.lengthening_step;
        let dw_like = model == RouteModel::Dynwave || model == RouteModel::Fv;

        if lengthening_step <= 0.0 || !dw_like {
            let cd = &mut ctx.link_subtypes.conduits;
            for r in 0..cd.count() {
                let ur = r as usize;
                cd.mod_length[ur] = cd.length[ur];
            }
        } else {
            let t_step = cmin(route_step, lengthening_step);

            for j in 0..n_links {
                let uj = j as usize;
                if ctx.links.type_[uj] != LinkType::Conduit {
                    continue;
                }
                let cr = ctx.link_subtypes.conduit_row(j);
                let ucr = cr as usize;

                let l = ctx.link_subtypes.conduits.true_length[ucr];
                if l <= 0.0 {
                    ctx.link_subtypes.conduits.mod_length[ucr] = ctx.link_subtypes.conduits.length[ucr];
                    continue;
                }

                let mut y_full = ctx.links.xsect_y_full[uj];
                let a_full = ctx.links.xsect_a_full[uj];
                let s_full = ctx.links.xsect_s_full[uj];
                let n_rough = conduit_manning_n(ctx, j);
                let slope_abs = ctx.link_subtypes.conduits.slope[ucr].abs();

                if xsect::is_open(ctx.links.xsect_batch_shape[uj]) {
                    let w_full = ctx.links.xsect_w_max[uj];
                    if w_full > 0.0 {
                        y_full = a_full / w_full;
                    }
                }

                let cd = &mut ctx.link_subtypes.conduits;
                if a_full > 0.0 && n_rough > 0.0 && slope_abs > 0.0 {
                    let v_full = PHI / n_rough * s_full * slope_abs.sqrt() / a_full;
                    let ratio = ((GRAVITY * y_full).sqrt() + v_full) * t_step / l;
                    let factor = if ratio > 1.0 { ratio } else { 1.0 };
                    courant_factor[ucr] = factor;
                    cd.mod_length[ucr] = if factor != 1.0 { factor * l } else { cd.length[ucr] };
                } else {
                    cd.mod_length[ucr] = cd.length[ucr];
                }
            }
        }
    }

    // Recompute conveyance properties accounting for lengthening.
    {
        use constants::{GRAVITY, PHI};
        for j in 0..n_links {
            let uj = j as usize;
            if ctx.links.type_[uj] != LinkType::Conduit {
                continue;
            }
            let cr = ctx.link_subtypes.conduit_row(j);
            let ucr = cr as usize;

            let factor = courant_factor[ucr];
            if factor == 1.0 {
                continue; // not lengthened
            }

            let slope_abs = ctx.link_subtypes.conduits.slope[ucr].abs() / factor;
            let roughness = conduit_manning_n(ctx, j) / factor.sqrt();

            let beta = PHI * slope_abs.sqrt() / roughness;
            let cd = &mut ctx.link_subtypes.conduits;
            cd.beta[ucr] = beta;
            cd.rough_factor[ucr] = GRAVITY * ((roughness / PHI) * (roughness / PHI));
            cd.q_full[ucr] = ctx.links.xsect_s_full[uj] * beta;
            cd.q_max[ucr] = ctx.links.xsect_s_max[uj] * beta;

            if ctx.links.xsect_shape[uj] == XsectShape::ForceMain && (model == RouteModel::Dynwave || model == RouteModel::Fv) {
                ctx.links.xsect_s_bot[uj] = rough_factor_for(ctx.options.force_main_eqn, ctx.links.xsect_r_bot[uj], factor);
            }
        }
    }
}

/// Per-node diversion-link index for toposort::sortLinks (legacy adjustAdjList).
fn build_divert_links(ctx: &SimulationContext) -> Vec<i32> {
    let mut divert = vec![-1i32; ctx.n_nodes() as usize];
    let divs = &ctx.node_subtypes.dividers;
    for i in 0..ctx.n_nodes() {
        if ctx.nodes.type_[i as usize] != NodeType::Divider {
            continue;
        }
        let r = ctx.node_subtypes.divider_row(i);
        if r < 0 {
            continue;
        }
        divert[i as usize] = divs.link[r as usize];
    }
    divert
}

/// The routing orchestrator: owns the XSectGroups, the KW and DW solvers and the FV state.
#[derive(Debug, Default)]
pub struct Router {
    pub checkpoint_context_attached_: bool,
    pub checkpoint_owned_session_: bool,
    pub checkpoint_aquifer_session_: bool,
    pub checkpoint_model_: Bytes,
    pub checkpoint_time_: f64,
    pub checkpoint_step_end_: f64,
    pub checkpoint_ready_: bool,
    pub checkpoint_prepared_: bool,
    pub checkpoint_step_finished_: bool,

    pub model_: RouteModel,
    pub groups_: XSectGroups,

    /// upstream `structures_ != nullptr` (see the file notes)
    pub structures_attached_: bool,

    pub steady_storage_updated_: Vec<i8>,
    pub steady_y_: Vec<f64>,

    // --- Explicit finite-volume solver (FLOW_ROUTING FV)
    pub fv_mesh_: NetworkMeshData,
    pub fv_state_: NetworkStateData,
    pub fv_opts_: FvOptions,
    pub fv_solver_: Option<ExplicitFvSolver>,
    pub fv_backend_: Bytes,
    pub fv_warnings_: Vec<Bytes>,
    pub fv_errors_: Vec<Bytes>,

    // --- Per-conduit time step summary ([REPORT] LINK_STEPS)
    pub lstep_enabled_: bool,
    pub lstep_rs_: f64,
    pub lstep_pmin_: Vec<f64>,
    pub lstep_pmax_: Vec<f64>,
    pub lstep_pdt_time_: Vec<f64>,
    pub lstep_ptime_: Vec<f64>,
    pub lstep_pbin_: Vec<f64>,

    pub fv_lateral_: Vec<f64>,
    pub fv_fixed_head_: Vec<f64>,
    pub fv_struct_flow_: Vec<f64>,
    pub fv_cond_loss_: Vec<f64>,
    pub fv_link_q_cap_: Vec<f64>,
    pub fv_struct_int_: Vec<f64>,
    pub fv_struct_t_prev_: f64,

    pub kw_solver_: KWSolver,
    pub dw_solver_: DWSolver,
    pub steady_sorted_links_: Vec<i32>,
    pub cycle_detected_: bool,
}

impl Router {
    /// Initialize routing groups, solver state and the selected model dispatch.
    pub fn init(&mut self, ctx: &mut SimulationContext, model: RouteModel) {
        self.checkpoint_context_attached_ = false;
        self.checkpoint_owned_session_ = false;
        self.checkpoint_model_.clear();
        self.checkpoint_step_end_ = 0.0;
        self.checkpoint_time_ = 0.0;
        self.checkpoint_step_finished_ = false;
        self.checkpoint_prepared_ = false;
        self.checkpoint_ready_ = false;
        self.model_ = model;

        let n_links = ctx.n_links();
        let n_nodes = ctx.n_nodes();

        // Per-link XSectParams from the link store (no tables attached: plain geometry).
        let mut xsect_params: Vec<XSectParams<'static>> = vec![XSectParams::default(); n_links as usize];
        for j in 0..n_links {
            let uj = j as usize;
            if ctx.links.type_[uj] != LinkType::Conduit {
                continue;
            }

            let xs = &mut xsect_params[uj];
            xs.type_ = link::translate_shape(ctx.links.xsect_shape[uj]);
            ctx.links.xsect_batch_shape[uj] = xs.type_;
            xs.y_full = ctx.links.xsect_y_full[uj];
            xs.a_full = ctx.links.xsect_a_full[uj];
            xs.w_max = ctx.links.xsect_w_max[uj];
            xs.r_full = ctx.links.xsect_r_full[uj];
            xs.s_full = ctx.links.xsect_s_full[uj];

            xs.y_bot = ctx.links.xsect_y_bot[uj];
            xs.a_bot = ctx.links.xsect_a_bot[uj];
            xs.s_bot = ctx.links.xsect_s_bot[uj];
            xs.r_bot = ctx.links.xsect_r_bot[uj];

            // Compute full-depth properties if not already set
            if xs.a_full == 0.0 && xs.y_full > 0.0 {
                let p = [xs.y_full, xs.w_max, 0.0, 0.0];
                let t = xs.type_;
                xsect::set_params(xs, t, &p, 1.0);
                ctx.links.xsect_a_full[uj] = xs.a_full;
                ctx.links.xsect_r_full[uj] = xs.r_full;
                ctx.links.xsect_s_full[uj] = xs.s_full;
                ctx.links.xsect_w_max[uj] = xs.w_max;
            }
        }

        apply_conduit_lengthening(ctx, model);

        // Build shape-grouped batch index
        self.groups_.build(&xsect_params, n_links);

        // Attach transect tables for IRREGULAR cross-sections
        self.groups_.attach_transect_tables(ctx);

        // Cache outfall -> connecting-conduit mapping.
        outfall::build_outfall_link_map(ctx);

        // Init solvers
        self.cycle_detected_ = false;
        match self.model_ {
            RouteModel::Kinwave => {
                self.kw_solver_.init(n_links);
                let mut sorted: Vec<i32> = Vec::new();
                let divert = build_divert_links(ctx);
                let n_sorted = topo_sort::sort_links(&ctx.links.node1, &ctx.links.node2, n_links, n_nodes, &mut sorted, Some(&divert));
                if n_sorted < n_links {
                    self.cycle_detected_ = true;
                }
                self.kw_solver_.set_link_order(&sorted);
            }
            RouteModel::Dynwave => {
                let ucf_len = crate::core::unit_conversion::UCF[crate::core::unit_conversion::LENGTH as usize]
                    [crate::core::unit_conversion::get_unit_system(ctx.options.flow_units as i32) as usize];
                self.dw_solver_.head_tol = if ctx.options.head_tol == 0.0 { constants::DEFAULT_HEAD_TOL } else { ctx.options.head_tol / ucf_len };
                self.dw_solver_.max_trials = if ctx.options.max_trials > 0 { ctx.options.max_trials } else { constants::DEFAULT_MAX_TRIALS };
                self.dw_solver_.surcharge_method = SurchargeMethod::from_i32(ctx.options.surcharge_method);
                self.dw_solver_.node_continuity = ctx.options.node_continuity;
                self.dw_solver_.anderson_accel = ctx.options.anderson_accel;
                self.dw_solver_.unsteady_friction = ctx.options.unsteady_friction;
                self.dw_solver_.uf_k3 = ctx.options.uf_k3;
                self.dw_solver_.tpa_celerity = ctx.options.tpa_celerity;

                self.dw_solver_.init(n_nodes, n_links, &self.groups_, ctx);
            }
            RouteModel::Fv => {
                self.init_fv(ctx);
            }
            RouteModel::Steady => {
                let divert = build_divert_links(ctx);
                let n_sorted = topo_sort::sort_links(&ctx.links.node1, &ctx.links.node2, n_links, n_nodes, &mut self.steady_sorted_links_, Some(&divert));
                if n_sorted < n_links {
                    self.cycle_detected_ = true;
                }
            }
        }
    }

    /// Enable conduit-step reporting when the selected solver can provide it.
    pub fn enable_link_step_stats(&mut self, ctx: &SimulationContext) -> bool {
        self.lstep_rs_ = ctx.options.routing_step;
        self.lstep_enabled_ = false;
        if self.model_ == RouteModel::Dynwave {
            self.lstep_enabled_ = true;
        } else if self.model_ == RouteModel::Fv {
            if self.fv_solver_.is_some() {
                let n = ctx.n_links() as usize;
                let nb = LinkData::N_LSTEP_BINS as usize;
                for v in [&mut self.lstep_pmin_, &mut self.lstep_pmax_, &mut self.lstep_pdt_time_, &mut self.lstep_ptime_] {
                    v.clear();
                    v.resize(n, 0.0);
                }
                self.lstep_pbin_.clear();
                self.lstep_pbin_.resize(n * nb, 0.0);
                // impl->setConduitStepSink(&Router::fvStepSink, this): the sink is the
                // step's (fvStepSink is translated with it); here it is only marked set.
                if let Some(s) = self.fv_solver_.as_mut() {
                    s.step_sink_set_ = true;
                }
                self.lstep_enabled_ = true;
            }
        }
        self.lstep_enabled_
    }

    /// Build the FV mesh, pick a backend and seed the solver state.
    pub fn init_fv(&mut self, ctx: &mut SimulationContext) {
        use crate::core::error_codes::{format_warning, WARN_DW_OPTION_UNDER_FV, WARN_FV_OPTION_INERT};
        use crate::core::options::fv::ScalarScheme;
        use crate::core::unit_conversion as ucf;
        self.fv_warnings_.clear();
        self.fv_errors_.clear();

        // FV_* keys are entered in project display units; convert to internal feet.
        let ucf_len = ucf::UCF[ucf::LENGTH as usize][ucf::get_unit_system(ctx.options.flow_units as i32) as usize];
        self.fv_opts_ = ctx.options.fv.clone();
        self.fv_opts_.cell_length = ctx.options.fv.cell_length / ucf_len;
        self.fv_opts_.slot_celerity = ctx.options.fv.slot_celerity / ucf_len;
        self.fv_opts_.dispersion = ctx.options.fv.dispersion / (ucf_len * ucf_len);
        self.fv_opts_.unsteady_friction = ctx.options.unsteady_friction;
        self.fv_opts_.uf_k3 = ctx.options.uf_k3;

        if ctx.options.fv.dispersion != 0.0 {
            self.fv_warnings_.push(format_warning(WARN_FV_OPTION_INERT, b"FV_DISPERSION"));
        }
        if ctx.options.fv.scalar_scheme != ScalarScheme::Muscl {
            self.fv_warnings_.push(format_warning(WARN_FV_OPTION_INERT, b"FV_SCALAR_SCHEME"));
        }
        if ctx.options.inertial_damping != 1 {
            self.fv_warnings_.push(format_warning(WARN_DW_OPTION_UNDER_FV, b"INERTIAL_DAMPING"));
        }
        if ctx.options.surcharge_method != 0 {
            self.fv_warnings_.push(format_warning(WARN_DW_OPTION_UNDER_FV, b"SURCHARGE_METHOD"));
        }
        if ctx.options.normal_flow_ltd != 2 {
            self.fv_warnings_.push(format_warning(WARN_DW_OPTION_UNDER_FV, b"NORMAL_FLOW_LIMITED"));
        }

        let rep = crate::hydraulics::fv::mesh_builder::build_network_mesh(ctx, &self.fv_opts_, &mut self.fv_mesh_);
        self.fv_warnings_.extend(rep.warnings.iter().cloned());
        self.fv_errors_ = rep.errors.clone();
        if !rep.errors.is_empty() {
            return;
        }

        let nc = self.fv_mesh_.n_cells();
        let nn = self.fv_mesh_.n_nodes();
        self.fv_state_.resize(nc, nn, 0);

        // Interpolated free surface at conduit-degree-2 virtual junctions (NaN = none).
        let mut vj_eta: Vec<f64> = Vec::new();
        {
            let n_nodes = ctx.n_nodes();
            let n_links = ctx.n_links();
            let mut any_vj = false;
            let mut n = 0;
            while n < n_nodes && !any_vj {
                any_vj = ctx.nodes.is_virtual[n as usize] != 0;
                n += 1;
            }
            if any_vj {
                vj_eta = vec![f64::NAN; n_nodes as usize]; // quiet_NaN
                let mut inc_a = vec![-1i32; n_nodes as usize];
                let mut inc_b = vec![-1i32; n_nodes as usize];
                let mut inc_n = vec![0i32; n_nodes as usize];
                for j in 0..n_links {
                    let ujj = j as usize;
                    if ctx.links.type_[ujj] != LinkType::Conduit {
                        continue;
                    }
                    for nd in [ctx.links.node1[ujj], ctx.links.node2[ujj]] {
                        if nd < 0 || nd >= n_nodes {
                            continue;
                        }
                        let und = nd as usize;
                        if ctx.nodes.is_virtual[und] == 0 {
                            continue;
                        }
                        if inc_n[und] == 0 {
                            inc_a[und] = j;
                        } else if inc_n[und] == 1 {
                            inc_b[und] = j;
                        }
                        inc_n[und] += 1;
                    }
                }
                // walk(start_node, via_link) -> Some((eta, dist, wet)) or None
                let walk = |start_node: i32, via_link: i32| -> Option<(f64, f64, bool)> {
                    let mut node = start_node;
                    let mut link = via_link;
                    let mut dist = 0.0;
                    for _guard in 0..=n_links {
                        let ul = link as usize;
                        let cr = ctx.link_subtypes.conduit_row(link);
                        dist += if cr >= 0 { ctx.link_subtypes.conduits.length[cr as usize] } else { 0.0 };
                        let other = if ctx.links.node1[ul] == node { ctx.links.node2[ul] } else { ctx.links.node1[ul] };
                        if other < 0 || other >= n_nodes {
                            return None;
                        }
                        let uo = other as usize;
                        if ctx.nodes.is_virtual[uo] == 0 {
                            return Some((ctx.nodes.invert_elev[uo] + ctx.nodes.depth[uo], dist, ctx.nodes.depth[uo] > 0.0));
                        }
                        if inc_n[uo] != 2 {
                            return None;
                        }
                        link = if inc_a[uo] == link { inc_b[uo] } else { inc_a[uo] };
                        if link < 0 {
                            return None;
                        }
                        node = other;
                    }
                    None // cycle of VJs (rejected elsewhere)
                };
                for i in 0..n_nodes {
                    let ui = i as usize;
                    if ctx.nodes.is_virtual[ui] == 0 || inc_n[ui] != 2 {
                        continue;
                    }
                    let ra = walk(i, inc_a[ui]);
                    let rb = walk(i, inc_b[ui]);
                    // upstream's out-parameters keep their 0.0 / false seeds on failure
                    let (ok_a, (e_a, d_a, wet_a)) = match ra {
                        Some(v) => (true, v),
                        None => (false, (0.0, 0.0, false)),
                    };
                    let (ok_b, (e_b, d_b, wet_b)) = match rb {
                        Some(v) => (true, v),
                        None => (false, (0.0, 0.0, false)),
                    };
                    if !(ok_a && wet_a) && !(ok_b && wet_b) {
                        continue;
                    }
                    if ok_a && ok_b {
                        let d_sum = d_a + d_b;
                        vj_eta[ui] = if d_sum > 0.0 { (e_a * d_b + e_b * d_a) / d_sum } else { 0.5 * (e_a + e_b) };
                    } else {
                        vj_eta[ui] = if ok_a { e_a } else { e_b };
                    }
                }
            }
        }

        // Seed cell state from the model's initial condition.
        for r in 0..self.fv_mesh_.n_conduits() {
            let ur = r as usize;
            let begin = self.fv_mesh_.conduit_cell_begin[ur];
            let count = self.fv_mesh_.conduit_cell_count[ur];
            if begin < 0 {
                continue;
            }
            let j = self.fv_mesh_.conduit_link[ur];
            let uj = j as usize;
            let g = &self.fv_mesh_.geom[self.fv_mesh_.conduit_section[ur] as usize];
            let q = ctx.links.flow[uj];

            let bank_eta = |nd: i32| -> Option<f64> {
                if nd < 0 {
                    return None;
                }
                let und = nd as usize;
                if und < vj_eta.len() && ctx.nodes.is_virtual[und] != 0 {
                    if !vj_eta[und].is_finite() || vj_eta[und] <= ctx.nodes.invert_elev[und] {
                        return None;
                    }
                    return Some(vj_eta[und]);
                }
                if ctx.nodes.depth[und] <= 0.0 {
                    return None;
                }
                Some(ctx.nodes.invert_elev[und] + ctx.nodes.depth[und])
            };
            let b1 = bank_eta(ctx.links.node1[uj]);
            let b2 = bank_eta(ctx.links.node2[uj]);
            let (wet1, eta1) = (b1.is_some(), b1.unwrap_or(0.0));
            let (wet2, eta2) = (b2.is_some(), b2.unwrap_or(0.0));
            let step = self.fv_mesh_.cell_dzdx[begin as usize] != 0.0;
            let shoreline = (wet1 != wet2) && step;
            let eta_wet = if wet1 { eta1 } else { eta2 };

            let mut init_depth = ctx.links.depth[uj];
            {
                let n1 = ctx.links.node1[uj];
                let n2 = ctx.links.node2[uj];
                let v1 = n1 >= 0 && (n1 as usize) < vj_eta.len() && ctx.nodes.is_virtual[n1 as usize] != 0;
                let v2 = n2 >= 0 && (n2 as usize) < vj_eta.len() && ctx.nodes.is_virtual[n2 as usize] != 0;
                if (v1 || v2) && n1 >= 0 && n2 >= 0 {
                    let end_depth = |nd: i32, v: bool| -> f64 {
                        let und = nd as usize;
                        if !v {
                            return ctx.nodes.depth[und];
                        }
                        if vj_eta[und].is_finite() {
                            cmax(0.0, vj_eta[und] - ctx.nodes.invert_elev[und])
                        } else {
                            0.0
                        }
                    };
                    init_depth = 0.5 * (end_depth(n1, v1) + end_depth(n2, v2));
                }
            }
            let uniform = crate::hydraulics::fv::kernels::area_of_depth(g, init_depth);

            let mut vol = 0.0;
            let mut a_len = 0.0;
            let mut sum_len = 0.0;
            let mut slot_vol = 0.0;
            for c in begin..begin + count {
                let uc = c as usize;
                let mut a_c = uniform;
                if shoreline {
                    a_c = crate::hydraulics::fv::kernels::area_of_depth(g, cmax(0.0, eta_wet - self.fv_mesh_.cell_zb[uc]));
                }
                self.fv_state_.cell_a[uc] = a_c;
                self.fv_state_.cell_q[uc] = if a_c > crate::hydraulics::fv::kernels::K_DRY_AREA { q } else { 0.0 };

                let dx = self.fv_mesh_.cell_dx[uc];
                sum_len += dx;
                a_len += a_c * dx;
                vol += a_c * dx;
                slot_vol += cmax(0.0, a_c - g.a_crown) * dx;
            }

            if sum_len > 0.0 {
                ctx.links.volume[uj] = vol;
                ctx.links.slot_volume[uj] = slot_vol;
                ctx.links.depth[uj] = cmin(crate::hydraulics::fv::kernels::depth_of_area(g, a_len / sum_len), g.y_full);
            }
        }
        for n in 0..nn {
            let un = n as usize;
            self.fv_state_.node_head[un] = ctx.nodes.invert_elev[un] + ctx.nodes.depth[un];
        }

        // makeNetworkSolver: a 1D reference build has no GPU plugin library, so every
        // backend request ends at the CPU explicit solver (NetworkSolverFactory.cpp).
        self.fv_backend_ = b"cpu (explicit finite volume)".to_vec();
        self.fv_solver_ = Some(ExplicitFvSolver::default());
        // StorageArea / SUBSTEP coupling need the CPU backend, which this is: no throw, no
        // clamp warning.

        let (mesh, state) = (&mut self.fv_mesh_, &mut self.fv_state_);
        let solver = self.fv_solver_.as_mut().unwrap();
        solver.initialize(mesh, state, &self.fv_opts_);
        for n in 0..nn {
            let un = n as usize;
            ctx.nodes.volume[un] = if ctx.nodes.type_[un] == NodeType::Outfall { 0.0 } else { state.node_volume[un] };
        }

        let legacy_quality = !ctx.options.ignore_quality
            && ctx.options.quality_solver == crate::core::options::QualitySolverKind::Legacy
            && (ctx.n_pollutants() > 0 || ctx.reactions.n_species() > 0 || ctx.options.water_age || ctx.options.heat_transport);
        if std::env::var_os("OPENSWMM_FV_RECEIPT_TRACE").is_some() || legacy_quality {
            solver.track_accepted_water(mesh, true);
            solver.require_conservative_node_balance(legacy_quality);
        }

        let unn = nn as usize;
        let nl = ctx.n_links() as usize;
        self.fv_lateral_.clear();
        self.fv_lateral_.resize(unn, 0.0);
        self.fv_fixed_head_.clear();
        self.fv_fixed_head_.resize(unn, f64::NAN); // quiet_NaN
        self.fv_struct_flow_.clear();
        self.fv_struct_flow_.resize(nl, 0.0);
        self.fv_cond_loss_.clear();
        self.fv_cond_loss_.resize(self.fv_mesh_.n_conduits() as usize, 0.0);
        self.fv_struct_int_.clear();
        self.fv_struct_int_.resize(nl, 0.0);
        self.fv_link_q_cap_.clear();
        self.fv_link_q_cap_.resize(nl, -1.0);
    }

    pub fn has_cycle(&self) -> bool {
        self.cycle_detected_
    }

    pub fn set_dw_num_threads(&mut self, n: i32, warnings: Option<&mut Vec<Bytes>>) {
        self.dw_solver_.set_num_threads(n, warnings);
    }

    pub fn set_structure_solver_attached(&mut self, attached: bool) {
        self.structures_attached_ = attached;
    }

    pub fn last_step_converged(&self) -> bool {
        if self.model_ == RouteModel::Dynwave {
            self.dw_solver_.last_converged()
        } else {
            true
        }
    }

    pub fn link_step_stats_enabled(&self) -> bool {
        self.lstep_enabled_
    }

    pub fn xsect_groups(&self) -> &XSectGroups {
        &self.groups_
    }

    pub fn fv_warnings(&self) -> &[Bytes] {
        &self.fv_warnings_
    }

    pub fn fv_errors(&self) -> &[Bytes] {
        &self.fv_errors_
    }

    pub fn fv_backend(&self) -> &[u8] {
        &self.fv_backend_
    }
}
