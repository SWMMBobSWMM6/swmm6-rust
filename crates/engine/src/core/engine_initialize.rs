// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/SWMMEngine.cpp (SWMMEngine::initialize,
// init_modules, reopenMassBalance, applyLegacyRunoffSection, computeFinalQualityMassBalance,
// and the file-local infilFromLegacy, buildupNorm and pollutantSurfaceMcf) and
// SWMMEngine.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// SWMMEngine::initialize (1D build: OPENSWMM_HAS_2D off): the state checks, the state-seeding
// phase (MIN_SURFAREA junction storage, LID nodes, initial node volumes, the pump wet-well
// override, the legacy DW volume convention, conduit q0 and normal depths, initNodeDepths,
// outfall depths and the backwater fill), [FILES] USE HOTSTART (legacy SWMM5 files), the
// clocks, init_modules (the four init legs in upstream's order, with the hot-start runoff and
// species blocks between them, then ctx_.forcing.resize), and the node inflow/outflow seed.
// With the `twod` feature (the 2D build, 12 November 2026) the 1D build's refusals of
// REACTION_SUBSURFACE and MESH_CONTACTS channel heat are absent, init_modules calls
// initializeChannelGround after init_geometry, and the reactive-LID check follows it
// (core/engine_twod.rs).
//
// Differences of form, never in result:
//   * initialize is noexcept upstream. Where an exception escapes it (U1/U65 inside
//     initHydraulics, U57 inside initQuality, a non-finite Richards pressure inside
//     lidnode::initialize), upstream terminates; here initialize returns Err(Terminated) at
//     that point. Ok(rc) is upstream's return code.
//   * perf:: timers only time the run; they are left out.
//   * A native OPENSWMM_HS_V1 hot-start file (HotStartManager::open/apply, the engine
//     checkpoint) is NOT YET TRANSLATED: initialize returns Err(stub) when it meets one. No
//     gate deck uses one.
//   * The warning callbacks that push onto ctx_.warnings receive the context as an argument.

use crate::core::cmath::{cclamp, cmax};
use crate::core::constants;
use crate::core::context::EngineState;
#[cfg(not(feature = "twod"))]
use crate::core::engine::SWMM_ERR_PLUGIN;
use crate::core::engine::{SwmmEngine, CFFI_ERR_BADPARAM, CFFI_ERR_HOTSTART, SWMM_ERR_BADPARAM, SWMM_ERR_PARSE, SWMM_ERR_WRONG_STATE, SWMM_OK};
use crate::core::engine_init_hydraulics::Terminated;
use crate::core::engine_step::stub;
use crate::core::hot_start_manager as hsm;
use crate::core::options::{Bytes, QualitySolverKind, RoutingModel};
use crate::core::unit_conversion as ucf;
#[cfg(not(feature = "twod"))]
use crate::data::heat_override::HeatElement;
use crate::data::link::LinkType;
use crate::data::node::NodeType;
use crate::data::pollutant::MassUnits;
use crate::data::table::TableType;
use crate::hydraulics::dynamic_wave::SurchargeMethod;
use crate::hydraulics::routing::{apply_conduit_lengthening, RouteModel};
use crate::hydraulics::{link, node, outfall, xsection};
use crate::hydrology::infiltration::InfilModel;
use crate::hydrology::lid_node_held as lidnode;
use crate::quality::landuse::LanduseSolver;

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

/// pollutantSurfaceMcf: per-pollutant surface buildup mass conversion multipliers.
pub(crate) fn pollutant_surface_mcf(ctx: &crate::core::context::SimulationContext) -> Vec<f64> {
    let np = ctx.n_pollutants();
    let mut mcf = vec![0.0; if np > 0 { np as usize } else { 0 }];
    let mass_ucf = ucf::ucf(ucf::MASS, &ctx.options);
    for p in 0..np {
        let mut m = mass_ucf;
        if (p as usize) < ctx.pollutants.units.len() {
            match ctx.pollutants.units[p as usize] {
                MassUnits::UgPerL => m = mass_ucf / 1000.0,
                MassUnits::CountsPerL => m = 1.0,
                _ => {}
            }
        }
        mcf[p as usize] = m;
    }
    mcf
}

/// infilFromLegacy: legacy runoff-section infiltration slots into the native six.
fn infil_from_legacy(model: i32, lg: &[f64], v6: &mut [f64; 6]) {
    *v6 = [0.0; 6];
    if model == InfilModel::Horton as i32 || model == InfilModel::ModHorton as i32 {
        v6[0] = lg[0];
        v6[1] = lg[1];
        v6[2] = lg[2];
    } else if model == InfilModel::GreenAmpt as i32 || model == InfilModel::ModGreenAmpt as i32 {
        // legacy grnampt_setState: Sat = (char)x[3]
        v6[0] = lg[0];
        v6[1] = lg[1];
        v6[2] = lg[2];
        v6[3] = lg[4];
        v6[4] = if crate::core::cmath::cvt_i32(lg[3]) as i8 != 0 { 1.0 } else { 0.0 };
    } else if model == InfilModel::CurveNum as i32 {
        v6[0] = lg[0];
        v6[1] = lg[4];
        v6[2] = lg[1];
        v6[3] = lg[2];
        v6[4] = lg[5];
        v6[5] = lg[3];
    }
}

/// buildupNorm: coverage fraction times area (PER_AREA) or curb length.
pub(crate) fn buildup_norm(ctx: &crate::core::context::SimulationContext, ls: &LanduseSolver, ui: usize, lu: i32, p: i32) -> f64 {
    let nlu = ctx.n_landuses();
    let np = ctx.n_pollutants();
    let cov = ui * nlu as usize + lu as usize;
    let frac = if cov < ctx.subcatches.coverage.len() { ctx.subcatches.coverage[cov] / 100.0 } else { 0.0 };
    if frac <= 0.0 {
        return 0.0;
    }
    let bp = &ls.buildup_params[(lu * np + p) as usize];
    if bp.normalizer == 0 {
        frac * ctx.subcatches.area[ui]
    } else {
        frac * ctx.subcatches.curb_length[ui]
    }
}

impl SwmmEngine {
    /// SWMMEngine::initialize.
    pub fn initialize(&mut self) -> Result<i32, Terminated> {
        if self.metadata_only_ {
            self.set_error(SWMM_ERR_WRONG_STATE, b"Metadata-only handles cannot initialize; close and open the model for simulation");
            return Ok(SWMM_ERR_WRONG_STATE);
        }
        // INITIALIZED is accepted so a caller can re-initialize before start().
        if self.ctx_.state != EngineState::Opened && self.ctx_.state != EngineState::Initialized {
            self.set_error(SWMM_ERR_WRONG_STATE, b"swmm_engine_initialize: must call open() first");
            return Ok(SWMM_ERR_WRONG_STATE);
        }
        // perf::begin_process_model: timing only.

        // without OPENSWMM_HAS_2D (with the twod feature, the 2D build's own initialize-time
        // checks decide these decks: core/engine_twod.rs)
        #[cfg(not(feature = "twod"))]
        if !self.ctx_.reactions.subsurface_laws.is_empty() {
            self.set_error(SWMM_ERR_PLUGIN, b"REACTION_SUBSURFACE requires a 2D-enabled groundwater engine.");
            return Ok(SWMM_ERR_PLUGIN);
        }
        #[cfg(not(feature = "twod"))]
        if self.ctx_.options.heat_transport && self.ctx_.heat_config.sediment_exchange {
            for l in 0..self.ctx_.n_links() {
                if crate::transport::heat_overrides::sediment_for(&self.ctx_, &HeatElement::link(l)).mesh_ground {
                    self.set_error(CFFI_ERR_BADPARAM, b"MESH_CONTACTS channel heat requires an engine with 2D routing support.");
                    return Ok(CFFI_ERR_BADPARAM);
                }
            }
        }

        if self.ctx_.options.quality_solver == QualitySolverKind::Lagrangian && self.ctx_.reactions.species_is_wall.iter().any(|&w| w != 0) {
            self.set_error(SWMM_ERR_BADPARAM, b"LAGRANGIAN reaction transport does not support attached WALL species; use QUALITY_SOLVER LEGACY.");
            return Ok(SWMM_ERR_BADPARAM);
        }
        if self.ctx_.reactions.has_subcatch_reactions() && self.ctx_.lid_usage.count() > 0 && self.ctx_.reactions.lid_controls.is_empty() {
            self.set_error(SWMM_ERR_BADPARAM, b"[REACTION_SUBCATCHMENTS] wet-pool reactions require non-LID subcatchments; porous LID reaction transport is unsupported.");
            return Ok(SWMM_ERR_BADPARAM);
        }
        self.ctx_.refresh_reported_species();

        // ---- the state-seeding phase ----
        self.ctx_.reset_state();
        // These clocks belong to the engine rather than the context.
        self.old_runoff_time_ = 0.0;
        self.new_runoff_time_ = 0.0;
        self.old_runoff_ms_ = 0.0;
        self.new_runoff_ms_ = 0.0;
        self.new_rule_time_ms_ = 0.0;
        self.prev_runoff_step_sec_ = 0.0;
        self.has_runoff_ = false;
        self.has_snow_ = false;

        // MIN_SURFAREA junction storage (full_volume) and the junction StorageArea checks.
        {
            let ucf_len = ucf::UCF[ucf::LENGTH as usize][ucf::get_unit_system(self.ctx_.options.flow_units as i32) as usize];
            let min_sa = if self.ctx_.options.min_surf_area > 0.0 { self.ctx_.options.min_surf_area / (ucf_len * ucf_len) } else { constants::MIN_SURFAREA };
            for i in 0..self.ctx_.n_nodes() {
                let ui = i as usize;
                let n = &self.ctx_.nodes;
                if !n.storage_area[ui].is_finite() || n.storage_area[ui] < 0.0 || !(n.storage_area[ui] * n.full_depth[ui]).is_finite() {
                    self.set_error(SWMM_ERR_BADPARAM, b"Junction StorageArea and chamber capacity must be finite and nonnegative");
                    return Ok(SWMM_ERR_BADPARAM);
                }
                if n.storage_area[ui] > 0.0 {
                    let rm = self.ctx_.options.routing_model;
                    if n.type_[ui] != NodeType::Junction || n.is_virtual[ui] != 0 || n.is_inlet[ui] != 0 || rm == RoutingModel::Kinwave || rm == RoutingModel::Steady {
                        self.set_error(SWMM_ERR_BADPARAM, b"Junction StorageArea requires an ordinary junction and DYNWAVE or supported FV routing");
                        return Ok(SWMM_ERR_BADPARAM);
                    }
                    if rm == RoutingModel::Dynwave && self.ctx_.options.surcharge_method == SurchargeMethod::DynamicSlot as i32 {
                        self.set_error(SWMM_ERR_BADPARAM, b"Junction StorageArea is not qualified with DYNAMIC_SLOT; use EXTRAN or SLOT");
                        return Ok(SWMM_ERR_BADPARAM);
                    }
                    if rm == RoutingModel::Fv && n.sur_depth[ui] > 0.0 {
                        self.set_error(SWMM_ERR_BADPARAM, b"FV junction StorageArea does not support above-rim pressure allowance; use DYNWAVE or set SurDepth to zero");
                        return Ok(SWMM_ERR_BADPARAM);
                    }
                    let n = &mut self.ctx_.nodes;
                    n.full_volume[ui] = n.storage_area[ui] * n.full_depth[ui];
                    continue;
                }
                if n.type_[ui] == NodeType::Storage {
                    continue;
                }
                let fd = n.full_depth[ui];
                if fd > 0.0 && !(n.full_volume[ui] > 0.0) {
                    self.ctx_.nodes.full_volume[ui] = min_sa * fd;
                }
            }
        }

        crate::hydrology::lid_node::validate(&mut self.ctx_);
        if let Some(e) = self.ctx_.errors.last().cloned() {
            self.set_error(CFFI_ERR_BADPARAM, &e);
            return Ok(CFFI_ERR_BADPARAM);
        }
        crate::hydrology::lid_node_init::initialize(&mut self.ctx_)?;
        if let Some(e) = self.ctx_.errors.last().cloned() {
            self.set_error(CFFI_ERR_BADPARAM, &e);
            return Ok(CFFI_ERR_BADPARAM);
        }

        let us = ucf::get_unit_system(self.ctx_.options.flow_units as i32);

        // Compute initial volumes from init_depth (matching legacy node_initState)
        for i in 0..self.ctx_.n_nodes() {
            let ui = i as usize;
            let c = &mut self.ctx_;
            let mut d = c.nodes.init_depth[ui];
            if lidnode::active(c, i) {
                d = cmax(d, node::get_depth(&c.nodes, i, c.nodes.volume[ui], Some(&c.tables), us, Some(&c.node_subtypes)));
                c.nodes.depth[ui] = d;
                c.nodes.old_depth[ui] = d;
                c.nodes.head[ui] = c.nodes.invert_elev[ui] + d;
            }
            if d > 0.0 {
                let vol = node::get_volume(&c.nodes, i, d, Some(&c.tables), us, Some(&c.node_subtypes));
                c.nodes.volume[ui] = vol;
                c.nodes.old_volume[ui] = vol;
            }
            let fd = c.nodes.full_depth[ui];
            if fd > 0.0 {
                c.nodes.full_volume[ui] = node::get_volume(&c.nodes, i, fd, Some(&c.tables), us, Some(&c.node_subtypes));
            }
        }

        // Legacy-convention full volume (rpt_full_volume).
        for i in 0..self.ctx_.n_nodes() {
            let ui = i as usize;
            let n = &mut self.ctx_.nodes;
            n.rpt_full_volume[ui] = if n.type_[ui] == NodeType::Storage || n.storage_area[ui] > 0.0 { n.full_volume[ui] } else { 0.0 };
        }

        // Type-1 (volume-controlled) pumps: the inlet junction is a wet well of the curve's
        // maximum volume (legacy pump_validate).
        {
            let ucf_vol = ucf::UCF[ucf::VOLUME as usize][us as usize];
            for j in 0..self.ctx_.n_links() {
                let uj = j as usize;
                let c = &mut self.ctx_;
                if c.links.type_[uj] != LinkType::Pump {
                    continue;
                }
                let pr = c.link_subtypes.pump_row(j);
                let ci = if pr >= 0 { c.link_subtypes.pumps.curve[pr as usize] } else { -1 };
                if ci < 0 || ci >= c.tables.tables.len() as i32 {
                    continue;
                }
                let tbl = &c.tables.tables[ci as usize];
                if tbl.type_ != TableType::CurvePump1 {
                    continue;
                }
                let n1 = c.links.node1[uj];
                if n1 < 0 {
                    continue;
                }
                let un1 = n1 as usize;
                if c.nodes.type_[un1] == NodeType::Storage || c.nodes.storage_area[un1] > 0.0 {
                    continue;
                }
                let mut xmax = tbl.x_max;
                if xmax <= 0.0 && !tbl.x.is_empty() {
                    // std::max_element: the first of the largest
                    let mut m = tbl.x[0];
                    for &x in &tbl.x[1..] {
                        if m < x {
                            m = x;
                        }
                    }
                    xmax = m;
                }
                let xmax_internal = xmax / ucf_vol;
                c.nodes.full_volume[un1] = cmax(c.nodes.full_volume[un1], xmax_internal);
                c.nodes.rpt_full_volume[un1] = cmax(c.nodes.rpt_full_volume[un1], xmax_internal);
                let d0 = c.nodes.init_depth[un1];
                if d0 > 0.0 {
                    let v0 = node::get_volume(&c.nodes, n1, d0, Some(&c.tables), us, Some(&c.node_subtypes));
                    c.nodes.volume[un1] = v0;
                    c.nodes.old_volume[un1] = v0;
                }
            }
        }

        // DYNWAVE books non-storage node volume in the legacy convention.
        if self.ctx_.options.routing_model == RoutingModel::Dynwave {
            let nn = self.ctx_.n_nodes() as usize;
            let n = &mut self.ctx_.nodes;
            for ui in 0..nn {
                if n.type_[ui] == NodeType::Storage {
                    continue;
                }
                let fd = n.full_depth[ui];
                let d0 = n.init_depth[ui];
                let v0 = if fd > 0.0 && d0 > 0.0 { n.rpt_full_volume[ui] * (d0 / fd) } else { 0.0 };
                n.volume[ui] = v0;
                n.old_volume[ui] = v0;
            }
        }

        // Cache the batch cross-section code of every conduit before the initial-state loops.
        for j in 0..self.ctx_.n_links() {
            let uj = j as usize;
            if self.ctx_.links.type_[uj] != LinkType::Conduit {
                continue;
            }
            self.ctx_.links.xsect_batch_shape[uj] = link::translate_shape(self.ctx_.links.xsect_shape[uj]);
        }

        // Lengthen the conduits and re-derive their conveyance before the initial state.
        {
            let mut rm0 = RouteModel::Dynwave;
            if self.ctx_.options.routing_model == RoutingModel::Kinwave {
                rm0 = RouteModel::Kinwave;
            } else if self.ctx_.options.routing_model == RoutingModel::Steady {
                rm0 = RouteModel::Steady;
            } else if self.ctx_.options.routing_model == RoutingModel::Fv {
                rm0 = RouteModel::Fv;
            }
            apply_conduit_lengthening(&mut self.ctx_, rm0);
        }

        // legacy link_initState / conduit_initState: q0 and the conduit's normal depth.
        let start_date0 = self.ctx_.options.start_date.floor();
        for j in 0..self.ctx_.n_links() {
            let uj = j as usize;
            let c = &mut self.ctx_;
            let q0 = c.links.q0[uj];
            c.links.flow[uj] = q0;
            c.links.old_flow[uj] = q0;
            c.links.time_last_set[uj] = start_date0;
            if c.links.type_[uj] == LinkType::Conduit && q0 != 0.0 {
                let cr = c.link_subtypes.conduit_row(j);
                let cd = &c.link_subtypes.conduits;
                let (y, vol) = {
                    let xs = link::build_xsect_params(&c.links, uj, Some(&c.transect_tables));
                    let barrels = if cr >= 0 { cd.barrels[cr as usize] } else { 1 };
                    let q_per_barrel = q0.abs() / std::cmp::max(barrels, 1) as f64;
                    let beta = if cr >= 0 { cd.beta[cr as usize] } else { 0.0 };
                    let q_max = if cr >= 0 { cd.q_max[cr as usize] } else { -1.0 };
                    let y = link::get_depth_from_flow(&xs, beta, q_per_barrel, q_max);
                    let len = if cr >= 0 { cd.true_length[cr as usize] } else { 0.0 };
                    let vol = xsection::get_a_of_y(&xs, y) * len * barrels as f64;
                    (y, vol)
                };
                c.links.depth[uj] = y;
                c.links.old_depth[uj] = y;
                c.links.volume[uj] = vol;
                c.links.old_volume[uj] = vol;
            }
        }

        // initNodeDepths (legacy flowrout.c): the average of the connecting links' depths.
        {
            let nn = self.ctx_.n_nodes();
            let mut acc = vec![0.0f64; nn as usize];
            let mut cnt = vec![0i32; nn as usize];
            for j in 0..self.ctx_.n_links() {
                let uj = j as usize;
                let l = &self.ctx_.links;
                let ld = l.depth[uj];
                let y = if ld > constants::FUDGE { ld + l.offset1[uj] } else { 0.0 };
                let n1 = l.node1[uj];
                let n2 = l.node2[uj];
                if n1 >= 0 {
                    acc[n1 as usize] += y;
                    cnt[n1 as usize] += 1;
                }
                if n2 >= 0 {
                    acc[n2 as usize] += y;
                    cnt[n2 as usize] += 1;
                }
            }
            for i in 0..nn {
                let ui = i as usize;
                let c = &mut self.ctx_;
                if c.nodes.type_[ui] == NodeType::Outfall {
                    continue;
                }
                if c.nodes.type_[ui] == NodeType::Storage {
                    continue;
                }
                if c.nodes.init_depth[ui] > 0.0 {
                    continue;
                }
                if cnt[ui] <= 0 {
                    continue;
                }
                let y = acc[ui] / cnt[ui] as f64;
                if y <= 0.0 {
                    continue;
                }
                c.nodes.depth[ui] = y;
                c.nodes.old_depth[ui] = y;
                c.nodes.head[ui] = c.nodes.invert_elev[ui] + y;
                let vol;
                if c.options.routing_model == RoutingModel::Dynwave {
                    let fd = c.nodes.full_depth[ui];
                    vol = if fd > 0.0 { c.nodes.rpt_full_volume[ui] * (y / fd) } else { 0.0 };
                } else {
                    vol = node::get_volume(&c.nodes, i, y, Some(&c.tables), us, None);
                }
                c.nodes.volume[ui] = vol;
                c.nodes.old_volume[ui] = vol;
            }
        }

        // Outfall boundary depths, then the backwater fill of the q0 = 0 conduits.
        outfall::build_outfall_link_map(&mut self.ctx_);
        outfall::set_all_outfall_depths(&mut self.ctx_, 0.0);
        for j in 0..self.ctx_.n_links() {
            let uj = j as usize;
            let c = &mut self.ctx_;
            if c.links.type_[uj] != LinkType::Conduit {
                continue;
            }
            if c.links.q0[uj] != 0.0 {
                continue;
            }
            let yfull = c.links.xsect_y_full[uj];
            let n1 = c.links.node1[uj];
            let n2 = c.links.node2[uj];
            if n1 < 0 || n2 < 0 {
                continue;
            }
            let y1 = cclamp(c.nodes.depth[n1 as usize] - c.links.offset1[uj], 0.0, yfull);
            let y2 = cclamp(c.nodes.depth[n2 as usize] - c.links.offset2[uj], 0.0, yfull);
            let y = cmax(0.5 * (y1 + y2), constants::FUDGE);
            c.links.depth[uj] = y;
            c.links.old_depth[uj] = y;
            let vol = {
                let xs = link::build_xsect_params(&c.links, uj, Some(&c.transect_tables));
                let cr = c.link_subtypes.conduit_row(j);
                let cd = &c.link_subtypes.conduits;
                let barrels = std::cmp::max(if cr >= 0 { cd.barrels[cr as usize] } else { 1 }, 1);
                let len = if cr >= 0 { cd.true_length[cr as usize] } else { 0.0 };
                xsection::get_a_of_y(&xs, y) * len * barrels as f64
            };
            c.links.volume[uj] = vol;
            c.links.old_volume[uj] = vol;
        }

        // ---- USE HOTSTART: load saved routing state from a hot-start file ----
        let hs_path: Bytes = if !self.ctx_.files.hotstart_use_path.absolute.is_empty() {
            self.ctx_.files.hotstart_use_path.absolute.clone()
        } else {
            self.ctx_.files.hotstart_use_path.original.clone()
        };
        let lid_restart = self.ctx_.node_subtypes.storages.lid.iter().any(|c| c.control >= 0);
        let hs_native = !self.ctx_.files.hotstart_use_path.is_empty() && hsm::is_native_file(&hs_path);
        let full_engine_restore = false;
        let deferred_hotstart = false;
        self.hotstart_deferred_ = false;
        if !self.ctx_.files.hotstart_use_path.is_empty() {
            let rc;
            if hs_native {
                // HotStartManager::open / engine checkpoint / deferred apply
                return Err(stub("HotStartManager native USE HOTSTART file"));
            } else {
                if lid_restart {
                    self.ctx_.warnings.push(cat(&[
                        b"USE HOTSTART: '",
                        &hs_path,
                        b"' is a legacy SWMM5 file, which carries no LID state; LID units start from their initial conditions.",
                    ]));
                }
                let mut runoff = std::mem::take(&mut self.hotstart_runoff_);
                let mut species = std::mem::take(&mut self.hotstart_species_);
                rc = hsm::apply_legacy_routing(&hs_path, &mut self.ctx_, &mut |c, m| c.warnings.push(m), Some(&mut runoff), Some(&mut species));
                self.hotstart_runoff_ = runoff;
                self.hotstart_species_ = species;
            }
            if rc != 0 {
                let m = cat(&[b"USE HOTSTART: ", &hsm::last_io_error()]);
                self.set_error(CFFI_ERR_HOTSTART, &m);
                return Ok(CFFI_ERR_HOTSTART);
            }
            if !full_engine_restore && !deferred_hotstart {
                // Recompute derived state from the applied depths/flows.
                for i in 0..self.ctx_.n_nodes() {
                    let ui = i as usize;
                    let c = &mut self.ctx_;
                    let y = c.nodes.depth[ui];
                    c.nodes.old_depth[ui] = y;
                    c.nodes.head[ui] = c.nodes.invert_elev[ui] + y;
                    let vol = node::get_volume(&c.nodes, i, y, Some(&c.tables), us, Some(&c.node_subtypes));
                    c.nodes.volume[ui] = vol;
                    c.nodes.old_volume[ui] = vol;
                }
                for j in 0..self.ctx_.n_links() {
                    let uj = j as usize;
                    let c = &mut self.ctx_;
                    c.links.old_flow[uj] = c.links.flow[uj];
                    if c.links.type_[uj] == LinkType::Conduit {
                        let y = c.links.depth[uj];
                        c.links.old_depth[uj] = y;
                        let vol = {
                            let xs = link::build_xsect_params(&c.links, uj, Some(&c.transect_tables));
                            let cr = c.link_subtypes.conduit_row(j);
                            let cd = &c.link_subtypes.conduits;
                            let barrels = std::cmp::max(if cr >= 0 { cd.barrels[cr as usize] } else { 1 }, 1);
                            let len = if cr >= 0 { cd.true_length[cr as usize] } else { 0.0 };
                            xsection::get_a_of_y(&xs, y) * len * barrels as f64
                        };
                        c.links.volume[uj] = vol;
                        c.links.old_volume[uj] = vol;
                    }
                }
            }
        }
        // Output timer + legacy-parity ms clocks
        self.ctx_.current_date = self.ctx_.options.start_date;
        self.ctx_.current_time = 0.0;
        self.ctx_.elapsed_ms = 0.0;
        self.ctx_.old_elapsed_ms = 0.0;
        self.ctx_.next_report_ms = 1000.0 * self.ctx_.options.report_step;
        self.new_rule_time_ms_ = 0.0;

        // Initialize all computational modules (batch SoA setup)
        self.init_modules()?;

        #[cfg(feature = "twod")]
        self.initialize_2d_reactive_lid_check();
        if let Some(e) = self.ctx_.errors.first().cloned() {
            self.set_error(SWMM_ERR_PARSE, &e);
            return Ok(SWMM_ERR_PARSE);
        }

        // Seed node inflow/outflow from the lateral flows and the initial link flows
        // (legacy initNodes, flowrout.c:440-473).
        for i in 0..self.ctx_.n_nodes() {
            let ui = i as usize;
            self.ctx_.nodes.inflow[ui] += self.ctx_.nodes.lat_flow[ui];
        }
        for j in 0..self.ctx_.n_links() {
            let uj = j as usize;
            let c = &mut self.ctx_;
            let q = c.links.flow[uj];
            let n1 = c.links.node1[uj];
            let n2 = c.links.node2[uj];
            if q >= 0.0 {
                if n1 >= 0 {
                    c.nodes.outflow[n1 as usize] += q;
                }
                if n2 >= 0 {
                    c.nodes.inflow[n2 as usize] += q;
                }
            } else {
                if n1 >= 0 {
                    c.nodes.inflow[n1 as usize] -= q;
                }
                if n2 >= 0 {
                    c.nodes.outflow[n2 as usize] -= q;
                }
            }
        }

        // (a deferred native hot start would be applied here; see the header note)

        self.ctx_.state = EngineState::Initialized;
        Ok(SWMM_OK)
    }

    /// SWMMEngine::init_modules: the init legs in upstream's order.
    pub(crate) fn init_modules(&mut self) -> Result<(), Terminated> {
        self.init_hydraulics()?;
        self.init_hydrology();
        self.init_quality()?;
        if !self.hotstart_runoff_.is_empty() {
            self.apply_legacy_runoff_section();
        }
        if self.hotstart_species_.n_landuses >= 0 {
            // after msxsurf::initState: the hot start replaces its store and initial buildup.
            let tr = std::mem::take(&mut self.hotstart_species_);
            if hsm::apply_species_surface(&tr, &mut self.ctx_, &mut |c, m| c.warnings.push(m)) {
                for m in 0..self.ctx_.reactions.surface.n_species {
                    let v = crate::transport::msx_surface::remaining_buildup(&self.ctx_, m);
                    self.ctx_.reactions.surface.led_init_buildup[m as usize] = v;
                }
            }
            self.hotstart_species_ = Default::default();
        }
        self.init_geometry();
        #[cfg(feature = "twod")]
        self.init_modules_channel_ground()?;
        self.init_mass_balance();

        // Allocate forcing arrays to match object counts
        let (nn, nl, ns, ng, np) = (self.ctx_.n_nodes(), self.ctx_.n_links(), self.ctx_.n_subcatches(), self.ctx_.n_gages(), self.ctx_.n_pollutants());
        self.ctx_.forcing.resize(nn, nl, ns, ng, np);
        Ok(())
    }

    /// SWMMEngine::reopenMassBalance (the deferred native hot start calls it; see above).
    #[allow(dead_code)]
    pub(crate) fn reopen_mass_balance(&mut self) {
        self.init_mass_balance();
    }

    /// SWMMEngine::applyLegacyRunoffSection: a legacy hot start's runoff section, once the
    /// hydrology and quality state it restores exists (legacy readRunoff).
    pub(crate) fn apply_legacy_runoff_section(&mut self) {
        let nsc = self.ctx_.n_subcatches();
        let np = self.ctx_.n_pollutants();
        let nlu = self.ctx_.n_landuses();
        if self.hotstart_runoff_.len() != hsm::legacy_runoff_size(&self.ctx_, nlu) {
            self.hotstart_runoff_.clear();
            return;
        }
        let xs = std::mem::take(&mut self.hotstart_runoff_);
        let mut x = 0usize;
        for i in 0..nsc {
            let ui = i as usize;
            {
                let soa = &mut self.runoff_.soa_;
                let sub = &mut self.ctx_.subcatches;
                soa.depth_imperv0[ui] = xs[x];
                x += 1;
                soa.depth_imperv1[ui] = xs[x];
                x += 1;
                soa.depth_perv[ui] = xs[x];
                x += 1;
                sub.runoff[ui] = xs[x];
                x += 1;
                sub.ponded_depth[ui] = soa.depth_imperv0[ui] * soa.frac_imperv0[ui] + soa.depth_imperv1[ui] * soa.frac_imperv1[ui] + soa.depth_perv[ui] * (1.0 - soa.imperv_pct[ui]);
            }
            let mut model = 0i32;
            let mut v6 = [0.0f64; 6];
            self.runoff_.infil_get_state(i, &mut model, &mut v6);
            infil_from_legacy(model, &xs[x..x + 6], &mut v6);
            self.runoff_.infil_set_state(i, model, &v6);
            x += 6;
            let sub = &mut self.ctx_.subcatches;
            if ui < sub.gw_aquifer.len() && sub.gw_aquifer[ui] >= 0 {
                // legacy gwater_setState
                let gw = &mut self.groundwater_.soa_;
                gw.theta[ui] = xs[x];
                gw.lower_depth[ui] = xs[x + 1] - gw.bottom_elev[ui];
                gw.old_flow[ui] = xs[x + 2];
                if xs[x + 3] != -1.0e10 {
                    gw.max_infil_vol[ui] = xs[x + 3];
                }
                x += 4;
            }
            if ui < sub.snowpack.len() && sub.snowpack[ui] >= 0 {
                let sn = &mut self.snow_.soa_;
                for j in 0..crate::hydrology::snow::N_SUBAREAS {
                    let pk = ui * crate::hydrology::snow::N_SUBAREAS as usize + j as usize;
                    if pk < sn.wsnow.len() {
                        sn.wsnow[pk] = xs[x];
                        sn.fw[pk] = xs[x + 1];
                        sn.coldc[pk] = xs[x + 2];
                        sn.ati[pk] = xs[x + 3];
                        sn.awe[pk] = xs[x + 4];
                    }
                    x += 5;
                }
            }
            if np > 0 {
                let unp = np as usize;
                for p in 0..unp {
                    sub.conc_old[ui * unp + p] = xs[x];
                    sub.conc[ui * unp + p] = xs[x];
                    x += 1;
                }
                for p in 0..unp {
                    sub.ponded_qual[ui * unp + p] = xs[x];
                    x += 1;
                }
                for lu in 0..nlu {
                    for p in 0..np {
                        let bu = self.surface_quality_.bu_idx(i, lu, p);
                        let norm = buildup_norm(&self.ctx_, &self.landuse_solver_, ui, lu, p);
                        if bu < self.surface_quality_.buildup.len() {
                            self.surface_quality_.buildup[bu] = if norm > 0.0 { xs[x] / norm } else { 0.0 };
                        }
                        x += 1;
                    }
                    let sub = &mut self.ctx_.subcatches;
                    let sw = ui * nlu as usize + lu as usize;
                    if sw < sub.sweep_last_swept.len() {
                        let n = sub.sweep_last_swept.len();
                        sub.sweep_last_date.resize(n, self.ctx_.options.start_date);
                        sub.sweep_last_date[sw] = xs[x];
                        sub.sweep_last_swept[sw] = self.ctx_.options.start_date - xs[x];
                    }
                    x += 1;
                }
            }
        }
        // legacy massbal_open runs after readRunoff: the initial load is what the hot
        // start left on the land, ponded mass included.
        if np > 0 {
            self.compute_final_quality_mass_balance();
            self.ctx_.mass_balance.qual_init_buildup = self.ctx_.mass_balance.qual_final_buildup.clone();
        }
        self.hotstart_runoff_.clear();
    }

    /// SWMMEngine::computeFinalQualityMassBalance: the buildup (and ponded mass, and ponds
    /// written off) still on the land, per pollutant (legacy massbal_getBuildup).
    pub(crate) fn compute_final_quality_mass_balance(&mut self) {
        if self.ctx_.n_pollutants() > 0 {
            let np = self.ctx_.n_pollutants();
            let nlu = std::cmp::max(self.ctx_.n_landuses(), 1);
            let mcf = pollutant_surface_mcf(&self.ctx_);
            for p in 0..np {
                let mut total = 0.0;
                for i in 0..self.ctx_.n_subcatches() {
                    let ui = i as usize;
                    let area_ac = self.ctx_.subcatches.area[ui];
                    for lu in 0..nlu {
                        let cov_idx = ui * nlu as usize + lu as usize;
                        let frac = if cov_idx < self.ctx_.subcatches.coverage.len() { self.ctx_.subcatches.coverage[cov_idx] / 100.0 } else { 0.0 };
                        if frac <= 0.0 {
                            continue;
                        }
                        let k = (lu * np + p) as usize;
                        let bp = &self.landuse_solver_.buildup_params[k];
                        let norm = if bp.normalizer == 0 { frac * area_ac } else { frac * self.ctx_.subcatches.curb_length[ui] };
                        let bu = self.surface_quality_.bu_idx(i, lu, p);
                        if bu < self.surface_quality_.buildup.len() {
                            total += self.surface_quality_.buildup[bu] * norm;
                        }
                    }
                    let sq = ui * np as usize + p as usize;
                    if sq < self.ctx_.subcatches.ponded_qual.len() {
                        total += self.ctx_.subcatches.ponded_qual[sq] * mcf[p as usize];
                    }
                }
                if (p as usize) < self.ctx_.mass_balance.qual_ponded_lost.len() {
                    total += self.ctx_.mass_balance.qual_ponded_lost[p as usize];
                }
                self.ctx_.mass_balance.qual_final_buildup[p as usize] = total;
            }
        }
    }
}
