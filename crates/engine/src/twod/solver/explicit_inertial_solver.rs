// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/solver/ExplicitInertialSolver.hpp and
// ExplicitInertialSolver.cpp (initialize, refreshSubsurfaceRows, reconstructAll,
// pinExternalSourceCells, publishMomentumVelocity) (HydroCouple/openswmm.engine, branch
// swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the Apache
// License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The explicit local-inertial LTS marcher (the serial CPU surface solver): its members, in
// upstream's order with upstream's names, and what SurfaceRouter2D::initialize runs of it. The
// step (advance and everything it calls) is wave 4.
//
// Differences in form, never in result:
//   * Upstream's non-owning pointers are flags here (true = attached): mesh_, state_ and
//     opts_ name the router's mesh_, state_ and options_ (the translation passes them
//     explicitly to each call); rk_reactions_, channel_ground_, common_stage_,
//     source_history_, engine_source_history_, engine_source_stream_ and gw_ belong to paths
//     phase 2a refuses (coupled-groundwater RK2, MESH_CONTACTS channel heat, the common source
//     stage, the source-result histories, the subsurface kernel) and stay false.
//   * rk_reaction_workspace_ (transport::RxWorkspace) is not modelled: only the refused
//     coupled-groundwater RK2 reaction path uses it, and initialize leaves it untouched.
//   * Rk2BoundarySignal carries no fields here: initialize only clears the list, and only the
//     refused RK2 path fills it.
//   * The exceptions initialize throws (std::invalid_argument, std::runtime_error, caught by
//     SWMMEngine::initHydraulics) are Err(InitError::Thrown(what)); the
//     OPENSWMM_2D_MARCHER_TELEMETRY development switch, read at the end of initialize, is
//     refused with the crate's stub when set (PLAN.md, the environment-switch policy).
//   * The fprintf(stderr, ...) notices are eprint! with the same text (core::cprintf formats
//     the numbers as glibc does).
//   * Upstream's OpenMP loops here are `schedule(static)` with each iteration writing only its
//     own cell, so the serial loops give the same bits at any thread count.

#![allow(clippy::all, non_snake_case)]

use crate::cformat;
use crate::core::cmath::cclamp;
use crate::core::engine_step::stub;
use crate::core::options::Bytes;
use crate::twod::coupling::node_coupling::CouplingPoint;
use crate::twod::data::boundary_data::{BoundaryData, BoundaryType};
use crate::twod::data::mesh_data::MeshData;
use crate::twod::data::solver_options_2d::{turbulence_name, Coriolis2D, Momentum2D, SolverOptions2D};
use crate::twod::data::surface_state_data::SurfaceStateData;
use crate::twod::mesh::vertex_reconstruction::build_vertex_stencils;
use crate::twod::solver::inertial_edges::InertialEdges;
use crate::twod::solver::inertial_kernels::cell_eta_depth;
use crate::twod::solver::momentum_physics_2d::MomentumPhysics2D;
use crate::twod::surface_router_2d::InitError;

/// Active-set / tier rebuild cadence in macro cycles (kRebuildEveryCycles).
pub const K_REBUILD_EVERY_CYCLES: i32 = 4;

/// ExplicitInertialSolver::Rk2BoundarySignal (no fields here; see the header).
#[derive(Debug, Clone, Default)]
pub struct Rk2BoundarySignal;

/// ExplicitInertialSolver::CellSourceStep.
#[derive(Debug, Clone, Default)]
pub struct CellSourceStep {
    pub area_dt: f64,
    pub rain: f64,
    pub coupling: f64,
    pub available: f64,
    pub infil: f64,
    pub evap: f64,
    pub out: f64,
    pub final_volume: f64,
    pub requested_infil: f64,
    /// `double* ledger` (a pointer; set or not).
    pub ledger: bool,
}

/// ExplicitInertialSolver::SourceLedger.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SourceLedger {
    Infiltration,
    CouplingOut,
    Rainfall,
    CouplingIn,
    Exfiltration,
    Boundary,
    ExternalIn,
    ExternalOut,
    BedExchange,
    Count,
}

#[derive(Debug, Clone)]
pub struct ExplicitInertialSolver {
    pub rk_reactions_: bool,
    pub rk_reaction_rows_: Vec<i32>,
    pub rk_decay_: Vec<f64>,
    pub rk_reacted_: Vec<f64>,
    pub rk_reacted_reported_: Vec<f64>,
    pub rk_reaction_identity_: Bytes,
    pub mesh_: bool,
    pub state_: bool,
    pub opts_: bool,
    pub edges_: InertialEdges,
    pub q_: Vec<f64>,
    pub qcx_: Vec<f64>,
    pub qcy_: Vec<f64>,
    pub mode_: Momentum2D,
    pub has_quads_: bool,
    pub wall_ptr_: Vec<i32>,
    pub wall_slot_: Vec<i32>,
    pub swe_bc_dqx_: f64,
    pub swe_bc_dqy_: f64,
    pub second_order_: bool,
    pub rk_transport_stage_: bool,
    pub rk_stage_mass_: Vec<f64>,
    pub rk_stage_volume_: Vec<f64>,
    pub gex_: Vec<f64>,
    pub gey_: Vec<f64>,
    pub gux_: Vec<f64>,
    pub guy_: Vec<f64>,
    pub gvx_: Vec<f64>,
    pub gvy_: Vec<f64>,
    pub rk_v0_: Vec<f64>,
    pub rk_qx0_: Vec<f64>,
    pub rk_qy0_: Vec<f64>,
    pub rk_inf0_: Vec<f64>,
    pub rk_bc0_: Vec<f64>,
    pub rk_ex0_: Vec<f64>,
    pub rk_drawn0_: Vec<f64>,
    pub rk_cpl0_: Vec<f64>,
    pub momentum_physics_: MomentumPhysics2D,
    pub wall_drag_slots_: Vec<i32>,
    pub wall_drag_energy_: f64,
    pub rk_boundary_signals_: Vec<Rk2BoundarySignal>,
    pub rk_boundary_anchor_: f64,
    pub front_rebuild_: bool,
    pub frontier_: Vec<u8>,
    pub front_breach_: bool,
    pub h_on_front_: f64,
    pub dw_slope_: Vec<f64>,
    pub cell_active_: Vec<u8>,
    pub active_cells_: Vec<i32>,
    pub active_faces_: Vec<i32>,
    pub pin_t0_: Vec<u8>,
    pub rebuild_seed_: Vec<u8>,
    pub rebuild_dt_cell_: Vec<f64>,
    pub rebuild_dt_partial_: Vec<f64>,
    pub tier_: Vec<u8>,
    pub face_tier_: Vec<u8>,
    pub face_accum_: Vec<f64>,
    pub face_accum_stride_: usize,
    pub sacc_L_: Vec<f64>,
    pub sacc_R_: Vec<f64>,
    pub species_on_: bool,
    pub plain_sources_: bool,
    pub source_ledgers_: Vec<f64>,
    pub source_ledger_stride_: usize,
    pub channel_ground_: bool,
    pub channel_sources_: Vec<CellSourceStep>,
    pub channel_source_ready_: Vec<u8>,
    pub channel_temperature_mixed_: Vec<u8>,
    pub accumulators_pending_: bool,
    pub cells_by_tier_: Vec<Vec<i32>>,
    pub edges_by_tier_: Vec<Vec<i32>>,
    pub dt0_: f64,
    pub bc_cell_: Vec<i32>,
    pub bc_slot_: Vec<i32>,
    pub gw_row_map_: Vec<i32>,
    pub bc_accum_: Vec<f64>,
    pub bc_q_: Vec<f64>,
    pub exch_: Vec<f64>,
    pub node_drawn_: Vec<f64>,
    pub spill_deficit_: f64,
    pub exch_head_slope_: Vec<f64>,
    pub exch_tau_: f64,
    pub t_last_sync_: f64,
    pub common_stage_: bool,
    pub source_history_: bool,
    pub engine_source_history_: bool,
    pub engine_source_stream_: bool,
    pub checkpoint_ready_: bool,
    pub joint_checkpoint_: bool,
    pub source_phase_end_: f64,
    pub cycles_since_rebuild_: i32,
    pub substeps_run_: i64,
    pub face_passes_: i64,
    pub last_steps_: i64,
    pub last_dt_: f64,
    pub initialized_: bool,
    pub telemetry_: Vec<(f64, i32)>,
    pub telemetry_path_: Bytes,
    pub active_frac_min_: f64,
    pub active_frac_max_: f64,
    pub active_frac_sum_: f64,
    pub active_samples_: i64,
    pub tier_occupancy_: [i64; 8],
    pub gw_: bool,
    pub gw_end_time_: f64,
    pub gw_quiescent_: bool,
}

impl Default for ExplicitInertialSolver {
    fn default() -> Self {
        ExplicitInertialSolver {
            rk_reactions_: false,
            rk_reaction_rows_: Vec::new(),
            rk_decay_: Vec::new(),
            rk_reacted_: Vec::new(),
            rk_reacted_reported_: Vec::new(),
            rk_reaction_identity_: Bytes::new(),
            mesh_: false,
            state_: false,
            opts_: false,
            edges_: InertialEdges::default(),
            q_: Vec::new(),
            qcx_: Vec::new(),
            qcy_: Vec::new(),
            mode_: Momentum2D::LocalInertial,
            has_quads_: false,
            wall_ptr_: Vec::new(),
            wall_slot_: Vec::new(),
            swe_bc_dqx_: 0.0,
            swe_bc_dqy_: 0.0,
            second_order_: false,
            rk_transport_stage_: false,
            rk_stage_mass_: Vec::new(),
            rk_stage_volume_: Vec::new(),
            gex_: Vec::new(),
            gey_: Vec::new(),
            gux_: Vec::new(),
            guy_: Vec::new(),
            gvx_: Vec::new(),
            gvy_: Vec::new(),
            rk_v0_: Vec::new(),
            rk_qx0_: Vec::new(),
            rk_qy0_: Vec::new(),
            rk_inf0_: Vec::new(),
            rk_bc0_: Vec::new(),
            rk_ex0_: Vec::new(),
            rk_drawn0_: Vec::new(),
            rk_cpl0_: Vec::new(),
            momentum_physics_: MomentumPhysics2D::default(),
            wall_drag_slots_: Vec::new(),
            wall_drag_energy_: 0.0,
            rk_boundary_signals_: Vec::new(),
            rk_boundary_anchor_: 0.0,
            front_rebuild_: false,
            frontier_: Vec::new(),
            front_breach_: false,
            h_on_front_: 0.0,
            dw_slope_: Vec::new(),
            cell_active_: Vec::new(),
            active_cells_: Vec::new(),
            active_faces_: Vec::new(),
            pin_t0_: Vec::new(),
            rebuild_seed_: Vec::new(),
            rebuild_dt_cell_: Vec::new(),
            rebuild_dt_partial_: Vec::new(),
            tier_: Vec::new(),
            face_tier_: Vec::new(),
            face_accum_: Vec::new(),
            face_accum_stride_: 2,
            sacc_L_: Vec::new(),
            sacc_R_: Vec::new(),
            species_on_: false,
            plain_sources_: false,
            source_ledgers_: Vec::new(),
            source_ledger_stride_: 0,
            channel_ground_: false,
            channel_sources_: Vec::new(),
            channel_source_ready_: Vec::new(),
            channel_temperature_mixed_: Vec::new(),
            accumulators_pending_: false,
            cells_by_tier_: Vec::new(),
            edges_by_tier_: Vec::new(),
            dt0_: 0.0,
            bc_cell_: Vec::new(),
            bc_slot_: Vec::new(),
            gw_row_map_: Vec::new(),
            bc_accum_: Vec::new(),
            bc_q_: Vec::new(),
            exch_: Vec::new(),
            node_drawn_: Vec::new(),
            spill_deficit_: 0.0,
            exch_head_slope_: Vec::new(),
            exch_tau_: 0.0,
            t_last_sync_: 0.0,
            common_stage_: false,
            source_history_: false,
            engine_source_history_: false,
            engine_source_stream_: false,
            checkpoint_ready_: false,
            joint_checkpoint_: false,
            source_phase_end_: 0.0,
            cycles_since_rebuild_: 1000,
            substeps_run_: 0,
            face_passes_: 0,
            last_steps_: 0,
            last_dt_: 0.0,
            initialized_: false,
            telemetry_: Vec::new(),
            telemetry_path_: Bytes::new(),
            active_frac_min_: 1.0e30,
            active_frac_max_: -1.0e30,
            active_frac_sum_: 0.0,
            active_samples_: 0,
            tier_occupancy_: [0; 8],
            gw_: false,
            gw_end_time_: -1.0,
            gw_quiescent_: false,
        }
    }
}

fn assign<T: Clone>(v: &mut Vec<T>, n: usize, x: T) {
    v.clear();
    v.resize(n, x);
}

fn invalid(what: &str) -> InitError {
    InitError::Thrown(what.as_bytes().to_vec())
}

/// What initialize reads of the surface state's attachments (state.boundary,
/// state.node_coupling, state.nodes_1d), passed explicitly.
pub struct SolverAttachments<'a> {
    /// *state.boundary (the router's boundary_), when attached.
    pub boundary: Option<&'a BoundaryData>,
    /// *state.node_coupling (the router's node_coupling_points_), when attached.
    pub node_coupling: Option<&'a [CouplingPoint]>,
    /// state.nodes_1d->volume.size(), when nodes_1d is attached.
    pub nodes_1d_len: Option<usize>,
}

impl ExplicitInertialSolver {
    pub fn is_initialized(&self) -> bool {
        self.initialized_
    }

    pub fn inertial_edges(&self) -> &InertialEdges {
        &self.edges_
    }

    /// initialize(mesh, state, opts).
    pub fn initialize(&mut self, mesh: &mut MeshData, state: &mut SurfaceStateData, att: &SolverAttachments, opts: &mut SolverOptions2D) -> Result<(), InitError> {
        self.wall_drag_slots_.clear();
        self.wall_drag_energy_ = 0.0;
        if let Some(b) = att.boundary {
            for slot in 0..b.edge_wall_drag.len() as i32 {
                let cw = b.edge_wall_drag[slot as usize];
                if !cw.is_finite() || cw < 0.0 {
                    return Err(invalid("WALL_DRAG requires finite nonnegative Cw."));
                }
                if cw == 0.0 {
                    continue;
                }
                let i = MeshData::slot_cell(slot);
                let e = MeshData::slot_local(slot);
                if i >= mesh.n_cells() || e >= mesh.cell_vertex_count(i) || mesh.cell_neighbour(i, e) >= 0 || b.edge_bc_type[slot as usize] != BoundaryType::Wall as i8 {
                    return Err(invalid("WALL_DRAG must name an exterior WALL edge."));
                }
                self.wall_drag_slots_.push(slot);
            }
        }
        opts.wall_drag_enabled = !self.wall_drag_slots_.is_empty();
        if !opts.eddy_viscosity_background.is_finite() || opts.eddy_viscosity_background < 0.0 {
            return Err(invalid("EDDY_VISCOSITY_BACKGROUND must be finite and nonnegative."));
        }
        if opts.momentum_physics() {
            if opts.momentum != Momentum2D::FullSwe || opts.lts_tiers != 1 {
                return Err(invalid("Coriolis/Smagorinsky/wall drag require FULL_SWE and LTS_TIERS 1."));
            }
            if opts.coriolis == Coriolis2D::Latitude && !opts.coriolis_latitude_set {
                return Err(invalid("CORIOLIS LATITUDE requires CORIOLIS_LATITUDE."));
            }
            if opts.coriolis == Coriolis2D::ConstantF && !opts.coriolis_f_set {
                return Err(invalid("CORIOLIS CONSTANT_F requires CORIOLIS_F."));
            }
            if !opts.coriolis_parameter().is_finite()
                || !opts.smagorinsky_cs.is_finite()
                || !opts.mixing_length_coefficient.is_finite()
                || opts.mixing_length_coefficient < 0.0
                || opts.smagorinsky_cs < 0.0
                || !opts.coriolis_latitude.is_finite()
                || opts.coriolis_latitude.abs() > 90.0
            {
                return Err(invalid("Invalid Coriolis/Smagorinsky parameters."));
            }
            if opts.groundwater == 1 || state.transport.active() || state.heat.enabled || state.snow.active() || att.node_coupling.is_some_and(|nc| !nc.is_empty()) {
                return Err(invalid("Coriolis/Smagorinsky/wall drag currently require uncoupled surface hydraulics."));
            }
            let m = cformat!(
                b"[openswmm 2D] optional momentum physics: f=%.17g s^-1, turbulence=%s Cs=%.17g, mixing-length c=%.17g, background nu=%.17g m2/s, CPU global stepping.\n",
                opts.coriolis_parameter(),
                turbulence_name(opts.turbulence),
                opts.smagorinsky_cs,
                opts.mixing_length_coefficient,
                opts.eddy_viscosity_background
            );
            eprint_bytes(&m);
        }
        self.momentum_physics_ = MomentumPhysics2D::default();
        self.mesh_ = true;
        self.state_ = true;
        self.opts_ = true;
        self.rk_boundary_signals_.clear();
        self.rk_boundary_anchor_ = 0.0;
        self.common_stage_ = false;
        self.rk_reactions_ = false;
        self.rk_reaction_rows_.clear();
        self.rk_decay_.clear();
        self.rk_reacted_.clear();
        self.rk_reacted_reported_.clear();
        self.rk_reaction_identity_.clear();
        self.channel_ground_ = false;
        self.channel_sources_.clear();
        self.channel_source_ready_.clear();
        self.channel_temperature_mixed_.clear();
        self.source_history_ = false;
        self.engine_source_history_ = false;
        self.engine_source_stream_ = false;
        self.checkpoint_ready_ = true;
        self.source_phase_end_ = 0.0;

        let nt = mesh.n_triangles();
        if nt <= 0 {
            return Ok(());
        }
        let unt = nt as usize;

        self.edges_.build(mesh);
        let ne = self.edges_.ne as usize;
        assign(&mut self.q_, ne, 0.0);
        self.face_accum_stride_ = if opts.momentum == Momentum2D::FullSwe { 6 } else { 2 };
        assign(&mut self.face_accum_, ne * self.face_accum_stride_, 0.0);
        assign(&mut self.face_tier_, ne, 0u8);

        // Momentum closure (MOMENTUM_EQUATION).
        self.mode_ = opts.momentum;
        self.second_order_ = false;

        self.dw_slope_.clear();
        if self.mode_ == Momentum2D::FullSwe {
            if opts.cfl_number > 0.5 {
                eprint_bytes(&cformat!(b"[openswmm 2D] MOMENTUM_EQUATION FULL_SWE: CFL_NUMBER %.3g reduced to 0.5 (positivity bound of the Godunov update).\n", opts.cfl_number));
                opts.cfl_number = 0.5;
            }
            if opts.reconstruction_order == 2 && opts.lts_tiers > 1 {
                eprint_bytes(&cformat!(b"[openswmm 2D] RECONSTRUCTION_ORDER 2 runs in global-dt mode: LTS_TIERS %d reduced to 1.\n", opts.lts_tiers));
                opts.lts_tiers = 1;
            }
            self.second_order_ = opts.reconstruction_order == 2;
            if self.second_order_ && state.transport.active() && opts.groundwater == 0 {
                eprint_bytes(b"[openswmm 2D] RECONSTRUCTION_ORDER 2 with overland transport is not supported yet; using first order.\n");
                self.second_order_ = false;
            }
            if self.second_order_ {
                if mesh.vert_stencil_ptr.len() != (mesh.n_vertices() + 1) as usize {
                    build_vertex_stencils(mesh);
                }
                assign(&mut self.gex_, 2 * unt, 0.0);
                assign(&mut self.gey_, 2 * unt, 0.0);
                assign(&mut self.gux_, unt, 0.0);
                assign(&mut self.guy_, unt, 0.0);
                assign(&mut self.gvx_, unt, 0.0);
                assign(&mut self.gvy_, unt, 0.0);
            } else {
                self.gex_.clear();
                self.gey_.clear();
                self.gux_.clear();
                self.guy_.clear();
                self.gvx_.clear();
                self.gvy_.clear();
            }
        } else if self.mode_ == Momentum2D::DiffusiveWave {
            assign(&mut self.dw_slope_, unt, 0.0);
        }
        // FRONT_REBUILD AUTO: on for FULL_SWE only.
        self.front_rebuild_ = if opts.front_rebuild < 0 { self.mode_ == Momentum2D::FullSwe } else { opts.front_rebuild != 0 };
        assign(&mut self.frontier_, if self.front_rebuild_ { unt } else { 0 }, 0u8);
        self.front_breach_ = false;
        self.has_quads_ = mesh.n_quads() > 0;
        // FULL_SWE: WALL boundary slots per cell, CSR.
        self.wall_ptr_.clear();
        self.wall_slot_.clear();
        if self.mode_ == Momentum2D::FullSwe {
            self.wall_ptr_.reserve(unt + 1);
            self.wall_ptr_.push(0);
            for i in 0..nt {
                let nvc = mesh.cell_vertex_count(i);
                for k2 in 0..nvc {
                    if mesh.cell_neighbour(i, k2) >= 0 {
                        continue;
                    }
                    let slot = MeshData::slot(i, k2);
                    if let Some(b) = att.boundary {
                        if b.edge_bc_type[slot as usize] != BoundaryType::Wall as i8 {
                            continue;
                        }
                    }
                    self.wall_slot_.push(slot);
                }
                self.wall_ptr_.push(self.wall_slot_.len() as i32);
            }
        }
        // S1: species accumulators, sized only when transport is live.
        self.species_on_ = state.transport.active();
        if self.second_order_ && opts.groundwater != 0 && state.transport.age_row >= 0 {
            let ns = state.transport.n_species as usize;
            assign(&mut state.transport.generated_age, ns, 0.0);
        }
        if self.species_on_ {
            let ns = state.transport.n_species as usize;
            assign(&mut self.sacc_L_, ns * ne, 0.0);
            assign(&mut self.sacc_R_, ns * ne, 0.0);
        } else {
            self.sacc_L_.clear();
            self.sacc_R_.clear();
        }
        let nledger = SourceLedger::Count as usize * state.transport.n_species as usize;
        self.source_ledger_stride_ = ((nledger + 1 + 7) / 8) * 8;
        assign(&mut self.source_ledgers_, self.source_ledger_stride_ * std::cmp::max(1, opts.num_threads) as usize, 0.0);
        // The Perot cell vectors (theta-blend, convective term, FULL_SWE cell momentum).
        self.qcx_.clear();
        self.qcy_.clear();
        if opts.theta < 1.0 || opts.advection || self.mode_ == Momentum2D::FullSwe {
            assign(&mut self.qcx_, unt, 0.0);
            assign(&mut self.qcy_, unt, 0.0);
        }
        assign(&mut self.cell_active_, unt, 0u8);
        assign(&mut self.tier_, unt, 0u8);

        let k = cclamp(opts.lts_tiers, 1, 8);
        assign(&mut self.cells_by_tier_, k as usize, Vec::new());
        assign(&mut self.edges_by_tier_, k as usize, Vec::new());

        // Non-WALL boundary entries, evaluated at their owning cell's firing.
        self.bc_cell_.clear();
        self.bc_slot_.clear();
        if let Some(b) = att.boundary {
            for i in 0..nt {
                let nvc = mesh.cell_vertex_count(i);
                for e in 0..nvc {
                    let idx = MeshData::slot(i, e);
                    let ty = b.edge_bc_type[idx as usize];
                    let interior = mesh.cell_neighbour(i, e) >= 0;
                    if !interior && ty != BoundaryType::Wall as i8 {
                        self.bc_cell_.push(i);
                        self.bc_slot_.push(idx);
                    }
                }
            }
        }
        assign(&mut self.bc_accum_, self.bc_cell_.len(), 0.0);
        assign(&mut self.bc_q_, self.bc_cell_.len(), 0.0);

        // S2: [2D_BOUNDARY_QUALITY] rows onto this solver's boundary-edge list.
        if state.transport.active() {
            let ns = state.transport.n_species as usize;
            assign(&mut state.transport.bc_conc, self.bc_cell_.len() * ns, 0.0);
            let rows = state.transport.bc_quality_rows.clone();
            for row in &rows {
                let mut kk = 0usize;
                let mut found = false;
                while kk < self.bc_slot_.len() {
                    if self.bc_slot_[kk] == row.slot {
                        found = true;
                        break;
                    }
                    kk += 1;
                }
                if !found {
                    return Err(InitError::Thrown(
                        format!("[2D_BOUNDARY_QUALITY] edge slot {} is not a non-WALL boundary edge, so its concentration could never enter the mesh.", row.slot).into_bytes(),
                    ));
                }
                if row.species < 0 || row.species as usize >= ns {
                    return Err(invalid("[2D_BOUNDARY_QUALITY] species index out of range."));
                }
                state.transport.bc_conc[kk * ns + row.species as usize] = row.conc;
            }
        }

        // Live junction exchange: one integral accumulator per point; spill budget per node.
        assign(&mut self.exch_, att.node_coupling.map_or(0, |nc| nc.len()), 0.0);
        assign(&mut self.node_drawn_, att.nodes_1d_len.unwrap_or(0), 0.0);
        assign(&mut self.pin_t0_, unt, 0u8);
        for &i in &self.bc_cell_ {
            self.pin_t0_[i as usize] = 1;
        }
        if let Some(nc) = att.node_coupling {
            for cp in nc {
                if cp.cell_idx >= 0 {
                    self.pin_t0_[cp.cell_idx as usize] = 1;
                }
            }
        }

        self.reconstruct_all(mesh, state, opts);

        // Seed face momentum from the optional [2D_INITIAL_VELOCITY] rows.
        {
            let mut any_uv = false;
            let mut i = 0usize;
            while i < unt && !any_uv {
                any_uv = mesh.tri_init_u[i] != 0.0 || mesh.tri_init_v[i] != 0.0;
                i += 1;
            }
            if any_uv && self.mode_ == Momentum2D::FullSwe {
                // Cell momentum seeds directly: q = h u.
                for i in 0..unt {
                    self.qcx_[i] = state.depth[i] * mesh.tri_init_u[i];
                    self.qcy_[i] = state.depth[i] * mesh.tri_init_v[i];
                }
            } else if any_uv {
                for e in 0..ne {
                    let a = self.edges_.cL[e] as usize;
                    let b = self.edges_.cR[e] as usize;
                    let qax = state.depth[a] * mesh.tri_init_u[a];
                    let qay = state.depth[a] * mesh.tri_init_v[a];
                    let qbx = state.depth[b] * mesh.tri_init_u[b];
                    let qby = state.depth[b] * mesh.tri_init_v[b];
                    self.q_[e] = if self.edges_.conveyance[e] > 0.0 { 0.5 * ((qax + qbx) * self.edges_.nx[e] + (qay + qby) * self.edges_.ny[e]) } else { 0.0 };
                }
                // Perot gather so the theta-blend and CFL speed see the seeded momentum.
                if !self.qcx_.is_empty() {
                    let ed = &self.edges_;
                    for i in 0..unt {
                        let mut sx = 0.0;
                        let mut sy = 0.0;
                        for p in ed.cell_ptr[i] as usize..ed.cell_ptr[i + 1] as usize {
                            let e = ed.cell_edge[p] as usize;
                            let mut fq = ed.cell_sign[p] as f64 * self.q_[e] * ed.xi[e];
                            if self.mode_ == Momentum2D::LocalInertial {
                                fq *= ed.conveyance[e];
                            }
                            sx += fq * (ed.mx[e] - mesh.tri_cx[i]);
                            sy += fq * (ed.my[e] - mesh.tri_cy[i]);
                        }
                        let inv_a = 1.0 / mesh.tri_area[i];
                        self.qcx_[i] = sx * inv_a;
                        self.qcy_[i] = sy * inv_a;
                    }
                }
            }
        }
        self.t_last_sync_ = 0.0;
        self.last_steps_ = 0;
        self.face_passes_ = 0;
        self.substeps_run_ = 0;
        self.last_dt_ = 0.0;
        self.telemetry_.clear();
        self.active_frac_min_ = 1.0e30;
        self.active_frac_max_ = -1.0e30;
        self.active_frac_sum_ = 0.0;
        self.active_samples_ = 0;
        if std::env::var_os("OPENSWMM_2D_MARCHER_TELEMETRY").is_some() {
            return Err(InitError::Stub(stub("the OPENSWMM_2D_MARCHER_TELEMETRY development switch (ExplicitInertialSolver::initialize)")));
        }
        self.refresh_subsurface_rows(state);
        if let Some(b) = att.boundary {
            for &slot in &self.bc_slot_ {
                b.ensure_receipt(slot, state.transport.n_species);
            }
        }
        self.initialized_ = true;
        Ok(())
    }

    /// refreshSubsurfaceRows(): the surface -> groundwater row map (no subsurface kernel in
    /// phase 2a: gw_ is never attached, so every row maps to -1).
    pub fn refresh_subsurface_rows(&mut self, state: &SurfaceStateData) {
        assign(&mut self.gw_row_map_, state.transport.n_species as usize, -1);
    }

    /// reconstructAll(): every cell's head and depth from its volume.
    pub fn reconstruct_all(&mut self, mesh: &MeshData, state: &mut SurfaceStateData, opts: &SolverOptions2D) {
        let nt = mesh.n_triangles();
        for i in 0..nt {
            let ui = i as usize;
            let (eta, depth) = cell_eta_depth(mesh, opts, i, state.volume[ui]);
            state.head[ui] = eta;
            state.depth[ui] = depth;
        }
    }

    /// pinExternalSourceCells(cells): pin the cells to tier 0 and force a rebuild.
    pub fn pin_external_source_cells(&mut self, cells: &[i32]) {
        for &cell in cells {
            self.pin_t0_[cell as usize] = 1;
        }
        self.cycles_since_rebuild_ = K_REBUILD_EVERY_CYCLES;
    }

    /// publishMomentumVelocity(): FULL_SWE's output cell velocity from the cell momentum.
    pub fn publish_momentum_velocity(&self, mesh: &MeshData, state: &mut SurfaceStateData, opts: &SolverOptions2D) {
        if !self.initialized_ || self.mode_ != Momentum2D::FullSwe {
            return;
        }
        for i in 0..mesh.n_cells() as usize {
            let h = state.depth[i];
            state.face_vx[i] = if h > opts.dry_depth { -self.qcx_[i] / h } else { 0.0 };
            state.face_vy[i] = if h > opts.dry_depth { -self.qcy_[i] / h } else { 0.0 };
        }
    }
}

/// fprintf(stderr, ...) of already formatted bytes.
pub fn eprint_bytes(b: &[u8]) {
    use std::io::Write;
    let _ = std::io::stderr().lock().write_all(b);
}
