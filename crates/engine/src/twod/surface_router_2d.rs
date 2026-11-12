// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/SurfaceRouter2D.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The 2D surface router: every data member, in upstream's order and with upstream's names
// (wave 1 carried the members open() fills; wave 3 added the rest), and the accessors.
// initialize() and what it calls are in twod/surface_router_2d_init.rs.
//
// Upstream's SWMMEngine wires SimulationContext::twod_io with raw pointers to some of these
// members (wire2DModelIO), and the section handlers capture references to them. Here each such
// member is a shared handle (twod::Shared, Rc<RefCell<T>>): the router, the context's twod_io
// and the registered handlers hold clones of the same handle.
//
// Differences in form, never in result:
//   * The non-owning pointers joint_hydraulic_, common_stage_ and source_context_ are flags
//     (set or not): the paths that set them (joint hydraulic pairing, the common source stage,
//     the RK2 checkpoint) are phase 2b and leave them false.
//   * solver_ (std::unique_ptr<ISurfaceSolver>) is an Option<ExplicitInertialSolver>: the
//     serial CPU marcher is the only solver this build constructs (the GPU/Kokkos plugin is not
//     built; twod/solver/surface_solver_factory.rs).
//   * mb_seen_nodes_ (std::unordered_set<int>) is a HashSet: initialize never fills it, and the
//     step only inserts and clears it (no iteration).

#![allow(clippy::all)]

use std::collections::HashSet;

use crate::core::engine_init_hydraulics::Terminated;
use crate::core::options::Bytes;
use crate::twod::coupling::node_coupling::CouplingPoint;
use crate::twod::coupling::outfall_exchange::OutfallExchange;
use crate::twod::data::boundary_data::BoundaryData;
use crate::twod::data::mesh_data::MeshData;
use crate::twod::data::pending_rows_2d::*;
use crate::twod::data::solver_options_2d::SolverOptions2D;
use crate::twod::data::surface_state_data::SurfaceStateData;
use crate::twod::gw::gw_transport_data::GwTransportData;
use crate::twod::heat::channel_ground::ChannelGroundState;
use crate::twod::infil::infil_2d::Infil2D;
use crate::twod::mesh::rainfall_interpolator::RainfallInterpolator;
use crate::twod::quality::surface_quality_2d::SurfaceQuality2D;
use crate::twod::solver::explicit_inertial_solver::ExplicitInertialSolver;
use crate::twod::subsurface::subsurface_data::SubsurfaceConfig;
use crate::twod::subsurface::subsurface_solver::SubsurfaceSolver;
use crate::twod::{shared, Shared};

/// An exception SurfaceRouter2D::initialize (or what it calls) lets escape, as the catch site
/// in SWMMEngine::initHydraulics / init_modules sees it: `Thrown(what())` for a
/// std::exception, or `Stub` where the translation refuses a path phase 2a does not cover
/// (the crate's stub, core::engine_step::stub).
#[derive(Debug, Clone)]
pub enum InitError {
    Thrown(Bytes),
    Stub(Terminated),
}

impl From<Terminated> for InitError {
    fn from(t: Terminated) -> Self {
        InitError::Stub(t)
    }
}

/// GwLinkShare (subsurface/SubsurfaceData.hpp): a conduit's length share over one aquifer
/// cell. Only the subsurface kernel (phase 2b) resolves them; the list stays empty in 2a.
#[derive(Debug, Clone, Default)]
pub struct GwLinkShare {
    pub link: i32,
    pub cell: i32,
    pub weight: f64,
}

#[derive(Debug)]
pub struct SurfaceRouter2D {
    pub joint_hydraulic_: bool,
    pub common_stage_: bool,
    pub checkpoint_owned_initialization_: bool,
    pub source_context_: bool,
    pub source_router_model_: Bytes,
    pub source_boundary_ready_: bool,
    pub mesh_: Shared<MeshData>,
    pub state_: SurfaceStateData,
    pub channel_ground_: ChannelGroundState,
    pub options_: Shared<SolverOptions2D>,
    pub boundary_: Shared<BoundaryData>,
    pub thread_warnings_: Vec<Bytes>,
    pub backend_name_: Bytes,
    pub infil_: Shared<Infil2D>,
    pub infil_elapsed_: f64,
    pub infil_cum_applied_: Vec<f64>,
    pub coupling_cum_applied_: Vec<f64>,
    pub cell_subcatch_: Vec<i32>,
    pub subcatch_recharge_: Vec<f64>,
    pub gw_: Shared<GwTransportData>,
    pub surface_quality_: Shared<SurfaceQuality2D>,
    pub aquifer_cfg_: Shared<SubsurfaceConfig>,
    pub aquifer_node_names_: Shared<Vec<Bytes>>,
    pub aquifer_link_names_: Shared<Vec<Bytes>>,
    pub subsurface_: SubsurfaceSolver,
    pub aquifer_link_shares_: Vec<GwLinkShare>,
    pub aquifer_link_ids_: Vec<i32>,
    pub rain_cum_: Vec<f64>,
    pub pending_bc_rows_: Shared<Vec<PendingBoundaryRow>>,
    pub pending_iq_rows_: Shared<Vec<PendingInitialQualityRow>>,
    pub pending_ibq_rows_: Shared<Vec<PendingInitialBedQualityRow>>,
    pub pending_bq_rows_: Shared<Vec<PendingBoundaryQualityRow>>,
    pub pending_edge_conveyance_rows_: Shared<Vec<PendingEdgeConveyanceRow>>,
    pub coupling_points_: Vec<CouplingPoint>,
    pub outfall_exchange_: OutfallExchange,
    pub outfall_state_changed_: bool,
    pub node_coupling_points_: Vec<CouplingPoint>,
    pub coupling_head_cells_: Vec<i32>,
    pub pond_track_ptr_: Vec<i32>,
    pub pond_track_cell_: Vec<i32>,
    pub pond_track_share_: Vec<f64>,
    pub pond_chamber_ft2_: Vec<f64>,
    pub node_pond_area_m2_: Vec<f64>,
    pub node_row_conc_: Vec<f64>,
    pub react_phi_: Vec<f64>,
    pub react_a_: Vec<f64>,
    pub react_dx_: Vec<f64>,
    pub co_refresh_elapsed_: f64,
    pub co_render_elapsed_: f64,
    pub render_dirty_: bool,
    pub render_eta_: Vec<f64>,
    pub bc_nonwall_slots_: Vec<i32>,
    pub bc_nonwall_dirty_: bool,
    pub mb_seen_nodes_: HashSet<i32>,
    pub co_forcing_elapsed_: f64,
    pub co_forcing_first_: bool,
    pub active_: bool,
    pub sim_time_: f64,
    pub pending_dt_: f64,
    pub report_old_depth_: Vec<f64>,
    pub report_old_head_: Vec<f64>,
    pub report_window_old_ms_: f64,
    pub report_window_new_ms_: f64,
    pub report_old_valid_: bool,
    pub ramp_prev_head_: Vec<f64>,
    pub ramp_prev_span_: f64,
    pub prev_boundary_cum_: f64,
    pub prev_boundary_in_: f64,
    pub prev_boundary_out_: f64,
    pub window_outfall_accum_: Vec<f64>,
    pub last_t_: f64,
    pub boundary_names_resolved_: bool,
    pub bc_stage_scale_: f64,
    pub bc_stage_ts_scale_: f64,
    pub bc_flow_scale_: f64,
    pub solver_: Option<ExplicitInertialSolver>,
    pub interp_: RainfallInterpolator,
    pub rain_gage_last_: Vec<f64>,
    pub gage_rate_: Vec<f64>,
    pub gage_cum_: Vec<f64>,
    pub gage_window_depth_: Vec<f64>,
    pub gage_window_rate_: Vec<f64>,
    pub gage_since_: f64,
    pub gage_report_depth_: Vec<f64>,
    pub report_rain_cum_: Vec<f64>,
    pub report_rain_ms_: f64,
}

impl Default for SurfaceRouter2D {
    fn default() -> Self {
        SurfaceRouter2D {
            joint_hydraulic_: false,
            common_stage_: false,
            checkpoint_owned_initialization_: false,
            source_context_: false,
            source_router_model_: Bytes::new(),
            source_boundary_ready_: false,
            mesh_: shared(MeshData::default()),
            state_: SurfaceStateData::default(),
            channel_ground_: ChannelGroundState::default(),
            options_: shared(SolverOptions2D::default()),
            boundary_: shared(BoundaryData::default()),
            thread_warnings_: Vec::new(),
            backend_name_: Bytes::new(),
            infil_: shared(Infil2D::default()),
            infil_elapsed_: 0.0,
            infil_cum_applied_: Vec::new(),
            coupling_cum_applied_: Vec::new(),
            cell_subcatch_: Vec::new(),
            subcatch_recharge_: Vec::new(),
            gw_: shared(GwTransportData::default()),
            surface_quality_: shared(SurfaceQuality2D::default()),
            aquifer_cfg_: shared(SubsurfaceConfig::default()),
            aquifer_node_names_: shared(Vec::new()),
            aquifer_link_names_: shared(Vec::new()),
            subsurface_: SubsurfaceSolver::default(),
            aquifer_link_shares_: Vec::new(),
            aquifer_link_ids_: Vec::new(),
            rain_cum_: Vec::new(),
            pending_bc_rows_: shared(Vec::new()),
            pending_iq_rows_: shared(Vec::new()),
            pending_ibq_rows_: shared(Vec::new()),
            pending_bq_rows_: shared(Vec::new()),
            pending_edge_conveyance_rows_: shared(Vec::new()),
            coupling_points_: Vec::new(),
            outfall_exchange_: OutfallExchange::default(),
            outfall_state_changed_: false,
            node_coupling_points_: Vec::new(),
            coupling_head_cells_: Vec::new(),
            pond_track_ptr_: Vec::new(),
            pond_track_cell_: Vec::new(),
            pond_track_share_: Vec::new(),
            pond_chamber_ft2_: Vec::new(),
            node_pond_area_m2_: Vec::new(),
            node_row_conc_: Vec::new(),
            react_phi_: Vec::new(),
            react_a_: Vec::new(),
            react_dx_: Vec::new(),
            co_refresh_elapsed_: 0.0,
            co_render_elapsed_: 0.0,
            render_dirty_: false,
            render_eta_: Vec::new(),
            bc_nonwall_slots_: Vec::new(),
            bc_nonwall_dirty_: true,
            mb_seen_nodes_: HashSet::new(),
            co_forcing_elapsed_: 0.0,
            co_forcing_first_: true,
            active_: false,
            sim_time_: 0.0,
            pending_dt_: 0.0,
            report_old_depth_: Vec::new(),
            report_old_head_: Vec::new(),
            report_window_old_ms_: 0.0,
            report_window_new_ms_: 0.0,
            report_old_valid_: false,
            ramp_prev_head_: Vec::new(),
            ramp_prev_span_: 0.0,
            prev_boundary_cum_: 0.0,
            prev_boundary_in_: 0.0,
            prev_boundary_out_: 0.0,
            window_outfall_accum_: Vec::new(),
            last_t_: 0.0,
            boundary_names_resolved_: false,
            bc_stage_scale_: 1.0,
            bc_stage_ts_scale_: 1.0,
            bc_flow_scale_: 1.0,
            solver_: None,
            interp_: RainfallInterpolator::default(),
            rain_gage_last_: Vec::new(),
            gage_rate_: Vec::new(),
            gage_cum_: Vec::new(),
            gage_window_depth_: Vec::new(),
            gage_window_rate_: Vec::new(),
            gage_since_: 0.0,
            gage_report_depth_: Vec::new(),
            report_rain_cum_: Vec::new(),
            report_rain_ms_: -1.0,
        }
    }
}

impl SurfaceRouter2D {
    pub fn mesh(&self) -> &Shared<MeshData> {
        &self.mesh_
    }
    pub fn options(&self) -> &Shared<SolverOptions2D> {
        &self.options_
    }
    pub fn boundary(&self) -> &Shared<BoundaryData> {
        &self.boundary_
    }
    pub fn infil(&self) -> &Shared<Infil2D> {
        &self.infil_
    }
    pub fn gw_transport(&self) -> &Shared<GwTransportData> {
        &self.gw_
    }
    pub fn surface_quality(&self) -> &Shared<SurfaceQuality2D> {
        &self.surface_quality_
    }
    pub fn aquifer_config(&self) -> &Shared<SubsurfaceConfig> {
        &self.aquifer_cfg_
    }
    pub fn aquifer_node_names(&self) -> &Shared<Vec<Bytes>> {
        &self.aquifer_node_names_
    }
    pub fn aquifer_link_names(&self) -> &Shared<Vec<Bytes>> {
        &self.aquifer_link_names_
    }
    pub fn pending_bc_rows(&self) -> &Shared<Vec<PendingBoundaryRow>> {
        &self.pending_bc_rows_
    }
    pub fn pending_edge_conveyance_rows(&self) -> &Shared<Vec<PendingEdgeConveyanceRow>> {
        &self.pending_edge_conveyance_rows_
    }
    pub fn pending_initial_quality_rows(&self) -> &Shared<Vec<PendingInitialQualityRow>> {
        &self.pending_iq_rows_
    }
    pub fn pending_initial_bed_quality_rows(&self) -> &Shared<Vec<PendingInitialBedQualityRow>> {
        &self.pending_ibq_rows_
    }
    pub fn pending_boundary_quality_rows(&self) -> &Shared<Vec<PendingBoundaryQualityRow>> {
        &self.pending_bq_rows_
    }
    /// isActive()
    pub fn is_active(&self) -> bool {
        self.active_
    }
    /// state()
    pub fn state(&self) -> &SurfaceStateData {
        &self.state_
    }
    /// subsurface()
    pub fn subsurface(&self) -> &SubsurfaceSolver {
        &self.subsurface_
    }
    /// threadWarnings()
    pub fn thread_warnings(&self) -> &Vec<Bytes> {
        &self.thread_warnings_
    }
    /// backendName()
    pub fn backend_name(&self) -> &Bytes {
        &self.backend_name_
    }

    /// The targets register2DSections writes through.
    pub fn sections_2d_targets(&self) -> crate::twod::input::section_handlers_2d::Sections2DTargets {
        crate::twod::input::section_handlers_2d::Sections2DTargets {
            mesh: self.mesh_.clone(),
            options: self.options_.clone(),
            pending_bc_rows: self.pending_bc_rows_.clone(),
            pending_ec_rows: self.pending_edge_conveyance_rows_.clone(),
            pending_iq_rows: self.pending_iq_rows_.clone(),
            pending_ibq_rows: self.pending_ibq_rows_.clone(),
            pending_bq_rows: self.pending_bq_rows_.clone(),
        }
    }

    /// resetWindowAccumulators()
    pub fn reset_window_accumulators(&mut self) {
        for x in self.window_outfall_accum_.iter_mut() {
            *x = 0.0;
        }
        self.outfall_exchange_.reset_window();
    }

    /// totalVolume(): the summed cell volumes (m3).
    pub fn total_volume(&self) -> f64 {
        let mut vol = 0.0;
        let nt = self.mesh_.borrow().n_triangles();
        for i in 0..nt as usize {
            vol += self.state_.volume[i];
        }
        vol
    }
}
