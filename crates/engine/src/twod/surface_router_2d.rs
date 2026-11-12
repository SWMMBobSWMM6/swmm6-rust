// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/SurfaceRouter2D.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The 2D surface router, as far as open() reaches it: the members the [2D_*], [GW_*] and
// [2D_AQUIFER*] handlers fill, in upstream's order and with upstream's names. Later waves add
// the rest of the members (state_, the solver, coupling, the rainfall interpolator, the
// clocks, ...) at their upstream positions, and the methods.
//
// Upstream's SWMMEngine wires SimulationContext::twod_io with raw pointers to these members
// (wire2DModelIO), and the section handlers capture references to them. Here each such member
// is a shared handle (twod::Shared, Rc<RefCell<T>>): the router, the context's twod_io and the
// registered handlers hold clones of the same handle.

#![allow(clippy::all)]

use crate::core::options::Bytes;
use crate::twod::data::boundary_data::BoundaryData;
use crate::twod::data::mesh_data::MeshData;
use crate::twod::data::pending_rows_2d::*;
use crate::twod::data::solver_options_2d::SolverOptions2D;
use crate::twod::gw::gw_transport_data::GwTransportData;
use crate::twod::infil::infil_2d::Infil2D;
use crate::twod::quality::surface_quality_2d::SurfaceQuality2D;
use crate::twod::subsurface::subsurface_data::SubsurfaceConfig;
use crate::twod::{shared, Shared};

#[derive(Debug)]
pub struct SurfaceRouter2D {
    pub mesh_: Shared<MeshData>,
    // state_ (SurfaceStateData), channel_ground_: later waves
    pub options_: Shared<SolverOptions2D>,
    pub boundary_: Shared<BoundaryData>,
    // thread_warnings_, backend_name_: later waves
    pub infil_: Shared<Infil2D>,
    // infil_elapsed_, infil_cum_applied_, coupling_cum_applied_, cell_subcatch_,
    // subcatch_recharge_: later waves
    pub gw_: Shared<GwTransportData>,
    pub surface_quality_: Shared<SurfaceQuality2D>,
    pub aquifer_cfg_: Shared<SubsurfaceConfig>,
    pub aquifer_node_names_: Shared<Vec<Bytes>>,
    pub aquifer_link_names_: Shared<Vec<Bytes>>,
    // subsurface_ (SubsurfaceSolver), aquifer_link_shares_, aquifer_link_ids_, rain_cum_:
    // later waves (the subsurface solver: phase 2b)
    pub pending_bc_rows_: Shared<Vec<PendingBoundaryRow>>,
    pub pending_iq_rows_: Shared<Vec<PendingInitialQualityRow>>,
    pub pending_ibq_rows_: Shared<Vec<PendingInitialBedQualityRow>>,
    pub pending_bq_rows_: Shared<Vec<PendingBoundaryQualityRow>>,
    pub pending_edge_conveyance_rows_: Shared<Vec<PendingEdgeConveyanceRow>>,
    // coupling_points_ and everything after: later waves
}

impl Default for SurfaceRouter2D {
    fn default() -> Self {
        SurfaceRouter2D {
            mesh_: shared(MeshData::default()),
            options_: shared(SolverOptions2D::default()),
            boundary_: shared(BoundaryData::default()),
            infil_: shared(Infil2D::default()),
            gw_: shared(GwTransportData::default()),
            surface_quality_: shared(SurfaceQuality2D::default()),
            aquifer_cfg_: shared(SubsurfaceConfig::default()),
            aquifer_node_names_: shared(Vec::new()),
            aquifer_link_names_: shared(Vec::new()),
            pending_bc_rows_: shared(Vec::new()),
            pending_iq_rows_: shared(Vec::new()),
            pending_ibq_rows_: shared(Vec::new()),
            pending_bq_rows_: shared(Vec::new()),
            pending_edge_conveyance_rows_: shared(Vec::new()),
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
}
