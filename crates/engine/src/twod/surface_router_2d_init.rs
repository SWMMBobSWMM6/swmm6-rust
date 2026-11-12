// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/SurfaceRouter2D.cpp (headFromMeanDepth,
// refreshOutputGradients, drainPendingRows, initialize, configureSnowGroundContacts,
// initializeChannelGround) (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// SurfaceRouter2D::initialize: the mesh to SI, topology and validation, the vertex stencils,
// the rainfall weights, the coupling maps and points (node coupling types, exchange structures,
// the outfall exchange, the coupled-junction ponds, the live exchange points and their head
// cells, the datum check), the thread count, the boundary drain, the solver (factory and
// ExplicitInertialSolver::initialize), 2D infiltration, the initial storage and the render
// fields. Upstream's std::runtime_error / std::invalid_argument throws are
// Err(InitError::Thrown(what)); SWMMEngine::initHydraulics turns them into "2D initialization
// failed: <what>", as its catch does (core/engine_twod.rs).
//
// Refused with the crate's stub where upstream would first act on them (phase 2b; PLAN.md):
//   * species transport on the surface (a model whose surface carries transported rows:
//     pollutants, reaction species, water age or temperature), before the rows are sized;
//   * 2D surface quality (authored [2D_LOADINGS] / coverages: SurfaceQuality2D::resolve);
//   * the two-zone groundwater kernel (an authored [2D_AQUIFER] that runs);
//   * snow on the mesh (MeshSnow::initialize);
//   * MESH_CONTACTS channel heat on an active mesh (ChannelGroundState::initialize);
//   * the development switches OPENSWMM_2D_RAINFALL_MODE and OPENSWMM_2D_FLUX_DH_EPS (and,
//     in the solver factory and the solver, OPENSWMM_2D_BACKEND, _MIN_PARALLEL_CELLS(_DEVICE),
//     _MARCHER_TELEMETRY and OPENSWMM_GPU_PLUGIN_PATH).
// Heat on the mesh, [2D_INITIAL_QUALITY] seeding, [2D_INITIAL_BED_QUALITY] seeding and the
// porewater and ground-heat owners all need a transported surface row or an active aquifer, so
// with those refused upstream's own checks decide them here (an error or an inactive path).
//
// Differences in form, never in result:
//   * The shared members (mesh_, options_, boundary_, infil_, the pending rows) are borrowed
//     for the whole call through clones of their handles; the functions upstream lets read
//     them through ctx.twod_io (transport::surface2DEnables) get them passed instead.
//   * snprintf into a fixed buffer keeps the first size - 1 bytes (snprintf_n).
//   * The erase(remove_if) passes are Vec::retain, which visits the points in order, once each,
//     as std::remove_if does (the predicates' side effects, warnings included, keep their order).
//   * drainPendingRows' unordered_map of edge keys is a HashMap: it is only filled, searched and
//     scanned to flag single-slot edges, which no order can change.

#![allow(clippy::all)]

use std::collections::{HashMap, HashSet};

use crate::cformat;
use crate::core::cmath::{cclamp, cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::engine_step::stub;
use crate::core::node_coupling_2d as nodecoupling2d;
use crate::core::options::{Bytes, RoutingModel};
use crate::core::polygon_rings as polygon;
use crate::core::thread_info;
use crate::core::unit_conversion as ucf;
use crate::data::node_coupling_2d::{NodeCoupling2DKey, NodeCoupling2DType};
use crate::hydrology::infil_bank::infiltration_factors;
use crate::hydrology::ueb_types::MeshSnowConfigMode;
use crate::twod::coupling::node_coupling::{build_coupling_points, resolve_exchange_structure, CouplingPoint};
use crate::twod::coupling::runtime_coupling::RUNTIME_WATER_HEAT_CAPACITY;
use crate::twod::data::boundary_data::{BoundaryData, BoundaryType};
use crate::twod::data::mesh_data::{MeshData, TriCouplingRow};
use crate::twod::data::pending_rows_2d::{PendingBoundaryRow, PendingEdgeConveyanceRow};
use crate::twod::data::solver_options_2d::{CellClosure2D, CouplingHead2D, CouplingLaw2D, CouplingPondMode2D, Momentum2D, RainfallMode, SolverOptions2D};
use crate::twod::data::surface_state_data::{PondBase, SurfaceStateData};
use crate::twod::heat::mesh_heat::{configure_ground_lid_owners, initialize_mesh_heat};
use crate::twod::infil::infil_2d::{infil_2d_method_token, parse_infil_2d_dest, Infil2DDest};
use crate::twod::mesh::mesh_builder::{build_mesh_topology, validate_mesh};
use crate::twod::mesh::quad_vfr::{quad_eta_from_mean_depth, K_QUAD_VFR_Z};
use crate::twod::mesh::rainfall_interpolator::Method as RainMethod;
use crate::twod::mesh::vertex_reconstruction::{build_vertex_stencils, reconstruct_vertex_heads, reconstruct_vertex_render_depths};
use crate::twod::mesh::vfr_closure::{vfr_eta_from_mean_depth, vfr_sort3};
use crate::twod::solver::explicit_inertial_solver::{eprint_bytes, SolverAttachments};
use crate::twod::solver::surface_flux_calculator::{compute_face_velocity, compute_limited_gradients, compute_unlimited_gradients};
use crate::twod::solver::surface_solver_factory::make_surface_solver;
use crate::twod::surface_router_2d::{InitError, SurfaceRouter2D};

fn thrown(m: Bytes) -> InitError {
    InitError::Thrown(m)
}

fn cat(parts: &[&[u8]]) -> Bytes {
    parts.concat()
}

/// snprintf(buf, n, ...): the formatted bytes cut to n - 1.
fn snprintf_n(mut b: Bytes, n: usize) -> Bytes {
    if b.len() > n - 1 {
        b.truncate(n - 1);
    }
    b
}

/// headFromMeanDepth(mesh, opts, i, d): cell i's free surface at mean depth d under the
/// configured closure.
fn head_from_mean_depth(mesh: &MeshData, opts: &SolverOptions2D, i: i32, d: f64) -> f64 {
    if opts.cell_closure == CellClosure2D::Vfr {
        if mesh.cell_vertex_count(i) == 4 {
            let u = i as usize;
            let kz = K_QUAD_VFR_Z as usize;
            return quad_eta_from_mean_depth(&mesh.quad_vfr_z[u * kz..u * kz + kz], mesh.quad_vfr_a[u * 2], mesh.quad_vfr_a[u * 2 + 1], d, opts.vfr_min_wet_frac);
        }
        let mut z1 = mesh.vz[mesh.cell_vertex(i, 0) as usize];
        let mut z2 = mesh.vz[mesh.cell_vertex(i, 1) as usize];
        let mut z3 = mesh.vz[mesh.cell_vertex(i, 2) as usize];
        vfr_sort3(&mut z1, &mut z2, &mut z3);
        return vfr_eta_from_mean_depth(z1, z2, z3, d, opts.vfr_min_wet_frac);
    }
    mesh.tri_cz[i as usize] + d
}

/// refreshOutputGradients(mesh, state, opts): the output gradients when 2D reporting is on.
fn refresh_output_gradients(mesh: &MeshData, state: &mut SurfaceStateData, opts: &SolverOptions2D) {
    if !opts.report_2d {
        return;
    }
    compute_unlimited_gradients(mesh, state, opts.num_threads);
    compute_limited_gradients(mesh, state, opts.limiter_epsilon, opts.num_threads);
}

/// drainPendingRows(): the authored [2D_BOUNDARY_CONDITIONS] / [2D_EDGE_CONVEYANCE] rows into
/// BoundaryData and the mesh edge slots.
fn drain_pending_rows(
    mesh: &mut MeshData,
    boundary: &mut BoundaryData,
    options: &mut SolverOptions2D,
    pending_bc_rows: &[PendingBoundaryRow],
    pending_edge_conveyance_rows: &[PendingEdgeConveyanceRow],
) -> Result<(), InitError> {
    if mesh.n_triangles() < 1 {
        return Ok(()); // nothing to drain into
    }
    boundary.resize(mesh.n_edge_slots());

    // V-E3: drain the [2D_BOUNDARY_CONDITIONS] rows; out-of-range rows are skipped.
    {
        let n_edges = boundary.size();
        for r in pending_bc_rows {
            if r.tri < 0 || r.tri >= mesh.n_triangles() || r.edge < 0 || r.edge >= mesh.cell_vertex_count(r.tri) {
                if r.bc_type == 0 && r.param1 > 0.0 {
                    return Err(thrown(b"WALL_DRAG names an invalid cell/edge.".to_vec()));
                }
                continue;
            }
            let idx = MeshData::slot(r.tri, r.edge);
            if idx < 0 || idx >= n_edges {
                continue;
            }
            let ui = idx as usize;
            boundary.edge_bc_type[ui] = r.bc_type as i8;
            if r.bc_type == BoundaryType::NormalFlow as i32 {
                boundary.edge_bed_slope[ui] = r.param1;
            } else if r.bc_type == BoundaryType::SpecifiedStage as i32 {
                if !r.name.is_empty() {
                    boundary.edge_bc_tseries_name[ui] = r.name.clone();
                    boundary.edge_bc_tseries[ui] = -2; // deferred resolve
                } else {
                    boundary.edge_bc_head[ui] = r.param1;
                }
            } else if r.bc_type == BoundaryType::SpecifiedFlow as i32 {
                if !r.name.is_empty() {
                    boundary.edge_bc_flow_tseries_name[ui] = r.name.clone();
                    boundary.edge_bc_flow_tseries[ui] = -2;
                } else {
                    boundary.edge_bc_flow[ui] = r.param1;
                }
            } else if r.bc_type == BoundaryType::RatingCurve as i32 {
                boundary.edge_bc_rating_curve_name[ui] = r.name.clone();
                boundary.edge_bc_rating_curve[ui] = -2;
            } else if r.bc_type == BoundaryType::Wall as i32 {
                boundary.edge_wall_drag[ui] = r.param1;
            }
        }
    }

    // [2D_EDGE_CONVEYANCE]: a vertex-pair -> slots map, then each factor into every slot.
    {
        let make_key = |va: i32, vb: i32| -> i64 {
            let lo = cmin(va, vb) as i64;
            let hi = cmax(va, vb) as i64;
            (lo << 32) | (hi & 0xFFFF_FFFF)
        };
        let mut edge_key_to_slots: HashMap<i64, Vec<i32>> = HashMap::with_capacity(mesh.n_edge_slots() as usize);
        let nt = mesh.n_triangles();
        for t in 0..nt {
            let nvc = mesh.cell_vertex_count(t);
            for e in 0..nvc {
                let (mut va, mut vb) = (0, 0);
                mesh.cell_edge_vertices(t, e, &mut va, &mut vb);
                edge_key_to_slots.entry(make_key(va, vb)).or_default().push(MeshData::slot(t, e));
            }
        }
        // Exterior faces: the endpoint pairs only one slot names.
        for slots in edge_key_to_slots.values() {
            if slots.len() == 1 {
                boundary.edge_is_exterior[slots[0] as usize] = 1;
            }
        }
        for r in pending_edge_conveyance_rows {
            if r.v_from < 0 || r.v_from >= mesh.n_vertices() || r.v_to < 0 || r.v_to >= mesh.n_vertices() {
                return Err(thrown(format!("[2D_EDGE_CONVEYANCE]: vertex index out of range ({}, {}); n_vertices = {}", r.v_from, r.v_to, mesh.n_vertices()).into_bytes()));
            }
            match edge_key_to_slots.get(&make_key(r.v_from, r.v_to)) {
                None => {
                    return Err(thrown(format!("[2D_EDGE_CONVEYANCE]: vertices ({}, {}) do not form a mesh edge", r.v_from, r.v_to).into_bytes()));
                }
                Some(slots) => {
                    for &slot in slots {
                        mesh.edge_conveyance[slot as usize] = r.conveyance;
                    }
                }
            }
        }
    }

    // The drained arrays are the live state from here on.
    options.wall_drag_enabled = boundary.edge_wall_drag.iter().any(|&c| c > 0.0);
    options.pending_rows_drained = true;
    Ok(())
}

/// The names of the first 8 nodes, joined (", ..." past 8).
fn first_names(ctx: &SimulationContext, nodes: &[i32]) -> Bytes {
    let mut s = Bytes::new();
    let mut i = 0;
    while i < nodes.len() && i < 8 {
        if i != 0 {
            s.extend_from_slice(b", ");
        }
        s.extend_from_slice(ctx.node_names.name_of(nodes[i]));
        i += 1;
    }
    if nodes.len() > 8 {
        s.extend_from_slice(b", ...");
    }
    s
}

impl SurfaceRouter2D {
    /// initialize(ctx).
    pub fn initialize(&mut self, ctx: &mut SimulationContext) -> Result<(), InitError> {
        self.channel_ground_.clear();
        ctx.heat_state.mesh_contact_link.clear();
        self.common_stage_ = false;
        self.source_context_ = false;
        self.joint_hydraulic_ = false;
        self.source_router_model_.clear();
        self.source_boundary_ready_ = false;
        if !self.checkpoint_owned_initialization_ && !self.aquifer_cfg_.borrow().surface_owners.is_empty() {
            return Err(thrown(b"[2D_SURFACE_OWNERSHIP] Reviewed authoring is available; runtime recharge is unavailable until the R4 completed-interval water/quality/ET adapter is qualified. Remove ownership records to run the existing formulation.".to_vec()));
        }
        ctx.surface_outfall_link_limit.clear();
        // Check whether 2D sections were parsed (vertices present)
        {
            let mesh = self.mesh_.borrow();
            if mesh.n_vertices() < 3 || mesh.n_triangles() < 1 {
                self.active_ = false;
                return Ok(());
            }
        }

        // IGNORE_2D: the mesh stays parsed, the solver never activates.
        if ctx.options.ignore_2d {
            eprint_bytes("[openswmm 2D] IGNORE_2D YES \u{2014} 2D surface routing disabled; running 1D-only.\n".as_bytes());
            self.active_ = false;
            return Ok(());
        }

        let mesh_rc = self.mesh_.clone();
        let options_rc = self.options_.clone();
        let boundary_rc = self.boundary_.clone();
        let mut mesh = mesh_rc.borrow_mut();
        let mut options = options_rc.borrow_mut();

        // The 1D side is always feet: the coupling converts feet <-> metres.
        const FT_TO_M: f64 = 0.3048;
        options.len_1d_to_2d = FT_TO_M;
        options.len_2d_to_1d = 1.0 / FT_TO_M;
        options.vol_1d_to_2d = FT_TO_M * FT_TO_M * FT_TO_M;
        options.flow_1d_to_2d = options.vol_1d_to_2d;
        options.flow_2d_to_1d = 1.0 / options.vol_1d_to_2d;

        // The mesh is authored in the project's display length units.
        let us = ucf::get_unit_system(ctx.options.flow_units as i32);
        let mesh_to_si = if us == 0 { FT_TO_M } else { 1.0 };

        self.bc_stage_scale_ = if !options.mesh_units_si && mesh_to_si != 1.0 { mesh_to_si } else { 1.0 };
        self.bc_stage_ts_scale_ = mesh_to_si;
        {
            const K_FLOW_TO_CMS: [f64; 6] = [
                0.028316846592,  // CFS
                6.30901964e-05,  // GPM
                0.043812636574,  // MGD
                1.0,             // CMS
                0.001,           // LPS
                0.011574074074,  // MLD
            ];
            let fu = ctx.options.flow_units as i32;
            self.bc_flow_scale_ = if fu >= 0 && fu < 6 { K_FLOW_TO_CMS[fu as usize] } else { 1.0 };
        }

        options.mesh_to_si_factor = if !options.mesh_units_si && mesh_to_si != 1.0 { mesh_to_si } else { 1.0 };

        // Mesh geometry from project length units to SI, before the topology.
        if !options.mesh_units_si && !options.mesh_scaled_to_si && mesh_to_si != 1.0 {
            let f = mesh_to_si;
            let f2 = f * f;
            for v in mesh.vx.iter_mut() {
                *v *= f;
            }
            for v in mesh.vy.iter_mut() {
                *v *= f;
            }
            for v in mesh.vz.iter_mut() {
                *v *= f;
            }
            for a in mesh.vert_coupling_area.iter_mut() {
                *a *= f2;
            }
            for a in mesh.tri_coupling_area.iter_mut() {
                *a *= f2;
            }
            for r in mesh.tri_couplings.iter_mut() {
                r.area *= f2;
            }
            for d in mesh.tri_init_depth.iter_mut() {
                *d *= f;
            }
            options.mesh_scaled_to_si = true;
        }

        // Build mesh topology (neighbours, edge geometry, areas)
        build_mesh_topology(&mut mesh);

        // Validate mesh
        let err = validate_mesh(&mesh);
        if !err.is_empty() {
            return Err(thrown(cat(&[b"2D mesh validation failed: ", &err])));
        }

        // Pseudo-Laplacian vertex reconstruction stencils
        build_vertex_stencils(&mut mesh);

        // Development switches (refused when set; PLAN.md).
        if std::env::var_os("OPENSWMM_2D_RAINFALL_MODE").is_some() {
            return Err(stub("the OPENSWMM_2D_RAINFALL_MODE development switch (SurfaceRouter2D::initialize)").into());
        }
        if std::env::var_os("OPENSWMM_2D_FLUX_DH_EPS").is_some() {
            return Err(stub("the OPENSWMM_2D_FLUX_DH_EPS development switch (SurfaceRouter2D::initialize)").into());
        }

        // The static rainfall-interpolation weights.
        {
            let method = if options.rainfall_mode == RainfallMode::NearestNeighbour { RainMethod::NearestNeighbour } else { RainMethod::NaturalNeighbour };
            let n_gages = ctx.n_gages();
            self.interp_.build(&mesh.tri_cx, &mesh.tri_cy, &ctx.spatial.gage_x, &ctx.spatial.gage_y, n_gages, mesh_to_si, method).map_err(thrown)?;
        }
        if options.rainfall_mode == RainfallMode::NaturalNeighbour || options.rainfall_mode == RainfallMode::NearestNeighbour {
            let d = self.interp_.diagnostics().clone();
            if d.unlocated != 0 || d.invalid != 0 || d.duplicates != 0 {
                ctx.warnings.push(
                    format!("2D rainfall gages: {} without a location (0,0), {} invalid coordinates, {} duplicate locations (first gage retained).", d.unlocated, d.invalid, d.duplicates).into_bytes(),
                );
            }
            if !self.interp_.ready() {
                ctx.warnings.push(b"2D rainfall: no located gages; using the uniform system mean.".to_vec());
            }
            if d.idwCells != 0 {
                ctx.warnings.push(format!("2D rainfall: inverse-distance fallback used for {} cells outside the gage hull or with degenerate weights.", d.idwCells).into_bytes());
            }
        }

        // Resolve deferred coupling node names -> indices
        for v in 0..mesh.n_vertices() {
            let uv = v as usize;
            if mesh.vert_coupled_node_name[uv].is_empty() {
                continue;
            }
            let node_idx = ctx.node_names.find(&mesh.vert_coupled_node_name[uv]);
            if node_idx >= 0 {
                mesh.vert_coupled_node[uv] = node_idx;
            } else {
                return Err(thrown(cat(&[format!("2D vertex {} coupled to unknown node '", v).as_bytes(), &mesh.vert_coupled_node_name[uv], b"'"])));
            }
        }

        // Cell couplings: the row vector is the source of truth; legacy per-triangle arrays are
        // folded in by synthesising rows first.
        if mesh.tri_couplings.is_empty() {
            for t in 0..mesh.n_triangles() {
                let ut = t as usize;
                let named = !mesh.tri_coupled_node_name[ut].is_empty();
                let indexed = mesh.tri_coupled_node[ut] >= 0;
                if !named && !indexed {
                    continue;
                }
                let mut row = TriCouplingRow::default();
                row.tri = t;
                row.node = if indexed { mesh.tri_coupled_node[ut] } else { -1 };
                row.node_name = mesh.tri_coupled_node_name[ut].clone();
                row.cd = mesh.tri_coupling_cd[ut];
                row.cd_set = true; // the arrays carry an explicit value
                row.area = mesh.tri_coupling_area[ut];
                mesh.tri_couplings.push(row);
            }
        }
        for k in 0..mesh.tri_couplings.len() {
            if mesh.tri_couplings[k].node < 0 {
                if mesh.tri_couplings[k].node_name.is_empty() {
                    continue;
                }
                let node_idx = ctx.node_names.find(&mesh.tri_couplings[k].node_name);
                if node_idx < 0 {
                    let r = &mesh.tri_couplings[k];
                    return Err(thrown(cat(&[format!("2D triangle {} coupled to unknown node '", r.tri).as_bytes(), &r.node_name, b"'"])));
                }
                mesh.tri_couplings[k].node = node_idx;
            }
            // Mirror into the legacy arrays (last row wins).
            let (tri, node, cd, area) = {
                let r = &mesh.tri_couplings[k];
                (r.tri as usize, r.node, r.cd, r.area)
            };
            mesh.tri_coupled_node[tri] = node;
            mesh.tri_coupling_cd[tri] = cd;
            mesh.tri_coupling_area[tri] = area;
        }

        // Initialize surface state
        self.state_.resize(mesh.n_triangles(), mesh.n_vertices());

        // Drain the pending boundary / conveyance rows (sizes boundary_, flips the flag).
        let mut boundary = boundary_rc.borrow_mut();
        {
            let pbc = self.pending_bc_rows_.borrow();
            let pec = self.pending_edge_conveyance_rows_.borrow();
            drain_pending_rows(&mut mesh, &mut boundary, &mut options, &pbc, &pec)?;
        }

        // Constant stage heads / flow discharges to SI (table-driven entries scale at lookup).
        for idx in 0..boundary.size() as usize {
            let bt = boundary.edge_bc_type[idx];
            if bt == BoundaryType::SpecifiedStage as i8 && boundary.edge_bc_tseries[idx] == -1 {
                boundary.edge_bc_head[idx] *= self.bc_stage_scale_;
            } else if bt == BoundaryType::SpecifiedFlow as i8 && boundary.edge_bc_flow_tseries[idx] == -1 {
                boundary.edge_bc_flow[idx] *= self.bc_flow_scale_;
            }
        }
        options.bc_flow_to_si_applied = self.bc_flow_scale_;

        // A NORMAL_FLOW edge with zero bed slope behaves as a wall: warn once with a count.
        {
            let mut n_zero_slope = 0;
            for i in 0..boundary.edge_bc_type.len() {
                if boundary.edge_bc_type[i] == BoundaryType::NormalFlow as i8 && boundary.edge_bed_slope[i] == 0.0 {
                    n_zero_slope += 1;
                }
            }
            if n_zero_slope > 0 {
                ctx.warnings.push(
                    format!(
                        "WARNING: {} 2D NORMAL_FLOW boundary edge(s) have zero bed slope (PARAM_1) \u{2014} Manning outflow is zero, so they behave as WALL edges. Author a non-zero slope to enable outflow.",
                        n_zero_slope
                    )
                    .into_bytes(),
                );
            }
        }

        // Initial heads (bed plus INIT_DEPTH under the closure) and volumes.
        for i in 0..mesh.n_triangles() {
            let ui = i as usize;
            self.state_.head[ui] = head_from_mean_depth(&mesh, &options, i, mesh.tri_init_depth[ui]);
            self.state_.volume[ui] = mesh.tri_init_depth[ui] * mesh.tri_area[ui];
        }
        reconstruct_vertex_heads(&mesh, &mut self.state_, options.num_threads);

        // Coupling point descriptors; a point without an authored CD takes COUPLING_CD.
        self.coupling_points_ = build_coupling_points(&mesh, ctx);
        let cd_default = if options.coupling_cd > 0.0 { options.coupling_cd } else { 0.65 };
        let mut cd_fallbacks: usize = 0;
        for cp in self.coupling_points_.iter_mut() {
            if !cp.cd_authored {
                cp.cd = cd_default;
                cd_fallbacks += 1;
            }
        }
        if !(options.coupling_cd > 0.0) && cd_fallbacks > 0 {
            let b = cformat!(b"2D_OPTIONS COUPLING_CD %g is not a discharge coefficient; %zu coupling point(s) without a CD use 0.65.", options.coupling_cd, cd_fallbacks);
            ctx.warnings.push(snprintf_n(b, 200));
        }

        // [2D_NODE_COUPLING]: a node's coupling type decides whether its map rows couple.
        {
            let nn = ctx.n_nodes() as usize;
            let mut row_of_node = vec![-1i32; nn];
            for r in 0..ctx.node_coupling_2d.count() {
                let n = ctx.node_coupling_2d.node[r as usize];
                let why = nodecoupling2d::check_type(ctx, n, ctx.node_coupling_2d.type_[r as usize]);
                if !why.is_empty() {
                    return Err(thrown(cat(&[b"[2D_NODE_COUPLING] node '", ctx.node_names.name_of(n), b"': ", &why])));
                }
                row_of_node[n as usize] = r;
            }
            let type_of = |ctx: &SimulationContext, n: i32| -> NodeCoupling2DType {
                let r = row_of_node[n as usize];
                if r >= 0 {
                    ctx.node_coupling_2d.type_[r as usize]
                } else {
                    nodecoupling2d::default_type(ctx, n)
                }
            };

            let mut ignored: Vec<i32> = Vec::new();
            let mut mapped = vec![0u8; nn];
            {
                let ctx_ref: &SimulationContext = ctx;
                self.coupling_points_.retain(|cp| {
                    if type_of(ctx_ref, cp.node_idx) != NodeCoupling2DType::None {
                        mapped[cp.node_idx as usize] = 1;
                        return true;
                    }
                    if !ignored.contains(&cp.node_idx) {
                        ignored.push(cp.node_idx);
                    }
                    false
                });
            }
            if !ignored.is_empty() {
                let m = cat(&[
                    format!(
                        "WARNING: 2D map rows ignored for {} node(s) that are not coupleable ([2D_NODE_COUPLING] type NONE, explicit or the inlet/virtual-junction default): ",
                        ignored.len()
                    )
                    .as_bytes(),
                    &first_names(ctx, &ignored),
                ]);
                ctx.warnings.push(m);
            }

            let mut unmapped: Vec<i32> = Vec::new();
            for r in 0..ctx.node_coupling_2d.count() as usize {
                let n = ctx.node_coupling_2d.node[r];
                if ctx.node_coupling_2d.type_[r] != NodeCoupling2DType::None && mapped[n as usize] == 0 {
                    unmapped.push(n);
                }
            }
            if !unmapped.is_empty() {
                let m = cat(&[
                    format!(
                        "WARNING: {} node(s) authored as coupleable in [2D_NODE_COUPLING] have no [2D_VERTEX_NODE_MAP] / [2D_TRIANGLE_NODE_MAP] row and do not couple: ",
                        unmapped.len()
                    )
                    .as_bytes(),
                    &first_names(ctx, &unmapped),
                ]);
                ctx.warnings.push(m);
            }

            // Exchange law, structure and coefficients of every junction point; an outfall's
            // SPREAD (project length -> m) and TAILWATER_RELAX (s).
            for k in 0..self.coupling_points_.len() {
                let node = self.coupling_points_[k].node_idx;
                let row = row_of_node[node as usize];
                if self.coupling_points_[k].is_outfall {
                    if row < 0 {
                        continue;
                    }
                    let p = &ctx.node_coupling_2d.params[row as usize];
                    let spread = p[NodeCoupling2DKey::Spread as usize];
                    let relax = p[NodeCoupling2DKey::TailwaterRelax as usize];
                    let cp = &mut self.coupling_points_[k];
                    if spread > 0.0 {
                        cp.spread = spread * mesh_to_si;
                    }
                    if relax > 0.0 {
                        cp.tw_relax = relax;
                    }
                    continue;
                }
                let t = type_of(ctx, node);
                resolve_exchange_structure(&mut self.coupling_points_[k], &mesh, ctx, &options, row, t, mesh_to_si);
            }
        }

        // Unused outfalls (no incident hydraulic link) cannot exchange: drop their mappings.
        {
            let mut has_link = vec![0u8; ctx.n_nodes() as usize];
            for j in 0..ctx.n_links() as usize {
                if ctx.links.node1[j] >= 0 {
                    has_link[ctx.links.node1[j] as usize] = 1;
                }
                if ctx.links.node2[j] >= 0 {
                    has_link[ctx.links.node2[j] as usize] = 1;
                }
            }
            let mut warned_outfalls: HashSet<i32> = HashSet::new();
            let mut warnings: Vec<Bytes> = Vec::new();
            {
                let ctx_ref: &SimulationContext = ctx;
                self.coupling_points_.retain(|cp| {
                    if !cp.is_outfall || has_link[cp.node_idx as usize] != 0 {
                        return true;
                    }
                    if warned_outfalls.insert(cp.node_idx) {
                        warnings.push(cat(&[
                            b"WARNING: 2D coupling ignored for outfall '",
                            ctx_ref.node_names.name_of(cp.node_idx),
                            b"': no incident hydraulic link; all vertex and cell mappings to this outfall are ignored.",
                        ]));
                    }
                    false
                });
            }
            ctx.warnings.extend(warnings);
        }

        // Auto-align ponded storage on 2D-coupled junctions with the 2D surface.
        ctx.coupled_node.clear();
        ctx.coupled_node.resize(ctx.n_nodes() as usize, 0u8);

        // 2D cell areas are SI m2; 1D ponded_area is ft2.
        let area_2d_to_1d = options.len_2d_to_1d * options.len_2d_to_1d;

        // Largest conduit area at each node, SI m2.
        let mut a_pipe_max_m2 = vec![0.0f64; ctx.n_nodes() as usize];
        for l in 0..ctx.n_links() as usize {
            for n in [ctx.links.node1[l], ctx.links.node2[l]] {
                if n < 0 {
                    continue;
                }
                let a = &mut a_pipe_max_m2[n as usize];
                *a = cmax(*a, ctx.links.xsect_a_full[l]);
            }
        }
        for a in a_pipe_max_m2.iter_mut() {
            *a *= options.len_1d_to_2d * options.len_1d_to_2d;
        }

        for k in 0..self.coupling_points_.len() {
            if self.coupling_points_[k].is_outfall {
                continue;
            }
            let ni = self.coupling_points_[k].node_idx as usize;
            let a_pipe_max = a_pipe_max_m2[ni];

            // COUPLING_AREA AUTO: the exchange area from the largest connected conduit.
            {
                let cp = &mut self.coupling_points_[k];
                if options.coupling_area_auto && a_pipe_max > 0.0 && cp.law == CouplingLaw2D::Orifice as i8 {
                    cp.area = cclamp(1.25 * a_pipe_max, 0.05, 2.0);
                }
            }
            let cp = self.coupling_points_[k].clone();

            // First time this node is touched: warn about overridden values and reset
            // ponded_area before accumulating the footprint.
            if ctx.coupled_node[ni] == 0 {
                let nname = ctx.node_names.name_of(cp.node_idx).to_vec();
                if ctx.nodes.sur_depth[ni] > 0.0 {
                    ctx.warnings.push(cat(&[
                        b"WARNING: 2D-coupled node '",
                        &nname,
                        "' has sur_depth > 0 \u{2014} ignored while the node is 2D-coupled: ponding and the exchange gate start at the rim (invert + full_depth).".as_bytes(),
                    ]));
                }
                if ctx.nodes.ponded_area[ni] > 0.0 {
                    ctx.warnings.push(cat(&[
                        b"WARNING: 2D-coupled node '",
                        &nname,
                        "' has ponded_area > 0 \u{2014} it is being overridden with the surrounding 2D-cell footprint so the 1D HGL stays aligned with the 2D surface.".as_bytes(),
                    ]));
                }
                ctx.nodes.ponded_area[ni] = 0.0;
                ctx.coupled_node[ni] = 1;

                // Exchange-area sanity (ORIFICE points only).
                if cp.law == CouplingLaw2D::Orifice as i8 && a_pipe_max > 0.0 && cp.area > 10.0 * a_pipe_max {
                    let b = cformat!(
                        "WARNING: 2D-coupled node '%s' has exchange AREA %.3f m\u{b2} \u{2014} %.0fx the largest connected conduit area (%.4f m\u{b2}). The orifice can inject far more than the pipe can convey; expect the node to fill and spill back each window. Consider AREA ~ the inlet/barrel area."
                            .as_bytes(),
                        &nname[..],
                        cp.area,
                        cp.area / a_pipe_max,
                        a_pipe_max
                    );
                    ctx.warnings.push(snprintf_n(b, 320));
                }
            }

            // Median-dual footprint of the cells around the coupling point (SI m2).
            let foot_m2;
            if cp.vertex_idx >= 0 {
                let s = mesh.vert_stencil_ptr[cp.vertex_idx as usize];
                let e = mesh.vert_stencil_ptr[cp.vertex_idx as usize + 1];
                let mut foot_tri = 0.0;
                let mut foot_quad = 0.0;
                for kk in s..e {
                    let c = mesh.vert_stencil_idx[kk as usize];
                    if mesh.cell_vertex_count(c) == 4 {
                        foot_quad += mesh.tri_area[c as usize];
                    } else {
                        foot_tri += mesh.tri_area[c as usize];
                    }
                }
                foot_m2 = foot_tri / 3.0 + foot_quad / 4.0;
            } else {
                foot_m2 = mesh.tri_area[cp.cell_idx as usize];
            }
            ctx.nodes.ponded_area[ni] += foot_m2 * area_2d_to_1d;
        }

        // Above-rim pond of the coupled junctions.
        ctx.coupled_pond_area.clear();
        ctx.coupled_pond_band.clear();
        self.pond_track_ptr_.clear();
        self.pond_track_cell_.clear();
        self.pond_track_share_.clear();
        self.pond_chamber_ft2_.clear();
        {
            let mode = options.coupling_pond_mode;
            let band_m = if options.coupling_rim_band >= 0.0 {
                options.coupling_rim_band
            } else if mode == CouplingPondMode2D::Footprint {
                0.0
            } else {
                0.05
            };
            let configured = mode != CouplingPondMode2D::Footprint || band_m > 0.0;
            if configured && ctx.options.routing_model != RoutingModel::Dynwave {
                ctx.warnings.push(b"WARNING: COUPLING_POND_MODE / COUPLING_RIM_BAND apply to dynamic wave routing; coupled junctions keep the footprint step.".to_vec());
            } else if configured {
                let n = ctx.n_nodes() as usize;
                let mut chamber_m2 = vec![0.0f64; n];
                for cp in &self.coupling_points_ {
                    if !cp.is_outfall {
                        let c = &mut chamber_m2[cp.node_idx as usize];
                        *c = cmax(*c, cp.chamber_area);
                    }
                }
                ctx.coupled_pond_area.clear();
                ctx.coupled_pond_area.resize(n, 0.0);
                ctx.coupled_pond_band.clear();
                ctx.coupled_pond_band.resize(n, band_m * options.len_2d_to_1d);
                self.pond_chamber_ft2_.clear();
                self.pond_chamber_ft2_.resize(n, 0.0);
                for i in 0..n {
                    if ctx.coupled_node[i] == 0 {
                        continue;
                    }
                    let chamber = if chamber_m2[i] > 0.0 { chamber_m2[i] * area_2d_to_1d } else { ctx.nodes.ponded_area[i] };
                    self.pond_chamber_ft2_[i] = chamber;
                    ctx.coupled_pond_area[i] = if mode == CouplingPondMode2D::Footprint { ctx.nodes.ponded_area[i] } else { chamber };
                }
                if mode == CouplingPondMode2D::Tracking {
                    // Each coupled node tracks the median-dual share of the cells around its
                    // coupling points.
                    let mut cells: Vec<Vec<(i32, f64)>> = vec![Vec::new(); n];
                    for cp in &self.coupling_points_ {
                        if cp.is_outfall {
                            continue;
                        }
                        let list = &mut cells[cp.node_idx as usize];
                        if cp.vertex_idx >= 0 {
                            let v = cp.vertex_idx as usize;
                            for kk in mesh.vert_stencil_ptr[v]..mesh.vert_stencil_ptr[v + 1] {
                                let c = mesh.vert_stencil_idx[kk as usize];
                                list.push((c, mesh.tri_area[c as usize] / mesh.cell_vertex_count(c) as f64));
                            }
                        } else {
                            list.push((cp.cell_idx, mesh.tri_area[cp.cell_idx as usize]));
                        }
                    }
                    self.pond_track_ptr_.clear();
                    self.pond_track_ptr_.resize(n + 1, 0);
                    for i in 0..n {
                        for &(c, share) in &cells[i] {
                            self.pond_track_cell_.push(c);
                            self.pond_track_share_.push(share);
                        }
                        self.pond_track_ptr_[i + 1] = self.pond_track_cell_.len() as i32;
                    }
                }
            }
        }

        // Live in-marcher exchange points: every non-outfall point becomes a single-cell point
        // (the lowest-bed incident cell), with its COUPLING_HEAD sampling cells.
        {
            self.node_coupling_points_.clear();
            self.coupling_head_cells_.clear();
            let mut seen = vec![0u8; mesh.n_cells() as usize];
            let coupling_head = options.coupling_head;
            let add_head_cells = |cp: &CouplingPoint, cell: i32, out: &mut Vec<i32>, seen: &mut Vec<u8>| {
                let mut start: Vec<i32> = Vec::new();
                if cp.vertex_idx >= 0 {
                    let v = cp.vertex_idx as usize;
                    for kk in mesh.vert_stencil_ptr[v]..mesh.vert_stencil_ptr[v + 1] {
                        start.push(mesh.vert_stencil_idx[kk as usize]);
                    }
                } else {
                    start.push(cell);
                }
                if coupling_head == CouplingHead2D::Vertex {
                    if cp.vertex_idx >= 0 {
                        out.extend_from_slice(&start);
                    }
                    return;
                }
                const K_PI: f64 = 3.14159265358979323846;
                let cx = if cp.vertex_idx >= 0 { mesh.vx[cp.vertex_idx as usize] } else { mesh.tri_cx[cell as usize] };
                let cy = if cp.vertex_idx >= 0 { mesh.vy[cp.vertex_idx as usize] } else { mesh.tri_cy[cell as usize] };
                let d_struct = if cp.law == CouplingLaw2D::WeirOrifice as i8 { cp.geo.l_weir / K_PI } else { 2.0 * (cp.area / K_PI).sqrt() };
                let r = cmax(d_struct, mesh.tri_area[cell as usize].sqrt());
                const K_MAX_CELLS: usize = 256;
                let mut queue: Vec<i32> = start.clone();
                let mut marked: Vec<i32> = start.clone();
                for &c in &start {
                    seen[c as usize] = 1;
                }
                let mut q = 0usize;
                while q < queue.len() && q < K_MAX_CELLS {
                    let c = queue[q];
                    for kk in 0..mesh.cell_vertex_count(c) {
                        let nb = mesh.cell_neighbour(c, kk);
                        if nb < 0 || seen[nb as usize] != 0 {
                            continue;
                        }
                        seen[nb as usize] = 1;
                        marked.push(nb);
                        let dx = mesh.tri_cx[nb as usize] - cx;
                        let dy = mesh.tri_cy[nb as usize] - cy;
                        if dx * dx + dy * dy <= r * r {
                            queue.push(nb);
                        }
                    }
                    q += 1;
                }
                for &c in &marked {
                    seen[c as usize] = 0;
                }
                if queue.len() > K_MAX_CELLS {
                    queue.truncate(K_MAX_CELLS);
                }
                out.extend_from_slice(&queue);
            };
            for k in 0..self.coupling_points_.len() {
                if self.coupling_points_[k].is_outfall {
                    // An outfall's tailwater samples the same way.
                    if coupling_head != CouplingHead2D::Cell {
                        let cp = self.coupling_points_[k].clone();
                        let begin = self.coupling_head_cells_.len() as i32;
                        add_head_cells(&cp, cp.cell_idx, &mut self.coupling_head_cells_, &mut seen);
                        let end = self.coupling_head_cells_.len() as i32;
                        self.coupling_points_[k].head_begin = begin;
                        self.coupling_points_[k].head_end = end;
                    }
                    continue;
                }
                let cp = self.coupling_points_[k].clone();
                let mut sc = cp.clone();
                if coupling_head != CouplingHead2D::Cell {
                    sc.head_begin = self.coupling_head_cells_.len() as i32;
                    add_head_cells(&cp, cp.cell_idx, &mut self.coupling_head_cells_, &mut seen);
                    sc.head_end = self.coupling_head_cells_.len() as i32;
                }
                if cp.vertex_idx >= 0 {
                    let v = cp.vertex_idx as usize;
                    let s = mesh.vert_stencil_ptr[v];
                    let e = mesh.vert_stencil_ptr[v + 1];
                    let mut lo = cp.cell_idx;
                    let mut zlo = if lo >= 0 { mesh.tri_cz[lo as usize] } else { 1.0e300 };
                    for kk in s..e {
                        let t = mesh.vert_stencil_idx[kk as usize];
                        if mesh.tri_cz[t as usize] < zlo {
                            zlo = mesh.tri_cz[t as usize];
                            lo = t;
                        }
                    }
                    sc.cell_idx = lo;
                    sc.vertex_idx = -1;
                }
                self.node_coupling_points_.push(sc);
            }
            self.state_.coupling_head_cells = !self.coupling_head_cells_.is_empty();
            if !self.node_coupling_points_.is_empty() {
                self.state_.node_coupling = true;
                self.state_.nodes_1d = true;
                self.state_.node_pond_base = if ctx.options.routing_model == RoutingModel::Dynwave { PondBase::RptFullVolume } else { PondBase::FullVolume };
            }
        }

        // C3b: vertical-datum consistency of the coupled nodes and the mesh.
        if !mesh.vz.is_empty() {
            let ft_to_m = options.len_1d_to_2d;
            let mut zmin = mesh.vz[0];
            let mut zmax = mesh.vz[0];
            for &z in &mesh.vz {
                zmin = cmin(zmin, z);
                zmax = cmax(zmax, z);
            }
            let margin = cmax(10.0, 10.0 * (zmax - zmin)); // metres
            for cp in &self.coupling_points_ {
                let ni = cp.node_idx as usize;
                let bed_z = if cp.vertex_idx >= 0 { mesh.vz[cp.vertex_idx as usize] } else { mesh.tri_cz[cp.cell_idx as usize] };
                let rim_m = (ctx.nodes.invert_elev[ni] + ctx.nodes.full_depth[ni]) * ft_to_m;
                if rim_m < zmin - margin || rim_m > zmax + margin {
                    let b = cformat!(
                        "WARNING: 2D-coupled node '%s' is on a different vertical datum than the 2D mesh: node ground (invert+MaxDepth) = %.2f m, but the mesh spans %.2f..%.2f m (bed = %.2f m at the coupling cell) \u{2014} a ~%.1f m offset. The coupling head dh = h_2d - h_1d is dominated by this offset and will drive a spurious exchange (it can drown an outfall or back the 1D network up). Check that the 1D inverts/MaxDepth and the 2D mesh elevations use the same datum and units."
                            .as_bytes(),
                        ctx.node_names.name_of(cp.node_idx),
                        rim_m,
                        zmin,
                        zmax,
                        bed_z,
                        rim_m - bed_z
                    );
                    ctx.warnings.push(snprintf_n(b, 512));
                }
            }
        }

        // The 2D solver's thread count (the 1D solvers' resolution rules, the triangle gate).
        self.thread_warnings_.clear();
        options.requested_threads = ctx.options.num_threads;
        options.num_threads = thread_info::two_d_threads(ctx.options.num_threads, mesh.n_triangles(), Some(&mut self.thread_warnings_));

        // Attach the per-edge boundary conditions to the state.
        self.state_.boundary = true;

        // Seed the output gradients.
        refresh_output_gradients(&mesh, &mut self.state_, &options);

        // ---- OPENSWMM_HAS_2D ----
        // Outfall batch ledger and per-routing-step donor reservations.
        self.window_outfall_accum_.clear();
        self.window_outfall_accum_.resize(self.coupling_points_.len(), 0.0);
        self.outfall_exchange_.initialize(&mesh, &self.coupling_points_, ctx).map_err(thrown)?;
        self.outfall_state_changed_ = false;

        // The surface's transported rows (TransportPolicy's canonical layout).
        {
            let en = crate::transport::transport_policy::surface_2d_enables(ctx, Some(&options));
            let l = crate::transport::transport_policy::canonical_rows(ctx, &en);
            if l.ns > 0 {
                return Err(stub("2D species transport on the surface (transported pollutant, reaction, age or temperature rows; phase 2b)").into());
            }
            if !ctx.options.ignore_quality && options.transport_msx && crate::transport::reaction_ard_binding::ard_reactions_active(ctx) && en.msx_has_wall {
                ctx.warnings.push(
                    "2D transport: WALL species have no transport semantics on the surface yet \u{2014} MSX rows are not carried on the mesh (the LEGACY element-local binding still runs them in 1D)."
                        .as_bytes()
                        .to_vec(),
                );
            }
            let ns = l.ns;
            self.state_.transport.resize(ns, mesh.n_triangles(), self.node_coupling_points_.len() as i32);
            let tr = &mut self.state_.transport;
            tr.n_pollut = l.n_pollut;
            tr.n_msx = l.n_msx;
            tr.age_row = l.age_row;
            tr.temp_row = l.temp_row;
            tr.row_names = l.names.clone();
            ctx.nodes.coupling_tuple_age = tr.age_row >= 0;
            ctx.nodes.coupling_tuple_temp = tr.temp_row >= 0;
            if !tr.active() && !self.pending_bq_rows_.borrow().is_empty() {
                return Err(thrown(
                    b"[2D_BOUNDARY_QUALITY] rows are present but the model carries no transported species ([POLLUTANTS] empty or IGNORE_QUALITY set), so they could carry nothing.".to_vec(),
                ));
            }
            // S4: the node row-value publisher.
            self.node_row_conc_.clear();
            self.node_row_conc_.resize(ns as usize * ctx.n_nodes() as usize, 0.0);
            self.state_.node_row_conc = true;
        }

        // S7: cell coverages -> buildup / washoff (2D surface quality: phase 2b).
        if self.surface_quality_.borrow().authored() {
            return Err(stub("2D surface quality (SurfaceQuality2D::resolve; phase 2b)").into());
        }

        let (porous_heat_requested, gw_wanted, groundwater_requested) = {
            let aq = self.aquifer_cfg_.borrow();
            let gw = self.gw_.borrow();
            let porous = ctx.options.heat_transport && options.groundwater != 0 && !aq.empty() && gw.options.thermal_model.as_slice() == b"POROUS_EQUILIBRIUM" && gw.options.transport_temperature;
            let wanted = if options.groundwater < 0 { !aq.empty() } else { options.groundwater > 0 };
            let requested = wanted && !aq.empty();
            (porous, wanted, requested)
        };
        if self.solver_.is_none() {
            let plugin_ok = mesh.n_quads() == 0 && options.momentum == Momentum2D::LocalInertial && !self.state_.transport.active() && !porous_heat_requested && !groundwater_requested;
            let solver = make_surface_solver(&options, &mut self.backend_name_, mesh.n_triangles(), plugin_ok)?;
            self.solver_ = Some(solver);
            let closure: &[u8] = match options.momentum {
                Momentum2D::FullSwe => b"FULL_SWE",
                Momentum2D::DiffusiveWave => b"DIFFUSIVE_WAVE",
                _ => b"LOCAL_INERTIAL",
            };
            let tail: &[u8] = if plugin_ok {
                b""
            } else if porous_heat_requested && !self.state_.transport.active() {
                " \u{2014} porous-ground heat requires the CPU marcher".as_bytes()
            } else if groundwater_requested {
                " \u{2014} 2D groundwater requires the CPU marcher".as_bytes()
            } else if self.state_.transport.active() {
                " \u{2014} the GPU/Kokkos plugin carries no transport rows (species / age / temperature); CPU marcher used".as_bytes()
            } else {
                " \u{2014} the GPU/Kokkos plugin serves all-triangle LOCAL_INERTIAL models only".as_bytes()
            };
            let b = cformat!(
                b"2D solver: %s; MOMENTUM_EQUATION %s; LTS_TIERS %d; %d cells (%d quads)%s",
                &self.backend_name_[..],
                closure,
                options.lts_tiers,
                mesh.n_triangles(),
                mesh.n_quads(),
                tail
            );
            self.thread_warnings_.push(snprintf_n(b, 320));
            if options.momentum_physics() {
                let b = cformat!(
                    b"2D optional momentum: f=%.17g s^-1; turbulence=%s; Cs=%.17g; mixing-length c=%.17g; background nu=%.17g m2/s; interior eddy-viscosity stress; selected quadratic wall traction; accepted cell velocities.",
                    options.coriolis_parameter(),
                    crate::twod::data::solver_options_2d::turbulence_name(options.turbulence),
                    options.smagorinsky_cs,
                    options.mixing_length_coefficient,
                    options.eddy_viscosity_background
                );
                self.thread_warnings_.push(snprintf_n(b, 320));
            }
        }
        {
            let aq = self.aquifer_cfg_.borrow();
            let gw = self.gw_.borrow();
            options.water_heat_capacity = if ctx.options.heat_transport
                && ((options.coupled_heat && self.state_.transport.temp_row >= 0)
                    || (options.groundwater != 0 && !aq.empty() && gw.options.thermal_model.as_slice() == b"POROUS_EQUILIBRIUM" && gw.options.transport_temperature))
            {
                ctx.options.water_density * ctx.options.water_specific_heat
            } else {
                RUNTIME_WATER_HEAT_CAPACITY
            };
        }
        if !(options.water_heat_capacity > 0.0) || !options.water_heat_capacity.is_finite() {
            return Err(thrown(b"[2D heat] Water density and specific heat must give a finite positive capacity.".to_vec()));
        }
        {
            let solver = self.solver_.as_mut().expect("solver_ is set above");
            let att = SolverAttachments {
                boundary: if self.state_.boundary { Some(&*boundary) } else { None },
                node_coupling: if self.state_.node_coupling { Some(&self.node_coupling_points_[..]) } else { None },
                nodes_1d_len: if self.state_.nodes_1d { Some(ctx.nodes.volume.len()) } else { None },
            };
            solver.initialize(&mut mesh, &mut self.state_, &att, &mut options)?;
            // The CPU marcher pins the outfall donor cells to tier 0.
            let mut cells: Vec<i32> = Vec::new();
            for cp in &self.coupling_points_ {
                if !cp.is_outfall {
                    continue;
                }
                if cp.vertex_idx < 0 {
                    cells.push(cp.cell_idx);
                } else {
                    let v = cp.vertex_idx as usize;
                    for p in mesh.vert_stencil_ptr[v]..mesh.vert_stencil_ptr[v + 1] {
                        cells.push(mesh.vert_stencil_idx[p as usize]);
                    }
                }
            }
            solver.pin_external_source_cells(&cells);
        }

        // G1: the two-zone groundwater kernel (phase 2b).
        if options.groundwater == 0 && !self.aquifer_cfg_.borrow().empty() {
            ctx.warnings.push(b"[2D_OPTIONS] GROUNDWATER NO \xe2\x80\x94 the [2D_AQUIFER*] rows are kept but no subsurface runs this simulation.".to_vec());
        }
        if gw_wanted && !self.aquifer_cfg_.borrow().empty() {
            return Err(stub("the 2D two-zone groundwater kernel ([2D_AQUIFER]: resolveSubsurface / SubsurfaceSolver::initialize; phase 2b)").into());
        }

        // Track I: resolve the [2D_INFILTRATION*] rows against the mesh.
        {
            let infil_rc = self.infil_.clone();
            let mut infil = infil_rc.borrow_mut();
            // INFIL_DESTINATION fills every row that did not spell its own DEST.
            if !options.infil_destination.is_empty() && options.infil_destination.as_slice() != b"LOST" {
                let mut dest = Infil2DDest::Lost;
                if parse_infil_2d_dest(&options.infil_destination, &mut dest) {
                    for d in infil.defaults_mut().iter_mut() {
                        if !d.row.dest_explicit {
                            d.row.dest = dest;
                        }
                    }
                    for o in infil.overrides_mut().iter_mut() {
                        if !o.row.dest_explicit {
                            o.row.dest = dest;
                        }
                    }
                }
            }
            if !options.infil_default_method.is_empty() {
                let star = infil.defaults().iter().position(|d| d.tag.as_slice() == b"*");
                if options.infil_default_method.as_slice() == b"NONE" {
                    if let Some(si) = star {
                        ctx.warnings.push(
                            b"[2D_OPTIONS] INFIL_DEFAULT_METHOD NONE \xe2\x80\x94 the '*' row of [2D_INFILTRATION_DEFAULTS] is not applied this run (tag and per-cell rows still are)."
                                .to_vec(),
                        );
                        infil.defaults_mut().remove(si);
                    }
                } else if star.is_none() {
                    ctx.warnings.push(cat(&[
                        b"[2D_OPTIONS] INFIL_DEFAULT_METHOD ",
                        &options.infil_default_method,
                        b" names no '*' row in [2D_INFILTRATION_DEFAULTS] \xe2\x80\x94 the method needs that row's parameters; nothing infiltrates by default.",
                    ]));
                } else {
                    let token = infil_2d_method_token(&infil.defaults()[star.unwrap()].row);
                    if options.infil_default_method.as_slice() != token {
                        return Err(thrown(cat(&[
                            b"[2D_OPTIONS] INFIL_DEFAULT_METHOD ",
                            &options.infil_default_method,
                            b" conflicts with the '*' row of [2D_INFILTRATION_DEFAULTS] (",
                            token,
                            b") \xe2\x80\x94 edit the row, whose parameters belong to its method.",
                        ])));
                    }
                }
            }
            // R2: coverage selects the aquifer owner.
            infil.set_aquifer_2d_available(self.subsurface_.active());
            {
                let mut owners = vec![if self.subsurface_.active() && options.infiltration != 0 { 1u8 } else { 0u8 }; mesh.n_triangles() as usize];
                let covered = self.subsurface_.coverage();
                if covered.len() == owners.len() {
                    for i in 0..owners.len() {
                        owners[i] &= covered[i];
                    }
                }
                infil.set_aquifer_owners(owners);
            }
            if options.infil_destination.as_slice() == b"AQUIFER_2D" {
                return Err(thrown(
                    b"[2D_OPTIONS] INFIL_DESTINATION AQUIFER_2D is obsolete; remove this option. Aquifer-owned cells compute their receiving capacity.".to_vec(),
                ));
            }
            let mut infil_err = Bytes::new();
            if !infil.resolve(&mesh, &ctx.options, &mut infil_err) {
                return Err(thrown(infil_err));
            }
            for message in infil.ownership_messages().clone() {
                ctx.warnings.push(message);
            }
            // U3: cell -> containing subcatchment for every SUBCATCH_AQUIFER cell.
            self.cell_subcatch_.clear();
            self.subcatch_recharge_.clear();
            {
                let mut any_aq = false;
                for r in infil.resolved_rows() {
                    if r.has_method && r.dest == Infil2DDest::SubcatchAquifer {
                        any_aq = true;
                        break;
                    }
                }
                if any_aq {
                    let nsub = ctx.n_subcatches();
                    if nsub <= 0 {
                        return Err(thrown(
                            b"[2D_OPTIONS] INFIL_DESTINATION SUBCATCH_AQUIFER needs subcatchments with aquifers to receive the recharge; this model has none.".to_vec(),
                        ));
                    }
                    for spec in &ctx.process_component_specs {
                        if spec.id.as_slice() == b"org.hydrocouple.openswmm.integrated2d" {
                            return Err(thrown(
                                b"[2D_OPTIONS] INFIL_DESTINATION SUBCATCH_AQUIFER conflicts with the registered org.hydrocouple.openswmm.integrated2d component: a subcatchment's aquifer has one owner. Use AQUIFER_2D with that component, or remove it."
                                    .to_vec(),
                            ));
                        }
                    }
                    if self.subsurface_.active() {
                        return Err(thrown(
                            b"[2D_OPTIONS] INFIL_DESTINATION SUBCATCH_AQUIFER conflicts with the [2D_AQUIFER] section: infiltration has one owner, and a [2D_AQUIFER] already receives it. Remove the explicit surface row or remove the [2D_AQUIFER]."
                                .to_vec(),
                        ));
                    }
                    let px = &ctx.spatial.subcatch_polygon_x;
                    let py = &ctx.spatial.subcatch_polygon_y;
                    let contains = |sub: usize, x: f64, y: f64| -> bool {
                        if sub >= px.len() || px[sub].len() < 3 {
                            return false;
                        }
                        polygon::contains(&px[sub], &py[sub], x, y)
                    };
                    // The mesh is SI here; the polygons are in the project's map units.
                    let inv = if options.mesh_to_si_factor > 0.0 { 1.0 / options.mesh_to_si_factor } else { 1.0 };
                    self.cell_subcatch_.resize(mesh.n_triangles() as usize, -1);
                    self.subcatch_recharge_.resize(nsub as usize, 0.0);
                    let mut unowned = 0;
                    let mut overlapped = 0;
                    for c in 0..mesh.n_triangles() as usize {
                        let row = &infil.resolved_rows()[c];
                        if !row.has_method || row.dest != Infil2DDest::SubcatchAquifer {
                            continue;
                        }
                        let x = mesh.tri_cx[c] * inv;
                        let y = mesh.tri_cy[c] * inv;
                        let mut owner = -1;
                        let mut hits = 0;
                        for sidx in 0..nsub {
                            if contains(sidx as usize, x, y) {
                                if owner < 0 {
                                    owner = sidx;
                                }
                                hits += 1;
                            }
                        }
                        self.cell_subcatch_[c] = owner;
                        if owner < 0 {
                            unowned += 1;
                        } else if hits > 1 {
                            overlapped += 1;
                        }
                    }
                    if unowned > 0 {
                        ctx.warnings.push(
                            format!(
                                "[2D_OPTIONS] INFIL_DESTINATION SUBCATCH_AQUIFER: {} cell(s) lie outside every subcatchment polygon \u{2014} their infiltration is LOST, as before.",
                                unowned
                            )
                            .into_bytes(),
                        );
                    }
                    if overlapped > 0 {
                        ctx.warnings.push(
                            format!(
                                "[2D_OPTIONS] INFIL_DESTINATION SUBCATCH_AQUIFER: {} cell(s) fall inside more than one subcatchment polygon \u{2014} the first by index receives the recharge.",
                                overlapped
                            )
                            .into_bytes(),
                        );
                    }
                }
            }

            // INFILTRATION NO keeps the rows but runs dry; YES with nothing resolved says so.
            if options.infiltration == 0 && infil.active() {
                ctx.warnings.push(b"[2D_OPTIONS] INFILTRATION NO \xe2\x80\x94 the [2D_INFILTRATION*] rows are kept but no cell infiltrates this run.".to_vec());
                infil.deactivate();
            } else if options.infiltration == 1 && !infil.active() {
                ctx.warnings.push(
                    b"[2D_OPTIONS] INFILTRATION YES but no [2D_INFILTRATION_DEFAULTS] / [2D_INFILTRATION] row resolved to a method \xe2\x80\x94 no cell infiltrates this run.".to_vec(),
                );
            }
        }

        // S1/S2: [2D_INITIAL_QUALITY] needs transported rows, which phase 2a refuses above, so a
        // row here is upstream's own refusal.
        if !self.pending_iq_rows_.borrow().is_empty() {
            if !self.state_.transport.active() {
                return Err(thrown(
                    b"[2D_INITIAL_QUALITY] rows are present but the model carries no transported species ([POLLUTANTS] empty or IGNORE_QUALITY set), so they could seed nothing.".to_vec(),
                ));
            }
            return Err(stub("[2D_INITIAL_QUALITY] seeding (phase 2b)").into());
        }
        // S7: initial buildup (surface quality is refused above, so never active here).
        self.state_.snow_context = true;
        {
            let cfg = &ctx.snowpacks.mesh_config;
            if !ctx.options.ignore_snow_melt && cfg.mode != MeshSnowConfigMode::No && (cfg.assigned() || cfg.mode == MeshSnowConfigMode::Yes) {
                return Err(stub("snow on the 2D mesh (MeshSnow::initialize; phase 2b)").into());
            }
        }
        {
            let e = initialize_mesh_heat(ctx, &mesh, &mut self.state_, &options)?;
            if !e.is_empty() {
                return Err(thrown(cat(&[b"[2D heat] ", &e])));
            }
        }
        // [2D_INITIAL_BED_QUALITY]: needs coupled heat with a LOCAL_STORAGE bed, which needs a
        // surface temperature row (refused above), so a row here is upstream's own refusal.
        if !self.pending_ibq_rows_.borrow().is_empty() {
            let heat = &self.state_.heat;
            if !heat.enabled || !heat.bed || heat.ground_column || heat.bed_pore_volume.is_empty() || heat.bed_species <= 0 {
                return Err(thrown(
                    b"[2D_INITIAL_BED_QUALITY] requires HEAT_MODEL COUPLED and BED_HEAT_MODEL LOCAL_STORAGE with finite bed porewater; GROUND_COLUMN bed species are authored through [GW_INITIAL_QUALITY]."
                        .to_vec(),
                ));
            }
            return Err(stub("[2D_INITIAL_BED_QUALITY] seeding (phase 2b)").into());
        }
        {
            let n_beds = self.subsurface_.node_beds().len();
            let e = configure_ground_lid_owners(ctx, &mesh, self.subsurface_.heat_mut(), n_beds)?;
            if !e.is_empty() {
                return Err(thrown(cat(&[b"[2D LID heat] ", &e])));
            }
        }
        // configureSnowGroundContacts(): mesh snow is never active here (refused above).
        if self.state_.snow.active() {
            return Err(stub("configureSnowGroundContacts (snow on the 2D mesh; phase 2b)").into());
        }
        // The porewater / ground-column / ground-heat checks need an active subsurface heat
        // owner or a ground-column bed, neither of which phase 2a has (see the header).
        if (self.state_.heat.ground_column || self.subsurface_.transport().n_layers > 0 || !ctx.reactions.subsurface_laws.is_empty()) && self.subsurface_.heat().active() {
            return Err(stub("the 2D porewater owner (SubsurfaceSolver::initPorewater; phase 2b)").into());
        }
        if self.state_.heat.ground_column && !self.subsurface_.transport().generated_age.is_empty() {
            return Err(stub("the 2D shared bed contacts (phase 2b)").into());
        }
        if self.subsurface_.heat().active() {
            return Err(stub("2D ground heat (phase 2b)").into());
        }
        // subsurface_.active() is false here (the kernel is refused above).

        self.infil_elapsed_ = 0.0;
        {
            let infil_rc = self.infil_.clone();
            let mut infil = infil_rc.borrow_mut();
            if infil.active() {
                // Publish an initial held rate so the very first substep infiltrates.
                let step = infil.step_seconds();
                let factors = infiltration_factors(ctx, 0.0);
                infil.update_rates(&mesh, &mut self.state_, step, factors, options.num_threads);
                self.infil_cum_applied_.clear();
                self.infil_cum_applied_.resize(mesh.n_triangles() as usize, 0.0);
            } else {
                // Left empty (not zero-filled), as Infil2D::cumulative() is.
                self.infil_cum_applied_.clear();
            }
        }
        self.rain_cum_.clear();
        self.rain_cum_.resize(mesh.n_triangles() as usize, 0.0);
        // C1: always sized.
        self.coupling_cum_applied_.clear();
        self.coupling_cum_applied_.resize(mesh.n_triangles() as usize, 0.0);

        self.subsurface_.prepare_clock();
        self.active_ = true;
        self.sim_time_ = 0.0;
        self.gage_rate_.clear();
        self.gage_cum_.clear();
        self.gage_window_depth_.clear();
        self.gage_window_rate_.clear();
        self.rain_gage_last_.clear();
        self.gage_since_ = 0.0;
        self.gage_report_depth_.clear();
        self.report_rain_cum_.clear();
        self.report_rain_ms_ = -1.0;
        self.pending_dt_ = 0.0;
        self.report_old_depth_.clear();
        self.report_old_depth_.resize(mesh.n_triangles() as usize, 0.0);
        self.report_old_head_.clear();
        self.report_old_head_.resize(mesh.n_triangles() as usize, 0.0);
        self.report_window_old_ms_ = 0.0;
        self.report_window_new_ms_ = 0.0;
        self.report_old_valid_ = false;

        self.reset_window_accumulators();

        // Seed the global 2D mass balance.
        ctx.mass_balance_2d.active = true;
        {
            let mut vol = 0.0;
            for i in 0..mesh.n_triangles() as usize {
                vol += self.state_.volume[i];
            }
            ctx.mass_balance_2d.init_storage = vol;
        }
        ctx.mass_balance_2d.final_storage = ctx.mass_balance_2d.init_storage;
        self.prev_boundary_cum_ = 0.0;
        self.prev_boundary_out_ = 0.0;
        self.prev_boundary_in_ = 0.0;

        // Seed the render-only vertex free surface.
        if options.momentum == Momentum2D::FullSwe {
            if let Some(cpu) = self.solver_.as_ref() {
                cpu.publish_momentum_velocity(&mesh, &mut self.state_, &options);
            }
        } else {
            let b = if self.state_.boundary { Some(&*boundary) } else { None };
            compute_face_velocity(&mesh, &mut self.state_, b, &options);
        }
        reconstruct_vertex_render_depths(&mesh, &mut self.state_, options.dry_depth, options.num_threads, &mut self.render_eta_);
        self.render_dirty_ = false;
        Ok(())
    }

    /// initializeChannelGround(ctx): MESH_CONTACTS channel heat. Requested only by heat on a
    /// link with a mesh-ground sediment; an inactive mesh refuses it (upstream's message), an
    /// active one would initialize the contacts, which is phase 2b.
    pub fn initialize_channel_ground(&mut self, ctx: &mut SimulationContext) -> Result<(), InitError> {
        let mut requested = false;
        if ctx.options.heat_transport && ctx.heat_config.sediment_exchange {
            for link in 0..ctx.n_links() {
                if crate::transport::heat_overrides::sediment_for(ctx, &crate::data::heat_override::HeatElement::link(link)).mesh_ground {
                    requested = true;
                }
            }
        }
        if !requested {
            return Ok(());
        }
        if !self.active_ {
            return Err(thrown(b"MESH_CONTACTS channel heat requires an enabled 2D mesh and porous thermal aquifer.".to_vec()));
        }
        Err(stub("MESH_CONTACTS channel heat on the 2D mesh (ChannelGroundState::initialize; phase 2b)").into())
    }
}
