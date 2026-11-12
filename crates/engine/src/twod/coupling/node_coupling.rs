// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/coupling/NodeCoupling.hpp, NodeCoupling.cpp
// (buildCouplingPoints) and ExchangeStructure.cpp (resolveExchangeStructure)
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The 1D-2D coupling point descriptor and its initialize-time construction: the points from
// the vertex and triangle node maps, and each junction point's exchange law, structure and
// coefficients. The run-time exchange (computeNodeCouplingQ and the rest of NodeCoupling.cpp)
// comes with the step (wave 4).
//
// Differences in form, never in result: CouplingPoint's members without an initializer upstream
// (cell_idx, vertex_idx, node_idx, cd, area, is_outfall, has_flap_gate) default to zero here;
// buildCouplingPoints sets every one of them before a point is kept. `static_cast<int8_t>(v)`
// of an authored LAW value is core::cmath::cvt_i32(v) narrowed to i8, as g++ converts.

#![allow(clippy::all)]

use crate::core::cmath::{cmax, cmin, cvt_i32};
use crate::core::context::SimulationContext;
use crate::core::unit_conversion as ucf;
use crate::data::node::NodeType;
use crate::data::node_coupling_2d::{NodeCoupling2DCoeffs, NodeCoupling2DCrest, NodeCoupling2DKey, NodeCoupling2DParams, NodeCoupling2DType};
use crate::hydraulics::node::get_surf_area;
use crate::twod::coupling::exchange_law::{Coefficients, Geometry};
use crate::twod::data::mesh_data::MeshData;
use crate::twod::data::solver_options_2d::{CouplingLaw2D, SolverOptions2D};

#[derive(Debug, Clone)]
pub struct CouplingPoint {
    pub cell_idx: i32,
    pub vertex_idx: i32,
    pub node_idx: i32,
    pub cd: f64,
    pub area: f64,
    pub is_outfall: bool,
    pub has_flap_gate: bool,
    pub area_authored: bool,
    pub cd_authored: bool,
    pub law: i8,
    pub geo: Geometry,
    pub coef: Coefficients,
    pub chamber_area: f64,
    pub equalize: bool,
    pub head_begin: i32,
    pub head_end: i32,
    pub spread: f64,
    pub tw_relax: f64,
}

impl Default for CouplingPoint {
    fn default() -> Self {
        CouplingPoint {
            cell_idx: 0,
            vertex_idx: 0,
            node_idx: 0,
            cd: 0.0,
            area: 0.0,
            is_outfall: false,
            has_flap_gate: false,
            area_authored: true,
            cd_authored: true,
            law: 0,
            geo: Geometry::default(),
            coef: Coefficients::default(),
            chamber_area: 0.0,
            equalize: false,
            head_begin: 0,
            head_end: 0,
            spread: 0.0,
            tw_relax: 0.0,
        }
    }
}

/// buildCouplingPoints(mesh, ctx): one point per coupled vertex, then one per triangle row.
pub fn build_coupling_points(mesh: &MeshData, ctx: &SimulationContext) -> Vec<CouplingPoint> {
    let mut cps: Vec<CouplingPoint> = Vec::new();

    // Vertex-to-node couplings: the first (lowest-index) cell touching each vertex.
    let nv = mesh.n_vertices();
    let mut first_cell_of_vertex = vec![-1i32; nv as usize];
    {
        let nt_all = mesh.n_triangles();
        for t in 0..nt_all {
            let nvc = mesh.cell_vertex_count(t);
            for k in 0..nvc {
                let vv = mesh.cell_vertex(t, k);
                if vv >= 0 && vv < nv && first_cell_of_vertex[vv as usize] < 0 {
                    first_cell_of_vertex[vv as usize] = t;
                }
            }
        }
    }
    for v in 0..nv {
        let uv = v as usize;
        let node_idx = mesh.vert_coupled_node[uv];
        if node_idx < 0 {
            continue;
        }
        let mut cp = CouplingPoint::default();
        cp.vertex_idx = v;
        cp.node_idx = node_idx;
        cp.cd = mesh.vert_coupling_cd[uv];
        cp.area = mesh.vert_coupling_area[uv];
        cp.area_authored = uv < mesh.vert_coupling_area_set.len() && mesh.vert_coupling_area_set[uv] != 0;
        cp.cd_authored = uv < mesh.vert_coupling_cd_set.len() && mesh.vert_coupling_cd_set[uv] != 0;

        cp.cell_idx = first_cell_of_vertex[uv];
        if cp.cell_idx < 0 {
            continue;
        }

        let ui = node_idx as usize;
        cp.is_outfall = ctx.nodes.type_[ui] == NodeType::Outfall;
        let ofr = ctx.node_subtypes.outfall_row(node_idx);
        cp.has_flap_gate = cp.is_outfall && ofr >= 0 && ctx.node_subtypes.outfalls.has_flap_gate[ofr as usize] != 0;

        cps.push(cp);
    }

    // Triangle-to-node couplings: one point per authored row.
    for row in &mesh.tri_couplings {
        let node_idx = row.node;
        if node_idx < 0 {
            continue;
        }
        if row.tri < 0 || row.tri >= mesh.n_triangles() {
            continue;
        }
        let mut cp = CouplingPoint::default();
        cp.cell_idx = row.tri;
        cp.vertex_idx = -1; // centroid coupling
        cp.node_idx = node_idx;
        cp.cd = row.cd;
        cp.area = row.area;
        cp.area_authored = row.area_set;
        cp.cd_authored = row.cd_set;

        let ui = node_idx as usize;
        cp.is_outfall = ctx.nodes.type_[ui] == NodeType::Outfall;
        let ofr = ctx.node_subtypes.outfall_row(node_idx);
        cp.has_flap_gate = cp.is_outfall && ofr >= 0 && ctx.node_subtypes.outfalls.has_flap_gate[ofr as usize] != 0;

        cps.push(cp);
    }
    cps
}

/// resolveExchangeStructure(cp, mesh, ctx, opts, row, type, project_len_to_m): a junction
/// point's exchange law, its WEIR_ORIFICE structure (SI) and coefficients.
pub fn resolve_exchange_structure(
    cp: &mut CouplingPoint,
    mesh: &MeshData,
    ctx: &SimulationContext,
    opts: &SolverOptions2D,
    row: i32,
    type_: NodeCoupling2DType,
    project_len_to_m: f64,
) {
    let ni = cp.node_idx as usize;
    let p: Option<&NodeCoupling2DParams> = if row >= 0 { Some(&ctx.node_coupling_2d.params[row as usize]) } else { None };
    let authored = |k: NodeCoupling2DKey, out: &mut f64| -> bool {
        let Some(p) = p else {
            return false;
        };
        let v = p[k as usize];
        if v.is_nan() {
            return false;
        }
        *out = v;
        true
    };
    let mut v = 0.0;
    let f = project_len_to_m;
    let f2 = f * f;
    cp.law = if authored(NodeCoupling2DKey::Law, &mut v) { cvt_i32(v) as i8 } else { opts.coupling_law as i8 };
    let weir_orifice = cp.law == CouplingLaw2D::WeirOrifice as i8;
    cp.equalize = if opts.coupling_equalize >= 0 { opts.coupling_equalize != 0 } else { weir_orifice };
    let chamber_node = authored(NodeCoupling2DKey::ChamberArea, &mut v);
    let chamber_val = v * f2;
    cp.chamber_area = if chamber_node { chamber_val } else { cp.area };
    if !weir_orifice {
        return;
    }

    // ---- Coefficients ----
    let coeffs_authored = authored(NodeCoupling2DKey::Coeffs, &mut v);
    let preset = if coeffs_authored { cvt_i32(v) } else { opts.coupling_coeffs };
    let mut c = Coefficients::default();
    if preset == NodeCoupling2DCoeffs::Rubinato2017 as i32 {
        c.c_w = 0.540;
        c.c_sw = 0.056;
        c.c_oi = 0.167;
        c.c_oo = 0.160;
    } else if preset == NodeCoupling2DCoeffs::Rubinato2017Lower as i32 {
        c.c_w = 0.493;
        c.c_sw = 0.055;
        c.c_oi = 0.160;
        c.c_oo = 0.159;
    } else if preset == NodeCoupling2DCoeffs::Rubinato2017Upper as i32 {
        c.c_w = 0.587;
        c.c_sw = 0.057;
        c.c_oi = 0.174;
        c.c_oo = 0.161;
    } else if preset == NodeCoupling2DCoeffs::Default as i32 {
        c.c_w = 0.540;
        c.c_oi = 0.167;
        c.c_oo = 0.160;
    }
    let open_pond = type_ == NodeCoupling2DType::Open && !coeffs_authored;
    if open_pond {
        c.c_oo = 1.0; // open pond: the heads track
        c.c_oi = 1.0;
        c.c_w = 1.0;
    }
    if authored(NodeCoupling2DKey::Cw, &mut v) {
        c.c_w = v;
    }
    // DEFAULT (and the open pond) keep the two weir forms continuous at h1 = crest:
    // C_sw = 2/3 C_w.
    if preset == NodeCoupling2DCoeffs::Default as i32 || open_pond {
        c.c_sw = (2.0 / 3.0) * c.c_w;
    }
    if authored(NodeCoupling2DKey::Csw, &mut v) {
        c.c_sw = v;
    }
    if authored(NodeCoupling2DKey::Coi, &mut v) {
        c.c_oi = v;
    }
    if authored(NodeCoupling2DKey::Coo, &mut v) {
        c.c_oo = v;
    }
    if cp.cd_authored {
        c.c_oo = cp.cd; // legacy map-row CD
        c.c_oi = c.c_oo;
    }

    // ---- Geometry (SI) ----
    const K_PI: f64 = 3.14159265358979323846;
    // Default manhole: 1.2 m, or 4 ft in US units.
    let d_default = if f == 1.0 { 1.2 } else { 4.0 * f };
    let d = if authored(NodeCoupling2DKey::D, &mut v) { v * f } else { d_default };
    let mut area = K_PI * d * d / 4.0;
    let mut perim = K_PI * d;
    let area_node = authored(NodeCoupling2DKey::Area, &mut v);
    let area_val = v * f2;
    let mut perim_val = 0.0;
    let perim_node = authored(NodeCoupling2DKey::Perimeter, &mut perim_val);
    perim_val *= f;
    match type_ {
        NodeCoupling2DType::Grate | NodeCoupling2DType::Curb => {
            if area_node {
                area = area_val;
            } else if cp.area_authored {
                area = cp.area;
            }
            perim = if perim_node { perim_val } else { 4.0 * area.sqrt() };
        }
        NodeCoupling2DType::Open => {
            if area_node {
                area = area_val;
            } else {
                let us = ucf::get_unit_system(ctx.options.flow_units as i32);
                let a_ft2 = get_surf_area(&ctx.nodes, cp.node_idx, ctx.nodes.full_depth[ni], Some(&ctx.tables), us, Some(&ctx.node_subtypes));
                if a_ft2 > 0.0 {
                    area = a_ft2 * opts.len_1d_to_2d * opts.len_1d_to_2d;
                }
            }
            perim = if perim_node { perim_val } else { (4.0 * K_PI * area).sqrt() };
        }
        _ => {
            // MANHOLE
            if area_node {
                area = area_val;
            }
            if perim_node {
                perim = perim_val;
            }
        }
    }
    if cp.area_authored {
        area = cp.area; // legacy map-row AREA
    }
    cp.chamber_area = if chamber_node { chamber_val } else { area };
    let lid = if authored(NodeCoupling2DKey::LidOpenFraction, &mut v) { v } else { 1.0 };
    cp.geo.a_orif = area * lid;
    cp.geo.l_weir = perim * lid;

    // ---- Crest ----
    let rim = (ctx.nodes.invert_elev[ni] + ctx.nodes.full_depth[ni]) * opts.len_1d_to_2d;
    let bed = if cp.vertex_idx >= 0 { mesh.vz[cp.vertex_idx as usize] } else { mesh.tri_cz[cp.cell_idx as usize] };
    let crest = if authored(NodeCoupling2DKey::Crest, &mut v) { cvt_i32(v) } else { 0 };
    if crest == NodeCoupling2DCrest::MeshBed as i32 {
        cp.geo.z_crest = bed;
    } else if crest == NodeCoupling2DCrest::Max as i32 {
        cp.geo.z_crest = cmax(rim, bed);
    } else if crest == NodeCoupling2DCrest::Min as i32 {
        cp.geo.z_crest = cmin(rim, bed);
    } else if crest == NodeCoupling2DCrest::NodeRim as i32 {
        cp.geo.z_crest = rim;
    }
    cp.coef = c;
}
