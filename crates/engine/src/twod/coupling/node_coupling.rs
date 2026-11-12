// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/coupling/NodeCoupling.hpp, NodeCoupling.cpp,
// ExchangeStructure.cpp (resolveExchangeStructure) and the vertexHeadAt helper of
// src/engine/2d/mesh/VertexReconstruction.hpp (HydroCouple/openswmm.engine, branch swmm6_rel,
// commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The 1D-2D coupling point descriptor, its initialize-time construction (the points from the
// vertex and triangle node maps, and each junction point's exchange law, structure and
// coefficients: wave 3), and the run-time exchange (wave 4b): the live junction exchange
// computeNodeCouplingQ (ORIFICE with the capped-pipe gate and the area ramp, or WEIR_ORIFICE;
// the CELL / VERTEX / STENCIL_WET_MEAN coupling heads; the source-side wet ramp; the
// COUPLING_EQUALIZE cap), its head sensitivity computeNodeCouplingDQdh1d, the outfall
// tailwater cache updateOutfallBoundaries, and the window helpers scatterCouplingToYdot,
// accumulateOutfallDischargeStep and injectAccumulatedExchange (the last three have no caller
// upstream; they are translated with the rest of the file and checked on real states).
//
// Differences in form, never in result:
//   * CouplingPoint's members without an initializer upstream (cell_idx, vertex_idx, node_idx,
//     cd, area, is_outfall, has_flap_gate) default to zero here; buildCouplingPoints sets every
//     one of them before a point is kept.
//   * `static_cast<int8_t>(v)` of an authored LAW value is core::cmath::cvt_i32(v) narrowed to
//     i8, as g++ converts.
//   * SurfaceStateData's two pointers the exchange reads (coupling_head_cells -> the router's
//     coupling_head_cells_, node_pond_area_m2 -> the router's node_pond_area_m2_) are flags in
//     the translated state; the targets are passed in a CouplingAttachments, None where
//     upstream's pointer is null (as wave 3's SolverAttachments does for the solver).
//   * Output pointers (ydot, sample_row, provisional_vol_m3, node_conc) are slices, None for a
//     null pointer.
#![allow(clippy::all)]

use crate::core::cmath::{cclamp, cmax, cmin, cvt_i32};
use crate::core::context::SimulationContext;
use crate::core::unit_conversion as ucf;
use crate::data::node::{NodeData, NodeType};
use crate::data::node_coupling_2d::{NodeCoupling2DCoeffs, NodeCoupling2DCrest, NodeCoupling2DKey, NodeCoupling2DParams, NodeCoupling2DType};
use crate::hydraulics::node::get_surf_area;
use crate::twod::coupling::exchange_law::{self, Coefficients, Geometry};
use crate::twod::coupling::pond_model as pond;
use crate::twod::data::mesh_data::MeshData;
use crate::twod::data::solver_options_2d::{CellClosure2D, CouplingLaw2D, SolverOptions2D};
use crate::twod::data::surface_state_data::SurfaceStateData;

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

// ---------------------------------------------------------------------------------------------
// The run-time exchange (NodeCoupling.cpp, wave 4b)
// ---------------------------------------------------------------------------------------------

/// The targets of the two SurfaceStateData pointers the exchange reads (upstream's
/// `state.coupling_head_cells` and `state.node_pond_area_m2`): None where the pointer is null.
/// The router attaches them (SurfaceRouter2D::initialize sets the first, publishEqualizeAreas
/// the second); a caller builds this from the state's flags and the router's members with
/// `CouplingAttachments::of`.
#[derive(Debug, Clone, Copy, Default)]
pub struct CouplingAttachments<'a> {
    /// *state.coupling_head_cells (the router's coupling_head_cells_), when attached.
    pub coupling_head_cells: Option<&'a [i32]>,
    /// *state.node_pond_area_m2 (the router's node_pond_area_m2_), when attached.
    pub node_pond_area_m2: Option<&'a [f64]>,
}

impl<'a> CouplingAttachments<'a> {
    /// The attachments of `state` given the router's two target members.
    pub fn of(state: &SurfaceStateData, coupling_head_cells: &'a [i32], node_pond_area_m2: &'a [f64]) -> Self {
        CouplingAttachments {
            coupling_head_cells: if state.coupling_head_cells { Some(coupling_head_cells) } else { None },
            node_pond_area_m2: if state.node_pond_area_m2 { Some(node_pond_area_m2) } else { None },
        }
    }
}

/// m/s2
const GRAVITY: f64 = 9.80665;

/// Regularization head (m) below which the orifice root law is replaced by a C1 quadratic.
const ORIFICE_H_EPS: f64 = 0.02;

/// vertexHeadAt (VertexReconstruction.hpp): the pseudo-Laplacian head at one vertex.
#[inline]
pub fn vertex_head_at(mesh: &MeshData, state: &SurfaceStateData, v: i32) -> f64 {
    let start = mesh.vert_stencil_ptr[v as usize];
    let end = mesh.vert_stencil_ptr[v as usize + 1];
    let mut h = 0.0;
    for k in start..end {
        let uk = k as usize;
        h += mesh.vert_stencil_wt[uk] * state.head[mesh.vert_stencil_idx[uk] as usize];
    }
    h
}

/// orificePhi: sqrt(a) for a >= eps, a C1 quadratic below.
#[inline]
fn orifice_phi(a: f64) -> f64 {
    if a >= ORIFICE_H_EPS {
        return a.sqrt();
    }
    let inv_sqrt_e = 1.0 / ORIFICE_H_EPS.sqrt();
    (1.5 * inv_sqrt_e) * a - (0.5 * inv_sqrt_e / ORIFICE_H_EPS) * a * a
}

/// orificeFlow: Q = Cd A sign(dh) sqrt(2g) phi(|dh|).
#[inline]
fn orifice_flow(dh: f64, cd: f64, area: f64) -> f64 {
    let a = dh.abs();
    if a < 1.0e-12 {
        return 0.0;
    }
    let sign = if dh > 0.0 { 1.0 } else { -1.0 };
    sign * cd * area * (2.0 * GRAVITY).sqrt() * orifice_phi(a)
}

/// effectiveArea: the inlet area below the crown, ramped to the manhole area over 5 cm above
/// it (full_depth is unused upstream too).
#[inline]
fn effective_area(h_max: f64, z_ground: f64, _full_depth: f64, a_inlet: f64, a_manhole: f64) -> f64 {
    if h_max < z_ground {
        return a_inlet;
    }
    let d_trans = 0.05; // 5 cm transition
    let frac = cmin((h_max - z_ground) / d_trans, 1.0);
    a_inlet + frac * (a_manhole - a_inlet)
}

/// scatterCouplingFlux: a signed exchange Q (m3/s) onto the point's cell(s) as a coupling_flux
/// rate (m/s), upwind-HGL weighted over a vertex stencil (geometric weights on a flat surface);
/// with `conc` and Q > 0 the species density rides along in transport.coupling_src.
fn scatter_coupling_flux(mesh: &MeshData, state: &mut SurfaceStateData, cp: &CouplingPoint, q: f64, conc: Option<&[f64]>) {
    let with_species = conc.is_some() && q > 0.0 && state.transport.active() && !state.transport.coupling_src.is_empty();
    let add = |state: &mut SurfaceStateData, cell: i32, q_share: f64| {
        let area = mesh.tri_area[cell as usize];
        if !(area > 1.0e-30) {
            return;
        }
        let dens = q_share / area;
        state.coupling_flux[cell as usize] += dens;
        if with_species {
            let conc = conc.unwrap();
            let tr = &mut state.transport;
            for s in 0..tr.n_species {
                let ix = tr.idx(s, cell);
                tr.coupling_src[ix] += dens * conc[s as usize];
            }
        }
    };
    // Triangle coupling: single cell, head and flux already co-located.
    if cp.vertex_idx < 0 {
        add(state, cp.cell_idx, q);
        return;
    }

    let v = cp.vertex_idx;
    let start = mesh.vert_stencil_ptr[v as usize];
    let end = mesh.vert_stencil_ptr[v as usize + 1];
    let hv = vertex_head_at(mesh, state, v);
    let vx = mesh.vx[v as usize];
    let vy = mesh.vy[v as usize];
    let sign = if q >= 0.0 { 1.0 } else { -1.0 }; // source -> downhill, sink -> uphill

    // Pass 1: sum the upwind-slope weights over the stencil.
    let mut wsum = 0.0;
    for k in start..end {
        let kc = mesh.vert_stencil_idx[k as usize];
        let dx = mesh.tri_cx[kc as usize] - vx;
        let dy = mesh.tri_cy[kc as usize] - vy;
        let d = (dx * dx + dy * dy).sqrt();
        if d < 1.0e-9 {
            continue;
        }
        let slope = sign * (hv - state.head[kc as usize]) / d;
        if slope > 0.0 {
            wsum += slope;
        }
    }

    // Pass 2: scatter Q, upwind when a usable gradient exists, else the geometric
    // partition-of-unity fallback (flat surface).
    if wsum > 1.0e-30 {
        for k in start..end {
            let kc = mesh.vert_stencil_idx[k as usize];
            let dx = mesh.tri_cx[kc as usize] - vx;
            let dy = mesh.tri_cy[kc as usize] - vy;
            let d = (dx * dx + dy * dy).sqrt();
            if d < 1.0e-9 {
                continue;
            }
            let slope = sign * (hv - state.head[kc as usize]) / d;
            if slope <= 0.0 {
                continue;
            }
            add(state, kc, q * (slope / wsum));
        }
    } else {
        for k in start..end {
            let kc = mesh.vert_stencil_idx[k as usize];
            add(state, kc, q * mesh.vert_stencil_wt[k as usize]);
        }
    }
}

/// wetVertexEta: the depth-weighted head of the vertex stencil's wet cells (depth >=
/// dry_depth), the vertex ground elevation when all are dry.
#[inline]
fn wet_vertex_eta(mesh: &MeshData, state: &SurfaceStateData, v: i32, dry_depth: f64) -> f64 {
    let s = mesh.vert_stencil_ptr[v as usize];
    let e = mesh.vert_stencil_ptr[v as usize + 1];
    let mut num = 0.0;
    let mut den = 0.0;
    for k in s..e {
        let t = mesh.vert_stencil_idx[k as usize] as usize;
        let h = state.depth[t];
        if !(h >= dry_depth) {
            continue; // dry skip (NaN-robust)
        }
        num += h * state.head[t];
        den += h;
    }
    if den > 0.0 {
        num / den
    } else {
        mesh.vz[v as usize]
    }
}

/// couplingHead2D: the point's sampling-cell wet-depth mean (COUPLING_HEAD VERTEX /
/// STENCIL_WET_MEAN), else the cell head (centroid coupling), the wet-masked vertex head (VFR)
/// or the pseudo-Laplacian vertex head (FLAT).
#[inline]
fn coupling_head_2d(cp: &CouplingPoint, mesh: &MeshData, state: &SurfaceStateData, att: &CouplingAttachments, opts: &SolverOptions2D) -> f64 {
    if cp.head_end > cp.head_begin {
        if let Some(cells) = att.coupling_head_cells {
            let mut num = 0.0;
            let mut den = 0.0;
            for k in cp.head_begin..cp.head_end {
                let t = cells[k as usize] as usize;
                let h = state.depth[t];
                if !(h >= opts.dry_depth) {
                    continue;
                }
                num += h * state.head[t];
                den += h;
            }
            return if den > 0.0 { num / den } else { state.head[cp.cell_idx as usize] };
        }
    }
    if cp.vertex_idx < 0 {
        return state.head[cp.cell_idx as usize];
    }
    if opts.cell_closure == CellClosure2D::Vfr {
        return wet_vertex_eta(mesh, state, cp.vertex_idx, opts.dry_depth);
    }
    vertex_head_at(mesh, state, cp.vertex_idx)
}

/// sourceWetRamp: the Hermite smoothstep of clamp(d / dry_depth, 0, 1).
#[inline]
fn source_wet_ramp(d: f64, dry_depth: f64) -> f64 {
    let t = cmin(1.0, cmax(0.0, d / dry_depth));
    t * t * (3.0 - 2.0 * t)
}

/// computeNodeCouplingQ(cp, mesh, state, nodes, opts, provisional_vol_m3, h1d_offset_m): the
/// live exchange of one non-outfall point, m3/s (> 0 drains 2D -> 1D, < 0 spills 1D -> 2D).
/// `provisional_vol_m3` (upstream nullable, None on the live path) re-derives the source
/// depths and heads from per-cell volumes under the flat closure.
pub fn compute_node_coupling_q(
    cp: &CouplingPoint,
    mesh: &MeshData,
    state: &SurfaceStateData,
    att: &CouplingAttachments,
    nodes: &NodeData,
    opts: &SolverOptions2D,
    provisional_vol_m3: Option<&[f64]>,
    h1d_offset_m: f64,
) -> f64 {
    let ni = cp.node_idx as usize;
    let ci = cp.cell_idx;

    let prov_depth = |i: i32| -> f64 {
        let a = mesh.tri_area[i as usize];
        let v = cmax(0.0, provisional_vol_m3.unwrap()[i as usize]);
        if a > 1.0e-30 {
            v / a
        } else {
            0.0
        }
    };

    // Live 2D head at the coupling point.
    let h_2d;
    if provisional_vol_m3.is_none() {
        h_2d = coupling_head_2d(cp, mesh, state, att, opts);
    } else if cp.vertex_idx < 0 {
        h_2d = mesh.tri_cz[ci as usize] + prov_depth(ci);
    } else {
        // Wet-depth-weighted eta over the stencil, on provisional depths.
        let v = cp.vertex_idx;
        let s = mesh.vert_stencil_ptr[v as usize];
        let e = mesh.vert_stencil_ptr[v as usize + 1];
        let mut num = 0.0;
        let mut den = 0.0;
        for k in s..e {
            let t = mesh.vert_stencil_idx[k as usize];
            let d = prov_depth(t);
            if !(d >= opts.dry_depth) {
                continue;
            }
            num += d * (mesh.tri_cz[t as usize] + d);
            den += d;
        }
        h_2d = if den > 0.0 { num / den } else { mesh.vz[v as usize] };
    }
    // 1D node head, frozen across the 2D advance window.
    let h_1d = nodes.head[ni] * opts.len_1d_to_2d + h1d_offset_m;

    // COUPLING_EQUALIZE: no single coupling window may carry the 1D head past the 2D head.
    let equalized = |q: f64| -> f64 {
        let Some(pond_area) = att.node_pond_area_m2 else {
            return q;
        };
        if !cp.equalize || !(state.coupling_window_s > 0.0) {
            return q;
        }
        let cap = pond::equalize_cap(h_2d - h_1d, state.coupling_window_s, pond_area[ni]);
        cclamp(q, -cap, cap)
    };

    // Live max-over-stencil source depth (the wet/dry self-limiter).
    let depth_2d_avail;
    if cp.vertex_idx >= 0 {
        let v = cp.vertex_idx;
        let s = mesh.vert_stencil_ptr[v as usize];
        let e = mesh.vert_stencil_ptr[v as usize + 1];
        let mut d_avail = 0.0;
        for k in s..e {
            let t = mesh.vert_stencil_idx[k as usize];
            d_avail = cmax(d_avail, if provisional_vol_m3.is_none() { state.depth[t as usize] } else { prov_depth(t) });
        }
        depth_2d_avail = d_avail;
    } else {
        depth_2d_avail = if provisional_vol_m3.is_none() { state.depth[ci as usize] } else { prov_depth(ci) };
    }
    let depth_1d_avail = nodes.depth[ni] * opts.len_1d_to_2d;

    // WEIR_ORIFICE: the regime-aware law on the point's structure; the source-side wet ramp
    // still applies.
    if cp.law == CouplingLaw2D::WeirOrifice as i8 {
        let q = exchange_law::exchange(&cp.geo, &cp.coef, h_2d, h_1d).q;
        return equalized(q * source_wet_ramp(if q > 0.0 { depth_2d_avail } else { depth_1d_avail }, opts.dry_depth));
    }

    // Capped-pipe gate at the crown (= surcharge / slot-engagement threshold).
    let crown = (nodes.invert_elev[ni] + nodes.full_depth[ni]) * opts.len_1d_to_2d;
    let z_top = crown;
    let h_max = cmax(h_1d, h_2d);
    let a_eff = effective_area(h_max, crown, nodes.full_depth[ni], cp.area, cp.area * 2.0);

    let mut q = orifice_flow(h_2d - h_1d, cp.cd, a_eff);

    // C1 Hermite gate over a 5 cm band above the crown.
    const CAP_BAND: f64 = 0.05;
    let ct = cmin(1.0, cmax(0.0, (h_max - z_top) / CAP_BAND));
    q *= ct * ct * (3.0 - 2.0 * ct);

    // Source-side wet/dry self-limit on the live depth.
    let wet_ramp = |d: f64| -> f64 {
        let t = cmin(1.0, cmax(0.0, d / opts.dry_depth));
        t * t * (3.0 - 2.0 * t)
    };
    q *= if q > 0.0 { wet_ramp(depth_2d_avail) } else { wet_ramp(depth_1d_avail) };
    equalized(q)
}

/// computeNodeCouplingDQdh1d(cp, mesh, state, nodes, opts): the head sensitivity G = -dQ/dh_1d
/// >= 0 (gate, ramp and cap derivatives dropped), m3/s per m of SI 1D head.
pub fn compute_node_coupling_dqdh1d(cp: &CouplingPoint, mesh: &MeshData, state: &SurfaceStateData, att: &CouplingAttachments, nodes: &NodeData, opts: &SolverOptions2D) -> f64 {
    let ni = cp.node_idx as usize;
    let h_2d = coupling_head_2d(cp, mesh, state, att, opts);
    let h_1d = nodes.head[ni] * opts.len_1d_to_2d;

    if cp.law == CouplingLaw2D::WeirOrifice as i8 {
        let depth_1d = nodes.depth[ni] * opts.len_1d_to_2d;
        let depth_2d = if cp.cell_idx >= 0 { state.depth[cp.cell_idx as usize] } else { 0.0 };
        let g = exchange_law::exchange(&cp.geo, &cp.coef, h_2d, h_1d).g * source_wet_ramp(cmax(depth_1d, depth_2d), opts.dry_depth);
        return if g > 0.0 { g } else { 0.0 };
    }

    let a = (h_2d - h_1d).abs();

    let crown = (nodes.invert_elev[ni] + nodes.full_depth[ni]) * opts.len_1d_to_2d;
    let h_max = cmax(h_1d, h_2d);
    let a_eff = effective_area(h_max, crown, nodes.full_depth[ni], cp.area, cp.area * 2.0);

    // phi' of the C1-regularized root.
    let phi_p;
    if a >= ORIFICE_H_EPS {
        phi_p = 0.5 / a.sqrt();
    } else {
        let inv_sqrt_e = 1.0 / ORIFICE_H_EPS.sqrt();
        phi_p = 1.5 * inv_sqrt_e - (inv_sqrt_e / ORIFICE_H_EPS) * a;
    }
    let mut g = cp.cd * a_eff * (2.0 * GRAVITY).sqrt() * phi_p;

    const CAP_BAND: f64 = 0.05;
    let ct = cmin(1.0, cmax(0.0, (h_max - crown) / CAP_BAND));
    g *= ct * ct * (3.0 - 2.0 * ct);

    let wet_ramp = |d: f64| -> f64 {
        let t = cmin(1.0, cmax(0.0, d / opts.dry_depth));
        t * t * (3.0 - 2.0 * t)
    };
    let depth_1d = nodes.depth[ni] * opts.len_1d_to_2d;
    let depth_2d = if cp.cell_idx >= 0 { state.depth[cp.cell_idx as usize] } else { 0.0 };
    g *= wet_ramp(cmax(depth_1d, depth_2d));
    if g > 0.0 {
        g
    } else {
        0.0
    }
}

/// scatterCouplingToYdot(mesh, state, cp, Q, ydot): a signed exchange (m3/s, > 0 into the
/// cells) added onto the cell derivatives, upwind-HGL weighted over a vertex stencil.
pub fn scatter_coupling_to_ydot(mesh: &MeshData, state: &SurfaceStateData, cp: &CouplingPoint, q: f64, ydot: &mut [f64]) {
    if cp.vertex_idx < 0 {
        ydot[cp.cell_idx as usize] += q;
        return;
    }

    let v = cp.vertex_idx;
    let start = mesh.vert_stencil_ptr[v as usize];
    let end = mesh.vert_stencil_ptr[v as usize + 1];
    let hv = vertex_head_at(mesh, state, v);
    let vx = mesh.vx[v as usize];
    let vy = mesh.vy[v as usize];
    let sign = if q >= 0.0 { 1.0 } else { -1.0 }; // source -> downhill, sink -> uphill

    let mut wsum = 0.0;
    for k in start..end {
        let kc = mesh.vert_stencil_idx[k as usize] as usize;
        let dx = mesh.tri_cx[kc] - vx;
        let dy = mesh.tri_cy[kc] - vy;
        let d = (dx * dx + dy * dy).sqrt();
        if d < 1.0e-9 {
            continue;
        }
        let slope = sign * (hv - state.head[kc]) / d;
        if slope > 0.0 {
            wsum += slope;
        }
    }

    if wsum > 1.0e-30 {
        for k in start..end {
            let kc = mesh.vert_stencil_idx[k as usize] as usize;
            let dx = mesh.tri_cx[kc] - vx;
            let dy = mesh.tri_cy[kc] - vy;
            let d = (dx * dx + dy * dy).sqrt();
            if d < 1.0e-9 {
                continue;
            }
            let slope = sign * (hv - state.head[kc]) / d;
            if slope <= 0.0 {
                continue;
            }
            ydot[kc] += q * (slope / wsum);
        }
    } else {
        for k in start..end {
            let kc = mesh.vert_stencil_idx[k as usize] as usize;
            ydot[kc] += q * mesh.vert_stencil_wt[k as usize];
        }
    }
}

/// updateOutfallBoundaries(cps, mesh, state, ctx, opts, dt): cache each coupled outfall's 2D
/// stage (1D length) and Hermite wet ramp in the outfall rows (head_2d, ramp_2d), with
/// TAILWATER_RELAX. Several points on one outfall: the last one's stage, ramp and relaxation.
pub fn update_outfall_boundaries(cps: &[CouplingPoint], mesh: &MeshData, state: &SurfaceStateData, att: &CouplingAttachments, ctx: &mut SimulationContext, opts: &SolverOptions2D, dt: f64) {
    let wet_ramp = |d: f64| -> f64 {
        let t = cmin(1.0, cmax(0.0, (d - opts.dry_depth) / opts.dry_depth));
        t * t * (3.0 - 2.0 * t) // smoothstep, C1 continuous
    };
    let mut target: Vec<f64> = Vec::new(); // per outfall row, 1D length; NaN = unmapped
    let mut relax: Vec<f64> = vec![0.0; ctx.node_subtypes.outfalls.head_2d.len()];
    for cp in cps {
        if !cp.is_outfall {
            continue;
        }

        // 2D stage and wetness at the outfall coupling point: the deepest stencil cell for a
        // vertex point (wetness and stage), the cell itself for a centroid point.
        let mut h_2d;
        let depth_2d;
        if cp.vertex_idx >= 0 {
            let v = cp.vertex_idx;
            let s = mesh.vert_stencil_ptr[v as usize];
            let e = mesh.vert_stencil_ptr[v as usize + 1];
            let mut deepest = -1;
            let mut dmax = 0.0;
            for k in s..e {
                let c = mesh.vert_stencil_idx[k as usize];
                if state.depth[c as usize] > dmax {
                    dmax = state.depth[c as usize];
                    deepest = c;
                }
            }
            depth_2d = dmax;
            h_2d = if deepest >= 0 { state.head[deepest as usize] } else { mesh.vz[v as usize] };
        } else {
            h_2d = state.head[cp.cell_idx as usize];
            depth_2d = state.depth[cp.cell_idx as usize];
        }
        // COUPLING_HEAD VERTEX / STENCIL_WET_MEAN: the wet-depth-weighted eta over the
        // sampling cells; wetness still from the deepest.
        if cp.head_end > cp.head_begin {
            if let Some(cells) = att.coupling_head_cells {
                let mut num = 0.0;
                let mut den = 0.0;
                for k in cp.head_begin..cp.head_end {
                    let t = cells[k as usize] as usize;
                    if !(state.depth[t] >= opts.dry_depth) {
                        continue;
                    }
                    num += state.depth[t] * state.head[t];
                    den += state.depth[t];
                }
                if den > 0.0 {
                    h_2d = num / den;
                }
            }
        }
        let orow = ctx.node_subtypes.outfall_row(cp.node_idx);
        if orow >= 0 {
            let ur = orow as usize;
            if target.is_empty() {
                target.clear();
                target.resize(ctx.node_subtypes.outfalls.head_2d.len(), f64::NAN);
            }
            target[ur] = h_2d * opts.len_2d_to_1d; // several points: the last wins
            relax[ur] = cp.tw_relax;
            ctx.node_subtypes.outfalls.ramp_2d[ur] = wet_ramp(depth_2d);
        }
    }
    // TAILWATER_RELAX tau: h_tw <- h_tw + (dt/tau)(h_2d - h_tw), once per outfall; the first
    // call (still the -1e30 sentinel) starts it at the surface.
    for ur in 0..target.len() {
        if target[ur].is_nan() {
            continue;
        }
        let h_tw = &mut ctx.node_subtypes.outfalls.head_2d[ur];
        if relax[ur] > 0.0 && dt > 0.0 && *h_tw > -1.0e29 {
            *h_tw += cmin(1.0, dt / relax[ur]) * (target[ur] - *h_tw);
        } else {
            *h_tw = target[ur];
        }
    }
}

/// budgetAvail: the remaining withdrawal budget (m3) a point sees.
#[inline]
fn budget_avail(mesh: &MeshData, cp: &CouplingPoint, cell_budget_m3: &[f64]) -> f64 {
    if cp.vertex_idx < 0 {
        return cmax(0.0, cell_budget_m3[cp.cell_idx as usize]);
    }
    let mut avail = 0.0;
    let v = cp.vertex_idx as usize;
    for k in mesh.vert_stencil_ptr[v]..mesh.vert_stencil_ptr[v + 1] {
        avail += cmax(0.0, cell_budget_m3[mesh.vert_stencil_idx[k as usize] as usize]);
    }
    avail
}

/// budgetDraw: draw `vol` (m3) down from the point's budget cells, weighted then greedy,
/// each cell floored at 0.
#[inline]
fn budget_draw(mesh: &MeshData, cp: &CouplingPoint, vol: f64, cell_budget_m3: &mut [f64]) {
    if vol <= 0.0 {
        return;
    }
    if cp.vertex_idx < 0 {
        let b = &mut cell_budget_m3[cp.cell_idx as usize];
        *b = cmax(0.0, *b - vol);
        return;
    }
    let v = cp.vertex_idx as usize;
    let s = mesh.vert_stencil_ptr[v];
    let e = mesh.vert_stencil_ptr[v + 1];
    // Pass 1: weighted draw, clamped per cell; collect the undrawn remainder.
    let mut remainder = 0.0;
    for k in s..e {
        let b = &mut cell_budget_m3[mesh.vert_stencil_idx[k as usize] as usize];
        let want = vol * mesh.vert_stencil_wt[k as usize];
        let take = cmin(want, cmax(0.0, *b));
        *b -= take;
        remainder += want - take;
    }
    // Pass 2: greedy spill of the remainder into cells with budget left.
    let mut k = s;
    while k < e && remainder > 0.0 {
        let b = &mut cell_budget_m3[mesh.vert_stencil_idx[k as usize] as usize];
        let take = cmin(remainder, cmax(0.0, *b));
        *b -= take;
        remainder -= take;
        k += 1;
    }
}

/// budgetCredit: credit `vol` (m3) to the cells the point feeds, by the geometric weights.
#[inline]
fn budget_credit(mesh: &MeshData, cp: &CouplingPoint, vol: f64, cell_budget_m3: &mut [f64]) {
    if vol <= 0.0 {
        return;
    }
    if cp.vertex_idx < 0 {
        cell_budget_m3[cp.cell_idx as usize] += vol;
        return;
    }
    let v = cp.vertex_idx as usize;
    for k in mesh.vert_stencil_ptr[v]..mesh.vert_stencil_ptr[v + 1] {
        cell_budget_m3[mesh.vert_stencil_idx[k as usize] as usize] += vol * mesh.vert_stencil_wt[k as usize];
    }
}

/// accumulateOutfallDischargeStep(cps, mesh, state, ctx, opts, dt, accum_m3, cell_budget_m3,
/// sample_row): one routing step's net outfall exchange accumulated per point (m3, + = onto
/// the surface), withdrawals capped by the shared cell budget. Returns the number of outfalls
/// whose withdrawal was clamped. `state` is unused upstream too.
pub fn accumulate_outfall_discharge_step(
    cps: &[CouplingPoint],
    mesh: &MeshData,
    _state: &SurfaceStateData,
    ctx: &SimulationContext,
    opts: &SolverOptions2D,
    dt: f64,
    accum_m3: &mut [f64],
    cell_budget_m3: &mut [f64],
    mut sample_row: Option<&mut [f64]>,
) -> i32 {
    if dt <= 0.0 {
        return 0;
    }
    let nodes = &ctx.nodes;
    let mut clamped = 0;

    for k in 0..cps.len() {
        let cp = &cps[k];
        if !cp.is_outfall {
            continue;
        }

        let ni = cp.node_idx as usize;

        // Live net signed 1D -> 2D exchange this routing step.
        let mut q_net = (nodes.inflow[ni] - nodes.outflow[ni]) * opts.flow_1d_to_2d;

        // Withdrawal capped by the shared frozen-2D budget (window-cumulative).
        if q_net < 0.0 {
            let q_min = -budget_avail(mesh, cp, cell_budget_m3) / dt;
            if q_net < q_min {
                q_net = q_min;
                clamped += 1;
            }
        }
        if q_net == 0.0 {
            continue;
        }

        let vol_m3 = q_net * dt; // m3 this routing step
        if q_net < 0.0 {
            budget_draw(mesh, cp, -vol_m3, cell_budget_m3);
        } else {
            budget_credit(mesh, cp, vol_m3, cell_budget_m3);
        }

        accum_m3[k] += vol_m3;

        // Series sample: Q_net is already net-2D-source-positive.
        if let Some(row) = sample_row.as_deref_mut() {
            row[k] = q_net;
        }
    }
    clamped
}

/// injectAccumulatedExchange(cps, mesh, state, accum_m3, window_dt, sign, node_conc): each
/// point's accumulated volume as the window's mean rate, sign-adjusted, scattered into
/// state.coupling_flux (and, with node_conc, transport.coupling_src). Nothing is cleared.
pub fn inject_accumulated_exchange(cps: &[CouplingPoint], mesh: &MeshData, state: &mut SurfaceStateData, accum_m3: &[f64], window_dt: f64, sign: f64, node_conc: Option<&[f64]>) {
    if window_dt <= 0.0 {
        return;
    }
    let ns = state.transport.n_species as usize;
    for k in 0..cps.len() {
        let vol = accum_m3[k];
        if vol == 0.0 {
            continue;
        }
        // S3: a 2D source carries the node's published concentration, when the row exists at
        // this stride.
        let mut conc: Option<&[f64]> = None;
        if let Some(nc) = node_conc {
            if ns > 0 && cps[k].node_idx >= 0 {
                let ni = cps[k].node_idx as usize;
                if (ni + 1) * ns <= nc.len() {
                    conc = Some(&nc[ni * ns..]);
                }
            }
        }
        scatter_coupling_flux(mesh, state, &cps[k], sign * vol / window_dt, conc);
    }
}
