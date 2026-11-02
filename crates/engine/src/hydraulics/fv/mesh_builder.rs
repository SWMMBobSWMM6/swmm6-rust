// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/fv/NetworkMeshBuilder.hpp and
// NetworkMeshBuilder.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit
// 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the Apache License,
// Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// FV network mesh construction: buildClosure, buildGeometry and buildNetworkMesh.
//
// Differences in form, never in result:
//   * The tabulation memo `std::map<GeomKey,int>` is core::std_map::StdMap (libstdc++'s
//     tree), so a key holding NaN finds what upstream's finds. Its three table pointers
//     are keyed by the transect-table index + 1 (0 for null): the tables sit in one vector
//     in index order, so the pointer order and this order agree.
//   * FvGeometry::xs keeps the scalars and the transect index; the borrowed tables are
//     used while the closure is built (see hydraulics::fv). The face-section averaging
//     tests `a.xs.area_tbl != nullptr` as `a.xs.transect >= 0`: buildXSectParams attaches
//     tables exactly when it sets the transect index.
//   * std::max over an initializer list keeps its first-largest-wins comparison order.
//   * snprintf("%.0f") is Rust's `{:.0}` (both round the exact binary value to nearest,
//     ties to even).
//   * The local lambdas (add_face, chain_end_cell, invert_at, attach_cell_face, step)
//     take the mesh explicitly.

use crate::core::cmath::{cmax, cmin, cvt_i32};
use crate::core::constants;
use crate::core::context::SimulationContext;
use crate::core::error_codes::{format_warning, WARN_FV_SLOT_CAP};
use crate::core::options::{fv::FvOptions, Bytes};
use crate::core::std_map::StdMap;
use crate::core::unit_conversion as ucf;
use crate::data::link::{LinkType, XsectShape};
use crate::data::node::NodeType;
use crate::hydraulics::fv::kernels::{secgeom, slot_ramp, slot_ramp_integral, K_GRAVITY};
use crate::hydraulics::fv::{
    FvClosure, FvGeometry, NetworkMeshData, K_CLOSURE_PANELS, K_CLOSURE_POLYNOMIAL, K_CLOSURE_TABULATED, K_NODE_JUNCTION, K_NODE_OUTFALL, K_NODE_STORAGE, K_NODE_VIRTUAL, K_NODE_VOL_SAMPLES,
};
use crate::hydraulics::hyd_closure_kernels::CulvertCurve;
use crate::hydraulics::xsect_batch::{XSectParams, XSectShape};
use crate::hydraulics::{culvert, link, node, xsection as xsect};

/// Diagnostics produced while building the mesh.
#[derive(Debug, Clone, Default)]
pub struct MeshBuildReport {
    pub n_cells: i32,
    pub n_faces: i32,
    pub n_conduits: i32,
    pub n_virtual: i32,
    pub min_dx: f64,
    pub max_dx: f64,
    pub warnings: Vec<Bytes>,
    pub errors: Vec<Bytes>,
}

/// The scalars of an XSectParams, without its borrowed tables.
fn strip(xs: &XSectParams) -> XSectParams<'static> {
    XSectParams {
        type_: xs.type_,
        culvert_code: xs.culvert_code,
        transect: xs.transect,
        y_full: xs.y_full,
        w_max: xs.w_max,
        yw_max: xs.yw_max,
        a_full: xs.a_full,
        r_full: xs.r_full,
        s_full: xs.s_full,
        s_max: xs.s_max,
        y_bot: xs.y_bot,
        a_bot: xs.a_bot,
        s_bot: xs.s_bot,
        r_bot: xs.r_bot,
        area_tbl: None,
        hrad_tbl: None,
        width_tbl: None,
        transect_tbl_size: xs.transect_tbl_size,
        area_lut: None,
    }
}

/// buildClosure: the closure tables of `g`, sampled from `xs` (g.xs with its tables).
fn build_closure(g: &mut FvGeometry, xs: &XSectParams) {
    let mut c = FvClosure::default();
    let n = K_CLOSURE_PANELS;
    if !(g.y_full > 0.0) {
        g.closure_tbl = Box::new(c);
        return;
    }

    c.y_full = g.y_full;
    c.dh = g.y_full / n as f64;
    c.inv_dh = n as f64 / g.y_full;
    c.t_slot = g.t_slot;
    c.r_full = g.r_full;
    let bs = g.barrel_scale;
    let band = g.y_full - g.y_crown;

    if g.is_open != 0 && secgeom::is_polynomial_open(xs.type_) {
        c.kind = K_CLOSURE_POLYNOMIAL;
        match XSectShape::from_code(xs.type_) {
            Some(XSectShape::RectOpen) => {
                c.c2 = 0.0;
                c.c1 = bs * xs.w_max;
                c.p0 = bs * xs.w_max;
                c.p1 = bs * (2.0 - xs.s_bot);
            }
            Some(XSectShape::Trapezoidal) => {
                c.c2 = bs * xs.s_bot;
                c.c1 = bs * xs.y_bot;
                c.p0 = bs * xs.y_bot;
                c.p1 = bs * xs.r_bot;
            }
            _ => {
                // TRIANGULAR
                c.c2 = bs * xs.s_bot;
                c.c1 = 0.0;
                c.p0 = 0.0;
                c.p1 = bs * 2.0 * xs.r_bot;
            }
        }
        c.a_crown = (c.c2 * g.y_full + c.c1) * g.y_full;
        c.i1_crown = (c.c2 * g.y_full / 3.0 + 0.5 * c.c1) * g.y_full * g.y_full;
        c.inv_a_crown_n = if c.a_crown > 0.0 { n as f64 / c.a_crown } else { 0.0 };
    } else {
        c.kind = K_CLOSURE_TABULATED;
        for i in 0..=n {
            let h = g.y_full * i as f64 / n as f64;
            let mut a = bs * secgeom::area_of_depth(xs, h);
            let mut w = bs * secgeom::width_of_depth(xs, h);
            let r = secgeom::hyd_rad_of_depth(xs, h);
            if band > 0.0 && h > g.y_crown {
                let s = (h - g.y_crown) / band;
                a += g.t_slot * band * slot_ramp_integral(s);
                w += g.t_slot * slot_ramp(s);
            }
            c.A[i] = a;
            c.M[i] = if w > 0.0 { w } else { 0.0 };
            c.R[i] = if r > 0.0 { r } else { 0.0 };
        }
        c.A[0] = 0.0;
        c.R[n] = g.r_full;

        if !secgeom::has_closed_form(xs.type_) {
            let mut d = [0.0f64; K_CLOSURE_PANELS];
            for i in 0..n {
                d[i] = (c.A[i + 1] - c.A[i]) * c.inv_dh;
            }
            c.M[0] = d[0];
            c.M[n] = d[n - 1];
            for i in 1..n {
                if d[i - 1] > 0.0 && d[i] > 0.0 {
                    c.M[i] = 2.0 / (1.0 / d[i - 1] + 1.0 / d[i]);
                } else {
                    c.M[i] = 0.0;
                }
            }
            if band > 0.0 {
                c.M[n] = cmax(c.M[n], g.t_slot);
            }
        }

        // Fritsch-Carlson
        for i in 0..n {
            let delta = (c.A[i + 1] - c.A[i]) * c.inv_dh;
            if !(delta > 0.0) {
                c.M[i] = 0.0;
                c.M[i + 1] = 0.0;
                continue;
            }
            let mut alpha = c.M[i] / delta;
            let mut beta = c.M[i + 1] / delta;
            if alpha < 0.0 {
                alpha = 0.0;
                c.M[i] = 0.0;
            }
            if beta < 0.0 {
                beta = 0.0;
                c.M[i + 1] = 0.0;
            }
            let s2 = alpha * alpha + beta * beta;
            if s2 > 9.0 {
                let tau = 3.0 / s2.sqrt();
                c.M[i] = tau * alpha * delta;
                c.M[i + 1] = tau * beta * delta;
            }
        }

        // I1 = running integral of the panel cubics.
        c.I1[0] = 0.0;
        for i in 0..n {
            c.I1[i + 1] = c.I1[i] + c.dh * (0.5 * (c.A[i] + c.A[i + 1]) + c.dh / 12.0 * (c.M[i] - c.M[i + 1]));
        }
        c.a_crown = c.A[n];
        c.i1_crown = c.I1[n];
        c.inv_a_crown_n = if c.a_crown > 0.0 { n as f64 / c.a_crown } else { 0.0 };

        // Area-uniform panel index for the inverse.
        for k in 0..=n {
            let ak = c.a_crown * k as f64 / n as f64;
            let mut lo = 0usize;
            let mut hi = n;
            while hi - lo > 1 {
                let mid = (lo + hi) / 2;
                if c.A[mid] <= ak {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            c.j_of_a[k] = (if lo > n - 1 { n - 1 } else { lo }) as i16;
        }
    }

    g.a_crown = c.a_crown;
    g.i1_crown = c.i1_crown;
    g.closure_tbl = Box::new(c);
}

/// buildGeometry: an FV section block from a scalar section view (tables borrowed).
pub fn build_geometry(xs: &XSectParams, is_open: bool, slot_celerity: f64, g: &mut FvGeometry, barrels: i32) {
    g.xs = strip(xs);
    g.barrels = cmax(1, barrels);
    g.barrel_scale = g.barrels as f64;
    g.y_full = xs.y_full;
    g.a_full = g.barrel_scale * xs.a_full;
    g.w_max = g.barrel_scale * xs.w_max;
    g.r_full = xs.r_full;
    g.is_open = if is_open { 1 } else { 0 };

    if is_open {
        g.y_crown = xs.y_full;
        g.t_slot = if g.w_max > 0.0 { g.w_max } else { 1.0 };
    } else {
        g.y_crown = constants::SLOT_CROWN_CUTOFF * xs.y_full;
        let c = if slot_celerity > 1.0 { slot_celerity } else { 1.0 };
        g.t_slot = K_GRAVITY * g.a_full / (c * c);
        let cap = 0.05 * if g.w_max > 0.0 { g.w_max } else { 1.0 };
        if g.t_slot > cap {
            g.t_slot = cap;
        }
        if g.t_slot <= 0.0 {
            g.t_slot = 1.0e-6;
        }
    }

    let band = g.y_full - g.y_crown;
    g.a_crown = g.a_full + g.t_slot * band * 0.5;

    build_closure(g, xs);
}

/// The tabulation memo key: every XSectParams field, barrels and openness.
#[derive(Clone, Copy)]
struct GeomKey {
    ints: [i32; 3],
    dbl: [f64; 11],
    ptrs: [u64; 3],
    tbl_size: i32,
    barrels: i32,
    is_open: bool,
}

/// std::tuple's operator<: element by element, `<` then `>`.
fn key_less(a: &GeomKey, b: &GeomKey) -> bool {
    for k in 0..3 {
        if a.ints[k] < b.ints[k] {
            return true;
        }
        if b.ints[k] < a.ints[k] {
            return false;
        }
    }
    for k in 0..11 {
        if a.dbl[k] < b.dbl[k] {
            return true;
        }
        if b.dbl[k] < a.dbl[k] {
            return false;
        }
    }
    for k in 0..3 {
        if a.ptrs[k] < b.ptrs[k] {
            return true;
        }
        if b.ptrs[k] < a.ptrs[k] {
            return false;
        }
    }
    if a.tbl_size < b.tbl_size {
        return true;
    }
    if b.tbl_size < a.tbl_size {
        return false;
    }
    if a.barrels < b.barrels {
        return true;
    }
    if b.barrels < a.barrels {
        return false;
    }
    !a.is_open & b.is_open
}

fn geom_key(x: &XSectParams, barrels: i32, is_open: bool) -> GeomKey {
    let tp = |t: Option<&[f64]>| -> u64 {
        match t {
            None => 0,
            Some(_) => (x.transect as i64 + 1) as u64,
        }
    };
    GeomKey {
        ints: [x.type_, x.culvert_code, x.transect],
        dbl: [x.y_full, x.w_max, x.yw_max, x.a_full, x.r_full, x.s_full, x.s_max, x.y_bot, x.a_bot, x.s_bot, x.r_bot],
        ptrs: [tp(x.area_tbl), tp(x.hrad_tbl), tp(x.width_tbl)],
        tbl_size: x.transect_tbl_size,
        barrels,
        is_open,
    }
}

#[derive(Clone, Copy)]
struct Attachment {
    /// ConduitData row
    conduit: i32,
    /// 0 = conduit's upstream (node1) end, 1 = downstream
    end: i32,
}

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

/// std::max({a, b, c}): std::max_element, first-largest wins, `<` comparisons.
fn max3_f(a: f64, b: f64, c: f64) -> f64 {
    let mut m = a;
    if m < b {
        m = b;
    }
    if m < c {
        m = c;
    }
    m
}

fn max3_i(a: i32, b: i32, c: i32) -> i32 {
    let mut m = a;
    if m < b {
        m = b;
    }
    if m < c {
        m = c;
    }
    m
}

#[allow(clippy::too_many_arguments)]
fn add_face(mesh: &mut NetworkMeshData, cl: i32, cr: i32, node: i32, zb: f64, dx: f64, dl: i8, dr: i8, is_vj: bool) {
    mesh.face_cl.push(cl);
    mesh.face_cr.push(cr);
    mesh.face_node.push(node);
    mesh.face_zb.push(zb);
    mesh.face_dx.push(dx);
    mesh.face_dir_l.push(dl);
    mesh.face_dir_r.push(dr);
    mesh.face_virtual.push(if is_vj { 1 } else { 0 });
    mesh.face_vj_node.push(-1);
    mesh.face_gate.push(0);
    mesh.face_culvert.push(-1);
}

fn chain_end_cell(mesh: &NetworkMeshData, a: &Attachment) -> i32 {
    let ur = a.conduit as usize;
    let begin = mesh.conduit_cell_begin[ur];
    let count = mesh.conduit_cell_count[ur];
    if a.end == 0 {
        begin
    } else {
        begin + count - 1
    }
}

fn invert_at(mesh: &NetworkMeshData, ctx: &SimulationContext, a: &Attachment) -> f64 {
    let uj = mesh.conduit_link[a.conduit as usize] as usize;
    let n = if a.end == 0 { ctx.links.node1[uj] } else { ctx.links.node2[uj] };
    let off = if a.end == 0 { ctx.links.offset1[uj] } else { ctx.links.offset2[uj] };
    ctx.nodes.invert_elev[n as usize] + off
}

fn assign<T: Clone>(v: &mut Vec<T>, n: usize, d: T) {
    v.clear();
    v.resize(n, d);
}

pub fn build_network_mesh(ctx: &mut SimulationContext, opts: &FvOptions, mesh: &mut NetworkMeshData) -> MeshBuildReport {
    let mut rep = MeshBuildReport::default();
    mesh.clear();

    let n_links = ctx.n_links();
    let n_nodes = ctx.n_nodes();
    let n_cond = ctx.link_subtypes.conduits.count();
    let nn = n_nodes as usize;

    // ---------------------------------------------------------------- Nodes
    mesh.node_invert.resize(nn, 0.0);
    mesh.node_full_depth.resize(nn, 0.0);
    mesh.node_ponded_area.resize(nn, 0.0);
    mesh.node_sur_depth.resize(nn, 0.0);
    mesh.node_can_pond.resize(nn, 0);
    mesh.node_kind.resize(nn, 0);
    assign(&mut mesh.node_vj_face, nn, -1);
    assign(&mut mesh.node_area, nn, constants::MIN_SURFAREA);
    for i in 0..n_nodes {
        let ui = i as usize;
        mesh.node_invert[ui] = ctx.nodes.invert_elev[ui];
        mesh.node_full_depth[ui] = ctx.nodes.full_depth[ui];
        mesh.node_ponded_area[ui] = ctx.nodes.ponded_area[ui];
        mesh.node_sur_depth[ui] = ctx.nodes.sur_depth[ui];
        let is_coupled = ui < ctx.coupled_node.len() && ctx.coupled_node[ui] != 0;
        mesh.node_can_pond[ui] = ((ctx.options.allow_ponding || is_coupled) && ctx.nodes.ponded_area[ui] > 0.0 && ctx.nodes.type_[ui] != NodeType::Outfall) as u8;
        let mut kind = K_NODE_JUNCTION;
        if ui < ctx.nodes.is_virtual.len() && ctx.nodes.is_virtual[ui] != 0 {
            kind = K_NODE_VIRTUAL;
        } else if ctx.nodes.type_[ui] == NodeType::Storage || ctx.nodes.storage_area[ui] > 0.0 {
            kind = K_NODE_STORAGE;
        } else if ctx.nodes.type_[ui] == NodeType::Outfall {
            kind = K_NODE_OUTFALL;
        }
        mesh.node_kind[ui] = kind;

        let fd = ctx.nodes.full_depth[ui];
        let fvol = ctx.nodes.full_volume[ui];
        mesh.node_area[ui] = if fd > 0.0 && fvol > 0.0 { fvol / fd } else { constants::MIN_SURFAREA };
    }

    // Storage depth -> volume tables
    assign(&mut mesh.node_vol_off, nn, -1);
    assign(&mut mesh.node_vol_dmax, nn, 0.0);
    assign(&mut mesh.node_vol_atop, nn, 0.0);
    {
        let us = ucf::get_unit_system(ctx.options.flow_units as i32);
        for i in 0..n_nodes {
            let ui = i as usize;
            if mesh.node_kind[ui] != K_NODE_STORAGE {
                continue;
            }
            let mut dmax = ctx.nodes.full_depth[ui] + if ctx.nodes.storage_area[ui] > 0.0 { 0.0 } else { ctx.nodes.sur_depth[ui] };
            if !(dmax > 0.0) {
                dmax = cmax(ctx.nodes.full_depth[ui], 1.0);
            }
            mesh.node_vol_off[ui] = mesh.node_vol_tbl.len() as i32;
            mesh.node_vol_dmax[ui] = dmax;
            let dd = dmax / (K_NODE_VOL_SAMPLES - 1) as f64;
            for k in 0..K_NODE_VOL_SAMPLES {
                let d = k as f64 * dd;
                mesh.node_vol_tbl.push(node::get_volume(&ctx.nodes, i, d, Some(&ctx.tables), us, Some(&ctx.node_subtypes)));
            }
            mesh.node_vol_atop[ui] = if ctx.nodes.storage_area[ui] > 0.0 {
                0.0
            } else {
                cmax(node::get_surf_area(&ctx.nodes, i, dmax, Some(&ctx.tables), us, Some(&ctx.node_subtypes)), constants::MIN_SURFAREA)
            };
        }
    }

    // ------------------------------------------ Per-conduit geometry + cells
    mesh.geom.clear();
    {
        let un_cond = n_cond as usize;
        assign(&mut mesh.conduit_section, un_cond, -1);
        assign(&mut mesh.conduit_roughness, un_cond, 0.01);
        assign(&mut mesh.conduit_rough_factor, un_cond, 0.0);
        assign(&mut mesh.conduit_loss_inlet, un_cond, 0.0);
        assign(&mut mesh.conduit_loss_outlet, un_cond, 0.0);
        assign(&mut mesh.conduit_slope, un_cond, 0.0);
        assign(&mut mesh.conduit_culvert_code, un_cond, 0);
        assign(&mut mesh.conduit_culvert_curve, un_cond, CulvertCurve::default());
        assign(&mut mesh.conduit_culvert_mitered, un_cond, 0);
    }
    assign(&mut mesh.conduit_cell_begin, n_cond as usize, -1);
    assign(&mut mesh.conduit_cell_count, n_cond as usize, 0);
    assign(&mut mesh.conduit_link, n_cond as usize, -1);

    rep.min_dx = 1.0e30;

    let mut geom_memo: StdMap<GeomKey, i32> = StdMap::new(key_less);

    let mut n_closed = 0i32;
    let mut n_capped = 0i32;
    let mut c_implied_max = 0.0f64;

    for r in 0..n_cond {
        let ur = r as usize;
        let j = ctx.link_subtypes.conduits.link_idx[ur];
        if j < 0 || j >= n_links {
            continue;
        }
        let uj = j as usize;

        if ctx.links.xsect_shape[uj] == XsectShape::Dummy {
            continue;
        }

        mesh.conduit_link[ur] = j;

        let xs = link::build_xsect_params(&ctx.links, uj, Some(&ctx.transect_tables));
        if xs.y_full <= 0.0 || xs.a_full <= 0.0 {
            rep.errors.push(cat(&[
                b"FV routing: conduit '",
                ctx.link_names.name_of(j),
                b"' has no usable cross-section (full depth/area is zero). FLOW_ROUTING FV requires every conduit to carry real geometry.",
            ]));
            continue;
        }

        let cd = &ctx.link_subtypes.conduits;
        let is_open = xsect::is_open(xs.type_);
        let key = geom_key(&xs, cd.barrels[ur], is_open);
        let sec = match geom_memo.find(&key) {
            Some(&s) => s,
            None => {
                let s = mesh.geom.len() as i32;
                let mut g = FvGeometry::default();
                build_geometry(&xs, is_open, opts.slot_celerity, &mut g, cd.barrels[ur]);
                mesh.geom.push(g);
                geom_memo.emplace(key, s);
                s
            }
        };
        mesh.conduit_section[ur] = sec;
        let g = &mesh.geom[sec as usize];
        if g.is_open == 0 {
            let c_req = cmax(opts.slot_celerity, 1.0);
            let uncapped = K_GRAVITY * g.a_full / (c_req * c_req);
            let cap = 0.05 * if g.w_max > 0.0 { g.w_max } else { 1.0 };
            n_closed += 1;
            if uncapped > cap && g.t_slot > 0.0 {
                n_capped += 1;
                c_implied_max = cmax(c_implied_max, (K_GRAVITY * g.a_full / g.t_slot).sqrt());
            }
        }
        mesh.conduit_roughness[ur] = cd.roughness[ur];
        mesh.conduit_rough_factor[ur] = cd.rough_factor[ur];
        mesh.conduit_loss_inlet[ur] = cd.loss_inlet[ur];
        mesh.conduit_loss_outlet[ur] = cd.loss_outlet[ur];
        mesh.conduit_culvert_code[ur] = cd.culvert_code[ur];
        mesh.conduit_slope[ur] = cd.slope[ur];
        if cd.culvert_code[ur] > 0 {
            let cc = culvert::get_coeffs(cd.culvert_code[ur]);
            mesh.conduit_culvert_curve[ur] = CulvertCurve { K: cc.K, M: cc.M, C: cc.C, Y: cc.Y };
            let code = cd.culvert_code[ur];
            mesh.conduit_culvert_mitered[ur] = (code == 5 || code == 37 || code == 46) as u8;
        }

        let len = cd.length[ur];
        let mlen = cd.mod_length[ur];
        let l = max3_f(len, mlen, constants::FUDGE);

        let mut ncell = 1;
        if opts.cell_length > 0.0 {
            ncell = cvt_i32((l / opts.cell_length).ceil());
        }
        ncell = max3_i(ncell, opts.min_cells, 1);
        mesh.conduit_cell_begin[ur] = mesh.n_cells();
        mesh.conduit_cell_count[ur] = ncell;

        let dx = l / ncell as f64;
        rep.min_dx = cmin(rep.min_dx, dx);
        rep.max_dx = cmax(rep.max_dx, dx);

        let n1 = ctx.links.node1[uj];
        let n2 = ctx.links.node2[uj];
        let z1 = ctx.nodes.invert_elev[n1 as usize] + ctx.links.offset1[uj];
        let z2 = ctx.nodes.invert_elev[n2 as usize] + ctx.links.offset2[uj];

        for i in 0..ncell {
            let t = (i as f64 + 0.5) / ncell as f64;
            mesh.cell_geom.push(sec);
            mesh.cell_conduit.push(r);
            mesh.cell_dx.push(dx);
            mesh.cell_zb.push(z1 + (z2 - z1) * t);
            mesh.cell_dzdx.push((z2 - z1) / l);
        }
    }
    if !rep.errors.is_empty() {
        return rep;
    }
    if rep.min_dx > 1.0e29 {
        rep.min_dx = 0.0;
    }

    if n_capped > 0 {
        let buf = format!("{} of {} closed conduits (cap-implied celerity up to {:.0} ft/s)", n_capped, n_closed, c_implied_max);
        rep.warnings.push(format_warning(WARN_FV_SLOT_CAP, buf.as_bytes()));
    }

    // ---------------------------------------- Interior faces inside each conduit
    for r in 0..n_cond {
        let ur = r as usize;
        let begin = mesh.conduit_cell_begin[ur];
        let count = mesh.conduit_cell_count[ur];
        if begin < 0 || count < 2 {
            continue;
        }
        let uj = mesh.conduit_link[ur] as usize;
        let n1 = ctx.links.node1[uj];
        let n2 = ctx.links.node2[uj];
        let z1 = ctx.nodes.invert_elev[n1 as usize] + ctx.links.offset1[uj];
        let z2 = ctx.nodes.invert_elev[n2 as usize] + ctx.links.offset2[uj];
        let dx = mesh.cell_dx[begin as usize];
        for i in 1..count {
            let t = i as f64 / count as f64;
            add_face(mesh, begin + i - 1, begin + i, -1, z1 + (z2 - z1) * t, dx, 1, 1, false);
        }
    }

    // ------------------------------------------------------- Node attachments
    let mut attach: Vec<Vec<Attachment>> = vec![Vec::new(); nn];
    for r in 0..n_cond {
        let ur = r as usize;
        if mesh.conduit_cell_begin[ur] < 0 {
            continue;
        }
        let uj = mesh.conduit_link[ur] as usize;
        let n1 = ctx.links.node1[uj];
        let n2 = ctx.links.node2[uj];
        if n1 >= 0 && n1 < n_nodes {
            attach[n1 as usize].push(Attachment { conduit: r, end: 0 });
        }
        if n2 >= 0 && n2 < n_nodes {
            attach[n2 as usize].push(Attachment { conduit: r, end: 1 });
        }
    }

    // ----------------------------------------------------- Flap gates, per conduit
    let mut conduit_gate = vec![0u8; n_cond as usize];
    for r in 0..n_cond {
        let ur = r as usize;
        if mesh.conduit_cell_begin[ur] < 0 {
            continue;
        }
        let uj = mesh.conduit_link[ur] as usize;
        let mut gate: u8 = 0;
        let mut block = |s: i32| {
            gate |= if s > 0 { 1 } else { 2 };
        };

        if ctx.links.has_flap_gate[uj] != 0 {
            block(-ctx.links.direction[uj]);
        }

        let ends = [ctx.links.node1[uj], ctx.links.node2[uj]];
        for (e, &nd) in ends.iter().enumerate() {
            if nd < 0 || ctx.nodes.type_[nd as usize] != NodeType::Outfall {
                continue;
            }
            let ofr = ctx.node_subtypes.outfall_row(nd);
            if ofr >= 0 && ctx.node_subtypes.outfalls.has_flap_gate[ofr as usize] != 0 {
                block(if e == 0 { 1 } else { -1 });
            }
        }
        conduit_gate[ur] = gate;
    }

    let mut nf_idx: Vec<Vec<i32>> = vec![Vec::new(); nn];
    let mut nf_sign: Vec<Vec<f64>> = vec![Vec::new(); nn];
    let mut nf_zb: Vec<Vec<f64>> = vec![Vec::new(); nn];

    for i in 0..n_nodes {
        let ui = i as usize;
        let at = &attach[ui];

        if mesh.node_kind[ui] == K_NODE_VIRTUAL {
            if at.len() != 2 {
                rep.errors.push(cat(&[
                    b"FV routing: virtual junction '",
                    ctx.node_names.name_of(i),
                    b"' must join exactly two conduits (found ",
                    at.len().to_string().as_bytes(),
                    b").",
                ]));
                continue;
            }
            let a0 = at[0];
            let a1 = at[1];
            let cl = chain_end_cell(mesh, &a0);
            let cr = chain_end_cell(mesh, &a1);
            let dl: i8 = if a0.end == 1 { 1 } else { -1 };
            let dr: i8 = if a1.end == 0 { 1 } else { -1 };
            let dxl = mesh.cell_dx[cl as usize];
            let dxr = mesh.cell_dx[cr as usize];
            let zb = mesh.node_invert[ui];
            add_face(mesh, cl, cr, -1, zb, 0.5 * (dxl + dxr), dl, dr, true);
            *mesh.face_vj_node.last_mut().unwrap() = i;
            mesh.node_vj_face[ui] = mesh.n_faces() - 1;
            rep.n_virtual += 1;
            continue;
        }

        for a in at {
            let cell = chain_end_cell(mesh, a);
            let dx = mesh.cell_dx[cell as usize];
            let zb = invert_at(mesh, ctx, a);
            let fidx = mesh.n_faces();
            if a.end == 0 {
                add_face(mesh, -1, cell, i, zb, 0.5 * dx, 1, 1, false);
                nf_sign[ui].push(-1.0);
            } else {
                add_face(mesh, cell, -1, i, zb, 0.5 * dx, 1, 1, false);
                nf_sign[ui].push(1.0);
            }
            nf_idx[ui].push(fidx);
            nf_zb[ui].push(zb);
            mesh.face_gate[fidx as usize] = conduit_gate[a.conduit as usize];

            if a.end == 0 && mesh.conduit_culvert_code[a.conduit as usize] > 0 {
                mesh.face_culvert[fidx as usize] = a.conduit;
            }
        }
    }
    if !rep.errors.is_empty() {
        return rep;
    }

    // Flatten to CSR.
    assign(&mut mesh.node_face_ptr, nn + 1, 0);
    for i in 0..nn {
        mesh.node_face_ptr[i + 1] = mesh.node_face_ptr[i] + nf_idx[i].len() as i32;
    }
    let total = mesh.node_face_ptr[nn] as usize;
    mesh.node_face_idx.reserve(total);
    mesh.node_face_sign.reserve(total);
    mesh.node_face_zb.reserve(total);
    for i in 0..nn {
        for k in 0..nf_idx[i].len() {
            mesh.node_face_idx.push(nf_idx[i][k]);
            mesh.node_face_sign.push(nf_sign[i][k]);
            mesh.node_face_zb.push(nf_zb[i][k]);
        }
    }

    // ------------------------------------------------------------- Face sections
    mesh.face_geom.clear();
    mesh.derive_face_geom();
    {
        // looked up only: a BTreeMap is as good as std::map on (int,int) keys
        let mut pair_geom: std::collections::BTreeMap<(i32, i32), i32> = std::collections::BTreeMap::new();
        for i in 0..nn {
            let nb = mesh.node_face_ptr[i];
            let ne = mesh.node_face_ptr[i + 1];
            if ne - nb != 2 {
                continue;
            }
            let f0 = mesh.node_face_idx[nb as usize];
            let f1 = mesh.node_face_idx[(nb + 1) as usize];
            let g0 = mesh.face_geom[f0 as usize];
            let g1 = mesh.face_geom[f1 as usize];
            if g0 < 0 || g1 < 0 || g0 == g1 {
                continue;
            }
            let (a_xs, a_barrels) = (mesh.geom[g0 as usize].xs, mesh.geom[g0 as usize].barrels);
            let (b_xs, b_barrels) = (mesh.geom[g1 as usize].xs, mesh.geom[g1 as usize].barrels);
            if a_xs.type_ != b_xs.type_ || a_barrels != b_barrels || a_xs.transect >= 0 || b_xs.transect >= 0 {
                continue;
            }
            let key = (cmin(g0, g1), cmax(g0, g1));
            let gf = match pair_geom.get(&key) {
                Some(&g) => g,
                None => {
                    let mut xs = a_xs;
                    let xb = b_xs;
                    let mid = |p: f64, q: f64| 0.5 * (p + q);
                    xs.y_full = mid(a_xs.y_full, xb.y_full);
                    xs.w_max = mid(a_xs.w_max, xb.w_max);
                    xs.yw_max = mid(a_xs.yw_max, xb.yw_max);
                    xs.a_full = mid(a_xs.a_full, xb.a_full);
                    xs.r_full = mid(a_xs.r_full, xb.r_full);
                    xs.s_full = mid(a_xs.s_full, xb.s_full);
                    xs.s_max = mid(a_xs.s_max, xb.s_max);
                    xs.y_bot = mid(a_xs.y_bot, xb.y_bot);
                    xs.a_bot = mid(a_xs.a_bot, xb.a_bot);
                    xs.s_bot = mid(a_xs.s_bot, xb.s_bot);
                    xs.r_bot = mid(a_xs.r_bot, xb.r_bot);
                    let mut g = FvGeometry::default();
                    build_geometry(&xs, xsect::is_open(xs.type_), opts.slot_celerity, &mut g, a_barrels);
                    let gf = mesh.geom.len() as i32;
                    mesh.geom.push(g);
                    pair_geom.insert(key, gf);
                    gf
                }
            };
            mesh.face_geom[f0 as usize] = gf;
            mesh.face_geom[f1 as usize] = gf;
        }
    }

    // ---------------------------------------------------- Cell -> face references
    let nc = mesh.n_cells();
    let unc = nc as usize;
    assign(&mut mesh.cell_face0, unc, -1);
    assign(&mut mesh.cell_face1, unc, -1);
    assign(&mut mesh.cell_side0, unc, 0);
    assign(&mut mesh.cell_side1, unc, 0);
    for f in 0..mesh.n_faces() {
        let uf = f as usize;
        for (cell, side) in [(mesh.face_cl[uf], 0i8), (mesh.face_cr[uf], 1i8)] {
            if cell < 0 {
                continue;
            }
            let uc = cell as usize;
            if mesh.cell_face0[uc] < 0 {
                mesh.cell_face0[uc] = f;
                mesh.cell_side0[uc] = side;
            } else {
                mesh.cell_face1[uc] = f;
                mesh.cell_side1[uc] = side;
            }
        }
    }
    for c in 0..nc {
        let uc = c as usize;
        if mesh.cell_face0[uc] < 0 || mesh.cell_face1[uc] < 0 {
            rep.errors.push(format!("FV routing: internal error \u{2014} cell {} is not bounded by two faces (dangling conduit end).", c).into_bytes());
            break;
        }
    }
    if !rep.errors.is_empty() {
        return rep;
    }

    // ----------------------------------------------------------------- Cell chains
    assign(&mut mesh.cell_chain, unc, -1);
    assign(&mut mesh.cell_chain_pos, unc, -1);
    mesh.chain_ptr.push(0);

    for c0 in 0..nc {
        let u0 = c0 as usize;
        if mesh.cell_chain[u0] >= 0 {
            continue;
        }
        let fa = mesh.cell_face0[u0];
        let fb = mesh.cell_face1[u0];
        let end_a = mesh.face_cl[fa as usize] < 0 || mesh.face_cr[fa as usize] < 0;
        let end_b = mesh.face_cl[fb as usize] < 0 || mesh.face_cr[fb as usize] < 0;
        if !end_a && !end_b {
            continue;
        }

        let chain = mesh.n_chains();
        let mut cur = c0;
        let mut dir: i8 = 1;
        let mut in_face = if end_a { fa } else { fb };
        let mut pos = 0;
        while cur >= 0 {
            let uc = cur as usize;
            mesh.cell_chain[uc] = chain;
            mesh.cell_chain_pos[uc] = pos;
            pos += 1;
            mesh.chain_cells.push(cur);
            mesh.chain_dir.push(dir);

            let out_face = if mesh.cell_face0[uc] == in_face { mesh.cell_face1[uc] } else { mesh.cell_face0[uc] };
            let out_side = if mesh.cell_face0[uc] == in_face { mesh.cell_side1[uc] } else { mesh.cell_side0[uc] };
            // step(): the neighbour across out_face and the axis flip
            let uf = out_face as usize;
            let flip = (mesh.face_dir_l[uf] as i32 * mesh.face_dir_r[uf] as i32) as i8;
            let next = if out_side == 0 { mesh.face_cr[uf] } else { mesh.face_cl[uf] };
            if next < 0 || mesh.cell_chain[next as usize] >= 0 {
                break;
            }
            dir = (dir as i32 * flip as i32) as i8;
            in_face = out_face;
            cur = next;
        }
        mesh.chain_ptr.push(mesh.chain_cells.len() as i32);
    }
    for c in 0..nc {
        if mesh.cell_chain[c as usize] < 0 {
            rep.errors.push(
                format!(
                    "FV routing: internal error \u{2014} cell {} is not reachable from any chain end (closed conduit loop). A cycle of conduits joined only by virtual junctions has no boundary and cannot be meshed.",
                    c
                )
                .into_bytes(),
            );
            break;
        }
    }
    if !rep.errors.is_empty() {
        return rep;
    }

    // ------------------------------------------------------------ Non-conduit links
    assign(&mut mesh.node_dummy_drain, nn, 0);
    for j in 0..n_links {
        let uj = j as usize;
        let is_dummy = ctx.links.type_[uj] == LinkType::Conduit && ctx.links.xsect_shape[uj] == XsectShape::Dummy;
        if ctx.links.type_[uj] == LinkType::Conduit && !is_dummy {
            continue;
        }
        let a = ctx.links.node1[uj];
        let b = ctx.links.node2[uj];
        if a < 0 || b < 0 {
            continue;
        }
        mesh.struct_link.push(j);
        mesh.struct_n1.push(a);
        mesh.struct_n2.push(b);
        mesh.struct_is_dummy.push(if is_dummy { 1 } else { 0 });
        if is_dummy {
            mesh.node_dummy_drain[a as usize] = 1;
        }
    }
    if !mesh.struct_is_dummy.contains(&1) {
        mesh.struct_is_dummy.clear();
        mesh.node_dummy_drain.clear();
    }

    rep.n_cells = nc;
    rep.n_faces = mesh.n_faces();
    rep.n_conduits = n_cond;
    rep
}
