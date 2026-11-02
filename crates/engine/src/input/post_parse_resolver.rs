// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/PostParseResolver.cpp and
// PostParseResolver.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Post-parse cross-reference resolution (namespace openswmm::input): unit conversion of
// the parsed inputs, name binding, external files, cross-section and conduit properties,
// virtual- and inlet-junction validation. The external-file loaders are in
// post_parse_resolver/external_files.rs. The 1D build's code path only: the
// OPENSWMM_HAS_2D block of resolve_external_file_slots is not compiled there.
//
// Differences in form, never in result:
//   * Upstream passes `SimulationContext&` to helpers that also read names out of it;
//     here those names are copied first where the borrow checker needs it.
//   * A pointer to "the crest" (offset1 or a weir/outlet row's crest_height) is the enum
//     `Crest`, read and written through `get`/`set`.
//   * The case-insensitive name maps (CiHash/CiEqual) are only looked up, and emplace
//     keeps the first index, so each is a search for the first ieq() match.
//   * The CUSTOM memo std::map<std::pair<int,double>,int> is core::std_map::StdMap, a
//     replica of libstdc++'s tree, because a NaN depth makes the key's order non-strict
//     and find() then depends on the tree's shape.
//   * perf:: timers are dropped (they only accumulate wall time).
//   * recompute_conduit_flow_properties reads CD.slope[conduit_row(j)] before its
//     early return; for a conduit with no side-table row upstream reads out of bounds
//     and then returns (n is 0). Here the read is skipped in that case; the result is
//     the same (only xsect_generation changes).

mod external_files;

pub use external_files::load_external_rain_files;

use crate::core::cmath::cvt_i32;
use crate::core::constants::{GRAVITY, MIN_DELTA_Z, MISSING, PHI, TINY};
use crate::core::context::SimulationContext;
use crate::core::cstdio::{printf_g, to_string_long, IfStream};
use crate::core::error_codes::{
    format_error, format_error_detail, format_warning, ERR_CURVE_SEQUENCE, ERR_DIVIDER_LINK, ERR_IJ_CAPTURE_NODE,
    ERR_IJ_DESIGN, ERR_IJ_NOT_STREET, ERR_IJ_NO_USAGE, ERR_IJ_USAGE_ON_PAIR, ERR_INLET_USAGE_SHAPE, ERR_NAME,
    ERR_NUMBER, ERR_RAIN_FILE_CONFLICT, ERR_TIMESERIES_SEQUENCE, ERR_TRANSECT_MANNING, ERR_VJ_ROUTING_MODEL,
    WARN_ELEV_DROP_EXCEEDS, WARN_MAX_DEPTH_INCREASED, WARN_MIN_ELEV_DROP, WARN_MIN_SLOPE, WARN_NEGATIVE_OFFSET,
    WARN_REGULATOR_CREST_LOW,
};
use crate::core::options::{Bytes, FilePathPair, RoutingModel};
use crate::core::path_resolver::{parent_dir, resolve_relative};
use crate::core::std_map::StdMap;
use crate::core::unit_conversion::{get_unit_system, LENGTH, QCF, RAINDEPTH, RAINFALL, UCF, UCF_INV};
use crate::data::gage::RainSource;
use crate::data::initial_quality::InitialQualityData;
use crate::data::link::{LinkType, XsectShape};
use crate::data::name_index::ieq;
use crate::data::node::{NodeType, OutfallType};
use crate::data::table::TableType;
use crate::edit::virtual_junction_ops as edit;
use crate::hydraulics::force_main::{self as forcemain, FrictionModel};
use crate::hydraulics::link::{apply_tabulated_xsect_params, translate_shape};
use crate::hydraulics::street::{self, StreetParams};
use crate::hydraulics::transect;
use crate::hydraulics::transect_types::TransectData;
use crate::hydraulics::xsect_batch::XSectParams;
use crate::hydraulics::xsection as xsect;
use crate::input::multi_column_series_file::MultiColumnFileCache;
use crate::input::tokenizer;

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut v = Bytes::new();
    for p in parts {
        v.extend_from_slice(p);
    }
    v
}

/// The first index whose name is ieq() to `name`, or -1 (an emplace-first CI map).
fn find_ci(names: &[Bytes], name: &[u8]) -> i32 {
    for (i, n) in names.iter().enumerate() {
        if ieq(n, name) {
            return i as i32;
        }
    }
    -1
}

/// `static_cast<forcemain::FrictionModel>(eqn)`; None for a value outside the enum.
fn friction_model(eqn: i32) -> Option<FrictionModel> {
    match eqn {
        0 => Some(FrictionModel::HazenWilliams),
        1 => Some(FrictionModel::DarcyWeisbach),
        _ => None,
    }
}

/// forcemain::getEquivN; its switch falls through to n_raw for any other model value.
fn get_equiv_n(eqn: i32, r_bot: f64, y_full: f64, slope: f64, n_raw: f64) -> f64 {
    if slope <= 0.0 || y_full <= 0.0 {
        return n_raw;
    }
    match friction_model(eqn) {
        Some(m) => forcemain::get_equiv_n(m, r_bot, y_full, slope, n_raw),
        None => n_raw,
    }
}

/// forcemain::getRoughFactor; 0 for any other model value.
fn get_rough_factor(eqn: i32, r_bot: f64, length_factor: f64) -> f64 {
    match friction_model(eqn) {
        Some(m) => forcemain::get_rough_factor(m, r_bot, length_factor),
        None => 0.0,
    }
}

// legacy link.c getOffsetHeight: an ELEVATION offset as a height above the node invert.
fn legacy_offset_height(offset: f64, invert_elev: f64, warnings: &mut Vec<Bytes>, link_name: &[u8]) -> f64 {
    if offset <= MISSING {
        return 0.0;
    }
    let offset = offset - invert_elev;
    if offset >= 0.0 {
        return offset;
    }
    if offset >= -MIN_DELTA_Z {
        return 0.0;
    }
    warnings.push(format_warning(WARN_NEGATIVE_OFFSET, link_name));
    0.0
}

/// Where a regulator's crest lives: links.offset1, or a weir / outlet row's crest_height.
#[derive(Clone, Copy)]
enum Crest {
    Offset1(usize),
    Weir(usize),
    Outlet(usize),
}

impl Crest {
    fn get(self, ctx: &SimulationContext) -> f64 {
        match self {
            Crest::Offset1(j) => ctx.links.offset1[j],
            Crest::Weir(r) => ctx.link_subtypes.weirs.crest_height[r],
            Crest::Outlet(r) => ctx.link_subtypes.outlets.crest_height[r],
        }
    }
    fn set(self, ctx: &mut SimulationContext, v: f64) {
        match self {
            Crest::Offset1(j) => ctx.links.offset1[j] = v,
            Crest::Weir(r) => ctx.link_subtypes.weirs.crest_height[r] = v,
            Crest::Outlet(r) => ctx.link_subtypes.outlets.crest_height[r] = v,
        }
    }
}

/// The weir row's crest, else the outlet row's, else none.
fn weir_or_outlet_crest(ctx: &SimulationContext, j: i32) -> Option<Crest> {
    let wr = ctx.link_subtypes.weir_row(j);
    let olr = if wr < 0 { ctx.link_subtypes.outlet_row(j) } else { -1 };
    if wr >= 0 {
        Some(Crest::Weir(wr as usize))
    } else if olr >= 0 {
        Some(Crest::Outlet(olr as usize))
    } else {
        None
    }
}

// -------------------------------------------------------------------------
// Slice IO-3: resolve every external-file slot's `original` against the .inp directory.
// -------------------------------------------------------------------------

pub fn resolve_external_file_slots(ctx: &mut SimulationContext, anchor_dir: &[u8]) {
    let resolve = |slot: &mut FilePathPair| {
        if slot.original.is_empty() {
            slot.absolute.clear();
            return;
        }
        slot.absolute = resolve_relative(&slot.original, anchor_dir);
    };

    let f = &mut ctx.files;
    resolve(&mut f.rainfall_path);
    resolve(&mut f.runoff_path);
    resolve(&mut f.rdii_path);
    resolve(&mut f.inflows_path);
    resolve(&mut f.outflows_path);
    resolve(&mut f.hotstart_use_path);
    for save in f.hotstart_saves.iter_mut() {
        resolve(&mut save.path);
    }

    for gfp in ctx.gages.file_path.iter_mut() {
        resolve(gfp);
    }

    resolve(&mut ctx.options.temp_file);
    resolve(&mut ctx.snowpacks.config.met_file);

    for tbl in ctx.tables.tables.iter_mut() {
        if tbl.type_ != TableType::Timeseries {
            continue;
        }
        resolve(&mut tbl.file_path);
    }

    for rpt in ctx.lid_usage.rpt_file.iter_mut() {
        resolve(rpt);
    }
    // (2D mesh_file / output_file: OPENSWMM_HAS_2D only.)
}

// -------------------------------------------------------------------------
// Conduit roughness, length and flow properties
// -------------------------------------------------------------------------

/// The Manning's n a conduit routes with (transect n, force-main equivalent n, meander
/// factor), or 0 for a non-conduit or a conduit without a row.
pub fn conduit_manning_n(ctx: &SimulationContext, j: i32) -> f64 {
    let uj = j as usize;
    if ctx.links.type_[uj] != LinkType::Conduit {
        return 0.0;
    }
    let cr = ctx.link_subtypes.conduit_row(j);
    if cr < 0 {
        return 0.0;
    }
    let cd = &ctx.link_subtypes.conduits;
    let ucr = cr as usize;
    let shape = ctx.links.xsect_shape[uj];

    let mut n_val = cd.roughness[ucr];

    let ti = if shape == XsectShape::Irregular { ctx.links.xsect_curve[uj] } else { -1 };
    if ti >= 0 && (ti as usize) < ctx.transects.n_channel.len() {
        let nch = ctx.transects.n_channel[ti as usize];
        if nch > 0.0 {
            n_val = nch;
        }
    }

    let is_dw = ctx.options.routing_model == RoutingModel::Dynwave || ctx.options.routing_model == RoutingModel::Fv;
    if is_dw && shape == XsectShape::ForceMain {
        let r_bot = ctx.links.xsect_r_bot[uj];
        let y_full = ctx.links.xsect_y_full[uj];
        let slope = cd.slope[ucr].abs();
        let eq = get_equiv_n(ctx.options.force_main_eqn, r_bot, y_full, slope, n_val);
        n_val = if eq > 0.0 { eq } else { cd.roughness[ucr] };
    }

    if ti >= 0 && (ti as usize) < ctx.transects.length_factor.len() {
        let lf = ctx.transects.length_factor[ti as usize];
        if lf > 0.0 {
            n_val *= lf.sqrt();
        }
    }
    n_val
}

/// The length a conduit routes with: the authored length, divided by the transect's
/// length factor for an IRREGULAR section; 0 for a non-conduit or missing row.
pub fn conduit_true_length(ctx: &SimulationContext, j: i32) -> f64 {
    if j < 0 || j >= ctx.n_links() {
        return 0.0;
    }
    let uj = j as usize;
    if ctx.links.type_[uj] != LinkType::Conduit {
        return 0.0;
    }
    let cr = ctx.link_subtypes.conduit_row(j);
    if cr < 0 {
        return 0.0;
    }
    let len = ctx.link_subtypes.conduits.length[cr as usize];
    if ctx.links.xsect_shape[uj] != XsectShape::Irregular {
        return len;
    }
    let ti = ctx.links.xsect_curve[uj];
    if ti < 0 || ti as usize >= ctx.transects.length_factor.len() {
        return len;
    }
    let lf = ctx.transects.length_factor[ti as usize];
    if lf > 0.0 {
        len / lf
    } else {
        len
    }
}

/// Recompute one conduit's rough_factor, beta, q_full, q_max, volume (and a dynamic-wave
/// force main's roughness factor in xsect_s_bot).
pub fn recompute_conduit_flow_properties(ctx: &mut SimulationContext, j: i32) {
    let uj = j as usize;
    if ctx.links.type_[uj] != LinkType::Conduit {
        return;
    }
    let cr = ctx.link_subtypes.conduit_row(j);

    ctx.xsect_generation = ctx.xsect_generation.wrapping_add(1);

    let n_val = conduit_manning_n(ctx, j);
    // (see the header: upstream reads slope[conduit_row] here even when the row is -1)
    let slope = if cr >= 0 { ctx.link_subtypes.conduits.slope[cr as usize].abs() } else { 0.0 };
    let a_full = ctx.links.xsect_a_full[uj];
    let s_full = ctx.links.xsect_s_full[uj];
    let s_max = ctx.links.xsect_s_max[uj];

    if n_val <= 0.0 || a_full <= 0.0 {
        return;
    }
    let ucr = cr as usize;
    let cd = &mut ctx.link_subtypes.conduits;

    cd.rough_factor[ucr] = GRAVITY * ((n_val / PHI) * (n_val / PHI));

    let beta = PHI * slope.sqrt() / n_val;
    cd.beta[ucr] = beta;
    cd.q_full[ucr] = s_full * beta;
    cd.q_max[ucr] = s_max * beta;

    let mut mod_len = cd.mod_length[ucr];
    if mod_len <= 0.0 {
        mod_len = cd.length[ucr];
    }
    cd.mod_length[ucr] = mod_len;
    let length = cd.length[ucr];
    ctx.links.volume[uj] = a_full * mod_len;

    let is_force_main = ctx.links.xsect_shape[uj] == XsectShape::ForceMain;
    let is_dw = ctx.options.routing_model == RoutingModel::Dynwave || ctx.options.routing_model == RoutingModel::Fv;
    if is_dw && is_force_main {
        let length_factor = if length > 0.0 { mod_len / length } else { 1.0 };
        ctx.links.xsect_s_bot[uj] = get_rough_factor(ctx.options.force_main_eqn, ctx.links.xsect_r_bot[uj], length_factor);
    }
}

// -------------------------------------------------------------------------
// Display -> internal unit conversion (and its inverse)
// -------------------------------------------------------------------------

fn convert_inputs_to_internal(ctx: &mut SimulationContext, n_nodes: i32, n_links: i32, n_subcatch: i32) {
    let us = get_unit_system(ctx.options.flow_units as i32);
    let usz = us as usize;

    {
        let rain0 = UCF[RAINFALL as usize][usz];
        if rain0 != 1.0 {
            for j in 0..n_links {
                if ctx.links.type_[j as usize] != LinkType::Conduit {
                    continue;
                }
                let cr = ctx.link_subtypes.conduit_row(j);
                if cr >= 0 {
                    ctx.link_subtypes.conduits.seep_rate[cr as usize] /= rain0;
                }
            }
        }
    }

    let inv_len = UCF_INV[LENGTH as usize][usz];
    if inv_len == 1.0 {
        let qcf_us = QCF[ctx.options.flow_units as usize];
        if qcf_us != 1.0 {
            for i in 0..n_nodes {
                if ctx.nodes.type_[i as usize] != NodeType::Divider {
                    continue;
                }
                let r = ctx.node_subtypes.divider_row(i);
                if r >= 0 {
                    ctx.node_subtypes.dividers.cutoff[r as usize] /= qcf_us;
                }
            }
            for j in 0..n_links {
                let uj = j as usize;
                if ctx.links.type_[uj] != LinkType::Conduit {
                    continue;
                }
                ctx.links.q0[uj] /= qcf_us;
                ctx.links.q_limit[uj] /= qcf_us;
            }
        }
        return;
    }

    let len = UCF[LENGTH as usize][usz];
    let len2 = len * len;
    let qcf = QCF[ctx.options.flow_units as usize];
    let _rain = UCF[RAINFALL as usize][usz];

    for i in 0..n_nodes {
        let ui = i as usize;
        ctx.nodes.invert_elev[ui] /= len;
        ctx.nodes.full_depth[ui] /= len;
        ctx.nodes.init_depth[ui] /= len;
        ctx.nodes.sur_depth[ui] /= len;
        ctx.nodes.ponded_area[ui] /= len2;
        ctx.nodes.storage_area[ui] /= len2;
        ctx.nodes.rim_depth[ui] /= len;
        if ctx.nodes.type_[ui] == NodeType::Outfall {
            let r = ctx.node_subtypes.outfall_row(i);
            if r >= 0 && ctx.node_subtypes.outfalls.bc_type[r as usize] == OutfallType::Fixed {
                ctx.node_subtypes.outfalls.param[r as usize] /= len;
            }
        }
        if ctx.nodes.type_[ui] == NodeType::Divider {
            let r = ctx.node_subtypes.divider_row(i);
            if r >= 0 {
                ctx.node_subtypes.dividers.cutoff[r as usize] /= qcf;
            }
        }
    }

    for j in 0..n_links {
        let uj = j as usize;
        let shp = ctx.links.xsect_shape[uj];
        ctx.links.xsect_y_full[uj] /= len;
        if shp != XsectShape::ForceMain {
            ctx.links.xsect_w_max[uj] /= len;
        }
        if shp == XsectShape::RectTriang || shp == XsectShape::RectRound || shp == XsectShape::Modbaskethandle {
            ctx.links.xsect_y_bot[uj] /= len;
        }

        if ctx.links.type_[uj] == LinkType::Conduit {
            let cr = ctx.link_subtypes.conduit_row(j);
            if cr >= 0 && ctx.link_subtypes.conduits.length[cr as usize] > 0.0 {
                ctx.link_subtypes.conduits.length[cr as usize] /= len;
            }
            ctx.links.q0[uj] /= qcf;
            ctx.links.q_limit[uj] /= qcf;
        }
        ctx.links.offset1[uj] /= len;
        ctx.links.offset2[uj] /= len;
        let wr = ctx.link_subtypes.weir_row(j);
        if wr >= 0 {
            ctx.link_subtypes.weirs.crest_height[wr as usize] /= len;
            ctx.link_subtypes.weirs.road_width[wr as usize] /= len;
        } else {
            let olr = ctx.link_subtypes.outlet_row(j);
            if olr >= 0 {
                ctx.link_subtypes.outlets.crest_height[olr as usize] /= len;
            }
        }
    }

    for s in 0..n_subcatch {
        ctx.subcatches.width[s as usize] /= len;
    }
}

/// Exact inverse of the parse-time display -> internal conversion, for the .inp writer.
pub fn convert_internal_to_display(ctx: &mut SimulationContext) {
    let us = get_unit_system(ctx.options.flow_units as i32);
    let usz = us as usize;

    {
        let rain0 = UCF[RAINFALL as usize][usz];
        if rain0 != 1.0 {
            for j in 0..ctx.n_links() {
                if ctx.links.type_[j as usize] != LinkType::Conduit {
                    continue;
                }
                let cr = ctx.link_subtypes.conduit_row(j);
                if cr >= 0 {
                    ctx.link_subtypes.conduits.seep_rate[cr as usize] *= rain0;
                }
            }
        }
    }

    let len = UCF[LENGTH as usize][usz];
    if len == 1.0 {
        let qcf_us = QCF[ctx.options.flow_units as usize];
        if qcf_us != 1.0 {
            for i in 0..ctx.n_nodes() {
                if ctx.nodes.type_[i as usize] != NodeType::Divider {
                    continue;
                }
                let r = ctx.node_subtypes.divider_row(i);
                if r >= 0 {
                    ctx.node_subtypes.dividers.cutoff[r as usize] *= qcf_us;
                }
            }
            for j in 0..ctx.n_links() {
                let uj = j as usize;
                if ctx.links.type_[uj] != LinkType::Conduit {
                    continue;
                }
                ctx.links.q0[uj] *= qcf_us;
                ctx.links.q_limit[uj] *= qcf_us;
            }
        }
        return;
    }

    let area = len * len;
    let flow = QCF[ctx.options.flow_units as usize];
    let _rain = UCF[RAINFALL as usize][usz];

    let n_nodes = ctx.n_nodes();
    let n_links = ctx.n_links();
    let n_subcatch = ctx.subcatch_names.size();

    for i in 0..n_nodes {
        let ui = i as usize;
        ctx.nodes.invert_elev[ui] *= len;
        ctx.nodes.full_depth[ui] *= len;
        ctx.nodes.init_depth[ui] *= len;
        ctx.nodes.sur_depth[ui] *= len;
        ctx.nodes.ponded_area[ui] *= area;
        ctx.nodes.storage_area[ui] *= area;
        ctx.nodes.rim_depth[ui] *= len;
        if ctx.nodes.type_[ui] == NodeType::Outfall {
            let r = ctx.node_subtypes.outfall_row(i);
            if r >= 0 && ctx.node_subtypes.outfalls.bc_type[r as usize] == OutfallType::Fixed {
                ctx.node_subtypes.outfalls.param[r as usize] *= len;
            }
        }
        if ctx.nodes.type_[ui] == NodeType::Divider {
            let r = ctx.node_subtypes.divider_row(i);
            if r >= 0 {
                ctx.node_subtypes.dividers.cutoff[r as usize] *= flow;
            }
        }
    }

    for j in 0..n_links {
        let uj = j as usize;
        let shp = ctx.links.xsect_shape[uj];
        ctx.links.xsect_y_full[uj] *= len;
        if shp != XsectShape::ForceMain {
            ctx.links.xsect_w_max[uj] *= len;
        }
        if shp == XsectShape::RectTriang || shp == XsectShape::RectRound || shp == XsectShape::Modbaskethandle {
            ctx.links.xsect_y_bot[uj] *= len;
        }

        if ctx.links.type_[uj] == LinkType::Conduit {
            let cr = ctx.link_subtypes.conduit_row(j);
            if cr >= 0 && ctx.link_subtypes.conduits.length[cr as usize] > 0.0 {
                ctx.link_subtypes.conduits.length[cr as usize] *= len;
            }
            ctx.links.q0[uj] *= flow;
            ctx.links.q_limit[uj] *= flow;
        }
        ctx.links.offset1[uj] *= len;
        ctx.links.offset2[uj] *= len;
        let wr = ctx.link_subtypes.weir_row(j);
        if wr >= 0 {
            ctx.link_subtypes.weirs.crest_height[wr as usize] *= len;
            ctx.link_subtypes.weirs.road_width[wr as usize] *= len;
        } else {
            let olr = ctx.link_subtypes.outlet_row(j);
            if olr >= 0 {
                ctx.link_subtypes.outlets.crest_height[olr as usize] *= len;
            }
        }
        let pr = ctx.link_subtypes.pump_row(j);
        if pr >= 0 {
            let upr = pr as usize;
            if ctx.link_subtypes.pumps.startup[upr] > 0.0 {
                ctx.link_subtypes.pumps.startup[upr] *= len;
            }
            if ctx.link_subtypes.pumps.shutoff[upr] > 0.0 {
                ctx.link_subtypes.pumps.shutoff[upr] *= len;
            }
        }
    }

    for s in 0..n_subcatch {
        ctx.subcatches.width[s as usize] *= len;
    }
}

/// True when the context holds parse-time normalisations the .inp writer must undo.
pub fn needs_authored_conversion(ctx: &SimulationContext) -> bool {
    if ctx.options.link_offsets == 1 {
        return true;
    }
    for j in 0..ctx.n_links() {
        if ctx.links.direction[j as usize] < 0 {
            return true;
        }
    }
    for j in 0..ctx.n_links() {
        let uj = j as usize;
        if ctx.links.type_[uj] == LinkType::Conduit
            && ctx.links.xsect_shape[uj] == XsectShape::FilledCircular
            && ctx.links.xsect_y_bot[uj] > 0.0
        {
            return true;
        }
    }
    false
}

/// Undo the adverse-slope reversal on every conduit with direction -1; the count.
pub fn restore_authored_orientation(ctx: &mut SimulationContext) -> i32 {
    let mut n = 0;
    for j in 0..ctx.n_links() {
        let uj = j as usize;
        if ctx.links.type_[uj] != LinkType::Conduit || ctx.links.direction[uj] >= 0 {
            continue;
        }
        let (a, b) = (ctx.links.node1[uj], ctx.links.node2[uj]);
        ctx.links.node1[uj] = b;
        ctx.links.node2[uj] = a;
        let (o1, o2) = (ctx.links.offset1[uj], ctx.links.offset2[uj]);
        ctx.links.offset1[uj] = o2;
        ctx.links.offset2[uj] = o1;
        ctx.links.q0[uj] = -ctx.links.q0[uj];
        ctx.links.direction[uj] = 1;
        let cr = ctx.link_subtypes.conduit_row(j);
        if cr >= 0 {
            let ucr = cr as usize;
            let cd = &mut ctx.link_subtypes.conduits;
            let (li, lo) = (cd.loss_inlet[ucr], cd.loss_outlet[ucr]);
            cd.loss_inlet[ucr] = lo;
            cd.loss_outlet[ucr] = li;
            cd.slope[ucr] = -cd.slope[ucr];
        }
        n += 1;
    }
    n
}

/// Restore authored link data for writing (orientation, sediment bump, ELEVATION offsets).
pub fn convert_internal_to_authored(ctx: &mut SimulationContext) {
    let n_links = ctx.n_links();
    let n_nodes = ctx.n_nodes();

    restore_authored_orientation(ctx);

    {
        let us = get_unit_system(ctx.options.flow_units as i32);
        let len = UCF[LENGTH as usize][us as usize];
        for j in 0..n_links {
            let uj = j as usize;
            if ctx.links.type_[uj] != LinkType::Conduit || ctx.links.xsect_shape[uj] != XsectShape::FilledCircular {
                continue;
            }
            let dz = ctx.links.xsect_y_bot[uj] * len;
            ctx.links.offset1[uj] -= dz;
            ctx.links.offset2[uj] -= dz;
        }
    }

    if ctx.options.link_offsets != 1 {
        return;
    }
    for j in 0..n_links {
        let uj = j as usize;
        let lt = ctx.links.type_[uj];
        if lt == LinkType::Pump {
            continue;
        }
        let n1 = ctx.links.node1[uj];
        let n2 = ctx.links.node2[uj];
        let ok1 = n1 >= 0 && n1 < n_nodes;
        let ok2 = n2 >= 0 && n2 < n_nodes;
        let inv1 = if ok1 { ctx.nodes.invert_elev[n1 as usize] } else { 0.0 };
        let inv2 = if ok2 { ctx.nodes.invert_elev[n2 as usize] } else { 0.0 };

        if lt == LinkType::Conduit {
            if ok1 {
                ctx.links.offset1[uj] += inv1;
            }
            if ok2 {
                ctx.links.offset2[uj] += inv2;
            }
        } else if lt == LinkType::Orifice {
            if ok1 {
                ctx.links.offset1[uj] += inv1;
            }
        } else if lt == LinkType::Weir || lt == LinkType::Outlet {
            if let Some(crest) = weir_or_outlet_crest(ctx, j) {
                if ok1 {
                    let v = crest.get(ctx) + inv1;
                    crest.set(ctx, v);
                }
            }
        }
    }
}

// -------------------------------------------------------------------------
// Virtual junctions and inlet junctions
// -------------------------------------------------------------------------

fn validate_virtual_junctions(ctx: &mut SimulationContext) {
    let n_nodes = ctx.nodes.count();

    let mut any = false;
    for i in 0..n_nodes {
        if ctx.nodes.is_virtual[i as usize] != 0 {
            any = true;
            break;
        }
    }
    if !any {
        return;
    }

    if ctx.options.routing_model != RoutingModel::Dynwave && ctx.options.routing_model != RoutingModel::Fv {
        for i in 0..n_nodes {
            if ctx.nodes.is_virtual[i as usize] != 0 {
                let e = format_error(ERR_VJ_ROUTING_MODEL, ctx.node_names.name_of(i));
                ctx.errors.push(e);
                break;
            }
        }
    }

    for i in 0..n_nodes {
        let ui = i as usize;
        if ctx.nodes.is_virtual[ui] == 0 {
            continue;
        }
        let code = edit::vj_rule_violation(ctx, i);
        if code != 0 {
            let e = format_error(code, ctx.node_names.name_of(i));
            ctx.errors.push(e);
            continue;
        }
        edit::vj_apply_derived_geometry(ctx, i);
    }
}

fn ij_attached_conduits(ctx: &SimulationContext, node_idx: i32, link_a: &mut i32, link_b: &mut i32) -> i32 {
    *link_a = -1;
    *link_b = -1;
    let mut touches = 0;
    for j in 0..ctx.n_links() {
        let uj = j as usize;
        if ctx.links.type_[uj] != LinkType::Conduit {
            continue;
        }
        let ends = [ctx.links.node1[uj], ctx.links.node2[uj]];
        for n in ends {
            if n != node_idx {
                continue;
            }
            if touches == 0 {
                *link_a = j;
            } else if touches == 1 {
                *link_b = j;
            }
            touches += 1;
        }
    }
    touches
}

fn ij_street_of_link(ctx: &SimulationContext, link_idx: i32) -> i32 {
    if link_idx < 0 || link_idx >= ctx.n_links() {
        return -1;
    }
    let uj = link_idx as usize;
    if ctx.links.xsect_shape[uj] != XsectShape::StreetXsect {
        return -1;
    }
    if uj >= ctx.links.pump_curve_name.len() {
        return -1;
    }
    let sname = &ctx.links.pump_curve_name[uj];
    if sname.is_empty() {
        return -1;
    }
    for s in 0..ctx.streets.count() {
        if ieq(&ctx.streets.names[s as usize], sname) {
            return s;
        }
    }
    -1
}

fn ij_is_drop_design(ctx: &SimulationContext, design_idx: i32) -> bool {
    if design_idx < 0 || design_idx >= ctx.inlets.count() {
        return false;
    }
    let t = &ctx.inlets.inlet_type[design_idx as usize];
    t.as_slice() == b"DROP_GRATE" || t.as_slice() == b"DROP_CURB"
}

fn resolve_inlet_references(ctx: &mut SimulationContext) {
    for i in 0..ctx.inlets.count() {
        let ui = i as usize;
        if ctx.inlets.inlet_type[ui].as_slice() != b"CUSTOM" {
            continue;
        }
        let c = ctx.find_curve(&ctx.inlets.curve_id[ui]);
        if c < 0 {
            let e = format_error(ERR_NAME, &ctx.inlets.curve_id[ui]);
            ctx.errors.push(e);
            continue;
        }
        let tt = ctx.tables.tables[c as usize].type_;
        if tt == TableType::CurveDiversion {
            ctx.inlets.curve_kind[ui] = 1;
        } else if tt == TableType::CurveRating {
            ctx.inlets.curve_kind[ui] = 2;
        } else {
            let e = format_error_detail(ERR_CURVE_SEQUENCE, &ctx.inlets.curve_id[ui], b"(a custom inlet needs a DIVERSION or RATING curve)");
            ctx.errors.push(e);
            continue;
        }
        ctx.inlets.curve_index[ui] = c;
    }

    for r in 0..ctx.inlet_usages.count() {
        let ur = r as usize;
        let host_node = ctx.inlet_usages.node_host[ur];

        if host_node >= 0 {
            let node_name = ctx.node_names.name_of(host_node).to_vec();
            if !ctx.inlet_usages.pending_design_name[ur].is_empty() {
                let mut di = -1;
                for i in 0..ctx.inlets.count() {
                    if ieq(&ctx.inlets.names[i as usize], &ctx.inlet_usages.pending_design_name[ur]) {
                        di = i;
                        break;
                    }
                }
                if di < 0 {
                    let e = format_error_detail(ERR_IJ_DESIGN, &node_name, &cat(&[b"('", &ctx.inlet_usages.pending_design_name[ur], b"')"]));
                    ctx.errors.push(e);
                }
                ctx.inlet_usages.design_index[ur] = di;
                ctx.inlet_usages.pending_design_name[ur].clear();
            }
            if !ctx.inlet_usages.pending_capture_name[ur].is_empty() {
                let ci = ctx.node_names.find(&ctx.inlet_usages.pending_capture_name[ur]);
                let bad = ci < 0 || ci == host_node || ctx.nodes.is_virtual[ci as usize] != 0;
                if bad {
                    let e = format_error_detail(ERR_IJ_CAPTURE_NODE, &node_name, &cat(&[b"('", &ctx.inlet_usages.pending_capture_name[ur], b"')"]));
                    ctx.errors.push(e);
                }
                ctx.inlet_usages.node_index[ur] = if bad { -1 } else { ci };
                ctx.inlet_usages.pending_capture_name[ur].clear();
            }
        } else if !ctx.inlet_usages.pending_capture_name[ur].is_empty() {
            let ci = ctx.node_names.find(&ctx.inlet_usages.pending_capture_name[ur]);
            if ci < 0 {
                let e = format_error(ERR_NAME, &ctx.inlet_usages.pending_capture_name[ur]);
                ctx.errors.push(e);
            }
            ctx.inlet_usages.node_index[ur] = ci;
            ctx.inlet_usages.pending_capture_name[ur].clear();
        }

        let host_link = if host_node >= 0 { edit::ij_host_conduit(ctx, 1, host_node) } else { ctx.inlet_usages.link_index[ur] };
        ctx.inlet_usages.street_index[ur] = ij_street_of_link(ctx, host_link);

        let design = ctx.inlet_usages.design_index[ur];
        if host_node < 0 && design >= 0 && !edit::ij_usage_shape_ok(ctx, design, host_link) {
            let host: Bytes = if host_link >= 0 { cat(&[b"Link ", ctx.link_names.name_of(host_link)]) } else { b"an unknown host".to_vec() };
            let e = format_error_detail(ERR_INLET_USAGE_SHAPE, &ctx.inlets.names[design as usize], &cat(&[b"(", &host, b"; the inlet is ignored)"]));
            ctx.warnings.push(e);
            ctx.inlet_usages.design_index[ur] = -1;
        }
    }
}

fn validate_inlet_junctions(ctx: &mut SimulationContext) {
    let n_nodes = ctx.nodes.count();
    for i in 0..n_nodes {
        let ui = i as usize;
        if ui >= ctx.nodes.is_inlet.len() || ctx.nodes.is_inlet[ui] == 0 {
            continue;
        }
        let name = ctx.node_names.name_of(i).to_vec();

        if ctx.options.routing_model == RoutingModel::Fv {
            ctx.errors.push(format_error_detail(
                ERR_VJ_ROUTING_MODEL,
                &name,
                b"(an inlet junction needs DYNWAVE: the finite-volume mesh splices virtual junctions out and would drop the inlet)",
            ));
            continue;
        }

        let row = ctx.inlet_usages.find_by_node_host(i);
        if row < 0 || ctx.inlet_usages.design_index[row as usize] < 0 {
            ctx.errors.push(format_error(ERR_IJ_NO_USAGE, &name));
            continue;
        }

        let mut j1 = -1;
        let mut j2 = -1;
        if ij_attached_conduits(ctx, i, &mut j1, &mut j2) != 2 {
            continue;
        }

        let design = ctx.inlet_usages.design_index[row as usize];
        if !edit::ij_usage_shape_ok(ctx, design, j1) || !edit::ij_usage_shape_ok(ctx, design, j2) {
            let detail: &[u8] = if ij_is_drop_design(ctx, design) { b"(a drop inlet needs RECT_OPEN or TRAPEZOIDAL conduits)" } else { b"" };
            ctx.errors.push(format_error_detail(ERR_IJ_NOT_STREET, &name, detail));
        }

        if ctx.inlet_usages.find_by_link(j1) >= 0 || ctx.inlet_usages.find_by_link(j2) >= 0 {
            ctx.errors.push(format_error(ERR_IJ_USAGE_ON_PAIR, &name));
        }
    }
}

/// The walk of the virtual-junction seeding: from `start_vj` along `via_link` through
/// virtual junctions to the first real node; (head, distance, wet) when one is reached.
#[allow(clippy::too_many_arguments)]
fn vj_walk(
    ctx: &SimulationContext,
    inc: &[[i32; 2]],
    inc_n: &[i32],
    n_nodes: i32,
    n_links: i32,
    start_vj: i32,
    via_link: i32,
    head_out: &mut f64,
    dist_out: &mut f64,
    wet_out: &mut bool,
) -> bool {
    let mut node = start_vj;
    let mut link = via_link;
    let mut dist = 0.0;
    let mut guard = 0;
    while guard <= n_links {
        let ul = link as usize;
        let cr = ctx.link_subtypes.conduit_row(link);
        dist += if cr >= 0 { ctx.link_subtypes.conduits.length[cr as usize] } else { 0.0 };
        let other = if ctx.links.node1[ul] == node { ctx.links.node2[ul] } else { ctx.links.node1[ul] };
        if other < 0 || other >= n_nodes {
            return false;
        }
        let uo = other as usize;
        if ctx.nodes.is_virtual[uo] == 0 {
            *head_out = ctx.nodes.head[uo];
            *dist_out = dist;
            *wet_out = ctx.nodes.depth[uo] > 0.0;
            return true;
        }
        if inc_n[uo] != 2 {
            return false;
        }
        link = if inc[uo][0] == link { inc[uo][1] } else { inc[uo][0] };
        if link < 0 {
            return false;
        }
        node = other;
        guard += 1;
    }
    false
}

// -------------------------------------------------------------------------
// resolve_cross_references
// -------------------------------------------------------------------------

/// Resolve every cross-reference of a parsed context, convert its inputs to internal
/// units and derive its geometry; problems are appended to ctx.errors / ctx.warnings.
pub fn resolve_cross_references(ctx: &mut SimulationContext) {
    let n_nodes = ctx.node_names.size();
    let n_links = ctx.link_names.size();
    let n_subcatch = ctx.subcatch_names.size();
    let n_gages = ctx.gage_names.size();
    let n_polluts = ctx.pollutant_names.size();

    if ctx.nodes.count() != n_nodes {
        ctx.nodes.resize(n_nodes);
    }
    if ctx.links.count() != n_links {
        ctx.links.resize(n_links);
    }
    if ctx.subcatches.count() != n_subcatch {
        ctx.subcatches.resize(n_subcatch);
    }
    if ctx.gages.count() != n_gages {
        ctx.gages.resize(n_gages);
    }

    if n_polluts > 0 {
        if ctx.pollutants.n_pollutants() != n_polluts {
            ctx.pollutants.resize_pollutants(n_polluts);
        }
        ctx.nodes.resize_quality(n_polluts);
        ctx.links.resize_quality(n_polluts);
        if ctx.subcatches.conc_n_pollutants != n_polluts || ctx.subcatches.conc.len() as i32 != n_subcatch.wrapping_mul(n_polluts) {
            ctx.subcatches.resize_quality(n_polluts);
        }
        ctx.nodes.resize_loads(n_polluts);
        ctx.links.resize_loads(n_polluts);
    }

    if !ctx.gpkg_units_internal {
        convert_inputs_to_internal(ctx, n_nodes, n_links, n_subcatch);
    }

    // Spatial coordinate arrays
    {
        let sp = &mut ctx.spatial;
        let un = n_nodes as usize;
        if sp.node_x.len() < un {
            sp.node_x.resize(un, 0.0);
        }
        if sp.node_y.len() < un {
            sp.node_y.resize(un, 0.0);
        }
        let ul = n_links as usize;
        if sp.link_x.len() < ul {
            sp.link_x.resize(ul, 0.0);
        }
        if sp.link_y.len() < ul {
            sp.link_y.resize(ul, 0.0);
        }
        let us = n_subcatch as usize;
        if sp.subcatch_x.len() < us {
            sp.subcatch_x.resize(us, 0.0);
        }
        if sp.subcatch_y.len() < us {
            sp.subcatch_y.resize(us, 0.0);
        }
        let ug = n_gages as usize;
        if sp.gage_x.len() < ug {
            sp.gage_x.resize(ug, 0.0);
        }
        if sp.gage_y.len() < ug {
            sp.gage_y.resize(ug, 0.0);
        }
        if sp.link_vertices_x.len() < ul {
            sp.link_vertices_x.resize(ul, Vec::new());
        }
        if sp.link_vertices_y.len() < ul {
            sp.link_vertices_y.resize(ul, Vec::new());
        }
        if sp.subcatch_polygon_x.len() < us {
            sp.subcatch_polygon_x.resize(us, Vec::new());
        }
        if sp.subcatch_polygon_y.len() < us {
            sp.subcatch_polygon_y.resize(us, Vec::new());
        }
    }

    // Slice IO-3
    let inp_dir = parent_dir(&ctx.inp_file_path);
    resolve_external_file_slots(ctx, &inp_dir);
    {
        let mut errs: Vec<Bytes> = Vec::new();
        for row in ctx.snowpacks.mesh_config.defaults.values() {
            if !ieq(&row.pack, b"NONE") && ctx.snowpack_names.find(&row.pack) < 0 {
                errs.push(cat(&[b"[2D_SNOW] Unknown snowpack: ", &row.pack]));
            }
        }
        for row in ctx.snowpacks.mesh_config.cells.values() {
            if !ieq(&row.pack, b"NONE") && ctx.snowpack_names.find(&row.pack) < 0 {
                errs.push(cat(&[b"[2D_SNOW] Unknown snowpack: ", &row.pack]));
            }
        }
        ctx.errors.extend(errs);
    }
    {
        let error = ctx.snowpacks.config.parameters.validate();
        if !error.is_empty() {
            ctx.errors.push(cat(&[b"[UEB_OPTIONS] ", &error]));
        }
    }
    if !ctx.snowpacks.config.met_file.is_empty() {
        let path = ctx.snowpacks.config.met_file.absolute.clone();
        let error = ctx.snowpacks.meteorology.load(&path);
        if !error.is_empty() {
            ctx.errors.push(error);
        }
    }

    // External FILE timeseries and rain gages, sharing one parse-once cache.
    {
        let mut series_file_cache = MultiColumnFileCache::default();
        external_files::load_external_timeseries_files(ctx, &inp_dir, &mut series_file_cache);
        external_files::load_external_rain_files_impl(ctx, &mut series_file_cache);
    }

    // Timeseries date offset resolution
    for tbl in ctx.tables.tables.iter_mut() {
        if tbl.type_ != TableType::Timeseries {
            continue;
        }
        if tbl.n_relative <= 0 {
            continue;
        }
        let start = ctx.options.start_date;
        let delta = start - tbl.rel_anchor;
        if delta == 0.0 {
            continue;
        }
        let n = tbl.x.len().min(tbl.n_relative as usize);
        for k in 0..n {
            tbl.x[k] += delta;
        }
        tbl.rel_anchor = start;
    }

    // Gage timeseries re-resolution
    for g in 0..n_gages {
        let ug = g as usize;
        if ctx.gages.ts_index[ug] < 0 && ctx.gages.source[ug] == RainSource::Timeseries && !ctx.gages.ts_name[ug].is_empty() {
            ctx.gages.ts_index[ug] = ctx.find_timeseries(&ctx.gages.ts_name[ug]);
        }
    }

    // Subcatchment outlet re-resolution
    for s in 0..n_subcatch {
        let us = s as usize;
        if us >= ctx.subcatches.outlet_name.len() {
            continue;
        }
        let name = ctx.subcatches.outlet_name[us].clone();
        if name.is_empty() {
            continue;
        }
        if ctx.node_names.find(&name) >= 0 && ctx.subcatch_names.find(&name) >= 0 {
            continue;
        }
        if ctx.subcatches.outlet_node[us] >= 0 && ctx.subcatches.outlet_node[us] < n_nodes {
            continue;
        }
        if ctx.subcatches.outlet_subcatch[us] >= 0 && ctx.subcatches.outlet_subcatch[us] < n_subcatch {
            continue;
        }
        let node_idx = ctx.node_names.find(&name);
        if node_idx >= 0 {
            ctx.subcatches.outlet_node[us] = node_idx;
            ctx.subcatches.outlet_subcatch[us] = -1;
        } else {
            let sub_idx = ctx.subcatch_names.find(&name);
            if sub_idx >= 0 {
                ctx.subcatches.outlet_subcatch[us] = sub_idx;
                ctx.subcatches.outlet_node[us] = -1;
            } else {
                ctx.errors.push(format_error(ERR_NAME, &name));
            }
        }
    }

    // Subcatchment groundwater receiving-node / aquifer re-resolution
    {
        let pend = std::mem::take(&mut ctx.pending_gw_aquifers);
        for (si, nm) in &pend {
            let si = *si;
            if si < 0 || si >= n_subcatch {
                continue;
            }
            let us = si as usize;
            if us >= ctx.subcatches.gw_aquifer.len() || ctx.subcatches.gw_aquifer[us] >= 0 {
                continue;
            }
            ctx.subcatches.gw_aquifer[us] = ctx.aquifer_names.find(nm);
            if ctx.subcatches.gw_aquifer[us] < 0 {
                ctx.errors.push(format_error(ERR_NAME, nm));
            }
        }
        ctx.pending_gw_aquifers.clear();
        let pend = std::mem::take(&mut ctx.pending_gw_nodes);
        for (si, nm) in &pend {
            let si = *si;
            if si < 0 || si >= n_subcatch {
                continue;
            }
            let us = si as usize;
            if us >= ctx.subcatches.gw_node.len() || ctx.subcatches.gw_node[us] >= 0 {
                continue;
            }
            ctx.subcatches.gw_node[us] = ctx.node_names.find(nm);
            if ctx.subcatches.gw_node[us] < 0 {
                ctx.errors.push(format_error(ERR_NAME, nm));
            }
        }
        ctx.pending_gw_nodes.clear();
    }

    // Subcatchment snow pack resolution
    for s in 0..n_subcatch {
        let us = s as usize;
        if us >= ctx.subcatches.snowpack_name.len() {
            continue;
        }
        if ctx.subcatches.snowpack_name[us].is_empty() {
            continue;
        }
        if ctx.subcatches.snowpack[us] >= 0 {
            continue;
        }
        ctx.subcatches.snowpack[us] = ctx.snowpack_names.find(&ctx.subcatches.snowpack_name[us]);
        if ctx.subcatches.snowpack[us] < 0 {
            let e = format_error(ERR_NAME, &ctx.subcatches.snowpack_name[us]);
            ctx.errors.push(e);
        }
    }

    // Subcatchment gage re-resolution
    for s in 0..n_subcatch {
        let us = s as usize;
        let gi = ctx.subcatches.gage[us];
        if gi >= 0 && gi < n_gages {
            continue;
        }
        let name: Bytes = if us < ctx.subcatches.gage_name.len() { ctx.subcatches.gage_name[us].clone() } else { Bytes::new() };
        if name.is_empty() {
            ctx.subcatches.gage[us] = -1;
            continue;
        }
        let idx = ctx.gage_names.find(&name);
        if idx >= 0 && idx < n_gages {
            ctx.subcatches.gage[us] = idx;
        } else {
            ctx.subcatches.gage[us] = -1;
            ctx.errors.push(format_error(ERR_NAME, &name));
        }
    }

    // Gap #53: co-gage detection (legacy gage_validate)
    {
        let mut used = vec![0u8; n_gages.max(0) as usize];
        for s2 in 0..n_subcatch {
            let gi = ctx.subcatches.gage[s2 as usize];
            if gi >= 0 && gi < n_gages {
                used[gi as usize] = 1;
            }
        }
        for gname in &ctx.unit_hyds.gage_names {
            let gi = ctx.gage_names.find(gname);
            if gi >= 0 && gi < n_gages {
                used[gi as usize] = 1;
            }
        }
        let mut first_on_ts: std::collections::HashMap<i32, i32> = std::collections::HashMap::new();
        for gj in 0..n_gages {
            let ugj = gj as usize;
            ctx.gages.co_gage_index[ugj] = -1;
            if used[ugj] == 0 {
                continue;
            }
            if ctx.gages.source[ugj] != RainSource::Timeseries {
                continue;
            }
            let ts_j = ctx.gages.ts_index[ugj];
            if ts_j < 0 {
                continue;
            }
            match first_on_ts.get(&ts_j) {
                Some(&first) => ctx.gages.co_gage_index[ugj] = first,
                None => {
                    first_on_ts.insert(ts_j, gj);
                }
            }
        }
    }

    // Storage curve name resolution
    for i in 0..n_nodes {
        if ctx.nodes.type_[i as usize] != NodeType::Storage {
            continue;
        }
        let r = ctx.node_subtypes.storage_row(i);
        if r < 0 {
            continue;
        }
        let ur = r as usize;
        if ctx.node_subtypes.storages.curve[ur] < 0 && !ctx.node_subtypes.storages.curve_name[ur].is_empty() {
            let c = ctx.find_curve(&ctx.node_subtypes.storages.curve_name[ur]);
            ctx.node_subtypes.storages.curve[ur] = c;
            if c < 0 {
                let e = format_error(ERR_NAME, &ctx.node_subtypes.storages.curve_name[ur]);
                ctx.errors.push(e);
            }
        }
    }

    // Outfall stage-data name resolution
    for i in 0..n_nodes {
        if ctx.nodes.type_[i as usize] != NodeType::Outfall {
            continue;
        }
        let r = ctx.node_subtypes.outfall_row(i);
        if r < 0 {
            continue;
        }
        let ur = r as usize;
        let bc = ctx.node_subtypes.outfalls.bc_type[ur];
        let is_tidal = bc == OutfallType::Tidal;
        if !is_tidal && bc != OutfallType::Timeseries {
            continue;
        }
        if ctx.node_subtypes.outfalls.param[ur] >= 0.0 {
            continue;
        }
        let nm = ctx.node_subtypes.outfalls.param_name[ur].clone();
        let t = if nm.is_empty() {
            -1
        } else if is_tidal {
            ctx.find_curve(&nm)
        } else {
            ctx.find_timeseries(&nm)
        };
        if t < 0 {
            ctx.errors.push(format_error(ERR_NAME, &nm));
            continue;
        }
        let is_ts = ctx.tables.tables[t as usize].type_ == TableType::Timeseries;
        if is_tidal == is_ts {
            ctx.errors.push(format_error(ERR_NAME, &nm));
            continue;
        }
        ctx.node_subtypes.outfalls.param[ur] = t as f64;
    }

    // Pump curve name resolution
    for j in 0..n_links {
        let uj = j as usize;
        if ctx.links.type_[uj] != LinkType::Pump {
            continue;
        }
        let pr = ctx.link_subtypes.pump_row(j);
        if pr < 0 {
            continue;
        }
        let upr = pr as usize;
        if ctx.link_subtypes.pumps.curve[upr] >= 0 {
            continue;
        }
        if ctx.links.pump_curve_name[uj].is_empty() || ctx.links.pump_curve_name[uj].as_slice() == b"*" {
            continue;
        }
        let c = ctx.find_curve(&ctx.links.pump_curve_name[uj]);
        ctx.link_subtypes.pumps.curve[upr] = c;
        if c < 0 {
            let e = format_error(ERR_NAME, &ctx.links.pump_curve_name[uj]);
            ctx.errors.push(e);
        }
    }

    // Outlet (TABULAR) rating curve name resolution
    for j in 0..n_links {
        let uj = j as usize;
        if ctx.links.type_[uj] != LinkType::Outlet {
            continue;
        }
        let olr = ctx.link_subtypes.outlet_row(j);
        if olr < 0 {
            continue;
        }
        let uolr = olr as usize;
        let outlet_type = cvt_i32(ctx.link_subtypes.outlets.outlet_type[uolr]);
        if outlet_type < 2 {
            continue;
        }
        if ctx.link_subtypes.outlets.curve[uolr] >= 0 {
            continue;
        }
        if ctx.links.pump_curve_name[uj].is_empty() {
            continue;
        }
        let c = ctx.find_curve(&ctx.links.pump_curve_name[uj]);
        ctx.link_subtypes.outlets.curve[uolr] = c;
        if c < 0 {
            let e = format_error(ERR_NAME, &ctx.links.pump_curve_name[uj]);
            ctx.errors.push(e);
        }
    }

    // Weir discharge-coefficient curve
    for j in 0..n_links {
        let uj = j as usize;
        if ctx.links.type_[uj] != LinkType::Weir {
            continue;
        }
        let wr = ctx.link_subtypes.weir_row(j);
        if wr < 0 {
            continue;
        }
        let uwr = wr as usize;
        if ctx.link_subtypes.weirs.cd_curve[uwr] >= 0 {
            continue;
        }
        if ctx.links.pump_curve_name[uj].is_empty() || ctx.links.pump_curve_name[uj].as_slice() == b"*" {
            continue;
        }
        let c = ctx.find_curve(&ctx.links.pump_curve_name[uj]);
        ctx.link_subtypes.weirs.cd_curve[uwr] = c;
        if c < 0 {
            let e = format_error(ERR_NAME, &ctx.links.pump_curve_name[uj]);
            ctx.errors.push(e);
        }
    }

    // Pump startup/shutoff depths: display -> internal (ft)
    {
        let us = get_unit_system(ctx.options.flow_units as i32);
        let ucf_len = UCF[LENGTH as usize][us as usize];
        for j in 0..n_links {
            let uj = j as usize;
            if ctx.links.type_[uj] != LinkType::Pump {
                continue;
            }
            let pr = ctx.link_subtypes.pump_row(j);
            if pr < 0 {
                continue;
            }
            let upr = pr as usize;
            if ctx.link_subtypes.pumps.startup[upr] > 0.0 {
                ctx.link_subtypes.pumps.startup[upr] /= ucf_len;
            }
            if ctx.link_subtypes.pumps.shutoff[upr] > 0.0 {
                ctx.link_subtypes.pumps.shutoff[upr] /= ucf_len;
            }
        }
    }

    // Node init_depth -> head initialisation
    for i in 0..n_nodes as usize {
        ctx.nodes.head[i] = ctx.nodes.invert_elev[i] + ctx.nodes.init_depth[i];
        ctx.nodes.depth[i] = ctx.nodes.init_depth[i];
    }

    // [INFLOWS] / [DWF] / [RDII] node re-resolution (backwards, erasing unresolved rows)
    {
        let mut i = ctx.ext_inflows.count() - 1;
        while i >= 0 {
            let ui = i as usize;
            if ctx.ext_inflows.node_idx[ui] < 0 {
                let name = ctx.ext_inflows.node_name[ui].clone();
                let mut keep = false;
                if !name.is_empty() {
                    ctx.ext_inflows.node_idx[ui] = ctx.node_names.find(&name);
                    if ctx.ext_inflows.node_idx[ui] >= 0 {
                        keep = true;
                    } else {
                        ctx.errors.push(format_error(ERR_NAME, &name));
                    }
                }
                if !keep {
                    ctx.ext_inflows.erase(i);
                }
            }
            i -= 1;
        }
        let mut i = ctx.dwf_inflows.count() - 1;
        while i >= 0 {
            let ui = i as usize;
            if ctx.dwf_inflows.node_idx[ui] < 0 {
                let name = ctx.dwf_inflows.node_name[ui].clone();
                let mut keep = false;
                if !name.is_empty() {
                    ctx.dwf_inflows.node_idx[ui] = ctx.node_names.find(&name);
                    if ctx.dwf_inflows.node_idx[ui] >= 0 {
                        keep = true;
                    } else {
                        ctx.errors.push(format_error(ERR_NAME, &name));
                    }
                }
                if !keep {
                    ctx.dwf_inflows.erase(i);
                }
            }
            i -= 1;
        }
        let mut i = ctx.rdii_assigns.count() - 1;
        while i >= 0 {
            let ui = i as usize;
            if ctx.rdii_assigns.node_idx[ui] < 0 {
                let name = ctx.rdii_assigns.node_name[ui].clone();
                let mut keep = false;
                if !name.is_empty() {
                    ctx.rdii_assigns.node_idx[ui] = ctx.node_names.find(&name);
                    if ctx.rdii_assigns.node_idx[ui] >= 0 {
                        keep = true;
                    } else {
                        ctx.errors.push(format_error(ERR_NAME, &name));
                    }
                }
                if !keep {
                    ctx.rdii_assigns.erase(i);
                }
            }
            i -= 1;
        }
    }

    // [INFLOWS] baseline / [DWF] time pattern names
    {
        let mut errs: Vec<Bytes> = Vec::new();
        let names = &ctx.patterns.names;
        let mut check = |name: &[u8]| {
            if !name.is_empty() && find_ci(names, name) < 0 {
                errs.push(format_error(ERR_NAME, name));
            }
        };
        for name in &ctx.ext_inflows.pattern_name {
            check(name);
        }
        for i in 0..ctx.dwf_inflows.count() as usize {
            check(&ctx.dwf_inflows.pat1[i]);
            check(&ctx.dwf_inflows.pat2[i]);
            check(&ctx.dwf_inflows.pat3[i]);
            check(&ctx.dwf_inflows.pat4[i]);
        }
        ctx.errors.extend(errs);
    }

    // [INITIAL_QUALITY] element + constituent resolution
    resolve_initial_quality(ctx);

    // Link end-node re-resolution
    {
        let pend = std::mem::take(&mut ctx.pending_link_nodes);
        for (li, names) in &pend {
            let li = *li;
            if li < 0 || li >= ctx.links.count() {
                continue;
            }
            let ul = li as usize;
            for (end, name) in [(0, &names.0), (1, &names.1)] {
                let idx = if end == 0 { ctx.links.node1[ul] } else { ctx.links.node2[ul] };
                if idx >= 0 && idx < n_nodes {
                    continue;
                }
                let found = ctx.node_names.find(name);
                if end == 0 {
                    ctx.links.node1[ul] = found;
                } else {
                    ctx.links.node2[ul] = found;
                }
                if found < 0 {
                    ctx.errors.push(format_error(ERR_NAME, name));
                }
            }
        }
        ctx.pending_link_nodes.clear();
    }

    // Quality data matrix sizing
    let n_landuses = ctx.landuse_names.size();
    if n_landuses > 0 && n_polluts > 0 {
        if ctx.buildup.n_landuses == 0 {
            ctx.buildup.resize(n_landuses, n_polluts);
        }
        if ctx.washoff.n_landuses == 0 {
            ctx.washoff.resize(n_landuses, n_polluts);
        }
    }
    if n_polluts > 0 && n_nodes > 0 && ctx.treatment.n_nodes == 0 {
        ctx.treatment.resize(n_nodes, n_polluts);
    }

    // Subcatchment coverage matrix sizing
    if n_landuses > 0 && n_subcatch > 0 {
        let total = n_subcatch.wrapping_mul(n_landuses) as usize;
        if ctx.subcatches.coverage.len() < total {
            ctx.subcatches.coverage.resize(total, 0.0);
        }
        ctx.subcatches.coverage_n_landuses = n_landuses;
    }

    // Divider link and curve re-resolution
    for i in 0..n_nodes {
        if ctx.nodes.type_[i as usize] != NodeType::Divider {
            continue;
        }
        let r = ctx.node_subtypes.divider_row(i);
        if r < 0 {
            continue;
        }
        let ur = r as usize;
        let mut bad_name = false;
        {
            let d = &ctx.node_subtypes.dividers;
            if d.link[ur] < 0 && !d.link_name[ur].is_empty() && d.link_name[ur].as_slice() != b"*" {
                let found = ctx.link_names.find(&d.link_name[ur]);
                let lname = d.link_name[ur].clone();
                ctx.node_subtypes.dividers.link[ur] = found;
                if found < 0 {
                    ctx.errors.push(format_error(ERR_NAME, &lname));
                    bad_name = true;
                }
            }
        }
        let dl = ctx.node_subtypes.dividers.link[ur];
        if !bad_name
            && (dl < 0 || dl >= n_links || (ctx.links.node1[dl as usize] != i && ctx.links.node2[dl as usize] != i))
        {
            let e = format_error(ERR_DIVIDER_LINK, ctx.node_names.name_of(i));
            ctx.errors.push(e);
        }
        if ctx.node_subtypes.dividers.curve[ur] < 0 && !ctx.node_subtypes.dividers.curve_name[ur].is_empty() {
            let c = ctx.find_curve(&ctx.node_subtypes.dividers.curve_name[ur]);
            ctx.node_subtypes.dividers.curve[ur] = c;
        }
    }

    // Node degree
    for j in 0..n_links as usize {
        let n1 = ctx.links.node1[j];
        let n2 = ctx.links.node2[j];
        if n1 >= 0 && n1 < n_nodes {
            ctx.nodes.degree[n1 as usize] += 1;
        }
        if n2 >= 0 && n2 < n_nodes {
            ctx.nodes.degree[n2 as usize] += 1;
        }
    }

    // Virtual-junction initial-state seeding (DYNWAVE only)
    if ctx.options.routing_model == RoutingModel::Dynwave {
        let nn = n_nodes.max(0) as usize;
        let mut inc: Vec<[i32; 2]> = vec![[-1, -1]; nn];
        let mut inc_n: Vec<i32> = vec![0; nn];
        for j in 0..n_links {
            let uj = j as usize;
            if ctx.links.type_[uj] != LinkType::Conduit {
                continue;
            }
            for nd in [ctx.links.node1[uj], ctx.links.node2[uj]] {
                if nd < 0 || nd >= n_nodes {
                    continue;
                }
                let und = nd as usize;
                if ctx.nodes.is_virtual[und] == 0 {
                    continue;
                }
                if inc_n[und] < 2 {
                    inc[und][inc_n[und] as usize] = j;
                }
                inc_n[und] += 1;
            }
        }
        for i in 0..n_nodes {
            let ui = i as usize;
            if ctx.nodes.is_virtual[ui] == 0 || inc_n[ui] != 2 {
                continue;
            }
            let (mut h_a, mut d_a, mut h_b, mut d_b) = (0.0, 0.0, 0.0, 0.0);
            let (mut wet_a, mut wet_b) = (false, false);
            let ok_a = vj_walk(ctx, &inc, &inc_n, n_nodes, n_links, i, inc[ui][0], &mut h_a, &mut d_a, &mut wet_a);
            let ok_b = vj_walk(ctx, &inc, &inc_n, n_nodes, n_links, i, inc[ui][1], &mut h_b, &mut d_b, &mut wet_b);
            if !(ok_a && wet_a) && !(ok_b && wet_b) {
                continue;
            }
            let head;
            if ok_a && ok_b {
                let dt_sum = d_a + d_b;
                head = if dt_sum > 0.0 { (h_a * d_b + h_b * d_a) / dt_sum } else { 0.5 * (h_a + h_b) };
            } else if ok_a {
                head = h_a;
            } else if ok_b {
                head = h_b;
            } else {
                continue;
            }
            let depth = head - ctx.nodes.invert_elev[ui];
            if depth > 0.0 {
                ctx.nodes.init_depth[ui] = depth;
                ctx.nodes.depth[ui] = depth;
                ctx.nodes.head[ui] = head;
            }
        }
    }

    // Non-conduit link offset conversion (ELEVATION -> DEPTH), or DEPTH-mode zeroing
    if ctx.options.link_offsets == 1 && !ctx.gpkg_units_internal {
        for j in 0..n_links {
            let uj = j as usize;
            let lt = ctx.links.type_[uj];
            if lt == LinkType::Conduit || lt == LinkType::Pump {
                continue;
            }
            let n1 = ctx.links.node1[uj];
            if lt == LinkType::Orifice {
                if n1 >= 0 && n1 < n_nodes {
                    let name = ctx.link_names.name_of(j).to_vec();
                    ctx.links.offset1[uj] =
                        legacy_offset_height(ctx.links.offset1[uj], ctx.nodes.invert_elev[n1 as usize], &mut ctx.warnings, &name);
                }
            } else if lt == LinkType::Weir || lt == LinkType::Outlet {
                if let Some(crest) = weir_or_outlet_crest(ctx, j) {
                    if n1 >= 0 && n1 < n_nodes {
                        let raw_crest = crest.get(ctx);
                        let name = ctx.link_names.name_of(j).to_vec();
                        let v = legacy_offset_height(raw_crest, ctx.nodes.invert_elev[n1 as usize], &mut ctx.warnings, &name);
                        crest.set(ctx, v);
                    }
                }
            }
        }
    } else {
        for j in 0..n_links {
            let uj = j as usize;
            let lt = ctx.links.type_[uj];
            if lt == LinkType::Conduit {
                if ctx.links.offset1[uj] < 0.0 {
                    let w = format_warning(WARN_NEGATIVE_OFFSET, ctx.link_names.name_of(j));
                    ctx.warnings.push(w);
                    ctx.links.offset1[uj] = 0.0;
                }
                if ctx.links.offset2[uj] < 0.0 {
                    let w = format_warning(WARN_NEGATIVE_OFFSET, ctx.link_names.name_of(j));
                    ctx.warnings.push(w);
                    ctx.links.offset2[uj] = 0.0;
                }
            } else if lt == LinkType::Weir && ctx.links.offset1[uj] < 0.0 {
                ctx.links.offset1[uj] = 0.0;
            }
        }
    }

    // Regulator crest below downstream invert (legacy WARNING 10)
    for j in 0..n_links {
        let uj = j as usize;
        let lt = ctx.links.type_[uj];
        if lt != LinkType::Orifice && lt != LinkType::Weir && lt != LinkType::Outlet {
            continue;
        }
        let n1 = ctx.links.node1[uj];
        let n2 = ctx.links.node2[uj];
        if n1 < 0 || n1 >= n_nodes || n2 < 0 || n2 >= n_nodes {
            continue;
        }
        let off = if lt == LinkType::Orifice { Some(Crest::Offset1(uj)) } else { weir_or_outlet_crest(ctx, j) };
        let off = match off {
            None => continue,
            Some(o) => o,
        };
        let inv1 = ctx.nodes.invert_elev[n1 as usize];
        let inv2 = ctx.nodes.invert_elev[n2 as usize];
        ctx.links.offset2[uj] = off.get(ctx);
        if (lt == LinkType::Orifice || lt == LinkType::Weir) && off.get(ctx) < 0.0 {
            off.set(ctx, 0.0);
        }
        if inv1 + off.get(ctx) < inv2 {
            if ctx.options.routing_model == RoutingModel::Dynwave || ctx.options.routing_model == RoutingModel::Fv {
                off.set(ctx, inv2 - inv1);
                let w = cat(&[b"  WARNING 10: crest elevation raised to downstream invert for regulator Link ", ctx.link_names.name_of(j)]);
                ctx.warnings.push(w);
            } else {
                let w = format_warning(WARN_REGULATOR_CREST_LOW, ctx.link_names.name_of(j));
                ctx.warnings.push(w);
            }
        }
    }

    // Transect geometry tables for IRREGULAR cross-sections
    {
        let nt = ctx.transects.count();
        ctx.transect_tables.resize(nt.max(0) as usize, TransectData::default());
        let t_us = get_unit_system(ctx.options.flow_units as i32);
        let t_ucf = UCF[LENGTH as usize][t_us as usize];
        for t in 0..nt {
            let td = &mut ctx.transect_tables[t as usize];
            if !transect::build_from_store(&ctx.transects, t, t_ucf, td) {
                let e = format_error(ERR_TRANSECT_MANNING, &td.name);
                ctx.errors.push(e);
                continue;
            }
        }

        for j in 0..n_links {
            let uj = j as usize;
            if ctx.links.xsect_shape[uj] != XsectShape::Irregular {
                continue;
            }
            let tname = ctx.links.pump_curve_name[uj].clone();
            if !tname.is_empty() {
                let hit = find_ci(&ctx.transects.names, &tname);
                if hit >= 0 {
                    ctx.links.xsect_curve[uj] = hit;
                } else if ctx.links.xsect_curve[uj] < 0 {
                    ctx.errors.push(format_error(ERR_NAME, &tname));
                }
            } else if ctx.links.xsect_curve[uj] < 0 {
                let e = format_error(ERR_NAME, ctx.link_names.name_of(j));
                ctx.errors.push(e);
            }
            let ci = ctx.links.xsect_curve[uj];
            if ci >= 0 && ci < nt {
                let td = &ctx.transect_tables[ci as usize];
                ctx.links.xsect_y_full[uj] = td.y_full;
                ctx.links.xsect_a_full[uj] = td.a_full;
                ctx.links.xsect_r_full[uj] = td.r_full;
                ctx.links.xsect_w_max[uj] = td.w_max;
            }
        }
    }

    // CUSTOM shape curves
    {
        let n_tables = ctx.tables.tables.len() as i32;
        let mut custom_memo: StdMap<(i32, f64), i32> =
            StdMap::new(|a: &(i32, f64), b: &(i32, f64)| a.0 < b.0 || (!(b.0 < a.0) && a.1 < b.1));
        for j in 0..n_links {
            let uj = j as usize;
            if ctx.links.xsect_shape[uj] != XsectShape::Custom {
                continue;
            }
            let cname = ctx.links.pump_curve_name[uj].clone();
            if cname.is_empty() {
                continue;
            }
            let ci = ctx.find_curve(&cname);
            if ci >= 0 && ci < n_tables {
                ctx.links.xsect_curve[uj] = ci;
                let y_full = ctx.links.xsect_y_full[uj];
                if y_full <= 0.0 || ctx.tables.tables[ci as usize].x.len() < 2 {
                    continue;
                }
                if let Some(&memo) = custom_memo.find(&(ci, y_full)) {
                    let shared = &ctx.transect_tables[memo as usize];
                    ctx.links.xsect_a_full[uj] = shared.a_full;
                    ctx.links.xsect_r_full[uj] = shared.r_full;
                    ctx.links.xsect_w_max[uj] = shared.w_max;
                    ctx.links.xsect_curve[uj] = memo;
                    continue;
                }
                let tbl = &ctx.tables.tables[ci as usize];
                let mut w_max_norm = 0.0;
                for &yv in &tbl.y {
                    if yv > w_max_norm {
                        w_max_norm = yv;
                    }
                }
                let w_max = w_max_norm * y_full;
                ctx.links.xsect_w_max[uj] = w_max;

                let mut ctd = TransectData::default();
                ctd.name = cname.clone();
                ctd.y_full = y_full;
                ctd.a_full = 0.0;
                ctd.r_full = 0.0;
                ctd.w_max = w_max;
                transect::build_custom_tables(&mut ctd, y_full, &tbl.x, &tbl.y, tbl.x.len() as i32);
                ctx.links.xsect_a_full[uj] = ctd.a_full;
                ctx.links.xsect_r_full[uj] = ctd.r_full;
                ctx.links.xsect_w_max[uj] = ctd.w_max;
                let custom_idx = ctx.transect_tables.len() as i32;
                ctx.transect_tables.push(ctd);
                ctx.links.xsect_curve[uj] = custom_idx;
                custom_memo.emplace((ci, y_full), custom_idx);
            } else {
                ctx.errors.push(format_error(ERR_NAME, &cname));
            }
        }
    }

    // STREET cross-sections
    {
        let us = get_unit_system(ctx.options.flow_units as i32);
        let ucf_len = UCF[LENGTH as usize][us as usize];
        let mut street_memo: std::collections::HashMap<i32, i32> = std::collections::HashMap::new();
        for j in 0..n_links {
            let uj = j as usize;
            if ctx.links.xsect_shape[uj] != XsectShape::StreetXsect {
                continue;
            }
            let sname = ctx.links.pump_curve_name[uj].clone();
            if sname.is_empty() {
                continue;
            }
            let si = find_ci(&ctx.streets.names, &sname);
            if si < 0 {
                continue;
            }
            let su = si as usize;

            if let Some(&memo) = street_memo.get(&si) {
                let shared = &ctx.transect_tables[memo as usize];
                ctx.links.xsect_curve[uj] = memo;
                ctx.links.xsect_y_full[uj] = shared.y_full;
                ctx.links.xsect_a_full[uj] = shared.a_full;
                ctx.links.xsect_r_full[uj] = shared.r_full;
                ctx.links.xsect_w_max[uj] = shared.w_max;
                continue;
            }

            let st = &ctx.streets;
            let mut sp = StreetParams::default();
            sp.width = st.t_crown[su] / ucf_len;
            sp.curb_height = st.h_curb[su] / ucf_len;
            sp.slope = st.sx[su] / 100.0;
            sp.roughness = st.n_road[su];
            sp.gutter_depression = st.gutter_depres[su] / ucf_len;
            sp.gutter_width = st.gutter_width[su] / ucf_len;
            sp.sides = st.sides[su];
            sp.back_width = st.back_width[su] / ucf_len;
            sp.back_slope = st.back_slope[su] / 100.0;
            sp.back_roughness = st.back_n[su];

            let mut td = TransectData::default();
            td.name = sname.clone();
            street::build_transect(&sp, &mut td);

            let idx_tbl = ctx.transect_tables.len() as i32;
            ctx.transect_tables.push(td);
            street_memo.insert(si, idx_tbl);
            let built = &ctx.transect_tables[idx_tbl as usize];
            ctx.links.xsect_curve[uj] = idx_tbl;
            ctx.links.xsect_y_full[uj] = built.y_full;
            ctx.links.xsect_a_full[uj] = built.a_full;
            ctx.links.xsect_r_full[uj] = built.r_full;
            ctx.links.xsect_w_max[uj] = built.w_max;
        }
    }

    // Full-flow properties from cross-section geometry (legacy xsect_setParams)
    for j in 0..n_links {
        let uj = j as usize;
        let lt = ctx.links.type_[uj];
        if lt != LinkType::Conduit && lt != LinkType::Orifice && lt != LinkType::Weir {
            continue;
        }

        let mut y_full = ctx.links.xsect_y_full[uj];
        let mut w_max = ctx.links.xsect_w_max[uj];
        let shape = ctx.links.xsect_shape[uj];

        let mut a_full;
        let mut r_full;
        let mut s_full;
        let mut s_max;
        let mut yw_max;

        match shape {
            XsectShape::Custom => {
                a_full = ctx.links.xsect_a_full[uj];
                r_full = ctx.links.xsect_r_full[uj];
                w_max = ctx.links.xsect_w_max[uj];
                s_full = a_full * r_full.powf(2.0 / 3.0); // LIBM
                s_max = s_full;
                yw_max = y_full;
                let ci = ctx.links.xsect_curve[uj];
                if ci >= 0 && (ci as usize) < ctx.transect_tables.len() {
                    let td = &ctx.transect_tables[ci as usize];
                    let mut xs = XSectParams::default();
                    xs.type_ = translate_shape(shape);
                    apply_tabulated_xsect_params(&mut xs, td, y_full, ci);
                    a_full = xs.a_full;
                    r_full = xs.r_full;
                    w_max = xs.w_max;
                    s_full = xs.s_full;
                    s_max = xs.s_max;
                    yw_max = xs.yw_max;
                    let a_bot = xs.a_bot;
                    ctx.links.xsect_a_bot[uj] = a_bot;
                }
            }
            XsectShape::Irregular | XsectShape::StreetXsect => {
                let ci = ctx.links.xsect_curve[uj];
                if ci >= 0 && (ci as usize) < ctx.transect_tables.len() {
                    let td = &ctx.transect_tables[ci as usize];
                    let mut xs = XSectParams::default();
                    xs.type_ = translate_shape(shape);
                    apply_tabulated_xsect_params(&mut xs, td, y_full, ci);
                    y_full = xs.y_full;
                    a_full = xs.a_full;
                    r_full = xs.r_full;
                    w_max = xs.w_max;
                    s_full = xs.s_full;
                    s_max = xs.s_max;
                    yw_max = xs.yw_max;
                    let a_bot = xs.a_bot;
                    if shape == XsectShape::Irregular {
                        ctx.links.xsect_a_bot[uj] = a_bot;
                    }
                } else {
                    a_full = w_max * y_full;
                    let p_def = 2.0 * y_full + w_max;
                    r_full = if p_def > 0.0 { a_full / p_def } else { 0.0 };
                    s_full = a_full * r_full.powf(2.0 / 3.0); // LIBM
                    s_max = s_full;
                    yw_max = y_full;
                }
            }
            _ => {
                let us = get_unit_system(ctx.options.flow_units as i32);
                let ucf_len = UCF[LENGTH as usize][us as usize];
                let mut xs = XSectParams::default();
                let p = [ctx.links.xsect_geom1[uj], ctx.links.xsect_geom2[uj], ctx.links.xsect_geom3[uj], ctx.links.xsect_geom4[uj]];
                let rc = xsect::set_params(&mut xs, translate_shape(shape), &p, ucf_len);
                if rc == 0 && xs.a_full > 0.0 {
                    y_full = xs.y_full;
                    w_max = xs.w_max;
                    yw_max = xs.yw_max;
                    a_full = xs.a_full;
                    r_full = xs.r_full;
                    s_full = xs.s_full;
                    s_max = xs.s_max;
                    ctx.links.xsect_y_bot[uj] = xs.y_bot;
                    ctx.links.xsect_a_bot[uj] = xs.a_bot;
                    ctx.links.xsect_s_bot[uj] = xs.s_bot;
                    ctx.links.xsect_r_bot[uj] = xs.r_bot;
                    if shape == XsectShape::ForceMain && ctx.options.force_main_eqn == FrictionModel::DarcyWeisbach as i32 {
                        ctx.links.xsect_r_bot[uj] /= UCF[RAINDEPTH as usize][us as usize];
                    }
                } else {
                    if rc != 0 {
                        let mut geom = b"(Geom1-4: ".to_vec();
                        for (k, v) in p.iter().enumerate() {
                            if k > 0 {
                                geom.push(b' ');
                            }
                            geom.extend_from_slice(&printf_g(*v, 10));
                        }
                        geom.push(b')');
                        let name = cat(&[b"in the cross-section geometry of link ", ctx.link_names.name_of(j)]);
                        ctx.errors.push(format_error_detail(ERR_NUMBER, &name, &geom));
                    }
                    a_full = w_max * y_full;
                    let p_def = 2.0 * y_full + w_max;
                    r_full = if p_def > 0.0 { a_full / p_def } else { 0.0 };
                    s_full = a_full * r_full.powf(2.0 / 3.0); // LIBM
                    s_max = s_full;
                    yw_max = y_full;
                }
            }
        }

        ctx.links.xsect_y_full[uj] = y_full;
        ctx.links.xsect_a_full[uj] = a_full;
        ctx.links.xsect_r_full[uj] = r_full;
        ctx.links.xsect_s_full[uj] = s_full;
        ctx.links.xsect_s_max[uj] = s_max;
        ctx.links.xsect_w_max[uj] = w_max;
        ctx.links.xsect_yw_max[uj] = yw_max;
    }

    // Conduit slope computation (legacy conduit_getSlope)
    for j in 0..n_links {
        let uj = j as usize;
        if ctx.links.type_[uj] != LinkType::Conduit {
            continue;
        }
        let bump_filled = |ctx: &mut SimulationContext| {
            if ctx.links.xsect_shape[uj] == XsectShape::FilledCircular {
                ctx.links.offset1[uj] += ctx.links.xsect_y_bot[uj];
                ctx.links.offset2[uj] += ctx.links.xsect_y_bot[uj];
            }
        };

        let n1 = ctx.links.node1[uj];
        let n2 = ctx.links.node2[uj];
        if n1 < 0 || n2 < 0 {
            bump_filled(ctx);
            continue;
        }

        let cr = ctx.link_subtypes.conduit_row(j);
        let ucr = cr as usize;
        let mut length = if cr >= 0 { ctx.link_subtypes.conduits.length[ucr] } else { 0.0 };
        if length <= 0.0 {
            bump_filled(ctx);
            continue;
        }
        length = conduit_true_length(ctx, j);
        if length <= 0.0 {
            bump_filled(ctx);
            continue;
        }

        if ctx.options.link_offsets == 1 && !ctx.gpkg_units_internal {
            let lname = ctx.link_names.name_of(j).to_vec();
            let inv1 = ctx.nodes.invert_elev[n1 as usize];
            let inv2 = ctx.nodes.invert_elev[n2 as usize];
            ctx.links.offset1[uj] = legacy_offset_height(ctx.links.offset1[uj], inv1, &mut ctx.warnings, &lname);
            ctx.links.offset2[uj] = legacy_offset_height(ctx.links.offset2[uj], inv2, &mut ctx.warnings, &lname);
        }

        bump_filled(ctx);

        let elev1 = ctx.links.offset1[uj] + ctx.nodes.invert_elev[n1 as usize];
        let elev2 = ctx.links.offset2[uj] + ctx.nodes.invert_elev[n2 as usize];
        let mut delta = (elev1 - elev2).abs();

        if delta < MIN_DELTA_Z {
            delta = MIN_DELTA_Z;
            let w = format_warning(WARN_MIN_ELEV_DROP, ctx.link_names.name_of(j));
            ctx.warnings.push(w);
        }

        let mut slope;
        if delta >= length {
            slope = delta / length;
            let w = format_warning(WARN_ELEV_DROP_EXCEEDS, ctx.link_names.name_of(j));
            ctx.warnings.push(w);
        } else {
            slope = delta / (length * length - delta * delta).sqrt();
        }

        let min_slope = ctx.options.min_slope / 100.0;
        let mut min_slope_kw_sf = false;
        if min_slope > 0.0 && slope < min_slope {
            slope = min_slope;
            let w = format_warning(WARN_MIN_SLOPE, ctx.link_names.name_of(j));
            ctx.warnings.push(w);
            min_slope_kw_sf = ctx.options.routing_model == RoutingModel::Steady || ctx.options.routing_model == RoutingModel::Kinwave;
        }

        if elev1 < elev2 && !min_slope_kw_sf {
            slope = -slope;
        }

        if cr >= 0 {
            ctx.link_subtypes.conduits.slope[ucr] = slope;
        }

        if (ctx.options.routing_model == RoutingModel::Dynwave || ctx.options.routing_model == RoutingModel::Fv)
            && slope < 0.0
            && ctx.links.xsect_shape[uj] != XsectShape::Dummy
        {
            let l = &mut ctx.links;
            let (a, b) = (l.node1[uj], l.node2[uj]);
            l.node1[uj] = b;
            l.node2[uj] = a;
            let (o1, o2) = (l.offset1[uj], l.offset2[uj]);
            l.offset1[uj] = o2;
            l.offset2[uj] = o1;
            l.direction[uj] = l.direction[uj].wrapping_mul(-1);
            l.q0[uj] = -l.q0[uj];
            if cr >= 0 {
                let cd = &mut ctx.link_subtypes.conduits;
                let (li, lo) = (cd.loss_inlet[ucr], cd.loss_outlet[ucr]);
                cd.loss_inlet[ucr] = lo;
                cd.loss_outlet[ucr] = li;
                cd.slope[ucr] = -slope;
            }
        }
    }

    // Conduit flow properties (legacy conduit_validate)
    for j in 0..n_links {
        recompute_conduit_flow_properties(ctx, j);
    }

    validate_virtual_junctions(ctx);

    resolve_inlet_references(ctx);
    validate_inlet_junctions(ctx);

    // Table sequence check (legacy table_validate's monotonicity test only)
    {
        let mut errs: Vec<Bytes> = Vec::new();
        for tbl in &ctx.tables.tables {
            let is_ts = tbl.type_ == TableType::Timeseries;
            for k in 1..tbl.x.len() {
                if tbl.x[k] - tbl.x[k - 1] <= 0.0 {
                    errs.push(format_error(if is_ts { ERR_TIMESERIES_SEQUENCE } else { ERR_CURVE_SEQUENCE }, &tbl.id));
                    break;
                }
            }
        }
        ctx.errors.extend(errs);
    }

    // Rain-file station conflict (legacy rainFileConflict)
    for i in 0..ctx.gages.count() {
        let ui = i as usize;
        if ctx.gages.source[ui] != RainSource::FileRain {
            continue;
        }
        for j in 0..i {
            let uj = j as usize;
            if ctx.gages.source[uj] != RainSource::FileRain {
                continue;
            }
            if ieq(&ctx.gages.station_id[ui], &ctx.gages.station_id[uj])
                && !ieq(&ctx.gages.file_path[ui].original, &ctx.gages.file_path[uj].original)
            {
                let e = format_error(ERR_RAIN_FILE_CONFLICT, ctx.gage_names.name_of(i));
                ctx.errors.push(e);
                break;
            }
        }
    }

    // Node fullDepth adjustment from connected link crowns
    let authored_full_depth: Vec<f64> = ctx.nodes.full_depth.clone();
    let extend_to_crown = |ctx: &mut SimulationContext, n: i32, crown: f64| {
        if n < 0 || n >= n_nodes {
            return;
        }
        let un = n as usize;
        if ctx.nodes.is_virtual[un] != 0 {
            return;
        }
        if ctx.nodes.type_[un] == NodeType::Storage && !(ctx.nodes.sur_depth[un] > 0.0) {
            return;
        }
        if crown <= ctx.nodes.full_depth[un] {
            return;
        }
        ctx.nodes.full_depth[un] = crown;
    };

    for j in 0..n_links {
        let uj = j as usize;
        if ctx.links.type_[uj] == LinkType::Pump {
            continue;
        }
        if ctx.links.type_[uj] == LinkType::Orifice {
            let orr = ctx.link_subtypes.orifice_row(j);
            if orr >= 0 && ctx.link_subtypes.orifices.orifice_type[orr as usize] == 0.0 {
                continue;
            }
        }

        let mut y_full = ctx.links.xsect_y_full[uj];
        if ctx.links.type_[uj] == LinkType::Outlet
            || (ctx.links.type_[uj] == LinkType::Conduit && ctx.links.xsect_shape[uj] == XsectShape::Dummy)
        {
            y_full = TINY;
        }

        let mut off1 = ctx.links.offset1[uj];
        if ctx.links.type_[uj] == LinkType::Weir {
            let wr = ctx.link_subtypes.weir_row(j);
            if wr >= 0 {
                off1 = ctx.link_subtypes.weirs.crest_height[wr as usize];
            }
        } else if ctx.links.type_[uj] == LinkType::Outlet {
            let olr = ctx.link_subtypes.outlet_row(j);
            if olr >= 0 {
                off1 = ctx.link_subtypes.outlets.crest_height[olr as usize];
            }
        }

        let node1 = ctx.links.node1[uj];
        extend_to_crown(ctx, node1, off1 + y_full);

        if ctx.links.type_[uj] == LinkType::Conduit {
            let node2 = ctx.links.node2[uj];
            let crown = ctx.links.offset2[uj] + y_full;
            extend_to_crown(ctx, node2, crown);
        }
    }

    for n in 0..n_nodes {
        let un = n as usize;
        if ctx.nodes.full_depth[un] > authored_full_depth[un] && authored_full_depth[un] > 0.0 {
            let w = format_warning(WARN_MAX_DEPTH_INCREASED, ctx.node_names.name_of(n));
            ctx.warnings.push(w);
        }
    }

    // Report flag propagation (legacy project_validate)
    {
        let opts = &ctx.options;
        if opts.rpt_subcatchments == 1 {
            ctx.subcatches.rpt_flag.iter_mut().for_each(|f| *f = 1);
        } else if opts.rpt_subcatchments == 2 {
            for name in &opts.rpt_subcatch_names {
                let idx = ctx.subcatch_names.find(name);
                if idx >= 0 {
                    ctx.subcatches.rpt_flag[idx as usize] = 1;
                }
            }
        }
        if opts.rpt_nodes == 1 {
            ctx.nodes.rpt_flag.iter_mut().for_each(|f| *f = 1);
        } else if opts.rpt_nodes == 2 {
            for name in &opts.rpt_node_names {
                let idx = ctx.node_names.find(name);
                if idx >= 0 {
                    ctx.nodes.rpt_flag[idx as usize] = 1;
                }
            }
        }
        if opts.rpt_links == 1 {
            ctx.links.rpt_flag.iter_mut().for_each(|f| *f = 1);
        } else if opts.rpt_links == 2 {
            for name in &opts.rpt_link_names {
                let idx = ctx.link_names.find(name);
                if idx >= 0 {
                    ctx.links.rpt_flag[idx as usize] = 1;
                }
            }
        }
    }

    // Ignore flags for absent object classes (legacy project_validate)
    if ctx.aquifers.count() == 0 {
        ctx.options.ignore_groundwater = true;
    }
    if ctx.snowpacks.count() == 0 {
        ctx.options.ignore_snow_melt = true;
    }

    ctx.shrink_all_to_fit();

    let nn = ctx.nodes.count();
    ctx.node_subtypes.rebuild_index(nn);
    crate::hydrology::lid_node::validate(ctx);
}

/// [INITIAL_QUALITY]: load the FILE sidecar rows, resolve element and constituent of
/// every row (erasing failures), refuse duplicates, mirror MSX rows.
fn resolve_initial_quality(ctx: &mut SimulationContext) {
    // Rows previously loaded from the file are dropped first.
    {
        let iq = &mut ctx.initial_quality;
        let mut i = iq.count() - 1;
        while i >= 0 {
            if (i as usize) < iq.from_file.len() && iq.from_file[i as usize] != 0 {
                iq.erase(i);
            }
            i -= 1;
        }
    }
    if !ctx.initial_quality.file.is_empty() {
        let dir = parent_dir(&ctx.inp_file_path);
        let file = ctx.initial_quality.file.clone();
        let path = resolve_relative(&file, &dir);
        match IfStream::open(&path) {
            None => {
                ctx.errors.push(cat(&[b"[INITIAL_QUALITY] FILE '", &file, b"' not found or unreadable (", &path, b")."]));
            }
            Some(mut input) => {
                let mut line = Bytes::new();
                let mut lineno: i32 = 0;
                while input.getline(&mut line) {
                    lineno += 1;
                    for ch in line.iter_mut() {
                        if *ch == b',' || *ch == b';' || *ch == b'\t' {
                            *ch = b' ';
                        }
                    }
                    let tok = tokenizer::tokenize(&line);
                    if tok.is_empty() {
                        continue;
                    }
                    let ln = to_string_long(lineno as i64);
                    if tok.len() < 4 {
                        ctx.errors.push(cat(&[b"[INITIAL_QUALITY] FILE '", &file, b"' line ", &ln, b": needs scope element constituent value."]));
                        continue;
                    }
                    let scope = tokenizer::to_upper(&tok[0]);
                    let (r, _) = crate::core::charconv::from_chars_double(&tok[3]);
                    let (numeric, v) = match r {
                        Ok(v) => (true, v),
                        Err(_) => (false, 0.0),
                    };
                    if !numeric {
                        if lineno == 1 {
                            continue;
                        }
                        ctx.errors.push(cat(&[b"[INITIAL_QUALITY] FILE '", &file, b"' line ", &ln, b": bad value '", &tok[3], b"'."]));
                        continue;
                    }
                    if scope.as_slice() != b"NODE" && scope.as_slice() != b"LINK" {
                        ctx.errors.push(cat(&[b"[INITIAL_QUALITY] FILE '", &file, b"' line ", &ln, b": scope must be NODE or LINK."]));
                        continue;
                    }
                    ctx.initial_quality.add_with(
                        scope.as_slice() == b"LINK",
                        &tok[1],
                        &tok[2],
                        v,
                        -1,
                        InitialQualityData::K_KIND_UNRESOLVED,
                        true,
                    );
                }
            }
        }
    }

    let mut i = ctx.initial_quality.count() - 1;
    while i >= 0 {
        let ui = i as usize;
        let link = ctx.initial_quality.is_link[ui] != 0;
        let en = ctx.initial_quality.elem_name[ui].clone();
        let cons = ctx.initial_quality.constituent[ui].clone();
        let scope_word: &[u8] = if link { b"link" } else { b"node" };

        let ei = if link { ctx.link_names.find(&en) } else { ctx.node_names.find(&en) };
        ctx.initial_quality.elem_idx[ui] = ei;
        if ei < 0 {
            ctx.errors.push(cat(&[b"[INITIAL_QUALITY] unknown ", scope_word, b" '", &en, b"'."]));
            ctx.initial_quality.erase(i);
            i -= 1;
            continue;
        }

        if cons.as_slice() == b"__WATER_AGE__" {
            ctx.initial_quality.kind[ui] = InitialQualityData::K_KIND_WATER_AGE;
            if !ctx.options.water_age {
                ctx.warnings.push(
                    "[INITIAL_QUALITY] __WATER_AGE__ rows are present but [OPTIONS] WATER_AGE is OFF \u{2014} the rows are inert this simulation."
                        .as_bytes()
                        .to_vec(),
                );
            }
        } else if cons.as_slice() == b"__TEMPERATURE__" {
            ctx.initial_quality.kind[ui] = InitialQualityData::K_KIND_TEMPERATURE;
            if !ctx.options.heat_transport {
                ctx.warnings.push(
                    "[INITIAL_QUALITY] __TEMPERATURE__ rows are present but [OPTIONS] HEAT_TRANSPORT is OFF \u{2014} the rows are inert this simulation."
                        .as_bytes()
                        .to_vec(),
                );
            }
        } else {
            ctx.initial_quality.kind[ui] = ctx.pollutant_names.find(&cons);
            let mut resolved = ctx.initial_quality.kind[ui] >= 0;
            if !resolved && ctx.reactions.configured {
                let m = ctx.reactions.find_species(&cons);
                if m >= 0 {
                    ctx.initial_quality.kind[ui] = InitialQualityData::msx_kind(m);
                    resolved = true;
                }
            }
            if !resolved {
                ctx.errors.push(cat(&[
                    b"[INITIAL_QUALITY] unknown constituent '",
                    &cons,
                    b"' at ",
                    scope_word,
                    b" '",
                    &en,
                    "' \u{2014} [POLLUTANTS] names, reactions-component species and __WATER_AGE__/__TEMPERATURE__ are accepted here."
                        .as_bytes(),
                ]));
                ctx.initial_quality.erase(i);
                i -= 1;
                continue;
            }
            if ctx.initial_quality.value[ui] < 0.0 {
                ctx.errors.push(cat(&[b"[INITIAL_QUALITY] negative value for '", &cons, b"' at ", scope_word, b" '", &en, b"'."]));
                ctx.initial_quality.erase(i);
                i -= 1;
                continue;
            }
        }
        i -= 1;
    }

    // Duplicate (scope, element, constituent) keys: the second occurrence is refused.
    let mut i: i32 = 0;
    while i < ctx.initial_quality.count() {
        let ui = i as usize;
        for j in 0..i {
            let uj = j as usize;
            let iq = &ctx.initial_quality;
            if iq.is_link[ui] == iq.is_link[uj] && iq.elem_idx[ui] == iq.elem_idx[uj] && iq.kind[ui] == iq.kind[uj] {
                let scope_word: &[u8] = if iq.is_link[ui] != 0 { b"link" } else { b"node" };
                let e = cat(&[b"[INITIAL_QUALITY] duplicate row for '", &iq.constituent[ui], b"' at ", scope_word, b" '", &iq.elem_name[ui], b"'."]);
                ctx.errors.push(e);
                ctx.initial_quality.erase(i);
                i -= 1;
                break;
            }
        }
        i += 1;
    }

    if ctx.reactions.configured {
        let errs = crate::transport::msx_initial_quality::mirror_initial_quality_msx_rows(ctx);
        for e in errs {
            ctx.errors.push(e);
        }
    }
}
