// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/input/SectionHandlers2D.hpp and
// SectionHandlers2D.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// (SectionHandlers2D.cpp has no SPDX/license comment block; its doc comment carries the
// copyright and "@license Apache-2.0", reproduced above.)
//
// The [2D_*] section parsers, their registration (register2DSections), the external mesh file
// ([2D_MESH_FILE], load2DMeshExternalFile with its own "mini" registry) and the SI-units
// prescan (prescan2DUnitsHeader).
//   * tryParseInt is glibc strtol on the token's c_str() (core::cstdlib::strtol): leading
//     white space and a sign accepted, ok when something converts and the C string is
//     consumed; the long narrows to int as static_cast does (modulo 2^32). tryParseDouble is
//     glibc strtod the same way (core::cstdlib::strtod).
//   * Handlers reach the router's members through shared handles (twod::Shared) where upstream
//     captures references; the [2D_INFILTRATION*] handlers reach the Infil2D through
//     ctx.twod_io.infil, as upstream does.
//   * load2DMeshExternalFile's mini registry needs 'static handlers too: the warnings sink
//     upstream passes as a pointer (the main context's warnings) is a shared buffer the
//     sidecar's [2D_OPTIONS] handler appends to, moved into the caller's vector after the read
//     (nothing else writes it during the read, so the order is upstream's).
//   * Not translated here (the .inp writer's and the C API's, off the open path):
//     format2DOptionValue, format2DOptionValueEx, formatReport2DVars' callers, is2DOptionKey's
//     callers. is2DOptionKey and current2DOptionKeys are translated.

#![allow(clippy::all)]

use std::cell::RefCell;
use std::rc::Rc;

use crate::core::context::SimulationContext;
use crate::core::cstdio::IfStream;
use crate::core::cstdlib::{c_isspace, c_str, strtod, strtol};
use crate::core::error_codes::{format_warning, WARN_2D_OPTION_RETIRED};
use crate::core::fs_path::FsPath;
use crate::core::options::Bytes;
use crate::data::node_coupling_2d::NodeCoupling2DKey;
use crate::input::input_reader::InputReader;
use crate::input::parse_utils::parse_time_seconds;
use crate::input::section_registry::SectionRegistry;
use crate::input::tokenizer;
use crate::twod::data::boundary_data::BoundaryType;
use crate::twod::data::mesh_data::{MeshData, TriCouplingRow, K_MAX_CELL_VERTS};
use crate::twod::data::pending_rows_2d::*;
use crate::twod::data::report_2d_vars as report2d;
use crate::twod::data::solver_options_2d::*;
use crate::twod::infil::infil_2d::*;
use crate::twod::Shared;

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut v = Bytes::new();
    for p in parts {
        v.extend_from_slice(p);
    }
    v
}

/// iequals(): equal lengths and equal bytes under C-locale toupper.
fn iequals(a: &[u8], b: &[u8]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| x.to_ascii_uppercase() == y.to_ascii_uppercase())
}

/// tryParseInt(): strtol base 10 narrowed to int; `ok` when it converted the whole C string.
fn try_parse_int(s: &[u8], ok: &mut bool) -> i32 {
    let (val, used) = strtol(s);
    *ok = used > 0 && used == c_str(s).len();
    val as i32
}

/// tryParseDouble(): strtod; `ok` when it converted the whole C string.
fn try_parse_double(s: &[u8], ok: &mut bool) -> f64 {
    let (val, used) = strtod(s);
    *ok = used > 0 && used == c_str(s).len();
    val
}

fn find_vertex_by_tag(mesh: &MeshData, tag: &[u8]) -> i32 {
    for i in 0..mesh.n_vertices() {
        if mesh.vtag[i as usize] == tag {
            return i;
        }
    }
    -1
}

fn find_triangle_by_tag(mesh: &MeshData, tag: &[u8]) -> i32 {
    for i in 0..mesh.n_triangles() {
        if mesh.tri_tag[i as usize] == tag {
            return i;
        }
    }
    -1
}

const RETIRED_SUFFIX: &[u8] = b" was retired with the CVODE/ARKODE 2D solvers: the explicit local-inertial marcher is the only 2D integrator. Remove the line; marcher settings are MOMENTUM_EQUATION, THETA, CFL_NUMBER, LTS_TIERS, H_MOVE, FROUDE_MAX, RECONSTRUCTION_ORDER, MAX_TIMESTEP, COUPLING_AREA.";

/// retiredOption(): WARNING 104 into the sink (file load), or a hard error without one.
fn retired_option(what: &[u8], warnings: Option<&mut Vec<Bytes>>) -> Bytes {
    match warnings {
        Some(w) => {
            w.push(format_warning(WARN_2D_OPTION_RETIRED, what));
            Bytes::new()
        }
        None => cat(&[what, RETIRED_SUFFIX]),
    }
}

/// parse2DOptionsLine(): one `KEY VALUE...` row of [2D_OPTIONS] into `opts`; "" on success.
pub fn parse_2d_options_line(tokens: &[Bytes], opts: &mut SolverOptions2D, mut warnings: Option<&mut Vec<Bytes>>) -> Bytes {
    if tokens.len() < 2 {
        return b"Expected PARAMETER VALUE".to_vec();
    }
    let key = &tokens[0];
    let val = &tokens[1];
    let mut ok = false;

    if iequals(key, b"MAX_TIMESTEP") {
        opts.max_timestep = try_parse_double(val, &mut ok);
        if !ok {
            return b"Invalid MAX_TIMESTEP value".to_vec();
        }
    } else if iequals(key, b"DRY_DEPTH") {
        opts.dry_depth = try_parse_double(val, &mut ok);
        if !ok {
            return b"Invalid DRY_DEPTH value".to_vec();
        }
    } else if iequals(key, b"COUPLING_CD") {
        opts.coupling_cd = try_parse_double(val, &mut ok);
        if !ok {
            return b"Invalid COUPLING_CD value".to_vec();
        }
    } else if iequals(key, b"DISPERSION") {
        opts.dispersion = try_parse_double(val, &mut ok);
        if !ok || !opts.dispersion.is_finite() || opts.dispersion < 0.0 {
            return b"Invalid DISPERSION value (m2/s, must be finite and >= 0)".to_vec();
        }
    } else if iequals(key, b"COUPLING_SYNC") {
        opts.coupling_sync = try_parse_double(val, &mut ok);
        if !ok || opts.coupling_sync < 0.0 {
            return b"Invalid COUPLING_SYNC value (seconds, >= 0)".to_vec();
        }
    } else if iequals(key, b"LIMITER_EPSILON") {
        opts.limiter_epsilon = try_parse_double(val, &mut ok);
        if !ok {
            return b"Invalid LIMITER_EPSILON value".to_vec();
        }
    } else if iequals(key, b"FLUX_DH_EPS") {
        opts.flux_dh_eps = try_parse_double(val, &mut ok);
        if !ok {
            return b"Invalid FLUX_DH_EPS value".to_vec();
        }
    } else if iequals(key, b"CELL_CLOSURE") {
        if iequals(val, b"FLAT") {
            opts.cell_closure = CellClosure2D::Flat;
        } else if iequals(val, b"VFR") {
            opts.cell_closure = CellClosure2D::Vfr;
        } else {
            return cat(&[b"Unknown CELL_CLOSURE: ", val]);
        }
    } else if iequals(key, b"FACE_RECONSTRUCTION") {
        if iequals(val, b"MEAN") {
            opts.face_reconstruction = FaceDepth2D::Mean;
        } else if iequals(val, b"VFR_FACE") {
            opts.face_reconstruction = FaceDepth2D::VfrFace;
        } else {
            return cat(&[b"Unknown FACE_RECONSTRUCTION: ", val]);
        }
    } else if iequals(key, b"VFR_MIN_WET_FRAC") {
        let frac = try_parse_double(val, &mut ok);
        if !ok || frac <= 0.0 || frac > 0.5 {
            return b"Invalid VFR_MIN_WET_FRAC value (expected (0, 0.5])".to_vec();
        }
        opts.vfr_min_wet_frac = frac;
    } else if iequals(key, b"RAINFALL_MODE") {
        if iequals(val, b"NATURAL_NEIGHBOUR") || iequals(val, b"NATURAL_NEIGHBOR") {
            opts.rainfall_mode = RainfallMode::NaturalNeighbour;
        } else if iequals(val, b"NEAREST_NEIGHBOUR") || iequals(val, b"NEAREST_NEIGHBOR") {
            opts.rainfall_mode = RainfallMode::NearestNeighbour;
        } else if iequals(val, b"SYSTEM") {
            opts.rainfall_mode = RainfallMode::System;
        } else if iequals(val, b"NONE") {
            opts.rainfall_mode = RainfallMode::None;
        } else {
            return cat(&[b"Unknown RAINFALL_MODE: ", val]);
        }
    } else if iequals(key, b"REPORT_2D") {
        if iequals(val, b"YES") || val.as_slice() == b"1" {
            opts.report_2d = true;
        } else if iequals(val, b"NO") || val.as_slice() == b"0" {
            opts.report_2d = false;
        } else {
            return b"Invalid REPORT_2D value (YES/NO)".to_vec();
        }
    } else if iequals(key, b"OUTPUT_FILE") {
        // Stored as the raw token; resolved against the .inp directory later.
        opts.output_file.set(val);
    } else if iequals(key, b"INTEGRATOR") {
        if !iequals(val, b"EXPLICIT") {
            return retired_option(&cat(&[b"INTEGRATOR ", val]), warnings);
        }
    } else if iequals(key, b"MOMENTUM_EQUATION") {
        if iequals(val, b"LOCAL_INERTIAL") || iequals(val, b"LI") {
            opts.momentum = Momentum2D::LocalInertial;
        } else if iequals(val, b"FULL_SWE") || iequals(val, b"SWE") || iequals(val, b"FULL") {
            opts.momentum = Momentum2D::FullSwe;
        } else if iequals(val, b"DIFFUSIVE_WAVE") || iequals(val, b"DW") || iequals(val, b"DIFFUSIVE") {
            opts.momentum = Momentum2D::DiffusiveWave;
        } else {
            return cat(&[b"Unknown MOMENTUM_EQUATION: ", val, b" (expected LOCAL_INERTIAL|FULL_SWE|DIFFUSIVE_WAVE)"]);
        }
    } else if iequals(key, b"CORIOLIS") {
        if iequals(val, b"OFF") {
            opts.coriolis = Coriolis2D::Off;
        } else if iequals(val, b"LATITUDE") {
            opts.coriolis = Coriolis2D::Latitude;
        } else if iequals(val, b"CONSTANT_F") {
            opts.coriolis = Coriolis2D::ConstantF;
        } else {
            return b"Invalid CORIOLIS (OFF|LATITUDE|CONSTANT_F)".to_vec();
        }
    } else if iequals(key, b"CORIOLIS_LATITUDE")
        || iequals(key, b"CORIOLIS_F")
        || iequals(key, b"SMAGORINSKY_CS")
        || iequals(key, b"EDDY_VISCOSITY_BACKGROUND")
        || iequals(key, b"MIXING_LENGTH_COEFFICIENT")
    {
        if iequals(val, b"UNSET") && (iequals(key, b"CORIOLIS_LATITUDE") || iequals(key, b"CORIOLIS_F")) {
            if iequals(key, b"CORIOLIS_LATITUDE") {
                opts.coriolis_latitude_set = false;
            } else {
                opts.coriolis_f_set = false;
            }
            return Bytes::new();
        }
        let x = try_parse_double(val, &mut ok);
        if !ok || !x.is_finite() {
            return cat(&[b"Invalid finite value for ", key]);
        }
        if iequals(key, b"CORIOLIS_LATITUDE") {
            if x < -90.0 || x > 90.0 {
                return b"CORIOLIS_LATITUDE must be in [-90,90] degrees north".to_vec();
            }
            opts.coriolis_latitude = x;
            opts.coriolis_latitude_set = true;
        } else if iequals(key, b"CORIOLIS_F") {
            opts.coriolis_f = x;
            opts.coriolis_f_set = true;
        } else if iequals(key, b"SMAGORINSKY_CS") {
            if x < 0.0 {
                return b"SMAGORINSKY_CS must be nonnegative".to_vec();
            }
            opts.smagorinsky_cs = x;
        } else if iequals(key, b"EDDY_VISCOSITY_BACKGROUND") {
            if x < 0.0 {
                return b"EDDY_VISCOSITY_BACKGROUND must be nonnegative".to_vec();
            }
            opts.eddy_viscosity_background = x;
        } else {
            if x < 0.0 {
                return b"MIXING_LENGTH_COEFFICIENT must be nonnegative".to_vec();
            }
            opts.mixing_length_coefficient = x;
        }
    } else if iequals(key, b"TURBULENCE_MODEL") {
        if iequals(val, b"NONE") {
            opts.turbulence = Turbulence2D::None;
        } else if iequals(val, b"SMAGORINSKY") {
            opts.turbulence = Turbulence2D::Smagorinsky;
        } else if iequals(val, b"MIXING_LENGTH") {
            opts.turbulence = Turbulence2D::MixingLength;
        } else {
            return b"Invalid TURBULENCE_MODEL (NONE|SMAGORINSKY|MIXING_LENGTH)".to_vec();
        }
    } else if iequals(key, b"FRONT_REBUILD") {
        if iequals(val, b"AUTO") {
            opts.front_rebuild = -1;
        } else if iequals(val, b"YES") || iequals(val, b"ON") {
            opts.front_rebuild = 1;
        } else if iequals(val, b"NO") || iequals(val, b"OFF") {
            opts.front_rebuild = 0;
        } else {
            return cat(&[b"Unknown FRONT_REBUILD: ", val, b" (expected AUTO|YES|NO)"]);
        }
    } else if iequals(key, b"RECONSTRUCTION_ORDER") {
        let k = try_parse_int(val, &mut ok);
        if !ok || (k != 1 && k != 2) {
            return b"Invalid RECONSTRUCTION_ORDER value (expected 1|2)".to_vec();
        }
        opts.reconstruction_order = k;
    } else if iequals(key, b"THETA") {
        let th = try_parse_double(val, &mut ok);
        if !ok || !th.is_finite() || th <= 0.0 || th > 1.0 {
            return b"Invalid THETA value (expected (0, 1])".to_vec();
        }
        opts.theta = th;
    } else if iequals(key, b"THETA_REFERENCE_STEP") {
        let step = try_parse_double(val, &mut ok);
        if !ok || !step.is_finite() || step < 0.0 {
            return b"Invalid THETA_REFERENCE_STEP value (expected finite seconds >= 0)".to_vec();
        }
        opts.theta_reference_step = step;
    } else if iequals(key, b"CFL_NUMBER") {
        let c = try_parse_double(val, &mut ok);
        if !ok || c <= 0.0 || c > 1.0 {
            return b"Invalid CFL_NUMBER value (expected (0, 1])".to_vec();
        }
        opts.cfl_number = c;
    } else if iequals(key, b"H_MOVE") {
        let h = try_parse_double(val, &mut ok);
        if !ok || h < 0.0 {
            return b"Invalid H_MOVE value (expected metres >= 0)".to_vec();
        }
        opts.h_move = h;
    } else if iequals(key, b"LTS_TIERS") {
        let k = try_parse_int(val, &mut ok);
        if !ok || k < 1 || k > 8 {
            return b"Invalid LTS_TIERS value (expected 1..8)".to_vec();
        }
        opts.lts_tiers = k;
    } else if iequals(key, b"FROUDE_MAX") {
        let f = try_parse_double(val, &mut ok);
        if !ok || f <= 0.0 {
            return b"Invalid FROUDE_MAX value (expected > 0)".to_vec();
        }
        opts.froude_max = f;
    } else if iequals(key, b"ADVECTION") {
        // DEPRECATED: a file load names the replacement.
        if iequals(val, b"YES") || iequals(val, b"ON") || iequals(val, b"TRUE") {
            if let Some(w) = warnings.as_deref_mut() {
                w.push(
                    "[2D_OPTIONS] ADVECTION YES is deprecated \u{2014} use MOMENTUM_EQUATION FULL_SWE for a conservative convective term with shock capturing."
                        .as_bytes()
                        .to_vec(),
                );
            }
        }
        if iequals(val, b"YES") || iequals(val, b"ON") || iequals(val, b"TRUE") {
            opts.advection = true;
        } else if iequals(val, b"NO") || iequals(val, b"OFF") || iequals(val, b"FALSE") {
            opts.advection = false;
        } else {
            return cat(&[b"Unknown ADVECTION: ", val, b" (expected YES|NO)"]);
        }
    } else if iequals(key, b"COUPLING_IN_FLOODING") {
        if iequals(val, b"YES") || iequals(val, b"TRUE") {
            opts.coupling_in_flooding = true;
        } else if iequals(val, b"NO") || iequals(val, b"FALSE") {
            opts.coupling_in_flooding = false;
        } else {
            return cat(&[b"Unknown COUPLING_IN_FLOODING: ", val, b" (expected YES|NO)"]);
        }
    } else if iequals(key, b"COUPLING_LAW") {
        if iequals(val, b"ORIFICE") {
            opts.coupling_law = CouplingLaw2D::Orifice;
        } else if iequals(val, b"WEIR_ORIFICE") {
            opts.coupling_law = CouplingLaw2D::WeirOrifice;
        } else {
            return cat(&[b"Unknown COUPLING_LAW: ", val, b" (expected ORIFICE|WEIR_ORIFICE)"]);
        }
    } else if iequals(key, b"COUPLING_HEAD") {
        if iequals(val, b"CELL") {
            opts.coupling_head = CouplingHead2D::Cell;
        } else if iequals(val, b"VERTEX") {
            opts.coupling_head = CouplingHead2D::Vertex;
        } else if iequals(val, b"STENCIL_WET_MEAN") {
            opts.coupling_head = CouplingHead2D::StencilWetMean;
        } else {
            return cat(&[b"Unknown COUPLING_HEAD: ", val, b" (expected CELL|VERTEX|STENCIL_WET_MEAN)"]);
        }
    } else if iequals(key, b"COUPLING_POND_MODE") {
        if iequals(val, b"FOOTPRINT") {
            opts.coupling_pond_mode = CouplingPondMode2D::Footprint;
        } else if iequals(val, b"CHAMBER") {
            opts.coupling_pond_mode = CouplingPondMode2D::Chamber;
        } else if iequals(val, b"TRACKING") {
            opts.coupling_pond_mode = CouplingPondMode2D::Tracking;
        } else {
            return cat(&[b"Unknown COUPLING_POND_MODE: ", val, b" (expected FOOTPRINT|CHAMBER|TRACKING)"]);
        }
    } else if iequals(key, b"COUPLING_POND_RELAX") {
        let v = try_parse_double(val, &mut ok);
        if !ok || v < 0.0 {
            return b"Invalid COUPLING_POND_RELAX value (expected seconds >= 0)".to_vec();
        }
        opts.coupling_pond_relax = v;
    } else if iequals(key, b"COUPLING_RIM_BAND") {
        if iequals(val, b"DEFAULT") {
            opts.coupling_rim_band = -1.0;
            return Bytes::new();
        }
        let v = try_parse_double(val, &mut ok);
        if !ok || v < 0.0 {
            return b"Invalid COUPLING_RIM_BAND value (expected metres >= 0)".to_vec();
        }
        opts.coupling_rim_band = v;
    } else if iequals(key, b"COUPLING_EQUALIZE") {
        if iequals(val, b"YES") {
            opts.coupling_equalize = 1;
        } else if iequals(val, b"NO") {
            opts.coupling_equalize = 0;
        } else if iequals(val, b"DEFAULT") {
            opts.coupling_equalize = -1;
        } else {
            return cat(&[b"Unknown COUPLING_EQUALIZE: ", val, b" (expected YES|NO|DEFAULT)"]);
        }
    } else if iequals(key, b"COUPLING_COEFFS") {
        let mut preset = 0.0;
        let err = crate::core::node_coupling_2d::parse_value(NodeCoupling2DKey::Coeffs, val, &mut preset);
        if !err.is_empty() {
            return cat(&[b"Unknown ", &err]);
        }
        opts.coupling_coeffs = preset as i32;
    } else if iequals(key, b"COUPLING_AREA") {
        if iequals(val, b"AUTO") {
            opts.coupling_area_auto = true;
        } else if iequals(val, b"DEFAULT") {
            opts.coupling_area_auto = false;
        } else {
            return cat(&[b"Unknown COUPLING_AREA: ", val, b" (expected AUTO|DEFAULT)"]);
        }
    } else if iequals(key, b"BACKEND") {
        if iequals(val, b"AUTO") {
            opts.backend = Backend2D::Auto;
        } else if iequals(val, b"CPU") {
            opts.backend = Backend2D::Cpu;
        } else if iequals(val, b"OMP") {
            opts.backend = Backend2D::Omp;
        } else if iequals(val, b"CUDA") {
            opts.backend = Backend2D::Cuda;
        } else if iequals(val, b"HIP") {
            opts.backend = Backend2D::Hip;
        } else if iequals(val, b"SYCL") {
            opts.backend = Backend2D::Sycl;
        } else {
            return cat(&[b"Unknown BACKEND: ", val, b" (expected AUTO|CPU|OMP|CUDA|HIP|SYCL)"]);
        }
    } else if iequals(key, b"OUTPUT_PRECISION") {
        if iequals(val, b"FLOAT32") || iequals(val, b"F32") || iequals(val, b"SINGLE") {
            opts.output_precision = OutputPrecision2D::Float32;
        } else if iequals(val, b"FLOAT64") || iequals(val, b"F64") || iequals(val, b"DOUBLE") {
            opts.output_precision = OutputPrecision2D::Float64;
        } else {
            return cat(&[b"Unknown OUTPUT_PRECISION: ", val, b" (expected FLOAT32|FLOAT64)"]);
        }
    } else if iequals(key, b"OUTPUT_COMPRESSION") {
        let lvl = try_parse_int(val, &mut ok);
        if !ok || lvl < 0 || lvl > 9 {
            return b"Invalid OUTPUT_COMPRESSION value (expected 0..9)".to_vec();
        }
        opts.output_compression = lvl;
    } else if iequals(key, b"REPORT_2D_VARIABLES") {
        let mut joined: Bytes = Bytes::new();
        for i in 1..tokens.len() {
            if i > 1 {
                joined.push(b' ');
            }
            joined.extend_from_slice(&tokens[i]);
        }
        let mut mask: u32 = 0;
        let err = parse_report_2d_vars(&joined, &mut mask);
        if !err.is_empty() {
            return err;
        }
        opts.report_2d_vars = mask;
    } else if iequals(key, b"REPORT_2D_SPECIES") {
        let mut joined: Bytes = Bytes::new();
        for i in 1..tokens.len() {
            if i > 1 {
                joined.push(b' ');
            }
            joined.extend_from_slice(&tokens[i]);
        }
        opts.report_2d_species = report2d::parse_species(&joined); // ALL -> empty
    } else if iequals(key, b"REPORT_2D_STEP") {
        let sec = parse_time_seconds(val);
        if !(sec >= 0.0) || !sec.is_finite() {
            return b"Invalid REPORT_2D_STEP value (expected HH:MM:SS or seconds >= 0)".to_vec();
        }
        opts.report_2d_step = sec;
    } else if iequals(key, b"INFILTRATION") {
        if iequals(val, b"YES") || iequals(val, b"ON") {
            opts.infiltration = 1;
        } else if iequals(val, b"NO") || iequals(val, b"OFF") {
            opts.infiltration = 0;
        } else if iequals(val, b"AUTO") {
            opts.infiltration = -1;
        } else {
            return cat(&[b"Unknown INFILTRATION: ", val, b" (expected YES|NO|AUTO)"]);
        }
    } else if iequals(key, b"INFIL_STEP") {
        let secs = parse_time_seconds(val);
        if !(secs >= 0.0) || !secs.is_finite() {
            return b"Invalid INFIL_STEP value (expected hh:mm:ss or seconds >= 0)".to_vec();
        }
        opts.infil_step = secs;
    } else if iequals(key, b"INFIL_DEFAULT_METHOD") {
        const K_METHODS: [&[u8]; 7] = [b"NONE", b"HORTON", b"MOD_HORTON", b"GREEN_AMPT", b"MOD_GREEN_AMPT", b"CURVE_NUMBER", b"CONSTANT"];
        let mut known = false;
        for m in K_METHODS.iter() {
            if iequals(val, m) {
                opts.infil_default_method = m.to_vec();
                known = true;
                break;
            }
        }
        // Legacy [INFILTRATION] spellings of the same methods.
        if !known && iequals(val, b"MODIFIED_HORTON") {
            opts.infil_default_method = b"MOD_HORTON".to_vec();
            known = true;
        }
        if !known && iequals(val, b"MODIFIED_GREEN_AMPT") {
            opts.infil_default_method = b"MOD_GREEN_AMPT".to_vec();
            known = true;
        }
        if !known && iequals(val, b"CURVE_NUM") {
            opts.infil_default_method = b"CURVE_NUMBER".to_vec();
            known = true;
        }
        if !known {
            return cat(&[b"Unknown INFIL_DEFAULT_METHOD: ", val, b" (expected NONE|HORTON|MOD_HORTON|GREEN_AMPT|MOD_GREEN_AMPT|CURVE_NUMBER|CONSTANT)"]);
        }
    } else if iequals(key, b"INFIL_DESTINATION") {
        if iequals(val, b"LOST") {
            opts.infil_destination = b"LOST".to_vec();
        } else if iequals(val, b"SUBCATCH_AQUIFER") {
            opts.infil_destination = b"SUBCATCH_AQUIFER".to_vec();
        } else if iequals(val, b"AQUIFER_2D") {
            opts.infil_destination = b"AQUIFER_2D".to_vec();
        } else {
            return cat(&[b"Unknown INFIL_DESTINATION: ", val, b" (expected LOST|SUBCATCH_AQUIFER|AQUIFER_2D)"]);
        }
    } else if iequals(key, b"HEAT_MODEL") {
        if iequals(val, b"COUPLED") {
            opts.coupled_heat = true;
        } else if iequals(val, b"TRANSPORT_ONLY") {
            opts.coupled_heat = false;
        } else {
            return b"HEAT_MODEL expects TRANSPORT_ONLY or COUPLED".to_vec();
        }
        opts.heat_model_authored = true;
    } else if iequals(key, b"HEAT_SURFACE_EXCHANGE") || iequals(key, b"HEAT_RADIATIVE_EXCHANGE") {
        let mut mode = -1;
        if iequals(val, b"YES") {
            mode = 1;
        } else if iequals(val, b"NO") {
            mode = 0;
        } else if !iequals(val, b"AUTO") {
            return cat(&[key, b" expects AUTO, YES or NO"]);
        }
        if iequals(key, b"HEAT_SURFACE_EXCHANGE") {
            opts.heat_surface = mode;
            opts.heat_surface_authored = true;
        } else {
            opts.heat_radiative = mode;
            opts.heat_radiative_authored = true;
        }
    } else if iequals(key, b"BED_HEAT_MODEL") {
        if iequals(val, b"NONE") {
            opts.bed_heat_model = 0;
        } else if iequals(val, b"LOCAL_STORAGE") {
            opts.bed_heat_model = 1;
        } else if iequals(val, b"GROUND_COLUMN") {
            opts.bed_heat_model = 2;
        } else {
            return b"BED_HEAT_MODEL expects NONE, LOCAL_STORAGE or GROUND_COLUMN".to_vec();
        }
        opts.bed_heat_authored = true;
    } else if iequals(key, b"EVAPORATION") {
        if iequals(val, b"YES") || iequals(val, b"ON") || iequals(val, b"FORCING") {
            opts.evaporation = 1;
        } else if iequals(val, b"NO") || iequals(val, b"OFF") {
            opts.evaporation = 0;
        } else if iequals(val, b"CLIMATE") {
            opts.evaporation = 2;
        } else {
            return cat(&[b"Unknown EVAPORATION: ", val, b" (expected NO|YES|CLIMATE)"]);
        }
    } else if iequals(key, b"TRANSPORT_POLLUTANTS") || iequals(key, b"TRANSPORT_MSX") || iequals(key, b"TRANSPORT_AGE") || iequals(key, b"TRANSPORT_TEMPERATURE") {
        let on;
        if iequals(val, b"YES") || iequals(val, b"ON") {
            on = true;
        } else if iequals(val, b"NO") || iequals(val, b"OFF") {
            on = false;
        } else {
            return cat(&[b"Unknown ", key, b": ", val, b" (expected YES|NO)"]);
        }
        if iequals(key, b"TRANSPORT_POLLUTANTS") {
            opts.transport_pollutants = on;
        } else if iequals(key, b"TRANSPORT_MSX") {
            opts.transport_msx = on;
        } else if iequals(key, b"TRANSPORT_AGE") {
            opts.transport_age = on;
        } else {
            opts.transport_temperature = on;
        }
    } else if iequals(key, b"GROUNDWATER") {
        if iequals(val, b"ON") || iequals(val, b"YES") {
            opts.groundwater = 1;
        } else if iequals(val, b"OFF") || iequals(val, b"NO") {
            opts.groundwater = 0;
        } else if iequals(val, b"AUTO") {
            opts.groundwater = -1;
        } else {
            return cat(&[b"Unknown GROUNDWATER: ", val, b" (expected YES|NO|AUTO)"]);
        }
    } else if iequals(key, b"GW_ET") {
        const K_ET: [&[u8]; 5] = [b"AUTO", b"NONE", b"CAPILLARY_RISE", b"BOUNDARY_ET", b"BOTH"];
        let mut known = false;
        for e in K_ET.iter() {
            if iequals(val, e) {
                opts.gw_et = e.to_vec();
                known = true;
                break;
            }
        }
        if !known {
            return cat(&[b"Unknown GW_ET: ", val, b" (expected NONE|CAPILLARY_RISE|BOUNDARY_ET|BOTH)"]);
        }
    } else if is_2d_retired_option_key(key) {
        return retired_option(key, warnings);
    } else {
        return cat(&[b"Unknown 2D_OPTIONS parameter: ", key]);
    }
    Bytes::new()
}

/// is2DRetiredOptionKey()
pub fn is_2d_retired_option_key(key: &[u8]) -> bool {
    const K_RETIRED: [&[u8]; 14] = [
        b"MIN_TIMESTEP",
        b"REL_TOLERANCE",
        b"ABS_TOLERANCE",
        b"MAX_CVODE_STEPS",
        b"MAX_KRYLOV_DIM",
        b"LINEAR_SOLVER",
        b"PRECONDITIONER",
        b"JACOBIAN",
        b"ATOL_AREA_REF",
        b"COUPLING_INTERVAL",
        b"COUPLING_WINDOW",
        b"ACTIVE_SET",
        b"ACTIVE_SET_HALO",
        b"MOMENTUM",
    ];
    K_RETIRED.iter().any(|k| iequals(key, k))
}

/// current2DOptionKeys(): every current key in documentation order.
pub fn current_2d_option_keys() -> &'static [&'static [u8]] {
    const K_KEYS: [&[u8]; 61] = [
        b"MAX_TIMESTEP",
        b"DRY_DEPTH",
        b"COUPLING_CD",
        b"COUPLING_SYNC",
        b"DISPERSION",
        b"LIMITER_EPSILON",
        b"FLUX_DH_EPS",
        b"RAINFALL_MODE",
        b"REPORT_2D",
        b"CELL_CLOSURE",
        b"FACE_RECONSTRUCTION",
        b"VFR_MIN_WET_FRAC",
        b"OUTPUT_FILE",
        b"INTEGRATOR",
        b"MOMENTUM_EQUATION",
        b"CORIOLIS",
        b"CORIOLIS_LATITUDE",
        b"CORIOLIS_F",
        b"TURBULENCE_MODEL",
        b"SMAGORINSKY_CS",
        b"EDDY_VISCOSITY_BACKGROUND",
        b"MIXING_LENGTH_COEFFICIENT",
        b"RECONSTRUCTION_ORDER",
        b"FRONT_REBUILD",
        b"THETA",
        b"THETA_REFERENCE_STEP",
        b"CFL_NUMBER",
        b"H_MOVE",
        b"LTS_TIERS",
        b"FROUDE_MAX",
        b"ADVECTION",
        b"COUPLING_AREA",
        b"COUPLING_IN_FLOODING",
        b"COUPLING_LAW",
        b"COUPLING_HEAD",
        b"COUPLING_COEFFS",
        b"COUPLING_POND_MODE",
        b"COUPLING_POND_RELAX",
        b"COUPLING_RIM_BAND",
        b"COUPLING_EQUALIZE",
        b"BACKEND",
        b"OUTPUT_PRECISION",
        b"OUTPUT_COMPRESSION",
        b"REPORT_2D_VARIABLES",
        b"REPORT_2D_SPECIES",
        b"REPORT_2D_STEP",
        b"INFILTRATION",
        b"INFIL_STEP",
        b"INFIL_DEFAULT_METHOD",
        b"INFIL_DESTINATION",
        b"EVAPORATION",
        b"TRANSPORT_POLLUTANTS",
        b"TRANSPORT_MSX",
        b"TRANSPORT_AGE",
        b"TRANSPORT_TEMPERATURE",
        b"GROUNDWATER",
        b"GW_ET",
        b"HEAT_MODEL",
        b"HEAT_SURFACE_EXCHANGE",
        b"HEAT_RADIATIVE_EXCHANGE",
        b"BED_HEAT_MODEL",
    ];
    &K_KEYS
}

/// is2DOptionKey()
pub fn is_2d_option_key(key: &[u8]) -> bool {
    current_2d_option_keys().iter().any(|k| iequals(key, k))
}

/// report2DVarTokens()
pub fn report_2d_var_tokens() -> &'static [&'static [u8]] {
    &report2d::TOKENS
}

/// parseReport2DVars()
pub fn parse_report_2d_vars(text: &[u8], mask: &mut u32) -> Bytes {
    report2d::parse_mask(text, mask)
}

/// formatReport2DVars()
pub fn format_report_2d_vars(mask: u32) -> Bytes {
    report2d::format_mask(mask)
}

/// parse2DVertexLine(): `X Y Z [TAG]` appended as the next vertex.
pub fn parse_2d_vertex_line(tokens: &[Bytes], mesh: &mut MeshData) -> Bytes {
    if tokens.len() < 3 {
        return b"Expected X Y Z [TAG]".to_vec();
    }
    let mut ok = false;
    let x = try_parse_double(&tokens[0], &mut ok);
    if !ok {
        return b"Invalid X coordinate".to_vec();
    }
    let y = try_parse_double(&tokens[1], &mut ok);
    if !ok {
        return b"Invalid Y coordinate".to_vec();
    }
    let z = try_parse_double(&tokens[2], &mut ok);
    if !ok {
        return b"Invalid Z coordinate".to_vec();
    }
    let mut tag: Bytes = Bytes::new();
    if tokens.len() >= 4 {
        tag = tokens[3].clone();
    }
    let idx = mesh.n_vertices();
    mesh.resize_vertices(idx + 1);
    let u = idx as usize;
    mesh.vx[u] = x;
    mesh.vy[u] = y;
    mesh.vz[u] = z;
    mesh.vtag[u] = tag;
    Bytes::new()
}

/// parse2DTriangleLine(): `V1 V2 V3 MANNINGS_N [INIT_DEPTH] [TAG]`.
pub fn parse_2d_triangle_line(tokens: &[Bytes], mesh: &mut MeshData) -> Bytes {
    if tokens.len() < 4 {
        return b"Expected V1 V2 V3 MANNINGS_N [INIT_DEPTH] [TAG]".to_vec();
    }
    let mut ok = false;
    let v0 = try_parse_int(&tokens[0], &mut ok);
    if !ok {
        return b"Invalid V1 index".to_vec();
    }
    let v1 = try_parse_int(&tokens[1], &mut ok);
    if !ok {
        return b"Invalid V2 index".to_vec();
    }
    let v2 = try_parse_int(&tokens[2], &mut ok);
    if !ok {
        return b"Invalid V3 index".to_vec();
    }
    let n = try_parse_double(&tokens[3], &mut ok);
    if !ok {
        return b"Invalid MANNINGS_N value".to_vec();
    }
    let mut init_depth = 0.0;
    let mut tag: Bytes = Bytes::new();
    if tokens.len() >= 5 {
        let mut num = false;
        let d = try_parse_double(&tokens[4], &mut num);
        if num {
            if d < 0.0 {
                return b"Invalid INIT_DEPTH (must be >= 0)".to_vec();
            }
            init_depth = d;
            if tokens.len() >= 6 {
                tag = tokens[5].clone();
            }
        } else {
            tag = tokens[4].clone();
        }
    }
    let idx = mesh.n_triangles();
    // Cell index contract: triangles first, then quads.
    if idx > 0 && mesh.cell_vertex_count(idx - 1) == 4 {
        return b"[2D_TRIANGLES] rows must precede all [2D_QUADS] rows (cells are numbered triangles first, then quads)".to_vec();
    }
    mesh.resize_triangles(idx + 1);
    mesh.set_triangle(idx, v0, v1, v2);
    let u = idx as usize;
    mesh.mannings_n[u] = n;
    mesh.tri_init_depth[u] = init_depth;
    mesh.tri_tag[u] = tag;
    Bytes::new()
}

/// parse2DQuadLine(): `V1 V2 V3 V4 MANNINGS_N [INIT_DEPTH] [TAG]`, appended after the
/// triangles.
pub fn parse_2d_quad_line(tokens: &[Bytes], mesh: &mut MeshData) -> Bytes {
    if tokens.len() < 5 {
        return b"Expected V1 V2 V3 V4 MANNINGS_N [INIT_DEPTH] [TAG]".to_vec();
    }
    let mut ok = false;
    let mut v = [0i32; 4];
    for k in 0..4 {
        v[k] = try_parse_int(&tokens[k], &mut ok);
        if !ok {
            return cat(&[b"Invalid V", (k + 1).to_string().as_bytes(), b" index"]);
        }
    }
    let n = try_parse_double(&tokens[4], &mut ok);
    if !ok {
        return b"Invalid MANNINGS_N value".to_vec();
    }
    let mut init_depth = 0.0;
    let mut tag: Bytes = Bytes::new();
    if tokens.len() >= 6 {
        let mut num = false;
        let d = try_parse_double(&tokens[5], &mut num);
        if num {
            if d < 0.0 {
                return b"Invalid INIT_DEPTH (must be >= 0)".to_vec();
            }
            init_depth = d;
            if tokens.len() >= 7 {
                tag = tokens[6].clone();
            }
        } else {
            tag = tokens[5].clone();
        }
    }
    let idx = mesh.n_triangles();
    mesh.resize_triangles(idx + 1);
    mesh.set_quad(idx, v[0], v[1], v[2], v[3]);
    let u = idx as usize;
    mesh.mannings_n[u] = n;
    mesh.tri_init_depth[u] = init_depth;
    mesh.tri_tag[u] = tag;
    Bytes::new()
}

/// parse2DInitialVelocityLine(): `TRI U V` on an existing cell.
pub fn parse_2d_initial_velocity_line(tokens: &[Bytes], mesh: &mut MeshData) -> Bytes {
    if tokens.len() < 3 {
        return b"Expected TRI U V".to_vec();
    }
    let mut ok = false;
    let tri = try_parse_int(&tokens[0], &mut ok);
    if !ok || tri < 0 || tri >= mesh.n_triangles() {
        return cat(&[b"Invalid triangle index: ", &tokens[0]]);
    }
    let u = try_parse_double(&tokens[1], &mut ok);
    if !ok {
        return b"Invalid U value".to_vec();
    }
    let v = try_parse_double(&tokens[2], &mut ok);
    if !ok {
        return b"Invalid V value".to_vec();
    }
    mesh.tri_init_u[tri as usize] = u;
    mesh.tri_init_v[tri as usize] = v;
    Bytes::new()
}

/// parse2DVertexNodeMapLine(): `VERTEX_INDEX_OR_TAG SWMM_NODE_NAME [CD] [AREA]`.
pub fn parse_2d_vertex_node_map_line(tokens: &[Bytes], mesh: &mut MeshData) -> Bytes {
    if tokens.len() < 2 {
        return b"Expected VERTEX_INDEX_OR_TAG SWMM_NODE_NAME [CD] [AREA]".to_vec();
    }
    let mut ok = false;
    let mut vidx = try_parse_int(&tokens[0], &mut ok);
    if ok && (vidx < 0 || vidx >= mesh.n_vertices()) {
        ok = false; // numeric but not a valid index: a tag
    }
    if !ok {
        vidx = find_vertex_by_tag(mesh, &tokens[0]);
        if vidx < 0 {
            return cat(&[b"Unknown vertex index or tag: ", &tokens[0]]);
        }
    }
    if vidx < 0 || vidx >= mesh.n_vertices() {
        return cat(&[b"Vertex index out of range: ", &tokens[0]]);
    }
    let u = vidx as usize;
    mesh.vert_coupled_node_name[u] = tokens[1].clone();
    if tokens.len() >= 3 {
        let cd = try_parse_double(&tokens[2], &mut ok);
        if ok {
            mesh.vert_coupling_cd[u] = cd;
            mesh.vert_coupling_cd_set[u] = 1;
        }
    }
    if tokens.len() >= 4 {
        let area = try_parse_double(&tokens[3], &mut ok);
        if ok {
            mesh.vert_coupling_area[u] = area;
            mesh.vert_coupling_area_set[u] = 1;
        }
    }
    Bytes::new()
}

/// parse2DTriangleNodeMapLine(): `TRIANGLE_INDEX_OR_TAG SWMM_NODE_NAME [CD] [AREA]`; every row
/// appends a coupling and refreshes the per-cell mirror (last row wins there).
pub fn parse_2d_triangle_node_map_line(tokens: &[Bytes], mesh: &mut MeshData) -> Bytes {
    if tokens.len() < 2 {
        return b"Expected TRIANGLE_INDEX_OR_TAG SWMM_NODE_NAME [CD] [AREA]".to_vec();
    }
    let mut ok = false;
    let mut tidx = try_parse_int(&tokens[0], &mut ok);
    if ok && (tidx < 0 || tidx >= mesh.n_triangles()) {
        ok = false; // numeric but not a valid index: a tag
    }
    if !ok {
        tidx = find_triangle_by_tag(mesh, &tokens[0]);
        if tidx < 0 {
            return cat(&[b"Unknown triangle index or tag: ", &tokens[0]]);
        }
    }
    if tidx < 0 || tidx >= mesh.n_triangles() {
        return cat(&[b"Triangle index out of range: ", &tokens[0]]);
    }
    let mut row = TriCouplingRow::default();
    row.tri = tidx;
    row.node_name = tokens[1].clone();
    if tokens.len() >= 3 {
        let cd = try_parse_double(&tokens[2], &mut ok);
        if ok {
            row.cd = cd;
            row.cd_set = true;
        }
    }
    if tokens.len() >= 4 {
        let area = try_parse_double(&tokens[3], &mut ok);
        if ok {
            row.area = area;
            row.area_set = true;
        }
    }
    let u = tidx as usize;
    mesh.tri_coupled_node_name[u] = row.node_name.clone();
    mesh.tri_coupling_cd[u] = row.cd;
    mesh.tri_coupling_area[u] = row.area;
    mesh.tri_couplings.push(row);
    Bytes::new()
}

/// makeSectionHandler(): tokenizes each non-blank line and calls the row parser, stopping at
/// the first diagnostic (error code 5, "[2D] <err> — line: <raw>").
fn make_section_handler<F: Fn(&[Bytes]) -> Bytes + 'static>(lp: F) -> impl Fn(&mut SimulationContext, &[Bytes]) + 'static {
    move |ctx: &mut SimulationContext, lines: &[Bytes]| {
        for raw in lines {
            let tokens = tokenizer::tokenize(raw);
            if tokens.is_empty() {
                continue;
            }
            let err = lp(&tokens);
            if !err.is_empty() {
                ctx.error_code = 5; // SWMM_ERR_PARSE
                ctx.error_message = cat(&[b"[2D] ", &err, " \u{2014} line: ".as_bytes(), raw]);
                return;
            }
        }
    }
}

/// parse2DBoundaryConditionsLine(): `TRI EDGE TYPE [PARAM_1 [PARAM_2 [GROUP]]]`.
pub fn parse_2d_boundary_conditions_line(tokens: &[Bytes], pending_rows: &mut Vec<PendingBoundaryRow>) -> Bytes {
    if tokens.is_empty() {
        return Bytes::new();
    }
    if tokens.len() < 3 {
        return b"[2D_BOUNDARY_CONDITIONS] needs TRI EDGE TYPE [PARAM_1 [PARAM_2 [GROUP]]]".to_vec();
    }
    let mut ok = false;
    let mut row = PendingBoundaryRow::default();
    row.tri = try_parse_int(&tokens[0], &mut ok);
    if !ok || row.tri < 0 {
        return b"[2D_BOUNDARY_CONDITIONS] invalid TRI index".to_vec();
    }
    row.edge = try_parse_int(&tokens[1], &mut ok);
    if !ok || row.edge < 0 || row.edge >= K_MAX_CELL_VERTS {
        return b"[2D_BOUNDARY_CONDITIONS] invalid EDGE (must be 0..2 for a triangle, 0..3 for a quad)".to_vec();
    }
    let type_tok = &tokens[2];
    if iequals(type_tok, b"WALL") || iequals(type_tok, b"WALL_DRAG") {
        row.bc_type = BoundaryType::Wall as i32;
    } else if iequals(type_tok, b"NORMAL_FLOW") {
        row.bc_type = BoundaryType::NormalFlow as i32;
    } else if iequals(type_tok, b"SPECIFIED_STAGE") {
        row.bc_type = BoundaryType::SpecifiedStage as i32;
    } else if iequals(type_tok, b"TS_STAGE") {
        row.bc_type = BoundaryType::SpecifiedStage as i32;
    } else if iequals(type_tok, b"SPECIFIED_FLOW") {
        row.bc_type = BoundaryType::SpecifiedFlow as i32;
    } else if iequals(type_tok, b"TS_FLOW") {
        row.bc_type = BoundaryType::SpecifiedFlow as i32;
    } else if iequals(type_tok, b"RATING_CURVE") {
        row.bc_type = BoundaryType::RatingCurve as i32;
    } else {
        return cat(&[b"[2D_BOUNDARY_CONDITIONS] unknown TYPE: ", type_tok]);
    }
    let mut p1: Bytes = if tokens.len() > 3 { tokens[3].clone() } else { b"*".to_vec() };
    if p1 == b"*" {
        p1.clear();
    }
    // WALL keeps its legacy ignored parameter. WALL_DRAG explicitly opts in.
    if iequals(type_tok, b"WALL_DRAG") {
        row.param1 = try_parse_double(&p1, &mut ok);
        if !ok || !row.param1.is_finite() || row.param1 < 0.0 {
            return b"WALL_DRAG requires a finite nonnegative dimensionless Cw".to_vec();
        }
    }
    if !p1.is_empty() && !iequals(type_tok, b"WALL") && !iequals(type_tok, b"WALL_DRAG") {
        if iequals(type_tok, b"NORMAL_FLOW") || iequals(type_tok, b"SPECIFIED_STAGE") || iequals(type_tok, b"SPECIFIED_FLOW") {
            row.param1 = try_parse_double(&p1, &mut ok);
            if !ok {
                return cat(&[b"[2D_BOUNDARY_CONDITIONS] non-numeric PARAM_1: ", &p1]);
            }
        } else {
            // TS or curve name.
            row.name = p1.clone();
        }
    }
    // PARAM_2 reserved; GROUP optional.
    if tokens.len() > 5 && tokens[5] != b"*" {
        row.group = tokens[5].clone();
    }
    pending_rows.push(row);
    Bytes::new()
}

/// parse2DEdgeConveyanceLine(): `FROM_VERTEX TO_VERTEX CONVEYANCE` with CONVEYANCE in [0, 1].
pub fn parse_2d_edge_conveyance_line(tokens: &[Bytes], pending_rows: &mut Vec<PendingEdgeConveyanceRow>) -> Bytes {
    if tokens.is_empty() {
        return Bytes::new();
    }
    if tokens.len() < 3 {
        return b"[2D_EDGE_CONVEYANCE] needs FROM_VERTEX TO_VERTEX CONVEYANCE".to_vec();
    }
    let mut row = PendingEdgeConveyanceRow::default();
    let mut ok = false;
    row.v_from = try_parse_int(&tokens[0], &mut ok);
    if !ok || row.v_from < 0 {
        return cat(&[b"[2D_EDGE_CONVEYANCE] invalid FROM_VERTEX: ", &tokens[0]]);
    }
    row.v_to = try_parse_int(&tokens[1], &mut ok);
    if !ok || row.v_to < 0 {
        return cat(&[b"[2D_EDGE_CONVEYANCE] invalid TO_VERTEX: ", &tokens[1]]);
    }
    if row.v_from == row.v_to {
        return b"[2D_EDGE_CONVEYANCE] FROM_VERTEX and TO_VERTEX must differ".to_vec();
    }
    row.conveyance = try_parse_double(&tokens[2], &mut ok);
    if !ok {
        return cat(&[b"[2D_EDGE_CONVEYANCE] invalid CONVEYANCE: ", &tokens[2]]);
    }
    if row.conveyance < 0.0 || row.conveyance > 1.0 {
        return cat(&[b"[2D_EDGE_CONVEYANCE] CONVEYANCE must be in [0, 1] (got ", &tokens[2], b")"]);
    }
    pending_rows.push(row);
    Bytes::new()
}

/// parseInfil2DRowTail(): `METHOD [P1 .. P5] [DEST]` from tokens[first]; '-' leaves a slot.
fn parse_infil_2d_row_tail(tokens: &[Bytes], first: usize, section: &[u8], row: &mut Infil2DRow) -> Bytes {
    if tokens.len() <= first {
        return cat(&[section, b" missing METHOD"]);
    }
    if !parse_infil_2d_method(&tokens[first], &mut row.method, &mut row.has_method) {
        return cat(&[section, b" unknown METHOD: ", &tokens[first]]);
    }
    let mut end = tokens.len();
    if end > first + 1 {
        let last = &tokens[end - 1];
        let mut numeric = false;
        let _ = try_parse_double(last, &mut numeric);
        if !numeric && last.as_slice() != b"-" {
            if !parse_infil_2d_dest(last, &mut row.dest) {
                return cat(&[section, b" unknown DEST: ", last]);
            }
            row.dest_explicit = true;
            end -= 1;
        }
    }
    for k in (first + 1)..end {
        let idx = k - (first + 1);
        if idx >= K_INFIL_2D_MAX_PARAMS as usize {
            return cat(&[section, b" too many parameter columns (max ", K_INFIL_2D_MAX_PARAMS.to_string().as_bytes(), b")"]);
        }
        if tokens[k] == b"-" {
            continue; // unset
        }
        let mut ok = false;
        row.p[idx] = try_parse_double(&tokens[k], &mut ok);
        if !ok {
            return cat(&[section, b" invalid parameter: ", &tokens[k]]);
        }
    }
    Bytes::new()
}

/// makeInfilSectionHandler(): the [2D_INFILTRATION*] handlers write into the Infil2D reached
/// through ctx.twod_io.infil; none (a context without 2D) skips the section.
fn make_infil_section_handler(lp: fn(&[Bytes], &mut Infil2D) -> Bytes) -> impl Fn(&mut SimulationContext, &[Bytes]) + 'static {
    move |ctx: &mut SimulationContext, lines: &[Bytes]| {
        let infil = match &ctx.twod_io.infil {
            Some(i) => i.clone(),
            None => return,
        };
        for raw in lines {
            let tokens = tokenizer::tokenize(raw);
            if tokens.is_empty() {
                continue;
            }
            let err = lp(&tokens, &mut infil.borrow_mut());
            if !err.is_empty() {
                ctx.error_code = 5; // SWMM_ERR_PARSE
                ctx.error_message = cat(&[b"[2D] ", &err, " \u{2014} line: ".as_bytes(), raw]);
                return;
            }
        }
    }
}

/// parse2DInfiltrationOptionsLine(): `INFIL_STEP duration`.
pub fn parse_2d_infiltration_options_line(tokens: &[Bytes], infil: &mut Infil2D) -> Bytes {
    if tokens.is_empty() {
        return Bytes::new();
    }
    if tokens.len() < 2 {
        return b"[2D_INFILTRATION_OPTIONS] needs PARAMETER VALUE".to_vec();
    }
    if !iequals(&tokens[0], b"INFIL_STEP") {
        return cat(&[b"Unknown 2D_INFILTRATION_OPTIONS parameter: ", &tokens[0]]);
    }
    let secs = parse_time_seconds(&tokens[1]);
    if secs < 0.0 {
        return cat(&[b"[2D_INFILTRATION_OPTIONS] invalid INFIL_STEP (expected hh:mm:ss >= 0): ", &tokens[1]]);
    }
    infil.options_mut().infil_step = secs;
    Bytes::new()
}

/// parse2DInfiltrationDefaultsLine(): `TAG|* METHOD [P1..P5] [DEST]`.
pub fn parse_2d_infiltration_defaults_line(tokens: &[Bytes], infil: &mut Infil2D) -> Bytes {
    if tokens.is_empty() {
        return Bytes::new();
    }
    if tokens.len() < 2 {
        return b"[2D_INFILTRATION_DEFAULTS] needs TAG METHOD [P1..P5] [DEST]".to_vec();
    }
    let mut entry = Infil2DDefault::default();
    entry.tag = tokens[0].clone();
    let err = parse_infil_2d_row_tail(tokens, 1, b"[2D_INFILTRATION_DEFAULTS]", &mut entry.row);
    if !err.is_empty() {
        return err;
    }
    infil.defaults_mut().push(entry);
    Bytes::new()
}

/// parse2DInfiltrationLine(): `CELL METHOD [P1..P5] [DEST]`, CELL 1-based.
pub fn parse_2d_infiltration_line(tokens: &[Bytes], infil: &mut Infil2D) -> Bytes {
    if tokens.is_empty() {
        return Bytes::new();
    }
    if tokens.len() < 2 {
        return b"[2D_INFILTRATION] needs CELL METHOD [P1..P5] [DEST]".to_vec();
    }
    let mut ok = false;
    let cell = try_parse_int(&tokens[0], &mut ok);
    if !ok || cell < 1 {
        return cat(&[b"[2D_INFILTRATION] invalid CELL index (1-based): ", &tokens[0]]);
    }
    let mut entry = Infil2DOverride::default();
    entry.tri = cell - 1;
    let err = parse_infil_2d_row_tail(tokens, 1, b"[2D_INFILTRATION]", &mut entry.row);
    if !err.is_empty() {
        return err;
    }
    infil.overrides_mut().push(entry);
    Bytes::new()
}

/// iequals for the case-insensitive __TEMPERATURE__ test of the quality rows.
fn is_temperature(species: &[u8]) -> bool {
    iequals(species, b"__TEMPERATURE__")
}

/// parse2DInitialQualityLine(): `CELL <n> SP C | TAG <t> SP C | * SP C`.
pub fn parse_2d_initial_quality_line(tokens: &[Bytes], rows: &mut Vec<PendingInitialQualityRow>) -> Bytes {
    if tokens.is_empty() {
        return Bytes::new();
    }
    let mut r = PendingInitialQualityRow::default();
    let at;
    if tokens[0] == b"*" {
        if tokens.len() != 3 {
            return b"[2D_INITIAL_QUALITY] '*' row needs SPECIES CONC".to_vec();
        }
        r.all = true;
        at = 1;
    } else if iequals(&tokens[0], b"CELL") {
        if tokens.len() != 4 {
            return b"[2D_INITIAL_QUALITY] CELL row needs CELL <n> SPECIES CONC".to_vec();
        }
        let mut ok = false;
        let cell = try_parse_int(&tokens[1], &mut ok);
        if !ok || cell < 1 {
            return cat(&[b"[2D_INITIAL_QUALITY] invalid CELL index (1-based): ", &tokens[1]]);
        }
        r.tri = cell - 1;
        at = 2;
    } else if iequals(&tokens[0], b"TAG") {
        if tokens.len() != 4 {
            return b"[2D_INITIAL_QUALITY] TAG row needs TAG <name> SPECIES CONC".to_vec();
        }
        r.tag = tokens[1].clone();
        at = 2;
    } else {
        return cat(&[b"[2D_INITIAL_QUALITY] row must start with CELL, TAG or '*': ", &tokens[0]]);
    }
    r.species = tokens[at].clone();
    let mut okc = false;
    r.conc = try_parse_double(&tokens[at + 1], &mut okc);
    if !okc || !r.conc.is_finite() || (r.conc < 0.0 && !is_temperature(&r.species)) {
        return cat(&[b"[2D_INITIAL_QUALITY] CONC must be finite and non-negative except for __TEMPERATURE__, got '", &tokens[at + 1], b"'"]);
    }
    rows.push(r);
    Bytes::new()
}

/// parse2DInitialBedQualityLine(): the same scopes; CONC finite and non-negative.
pub fn parse_2d_initial_bed_quality_line(tokens: &[Bytes], rows: &mut Vec<PendingInitialBedQualityRow>) -> Bytes {
    if tokens.is_empty() {
        return Bytes::new();
    }
    let mut r = PendingInitialBedQualityRow::default();
    let at;
    if tokens[0] == b"*" {
        if tokens.len() != 3 {
            return b"[2D_INITIAL_BED_QUALITY] '*' row needs SPECIES CONC".to_vec();
        }
        r.all = true;
        at = 1;
    } else if iequals(&tokens[0], b"CELL") {
        if tokens.len() != 4 {
            return b"[2D_INITIAL_BED_QUALITY] CELL row needs CELL <n> SPECIES CONC".to_vec();
        }
        let mut ok = false;
        let cell = try_parse_int(&tokens[1], &mut ok);
        if !ok || cell < 1 {
            return cat(&[b"[2D_INITIAL_BED_QUALITY] invalid CELL index (1-based): ", &tokens[1]]);
        }
        r.tri = cell - 1;
        at = 2;
    } else if iequals(&tokens[0], b"TAG") {
        if tokens.len() != 4 {
            return b"[2D_INITIAL_BED_QUALITY] TAG row needs TAG <name> SPECIES CONC".to_vec();
        }
        r.tag = tokens[1].clone();
        at = 2;
    } else {
        return cat(&[b"[2D_INITIAL_BED_QUALITY] row must start with CELL, TAG or '*': ", &tokens[0]]);
    }
    r.species = tokens[at].clone();
    let mut okc = false;
    r.conc = try_parse_double(&tokens[at + 1], &mut okc);
    if !okc || !r.conc.is_finite() || r.conc < 0.0 {
        return cat(&[b"[2D_INITIAL_BED_QUALITY] CONC must be a finite non-negative number, got '", &tokens[at + 1], b"'"]);
    }
    rows.push(r);
    Bytes::new()
}

/// parse2DBoundaryQualityLine(): `TRI EDGE SPECIES CONC`.
pub fn parse_2d_boundary_quality_line(tokens: &[Bytes], rows: &mut Vec<PendingBoundaryQualityRow>) -> Bytes {
    if tokens.is_empty() {
        return Bytes::new();
    }
    if tokens.len() != 4 {
        return b"[2D_BOUNDARY_QUALITY] needs TRI EDGE SPECIES CONC".to_vec();
    }
    let mut ok = false;
    let mut r = PendingBoundaryQualityRow::default();
    r.tri = try_parse_int(&tokens[0], &mut ok);
    if !ok || r.tri < 0 {
        return cat(&[b"[2D_BOUNDARY_QUALITY] invalid TRI index: ", &tokens[0]]);
    }
    r.edge = try_parse_int(&tokens[1], &mut ok);
    if !ok || r.edge < 0 || r.edge >= K_MAX_CELL_VERTS {
        return cat(&[b"[2D_BOUNDARY_QUALITY] invalid EDGE (must be 0..2 for a triangle, 0..3 for a quad): ", &tokens[1]]);
    }
    r.species = tokens[2].clone();
    let mut okc = false;
    r.conc = try_parse_double(&tokens[3], &mut okc);
    if !okc || !r.conc.is_finite() || (r.conc < 0.0 && !is_temperature(&r.species)) {
        return cat(&[b"[2D_BOUNDARY_QUALITY] CONC must be finite and non-negative except for __TEMPERATURE__, got '", &tokens[3], b"'"]);
    }
    rows.push(r);
    Bytes::new()
}

/// The router members register2DSections writes through (upstream passes each by reference).
#[derive(Clone)]
pub struct Sections2DTargets {
    pub mesh: Shared<MeshData>,
    pub options: Shared<SolverOptions2D>,
    pub pending_bc_rows: Shared<Vec<PendingBoundaryRow>>,
    pub pending_ec_rows: Shared<Vec<PendingEdgeConveyanceRow>>,
    pub pending_iq_rows: Shared<Vec<PendingInitialQualityRow>>,
    pub pending_ibq_rows: Shared<Vec<PendingInitialBedQualityRow>>,
    pub pending_bq_rows: Shared<Vec<PendingBoundaryQualityRow>>,
}

/// register2DSections(): the mesh, options, native-quality, boundary and infiltration handlers.
pub fn register_2d_sections(t: &Sections2DTargets, registry: &mut SectionRegistry) {
    {
        let rows = t.pending_iq_rows.clone();
        registry.register_custom(b"2D_INITIAL_QUALITY", make_section_handler(move |tokens: &[Bytes]| parse_2d_initial_quality_line(tokens, &mut rows.borrow_mut())));
    }
    {
        let rows = t.pending_ibq_rows.clone();
        registry.register_custom(b"2D_INITIAL_BED_QUALITY", make_section_handler(move |tokens: &[Bytes]| parse_2d_initial_bed_quality_line(tokens, &mut rows.borrow_mut())));
    }
    {
        let rows = t.pending_bq_rows.clone();
        registry.register_custom(b"2D_BOUNDARY_QUALITY", make_section_handler(move |tokens: &[Bytes]| parse_2d_boundary_quality_line(tokens, &mut rows.borrow_mut())));
    }
    {
        // Full-form handler: parse2DOptionsLine needs ctx.warnings for the retired keys.
        let options = t.options.clone();
        registry.register_custom(b"2D_OPTIONS", move |ctx: &mut SimulationContext, lines: &[Bytes]| {
            for raw in lines {
                let tokens = tokenizer::tokenize(raw);
                if tokens.is_empty() {
                    continue;
                }
                let err = parse_2d_options_line(&tokens, &mut options.borrow_mut(), Some(&mut ctx.warnings));
                if !err.is_empty() {
                    ctx.error_code = 5; // SWMM_ERR_PARSE
                    ctx.error_message = cat(&[b"[2D] ", &err, " \u{2014} line: ".as_bytes(), raw]);
                    return;
                }
            }
        });
    }
    {
        let mesh = t.mesh.clone();
        registry.register_custom(b"2D_VERTICES", make_section_handler(move |tokens: &[Bytes]| parse_2d_vertex_line(tokens, &mut mesh.borrow_mut())));
    }
    {
        let mesh = t.mesh.clone();
        registry.register_custom(b"2D_TRIANGLES", make_section_handler(move |tokens: &[Bytes]| parse_2d_triangle_line(tokens, &mut mesh.borrow_mut())));
    }
    {
        let mesh = t.mesh.clone();
        registry.register_custom(b"2D_QUADS", make_section_handler(move |tokens: &[Bytes]| parse_2d_quad_line(tokens, &mut mesh.borrow_mut())));
    }
    {
        let mesh = t.mesh.clone();
        registry.register_custom(b"2D_INITIAL_VELOCITY", make_section_handler(move |tokens: &[Bytes]| parse_2d_initial_velocity_line(tokens, &mut mesh.borrow_mut())));
    }
    {
        let mesh = t.mesh.clone();
        registry.register_custom(b"2D_VERTEX_NODE_MAP", make_section_handler(move |tokens: &[Bytes]| parse_2d_vertex_node_map_line(tokens, &mut mesh.borrow_mut())));
    }
    {
        let mesh = t.mesh.clone();
        registry.register_custom(b"2D_TRIANGLE_NODE_MAP", make_section_handler(move |tokens: &[Bytes]| parse_2d_triangle_node_map_line(tokens, &mut mesh.borrow_mut())));
    }
    {
        let rows = t.pending_bc_rows.clone();
        registry.register_custom(b"2D_BOUNDARY_CONDITIONS", make_section_handler(move |tokens: &[Bytes]| parse_2d_boundary_conditions_line(tokens, &mut rows.borrow_mut())));
    }
    {
        let rows = t.pending_ec_rows.clone();
        registry.register_custom(b"2D_EDGE_CONVEYANCE", make_section_handler(move |tokens: &[Bytes]| parse_2d_edge_conveyance_line(tokens, &mut rows.borrow_mut())));
    }
    registry.register_custom(b"2D_INFILTRATION_OPTIONS", make_infil_section_handler(parse_2d_infiltration_options_line));
    registry.register_custom(b"2D_INFILTRATION_DEFAULTS", make_infil_section_handler(parse_2d_infiltration_defaults_line));
    registry.register_custom(b"2D_INFILTRATION", make_infil_section_handler(parse_2d_infiltration_line));
    {
        // [2D_MESH_FILE]: the first FILE row only; the mesh is loaded after the .inp is read.
        let options = t.options.clone();
        registry.register_custom(b"2D_MESH_FILE", move |_ctx: &mut SimulationContext, lines: &[Bytes]| {
            for raw in lines {
                let tokens = tokenizer::tokenize(raw);
                if tokens.len() >= 2 && iequals(&tokens[0], b"FILE") {
                    // FILE rows consume the raw tail so unquoted names keep repeated spaces.
                    let line = tokenizer::trim(tokenizer::strip_comment(raw));
                    let mut path: Bytes = tokens[1].clone();
                    for k in 2..tokens.len() {
                        path.push(b' ');
                        path.extend_from_slice(&tokens[k]);
                    }
                    if line.len() > 4 && c_isspace(line[4]) {
                        let tail = tokenizer::trim(&line[4..]);
                        if !tail.is_empty() && tail[0] != b'"' && tail[0] != b',' {
                            path = tail.to_vec();
                        }
                    }
                    options.borrow_mut().mesh_file.set(&path);
                    return; // only the first FILE line
                }
            }
        });
    }
}

/// load2DMeshExternalFile(): reads a mesh sidecar through its own registry into the targets.
/// Returns "" on success, otherwise the diagnostic. `warnings` (the caller's warnings) receives
/// the retired-option warnings of the sidecar's [2D_OPTIONS]; none makes them hard errors.
pub fn load_2d_mesh_external_file(
    mesh: &Shared<MeshData>,
    opts: &Shared<SolverOptions2D>,
    pending_bc_rows: &Shared<Vec<PendingBoundaryRow>>,
    pending_ec_rows: &Shared<Vec<PendingEdgeConveyanceRow>>,
    infil: Option<&Shared<Infil2D>>,
    mesh_file: &[u8],
    inp_base_dir: &[u8],
    warnings: Option<&mut Vec<Bytes>>,
    pending_iq_rows: Option<&Shared<Vec<PendingInitialQualityRow>>>,
    pending_ibq_rows: Option<&Shared<Vec<PendingInitialBedQualityRow>>>,
    pending_bq_rows: Option<&Shared<Vec<PendingBoundaryQualityRow>>>,
) -> Bytes {
    // Resolve the path.
    let mut p = FsPath::new(mesh_file);
    if !p.is_absolute() && !inp_base_dir.is_empty() {
        p = FsPath::new(inp_base_dir).join(&p);
    }

    // A minimal registry (no 2D_MESH_FILE: no recursion).
    let mut mini = SectionRegistry::new();
    let sink: Option<Shared<Vec<Bytes>>> = if warnings.is_some() { Some(Rc::new(RefCell::new(Vec::new()))) } else { None };
    {
        let opts = opts.clone();
        let sink = sink.clone();
        mini.register_custom(
            b"2D_OPTIONS",
            make_section_handler(move |tokens: &[Bytes]| match &sink {
                Some(s) => parse_2d_options_line(tokens, &mut opts.borrow_mut(), Some(&mut s.borrow_mut())),
                None => parse_2d_options_line(tokens, &mut opts.borrow_mut(), None),
            }),
        );
    }
    {
        let mesh = mesh.clone();
        mini.register_custom(b"2D_VERTICES", make_section_handler(move |tokens: &[Bytes]| parse_2d_vertex_line(tokens, &mut mesh.borrow_mut())));
    }
    {
        let mesh = mesh.clone();
        mini.register_custom(b"2D_TRIANGLES", make_section_handler(move |tokens: &[Bytes]| parse_2d_triangle_line(tokens, &mut mesh.borrow_mut())));
    }
    {
        let mesh = mesh.clone();
        mini.register_custom(b"2D_QUADS", make_section_handler(move |tokens: &[Bytes]| parse_2d_quad_line(tokens, &mut mesh.borrow_mut())));
    }
    {
        let mesh = mesh.clone();
        mini.register_custom(b"2D_INITIAL_VELOCITY", make_section_handler(move |tokens: &[Bytes]| parse_2d_initial_velocity_line(tokens, &mut mesh.borrow_mut())));
    }
    {
        let mesh = mesh.clone();
        mini.register_custom(b"2D_VERTEX_NODE_MAP", make_section_handler(move |tokens: &[Bytes]| parse_2d_vertex_node_map_line(tokens, &mut mesh.borrow_mut())));
    }
    {
        let mesh = mesh.clone();
        mini.register_custom(b"2D_TRIANGLE_NODE_MAP", make_section_handler(move |tokens: &[Bytes]| parse_2d_triangle_node_map_line(tokens, &mut mesh.borrow_mut())));
    }
    {
        let rows = pending_bc_rows.clone();
        mini.register_custom(b"2D_BOUNDARY_CONDITIONS", make_section_handler(move |tokens: &[Bytes]| parse_2d_boundary_conditions_line(tokens, &mut rows.borrow_mut())));
    }
    // The surface quality sections travel with the mesh they address, when the caller passes
    // their stores.
    if let Some(rows) = pending_iq_rows {
        let rows = rows.clone();
        mini.register_custom(b"2D_INITIAL_QUALITY", make_section_handler(move |tokens: &[Bytes]| parse_2d_initial_quality_line(tokens, &mut rows.borrow_mut())));
    }
    if let Some(rows) = pending_ibq_rows {
        let rows = rows.clone();
        mini.register_custom(b"2D_INITIAL_BED_QUALITY", make_section_handler(move |tokens: &[Bytes]| parse_2d_initial_bed_quality_line(tokens, &mut rows.borrow_mut())));
    }
    if let Some(rows) = pending_bq_rows {
        let rows = rows.clone();
        mini.register_custom(b"2D_BOUNDARY_QUALITY", make_section_handler(move |tokens: &[Bytes]| parse_2d_boundary_quality_line(tokens, &mut rows.borrow_mut())));
    }
    {
        let rows = pending_ec_rows.clone();
        mini.register_custom(b"2D_EDGE_CONVEYANCE", make_section_handler(move |tokens: &[Bytes]| parse_2d_edge_conveyance_line(tokens, &mut rows.borrow_mut())));
    }
    // [2D_INFILTRATION*] are parsed into a scratch object and transplanted after the read: a
    // section the sidecar carries replaces the inline one.
    let sidecar_infil: Shared<Infil2D> = Rc::new(RefCell::new(Infil2D::default()));
    {
        let s = sidecar_infil.clone();
        mini.register_custom(b"2D_INFILTRATION_OPTIONS", make_section_handler(move |tokens: &[Bytes]| parse_2d_infiltration_options_line(tokens, &mut s.borrow_mut())));
    }
    {
        let s = sidecar_infil.clone();
        mini.register_custom(b"2D_INFILTRATION_DEFAULTS", make_section_handler(move |tokens: &[Bytes]| parse_2d_infiltration_defaults_line(tokens, &mut s.borrow_mut())));
    }
    {
        let s = sidecar_infil.clone();
        mini.register_custom(b"2D_INFILTRATION", make_section_handler(move |tokens: &[Bytes]| parse_2d_infiltration_line(tokens, &mut s.borrow_mut())));
    }

    // The sidecar may carry its own `;; UNITS:` header.
    let p_utf8: Bytes = p.native().to_vec();
    prescan_2d_units_header(&p_utf8, &mut opts.borrow_mut());

    let ok = {
        let mut reader = InputReader::new(&mini);
        let mut dummy = SimulationContext::default();
        let ok = reader.read(&p_utf8, &mut dummy);
        if let (Some(w), Some(s)) = (warnings, &sink) {
            w.append(&mut s.borrow_mut());
        }
        if !ok {
            Err(dummy.error_message.clone())
        } else {
            Ok(())
        }
    };
    if let Err(msg) = ok {
        return cat(&[b"2D_MESH_FILE: error reading '", &p_utf8, b"': ", &msg]);
    }

    if let Some(infil) = infil {
        let mut side = sidecar_infil.borrow_mut();
        let mut target = infil.borrow_mut();
        if !side.defaults().is_empty() {
            *target.defaults_mut() = std::mem::take(side.defaults_mut());
        }
        if !side.overrides().is_empty() {
            *target.overrides_mut() = std::mem::take(side.overrides_mut());
        }
        if side.options().infil_step > 0.0 {
            *target.options_mut() = side.options().clone();
        }
    }
    Bytes::new()
}

/// prescan2DUnitsHeader(): sets mesh_units_si on a `;; UNITS: SI (m)` (or metric synonym)
/// comment anywhere in the file; never clears it. A missing file returns silently.
pub fn prescan_2d_units_header(inp_path: &[u8], opts: &mut SolverOptions2D) {
    let mut input = match IfStream::open(inp_path) {
        Some(i) => i,
        None => return, // file missing: the caller surfaces the error
    };
    fn trim(s: &[u8]) -> &[u8] {
        let mut a = 0;
        let mut b = s.len();
        while b > a && c_isspace(s[b - 1]) {
            b -= 1;
        }
        while a < b && c_isspace(s[a]) {
            a += 1;
        }
        &s[a..b]
    }
    let mut line: Bytes = Bytes::new();
    while input.getline(&mut line) {
        let t = trim(&line);
        if t.len() < 2 || t[0] != b';' || t[1] != b';' {
            continue;
        }
        let rest = trim(&t[2..]);
        const K_KEY: &[u8] = b"UNITS:";
        if rest.len() < K_KEY.len() {
            continue;
        }
        let mut matched = true;
        for i in 0..K_KEY.len() {
            if rest[i].to_ascii_uppercase() != K_KEY[i] {
                matched = false;
                break;
            }
        }
        if !matched {
            continue;
        }
        let value = trim(&rest[K_KEY.len()..]);
        let si = iequals(value, b"SI (m)") || iequals(value, b"m") || iequals(value, b"metre") || iequals(value, b"metres") || iequals(value, b"meter") || iequals(value, b"meters");
        if si {
            opts.mesh_units_si = true;
        }
    }
}
