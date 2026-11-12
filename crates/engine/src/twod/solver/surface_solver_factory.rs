// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/solver/SurfaceSolverFactory.hpp and
// SurfaceSolverFactory.cpp (makeSurfaceSolver) (HydroCouple/openswmm.engine, branch swmm6_rel,
// commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version
// 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The surface solver selection: the backend policy ([2D_OPTIONS] BACKEND, the AUTO mesh-size
// floors) and its notices, for a build without the GPU/Kokkos plugin. Upstream's plugin
// loader searches OPENSWMM_GPU_PLUGIN_PATH, the engine library's own folder and its gpu/
// subfolder for libopenswmm_gpu_<backend>.so; the reference build (OPENSWMM_BUILD_GPU_PLUGIN
// OFF) installs none, so every try_plugin finds nothing and the serial CPU marcher
// (ExplicitInertialSolver) serves every model, as here: try_plugin is the empty search.
// Every notice goes to stderr, as upstream's fprintf(stderr, ...).
//
// Refused (the crate's stub, where upstream reads them; PLAN.md, the environment-switch
// policy): the development switches OPENSWMM_2D_BACKEND, OPENSWMM_2D_MIN_PARALLEL_CELLS and
// OPENSWMM_2D_MIN_PARALLEL_CELLS_DEVICE, and OPENSWMM_GPU_PLUGIN_PATH (a plugin folder the
// translation cannot load from).

#![allow(clippy::all)]

use crate::cformat;
use crate::core::engine_step::stub;
use crate::core::options::Bytes;
use crate::twod::data::solver_options_2d::{Backend2D, Momentum2D, SolverOptions2D};
use crate::twod::solver::explicit_inertial_solver::{eprint_bytes, ExplicitInertialSolver};
use crate::twod::surface_router_2d::InitError;

fn backend_name(b: Backend2D) -> &'static [u8] {
    match b {
        Backend2D::Cpu => b"cpu",
        Backend2D::Omp => b"omp",
        Backend2D::Cuda => b"cuda",
        Backend2D::Hip => b"hip",
        Backend2D::Sycl => b"sycl",
        Backend2D::Auto => b"auto",
    }
}

fn env_set(name: &str) -> bool {
    std::env::var_os(name).is_some()
}

/// try_plugin(backend, chosen, threads): no plugin is installed (see the header).
fn try_plugin(_backend: &[u8]) -> Result<Option<ExplicitInertialSolver>, InitError> {
    if env_set("OPENSWMM_GPU_PLUGIN_PATH") {
        return Err(InitError::Stub(stub("OPENSWMM_GPU_PLUGIN_PATH (the 2D GPU/Kokkos plugin search, SurfaceSolverFactory)")));
    }
    Ok(None)
}

/// makeSurfaceSolver(opts, chosen, n_cells, plugin_capable).
pub fn make_surface_solver(opts: &SolverOptions2D, chosen: &mut Bytes, n_cells: i32, plugin_capable: bool) -> Result<ExplicitInertialSolver, InitError> {
    let mut plugin_capable = plugin_capable;
    if opts.momentum_physics() {
        if opts.momentum != Momentum2D::FullSwe || opts.lts_tiers != 1 {
            return Err(InitError::Thrown(b"Coriolis/Smagorinsky/wall drag require FULL_SWE and LTS_TIERS 1.".to_vec()));
        }
        // Never hand an enabled process to a plugin that has no such kernel.
        plugin_capable = false;
    }
    let serial_marcher = |chosen: &mut Bytes, label: &[u8]| -> ExplicitInertialSolver {
        *chosen = label.to_vec();
        ExplicitInertialSolver::default()
    };
    // Policy: the environment variable wins; otherwise the model's [2D_OPTIONS] BACKEND.
    if env_set("OPENSWMM_2D_BACKEND") {
        return Err(InitError::Stub(stub("the OPENSWMM_2D_BACKEND development switch (makeSurfaceSolver)")));
    }
    let mode: &[u8] = if opts.backend == Backend2D::Auto && opts.auto_backend_cpu { b"cpu" } else { backend_name(opts.backend) };
    if mode == b"cpu" {
        return Ok(serial_marcher(chosen, b"cpu (explicit marcher)"));
    }
    // The Kokkos plugin implements the triangle-only LOCAL_INERTIAL scheme without transport
    // rows: anything else is refused loudly and served by the CPU marcher.
    if !plugin_capable {
        let m = cformat!(
            b"[openswmm 2D] backend '%s': the GPU/Kokkos plugin serves only all-triangle meshes under MOMENTUM_EQUATION LOCAL_INERTIAL with no transport rows (species / age / temperature) or coupled groundwater/porous heat; using the CPU marcher for this model.\n",
            if mode.is_empty() { &b"auto"[..] } else { mode }
        );
        eprint_bytes(&m);
        return Ok(serial_marcher(chosen, b"cpu (explicit marcher; plugin not applicable)"));
    }
    // An explicit backend request bypasses every mesh-size gate.
    if mode == b"omp" || mode == b"cuda" || mode == b"hip" || mode == b"sycl" {
        if let Some(s) = try_plugin(mode)? {
            return Ok(s);
        }
        eprint_bytes(&cformat!(b"[openswmm 2D] backend '%s' requested but no usable plugin/device found; falling back to the CPU marcher.\n", mode));
        return Ok(serial_marcher(chosen, b"cpu (explicit marcher)"));
    }
    // auto: separate mesh-size floors for the device and the host (OpenMP) plugins.
    if env_set("OPENSWMM_2D_MIN_PARALLEL_CELLS_DEVICE") || env_set("OPENSWMM_2D_MIN_PARALLEL_CELLS") {
        return Err(InitError::Stub(stub("the OPENSWMM_2D_MIN_PARALLEL_CELLS(_DEVICE) development switches (makeSurfaceSolver)")));
    }
    let dev_floor: i64 = 10000;
    let omp_floor: i64 = 50000;
    if n_cells <= 0 || n_cells as i64 >= dev_floor {
        for b in [&b"cuda"[..], b"hip", b"sycl"] {
            if let Some(s) = try_plugin(b)? {
                return Ok(s);
            }
        }
    }
    if n_cells <= 0 || n_cells as i64 >= omp_floor {
        if let Some(s) = try_plugin(b"omp")? {
            return Ok(s);
        }
    }
    eprint_bytes(&cformat!(
        b"[openswmm 2D] using CPU marcher (no GPU device selected; n_cells=%d, device floor=%ld cells). Set [2D_OPTIONS] BACKEND CUDA (or OPENSWMM_2D_BACKEND=cuda) to force the GPU, or lower OPENSWMM_2D_MIN_PARALLEL_CELLS_DEVICE.\n",
        n_cells,
        dev_floor
    ));
    Ok(serial_marcher(chosen, b"cpu (explicit marcher)"))
}
