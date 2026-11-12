// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/SWMMEngine.cpp (the OPENSWMM_HAS_2D blocks of the
// constructor, wire2DModelIO, open(), initialize(), init_modules(), initHydraulics() and
// start()) (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Compiled only with the `twod` feature. open() (core/engine.rs) calls these at upstream's
// points:
//   * open_2d_clear: `surface_router_.gwTransport().clear()` after ctx_.reset();
//   * open_2d_register_and_prescan: the 2D, [GW_*], surface-quality and subsurface section
//     registration into the default input plugin, then the `;; UNITS: SI` prescan of the .inp;
//   * open_2d_after_read: the external mesh file ([2D_MESH_FILE]) and the INFIL_STEP / GW_ET
//     alias folding;
//   * open_2d_resolve: resolveGwTransport and its two warnings, after msxsurf::resolve;
//   * open_2d_output_file_warning: [2D_OPTIONS] OUTPUT_FILE without the HDF5 writer (the
//     reference build is OPENSWMM_2D_HDF5_OUTPUT=OFF, as this translation is).
// and the later lifecycle (wave 3, 12 November 2026):
//   * init_hydraulics_2d: initHydraulics' 1a block, surface_router_.initialize(ctx_) and the two
//     follow-up checks inside upstream's try; its catch's "2D initialization failed: <what>";
//   * init_hydraulics_2d_has_domain: the ERR_NO_OUTLETS exemption of a model with a mesh;
//   * the router's thread warnings, appended to initHydraulics' (core/engine_init_hydraulics.rs);
//   * init_modules_channel_ground: init_modules' initializeChannelGround and its catch;
//   * initialize_2d_reactive_lid_check: initialize's check after init_modules;
//   * start's do_routing_ for a surface-only engine (core/engine_start.rs).
// A path the translation refuses inside SurfaceRouter2D::initialize (phase 2b) returns
// Err(Terminated) with the crate's STUB_EXCEPTION from initialize().
// The get_snow_state / set_snow_state accessors' mesh snow pointer is not here: the state
// accessors are not translated (core/context.rs), and the mesh snow state is a later wave's.

#![allow(clippy::all)]

use crate::core::engine::{SwmmEngine, CFFI_ERR_BADPARAM};
use crate::core::engine_init_hydraulics::Terminated;
use crate::data::lid_node::LidNodeBottom;
use crate::twod::surface_router_2d::InitError;
use crate::core::fs_path::FsPath;
use crate::core::options::Bytes;
use crate::twod::gw::gw_transport_sections::{gw_transport_inert_warning, register_gw_transport_sections, resolve_gw_transport};
use crate::twod::input::section_handlers_2d::{load_2d_mesh_external_file, prescan_2d_units_header, register_2d_sections};
use crate::twod::quality::surface_quality_2d::register_surface_quality_sections;
use crate::twod::subsurface::subsurface_sections::register_subsurface_sections;
use crate::twod::surface_router_2d::SurfaceRouter2D;

/// SWMM_ERR_PARSE
const SWMM_ERR_PARSE: i32 = 5;

impl SwmmEngine {
    /// surfaceRouter2D()
    pub fn surface_router_2d(&self) -> &SurfaceRouter2D {
        &self.surface_router_
    }

    /// wire2DModelIO(): the context's 2D handles point at the router's members.
    pub(crate) fn wire_2d_model_io(&mut self) {
        let r = &self.surface_router_;
        let io = &mut self.ctx_.twod_io;
        io.mesh = Some(r.mesh().clone());
        io.options = Some(r.options().clone());
        io.boundary = Some(r.boundary().clone());
        io.pending_bc = Some(r.pending_bc_rows().clone());
        io.pending_ec = Some(r.pending_edge_conveyance_rows().clone());
        io.pending_iq = Some(r.pending_initial_quality_rows().clone());
        io.pending_ibq = Some(r.pending_initial_bed_quality_rows().clone());
        io.pending_bq = Some(r.pending_boundary_quality_rows().clone());
        io.infil = Some(r.infil().clone());
        io.gw = Some(r.gw_transport().clone());
        io.surface_quality = Some(r.surface_quality().clone());
        io.aquifer = Some(r.aquifer_config().clone());
        io.aquifer_nodes = Some(r.aquifer_node_names().clone());
        io.aquifer_links = Some(r.aquifer_link_names().clone());
        io.surface_router = true;
    }

    /// Router-owned groundwater rows outlive ctx_.reset(): a fresh open replaces them.
    pub(crate) fn open_2d_clear(&mut self) {
        self.surface_router_.gw_transport().borrow_mut().clear();
    }

    /// Register the 2D section handlers before reading, then scan the .inp for the
    /// `;; UNITS: SI (m)` header.
    pub(crate) fn open_2d_register_and_prescan(&mut self, inp_path: Option<&[u8]>) {
        {
            let r = &self.surface_router_;
            let reg = self.plugins_.input_plugins_mut()[0].registry();
            register_2d_sections(&r.sections_2d_targets(), reg);
            register_gw_transport_sections(r.gw_transport(), reg);
            register_surface_quality_sections(r.surface_quality(), reg);
            register_subsurface_sections(r.aquifer_config(), r.aquifer_node_names(), r.aquifer_link_names(), reg);
        }
        if let Some(p) = inp_path {
            if !p.is_empty() {
                prescan_2d_units_header(p, &mut self.surface_router_.options().borrow_mut());
            }
        }
    }

    /// After the read: load the [2D_MESH_FILE] sidecar, then fold the [2D_OPTIONS] INFIL_STEP and
    /// GW_ET aliases into their owners. Some(code) when open() must return it.
    pub(crate) fn open_2d_after_read(&mut self, inp_path: Option<&[u8]>) -> Option<i32> {
        let mf: Bytes = self.surface_router_.options().borrow().mesh_file.original.clone();
        if !mf.is_empty() {
            let mut base_dir: Bytes = Bytes::new();
            if let Some(p) = inp_path {
                if !p.is_empty() {
                    base_dir = FsPath::new(p).parent_path().native().to_vec();
                }
            }
            let err = {
                let r = &self.surface_router_;
                load_2d_mesh_external_file(
                    r.mesh(),
                    r.options(),
                    r.pending_bc_rows(),
                    r.pending_edge_conveyance_rows(),
                    Some(r.infil()),
                    &mf,
                    &base_dir,
                    Some(&mut self.ctx_.warnings),
                    Some(r.pending_initial_quality_rows()),
                    Some(r.pending_initial_bed_quality_rows()),
                    Some(r.pending_boundary_quality_rows()),
                )
            };
            if !err.is_empty() {
                if self.lenient_open() && !self.metadata_only_ {
                    // Discard the partly read mesh, keep the diagnostic, continue 1D-only.
                    *self.surface_router_.mesh().borrow_mut() = Default::default();
                    self.surface_router_.pending_bc_rows().borrow_mut().clear();
                    self.surface_router_.pending_edge_conveyance_rows().borrow_mut().clear();
                    self.ctx_.errors.push(err);
                } else {
                    self.ctx_.error_code = SWMM_ERR_PARSE;
                    self.ctx_.error_message = err.clone();
                    self.ctx_.errors.push(err);
                    self.write_open_failure_report();
                    return Some(SWMM_ERR_PARSE);
                }
            }
        }
        // E2: [2D_OPTIONS] INFIL_STEP is an alias of [2D_INFILTRATION_OPTIONS] INFIL_STEP; the
        // alias wins when both were spelled.
        {
            let r = &self.surface_router_;
            let mut o2d = r.options().borrow_mut();
            if o2d.infil_step > 0.0 {
                r.infil().borrow_mut().options_mut().infil_step = o2d.infil_step;
                o2d.infil_step = 0.0;
            }
            // U5: [2D_OPTIONS] GW_ET is the same kind of alias, for [2D_AQUIFER_OPTIONS] GW_ET;
            // it counts as authoring the section.
            if !o2d.gw_et.is_empty() {
                let mut gw = r.aquifer_config().borrow_mut();
                gw.options.gw_et = o2d.gw_et.clone();
                gw.options.authored = true;
                o2d.gw_et.clear();
            }
        }
        None
    }

    /// U4: validate the [GW_*] rows once the mesh, species and series are final, then the
    /// "authored but inert" and GROUNDWATER YES warnings.
    pub(crate) fn open_2d_resolve(&mut self) {
        let gw_errs = {
            let r = &self.surface_router_;
            let mesh = r.mesh().borrow();
            let mut gw = r.gw_transport().borrow_mut();
            resolve_gw_transport(&self.ctx_, &mesh, &mut gw)
        };
        for e in gw_errs {
            self.ctx_.errors.push(e);
        }
        if self.surface_router_.aquifer_config().borrow().empty() {
            let inert = gw_transport_inert_warning(&self.surface_router_.gw_transport().borrow());
            if !inert.is_empty() {
                self.push_report_warning(&inert, 0);
            }
        }
        let groundwater = self.surface_router_.options().borrow().groundwater;
        if groundwater > 0 && self.surface_router_.aquifer_config().borrow().empty() {
            self.push_report_warning(
                "[2D_OPTIONS] GROUNDWATER YES but no [2D_AQUIFER_OPTIONS] / [2D_AQUIFER] row was authored \u{2014} there is no subsurface to run this simulation.".as_bytes(),
                0,
            );
        }
    }

    /// Built without the HDF5 results writer: [2D_OPTIONS] OUTPUT_FILE draws a warning and the
    /// 2D results stay in the engine.
    pub(crate) fn open_2d_output_file_warning(&mut self) {
        let of: Bytes = self.surface_router_.options().borrow().output_file.original.clone();
        if !of.is_empty() {
            let mut m = b"[2D_OPTIONS] OUTPUT_FILE ".to_vec();
            m.extend_from_slice(&of);
            m.extend_from_slice(b" ignored: this engine was built without the 2D HDF5 results writer; read 2D results through the swmm_2d_* API");
            self.push_report_warning(&m, 0);
        }
    }
}

impl SwmmEngine {
    /// initHydraulics' has_2d_domain: the mesh is parsed at open, so its counts are valid even
    /// though surface_router_.initialize() (which sets isActive()) runs later.
    pub(crate) fn init_hydraulics_2d_has_domain(&self) -> bool {
        let m = self.surface_router_.mesh().borrow();
        m.n_triangles() >= 1 && m.n_vertices() >= 3
    }

    /// initHydraulics 1a: the 2D surface routing module. initialize() and the two checks after
    /// it run inside upstream's try; its catch reports "2D initialization failed: <what>"
    /// through the error path. Err only for a path the translation refuses (a stub).
    pub(crate) fn init_hydraulics_2d(&mut self) -> Result<(), Terminated> {
        let r = self.init_hydraulics_2d_try();
        match r {
            Ok(()) => Ok(()),
            Err(InitError::Thrown(what)) => {
                let mut m = b"2D initialization failed: ".to_vec();
                m.extend_from_slice(&what);
                self.ctx_.errors.push(m);
                let e = self.ctx_.errors.last().cloned().unwrap_or_default();
                self.set_error(SWMM_ERR_PARSE, &e);
                Ok(())
            }
            Err(InitError::Stub(t)) => Err(t),
        }
    }

    fn init_hydraulics_2d_try(&mut self) -> Result<(), InitError> {
        self.surface_router_.initialize(&mut self.ctx_)?;
        let groundwater = self.surface_router_.options().borrow().groundwater;
        if !self.ctx_.options.ignore_2d
            && groundwater != 0
            && !self.surface_router_.subsurface().active()
            && (!self.ctx_.reactions.subsurface_laws.is_empty() || !self.surface_router_.aquifer_config().borrow().boundaries.is_empty())
        {
            return Err(InitError::Thrown(b"Groundwater edge forcing and REACTION_SUBSURFACE require an active resolved mesh aquifer.".to_vec()));
        }
        for n in 0..self.ctx_.n_nodes() {
            if crate::hydrology::lid_node_step::bottom_mode(&self.ctx_, n) == LidNodeBottom::MeshAquifer
                && (!self.surface_router_.is_active()
                    || !self.surface_router_.subsurface().active()
                    || n as usize >= self.ctx_.nodes.aquifer2d_bed.len()
                    || self.ctx_.nodes.aquifer2d_bed[n as usize] == 0)
            {
                let mut m = b"Mesh LID '".to_vec();
                m.extend_from_slice(self.ctx_.node_names.name_of(n));
                m.extend_from_slice(b"' requires an active resolved 2D aquifer bottom connection.");
                return Err(InitError::Thrown(m));
            }
        }
        Ok(())
    }

    /// init_modules: surface_router_.initializeChannelGround(ctx_) and its catch. Err only for
    /// a path the translation refuses (a stub).
    pub(crate) fn init_modules_channel_ground(&mut self) -> Result<(), Terminated> {
        match self.surface_router_.initialize_channel_ground(&mut self.ctx_) {
            Ok(()) => Ok(()),
            Err(InitError::Thrown(what)) => {
                let mut m = b"Channel ground heat initialization failed: ".to_vec();
                m.extend_from_slice(&what);
                self.ctx_.errors.push(m);
                let e = self.ctx_.errors.last().cloned().unwrap_or_default();
                self.set_error(CFFI_ERR_BADPARAM, &e);
                Ok(())
            }
            Err(InitError::Stub(t)) => Err(t),
        }
    }

    /// initialize(), after init_modules(): reactive conventional LIDs need ordinary runoff
    /// ownership, which a 2D source owner does not give.
    pub(crate) fn initialize_2d_reactive_lid_check(&mut self) {
        if self.surface_router_.is_active() && self.ctx_.reactions.lid.subcatchments.iter().any(|&v| v != 0) {
            self.ctx_.errors.push(b"Reactive conventional LIDs require ordinary runoff ownership; 2D source-owner species receipts are not implemented.".to_vec());
        }
    }
}
