// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/SWMMEngine.cpp (the OPENSWMM_HAS_2D blocks of the
// constructor, wire2DModelIO and open()) (HydroCouple/openswmm.engine, branch swmm6_rel, commit
// 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
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
// The get_snow_state / set_snow_state accessors' mesh snow pointer is not here: the state
// accessors are not translated (core/context.rs), and the mesh snow state is a later wave's.

#![allow(clippy::all)]

use crate::core::engine::SwmmEngine;
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
