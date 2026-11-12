// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/SWMMEngine.cpp (fillSurfaceSnapshot and the
// OPENSWMM_HAS_2D blocks of postOutputSnapshot, end, writeHotstartSave, report and close) and
// src/engine/core/SWMMEngineCheckpoint.cpp (detachEngineSourceHistory)
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin,
// licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// Compiled only with the `twod` feature (phase 2a wave 5). core/engine_output.rs calls these at
// upstream's points:
//   * fill_surface_snapshot: the snapshot's 2D surface block, at each report period
//     (postOutputSnapshot, after the batch flush and the render refresh);
//   * detach_engine_source_history: end() (before the groundwater storage) and close();
//   * the router's finalize (twod/surface_router_2d_end.rs): end(), before the max-stats
//     finalize;
//   * report_mesh_snow_rows: report(), the mesh snow rows of the snow balance;
//   * hotstart_independent_clock: writeHotstartSave's full-continuation test;
//   * close_2d: close()'s resets.
//
// Differences in form, never in result:
//   * Without the HDF5 2D writer (the reference build is OPENSWMM_2D_HDF5_OUTPUT=OFF) there is
//     no surface_output_plugin_, so postOutputSnapshot's surface-only snapshot branch (2D output
//     without a .out) cannot run; the surface block is filled into the snapshot the 1D writers
//     receive, which none of them reads.
//   * The engine members only the checkpoint and source-delivery paths set
//     (checkpoint_surface_output_, runtime_water_driver_, surface_output_plugin_,
//     engine_source_history_, engine_source_stream_, checkpoint_file_outputs_,
//     source_delivery_, source_stage_, source_totals_) are not carried: no CLI run attaches
//     them (engine checkpointing and the runtime source API are not translated), so they stay
//     null/false/zero, which is what close() resets them to.
//   * A phase-2b owner the snapshot or the report would read (mesh snow, surface species, the
//     surface quality buildup, the two-zone aquifer) is refused by SurfaceRouter2D::initialize
//     before a run can reach these; each is guarded here with the crate's stub.

#![allow(clippy::all)]

use crate::core::engine::SwmmEngine;
use crate::core::engine_init_hydraulics::Terminated;
use crate::core::engine_step::stub;
use crate::core::output_plugin::SimulationSnapshot;
use crate::twod::data::mesh_data::MeshData;

impl SwmmEngine {
    /// fillSurfaceSnapshot(snap): the 2D surface state into the snapshot (SI-native), the
    /// depth and head blended to the report instant. Upstream's is noexcept with no error path;
    /// a phase-2b owner (unreachable after initialize's refusals) is a "STUB:" panic.
    pub(crate) fn fill_surface_snapshot(&self, snap: &mut SimulationSnapshot) {
        let r = &self.surface_router_;
        if !r.is_active() {
            return;
        }
        let st = r.state();
        let msh_rc = r.mesh().clone();
        let msh = msh_rc.borrow();
        let s = &mut snap.surface;
        s.tri_count = msh.n_triangles();
        s.vert_count = msh.n_vertices();
        s.quad_count = msh.n_quads();
        s.edge_stride = msh.edge_stride();
        if s.quad_count > 0 {
            s.cell_nv = msh.cell_nv.clone();
        } else {
            s.cell_nv.clear();
        }
        s.depth = st.depth.clone();
        s.head = st.head.clone();
        // Report-time blend: depth and head to the exact report instant.
        {
            let f = r.report_weight(self.ctx_.next_report_ms);
            if f < 1.0 {
                let f1 = 1.0 - f;
                let od = r.report_old_depth();
                let oh = r.report_old_head();
                let n = std::cmp::min(s.depth.len(), od.len());
                for i in 0..n {
                    s.depth[i] = f1 * od[i] + f * s.depth[i];
                    s.head[i] = f1 * oh[i] + f * s.head[i];
                }
            }
        }
        s.grad_hx = st.grad_hx.clone();
        s.grad_hy = st.grad_hy.clone();
        s.grad_hx_lim = st.grad_hx_lim.clone();
        s.grad_hy_lim = st.grad_hy_lim.clone();
        // Rainfall as the .out reports subcatchment rainfall: the gage rate at the report instant.
        r.report_rainfall(&self.ctx_, snap.sim_time, &mut s.rainfall);
        for field in s.snow.iter_mut() {
            field.clear();
        }
        if st.snow.active() {
            panic!("STUB: snow on the 2D mesh in the report snapshot (SWMMEngine::fillSurfaceSnapshot; phase 2b)");
        }

        s.coupling_flux = st.coupling_flux.clone();
        s.net_source = st.net_source.clone();
        s.infil_rate = st.infil_rate.clone();
        {
            let tr = &st.transport;
            s.species_count = tr.n_species;
            if tr.active() {
                panic!("STUB: 2D surface species in the report snapshot (SWMMEngine::fillSurfaceSnapshot; phase 2b)");
            } else {
                s.species_conc.clear();
            }
            let sq_rc = r.surface_quality().clone();
            let sq = sq_rc.borrow();
            if sq.active() {
                panic!("STUB: 2D surface quality buildup in the report snapshot (SWMMEngine::fillSurfaceSnapshot; phase 2b)");
            } else {
                s.buildup_count = 0;
                s.buildup.clear();
            }
        }
        let infil_cum = r.infil_cumulative();
        if infil_cum.len() == st.infil_rate.len() {
            s.infil_cum = infil_cum.clone();
        } else {
            s.infil_cum.clear();
            s.infil_cum.resize(st.infil_rate.len(), 0.0);
        }
        s.rain_cum = r.report_rain_cumulative(self.ctx_.next_report_ms).clone();
        // Output sign convention: outward-positive edge flux and down-gradient velocity.
        {
            let stride = s.edge_stride;
            let nc = s.tri_count;
            s.edge_flux.clear();
            s.edge_flux.resize(nc as usize * stride as usize, 0.0);
            for c in 0..nc {
                for k in 0..stride {
                    s.edge_flux[c as usize * stride as usize + k as usize] = -st.edge_flux[MeshData::slot(c, k) as usize];
                }
            }
        }
        s.vert_head = st.vert_head.clone();
        s.vert_depth = st.vert_depth_signed.clone();
        s.face_vx = st.face_vx.clone();
        s.face_vy = st.face_vy.clone();
        for v in s.face_vx.iter_mut() {
            *v = -*v;
        }
        for v in s.face_vy.iter_mut() {
            *v = -*v;
        }
        s.continuity_err = st.cell_continuity_err.clone();
        s.stat_max_depth = st.stat_max_depth.clone();
        s.stat_max_velocity = st.stat_max_velocity.clone();
        s.stat_max_cont_err = st.stat_max_cont_err.clone();
        // G-O: the two-zone aquifer's per-cell state (phase 2b).
        {
            let gw = &r.subsurface_;
            s.gw2d_active = gw.active();
            if gw.active() {
                panic!("STUB: the 2D two-zone aquifer in the report snapshot (SWMMEngine::fillSurfaceSnapshot; phase 2b)");
            }
        }
    }

    /// detachEngineSourceHistory(): the CPU marcher's engine source history and stream are
    /// detached (the engine's own pointers are not carried; see the file notes).
    pub(crate) fn detach_engine_source_history(&mut self) {
        if let Some(cpu) = self.surface_router_.solver_.as_mut() {
            cpu.engine_source_history_ = false;
            cpu.engine_source_stream_ = false;
        }
    }

    /// end()'s finalize of the 2D module, before the max-stats finalize. end() is noexcept: an
    /// exception from the flush becomes a termination.
    pub(crate) fn end_2d_finalize(&mut self) -> Result<(), Terminated> {
        self.surface_router_.finalize(&mut self.ctx_).map_err(|t| match t.exception.strip_prefix(crate::core::engine_step::THROWN_PREFIX) {
            Some(ty) => Terminated { exception: ty, what: t.what },
            None => t,
        })
    }

    /// report()'s mesh snow rows of the snow balance (phase 2b: snow on the mesh is refused at
    /// initialize).
    pub(crate) fn report_mesh_snow_rows(&mut self) -> Result<(), Terminated> {
        if self.surface_router_.state().snow.active() {
            return Err(stub("the mesh snow rows of the snow balance (SWMMEngine::report; phase 2b)"));
        }
        Ok(())
    }

    /// writeHotstartSave's `surface_router_.subsurface().independentClock()`.
    pub(crate) fn hotstart_independent_clock(&self) -> bool {
        self.surface_router_.subsurface_.independent_clock()
    }

    /// close()'s OPENSWMM_HAS_2D resets.
    pub(crate) fn close_2d(&mut self) {
        // checkpoint_surface_output_, runtime_water_driver_, surface_output_plugin_: not carried
        self.detach_engine_source_history();
        // source_delivery_, source_stage_, source_totals_: not carried
    }
}

// Gate only (cargo features `twod` and `gate`; examples/snap2d.rs): a hook that sees each
// report period's snapshot right after fill_surface_snapshot, as gate/diff/snap2d's
// interposed fillSurfaceSnapshot does upstream's. Not part of the translation.
#[cfg(feature = "gate")]
thread_local! {
    static GATE_SNAPSHOT_HOOK: std::cell::RefCell<Option<Box<dyn FnMut(&SimulationSnapshot)>>> = std::cell::RefCell::new(None);
}

/// Gate only: install (or remove) the snapshot hook.
#[cfg(feature = "gate")]
pub fn gate_set_snapshot_hook(hook: Option<Box<dyn FnMut(&SimulationSnapshot)>>) {
    GATE_SNAPSHOT_HOOK.with(|h| *h.borrow_mut() = hook);
}

/// Gate only: call the snapshot hook, when one is installed.
#[cfg(feature = "gate")]
pub(crate) fn gate_snapshot_hook(snap: &SimulationSnapshot) {
    GATE_SNAPSHOT_HOOK.with(|h| {
        if let Some(f) = h.borrow_mut().as_mut() {
            f(snap);
        }
    });
}
