// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/SWMMEngine.cpp (SWMMEngine::start,
// openRunoffIfaceWrite, openRunoffIfaceRead, closeRunoffIface and the file-local
// plugin_io_error) (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// SWMMEngine::start (1D build: OPENSWMM_HAS_2D and OPENSWMM_HAS_2D_HDF5_OUTPUT off): the
// one-report-period clamp, DoRouting, the [FILES] interface files (routing USE INFLOWS / SAVE
// OUTFLOWS, runoff SAVE/USE, RDII USE/SAVE, the RAINFALL slot warning), the plugins'
// initialize / validate / prepare, the averaging accumulators, the mass-balance NodeInflow
// seed, the dated SAVE HOTSTART flags and the RUNNING state. With the `twod` feature (the 2D
// build, 12 November 2026) an active 2D router also turns DoRouting on, and an independent
// subsurface clock asks for the dated hot-start snapshot.
//
// Differences of form, never in result:
//   * io_thread_.start(): the writer thread only moves output writes off the stepping thread;
//     there is none here.
//   * perf:: timers only time the run; they are left out.
//   * prepareNativeHotstartCheckpointing and loadCheckpointBytes (a native USE HOTSTART
//     checkpoint, or SAVE HOTSTART rows under coupled heat) are NOT YET TRANSLATED: start
//     returns Err(stub) there. start is noexcept upstream; Ok(rc) is its return code.

use crate::core::cmath::cvt_i32;
use crate::core::context::{EngineState, FileMode};
use crate::core::engine::{SwmmEngine, SWMM_ERR_IO, SWMM_ERR_PLUGIN, SWMM_ERR_WRONG_STATE, SWMM_OK};
use crate::core::engine_init_hydraulics::Terminated;
use crate::core::engine_step::stub;
use crate::core::error_codes::{
    format_error, format_warning, ERR_RDII_IFACE, ERR_RDII_IFACE_FORMAT, ERR_ROUTING_IFACE, ERR_ROUTING_IFACE_FORMAT, ERR_ROUTING_IFACE_SAME, ERR_RUNOFF_IFACE, ERR_RUNOFF_IFACE_COMPAT,
    WARN_FILES_SLOT_UNSUPPORTED,
};
use crate::core::options::{Bytes, FilePathPair};
use crate::core::plugins::PluginFactory;
use crate::hydrology::rdii_interface::rdii_node_list;
use crate::hydrology::runoff_interface::RunoffInterfaceFile;

/// The resolved absolute path, else the verbatim token (IO-3).
fn chosen(p: &FilePathPair) -> Bytes {
    if !p.absolute.is_empty() {
        p.absolute.clone()
    } else {
        p.original.clone()
    }
}

/// plugin_io_error: the first plugin's last error after the fallback, or the fallback.
fn plugin_io_error(plugins: &PluginFactory, fallback: &[u8]) -> Bytes {
    let mut details: Vec<Bytes> = Vec::new();
    details.push(plugins.output_plugins().iter().map(|p| crate::core::cstdlib::c_str(p.last_error_message()).to_vec()).find(|m| !m.is_empty()).unwrap_or_default());
    details.push(plugins.report_plugins().iter().map(|p| crate::core::cstdlib::c_str(p.last_error_message()).to_vec()).find(|m| !m.is_empty()).unwrap_or_default());
    details.push(plugins.state_io_plugins().iter().map(|p| crate::core::cstdlib::c_str(p.last_error_message()).to_vec()).find(|m| !m.is_empty()).unwrap_or_default());
    for m in details {
        if !m.is_empty() {
            let mut s = fallback.to_vec();
            s.extend_from_slice(b": ");
            s.extend_from_slice(&m);
            return s;
        }
    }
    fallback.to_vec()
}

impl SwmmEngine {
    /// SWMMEngine::openRunoffIfaceWrite.
    pub(crate) fn open_runoff_iface_write(&mut self, path: &[u8]) -> i32 {
        if self.runoff_iface_file_.as_ref().map(|f| f.is_open()).unwrap_or(false) {
            return -10;
        }
        let mut f = RunoffInterfaceFile::default();
        let rc = f.open_for_write(path, self.ctx_.n_subcatches(), self.ctx_.n_pollutants(), self.ctx_.options.flow_units as i32);
        self.runoff_iface_file_ = if rc != 0 { None } else { Some(f) };
        rc
    }

    /// SWMMEngine::openRunoffIfaceRead.
    pub(crate) fn open_runoff_iface_read(&mut self, path: &[u8]) -> i32 {
        if self.runoff_iface_file_.as_ref().map(|f| f.is_open()).unwrap_or(false) {
            return -10;
        }
        let mut f = RunoffInterfaceFile::default();
        let rc = f.open_for_read(path, self.ctx_.n_subcatches(), self.ctx_.n_pollutants(), self.ctx_.options.flow_units as i32);
        self.runoff_iface_file_ = if rc != 0 { None } else { Some(f) };
        rc
    }

    /// SWMMEngine::closeRunoffIface (end() and the C API call it).
    #[allow(dead_code)]
    pub(crate) fn close_runoff_iface(&mut self) {
        if let Some(mut f) = self.runoff_iface_file_.take() {
            f.close();
        }
    }

    /// SWMMEngine::start(save_results).
    pub fn start(&mut self, save_results: i32) -> Result<i32, Terminated> {
        if self.ctx_.state != EngineState::Initialized {
            self.set_error(SWMM_ERR_WRONG_STATE, b"swmm_engine_start: must call initialize() first");
            return Ok(SWMM_ERR_WRONG_STATE);
        }
        self.save_results_ = save_results;

        // Short runs still need one reporting period (legacy project_validate).
        let duration = self.ctx_.options.total_duration_ms() / 1000.0;
        if duration > 0.0 && self.ctx_.options.report_step > duration {
            self.ctx_.options.report_step = duration;
            self.ctx_.next_report_ms = 1000.0 * duration;
        }

        // Legacy DoRouting (swmm5.c:748)
        self.do_routing_ = self.ctx_.n_nodes() > 0 && !self.ctx_.options.ignore_routing;
        // Surface-only engines still need the routing clock and surface advance.
        #[cfg(feature = "twod")]
        {
            self.do_routing_ = self.do_routing_ || (self.surface_router_.is_active() && !self.ctx_.options.ignore_routing);
        }

        // Routing interface files ([FILES] USE INFLOWS / SAVE OUTFLOWS).
        if !self.ctx_.files.inflows_path.is_empty() || !self.ctx_.files.outflows_path.is_empty() {
            let in_path = chosen(&self.ctx_.files.inflows_path);
            let out_path = chosen(&self.ctx_.files.outflows_path);
            let rc = self.iface_.open_files(&in_path, &out_path);
            if rc == -1 {
                self.ctx_.errors.push(format_error(ERR_ROUTING_IFACE_SAME, b""));
                self.set_error(SWMM_ERR_IO, b"swmm_engine_start: inflows and outflows routing interface files have the same name");
                return Ok(SWMM_ERR_IO);
            }
            if rc == -2 {
                self.ctx_.errors.push(format_error(ERR_ROUTING_IFACE, &out_path));
                self.set_error(SWMM_ERR_IO, b"swmm_engine_start: cannot open routing interface file for writing");
                return Ok(SWMM_ERR_IO);
            }
            if rc == -3 {
                self.iface_.close_files();
                self.ctx_.errors.push(format_error(ERR_ROUTING_IFACE, &in_path));
                self.set_error(SWMM_ERR_IO, b"swmm_engine_start: cannot open routing interface file for reading");
                return Ok(SWMM_ERR_IO);
            }
            if !in_path.is_empty() && self.iface_.read_file_header(&mut self.ctx_) != 0 {
                self.iface_.close_files();
                self.ctx_.errors.push(format_error(ERR_ROUTING_IFACE_FORMAT, &in_path));
                self.set_error(SWMM_ERR_IO, b"swmm_engine_start: invalid format for routing interface file");
                return Ok(SWMM_ERR_IO);
            }
            if !out_path.is_empty() {
                self.iface_.write_file_header(&self.ctx_);
            }
        }

        // [FILES] SAVE/USE RUNOFF (legacy runoff_open -> openRunoffFile).
        if self.ctx_.n_subcatches() > 0
            && self.ctx_.files.runoff_mode != FileMode::None
            && !self.ctx_.files.runoff_path.is_empty()
            && !self.runoff_iface_file_.as_ref().map(|f| f.is_open()).unwrap_or(false)
        {
            let rp = chosen(&self.ctx_.files.runoff_path);
            let rc = if self.ctx_.files.runoff_mode == FileMode::Save { self.open_runoff_iface_write(&rp) } else { self.open_runoff_iface_read(&rp) };
            if rc != 0 {
                self.ctx_.errors.push(format_error(if rc == -3 { ERR_RUNOFF_IFACE_COMPAT } else { ERR_RUNOFF_IFACE }, &rp));
                self.set_error(SWMM_ERR_IO, b"swmm_engine_start: cannot open runoff interface file");
                return Ok(SWMM_ERR_IO);
            }
        }

        // [FILES] USE/SAVE RDII (legacy: rain_open() gates rdii_openRdii(), swmm5.c:735).
        if self.ctx_.files.rdii_mode != FileMode::None && !self.ctx_.files.rdii_path.is_empty() && !self.ctx_.options.ignore_rdii && !self.ctx_.options.ignore_rainfall {
            let dp = chosen(&self.ctx_.files.rdii_path);
            if self.ctx_.files.rdii_mode == FileMode::Use {
                let rc = self.rdii_iface_file_.open_for_read(&mut self.ctx_, &dp);
                if rc != 0 {
                    self.ctx_.errors.push(format_error(if rc == -2 { ERR_RDII_IFACE_FORMAT } else { ERR_RDII_IFACE }, &dp));
                    self.set_error(SWMM_ERR_IO, b"swmm_engine_start: cannot open RDII interface file");
                    return Ok(SWMM_ERR_IO);
                }
            } else {
                let rdii_nodes = rdii_node_list(&self.rdii_);
                let rc = self.rdii_iface_file_.open_for_write(&dp, cvt_i32(self.ctx_.options.wet_step), &rdii_nodes);
                if rc == -1 {
                    self.ctx_.errors.push(format_error(ERR_RDII_IFACE, &dp));
                    self.set_error(SWMM_ERR_IO, b"swmm_engine_start: cannot open RDII interface file for writing");
                    return Ok(SWMM_ERR_IO);
                }
                // rc == -2: no RDII inflows; the slot stays inert, as legacy's does.
            }
        }

        // [FILES] USE/SAVE RAINFALL: not implemented upstream; a warning.
        if self.ctx_.files.rainfall_mode != FileMode::None && !self.ctx_.files.rainfall_path.is_empty() {
            let what: &[u8] = if self.ctx_.files.rainfall_mode == FileMode::Use { b"USE RAINFALL" } else { b"SAVE RAINFALL" };
            self.ctx_.warnings.push(format_warning(WARN_FILES_SLOT_UNSUPPORTED, what));
            self.emit_warning(SWMM_ERR_IO, b"[FILES] RAINFALL interface file is not supported by this engine and was ignored");
        }

        // Phase 4: prepare() on all plugins (opens output files/headers)
        if !self.plugins_.is_empty() {
            self.plugins_prepare_attempted_ = true;
            if self.plugins_.initialize_all(&mut self.ctx_) != 0 || self.plugins_.validate_all(&mut self.ctx_) != 0 {
                let m = plugin_io_error(&self.plugins_, b"Selected plugin initialization or validation failed");
                self.set_error(SWMM_ERR_PLUGIN, &m);
                return Ok(SWMM_ERR_PLUGIN);
            }
            let rc = self.plugins_.prepare_all(&mut self.ctx_);
            if rc != 0 {
                let m = plugin_io_error(&self.plugins_, b"swmm_engine_start: plugin prepare() failed");
                self.set_error(SWMM_ERR_PLUGIN, &m);
                return Ok(SWMM_ERR_PLUGIN);
            }
        }

        // Phase 5: io_thread_.start() (no writer thread here)

        if self.ctx_.options.rpt_averages {
            let (nn, nl) = (self.ctx_.n_nodes(), self.ctx_.n_links());
            self.avg_.resize(nn, nl);
        }

        // legacy massbal_open: NodeInflow starts at each node's stored volume.
        for j in 0..self.ctx_.n_nodes() as usize {
            self.ctx_.nodes.stat_total_inflow_vol[j] = self.ctx_.nodes.volume[j];
        }

        // No dated SAVE HOTSTART row has been written in this run yet.
        self.hotstart_save_done_.clear();
        self.hotstart_save_done_.resize(self.ctx_.files.hotstart_saves.len(), 0);

        self.ctx_.state = EngineState::Running;
        #[allow(unused_mut)]
        let mut date_native_snapshot = self.ctx_.options.heat_transport;
        #[cfg(feature = "twod")]
        {
            date_native_snapshot |= self.surface_router_.subsurface().independent_clock();
        }
        let wants_native_hotstart = !self.pending_hotstart_checkpoint_.is_empty() || (date_native_snapshot && !self.ctx_.files.hotstart_saves.is_empty());
        if wants_native_hotstart {
            return Err(stub("SWMMEngine::prepareNativeHotstartCheckpointing"));
        }
        Ok(SWMM_OK)
    }
}
