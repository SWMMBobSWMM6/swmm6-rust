// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/InterfaceFile.hpp and InterfaceFile.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Routing interface files: the manager's state and InterfaceManager::init.
//
// NOT YET TRANSLATED (run-time functions, with start() and the routing step):
// InterfaceFileSoA::resize, openFiles, readFileHeader, writeFileHeader, readInflows,
// writeOutfallResults, closeFiles (beyond what init needs), getFlow, getQuality,
// getIfaceNode, readNextPeriod, setOldValues, isOutletNode.
//
// Differences in form, never in result:
//   * The two `FILE*` handles are `Option<std::fs::File>`; closeFiles drops them.
//   * Members keep upstream's names and are `pub` (upstream's private ones too) so the
//     gate's state dumper can read them.

use crate::core::context::SimulationContext;

#[derive(Debug, Clone, Default)]
pub struct InterfaceFileSoA {
    pub count: i32,
    pub node_idx: Vec<i32>,
    pub old_values: Vec<f64>,
    pub new_values: Vec<f64>,
}

#[derive(Debug)]
pub struct InterfaceManager {
    pub data_: InterfaceFileSoA,
    pub infile_: Option<std::fs::File>,
    pub outfile_: Option<std::fs::File>,
    pub old_time_: f64,
    pub new_time_: f64,
    pub n_iface_polluts_: i32,
    pub n_values_per_node_: i32,
    pub iface_flow_units_: i32,
    pub iface_step_: i32,
    pub iface_frac_: f64,
    pub pollut_map_: Vec<i32>,
}

impl Default for InterfaceManager {
    fn default() -> Self {
        InterfaceManager {
            data_: InterfaceFileSoA::default(),
            infile_: None,
            outfile_: None,
            old_time_: 0.0,
            new_time_: 0.0,
            n_iface_polluts_: 0,
            n_values_per_node_: 1,
            iface_flow_units_: 0,
            iface_step_: 0,
            iface_frac_: 0.0,
            pollut_map_: Vec::new(),
        }
    }
}

impl InterfaceManager {
    pub fn init(&mut self, _ctx: &mut SimulationContext) {
        // Re-init safety: release any files left open by a previous run.
        self.close_files();
        self.data_.count = 0;
        self.infile_ = None;
        self.outfile_ = None;
        self.old_time_ = 0.0;
        self.new_time_ = 0.0;
        self.n_iface_polluts_ = 0;
        self.n_values_per_node_ = 1;
        self.iface_frac_ = 0.0;
        self.pollut_map_.clear();
    }

    /// closeFiles: close and forget both handles.
    pub fn close_files(&mut self) {
        self.infile_ = None;
        self.outfile_ = None;
    }

    pub fn get_num_iface_nodes(&self) -> i32 {
        self.data_.count
    }

    pub fn get_frac(&self) -> f64 {
        self.iface_frac_
    }
}
