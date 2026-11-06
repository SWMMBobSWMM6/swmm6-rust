// SPDX-License-Identifier: Apache-2.0
//
// STUB, written 5 November 2026 by Claude (Anthropic) at the request of Robert Dickinson, so
// that SWMMEngine::start and step (core/engine_start.rs, engine_step.rs) compile and the
// state gates run. The translation of the output and report side (OpenSWMM Engine
// src/engine/core/SWMMEngine.cpp postOutputSnapshot, accumulateAvgResults,
// writeDueHotstartSaves; plugins/DefaultOutputPlugin.cpp writeHeader;
// plugins/DefaultReportPlugin.cpp write_preamble; HydroCouple/openswmm.engine, branch
// swmm6_rel, commit 88329d6, Copyright 2026 Caleb Buahin, Apache-2.0) REPLACES THIS WHOLE FILE.
//
// The signatures start() and step() call (see PLAN.md):
//   SwmmEngine::post_output_snapshot(&mut self, dt_step: f64) -> Result<(), Terminated>
//   SwmmEngine::accumulate_avg_results(&mut self) -> Result<(), Terminated>
//   SwmmEngine::write_due_hotstart_saves(&mut self) -> Result<(), Terminated>
//   DefaultOutputPlugin::write_header(&mut self, ctx: &SimulationContext)
//   DefaultReportPlugin::write_preamble(&mut self, ctx: &SimulationContext)
// write_header and write_preamble are called from DefaultOutputPlugin::prepare and
// DefaultReportPlugin::prepare (core/plugins.rs) inside start(). They write file content
// only (the .out header, the .rpt preamble), which no state gate sees; so these two stubs
// do nothing, and start() runs through, leaving the two files without their headers. The
// others return Err(stub) so the step harness reports where it stopped.

use crate::core::context::SimulationContext;
use crate::core::engine::SwmmEngine;
use crate::core::engine_init_hydraulics::Terminated;
use crate::core::engine_step::stub;
use crate::core::plugins::{DefaultOutputPlugin, DefaultReportPlugin};

impl SwmmEngine {
    /// STUB of SWMMEngine::postOutputSnapshot.
    pub(crate) fn post_output_snapshot(&mut self, _dt_step: f64) -> Result<(), Terminated> {
        Err(stub("SWMMEngine::postOutputSnapshot"))
    }
    /// STUB of SWMMEngine::accumulateAvgResults.
    pub(crate) fn accumulate_avg_results(&mut self) -> Result<(), Terminated> {
        Err(stub("SWMMEngine::accumulateAvgResults"))
    }
    /// STUB of SWMMEngine::writeDueHotstartSaves.
    pub(crate) fn write_due_hotstart_saves(&mut self) -> Result<(), Terminated> {
        Err(stub("SWMMEngine::writeDueHotstartSaves"))
    }
}

impl DefaultOutputPlugin {
    /// STUB of DefaultOutputPlugin::writeHeader: writes nothing (see above).
    pub(crate) fn write_header(&mut self, _ctx: &SimulationContext) {}
}

impl DefaultReportPlugin {
    /// STUB of DefaultReportPlugin::write_preamble: writes nothing (see above).
    pub(crate) fn write_preamble(&mut self, _ctx: &SimulationContext) {}
}
