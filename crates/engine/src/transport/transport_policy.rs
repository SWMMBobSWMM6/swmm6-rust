// SPDX-License-Identifier: Apache-2.0
//
// STUB (5 November 2026, Claude for Robert Dickinson): transport/TransportPolicy.cpp
// (resolve, formatReportBlock) is not translated yet. The .rpt preamble prints its block only
// when [OPTIONS] WATER_AGE or HEAT_TRANSPORT is on or a reactions component is compiled.

use crate::core::context::SimulationContext;
use crate::core::options::Bytes;

/// formatReportBlock(resolve(ctx)): STUB.
pub fn report_block(_ctx: &SimulationContext) -> Bytes {
    panic!("STUB: transport::formatReportBlock(transport::resolve(ctx))")
}
