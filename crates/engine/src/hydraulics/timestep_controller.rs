// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/TimestepController.hpp and
// TimestepController.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The routing-step clock (static methods of TimestepController upstream; free functions
// here): compute_next, advance, reset_output_timer, output_due, simulation_complete.

use crate::core::cmath::{cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::datetime;

pub const SEC_PER_DAY: f64 = 86400.0;

/// compute_next(ctx, dt_cfl): the next routing step (legacy execRouting +
/// routing_getRoutingStep): the CFL step capped at ROUTING_STEP, the duration clamp
/// (1 ms floor), the control-rule clamp; never zero or negative.
pub fn compute_next(ctx: &SimulationContext, dt_cfl: f64) -> f64 {
    let opt = &ctx.options;
    let min_step = opt.min_routing_step;
    let mut dt = cmin(opt.routing_step, dt_cfl);
    let total_msec = opt.total_duration_ms();
    if ctx.elapsed_ms + 1000.0 * dt > total_msec {
        dt = (total_msec - ctx.elapsed_ms) / 1000.0;
        dt = cmax(dt, 0.001);
    }
    if ctx.dt_controls_remaining > 0.0 && ctx.dt_controls_remaining < dt {
        dt = ctx.dt_controls_remaining;
    }
    if dt <= 0.0 {
        dt = min_step;
    }
    dt
}

/// advance(ctx, dt_taken): the clocks after a step (routing.c:308-309 for the ms clock).
pub fn advance(ctx: &mut SimulationContext, dt_taken: f64) {
    ctx.current_time += dt_taken;
    ctx.current_date = datetime::add_seconds(ctx.options.start_date, ctx.current_time);
    ctx.old_elapsed_ms = ctx.elapsed_ms;
    ctx.elapsed_ms = ctx.elapsed_ms + 1000.0 * dt_taken;
    if ctx.dt_controls_remaining > 0.0 {
        ctx.dt_controls_remaining -= dt_taken;
        if ctx.dt_controls_remaining < 0.0 {
            ctx.dt_controls_remaining = 0.0;
        }
    }
}

/// reset_output_timer: ReportTime += 1000 * ReportStep (swmm5.c:1044).
pub fn reset_output_timer(ctx: &mut SimulationContext) {
    ctx.next_report_ms = ctx.next_report_ms + 1000.0 * ctx.options.report_step;
}

/// output_due: NewRoutingTime >= ReportTime (swmm5.c:1019).
pub fn output_due(ctx: &SimulationContext) -> bool {
    ctx.elapsed_ms >= ctx.next_report_ms
}

/// simulation_complete: NewRoutingTime >= RoutingDuration (ms).
pub fn simulation_complete(ctx: &SimulationContext) -> bool {
    ctx.elapsed_ms >= ctx.options.total_duration_ms()
}
