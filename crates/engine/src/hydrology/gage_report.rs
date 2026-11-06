// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/Gage.cpp (getReportRainfall,
// getReportRainfallFromSeries and the file-local reportRainfall, gageRainSeries,
// gageSeriesRainType, gageUnitsFactor, recordRate, currentInterval) (HydroCouple/
// openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin, licensed
// under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// A gage's rainfall at a report instant (legacy gage_setReportRainfall), read-only on the gage
// state the last runoff step left. The report writer's (SWMMEngine::postOutputSnapshot) use of
// Gage.cpp. The helpers are private copies of the Gage.cpp functions the runoff step's
// translation (hydrology/gage.rs) also carries; they translate the same C++ lines.

use crate::core::context::SimulationContext;
use crate::core::datetime;
use crate::data::gage::{RainFileFormat, RainSource};
use crate::data::table::Table;

const K_NO_STATE: i32 = -2;
const K_PRE_RECORD: i32 = -1;

fn gage_rain_series(ctx: &SimulationContext, gage_idx: i32) -> Option<&Table> {
    if gage_idx < 0 || gage_idx >= ctx.n_gages() {
        return None;
    }
    let ug = gage_idx as usize;
    if ctx.gages.source[ug] == RainSource::FileRain {
        if ug < ctx.gages.rain_series.len() && !ctx.gages.rain_series[ug].x.is_empty() {
            return Some(&ctx.gages.rain_series[ug]);
        }
        return None;
    }
    let ts_idx = ctx.gages.ts_index[ug];
    if ts_idx >= 0 && ts_idx < ctx.tables.tables.len() as i32 {
        return Some(&ctx.tables.tables[ts_idx as usize]);
    }
    None
}

fn gage_series_rain_type(ctx: &SimulationContext, gage_idx: i32) -> i32 {
    let ug = gage_idx as usize;
    if ctx.gages.source[ug] == RainSource::FileRain && ctx.gages.file_format[ug] == RainFileFormat::StanPrcp {
        return 1;
    }
    ctx.gages.rain_type[ug]
}

fn gage_units_factor(ctx: &SimulationContext, gage_idx: i32) -> f64 {
    if gage_idx < 0 || gage_idx >= ctx.n_gages() {
        return 1.0;
    }
    let ug = gage_idx as usize;
    if ctx.gages.source[ug] == RainSource::FileRain && ctx.gages.file_format[ug] == RainFileFormat::StanPrcp && ctx.options.flow_units as i32 >= 3 {
        return 25.40;
    }
    1.0
}

fn record_rate(ctx: &SimulationContext, j: i32, tbl: &Table, k: i32) -> f64 {
    let uk = k as usize;
    let uj = j as usize;
    let rain_type = gage_series_rain_type(ctx, j);
    let interval = ctx.gages.interval_sec[uj] as f64;
    let mut r = tbl.y[uk];
    if rain_type == 1 && interval > 0.0 {
        r = r / interval * 3600.0;
    } else if rain_type == 2 && interval > 0.0 {
        let prev = if k > 0 { tbl.y[uk - 1] } else { 0.0 };
        r = if r < prev { r / interval * 3600.0 } else { (r - prev) / interval * 3600.0 };
    }
    r * gage_units_factor(ctx, j) * ctx.gages.scale_factor[uj]
}

fn current_interval(ctx: &SimulationContext, j: i32, tbl: &Table) -> (f64, f64) {
    let uj = j as usize;
    let cur = ctx.gages.st_cur[uj];
    if cur == K_PRE_RECORD {
        (ctx.options.start_date, tbl.x[0])
    } else {
        let start = tbl.x[cur as usize];
        (start, datetime::add_seconds(start, ctx.gages.interval_sec[uj] as f64))
    }
}

fn report_rainfall(ctx: &SimulationContext, gage_idx: i32, report_date: f64, from_series: bool) -> f64 {
    if ctx.options.ignore_rainfall {
        return 0.0;
    }
    let ug = gage_idx as usize;

    // Co-gage: the primary's report rainfall re-scaled by this gage's factor
    let co = if ug < ctx.gages.co_gage_index.len() { ctx.gages.co_gage_index[ug] } else { -1 };
    if co >= 0 && co < gage_idx {
        let uco = co as usize;
        let primary_sf = ctx.gages.scale_factor[uco];
        let this_sf = ctx.gages.scale_factor[ug];
        let ratio = if primary_sf > 0.0 { this_sf / primary_sf } else { 1.0 };
        return report_rainfall(ctx, co, report_date, from_series) * ratio;
    }

    // API override
    if ctx.gages.api_rainfall[ug] >= 0.0 {
        return ctx.gages.api_rainfall[ug];
    }

    let rtbl = match gage_rain_series(ctx, gage_idx) {
        Some(t) if !t.x.is_empty() => t,
        _ => return 0.0,
    };
    if ug >= ctx.gages.st_cur.len() || ctx.gages.st_cur[ug] == K_NO_STATE {
        return 0.0;
    }
    let t = report_date + datetime::ONE_SECOND;

    let mon = datetime::month_of_year(report_date) - 1;
    let radj = if (0..12).contains(&mon) { ctx.adjust_rain[mon as usize] } else { 1.0 };

    if from_series {
        // std::upper_bound(x, t): the first element greater than t
        let it = rtbl.x.partition_point(|&v| !(t < v));
        if it == 0 {
            return 0.0;
        }
        let k = it as i32 - 1;
        let end = datetime::add_seconds(rtbl.x[k as usize], ctx.gages.interval_sec[ug] as f64);
        return if t < end { record_rate(ctx, gage_idx, rtbl, k) * radj } else { 0.0 };
    }

    let (_start, end) = current_interval(ctx, gage_idx, rtbl);
    if t < end {
        return ctx.gages.st_rain[ug];
    }
    let nxt = ctx.gages.st_next[ug];
    if nxt < 0 {
        return 0.0;
    }
    if t < rtbl.x[nxt as usize] {
        return 0.0;
    }
    ctx.gages.next_rainfall[ug]
}

/// gage::getReportRainfall: the rainfall a gage reports at `report_date`.
pub fn get_report_rainfall(ctx: &SimulationContext, gage_idx: i32, report_date: f64) -> f64 {
    report_rainfall(ctx, gage_idx, report_date, false)
}

/// gage::getReportRainfallFromSeries: the same, read from the series rather than the cursor.
pub fn get_report_rainfall_from_series(ctx: &SimulationContext, gage_idx: i32, report_date: f64) -> f64 {
    report_rainfall(ctx, gage_idx, report_date, true)
}
