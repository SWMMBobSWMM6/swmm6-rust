// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/Gage.hpp and Gage.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The parts of the rain-gage module translated so far: gageIsUsed (SWMMEngine::open
// reaches it), and the constants, RainType and GageState, which SWMMEngine::initHydrology
// sizes (gage_states_); and (5 November 2026, for the runoff step) splitPrecip,
// updatePastRain (the context form), gageRainSeries, gageSeriesRainType, gageUnitsFactor,
// the legacy rain-gage state machine (recordRate, nextNonzeroRecord, initGageState,
// currentInterval, setGageState), gageNextRainDate, gageNextRateChange (UEB snow, 6 November
// 2026) and updateAllGages.
// NOT YET TRANSLATED: convertRainfall, the GageState overloads of updatePastRain and
// getPastRain, convertGageValue and the report-rainfall functions.
//
// Differences of form, never in result:
//   * gageRainSeries returns a pointer into either ctx.gages.rain_series or ctx.tables.
//     Here gage_rain_series returns the reference for read-only callers, and
//     update_all_gages moves the series out of the context for the duration of one gage's
//     update (SeriesLoc, take_series/put_series) so the state machine can write the
//     gage arrays beside it; nothing in that update reads the series through the context.
// `double past_rain[MAXPASTRAIN]` is a fixed array; member initializers carry into a
// hand-written Default for GageState.
//
// 1D build: the 2D-mesh clause (every gage is used under a rain-on-grid RAINFALL_MODE)
// tests ctx.twod_io.mesh, which a 1D build never sets, so it is always false here.

use crate::core::cmath::{cclamp, cmin};
use crate::core::context::SimulationContext;
use crate::core::datetime;
use crate::core::unit_conversion as ucf;
use crate::data::forcing::ForcingMode;
use crate::data::gage::{RainFileFormat, RainSource};
use crate::data::table::{Table, TableType};
use crate::hydrology::ueb_types::SnowKind;

/// One second in days
pub const ONE_SECOND: f64 = 1.1574074e-5;
/// Max past hours tracked
pub const MAXPASTRAIN: usize = 48;

/// Rain type codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum RainType {
    Intensity = 0,
    Volume = 1,
    Cumulative = 2,
}

/// Per-gage state.
#[derive(Debug, Clone, Copy)]
pub struct GageState {
    /// Current intensity (project in/hr or mm/hr)
    pub rainfall: f64,
    /// Current snow intensity (project in/hr or mm/hr)
    pub snowfall: f64,
    /// rainfall + snowfall
    pub total_precip: f64,
    /// Cumulative rain accumulator (for CUMULATIVE type)
    pub rain_accum: f64,
    /// API-overridden rainfall (-1 = no override)
    pub api_rainfall: f64,
    /// Snow catch factor
    pub snow_factor: f64,
    /// Rainfall scaling factor (1.0 = no scaling)
    pub scale_factor: f64,
    /// Unit conversion factor
    pub units_factor: f64,
    /// Monthly/seasonal adjustment factor
    pub adjust_factor: f64,
    /// Rain recording interval (seconds)
    pub rain_interval: f64,
    pub rain_type: RainType,

    /// Past hourly rainfall totals
    pub past_rain: [f64; MAXPASTRAIN],
    /// Accumulator for current hour
    pub past_rain_accum: f64,
    /// Seconds in the current partial hour (0..3600)
    pub past_rain_interval: f64,
}

impl Default for GageState {
    fn default() -> Self {
        GageState {
            rainfall: 0.0,
            snowfall: 0.0,
            total_precip: 0.0,
            rain_accum: 0.0,
            api_rainfall: -1.0,
            snow_factor: 1.0,
            scale_factor: 1.0,
            units_factor: 1.0,
            adjust_factor: 1.0,
            rain_interval: 0.0,
            rain_type: RainType::Intensity,
            past_rain: [0.0; MAXPASTRAIN],
            past_rain_accum: 0.0,
            past_rain_interval: 0.0,
        }
    }
}

/// True when a subcatchment or an RDII unit hydrograph reads the gage.
pub fn gage_is_used(ctx: &SimulationContext, gage_idx: i32) -> bool {
    for i in 0..ctx.n_subcatches() {
        if ctx.subcatches.gage[i as usize] == gage_idx {
            return true;
        }
    }
    for gname in &ctx.unit_hyds.gage_names {
        if ctx.gage_names.find(gname) == gage_idx {
            return true;
        }
    }
    false
}

// ============================================================================
// Run-time functions (Gage.cpp), translated 5 November 2026 for the runoff step
// ============================================================================

/// PrecipSplit: one subcatchment's precipitation, ft/s.
#[derive(Debug, Clone, Copy, Default)]
pub struct PrecipSplit {
    pub rainfall: f64,
    pub snowfall: f64,
}

/// splitPrecip: the gage intensity split into rain and snow by the air temperature, with
/// the snow catch and the subcatchment scale factors (legacy gage.c:513-523,
/// subcatch.c:789-790).
pub fn split_precip(ctx: &SimulationContext, sub: usize) -> PrecipSplit {
    let mut out = PrecipSplit::default();

    let gi = ctx.subcatches.gage[sub];
    if gi < 0 || gi >= ctx.n_gages() {
        return out;
    }
    let ug = gi as usize;

    let intensity = ctx.gages.rainfall[ug];

    let is_snowing = !ctx.options.ignore_snow_melt && (ctx.climate_state.temperature <= ctx.climate_state.snow_divt);

    let ucf_rain = ucf::ucf(ucf::RAINFALL, &ctx.options);
    let pack = if sub < ctx.subcatches.snowpack.len() { ctx.subcatches.snowpack[sub] } else { -1 };
    let snow_config = &ctx.snowpacks.config;
    if !ctx.options.ignore_snow_melt && pack >= 0 && snow_config.model_for(pack) == SnowKind::Ueb {
        let p = &snow_config.parameters;
        let temperature = (ctx.climate_state.temperature - 32.) / 1.8;
        let fraction = if snow_config.division_temperature {
            if is_snowing {
                1.
            } else {
                0.
            }
        } else {
            cclamp((p.rain_temperature - temperature) / (p.rain_temperature - p.snow_temperature), 0., 1.)
        };
        let drift = match snow_config.packs.get(&pack) {
            None => 1.0,
            Some(site) => site.site.drift,
        };
        out.rainfall = intensity / ucf_rain * (1.0 - fraction) * ctx.subcatches.rain_scale_factor[sub];
        out.snowfall = intensity / ucf_rain * fraction * ctx.gages.snow_factor[ug] * drift * ctx.subcatches.snow_scale_factor[sub];
        return out;
    }
    if is_snowing {
        out.snowfall = intensity * ctx.gages.snow_factor[ug] / ucf_rain;
        out.snowfall *= ctx.subcatches.snow_scale_factor[sub];
    } else {
        out.rainfall = intensity / ucf_rain;
        out.rainfall *= ctx.subcatches.rain_scale_factor[sub];
    }
    out
}

/// accumulatePastRain (file-static): integrate a held intensity into the hourly history.
fn accumulate_past_rain(completed: &mut [f64], partial: &mut f64, interval: &mut f64, rainfall: f64, mut time_step: i32) {
    let rate = rainfall / 3600.0;
    while time_step > 0 {
        let remaining = 3600 - (*interval as i32);
        if time_step > remaining {
            *partial += remaining as f64 * rate;
            let mut h = MAXPASTRAIN - 1;
            while h > 0 {
                completed[h] = completed[h - 1];
                h -= 1;
            }
            completed[0] = *partial;
            *partial = 0.0;
            *interval = 0.0;
            time_step -= remaining;
        } else {
            *partial += time_step as f64 * rate;
            *interval += time_step as f64;
            time_step = 0;
        }
    }
}

/// updatePastRain(ctx, time_step): every gage's hourly history.
pub fn update_past_rain(ctx: &mut SimulationContext, time_step: i32) {
    let g = &mut ctx.gages;
    for j in 0..g.count() {
        let uj = j as usize;
        let (a, b) = (uj * MAXPASTRAIN, uj * MAXPASTRAIN + MAXPASTRAIN);
        accumulate_past_rain(&mut g.past_rain[a..b], &mut g.past_rain_accum[uj], &mut g.past_rain_interval[uj], g.rainfall[uj], time_step);
    }
}

/// Where a gage's rain series lives (gageRainSeries's pointer).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeriesLoc {
    None,
    File(usize),
    Ts(usize),
}

pub fn gage_rain_series_loc(ctx: &SimulationContext, gage_idx: i32) -> SeriesLoc {
    if gage_idx < 0 || gage_idx >= ctx.n_gages() {
        return SeriesLoc::None;
    }
    let ug = gage_idx as usize;
    if ctx.gages.source[ug] == RainSource::FileRain {
        if ug < ctx.gages.rain_series.len() && !ctx.gages.rain_series[ug].is_empty() {
            return SeriesLoc::File(ug);
        }
        return SeriesLoc::None;
    }
    let ts_idx = ctx.gages.ts_index[ug];
    if ts_idx >= 0 && (ts_idx as usize) < ctx.tables.tables.len() {
        return SeriesLoc::Ts(ts_idx as usize);
    }
    SeriesLoc::None
}

/// gageRainSeries: the gage's series, read only.
pub fn gage_rain_series(ctx: &SimulationContext, gage_idx: i32) -> Option<&Table> {
    match gage_rain_series_loc(ctx, gage_idx) {
        SeriesLoc::None => None,
        SeriesLoc::File(u) => Some(&ctx.gages.rain_series[u]),
        SeriesLoc::Ts(u) => Some(&ctx.tables.tables[u]),
    }
}

fn take_series(ctx: &mut SimulationContext, loc: SeriesLoc) -> Option<Table> {
    let slot = match loc {
        SeriesLoc::None => return None,
        SeriesLoc::File(u) => &mut ctx.gages.rain_series[u],
        SeriesLoc::Ts(u) => &mut ctx.tables.tables[u],
    };
    Some(std::mem::replace(slot, Table::new(b"", TableType::Timeseries)))
}

fn put_series(ctx: &mut SimulationContext, loc: SeriesLoc, tbl: Option<Table>) {
    if let Some(t) = tbl {
        match loc {
            SeriesLoc::None => {}
            SeriesLoc::File(u) => ctx.gages.rain_series[u] = t,
            SeriesLoc::Ts(u) => ctx.tables.tables[u] = t,
        }
    }
}

pub fn gage_series_rain_type(ctx: &SimulationContext, gage_idx: i32) -> i32 {
    let ug = gage_idx as usize;
    if ctx.gages.source[ug] == RainSource::FileRain && ctx.gages.file_format[ug] == RainFileFormat::StanPrcp {
        return 1;
    }
    ctx.gages.rain_type[ug]
}

pub fn gage_units_factor(ctx: &SimulationContext, gage_idx: i32) -> f64 {
    if gage_idx < 0 || gage_idx >= ctx.n_gages() {
        return 1.0;
    }
    let ug = gage_idx as usize;
    if ctx.gages.source[ug] == RainSource::FileRain && ctx.gages.file_format[ug] == RainFileFormat::StanPrcp && ctx.options.flow_units as i32 >= 3 {
        return 25.40;
    }
    1.0
}

const K_NO_STATE: i32 = -2;
const K_PRE_RECORD: i32 = -1;
const K_NO_DATE: i32 = -1;
const K_NO_DATE_VALUE: f64 = -693594.0;

/// recordRate: one record's rate (project units per hour), unit then scale factor.
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

fn next_nonzero_record(ctx: &SimulationContext, j: i32, tbl: &Table, from: i32, adjustment: f64) -> i32 {
    let n = tbl.x.len() as i32;
    let mut k = from + 1;
    while k < n {
        if record_rate(ctx, j, tbl, k) * adjustment != 0.0 {
            return k;
        }
        k += 1;
    }
    K_NO_DATE
}

fn init_gage_state(ctx: &mut SimulationContext, j: i32, tbl: Option<&Table>) {
    let uj = j as usize;
    ctx.gages.st_init[uj] = 1;
    ctx.gages.st_used[uj] = if gage_is_used(ctx, j) { 1 } else { 0 };
    ctx.gages.st_cur[uj] = K_NO_STATE;
    ctx.gages.st_next[uj] = K_NO_DATE;
    ctx.gages.st_rain[uj] = 0.0;
    ctx.gages.next_rainfall[uj] = 0.0;
    if ctx.options.ignore_rainfall {
        return;
    }
    let tbl = match tbl {
        Some(t) if !t.x.is_empty() => t,
        _ => return,
    };
    if tbl.x[0] > ctx.options.start_date {
        ctx.gages.st_cur[uj] = K_PRE_RECORD;
        ctx.gages.st_next[uj] = 0;
    } else {
        ctx.gages.st_cur[uj] = 0;
        ctx.gages.st_rain[uj] = record_rate(ctx, j, tbl, 0);
        ctx.gages.st_next[uj] = next_nonzero_record(ctx, j, tbl, 0, 1.0);
    }
    let next = ctx.gages.st_next[uj];
    if next >= 0 {
        ctx.gages.next_rainfall[uj] = record_rate(ctx, j, tbl, next);
    }
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

fn set_gage_state(ctx: &mut SimulationContext, j: i32, tbl: &Table, t: f64, adjustment: f64) {
    let uj = j as usize;
    loop {
        if ctx.gages.st_cur[uj] == K_NO_STATE {
            ctx.gages.st_rain[uj] = 0.0;
            return;
        }
        let (start, end) = current_interval(ctx, j, tbl);
        if t < start {
            ctx.gages.st_rain[uj] = 0.0;
            return;
        }
        if t < end {
            return;
        }
        let nxt = ctx.gages.st_next[uj];
        if nxt < 0 {
            ctx.gages.st_rain[uj] = 0.0;
            return;
        }
        if t < tbl.x[nxt as usize] {
            ctx.gages.st_rain[uj] = 0.0;
            return;
        }
        ctx.gages.st_cur[uj] = nxt;
        ctx.gages.st_rain[uj] = ctx.gages.next_rainfall[uj];
        let next = next_nonzero_record(ctx, j, tbl, nxt, adjustment);
        ctx.gages.st_next[uj] = next;
        ctx.gages.next_rainfall[uj] = if next >= 0 { record_rate(ctx, j, tbl, next) * adjustment } else { 0.0 };
    }
}

/// gageNextRainDate (legacy gage_getNextRainDate).
pub fn gage_next_rain_date(ctx: &SimulationContext, gage_idx: i32, t: f64) -> f64 {
    let ug = gage_idx as usize;
    if ug >= ctx.gages.st_cur.len() || ctx.gages.st_used[ug] == 0 {
        return t;
    }
    let cur = ctx.gages.st_cur[ug];
    let tbl = match gage_rain_series(ctx, gage_idx) {
        Some(t) if cur != K_NO_STATE && !t.x.is_empty() => t,
        _ => return K_NO_DATE_VALUE,
    };
    let t1 = t + datetime::ONE_SECOND;
    let (start, end) = current_interval(ctx, gage_idx, tbl);
    if t1 < start {
        return start;
    }
    if t1 < end {
        return end;
    }
    let nxt = ctx.gages.st_next[ug];
    if nxt >= 0 {
        tbl.x[nxt as usize]
    } else {
        K_NO_DATE_VALUE
    }
}

/// gageNextRateChange: the first date in (date, limit) at which gage j's rain rate changes,
/// or limit (the UEB snow step's horizon).
pub fn gage_next_rate_change(ctx: &SimulationContext, j: i32, date: f64, limit: f64) -> f64 {
    let table = match gage_rain_series(ctx, j) {
        Some(t) if !t.x.is_empty() && !ctx.options.ignore_rainfall => t,
        _ => return limit,
    };
    const TOLERANCE: f64 = 1e-5 / 86400.0;
    let uj = j as usize;
    let duration = ctx.gages.interval_sec[uj] as f64 / 86400.0;
    let mut current = ctx.gages.st_cur[uj];
    let mut end = date;
    let mut rate = 0.0;
    if current >= 0 {
        end = table.x[current as usize] + duration;
        if date < end - TOLERANCE {
            rate = record_rate(ctx, j, table, current);
        }
    }
    if rate == 0.0 {
        let mut next = ctx.gages.st_next[uj];
        while next >= 0 && table.x[next as usize] <= date + TOLERANCE {
            next = next_nonzero_record(ctx, j, table, next, 1.0);
        }
        return if next >= 0 { cmin(limit, table.x[next as usize]) } else { limit };
    }
    while end < limit - TOLERANCE {
        let next = next_nonzero_record(ctx, j, table, current, 1.0);
        if next < 0 || table.x[next as usize] > end + TOLERANCE || record_rate(ctx, j, table, next) != rate {
            return end;
        }
        current = next;
        end = table.x[current as usize] + duration;
    }
    limit
}

/// updateAllGages: every gage's rainfall at `current_time` (absolute days).
pub fn update_all_gages(ctx: &mut SimulationContext, current_time: f64) {
    for j in 0..ctx.n_gages() {
        let uj = j as usize;
        let loc = gage_rain_series_loc(ctx, j);
        let rtbl = take_series(ctx, loc);
        update_one_gage(ctx, j, uj, rtbl.as_ref(), current_time);
        put_series(ctx, loc, rtbl);
    }
    if !ctx.options.ignore_rainfall {
        for j in 0..ctx.n_gages() {
            let uj = j as usize;
            if ctx.gages.st_used[uj] == 0 || uj >= ctx.forcing.gage_rainfall_mode.len() {
                continue;
            }
            if ctx.forcing.gage_rainfall_mode[uj] == ForcingMode::Override {
                ctx.gages.rainfall[uj] = ctx.forcing.gage_rainfall_value[uj];
            } else if ctx.forcing.gage_rainfall_mode[uj] == ForcingMode::Add {
                ctx.gages.rainfall[uj] += ctx.forcing.gage_rainfall_value[uj];
            }
        }
    }
}

/// The body of updateAllGages's per-gage loop (a `continue` there is a return here).
fn update_one_gage(ctx: &mut SimulationContext, j: i32, uj: usize, rtbl: Option<&Table>, current_time: f64) {
    if uj < ctx.gages.st_init.len() && ctx.gages.st_init[uj] == 0 {
        init_gage_state(ctx, j, rtbl);
    }
    if uj < ctx.gages.st_used.len() && ctx.gages.st_used[uj] == 0 {
        ctx.gages.rainfall[uj] = ctx.gages.st_rain[uj];
        return;
    }
    if ctx.options.ignore_rainfall {
        ctx.gages.rainfall[uj] = 0.0;
        return;
    }
    let co = if uj < ctx.gages.co_gage_index.len() { ctx.gages.co_gage_index[uj] } else { -1 };
    if co >= 0 && co < j {
        let uco = co as usize;
        let primary_sf = ctx.gages.scale_factor[uco];
        let this_sf = ctx.gages.scale_factor[uj];
        let ratio = if primary_sf > 0.0 { this_sf / primary_sf } else { 1.0 };
        ctx.gages.rainfall[uj] = ctx.gages.rainfall[uco] * ratio;
        return;
    }
    if ctx.gages.api_rainfall[uj] >= 0.0 {
        ctx.gages.rainfall[uj] = ctx.gages.api_rainfall[uj];
        return;
    }
    let mut raw_value = 0.0;
    if let Some(tbl) = rtbl {
        let month = datetime::month_of_year(current_time) - 1;
        let adjustment = if (0..12).contains(&month) { ctx.adjust_rain[month as usize] } else { 1.0 };
        set_gage_state(ctx, j, tbl, current_time + datetime::ONE_SECOND, adjustment);
        raw_value = ctx.gages.st_rain[uj];
    }
    ctx.gages.rainfall[uj] = raw_value;
}
