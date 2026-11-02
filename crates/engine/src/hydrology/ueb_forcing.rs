// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/UebForcing.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Meteorological forcing preparation for the UEB snow models, the part
// SWMMEngine::initHydrology reaches: climate::initializeRadiation (validate the UEB
// parameters, apply the [TEMPERATURE] SITE, load the hourly MET_FILE, build the daily
// snow-weather cache from a climate file, check radiation time series). The error is a
// std::string upstream (empty = success); here Bytes, empty on success.
//
// NOT YET TRANSLATED (run-time step functions, phase 1c): resolveRadiation,
// resolveUebForcing, and the file-static helpers seriesAverage, dailyRange,
// meanDailyTemperature, saturatedVapourPressure.

use crate::core::context::SimulationContext;
use crate::core::options::Bytes;
use crate::core::path_resolver::{parent_dir, resolve_relative};
use crate::core::unit_conversion as ucf;
use crate::hydrology::climate_file::{ClimateFileReader, DailyClimateRecord};
use crate::hydrology::hourly_meteorology::{HourlyMeteorology, SnowDailyWeather};
use crate::hydrology::ueb_types::{RadiationSourceMode, SnowKind};

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

pub fn initialize_radiation(ctx: &mut SimulationContext) -> Bytes {
    {
        let error = ctx.snowpacks.config.parameters.validate();
        if !error.is_empty() {
            return error;
        }
    }
    if let Some(site) = ctx.snowpacks.config.site {
        let solar = &ctx.heat_config.solar;
        let state = &ctx.climate_state;
        let conflict = (solar.has_latitude && solar.latitude_deg != site.latitude)
            || (solar.has_longitude && solar.longitude_deg != site.longitude)
            || (solar.has_timezone && solar.timezone_hours != site.timezone)
            || (solar.has_elevation && solar.elevation_m != site.elevation)
            || (ctx.options.snow_lat != 0.0 && ctx.options.snow_lat != site.latitude)
            || (ctx.options.snow_elev != 0.0 && state.elev * 0.3048 != site.elevation);
        let warning: &[u8] = b"[TEMPERATURE] SITE overrides conflicting legacy snow/solar site values.";
        if conflict && !ctx.warnings.iter().any(|w| w.as_slice() == warning) {
            ctx.warnings.push(warning.to_vec());
        }
        let state = &mut ctx.climate_state;
        state.latitude = site.latitude;
        state.elev = site.elevation / 0.3048;
        state.dtlong = site.timezone - site.longitude / 15.0;
        let z = state.elev / 1000.0;
        let pressure_inhg = if z <= 0.0 { 29.9 } else { 29.9 - 1.02 * z + 0.0032 * z.powf(2.4) }; // LIBM
        state.gamma = 0.000359 * pressure_inhg;
    }
    if !ctx.snowpacks.config.met_file.is_empty() {
        ctx.snowpacks.config.met_file.absolute = resolve_relative(&ctx.snowpacks.config.met_file.original, &parent_dir(&ctx.inp_file_path));
        let path = ctx.snowpacks.config.met_file.absolute.clone();
        let error = ctx.snowpacks.meteorology.load(&path);
        if !error.is_empty() {
            return error;
        }
        if ctx.options.end_date > ctx.options.start_date {
            let mut values = [0.0f64; 5];
            let error = ctx.snowpacks.meteorology.average(ctx.options.start_date, (ctx.options.end_date - ctx.options.start_date) * 86400.0, &mut values);
            if !error.is_empty() {
                return error;
            }
        }
    } else {
        ctx.snowpacks.meteorology = HourlyMeteorology::default();
    }
    ctx.snowpacks.daily_weather.clear();
    let cfg = &ctx.snowpacks.config;
    let needed = cfg.model == SnowKind::Ueb || !cfg.packs.is_empty() || ctx.snowpacks.mesh_config.assigned() || cfg.shortwave.authored || cfg.longwave.authored;
    if needed && ctx.snowpacks.meteorology.is_empty() && !ctx.options.temp_file.is_empty() && (ctx.options.temp_source == 2 || ctx.options.wind_type == 1) {
        let mut reader = ClimateFileReader::default();
        let path = resolve_relative(&ctx.options.temp_file.original, &parent_dir(&ctx.inp_file_path));
        if !reader.open(&path, ctx.options.temp_file_start, ucf::get_unit_system(ctx.options.flow_units as i32), ctx.options.temp_units) {
            return b"Snow weather could not read the configured daily climate file.".to_vec();
        }
        reader.set_day_offset(if ctx.options.temp_file_start > 0.0 { ctx.options.temp_file_start.floor() - ctx.options.start_date.floor() } else { 0.0 });
        let (mut low, mut high, mut wind, mut previous) = (70.0, 70.0, 0.0, 70.0);
        let mut day = ctx.options.start_date.floor();
        while day < ctx.options.end_date {
            let mut rec = DailyClimateRecord::default();
            if reader.get_record(day, &mut rec) {
                if rec.tmin.is_finite() {
                    low = rec.tmin;
                }
                if rec.tmax.is_finite() {
                    high = rec.tmax;
                }
                if rec.wind.is_finite() {
                    wind = rec.wind;
                }
            }
            if ctx.snowpacks.daily_weather.is_empty() {
                previous = high;
            }
            ctx.snowpacks.daily_weather.push(SnowDailyWeather { date: day, tmin: low, tmax: high, previous_max: previous, wind });
            previous = high;
            day += 1.0;
        }
    }
    let sources = [ctx.snowpacks.config.shortwave.clone(), ctx.snowpacks.config.longwave.clone()];
    for source in sources.iter() {
        if source.mode == RadiationSourceMode::File && ctx.snowpacks.meteorology.is_empty() {
            return b"Radiation FILE source requires MET_FILE.".to_vec();
        }
        if source.mode != RadiationSourceMode::Timeseries || !ctx.snowpacks.meteorology.is_empty() {
            continue;
        }
        let index = ctx.find_timeseries(&source.series);
        if index < 0 {
            return cat(&[b"Radiation time series not found: ", &source.series]);
        }
        let table = &ctx.tables.tables[index as usize];
        if table.x.is_empty() || table.x[0] > ctx.options.start_date || *table.x.last().unwrap() < ctx.options.end_date {
            return cat(&[b"Radiation time series does not cover the simulation: ", &source.series]);
        }
        for &value in &table.y {
            if !value.is_finite() || value < 0.0 {
                return cat(&["Radiation time series must contain finite nonnegative W/m\u{b2}: ".as_bytes(), &source.series]);
            }
        }
    }
    Bytes::new()
}
