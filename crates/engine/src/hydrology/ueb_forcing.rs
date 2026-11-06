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
// The run-time step functions (6 November 2026): resolveRadiation, resolveUebForcing, and
// the file-static helpers seriesAverage, dailyRange, meanDailyTemperature,
// saturatedVapourPressure. surface::solarIntegral throws std::invalid_argument; inside
// resolveRadiation nothing catches it, so resolve_radiation returns Result<Bytes, String>
// with Err(what) for that exception on its way up (resolveUebForcing lets it pass, and
// catches only its own solarIntegral call, as upstream).

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

// ---------------------------------------------------------------------------
// The run-time step (6 November 2026): resolveRadiation, resolveUebForcing and the
// file-static helpers seriesAverage, dailyRange, meanDailyTemperature, saturatedVapourPressure.
// ---------------------------------------------------------------------------

/// seriesAverage: the trapezoidal mean of a series over [date, date + seconds].
fn series_average(ctx: &mut SimulationContext, index: i32, date: f64, seconds: f64, output: &mut f64, coverage: bool) -> Bytes {
    if index < 0 || index >= ctx.tables.count() as i32 {
        return b"Radiation time series was not found.".to_vec();
    }
    let table = &mut ctx.tables[index];
    let end = date + seconds / 86400.0;
    if table.x.is_empty() || table.y.len() != table.x.len() || (coverage && (date < table.x[0] || end > *table.x.last().unwrap())) {
        return cat(&[b"Time series ", &table.id, b" does not cover the forcing interval."]);
    }
    let mut a = date;
    let mut integral = 0.0;
    // std::upper_bound(x, a)
    let mut it = table.x.partition_point(|&v| v <= a);
    while a < end {
        let b = if it == table.x.len() { end } else { crate::core::cmath::cmin(end, table.x[it]) };
        let va = crate::data::table::table_lookup_cursor(table, a);
        let vb = crate::data::table::table_lookup_cursor(table, b);
        if !va.is_finite() || !vb.is_finite() {
            return cat(&[b"Nonfinite forcing in time series ", &table.id, b"."]);
        }
        integral += 0.5 * (va + vb) * (b - a);
        a = b;
        if it != table.x.len() {
            it += 1;
        }
    }
    *output = integral / (end - date);
    Bytes::new()
}

/// dailyRange: the day's temperature range (C), or -1.
fn daily_range(ctx: &mut SimulationContext, date: f64) -> f64 {
    let mut low = f64::INFINITY;
    let mut high = -low;
    let day = date.floor();
    let accumulate = |t: f64, low: &mut f64, high: &mut f64| {
        *low = crate::core::cmath::cmin(*low, t);
        *high = crate::core::cmath::cmax(*high, t);
    };
    if !ctx.snowpacks.meteorology.is_empty() {
        let records = ctx.snowpacks.meteorology.records();
        let mut i = records.partition_point(|r| r.date < day);
        while i < records.len() && records[i].date < day + 1.0 {
            accumulate(records[i].values[0], &mut low, &mut high);
            i += 1;
        }
    } else if !ctx.snowpacks.daily_weather.is_empty() {
        let records = &ctx.snowpacks.daily_weather;
        let i = records.partition_point(|r| r.date < day);
        if i != records.len() && records[i].date == day {
            return crate::core::cmath::cmax(0.0, (records[i].tmax - records[i].tmin) * 5.0 / 9.0);
        }
    } else if ctx.climate_state.has_minmax {
        return crate::core::cmath::cmax(0.0, (ctx.climate_state.tmax_daily - ctx.climate_state.tmin_daily) * 5.0 / 9.0);
    } else if ctx.climate_state.temp_ts_index >= 0 {
        let scale = if ucf::get_unit_system(ctx.options.flow_units as i32) == 0 { 5.0 / 9.0 } else { 1.0 };
        let ti = ctx.climate_state.temp_ts_index;
        let table = &mut ctx.tables[ti];
        if !table.x.is_empty() {
            let v = crate::data::table::table_lookup_cursor(table, day) * scale;
            accumulate(v, &mut low, &mut high);
            let v = crate::data::table::table_lookup_cursor(table, day + 1.0) * scale;
            accumulate(v, &mut low, &mut high);
            let start = table.x.partition_point(|&v| v < day);
            let mut k = start;
            while k < table.x.len() && table.x[k] < day + 1.0 {
                accumulate(table.y[k] * scale, &mut low, &mut high);
                k += 1;
            }
        }
    }
    if low.is_finite() && high.is_finite() {
        high - low
    } else {
        -1.0
    }
}

/// meanDailyTemperature: the mean of the legacy sub-daily temperature curve over hours [a, b].
fn mean_daily_temperature(s: &crate::hydrology::climate_types::ClimateState, a: f64, b: f64) -> f64 {
    let pi = 3.14159265358979323846;
    let mut sum = 0.0;
    let mut integrate = |lo: f64, hi: f64, base: f64, amplitude: f64, phase: f64, frequency: f64| {
        let lo = crate::core::cmath::cmax(a, lo);
        let hi = crate::core::cmath::cmin(b, hi);
        if hi > lo {
            sum += base * (hi - lo) + amplitude * ((phase + frequency * lo).cos() - (phase + frequency * hi).cos()) / frequency; // LIBM
        }
    };
    let night = pi / (if s.dydif.abs() > 1e-10 { s.dydif } else { 1.0 });
    let daylight = pi / (if s.dhrdy.abs() > 1e-10 { s.dhrdy } else { -1.0 });
    integrate(0.0, s.hrsr, s.tmin_daily, 0.5 * (s.prev_tmax - s.tmin_daily), night * s.hrsr, -night);
    integrate(s.hrsr, s.hrss, 0.5 * (s.tmin_daily + s.tmax_daily), 0.5 * (s.tmax_daily - s.tmin_daily), daylight * s.hrday, -daylight);
    integrate(s.hrss, 24.0, s.tmax_daily, -0.5 * (s.tmax_daily - s.tmin_daily), -night * s.hrss, night);
    sum / (b - a)
}

fn saturated_vapour_pressure(t: f64) -> f64 {
    100.0 * (6.107799961 + t * (0.4436518521 + t * (0.01428945805 + t * (0.0002650648471 + t * (0.000003031240396 + t * (0.00000002034080948 + t * 0.00000000006136820929))))))
}

/// climate::resolveRadiation: the snow model's meteorology and radiation over one interval.
pub fn resolve_radiation(ctx: &mut SimulationContext, date: f64, seconds: f64, output: &mut crate::hydrology::radiation::RadiationState) -> Result<Bytes, String> {
    use crate::data::heat::ShortwaveMode;
    use crate::data::forcing::ForcingMode;
    use crate::hydrology::radiation::{RadiationOrigin, RadiationState};
    use crate::transport::heat_fluxes::{bird_clear_sky_ghi, cloud_shortwave_factor, pressure_from_elevation, solar_position, K_SOLAR_CONSTANT};
    if !date.is_finite() || !seconds.is_finite() || seconds <= 0.0 || seconds > 86400.0 {
        return Ok(b"Invalid radiation interval.".to_vec());
    }
    let mut next = RadiationState::default();
    let mut longitude;
    let mut timezone;
    let mut elevation;
    {
        let cfg = &ctx.snowpacks.config;
        let climate = &ctx.climate_state;
        let solar = &ctx.heat_config.solar;
        next.date = date;
        next.seconds = seconds;
        next.temperature = (climate.temperature - 32.0) * 5.0 / 9.0;
        next.humidity = climate.humidity * 0.01;
        next.wind = climate.wind_speed * 0.44704;
        longitude = if solar.has_longitude { solar.longitude_deg } else { -climate.dtlong * 15.0 };
        timezone = if solar.has_longitude { solar.timezone_hours } else { 0.0 };
        elevation = if solar.has_elevation { solar.elevation_m } else { climate.elev * 0.3048 };
        next.latitude = if solar.has_latitude { solar.latitude_deg } else { climate.latitude };
        if let Some(site) = cfg.site.as_ref() {
            next.latitude = site.latitude;
            longitude = site.longitude;
            timezone = site.timezone;
            elevation = site.elevation;
        }
    }
    next.pressure = pressure_from_elevation(elevation) * 100.0;
    let pos = solar_position(crate::core::datetime::day_of_year(date), (date - date.floor()) * 24.0, next.latitude, longitude, timezone);
    next.declination = pos.declination_deg;
    next.hour_angle = pos.hour_angle_deg;
    let integral = match crate::hydrology::solar_geometry::solar_integral(next.latitude, next.declination, next.hour_angle, seconds / 3600.0, 0.0, 0.0) {
        Ok(v) => v,
        Err(w) => return Err(w.to_string()),
    };
    next.cosine = integral.horizontal;
    next.horizontal = integral.horizontal;
    let (_year, month, _day) = crate::core::datetime::decode_date(date);
    let mut range = daily_range(ctx, date);
    let monthly = ctx.snowpacks.config.monthly_range[crate::core::cmath::cclamp(month - 1, 0, 11) as usize];
    if range < 0.0 {
        range = monthly.unwrap_or(0.0);
    }
    let mean_range = monthly.unwrap_or(range);
    let (bc_a, bc_c) = (ctx.snowpacks.config.bc_a, ctx.snowpacks.config.bc_c);
    let b = 0.036 * (-0.154 * mean_range).exp(); // LIBM
    next.transmittance = bc_a * (1.0 - (-b * range.powf(bc_c)).exp()); // LIBM
    next.cloud = 1.0 - next.transmittance / bc_a;
    let cloud = ctx.heat_config.cloud;
    if cloud.configured {
        next.cloud = cloud.fraction;
        if cloud.use_timeseries {
            let mut v = next.cloud;
            let error = series_average(ctx, cloud.ts_index, date, seconds, &mut v, true);
            next.cloud = v;
            if !error.is_empty() {
                return Ok(error);
            }
        }
        if !next.cloud.is_finite() || next.cloud < 0.0 || next.cloud > 1.0 {
            return Ok(b"Cloud forcing must be a fraction in [0,1].".to_vec());
        }
    }
    next.shortwave = next.transmittance * K_SOLAR_CONSTANT * pos.eccentricity * integral.horizontal;
    if !ctx.snowpacks.config.shortwave.authored && ctx.snowpacks.meteorology.is_empty() && ctx.options.heat_transport && ctx.heat_config.radiative_exchange {
        let legacy = ctx.heat_config.radiative.clone();
        if legacy.sw_mode == ShortwaveMode::Timeseries {
            let mut v = next.shortwave;
            let error = series_average(ctx, legacy.sw_ts_index, date, seconds, &mut v, true);
            next.shortwave = v;
            if !error.is_empty() {
                return Ok(error);
            }
        } else if legacy.sw_mode == ShortwaveMode::Constant {
            next.shortwave = legacy.shortwave_wm2;
        } else {
            next.shortwave = bird_clear_sky_ghi(&pos, next.pressure / 100.0, &ctx.heat_config.solar);
        }
        next.shortwave *= cloud_shortwave_factor(if cloud.configured { next.cloud } else { 0.0 }, cloud.sw_atten_k, cloud.sw_atten_n);
        next.shortwave_source = RadiationOrigin::LegacyHeat;
    }
    for channel in 0..2 {
        let source = if channel == 1 { ctx.snowpacks.config.longwave.clone() } else { ctx.snowpacks.config.shortwave.clone() };
        if source.mode == RadiationSourceMode::File && ctx.snowpacks.meteorology.is_empty() {
            return Ok(b"Radiation FILE source requires MET_FILE.".to_vec());
        }
        if source.mode != RadiationSourceMode::Timeseries || !ctx.snowpacks.meteorology.is_empty() {
            continue;
        }
        let mut value = 0.0;
        let idx = ctx.find_timeseries(&source.series);
        let error = series_average(ctx, idx, date, seconds, &mut value, true);
        if !error.is_empty() {
            return Ok(error);
        }
        if channel == 1 {
            next.longwave = value;
            next.longwave_source = RadiationOrigin::Timeseries;
        } else {
            next.shortwave = value;
            next.shortwave_source = RadiationOrigin::Timeseries;
        }
    }
    if !ctx.snowpacks.meteorology.is_empty() {
        let mut values = [0.0f64; 5];
        let error = ctx.snowpacks.meteorology.average(date, seconds, &mut values);
        if !error.is_empty() {
            return Ok(error);
        }
        next.temperature = values[0];
        next.humidity = values[1] * 0.01;
        next.wind = values[2];
        next.shortwave = values[3];
        next.longwave = values[4];
        next.shortwave_source = RadiationOrigin::MetFile;
        next.longwave_source = RadiationOrigin::MetFile;
    }
    if ctx.snowpacks.meteorology.is_empty() {
        let mon = (crate::core::datetime::month_of_year(date) - 1) as usize;
        let si = ucf::get_unit_system(ctx.options.flow_units as i32) == 1;
        let mut temperature = ctx.climate_state.temperature_src;
        let mut wind = ctx.climate_state.wind_speed_src;
        if ctx.climate_state.temp_ts_index >= 0 {
            let ti = ctx.climate_state.temp_ts_index;
            let error = series_average(ctx, ti, date, seconds, &mut temperature, false);
            if !error.is_empty() {
                return Ok(error);
            }
            if si {
                temperature = temperature * 1.8 + 32.0;
            }
            temperature += ctx.climate_state.adjust_temp[mon];
        } else if ctx.options.temp_source == 2 && !ctx.snowpacks.daily_weather.is_empty() {
            let mut sum = 0.0;
            let mut span = 0.0;
            let end = date + seconds / 86400.0;
            let mut a = date;
            while a < end {
                let day = a.floor();
                let b = crate::core::cmath::cmin(end, day + 1.0);
                let records = &ctx.snowpacks.daily_weather;
                let ri = records.partition_point(|r| r.date < day);
                if ri == records.len() || records[ri].date != day {
                    return Ok(b"Daily climate cache does not cover the snow interval.".to_vec());
                }
                let row = &records[ri];
                let mut daily = ctx.climate_state;
                let month = (crate::core::datetime::month_of_year(a) - 1) as usize;
                daily.tmin_daily = row.tmin + ctx.climate_state.adjust_temp[month];
                daily.tmax_daily = row.tmax + ctx.climate_state.adjust_temp[month];
                daily.prev_tmax = row.previous_max + ctx.climate_state.adjust_temp[month];
                crate::hydrology::climate::update_temp_times(&mut daily, crate::core::datetime::day_of_year(a));
                sum += mean_daily_temperature(&daily, (a - day) * 24.0, (b - day) * 24.0) * (b - a);
                span += b - a;
                a = b;
            }
            temperature = sum / span;
        } else if ctx.forcing.climate_temperature_mode == ForcingMode::None {
            temperature = ctx.climate_state.temperature;
        }
        if ctx.options.wind_type == 0 {
            wind = ctx.options.wind_speed[mon] / ucf::ucf(ucf::WINDSPEED, &ctx.options);
        } else if !ctx.snowpacks.daily_weather.is_empty() {
            let rows = &ctx.snowpacks.daily_weather;
            let ri = rows.partition_point(|r| r.date < date.floor());
            if ri != rows.len() {
                wind = rows[ri].wind;
            }
        }
        next.temperature = (ctx.forcing.effective_temperature(temperature) - 32.0) * 5.0 / 9.0;
        next.wind = ctx.forcing.effective_wind(wind) * 0.44704;
        let mut humidity = ctx.options.humidity[mon];
        if ctx.climate_state.humidity_ts_index >= 0 {
            let hi = ctx.climate_state.humidity_ts_index;
            let error = series_average(ctx, hi, date, seconds, &mut humidity, false);
            if !error.is_empty() {
                return Ok(error);
            }
        }
        if ctx.options.humidity_var == 1 {
            let dew = if si { humidity } else { (humidity - 32.0) * 5.0 / 9.0 };
            humidity = 100.0 * saturated_vapour_pressure(dew) / saturated_vapour_pressure(next.temperature);
        }
        next.humidity = crate::core::cmath::cclamp(humidity * 0.01, 0.0, 1.0);
    } else {
        next.temperature = (ctx.forcing.effective_temperature(next.temperature * 1.8 + 32.0) - 32.0) * 5.0 / 9.0;
        next.wind = ctx.forcing.effective_wind(next.wind / 0.44704) * 0.44704;
    }
    if next.longwave_source == RadiationOrigin::Computed {
        let kelvin = next.temperature + 273.15;
        let vapour = crate::core::cmath::cmax(0.0, saturated_vapour_pressure(next.temperature) * next.humidity);
        let emissivity = 1.08 * (1.0 - (-(vapour / 100.0).powf(kelvin / 2016.0)).exp()); // LIBM
        next.longwave = (next.cloud + (1.0 - next.cloud) * emissivity) * (2.0747e-7 / 3.6) * kelvin.powf(4.0); // LIBM
    }
    if !next.shortwave.is_finite() || next.shortwave < 0.0 || !next.longwave.is_finite() || next.longwave < 0.0 || !next.temperature.is_finite() || next.temperature <= -100.0 || next.temperature > 70.0 || !next.humidity.is_finite() || next.humidity < 0.0 || next.humidity > 1.0 || !next.wind.is_finite() || next.wind < 0.0 {
        return Ok(b"Invalid resolved UEB meteorological forcing.".to_vec());
    }
    *output = next;
    Ok(Bytes::new())
}

/// climate::resolveUebForcing: the point forcing for one UEB interval (and the site pressure).
pub fn resolve_ueb_forcing(ctx: &mut SimulationContext, date: f64, seconds: f64, slope: f64, aspect: f64, forcing: &mut crate::hydrology::ueb_types::Forcing, pressure: &mut f64) -> Result<Bytes, String> {
    let mut radiation = crate::hydrology::radiation::RadiationState::default();
    let error = resolve_radiation(ctx, date, seconds, &mut radiation)?;
    if !error.is_empty() {
        return Ok(error);
    }
    let geometry = match crate::hydrology::solar_geometry::solar_integral(radiation.latitude, radiation.declination, radiation.hour_angle, seconds / 3600.0, slope, aspect) {
        Ok(v) => v,
        Err(w) => return Ok(w.as_bytes().to_vec()),
    };
    let mut next = *forcing;
    next.seconds = seconds;
    next.temperature = radiation.temperature;
    next.humidity = radiation.humidity;
    next.wind = radiation.wind;
    next.longwave = radiation.longwave;
    next.cosine = geometry.cosine;
    if ctx.snowpacks.config.division_temperature {
        next.snow_fraction = if next.temperature <= (ctx.climate_state.snow_divt - 32.0) * 5.0 / 9.0 { 1.0 } else { 0.0 };
    }
    next.shortwave = if slope == 0.0 || geometry.horizontal <= 0.0 { radiation.shortwave } else { radiation.shortwave * geometry.hri / geometry.horizontal };
    *forcing = next;
    *pressure = radiation.pressure;
    ctx.climate_state.radiation = radiation;
    Ok(Bytes::new())
}
