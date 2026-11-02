// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/HourlyMeteorology.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The hourly meteorology record, the daily snow-weather cache row and the hourly cache
// class of namespace openswmm::climate. The class keeps its one private member
// (`records_`, here `records`) private. Its two inline accessors are here; `empty()` is
// `is_empty()`, as elsewhere in this crate.
//
// Also derived from src/engine/hydrology/HourlyMeteorology.cpp (same repository, branch
// and commit, Copyright (c) 2026 Caleb Buahin, Apache License 2.0), translated 2 November
// 2026: load, read and average. `std::istream` is core::cstdio::IfStream; read() takes
// one; load() opens it. Messages are byte strings, empty for success.

#![allow(clippy::all)]

/// Date/time marks the START of the hourly interval. Temperatures are °C, humidity is
/// percent, wind m/s and radiation W/m². Every initializer is zero (`values{}`), so
/// Default is derived.
#[derive(Debug, Clone, Copy, Default)]
pub struct HourlyMetRecord {
    pub date: f64,
    /// temperature, RH, wind, Qsi, Qli
    pub values: [f64; 5],
}

/// Immutable daily source cache for independent surface clocks.
#[derive(Debug, Clone, Copy)]
pub struct SnowDailyWeather {
    pub date: f64,
    pub tmin: f64,
    pub tmax: f64,
    pub previous_max: f64,
    pub wind: f64,
}

impl Default for SnowDailyWeather {
    fn default() -> Self {
        SnowDailyWeather {
            date: 0.0,
            tmin: 70.0,
            tmax: 70.0,
            previous_max: 70.0,
            wind: 0.0,
        }
    }
}

/// Owned consecutive hourly records with interval-weighted means.
#[derive(Debug, Clone, Default)]
pub struct HourlyMeteorology {
    records: Vec<HourlyMetRecord>,
}

impl HourlyMeteorology {
    /// Open `path` and read() it; "Cannot open ..." when it does not open.
    pub fn load(&mut self, path: &[u8]) -> Vec<u8> {
        match crate::core::cstdio::IfStream::open(path) {
            None => {
                let mut m = b"Cannot open hourly meteorology file: ".to_vec();
                m.extend_from_slice(path);
                m
            }
            Some(mut stream) => self.read(&mut stream),
        }
    }

    /// Replace the records with the stream's consecutive hourly rows, or keep them and
    /// return why not.
    pub fn read(&mut self, stream: &mut crate::core::cstdio::IfStream) -> Vec<u8> {
        use crate::data::table::{sscanf, ScanArg};
        let mut candidate: Vec<HourlyMetRecord> = Vec::new();
        let mut line: Vec<u8> = Vec::new();
        let mut line_number: i32 = 0;
        while stream.getline(&mut line) {
            line_number += 1;
            let tok = crate::input::tokenizer::tokenize(&line);
            if tok.is_empty() {
                continue;
            }
            let error = || -> Vec<u8> {
                format!(
                    "Invalid hourly meteorology at line {}: expected MM/DD/YYYY HH:00 temperature_C RH_percent wind_m/s shortwave_W/m2 longwave_W/m2.",
                    line_number
                )
                .into_bytes()
            };
            if tok.len() != 7 {
                return error();
            }
            let (mut m, mut d, mut y, mut h, mut minute, mut end) = (0i32, 0i32, 0i32, 0i32, 0i32, 0i32);
            if sscanf(&tok[0], b"%d/%d/%d%n", &mut [ScanArg::Int(&mut m), ScanArg::Int(&mut d), ScanArg::Int(&mut y), ScanArg::Count(&mut end)]) != 3
                || end != tok[0].len() as i32
            {
                return error();
            }
            end = 0;
            if sscanf(&tok[1], b"%d:%d%n", &mut [ScanArg::Int(&mut h), ScanArg::Int(&mut minute), ScanArg::Count(&mut end)]) != 2
                || end != tok[1].len() as i32
                || h < 0
                || h > 23
                || minute != 0
            {
                return error();
            }
            let day = crate::core::datetime::encode_date(y, m, d);
            if day <= 0.0 {
                return error();
            }
            let mut record = HourlyMetRecord::default();
            record.date = day + h as f64 / 24.0;
            for i in 0..5 {
                if !crate::input::parse_utils::parse_double_strict(&tok[i + 2], &mut record.values[i]) || !record.values[i].is_finite() {
                    return error();
                }
            }
            let v = &record.values;
            if v[0] <= -100.0 || v[0] > 70.0 || v[1] < 0.0 || v[1] > 100.0 || v[2] < 0.0 || v[3] < 0.0 || v[4] < 0.0 {
                return error();
            }
            if !candidate.is_empty() && ((record.date - candidate.last().unwrap().date) * 86400.0 - 3600.0).abs() > 1e-5 {
                return format!(
                    "Hourly meteorology must have consecutive, unique hourly records; gap or duplicate at line {}.",
                    line_number
                )
                .into_bytes();
            }
            candidate.push(record);
        }
        if stream.bad() {
            return b"Could not read hourly meteorology.".to_vec();
        }
        if candidate.is_empty() {
            return b"Hourly meteorology contains no records.".to_vec();
        }
        self.records = candidate;
        Vec::new()
    }

    /// Interval-weighted mean of the five values over [start, start + seconds), or why
    /// not (result untouched).
    pub fn average(&self, start: f64, seconds: f64, result: &mut [f64; 5]) -> Vec<u8> {
        use crate::core::cmath::{cmax, cmin};
        if !start.is_finite() || !seconds.is_finite() || seconds <= 0.0 || self.records.is_empty() {
            return b"Invalid hourly meteorology request.".to_vec();
        }
        const TOLERANCE: f64 = 1e-5;
        let mut offset = (start - self.records[0].date) * 86400.0;
        let nearest_hour = (offset / 3600.0).round() * 3600.0;
        if (offset - nearest_hour).abs() <= TOLERANCE {
            offset = nearest_hour;
        }
        let available = self.records.len() as f64 * 3600.0;
        if !offset.is_finite() || offset < -TOLERANCE || offset >= available || seconds > available - offset + TOLERANCE {
            return b"Hourly meteorology does not cover the requested simulation interval.".to_vec();
        }
        offset = cmax(0.0, offset);
        let mut index = (offset / 3600.0) as usize;
        let mut within_hour = offset - index as f64 * 3600.0;
        let mut remaining = seconds;
        let mut sum = [0.0f64; 5];
        let mut elapsed = 0.0;
        while index < self.records.len() && remaining > 0.0 {
            let duration = cmin(remaining, 3600.0 - within_hour);
            for i in 0..5 {
                sum[i] += duration * self.records[index].values[i];
            }
            elapsed += duration;
            remaining -= duration;
            within_hour = 0.0;
            index += 1;
        }
        if !(elapsed > 0.0) || remaining > TOLERANCE {
            return b"Hourly meteorology interval has incomplete coverage.".to_vec();
        }
        for v in sum.iter_mut() {
            *v /= elapsed;
        }
        *result = sum;
        Vec::new()
    }

    /// True when no validated records are stored (upstream empty()).
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    /// Borrow the immutable hourly starts and their complete-hour weather values.
    #[inline]
    pub fn records(&self) -> &Vec<HourlyMetRecord> {
        &self.records
    }
}

/// Gate only (cargo feature `gate`): replace the records, for the context loader in
/// examples/common/ctxload.rs. Not part of the translation.
#[cfg(feature = "gate")]
impl HourlyMeteorology {
    pub fn gate_set_records(&mut self, records: Vec<HourlyMetRecord>) {
        self.records = records;
    }
}
