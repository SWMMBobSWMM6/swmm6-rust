// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/GageData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson (first pass by tools/soa2rs.py, finished by hand). Not an
// OpenSWMM release.
//
// The rain gage store, one Vec per attribute (structure of arrays), sized and cleared
// exactly as upstream sizes and clears each one: resize() and grow_to() set every field,
// erase_at() and shrink_to_fit() only the ones upstream names (erase_at leaves
// station_id, rain_units, file_first_date, file_last_date, file_periods_precip and
// rain_series alone, as upstream does). `past_rain` is flat [gage * MAXPASTRAIN + hour].
//
// Differences in form only: upstream's `Table{}` is `default_table()` below (an empty
// TIMESERIES table, value-initialized `type` being 0), and the Table rows are filled with
// resize_with because Table holds a file handle and is not Clone. For the same reason
// GageData is not Clone.

#![allow(clippy::all)]
use crate::core::options::{Bytes, FilePathPair};
use crate::data::table::{Table, TableType};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i8)]
pub enum RainSource {
    Timeseries = 0,
    FileRain = 1,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i8)]
pub enum RainFileFormat {
    Unknown = -1,
    Nws15 = 0,
    NwsHourly = 1,
    Dsi3240 = 2,
    Dsi3260 = 3,
    HlyPrcp = 4,
    StanPrcp = 5,
    UserCsv = 6,
}

/// Upstream's `Table{}`: every member at its initializer, `type` value-initialized to 0
/// (TIMESERIES), `id` empty.
fn default_table() -> Table {
    Table::new(b"", TableType::Timeseries)
}

#[derive(Debug, Default)]
pub struct GageData {
    pub rain_type: Vec<i32>,
    pub source: Vec<RainSource>,
    pub ts_index: Vec<i32>,
    pub file_path: Vec<FilePathPair>,
    pub ts_name: Vec<Bytes>,
    pub col_name: Vec<Bytes>,
    pub station_id: Vec<Bytes>,
    pub rain_units: Vec<i32>,
    pub file_first_date: Vec<f64>,
    pub file_last_date: Vec<f64>,
    pub file_periods_precip: Vec<i64>,
    pub file_start_date: Vec<f64>,
    pub rain_series: Vec<Table>,
    pub file_format: Vec<RainFileFormat>,
    pub interval_sec: Vec<i32>,
    pub snow_factor: Vec<f64>,
    pub scale_factor: Vec<f64>,
    pub rainfall: Vec<f64>,
    pub next_rainfall: Vec<f64>,
    pub api_rainfall: Vec<f64>,
    pub next_rain_date: Vec<f64>,
    pub is_raining: Vec<bool>,
    pub st_cur: Vec<i32>,
    pub st_next: Vec<i32>,
    pub st_rain: Vec<f64>,
    pub st_init: Vec<u8>,
    pub st_used: Vec<u8>,
    pub past_rain: Vec<f64>,
    pub past_rain_accum: Vec<f64>,
    pub past_rain_interval: Vec<f64>,
    pub cumul_rain_accum: Vec<f64>,
    pub comments: Vec<Bytes>,
    pub co_gage_index: Vec<i32>,
}

// Upstream has no member initializers here (MAXPASTRAIN is a static constant), so
// Default is derived.

impl GageData {
    pub const MAXPASTRAIN: i32 = 48;

    pub fn count(&self) -> i32 {
        self.source.len() as i32
    }

    pub fn resize(&mut self, n: i32) {
        let un = usize::try_from(n).unwrap_or(0);
        self.rain_type.clear();
        self.rain_type.resize(un, 0);
        self.source.clear();
        self.source.resize(un, RainSource::Timeseries);
        self.ts_index.clear();
        self.ts_index.resize(un, -1);
        self.ts_name.clear();
        self.ts_name.resize(un, Vec::new());
        self.file_path.clear();
        self.file_path.resize(un, FilePathPair::default());
        self.col_name.clear();
        self.col_name.resize(un, Vec::new());
        self.station_id.clear();
        self.station_id.resize(un, Vec::new());
        self.rain_units.clear();
        self.rain_units.resize(un, 0);
        self.file_first_date.clear();
        self.file_first_date.resize(un, 0.0);
        self.file_last_date.clear();
        self.file_last_date.resize(un, 0.0);
        self.file_periods_precip.clear();
        self.file_periods_precip.resize(un, 0);
        self.file_start_date.clear();
        self.file_start_date.resize(un, 0.0);
        self.rain_series.clear();
        self.rain_series.resize_with(un, default_table);
        self.file_format.clear();
        self.file_format.resize(un, RainFileFormat::Unknown);
        self.interval_sec.clear();
        self.interval_sec.resize(un, 3600);
        self.snow_factor.clear();
        self.snow_factor.resize(un, 1.0);
        self.scale_factor.clear();
        self.scale_factor.resize(un, 1.0);

        self.rainfall.clear();
        self.rainfall.resize(un, 0.0);
        self.next_rainfall.clear();
        self.next_rainfall.resize(un, 0.0);
        self.api_rainfall.clear();
        self.api_rainfall.resize(un, -1.0); // -1.0 means no API override
        self.next_rain_date.clear();
        self.next_rain_date.resize(un, 0.0);
        self.is_raining.clear();
        self.is_raining.resize(un, false);
        self.st_cur.clear();
        self.st_cur.resize(un, -2);
        self.st_next.clear();
        self.st_next.resize(un, -1);
        self.st_rain.clear();
        self.st_rain.resize(un, 0.0);
        self.st_init.clear();
        self.st_init.resize(un, 0);
        self.st_used.clear();
        self.st_used.resize(un, 0);

        self.past_rain.clear();
        self.past_rain.resize(un * Self::MAXPASTRAIN as usize, 0.0);
        self.past_rain_accum.clear();
        self.past_rain_accum.resize(un, 0.0);
        self.past_rain_interval.clear();
        self.past_rain_interval.resize(un, 0.0);
        self.cumul_rain_accum.clear();
        self.cumul_rain_accum.resize(un, 0.0);
        self.co_gage_index.clear();
        self.co_gage_index.resize(un, -1);
        self.comments.clear();
        self.comments.resize(un, Vec::new());
    }

    pub fn grow_to(&mut self, n: i32) {
        if n <= self.count() {
            return;
        }
        let un = usize::try_from(n).unwrap_or(0);
        self.rain_type.resize(un, 0);
        self.source.resize(un, RainSource::Timeseries);
        self.ts_index.resize(un, -1);
        self.ts_name.resize(un, Vec::new());
        self.file_path.resize(un, FilePathPair::default());
        self.col_name.resize(un, Vec::new());
        self.station_id.resize(un, Vec::new());
        self.rain_units.resize(un, 0);
        self.file_first_date.resize(un, 0.0);
        self.file_last_date.resize(un, 0.0);
        self.file_periods_precip.resize(un, 0);
        self.file_start_date.resize(un, 0.0);
        self.rain_series.resize_with(un, default_table);
        self.file_format.resize(un, RainFileFormat::Unknown);
        self.interval_sec.resize(un, 3600);
        self.snow_factor.resize(un, 1.0);
        self.scale_factor.resize(un, 1.0);
        self.rainfall.resize(un, 0.0);
        self.next_rainfall.resize(un, 0.0);
        self.api_rainfall.resize(un, -1.0);
        self.next_rain_date.resize(un, 0.0);
        self.is_raining.resize(un, false);
        self.st_cur.resize(un, -2);
        self.st_next.resize(un, -1);
        self.st_rain.resize(un, 0.0);
        self.st_init.resize(un, 0);
        self.st_used.resize(un, 0);
        // Flat 2D: [gage * MAXPASTRAIN + hour]
        self.past_rain.resize(un * Self::MAXPASTRAIN as usize, 0.0);
        self.past_rain_accum.resize(un, 0.0);
        self.past_rain_interval.resize(un, 0.0);
        self.cumul_rain_accum.resize(un, 0.0);
        self.co_gage_index.resize(un, -1);
        self.comments.resize(un, Vec::new());
    }

    pub fn erase_at(&mut self, idx: i32) {
        let ui = idx as usize;
        if ui < self.rain_type.len() { self.rain_type.remove(ui); }
        if ui < self.source.len() { self.source.remove(ui); }
        if ui < self.ts_index.len() { self.ts_index.remove(ui); }
        if ui < self.ts_name.len() { self.ts_name.remove(ui); }
        if ui < self.file_path.len() { self.file_path.remove(ui); }
        if ui < self.col_name.len() { self.col_name.remove(ui); }
        if ui < self.file_format.len() { self.file_format.remove(ui); }
        if ui < self.interval_sec.len() { self.interval_sec.remove(ui); }
        if ui < self.snow_factor.len() { self.snow_factor.remove(ui); }
        if ui < self.scale_factor.len() { self.scale_factor.remove(ui); }
        if ui < self.file_start_date.len() { self.file_start_date.remove(ui); }
        if ui < self.rainfall.len() { self.rainfall.remove(ui); }
        if ui < self.next_rainfall.len() { self.next_rainfall.remove(ui); }
        if ui < self.api_rainfall.len() { self.api_rainfall.remove(ui); }
        if ui < self.next_rain_date.len() { self.next_rain_date.remove(ui); }
        if ui < self.is_raining.len() { self.is_raining.remove(ui); }
        if ui < self.st_cur.len() { self.st_cur.remove(ui); }
        if ui < self.st_next.len() { self.st_next.remove(ui); }
        if ui < self.st_rain.len() { self.st_rain.remove(ui); }
        if ui < self.st_init.len() { self.st_init.remove(ui); }
        if ui < self.st_used.len() { self.st_used.remove(ui); }
        if ui < self.past_rain_accum.len() { self.past_rain_accum.remove(ui); }
        if ui < self.past_rain_interval.len() { self.past_rain_interval.remove(ui); }
        if ui < self.cumul_rain_accum.len() { self.cumul_rain_accum.remove(ui); }
        if ui < self.co_gage_index.len() { self.co_gage_index.remove(ui); }
        if ui < self.comments.len() { self.comments.remove(ui); }

        // Flat 2D past_rain: [gage * MAXPASTRAIN + hour]
        let base = ui * Self::MAXPASTRAIN as usize;
        let end = base + Self::MAXPASTRAIN as usize;
        if end <= self.past_rain.len() {
            self.past_rain.drain(base..end);
        }
    }

    pub fn shrink_to_fit(&mut self) {
        self.rain_type.shrink_to_fit();
        self.source.shrink_to_fit();
        self.ts_index.shrink_to_fit();
        self.ts_name.shrink_to_fit();
        self.file_path.shrink_to_fit();
        self.col_name.shrink_to_fit();
        self.file_format.shrink_to_fit();
        self.interval_sec.shrink_to_fit();
        self.snow_factor.shrink_to_fit();
        self.scale_factor.shrink_to_fit();

        self.rainfall.shrink_to_fit();
        self.next_rainfall.shrink_to_fit();
        self.api_rainfall.shrink_to_fit();
        self.next_rain_date.shrink_to_fit();
        self.is_raining.shrink_to_fit();
        self.st_cur.shrink_to_fit();
        self.st_next.shrink_to_fit();
        self.st_rain.shrink_to_fit();
        self.st_init.shrink_to_fit();
        self.st_used.shrink_to_fit();

        self.past_rain.shrink_to_fit();
        self.past_rain_accum.shrink_to_fit();
        self.past_rain_interval.shrink_to_fit();
        self.cumul_rain_accum.shrink_to_fit();
        self.co_gage_index.shrink_to_fit();
        self.comments.shrink_to_fit();
    }

    pub fn reset_state(&mut self) {
        self.rainfall.fill(0.0);
        self.next_rainfall.fill(0.0);
        self.api_rainfall.fill(-1.0); // -1 = no override
        self.is_raining.fill(false);
        self.st_cur.fill(-2);
        self.st_next.fill(-1);
        self.st_rain.fill(0.0);
        self.st_init.fill(0);
        self.st_used.fill(0);
        self.past_rain.fill(0.0);
        self.past_rain_accum.fill(0.0);
        self.past_rain_interval.fill(0.0);
        self.cumul_rain_accum.fill(0.0);
        // co_gage_index is static (set at parse time), not reset between runs
    }
}
