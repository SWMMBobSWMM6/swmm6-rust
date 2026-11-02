// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/TableData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The in-memory tables and every lookup on them, operation for operation. The functions
// that page an external time-series file (TableData.cpp: validate_table, table_open_file,
// table_load_cache, table_lookup_column, table_step_column) are translated with the rest
// of TableData.cpp; the fields they use are here already.

use crate::core::cmath::{cclamp, cmin};
use crate::core::options::{Bytes, FilePathPair};
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum TableType {
    Timeseries = 0,
    CurveStorage = 1,
    CurveDiversion = 2,
    CurveRating = 3,
    CurveShape = 4,
    CurveControl = 5,
    CurveTidal = 6,
    CurvePump1 = 7,
    CurvePump2 = 8,
    CurvePump3 = 9,
    CurvePump4 = 10,
    CurvePump5 = 11,
    CurveWeir = 12,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TableCursor {
    pub index: i32,
    pub direction: i32,
}

impl Default for TableCursor {
    fn default() -> Self {
        TableCursor { index: 0, direction: 1 }
    }
}

impl TableCursor {
    pub fn reset(&mut self) {
        self.index = 0;
        self.direction = 1;
    }
}

#[derive(Debug, Clone)]
pub struct TableBlock {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub num_cols: i32,
    pub file_row_start: usize,
}

impl Default for TableBlock {
    fn default() -> Self {
        TableBlock { x: Vec::new(), y: Vec::new(), num_cols: 1, file_row_start: 0 }
    }
}

impl TableBlock {
    pub fn num_rows(&self) -> usize {
        self.x.len()
    }
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }
}

#[derive(Debug)]
pub struct Table {
    pub id: Bytes,
    pub type_: TableType,
    pub comment: Bytes,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub cursor: TableCursor,
    pub n_relative: i32,
    pub rel_anchor: f64,
    // ---- external file (paged)
    pub is_file_based: bool,
    pub file_handle: Option<std::io::BufReader<std::fs::File>>,
    pub file_path: FilePathPair,
    pub first_boundary: TableBlock,
    pub last_boundary: TableBlock,
    pub cache: TableBlock,
    pub dx_min: f64,
    pub x_min: f64,
    pub x_max: f64,
    pub num_cols: i32,
    pub total_rows: usize,
    pub num_cache_rows: usize,
    pub data_start_offset: i64,
    pub row_offsets: Vec<i64>,
    pub column_ids: Vec<Bytes>,
    pub column_map: HashMap<Bytes, i32>,
}

impl Table {
    pub const INDEX_STRIDE: usize = 4096;

    pub fn new(id: &[u8], type_: TableType) -> Self {
        Table {
            id: id.to_vec(),
            type_,
            comment: Vec::new(),
            x: Vec::new(),
            y: Vec::new(),
            cursor: TableCursor::default(),
            n_relative: 0,
            rel_anchor: 0.0,
            is_file_based: false,
            file_handle: None,
            file_path: FilePathPair::default(),
            first_boundary: TableBlock::default(),
            last_boundary: TableBlock::default(),
            cache: TableBlock::default(),
            dx_min: 0.0,
            x_min: 0.0,
            x_max: 0.0,
            num_cols: 1,
            total_rows: 0,
            num_cache_rows: 8192,
            data_start_offset: 0,
            row_offsets: Vec::new(),
            column_ids: Vec::new(),
            column_map: HashMap::new(),
        }
    }
    pub fn size(&self) -> usize {
        self.x.len()
    }
    pub fn is_empty(&self) -> bool {
        self.x.is_empty()
    }
}

/// Linear interpolation with a moving cursor; constant beyond both ends.
pub fn table_lookup_cursor(tbl: &mut Table, x_query: f64) -> f64 {
    let n = tbl.x.len() as i32;
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return tbl.y[0];
    }
    if x_query <= tbl.x[0] {
        tbl.cursor.index = 0;
        tbl.cursor.direction = 1;
        return tbl.y[0];
    }
    if x_query > tbl.x[(n - 1) as usize] {
        tbl.cursor.index = n - 1;
        tbl.cursor.direction = -1;
        return tbl.y[(n - 1) as usize];
    }
    let mut idx = cclamp(tbl.cursor.index, 0, n - 2);
    while idx < n - 1 && tbl.x[(idx + 1) as usize] < x_query {
        idx += 1;
        tbl.cursor.direction = 1;
    }
    while idx > 0 && tbl.x[idx as usize] >= x_query {
        idx -= 1;
        tbl.cursor.direction = -1;
    }
    tbl.cursor.index = idx;
    let (x1, y1) = (tbl.x[idx as usize], tbl.y[idx as usize]);
    let (x2, y2) = (tbl.x[(idx + 1) as usize], tbl.y[(idx + 1) as usize]);
    let dx = x2 - x1;
    if dx.abs() < 1.0e-20 {
        return (y1 + y2) / 2.;
    }
    y1 + (x_query - x1) * (y2 - y1) / dx
}

/// Time-series interpolation: zero outside the series, the end value exactly at an end.
pub fn table_tseries_lookup_cursor(tbl: &mut Table, x_query: f64) -> f64 {
    let n = tbl.x.len() as i32;
    if n == 0 {
        return 0.0;
    }
    if x_query <= tbl.x[0] {
        tbl.cursor.index = 0;
        tbl.cursor.direction = 1;
        return if x_query < tbl.x[0] { 0.0 } else { tbl.y[0] };
    }
    let last = (n - 1) as usize;
    if x_query >= tbl.x[last] {
        tbl.cursor.index = n - 1;
        tbl.cursor.direction = -1;
        return if x_query > tbl.x[last] { 0.0 } else { tbl.y[last] };
    }
    table_lookup_cursor(tbl, x_query)
}

/// Step (zero-order hold) lookup.
pub fn table_step_cursor(tbl: &mut Table, x_query: f64) -> f64 {
    let n = tbl.x.len() as i32;
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return tbl.y[0];
    }
    if x_query < tbl.x[0] {
        tbl.cursor.index = 0;
        tbl.cursor.direction = 1;
        return 0.0;
    }
    if x_query >= tbl.x[(n - 1) as usize] {
        tbl.cursor.index = n - 1;
        tbl.cursor.direction = -1;
        return tbl.y[(n - 1) as usize];
    }
    let mut idx = cclamp(tbl.cursor.index, 0, n - 2);
    while idx < n - 1 && tbl.x[(idx + 1) as usize] <= x_query {
        idx += 1;
        tbl.cursor.direction = 1;
    }
    while idx > 0 && tbl.x[idx as usize] > x_query {
        idx -= 1;
        tbl.cursor.direction = -1;
    }
    tbl.cursor.index = idx;
    tbl.y[idx as usize]
}

/// Volume under a storage area curve up to `depth` (trapezoids, extrapolated past the end).
pub fn table_get_storage_volume(tbl: &Table, depth: f64) -> f64 {
    let n = tbl.x.len();
    if n == 0 || depth <= 0.0 {
        return 0.0;
    }
    let mut x1 = tbl.x[0];
    let mut a1 = tbl.y[0];
    if depth <= x1 {
        if x1 < 1.0e-6 {
            return 0.0;
        }
        return (a1 / x1) * depth * depth / 2.0;
    }
    let mut v = 0.0;
    let (mut dx, mut dy) = (0.0, 0.0);
    for i in 1..n {
        let x2 = tbl.x[i];
        let a2 = tbl.y[i];
        if x2 >= depth {
            let dxi = x2 - x1;
            let a = if dxi.abs() < 1.0e-20 { (a1 + a2) / 2. } else { a1 + (depth - x1) * (a2 - a1) / dxi };
            return v + (a1 + a) / 2.0 * (depth - x1);
        }
        dx = x2 - x1;
        dy = a2 - a1;
        v += (a1 + a2) / 2.0 * dx;
        x1 = x2;
        a1 = a2;
    }
    if dx > 1.0e-6 {
        let s = dy / dx;
        let a = a1 + s * (depth - x1);
        if a < 0.0 {
            v -= a1 * a1 / s / 2.0;
        } else {
            v += (a1 + a) / 2.0 * (depth - x1);
        }
    }
    v
}

/// Depth at which a storage area curve holds `volume`.
pub fn table_get_storage_depth(tbl: &Table, volume: f64) -> f64 {
    let n = tbl.x.len();
    if n == 0 || volume <= 0.0 {
        return 0.0;
    }
    let mut x1 = tbl.x[0];
    let mut a1 = tbl.y[0];
    let mut v_accum = 0.0;
    for i in 1..n {
        let x2 = tbl.x[i];
        let a2 = tbl.y[i];
        let dv = (a1 + a2) / 2.0 * (x2 - x1);
        if v_accum + dv >= volume {
            let dx = x2 - x1;
            let dv_need = volume - v_accum;
            if dx < 1.0e-10 {
                return x1;
            }
            let s = (a2 - a1) / dx;
            let depth = if s.abs() < 1.0e-10 {
                if a1 > 0.0 {
                    x1 + dv_need / a1
                } else {
                    x1
                }
            } else {
                let mut disc = a1 * a1 + 2.0 * s * dv_need;
                if disc < 0.0 {
                    disc = 0.0;
                }
                x1 + (-a1 + disc.sqrt()) / s
            };
            return cmin(depth, x2);
        }
        v_accum += dv;
        x1 = x2;
        a1 = a2;
    }
    if n >= 2 {
        let last = n - 1;
        let dx = tbl.x[last] - tbl.x[last - 1];
        let da = tbl.y[last] - tbl.y[last - 1];
        let s = if dx > 1.0e-10 { da / dx } else { 0.0 };
        let dv_need = volume - v_accum;
        let a1_ext = tbl.y[last];
        let depth = if s.abs() < 1.0e-10 {
            if a1_ext > 0.0 {
                tbl.x[last] + dv_need / a1_ext
            } else {
                tbl.x[last]
            }
        } else {
            let mut disc = a1_ext * a1_ext + 2.0 * s * dv_need;
            if disc < 0.0 {
                disc = 0.0;
            }
            tbl.x[last] + (-a1_ext + disc.sqrt()) / s
        };
        return depth;
    }
    tbl.x[n - 1]
}

/// x at which y first reaches `y_query` (y assumed non-decreasing).
pub fn table_inverse_lookup(tbl: &Table, y_query: f64) -> f64 {
    let n = tbl.x.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return tbl.x[0];
    }
    if y_query <= tbl.y[0] {
        return tbl.x[0];
    }
    if y_query >= tbl.y[n - 1] {
        return tbl.x[n - 1];
    }
    for i in 1..n {
        if tbl.y[i] >= y_query {
            let dy = tbl.y[i] - tbl.y[i - 1];
            if dy <= 0.0 {
                return tbl.x[i - 1];
            }
            let t = (y_query - tbl.y[i - 1]) / dy;
            return tbl.x[i - 1] + t * (tbl.x[i] - tbl.x[i - 1]);
        }
    }
    tbl.x[n - 1]
}

/// Interpolation with linear extrapolation past the end (never with a negative slope) and
/// proportional scaling below the first point.
pub fn table_lookup_ex(tbl: &Table, x_query: f64) -> f64 {
    let n = tbl.x.len();
    if n == 0 {
        return 0.0;
    }
    let mut x1 = tbl.x[0];
    let mut y1 = tbl.y[0];
    let mut s = 0.0;
    if x_query <= x1 {
        if x1 > 0.0 {
            return x_query / x1 * y1;
        } else {
            return y1;
        }
    }
    for i in 1..n {
        let x2 = tbl.x[i];
        let y2 = tbl.y[i];
        if x2 != x1 {
            s = (y2 - y1) / (x2 - x1);
        }
        if x_query <= x2 {
            let dx = x2 - x1;
            if dx.abs() < 1.0e-20 {
                return (y1 + y2) / 2.;
            }
            return y1 + (x_query - x1) * (y2 - y1) / dx;
        }
        x1 = x2;
        y1 = y2;
    }
    if s < 0.0 {
        s = 0.0;
    }
    y1 + s * (x_query - x1)
}

/// The y of the first entry whose x exceeds `x_query`, else the last y.
pub fn table_interval_lookup(tbl: &Table, x_query: f64) -> f64 {
    let n = tbl.x.len();
    if n == 0 {
        return 0.0;
    }
    for i in 0..n {
        if tbl.x[i] > x_query {
            return tbl.y[i];
        }
    }
    tbl.y[n - 1]
}

pub fn table_get_slope(tbl: &Table, x_query: f64) -> f64 {
    let n = tbl.x.len();
    if n < 2 {
        return 0.0;
    }
    if x_query <= tbl.x[0] {
        let dx = tbl.x[1] - tbl.x[0];
        return if dx > 0.0 { (tbl.y[1] - tbl.y[0]) / dx } else { 0.0 };
    }
    if x_query > tbl.x[n - 1] {
        return 0.0;
    }
    for i in 1..n {
        if tbl.x[i] >= x_query {
            let dx = tbl.x[i] - tbl.x[i - 1];
            return if dx > 0.0 { (tbl.y[i] - tbl.y[i - 1]) / dx } else { 0.0 };
        }
    }
    0.0
}

/// The largest y reached before y first decreases.
pub fn table_get_max_y(tbl: &Table) -> f64 {
    let n = tbl.x.len();
    if n == 0 {
        return 0.0;
    }
    let mut ymax = tbl.y[0];
    for i in 1..n {
        if tbl.y[i] < ymax {
            break;
        }
        ymax = tbl.y[i];
    }
    ymax
}

/// All tables in index order, findable by name (case-insensitive). A time series and a
/// curve may share a name, so a name maps to every index that carries it.
#[derive(Debug, Default)]
pub struct TableData {
    pub tables: Vec<Table>,
    pub by_name: HashMap<Bytes, Vec<i32>>,
}

fn fold(s: &[u8]) -> Bytes {
    s.iter().map(|&c| crate::data::name_index::ci_upper(c)).collect()
}

impl TableData {
    pub fn count(&self) -> usize {
        self.tables.len()
    }

    pub fn add(&mut self, id: &[u8], type_: TableType) -> i32 {
        self.tables.push(Table::new(id, type_));
        let idx = self.tables.len() as i32 - 1;
        self.by_name.entry(fold(id)).or_default().push(idx);
        idx
    }

    pub fn rebuild_index(&mut self) {
        self.by_name.clear();
        for (i, t) in self.tables.iter().enumerate() {
            self.by_name.entry(fold(&t.id)).or_default().push(i as i32);
        }
    }

    pub fn find_by_kind(&self, name: &[u8], want_timeseries: bool) -> i32 {
        match self.by_name.get(&fold(name)) {
            None => -1,
            Some(v) => {
                for &i in v {
                    let is_ts = self.tables[i as usize].type_ == TableType::Timeseries;
                    if is_ts == want_timeseries {
                        return i;
                    }
                }
                -1
            }
        }
    }

    pub fn reset_cursors(&mut self) {
        for t in &mut self.tables {
            t.cursor.reset();
        }
    }
}

impl std::ops::Index<i32> for TableData {
    type Output = Table;
    fn index(&self, i: i32) -> &Table {
        &self.tables[i as usize]
    }
}

impl std::ops::IndexMut<i32> for TableData {
    fn index_mut(&mut self, i: i32) -> &mut Table {
        &mut self.tables[i as usize]
    }
}

#[derive(Debug, Clone)]
pub struct TableValidation {
    pub valid: bool,
    pub errors: Vec<Bytes>,
    pub warnings: Vec<Bytes>,
}

impl Default for TableValidation {
    fn default() -> Self {
        TableValidation { valid: true, errors: Vec::new(), warnings: Vec::new() }
    }
}
