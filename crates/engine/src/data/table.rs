// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/TableData.hpp and
// src/engine/data/TableData.cpp (HydroCouple/openswmm.engine, branch swmm6_rel,
// commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the Apache License,
// Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The in-memory tables and every lookup on them, operation for operation, then
// TableData.cpp: validate_table and the paged external-file reader (table_open_file,
// table_load_cache, table_lookup_column, table_step_column).
//
// Differences in form, never in result:
//   * `std::FILE*` is `CFile`, a buffered file with the four stdio calls upstream makes
//     (fread, fgets with upstream's buffer sizes, fseek, ftell). A read error reads as end
//     of file, as glibc's fgets returns NULL for both.
//   * `std::sscanf` and `strtol` are emulated (glibc, C locale) by `sscanf` below, which
//     covers exactly the conversions upstream's format strings use (%d, %Ns, %n, %*[..]).
//   * The four identical "skip N data rows" loops are one helper, `skip_data_rows`.
//   * `prepare_lookup` returns (n, stride, cursor) and the callers index the cache or
//     the in-memory arrays directly instead of through raw pointers.

use crate::core::charconv::{from_chars_double, from_chars_i32};
use crate::core::cmath::{cclamp, cmin};
use crate::core::datetime;
use crate::core::error_codes::{
    format_error, format_error_detail, format_warning, ERR_CURVE_SEQUENCE, ERR_TABLE_COL_MISMATCH,
    ERR_TABLE_FILE_OPEN, ERR_TABLE_FILE_READ, ERR_TIMESERIES_EMPTY, ERR_TIMESERIES_NAN,
    ERR_TIMESERIES_SEQUENCE, WARN_BOUNDARY_OVERLAP, WARN_TIMESERIES_DUPLICATE_X,
};
use crate::core::options::{Bytes, FilePathPair};
use std::collections::HashMap;
use std::io::{BufRead, Seek, SeekFrom};

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
    pub file_handle: Option<CFile>,
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

// ============================================================================
// TableData.cpp
// ============================================================================

/// `boundary_rows`' default in `table_open_file(Table&, std::size_t boundary_rows = 128)`.
pub const TABLE_OPEN_FILE_DEFAULT_BOUNDARY_ROWS: usize = 128;

// ---- the C library pieces upstream calls -----------------------------------

/// `isspace` in the C locale.
#[inline]
pub(crate) fn c_isspace(c: u8) -> bool {
    matches!(c, b' ' | b'\t' | b'\n' | 0x0b | 0x0c | b'\r')
}

/// What a C string function sees of a buffer: everything before the first NUL.
#[inline]
pub(crate) fn c_str(b: &[u8]) -> &[u8] {
    match b.iter().position(|&c| c == 0) {
        Some(p) => &b[..p],
        None => b,
    }
}

/// A `std::FILE*` opened with mode "r" (no newline translation on the reference platform),
/// with the four stdio calls TableData.cpp makes on it.
#[derive(Debug)]
pub struct CFile {
    r: std::io::BufReader<std::fs::File>,
}

impl CFile {
    /// `fopen(path, "r")`, or None where fopen returns NULL. As with glibc, a directory
    /// opens and then reads nothing.
    pub fn open(path: &[u8]) -> Option<CFile> {
        std::fs::File::open(crate::core::file_io::path_from_bytes(path))
            .ok()
            .map(|f| CFile { r: std::io::BufReader::new(f) })
    }

    /// `fread(buf, 1, n, fp)`: the bytes read, fewer than `n` at end of file.
    pub fn fread(&mut self, n: usize) -> Vec<u8> {
        let mut out = Vec::with_capacity(n);
        while out.len() < n {
            let avail = match self.r.fill_buf() {
                Ok(b) => b,
                Err(_) => break,
            };
            if avail.is_empty() {
                break;
            }
            let take = avail.len().min(n - out.len());
            out.extend_from_slice(&avail[..take]);
            self.r.consume(take);
        }
        out
    }

    /// `fseek(fp, off, SEEK_SET)`: 0, or -1 for a negative offset or a failed seek.
    pub fn fseek(&mut self, off: i64) -> i32 {
        if off < 0 {
            return -1;
        }
        match self.r.seek(SeekFrom::Start(off as u64)) {
            Ok(_) => 0,
            Err(_) => -1,
        }
    }

    /// `ftell(fp)`, -1 on failure.
    pub fn ftell(&mut self) -> i64 {
        match self.r.stream_position() {
            Ok(p) => p as i64,
            Err(_) => -1,
        }
    }

    /// `fgets(buf, size, fp)`: up to `size - 1` bytes, through the first '\n'. On success
    /// `buf` holds exactly those bytes (the C buffer's NUL terminator is implied) and the
    /// result is true. At end of file with nothing read it is false and `buf` is left as
    /// it was, as fgets leaves the array alone.
    pub fn fgets(&mut self, size: usize, buf: &mut Bytes) -> bool {
        let max = size - 1;
        let mut chunk = Bytes::new();
        while chunk.len() < max {
            let avail = match self.r.fill_buf() {
                Ok(b) => b,
                Err(_) => break,
            };
            if avail.is_empty() {
                break;
            }
            let take = avail.len().min(max - chunk.len());
            if let Some(p) = avail[..take].iter().position(|&c| c == b'\n') {
                chunk.extend_from_slice(&avail[..=p]);
                self.r.consume(p + 1);
                break;
            }
            chunk.extend_from_slice(&avail[..take]);
            self.r.consume(take);
        }
        if chunk.is_empty() {
            return false;
        }
        *buf = chunk;
        true
    }
}

/// One destination of a `sscanf` call.
pub(crate) enum ScanArg<'a> {
    /// `%d` into an `int`.
    Int(&'a mut i32),
    /// `%Ns` into a `char[N+1]`.
    Str(&'a mut Bytes),
    /// `%n` into an `int`.
    Count(&'a mut i32),
}

/// glibc `strtol(s, NULL, 10)` on a buffer `%d` collected (a sign, then digits):
/// LONG_MAX / LONG_MIN when out of range.
fn strtol_digits(b: &[u8]) -> i64 {
    let (neg, digits) = match b.first() {
        Some(b'-') => (true, &b[1..]),
        Some(b'+') => (false, &b[1..]),
        _ => (false, b),
    };
    let mut v: i128 = 0;
    for &c in digits {
        v = v * 10 + (c - b'0') as i128;
        if v > i64::MAX as i128 + 1 {
            v = i64::MAX as i128 + 2; // saturated; stays out of range either way
        }
    }
    if neg {
        v = -v;
    }
    if v > i64::MAX as i128 {
        i64::MAX
    } else if v < i64::MIN as i128 {
        i64::MIN
    } else {
        v as i64
    }
}

/// glibc `sscanf` in the C locale, for the conversions upstream uses: `%d` (strtol, then
/// truncated to int), `%Ns`, `%n`, and a suppressed scan set `%*[...]`. White space in
/// the format matches any amount of input white space; other bytes match themselves.
/// Returns the number of assignments, or -1 (EOF) when input ran out before the first.
/// The input is a C string: it ends at the first NUL.
pub(crate) fn sscanf(input: &[u8], fmt: &[u8], args: &mut [ScanArg<'_>]) -> i32 {
    let s = c_str(input);
    let mut pos = 0usize; // glibc's read_in
    let mut done: i32 = 0;
    let mut ai = 0usize;
    let mut fi = 0usize;
    let mut skip_space = false;
    macro_rules! input_error {
        () => {
            return if done == 0 { -1 } else { done }
        };
    }
    while fi < fmt.len() {
        let fc = fmt[fi];
        fi += 1;
        if fc != b'%' {
            if c_isspace(fc) {
                skip_space = true;
                continue;
            }
            if pos >= s.len() {
                input_error!();
            }
            let mut c = s[pos];
            pos += 1;
            if skip_space {
                while c_isspace(c) {
                    if pos >= s.len() {
                        input_error!();
                    }
                    c = s[pos];
                    pos += 1;
                }
                skip_space = false;
            }
            if c != fc {
                return done; // conv_error (the byte is pushed back)
            }
            continue;
        }
        let mut suppress = false;
        if fi < fmt.len() && fmt[fi] == b'*' {
            suppress = true;
            fi += 1;
        }
        let mut width: i32 = -1;
        if fi < fmt.len() && fmt[fi].is_ascii_digit() {
            width = 0;
            while fi < fmt.len() && fmt[fi].is_ascii_digit() {
                width = width * 10 + (fmt[fi] - b'0') as i32;
                fi += 1;
            }
        }
        let conv = fmt[fi];
        fi += 1;
        if skip_space || !matches!(conv, b'[' | b'c' | b'n') {
            while pos < s.len() && c_isspace(s[pos]) {
                pos += 1;
            }
            skip_space = false;
        }
        match conv {
            b'n' => {
                if !suppress {
                    if let ScanArg::Count(p) = &mut args[ai] {
                        **p = pos as i32;
                    } else {
                        unreachable!("%n needs ScanArg::Count");
                    }
                    ai += 1;
                }
            }
            b'd' => {
                if pos >= s.len() {
                    input_error!();
                }
                let mut buf = Bytes::new();
                let mut w = width;
                if s[pos] == b'-' || s[pos] == b'+' {
                    buf.push(s[pos]);
                    pos += 1;
                    if w > 0 {
                        w -= 1;
                    }
                }
                while pos < s.len() && w != 0 && s[pos].is_ascii_digit() {
                    buf.push(s[pos]);
                    pos += 1;
                    if w > 0 {
                        w -= 1;
                    }
                }
                if buf.is_empty() || (buf.len() == 1 && (buf[0] == b'+' || buf[0] == b'-')) {
                    return done; // conv_error
                }
                let v = strtol_digits(&buf) as i32;
                if !suppress {
                    if let ScanArg::Int(p) = &mut args[ai] {
                        **p = v;
                    } else {
                        unreachable!("%d needs ScanArg::Int");
                    }
                    ai += 1;
                    done += 1;
                }
            }
            b's' => {
                if pos >= s.len() {
                    input_error!();
                }
                let mut w = width;
                let mut out = Bytes::new();
                loop {
                    let c = s[pos];
                    if c_isspace(c) {
                        break;
                    }
                    out.push(c);
                    pos += 1;
                    let more = if w <= 0 {
                        true
                    } else {
                        w -= 1;
                        w > 0
                    };
                    if !more || pos >= s.len() {
                        break;
                    }
                }
                if !suppress {
                    if let ScanArg::Str(p) = &mut args[ai] {
                        **p = out;
                    } else {
                        unreachable!("%s needs ScanArg::Str");
                    }
                    ai += 1;
                    done += 1;
                }
            }
            b'[' => {
                let mut not_in = false;
                if fi < fmt.len() && fmt[fi] == b'^' {
                    not_in = true;
                    fi += 1;
                }
                let set_start = fi;
                if fi < fmt.len() && fmt[fi] == b']' {
                    fi += 1;
                }
                while fi < fmt.len() && fmt[fi] != b']' {
                    fi += 1;
                }
                let set = &fmt[set_start..fi];
                fi += 1; // the ']'
                if pos >= s.len() {
                    input_error!();
                }
                let now = pos;
                let mut w: i32 = if width < 0 { i32::MAX } else { width };
                loop {
                    let c = s[pos];
                    if set.contains(&c) == not_in {
                        break;
                    }
                    pos += 1;
                    w -= 1;
                    if w <= 0 || pos >= s.len() {
                        break;
                    }
                }
                if pos == now {
                    return done; // conv_error
                }
                assert!(suppress, "only the suppressed scan set is used");
            }
            _ => unreachable!("conversion not used by upstream"),
        }
    }
    done
}

// ---- internal helpers (anonymous namespace) ---------------------------------

/// The comment/blank test every row loop makes on the first byte of an fgets chunk.
#[inline]
fn is_skipped_row(b0: u8) -> bool {
    b0 == b';' || b0 == b'\n' || b0 == b'\r' || b0 == b'#'
}

/// Heuristically identify a delimited text header: a comma or tab and an ASCII letter
/// both before the line end (CR, LF or NUL).
fn looks_like_csv_header(line: &[u8]) -> bool {
    let mut has_delim = false;
    let mut has_alpha = false;
    for &c in c_str(line) {
        if c == b'\n' || c == b'\r' {
            break;
        }
        if c == b',' || c == b'\t' {
            has_delim = true;
        }
        if c.is_ascii_uppercase() || c.is_ascii_lowercase() {
            has_alpha = true;
        }
    }
    has_delim && has_alpha
}

/// Split one unquoted comma-or-tab line (any comma anywhere selects comma) and trim
/// ASCII spaces. Interior empty fields are kept; a final delimiter adds no field.
fn split_csv_line(line: &[u8]) -> Vec<Bytes> {
    let line = c_str(line);
    let mut tokens = Vec::new();
    let mut delim = b',';
    if !line.contains(&b',') {
        delim = b'\t';
    }
    let at = |p: usize| if p < line.len() { line[p] } else { 0 };
    let mut p = 0usize;
    while at(p) != 0 && at(p) != b'\n' && at(p) != b'\r' {
        while at(p) == b' ' {
            p += 1;
        }
        let start = p;
        while at(p) != 0 && at(p) != delim && at(p) != b'\n' && at(p) != b'\r' {
            p += 1;
        }
        let mut end = p;
        while end > start && line[end - 1] == b' ' {
            end -= 1;
        }
        tokens.push(line[start..end].to_vec());
        if at(p) == delim {
            p += 1;
        }
    }
    tokens
}

/// Parse YYYY-MM-DD or M/D/Y, optionally followed by H:M[:S], into an OADate; -1 when
/// no supported date prefix scans.
fn parse_csv_datetime(tok: &[u8]) -> f64 {
    use ScanArg::Int;
    let (mut month, mut day, mut year) = (0i32, 0i32, 0i32);
    let (mut hour, mut minute, mut second) = (0i32, 0i32, 0i32);

    if sscanf(
        tok,
        b"%d-%d-%d %d:%d:%d",
        &mut [Int(&mut year), Int(&mut month), Int(&mut day), Int(&mut hour), Int(&mut minute), Int(&mut second)],
    ) >= 5
    {
        return datetime::encode_date(year, month, day) + datetime::encode_time(hour, minute, second);
    }
    if sscanf(
        tok,
        b"%d/%d/%d %d:%d:%d",
        &mut [Int(&mut month), Int(&mut day), Int(&mut year), Int(&mut hour), Int(&mut minute), Int(&mut second)],
    ) >= 5
    {
        return datetime::encode_date(year, month, day) + datetime::encode_time(hour, minute, second);
    }
    if sscanf(tok, b"%d-%d-%d", &mut [Int(&mut year), Int(&mut month), Int(&mut day)]) == 3 {
        return datetime::encode_date(year, month, day);
    }
    if sscanf(tok, b"%d/%d/%d", &mut [Int(&mut month), Int(&mut day), Int(&mut year)]) == 3 {
        return datetime::encode_date(year, month, day);
    }
    -1.0
}

/// Parse a "M/D/Y H:M[:S] v1 v2 ..." line: the number of values appended to `vals_out`
/// (cleared first), 0 for a malformed date or time. `x_out` is set once the date and time
/// scan, even when no value follows.
fn parse_dat_line(line: &[u8], x_out: &mut f64, vals_out: &mut Vec<f64>) -> i32 {
    use ScanArg::{Count, Int, Str};
    vals_out.clear();

    let mut date_str = Bytes::new();
    let mut time_str = Bytes::new();
    let mut offset: i32 = 0;
    if sscanf(line, b"%31s %31s%n", &mut [Str(&mut date_str), Str(&mut time_str), Count(&mut offset)]) < 2 {
        return 0;
    }

    let (mut month, mut day, mut year) = (0i32, 0i32, 0i32);
    if sscanf(&date_str, b"%d/%d/%d", &mut [Int(&mut month), Int(&mut day), Int(&mut year)]) != 3 {
        return 0;
    }

    let (mut hour, mut minute, mut second) = (0i32, 0i32, 0i32);
    if sscanf(&time_str, b"%d:%d:%d", &mut [Int(&mut hour), Int(&mut minute), Int(&mut second)]) < 2 {
        return 0;
    }

    *x_out = datetime::encode_date(year, month, day) + datetime::encode_time(hour, minute, second);

    let line = c_str(line);
    let at = |p: usize| if p < line.len() { line[p] } else { 0 };
    let mut p = offset as usize;
    while at(p) != 0 {
        while at(p) == b' ' || at(p) == b'\t' {
            p += 1;
        }
        if at(p) == b'\n' || at(p) == b'\r' || at(p) == 0 {
            break;
        }
        let (r, used) = from_chars_double(&line[p..]);
        match r {
            Ok(v) => {
                vals_out.push(v);
                p += used;
            }
            Err(_) => break,
        }
    }
    vals_out.len() as i32
}

/// Count data rows from `data_start` (comment and blank fgets chunks excluded, 1023 bytes
/// at most per chunk) and record the byte offset after every INDEX_STRIDE-th row. Leaves
/// the file at `data_start`.
fn count_and_index_rows(fp: &mut CFile, data_start: i64, row_offsets: &mut Vec<i64>) -> usize {
    fp.fseek(data_start);
    row_offsets.clear();
    row_offsets.push(data_start); // row 0 starts at data_start

    let mut count = 0usize;
    let mut buf = Bytes::new(); // char buf[1024]
    loop {
        let _line_start = fp.ftell(); // read and never used upstream
        if !fp.fgets(1024, &mut buf) {
            break;
        }
        if is_skipped_row(buf[0]) {
            continue;
        }
        count += 1;
        // Record offset at every INDEX_STRIDE boundary
        if count.is_multiple_of(Table::INDEX_STRIDE) {
            row_offsets.push(fp.ftell());
        }
    }
    fp.fseek(data_start);
    count
}

/// Read up to `max_rows` data rows (2047-byte fgets chunks) into flat row-major x / y.
fn read_rows_flat(
    fp: &mut CFile,
    is_csv: bool,
    ncols: i32,
    max_rows: usize,
    x_out: &mut Vec<f64>,
    y_out: &mut Vec<f64>,
) -> usize {
    x_out.clear();
    y_out.clear();

    let mut buf = Bytes::new(); // char buf[2048]
    let mut rows_read = 0usize;

    while rows_read < max_rows && fp.fgets(2048, &mut buf) {
        if is_skipped_row(buf[0]) {
            continue;
        }

        if is_csv {
            let tokens = split_csv_line(&buf);
            if tokens.is_empty() {
                continue;
            }

            let dt = parse_csv_datetime(&tokens[0]);
            if dt < 0.0 {
                continue;
            }

            x_out.push(dt);
            for c in 0..ncols {
                let mut v = 0.0;
                if ((c + 1) as usize) < tokens.len() {
                    if let Ok(parsed) = from_chars_double(&tokens[(c + 1) as usize]).0 {
                        v = parsed;
                    }
                }
                y_out.push(v);
            }
            rows_read += 1;
        } else {
            let mut xv = 0.0;
            let mut vals = Vec::new();
            let nv = parse_dat_line(&buf, &mut xv, &mut vals);
            if nv == 0 {
                continue;
            }

            x_out.push(xv);
            for c in 0..ncols {
                y_out.push(if c < nv { vals[c as usize] } else { 0.0 });
            }
            rows_read += 1;
        }
    }
    rows_read
}

/// Skip `n` data rows (2047-byte fgets chunks; comment and blank chunks do not count),
/// stopping early at end of file. Upstream writes this loop out four times.
fn skip_data_rows(fp: &mut CFile, n: usize) {
    let mut buf = Bytes::new(); // char buf[2048]
    let mut skipped = 0usize;
    while skipped < n && fp.fgets(2048, &mut buf) {
        if buf[0] == b';' || buf[0] == b'#' || buf[0] == b'\n' || buf[0] == b'\r' {
            continue;
        }
        skipped += 1;
    }
}

/// Column-selector mode from the first ID alone: false for none or a fully parsed int,
/// true otherwise.
fn has_named_columns(column_ids: &[Bytes]) -> bool {
    if column_ids.is_empty() {
        return false;
    }
    let (r, used) = from_chars_i32(&column_ids[0]);
    r.is_err() || used != column_ids[0].len()
}

// ---- validate_table ----------------------------------------------------------

/// An OADate as "MM/DD/YYYY HH:MM:SS".
fn format_datetime(dt: f64) -> Bytes {
    let (y, m, d) = datetime::decode_date(dt);
    let (h, mi, s) = datetime::decode_time(dt);
    format!("{:02}/{:02}/{:04} {:02}:{:02}:{:02}", m, d, y, h, mi, s).into_bytes()
}

/// dx_min (smallest positive interval, 1e30 if none) and the index of the first decrease
/// (-1 if none); true when non-decreasing.
fn check_monotonicity(x: &[f64], dx_min: &mut f64, bad_idx: &mut i32) -> bool {
    *dx_min = 1.0e30;
    *bad_idx = -1;
    for i in 1..x.len() {
        let dx = x[i] - x[i - 1];
        if dx < 0.0 {
            *bad_idx = i as i32;
            return false;
        }
        if dx > 0.0 && dx < *dx_min {
            *dx_min = dx;
        }
    }
    true
}

fn sequence_error(is_ts: bool, seq_err: i32, id: &[u8], at_x: f64) -> Bytes {
    if is_ts {
        let mut detail = b"at ".to_vec();
        detail.extend_from_slice(&format_datetime(at_x));
        format_error_detail(seq_err, id, &detail)
    } else {
        format_error(seq_err, id)
    }
}

/// Validate a table: emptiness, then either the file-backed boundary checks or the
/// in-memory monotonicity / finiteness / y-size / duplicate checks. Sets `tbl.dx_min`.
pub fn validate_table(tbl: &mut Table) -> TableValidation {
    let mut result = TableValidation::default();
    let is_ts = tbl.type_ == TableType::Timeseries;
    let seq_err = if is_ts { ERR_TIMESERIES_SEQUENCE } else { ERR_CURVE_SEQUENCE };

    // ---- Empty check ----
    if tbl.is_empty() {
        result.valid = false;
        result.errors.push(format_error(ERR_TIMESERIES_EMPTY, &tbl.id));
        return result;
    }

    // ---- File-backed checks ----
    if tbl.is_file_based {
        if tbl.file_handle.is_none() {
            result.valid = false;
            result.errors.push(format_error(ERR_TABLE_FILE_OPEN, &tbl.id));
        }

        if tbl.first_boundary.is_empty() {
            result.valid = false;
            result.errors.push(format_error(ERR_TABLE_FILE_READ, &tbl.id));
        }

        // Boundary monotonicity
        let mut bdx_min = 1.0e30;
        let mut bad_idx: i32 = -1;
        if !tbl.first_boundary.is_empty() {
            check_monotonicity(&tbl.first_boundary.x, &mut bdx_min, &mut bad_idx);
            if bad_idx >= 0 {
                result.valid = false;
                result.errors.push(sequence_error(is_ts, seq_err, &tbl.id, tbl.first_boundary.x[bad_idx as usize]));
            }
            if bdx_min < tbl.dx_min || tbl.dx_min == 0.0 {
                tbl.dx_min = bdx_min;
            }
        }
        if !tbl.last_boundary.is_empty() {
            check_monotonicity(&tbl.last_boundary.x, &mut bdx_min, &mut bad_idx);
            if bad_idx >= 0 {
                result.valid = false;
                result.errors.push(sequence_error(is_ts, seq_err, &tbl.id, tbl.last_boundary.x[bad_idx as usize]));
            }
            if bdx_min < tbl.dx_min || tbl.dx_min == 0.0 {
                tbl.dx_min = bdx_min;
            }
        }

        // Cross-boundary ordering
        if !tbl.first_boundary.is_empty()
            && !tbl.last_boundary.is_empty()
            && tbl.first_boundary.x[tbl.first_boundary.x.len() - 1] > tbl.last_boundary.x[0]
        {
            result.warnings.push(format_warning(WARN_BOUNDARY_OVERLAP, &tbl.id));
        }

        // Boundary y-size consistency
        let expected = tbl.first_boundary.num_rows().wrapping_mul(tbl.first_boundary.num_cols as usize);
        if !tbl.first_boundary.is_empty() && tbl.first_boundary.y.len() != expected {
            result.valid = false;
            result.errors.push(format_error(ERR_TABLE_COL_MISMATCH, &tbl.id));
        }

        return result;
    }

    // ---- In-memory checks ----

    // Monotonicity + dx_min (matches legacy table_validate)
    {
        let mut dx_min_val = 1.0e30;
        let mut bad_idx: i32 = -1;
        check_monotonicity(&tbl.x, &mut dx_min_val, &mut bad_idx);
        tbl.dx_min = dx_min_val;

        if bad_idx >= 0 {
            result.valid = false;
            result.errors.push(sequence_error(is_ts, seq_err, &tbl.id, tbl.x[bad_idx as usize]));
        }
    }

    // NaN/Inf check
    for i in 0..tbl.x.len() {
        if !tbl.x[i].is_finite() {
            result.valid = false;
            result.errors.push(format_error(ERR_TIMESERIES_NAN, &tbl.id));
            break;
        }
    }

    // Flat y size consistency
    let expected_y = tbl.x.len().wrapping_mul(tbl.num_cols as usize);
    if tbl.y.len() != expected_y {
        result.valid = false;
        result.errors.push(format_error(ERR_TABLE_COL_MISMATCH, &tbl.id));
    }

    // Duplicate x values (warning only)
    for i in 1..tbl.x.len() {
        if tbl.x[i] == tbl.x[i - 1] {
            result.warnings.push(format_warning(WARN_TIMESERIES_DUPLICATE_X, &tbl.id));
            break;
        }
    }

    result
}

// ---- table_open_file ---------------------------------------------------------
// Upstream's note (2026-08-17): this streaming, sliding-cache reader and its
// column_ids / column_map / table_lookup_column / table_step_column companions have no
// callers; the multi-column series feature uses input::multi_column_series_file instead.

/// Open `tbl.file_path` (its authored token) as a paged external series: detect CSV
/// (a header with a delimiter and a letter) or SWMM .dat (columns counted from the first
/// data line), count rows, build the sparse offset index, and read the first and last
/// `boundary_rows` rows. False if there is no path or it does not open.
/// `boundary_rows` is 128 upstream by default (`TABLE_OPEN_FILE_DEFAULT_BOUNDARY_ROWS`).
pub fn table_open_file(tbl: &mut Table, boundary_rows: usize) -> bool {
    if tbl.file_path.is_empty() {
        return false;
    }

    let fp = match CFile::open(&tbl.file_path.original) {
        Some(f) => f,
        None => return false,
    };

    let fp = tbl.file_handle.insert(fp);
    tbl.is_file_based = true;

    // Skip UTF-8 BOM if present
    let bom = fp.fread(3);
    if bom.len() == 3 {
        if !(bom[0] == 0xEF && bom[1] == 0xBB && bom[2] == 0xBF) {
            fp.fseek(0);
        }
    } else {
        fp.fseek(0);
    }

    // Read first non-comment, non-empty line to detect format
    let mut first_line = Bytes::new(); // char first_line[2048] = {}
    let mut pre_data_pos = fp.ftell();
    let mut is_csv = false;

    while fp.fgets(2048, &mut first_line) {
        if first_line[0] == b';' || first_line[0] == b'#' || first_line[0] == b'\n' || first_line[0] == b'\r' {
            pre_data_pos = fp.ftell();
            continue;
        }
        break;
    }

    // Detect CSV with header vs SWMM .dat
    if looks_like_csv_header(&first_line) {
        is_csv = true;
        let headers = split_csv_line(&first_line);

        tbl.column_ids.clear();
        tbl.column_map.clear();
        for header in headers.iter().skip(1) {
            let col_idx = tbl.column_ids.len() as i32;
            tbl.column_ids.push(header.clone());
            tbl.column_map.insert(header.clone(), col_idx);
        }
        tbl.data_start_offset = fp.ftell();
    } else {
        // SWMM .dat format — detect column count from first data line
        let mut xv = 0.0;
        let mut vals = Vec::new();
        let mut ncols = parse_dat_line(&first_line, &mut xv, &mut vals);
        if ncols == 0 {
            ncols = 1;
        }

        tbl.column_ids.clear();
        tbl.column_map.clear();
        for c in 0..ncols {
            let cid = (c + 1).to_string().into_bytes();
            tbl.column_ids.push(cid.clone());
            tbl.column_map.insert(cid, c);
        }
        tbl.data_start_offset = pre_data_pos;
    }

    let mut ncols = tbl.column_ids.len() as i32;
    if ncols == 0 {
        ncols = 1;
    }
    tbl.num_cols = ncols;

    // Count total data rows and build sparse byte-offset index
    tbl.total_rows = count_and_index_rows(fp, tbl.data_start_offset, &mut tbl.row_offsets);

    // ---- Read first boundary ----
    fp.fseek(tbl.data_start_offset);
    let first_n = cmin(boundary_rows, tbl.total_rows);
    tbl.first_boundary.num_cols = ncols;
    read_rows_flat(fp, is_csv, ncols, first_n, &mut tbl.first_boundary.x, &mut tbl.first_boundary.y);

    // ---- Read last boundary ----
    if tbl.total_rows > boundary_rows {
        // Use row-offset index for O(1) seek to the last boundary region
        let skip = tbl.total_rows - boundary_rows;
        let chunk = skip / Table::INDEX_STRIDE;
        if chunk < tbl.row_offsets.len() {
            fp.fseek(tbl.row_offsets[chunk]);
            // Skip remaining rows within the chunk
            let remaining = skip - chunk * Table::INDEX_STRIDE;
            skip_data_rows(fp, remaining);
        } else {
            fp.fseek(tbl.data_start_offset);
            skip_data_rows(fp, skip);
        }
        tbl.last_boundary.num_cols = ncols;
        read_rows_flat(fp, is_csv, ncols, boundary_rows, &mut tbl.last_boundary.x, &mut tbl.last_boundary.y);
    } else {
        tbl.last_boundary = tbl.first_boundary.clone();
    }

    // Set x_min / x_max from boundaries
    if !tbl.first_boundary.is_empty() {
        tbl.x_min = tbl.first_boundary.x[0];
    }
    if !tbl.last_boundary.is_empty() {
        tbl.x_max = tbl.last_boundary.x[tbl.last_boundary.x.len() - 1];
    }

    // Rewind to data start for future cache loads
    let fp = tbl.file_handle.as_mut().expect("opened above");
    fp.fseek(tbl.data_start_offset);

    true
}

// ---- table_load_cache --------------------------------------------------------

/// Fill the sliding cache with up to `num_cache_rows` rows starting at data row
/// `start_row` (sparse index, then a linear skip). CSV-vs-.dat follows
/// `has_named_columns(column_ids)`. Returns the rows read; 0 if no file is open.
pub fn table_load_cache(tbl: &mut Table, start_row: usize) -> usize {
    if tbl.file_handle.is_none() || !tbl.is_file_based {
        return 0;
    }

    let is_csv = has_named_columns(&tbl.column_ids);

    let mut ncols = tbl.num_cols;
    if ncols == 0 {
        ncols = 1;
    }

    let fp = tbl.file_handle.as_mut().expect("checked above");

    // Seek to start_row using the row-offset index for O(1) positioning
    {
        let chunk = start_row / Table::INDEX_STRIDE;
        if !tbl.row_offsets.is_empty() && chunk < tbl.row_offsets.len() {
            fp.fseek(tbl.row_offsets[chunk]);
            // Skip remaining rows within the chunk
            let remaining = start_row - chunk * Table::INDEX_STRIDE;
            skip_data_rows(fp, remaining);
        } else {
            // Fallback: linear scan from data start
            fp.fseek(tbl.data_start_offset);
            skip_data_rows(fp, start_row);
        }
    }

    tbl.cache.num_cols = ncols;
    let rows = read_rows_flat(fp, is_csv, ncols, tbl.num_cache_rows, &mut tbl.cache.x, &mut tbl.cache.y);
    tbl.cache.file_row_start = start_row;

    rows
}

// ---- multicolumn lookup ------------------------------------------------------

/// Make the data for a column lookup available: for a file-backed table, load or slide
/// the cache so it covers `x_query` (a cache load centered on the row `x_query` would
/// occupy if the rows were evenly spaced). Returns (rows, stride = columns, cursor row
/// relative to the data), or None if there is no data. The data are `tbl.cache.x/y`
/// when file-backed, else `tbl.x/y`.
fn prepare_lookup(tbl: &mut Table, x_query: f64) -> Option<(i32, i32, i32)> {
    if tbl.is_file_based {
        if tbl.cache.is_empty() || tbl.total_rows == 0 {
            table_load_cache(tbl, 0);
        }
        if tbl.cache.is_empty() {
            return None;
        }

        // Check if x_query is outside cache window
        if x_query < tbl.cache.x[0] || x_query > tbl.cache.x[tbl.cache.x.len() - 1] {
            let mut target_row = 0usize;
            if tbl.x_max > tbl.x_min && tbl.total_rows > 1 {
                let frac = cclamp((x_query - tbl.x_min) / (tbl.x_max - tbl.x_min), 0.0, 1.0);
                target_row = (frac * (tbl.total_rows - 1) as f64) as usize;
            }
            let half = tbl.num_cache_rows / 2;
            let start = target_row.saturating_sub(half); // (target_row > half) ? target_row - half : 0
            table_load_cache(tbl, start);
            if tbl.cache.is_empty() {
                return None;
            }
        }

        let n = tbl.cache.num_rows() as i32;
        let stride = tbl.cache.num_cols;

        let row_offset = tbl.cache.file_row_start as i32;
        let cursor_idx = cclamp(tbl.cursor.index - row_offset, 0, n - 1);
        return Some((n, stride, cursor_idx));
    }

    // In-memory
    let n = tbl.x.len() as i32;
    if n == 0 {
        return None;
    }

    let stride = tbl.num_cols;
    let cursor_idx = cclamp(tbl.cursor.index, 0, n - 1);
    Some((n, stride, cursor_idx))
}

/// Linear interpolation in column `col_idx` (0 if out of range); the end values beyond
/// both ends. Moves the cursor (an absolute file row when file-backed).
pub fn table_lookup_column(tbl: &mut Table, col_idx: i32, x_query: f64) -> f64 {
    let (n, stride, mut idx) = match prepare_lookup(tbl, x_query) {
        Some(t) => t,
        None => return 0.0,
    };
    let mut col_idx = col_idx;
    let fb = tbl.is_file_based;
    let frs = tbl.cache.file_row_start as i32;
    let (x_data, y_base) = if fb { (&tbl.cache.x, &tbl.cache.y) } else { (&tbl.x, &tbl.y) };
    let su = stride as usize;

    if col_idx < 0 || col_idx >= stride {
        col_idx = 0;
    }
    let cu = col_idx as usize;
    if n == 1 {
        return y_base[cu];
    }

    // Clamp below
    if x_query <= x_data[0] {
        tbl.cursor.index = if fb { frs } else { 0 };
        tbl.cursor.direction = 1;
        return y_base[cu];
    }

    // Clamp above
    if x_query >= x_data[(n - 1) as usize] {
        tbl.cursor.index = if fb { frs + n - 1 } else { n - 1 };
        tbl.cursor.direction = -1;
        return y_base[(n - 1) as usize * su + cu];
    }

    idx = cclamp(idx, 0, n - 2);

    while idx < n - 1 && x_data[(idx + 1) as usize] < x_query {
        idx += 1;
        tbl.cursor.direction = 1;
    }
    while idx > 0 && x_data[idx as usize] > x_query {
        idx -= 1;
        tbl.cursor.direction = -1;
    }

    tbl.cursor.index = if fb { frs + idx } else { idx };

    let dx = x_data[(idx + 1) as usize] - x_data[idx as usize];
    let y0 = y_base[idx as usize * su + cu];
    let y1 = y_base[(idx + 1) as usize * su + cu];
    if dx <= 0.0 {
        return y0;
    }
    let t = (x_query - x_data[idx as usize]) / dx;
    y0 + t * (y1 - y0)
}

/// Step (zero-order hold) lookup in column `col_idx`: 0 before the first row, the last
/// value from the last row on.
pub fn table_step_column(tbl: &mut Table, col_idx: i32, x_query: f64) -> f64 {
    let (n, stride, mut idx) = match prepare_lookup(tbl, x_query) {
        Some(t) => t,
        None => return 0.0,
    };
    let mut col_idx = col_idx;
    let fb = tbl.is_file_based;
    let frs = tbl.cache.file_row_start as i32;
    let (x_data, y_base) = if fb { (&tbl.cache.x, &tbl.cache.y) } else { (&tbl.x, &tbl.y) };
    let su = stride as usize;

    if col_idx < 0 || col_idx >= stride {
        col_idx = 0;
    }
    let cu = col_idx as usize;
    if n == 1 {
        return y_base[cu];
    }

    if x_query < x_data[0] {
        tbl.cursor.index = if fb { frs } else { 0 };
        tbl.cursor.direction = 1;
        return 0.0;
    }

    if x_query >= x_data[(n - 1) as usize] {
        tbl.cursor.index = if fb { frs + n - 1 } else { n - 1 };
        tbl.cursor.direction = -1;
        return y_base[(n - 1) as usize * su + cu];
    }

    idx = cclamp(idx, 0, n - 2);

    while idx < n - 1 && x_data[(idx + 1) as usize] <= x_query {
        idx += 1;
        tbl.cursor.direction = 1;
    }
    while idx > 0 && x_data[idx as usize] > x_query {
        idx -= 1;
        tbl.cursor.direction = -1;
    }

    tbl.cursor.index = if fb { frs + idx } else { idx };

    y_base[idx as usize * su + cu]
}
