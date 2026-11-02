// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/MultiColumnSeriesFile.hpp and
// src/engine/input/MultiColumnSeriesFile.cpp (HydroCouple/openswmm.engine, branch
// swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the Apache
// License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The parse-once cache for multi-column external series files (CSV / TSV / PCSWMM TSF).
// Values are stored raw (file units); a missing or unparseable cell is a quiet NaN.
//
// Differences in form, never in result:
//   * The file is read whole and split as std::getline splits it (on '\n', the '\n'
//     dropped, a final empty piece not returned). A read error ends the input, as it
//     ends getline; a directory opens and reads nothing, as with libstdc++.
//   * Cells go through C functions upstream calls (sscanf, strtod), which see a C
//     string: the bytes before the first NUL. Both are emulated with glibc's C-locale
//     behavior; `strtod` here is glibc's, including hex floats and NaN payloads.
//   * The row sort is libstdc++'s std::sort (introsort, GCC 13), ported step for step,
//     because rows with equal datetimes come out in the order that algorithm leaves them.
//   * `get_or_parse` returns Option<&ParsedSeriesFile> for the pointer-or-nullptr result.

use crate::core::datetime;
use crate::core::file_io::path_from_bytes;
use crate::core::options::Bytes;
use crate::data::table::{sscanf, ScanArg};
use std::collections::HashMap;
use std::io::Read;
use std::sync::atomic::{AtomicI64, Ordering};

/// One fully parsed multi-column series file (raw file-unit values). Rows are sorted
/// ascending by datetime; a missing or unparseable cell is a quiet NaN so all columns
/// stay aligned to `dates`.
#[derive(Debug, Clone, Default)]
pub struct ParsedSeriesFile {
    /// Header names; index 0 is the time-column label and never matches.
    pub headers: Vec<Bytes>,
    /// Row datetimes (SWMM OADate, whole file, ascending).
    pub dates: Vec<f64>,
    /// columns[i] (i >= 1) aligned to `dates`; columns[0] is unused and empty.
    pub columns: Vec<Vec<f64>>,
    /// Whole-file first / last datetime (0 when the file has no data rows).
    pub first_date: f64,
    pub last_date: f64,
    /// Per-column statistics over the whole file (index-aligned with `headers`).
    pub col_first_date: Vec<f64>,
    pub col_last_date: Vec<f64>,
    /// Cells with value > 0.
    pub col_periods_precip: Vec<i64>,
    /// Value-parse failures.
    pub col_unparsed_cells: Vec<i64>,
    /// Data rows whose datetime cell failed to parse (skipped entirely).
    pub unparsed_rows: i64,
}

impl ParsedSeriesFile {
    /// The first non-time header equal to `name` under C-locale case folding, or -1.
    pub fn find_column(&self, name: &[u8]) -> i32 {
        // Column 0 is the time stamp, so a match there is as wrong as no match.
        for i in 1..self.headers.len() {
            if iequals(&self.headers[i], name) {
                return i as i32;
            }
        }
        -1
    }

    /// 1 if there is at least one data column, else -1.
    pub fn first_data_column(&self) -> i32 {
        if self.headers.len() > 1 {
            1
        } else {
            -1
        }
    }
}

/// Why a parse failed (or that it did not).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum SeriesFileStatus {
    /// parsed (possibly zero data rows)
    Ok = 0,
    /// file could not be opened
    OpenFailed = 1,
    /// no header / no data columns
    FormatFailed = 2,
}

static G_PARSE_COUNT_TOTAL: AtomicI64 = AtomicI64::new(0);

// ---------------------------------------------------------------------------
// C library pieces
// ---------------------------------------------------------------------------

/// `std::tolower` / `std::toupper` in the C locale.
#[inline]
fn c_tolower(c: u8) -> u8 {
    c.to_ascii_lowercase()
}
#[inline]
fn c_toupper(c: u8) -> u8 {
    c.to_ascii_uppercase()
}

/// `std::strtod`: the crate's glibc-checked emulation (core::cstdlib).
pub use crate::core::cstdlib::strtod;

/// libstdc++ (GCC 13) `std::sort(first, last, comp)` on a slice: introsort with a
/// depth limit of 2*floor(log2 n), heap sort past it, then a final insertion sort.
pub(crate) fn std_sort<T: Copy>(v: &mut [T], comp: &impl Fn(&T, &T) -> bool) {
    let n = v.len();
    if n != 0 {
        let lg = (usize::BITS - 1 - n.leading_zeros()) as i64;
        introsort_loop(v, 0, n, lg * 2, comp);
        final_insertion_sort(v, 0, n, comp);
    }
}

const S_THRESHOLD: usize = 16;

fn introsort_loop<T: Copy>(v: &mut [T], first: usize, mut last: usize, mut depth_limit: i64, comp: &impl Fn(&T, &T) -> bool) {
    while last - first > S_THRESHOLD {
        if depth_limit == 0 {
            partial_sort(v, first, last, last, comp);
            return;
        }
        depth_limit -= 1;
        let cut = unguarded_partition_pivot(v, first, last, comp);
        introsort_loop(v, cut, last, depth_limit, comp);
        last = cut;
    }
}

fn move_median_to_first<T: Copy>(v: &mut [T], result: usize, a: usize, b: usize, c: usize, comp: &impl Fn(&T, &T) -> bool) {
    if comp(&v[a], &v[b]) {
        if comp(&v[b], &v[c]) {
            v.swap(result, b);
        } else if comp(&v[a], &v[c]) {
            v.swap(result, c);
        } else {
            v.swap(result, a);
        }
    } else if comp(&v[a], &v[c]) {
        v.swap(result, a);
    } else if comp(&v[b], &v[c]) {
        v.swap(result, c);
    } else {
        v.swap(result, b);
    }
}

fn unguarded_partition<T: Copy>(v: &mut [T], mut first: usize, mut last: usize, pivot: usize, comp: &impl Fn(&T, &T) -> bool) -> usize {
    loop {
        while comp(&v[first], &v[pivot]) {
            first += 1;
        }
        last -= 1;
        while comp(&v[pivot], &v[last]) {
            last -= 1;
        }
        if !(first < last) {
            return first;
        }
        v.swap(first, last);
        first += 1;
    }
}

fn unguarded_partition_pivot<T: Copy>(v: &mut [T], first: usize, last: usize, comp: &impl Fn(&T, &T) -> bool) -> usize {
    let mid = first + (last - first) / 2;
    move_median_to_first(v, first, first + 1, mid, last - 1, comp);
    unguarded_partition(v, first + 1, last, first, comp)
}

fn unguarded_linear_insert<T: Copy>(v: &mut [T], mut last: usize, comp: &impl Fn(&T, &T) -> bool) {
    let val = v[last];
    let mut next = last - 1;
    while comp(&val, &v[next]) {
        v[last] = v[next];
        last = next;
        next -= 1;
    }
    v[last] = val;
}

fn insertion_sort<T: Copy>(v: &mut [T], first: usize, last: usize, comp: &impl Fn(&T, &T) -> bool) {
    if first == last {
        return;
    }
    for i in first + 1..last {
        if comp(&v[i], &v[first]) {
            let val = v[i];
            v.copy_within(first..i, first + 1);
            v[first] = val;
        } else {
            unguarded_linear_insert(v, i, comp);
        }
    }
}

fn unguarded_insertion_sort<T: Copy>(v: &mut [T], first: usize, last: usize, comp: &impl Fn(&T, &T) -> bool) {
    for i in first..last {
        unguarded_linear_insert(v, i, comp);
    }
}

fn final_insertion_sort<T: Copy>(v: &mut [T], first: usize, last: usize, comp: &impl Fn(&T, &T) -> bool) {
    if last - first > S_THRESHOLD {
        insertion_sort(v, first, first + S_THRESHOLD, comp);
        unguarded_insertion_sort(v, first + S_THRESHOLD, last, comp);
    } else {
        insertion_sort(v, first, last, comp);
    }
}

// Heap pieces, on the subrange starting at `first` with signed distances as upstream.
fn push_heap<T: Copy>(v: &mut [T], first: usize, mut hole: isize, top: isize, value: T, comp: &impl Fn(&T, &T) -> bool) {
    let mut parent = (hole - 1) / 2;
    while hole > top && comp(&v[first + parent as usize], &value) {
        v[first + hole as usize] = v[first + parent as usize];
        hole = parent;
        parent = (hole - 1) / 2;
    }
    v[first + hole as usize] = value;
}

fn adjust_heap<T: Copy>(v: &mut [T], first: usize, mut hole: isize, len: isize, value: T, comp: &impl Fn(&T, &T) -> bool) {
    let top = hole;
    let mut second = hole;
    while second < (len - 1) / 2 {
        second = 2 * (second + 1);
        if comp(&v[first + second as usize], &v[first + (second - 1) as usize]) {
            second -= 1;
        }
        v[first + hole as usize] = v[first + second as usize];
        hole = second;
    }
    if (len & 1) == 0 && second == (len - 2) / 2 {
        second = 2 * (second + 1);
        v[first + hole as usize] = v[first + (second - 1) as usize];
        hole = second - 1;
    }
    push_heap(v, first, hole, top, value, comp);
}

fn make_heap<T: Copy>(v: &mut [T], first: usize, last: usize, comp: &impl Fn(&T, &T) -> bool) {
    if last - first < 2 {
        return;
    }
    let len = (last - first) as isize;
    let mut parent = (len - 2) / 2;
    loop {
        let value = v[first + parent as usize];
        adjust_heap(v, first, parent, len, value, comp);
        if parent == 0 {
            return;
        }
        parent -= 1;
    }
}

fn pop_heap<T: Copy>(v: &mut [T], first: usize, last: usize, result: usize, comp: &impl Fn(&T, &T) -> bool) {
    let value = v[result];
    v[result] = v[first];
    adjust_heap(v, first, 0, (last - first) as isize, value, comp);
}

fn heap_select<T: Copy>(v: &mut [T], first: usize, middle: usize, last: usize, comp: &impl Fn(&T, &T) -> bool) {
    make_heap(v, first, middle, comp);
    let mut i = middle;
    while i < last {
        if comp(&v[i], &v[first]) {
            pop_heap(v, first, middle, i, comp);
        }
        i += 1;
    }
}

fn sort_heap<T: Copy>(v: &mut [T], first: usize, mut last: usize, comp: &impl Fn(&T, &T) -> bool) {
    while last - first > 1 {
        last -= 1;
        pop_heap(v, first, last, last, comp);
    }
}

fn partial_sort<T: Copy>(v: &mut [T], first: usize, middle: usize, last: usize, comp: &impl Fn(&T, &T) -> bool) {
    heap_select(v, first, middle, last, comp);
    sort_heap(v, first, middle, comp);
}

// ---------------------------------------------------------------------------
// Internal helpers (anonymous namespace)
// ---------------------------------------------------------------------------

/// A cell without ASCII edge white space (space, tab, CR, LF) and one surrounding pair of
/// double quotes.
fn cell_trim(s: &[u8]) -> Bytes {
    let (mut a, mut b) = (0usize, s.len());
    let space = |c: u8| c == b' ' || c == b'\t' || c == b'\r' || c == b'\n';
    while a < b && space(s[a]) {
        a += 1;
    }
    while b > a && space(s[b - 1]) {
        b -= 1;
    }
    if b - a >= 2 && s[a] == b'"' && s[b - 1] == b'"' {
        a += 1;
        b -= 1;
    }
    s[a..b].to_vec()
}

/// Split one physical line on unquoted delimiter bytes (a double quote toggles grouping
/// and stays in the cell until the per-cell trim).
fn cell_split(line: &[u8], delim: u8) -> Vec<Bytes> {
    let mut out = Vec::new();
    let mut cur = Bytes::new();
    let mut in_quotes = false;
    for &c in line {
        if c == b'"' {
            in_quotes = !in_quotes;
            cur.push(c);
        } else if c == delim && !in_quotes {
            out.push(cell_trim(&cur));
            cur.clear();
        } else {
            cur.push(c);
        }
    }
    out.push(cell_trim(&cur));
    out
}

/// Delimiter bytes outside toggled double-quote groups.
fn count_unquoted(line: &[u8], delim: u8) -> i64 {
    let mut n = 0i64;
    let mut in_quotes = false;
    for &c in line {
        if c == b'"' {
            in_quotes = !in_quotes;
        } else if c == delim && !in_quotes {
            n += 1;
        }
    }
    n
}

/// Equal-length byte strings equal under C-locale lowercase mapping.
fn iequals(a: &[u8], b: &[u8]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    for i in 0..a.len() {
        if c_tolower(a[i]) != c_tolower(b[i]) {
            return false;
        }
    }
    true
}

/// 1 for a trailing "AM", 2 for "PM" (after ASCII edge white space), else 0.
fn trailing_meridiem(cell: &[u8]) -> i32 {
    let e = match cell.iter().rposition(|&c| !(c == b' ' || c == b'\t' || c == b'\r' || c == b'\n')) {
        Some(e) => e,
        None => return 0,
    };
    if e < 1 {
        return 0;
    }
    let m1 = c_toupper(cell[e]);
    let m0 = c_toupper(cell[e - 1]);
    if m1 != b'M' || (m0 != b'A' && m0 != b'P') {
        return 0;
    }
    if m0 == b'A' {
        1
    } else {
        2
    }
}

/// Remove a first-line UTF-8 BOM and one trailing CR.
fn normalize_line(line: &mut Bytes, first_line: bool) {
    if first_line && line.len() >= 3 && line[0] == 0xEF && line[1] == 0xBB && line[2] == 0xBF {
        line.drain(0..3);
    }
    if line.last() == Some(&b'\r') {
        line.pop();
    }
}

/// Empty, all space/tab, or a ';' / '#' comment after leading space/tab.
fn is_blank_or_comment(line: &[u8]) -> bool {
    for &c in line {
        if c == b' ' || c == b'\t' {
            continue;
        }
        return c == b';' || c == b'#';
    }
    true // all-whitespace / empty
}

/// `line` begins with `prefix` under C-locale case folding.
fn istarts_with(line: &[u8], prefix: &[u8]) -> bool {
    if line.len() < prefix.len() {
        return false;
    }
    iequals(&line[..prefix.len()], prefix)
}

/// An `std::ifstream` opened on `abs_path` and read line by line with `std::getline`.
struct LineStream {
    data: Bytes,
    pos: usize,
}

impl LineStream {
    /// None where `is_open()` is false.
    fn open(abs_path: &[u8]) -> Option<LineStream> {
        let mut f = std::fs::File::open(path_from_bytes(abs_path)).ok()?;
        let mut data = Bytes::new();
        let mut chunk = [0u8; 65536];
        loop {
            match f.read(&mut chunk) {
                Ok(0) => break,
                Ok(k) => data.extend_from_slice(&chunk[..k]),
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => continue,
                Err(_) => break,
            }
        }
        Some(LineStream { data, pos: 0 })
    }

    /// `std::getline(in, line)`: false once nothing is left.
    fn getline(&mut self, line: &mut Bytes) -> bool {
        if self.pos >= self.data.len() {
            return false;
        }
        let rest = &self.data[self.pos..];
        match rest.iter().position(|&c| c == b'\n') {
            Some(k) => {
                *line = rest[..k].to_vec();
                self.pos += k + 1;
            }
            None => {
                *line = rest.to_vec();
                self.pos = self.data.len();
            }
        }
        true
    }
}

// ---------------------------------------------------------------------------
// Datetime parsing (relocated csv_parse_datetime + 12-hour AM/PM)
// ---------------------------------------------------------------------------

/// Parse "YYYY-MM-DD[ T]HH:MM[:SS]" or "MM/DD/YYYY[ HH:MM[:SS]][ AM|PM]" into SWMM decimal
/// days. `out` is written when the three date fields scan; false otherwise.
pub fn parse_series_datetime(cell: &[u8], out: &mut f64) -> bool {
    use ScanArg::Int;
    let (mut y, mut mo, mut d, mut h, mut mi, mut s) = (0i32, 0i32, 0i32, 0i32, 0i32, 0i32);

    // ISO-8601: YYYY-MM-DD[ T]HH:MM[:SS]
    let mut n = sscanf(
        cell,
        b"%d-%d-%d%*[ T]%d:%d:%d",
        &mut [Int(&mut y), Int(&mut mo), Int(&mut d), Int(&mut h), Int(&mut mi), Int(&mut s)],
    );
    if n >= 3 {
        if n < 4 {
            h = 0;
            mi = 0;
            s = 0;
        } else if n < 6 {
            s = 0;
        }
        *out = datetime::encode_date(y, mo, d) + datetime::encode_time(h, mi, s);
        return true;
    }

    // US: MM/DD/YYYY[ HH:MM[:SS]], optionally 12-hour with a trailing AM/PM token.
    n = sscanf(
        cell,
        b"%d/%d/%d %d:%d:%d",
        &mut [Int(&mut mo), Int(&mut d), Int(&mut y), Int(&mut h), Int(&mut mi), Int(&mut s)],
    );
    if n >= 3 {
        if n < 4 {
            h = 0;
            mi = 0;
            s = 0;
        } else if n < 6 {
            s = 0;
        }
        if n >= 4 {
            let meridiem = trailing_meridiem(cell);
            // Some TSF exporters append AM/PM to a 24-hour clock (00 AM, 13..23 PM).
            // Those hours are already absolute; do not add 12.
            if meridiem != 0 && (1..=12).contains(&h) {
                if h == 12 {
                    h = 0; // 12:xx AM -> 00:xx, 12:xx PM -> 12:xx
                }
                if meridiem == 2 {
                    h += 12; // PM
                }
            }
        }
        *out = datetime::encode_date(y, mo, d) + datetime::encode_time(h, mi, s);
        return true;
    }
    false
}

// ---------------------------------------------------------------------------
// Sniff
// ---------------------------------------------------------------------------

/// True when the first content row starts with "IDs:", has a comma, or has tabs and a
/// first cell that is not a datetime.
pub fn looks_like_multicolumn_series_file(abs_path: &[u8]) -> bool {
    let mut input = match LineStream::open(abs_path) {
        Some(s) => s,
        None => return false,
    };

    let mut line = Bytes::new();
    let mut first = true;
    while input.getline(&mut line) {
        normalize_line(&mut line, first);
        first = false;
        if is_blank_or_comment(&line) {
            continue;
        }

        if istarts_with(&line, b"IDs:") {
            return true; // PCSWMM TSF
        }
        if line.contains(&b',') {
            return true; // comma CSV
        }
        if line.contains(&b'\t') {
            // Tab-delimited: a legacy .dat data row starts with a date; a TSV header row
            // does not.
            let cells = cell_split(&line, b'\t');
            let mut dt = 0.0;
            return !cells.is_empty() && !parse_series_datetime(&cells[0], &mut dt);
        }
        return false;
    }
    false
}

// ---------------------------------------------------------------------------
// `path:column` token split (one rule for gages AND timeseries)
// ---------------------------------------------------------------------------

/// Split "path:column" at the last colon, unless that colon is a drive letter's (index 1
/// after a letter) or the suffix holds a slash or backslash.
pub fn split_series_file_token(token: &[u8], path_out: &mut Bytes, column_out: &mut Bytes) -> bool {
    *path_out = token.to_vec();
    column_out.clear();

    let colon = match token.iter().rposition(|&c| c == b':') {
        Some(c) => c,
        None => return false,
    };

    // Windows drive letter ("C:\dir\f.csv"): not a column separator.
    if colon == 1 && token[0].is_ascii_alphabetic() {
        return false;
    }

    // A colon inside a directory name ("/data/Rich:2024/f.csv"): the text after it is
    // still path, not a column name.
    let suffix = &token[colon + 1..];
    if suffix.contains(&b'/') || suffix.contains(&b'\\') {
        return false;
    }

    *path_out = token[..colon].to_vec();
    *column_out = suffix.to_vec();
    true
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse a CSV / TSV / TSF file into `out` (reset first): headers, ascending dates and
/// raw column values. False (with a message appended to `errors`) if the file does not
/// open, has no header row, or has no data column.
pub fn parse_multicolumn_series_file(
    abs_path: &[u8],
    out: &mut ParsedSeriesFile,
    errors: &mut Vec<Bytes>,
    mut status: Option<&mut SeriesFileStatus>,
) -> bool {
    let mut fail = |status: Option<&mut SeriesFileStatus>, st: SeriesFileStatus, msg: &[u8]| -> bool {
        if let Some(s) = status {
            *s = st;
        }
        let mut m = msg.to_vec();
        m.extend_from_slice(abs_path);
        errors.push(m);
        false
    };

    *out = ParsedSeriesFile::default();
    G_PARSE_COUNT_TOTAL.fetch_add(1, Ordering::Relaxed);

    let mut input = match LineStream::open(abs_path) {
        Some(s) => s,
        None => return fail(status, SeriesFileStatus::OpenFailed, b"cannot open series file "),
    };

    // --- Header ---------------------------------------------------------
    let mut line = Bytes::new();
    let mut first = true;
    let mut have_header = false;
    while input.getline(&mut line) {
        normalize_line(&mut line, first);
        first = false;
        if is_blank_or_comment(&line) {
            continue;
        }
        have_header = true;
        break;
    }
    if !have_header {
        return fail(status, SeriesFileStatus::FormatFailed, b"no header row in series file ");
    }

    let delim: u8;
    if istarts_with(&line, b"IDs:") {
        // PCSWMM TSF: tab-delimited; line 1 tokens 1..N name the columns, line 2
        // (parameter row) and line 3 (units row) carry no data; skip both.
        delim = b'\t';
        out.headers = cell_split(&line, delim);
        let mut skip = Bytes::new();
        let mut i = 0;
        while i < 2 && input.getline(&mut skip) {
            i += 1;
        }
    } else {
        // CSV/TSV: sniff the delimiter from the header row, counting only delimiters
        // outside quotes.
        let tabs = count_unquoted(&line, b'\t');
        let commas = count_unquoted(&line, b',');
        delim = if tabs > 0 && tabs >= commas { b'\t' } else { b',' };
        out.headers = cell_split(&line, delim);
    }

    let ncols = out.headers.len();
    if ncols < 2 {
        return fail(status.take(), SeriesFileStatus::FormatFailed, b"no data columns in series file ");
    }

    out.columns.clear();
    out.columns.resize(ncols, Vec::new());
    out.col_first_date.clear();
    out.col_first_date.resize(ncols, 0.0);
    out.col_last_date.clear();
    out.col_last_date.resize(ncols, 0.0);
    out.col_periods_precip.clear();
    out.col_periods_precip.resize(ncols, 0);
    out.col_unparsed_cells.clear();
    out.col_unparsed_cells.resize(ncols, 0);

    // --- Data rows ------------------------------------------------------
    let k_nan = f64::NAN; // std::numeric_limits<double>::quiet_NaN()
    while input.getline(&mut line) {
        normalize_line(&mut line, false);
        if is_blank_or_comment(&line) {
            continue;
        }

        let cells = cell_split(&line, delim);
        let mut dt = 0.0;
        if cells.is_empty() || !parse_series_datetime(&cells[0], &mut dt) {
            out.unparsed_rows += 1;
            continue;
        }

        out.dates.push(dt);
        if out.first_date == 0.0 || dt < out.first_date {
            out.first_date = dt;
        }
        if dt > out.last_date {
            out.last_date = dt;
        }

        for i in 1..ncols {
            let mut v = k_nan;
            if i < cells.len() {
                let cell = &cells[i];
                let (parsed, used) = strtod(cell);
                if used != 0 {
                    v = parsed;
                    if out.col_first_date[i] == 0.0 || dt < out.col_first_date[i] {
                        out.col_first_date[i] = dt;
                    }
                    if dt > out.col_last_date[i] {
                        out.col_last_date[i] = dt;
                    }
                    if v > 0.0 {
                        out.col_periods_precip[i] += 1;
                    }
                } else {
                    // Includes present-but-empty cells (counted as unreadable).
                    out.col_unparsed_cells[i] += 1;
                }
            }
            // A short row (cell absent for this column) stays NaN without counting.
            out.columns[i].push(v);
        }
    }

    // --- Sort rows ascending (records may be interleaved/out of order) --
    let is_sorted = !(1..out.dates.len()).any(|k| out.dates[k] < out.dates[k - 1]);
    if !is_sorted {
        let mut order: Vec<usize> = (0..out.dates.len()).collect();
        {
            let dates = &out.dates;
            std_sort(&mut order, &|a: &usize, b: &usize| dates[*a] < dates[*b]);
        }
        let permute = |v: &mut Vec<f64>| {
            let mut tmp = vec![0.0; v.len()];
            for i in 0..order.len() {
                tmp[i] = v[order[i]];
            }
            std::mem::swap(v, &mut tmp);
        };
        permute(&mut out.dates);
        for i in 1..ncols {
            permute(&mut out.columns[i]);
        }
    }

    if let Some(s) = status {
        *s = SeriesFileStatus::Ok;
    }
    true
}

// ---------------------------------------------------------------------------
// Cache
// ---------------------------------------------------------------------------

#[derive(Debug)]
struct Entry {
    status: SeriesFileStatus,
    file: ParsedSeriesFile,
}

impl Default for Entry {
    fn default() -> Self {
        Entry { status: SeriesFileStatus::OpenFailed, file: ParsedSeriesFile::default() }
    }
}

/// Per-resolve-pass cache keyed by resolved path (exact bytes). Failures are cached too.
#[derive(Debug, Default)]
pub struct MultiColumnFileCache {
    files: HashMap<Bytes, Entry>,
    parse_count: i32,
}

impl MultiColumnFileCache {
    /// The parsed file for `abs_path` (parsing it on first request), or None for a cached
    /// failure. `errors` receives diagnostics only on the first parse of a path.
    pub fn get_or_parse(
        &mut self,
        abs_path: &[u8],
        errors: &mut Vec<Bytes>,
        status: Option<&mut SeriesFileStatus>,
    ) -> Option<&ParsedSeriesFile> {
        if !self.files.contains_key(abs_path) {
            let mut entry = Entry::default();
            parse_multicolumn_series_file(abs_path, &mut entry.file, errors, Some(&mut entry.status));
            self.parse_count += 1;
            self.files.insert(abs_path.to_vec(), entry);
        }
        let e = &self.files[abs_path];
        if let Some(s) = status {
            *s = e.status;
        }
        if e.status == SeriesFileStatus::Ok {
            Some(&e.file)
        } else {
            None
        }
    }

    /// Parse attempts made by this cache (failures count, hits do not).
    pub fn parse_count(&self) -> i32 {
        self.parse_count
    }
}

/// Calls made to `parse_multicolumn_series_file` in this process, failed opens included.
pub fn multicolumn_parse_count_total() -> i64 {
    G_PARSE_COUNT_TOTAL.load(Ordering::Relaxed)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Recorded from glibc strtod (2.39) on the reference machine.
    #[test]
    fn strtod_matches_glibc() {
        let cases: &[(&str, u64, usize)] = &[
            ("-nan", 0xfff8000000000000, 4),
            ("nan(123)", 0x7ff800000000007b, 8),
            ("nan(0b1)", 0x7ff8000000000000, 8),
            ("nan(0x1f)", 0x7ff800000000001f, 9),
            ("nan(012)", 0x7ff800000000000a, 8),
            ("nan(08)", 0x7ff8000000000000, 7),
            ("nan()", 0x7ff8000000000000, 5),
            ("nan(_)", 0x7ff8000000000000, 6),
            ("-nan(0xffffffffffffffffff)", 0xffffffffffffffff, 26),
            ("0x", 0, 1),
            ("-0x", 0x8000000000000000, 2),
            ("0xg", 0, 1),
            ("1.e5", 0x40f86a0000000000, 4),
            ("0x.8", 0x3fe0000000000000, 4),
            ("0x1.", 0x3ff0000000000000, 4),
            ("0x1p", 0x3ff0000000000000, 3),
            ("0x1p-1074", 1, 9),
            ("0x1p-1075", 0, 9),
            ("0x1.8p-1075", 1, 11),
            ("0x3p-1076", 1, 9),
            ("0x1.fffffffffffff8p1023", 0x7ff0000000000000, 23),
            (" \x0b\x0c+1", 0x3ff0000000000000, 5),
            ("+-1", 0, 0),
            ("0x00g", 0, 4),
            ("-0e", 0x8000000000000000, 2),
            ("1e99999999999999999999", 0x7ff0000000000000, 22),
            ("0x1p99999999999999999999", 0x7ff0000000000000, 24),
            (".e1", 0, 0),
            ("0x.p1", 0, 1),
        ];
        for &(s, bits, used) in cases {
            let (v, u) = strtod(s.as_bytes());
            assert_eq!((v.to_bits(), u), (bits, used), "{s}");
        }
    }
}
