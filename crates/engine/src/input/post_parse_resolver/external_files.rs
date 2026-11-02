// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/input/PostParseResolver.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The external-file loaders of PostParseResolver.cpp: FILE timeseries (legacy
// "date time value" text and the shared multi-column cache), and FILE rain gages
// (multi-column USER_CSV, the NWS/NCDC layouts, and the standard station file).
//
// Differences in form, never in result:
//   * `std::FILE*` is data::table::CFile (fgets with upstream's buffer sizes, rewind as
//     fseek 0). A path is a C string, so it ends at its first NUL, as fopen sees it.
//   * `std::sscanf`/`std::strtod` are data::table::sscanf and core::cstdlib::strtod.
//   * nws_read_file reads fields at fixed columns of its `char line[4096]`, sometimes
//     past the line's NUL (ONLINE value column); the buffer is therefore kept as a
//     persistent 4096-byte array that fgets overwrites, so such a read sees the bytes an
//     earlier, longer line left there, as in C. Bytes never written read as 0 here (stack
//     garbage in C).
//   * `static_cast<long>(float)` is cvt_i64_f32 (cvttss2si: NaN or out of range gives
//     LONG_MIN).
//   * std::sort of the station-file rows is libstdc++'s introsort
//     (multi_column_series_file::std_sort), so equal timestamps keep upstream's order.

use crate::core::context::SimulationContext;
use crate::core::cstdlib::{c_str, strtod};
use crate::core::datetime;
use crate::core::error_codes::{
    format_error, format_error_detail, format_warning_detail, ERR_RAIN_FILE_FORMAT, ERR_RAIN_FILE_OPEN,
    ERR_TABLE_FILE_OPEN, ERR_TABLE_FILE_READ, WARN_RAIN_CSV_ROWS_SKIPPED,
};
use crate::core::options::Bytes;
use crate::data::gage::{RainFileFormat, RainSource};
use crate::data::table::{sscanf, CFile, ScanArg, Table, TableType};
use crate::input::multi_column_series_file::{
    looks_like_multicolumn_series_file, split_series_file_token, std_sort, MultiColumnFileCache, SeriesFileStatus,
};

/// `fopen_utf8(path, "r")` on POSIX.
fn fopen(path: &[u8]) -> Option<CFile> {
    CFile::open(c_str(path))
}

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut v = Bytes::new();
    for p in parts {
        v.extend_from_slice(p);
    }
    v
}

/// `static_cast<long>(f)` as x86-64 cvttss2si does it.
fn cvt_i64_f32(x: f32) -> i64 {
    if x >= -9.223_372_036_854_775_808e18_f32 && x < 9.223_372_036_854_775_808e18_f32 {
        x as i64
    } else {
        i64::MIN
    }
}

// ---------------------------------------------------------------------------
// FILE timeseries
// ---------------------------------------------------------------------------

/// The tokens `next_tok` yields from a line: maximal runs of bytes other than space,
/// tab, CR, LF and comma, up to the first three.
fn line_tokens(line: &[u8]) -> Vec<&[u8]> {
    let is_sep = |c: u8| c == b' ' || c == b'\t' || c == b'\n' || c == b'\r' || c == b',';
    let mut out = Vec::new();
    let mut i = 0;
    while i < line.len() && out.len() < 3 {
        while i < line.len() && is_sep(line[i]) {
            i += 1;
        }
        if i >= line.len() {
            break;
        }
        let start = i;
        while i < line.len() && !is_sep(line[i]) {
            i += 1;
        }
        out.push(&line[start..i]);
    }
    out
}

fn parse_file_date(s: &[u8], d_out: &mut f64) -> bool {
    if !s.contains(&b'/') && !s.contains(&b'-') {
        return false;
    }
    let (mut m, mut d, mut y): (u32, u32, u32) = (0, 0, 0);
    let (mut sep1, mut sep2): (u8, u8) = (0, 0);
    if sscanf(
        s,
        b"%u%c%u%c%u",
        &mut [ScanArg::UInt(&mut m), ScanArg::Char(&mut sep1), ScanArg::UInt(&mut d), ScanArg::Char(&mut sep2), ScanArg::UInt(&mut y)],
    ) < 5
    {
        let mut mon = Bytes::new();
        if sscanf(
            s,
            b"%3[A-Za-z]%c%u%c%u",
            &mut [ScanArg::Str(&mut mon), ScanArg::Char(&mut sep1), ScanArg::UInt(&mut d), ScanArg::Char(&mut sep2), ScanArg::UInt(&mut y)],
        ) < 5
        {
            return false;
        }
        const K_MONTHS: [&[u8; 3]; 12] =
            [b"JAN", b"FEB", b"MAR", b"APR", b"MAY", b"JUN", b"JUL", b"AUG", b"SEP", b"OCT", b"NOV", b"DEC"];
        let mc = |k: usize| -> u8 { mon.get(k).copied().unwrap_or(0) };
        m = 0;
        for i in 0..12u32 {
            let km = K_MONTHS[i as usize];
            if mc(0).to_ascii_uppercase() == km[0] && mc(1).to_ascii_uppercase() == km[1] && mc(2).to_ascii_uppercase() == km[2] {
                m = i + 1;
                break;
            }
        }
        if m == 0 {
            return false;
        }
    }
    let enc = datetime::encode_date(y as i32, m as i32, d as i32);
    if enc == -(datetime::DATE_DELTA as f64) {
        return false;
    }
    *d_out = enc;
    true
}

fn parse_file_time(s: &[u8], t_out: &mut f64) -> bool {
    let (hrs, used) = strtod(s);
    if used == s.len() {
        // decimal hours
        *t_out = hrs / 24.0;
        return true;
    }
    let (mut hr, mut min, mut sec) = (0i32, 0i32, 0i32);
    if sscanf(s, b"%d:%d:%d", &mut [ScanArg::Int(&mut hr), ScanArg::Int(&mut min), ScanArg::Int(&mut sec)]) < 1 {
        return false;
    }
    if hr < 0 || min < 0 || sec < 0 {
        return false;
    }
    *t_out = datetime::encode_time(hr, min, sec);
    true
}

/// Append each FILE timeseries' rows to its in-memory x/y.
pub(super) fn load_external_timeseries_files(ctx: &mut SimulationContext, inp_dir: &[u8], file_cache: &mut MultiColumnFileCache) {
    let start_date = ctx.options.start_date;
    let errors = &mut ctx.errors;
    for t in 0..ctx.tables.tables.len() {
        let tbl = &mut ctx.tables.tables[t];
        if tbl.type_ != TableType::Timeseries {
            continue;
        }
        if tbl.file_path.is_empty() {
            continue;
        }

        let mut file_path: Bytes;
        if !tbl.file_path.absolute.is_empty() {
            file_path = tbl.file_path.absolute.clone();
        } else {
            file_path = tbl.file_path.original.clone();
            if !file_path.is_empty() && file_path[0] != b'/' && file_path[0] != b'\\' && !inp_dir.is_empty() {
                file_path = cat(&[inp_dir, b"/", &file_path]);
            }
        }

        let mut col_name = Bytes::new();
        {
            let mut path_only = Bytes::new();
            split_series_file_token(&file_path, &mut path_only, &mut col_name);
            file_path = path_only;
        }

        if !col_name.is_empty() || looks_like_multicolumn_series_file(&file_path) {
            let mut file_errors: Vec<Bytes> = Vec::new();
            let mut st = SeriesFileStatus::Ok;
            let mut ok = file_cache.get_or_parse(&file_path, &mut file_errors, Some(&mut st)).is_some();
            let mut used_path = file_path.clone();
            if !ok && st == SeriesFileStatus::OpenFailed {
                let mut verbatim = Bytes::new();
                let mut vcol = Bytes::new();
                split_series_file_token(&tbl.file_path.original, &mut verbatim, &mut vcol);
                if verbatim != file_path {
                    ok = file_cache.get_or_parse(&verbatim, &mut file_errors, Some(&mut st)).is_some();
                    used_path = verbatim;
                }
            }
            if !ok {
                errors.push(format_error(
                    if st == SeriesFileStatus::OpenFailed { ERR_TABLE_FILE_OPEN } else { ERR_TABLE_FILE_READ },
                    &tbl.id,
                ));
                continue;
            }
            let pf = file_cache.get_or_parse(&used_path, &mut file_errors, None).unwrap();
            let col = if col_name.is_empty() { pf.first_data_column() } else { pf.find_column(&col_name) };
            if col < 0 {
                errors.push(format_error_detail(
                    ERR_TABLE_FILE_READ,
                    &tbl.id,
                    &cat(&[b"column \"", &col_name, b"\" not found in ", &file_path]),
                ));
                continue;
            }
            let vals = &pf.columns[col as usize];
            tbl.x.reserve(pf.dates.len().saturating_sub(tbl.x.len()));
            tbl.y.reserve(pf.dates.len().saturating_sub(tbl.y.len()));
            for i in 0..pf.dates.len() {
                if vals[i].is_nan() {
                    continue;
                }
                tbl.x.push(pf.dates[i]);
                tbl.y.push(vals[i]);
            }
            if tbl.x.is_empty() {
                errors.push(format_error_detail(ERR_TABLE_FILE_READ, &tbl.id, &file_path));
                continue;
            }
            tbl.x.shrink_to_fit();
            tbl.y.shrink_to_fit();
            continue;
        }

        let mut fp = match fopen(&file_path) {
            Some(f) => f,
            None => match fopen(&tbl.file_path.original) {
                Some(f) => f,
                None => {
                    errors.push(format_error(ERR_TABLE_FILE_OPEN, &tbl.id));
                    continue;
                }
            },
        };

        {
            // Reserve from the file's size (capacity only).
            let rows = match std::fs::metadata(crate::core::file_io::path_from_bytes(c_str(&file_path))) {
                Ok(md) if !md.is_dir() => (md.len() as usize) / 24 + 16,
                _ => 1024,
            };
            let rows = rows.min(2_000_000);
            tbl.x.reserve(rows.saturating_sub(tbl.x.len()));
            tbl.y.reserve(rows.saturating_sub(tbl.y.len()));
        }

        let mut last_date = start_date;
        let mut line = Bytes::new();
        while fp.fgets(256, &mut line) {
            let toks = line_tokens(c_str(&line));
            if toks.is_empty() || toks[0][0] == b';' {
                continue;
            }
            let time_tok: &[u8];
            let value_tok: &[u8];
            let mut d = 0.0;
            if toks.len() >= 3 {
                if !parse_file_date(toks[0], &mut d) {
                    continue;
                }
                last_date = d;
                time_tok = toks[1];
                value_tok = toks[2];
            } else if toks.len() == 2 {
                d = last_date;
                time_tok = toks[0];
                value_tok = toks[1];
            } else {
                continue;
            }

            let mut t = 0.0;
            if !parse_file_time(time_tok, &mut t) {
                continue;
            }
            let (value, used) = strtod(value_tok);
            if used == 0 || used != value_tok.len() {
                continue;
            }
            tbl.x.push(d + t);
            tbl.y.push(value);
        }
        drop(fp);

        if tbl.x.is_empty() {
            errors.push(format_error_detail(ERR_TABLE_FILE_READ, &tbl.id, &file_path));
            continue;
        }
        tbl.x.shrink_to_fit();
        tbl.y.shrink_to_fit();
    }
}

// ---------------------------------------------------------------------------
// USER_CSV rain files
// ---------------------------------------------------------------------------

fn load_rain_file_user_csv(ctx: &mut SimulationContext, g: i32, path: &[u8], file_cache: &mut MultiColumnFileCache) {
    let ug = g as usize;
    let mut file_errors: Vec<Bytes> = Vec::new();
    let mut st = SeriesFileStatus::Ok;
    let mut used_path = path.to_vec();
    let mut ok = file_cache.get_or_parse(path, &mut file_errors, Some(&mut st)).is_some();
    if !ok && st == SeriesFileStatus::OpenFailed && ctx.gages.file_path[ug].original.as_slice() != path {
        used_path = ctx.gages.file_path[ug].original.clone();
        ok = file_cache.get_or_parse(&used_path, &mut file_errors, Some(&mut st)).is_some();
    }
    if !ok {
        ctx.errors.push(format_error(
            if st == SeriesFileStatus::OpenFailed { ERR_RAIN_FILE_OPEN } else { ERR_RAIN_FILE_FORMAT },
            path,
        ));
        return;
    }
    let pf = file_cache.get_or_parse(&used_path, &mut file_errors, None).unwrap();

    let want = ctx.gages.col_name[ug].clone();
    let col = if want.is_empty() { pf.first_data_column() } else { pf.find_column(&want) };
    if col < 0 {
        ctx.errors.push(format_error(ERR_RAIN_FILE_FORMAT, &cat(&[path, b" (column \"", &want, b"\")"])));
        return;
    }
    let uc = col as usize;

    let mut series = Table::new(ctx.gage_names.name_of(g), TableType::Timeseries);

    let vals = &pf.columns[uc];
    for i in 0..pf.dates.len() {
        let val = vals[i];
        if val.is_nan() {
            continue;
        }
        let dt = pf.dates[i];
        series.x.push(dt);
        series.y.push(val);
    }

    let unparsed_rows: i64 = pf.unparsed_rows + pf.col_unparsed_cells[uc];
    if unparsed_rows > 0 {
        let detail = cat(&[unparsed_rows.to_string().as_bytes(), b" row(s), ", path]);
        let name = ctx.gage_names.name_of(g).to_vec();
        ctx.warnings.push(format_warning_detail(WARN_RAIN_CSV_ROWS_SKIPPED, &name, &detail));
    }

    series.x.shrink_to_fit();
    series.y.shrink_to_fit();

    ctx.gages.file_first_date[ug] = pf.col_first_date[uc];
    ctx.gages.file_last_date[ug] = pf.col_last_date[uc];
    ctx.gages.file_periods_precip[ug] = pf.col_periods_precip[uc];
    ctx.gages.rain_series[ug] = series;
}

// ---------------------------------------------------------------------------
// NWS / NCDC external rain file formats (legacy rain.c)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum NwsFormat {
    None,
    Tape,
    SpaceDelimited,
    CommaDelimited,
    Online60,
    Online15,
}

#[derive(Debug, Clone, Copy)]
struct NwsFileSpec {
    fmt: NwsFormat,
    interval_sec: i32,
    time_offset: i32,
    hdr_lines: i32,
    has_station_name: bool,
    value_offset: i32,
    data_offset: i32,
}

impl Default for NwsFileSpec {
    fn default() -> Self {
        NwsFileSpec {
            fmt: NwsFormat::None,
            interval_sec: 0,
            time_offset: 0,
            hdr_lines: 0,
            has_station_name: false,
            value_offset: 0,
            data_offset: 0,
        }
    }
}

fn nws_interval(elem_type: &[u8]) -> i32 {
    if elem_type == b"HPCP" {
        return 3600;
    }
    if elem_type == b"QPCP" {
        return 900;
    }
    if elem_type == b"QGAG" {
        return 900;
    }
    0
}

fn find_sub(hay: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.len() > hay.len() {
        return None;
    }
    (0..=hay.len() - needle.len()).find(|&i| &hay[i..i + needle.len()] == needle)
}

fn nws_detect_online(f: &mut CFile, spec: &mut NwsFileSpec) -> bool {
    let mut buf = Bytes::new();
    f.fseek(0);
    if !f.fgets(1024, &mut buf) {
        return false;
    }
    let line = c_str(&buf);
    if let Some(s) = find_sub(line, b"HPCP") {
        spec.interval_sec = 3600;
        spec.value_offset = s as i32;
        spec.fmt = NwsFormat::Online60;
    } else if let Some(q) = find_sub(line, b"QPCP") {
        spec.interval_sec = 900;
        spec.value_offset = q as i32;
        spec.fmt = NwsFormat::Online15;
    } else {
        return false;
    }
    spec.time_offset = spec.interval_sec;

    for _n in 1..=5 {
        if !f.fgets(1024, &mut buf) {
            return false;
        }
        let line = c_str(&buf);
        if find_sub(line, b"COOP:").is_none() {
            continue;
        }
        let last = match line.iter().rposition(|&c| c == b':') {
            Some(l) => l,
            None => return false,
        };
        spec.data_offset = last as i32 - 11;
        return true;
    }
    false
}

fn nws_detect_format(f: &mut CFile) -> NwsFileSpec {
    let mut spec = NwsFileSpec::default();
    let mut buf = Bytes::new();
    let mut elem: Bytes = Bytes::new();
    let mut rec: Bytes = Bytes::new();
    let mut coop: Bytes = Bytes::new();
    let mut sn2: i64 = 0;
    let mut div: i32 = 0;
    let mut year: i32 = 0;

    f.fseek(0);
    for _line_count in 1..=5 {
        if !f.fgets(1024, &mut buf) {
            return spec;
        }
        let line = c_str(&buf).to_vec();
        let len = line.len();

        let accept = |spec: &mut NwsFileSpec, elem: &[u8], fmt: NwsFormat, named: bool| -> bool {
            spec.interval_sec = nws_interval(elem);
            spec.time_offset = spec.interval_sec;
            spec.has_station_name = named;
            if spec.interval_sec > 0 {
                spec.fmt = fmt;
            }
            spec.interval_sec > 0
        };

        if sscanf(&line, b"%6ld %2d %4s", &mut [ScanArg::Long(&mut sn2), ScanArg::Int(&mut div), ScanArg::Str(&mut elem)]) == 3
            && accept(&mut spec, &elem, NwsFormat::SpaceDelimited, false)
        {
            return spec;
        }

        if len > 37
            && sscanf(
                &line[37..],
                b"%2d %4s %2s %4d",
                &mut [ScanArg::Int(&mut div), ScanArg::Str(&mut elem), ScanArg::Str(&mut rec), ScanArg::Int(&mut year)],
            ) == 4
            && accept(&mut spec, &elem, NwsFormat::SpaceDelimited, true)
        {
            return spec;
        }

        if sscanf(&line, b"%6ld,%2d,%4s", &mut [ScanArg::Long(&mut sn2), ScanArg::Int(&mut div), ScanArg::Str(&mut elem)]) == 3
            && accept(&mut spec, &elem, NwsFormat::CommaDelimited, false)
        {
            return spec;
        }

        if len > 37
            && sscanf(
                &line[37..],
                b"%2d,%4s,%2s,%4d",
                &mut [ScanArg::Int(&mut div), ScanArg::Str(&mut elem), ScanArg::Str(&mut rec), ScanArg::Int(&mut year)],
            ) == 4
            && accept(&mut spec, &elem, NwsFormat::CommaDelimited, true)
        {
            return spec;
        }

        if sscanf(
            &line,
            b"%3s%6ld%2d%4s",
            &mut [ScanArg::Str(&mut rec), ScanArg::Long(&mut sn2), ScanArg::Int(&mut div), ScanArg::Str(&mut elem)],
        ) == 4
            && accept(&mut spec, &elem, NwsFormat::Tape, false)
        {
            return spec;
        }

        if sscanf(&line, b"%5s%6ld", &mut [ScanArg::Str(&mut coop), ScanArg::Long(&mut sn2)]) == 2 && coop.as_slice() == b"COOP:" {
            if !nws_detect_online(f, &mut spec) {
                spec.fmt = NwsFormat::None;
            }
            return spec;
        }

        spec.hdr_lines += 1;
    }
    spec
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum NwsCondition {
    None,
    Accumulated,
    Deleted,
    Missing,
}

fn nws_condition(flag: u8) -> NwsCondition {
    match flag {
        b'a' | b'A' => NwsCondition::Accumulated,
        b'{' | b'}' => NwsCondition::Deleted,
        b'[' | b']' => NwsCondition::Missing,
        _ => NwsCondition::None,
    }
}

fn nws_read_online_value(s: &[u8], v: &mut i64, flag: &mut u8) -> i32 {
    if c_str(s).contains(&b'.') {
        let mut x: f32 = 99.99;
        let n = sscanf(s, b"%f %c", &mut [ScanArg::Float(&mut x), ScanArg::Char(flag)]);
        *v = cvt_i64_f32(100.0f32 * x + 0.5f32);
        return n;
    }
    sscanf(s, b"%ld %c", &mut [ScanArg::Long(v), ScanArg::Char(flag)])
}

/// The bytes `&line[k]` points at, as a C string (a k past the buffer reads as empty).
fn at(line: &[u8], k: i32) -> &[u8] {
    if k < 0 || k as usize >= line.len() {
        &[]
    } else {
        &line[k as usize..]
    }
}

fn nws_read_file(
    f: &mut CFile,
    spec: &NwsFileSpec,
    file_start: f64,
    series: &mut Table,
    first_date: &mut f64,
    last_date: &mut f64,
    periods_precip: &mut i64,
) {
    // char line[4096], kept between fgets calls (see the header).
    let mut line = vec![0u8; 4096];
    let mut chunk = Bytes::new();
    let mut fgets = |f: &mut CFile, line: &mut Vec<u8>| -> bool {
        if !f.fgets(4096, &mut chunk) {
            return false;
        }
        line[..chunk.len()].copy_from_slice(&chunk);
        line[chunk.len()] = 0;
        true
    };
    #[allow(unused_assignments)]
    let mut condition = NwsCondition::None;
    let mut accum_start = 0.0f64;

    let emit = |series: &mut Table, first_date: &mut f64, last_date: &mut f64, periods_precip: &mut i64, day: f64, hour: i32, minute: i32, x: f32| {
        let when = datetime::add_seconds(day, 3600.0 * hour as f64 + 60.0 * minute as f64 - spec.time_offset as f64);
        if *first_date == 0.0 || when < *first_date {
            *first_date = when;
        }
        if when > *last_date {
            *last_date = when;
        }
        if x > 0.0f32 {
            *periods_precip += 1;
        }
        series.x.push(when);
        series.y.push(x as f64);
    };

    f.fseek(0);
    for _i in 1..=spec.hdr_lines {
        if !fgets(f, &mut line) {
            return;
        }
    }

    while fgets(f, &mut line) {
        let line_length = c_str(&line).len() as i32 - 1;
        let (mut y, mut m, mut d, mut n) = (0i32, 0i32, 0i32, 0i32);
        let mut k: i32;

        match spec.fmt {
            NwsFormat::Tape => {
                if line_length <= 30 {
                    continue;
                }
                if sscanf(at(&line, 17), b"%4d%2d%4d%3d", &mut [ScanArg::Int(&mut y), ScanArg::Int(&mut m), ScanArg::Int(&mut d), ScanArg::Int(&mut n)]) < 4 {
                    continue;
                }
                k = 30;
            }
            NwsFormat::SpaceDelimited => {
                let name_length = if spec.has_station_name { 31 } else { 0 };
                if line_length <= 28 + name_length {
                    continue;
                }
                k = 18 + name_length;
                if sscanf(at(&line, k), b"%4d %2d %2d", &mut [ScanArg::Int(&mut y), ScanArg::Int(&mut m), ScanArg::Int(&mut d)]) < 3 {
                    continue;
                }
                k += 10;
            }
            NwsFormat::CommaDelimited => {
                if line_length <= 28 {
                    continue;
                }
                if sscanf(at(&line, 18), b"%4d,%2d,%2d", &mut [ScanArg::Int(&mut y), ScanArg::Int(&mut m), ScanArg::Int(&mut d)]) < 3 {
                    continue;
                }
                k = 28;
            }
            NwsFormat::Online60 | NwsFormat::Online15 => {
                if line_length <= spec.data_offset + 23 {
                    continue;
                }
                if sscanf(at(&line, spec.data_offset), b"%4d%2d%2d", &mut [ScanArg::Int(&mut y), ScanArg::Int(&mut m), ScanArg::Int(&mut d)]) < 3 {
                    continue;
                }
                k = spec.data_offset + 8;
            }
            NwsFormat::None => return,
        }

        let mut date1 = datetime::encode_date(y, m, d);
        if date1 < file_start {
            continue;
        }

        while k < line_length {
            let mut flag1: u8 = 0;
            let mut flag2: u8 = 0;
            let mut v: i64 = 99999;
            let mut hour: i32 = 25;
            let mut minute: i32 = 0;

            match spec.fmt {
                NwsFormat::Tape => {
                    n = sscanf(
                        at(&line, k),
                        b"%2d%2d%6ld%c%c",
                        &mut [ScanArg::Int(&mut hour), ScanArg::Int(&mut minute), ScanArg::Long(&mut v), ScanArg::Char(&mut flag1), ScanArg::Char(&mut flag2)],
                    );
                    k += 12;
                }
                NwsFormat::SpaceDelimited => {
                    n = sscanf(
                        at(&line, k),
                        b" %2d%2d %6ld %c %c",
                        &mut [ScanArg::Int(&mut hour), ScanArg::Int(&mut minute), ScanArg::Long(&mut v), ScanArg::Char(&mut flag1), ScanArg::Char(&mut flag2)],
                    );
                    k += 16;
                }
                NwsFormat::CommaDelimited => {
                    n = sscanf(
                        at(&line, k),
                        b",%2d%2d,%6ld,%c,%c",
                        &mut [ScanArg::Int(&mut hour), ScanArg::Int(&mut minute), ScanArg::Long(&mut v), ScanArg::Char(&mut flag1), ScanArg::Char(&mut flag2)],
                    );
                    k += 16;
                }
                _ => {
                    n = sscanf(at(&line, k), b" %2d:%2d", &mut [ScanArg::Int(&mut hour), ScanArg::Int(&mut minute)]);
                    n += nws_read_online_value(at(&line, spec.value_offset), &mut v, &mut flag1);
                    if hour == 0 {
                        hour = 24;
                        date1 -= 1.0;
                    }
                    k += line_length;
                }
            }

            if n < 3 || hour >= 25 {
                break;
            }

            condition = nws_condition(flag1);
            let is_missing = condition == NwsCondition::Deleted || condition == NwsCondition::Missing || flag1 == b'M' || v >= 9999;

            if flag1 == b'a' {
                accum_start = date1 + datetime::encode_time(hour, minute, 0);
            } else if flag1 == b'A' {
                if accum_start != 0.0 && v != 99999 && spec.interval_sec > 0 {
                    let date2 = date1 + datetime::encode_time(hour, minute, 0);
                    let np = (datetime::time_diff(date2, accum_start) / spec.interval_sec as i64) as i32 + 1;
                    let x = v as f32 / np as f32 / 100.0f32;
                    if x > 0.0f32 {
                        let mut when = datetime::add_seconds(accum_start, -spec.time_offset as f64);
                        for _j in 0..np {
                            if *first_date == 0.0 || when < *first_date {
                                *first_date = when;
                            }
                            if when > *last_date {
                                *last_date = when;
                            }
                            *periods_precip += 1;
                            series.x.push(when);
                            series.y.push(x as f64);
                            when = datetime::add_seconds(when, spec.interval_sec as f64);
                        }
                    }
                }
                accum_start = 0.0;
            } else if !is_missing {
                emit(series, first_date, last_date, periods_precip, date1, hour, minute, v as f32 / 100.0f32);
            }

            if flag1 == b'A' || flag1 == b'}' || flag1 == b']' {
                condition = NwsCondition::None;
            }
        }
    }
    let _ = condition;
}

// ---------------------------------------------------------------------------
// FILE rain gages
// ---------------------------------------------------------------------------

pub(super) fn load_external_rain_files_impl(ctx: &mut SimulationContext, file_cache: &mut MultiColumnFileCache) {
    let n_gages = ctx.gages.count();
    if n_gages == 0 {
        return;
    }

    for g in 0..n_gages {
        let ug = g as usize;
        if ctx.gages.source[ug] != RainSource::FileRain {
            continue;
        }

        let mut is_csv = ctx.gages.file_format[ug] == RainFileFormat::UserCsv;
        if !is_csv && ctx.gages.file_format[ug] != RainFileFormat::StanPrcp {
            continue;
        }

        let path: Bytes = if !ctx.gages.file_path[ug].absolute.is_empty() {
            ctx.gages.file_path[ug].absolute.clone()
        } else {
            ctx.gages.file_path[ug].original.clone()
        };

        if !is_csv && ctx.gages.station_id[ug].is_empty() && looks_like_multicolumn_series_file(&path) {
            ctx.gages.file_format[ug] = RainFileFormat::UserCsv;
            is_csv = true;
        }

        if is_csv {
            load_rain_file_user_csv(ctx, g, &path, file_cache);
            continue;
        }

        let mut fp = match fopen(&path) {
            Some(f) => f,
            None => match fopen(&ctx.gages.file_path[ug].original) {
                Some(f) => f,
                None => {
                    ctx.errors.push(format_error(ERR_RAIN_FILE_OPEN, &path));
                    continue;
                }
            },
        };

        let nws = nws_detect_format(&mut fp);
        if nws.fmt != NwsFormat::None {
            let mut nws_series = Table::new(ctx.gage_names.name_of(g), TableType::Timeseries);
            let mut nws_first = 0.0;
            let mut nws_last = 0.0;
            let mut nws_periods: i64 = 0;
            nws_read_file(&mut fp, &nws, ctx.gages.file_start_date[ug], &mut nws_series, &mut nws_first, &mut nws_last, &mut nws_periods);
            drop(fp);

            ctx.gages.interval_sec[ug] = nws.interval_sec;
            ctx.gages.file_first_date[ug] = nws_first;
            ctx.gages.file_last_date[ug] = nws_last;
            ctx.gages.file_periods_precip[ug] = nws_periods;
            ctx.gages.rain_series[ug] = nws_series;
            continue;
        }
        fp.fseek(0);

        let file_mm = ctx.gages.rain_units[ug] == 1;
        let units_factor_f: f32 = if file_mm { (1.0f64 / 25.40) as f32 } else { 1.0f32 };
        let rain_type_w = ctx.gages.rain_type[ug];
        let interval_f = ctx.gages.interval_sec[ug] as f32;
        let mut rain_accum_f: f32 = 0.0;

        let sta = ctx.gages.station_id[ug].clone();
        let file_start = ctx.gages.file_start_date[ug];

        let mut series = Table::new(ctx.gage_names.name_of(g), TableType::Timeseries);
        series.x.reserve(8192);
        series.y.reserve(8192);

        let mut first_date = 0.0f64;
        let mut last_date = 0.0f64;
        let mut periods_precip: i64 = 0;

        let mut line = Bytes::new();
        let mut tok = Bytes::new();
        let (mut yr, mut mo, mut dy, mut hr, mut mn) = (0i32, 0i32, 0i32, 0i32, 0i32);
        let mut val: f32 = 0.0;
        while fp.fgets(256, &mut line) {
            if line[0] == b';' || line[0] == b'\n' || line[0] == b'\r' {
                continue;
            }
            if sscanf(
                &line,
                b"%63s %d %d %d %d %d %f",
                &mut [
                    ScanArg::Str(&mut tok),
                    ScanArg::Int(&mut yr),
                    ScanArg::Int(&mut mo),
                    ScanArg::Int(&mut dy),
                    ScanArg::Int(&mut hr),
                    ScanArg::Int(&mut mn),
                    ScanArg::Float(&mut val),
                ],
            ) != 7
            {
                continue;
            }
            if !sta.is_empty() && sta != tok {
                continue;
            }

            let day = datetime::encode_date(yr, mo, dy);
            if day < file_start {
                continue;
            }
            let dt = day + datetime::encode_time(hr, mn, 0);

            if first_date == 0.0 || dt < first_date {
                first_date = dt;
            }
            if dt > last_date {
                last_date = dt;
            }
            if val > 0.0f32 {
                periods_precip += 1;
            }

            let mut x: f32 = val;
            if rain_type_w == 0 {
                x = x * interval_f / 3600.0f32;
            } else if rain_type_w == 2 {
                if x >= rain_accum_f {
                    x = x - rain_accum_f;
                    rain_accum_f += x;
                } else {
                    rain_accum_f = x;
                }
            }
            x *= units_factor_f;

            series.x.push(dt);
            series.y.push(x as f64);
        }
        drop(fp);

        // std::is_sorted, then std::sort of the row order by x.
        let sorted = !(1..series.x.len()).any(|i| series.x[i] < series.x[i - 1]);
        if !sorted {
            let mut order: Vec<usize> = (0..series.x.len()).collect();
            let xs = &series.x;
            std_sort(&mut order, &|a: &usize, b: &usize| xs[*a] < xs[*b]);
            let mut sx = vec![0.0f64; series.x.len()];
            let mut sy = vec![0.0f64; series.y.len()];
            for i in 0..order.len() {
                sx[i] = series.x[order[i]];
                sy[i] = series.y[order[i]];
            }
            std::mem::swap(&mut series.x, &mut sx);
            std::mem::swap(&mut series.y, &mut sy);
        }
        series.x.shrink_to_fit();
        series.y.shrink_to_fit();

        ctx.gages.file_first_date[ug] = first_date;
        ctx.gages.file_last_date[ug] = last_date;
        ctx.gages.file_periods_precip[ug] = periods_precip;
        ctx.gages.rain_series[ug] = series;
    }
}

/// (Re)load every FILE-source rain gage's data with a fresh multi-column cache.
pub fn load_external_rain_files(ctx: &mut SimulationContext) {
    let mut file_cache = MultiColumnFileCache::default();
    load_external_rain_files_impl(ctx, &mut file_cache);
}
