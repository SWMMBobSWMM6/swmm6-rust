// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/ClimateFile.hpp and ClimateFile.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The daily climate file reader (USER_PREPARED, GHCND, TD3200, DLY0204 formats).
// SWMMEngine::initHydrology opens it (open, setDayOffset); climate::initializeRadiation
// opens a private one and reads it day by day (getRecord), so the whole reader is here.
//
// Differences of form:
//   * The `std::FILE*` is a data::table::CFile (fopen "r", fgets, ftell with glibc's
//     behavior: a directory opens and then reads nothing). The destructor's fclose is
//     the drop of the CFile.
//   * Line buffers are C strings: `char line[MAX_LINE + 1]` filled by fgets is the Bytes
//     fgets returned, read through c_str (a NUL byte ends it, as strlen sees it).
//     saved_line_ holds strncpy's copy (the bytes before the first NUL; strncpy zero-pads
//     the rest of the array, so nothing else is held there).
//   * `strncpy(buf, &line[k], n); buf[n] = '\0'` is strncpy_at below.
//   * convertTemp, convertEvap and convertWind are declared in the header but defined
//     nowhere upstream; they are not here either.
//   * The members are public so the gate dumper can read them; upstream's are private.

use crate::core::cstdlib::{atof, atoi, c_str};
use crate::core::datetime;
use crate::core::options::Bytes;
use crate::data::table::{sscanf, CFile, ScanArg};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum ClimateFileFormat {
    Unknown = 0,
    /// StationID YYYY MM DD TMAX TMIN EVAP WIND
    UserPrepared = 1,
    /// NCDC Global Historical Climatology Network Daily
    Ghcnd = 2,
    /// NCDC TD3200 (NWS cooperative observer)
    Td3200 = 3,
    /// Canadian DLY02/DLY04
    Dly0204 = 4,
}

/// Climate variable indices (legacy ClimateVarType)
pub const TMIN: usize = 0;
pub const TMAX: usize = 1;
pub const EVAP: usize = 2;
pub const WIND: usize = 3;

/// Wind field type in GHCND
pub const WDMV: i32 = 0;
pub const AWND: i32 = 1;

/// Temperature unit encoding in GHCND
pub const DEG_C10: i32 = 0;
pub const DEG_C: i32 = 1;
pub const DEG_F: i32 = 2;

/// Daily climate record (output of a lookup). Members default to NAN.
#[derive(Debug, Clone, Copy)]
pub struct DailyClimateRecord {
    /// Minimum temperature (deg F)
    pub tmin: f64,
    /// Maximum temperature (deg F)
    pub tmax: f64,
    /// Pan evaporation (in/day US, mm/day SI; NAN if missing)
    pub evap: f64,
    /// Wind speed (mph; NAN if missing)
    pub wind: f64,
}

impl Default for DailyClimateRecord {
    fn default() -> Self {
        DailyClimateRecord { tmin: f64::NAN, tmax: f64::NAN, evap: f64::NAN, wind: f64::NAN }
    }
}

pub const MISSING: f64 = -1.0e10;
pub const MM_PER_INCH: f64 = 25.4;
pub const MAX_LINE: usize = 512;

#[derive(Debug)]
pub struct ClimateFileReader {
    // File state
    pub file_: Option<CFile>,
    pub format_: ClimateFileFormat,
    /// 0=US, 1=SI
    pub unit_system_: i32,

    /// Current month buffer: [variable][day], day 0 unused, 1-31 valid
    pub file_data_: [[f64; 32]; 4],
    pub buf_year_: i32,
    pub buf_month_: i32,
    pub day_offset_: f64,

    /// GHCND header-derived field positions (column start for TMIN,TMAX,EVAP,WIND)
    pub field_pos_: [i32; 4],
    pub date_field_pos_: i32,
    pub wind_type_: i32,
    pub temp_units_: i32,

    /// Saved line from month-boundary read-ahead
    pub saved_line_: Bytes,
    pub has_saved_line_: bool,
}

impl Default for ClimateFileReader {
    fn default() -> Self {
        ClimateFileReader {
            file_: None,
            format_: ClimateFileFormat::Unknown,
            unit_system_: 0,
            file_data_: [[0.0; 32]; 4],
            buf_year_: -1,
            buf_month_: -1,
            day_offset_: 0.0,
            field_pos_: [-1, -1, -1, -1],
            date_field_pos_: 0,
            wind_type_: WDMV,
            temp_units_: DEG_C10,
            saved_line_: Bytes::new(),
            has_saved_line_: false,
        }
    }
}

/// `char buf[n + 1] = {}; strncpy(buf, &line[k], n); buf[n] = '\0';` read as a C string.
/// `k` is at most strlen(line) wherever upstream calls it.
fn strncpy_at(line: &[u8], k: usize, n: usize) -> &[u8] {
    let l = c_str(line);
    let s = if k <= l.len() { &l[k..] } else { &l[l.len()..] };
    &s[..s.len().min(n)]
}

/// C `strstr(line, word)` offset in the C string `line`.
fn strstr(line: &[u8], word: &[u8]) -> Option<usize> {
    let l = c_str(line);
    if word.len() > l.len() {
        return None;
    }
    (0..=l.len() - word.len()).find(|&i| &l[i..i + word.len()] == word)
}

impl ClimateFileReader {
    pub fn close(&mut self) {
        self.file_ = None;
        self.format_ = ClimateFileFormat::Unknown;
        self.buf_year_ = -1;
        self.buf_month_ = -1;
        self.has_saved_line_ = false;
    }

    /// Open a climate file and detect its format. The start date is unused upstream.
    pub fn open(&mut self, path: &[u8], _start_oa_date: f64, unit_system: i32, temp_units: i32) -> bool {
        self.close();
        self.unit_system_ = unit_system;

        self.file_ = CFile::open(c_str(path));
        let f = match self.file_.as_mut() {
            None => return false,
            Some(f) => f,
        };

        // Read first line for format detection
        let mut line: Bytes = Bytes::new();
        if !f.fgets(MAX_LINE, &mut line) {
            self.close();
            return false;
        }

        self.format_ = self.detect_format(&line);
        if self.format_ == ClimateFileFormat::Unknown {
            self.close();
            return false;
        }

        // Config-specified units (legacy [TEMPERATURE] FILE keyword) win over the
        // per-format default/header auto-detection. -1 keeps the existing behavior.
        if temp_units >= 0 {
            self.temp_units_ = temp_units;
        }

        if self.format_ == ClimateFileFormat::Ghcnd {
            // Header already parsed by isGhcndFormat; do not save as data line
        } else {
            // For other formats, the first line IS data — save it
            self.saved_line_ = strncpy_at(&line, 0, MAX_LINE).to_vec();
            self.has_saved_line_ = true;
        }
        true
    }

    pub fn set_day_offset(&mut self, days: f64) {
        self.day_offset_ = days;
    }

    pub fn format(&self) -> ClimateFileFormat {
        self.format_
    }

    pub fn is_open(&self) -> bool {
        self.file_.is_some()
    }

    /// legacy getFileFormat
    fn detect_format(&mut self, line: &[u8]) -> ClimateFileFormat {
        let l = c_str(line);
        // 1. TD3200: starts with "DLY" and has "9999" at positions 23-26
        if l.len() > 26 && l[0] == b'D' && l[1] == b'L' && l[2] == b'Y' && l[23] == b'9' && l[24] == b'9' && l[25] == b'9' && l[26] == b'9' {
            return ClimateFileFormat::Td3200;
        }

        // 2. DLY0204: line >= 233 chars, positions 13-15 are param code 1, 2, or 151
        if l.len() >= 233 {
            let p = atoi(strncpy_at(l, 13, 3));
            if p == 1 || p == 2 || p == 151 {
                return ClimateFileFormat::Dly0204;
            }
        }

        // 3. USER_PREPARED: "StationID YYYY MM DD value..."
        {
            let (mut sta, mut s) = (Bytes::new(), Bytes::new());
            let (mut y, mut m, mut d) = (0i32, 0i32, 0i32);
            if sscanf(l, b"%79s %d %d %d %79s", &mut [ScanArg::Str(&mut sta), ScanArg::Int(&mut y), ScanArg::Int(&mut m), ScanArg::Int(&mut d), ScanArg::Str(&mut s)]) == 5 {
                return ClimateFileFormat::UserPrepared;
            }
        }

        // 4. GHCND: header contains "DATE" label
        if self.is_ghcnd_format(l) {
            return ClimateFileFormat::Ghcnd;
        }

        ClimateFileFormat::Unknown
    }

    fn is_ghcnd_format(&mut self, line: &[u8]) -> bool {
        let p = match strstr(line, b"DATE") {
            None => return false,
            Some(p) => p,
        };
        self.date_field_pos_ = p as i32;

        // Reset field positions
        for i in 0..4 {
            self.field_pos_[i] = -1;
        }

        if let Some(p) = strstr(line, b"TMIN") {
            self.field_pos_[TMIN] = p as i32;
        }
        if let Some(p) = strstr(line, b"TMAX") {
            self.field_pos_[TMAX] = p as i32;
        }
        if let Some(p) = strstr(line, b"EVAP") {
            self.field_pos_[EVAP] = p as i32;
        }

        // Wind: try WDMV first, then AWND
        if let Some(p) = strstr(line, b"WDMV") {
            self.field_pos_[WIND] = p as i32;
            self.wind_type_ = WDMV;
        } else if let Some(p) = strstr(line, b"AWND") {
            self.field_pos_[WIND] = p as i32;
            self.wind_type_ = AWND;
        }

        // Detect temperature units from column labels (e.g., "TMIN_C10"); the default
        // for GHCND is tenths of degrees Celsius either way.
        if strstr(line, b"C10").is_some() {
            self.temp_units_ = DEG_C10;
        } else {
            self.temp_units_ = DEG_C10;
        }
        true
    }

    fn oa_date_to_ymd(oa_date: f64) -> (i32, i32, i32) {
        datetime::decode_date(oa_date)
    }

    fn clear_buffer(&mut self) {
        for v in 0..4 {
            for d in 0..32 {
                self.file_data_[v][d] = MISSING;
            }
        }
    }

    fn process_line(&mut self, ln: &[u8]) {
        match self.format_ {
            ClimateFileFormat::UserPrepared => self.parse_user_line(ln),
            ClimateFileFormat::Ghcnd => self.parse_ghcnd_line(ln),
            ClimateFileFormat::Td3200 => self.parse_td3200_line(ln),
            ClimateFileFormat::Dly0204 => self.parse_dly0204_line(ln),
            _ => {}
        }
    }

    fn buffer_month(&mut self, year: i32, month: i32) -> bool {
        if self.file_.is_none() {
            return false;
        }
        self.clear_buffer();
        self.buf_year_ = year;
        self.buf_month_ = month;

        // Process saved line from previous read-ahead
        if self.has_saved_line_ {
            let saved = self.saved_line_.clone();
            self.process_line(&saved);
            self.has_saved_line_ = false;
        }

        // Read remaining lines for this month
        let mut raw: Bytes = Bytes::new();
        loop {
            if !self.file_.as_mut().unwrap().fgets(MAX_LINE, &mut raw) {
                break;
            }
            // The C buffer: the bytes fgets stored, then NUL.
            let mut line = raw.clone();
            // Remove trailing newline (len is strlen before either edit)
            let len = c_str(&line).len();
            if len > 0 && line[len - 1] == b'\n' {
                line[len - 1] = 0;
            }
            if len > 1 && line[len - 2] == b'\r' {
                line[len - 2] = 0;
            }

            // Peek at the date to check if we've crossed into the next month
            let (mut ly, mut lm, mut ld) = (0i32, 0i32, 0i32);
            let mut got_date = false;

            if self.format_ == ClimateFileFormat::UserPrepared {
                let mut sta = Bytes::new();
                if sscanf(&line, b"%79s %d %d %d", &mut [ScanArg::Str(&mut sta), ScanArg::Int(&mut ly), ScanArg::Int(&mut lm), ScanArg::Int(&mut ld)]) >= 4 {
                    got_date = true;
                }
            } else if self.format_ == ClimateFileFormat::Ghcnd {
                if c_str(&line).len() as i32 > self.date_field_pos_ + 8 {
                    let at = &line[self.date_field_pos_ as usize..];
                    if sscanf(at, b"%4d%2d%2d", &mut [ScanArg::Int(&mut ly), ScanArg::Int(&mut lm), ScanArg::Int(&mut ld)]) == 3 {
                        got_date = true;
                    }
                }
            } else if self.format_ == ClimateFileFormat::Td3200 {
                if c_str(&line).len() > 22 {
                    ly = atoi(strncpy_at(&line, 17, 4));
                    lm = atoi(strncpy_at(&line, 21, 2));
                    got_date = true;
                }
            } else if self.format_ == ClimateFileFormat::Dly0204 {
                if c_str(&line).len() >= 13 {
                    ly = atoi(strncpy_at(&line, 7, 4));
                    lm = atoi(strncpy_at(&line, 11, 2));
                    got_date = true;
                }
            }
            let _ = ld;

            // If this line belongs to a future month, save it and stop
            if got_date && (ly > year || (ly == year && lm > month)) {
                self.saved_line_ = strncpy_at(&line, 0, MAX_LINE).to_vec();
                self.has_saved_line_ = true;
                break;
            }

            self.process_line(&line);
        }
        true
    }

    /// USER_PREPARED: StationID YYYY MM DD TMAX TMIN EVAP WIND
    fn parse_user_line(&mut self, line: &[u8]) {
        let mut sta = Bytes::new();
        let (mut s0, mut s1, mut s2, mut s3) = (Bytes::new(), Bytes::new(), Bytes::new(), Bytes::new());
        let (mut y, mut m, mut d) = (0i32, 0i32, 0i32);

        if sscanf(line, b"%79s %d %d %d %79s %79s %79s %79s", &mut [ScanArg::Str(&mut sta), ScanArg::Int(&mut y), ScanArg::Int(&mut m), ScanArg::Int(&mut d), ScanArg::Str(&mut s0), ScanArg::Str(&mut s1), ScanArg::Str(&mut s2), ScanArg::Str(&mut s3)]) < 4 {
            return;
        }
        if d < 1 || d > 31 {
            return;
        }
        let d = d as usize;
        let first = |s: &Bytes| -> u8 { if s.is_empty() { 0 } else { s[0] } };

        // TMAX
        if first(&s0) != 0 && first(&s0) != b'*' {
            let mut v = atof(&s0);
            // User-prepared files: temps in user units, convert if SI
            if self.unit_system_ == 1 {
                v = v * 9.0 / 5.0 + 32.0; // degC -> degF
            }
            self.file_data_[TMAX][d] = v;
        }
        // TMIN
        if first(&s1) != 0 && first(&s1) != b'*' {
            let mut v = atof(&s1);
            if self.unit_system_ == 1 {
                v = v * 9.0 / 5.0 + 32.0;
            }
            self.file_data_[TMIN][d] = v;
        }
        // EVAP
        if first(&s2) != 0 && first(&s2) != b'*' {
            self.file_data_[EVAP][d] = atof(&s2);
        }
        // WIND
        if first(&s3) != 0 && first(&s3) != b'*' {
            self.file_data_[WIND][d] = atof(&s3);
        }
    }

    /// GHCND: date YYYYMMDD at date_field_pos_; data fields at field_pos_[i], 8 wide.
    fn parse_ghcnd_line(&mut self, line: &[u8]) {
        let len = c_str(line).len() as i32;
        let (mut y, mut m, mut d) = (0i32, 0i32, 0i32);

        if len <= self.date_field_pos_ + 8 {
            return;
        }
        if sscanf(&line[self.date_field_pos_ as usize..], b"%4d%2d%2d", &mut [ScanArg::Int(&mut y), ScanArg::Int(&mut m), ScanArg::Int(&mut d)]) != 3 {
            return;
        }
        if d < 1 || d > 31 {
            return;
        }

        for i in 0..4 {
            if self.field_pos_[i] < 0 || self.field_pos_[i] + 8 > len {
                continue;
            }
            let mut v = 0.0f64;
            if sscanf(&line[self.field_pos_[i] as usize..], b"%8lf", &mut [ScanArg::Double(&mut v)]) != 1 {
                continue;
            }

            // Convert based on variable type and units
            if i == TMIN || i == TMAX {
                if self.temp_units_ == DEG_C10 {
                    v = v / 10.0 * 9.0 / 5.0 + 32.0; // tenths degC -> degF
                } else if self.temp_units_ == DEG_C {
                    v = v * 9.0 / 5.0 + 32.0; // degC -> degF
                }
                // DEG_F: no conversion
            } else if i == EVAP {
                if self.temp_units_ == DEG_C10 {
                    v /= 10.0; // tenths mm -> mm
                }
                if self.unit_system_ == 0 {
                    v /= MM_PER_INCH; // mm -> inches for US
                }
            } else if i == WIND {
                if self.wind_type_ == WDMV {
                    v = v / 10.0; // tenths km -> km
                    v = v * 0.621371 / 24.0; // km/day -> mph
                } else {
                    v = v / 10.0 * 2.23694; // tenths m/s -> m/s -> mph
                }
            }

            self.file_data_[i][d as usize] = v;
        }
    }

    /// TD3200: 0-2="DLY", 11-14=param, 17-20=year, 21-22=month; data at 30+, 12 chars a day.
    fn parse_td3200_line(&mut self, line: &[u8]) {
        let len = c_str(line).len() as i32;
        if len < 30 {
            return;
        }

        // Identify parameter: char param[5], trailing spaces trimmed from index 3 down.
        let mut param = [0u8; 5];
        let src = strncpy_at(line, 11, 4);
        param[..src.len()].copy_from_slice(src);
        for i in (0..=3).rev() {
            if param[i] == b' ' {
                param[i] = 0;
            } else {
                break;
            }
        }
        let param = c_str(&param);

        let var: usize;
        if param == b"TMIN" {
            var = TMIN;
        } else if param == b"TMAX" {
            var = TMAX;
        } else if param == b"EVAP" {
            var = EVAP;
        } else if param == b"WDMV" {
            var = WIND;
            self.wind_type_ = WDMV;
        } else if param == b"AWND" {
            var = WIND;
            self.wind_type_ = AWND;
        } else {
            return;
        }

        // Number of values
        let nvals = atoi(strncpy_at(line, 27, 3));

        for j in 0..nvals {
            let k = 30 + j * 12;
            if k + 12 > len {
                break;
            }
            let k = k as usize;

            // Day of month (2 chars)
            let d = atoi(strncpy_at(line, k, 2));
            if d < 1 || d > 31 {
                continue;
            }

            // Sign + value (1 + 5 chars)
            let sign = line[k + 4];
            let vstr = strncpy_at(line, k + 5, 5);

            // Quality flag
            let flag = line[k + 11];
            if vstr == b"99999" {
                continue; // missing
            }
            if flag != b'0' && flag != b'1' && flag != b' ' {
                continue; // bad quality
            }

            let mut v = atof(vstr);
            if sign == b'-' {
                v = -v;
            }

            // Unit conversions
            if var == TMIN || var == TMAX {
                // Already in degF for TD3200
            } else if var == EVAP {
                v /= 100.0; // hundredths of inches -> inches
                if self.unit_system_ == 1 {
                    v *= MM_PER_INCH; // -> mm for SI
                }
            } else if var == WIND {
                v /= 24.0; // miles/day -> mph
            }

            self.file_data_[var][d as usize] = v;
        }
    }

    /// DLY0204: 7-10=year, 11-12=month, 13-15=param code; data at 16, 31 x 7-char blocks.
    fn parse_dly0204_line(&mut self, line: &[u8]) {
        let len = c_str(line).len() as i32;
        if len < 233 {
            return;
        }

        // Parameter code
        let pcode = atoi(strncpy_at(line, 13, 3));

        let var = if pcode == 1 {
            TMAX
        } else if pcode == 2 {
            TMIN
        } else if pcode == 151 {
            EVAP
        } else {
            return; // skip other parameters
        };

        let mut k = 16usize;
        for d in 1..=31usize {
            let sign = line[k];
            let vstr = strncpy_at(line, k + 1, 5);
            k += 7;

            // Skip missing values
            if vstr == b"99999" || vstr == b"     " {
                continue;
            }

            let mut v = atof(vstr);
            if sign == b'-' {
                v = -v;
            }

            // Unit conversions
            if var == TMIN || var == TMAX {
                // Tenths of degrees Celsius -> degF
                v = v / 10.0 * 9.0 / 5.0 + 32.0;
            } else if var == EVAP {
                v /= 10.0; // tenths mm -> mm
                if self.unit_system_ == 0 {
                    v /= MM_PER_INCH; // -> inches for US
                }
            }

            self.file_data_[var][d] = v;
        }
    }

    /// Daily record for `oa_date` (shifted by the day offset), buffering its month first.
    pub fn get_record(&mut self, oa_date: f64, rec: &mut DailyClimateRecord) -> bool {
        let (y, m, d) = Self::oa_date_to_ymd(oa_date + self.day_offset_);

        // Buffer the month if needed
        if y != self.buf_year_ || m != self.buf_month_ {
            if !self.buffer_month(y, m) {
                return false;
            }
        }

        if d < 1 || d > 31 {
            return false;
        }
        let d = d as usize;

        let pick = |x: f64| if x > MISSING + 1.0 { x } else { f64::NAN };
        rec.tmin = pick(self.file_data_[TMIN][d]);
        rec.tmax = pick(self.file_data_[TMAX][d]);
        rec.evap = pick(self.file_data_[EVAP][d]);
        rec.wind = pick(self.file_data_[WIND][d]);
        true
    }
}
