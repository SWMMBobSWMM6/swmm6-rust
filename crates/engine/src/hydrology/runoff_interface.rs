// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/RunoffInterface.hpp and
// RunoffInterface.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The runoff interface file ([FILES] SAVE/USE RUNOFF, legacy "Frunoff"): the part
// SWMMEngine::start reaches, openForWrite (header with a maxSteps placeholder), openForRead
// (header check), close (maxSteps patched in) and isOpen.
//
// NOT YET TRANSLATED (run-time, with stepRunoff): saveResults, readResults.
//
// Differences of form, never in result:
//   * The FILE* is a std::fs::File, read through data::table::CFile in read mode. Writes go
//     straight to the file rather than through a stdio buffer; nothing reads the file while
//     it is open for writing.
//   * Members keep upstream's names and are pub so the gate's state dumper can read them.

use std::io::{Seek, SeekFrom, Write};

use crate::core::unit_conversion as ucf;
use crate::data::table::CFile;

/// Fixed result fields per subcatchment
pub const N_FIXED_RESULTS: i32 = 8;
pub const FILE_STAMP: &[u8; 12] = b"SWMM5-RUNOFF";

#[derive(Debug)]
pub enum RunoffFile {
    Write(std::fs::File),
    Read(CFile),
}

#[derive(Debug, Default)]
pub struct RunoffInterfaceFile {
    pub fp_: Option<RunoffFile>,
    pub writing_: bool,
    pub n_subcatch_: i32,
    pub n_pollut_: i32,
    /// N_FIXED_RESULTS + n_pollut
    pub n_results_: i32,
    /// Project flow units (for UCF(FLOW) = Qcf)
    pub flow_units_: i32,
    /// 0 = US, 1 = SI (for Ucf lookups)
    pub unit_system_: i32,
    pub step_count_: i32,
    /// File position of maxSteps field
    pub max_steps_pos_: i64,
    /// Reusable buffer for read/write
    pub buf_: Vec<f32>,
}

/// fopen(path, "wb")
pub(crate) fn fopen_wb(path: &[u8]) -> Option<std::fs::File> {
    std::fs::OpenOptions::new().write(true).create(true).truncate(true).open(crate::core::file_io::path_from_bytes(path)).ok()
}

impl RunoffInterfaceFile {
    /// openForWrite: 0, or -1 when the file cannot be created.
    pub fn open_for_write(&mut self, path: &[u8], n_subcatch: i32, n_pollut: i32, flow_units: i32) -> i32 {
        let mut f = match fopen_wb(path) {
            Some(f) => f,
            None => return -1,
        };
        self.writing_ = true;
        self.n_subcatch_ = n_subcatch;
        self.n_pollut_ = n_pollut;
        self.n_results_ = N_FIXED_RESULTS + n_pollut;
        self.flow_units_ = flow_units;
        self.unit_system_ = ucf::get_unit_system(flow_units);
        self.step_count_ = 0;
        self.buf_.resize(self.n_results_ as usize, 0.0);

        // Write header (write errors go unchecked, as upstream's fwrite calls do)
        let _ = f.write_all(FILE_STAMP);
        let _ = f.write_all(&n_subcatch.to_le_bytes());
        let _ = f.write_all(&n_pollut.to_le_bytes());
        let _ = f.write_all(&flow_units.to_le_bytes());
        self.max_steps_pos_ = f.stream_position().map(|p| p as i64).unwrap_or(-1);
        let _ = f.write_all(&0i32.to_le_bytes());
        self.fp_ = Some(RunoffFile::Write(f));
        0
    }

    /// openForRead: 0; -1 cannot open, -2 not a runoff file, -3 incompatible with the model.
    pub fn open_for_read(&mut self, path: &[u8], n_subcatch: i32, n_pollut: i32, flow_units: i32) -> i32 {
        let mut f = match CFile::open(path) {
            Some(f) => f,
            None => return -1,
        };
        self.writing_ = false;
        let mut stamp = f.fread(12);
        stamp.resize(12, 0);
        // strncmp(stamp, FILE_STAMP, 12) over a char[12] without a terminator
        let n = stamp.iter().position(|&c| c == 0).map(|p| p + 1).unwrap_or(12).min(12);
        if stamp[..n] != FILE_STAMP[..n] {
            return -2;
        }
        let rd = |f: &mut CFile| -> i32 {
            let b = f.fread(4);
            if b.len() == 4 {
                i32::from_le_bytes([b[0], b[1], b[2], b[3]])
            } else {
                // an uninitialized int32 upstream; a short header fails the checks below
                // unless it happens to match (not reproducible): treat as a mismatch
                i32::MIN
            }
        };
        let ns = rd(&mut f);
        let np = rd(&mut f);
        let fu = rd(&mut f);
        let _ms = rd(&mut f);
        if ns != n_subcatch || np != n_pollut || fu != flow_units {
            return -3;
        }
        self.n_subcatch_ = n_subcatch;
        self.n_pollut_ = n_pollut;
        self.n_results_ = N_FIXED_RESULTS + n_pollut;
        self.flow_units_ = flow_units;
        self.unit_system_ = ucf::get_unit_system(flow_units);
        self.step_count_ = 0;
        self.buf_.resize(self.n_results_ as usize, 0.0);
        self.fp_ = Some(RunoffFile::Read(f));
        0
    }

    /// close: patch maxSteps into a written header, then close.
    pub fn close(&mut self) {
        let fp = match self.fp_.take() {
            Some(f) => f,
            None => return,
        };
        if let RunoffFile::Write(mut f) = fp {
            if self.writing_ && self.max_steps_pos_ > 0 {
                let _ = f.seek(SeekFrom::Start(self.max_steps_pos_ as u64));
                let _ = f.write_all(&self.step_count_.to_le_bytes());
            }
        }
    }

    pub fn is_open(&self) -> bool {
        self.fp_.is_some()
    }

    /// The open file's position (the gate dumps it): ftell, or -1 when closed.
    pub fn tell(&mut self) -> i64 {
        match self.fp_.as_mut() {
            Some(RunoffFile::Write(f)) => f.stream_position().map(|p| p as i64).unwrap_or(-1),
            Some(RunoffFile::Read(f)) => f.ftell(),
            None => -1,
        }
    }
}
