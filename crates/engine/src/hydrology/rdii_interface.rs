// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/RdiiInterface.hpp and RdiiInterface.cpp,
// and RDIISolver::rdiiNodeList in RDII.cpp (HydroCouple/openswmm.engine, branch swmm6_rel,
// commit 88329d6), Copyright 2026 Caleb Buahin, licensed under the Apache License,
// Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The RDII interface file ([FILES] USE/SAVE RDII, legacy "Frdii"): openForRead (binary or
// text header, then the first record primed), openForWrite (binary header), applyFlows,
// saveFlows, close.
//
// Differences of form, never in result:
//   * The FILE* is data::table::CFile when reading (fread/fgets as glibc does them) and a
//     std::fs::File when writing (written straight through, not buffered; nothing reads it
//     while it is open). Upstream reopens the file in text mode ("rt"), which is the same
//     byte stream on the reference platform.
//   * Members keep upstream's names and are pub so the gate's state dumper can read them.

use std::io::Write;

use crate::core::constants;
use crate::core::context::SimulationContext;
use crate::core::datetime;
use crate::core::options::Bytes;
use crate::core::unit_conversion as ucf;
use crate::data::table::{sscanf, CFile, ScanArg};
use crate::hydrology::rdii::RDIISolver;
use crate::hydrology::runoff_interface::fopen_wb;

pub const RDII_FILE_STAMP: &[u8] = b"SWMM5-RDII";
const MAXLINE: usize = 1024;
const NO_DATE: f64 = -1.0e10;

#[derive(Debug)]
pub enum RdiiFile {
    Write(std::fs::File),
    Read(CFile),
}

#[derive(Debug)]
pub struct RdiiInterfaceFile {
    pub fp_: Option<RdiiFile>,
    pub writing_: bool,
    pub binary_: bool,
    /// Record validity interval (sec)
    pub rdii_step_: i32,
    /// Text format: declared flow units
    pub flow_units_: i32,
    /// Interface entry -> project node
    pub node_idx_: Vec<i32>,
    /// Current record flows (cfs)
    pub flows_: Vec<f32>,
    /// Current record start (OADate)
    pub start_date_: f64,
    /// Current record end (OADate)
    pub end_date_: f64,
    pub have_record_: bool,
}

impl Default for RdiiInterfaceFile {
    fn default() -> Self {
        RdiiInterfaceFile {
            fp_: None,
            writing_: false,
            binary_: true,
            rdii_step_: 0,
            flow_units_: 0,
            node_idx_: Vec::new(),
            flows_: Vec::new(),
            start_date_: 0.0,
            end_date_: 0.0,
            have_record_: false,
        }
    }
}

/// strncmp(a, b, n) == 0 over C strings
fn strncmp_eq(a: &[u8], b: &[u8], n: usize) -> bool {
    for i in 0..n {
        let x = a.get(i).copied().unwrap_or(0);
        let y = b.get(i).copied().unwrap_or(0);
        if x != y {
            return false;
        }
        if x == 0 {
            return true;
        }
    }
    true
}

impl RdiiInterfaceFile {
    fn reader(&mut self) -> Option<&mut CFile> {
        match self.fp_.as_mut() {
            Some(RdiiFile::Read(f)) => Some(f),
            _ => None,
        }
    }

    fn fgets(&mut self, line: &mut Bytes) -> bool {
        match self.reader() {
            Some(f) => f.fgets(MAXLINE, line),
            None => false,
        }
    }

    /// openForRead: 0; -1 cannot open, -2 bad format (the file is closed again).
    pub fn open_for_read(&mut self, ctx: &mut SimulationContext, path: &[u8]) -> i32 {
        self.close();
        self.writing_ = false;
        self.have_record_ = false;
        self.start_date_ = NO_DATE;
        self.end_date_ = NO_DATE;

        let mut f = match CFile::open(path) {
            Some(f) => f,
            None => return -1,
        };
        let stamp_len = RDII_FILE_STAMP.len();
        let stamp = f.fread(stamp_len);
        if stamp.len() == stamp_len && strncmp_eq(&stamp, RDII_FILE_STAMP, stamp_len) {
            self.binary_ = true;
            self.fp_ = Some(RdiiFile::Read(f));
            let rc = self.read_binary_header(ctx);
            if rc != 0 {
                self.close();
                return rc;
            }
        } else {
            // Not binary: reopen in text mode (legacy openRdiiTextFile).
            drop(f);
            let f = match CFile::open(path) {
                Some(f) => f,
                None => return -1,
            };
            self.fp_ = Some(RdiiFile::Read(f));
            self.binary_ = false;
            let rc = self.read_text_header(ctx);
            if rc != 0 {
                self.close();
                return rc;
            }
        }
        // Prime the first record (legacy rdii_openRdii -> readRdiiFlows()).
        self.read_next_record();
        0
    }

    fn read_binary_header(&mut self, ctx: &SimulationContext) -> i32 {
        let f = self.reader().unwrap();
        let b = f.fread(4);
        if b.len() != 4 {
            return -2;
        }
        let step = i32::from_le_bytes([b[0], b[1], b[2], b[3]]);
        if step <= 0 {
            return -2;
        }
        let b = f.fread(4);
        if b.len() != 4 {
            return -2;
        }
        let n_nodes = i32::from_le_bytes([b[0], b[1], b[2], b[3]]);
        if n_nodes <= 0 {
            return -2;
        }
        let un = n_nodes as usize;
        let raw = f.fread(un * 4);
        self.rdii_step_ = step;
        self.node_idx_.clear();
        self.node_idx_.resize(un, -1);
        self.flows_.clear();
        self.flows_.resize(un, 0.0);
        if raw.len() / 4 != un {
            return -2;
        }
        for i in 0..un {
            let j = i32::from_le_bytes([raw[4 * i], raw[4 * i + 1], raw[4 * i + 2], raw[4 * i + 3]]);
            if j < 0 || j >= ctx.n_nodes() {
                return -2; // legacy: node must exist
            }
            self.node_idx_[i] = j;
        }
        0
    }

    fn read_text_header(&mut self, ctx: &SimulationContext) -> i32 {
        let mut line = Bytes::new();
        let mut s1 = Bytes::new();
        let mut s2 = Bytes::new();

        // Line 1: "SWMM5 ..." magic
        if !self.fgets(&mut line) {
            return -2;
        }
        if sscanf(&line, b"%s", &mut [ScanArg::Str(&mut s1)]) != 1 || !strncmp_eq(&s1, b"SWMM5", 5) {
            return -2;
        }
        // Line 2: title (skipped)
        if !self.fgets(&mut line) {
            return -2;
        }
        // Line 3: RDII time step (sec)
        if !self.fgets(&mut line) {
            return -2;
        }
        let mut step = self.rdii_step_;
        let n = sscanf(&line, b"%d", &mut [ScanArg::Int(&mut step)]);
        self.rdii_step_ = step;
        if n != 1 || self.rdii_step_ <= 0 {
            return -2;
        }
        // Line 4: number of constituents (skipped)
        if !self.fgets(&mut line) {
            return -2;
        }
        // Line 5: "FLOW <units>"
        if !self.fgets(&mut line) {
            return -2;
        }
        if sscanf(&line, b"%s %s", &mut [ScanArg::Str(&mut s1), ScanArg::Str(&mut s2)]) < 2 {
            return -2;
        }
        const FLOW_WORDS: [&[u8]; 6] = [b"CFS", b"GPM", b"MGD", b"CMS", b"LPS", b"MLD"];
        self.flow_units_ = -1;
        for (fi, w) in FLOW_WORDS.iter().enumerate() {
            if strncmp_eq(&s2, w, 3) {
                self.flow_units_ = fi as i32;
                break;
            }
        }
        if self.flow_units_ < 0 {
            return -2;
        }
        // Line 6: number of RDII nodes
        let mut n_nodes = 0i32;
        if !self.fgets(&mut line) {
            return -2;
        }
        if sscanf(&line, b"%d", &mut [ScanArg::Int(&mut n_nodes)]) != 1 || n_nodes <= 0 {
            return -2;
        }
        self.node_idx_.clear();
        self.node_idx_.resize(n_nodes as usize, -1);
        self.flows_.clear();
        self.flows_.resize(n_nodes as usize, 0.0);
        for i in 0..n_nodes as usize {
            if !self.fgets(&mut line) {
                return -2;
            }
            if sscanf(&line, b"%s", &mut [ScanArg::Str(&mut s1)]) != 1 {
                return -2;
            }
            let j = ctx.node_names.find(&s1);
            if j < 0 {
                return -2;
            }
            self.node_idx_[i] = j;
        }
        // Column headings line (skipped)
        if !self.fgets(&mut line) {
            return -2;
        }
        0
    }

    /// readNextRecord: advance to the next record (either format).
    fn read_next_record(&mut self) -> bool {
        self.have_record_ = false;
        self.start_date_ = NO_DATE;
        self.end_date_ = NO_DATE;
        if self.fp_.is_none() || self.writing_ {
            return false;
        }
        if self.binary_ {
            let nf = self.flows_.len();
            let f = self.reader().unwrap();
            let b = f.fread(8);
            if b.len() != 8 {
                return false;
            }
            let date = f64::from_le_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]]);
            let raw = f.fread(nf * 4);
            // fread fills what it read before reporting the short count
            for (k, c) in raw.chunks(4).enumerate() {
                if c.len() == 4 {
                    self.flows_[k] = f32::from_le_bytes([c[0], c[1], c[2], c[3]]);
                } else {
                    let mut cur = self.flows_[k].to_le_bytes();
                    cur[..c.len()].copy_from_slice(c);
                    self.flows_[k] = f32::from_le_bytes(cur);
                }
            }
            if raw.len() / 4 != nf {
                return false;
            }
            self.start_date_ = date;
        } else {
            let mut line = Bytes::new();
            let mut s = Bytes::new();
            let (mut yr, mut mon, mut day, mut hr, mut mn, mut sec) = (0i32, 0i32, 0i32, 0i32, 0i32, 0i32);
            let mut x = 0.0f64;
            for i in 0..self.flows_.len() {
                if !self.fgets(&mut line) {
                    return false;
                }
                if sscanf(
                    &line,
                    b"%s %d %d %d %d %d %d %lf",
                    &mut [
                        ScanArg::Str(&mut s),
                        ScanArg::Int(&mut yr),
                        ScanArg::Int(&mut mon),
                        ScanArg::Int(&mut day),
                        ScanArg::Int(&mut hr),
                        ScanArg::Int(&mut mn),
                        ScanArg::Int(&mut sec),
                        ScanArg::Double(&mut x),
                    ],
                ) < 8
                {
                    return false;
                }
                // Convert declared units -> cfs (legacy: x / Qcf[RdiiFlowUnits])
                self.flows_[i] = (x / ucf::QCF[self.flow_units_ as usize]) as f32;
            }
            self.start_date_ = datetime::encode_date(yr, mon, day) + datetime::encode_time(hr, mn, sec);
        }
        self.end_date_ = datetime::add_seconds(self.start_date_, self.rdii_step_ as f64);
        self.have_record_ = true;
        true
    }

    /// applyFlows: add the record bracketing current_date to the nodes' RDII inflow.
    pub fn apply_flows(&mut self, ctx: &mut SimulationContext, current_date: f64) {
        if self.fp_.is_none() || self.writing_ {
            return;
        }
        while self.have_record_ {
            if current_date < self.start_date_ {
                return;
            }
            if current_date < self.end_date_ {
                break;
            }
            self.read_next_record();
        }
        if !self.have_record_ {
            return;
        }
        for i in 0..self.node_idx_.len() {
            let j = self.node_idx_[i];
            if j < 0 || j >= ctx.n_nodes() {
                continue;
            }
            let q = self.flows_[i] as f64;
            if q.abs() < constants::FLOW_TOL {
                continue;
            }
            ctx.nodes.rdii_inflow[j as usize] += q;
        }
    }

    /// openForWrite: 0; -2 no RDII nodes (nothing to save), -1 cannot create the file.
    pub fn open_for_write(&mut self, path: &[u8], rdii_step: i32, node_idx: &[i32]) -> i32 {
        self.close();
        if node_idx.is_empty() {
            return -2;
        }
        let mut f = match fopen_wb(path) {
            Some(f) => f,
            None => return -1,
        };
        self.writing_ = true;
        self.binary_ = true;
        self.rdii_step_ = if rdii_step > 0 { rdii_step } else { 300 };
        self.node_idx_ = node_idx.to_vec();
        self.flows_.clear();
        self.flows_.resize(self.node_idx_.len(), 0.0);
        let _ = f.write_all(RDII_FILE_STAMP);
        let _ = f.write_all(&self.rdii_step_.to_le_bytes());
        let _ = f.write_all(&(self.node_idx_.len() as i32).to_le_bytes());
        for &j in &self.node_idx_ {
            let _ = f.write_all(&j.to_le_bytes());
        }
        self.fp_ = Some(RdiiFile::Write(f));
        0
    }

    /// saveFlows: one record (date, per-entry flows) from the per-node flow vector.
    pub fn save_flows(&mut self, current_date: f64, node_flows: &[f64]) {
        if !self.writing_ {
            return;
        }
        for i in 0..self.node_idx_.len() {
            let j = self.node_idx_[i] as usize;
            self.flows_[i] = if j < node_flows.len() { node_flows[j] as f32 } else { 0.0 };
        }
        if let Some(RdiiFile::Write(f)) = self.fp_.as_mut() {
            let _ = f.write_all(&current_date.to_le_bytes());
            for q in &self.flows_ {
                let _ = f.write_all(&q.to_le_bytes());
            }
        }
    }

    pub fn close(&mut self) {
        self.fp_ = None;
        self.writing_ = false;
        self.have_record_ = false;
    }

    pub fn is_open(&self) -> bool {
        self.fp_.is_some()
    }
    pub fn is_writing(&self) -> bool {
        self.writing_
    }
    pub fn step(&self) -> i32 {
        self.rdii_step_
    }

    /// The open file's position (the gate dumps it): ftell, or -1 when closed.
    pub fn tell(&mut self) -> i64 {
        match self.fp_.as_mut() {
            Some(RdiiFile::Write(f)) => std::io::Seek::stream_position(f).map(|p| p as i64).unwrap_or(-1),
            Some(RdiiFile::Read(f)) => f.ftell(),
            None => -1,
        }
    }
}

/// RDIISolver::rdiiNodeList: sorted unique node indices carrying an RDII UH inflow.
pub fn rdii_node_list(rdii: &RDIISolver) -> Vec<i32> {
    let mut nodes: Vec<i32> = rdii.groups_.node_idx.clone();
    nodes.sort();
    nodes.dedup();
    nodes
}
