// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/InterfaceFile.hpp and InterfaceFile.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Routing interface files: the manager's state and InterfaceManager::init; and (5 November
// 2026, for SWMMEngine::start) InterfaceFileSoA::resize, openFiles, readFileHeader,
// readNextPeriod, setOldValues, writeFileHeader and isOutletNode; and (5 November 2026, for
// SWMMEngine::stepRouting) readInflows, getFlow, getQuality and getIfaceNode.
//
// NOT YET TRANSLATED (run-time, with the output): writeOutfallResults.
//
// Differences in form, never in result:
//   * The input `FILE*` is data::table::CFile (fgets as glibc does it), the output one a
//     std::fs::File written straight through (upstream's stdio buffer is flushed at close;
//     nothing reads the file while it is open). closeFiles drops them.
//   * readNextPeriod's strtok_r over " \t\n\r" is a split on those bytes that skips empty
//     tokens, as strtok does.
//   * Members keep upstream's names and are `pub` (upstream's private ones too) so the
//     gate's state dumper can read them.

use crate::core::context::SimulationContext;
use crate::core::options::Bytes;

#[derive(Debug, Clone, Default)]
pub struct InterfaceFileSoA {
    pub count: i32,
    pub node_idx: Vec<i32>,
    pub old_values: Vec<f64>,
    pub new_values: Vec<f64>,
}

#[derive(Debug)]
pub struct InterfaceManager {
    pub data_: InterfaceFileSoA,
    pub infile_: Option<crate::data::table::CFile>,
    pub outfile_: Option<std::fs::File>,
    pub old_time_: f64,
    pub new_time_: f64,
    pub n_iface_polluts_: i32,
    pub n_values_per_node_: i32,
    pub iface_flow_units_: i32,
    pub iface_step_: i32,
    pub iface_frac_: f64,
    pub pollut_map_: Vec<i32>,
}

impl Default for InterfaceManager {
    fn default() -> Self {
        InterfaceManager {
            data_: InterfaceFileSoA::default(),
            infile_: None,
            outfile_: None,
            old_time_: 0.0,
            new_time_: 0.0,
            n_iface_polluts_: 0,
            n_values_per_node_: 1,
            iface_flow_units_: 0,
            iface_step_: 0,
            iface_frac_: 0.0,
            pollut_map_: Vec::new(),
        }
    }
}

impl InterfaceManager {
    pub fn init(&mut self, _ctx: &mut SimulationContext) {
        // Re-init safety: release any files left open by a previous run.
        self.close_files();
        self.data_.count = 0;
        self.infile_ = None;
        self.outfile_ = None;
        self.old_time_ = 0.0;
        self.new_time_ = 0.0;
        self.n_iface_polluts_ = 0;
        self.n_values_per_node_ = 1;
        self.iface_frac_ = 0.0;
        self.pollut_map_.clear();
    }

    /// closeFiles: close and forget both handles.
    pub fn close_files(&mut self) {
        self.infile_ = None;
        self.outfile_ = None;
    }

    pub fn get_num_iface_nodes(&self) -> i32 {
        self.data_.count
    }

    pub fn get_frac(&self) -> f64 {
        self.iface_frac_
    }
}

const MAXLINE: usize = 1024;
const NO_DATE: f64 = -1.0e10;

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

impl InterfaceFileSoA {
    /// resize(n): count and node indices (-1); the value arrays are sized by the caller.
    pub fn resize(&mut self, n: i32) {
        self.count = n;
        self.node_idx.clear();
        self.node_idx.resize(n as usize, -1);
    }
}

impl InterfaceManager {
    /// openFiles(infile, outfile): 0; -1 same file, -2 cannot write, -3 cannot read.
    pub fn open_files(&mut self, infile_path: &[u8], outfile_path: &[u8]) -> i32 {
        if !infile_path.is_empty() && !outfile_path.is_empty() && infile_path == outfile_path {
            return -1;
        }
        if !outfile_path.is_empty() {
            self.outfile_ = crate::hydrology::runoff_interface::fopen_wb(outfile_path);
            if self.outfile_.is_none() {
                return -2;
            }
        }
        if !infile_path.is_empty() {
            self.infile_ = crate::data::table::CFile::open(infile_path);
            if self.infile_.is_none() {
                return -3;
            }
        }
        0
    }

    fn fgets(&mut self, line: &mut Bytes) -> bool {
        match self.infile_.as_mut() {
            Some(f) => f.fgets(MAXLINE, line),
            None => false,
        }
    }

    /// readFileHeader: parse an inflow interface file's header and read its first period.
    /// 0, or -1 on any format error.
    pub fn read_file_header(&mut self, ctx: &mut SimulationContext) -> i32 {
        use crate::data::table::{sscanf, ScanArg};
        if self.infile_.is_none() {
            return -1;
        }
        let mut line = Bytes::new();
        let mut s1 = Bytes::new();
        let mut s2 = Bytes::new();

        // Line 1: check "SWMM5 Interface File"
        if !self.fgets(&mut line) {
            return -1;
        }
        if sscanf(&line, b"%s", &mut [ScanArg::Str(&mut s1)]) != 1 {
            return -1;
        }
        if !strncmp_eq(&s1, b"SWMM5", 5) {
            return -1;
        }
        // Line 2: skip title
        if !self.fgets(&mut line) {
            return -1;
        }
        // Line 3: reporting time step
        if !self.fgets(&mut line) {
            return -1;
        }
        let mut st = self.iface_step_;
        let n = sscanf(&line, b"%d", &mut [ScanArg::Int(&mut st)]);
        self.iface_step_ = st;
        if n != 1 || self.iface_step_ <= 0 {
            return -1;
        }
        // Line 4: number of constituents (including FLOW)
        let mut n_constit = 0i32;
        if !self.fgets(&mut line) {
            return -1;
        }
        if sscanf(&line, b"%d", &mut [ScanArg::Int(&mut n_constit)]) != 1 {
            return -1;
        }
        self.n_iface_polluts_ = n_constit.wrapping_sub(1);
        if self.n_iface_polluts_ < 0 {
            return -1;
        }
        self.n_values_per_node_ = 1 + self.n_iface_polluts_;
        // Line 5: FLOW units
        if !self.fgets(&mut line) {
            return -1;
        }
        if sscanf(&line, b"%s %s", &mut [ScanArg::Str(&mut s1), ScanArg::Str(&mut s2)]) < 2 {
            return -1;
        }
        if !strncmp_eq(&s1, b"FLOW", 4) {
            return -1;
        }
        const FLOW_WORDS: [&[u8]; 6] = [b"CFS", b"GPM", b"MGD", b"CMS", b"LPS", b"MLD"];
        self.iface_flow_units_ = 0;
        for (fi, w) in FLOW_WORDS.iter().enumerate() {
            if strncmp_eq(&s2, w, 3) {
                self.iface_flow_units_ = fi as i32;
                break;
            }
        }
        // Pollutant mapping: project pollutant -> interface file column
        let n_polluts = ctx.n_pollutants();
        if n_polluts > 0 {
            self.pollut_map_.clear();
            self.pollut_map_.resize(n_polluts as usize, -1);
        }
        for i in 0..self.n_iface_polluts_ {
            if !self.fgets(&mut line) {
                return -1;
            }
            if sscanf(&line, b"%s %s", &mut [ScanArg::Str(&mut s1), ScanArg::Str(&mut s2)]) < 2 {
                return -1;
            }
            let j = ctx.pollutant_names.find(&s1);
            if j >= 0 && j < n_polluts {
                self.pollut_map_[j as usize] = i;
            }
        }
        // Number of interface nodes
        let mut n_nodes = 0i32;
        if !self.fgets(&mut line) {
            return -1;
        }
        if sscanf(&line, b"%d", &mut [ScanArg::Int(&mut n_nodes)]) != 1 || n_nodes <= 0 {
            return -1;
        }
        self.data_.resize(n_nodes);
        let total_values = n_nodes as usize * self.n_values_per_node_ as usize;
        self.data_.old_values.clear();
        self.data_.old_values.resize(total_values, 0.0);
        self.data_.new_values.clear();
        self.data_.new_values.resize(total_values, 0.0);
        for i in 0..n_nodes as usize {
            if !self.fgets(&mut line) {
                return -1;
            }
            let mut name = Bytes::new();
            if sscanf(&line, b"%s", &mut [ScanArg::Str(&mut name)]) != 1 {
                return -1;
            }
            self.data_.node_idx[i] = ctx.node_names.find(&name);
        }
        // Skip column headings line
        if !self.fgets(&mut line) {
            return -1;
        }
        // Read first set of interface values
        self.read_next_period();
        self.old_time_ = self.new_time_;
        0
    }

    /// readNextPeriod: read the next period's values; false at end of file or bad data.
    pub fn read_next_period(&mut self) -> bool {
        use crate::core::cstdlib::{atof, atoi};
        if self.infile_.is_none() {
            return false;
        }
        self.new_time_ = NO_DATE;
        for i in 0..self.data_.count {
            let ui = i as usize;
            // feof: fgets below returns false at end of file in the same cases
            let mut line = Bytes::new();
            if !self.fgets(&mut line) {
                return false;
            }
            let line = crate::core::cstdlib::c_str(&line).to_vec();
            let mut toks = line.split(|c| matches!(*c, b' ' | b'\t' | b'\n' | b'\r')).filter(|t| !t.is_empty());
            if toks.next().is_none() {
                return false;
            }
            let mut f = [0i32; 6];
            for v in f.iter_mut() {
                match toks.next() {
                    Some(t) => *v = atoi(t),
                    None => return false,
                }
            }
            let (yr, mon, day, hr, mn, sec) = (f[0], f[1], f[2], f[3], f[4], f[5]);
            let tok = match toks.next() {
                Some(t) => t,
                None => return false,
            };
            let base = ui * self.n_values_per_node_ as usize;
            self.data_.new_values[base] = atof(tok);
            if self.iface_flow_units_ >= 0 && self.iface_flow_units_ <= 5 {
                self.data_.new_values[base] /= crate::core::unit_conversion::QCF[self.iface_flow_units_ as usize];
            }
            for j in 1..=self.n_iface_polluts_ {
                let tok = match toks.next() {
                    Some(t) => t,
                    None => return false,
                };
                self.data_.new_values[base + j as usize] = atof(tok);
            }
            self.new_time_ = crate::core::datetime::encode_date(yr, mon, day) + crate::core::datetime::encode_time(hr, mn, sec);
        }
        true
    }

    /// setOldValues: shift new values to old for the next interpolation.
    pub fn set_old_values(&mut self) {
        self.old_time_ = self.new_time_;
        let total = self.data_.count as usize * self.n_values_per_node_ as usize;
        for i in 0..total {
            self.data_.old_values[i] = self.data_.new_values[i];
        }
    }

    /// writeFileHeader: the outflow interface file's header.
    pub fn write_file_header(&mut self, ctx: &SimulationContext) {
        use std::io::Write;
        let f = match self.outfile_.as_mut() {
            Some(f) => f,
            None => return,
        };
        let mut o: Vec<u8> = Vec::new();
        o.extend_from_slice(b"SWMM5 Interface File");
        o.extend_from_slice(b"\nOpenSWMM Engine Output");
        let report_step = crate::core::cmath::cvt_i32(ctx.options.report_step);
        o.extend_from_slice(format!("\n{:<4} - reporting time step in sec", report_step).as_bytes());
        let n_polluts = ctx.n_pollutants();
        o.extend_from_slice(format!("\n{:<4} - number of constituents as listed below:", n_polluts + 1).as_bytes());
        const FLOW_WORDS: [&[u8]; 6] = [b"CFS", b"GPM", b"MGD", b"CMS", b"LPS", b"MLD"];
        let mut fu = ctx.options.flow_units as i32;
        if !(0..=5).contains(&fu) {
            fu = 0;
        }
        o.extend_from_slice(b"\nFLOW ");
        o.extend_from_slice(FLOW_WORDS[fu as usize]);
        for i in 0..n_polluts {
            o.extend_from_slice(b"\n");
            o.extend_from_slice(crate::core::cstdlib::c_str(&ctx.pollutant_names.names()[i as usize]));
            o.extend_from_slice(b" MG/L");
        }
        let mut n_outlets = 0;
        for i in 0..ctx.n_nodes() {
            if Self::is_outlet_node(ctx, i) {
                n_outlets += 1;
            }
        }
        o.extend_from_slice(format!("\n{:<4} - number of nodes as listed below:", n_outlets).as_bytes());
        for i in 0..ctx.n_nodes() {
            if !Self::is_outlet_node(ctx, i) {
                continue;
            }
            o.extend_from_slice(b"\n");
            o.extend_from_slice(crate::core::cstdlib::c_str(&ctx.node_names.names()[i as usize]));
        }
        o.extend_from_slice(b"\nNode             Year Mon Day Hr  Min Sec FLOW      ");
        for i in 0..n_polluts {
            let pname = crate::core::cstdlib::c_str(&ctx.pollutant_names.names()[i as usize]);
            o.push(b' ');
            o.extend_from_slice(pname);
            for _ in pname.len()..10 {
                o.push(b' ');
            }
        }
        let _ = f.write_all(&o);
    }

    /// isOutletNode: an OUTFALL under DYNWAVE/FV; a node with no outflow links otherwise.
    pub fn is_outlet_node(ctx: &SimulationContext, node_idx: i32) -> bool {
        use crate::core::options::RoutingModel;
        use crate::data::node::NodeType;
        let ui = node_idx as usize;
        if ctx.options.routing_model == RoutingModel::Dynwave || ctx.options.routing_model == RoutingModel::Fv {
            return ctx.nodes.type_[ui] == NodeType::Outfall;
        }
        ctx.nodes.degree[ui] == 0
    }
}

// ============================================================================
// Run time (5 November 2026, for SWMMEngine::stepRouting): readInflows, getFlow, getQuality,
// getIfaceNode
// ============================================================================

impl InterfaceManager {
    /// readInflows: add the interface file's interpolated flows (and positive-flow pollutant
    /// loads) at current_time to the nodes' interface inflow.
    pub fn read_inflows(&mut self, ctx: &mut SimulationContext, current_time: f64) {
        if self.infile_.is_none() {
            return;
        }
        // Fallback: read the header if start() has not done so eagerly.
        if self.data_.old_values.is_empty() && self.read_file_header(ctx) != 0 {
            self.close_files();
            return;
        }
        if self.data_.count == 0 {
            return;
        }
        if self.old_time_ > current_time {
            return;
        }
        while self.new_time_ < current_time && self.new_time_ != NO_DATE {
            self.set_old_values();
            self.read_next_period();
        }
        if self.new_time_ == NO_DATE {
            return;
        }
        let dt_iface = self.new_time_ - self.old_time_;
        if dt_iface > 0.0 {
            self.iface_frac_ = (current_time - self.old_time_) / dt_iface;
            self.iface_frac_ = crate::core::cmath::cmax(0.0, crate::core::cmath::cmin(1.0, self.iface_frac_));
        } else {
            self.iface_frac_ = 0.0;
        }
        for i in 0..self.data_.count {
            let node = self.data_.node_idx[i as usize];
            if node < 0 || node >= ctx.n_nodes() {
                continue;
            }
            let flow = self.get_flow(i, self.iface_frac_);
            // legacy addIfaceInflows: the flow and its loads are skipped below FLOW_TOL
            if flow.abs() < crate::core::constants::FLOW_TOL {
                continue;
            }
            ctx.nodes.iface_inflow[node as usize] += flow;
            let n_polluts = ctx.n_pollutants();
            if flow > 0.0 && n_polluts > 0 && !self.pollut_map_.is_empty() && !ctx.nodes.iface_qual_mass.is_empty() {
                let base_q = node as usize * n_polluts as usize;
                for p in 0..n_polluts {
                    let qual = self.get_quality(i, p, self.iface_frac_);
                    if qual > 0.0 {
                        ctx.nodes.iface_qual_mass[base_q + p as usize] += flow * qual;
                    }
                }
            }
        }
    }

    /// getFlow: the interpolated flow of interface entry `index`.
    pub fn get_flow(&self, index: i32, frac: f64) -> f64 {
        if index < 0 || index >= self.data_.count {
            return 0.0;
        }
        let base = index as usize * self.n_values_per_node_ as usize;
        let q1 = self.data_.old_values[base];
        let q2 = self.data_.new_values[base];
        (1.0 - frac) * q1 + frac * q2
    }

    /// getQuality: the interpolated concentration of pollutant `pollut_idx` at entry `index`.
    pub fn get_quality(&self, index: i32, pollut_idx: i32, frac: f64) -> f64 {
        if index < 0 || index >= self.data_.count {
            return 0.0;
        }
        if pollut_idx < 0 || pollut_idx >= self.pollut_map_.len() as i32 {
            return 0.0;
        }
        let col = self.pollut_map_[pollut_idx as usize];
        if col < 0 {
            return 0.0;
        }
        let base = index as usize * self.n_values_per_node_ as usize;
        let offset = base + (col + 1) as usize;
        let c1 = self.data_.old_values[offset];
        let c2 = self.data_.new_values[offset];
        (1.0 - frac) * c1 + frac * c2
    }

    /// getIfaceNode: the project node of interface entry `index` (-1 outside the table).
    pub fn get_iface_node(&self, index: i32) -> i32 {
        if index >= 0 && index < self.data_.count {
            self.data_.node_idx[index as usize]
        } else {
            -1
        }
    }
}
