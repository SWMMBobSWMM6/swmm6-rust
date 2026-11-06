// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/core/HotStartManager.hpp and HotStartManager.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The part of the hot-start manager SWMMEngine::initialize reaches with [FILES] USE HOTSTART:
// last_io_error, is_native_file, legacy_runoff_size, apply_legacy_routing (a legacy EPA SWMM5
// "SWMM5-HOTSTART<n>" file: the runoff section is handed back for applyLegacyRunoffSection,
// the node and link routing state is applied, and the optional OPENSWMM-HSMSX1 reaction
// species trailer is read) and apply_species_surface.
//
// NOT YET TRANSLATED: the native OPENSWMM_HS_V1 format (open, read_file, apply, save, the
// engine-checkpoint blocks, write_file, CRC32) and save_legacy_routing. No gate deck uses a
// native file; initialize() reports reaching that path as a stub (see engine_initialize.rs).
//
// Differences of form, never in result:
//   * The thread_local `tl_last_io_error` is a thread_local RefCell<Bytes>.
//   * The std::ifstream is a byte buffer read whole at open with a cursor that follows
//     libstdc++'s rules for what this reader does: read() of n bytes past the end copies what
//     is left and sets fail and eof; a failed stream reads nothing; seekg(off, cur) on a good
//     stream may move past the end without failing (a filebuf lseek); clear() resets the
//     state. A directory opens and reads as an empty file, as std::ifstream does (U53).
//   * The 2D-only refusals (layeredGroundwaterRestart, coupledThermalRestartUnsupported) are
//     false without OPENSWMM_HAS_2D, as upstream's are.

use std::cell::RefCell;

use crate::core::context::SimulationContext;
use crate::core::file_io::path_from_bytes;
use crate::core::options::Bytes;
use crate::data::lid_node::LidNodeBottom;
use crate::data::node::NodeType;

thread_local! {
    static TL_LAST_IO_ERROR: RefCell<Bytes> = const { RefCell::new(Vec::new()) };
}

/// HotStartManager::last_io_error(): the message of the last failed hot-start I/O call.
pub fn last_io_error() -> Bytes {
    TL_LAST_IO_ERROR.with(|e| e.borrow().clone())
}

pub(crate) fn set_last_io_error(m: Bytes) {
    set_io_error(m)
}

fn set_io_error(m: Bytes) {
    TL_LAST_IO_ERROR.with(|e| *e.borrow_mut() = m);
}

fn cat(parts: &[&[u8]]) -> Bytes {
    let mut m = Vec::new();
    for p in parts {
        m.extend_from_slice(p);
    }
    m
}

/// LegacySpeciesTrailer: the reaction species block a v6 engine appends to a legacy file.
#[derive(Debug, Clone, PartialEq)]
pub struct LegacySpeciesTrailer {
    pub species: Vec<Bytes>,
    pub n_subcatch: i32,
    /// -1: the file has no surface block
    pub n_landuses: i32,
    /// [sc * ns + s] washoff concentration
    pub washoff: Vec<f64>,
    /// [sc * ns + s] ponded mass
    pub ponded: Vec<f64>,
    /// [(sc * nlu + lu) * ns + s], per normalizer unit
    pub buildup: Vec<f64>,
}

impl Default for LegacySpeciesTrailer {
    fn default() -> Self {
        LegacySpeciesTrailer { species: Vec::new(), n_subcatch: 0, n_landuses: -1, washoff: Vec::new(), ponded: Vec::new(), buildup: Vec::new() }
    }
}

/// The std::ifstream this file reads through (see the header note).
struct Stream {
    data: Vec<u8>,
    pos: usize,
    fail: bool,
    eof: bool,
}

impl Stream {
    /// std::ifstream(path, binary): None when the open fails.
    fn open(path: &[u8]) -> Option<Stream> {
        let p = path_from_bytes(path);
        let mut f = std::fs::File::open(&p).ok()?;
        let mut data = Vec::new();
        // a directory opens; its reads fail (U53): an empty stream
        if std::io::Read::read_to_end(&mut f, &mut data).is_err() {
            data.clear();
        }
        Some(Stream { data, pos: 0, fail: false, eof: false })
    }
    fn good(&self) -> bool {
        !self.fail && !self.eof
    }
    fn ok(&self) -> bool {
        !self.fail
    }
    /// istream::read: Some(bytes) when all n were extracted.
    fn read(&mut self, n: usize) -> Option<Vec<u8>> {
        if !self.good() {
            self.fail = true;
            return None;
        }
        let avail = self.data.len().saturating_sub(self.pos);
        if avail < n {
            self.pos = self.data.len();
            self.fail = true;
            self.eof = true;
            return None;
        }
        let v = self.data[self.pos..self.pos + n].to_vec();
        self.pos += n;
        Some(v)
    }
    /// istream::read into a prefix: the bytes extracted (all n, or what was left).
    fn read_partial(&mut self, n: usize) -> Vec<u8> {
        if !self.good() {
            self.fail = true;
            return Vec::new();
        }
        let avail = self.data.len().saturating_sub(self.pos);
        let k = avail.min(n);
        let v = self.data[self.pos..self.pos + k].to_vec();
        self.pos += k;
        if k < n {
            self.fail = true;
            self.eof = true;
        }
        v
    }
    /// seekg(off, ios::cur) (C++11: clears eofbit first; does nothing on a failed stream).
    fn seek_cur(&mut self, off: i64) {
        self.eof = false;
        if self.fail {
            return;
        }
        let np = self.pos as i64 + off;
        if np < 0 {
            self.fail = true;
        } else {
            self.pos = np as usize;
        }
    }
    fn clear(&mut self) {
        self.fail = false;
        self.eof = false;
    }
    /// read_pod<T>: is.read(sizeof T); return is.good()
    fn read_i32(&mut self, v: &mut i32) -> bool {
        if let Some(b) = self.read(4) {
            *v = i32::from_le_bytes([b[0], b[1], b[2], b[3]]);
        }
        self.good()
    }
    fn read_u32(&mut self, v: &mut u32) -> bool {
        if let Some(b) = self.read(4) {
            *v = u32::from_le_bytes([b[0], b[1], b[2], b[3]]);
        }
        self.good()
    }
    fn read_f32(&mut self, v: &mut f32) -> bool {
        if let Some(b) = self.read(4) {
            *v = f32::from_le_bytes([b[0], b[1], b[2], b[3]]);
        }
        self.good()
    }
    /// is.read(reinterpret_cast<char*>(v.data()), n * 8): the doubles read in place (a
    /// short read fills the prefix, as the raw read does).
    fn read_f64s(&mut self, v: &mut [f64]) {
        let b = self.read_partial(v.len() * 8);
        for (k, c) in b.chunks(8).enumerate() {
            if c.len() == 8 {
                v[k] = f64::from_le_bytes([c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7]]);
            } else {
                // a partial double: its low bytes are overwritten, the rest kept
                let mut cur = v[k].to_le_bytes();
                cur[..c.len()].copy_from_slice(c);
                v[k] = f64::from_le_bytes(cur);
            }
        }
    }
    /// read_string: a uint32 length (1..4096) then the bytes, the last one dropped.
    fn read_string(&mut self, s: &mut Bytes) -> bool {
        let mut len: u32 = 0;
        if !self.read_u32(&mut len) {
            return false;
        }
        if len == 0 || len > 4096 {
            return false;
        }
        let b = self.read_partial(len as usize);
        if !self.good() {
            return false;
        }
        *s = b[..(len - 1) as usize].to_vec();
        true
    }
}

/// HotStartManager::is_native_file: the file begins with the 16-byte native magic.
pub fn is_native_file(path: &[u8]) -> bool {
    let mut s = match Stream::open(path) {
        Some(s) => s,
        None => return false,
    };
    let magic = s.read_partial(16);
    magic.len() == 16 && &magic[..15] == b"OPENSWMM_HS_V1\0"
}

/// meshLidRestartUnsupported: a storage LID on the mesh aquifer.
pub(crate) fn mesh_lid_restart_unsupported(ctx: &SimulationContext) -> bool {
    ctx.node_subtypes.storages.lid.iter().any(|c| c.control >= 0 && c.bottom == LidNodeBottom::MeshAquifer)
}

/// HotStartManager::legacy_runoff_size: doubles in the legacy runoff section for this model.
pub fn legacy_runoff_size(ctx: &SimulationContext, n_land: i32) -> usize {
    let sub = &ctx.subcatches;
    let np = ctx.n_pollutants().max(0) as usize;
    let mut n: usize = 0;
    for i in 0..ctx.n_subcatches() {
        let ui = i as usize;
        n += 4 + 6;
        if ui < sub.gw_aquifer.len() && sub.gw_aquifer[ui] >= 0 {
            n += 4;
        }
        if ui < sub.snowpack.len() && sub.snowpack[ui] >= 0 {
            n += 15;
        }
        if np > 0 {
            n += 2 * np + (n_land.max(0) as usize) * (np + 1);
        }
    }
    n
}

const K_SPECIES_TRAILER_TAG: [u8; 15] = *b"OPENSWMM-HSMSX1";

/// HotStartManager::apply_legacy_routing: read a legacy hot-start file into the context
/// (node depths and lateral flows, storage HRT, link flows, depths and settings, qualities);
/// the runoff section goes to `runoff` and the species trailer to `species`.
/// 0 on success; 1 I/O or format, 2 not a legacy file, 3 count mismatch, -1 thermal.
pub fn apply_legacy_routing(
    path: &[u8],
    ctx: &mut SimulationContext,
    warn_cb: &mut dyn FnMut(&mut SimulationContext, Bytes),
    runoff: Option<&mut Vec<f64>>,
    species: Option<&mut LegacySpeciesTrailer>,
) -> i32 {
    if !ctx.reactions.lid_controls.is_empty() {
        set_io_error(b"LID reaction inventories require a full native engine checkpoint".to_vec());
        return 1;
    }
    // layeredGroundwaterRestart / coupledThermalRestartUnsupported: false without 2D
    if mesh_lid_restart_unsupported(ctx) {
        set_io_error(b"Mesh LID hotstart requires the R5 joint owner restart extension".to_vec());
        return 1;
    }
    let mut file = match Stream::open(path) {
        Some(f) => f,
        None => {
            set_io_error(cat(&[b"Cannot open hotstart file '", path, b"'"]));
            return 1;
        }
    };

    // --- detect version from the file stamp (hotstart_is_valid) ---
    let stamp = file.read_partial(15);
    if !file.ok() {
        set_io_error(cat(&[b"Hotstart file '", path, b"' too small for header"]));
        return 1;
    }
    let version;
    if stamp.as_slice() == b"SWMM5-HOTSTART4" {
        version = 4;
    } else if stamp.as_slice() == b"SWMM5-HOTSTART3" {
        version = 3;
    } else if stamp.as_slice() == b"SWMM5-HOTSTART2" {
        version = 2;
    } else {
        file.clear();
        file.pos = 0; // seekg(0, beg)
        let s14 = file.read_partial(14);
        if file.ok() && s14.as_slice() == b"SWMM5-HOTSTART" {
            version = 1;
        } else {
            set_io_error(cat(&[b"'", path, b"' is not a legacy EPA SWMM5 hotstart file"]));
            return 2;
        }
    }

    // --- header object counts (legacy initializeFromHotstartFile) ---
    let (mut n_sub, mut n_land, mut n_nodes, mut n_links, mut n_pollut, mut flow_units) = (0i32, 0i32, 0i32, 0i32, 0i32, 0i32);
    if version >= 2 {
        if !file.read_i32(&mut n_sub) {
            return 1;
        }
    } else {
        n_sub = ctx.n_subcatches();
    }
    if version >= 3 && !file.read_i32(&mut n_land) {
        return 1;
    }
    if !file.read_i32(&mut n_nodes) {
        return 1;
    }
    if !file.read_i32(&mut n_links) {
        return 1;
    }
    if !file.read_i32(&mut n_pollut) {
        return 1;
    }
    if !file.read_i32(&mut flow_units) {
        return 1;
    }

    if n_nodes != ctx.n_nodes() || n_links != ctx.n_links() {
        set_io_error(cat(&[
            b"Hotstart node/link count mismatch (file ",
            n_nodes.to_string().as_bytes(),
            b"/",
            n_links.to_string().as_bytes(),
            b" vs model ",
            ctx.n_nodes().to_string().as_bytes(),
            b"/",
            ctx.n_links().to_string().as_bytes(),
            b")",
        ]));
        return 3;
    }
    if n_pollut != ctx.n_pollutants() {
        set_io_error(b"Hotstart pollutant count mismatch".to_vec());
        return 3;
    }
    if n_sub != ctx.n_subcatches() {
        set_io_error(cat(&[
            b"Hotstart subcatchment count mismatch (file ",
            n_sub.to_string().as_bytes(),
            b" vs model ",
            ctx.n_subcatches().to_string().as_bytes(),
            b")",
        ]));
        return 3;
    }
    let has_runoff = runoff.is_some();
    if has_runoff && version >= 3 && n_land != ctx.n_landuses() {
        set_io_error(cat(&[
            b"Hotstart land use count mismatch (file ",
            n_land.to_string().as_bytes(),
            b" vs model ",
            ctx.n_landuses().to_string().as_bytes(),
            b")",
        ]));
        return 3;
    }

    if let Some(runoff) = runoff {
        runoff.clear();
        if version >= 3 && n_sub > 0 {
            runoff.resize(legacy_runoff_size(ctx, n_land), 0.0);
            file.read_f64s(runoff);
            if !file.ok() {
                set_io_error(cat(&[b"Hotstart runoff section truncated in '", path, b"'"]));
                return 1;
            }
        }
    }

    if n_sub > 0 && (!has_runoff || version < 3) {
        warn_cb(
            ctx,
            b"USE HOTSTART: subcatchment runoff/infiltration/groundwater/snowpack state was skipped (routing state applied only)".to_vec(),
        );
    }

    if !has_runoff && version >= 3 && n_sub > 0 {
        let sub = &ctx.subcatches;
        let mut skip_doubles: i64 = 0;
        for i in 0..n_sub {
            let ui = i as usize;
            let mut nd: i64 = 4 + 6;
            if ui < sub.gw_aquifer.len() && sub.gw_aquifer[ui] >= 0 {
                nd += 4;
            }
            if ui < sub.snowpack.len() && sub.snowpack[ui] >= 0 {
                nd += 15;
            }
            if n_pollut > 0 {
                nd += 2 * n_pollut as i64;
                nd += n_land as i64 * (n_pollut as i64 + 1);
            }
            skip_doubles += nd;
        }
        file.seek_cur(skip_doubles * 8);
        if !file.ok() {
            set_io_error(cat(&[b"Hotstart runoff section truncated in '", path, b"'"]));
            return 1;
        }
    } else if version == 2 && n_sub > 0 {
        file.seek_cur(n_sub as i64 * 2 * 4);
        if !file.ok() {
            set_io_error(cat(&[b"Hotstart v2 groundwater block truncated in '", path, b"'"]));
            return 1;
        }
    }

    // --- node states (legacy readRouting) ---
    for i in 0..n_nodes {
        let ui = i as usize;
        let (mut depth, mut lat) = (0.0f32, 0.0f32);
        if !file.read_f32(&mut depth) {
            return 1;
        }
        if !file.read_f32(&mut lat) {
            return 1;
        }
        ctx.nodes.depth[ui] = depth as f64;
        ctx.nodes.lat_flow[ui] = lat as f64;
        if version >= 4 && ctx.nodes.type_[ui] == NodeType::Storage {
            let mut hrt = 0.0f32;
            if !file.read_f32(&mut hrt) {
                return 1;
            }
            if ui < ctx.nodes.hrt.len() {
                ctx.nodes.hrt[ui] = hrt as f64;
            }
        }
        for j in 0..n_pollut {
            let mut q = 0.0f32;
            if !file.read_f32(&mut q) {
                return 1;
            }
            let qi = ui * n_pollut as usize + j as usize;
            ctx.nodes.conc_old[qi] = q as f64;
            ctx.nodes.conc[qi] = q as f64;
        }
        if version <= 2 {
            for _ in 0..n_pollut {
                let mut q = 0.0f32;
                if !file.read_f32(&mut q) {
                    return 1;
                }
            }
        }
    }

    // --- link states (legacy readRouting) ---
    for i in 0..n_links {
        let ui = i as usize;
        let (mut flow, mut depth, mut setting) = (0.0f32, 0.0f32, 0.0f32);
        if !file.read_f32(&mut flow) {
            return 1;
        }
        if !file.read_f32(&mut depth) {
            return 1;
        }
        if !file.read_f32(&mut setting) {
            return 1;
        }
        ctx.links.flow[ui] = flow as f64;
        ctx.links.depth[ui] = depth as f64;
        ctx.links.setting[ui] = setting as f64;
        ctx.links.target_setting[ui] = setting as f64;
        for j in 0..n_pollut {
            let mut q = 0.0f32;
            if !file.read_f32(&mut q) {
                return 1;
            }
            let qi = ui * n_pollut as usize + j as usize;
            ctx.links.conc_old[qi] = q as f64;
            ctx.links.conc[qi] = q as f64;
        }
    }

    // --- reaction species trailer (v6 only; legacy stops reading above) ---
    if let Some(species) = species {
        *species = LegacySpeciesTrailer::default();
        let tag = file.read_partial(15);
        if file.ok() && tag.as_slice() == K_SPECIES_TRAILER_TAG {
            let (mut ns, mut nn, mut nl, mut surf) = (0i32, 0i32, 0i32, 0i32);
            let mut ok = file.read_i32(&mut ns) && ns >= 0 && ns < 100000;
            let mut names: Vec<Bytes> = vec![Vec::new(); if ok { ns as usize } else { 0 }];
            for n in names.iter_mut() {
                ok = ok && file.read_string(n);
            }
            ok = ok && file.read_i32(&mut nn) && file.read_i32(&mut nl);
            ok = ok && nn == n_nodes && nl == n_links;
            let uns = ns.max(0) as usize;
            let mut node = vec![0.0f64; if ok { nn as usize * uns } else { 0 }];
            let mut link = vec![0.0f64; if ok { nl as usize * uns } else { 0 }];
            if ok {
                file.read_f64s(&mut node);
                file.read_f64s(&mut link);
                ok = file.ok() && file.read_i32(&mut surf);
            }
            if ok && surf == 1 {
                let (mut nsc, mut nlu) = (0i32, 0i32);
                ok = file.read_i32(&mut nsc) && file.read_i32(&mut nlu) && nsc >= 0 && nlu >= 0;
                if ok {
                    species.n_subcatch = nsc;
                    species.n_landuses = nlu;
                    let a = nsc as usize * uns;
                    species.washoff.resize(a, 0.0);
                    species.ponded.resize(a, 0.0);
                    species.buildup.resize(a * nlu as usize, 0.0);
                    file.read_f64s(&mut species.washoff);
                    file.read_f64s(&mut species.ponded);
                    file.read_f64s(&mut species.buildup);
                    ok = file.ok();
                }
            }
            if !ok {
                set_io_error(cat(&[b"Hotstart reaction species block is malformed in '", path, b"'"]));
                return 1;
            }
            species.species = names.clone();

            // Network concentrations, by name, before init_modules().
            let nm = if ctx.reactions.configured { ctx.reactions.n_species() } else { 0 };
            let unm = nm as usize;
            let mut col = vec![-1i32; uns];
            for s in 0..uns {
                col[s] = if nm > 0 { ctx.reactions.find_species(&names[s]) } else { -1 };
                if col[s] < 0 {
                    warn_cb(ctx, cat(&[b"USE HOTSTART: reaction species '", &names[s], b"' is not in the current model; its state is skipped"]));
                }
            }
            if nm > 0 {
                crate::transport::reaction_legacy_binding::ensure_msx_state(ctx);
                let rx = &mut ctx.reactions;
                for s in 0..uns {
                    if col[s] < 0 {
                        continue;
                    }
                    let m = col[s] as usize;
                    for e in 0..nn as usize {
                        if e * unm + m < rx.msx_node_conc.len() {
                            rx.msx_node_conc[e * unm + m] = node[e * uns + s];
                        }
                    }
                    for e in 0..nl as usize {
                        if e * unm + m < rx.msx_link_conc.len() {
                            rx.msx_link_conc[e * unm + m] = link[e * uns + s];
                        }
                    }
                }
            }
        } else if file.ok() {
            warn_cb(ctx, b"USE HOTSTART: unrecognised data after the routing records was ignored".to_vec());
        }
    }

    let _ = flow_units;
    let _ = n_land;
    set_io_error(Vec::new());
    0
}

/// HotStartManager::apply_species_surface: install the trailer's surface block in the
/// species store; false when it has none or the model's counts differ.
pub fn apply_species_surface(tr: &LegacySpeciesTrailer, ctx: &mut SimulationContext, warn: &mut dyn FnMut(&mut SimulationContext, Bytes)) -> bool {
    if tr.n_landuses < 0 || !ctx.reactions.surface.active() {
        return false;
    }
    if tr.n_subcatch != ctx.n_subcatches() || tr.n_landuses != ctx.reactions.surface.n_landuses {
        warn(ctx, b"USE HOTSTART: reaction species surface state skipped (subcatchment or land use count differs from the model)".to_vec());
        return false;
    }
    let uns = tr.species.len();
    let nlu = tr.n_landuses;
    for k in 0..uns {
        let m = ctx.reactions.find_species(&tr.species[k]);
        let s = &mut ctx.reactions.surface;
        if m < 0 || m >= s.n_species {
            continue;
        }
        for sc in 0..tr.n_subcatch {
            let src = sc as usize * uns + k;
            let dst = s.sidx(sc, m);
            if dst < s.washoff_conc.len() {
                s.washoff_conc[dst] = tr.washoff[src];
            }
            if dst < s.washoff_conc_old.len() {
                s.washoff_conc_old[dst] = tr.washoff[src];
            }
            if dst < s.ponded_qual.len() {
                s.ponded_qual[dst] = tr.ponded[src];
            }
            for lu in 0..nlu {
                let b = (sc as usize * nlu as usize + lu as usize) * uns + k;
                let d = s.bu_idx(sc, lu, m);
                if d < s.buildup.len() {
                    s.buildup[d] = tr.buildup[b];
                }
            }
        }
    }
    true
}
