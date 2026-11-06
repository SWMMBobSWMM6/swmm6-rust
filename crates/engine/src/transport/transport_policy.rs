// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/TransportPolicy.hpp and TransportPolicy.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The Domain x Species transport matrix (resolve) and its .rpt block (formatReportBlock), as
// the report preamble prints them.
//
// Differences of form, never in result:
//   * 1D build: ctx.twod_io's mesh, options, aquifer, aquifer_state and gw pointers are null
//     (no 2D router attaches them), so hasMesh2D is false, opts2D null, the groundwater cell
//     takes its legacy branch and the 2D surface cells read "no 2D mesh". The Draft argument
//     (a planning pass of the C API) is not used by the report: resolve(ctx) takes none.
//     groundwaterEnables and surface2DEnables are not on the report path; not translated.

use crate::core::context::SimulationContext;
use crate::core::options::{Bytes, QualitySolverKind};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum Domain {
    Runoff = 0,
    Groundwater = 1,
    Network1D = 2,
    Surface2D = 3,
}
const DOMAIN_COUNT: usize = 4;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(i32)]
pub enum SpeciesClass {
    Pollutants = 0,
    Msx = 1,
    Age = 2,
    Temperature = 3,
}
const SPECIES_CLASS_COUNT: usize = 4;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[repr(i32)]
pub enum CellState {
    Enabled = 0,
    DisabledByUser = 1,
    #[default]
    Unavailable = 2,
}

#[derive(Debug, Clone, Default)]
pub struct Cell {
    pub state: CellState,
    pub count: i32,
    pub reason: Bytes,
}

#[derive(Debug, Clone, Default)]
pub struct Matrix {
    pub cells: [[Cell; SPECIES_CLASS_COUNT]; DOMAIN_COUNT],
}

impl Matrix {
    fn at(&mut self, d: Domain, c: SpeciesClass) -> &mut Cell {
        &mut self.cells[d as usize][c as usize]
    }
}

fn reactions_active(ctx: &SimulationContext) -> bool {
    ctx.reactions.configured && ctx.reactions.compiled
}

fn reactions_have_wall(ctx: &SimulationContext) -> bool {
    ctx.reactions.species_is_wall.iter().any(|&w| w != 0)
}

/// ClassEnables: the species rows a 1D or 2D transport owner carries.
#[derive(Debug, Clone, Copy, Default)]
pub struct ClassEnables {
    pub n_pollut: i32,
    pub n_msx: i32,
    pub age: bool,
    pub temperature: bool,
    pub msx_has_wall: bool,
    pub msx_bulk_only: bool,
}

impl ClassEnables {
    pub fn total(&self) -> i32 {
        self.n_pollut + self.n_msx + (if self.age { 1 } else { 0 }) + (if self.temperature { 1 } else { 0 })
    }
}

/// network1DEnables: the 1D network enables (LEGACY, ARD and LARD size from this).
pub fn network_1d_enables(ctx: &SimulationContext) -> ClassEnables {
    let q = !ctx.options.ignore_quality;
    ClassEnables {
        n_pollut: if q { ctx.n_pollutants() } else { 0 },
        msx_has_wall: reactions_have_wall(ctx),
        n_msx: if q && reactions_active(ctx) { ctx.reactions.n_species() } else { 0 },
        age: q && ctx.options.water_age,
        temperature: q && ctx.options.heat_transport,
        msx_bulk_only: false,
    }
}

fn enabled(c: &mut Cell, count: i32) {
    c.state = CellState::Enabled;
    c.count = count;
    c.reason.clear();
}
fn disabled(c: &mut Cell, key: &[u8]) {
    c.state = CellState::DisabledByUser;
    c.count = 0;
    c.reason = key.to_vec();
}
fn unavailable(c: &mut Cell, why: &[u8]) {
    c.state = CellState::Unavailable;
    c.count = 0;
    c.reason = why.to_vec();
}

/// resolve(ctx): the matrix for the project as configured (no draft).
pub fn resolve(ctx: &SimulationContext) -> Matrix {
    let mut m = Matrix::default();
    let np = ctx.n_pollutants();
    let ignore_q = ctx.options.ignore_quality;
    let msx_on = reactions_active(ctx);
    let nm = if msx_on { ctx.reactions.n_species() } else { 0 };
    let msx_wall = msx_on && reactions_have_wall(ctx);
    let age_on = ctx.options.water_age;
    let heat_on = ctx.options.heat_transport;

    // ---- Runoff / LID ----
    {
        let has_sub = ctx.n_subcatches() > 0;
        let pol = m.at(Domain::Runoff, SpeciesClass::Pollutants);
        if !has_sub {
            unavailable(pol, b"no subcatchments");
        } else if np == 0 {
            unavailable(pol, b"no [POLLUTANTS]");
        } else if ignore_q {
            disabled(pol, b"IGNORE_QUALITY");
        } else {
            enabled(pol, np);
        }

        let msx = m.at(Domain::Runoff, SpeciesClass::Msx);
        if !msx_on {
            unavailable(msx, b"no reactions component");
        } else if !has_sub {
            unavailable(msx, b"no subcatchments");
        } else if ignore_q {
            disabled(msx, b"IGNORE_QUALITY");
        } else if !ctx.reactions.surface.active() {
            unavailable(msx, b"no species buildup, washoff or rainfall rows");
        } else {
            enabled(msx, nm);
            msx.reason = b"loading and transport; no surface reactions".to_vec();
        }

        let age = m.at(Domain::Runoff, SpeciesClass::Age);
        if !has_sub {
            unavailable(age, b"no subcatchments");
        } else if !age_on {
            unavailable(age, b"WATER_AGE OFF");
        } else if ignore_q {
            disabled(age, b"IGNORE_QUALITY");
        } else {
            enabled(age, 1);
        }

        let tmp = m.at(Domain::Runoff, SpeciesClass::Temperature);
        if !has_sub {
            unavailable(tmp, b"no subcatchments");
        } else if !heat_on {
            unavailable(tmp, b"HEAT_TRANSPORT OFF");
        } else if ignore_q {
            disabled(tmp, b"IGNORE_QUALITY");
        } else {
            enabled(tmp, 1);
        }
    }

    // ---- Groundwater (no two-zone [2D_AQUIFER] kernel in a 1D build) ----
    {
        let has_aq = ctx.aquifer_names.size() > 0;
        let no_aq: &[u8] = b"no [AQUIFERS] or [2D_AQUIFER]";
        unavailable(m.at(Domain::Groundwater, SpeciesClass::Pollutants), if has_aq { b"no transported quality \xE2\x80\x94 the aquifer supplies the [POLLUTANTS] GW concentration at the node seam (legacy)" as &[u8] } else { no_aq });
        unavailable(m.at(Domain::Groundwater, SpeciesClass::Msx), if has_aq { b"no transported quality in groundwater" as &[u8] } else { no_aq });
        unavailable(m.at(Domain::Groundwater, SpeciesClass::Age), if has_aq { b"source volume only (new water enters at age 0)" as &[u8] } else { no_aq });
        unavailable(m.at(Domain::Groundwater, SpeciesClass::Temperature), if has_aq { b"source volume only (GROUNDWATER source temperature)" as &[u8] } else { no_aq });
    }

    // ---- 1D network ----
    {
        let pol = m.at(Domain::Network1D, SpeciesClass::Pollutants);
        if np == 0 {
            unavailable(pol, b"no [POLLUTANTS]");
        } else if ignore_q {
            disabled(pol, b"IGNORE_QUALITY");
        } else {
            enabled(pol, np);
        }

        let msx = m.at(Domain::Network1D, SpeciesClass::Msx);
        if !msx_on {
            unavailable(msx, b"no reactions component");
        } else if ignore_q {
            disabled(msx, b"IGNORE_QUALITY");
        } else {
            enabled(msx, nm);
        }
        if msx.state == CellState::Enabled && msx_wall && ctx.options.quality_solver == QualitySolverKind::EulerianArd {
            msx.reason = b"WALL species: EULERIAN_ARD falls back to LEGACY".to_vec();
        }

        let age = m.at(Domain::Network1D, SpeciesClass::Age);
        if !age_on {
            unavailable(age, b"WATER_AGE OFF");
        } else if ignore_q {
            disabled(age, b"IGNORE_QUALITY");
        } else {
            enabled(age, 1);
        }

        let tmp = m.at(Domain::Network1D, SpeciesClass::Temperature);
        if !heat_on {
            unavailable(tmp, b"HEAT_TRANSPORT OFF");
        } else if ignore_q {
            disabled(tmp, b"IGNORE_QUALITY");
        } else {
            enabled(tmp, 1);
        }
    }

    // ---- 2D surface (no mesh in a 1D build) ----
    for c in [SpeciesClass::Pollutants, SpeciesClass::Msx, SpeciesClass::Age, SpeciesClass::Temperature] {
        unavailable(m.at(Domain::Surface2D, c), b"no 2D mesh");
    }

    m
}

fn domain_name(d: usize) -> &'static [u8] {
    match d {
        0 => b"Runoff / LID",
        1 => b"Groundwater",
        2 => b"1D network",
        3 => b"2D surface",
        _ => b"?",
    }
}

fn species_class_name(c: usize) -> &'static [u8] {
    match c {
        0 => b"Pollutants",
        1 => b"MSX",
        2 => b"Age",
        3 => b"Temperature",
        _ => b"?",
    }
}

/// snprintf into a char[256]: at most 255 bytes.
fn buf256(mut b: Bytes) -> Bytes {
    b.truncate(255);
    b
}

/// formatReportBlock(m).
pub fn format_report_block(m: &Matrix) -> Bytes {
    let mut out: Bytes = Vec::new();
    out.extend_from_slice(b"  Transport by domain (Domain x Species; on(n) = rows carried,\n");
    out.extend_from_slice(b"  off:KEY = switched off by that option, n/a = not carried):\n");
    out.extend(buf256(crate::cformat!(b"  %-14s %-14s %-14s %-14s %-14s\n", b"", b"Pollutants", b"MSX", b"Age", b"Temperature")));
    for d in 0..DOMAIN_COUNT {
        let mut line: Bytes = b"  ".to_vec();
        line.extend(buf256(crate::cformat!(b"%-14s ", domain_name(d))));
        for c in 0..SPECIES_CLASS_COUNT {
            let cell = &m.cells[d][c];
            let tok: Bytes = match cell.state {
                CellState::Enabled => {
                    let mut t = b"on(".to_vec();
                    t.extend_from_slice(cell.count.to_string().as_bytes());
                    t.push(b')');
                    t
                }
                CellState::DisabledByUser => {
                    let mut t = b"off:".to_vec();
                    t.extend_from_slice(&cell.reason);
                    t
                }
                _ => b"n/a".to_vec(),
            };
            line.extend(buf256(crate::cformat!(b"%-14s ", tok)));
        }
        out.extend_from_slice(&line);
        out.push(b'\n');
    }
    let mut any = false;
    for d in 0..DOMAIN_COUNT {
        for c in 0..SPECIES_CLASS_COUNT {
            let cell = &m.cells[d][c];
            if cell.state != CellState::Unavailable || cell.reason.is_empty() {
                continue;
            }
            if !any {
                out.extend_from_slice(b"  n/a because:\n");
                any = true;
            }
            out.extend(buf256(crate::cformat!(b"    %-14s %-12s %s\n", domain_name(d), species_class_name(c), cell.reason)));
        }
    }
    out
}

/// formatReportBlock(resolve(ctx)), as the .rpt preamble prints it.
pub fn report_block(ctx: &SimulationContext) -> Bytes {
    format_report_block(&resolve(ctx))
}
