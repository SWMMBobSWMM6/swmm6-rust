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
//     groundwaterEnables is not on the report path; not translated.
//   * With the `twod` feature (12 November 2026) the groundwater and 2D surface cells read the
//     router's mesh, options, [2D_AQUIFER] rows and [GW_*] options through ctx.twod_io, as
//     upstream's 2D build does; see groundwater_two_zone for the kernel state. surface2DEnables
//     and the canonical row layout (canonicalRows) are translated for the 2D router's
//     initialize.

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

/// surface2DEnables(ctx): the 2D surface's class enables (twod feature). Upstream reads the
/// [2D_OPTIONS] TRANSPORT_* keys through ctx.twod_io.options; SurfaceRouter2D::initialize holds
/// that store borrowed, so it is passed here (None reads as upstream's null options: every class
/// inherits on).
#[cfg(feature = "twod")]
pub fn surface_2d_enables(ctx: &SimulationContext, o: Option<&crate::twod::data::solver_options_2d::SolverOptions2D>) -> ClassEnables {
    let q = !ctx.options.ignore_quality;
    let tp = o.map_or(true, |o| o.transport_pollutants);
    let tm = o.map_or(true, |o| o.transport_msx);
    let ta = o.map_or(true, |o| o.transport_age);
    let tt = o.map_or(true, |o| o.transport_temperature);
    let msx_has_wall = reactions_have_wall(ctx);
    ClassEnables {
        n_pollut: if q && tp { ctx.n_pollutants() } else { 0 },
        msx_has_wall,
        // WALL species have no surface transport semantics.
        n_msx: if q && tm && reactions_active(ctx) && !msx_has_wall { ctx.reactions.n_species() } else { 0 },
        age: ta && ctx.options.water_age,
        temperature: tt && ctx.options.heat_transport,
        msx_bulk_only: false,
    }
}

/// RowLayout: a canonical transported-row layout.
#[derive(Debug, Clone, Default)]
pub struct RowLayout {
    pub n_pollut: i32,
    pub n_msx: i32,
    pub age_row: i32,
    pub temp_row: i32,
    pub ns: i32,
    pub names: Vec<crate::core::options::Bytes>,
}

/// canonicalRows(ctx, e): pollutants, MSX species (BULK only when msx_bulk_only), age,
/// temperature.
pub fn canonical_rows(ctx: &SimulationContext, e: &ClassEnables) -> RowLayout {
    let mut l = RowLayout::default();
    l.n_pollut = e.n_pollut;
    l.n_msx = e.n_msx;
    l.age_row = if e.age { l.n_pollut + l.n_msx } else { -1 };
    l.temp_row = if e.temperature { l.n_pollut + l.n_msx + (if e.age { 1 } else { 0 }) } else { -1 };
    l.ns = e.total();
    l.names.reserve(l.ns.max(0) as usize);
    for p in 0..l.n_pollut {
        l.names.push(ctx.pollutant_names.name_of(p).to_vec());
    }
    let mut added = 0;
    let mut m = 0;
    while added < l.n_msx && m < ctx.reactions.n_species() {
        if e.msx_bulk_only && ctx.reactions.species_is_wall[m as usize] != 0 {
            m += 1;
            continue;
        }
        l.names.push(ctx.reactions.species_name[m as usize].clone());
        added += 1;
        m += 1;
    }
    if e.age {
        l.names.push(b"__WATER_AGE__".to_vec());
    }
    if e.temperature {
        l.names.push(b"__TEMPERATURE__".to_vec());
    }
    l
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

/// hasMesh2D(): a mesh with at least one cell (twod feature).
#[cfg(feature = "twod")]
fn has_mesh_2d(ctx: &SimulationContext) -> bool {
    ctx.twod_io.mesh.as_ref().is_some_and(|m| m.borrow().n_cells() > 0)
}

/// The groundwater cell's two-zone test (twod feature). Upstream: `authored = (!aq ||
/// aq->n_cells == 0) && hasMesh2D && aquifer && !aquifer->empty() && !ignore_2d &&
/// (!options || options->groundwater != 0)`, `two_zone = authored || (aq && aq->active &&
/// aq->n_cells > 0)`. The kernel state `aq` (twod_io.aquifer_state) belongs to phase 2b and is
/// never active here, as it is not upstream before SurfaceRouter2D::initialize resolves an
/// aquifer; so two_zone is `authored` and the counts are the planned ones.
#[cfg(feature = "twod")]
fn groundwater_two_zone(ctx: &SimulationContext) -> bool {
    let io = &ctx.twod_io;
    has_mesh_2d(ctx)
        && io.aquifer.as_ref().is_some_and(|a| !a.borrow().empty())
        && !ctx.options.ignore_2d
        && io.options.as_ref().map_or(true, |o| o.borrow().groundwater != 0)
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

    // ---- Groundwater: the two-zone [2D_AQUIFER] kernel (twod feature, below) ----
    #[cfg(feature = "twod")]
    let two_zone = groundwater_two_zone(ctx);
    #[cfg(feature = "twod")]
    if two_zone {
        // planned (authored) counts: groundwaterEnables would read an active kernel, which
        // the translation does not run yet (see groundwater_two_zone)
        let g = ctx.twod_io.gw.as_ref().map(|g| g.borrow());
        let pol = m.at(Domain::Groundwater, SpeciesClass::Pollutants);
        if np == 0 {
            unavailable(pol, b"no [POLLUTANTS]");
        } else if ignore_q {
            disabled(pol, b"IGNORE_QUALITY");
        } else if g.as_ref().is_some_and(|g| !g.options.transport_pollutants) {
            disabled(pol, b"TRANSPORT_POLLUTANTS NO");
        } else {
            enabled(pol, np);
        }
        let bulk = ctx.reactions.species_is_wall.iter().filter(|&&w| w == 0).count() as i32;
        let msx = m.at(Domain::Groundwater, SpeciesClass::Msx);
        if !reactions_active(ctx) {
            unavailable(msx, b"no reactions component");
        } else if bulk == 0 {
            unavailable(msx, b"no BULK species for subsurface transport");
        } else if ignore_q {
            disabled(msx, b"IGNORE_QUALITY");
        } else if g.as_ref().is_some_and(|g| !g.options.transport_msx) {
            disabled(msx, b"TRANSPORT_MSX NO");
        } else {
            enabled(msx, bulk);
        }
        let age = m.at(Domain::Groundwater, SpeciesClass::Age);
        if !age_on {
            unavailable(age, b"WATER_AGE OFF");
        } else if g.as_ref().is_some_and(|g| !g.options.transport_age) {
            disabled(age, b"TRANSPORT_AGE NO");
        } else {
            enabled(age, 1);
        }
        let tmp = m.at(Domain::Groundwater, SpeciesClass::Temperature);
        if !heat_on {
            unavailable(tmp, b"HEAT_TRANSPORT OFF");
        } else if g.as_ref().is_some_and(|g| !g.options.transport_temperature) {
            disabled(tmp, b"TRANSPORT_TEMPERATURE NO");
        } else {
            enabled(tmp, 1);
        }
    }
    #[cfg(not(feature = "twod"))]
    let two_zone = false;
    if !two_zone {
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
    #[cfg(not(feature = "twod"))]
    for c in [SpeciesClass::Pollutants, SpeciesClass::Msx, SpeciesClass::Age, SpeciesClass::Temperature] {
        unavailable(m.at(Domain::Surface2D, c), b"no 2D mesh");
    }
    #[cfg(feature = "twod")]
    {
        let mesh = has_mesh_2d(ctx);
        let o = ctx.twod_io.options.as_ref().map(|o| o.borrow());
        let ig2d = ctx.options.ignore_2d;
        let gate = |c: &mut Cell, project_on: bool, off_why: &[u8], key_on: bool, key: &[u8], quality_gated: bool, count: i32, unavailable_why: Option<&[u8]>| {
            if !mesh {
                unavailable(c, b"no 2D mesh");
                return;
            }
            if !project_on {
                unavailable(c, off_why);
                return;
            }
            if ig2d {
                disabled(c, b"IGNORE_2D");
                return;
            }
            if quality_gated && ignore_q {
                disabled(c, b"IGNORE_QUALITY");
                return;
            }
            if let Some(why) = unavailable_why {
                unavailable(c, why);
                return;
            }
            if !key_on {
                disabled(c, key);
                return;
            }
            enabled(c, count);
        };
        let tp = o.as_ref().map_or(true, |o| o.transport_pollutants);
        let tm = o.as_ref().map_or(true, |o| o.transport_msx);
        let ta = o.as_ref().map_or(true, |o| o.transport_age);
        let tt = o.as_ref().map_or(true, |o| o.transport_temperature);
        gate(m.at(Domain::Surface2D, SpeciesClass::Pollutants), np > 0, b"no [POLLUTANTS]", tp, b"TRANSPORT_POLLUTANTS", true, np, None);
        gate(
            m.at(Domain::Surface2D, SpeciesClass::Msx),
            msx_on,
            b"no reactions component",
            tm,
            b"TRANSPORT_MSX",
            true,
            nm,
            if msx_wall { Some(b"WALL species have no surface transport semantics" as &[u8]) } else { None },
        );
        gate(m.at(Domain::Surface2D, SpeciesClass::Age), age_on, b"WATER_AGE OFF", ta, b"TRANSPORT_AGE", false, 1, None);
        gate(m.at(Domain::Surface2D, SpeciesClass::Temperature), heat_on, b"HEAT_TRANSPORT OFF", tt, b"TRANSPORT_TEMPERATURE", false, 1, None);
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
