// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/coupling/OutfallExchange.hpp and
// OutfallExchange.cpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The outfall exchange: the donor cells, SPREAD receivers and per-node ledgers built at
// initialize (wave 3), and the routing step's reservation (prepare), its application of the
// routed outfall water to the surface cells (apply), the tailwater record (recordTailwater,
// tailwater) (wave 4b). checkpoint and validateCheckpoint are out of scope (checkpointing).
//
// Differences in form, never in result:
//   * initialize's two std::runtime_error throws are Err(what), caught by
//     SWMMEngine::initHydraulics' handler. `users.at(c)` cannot go out of range: every cell
//     comes from the mesh's own stencils.
//   * apply's three std::runtime_error throws travel up through advancePostRouting and
//     stepRouting (neither noexcept) to step()'s catch: they are Err(Terminated) whose
//     exception is "throw std::runtime_error" (core::engine_step::THROWN_PREFIX). The
//     ResetReceivers guard (recv_ cleared on every exit) is an explicit clear on each return.
//   * apply's perf::ProcessTimer (timing only) is left out.
//   * Point's `index` and `node` (no initializer upstream) default to 0; initialize sets both.
#![allow(clippy::all)]

use crate::core::cmath::{cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::engine_init_hydraulics::Terminated;
use crate::core::options::Bytes;
use crate::twod::coupling::node_coupling::CouplingPoint;
use crate::twod::data::mesh_data::MeshData;
use crate::twod::data::solver_options_2d::SolverOptions2D;
use crate::twod::data::surface_state_data::SurfaceStateData;
use crate::twod::solver::inertial_kernels::cell_eta_depth;

/// apply's std::runtime_error, on its way to step()'s catch.
fn thrown(what: &'static str) -> Terminated {
    Terminated { what, exception: "throw std::runtime_error" }
}

/// OutfallExchange::Point.
#[derive(Debug, Clone, Default)]
pub struct Point {
    pub index: i32,
    pub node: i32,
    pub cells: Vec<i32>,
    pub fraction: Vec<f64>,
    pub reserved: Vec<f64>,
    pub total: f64,
    pub area: f64,
    pub spread: Vec<i32>,
    pub bed_ref: f64,
}

#[derive(Debug, Clone, Default)]
pub struct OutfallExchange {
    pub points_: Vec<Point>,
    pub recv_: Vec<i32>,
    pub nodes_: Vec<i32>,
    pub degree_: Vec<i32>,
    pub node_area_: Vec<f64>,
    pub node_reserved_: Vec<f64>,
    pub tw_engaged_s_: Vec<f64>,
    pub tw_submerged_s_: Vec<f64>,
    pub tw_max_: Vec<f64>,
    pub pending_: Vec<f64>,
    pub prepared_dt_: f64,
    pub prepared_: bool,
}

fn assign<T: Clone>(v: &mut Vec<T>, n: usize, x: T) {
    v.clear();
    v.resize(n, x);
}

impl OutfallExchange {
    pub fn reset_window(&mut self) {
        for x in self.pending_.iter_mut() {
            *x = 0.0;
        }
    }

    /// initialize(mesh, cps, ctx).
    pub fn initialize(&mut self, mesh: &MeshData, cps: &[CouplingPoint], ctx: &mut SimulationContext) -> Result<(), Bytes> {
        self.points_.clear();
        self.recv_.clear();
        self.nodes_.clear();
        self.pending_.clear();
        self.prepared_dt_ = 0.0;
        self.prepared_ = false;
        ctx.surface_outfall_link_limit.clear();
        let nn = ctx.n_nodes() as usize;
        assign(&mut self.degree_, nn, 0);
        assign(&mut self.node_area_, nn, 0.0);
        assign(&mut self.node_reserved_, nn, 0.0);
        assign(&mut self.tw_engaged_s_, nn, 0.0);
        assign(&mut self.tw_submerged_s_, nn, 0.0);
        assign(&mut self.tw_max_, nn, -1.0e30);
        let mut users = vec![0i32; mesh.n_cells() as usize];
        let mut seen: Vec<u8> = Vec::new();
        for (k, cp) in cps.iter().enumerate() {
            if !cp.is_outfall {
                continue;
            }
            let mut p = Point::default();
            p.index = k as i32;
            p.node = cp.node_idx;
            if cp.vertex_idx < 0 {
                p.cells.push(cp.cell_idx);
            } else {
                let v = cp.vertex_idx as usize;
                for j in mesh.vert_stencil_ptr[v]..mesh.vert_stencil_ptr[v + 1] {
                    p.cells.push(mesh.vert_stencil_idx[j as usize]);
                }
            }
            if p.cells.is_empty() {
                return Err(b"2D outfall has no donor cells.".to_vec());
            }
            for &c in &p.cells {
                users[c as usize] += 1;
                p.area += mesh.tri_area[c as usize];
            }
            self.node_area_[p.node as usize] += p.area;
            p.reserved.resize(p.cells.len(), 0.0);
            // SPREAD r: the cells whose centroid lies within r of the point, gathered outward
            // from its own cells, receive its discharge.
            if cp.spread > 0.0 {
                if seen.is_empty() {
                    assign(&mut seen, mesh.n_cells() as usize, 0u8);
                }
                let cx = if cp.vertex_idx >= 0 { mesh.vx[cp.vertex_idx as usize] } else { mesh.tri_cx[cp.cell_idx as usize] };
                let cy = if cp.vertex_idx >= 0 { mesh.vy[cp.vertex_idx as usize] } else { mesh.tri_cy[cp.cell_idx as usize] };
                p.spread = p.cells.clone();
                p.bed_ref = f64::NEG_INFINITY;
                let mut marked: Vec<i32> = p.cells.clone();
                for &c in &p.cells {
                    seen[c as usize] = 1;
                    p.bed_ref = cmax(p.bed_ref, mesh.tri_cz[c as usize]);
                }
                let mut q = 0usize;
                while q < p.spread.len() {
                    let c = p.spread[q];
                    for e in 0..mesh.cell_vertex_count(c) {
                        let nb = mesh.cell_neighbour(c, e);
                        if nb < 0 || seen[nb as usize] != 0 {
                            continue;
                        }
                        seen[nb as usize] = 1;
                        marked.push(nb);
                        let dx = mesh.tri_cx[nb as usize] - cx;
                        let dy = mesh.tri_cy[nb as usize] - cy;
                        if dx * dx + dy * dy <= cp.spread * cp.spread {
                            p.spread.push(nb);
                        }
                    }
                    q += 1;
                }
                for &c in &marked {
                    seen[c as usize] = 0;
                }
            }
            if !self.nodes_.contains(&p.node) {
                self.nodes_.push(p.node);
            }
            self.points_.push(p);
        }
        if self.points_.is_empty() {
            return Ok(());
        }
        assign(&mut self.pending_, mesh.n_cells() as usize, 0.0);
        for j in 0..ctx.n_links() as usize {
            if ctx.links.node1[j] >= 0 {
                self.degree_[ctx.links.node1[j] as usize] += 1;
            }
            if ctx.links.node2[j] >= 0 {
                self.degree_[ctx.links.node2[j] as usize] += 1;
            }
        }
        for p in self.points_.iter_mut() {
            for &c in &p.cells {
                p.fraction.push(1.0 / users[c as usize] as f64);
            }
            if self.degree_[p.node as usize] == 0 {
                return Err(b"2D outfall has no incident hydraulic link.".to_vec());
            }
        }
        assign(&mut ctx.surface_outfall_link_limit, nn, f64::INFINITY);
        Ok(())
    }

    /// prepare(mesh, state, opts, ctx, dt): reserve each outfall's share of the frozen surface
    /// water for one routing step, and the per-link native flow limit of its node. The mesh is
    /// unused upstream too.
    pub fn prepare(&mut self, _mesh: &MeshData, state: &SurfaceStateData, opts: &SolverOptions2D, ctx: &mut SimulationContext, dt: f64) {
        if !(dt > 0.0) {
            return;
        }
        if !self.points_.is_empty() {
            self.prepared_dt_ = dt;
            self.prepared_ = true;
        }
        for p in &self.points_ {
            ctx.surface_outfall_link_limit[p.node as usize] = 0.0;
            self.node_reserved_[p.node as usize] = 0.0;
        }
        for p in self.points_.iter_mut() {
            p.total = 0.0;
            for j in 0..p.cells.len() {
                p.reserved[j] = cmax(0.0, state.volume[p.cells[j] as usize]) * p.fraction[j];
                p.total += p.reserved[j];
            }
            self.node_reserved_[p.node as usize] += p.total;
            // Conservative equal shares across incident links (aggregate-barrel rate).
            ctx.surface_outfall_link_limit[p.node as usize] += p.total / dt / self.degree_[p.node as usize] as f64 / opts.flow_1d_to_2d;
        }
    }

    /// apply(mesh, state, opts, ctx, dt, accumulated, node_conc, cell_accumulated): the routed
    /// outfall water (and its species) onto the surface cells, per point; true when any point
    /// moved water. Withdrawals draw on the prepared reservations, discharges spread by area.
    pub fn apply(
        &mut self,
        mesh: &MeshData,
        state: &mut SurfaceStateData,
        opts: &SolverOptions2D,
        ctx: &mut SimulationContext,
        dt: f64,
        accumulated: &mut [f64],
        node_conc: Option<&[f64]>,
        mut cell_accumulated: Option<&mut [f64]>,
    ) -> Result<bool, Terminated> {
        let r = self.apply_body(mesh, state, opts, ctx, dt, accumulated, node_conc, &mut cell_accumulated);
        self.recv_.clear(); // ResetReceivers
        r
    }

    fn apply_body(
        &mut self,
        mesh: &MeshData,
        state: &mut SurfaceStateData,
        opts: &SolverOptions2D,
        ctx: &mut SimulationContext,
        dt: f64,
        accumulated: &mut [f64],
        node_conc: Option<&[f64]>,
        cell_accumulated: &mut Option<&mut [f64]>,
    ) -> Result<bool, Terminated> {
        let mut changed = false;
        for pi in 0..self.points_.len() {
            let p = &self.points_[pi];
            let un = p.node as usize;
            let node_volume = (ctx.nodes.inflow[un] - ctx.nodes.outflow[un]) * opts.flow_1d_to_2d * dt;
            let node_weight = if node_volume < 0.0 { self.node_reserved_[un] } else { self.node_area_[un] };
            let point_weight = if node_volume < 0.0 { p.total } else { p.area };
            if node_volume < 0.0 && -node_volume > node_weight + 1e-10 * cmax(1.0, node_weight) {
                return Err(thrown("1D outfall exceeded its reserved surface water budget."));
            }
            let volume = if node_weight > 0.0 { node_volume * (point_weight / node_weight) } else { 0.0 };
            if volume == 0.0 {
                continue;
            }
            if -volume > p.total + 1e-10 * cmax(1.0, p.total) {
                return Err(thrown("1D outfall exceeded its reserved surface water budget."));
            }
            // SPREAD: a discharge lands on the own cells plus every cell in the radius that is
            // wet or no higher than the highest own cell (area-weighted); a withdrawal draws on
            // the own cells' reservations.
            let recv = &mut self.recv_;
            recv.clear();
            if volume > 0.0 && !p.spread.is_empty() {
                for j in 0..p.spread.len() {
                    let c = p.spread[j];
                    if j < p.cells.len() || state.depth[c as usize] >= opts.dry_depth || mesh.tri_cz[c as usize] <= p.bed_ref {
                        recv.push(c);
                    }
                }
            }
            let cells: &[i32] = if recv.is_empty() { &p.cells } else { &recv[..] };
            let mut total_weight = if volume < 0.0 { p.total } else { 0.0 };
            if volume > 0.0 {
                for &c in cells {
                    total_weight += mesh.tri_area[c as usize];
                }
            }
            if !(total_weight > 0.0) {
                return Err(thrown("Nonzero outfall exchange without a usable surface budget."));
            }
            let mut remaining = volume.abs();
            let mut applied = 0.0;
            let mut last = 0usize;
            for j in 0..cells.len() {
                if (if volume < 0.0 { p.reserved[j] } else { mesh.tri_area[cells[j] as usize] }) > 0.0 {
                    last = j;
                }
            }
            for j in 0..cells.len() {
                let c = cells[j];
                let uc = c as usize;
                let weight = if volume < 0.0 { p.reserved[j] } else { mesh.tri_area[uc] };
                if !(weight > 0.0) {
                    continue;
                }
                let mut take = if j == last { remaining } else { cmin(remaining, volume.abs() * (weight / total_weight)) };
                let old = state.volume[uc];
                if volume < 0.0 {
                    if take > old + 1e-10 * cmax(1.0, old) {
                        return Err(thrown("Overlapping outfalls overspent a donor cell."));
                    }
                    take = cmin(take, cmax(0.0, old));
                }
                remaining -= take;
                let dv = if volume < 0.0 { -take } else { take };
                // Review D2: mesh water entering the network through this outfall carries the
                // cell's quality into the outfall's next quality step.
                let to_1d = if opts.flow_1d_to_2d > 0.0 { 1.0 / opts.flow_1d_to_2d } else { 1.0 };
                let np_1d = ctx.n_pollutants();
                let ns_1d = ctx.reactions.n_species();
                let n_nodes = ctx.n_nodes();
                let tr = &mut state.transport;
                let q1d = &mut ctx.surface_outfall_quality;
                if volume < 0.0 && take > 0.0 && tr.n_species > 0 {
                    q1d.ensure(n_nodes, np_1d, ns_1d);
                    q1d.vol[un] += take * to_1d;
                    q1d.any = true;
                }
                for sp in 0..tr.n_species {
                    let ix = tr.idx(sp, c);
                    if volume < 0.0 {
                        let mass = tr.cell_mass[ix];
                        let dm = if old > 0.0 { mass * (take / old) } else { 0.0 };
                        tr.cell_mass[ix] -= dm;
                        tr.lost_coupling[sp as usize] += dm;
                        let m1 = dm * to_1d;
                        if sp < tr.n_pollut && sp < np_1d {
                            q1d.mass[un * np_1d as usize + sp as usize] += m1;
                        } else if sp >= tr.n_pollut && sp < tr.n_pollut + tr.n_msx && sp - tr.n_pollut < ns_1d {
                            q1d.msx[un * ns_1d as usize + (sp - tr.n_pollut) as usize] += m1;
                        } else if sp == tr.age_row {
                            q1d.age_vol[un] += m1;
                        } else if sp == tr.temp_row {
                            q1d.temp_vol[un] += m1;
                        }
                    } else if let Some(nc) = node_conc {
                        let ixc = un * tr.n_species as usize + sp as usize;
                        let dm = if ixc < nc.len() { take * nc[ixc] } else { 0.0 };
                        tr.cell_mass[ix] += dm;
                        tr.gained_coupling[sp as usize] += dm;
                    }
                }
                state.volume[uc] = old + dv;
                state.coupling_applied[uc] += dv;
                self.pending_[uc] += dv;
                if let Some(ca) = cell_accumulated.as_deref_mut() {
                    ca[uc] += dv;
                }
                if !state.transport.cell_runoff_vol.is_empty() {
                    state.transport.cell_runoff_vol[uc] -= dv;
                }
                let (eta, depth) = cell_eta_depth(mesh, opts, c, state.volume[uc]);
                state.head[uc] = eta;
                state.depth[uc] = depth;
                applied += dv;
            }
            accumulated[p.index as usize] += applied;
            // Book each routing-step direction before a surface batch can net a withdrawal
            // against a later discharge.
            let mb = &mut ctx.mass_balance_2d;
            if applied > 0.0 {
                mb.outfall_in += applied;
            } else {
                mb.outfall_out -= applied;
            }
            mb.final_storage += applied;
            changed = true;
        }
        self.prepared_ = false;
        Ok(changed)
    }

    /// recordTailwater(ctx, dt): book one routing step of each coupled outfall's tailwater:
    /// time engaged, time above the outlet crown, the highest engaged stage.
    pub fn record_tailwater(&mut self, ctx: &SimulationContext, dt: f64) {
        if !(dt > 0.0) {
            return;
        }
        for &n in &self.nodes_ {
            let row = ctx.node_subtypes.outfall_row(n);
            if row < 0 {
                continue;
            }
            let ur = row as usize;
            let un = n as usize;
            let h = ctx.node_subtypes.outfalls.head_2d[ur];
            let z = ctx.nodes.invert_elev[un];
            // The predicate Outfall::setAllOutfallDepths applies the override on.
            if !(ctx.node_subtypes.outfalls.ramp_2d[ur] > 0.0 && h > z) {
                continue;
            }
            self.tw_engaged_s_[un] += dt;
            if h >= z + ctx.nodes.full_depth[un] {
                self.tw_submerged_s_[un] += dt;
            }
            self.tw_max_[un] = cmax(self.tw_max_[un], h);
        }
    }

    /// tailwater(node, engaged_s, submerged_s, max_stage): a coupled outfall's tailwater
    /// record, false (outputs unchanged) for an uncoupled node.
    pub fn tailwater(&self, node: i32, engaged_s: &mut f64, submerged_s: &mut f64, max_stage: &mut f64) -> bool {
        if !self.nodes_.contains(&node) {
            return false;
        }
        let un = node as usize;
        *engaged_s = self.tw_engaged_s_[un];
        *submerged_s = self.tw_submerged_s_[un];
        *max_stage = self.tw_max_[un];
        true
    }
}
