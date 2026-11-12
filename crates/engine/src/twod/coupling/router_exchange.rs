// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/SurfaceRouter2D.cpp (updateTrackingPonds,
// publishEqualizeAreas, computeCouplingConductances, and the junction ledger of coAdvanceStep)
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb Buahin,
// licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The router's own run-time exchange members (wave 4b), for the router's step (wave 4a's
// SurfaceRouter2D::updateOutfallsPreRouting / coAdvanceStep) and the routing step to call:
//   * update_tracking_ponds: COUPLING_POND_MODE TRACKING, the pond area relaxed toward the wet
//     footprint, the pond volume kept (updateOutfallsPreRouting's last call);
//   * publish_equalize_areas: COUPLING_EQUALIZE's per-node 1D surface area and the held
//     window (coAdvanceStep, just before the solver's advance);
//   * compute_coupling_conductances: the exchange's head sensitivities for the dynamic-wave
//     node continuity (stepRouting's non-conduit callback);
//   * co_advance_junction_ledger and co_advance_junction_queue: the two halves of
//     coAdvanceStep's junction ledger (the marcher's per-point integral of the exchange into
//     nodes.coupling_volume and the quality queues, then into nodes.coupling_queue), which
//     upstream writes inline; the first runs right after applyAgingAndReactions, the second
//     after the aquifer's node exchange (accumulateMassBalance and the subsurface block run
//     between them). `exch` is the solver's last_coupling_exchange().
//
// Differences in form, never in result: the router's options and mesh are shared handles
// (twod::Shared), borrowed for the call; the solver's exchange vector is an argument (the
// caller passes the solver's own, cloned, since it lives inside the router).

#![allow(clippy::all)]

use crate::core::cmath::{cmax, cmin};
use crate::core::constants;
use crate::core::context::SimulationContext;
use crate::core::engine::SwmmEngine;
use crate::core::unit_conversion as ucf;
use crate::twod::coupling::node_coupling::{compute_node_coupling_dqdh1d, CouplingAttachments};
use crate::twod::coupling::pond_model as pond;
use crate::twod::surface_router_2d::SurfaceRouter2D;

impl SurfaceRouter2D {
    /// updateTrackingPonds(ctx, dt): each tracking pond's area A_eq relaxes toward the wet
    /// footprint around its node (floored at its chamber area), alpha = min(1, dt / tau) with
    /// tau = COUPLING_POND_RELAX or 5 dt; a ponded node's depth is re-derived so it holds the
    /// same pond volume.
    pub fn update_tracking_ponds(&mut self, ctx: &mut SimulationContext, dt: f64) {
        if self.pond_track_ptr_.is_empty() || !(dt > 0.0) {
            return;
        }
        let opts = self.options_.borrow();
        let tau = if opts.coupling_pond_relax > 0.0 { opts.coupling_pond_relax } else { 5.0 * dt };
        let alpha = cmin(1.0, dt / tau);
        let m2_to_ft2 = opts.len_2d_to_1d * opts.len_2d_to_1d;
        let n = self.pond_track_ptr_.len() - 1;
        for i in 0..n {
            if self.pond_track_ptr_[i] == self.pond_track_ptr_[i + 1] {
                continue;
            }
            let mut wet = 0.0;
            for k in self.pond_track_ptr_[i]..self.pond_track_ptr_[i + 1] {
                if self.state_.depth[self.pond_track_cell_[k as usize] as usize] >= opts.dry_depth {
                    wet += self.pond_track_share_[k as usize];
                }
            }
            let target = cmax(self.pond_chamber_ft2_[i], wet * m2_to_ft2);
            let a_old = ctx.coupled_pond_area[i];
            let a_new = a_old + alpha * (target - a_old);
            if a_new == a_old {
                continue;
            }
            let band = ctx.coupled_pond_band[i];
            let x = ctx.nodes.depth[i] - ctx.nodes.full_depth[i];
            if x > 0.0 {
                let v = pond::volume(x, a_old, band);
                ctx.nodes.depth[i] = ctx.nodes.full_depth[i] + pond::depth_for_volume(v, a_new, band);
                ctx.nodes.head[i] = ctx.nodes.invert_elev[i] + ctx.nodes.depth[i];
            }
            ctx.coupled_pond_area[i] = a_new;
        }
    }

    /// publishEqualizeAreas(ctx, window): each coupled node's 1D surface area (m2) for the
    /// advance about to run (its pond area at the current depth, floored at MIN_SURFAREA) and
    /// the held window; the state's node_pond_area_m2 stays detached when no point equalizes.
    pub fn publish_equalize_areas(&mut self, ctx: &SimulationContext, window: f64) {
        self.state_.node_pond_area_m2 = false;
        self.state_.coupling_window_s = window;
        if !self.node_coupling_points_.iter().any(|cp| cp.equalize) {
            return;
        }
        let us = ucf::get_unit_system(ctx.options.flow_units as i32);
        let ucf_len = ucf::UCF[ucf::LENGTH as usize][us as usize];
        let min_sa = if ctx.options.min_surf_area > 0.0 { ctx.options.min_surf_area / (ucf_len * ucf_len) } else { constants::MIN_SURFAREA };
        let opts = self.options_.borrow();
        let ft2_to_m2 = opts.len_1d_to_2d * opts.len_1d_to_2d;
        self.node_pond_area_m2_.clear();
        self.node_pond_area_m2_.resize(ctx.n_nodes() as usize, 0.0);
        for cp in &self.node_coupling_points_ {
            let i = cp.node_idx as usize;
            let x = ctx.nodes.depth[i] - ctx.nodes.full_depth[i];
            let pond_ft2 = if i < ctx.coupled_pond_area.len() {
                pond::area(x, ctx.coupled_pond_area[i], ctx.coupled_pond_band[i])
            } else if x > 0.0 {
                ctx.nodes.ponded_area[i]
            } else {
                0.0
            };
            self.node_pond_area_m2_[i] = cmax(min_sa, pond_ft2) * ft2_to_m2;
        }
        self.state_.node_pond_area_m2 = true;
    }

    /// computeCouplingConductances(ctx, out): appends (node, G) for each junction point with a
    /// positive head sensitivity, G in ft3/s per ft of the node's own head.
    pub fn compute_coupling_conductances(&self, ctx: &SimulationContext, out: &mut Vec<(i32, f64)>) {
        if !self.active_ {
            return;
        }
        let opts = self.options_.borrow();
        let mesh = self.mesh_.borrow();
        let att = CouplingAttachments::of(&self.state_, &self.coupling_head_cells_, &self.node_pond_area_m2_);
        // G in SI is m3/s per metre of (2D-frame) 1D head; the 1D consumes ft3/s per ft.
        let to_1d = opts.flow_2d_to_1d * opts.len_1d_to_2d;
        for cp in &self.node_coupling_points_ {
            let g = compute_node_coupling_dqdh1d(cp, &mesh, &self.state_, &att, &ctx.nodes, &opts);
            if g > 0.0 {
                out.push((cp.node_idx, g * to_1d));
            }
        }
    }

    /// coAdvanceStep's junction ledger, first half (right after applyAgingAndReactions): the
    /// batch's per-point exchange integral `exch` (m3, + = 2D -> 1D drain) into
    /// nodes.coupling_volume (1D ft3; cleared, then summed per node) and, for transported
    /// rows, the drain's species mass (net of this window's paired spill) into the nodes'
    /// quality queues. Under a common source stage the integral also leaves the coupling
    /// cell's cumulative receipt.
    pub fn co_advance_junction_ledger(&mut self, ctx: &mut SimulationContext, exch: &[f64]) {
        if self.common_stage_ {
            for k in 0..exch.len() {
                self.coupling_cum_applied_[self.node_coupling_points_[k].cell_idx as usize] -= exch[k];
            }
        }
        if exch.is_empty() {
            return;
        }
        let opts = self.options_.borrow();
        for cp in &self.node_coupling_points_ {
            ctx.nodes.coupling_volume[cp.node_idx as usize] = 0.0;
        }
        let n = cmin(exch.len(), self.node_coupling_points_.len());
        for k in 0..n {
            let ni = self.node_coupling_points_[k].node_idx as usize;
            ctx.nodes.coupling_volume[ni] += exch[k] * opts.flow_2d_to_1d;
        }
        // S3: the drain's species mass, booked per point by the marcher at the cell's
        // concentration, to the node's mass queues in 1D mass units.
        let tr = &self.state_.transport;
        let ns = tr.n_species as usize;
        if ns > 0 && !tr.exch_mass.is_empty() {
            for k in 0..n {
                let ni = self.node_coupling_points_[k].node_idx as usize;
                if (ni + 1) * ns > ctx.nodes.coupling_qual_queue.len() || (k + 1) * ns > tr.exch_mass.len() {
                    continue;
                }
                // The spill pairs against this window's own drain first.
                let s_ = if k < tr.exch_spill.len() { tr.exch_spill[k] } else { 0.0 };
                let d_ = exch[k] + s_; // gross drain volume
                let paired = cmin(s_, if d_ > 0.0 { d_ } else { 0.0 });
                let npol = tr.n_pollut as usize;
                let nmsx = tr.n_msx as usize;
                let nsp_all = ctx.reactions.n_species() as usize;
                if nmsx > 0 && ctx.reactions.msx_coupling_queue.len() != ctx.n_nodes() as usize * nsp_all {
                    let len = ctx.n_nodes() as usize * nsp_all;
                    ctx.reactions.msx_coupling_queue.clear();
                    ctx.reactions.msx_coupling_queue.resize(len, 0.0);
                }
                for sp in 0..ns {
                    let mut inc = tr.exch_mass[k * ns + sp];
                    if paired > 0.0 && (ni + 1) * ns <= self.node_row_conc_.len() {
                        inc -= paired * self.node_row_conc_[ni * ns + sp];
                    }
                    let row: Option<&mut f64> = if sp < npol {
                        if (ni + 1) * npol <= ctx.nodes.coupling_qual_queue.len() {
                            Some(&mut ctx.nodes.coupling_qual_queue[ni * npol + sp])
                        } else {
                            None
                        }
                    } else if sp as i32 == tr.age_row {
                        if ni < ctx.nodes.coupling_age_vol_queue.len() {
                            Some(&mut ctx.nodes.coupling_age_vol_queue[ni])
                        } else {
                            None
                        }
                    } else if sp as i32 == tr.temp_row {
                        if ni < ctx.nodes.coupling_temp_vol_queue.len() {
                            Some(&mut ctx.nodes.coupling_temp_vol_queue[ni])
                        } else {
                            None
                        }
                    } else if sp >= npol && sp < npol + nmsx {
                        let kk = ni * nsp_all + (sp - npol);
                        if kk < ctx.reactions.msx_coupling_queue.len() {
                            Some(&mut ctx.reactions.msx_coupling_queue[kk])
                        } else {
                            None
                        }
                    } else {
                        None
                    };
                    let Some(row) = row else {
                        continue;
                    };
                    *row += inc * opts.flow_2d_to_1d;
                    // Mass cannot be owed (not for the signed temperature row).
                    if *row < 0.0 && !tr.signed_row(sp as i32) {
                        *row = 0.0;
                    }
                }
            }
        }
    }

    /// coAdvanceStep's junction ledger, second half (after accumulateMassBalance and the
    /// aquifer's node exchange): each coupled node's coupling_volume moves to its delivery
    /// queue, which assembleLateralInflows drains at a uniform rate over the batch span `dt`.
    pub fn co_advance_junction_queue(&mut self, ctx: &mut SimulationContext, exch: &[f64], dt: f64) {
        if exch.is_empty() {
            return;
        }
        let mut queued = false;
        for cp in &self.node_coupling_points_ {
            let ni = cp.node_idx as usize;
            if ctx.nodes.coupling_volume[ni] != 0.0 {
                ctx.nodes.coupling_queue[ni] += ctx.nodes.coupling_volume[ni];
                ctx.nodes.coupling_volume[ni] = 0.0;
                queued = true;
            }
        }
        if queued {
            ctx.coupling_delivery_remaining = cmax(ctx.coupling_delivery_remaining, dt);
        }
    }
}

impl SwmmEngine {
    /// The 2D router and the context together, mutable (a gate's access to the run-time
    /// exchange members; upstream's harnesses reach them through `#define private public`).
    #[doc(hidden)]
    pub fn gate_surface_router_2d_and_context(&mut self) -> (&mut SurfaceRouter2D, &mut SimulationContext) {
        (&mut self.surface_router_, &mut self.ctx_)
    }

    /// The 1D router and the context together, mutable (the exchange gate reaches the
    /// dynamic-wave solver's initNodeStates through it).
    #[doc(hidden)]
    pub fn gate_router_and_context(&mut self) -> (&mut crate::hydraulics::routing::Router, &mut SimulationContext) {
        (&mut self.router_, &mut self.ctx_)
    }
}
