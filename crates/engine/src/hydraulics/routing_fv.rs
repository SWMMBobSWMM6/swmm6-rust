// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Routing.cpp (stepFv,
// accumulateStructureFlows, refreshFvBoundaryFlows, publishFv, acceptedQualityWater,
// fvStepSink, the FV branches of getAdaptiveStep and accumulateLinkStepStats)
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5-6 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The Router's FV routing step: the forcing it assembles (laterals net of losses, the
// withdrawal rates, outfall stages, structure flows and their caps, distributed conduit
// losses), the per-substep structure refresh (FV_STRUCTURE_COUPLING SUBSTEP), the checked
// advance, the withdrawal refund, the conduit-loss acceptance, and the publish of the cell
// and node state onto the link and node reporting arrays.
//
// Differences in form, never in result:
//   * The forcing's arrays (fv_lateral_, fv_fixed_head_, fv_struct_flow_, fv_link_q_cap_,
//     fv_cond_loss_) move into the FvStepForcing for the advance and back after it; upstream
//     hands the solver pointers to them. The refresh hook's evaluate_structures writes the
//     same arrays, here through the forcing it is handed.
//   * The non-conduit callback gets the DW solver's node arrays as under DYNWAVE; under FV
//     that solver is never initialized, so its arrays are empty (upstream: data() of an
//     empty vector, a null pointer).
//   * OPENSWMM_FV_RECEIPT_TRACE (a CSV trace written to a file named by the environment) is
//     left out: it prints only.

use crate::core::cmath::{cmax, cmin};
use crate::core::constants;
use crate::core::context::SimulationContext;
use crate::core::engine_step_routing::StepError;
use crate::core::options::fv::{StructureCoupling, TimeIntegration};
use crate::data::forcing::ForcingMode;
use crate::core::unit_conversion as ucf;
use crate::data::link::{FlowClass, LinkData, LinkType, XsectShape};
use crate::data::node::NodeType;
use crate::hydraulics::dynamic_wave_step::{NcView, NonConduitFn};
use crate::hydraulics::fv::explicit_step::{advance_network_checked, FvStepForcing};
use crate::hydraulics::fv::kernels as fvk;
use crate::hydraulics::fv::K_NODE_VIRTUAL;
use crate::hydraulics::hyd_structures::StructureSolver;
use crate::hydraulics::node;
use crate::hydraulics::outfall;
use crate::hydraulics::routing::Router;
use crate::quality::quality_routing::{BoundaryWaterImport, FvTransportWater, WaterTransfer};

/// evaluate_structures: the non-conduit callback at the current heads, then the structure
/// flows and the DUMMY caps the solver reads.
#[allow(clippy::too_many_arguments)]
fn evaluate_structures(ctx: &mut SimulationContext, dt: f64, structures: &mut StructureSolver, non_conduit_fn: &mut Option<&mut NonConduitFn>, view: &mut NcView, struct_flow: &mut [f64], link_q_cap: &mut [f64]) -> Result<(), StepError> {
    if let Some(f) = non_conduit_fn.as_mut() {
        f(ctx, dt, 0, view, structures)?;
    }
    let nl = ctx.n_links();
    for j in 0..nl {
        let uj = j as usize;
        struct_flow[uj] = if ctx.links.type_[uj] == LinkType::Conduit { 0.0 } else { ctx.links.flow[uj] };
        if ctx.links.xsect_shape[uj] == XsectShape::Dummy {
            let lim = ctx.links.q_limit[uj];
            link_q_cap[uj] = if ctx.links.setting[uj] == 0.0 { 0.0 } else if lim > 0.0 { lim } else { -1.0 };
        }
    }
    Ok(())
}

/// accumulateStructureFlows(ctx, t_now): bank the structure discharge in force since the
/// last refresh.
fn accumulate_structure_flows(ctx: &SimulationContext, t_now: f64, struct_flow: &[f64], struct_int: &mut [f64], t_prev: &mut f64) {
    let span = t_now - *t_prev;
    if span <= 0.0 {
        return;
    }
    for j in 0..ctx.n_links() {
        let uj = j as usize;
        if ctx.links.type_[uj] == LinkType::Conduit {
            continue;
        }
        struct_int[uj] += struct_flow[uj] * span;
    }
    *t_prev = t_now;
}

impl Router {
    /// The FV paths the translation does not cover yet (checked before the step changes
    /// anything): the solver options that select them.
    pub fn fv_not_translated(&self, ctx: &SimulationContext) -> Option<&'static str> {
        let s = self.fv_solver_.as_ref()?;
        let o = &s.opts_;
        if o.pressurized_implicit {
            return Some("FV_PRESSURIZED_IMPLICIT (PressurizedHeadSolver)");
        }
        if s.tpa_ {
            return Some("FV TPA pressure closure");
        }
        if o.time_integration == TimeIntegration::Rk2 {
            return Some("FV RK2 time integration");
        }
        if o.unsteady_friction != 0 {
            return Some("FV unsteady friction");
        }
        if self.fv_state_.n_species > 0 {
            return Some("FV species transport");
        }
        let _ = ctx;
        None
    }

    /// stepFv(ctx, dt, non_conduit_fn): one FV routing step (assemble forcing, advance,
    /// publish). Returns the substep count.
    pub fn step_fv(&mut self, ctx: &mut SimulationContext, dt: f64, structures: &mut StructureSolver, mut non_conduit_fn: Option<&mut NonConduitFn>) -> Result<i32, StepError> {
        if self.fv_solver_.is_none() {
            return Ok(0);
        }
        let nn = ctx.n_nodes();
        let nl = ctx.n_links();
        let mut has_withdrawals = false;
        let user_rate = |ctx: &SimulationContext, n: usize| -> f64 {
            let mut q = ctx.nodes.user_lat_flow[n];
            if n < ctx.forcing.node_lat_inflow_mode.len() {
                if ctx.forcing.node_lat_inflow_mode[n] == ForcingMode::Override {
                    q = ctx.forcing.node_lat_inflow_value[n];
                } else if ctx.forcing.node_lat_inflow_mode[n] == ForcingMode::Add {
                    q += ctx.forcing.node_lat_inflow_value[n];
                }
            }
            q
        };

        for n in 0..nn {
            let un = n as usize;
            self.fv_lateral_[un] = ctx.nodes.lat_flow[un] - ctx.nodes.losses[un];
            {
                let mut withdrawal = 0.0;
                for q in [ctx.nodes.ext_inflow[un], ctx.nodes.dwf_inflow[un], ctx.nodes.gw_inflow[un], ctx.nodes.rdii_inflow[un], ctx.nodes.iface_inflow[un], user_rate(ctx, un)] {
                    withdrawal += cmax(0.0, -q);
                }
                let solver = self.fv_solver_.as_mut().unwrap();
                solver.withdrawal_rates_[un] = withdrawal;
                has_withdrawals |= withdrawal > 0.0;
                if withdrawal > 0.0 && solver.water_.face_forward.is_empty() {
                    solver.track_accepted_water(&self.fv_mesh_, true);
                }
            }
            self.fv_fixed_head_[un] = if ctx.nodes.type_[un] == NodeType::Outfall { ctx.nodes.invert_elev[un] + ctx.nodes.depth[un] } else { f64::NAN };
        }

        // The forcing takes the Router's arrays for the advance.
        let mut forcing = FvStepForcing {
            node_lateral: std::mem::take(&mut self.fv_lateral_),
            node_fixed_head: std::mem::take(&mut self.fv_fixed_head_),
            node_source_max: if ctx.surface_outfall_link_limit.is_empty() { None } else { Some(ctx.surface_outfall_link_limit.clone()) },
            structure_flow: std::mem::take(&mut self.fv_struct_flow_),
            link_q_cap: std::mem::take(&mut self.fv_link_q_cap_),
            conduit_loss: std::mem::take(&mut self.fv_cond_loss_),
            has_refresh: false,
            n_nodes: nn,
            n_links: nl,
            node_withdrawal: has_withdrawals,
        };

        let result = self.step_fv_advance(ctx, dt, structures, &mut non_conduit_fn, &mut forcing);

        self.fv_lateral_ = forcing.node_lateral;
        self.fv_fixed_head_ = forcing.node_fixed_head;
        self.fv_struct_flow_ = forcing.structure_flow;
        self.fv_link_q_cap_ = forcing.link_q_cap;
        self.fv_cond_loss_ = forcing.conduit_loss;
        result?;

        // Withdrawals the solver could not honour are refunded.
        if has_withdrawals && dt > 0.0 {
            let solver = self.fv_solver_.as_mut().unwrap();
            let mut refused = solver.unmet_withdrawal_.clone();
            if !solver.water_.face_forward.is_empty() {
                let water = &solver.water_;
                let mut source_port = vec![0.0f64; ctx.n_nodes() as usize];
                for f in 0..self.fv_mesh_.n_faces() as usize {
                    let mut n = self.fv_mesh_.face_node[f];
                    if n < 0 {
                        n = self.fv_mesh_.face_vj_node[f];
                    }
                    if n >= 0 {
                        source_port[n as usize] += water.face_lateral_left[f] + water.face_lateral_right[f];
                    }
                }
                for n in 0..ctx.n_nodes() as usize {
                    if self.fv_mesh_.node_kind[n] == K_NODE_VIRTUAL && self.fv_lateral_[n] < 0.0 {
                        refused[n] += cmax(0.0, -self.fv_lateral_[n] * dt + source_port[n]);
                    }
                }
            }
            for n in 0..ctx.n_nodes() {
                let un = n as usize;
                let requested = solver.withdrawal_rates_[un];
                let refund = cmin(requested, cmax(0.0, refused[un] / dt));
                let fraction = if requested > 0.0 { refund / requested } else { 0.0 };
                let accepted_user = cmax(0.0, -user_rate(ctx, un)) * (1.0 - fraction);
                ctx.mass_balance.step_ext_withdrawal += accepted_user;
                solver.withdrawal_rates_[un] = cmax(0.0, requested - refund);
                if refund <= 0.0 {
                    continue;
                }
                fn reduce(q: &mut f64, category: &mut f64, sign: f64, fraction: f64) {
                    if *q >= 0.0 {
                        return;
                    }
                    let accepted_refund = -*q * fraction;
                    *q += accepted_refund;
                    *category += sign * accepted_refund;
                }
                reduce(&mut ctx.nodes.ext_inflow[un], &mut ctx.mass_balance.step_ext_withdrawal, -1.0, fraction);
                reduce(&mut ctx.nodes.dwf_inflow[un], &mut ctx.mass_balance.step_dw_inflow, 1.0, fraction);
                reduce(&mut ctx.nodes.gw_inflow[un], &mut ctx.mass_balance.step_gw_inflow, 1.0, fraction);
                reduce(&mut ctx.nodes.rdii_inflow[un], &mut ctx.mass_balance.step_rdii_inflow, 1.0, fraction);
                reduce(&mut ctx.nodes.iface_inflow[un], &mut ctx.mass_balance.step_ext_inflow, 1.0, fraction);
                ctx.nodes.lat_flow[un] += refund;
            }
        }
        // Distributed conduit losses the cells could not supply are scaled back.
        {
            let solver = self.fv_solver_.as_ref().unwrap();
            if dt > 0.0 && !solver.water_.face_forward.is_empty() {
                let cd = &mut ctx.link_subtypes.conduits;
                let accepted = &solver.water_;
                for row in 0..self.fv_mesh_.n_conduits() as usize {
                    let begin = self.fv_mesh_.conduit_cell_begin[row];
                    let count = self.fv_mesh_.conduit_cell_count[row];
                    if begin < 0 || count <= 0 {
                        continue;
                    }
                    let requested = (cd.seep_loss_rate[row] + cd.evap_loss_rate[row]) * cmax(1, cd.barrels[row]) as f64 * dt;
                    let mut actual = 0.0;
                    for c in begin..begin + count {
                        actual += accepted.cell_loss[c as usize];
                    }
                    if requested != 0.0 {
                        let ratio = actual / requested;
                        cd.seep_loss_rate[row] *= ratio;
                        cd.evap_loss_rate[row] *= ratio;
                    }
                }
            }
        }
        self.publish_fv(ctx, dt);
        Ok(self.fv_solver_.as_ref().unwrap().last_num_steps() as i32)
    }

    /// stepFv from evaluate_structures to the checked advance, with the forcing in hand.
    fn step_fv_advance(&mut self, ctx: &mut SimulationContext, dt: f64, structures: &mut StructureSolver, non_conduit_fn: &mut Option<&mut NonConduitFn>, forcing: &mut FvStepForcing) -> Result<(), StepError> {
        let Router { fv_solver_, fv_mesh_, fv_state_, fv_opts_, fv_struct_int_, fv_struct_t_prev_, dw_solver_, lstep_rs_, lstep_pmin_, lstep_pmax_, lstep_pdt_time_, lstep_ptime_, lstep_pbin_, .. } = self;
        let initialized = dw_solver_.is_initialized();
        let mut view = NcView { new_surf_area: &mut dw_solver_.xnode_.new_surf_area, sumdqdh: &mut dw_solver_.xnode_.sumdqdh, bypassed: &dw_solver_.bypassed_, initialized };

        fv_struct_int_.iter_mut().for_each(|x| *x = 0.0);
        *fv_struct_t_prev_ = 0.0;
        evaluate_structures(ctx, dt, structures, non_conduit_fn, &mut view, &mut forcing.structure_flow, &mut forcing.link_q_cap)?;

        // Distributed conduit losses as a rate per unit length (ft2/s).
        {
            let solver = fv_solver_.as_mut().unwrap();
            for r in 0..fv_mesh_.n_conduits() as usize {
                let begin = fv_mesh_.conduit_cell_begin[r];
                if begin < 0 {
                    forcing.conduit_loss[r] = 0.0;
                    continue;
                }
                let cd = &ctx.link_subtypes.conduits;
                let rate = (cd.evap_loss_rate[r] + cd.seep_loss_rate[r]) * cmax(1, cd.barrels[r]) as f64;
                let len = fv_mesh_.cell_dx[begin as usize] * fv_mesh_.conduit_cell_count[r] as f64;
                forcing.conduit_loss[r] = if len > 0.0 { rate / len } else { 0.0 };
                if rate != 0.0 && solver.water_.face_forward.is_empty() {
                    solver.track_accepted_water(fv_mesh_, true);
                }
            }
        }

        let substep = fv_opts_.structure_coupling == StructureCoupling::Substep;
        forcing.has_refresh = substep;
        let mesh = &*fv_mesh_;
        let mut refresh = |state: &crate::hydraulics::fv::NetworkStateData, fc: &mut FvStepForcing, t_elapsed: f64| -> Result<(), StepError> {
            // refreshFvBoundaryFlows(ctx, t_elapsed, eval_structures)
            accumulate_structure_flows(ctx, t_elapsed, &fc.structure_flow, fv_struct_int_, fv_struct_t_prev_);
            let nn = ctx.n_nodes();
            for n in 0..nn {
                let un = n as usize;
                if mesh.node_kind[un] == K_NODE_VIRTUAL {
                    continue;
                }
                let depth = state.node_head[un] - mesh.node_invert[un];
                ctx.nodes.depth[un] = cmax(0.0, depth);
                ctx.nodes.head[un] = state.node_head[un];
            }
            outfall::set_all_outfall_depths(ctx, t_elapsed);
            for n in 0..nn {
                let un = n as usize;
                if ctx.nodes.type_[un] != NodeType::Outfall {
                    continue;
                }
                fc.node_fixed_head[un] = ctx.nodes.invert_elev[un] + ctx.nodes.depth[un];
            }
            evaluate_structures(ctx, dt, structures, non_conduit_fn, &mut view, &mut fc.structure_flow, &mut fc.link_q_cap)
        };
        let rs = *lstep_rs_;
        let mut sink = |conduit: i32, dts: f64, dur: f64| {
            // fvStepSink(user, conduit, dt, dur)
            let j = mesh.conduit_link[conduit as usize];
            if j < 0 {
                return;
            }
            let dts = cmin(dts, rs);
            let u = j as usize;
            if lstep_ptime_[u] <= 0.0 || dts < lstep_pmin_[u] {
                lstep_pmin_[u] = dts;
            }
            if dts > lstep_pmax_[u] {
                lstep_pmax_[u] = dts;
            }
            lstep_pdt_time_[u] += dts * dur;
            lstep_ptime_[u] += dur;
            lstep_pbin_[u * LinkData::N_LSTEP_BINS as usize + LinkData::lstep_bin(dts, rs) as usize] += dur;
        };
        let solver = fv_solver_.as_mut().unwrap();
        let sink_on = solver.step_sink_set_;
        advance_network_checked(solver, mesh, fv_state_, 0.0, dt, forcing, if substep { Some(&mut refresh) } else { None }, if sink_on { Some(&mut sink) } else { None })?;
        Ok(())
    }

    /// publishFv(ctx, dt): the cell and node state onto the link and node reporting arrays.
    pub fn publish_fv(&mut self, ctx: &mut SimulationContext, dt: f64) {
        let us = ucf::get_unit_system(ctx.options.flow_units as i32);
        let impl_ = self.fv_solver_.as_ref();
        let mesh = &self.fv_mesh_;
        let st = &self.fv_state_;

        // ---- links
        for r in 0..mesh.n_conduits() {
            let ur = r as usize;
            let begin = mesh.conduit_cell_begin[ur];
            let count = mesh.conduit_cell_count[ur];
            if begin < 0 || count <= 0 {
                continue;
            }
            let j = mesh.conduit_link[ur];
            let uj = j as usize;
            let g = &mesh.geom[mesh.conduit_section[ur] as usize];
            let mut sum_len = 0.0;
            let mut q_len = 0.0;
            let mut a_len = 0.0;
            let mut vol = 0.0;
            let mut slot_vol = 0.0;
            for c in begin..begin + count {
                let uc = c as usize;
                let dx = mesh.cell_dx[uc];
                let q = match impl_ {
                    Some(s) if dt > 0.0 => s.cell_q_int_[uc] / dt,
                    _ => st.cell_q[uc],
                };
                sum_len += dx;
                q_len += q * dx;
                a_len += st.cell_a[uc] * dx;
                vol += st.cell_a[uc] * dx;
                slot_vol += cmax(0.0, st.cell_a[uc] - g.a_crown) * dx;
            }
            let q_mean = if sum_len > 0.0 { q_len / sum_len } else { 0.0 };
            let a_mean = if sum_len > 0.0 { a_len / sum_len } else { 0.0 };
            ctx.links.flow[uj] = q_mean;
            ctx.links.depth[uj] = cmin(fvk::depth_of_area(g, a_mean), g.y_full);
            if let Some(s) = impl_ {
                let cr = ctx.link_subtypes.conduit_row(j);
                if cr >= 0 {
                    ctx.link_subtypes.conduits.inlet_control[cr as usize] = s.inlet_control_[ur];
                }
            }
            ctx.links.volume[uj] = vol;
            ctx.links.slot_volume[uj] = slot_vol;
            let a_conv = if g.is_open == 0 && a_mean > g.a_crown { g.a_crown } else { a_mean };
            let v = if a_conv > 0.0 { q_mean / a_conv } else { 0.0 };
            let w = fvk::width_of_depth(g, ctx.links.depth[uj]);
            let hyd_depth = if w > 0.0 { a_conv / fvk::width_of_depth(g, ctx.links.depth[uj]) } else { 0.0 };
            ctx.links.froude[uj] = if hyd_depth > 0.0 { v.abs() * constants::INV_SQRT_GRAVITY / hyd_depth.sqrt() } else { 0.0 };

            let uc_first = begin as usize;
            let uc_last = (begin + count - 1) as usize;
            let up_dry = st.cell_h[uc_first] <= fvk::K_DRY_DEPTH;
            let dn_dry = st.cell_h[uc_last] <= fvk::K_DRY_DEPTH;
            let fc = if up_dry && dn_dry {
                FlowClass::Dry
            } else if up_dry {
                FlowClass::UpDry
            } else if dn_dry {
                FlowClass::DnDry
            } else if ctx.links.froude[uj] > 1.0 {
                FlowClass::Supercritical
            } else {
                FlowClass::Subcritical
            };
            ctx.links.flow_class[uj] = fc;
            let cr_fs = ctx.link_subtypes.conduit_row(j);
            if cr_fs >= 0 {
                let yf = g.y_full;
                let fs = (if st.cell_h[uc_first] >= yf { 1 } else { 0 }) | (if st.cell_h[uc_last] >= yf { 2 } else { 0 });
                ctx.link_subtypes.conduits.full_state[cr_fs as usize] = fs as i8;
            }
        }

        // ---- nodes
        for n in 0..ctx.n_nodes() {
            let un = n as usize;
            if mesh.node_kind[un] == K_NODE_VIRTUAL {
                continue;
            }
            let mut head = st.node_head[un];
            if let Some(s) = impl_ {
                if un < s.node_pub_.len() && s.node_pub_[un] != 0 {
                    let b0 = mesh.node_face_ptr[un];
                    let e0 = mesh.node_face_ptr[un + 1];
                    let mut eta = -1.0e30f64;
                    let mut wet = false;
                    for p in b0..e0 {
                        let uf = mesh.node_face_idx[p as usize] as usize;
                        let c = if mesh.face_cl[uf] >= 0 { mesh.face_cl[uf] } else { mesh.face_cr[uf] };
                        if c < 0 {
                            continue;
                        }
                        let uc = c as usize;
                        if st.cell_h[uc] <= fvk::K_DRY_DEPTH {
                            continue;
                        }
                        if mesh.face_zb[uf] > head {
                            continue;
                        }
                        eta = cmax(eta, cmin(mesh.cell_zb[uc], mesh.face_zb[uf]) + st.cell_h[uc]);
                        wet = true;
                    }
                    if wet {
                        head = eta;
                    }
                }
            }
            let depth = head - mesh.node_invert[un];
            ctx.nodes.depth[un] = cmax(0.0, depth);
            ctx.nodes.head[un] = head;
            let full_d = mesh.node_full_depth[un];
            if full_d > 0.0 && ctx.nodes.depth[un] > full_d {
                ctx.nodes.volume[un] = if mesh.node_can_pond[un] != 0 { cmax(st.node_volume[un], ctx.nodes.full_volume[un]) } else { ctx.nodes.full_volume[un] };
            } else {
                ctx.nodes.volume[un] = node::get_volume(&ctx.nodes, n, ctx.nodes.depth[un], Some(&ctx.tables), us, Some(&ctx.node_subtypes));
            }
            if impl_.is_some() {
                ctx.nodes.volume[un] = if ctx.nodes.type_[un] == NodeType::Outfall { 0.0 } else { st.node_volume[un] };
            }
            if let Some(s) = impl_ {
                if dt > 0.0 {
                    ctx.nodes.overflow[un] = s.flood_vol_[un] / dt;
                }
            }
            if let Some(s) = impl_ {
                if dt > 0.0 {
                    let lat = ctx.nodes.lat_flow[un];
                    ctx.nodes.inflow[un] = s.node_in_[un] / dt + if lat > 0.0 { lat } else { 0.0 };
                    ctx.nodes.outflow[un] = s.node_out_[un] / dt + (if lat < 0.0 { -lat } else { 0.0 }) + ctx.nodes.losses[un];
                }
            }
        }

        // Structure flows: the time-weighted mean over the routing step.
        if dt > 0.0 {
            let fv_struct_flow = &self.fv_struct_flow_;
            accumulate_structure_flows(ctx, dt, fv_struct_flow, &mut self.fv_struct_int_, &mut self.fv_struct_t_prev_);
            for j in 0..ctx.n_links() {
                let uj = j as usize;
                if ctx.links.type_[uj] == LinkType::Conduit {
                    continue;
                }
                ctx.links.flow[uj] = self.fv_struct_int_[uj] / dt;
            }
        }
        let impl_ = self.fv_solver_.as_ref();
        let mesh = &self.fv_mesh_;
        let st = &self.fv_state_;
        if let Some(s) = impl_ {
            if dt > 0.0 {
                let dvol = &s.dummy_vol_;
                for si in 0..mesh.struct_link.len() {
                    if mesh.struct_is_dummy.is_empty() || mesh.struct_is_dummy[si] == 0 {
                        continue;
                    }
                    let uj = mesh.struct_link[si] as usize;
                    ctx.links.flow[uj] = dvol[si] / dt;
                    ctx.links.depth[uj] = 0.0;
                    ctx.links.volume[uj] = 0.0;
                }
            }
        }
        for j in 0..ctx.n_links() {
            let uj = j as usize;
            if ctx.links.type_[uj] == LinkType::Conduit && ctx.links.xsect_shape[uj] != XsectShape::Dummy {
                continue;
            }
            let n1 = ctx.links.node1[uj];
            let n2 = ctx.links.node2[uj];
            if n1 < 0 || n2 < 0 {
                continue;
            }
            let q = ctx.links.flow[uj];
            let u1 = n1 as usize;
            let u2 = n2 as usize;
            if q > 0.0 {
                ctx.nodes.outflow[u1] += q;
                ctx.nodes.inflow[u2] += q;
            } else {
                ctx.nodes.inflow[u1] -= q;
                ctx.nodes.outflow[u2] -= q;
            }
        }
        // Virtual junctions: the shared face state.
        for f in 0..mesh.n_faces() {
            let uf = f as usize;
            if mesh.face_virtual[uf] == 0 {
                continue;
            }
            let n = mesh.face_vj_node[uf];
            if n < 0 {
                continue;
            }
            let un = n as usize;
            let zf = mesh.face_zb[uf];
            let mut eta = zf;
            let mut any_tpa = false;
            let mut eta_tpa = 0.0;
            for c in [mesh.face_cl[uf], mesh.face_cr[uf]] {
                if c < 0 {
                    continue;
                }
                let uc = c as usize;
                if !st.cell_tpa.is_empty() && st.cell_tpa[uc] != 0 {
                    let e_t = mesh.cell_zb[uc] + st.cell_h[uc];
                    eta_tpa = if any_tpa { cmax(eta_tpa, e_t) } else { e_t };
                    any_tpa = true;
                    continue;
                }
                eta = cmax(eta, cmin(mesh.cell_zb[uc], zf) + st.cell_h[uc]);
            }
            if any_tpa {
                eta = if eta > zf { cmax(eta, eta_tpa) } else { eta_tpa };
            }
            ctx.nodes.head[un] = eta;
            ctx.nodes.depth[un] = cmax(0.0, eta - mesh.node_invert[un]);
            ctx.nodes.volume[un] = 0.0;
            ctx.nodes.overflow[un] = 0.0;
            let (q_in, q_out) = match impl_ {
                Some(s) if dt > 0.0 => (s.node_in_[un] / dt, s.node_out_[un] / dt),
                _ => (0.0, 0.0),
            };
            let lat = ctx.nodes.lat_flow[un];
            ctx.nodes.inflow[un] = q_in + if lat > 0.0 { lat } else { 0.0 };
            ctx.nodes.outflow[un] = q_out + if lat < 0.0 { -lat } else { 0.0 };
            if lat > 0.0 {
                ctx.nodes.outflow[un] += lat;
            }
        }
    }

    /// acceptedQualityWater(ctx, dt, out): the accepted water of the step, as the quality
    /// transport's mixing input.
    pub fn accepted_quality_water(&self, ctx: &SimulationContext, dt: f64, out: &mut FvTransportWater) {
        out.ready = false;
        let native = match self.fv_solver_.as_ref() {
            Some(s) => s,
            None => return,
        };
        if dt <= 0.0 || native.water_.face_forward.is_empty() {
            return;
        }
        let mesh = &self.fv_mesh_;
        let nn = ctx.n_nodes() as usize;
        let nl = ctx.n_links() as usize;
        let n = nn + nl;
        let r = &native.water_;
        {
            let w = &mut out.water;
            for v in [&mut w.old_volume, &mut w.new_volume, &mut w.source_volume, &mut w.exported, &mut w.seeped, &mut w.evaporated] {
                v.clear();
                v.resize(n, 0.0);
            }
            w.transfers.clear();
        }
        out.boundary_imports.clear();
        out.withdrawal.clear();
        out.withdrawal.resize(nn, 0.0);
        out.flooding.clear();
        out.flooding.resize(nn, 0.0);
        for i in 0..nn {
            if ctx.nodes.type_[i] != NodeType::Outfall {
                out.water.old_volume[i] = ctx.nodes.old_volume[i];
                out.water.new_volume[i] = ctx.nodes.volume[i];
            }
            let row = ctx.node_subtypes.storage_row(i as i32);
            if row >= 0 {
                let stg = &ctx.node_subtypes.storages;
                out.water.seeped[i] = stg.exfil_loss[row as usize];
                out.water.evaporated[i] = stg.evap_loss[row as usize];
            }
            if mesh.node_can_pond[i] == 0 {
                out.flooding[i] = native.flood_vol_[i];
            }
            out.withdrawal[i] = native.withdrawal_rates_[i] * dt;
            if ctx.nodes.type_[i] == NodeType::Outfall && out.withdrawal[i] > 0.0 {
                out.boundary_imports.push(BoundaryWaterImport { node: i, receiver: i, volume: out.withdrawal[i] });
                out.water.source_volume[i] += out.withdrawal[i];
            }
            out.water.exported[i] = out.withdrawal[i] + out.flooding[i];
        }
        for j in 0..nl {
            if ctx.links.type_[j] != LinkType::Conduit || ctx.links.xsect_shape[j] == XsectShape::Dummy {
                continue;
            }
            out.water.old_volume[nn + j] = ctx.links.old_volume[j];
            out.water.new_volume[nn + j] = ctx.links.volume[j];
            let cr = ctx.link_subtypes.conduit_row(j as i32);
            if cr >= 0 {
                let cd = &ctx.link_subtypes.conduits;
                let cr = cr as usize;
                out.water.seeped[nn + j] = cmax(0.0, cd.seep_loss_rate[cr]) * cmax(1, cd.barrels[cr]) as f64 * dt;
                out.water.evaporated[nn + j] = cmax(0.0, cd.evap_loss_rate[cr]) * cmax(1, cd.barrels[cr]) as f64 * dt;
            }
        }
        let transfer = |out: &mut FvTransportWater, donor: usize, receiver: usize, v: f64| {
            if v <= 0.0 || donor == receiver {
                return;
            }
            if donor < nn && ctx.nodes.type_[donor] == NodeType::Outfall {
                out.boundary_imports.push(BoundaryWaterImport { node: donor, receiver, volume: v });
                out.water.source_volume[receiver] += v;
            } else {
                out.water.transfers.push(WaterTransfer { donor, receiver, volume: v });
                if receiver < nn && ctx.nodes.type_[receiver] == NodeType::Outfall {
                    out.water.exported[receiver] += v;
                }
            }
        };
        let cellcv = |c: i32| -> usize { nn + mesh.conduit_link[mesh.cell_conduit[c as usize] as usize] as usize };
        for f in 0..mesh.n_faces() as usize {
            let cl = mesh.face_cl[f];
            let cr = mesh.face_cr[f];
            let mut nd = mesh.face_node[f];
            let vj = mesh.face_vj_node[f];
            if cl >= 0 && cr >= 0 && mesh.cell_conduit[cl as usize] == mesh.cell_conduit[cr as usize] {
                continue;
            }
            if vj >= 0 {
                let vju = vj as usize;
                transfer(out, cellcv(cl), vju, r.face_forward[f]);
                transfer(out, vju, cellcv(cr), r.face_forward[f]);
                transfer(out, cellcv(cr), vju, r.face_reverse[f]);
                transfer(out, vju, cellcv(cl), r.face_reverse[f]);
                nd = vj;
            } else if nd >= 0 {
                let left = if cl >= 0 { cellcv(cl) } else { nd as usize };
                let right = if cr >= 0 { cellcv(cr) } else { nd as usize };
                transfer(out, left, right, r.face_forward[f]);
                transfer(out, right, left, r.face_reverse[f]);
            }
            if nd >= 0 {
                for side in 0..2 {
                    let c = if side != 0 { cr } else { cl };
                    if c < 0 {
                        continue;
                    }
                    let v = if side != 0 { r.face_lateral_right[f] } else { r.face_lateral_left[f] };
                    if v >= 0.0 {
                        transfer(out, nd as usize, cellcv(c), v);
                    } else {
                        transfer(out, cellcv(c), nd as usize, -v);
                    }
                }
                if r.face_spill[f] > 0.0 {
                    transfer(out, cellcv(if cl >= 0 { cl } else { cr }), nd as usize, r.face_spill[f]);
                }
            }
        }
        for i in 0..mesh.struct_link.len() {
            transfer(out, mesh.struct_n1[i] as usize, mesh.struct_n2[i] as usize, r.structure_forward[i]);
            transfer(out, mesh.struct_n2[i] as usize, mesh.struct_n1[i] as usize, r.structure_reverse[i]);
        }
        out.ready = true;
    }
}
