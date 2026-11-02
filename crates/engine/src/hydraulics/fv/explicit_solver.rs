// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/fv/ExplicitFvSolver.hpp,
// ExplicitFvSolver.cpp, AcceptedWaterVolumes.hpp and FvKernels.hpp (FaceState, FaceFlux,
// tpaDepthOfArea) (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The CPU explicit finite-volume backend, initialization half: the state initialize()
// builds, initialize() itself, refreshDepths, cacheClosure (without the radius output),
// refreshNodeAreas, nodeVolumeFromDepth, algebraicActive, and the setters the Router
// calls at init (setConduitStepSink, trackAcceptedWater, requireConservativeNodeBalance).
//
// NOT YET TRANSLATED (the FV step, phase 1e): advance and everything it reaches
// (reinitialize, finalize, run_stats, nodeDepthFromVolume, rebuildActiveLists, the face
// flux / reconstruction / census / LTS / RK2 / TPA / node-solve machinery), and the
// members only that machinery reads or writes (save_*/rk_* snapshots, lo_flux_ and the
// other FCT scratch, cell_qlat_, node_lat_div_, node_qstruct_, live_stamp_, the due_*
// tier lists, div_link_/div_value_, last_nsteps_/last_h_/suggested_h_, the dt_argmin_*
// counters, press_step_ and the other step-local flags).
//
// Differences in form, never in result:
//   * The solver keeps no pointers to the Router's mesh and state (`mesh_`, `state_`):
//     initialize and the refresh helpers take them as arguments, and `attached_` records
//     that initialize ran (upstream: mesh_ != nullptr).
//   * The conduit-step sink (a C function pointer and its user pointer) is recorded as
//     `step_sink_set_`; the sink itself belongs to the step.
//   * The OpenMP-parallel refreshDepths loop is the serial loop: each iteration touches
//     only its own cell, as upstream's comment states, so the result is the same.
//   * The OPENSWMM_FV_* development switches are read from the environment as upstream
//     reads them (first byte tests, atof).

use crate::core::cmath::cmax;
use crate::core::constants;
use crate::core::options::fv::{FvOptions, RiemannSolver};
use crate::hydraulics::fv::kernels::{closure_eval, depth_of_area, K_DRY_AREA, K_DRY_DEPTH};
use crate::hydraulics::fv::{FvGeometry, NetworkMeshData, NetworkStateData, K_NODE_JUNCTION, K_NODE_OUTFALL, K_NODE_VIRTUAL, K_NODE_VOL_SAMPLES};

pub const K_MAX_LTS_TIERS: usize = 8;

/// Reconstructed state on one side of a face (FvKernels.hpp).
#[derive(Debug, Clone, Copy, Default)]
pub struct FaceState {
    pub a: f64,
    pub q: f64,
    pub u: f64,
    pub c: f64,
    pub i1: f64,
    pub press: u8,
}

/// Result of a face flux evaluation (FvKernels.hpp).
#[derive(Debug, Clone, Copy, Default)]
pub struct FaceFlux {
    pub mass: f64,
    pub mom: f64,
    pub sstar: f64,
    pub sl: f64,
    pub sr: f64,
}

/// Optional accepted-window water observations (AcceptedWaterVolumes.hpp).
#[derive(Debug, Clone, Default)]
pub struct AcceptedWaterVolumes {
    pub face_forward: Vec<f64>,
    pub face_reverse: Vec<f64>,
    pub structure_forward: Vec<f64>,
    pub structure_reverse: Vec<f64>,
    pub cell_lateral: Vec<f64>,
    pub cell_loss: Vec<f64>,
    pub node_deficit: Vec<f64>,
    pub face_lateral_left: Vec<f64>,
    pub face_lateral_right: Vec<f64>,
    pub face_spill: Vec<f64>,
}

fn assign<T: Clone>(v: &mut Vec<T>, n: usize, d: T) {
    v.clear();
    v.resize(n, d);
}

impl AcceptedWaterVolumes {
    pub fn resize(&mut self, nf: usize, ns: usize, nc: usize, nn: usize) {
        assign(&mut self.face_forward, nf, 0.0);
        assign(&mut self.face_reverse, nf, 0.0);
        assign(&mut self.structure_forward, ns, 0.0);
        assign(&mut self.structure_reverse, ns, 0.0);
        assign(&mut self.cell_lateral, nc, 0.0);
        assign(&mut self.cell_loss, nc, 0.0);
        assign(&mut self.node_deficit, nn, 0.0);
        assign(&mut self.face_lateral_left, nf, 0.0);
        assign(&mut self.face_lateral_right, nf, 0.0);
        assign(&mut self.face_spill, nf, 0.0);
    }
}

fn tpa_depth_of_area(g: &FvGeometry, a: f64) -> f64 {
    g.y_full + (a - g.a_crown) / g.t_slot
}

/// The CPU explicit FV backend (init half; see the module notes).
#[derive(Debug, Clone)]
pub struct ExplicitFvSolver {
    /// initialize() ran (upstream: mesh_ != nullptr)
    pub attached_: bool,
    pub opts_: FvOptions,

    pub f_mass_: Vec<f64>,
    pub f_mom_: Vec<f64>,
    pub f_sstar_: Vec<f64>,
    pub f_corr_l_: Vec<f64>,
    pub f_corr_r_: Vec<f64>,
    pub f_phi_l_: Vec<f64>,
    pub f_phi_r_: Vec<f64>,
    pub f_phi_flux_: Vec<f64>,
    pub f_state_l_: Vec<FaceState>,
    pub f_state_r_: Vec<FaceState>,
    pub f_flux_: Vec<FaceFlux>,
    pub hllc_: bool,
    pub cell_slope_: Vec<f64>,
    pub cell_eta_slope_: Vec<f64>,
    pub cell_u_slope_: Vec<f64>,
    pub cell_ho2_: Vec<i8>,
    pub cell_eta_: Vec<f64>,
    pub cell_u_: Vec<f64>,
    pub cell_ah_: Vec<f64>,
    pub cell_t_: Vec<f64>,
    pub cell_i1_: Vec<f64>,
    pub cell_q_int_: Vec<f64>,
    pub uf_grad_: Vec<f64>,

    pub tpa_: bool,
    pub mixed_wave_: bool,
    pub mw_rh_: bool,
    pub deg1_guard_slot_: bool,
    pub census_pass_far_: bool,
    pub lts_fit_: bool,
    pub tpa_scratch_: Vec<u8>,

    pub water_tracking_: bool,
    pub conservative_node_balance_: bool,
    pub withdrawal_rates_: Vec<f64>,
    pub unmet_withdrawal_: Vec<f64>,
    pub save_unmet_: Vec<f64>,
    pub rk_unmet_: Vec<f64>,
    pub water_: AcceptedWaterVolumes,
    pub save_water_: AcceptedWaterVolumes,
    pub rk_water_: AcceptedWaterVolumes,
    pub node_exch_: Vec<f64>,
    pub node_in_: Vec<f64>,
    pub node_out_: Vec<f64>,
    pub flood_vol_: Vec<f64>,
    pub node_qdummy_: Vec<f64>,
    pub dummy_vol_: Vec<f64>,
    pub inlet_control_: Vec<u8>,

    pub cell_tier_: Vec<u8>,
    pub node_tier_: Vec<u8>,
    pub face_tier_: Vec<u8>,
    pub relax_beta_: Vec<f64>,
    pub acc_a_: Vec<f64>,
    pub acc_q_: Vec<f64>,
    pub acc_nvol_: Vec<f64>,
    pub cells_by_tier_: Vec<Vec<i32>>,
    pub faces_by_tier_: Vec<Vec<i32>>,
    pub nodes_by_tier_: Vec<Vec<i32>>,
    pub lts_tiers_: i32,
    pub tier_occupancy_: [i64; K_MAX_LTS_TIERS],

    pub cell_active_: Vec<i8>,
    pub active_faces_: Vec<i32>,
    pub lists_valid_: bool,
    pub since_rebuild_: i32,
    pub census_count_: i32,
    pub dt_cache_: f64,
    pub dt_census_: f64,
    pub lts_valid_: bool,
    pub lts_countdown_: i32,
    pub lts_dt0_: f64,

    pub total_steps_: i64,
    pub total_flux_: i64,
    pub min_h_: f64,
    pub sim_time_: f64,
    pub active_sum_: f64,
    pub active_min_: f64,
    pub active_max_: f64,
    pub active_n_: i64,
    pub n_macro_cycles_: i64,
    pub n_macro_rejected_: i64,
    pub state_clean_: bool,
    pub acc_dirty_: bool,

    pub node_alg_: Vec<u8>,
    pub node_carry_: Vec<f64>,
    pub node_vfull_: Vec<f64>,
    pub node_pass_static_: Vec<u8>,
    pub node_pass_: Vec<u8>,
    pub node_pub_static_: Vec<u8>,
    pub node_pub_: Vec<u8>,
    pub face_a_rim_: Vec<f64>,
    pub face_lat_spill_: Vec<f64>,

    /// setConduitStepSink was called (the sink itself is the step's)
    pub step_sink_set_: bool,
}

impl Default for ExplicitFvSolver {
    fn default() -> Self {
        ExplicitFvSolver {
            attached_: false,
            opts_: FvOptions::default(),
            f_mass_: Vec::new(),
            f_mom_: Vec::new(),
            f_sstar_: Vec::new(),
            f_corr_l_: Vec::new(),
            f_corr_r_: Vec::new(),
            f_phi_l_: Vec::new(),
            f_phi_r_: Vec::new(),
            f_phi_flux_: Vec::new(),
            f_state_l_: Vec::new(),
            f_state_r_: Vec::new(),
            f_flux_: Vec::new(),
            hllc_: true,
            cell_slope_: Vec::new(),
            cell_eta_slope_: Vec::new(),
            cell_u_slope_: Vec::new(),
            cell_ho2_: Vec::new(),
            cell_eta_: Vec::new(),
            cell_u_: Vec::new(),
            cell_ah_: Vec::new(),
            cell_t_: Vec::new(),
            cell_i1_: Vec::new(),
            cell_q_int_: Vec::new(),
            uf_grad_: Vec::new(),
            tpa_: false,
            mixed_wave_: true,
            mw_rh_: true,
            deg1_guard_slot_: true,
            census_pass_far_: true,
            lts_fit_: false,
            tpa_scratch_: Vec::new(),
            water_tracking_: false,
            conservative_node_balance_: false,
            withdrawal_rates_: Vec::new(),
            unmet_withdrawal_: Vec::new(),
            save_unmet_: Vec::new(),
            rk_unmet_: Vec::new(),
            water_: AcceptedWaterVolumes::default(),
            save_water_: AcceptedWaterVolumes::default(),
            rk_water_: AcceptedWaterVolumes::default(),
            node_exch_: Vec::new(),
            node_in_: Vec::new(),
            node_out_: Vec::new(),
            flood_vol_: Vec::new(),
            node_qdummy_: Vec::new(),
            dummy_vol_: Vec::new(),
            inlet_control_: Vec::new(),
            cell_tier_: Vec::new(),
            node_tier_: Vec::new(),
            face_tier_: Vec::new(),
            relax_beta_: Vec::new(),
            acc_a_: Vec::new(),
            acc_q_: Vec::new(),
            acc_nvol_: Vec::new(),
            cells_by_tier_: Vec::new(),
            faces_by_tier_: Vec::new(),
            nodes_by_tier_: Vec::new(),
            lts_tiers_: 1,
            tier_occupancy_: [0; K_MAX_LTS_TIERS],
            cell_active_: Vec::new(),
            active_faces_: Vec::new(),
            lists_valid_: false,
            since_rebuild_: 0,
            census_count_: 0,
            dt_cache_: 0.0,
            dt_census_: 0.0,
            lts_valid_: false,
            lts_countdown_: 0,
            lts_dt0_: 0.0,
            total_steps_: 0,
            total_flux_: 0,
            min_h_: 0.0,
            sim_time_: 0.0,
            active_sum_: 0.0,
            active_min_: -1.0,
            active_max_: -1.0,
            active_n_: 0,
            n_macro_cycles_: 0,
            n_macro_rejected_: 0,
            state_clean_: false,
            acc_dirty_: false,
            node_alg_: Vec::new(),
            node_carry_: Vec::new(),
            node_vfull_: Vec::new(),
            node_pass_static_: Vec::new(),
            node_pass_: Vec::new(),
            node_pub_static_: Vec::new(),
            node_pub_: Vec::new(),
            face_a_rim_: Vec::new(),
            face_lat_spill_: Vec::new(),
            step_sink_set_: false,
        }
    }
}

fn env_first_byte(name: &str) -> Option<u8> {
    std::env::var_os(name).map(|v| v.as_encoded_bytes().first().copied().unwrap_or(0))
}

impl ExplicitFvSolver {
    fn tpa_cell(&self, state: &NetworkStateData, uc: usize) -> bool {
        self.tpa_ && state.cell_tpa[uc] != 0
    }

    /// cacheClosure(uc, g, h, nullptr)
    fn cache_closure(&mut self, state: &NetworkStateData, uc: usize, g: &FvGeometry, h: f64) {
        if h > K_DRY_DEPTH && !self.tpa_cell(state, uc) {
            let e = closure_eval(&g.closure_tbl, h);
            self.cell_ah_[uc] = e.a;
            self.cell_t_[uc] = e.t;
            self.cell_i1_[uc] = e.i1;
        } else {
            self.cell_ah_[uc] = 0.0;
            self.cell_t_[uc] = 0.0;
            self.cell_i1_[uc] = 0.0;
        }
    }

    pub fn refresh_depths(&mut self, mesh: &NetworkMeshData, state: &mut NetworkStateData) {
        let nc = mesh.n_cells();
        for c in 0..nc {
            let uc = c as usize;
            let g = &mesh.geom[mesh.cell_geom[uc] as usize];
            let mut a = state.cell_a[uc];
            if a < 0.0 {
                a = 0.0;
            }
            let tpa_cell = self.tpa_cell(state, uc);
            let h = if tpa_cell { tpa_depth_of_area(g, a) } else { depth_of_area(g, a) };
            state.cell_a[uc] = a;
            state.cell_h[uc] = h;
            self.cell_eta_[uc] = mesh.cell_zb[uc] + h;
            self.cache_closure(state, uc, g, h);
            if !tpa_cell && h <= K_DRY_DEPTH {
                state.cell_q[uc] = 0.0;
                self.cell_u_[uc] = 0.0;
            } else if tpa_cell {
                self.cell_u_[uc] = state.cell_q[uc] / cmax(a, K_DRY_AREA);
            } else {
                self.cell_u_[uc] = state.cell_q[uc] / a;
            }
        }
    }

    pub fn node_volume_from_depth(&self, mesh: &NetworkMeshData, state: &NetworkStateData, node: i32, depth: f64) -> f64 {
        let un = node as usize;
        if depth <= 0.0 {
            return 0.0;
        }
        let full = mesh.node_full_depth[un];
        if mesh.node_can_pond[un] != 0 && full > 0.0 && depth > full {
            return self.node_volume_from_depth(mesh, state, node, full) + (depth - full) * mesh.node_ponded_area[un];
        }
        let off = mesh.node_vol_off[un];
        if off < 0 {
            return state.node_surf_area[un] * depth;
        }

        let dmax = mesh.node_vol_dmax[un];
        let n = K_NODE_VOL_SAMPLES;
        let vmax = mesh.node_vol_tbl[(off + n - 1) as usize];
        if depth >= dmax {
            return vmax + mesh.node_vol_atop[un] * (depth - dmax);
        }

        let dd = dmax / (n - 1) as f64;
        let mut i = crate::core::cmath::cvt_i32(depth / dd);
        if i > n - 2 {
            i = n - 2;
        }
        if i < 0 {
            i = 0;
        }
        let t = (depth - i as f64 * dd) / dd;
        let v0 = mesh.node_vol_tbl[(off + i) as usize];
        let v1 = mesh.node_vol_tbl[(off + i + 1) as usize];
        v0 + (v1 - v0) * t
    }

    pub fn refresh_node_areas(&self, mesh: &NetworkMeshData, state: &mut NetworkStateData) {
        let nn = mesh.n_nodes();
        for n in 0..nn {
            let un = n as usize;
            if mesh.node_kind[un] == K_NODE_VIRTUAL {
                state.node_surf_area[un] = 0.0;
                continue;
            }
            if mesh.node_vol_off[un] >= 0 {
                let d = cmax(0.0, state.node_head[un] - mesh.node_invert[un]);
                let dd = mesh.node_vol_dmax[un] / (K_NODE_VOL_SAMPLES - 1) as f64;
                let v0 = self.node_volume_from_depth(mesh, state, n, d);
                let v1 = self.node_volume_from_depth(mesh, state, n, d + dd);
                state.node_surf_area[un] = cmax((v1 - v0) / dd, constants::MIN_SURFAREA);
                continue;
            }
            state.node_surf_area[un] = mesh.node_area[un];
        }
    }

    pub fn algebraic_active(&self, mesh: &NetworkMeshData, state: &NetworkStateData, n: i32) -> bool {
        let un = n as usize;
        if self.node_alg_[un] == 0 {
            return false;
        }
        if mesh.node_can_pond[un] != 0 && state.node_volume[un] > self.node_vfull_[un] {
            return false;
        }
        true
    }

    pub fn set_conduit_step_sink(&mut self) {
        self.step_sink_set_ = true;
    }

    /// trackAcceptedWater(enabled): (re)size the accepted-water observations.
    pub fn track_accepted_water(&mut self, mesh: &NetworkMeshData, enabled: bool) {
        self.water_tracking_ = enabled;
        if enabled && self.attached_ {
            self.water_.resize(mesh.n_faces() as usize, mesh.struct_link.len(), mesh.n_cells() as usize, mesh.n_nodes() as usize);
        } else {
            self.water_ = AcceptedWaterVolumes::default();
        }
        self.save_water_ = self.water_.clone();
        self.rk_water_ = self.water_.clone();
    }

    pub fn require_conservative_node_balance(&mut self, enabled: bool) {
        self.conservative_node_balance_ = enabled;
    }

    /// initialize(mesh, state, opts)
    pub fn initialize(&mut self, mesh: &mut NetworkMeshData, state: &mut NetworkStateData, opts: &FvOptions) {
        self.attached_ = true;
        self.opts_ = opts.clone();

        let nf = mesh.n_faces() as usize;
        let nc = mesh.n_cells() as usize;
        let nn = mesh.n_nodes() as usize;
        let ns = cmax(0, state.n_species) as usize;

        if self.water_tracking_ {
            self.water_.resize(nf, mesh.struct_link.len(), nc, nn);
            self.save_water_ = self.water_.clone();
            self.rk_water_ = self.water_.clone();
        }

        if mesh.node_sur_depth.len() < nn {
            mesh.node_sur_depth.resize(nn, 0.0);
        }
        if mesh.node_can_pond.len() < nn {
            mesh.node_can_pond.resize(nn, 0);
        }
        mesh.derive_face_geom();
        if mesh.face_gate.len() < nf {
            mesh.face_gate.resize(nf, 0);
        }
        if mesh.face_culvert.len() < nf {
            mesh.face_culvert.resize(nf, -1);
        }

        assign(&mut self.f_mass_, nf, 0.0);
        assign(&mut self.f_mom_, nf, 0.0);
        assign(&mut self.f_sstar_, nf, 0.0);
        assign(&mut self.f_corr_l_, nf, 0.0);
        assign(&mut self.f_corr_r_, nf, 0.0);
        assign(&mut self.f_phi_l_, ns * nf, 0.0);
        assign(&mut self.f_phi_r_, ns * nf, 0.0);
        assign(&mut self.f_phi_flux_, ns * nf, 0.0);
        assign(&mut self.f_state_l_, nf, FaceState::default());
        assign(&mut self.f_state_r_, nf, FaceState::default());
        assign(&mut self.f_flux_, nf, FaceFlux::default());
        self.hllc_ = opts.riemann == RiemannSolver::Hllc;
        assign(&mut self.cell_slope_, nc, 0.0);
        assign(&mut self.cell_eta_slope_, nc, 0.0);
        assign(&mut self.cell_u_slope_, nc, 0.0);
        assign(&mut self.cell_ho2_, nc, 0);

        assign(&mut self.cell_eta_, nc, 0.0);
        assign(&mut self.cell_u_, nc, 0.0);
        assign(&mut self.cell_ah_, nc, 0.0);
        assign(&mut self.cell_t_, nc, 0.0);
        assign(&mut self.cell_i1_, nc, 0.0);
        assign(&mut self.cell_q_int_, nc, 0.0);
        assign(&mut self.uf_grad_, if opts.unsteady_friction != 0 { nc } else { 0 }, 0.0);

        self.tpa_ = opts.pressure_closure == 1;
        {
            let mw = env_first_byte("OPENSWMM_FV_MIXED_WAVE");
            self.mixed_wave_ = mw != Some(b'0');
            self.mw_rh_ = self.mixed_wave_;
            let d1 = env_first_byte("OPENSWMM_FV_DEG1");
            self.deg1_guard_slot_ = d1 != Some(b'1');
            let cp = env_first_byte("OPENSWMM_FV_CENSUS_PASS");
            self.census_pass_far_ = cp != Some(b'0');
            let lf = env_first_byte("OPENSWMM_FV_LTS_FIT");
            self.lts_fit_ = lf == Some(b'1');
            if let Some(cf) = std::env::var_os("OPENSWMM_FV_CFL") {
                let v = crate::core::cstdlib::atof(cf.as_encoded_bytes());
                if v > 0.0 && v <= 1.0 {
                    self.opts_.cfl = v;
                }
            }
        }
        if self.tpa_ {
            assign(&mut state.cell_tpa, nc, 0);
            assign(&mut self.tpa_scratch_, nc, 0);
        } else {
            state.cell_tpa.clear();
            self.tpa_scratch_.clear();
        }

        assign(&mut self.withdrawal_rates_, nn, 0.0);
        assign(&mut self.unmet_withdrawal_, nn, 0.0);
        self.save_unmet_ = self.unmet_withdrawal_.clone();
        self.rk_unmet_ = self.unmet_withdrawal_.clone();
        assign(&mut self.node_exch_, nn, 0.0);
        assign(&mut self.node_in_, nn, 0.0);
        assign(&mut self.node_out_, nn, 0.0);
        assign(&mut self.flood_vol_, nn, 0.0);
        assign(&mut self.node_qdummy_, if mesh.struct_is_dummy.is_empty() { 0 } else { nn }, 0.0);
        assign(&mut self.dummy_vol_, mesh.struct_link.len(), 0.0);
        assign(&mut self.inlet_control_, mesh.n_conduits() as usize, 0);

        assign(&mut self.cell_tier_, nc, 0);
        assign(&mut self.node_tier_, nn, 0);
        assign(&mut self.face_tier_, nf, 0);
        assign(&mut self.relax_beta_, nf, 0.0);
        self.acc_a_.clear();
        self.acc_q_.clear();
        self.acc_nvol_.clear();
        self.cells_by_tier_.clear();
        self.faces_by_tier_.clear();
        self.nodes_by_tier_.clear();
        self.lts_tiers_ = 1;
        self.tier_occupancy_ = [0; K_MAX_LTS_TIERS];

        assign(&mut self.cell_active_, nc, 1);
        self.active_faces_.clear();
        self.lists_valid_ = false;
        self.since_rebuild_ = 0;
        self.census_count_ = 0;
        self.dt_cache_ = 0.0;
        self.dt_census_ = 0.0;
        self.lts_valid_ = false;
        self.lts_countdown_ = 0;
        self.lts_dt0_ = 0.0;

        self.total_steps_ = 0;
        self.total_flux_ = 0;
        self.min_h_ = 0.0;
        self.sim_time_ = 0.0;
        self.active_sum_ = 0.0;
        self.active_min_ = -1.0;
        self.active_max_ = -1.0;
        self.active_n_ = 0;
        self.n_macro_cycles_ = 0;
        self.n_macro_rejected_ = 0;
        self.state_clean_ = false;
        self.acc_dirty_ = false;

        self.refresh_depths(mesh, state);
        self.refresh_node_areas(mesh, state);

        // Algebraic-junction eligibility.
        assign(&mut self.node_alg_, nn, 0);
        assign(&mut self.node_carry_, nn, 0.0);
        assign(&mut self.node_vfull_, nn, 0.0);
        assign(&mut self.node_pass_static_, nn, 0);
        assign(&mut self.node_pass_, nn, 0);
        assign(&mut self.node_pub_static_, nn, 0);
        assign(&mut self.node_pub_, nn, 0);
        for un in 0..nn {
            self.node_vfull_[un] = self.node_volume_from_depth(mesh, state, un as i32, mesh.node_full_depth[un]);
            self.node_alg_[un] = if mesh.node_kind[un] == K_NODE_JUNCTION && mesh.node_vol_off[un] < 0 && mesh.node_face_ptr[un + 1] > mesh.node_face_ptr[un] { 1 } else { 0 };
            let mut clean = self.node_alg_[un] != 0;
            if clean {
                for p in mesh.node_face_ptr[un]..mesh.node_face_ptr[un + 1] {
                    let uf = mesh.node_face_idx[p as usize] as usize;
                    if mesh.face_culvert[uf] >= 0 || mesh.face_gate[uf] != 0 {
                        clean = false;
                    }
                }
            }
            self.node_pub_static_[un] = if clean { 1 } else { 0 };
            if clean && mesh.node_face_ptr[un + 1] - mesh.node_face_ptr[un] == 2 {
                self.node_pass_static_[un] = 1;
            }
        }

        // Rim cap for the diverted lateral.
        assign(&mut self.face_a_rim_, nf, f64::INFINITY);
        assign(&mut self.face_lat_spill_, nf, 0.0);
        for un in 0..nn {
            if self.node_pass_static_[un] == 0 {
                continue;
            }
            let full = mesh.node_full_depth[un];
            if !(full > 0.0) {
                continue;
            }
            let rim = mesh.node_invert[un] + full + if mesh.node_can_pond[un] != 0 { 0.0 } else { mesh.node_sur_depth[un] };
            for p in mesh.node_face_ptr[un]..mesh.node_face_ptr[un + 1] {
                let uf = mesh.node_face_idx[p as usize] as usize;
                let c = if mesh.face_cl[uf] >= 0 { mesh.face_cl[uf] } else { mesh.face_cr[uf] };
                if c < 0 {
                    continue;
                }
                let uc = c as usize;
                let g = &mesh.geom[mesh.cell_geom[uc] as usize];
                if g.is_open != 0 {
                    continue;
                }
                let h_rim = rim - mesh.cell_zb[uc];
                self.face_a_rim_[uf] = if h_rim > 0.0 { crate::hydraulics::fv::kernels::area_of_depth(g, h_rim) } else { 0.0 };
            }
        }
        // Seed the conserved node water.
        for n in 0..nn {
            if mesh.node_kind[n] == K_NODE_VIRTUAL || mesh.node_kind[n] == K_NODE_OUTFALL || self.algebraic_active(mesh, state, n as i32) {
                state.node_volume[n] = 0.0;
            } else {
                let v = self.node_volume_from_depth(mesh, state, n as i32, cmax(0.0, state.node_head[n] - mesh.node_invert[n]));
                state.node_volume[n] = v;
            }
        }
    }
}
