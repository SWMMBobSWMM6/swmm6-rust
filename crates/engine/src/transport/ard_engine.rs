// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine
// src/engine/transport/components/EulerianArdComponent/ArdEngine.cpp and ArdEngine.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The EULERIAN_ARD quality engine (QUALITY_SOLVER EULERIAN_ARD): species advected on the FV
// transport mesh under any hydraulic routing model. Each routing step projects the link and
// node hydraulics onto the mesh faces, subcycles advection (MUSCL / QUICKEST-ULTIMATE + FCT,
// transport/fv_species_kernels.rs) with node CSTR stores and the optional dispersion solve,
// then ages, reacts and publishes into links.conc / nodes.conc (and the age, temperature and
// reactions-species element arrays).
//
// Differences in form, never in result:
//   * step() returns Result: the reaction stage and the booster sources can meet code not
//     translated yet (StepError::NotTranslated).
//   * The DETAILED_OUTPUT sidecar is a C stdio stream (core::cprintf::CFile); the doubles are
//     written with "%g" at precision 6, which is what std::ostream's default formatting does.
//   * thermalRegions, speciesRegions, applyThermalIncrements, applySpeciesIncrements,
//     checkpointState and restoreCheckpointState serve the C API's state accessors and engine
//     checkpoints, which are not on the 1D run path; they are not translated.
//
// The reaction stage and the booster sources can still meet code not translated yet (the
// ROS2/BDF2 solvers, EQUIL species, booster rows), which they report as
// StepError::NotTranslated.

use crate::core::cmath::{cclamp, cmax, cmin, cvt_i32};
use crate::core::context::SimulationContext;
use crate::core::cprintf::CFile;
use crate::core::cstdio::printf_g;
use crate::core::engine_step_routing::StepError;
use crate::core::options::{Bytes, RoutingModel};
use crate::core::unit_conversion::{get_unit_system, LENGTH, UCF};
use crate::data::ard_config::ArdDispersionMode;
use crate::data::heat::HeatSource;
use crate::data::node::NodeType;
use crate::data::table::table_tseries_lookup_cursor;
use crate::data::water_age::WaterAgeSource;
use crate::hydraulics::fv::explicit_solver::{FaceFlux, FaceState};
use crate::hydraulics::fv::kernels::K_DRY_AREA;
use crate::hydraulics::fv::mesh_builder::build_network_mesh;
use crate::hydraulics::fv::{NetworkMeshData, NetworkStateData, K_NODE_VIRTUAL};
use crate::hydraulics::node::store_volume;
use crate::quality::quality_routing::{book_negative_cell_source_clamp, book_negative_source_clamp, update_untreated_storage_hrt};
use crate::transport::fv_species_kernels::{dispersion_solve, reconstruct_scalars, DispersionView, SpeciesKernelView};
use crate::data::heat_override::HeatElement;
use crate::hydraulics::fv::kernels::{depth_of_area, hyd_rad_of_depth, width_of_depth};
use crate::hydraulics::node::get_surf_area;
use crate::transport::heat_bed::{bed_coupling_from_contact, bed_exchange_enabled, bed_exchange_q, exchange_pair, ground_temp_for, ground_temperature, relax_pair};
use crate::transport::heat_fluxes::{net_flux_out, relax_t, update_solar_forcing, K_PROBE_C};
use crate::transport::heat_legacy::apply_initial_temp_overrides;
use crate::transport::heat_overrides::sediment_for;
use crate::transport::initial_quality_seeds::{initial_age_seconds_for, initial_temp_for};
use crate::transport::reaction_ard_binding::{ard_reactions_active, react_ard_stage, ArdArrays};
use crate::transport::reaction_legacy_binding::ensure_msx_state;
use crate::transport::reaction_legacy_step::reaction_booster_mass;
use crate::transport::transport_policy::network_1d_enables;

const K_CFL: f64 = 0.9;
const K_MAX_SUBSTEPS: i32 = 512;
const K_MIN_STORE_VOL: f64 = 1.0e-9;
const K_RESYNC_NOISE_REL: f64 = 1.0e-4;
const K_TEMP_FLOOR_C: f64 = -50.0;
const K_TEMP_CEIL_C: f64 = 100.0;
const K_ZERO_VOLUME_FT3: f64 = 0.0353147;
const K_GRAVITY: f64 = 32.174;
const K_MIN_SLOPE: f64 = 1.0e-5;
const K_MIN_FISCHER_DEPTH: f64 = 0.01;
use crate::data::ard_config::K_LITERS_PER_FT3;

/// transport::ArdEngine
#[derive(Default)]
pub struct ArdEngine {
    pub mesh_: NetworkMeshData,
    pub state_: NetworkStateData,
    pub face_q_: Vec<f64>,
    pub node_mass_: Vec<f64>,
    pub node_vol_: Vec<f64>,
    pub cell_src_vol_: Vec<f64>,
    pub cell_src_mass_: Vec<f64>,
    pub has_vj_src_: bool,
    pub f_mass_: Vec<f64>,
    pub f_sstar_: Vec<f64>,
    pub f_phi_l_: Vec<f64>,
    pub f_phi_r_: Vec<f64>,
    pub f_phi_flux_: Vec<f64>,
    pub cell_slope_: Vec<f64>,
    pub lo_flux_: Vec<f64>,
    pub anti_flux_: Vec<f64>,
    pub td_: Vec<f64>,
    pub anew_: Vec<f64>,
    pub rplus_: Vec<f64>,
    pub rminus_: Vec<f64>,
    pub cell_u_: Vec<f64>,
    pub f_state_l_: Vec<FaceState>,
    pub f_state_r_: Vec<FaceState>,
    pub f_flux_: Vec<FaceFlux>,
    pub cell_active_: Vec<u8>,
    pub active_faces_: Vec<i32>,
    pub bc_node_: Vec<i32>,
    pub bc_srow_: Vec<i32>,
    pub bc_ts_: Vec<i32>,
    pub bc_value_: Vec<f64>,
    pub bc_now_: Vec<f64>,
    pub src_crow_: Vec<i32>,
    pub src_srow_: Vec<i32>,
    pub src_ts_: Vec<i32>,
    pub src_value_: Vec<f64>,
    pub src_len_: Vec<f64>,
    pub src_now_: Vec<f64>,
    pub disp_active_: bool,
    pub disp_mode_: i32,
    pub disp_global_ft_: f64,
    pub conduit_disp_ft_: Vec<f64>,
    pub cell_disp_: Vec<f64>,
    pub detail_out_: Option<CFile>,
    pub detail_active_: bool,
    pub detail_time_s_: f64,
    pub age_row_: i32,
    pub nm_: i32,
    pub temp_row_: i32,
    pub warnings_: Vec<Bytes>,
    pub initialized_: bool,
    pub warned_cfl_clamp_: bool,
}

impl std::fmt::Debug for ArdEngine {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ArdEngine").field("initialized_", &self.initialized_).finish()
    }
}

fn assign<T: Clone>(v: &mut Vec<T>, n: usize, x: T) {
    v.clear();
    v.resize(n, x);
}

impl ArdEngine {
    /// A fresh engine (the member initializers: age_row_ and temp_row_ -1).
    pub fn new() -> ArdEngine {
        ArdEngine { age_row_: -1, temp_row_: -1, ..Default::default() }
    }

    pub fn initialized(&self) -> bool {
        self.initialized_
    }
    pub fn warnings(&self) -> &[Bytes] {
        &self.warnings_
    }

    /// ArdEngine::init: build the transport mesh and size and seed the state.
    pub fn init(&mut self, ctx: &mut SimulationContext) -> bool {
        self.warnings_.clear();
        self.initialized_ = false;

        let en = network_1d_enables(ctx);
        let np = en.n_pollut;
        let mut nm = 0;
        if ard_reactions_active(ctx) {
            if en.msx_has_wall {
                self.warnings_.push(b"QUALITY_SOLVER EULERIAN_ARD: WALL species have no transport semantics under this engine yet \xe2\x80\x94 falling back to LEGACY, which runs them per element (R4).".to_vec());
                return false;
            }
            nm = en.n_msx;
        }
        let na = if en.age { 1 } else { 0 };
        let nt = if en.temperature { 1 } else { 0 };
        let ns = np + nm + na + nt;
        self.nm_ = nm;
        self.age_row_ = if na > 0 { np + nm } else { -1 };
        self.temp_row_ = if nt > 0 { np + nm + na } else { -1 };
        if ns <= 0 {
            return false;
        }

        let mut mesh_opts = ctx.options.fv.clone();
        {
            let ucf_len0 = UCF[LENGTH as usize][get_unit_system(ctx.options.flow_units as i32) as usize];
            mesh_opts.cell_length /= ucf_len0;
            mesh_opts.slot_celerity /= ucf_len0;
            mesh_opts.dispersion /= ucf_len0 * ucf_len0;
            if ctx.ard_config.configured && ctx.ard_config.target_dx > 0.0 {
                if ctx.options.routing_model == RoutingModel::Fv {
                    self.warnings_.push(b"[TRANSPORT_OPTIONS] TARGET_DX is ignored under FLOW_ROUTING FV \xe2\x80\x94 the hydraulic solver's mesh governs transport (FV_CELL_LENGTH).".to_vec());
                } else {
                    mesh_opts.cell_length = ctx.ard_config.target_dx / ucf_len0;
                }
            }
        }
        let report = build_network_mesh(ctx, &mesh_opts, &mut self.mesh_);
        for w in report.warnings.iter() {
            self.warnings_.push(w.clone());
        }
        if !report.errors.is_empty() {
            for e in report.errors.iter() {
                let mut m = b"EULERIAN_ARD mesh error: ".to_vec();
                m.extend_from_slice(e);
                self.warnings_.push(m);
            }
            return false;
        }
        let nc = self.mesh_.n_cells();
        let nf = self.mesh_.n_faces();
        let nn = if self.mesh_.node_face_ptr.is_empty() { 0 } else { self.mesh_.node_face_ptr.len() as i32 - 1 };
        if nc <= 0 {
            return false;
        }

        self.state_.resize(nc, nn, ns);

        let unc = nc as usize;
        let unf = nf as usize;
        let uns = ns as usize;

        assign(&mut self.face_q_, unf, 0.0);
        assign(&mut self.f_mass_, unf, 0.0);
        assign(&mut self.f_sstar_, unf, 0.0);
        assign(&mut self.f_state_l_, unf, FaceState::default());
        assign(&mut self.f_state_r_, unf, FaceState::default());
        assign(&mut self.f_flux_, unf, FaceFlux::default());
        assign(&mut self.f_phi_l_, uns * unf, 0.0);
        assign(&mut self.f_phi_r_, uns * unf, 0.0);
        assign(&mut self.f_phi_flux_, uns * unf, 0.0);
        assign(&mut self.cell_slope_, unc, 0.0);
        assign(&mut self.lo_flux_, unf, 0.0);
        assign(&mut self.anti_flux_, unf, 0.0);
        assign(&mut self.td_, unc, 0.0);
        assign(&mut self.anew_, unc, 0.0);
        assign(&mut self.rplus_, unc, 1.0);
        assign(&mut self.rminus_, unc, 1.0);
        assign(&mut self.cell_u_, unc, 0.0);
        assign(&mut self.cell_active_, unc, 1u8);
        self.active_faces_.resize(unf, 0);
        for f in 0..nf {
            self.active_faces_[f as usize] = f;
        }

        assign(&mut self.node_mass_, nn as usize * uns, 0.0);
        assign(&mut self.node_vol_, nn as usize, 0.0);

        let unp = np as usize;
        let unm_msx = nm as usize;
        if nm > 0 {
            ensure_msx_state(ctx);
        }

        // A2a: hotstart-loaded ages win over INITIAL_STATE
        let mut hs_node_age: Vec<f64> = Vec::new();
        let mut hs_link_age: Vec<f64> = Vec::new();
        let mut age_from_hs = false;
        if self.age_row_ >= 0 && ctx.water_age_state.hotstart_loaded {
            age_from_hs = true;
            hs_node_age = ctx.water_age_state.node_age.clone();
            hs_link_age = ctx.water_age_state.link_age.clone();
        }
        let age_global = ctx.water_age_config.global_age[WaterAgeSource::InitialState as usize];
        let age_seed_link = |ctx: &SimulationContext, link: i32| -> f64 {
            if age_from_hs && (link as usize) < hs_link_age.len() {
                return hs_link_age[link as usize];
            }
            initial_age_seconds_for(ctx, true, link, age_global)
        };
        let temp_seed = ctx.heat_config.global_temp[HeatSource::InitialState as usize];
        let age_seed_node = |ctx: &SimulationContext, nd: usize| -> f64 {
            if age_from_hs && nd < hs_node_age.len() {
                return hs_node_age[nd];
            }
            initial_age_seconds_for(ctx, false, nd as i32, age_global)
        };

        let n_links = ctx.n_links();
        for r in 0..self.mesh_.conduit_link.len() {
            let link = self.mesh_.conduit_link[r];
            if link < 0 || link >= n_links {
                continue;
            }
            let b = self.mesh_.conduit_cell_begin[r];
            let n = self.mesh_.conduit_cell_count[r];
            let mut length = 0.0;
            for i in 0..n {
                length += self.mesh_.cell_dx[(b + i) as usize];
            }
            let vol = ctx.links.volume[link as usize];
            let a0 = if length > 0.0 { vol / length } else { 0.0 };
            for i in 0..n {
                let uc = (b + i) as usize;
                self.state_.cell_a[uc] = a0;
                for s in 0..unp {
                    self.state_.cell_phi[s * unc + uc] = ctx.links.conc[link as usize * unp + s];
                }
                for m in 0..unm_msx {
                    let xi = link as usize * unm_msx + m;
                    self.state_.cell_phi[(unp + m) * unc + uc] = if xi < ctx.reactions.msx_link_conc.len() { ctx.reactions.msx_link_conc[xi] } else { ctx.reactions.init_global[m] };
                }
                if self.age_row_ >= 0 {
                    self.state_.cell_phi[self.age_row_ as usize * unc + uc] = age_seed_link(ctx, link);
                }
                if self.temp_row_ >= 0 {
                    self.state_.cell_phi[self.temp_row_ as usize * unc + uc] = initial_temp_for(ctx, true, link, temp_seed);
                }
            }
        }
        let mut nd = 0;
        while nd < nn && nd < ctx.n_nodes() {
            let und = nd as usize;
            self.node_vol_[und] = store_volume(&ctx.nodes, nd);
            for s in 0..unp {
                self.node_mass_[und * uns + s] = ctx.nodes.conc[und * unp + s] * self.node_vol_[und];
            }
            for m in 0..unm_msx {
                let xi = und * unm_msx + m;
                self.node_mass_[und * uns + unp + m] = (if xi < ctx.reactions.msx_node_conc.len() { ctx.reactions.msx_node_conc[xi] } else { ctx.reactions.init_global[m] }) * self.node_vol_[und];
            }
            if self.age_row_ >= 0 {
                self.node_mass_[und * uns + self.age_row_ as usize] = age_seed_node(ctx, und) * self.node_vol_[und];
            }
            if self.temp_row_ >= 0 {
                self.node_mass_[und * uns + self.temp_row_ as usize] = initial_temp_for(ctx, false, nd, temp_seed) * self.node_vol_[und];
            }
            nd += 1;
        }
        if self.age_row_ >= 0 {
            let (a, b, c) = (ctx.n_nodes(), ctx.n_links(), ctx.n_subcatches());
            ctx.water_age_state.resize(a, b, c);
        }
        if self.temp_row_ >= 0 {
            let (a, b) = (ctx.n_nodes(), ctx.n_links());
            ctx.heat_state.resize(a, b, temp_seed);
            apply_initial_temp_overrides(ctx);
        }

        // E3: dispersion configuration
        self.disp_active_ = false;
        self.disp_mode_ = ArdDispersionMode::Off as i32;
        self.disp_global_ft_ = 0.0;
        assign(&mut self.conduit_disp_ft_, self.mesh_.conduit_link.len(), -1.0);
        assign(&mut self.cell_disp_, unc, 0.0);
        if ctx.ard_config.configured && ctx.ard_config.any_dispersion() {
            let ucf_len = UCF[LENGTH as usize][get_unit_system(ctx.options.flow_units as i32) as usize];
            self.disp_mode_ = ctx.ard_config.dispersion_mode as i32;
            self.disp_global_ft_ = ctx.ard_config.dispersion_value / (ucf_len * ucf_len);
            for i in 0..ctx.ard_config.conduit_disp_link.len() {
                let link = ctx.ard_config.conduit_disp_link[i];
                let d_ft = ctx.ard_config.conduit_disp_value[i] / (ucf_len * ucf_len);
                for r in 0..self.mesh_.conduit_link.len() {
                    if self.mesh_.conduit_link[r] == link {
                        self.conduit_disp_ft_[r] = d_ft;
                    }
                }
            }
            self.disp_active_ = true;
        }

        // E5a: boundary/source rows onto the mesh
        self.bc_node_.clear();
        self.bc_srow_.clear();
        self.bc_ts_.clear();
        self.bc_value_.clear();
        self.bc_now_.clear();
        self.src_crow_.clear();
        self.src_srow_.clear();
        self.src_ts_.clear();
        self.src_value_.clear();
        self.src_len_.clear();
        self.src_now_.clear();
        if nm > 0 {
            let cfg = &ctx.ard_config;
            for i in 0..cfg.bc_node.len() {
                if cfg.bc_msx[i] < 0 || cfg.bc_msx[i] >= nm {
                    continue;
                }
                if cfg.bc_node[i] < 0 || cfg.bc_node[i] >= nn {
                    continue;
                }
                self.bc_node_.push(cfg.bc_node[i]);
                self.bc_srow_.push(np + cfg.bc_msx[i]);
                self.bc_value_.push(cfg.bc_value[i]);
                self.bc_ts_.push(cfg.bc_ts[i]);
            }
            assign(&mut self.bc_now_, self.bc_node_.len(), 0.0);
            for i in 0..cfg.src_link.len() {
                if cfg.src_msx[i] < 0 || cfg.src_msx[i] >= nm {
                    continue;
                }
                let mut crow = -1;
                for r2 in 0..self.mesh_.conduit_link.len() {
                    if self.mesh_.conduit_link[r2] == cfg.src_link[i] {
                        crow = r2 as i32;
                        break;
                    }
                }
                if crow < 0 {
                    continue;
                }
                let mut length = 0.0;
                let b2 = self.mesh_.conduit_cell_begin[crow as usize];
                let n2 = self.mesh_.conduit_cell_count[crow as usize];
                for i2 in 0..n2 {
                    length += self.mesh_.cell_dx[(b2 + i2) as usize];
                }
                self.src_crow_.push(crow);
                self.src_srow_.push(np + cfg.src_msx[i]);
                self.src_value_.push(cfg.src_value[i]);
                self.src_ts_.push(cfg.src_ts[i]);
                self.src_len_.push(cmax(length, 1.0e-6));
            }
            assign(&mut self.src_now_, self.src_crow_.len(), 0.0);
        }

        // E5b: per-cell CSV sidecar
        self.detail_active_ = false;
        self.detail_time_s_ = 0.0;
        if let Some(mut f) = self.detail_out_.take() {
            f.flush();
        }
        if ctx.ard_config.configured && !ctx.ard_config.detailed_output_path.is_empty() {
            match CFile::create(&ctx.ard_config.detailed_output_path, false) {
                Some(mut f) => {
                    f.write_bytes(b"time_s,element,kind,cell,species,conc\n");
                    self.detail_out_ = Some(f);
                    self.detail_active_ = true;
                }
                None => {
                    let mut m = b"[TRANSPORT_OPTIONS] DETAILED_OUTPUT: could not open '".to_vec();
                    m.extend_from_slice(&ctx.ard_config.detailed_output_path);
                    m.extend_from_slice(b"' for writing \xe2\x80\x94 detail output disabled this run.");
                    self.warnings_.push(m);
                }
            }
        }

        self.initialized_ = true;
        true
    }

    /// updateTransportRows: the boundary and source values at the current date.
    fn update_transport_rows(&mut self, ctx: &mut SimulationContext) {
        for i in 0..self.bc_node_.len() {
            let mut c = self.bc_value_[i];
            let ts = self.bc_ts_[i];
            if ts >= 0 && ts < ctx.tables.count() as i32 {
                let d = ctx.current_date;
                c = table_tseries_lookup_cursor(&mut ctx.tables[ts], d);
            }
            self.bc_now_[i] = cmax(0.0, c);
        }
        for i in 0..self.src_crow_.len() {
            let mut r = self.src_value_[i];
            let ts = self.src_ts_[i];
            if ts >= 0 && ts < ctx.tables.count() as i32 {
                let d = ctx.current_date;
                r = table_tseries_lookup_cursor(&mut ctx.tables[ts], d) / K_LITERS_PER_FT3;
            }
            self.src_now_[i] = r;
        }
    }

    /// projectHydraulics: continuity-consistent face fluxes from the link flows and volumes.
    fn project_hydraulics(&mut self, ctx: &SimulationContext, dt: f64) {
        let mesh = &self.mesh_;
        let nf = mesh.n_faces();
        self.face_q_.iter_mut().for_each(|x| *x = 0.0);
        let mut face_set = vec![0u8; nf as usize];

        for r in 0..mesh.conduit_link.len() {
            let link = mesh.conduit_link[r];
            if link < 0 {
                continue;
            }
            let ul = link as usize;
            let b = mesh.conduit_cell_begin[r];
            let n = mesh.conduit_cell_count[r];
            if n <= 0 {
                continue;
            }
            let q = ctx.links.flow[ul];
            let dv = if dt > 0.0 { (ctx.links.volume[ul] - ctx.links.old_volume[ul]) / dt } else { 0.0 };
            for i in 0..=n {
                let face;
                if i < n {
                    let uc = (b + i) as usize;
                    face = if mesh.cell_side0[uc] == 1 { mesh.cell_face0[uc] } else { mesh.cell_face1[uc] };
                } else {
                    let uc = (b + n - 1) as usize;
                    face = if mesh.cell_side0[uc] == 0 { mesh.cell_face0[uc] } else { mesh.cell_face1[uc] };
                }
                if face < 0 || face >= nf {
                    continue;
                }
                let uf = face as usize;
                let frac = if n > 0 { 0.5 - i as f64 / n as f64 } else { 0.0 };
                let fq = q + dv * frac;
                if face_set[uf] != 0 {
                    self.face_q_[uf] = 0.5 * (self.face_q_[uf] + fq);
                } else {
                    self.face_q_[uf] = fq;
                    face_set[uf] = 1;
                }
            }
        }

        for f in 0..nf {
            let uf = f as usize;
            let fq = self.face_q_[uf];
            self.f_mass_[uf] = fq;
            self.f_sstar_[uf] = fq;
            self.f_flux_[uf] = FaceFlux { mass: fq, sstar: fq, ..FaceFlux::default() };
            let cl = mesh.face_cl[uf];
            let cr = mesh.face_cr[uf];
            let mut sl = FaceState::default();
            let mut sr = FaceState::default();
            if cl >= 0 {
                sl.a = self.state_.cell_a[cl as usize];
                sl.q = fq;
            }
            if cr >= 0 {
                sr.a = self.state_.cell_a[cr as usize];
                sr.q = fq;
            }
            self.f_state_l_[uf] = sl;
            self.f_state_r_[uf] = sr;
        }

        let nc = mesh.n_cells();
        for c in 0..nc {
            let uc = c as usize;
            let a = self.state_.cell_a[uc];
            self.cell_u_[uc] = if a > K_DRY_AREA { 0.5 * (self.face_q_[mesh.cell_face0[uc] as usize] + self.face_q_[mesh.cell_face1[uc] as usize]) / a } else { 0.0 };
        }
    }

    /// updateDispersion: the per-cell coefficients for this routing step.
    fn update_dispersion(&mut self, ctx: &SimulationContext) {
        let mode = self.disp_mode_;
        for r in 0..self.mesh_.conduit_link.len() {
            let mut d = 0.0;
            if self.conduit_disp_ft_[r] >= 0.0 {
                d = self.conduit_disp_ft_[r];
            } else if mode == ArdDispersionMode::Value as i32 {
                d = self.disp_global_ft_;
            } else if mode == ArdDispersionMode::Fischer as i32 {
                let link = self.mesh_.conduit_link[r];
                let cr2 = if link >= 0 { ctx.link_subtypes.conduit_row(link) } else { -1 };
                if link >= 0 && cr2 >= 0 {
                    let ul = link as usize;
                    let ucr = cr2 as usize;
                    let len = ctx.link_subtypes.conduits.length[ucr];
                    let vol = ctx.links.volume[ul];
                    let y = ctx.links.depth[ul];
                    let a = if len > 0.0 { vol / len } else { 0.0 };
                    if a > K_DRY_AREA && y > K_MIN_FISCHER_DEPTH {
                        let vel = ctx.links.flow[ul].abs() / a;
                        if vel > 0.0 {
                            let b_w = a / y;
                            let s = cmax(ctx.link_subtypes.conduits.slope[ucr].abs(), K_MIN_SLOPE);
                            let ustar = (K_GRAVITY * y * s).sqrt();
                            d = 0.011 * vel * vel * b_w * b_w / (y * ustar);
                        }
                    }
                }
            }
            let b = self.mesh_.conduit_cell_begin[r];
            let n = self.mesh_.conduit_cell_count[r];
            for i in 0..n {
                self.cell_disp_[(b + i) as usize] = d;
            }
        }
    }

    /// substep: advection, node mixing, sources and dispersion over dt_sub.
    fn substep(&mut self, ctx: &mut SimulationContext, dt_sub: f64, load_frac: f64) -> Result<(), StepError> {
        let ns = self.state_.n_species;
        let nc = self.mesh_.n_cells();
        let nf = self.mesh_.n_faces();
        let uns = ns as usize;
        let unc = nc as usize;
        let unf = nf as usize;

        // 1. Interior reconstruction + FCT through the shared kernels.
        {
            let mut v = SpeciesKernelView {
                mesh: &self.mesh_,
                state: &mut self.state_,
                scalar_scheme: ctx.options.fv.scalar_scheme,
                limiter: ctx.options.fv.limiter,
                hllc: true,
                f_mass: &self.f_mass_,
                f_sstar: &self.f_sstar_,
                f_state_l: &self.f_state_l_,
                f_state_r: &self.f_state_r_,
                f_flux: &self.f_flux_,
                cell_u: &self.cell_u_,
                cell_active: &self.cell_active_,
                active_faces: &self.active_faces_,
                f_phi_l: &mut self.f_phi_l_,
                f_phi_r: &mut self.f_phi_r_,
                f_phi_flux: &mut self.f_phi_flux_,
                cell_slope: &mut self.cell_slope_,
                lo_flux: &mut self.lo_flux_,
                anti_flux: &mut self.anti_flux_,
                td: &mut self.td_,
                anew: &mut self.anew_,
                rplus: &mut self.rplus_,
                rminus: &mut self.rminus_,
            };
            reconstruct_scalars(&mut v, dt_sub);
        }

        let nn = (self.node_mass_.len() / cmax(uns, 1)) as i32;
        let mesh = &self.mesh_;

        // 1a'. Loads at a virtual junction go to the two cells adjoining its spliced face.
        self.has_vj_src_ = false;
        if !mesh.node_vj_face.is_empty() {
            let np_v = ctx.n_pollutants();
            let unp_v = np_v as usize;
            let cell_src_vol = &mut self.cell_src_vol_;
            let cell_src_mass = &mut self.cell_src_mass_;
            let has_vj_src = &mut self.has_vj_src_;
            let ensure = |has: &mut bool, sv: &mut Vec<f64>, sm: &mut Vec<f64>| {
                if *has {
                    return;
                }
                assign(sv, unc, 0.0);
                assign(sm, uns * unc, 0.0);
                *has = true;
            };
            let mut nd = 0;
            while nd < nn && nd < ctx.n_nodes() {
                let und = nd as usize;
                if und >= mesh.node_kind.len() || mesh.node_kind[und] != K_NODE_VIRTUAL {
                    nd += 1;
                    continue;
                }
                let f = mesh.node_vj_face[und];
                if f < 0 {
                    nd += 1;
                    continue;
                }
                let uf = f as usize;
                let cells = [mesh.face_cl[uf], mesh.face_cr[uf]];
                let vol = load_frac * ctx.nodes.qual_vol_in[und];
                if vol != 0.0 {
                    ensure(has_vj_src, cell_src_vol, cell_src_mass);
                    for &c in cells.iter() {
                        if c >= 0 {
                            cell_src_vol[c as usize] += 0.5 * vol;
                        }
                    }
                }
                let credit_mass = |s_row: usize, amount: f64, has: &mut bool, sv: &mut Vec<f64>, sm: &mut Vec<f64>| {
                    if amount == 0.0 {
                        return;
                    }
                    ensure(has, sv, sm);
                    for &c in cells.iter() {
                        if c >= 0 {
                            sm[s_row * unc + c as usize] += 0.5 * amount;
                        }
                    }
                };
                for s in 0..unp_v {
                    credit_mass(s, dt_sub * ctx.nodes.qual_mass_in[und * unp_v + s], has_vj_src, cell_src_vol, cell_src_mass);
                }
                if self.nm_ > 0 && ctx.reactions.msx_ext_mass_in.len() >= ctx.n_nodes() as usize * self.nm_ as usize {
                    for m in 0..self.nm_ as usize {
                        credit_mass(unp_v + m, dt_sub * ctx.reactions.msx_ext_mass_in[und * self.nm_ as usize + m], has_vj_src, cell_src_vol, cell_src_mass);
                    }
                }
                if self.age_row_ >= 0 && und < ctx.water_age_state.node_age_vol_in.len() {
                    credit_mass(self.age_row_ as usize, dt_sub * ctx.water_age_state.node_age_vol_in[und], has_vj_src, cell_src_vol, cell_src_mass);
                }
                if self.temp_row_ >= 0 && und < ctx.heat_state.node_temp_vol_in.len() {
                    credit_mass(self.temp_row_ as usize, dt_sub * ctx.heat_state.node_temp_vol_in[und], has_vj_src, cell_src_vol, cell_src_mass);
                }
                nd += 1;
            }
            for i in 0..self.bc_node_.len() {
                let und = self.bc_node_[i] as usize;
                if und as i32 >= nn || und >= mesh.node_kind.len() || mesh.node_kind[und] != K_NODE_VIRTUAL {
                    continue;
                }
                let f = mesh.node_vj_face[und];
                let vol_ext = load_frac * ctx.nodes.qual_vol_in[und];
                if f >= 0 && vol_ext > 0.0 && self.bc_now_[i] > 0.0 {
                    let uf = f as usize;
                    let amount = vol_ext * self.bc_now_[i];
                    let s_row = self.bc_srow_[i] as usize;
                    if amount != 0.0 {
                        ensure(has_vj_src, cell_src_vol, cell_src_mass);
                        for &c in [mesh.face_cl[uf], mesh.face_cr[uf]].iter() {
                            if c >= 0 {
                                cell_src_mass[s_row * unc + c as usize] += 0.5 * amount;
                            }
                        }
                    }
                }
            }
        }

        // 1b. External loads mixed into the stores before the faces read them.
        {
            let np_l = ctx.n_pollutants();
            let unp_l = np_l as usize;
            let mut nd = 0;
            while nd < nn && nd < ctx.n_nodes() {
                let und = nd as usize;
                if und < mesh.node_kind.len() && mesh.node_kind[und] == K_NODE_VIRTUAL {
                    nd += 1;
                    continue;
                }
                self.node_vol_[und] += load_frac * ctx.nodes.qual_vol_in[und];
                for s in 0..np_l {
                    let mut delta = dt_sub * ctx.nodes.qual_mass_in[und * unp_l + s as usize];
                    let k = und * uns + s as usize;
                    if delta < 0.0 && self.node_mass_[k] + delta < 0.0 {
                        let short = -(self.node_mass_[k] + delta);
                        book_negative_source_clamp(ctx, nd, s, short);
                        delta = -self.node_mass_[k];
                    }
                    self.node_mass_[k] += delta;
                }
                if self.nm_ > 0 && ctx.reactions.msx_ext_mass_in.len() >= ctx.n_nodes() as usize * self.nm_ as usize {
                    for m in 0..self.nm_ as usize {
                        let mut delta = dt_sub * ctx.reactions.msx_ext_mass_in[und * self.nm_ as usize + m];
                        let k = und * uns + unp_l + m;
                        if delta < 0.0 && self.node_mass_[k] + delta < 0.0 {
                            delta = -self.node_mass_[k];
                        }
                        self.node_mass_[k] += delta;
                    }
                }
                if self.age_row_ >= 0 && und < ctx.water_age_state.node_age_vol_in.len() {
                    self.node_mass_[und * uns + self.age_row_ as usize] += dt_sub * ctx.water_age_state.node_age_vol_in[und];
                }
                if self.temp_row_ >= 0 && und < ctx.heat_state.node_temp_vol_in.len() {
                    self.node_mass_[und * uns + self.temp_row_ as usize] += dt_sub * ctx.heat_state.node_temp_vol_in[und];
                }
                nd += 1;
            }
        }

        // 1b(ii). Transport boundaries ride the same external inflow water.
        for i in 0..self.bc_node_.len() {
            let und = self.bc_node_[i] as usize;
            if und as i32 >= nn {
                continue;
            }
            if und < mesh.node_kind.len() && mesh.node_kind[und] == K_NODE_VIRTUAL {
                continue;
            }
            let vol_ext = load_frac * ctx.nodes.qual_vol_in[und];
            if vol_ext > 0.0 && self.bc_now_[i] > 0.0 {
                self.node_mass_[und * uns + self.bc_srow_[i] as usize] += vol_ext * self.bc_now_[i];
            }
        }

        // 1c. Face inflows mixed in before the donor is read.
        for nd in 0..nn {
            let und = nd as usize;
            let b = mesh.node_face_ptr[und];
            let e = mesh.node_face_ptr[und + 1];
            for p in b..e {
                let uf = mesh.node_face_idx[p as usize] as usize;
                let sign = mesh.node_face_sign[p as usize];
                let fq = self.face_q_[uf];
                if sign * fq <= 0.0 {
                    continue;
                }
                self.node_vol_[und] += dt_sub * sign * fq;
                for s in 0..uns {
                    self.node_mass_[und * uns + s] += dt_sub * sign * self.f_phi_flux_[s * unf + uf];
                }
            }
        }

        // 2. Node boundary faces: the node store is the inflow donor.
        for nd in 0..nn {
            let und = nd as usize;
            let b = mesh.node_face_ptr[und];
            let e = mesh.node_face_ptr[und + 1];
            for p in b..e {
                let uf = mesh.node_face_idx[p as usize] as usize;
                let sign = mesh.node_face_sign[p as usize];
                let fq = self.face_q_[uf];
                if sign * fq < 0.0 {
                    let vol = self.node_vol_[und];
                    let fresh_outfall = ctx.options.outfall_backflow_zero && ctx.nodes.type_[und] == NodeType::Outfall;
                    for s in 0..uns {
                        let cnode = if !fresh_outfall && vol > K_MIN_STORE_VOL { self.node_mass_[und * uns + s] / vol } else { 0.0 };
                        self.f_phi_flux_[s * unf + uf] = fq * cnode;
                    }
                }
            }
        }

        // 3. Cell update: species mass and area from the same fluxes.
        for c in 0..nc {
            let uc = c as usize;
            let faces = [mesh.cell_face0[uc], mesh.cell_face1[uc]];
            let sides = [mesh.cell_side0[uc], mesh.cell_side1[uc]];
            let inv_dx = 1.0 / mesh.cell_dx[uc];
            let mut d_a = 0.0;
            for e2 in 0..2 {
                let uf = faces[e2] as usize;
                let sg = if sides[e2] == 0 { -1.0 } else { 1.0 };
                d_a += sg * self.f_mass_[uf];
            }
            let a_old = self.state_.cell_a[uc];
            let mut a_new = cmax(0.0, a_old + dt_sub * d_a * inv_dx);
            if self.has_vj_src_ {
                a_new = cmax(0.0, a_new + self.cell_src_vol_[uc] * inv_dx);
            }
            for s in 0..ns {
                let sb = s as usize;
                let mut dm = 0.0;
                for e2 in 0..2 {
                    let uf = faces[e2] as usize;
                    let sg = if sides[e2] == 0 { -1.0 } else { 1.0 };
                    dm += sg * self.f_phi_flux_[sb * unf + uf];
                }
                let mut m = a_old * self.state_.cell_phi[sb * unc + uc] + dt_sub * dm * inv_dx;
                if self.has_vj_src_ {
                    m += self.cell_src_mass_[sb * unc + uc] * inv_dx;
                }
                let m_f = if s == self.temp_row_ { cclamp(m, K_TEMP_FLOOR_C * a_new, K_TEMP_CEIL_C * a_new) } else { cmax(0.0, m) };
                self.state_.cell_phi[sb * unc + uc] = if a_new > K_DRY_AREA { m_f / a_new } else { self.state_.cell_phi[sb * unc + uc] };
            }
            self.state_.cell_a[uc] = a_new;
        }

        // 4. Node stores: the outflow half of the face exchange.
        let mut nd = 0;
        while nd < nn && nd < ctx.n_nodes() {
            let und = nd as usize;
            let b = mesh.node_face_ptr[und];
            let e = mesh.node_face_ptr[und + 1];
            let mut dvol = 0.0;
            for p in b..e {
                let uf = mesh.node_face_idx[p as usize] as usize;
                let sign = mesh.node_face_sign[p as usize];
                let fq = self.face_q_[uf];
                if sign * fq > 0.0 {
                    continue;
                }
                dvol += sign * fq;
                for s in 0..uns {
                    self.node_mass_[und * uns + s] += dt_sub * sign * self.f_phi_flux_[s * unf + uf];
                }
            }
            self.node_vol_[und] = cmax(0.0, self.node_vol_[und] + dt_sub * dvol);
            for s in 0..ns {
                let k = und * uns + s as usize;
                if s == self.temp_row_ {
                    let v = self.node_vol_[und];
                    self.node_mass_[k] = cclamp(self.node_mass_[k], K_TEMP_FLOOR_C * v, K_TEMP_CEIL_C * v);
                } else {
                    self.node_mass_[k] = cmax(0.0, self.node_mass_[k]);
                }
            }
            nd += 1;
        }

        // 5. Structures: zero-volume passthrough of the donor node store.
        let n_struct = mesh.struct_link.len() as i32;
        for i in 0..n_struct {
            let ui = i as usize;
            let link = mesh.struct_link[ui];
            if link < 0 || link >= ctx.n_links() {
                continue;
            }
            let q = ctx.links.flow[link as usize];
            if q == 0.0 {
                continue;
            }
            let donor_n = if q >= 0.0 { mesh.struct_n1[ui] } else { mesh.struct_n2[ui] };
            let receiver_n = if q >= 0.0 { mesh.struct_n2[ui] } else { mesh.struct_n1[ui] };
            if donor_n < 0 || receiver_n < 0 || donor_n >= nn || receiver_n >= nn {
                continue;
            }
            let ud = donor_n as usize;
            let ur = receiver_n as usize;
            let dv = q.abs() * dt_sub;
            let donor_vol = self.node_vol_[ud];
            if donor_vol > K_MIN_STORE_VOL {
                for sb in 0..uns {
                    let c_donor = self.node_mass_[ud * uns + sb] / donor_vol;
                    let dm = cmin(dv * c_donor, self.node_mass_[ud * uns + sb]);
                    self.node_mass_[ud * uns + sb] -= dm;
                    self.node_mass_[ur * uns + sb] += dm;
                }
            }
            self.node_vol_[ud] = cmax(0.0, self.node_vol_[ud] - dv);
            self.node_vol_[ur] += dv;
        }

        // Booster stage
        if !ctx.reactions.sources.is_empty() {
            for nd in 0..nn {
                if mesh.node_kind[nd as usize] == K_NODE_VIRTUAL {
                    continue;
                }
                let mut arriving = cmax(0., load_frac * ctx.nodes.qual_vol_in[nd as usize]);
                for p in mesh.node_face_ptr[nd as usize]..mesh.node_face_ptr[nd as usize + 1] {
                    let f = mesh.node_face_idx[p as usize] as usize;
                    arriving += dt_sub * cmax(0., mesh.node_face_sign[p as usize] * self.face_q_[f]);
                }
                for j in 0..n_struct as usize {
                    let link = mesh.struct_link[j];
                    let q = ctx.links.flow[link as usize];
                    if (if q >= 0.0 { mesh.struct_n2[j] } else { mesh.struct_n1[j] }) == nd {
                        arriving += q.abs() * dt_sub;
                    }
                }
                for m in 0..self.nm_ {
                    let k = nd as usize * uns + ctx.n_pollutants() as usize + m as usize;
                    let add = reaction_booster_mass(ctx, nd, m, self.node_vol_[nd as usize], arriving)?;
                    self.node_mass_[k] += add;
                }
            }
        }

        // 5c. Distributed sources over the conduit's wet cells.
        for i in 0..self.src_crow_.len() {
            let r = self.src_now_[i];
            if r == 0.0 {
                continue;
            }
            let ucrow = self.src_crow_[i] as usize;
            let b2 = mesh.conduit_cell_begin[ucrow];
            let n2 = mesh.conduit_cell_count[ucrow];
            let sb = self.src_srow_[i] as usize;
            let mut shortfall = 0.0;
            for i2 in 0..n2 {
                let uc = (b2 + i2) as usize;
                let a = self.state_.cell_a[uc];
                if a <= K_DRY_AREA {
                    continue;
                }
                let mut dphi = dt_sub * r / (self.src_len_[i] * a);
                let phi = &mut self.state_.cell_phi[sb * unc + uc];
                if dphi < 0.0 && *phi + dphi < 0.0 {
                    shortfall += -(*phi + dphi) * a * mesh.cell_dx[uc];
                    dphi = -*phi;
                }
                *phi += dphi;
            }
            if shortfall > 0.0 {
                book_negative_cell_source_clamp(ctx, self.src_crow_[i], shortfall);
            }
        }

        // 6. Dispersion: implicit per-chain solve over the updated cells.
        let mut dv = DispersionView {
            mesh: &self.mesh_,
            state: &mut self.state_,
            dispersion: 0.0,
            cell_dispersion: if self.disp_active_ { Some(self.cell_disp_.as_slice()) } else { None },
        };
        dispersion_solve(&mut dv, dt_sub);
        Ok(())
    }

    /// ArdEngine::step: one routing step of transport, aging, heat, reactions and the bed,
    /// then publish.
    pub fn step(&mut self, ctx: &mut SimulationContext, dt: f64) -> Result<(), StepError> {
        if !self.initialized_ || dt <= 0.0 {
            return Ok(());
        }
        ctx.links.conc_old.clone_from(&ctx.links.conc);
        ctx.nodes.conc_old.clone_from(&ctx.nodes.conc);

        self.project_hydraulics(ctx, dt);
        if self.disp_active_ {
            self.update_dispersion(ctx);
        }
        if !self.bc_node_.is_empty() || !self.src_crow_.is_empty() {
            self.update_transport_rows(ctx);
        }

        // Node volumes resync to the solver's each routing step.
        let ns = self.state_.n_species;
        let nn = cmin(ctx.n_nodes(), self.node_vol_.len() as i32);
        let uns = ns as usize;
        for nd in 0..nn {
            let und = nd as usize;
            let v_old = self.node_vol_[und];
            let v_new = store_volume(&ctx.nodes, nd);
            if v_old > K_MIN_STORE_VOL {
                let ratio = v_new / v_old;
                if ratio < 1.0 || ratio - 1.0 <= K_RESYNC_NOISE_REL {
                    for s in 0..uns {
                        self.node_mass_[und * uns + s] *= ratio;
                    }
                }
            } else if v_new <= K_MIN_STORE_VOL {
                for s in 0..uns {
                    self.node_mass_[und * uns + s] = 0.0;
                }
            }
            self.node_vol_[und] = v_new;
        }

        // CFL subcycling on the projected velocities.
        let mut dt_cfl = dt;
        let nc = self.mesh_.n_cells();
        for c in 0..nc {
            let uc = c as usize;
            let u = self.cell_u_[uc].abs();
            if u > 1.0e-12 {
                dt_cfl = cmin(dt_cfl, K_CFL * self.mesh_.cell_dx[uc] / u);
            }
        }
        let mut nsub = cvt_i32((dt / cmax(dt_cfl, 1.0e-12)).ceil());
        if nsub > K_MAX_SUBSTEPS {
            if !self.warned_cfl_clamp_ {
                self.warned_cfl_clamp_ = true;
                ctx.warnings.push(format!("QUALITY_SOLVER EULERIAN_ARD: transport subcycling clamped at {} substeps (needed {}) \u{2014} effective CFL exceeds the stability target; shorten ROUTING_STEP or increase FV_CELL_LENGTH.", K_MAX_SUBSTEPS, nsub).into_bytes());
            }
            nsub = K_MAX_SUBSTEPS;
        }
        let dt_sub = dt / nsub as f64;
        let load_frac = 1.0 / nsub as f64;

        for _ in 0..nsub {
            self.substep(ctx, dt_sub, load_frac)?;
        }

        // A1a: aging
        if self.age_row_ >= 0 {
            let ar = self.age_row_ as usize;
            let unc0 = self.mesh_.n_cells() as usize;
            for c in 0..unc0 {
                self.state_.cell_phi[ar * unc0 + c] += dt;
            }
            let uns0 = self.state_.n_species as usize;
            for nd2 in 0..self.node_vol_.len() {
                self.node_mass_[nd2 * uns0 + ar] += dt * self.node_vol_[nd2];
            }
        }

        // H4: surface heat fluxes
        self.apply_heat_fluxes(ctx, dt)?;

        // E4: reaction stage
        update_untreated_storage_hrt(ctx, dt);
        {
            let np = ctx.n_pollutants();
            let arrays = ArdArrays {
                cell_phi: &mut self.state_.cell_phi,
                cell_a: &self.state_.cell_a,
                cell_dx: &self.mesh_.cell_dx,
                n_cells: self.mesh_.n_cells(),
                node_mass: &mut self.node_mass_,
                node_vol: &self.node_vol_,
                n_nodes: self.node_vol_.len() as i32,
                cell_conduit: &self.mesh_.cell_conduit,
                conduit_link: &self.mesh_.conduit_link,
            };
            react_ard_stage(ctx, dt, arrays, np, self.state_.n_species, K_MIN_STORE_VOL, self.temp_row_)?;
        }

        // H6b-ARD: bed solute exchange
        self.apply_bed_solute_exchange(ctx, dt)?;

        self.publish(ctx);

        if self.detail_active_ {
            self.detail_time_s_ += dt;
            self.write_detail_rows(ctx);
        }
        Ok(())
    }

    /// absorbTreatedNodeConc: the treated node concentrations back into the node stores.
    pub fn absorb_treated_node_conc(&mut self, ctx: &SimulationContext) {
        if !self.initialized_ {
            return;
        }
        let np = ctx.n_pollutants();
        if np <= 0 || !ctx.treatment.has_any() {
            return;
        }
        let uns = self.state_.n_species as usize;
        let unp = np as usize;
        let nn = cmin(ctx.n_nodes(), self.node_vol_.len() as i32);
        for nd in 0..nn {
            let und = nd as usize;
            if !ctx.treatment.has_treatment[und] {
                continue;
            }
            let vol = self.node_vol_[und];
            if vol <= K_MIN_STORE_VOL {
                continue;
            }
            for p in 0..unp {
                self.node_mass_[und * uns + p] = ctx.nodes.conc[und * unp + p] * vol;
            }
        }
    }

    /// writeDetailRows: every cell and node concentration to the CSV sidecar.
    fn write_detail_rows(&mut self, ctx: &SimulationContext) {
        let ns = self.state_.n_species;
        let np = ctx.n_pollutants();
        let unc = self.mesh_.n_cells() as usize;
        let uns = ns as usize;
        let (age_row, temp_row) = (self.age_row_, self.temp_row_);
        let species_name = |s: i32| -> Bytes {
            if s == age_row {
                return b"__WATER_AGE__".to_vec();
            }
            if s == temp_row {
                return b"__TEMPERATURE__".to_vec();
            }
            if s < np {
                ctx.pollutant_names.name_of(s).to_vec()
            } else {
                ctx.reactions.species_name[(s - np) as usize].clone()
            }
        };
        let t = printf_g(self.detail_time_s_, 6);
        let Some(out) = self.detail_out_.as_mut() else { return };
        for r in 0..self.mesh_.conduit_link.len() {
            let link = self.mesh_.conduit_link[r];
            if link < 0 || link >= ctx.n_links() {
                continue;
            }
            let lname = &ctx.link_names.names()[link as usize];
            let b = self.mesh_.conduit_cell_begin[r];
            let n = self.mesh_.conduit_cell_count[r];
            for i in 0..n {
                let uc = (b + i) as usize;
                for s in 0..ns {
                    let mut line = t.clone();
                    line.push(b',');
                    line.extend_from_slice(lname);
                    line.extend_from_slice(b",L,");
                    line.extend_from_slice(i.to_string().as_bytes());
                    line.push(b',');
                    line.extend_from_slice(&species_name(s));
                    line.push(b',');
                    line.extend_from_slice(&printf_g(self.state_.cell_phi[s as usize * unc + uc], 6));
                    line.push(b'\n');
                    out.write_bytes(&line);
                }
            }
        }
        let nn = cmin(ctx.n_nodes(), self.node_vol_.len() as i32);
        for nd in 0..nn {
            let und = nd as usize;
            let vol = self.node_vol_[und];
            let nname = &ctx.node_names.names()[und];
            for s in 0..ns {
                let c = if vol > K_MIN_STORE_VOL { self.node_mass_[und * uns + s as usize] / vol } else { 0.0 };
                let mut line = t.clone();
                line.push(b',');
                line.extend_from_slice(nname);
                line.extend_from_slice(b",N,0,");
                line.extend_from_slice(&species_name(s));
                line.push(b',');
                line.extend_from_slice(&printf_g(c, 6));
                line.push(b'\n');
                out.write_bytes(&line);
            }
        }
    }

    /// cellLink: the link a mesh cell belongs to, or -1.
    pub fn cell_link(&self, cell: usize) -> i32 {
        if cell >= self.mesh_.cell_conduit.len() {
            return -1;
        }
        let cr = self.mesh_.cell_conduit[cell];
        if cr < 0 || cr as usize >= self.mesh_.conduit_link.len() {
            return -1;
        }
        self.mesh_.conduit_link[cr as usize]
    }

    /// applyBedSoluteExchange: per-cell bed solute exchange (pollutant and reactions rows).
    fn apply_bed_solute_exchange(&mut self, ctx: &mut SimulationContext, dt: f64) -> Result<(), StepError> {
        if !bed_exchange_enabled(ctx) || !(dt > 0.0) {
            return Ok(());
        }
        let mut n_exch = self.state_.n_species;
        if self.age_row_ >= 0 {
            n_exch -= 1;
        }
        if self.temp_row_ >= 0 {
            n_exch -= 1;
        }
        if n_exch <= 0 {
            return Ok(());
        }
        let ncells = self.mesh_.n_cells();
        {
            let bed = &mut ctx.bed_state;
            if bed.cell_temp.len() != ncells as usize || bed.cell_n_species != n_exch {
                let sd = &ctx.heat_config.sediment;
                let t_keep = if bed.cell_temp.len() == ncells as usize && !bed.cell_temp.is_empty() {
                    bed.cell_temp[0]
                } else if sd.has_initial_temp {
                    sd.initial_temp
                } else {
                    sd.ground_temp
                };
                let mut keep_temp: Vec<f64> = Vec::new();
                if bed.cell_temp.len() == ncells as usize {
                    keep_temp = bed.cell_temp.clone();
                }
                bed.resize_cells(ncells, n_exch, t_keep);
                if !keep_temp.is_empty() {
                    bed.cell_temp = keep_temp;
                }
            }
        }
        const K_SQ_FT_TO_SQ_M: f64 = 0.09290304;
        const K_CU_FT_TO_CU_M: f64 = 0.028316846592;
        let unc = ncells as usize;
        let une = n_exch as usize;
        for c in 0..unc {
            let mesh_link = self.cell_link(c);
            if mesh_link >= 0 && (mesh_link as usize) < ctx.heat_state.mesh_contact_link.len() && ctx.heat_state.mesh_contact_link[mesh_link as usize] != 0 {
                continue;
            }
            let gi = self.mesh_.cell_geom[c] as usize;
            if gi >= self.mesh_.geom.len() {
                continue;
            }
            let g = &self.mesh_.geom[gi];
            let area_x = self.state_.cell_a[c];
            if !(area_x > 0.0) {
                continue;
            }
            let vol_ft3 = area_x * self.mesh_.cell_dx[c];
            if !(vol_ft3 > 0.0) {
                continue;
            }
            let h = depth_of_area(g, area_x);
            if !(h > 0.0) {
                continue;
            }
            let rh = hyd_rad_of_depth(g, h);
            if !(rh > 0.0) {
                continue;
            }
            let bed_m2 = (area_x / rh) * self.mesh_.cell_dx[c] * K_SQ_FT_TO_SQ_M;
            let sd_c = sediment_for(ctx, &HeatElement::link(self.cell_link(c))).clone();
            let vol_w = vol_ft3 * K_CU_FT_TO_CU_M;
            let vol_b = bed_m2 * sd_c.bed_thickness;
            let q_exch = bed_exchange_q(&sd_c, bed_m2);
            if !(q_exch > 0.0) || !(vol_b > 0.0) {
                continue;
            }
            for sp in 0..une {
                let cw = self.state_.cell_phi[sp * unc + c];
                let cb = ctx.bed_state.cell_conc[sp * unc + c];
                let st = exchange_pair(cw, cb, vol_w, vol_b, q_exch, dt);
                self.state_.cell_phi[sp * unc + c] += st.dc_w;
                ctx.bed_state.cell_conc[sp * unc + c] += st.dc_b;
            }
        }
        Ok(())
    }

    /// applyHeatFluxes: per-cell surface exchange and bed pair, and the storage-node stores.
    fn apply_heat_fluxes(&mut self, ctx: &mut SimulationContext, dt: f64) -> Result<(), StepError> {
        if self.temp_row_ < 0 || !(dt > 0.0) {
            return Ok(());
        }
        let bed_on = bed_exchange_enabled(ctx);
        let any = ctx.heat_config.surface_exchange || ctx.heat_config.radiative_exchange || bed_on;
        if !any {
            return Ok(());
        }
        let mut t_gr = 0.0;
        if bed_on {
            let ncells = self.mesh_.n_cells();
            let bed = &mut ctx.bed_state;
            if bed.cell_temp.len() != ncells as usize || !bed.cells_seeded {
                let sd = &ctx.heat_config.sediment;
                let n = bed.cell_n_species;
                bed.resize_cells(ncells, n, if sd.has_initial_temp { sd.initial_temp } else { sd.ground_temp });
            }
            t_gr = ground_temperature(ctx)?;
        }
        update_solar_forcing(ctx, dt)?;
        const K_SQ_FT_TO_SQ_M: f64 = 0.09290304;
        const K_CU_FT_TO_CU_M: f64 = 0.028316846592;
        let rho = ctx.options.water_density;
        let cp = ctx.options.water_specific_heat;
        let tr = self.temp_row_ as usize;
        let unc = self.mesh_.n_cells() as usize;
        for c in 0..unc {
            let mesh_link = self.cell_link(c);
            if mesh_link >= 0 && (mesh_link as usize) < ctx.heat_state.mesh_contact_link.len() && ctx.heat_state.mesh_contact_link[mesh_link as usize] != 0 {
                continue;
            }
            let gi = self.mesh_.cell_geom[c] as usize;
            if gi >= self.mesh_.geom.len() {
                continue;
            }
            let g = &self.mesh_.geom[gi];
            let area_x = self.state_.cell_a[c];
            if !(area_x > 0.0) {
                continue;
            }
            let vol_ft3 = area_x * self.mesh_.cell_dx[c];
            if !(vol_ft3 > 0.0) {
                continue;
            }
            let h = depth_of_area(g, area_x);
            if !(h > 0.0) {
                continue;
            }
            let mut surf_m2 = 0.0;
            if g.is_open != 0 {
                let width = width_of_depth(g, h);
                if width > 0.0 {
                    surf_m2 = width * self.mesh_.cell_dx[c] * K_SQ_FT_TO_SQ_M;
                }
            }
            let t_w = self.state_.cell_phi[tr * unc + c];
            let ce = HeatElement::link(self.cell_link(c));
            if bed_on && c < ctx.bed_state.cell_temp.len() {
                let rh = hyd_rad_of_depth(g, h);
                let bed_m2 = if rh > 0.0 { (area_x / rh) * self.mesh_.cell_dx[c] * K_SQ_FT_TO_SQ_M } else { 0.0 };
                let sd_c = sediment_for(ctx, &ce);
                let bc = bed_coupling_from_contact(ctx, sd_c, bed_m2, vol_ft3, ground_temp_for(ctx, sd_c, t_gr));
                if bc.viable() {
                    let ps = relax_pair(&bc, t_w, ctx.bed_state.cell_temp[c], net_flux_out(ctx, &ce, t_w), net_flux_out(ctx, &ce, t_w + K_PROBE_C), K_PROBE_C, surf_m2, dt);
                    self.state_.cell_phi[tr * unc + c] += ps.dt_w;
                    ctx.bed_state.cell_temp[c] += ps.dt_b;
                    continue;
                }
            }
            if !(surf_m2 > 0.0) {
                continue;
            }
            self.state_.cell_phi[tr * unc + c] += relax_t(net_flux_out(ctx, &ce, t_w), net_flux_out(ctx, &ce, t_w + K_PROBE_C), K_PROBE_C, surf_m2, vol_ft3 * K_CU_FT_TO_CU_M, dt, rho, cp);
        }
        let unit_sys = get_unit_system(ctx.options.flow_units as i32);
        let uns = self.state_.n_species as usize;
        let nn = cmin(ctx.n_nodes(), self.node_vol_.len() as i32);
        for nd in 0..nn {
            let und = nd as usize;
            let vol_ft3 = self.node_vol_[und];
            if !(vol_ft3 > 0.0) {
                continue;
            }
            let area_ft2 = get_surf_area(&ctx.nodes, nd, ctx.nodes.depth[und], Some(&ctx.tables), unit_sys, Some(&ctx.node_subtypes));
            if !(area_ft2 > 0.0) {
                continue;
            }
            let t_w = self.node_mass_[und * uns + tr] / vol_ft3;
            let ne = HeatElement::node(nd);
            let d_t = relax_t(net_flux_out(ctx, &ne, t_w), net_flux_out(ctx, &ne, t_w + K_PROBE_C), K_PROBE_C, area_ft2 * K_SQ_FT_TO_SQ_M, vol_ft3 * K_CU_FT_TO_CU_M, dt, rho, cp);
            self.node_mass_[und * uns + tr] += d_t * vol_ft3;
        }
        Ok(())
    }

    /// publish: the conduit means into links.conc and the stores into nodes.conc.
    fn publish(&mut self, ctx: &mut SimulationContext) {
        let ns = self.state_.n_species;
        let np = ctx.n_pollutants();
        let na = if self.age_row_ >= 0 { 1 } else { 0 };
        let nt = if self.temp_row_ >= 0 { 1 } else { 0 };
        let nm = ns - np - na - nt;
        let uns = ns as usize;
        let unp = np as usize;
        let unm = if nm > 0 { nm as usize } else { 0 };
        let unc = self.mesh_.n_cells() as usize;

        if nm > 0 {
            let want_l = ctx.n_links() as usize * unm;
            let want_n = ctx.n_nodes() as usize * unm;
            let rx = &mut ctx.reactions;
            if rx.msx_link_conc.len() != want_l {
                assign(&mut rx.msx_link_conc, want_l, 0.0);
            }
            if rx.msx_node_conc.len() != want_n {
                assign(&mut rx.msx_node_conc, want_n, 0.0);
            }
        }

        for r in 0..self.mesh_.conduit_link.len() {
            let link = self.mesh_.conduit_link[r];
            if link < 0 || link >= ctx.n_links() {
                continue;
            }
            let ul = link as usize;
            let b = self.mesh_.conduit_cell_begin[r];
            let n = self.mesh_.conduit_cell_count[r];
            for s in 0..ns {
                let sb = s as usize;
                let mut m = 0.0;
                let mut vol = 0.0;
                for i in 0..n {
                    let uc = (b + i) as usize;
                    let dv = self.state_.cell_a[uc] * self.mesh_.cell_dx[uc];
                    m += dv * self.state_.cell_phi[sb * unc + uc];
                    vol += dv;
                }
                let conc = if vol > 1.0e-12 { m / vol } else { 0.0 };
                if s == self.age_row_ {
                    if ul < ctx.water_age_state.link_age.len() {
                        ctx.water_age_state.link_age[ul] = conc;
                    }
                } else if s == self.temp_row_ {
                    if ul < ctx.heat_state.link_temp.len() && vol > K_ZERO_VOLUME_FT3 {
                        ctx.heat_state.link_temp[ul] = conc;
                    }
                } else if s < np {
                    ctx.links.conc[ul * unp + sb] = conc;
                } else {
                    ctx.reactions.msx_link_conc[ul * unm + (s - np) as usize] = conc;
                }
            }
        }

        let nn = cmin(ctx.n_nodes(), self.node_vol_.len() as i32);
        for nd in 0..nn {
            let und = nd as usize;
            let vol = self.node_vol_[und];
            for s in 0..ns {
                let conc = if vol > 1.0e-12 { self.node_mass_[und * uns + s as usize] / vol } else { 0.0 };
                if s == self.age_row_ {
                    if und < ctx.water_age_state.node_age.len() {
                        ctx.water_age_state.node_age[und] = conc;
                    }
                } else if s == self.temp_row_ {
                    if und < ctx.heat_state.node_temp.len() && vol > K_ZERO_VOLUME_FT3 {
                        ctx.heat_state.node_temp[und] = conc;
                    }
                } else if s < np {
                    ctx.nodes.conc[und * unp + s as usize] = conc;
                } else {
                    ctx.reactions.msx_node_conc[und * unm + (s - np) as usize] = conc;
                }
            }
        }

        let n_struct = self.mesh_.struct_link.len();
        for ui in 0..n_struct {
            let link = self.mesh_.struct_link[ui];
            if link < 0 || link >= ctx.n_links() {
                continue;
            }
            let ul = link as usize;
            let q = ctx.links.flow[ul];
            let donor_n = if q >= 0.0 { self.mesh_.struct_n1[ui] } else { self.mesh_.struct_n2[ui] };
            if donor_n < 0 || donor_n >= nn {
                continue;
            }
            let ud = donor_n as usize;
            for s in 0..unp {
                ctx.links.conc[ul * unp + s] = ctx.nodes.conc[ud * unp + s];
            }
            for m2 in 0..(if nm > 0 { nm as usize } else { 0 }) {
                ctx.reactions.msx_link_conc[ul * unm + m2] = ctx.reactions.msx_node_conc[ud * unm + m2];
            }
            if self.age_row_ >= 0 && ul < ctx.water_age_state.link_age.len() && ud < ctx.water_age_state.node_age.len() {
                ctx.water_age_state.link_age[ul] = ctx.water_age_state.node_age[ud];
            }
            if self.temp_row_ >= 0 && ul < ctx.heat_state.link_temp.len() && ud < ctx.heat_state.node_temp.len() {
                ctx.heat_state.link_temp[ul] = ctx.heat_state.node_temp[ud];
            }
        }
    }

    /// totalMass: species inventory over the cells and node stores.
    pub fn total_mass(&self) -> Vec<f64> {
        let ns = self.state_.n_species;
        let uns = ns as usize;
        let unc = self.mesh_.n_cells() as usize;
        let mut total = vec![0.0; uns];
        for sb in 0..uns {
            let mut m = 0.0;
            for c in 0..unc {
                m += self.state_.cell_a[c] * self.mesh_.cell_dx[c] * self.state_.cell_phi[sb * unc + c];
            }
            for nd in 0..self.node_vol_.len() {
                m += self.node_mass_[nd * uns + sb];
            }
            total[sb] = m;
        }
        total
    }
}
