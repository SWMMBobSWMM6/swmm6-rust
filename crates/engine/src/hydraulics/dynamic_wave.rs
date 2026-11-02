// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/DynamicWave.hpp and DynamicWave.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The dynamic-wave solver's state (DWSolver and its SoA blocks) and its initialization:
// DWSolver::init, refreshConduitTile, buildConduitNodeCSR, buildVirtualJunctionPairs,
// setNumThreads, and the helpers they call (isOpenShape, the file-local buildXSP).
//
// NOT YET TRANSLATED (run-time step functions, with the routing step): execute,
// getRoutingStep, accumulateLinkStepStats, getCrownCutoff, getSlotWidth, applyDPSGeometry,
// updateTpaLatch, applyTpaGeometry, updateDPSState, spatialSmoothP, vjPrepareIteration,
// vjAccumulateResiduals, initNodeStates, findBypassedLinks, computeLinkGeometry,
// recomputeConduitLossOne, limitConduitLosses, momentumKernels, processDryLink,
// applyFlowLimits, processManningLink, traceLinkTerms, processForceMainLink,
// updateNodeFlows, gatherConduitNodeFlows, computeAASkipFlags, updateNodeDepthsTeam,
// setNodeDepth, commitNodeDepthState, getLinkStep, and the file-local conveyArea,
// chamberEndArea, chamberPipeVolumeChange and computeYnorm.
//
// Differences in form, never in result:
//   * `const XSectGroups* groups_` (a pointer into the owning Router) is not stored: the
//     groups are owned by the Router next to this solver, and the step functions will
//     take them as an argument. init's `groups` argument is therefore only recorded as
//     attached (`groups_attached_`).
//   * The members keep upstream's names, trailing underscores included, and are `pub`
//     (upstream's private ones too) so the gate's state dumper can read them, as
//     upstream's harness does through `#define private public`.
//   * `std::vector<int8_t>` is Vec<i8>, `std::vector<uint8_t>` Vec<u8>, `XsectShape` and
//     `MomentumCategory` vectors hold the enums.
//   * setNumThreads reads SWMM_DW_THREADS through atoi, as upstream's std::atoi.

use crate::core::cmath::cmax;
use crate::core::constants;
use crate::core::context::SimulationContext;
use crate::core::options::{Bytes, NodeContinuity};
use crate::core::thread_info;
use crate::core::unit_conversion as ucf;
use crate::data::link::{LinkType, XsectShape};
use crate::data::node::NodeType;
use crate::hydraulics::link;
use crate::hydraulics::transect_types::N_TRANSECT_TBL;
use crate::hydraulics::xsect_batch::{XSectGroups, XSectParams};
use crate::hydraulics::xsection as xsect;

use constants::GRAVITY;

/// Configuration for the Dynamic Preissmann Slot (Sharior et al. 2023).
#[allow(non_snake_case)]
#[derive(Debug, Clone, Copy)]
pub struct DPSConfig {
    /// Target pressure celerity (ft/s, internal units)
    pub c_pT: f64,
    /// Surcharge shock parameter (>= 2)
    pub alpha: f64,
    /// Decay time scale for P -> 1 (seconds)
    pub r: f64,
    /// c_pT^2 (pre-computed)
    pub c_pT_sq: f64,
}

impl Default for DPSConfig {
    fn default() -> Self {
        DPSConfig { c_pT: 25.0, alpha: 3.0, r: 0.5, c_pT_sq: 625.0 }
    }
}

/// Per-conduit DPS state (SoA).
#[allow(non_snake_case)]
#[derive(Debug, Clone, Default)]
pub struct DPSLinkArrays {
    pub As: Vec<f64>,
    pub hs: Vec<f64>,
    pub hs_prev_iter: Vec<f64>,
    pub P: Vec<f64>,
    pub P_hat: Vec<f64>,
    pub P_hat_0: Vec<f64>,
    pub T_s_target: Vec<f64>,
    pub t_s: Vec<f64>,
    pub surcharged: Vec<u8>,
}

impl DPSLinkArrays {
    pub fn resize(&mut self, n: usize) {
        fn a<T: Clone>(v: &mut Vec<T>, n: usize, d: T) {
            v.clear();
            v.resize(n, d);
        }
        a(&mut self.As, n, 0.0);
        a(&mut self.hs, n, 0.0);
        a(&mut self.hs_prev_iter, n, 0.0);
        a(&mut self.P, n, 1.0);
        a(&mut self.P_hat, n, 1.0);
        a(&mut self.P_hat_0, n, 1.0);
        a(&mut self.T_s_target, n, 0.0);
        a(&mut self.t_s, n, 0.0);
        a(&mut self.surcharged, n, 0);
    }
}

/// Per-node dynamic-wave working arrays (SoA).
#[allow(non_snake_case)]
#[derive(Debug, Clone, Default)]
pub struct DWNodeArrays {
    pub new_surf_area: Vec<f64>,
    pub old_surf_area: Vec<f64>,
    pub sumdqdh: Vec<f64>,
    pub dYdT: Vec<f64>,
    pub conduit_loss: Vec<f64>,
    pub loss_factor: Vec<f64>,
    pub converged: Vec<u8>,
    pub is_surcharged: Vec<u8>,
}

impl DWNodeArrays {
    pub fn resize(&mut self, n: usize) {
        fn a<T: Clone>(v: &mut Vec<T>, n: usize, d: T) {
            v.clear();
            v.resize(n, d);
        }
        a(&mut self.new_surf_area, n, 0.0);
        a(&mut self.old_surf_area, n, 0.0);
        a(&mut self.sumdqdh, n, 0.0);
        a(&mut self.dYdT, n, 0.0);
        a(&mut self.conduit_loss, n, 0.0);
        a(&mut self.loss_factor, n, 1.0);
        a(&mut self.converged, n, 0);
        a(&mut self.is_surcharged, n, 0);
    }
}

#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SurchargeMethod {
    #[default]
    Extran = 0,
    Slot = 1,
    DynamicSlot = 2,
    Tpa = 3,
}

impl SurchargeMethod {
    /// static_cast<SurchargeMethod>(int) for the values the [OPTIONS] parser can store
    /// (0..3). Upstream would keep any other number and match no `==` test; this maps it to
    /// Extran, which only differs from upstream in the getCrownCutoff branch of the step.
    pub fn from_i32(v: i32) -> SurchargeMethod {
        match v {
            1 => SurchargeMethod::Slot,
            2 => SurchargeMethod::DynamicSlot,
            3 => SurchargeMethod::Tpa,
            _ => SurchargeMethod::Extran,
        }
    }
}

#[repr(u8)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum MomentumCategory {
    #[default]
    SkipDry = 0,
    ManningOpen = 1,
    ManningClosedFs = 2,
    ManningClosedFull = 3,
    ForceMainHw = 4,
    ForceMainDw = 5,
    NCategories = 6,
}

/// One virtual-junction pair (DWSolver::VJuncPair).
#[derive(Debug, Clone, Copy)]
pub struct VJuncPair {
    pub node: i32,
    pub link_a: i32,
    pub link_b: i32,
    pub up_link: i32,
    pub dn_link: i32,
    pub through: u8,
    pub lambda: f64,
    pub sigma_j: f64,
    pub a_up_mid: f64,
    pub r_up_mid: f64,
    pub dq4j: f64,
    pub active: u8,
    pub dirn: i8,
    pub resid_max: f64,
    pub resid_sum: f64,
    pub resid_n: i64,
}

impl Default for VJuncPair {
    fn default() -> Self {
        VJuncPair {
            node: -1,
            link_a: -1,
            link_b: -1,
            up_link: -1,
            dn_link: -1,
            through: 0,
            lambda: 0.0,
            sigma_j: 1.0,
            a_up_mid: 0.0,
            r_up_mid: 0.0,
            dq4j: 0.0,
            active: 0,
            dirn: 0,
            resid_max: 0.0,
            resid_sum: 0.0,
            resid_n: 0,
        }
    }
}

/// Per-node invariant tile (DWSolver::NodeTile; built by initNodeStates, not at init).
#[derive(Debug, Clone, Copy, Default)]
pub struct NodeTile {
    pub full_depth: f64,
    pub y_crown: f64,
    pub invert_elev: f64,
    pub ponded_area: f64,
    pub sur_depth: f64,
    pub full_volume: f64,
    pub rpt_full_volume: f64,
    pub degree: i32,
    pub is_storage: u8,
    pub is_outfall: u8,
    pub is_virtual: u8,
}

/// Fraction of y_full at which a fed virtual junction's wetting floor is evaluated.
pub const K_VJ_WET_SEED_FRAC: f64 = 0.02;

/// The dynamic-wave solver.
#[derive(Debug, Clone)]
pub struct DWSolver {
    // --- public configuration
    pub head_tol: f64,
    pub max_trials: i32,
    pub omega: f64,
    pub surcharge_method: SurchargeMethod,
    pub node_continuity: NodeContinuity,
    pub anderson_accel: bool,
    pub unsteady_friction: i32,
    pub uf_k3: f64,
    pub tpa_celerity: f64,
    pub uf_nb_up_: Vec<i32>,
    pub uf_nb_dn_: Vec<i32>,
    pub uf_sg_up_: Vec<i8>,
    pub uf_sg_dn_: Vec<i8>,
    pub evap_rate: f64,

    // --- private state
    pub n_nodes_: i32,
    pub n_links_: i32,
    pub n_conduits_: i32,
    pub num_threads_: i32,
    /// upstream `groups_ != nullptr` (see the file notes)
    pub groups_attached_: bool,
    pub any_conduit_seep_: bool,
    pub losses_all_zero_: bool,
    pub conduit_idx_: Vec<i32>,
    pub is_open_: Vec<u8>,
    pub is_force_main_: Vec<u8>,
    pub has_losses_: Vec<u8>,
    pub barrels_d_: Vec<f64>,
    pub cached_length_: Vec<f64>,
    pub inv_length_: Vec<f64>,

    pub tile_uj_: Vec<i32>,
    pub tile_crow_: Vec<i32>,
    pub tile_n1_: Vec<i32>,
    pub tile_n2_: Vec<i32>,
    pub tile_inv1_elev_: Vec<f64>,
    pub tile_inv2_elev_: Vec<f64>,
    pub tile_z1_off_: Vec<f64>,
    pub tile_z2_off_: Vec<f64>,
    pub tile_y_full_: Vec<f64>,
    pub tile_a_full_: Vec<f64>,
    pub tile_r_full_: Vec<f64>,
    pub tile_w_max_: Vec<f64>,
    pub tile_length_: Vec<f64>,
    pub tile_inv_length_: Vec<f64>,
    pub tile_links_length_: Vec<f64>,
    pub tile_beta_: Vec<f64>,
    pub tile_q_max_: Vec<f64>,
    pub tile_rough_factor_: Vec<f64>,
    pub tile_barrels_d_: Vec<f64>,
    pub tile_is_open_: Vec<u8>,
    pub tile_is_force_main_: Vec<u8>,
    pub tile_is_closed_: Vec<u8>,
    pub tile_has_losses_: Vec<u8>,
    pub tile_xsect_batch_shape_: Vec<i32>,
    pub tile_shape_: Vec<XsectShape>,
    pub tile_has_offset_: Vec<u8>,
    pub tile_uj_to_ci_: Vec<i32>,
    pub tile_culvert_code_: Vec<i32>,
    pub tile_slope_: Vec<f64>,
    pub tile_q_limit_: Vec<f64>,
    pub tile_loss_inlet_: Vec<f64>,
    pub tile_loss_outlet_: Vec<f64>,
    pub tile_loss_avg_: Vec<f64>,
    pub tile_roughness_: Vec<f64>,
    pub tile_fm_sbot_: Vec<f64>,
    pub tile_fm_rbot_: Vec<f64>,
    pub tile_has_flap_gate_: Vec<u8>,
    pub tile_direction_: Vec<i8>,

    pub dt_gravity_: f64,
    pub last_converged_: bool,
    pub min_surf_area_: f64,
    pub wcap_d1_: Vec<f64>,
    pub wcap_d2_: Vec<f64>,
    pub wcap_dm_: Vec<f64>,
    pub variable_step_: f64,
    pub xnode_: DWNodeArrays,

    pub area1_: Vec<f64>,
    pub area2_: Vec<f64>,
    pub area_mid_: Vec<f64>,
    pub hrad_mid_: Vec<f64>,
    pub width_mid_: Vec<f64>,
    pub depth1_: Vec<f64>,
    pub depth2_: Vec<f64>,
    pub depth_mid_: Vec<f64>,

    pub velocity_: Vec<f64>,
    pub froude_: Vec<f64>,
    pub sigma_: Vec<f64>,
    pub dqdh_: Vec<f64>,
    pub new_flow_: Vec<f64>,
    pub q1_: Vec<f64>,
    pub area_old_: Vec<f64>,
    pub bypassed_: Vec<u8>,

    pub csr_row_: Vec<i32>,
    pub csr_link_: Vec<i32>,
    pub csr_is_n2_: Vec<u8>,
    pub csr_other_outfall_: Vec<u8>,

    pub vjunc_: Vec<VJuncPair>,
    pub vj_pair_n1_: Vec<i32>,
    pub vj_pair_n2_: Vec<i32>,
    pub vj_wet_floor_: Vec<f64>,

    pub node_tile_: Vec<NodeTile>,
    pub unit_sys_: i32,

    pub surf_area1_: Vec<f64>,
    pub surf_area2_: Vec<f64>,
    pub hrad1_: Vec<f64>,
    pub width1_: Vec<f64>,
    pub width2_: Vec<f64>,
    pub h1_: Vec<f64>,
    pub h2_: Vec<f64>,
    pub fasnh_: Vec<f64>,

    pub aa_y_prev_: Vec<f64>,
    pub aa_g_prev_: Vec<f64>,
    pub aa_r_prev_: Vec<f64>,
    pub aa_skip_: Vec<u8>,

    pub category_: Vec<MomentumCategory>,

    pub dps_config_: DPSConfig,
    pub dps_: DPSLinkArrays,
    pub sim_time_: f64,

    pub tpa_w_: Vec<f64>,
    pub tpa_latch_: Vec<u8>,
    pub tpa_latch_changed_: Vec<u8>,
}

impl Default for DWSolver {
    fn default() -> Self {
        DWSolver {
            head_tol: constants::DEFAULT_HEAD_TOL,
            max_trials: constants::DEFAULT_MAX_TRIALS,
            omega: constants::OMEGA,
            surcharge_method: SurchargeMethod::Extran,
            node_continuity: NodeContinuity::Explicit,
            anderson_accel: false,
            unsteady_friction: 0,
            uf_k3: 0.015,
            tpa_celerity: 100.0,
            uf_nb_up_: Vec::new(),
            uf_nb_dn_: Vec::new(),
            uf_sg_up_: Vec::new(),
            uf_sg_dn_: Vec::new(),
            evap_rate: 0.0,
            n_nodes_: 0,
            n_links_: 0,
            n_conduits_: 0,
            num_threads_: 1,
            groups_attached_: false,
            any_conduit_seep_: false,
            losses_all_zero_: true,
            conduit_idx_: Vec::new(),
            is_open_: Vec::new(),
            is_force_main_: Vec::new(),
            has_losses_: Vec::new(),
            barrels_d_: Vec::new(),
            cached_length_: Vec::new(),
            inv_length_: Vec::new(),
            tile_uj_: Vec::new(),
            tile_crow_: Vec::new(),
            tile_n1_: Vec::new(),
            tile_n2_: Vec::new(),
            tile_inv1_elev_: Vec::new(),
            tile_inv2_elev_: Vec::new(),
            tile_z1_off_: Vec::new(),
            tile_z2_off_: Vec::new(),
            tile_y_full_: Vec::new(),
            tile_a_full_: Vec::new(),
            tile_r_full_: Vec::new(),
            tile_w_max_: Vec::new(),
            tile_length_: Vec::new(),
            tile_inv_length_: Vec::new(),
            tile_links_length_: Vec::new(),
            tile_beta_: Vec::new(),
            tile_q_max_: Vec::new(),
            tile_rough_factor_: Vec::new(),
            tile_barrels_d_: Vec::new(),
            tile_is_open_: Vec::new(),
            tile_is_force_main_: Vec::new(),
            tile_is_closed_: Vec::new(),
            tile_has_losses_: Vec::new(),
            tile_xsect_batch_shape_: Vec::new(),
            tile_shape_: Vec::new(),
            tile_has_offset_: Vec::new(),
            tile_uj_to_ci_: Vec::new(),
            tile_culvert_code_: Vec::new(),
            tile_slope_: Vec::new(),
            tile_q_limit_: Vec::new(),
            tile_loss_inlet_: Vec::new(),
            tile_loss_outlet_: Vec::new(),
            tile_loss_avg_: Vec::new(),
            tile_roughness_: Vec::new(),
            tile_fm_sbot_: Vec::new(),
            tile_fm_rbot_: Vec::new(),
            tile_has_flap_gate_: Vec::new(),
            tile_direction_: Vec::new(),
            dt_gravity_: 0.0,
            last_converged_: false,
            min_surf_area_: constants::MIN_SURFAREA,
            wcap_d1_: Vec::new(),
            wcap_d2_: Vec::new(),
            wcap_dm_: Vec::new(),
            variable_step_: 0.0,
            xnode_: DWNodeArrays::default(),
            area1_: Vec::new(),
            area2_: Vec::new(),
            area_mid_: Vec::new(),
            hrad_mid_: Vec::new(),
            width_mid_: Vec::new(),
            depth1_: Vec::new(),
            depth2_: Vec::new(),
            depth_mid_: Vec::new(),
            velocity_: Vec::new(),
            froude_: Vec::new(),
            sigma_: Vec::new(),
            dqdh_: Vec::new(),
            new_flow_: Vec::new(),
            q1_: Vec::new(),
            area_old_: Vec::new(),
            bypassed_: Vec::new(),
            csr_row_: Vec::new(),
            csr_link_: Vec::new(),
            csr_is_n2_: Vec::new(),
            csr_other_outfall_: Vec::new(),
            vjunc_: Vec::new(),
            vj_pair_n1_: Vec::new(),
            vj_pair_n2_: Vec::new(),
            vj_wet_floor_: Vec::new(),
            node_tile_: Vec::new(),
            unit_sys_: 0,
            surf_area1_: Vec::new(),
            surf_area2_: Vec::new(),
            hrad1_: Vec::new(),
            width1_: Vec::new(),
            width2_: Vec::new(),
            h1_: Vec::new(),
            h2_: Vec::new(),
            fasnh_: Vec::new(),
            aa_y_prev_: Vec::new(),
            aa_g_prev_: Vec::new(),
            aa_r_prev_: Vec::new(),
            aa_skip_: Vec::new(),
            category_: Vec::new(),
            dps_config_: DPSConfig::default(),
            dps_: DPSLinkArrays::default(),
            sim_time_: 0.0,
            tpa_w_: Vec::new(),
            tpa_latch_: Vec::new(),
            tpa_latch_changed_: Vec::new(),
        }
    }
}

/// Legacy open-profile shapes (xsect_isOpen), in the link-store shape numbering.
fn is_open_shape(shape: XsectShape) -> bool {
    matches!(
        shape,
        XsectShape::Dummy
            | XsectShape::RectOpen
            | XsectShape::Trapezoidal
            | XsectShape::Triangular
            | XsectShape::Parabolic
            | XsectShape::Power
            | XsectShape::Irregular
            | XsectShape::StreetXsect
    )
}

/// DynamicWave.cpp's file-local buildXSP: the scalar section view of one link, with yw_max
/// and, for tabulated shapes, the link's transect tables attached.
pub(crate) fn build_xsp(ctx: &SimulationContext, uk: usize) -> XSectParams<'_> {
    let links = &ctx.links;
    let mut xs = XSectParams::default();
    let ls = links.xsect_shape[uk];
    xs.type_ = link::translate_shape(ls);
    xs.y_full = links.xsect_y_full[uk];
    xs.a_full = links.xsect_a_full[uk];
    xs.w_max = links.xsect_w_max[uk];
    xs.yw_max = links.xsect_yw_max[uk];
    xs.r_full = links.xsect_r_full[uk];
    xs.s_full = links.xsect_s_full[uk];
    xs.s_max = links.xsect_s_max[uk];
    xs.y_bot = links.xsect_y_bot[uk];
    xs.a_bot = links.xsect_a_bot[uk];
    xs.s_bot = links.xsect_s_bot[uk];
    xs.r_bot = links.xsect_r_bot[uk];
    if ls == XsectShape::Irregular || ls == XsectShape::Custom || ls == XsectShape::StreetXsect {
        let ci = links.xsect_curve[uk];
        if ci >= 0 && (ci as usize) < ctx.transect_tables.len() {
            let td = &ctx.transect_tables[ci as usize];
            xs.transect = ci;
            xs.area_tbl = Some(&td.area_tbl[..]);
            xs.hrad_tbl = Some(&td.hrad_tbl[..]);
            xs.width_tbl = Some(&td.width_tbl[..]);
            xs.area_lut = Some(&td.area_lut);
            xs.transect_tbl_size = N_TRANSECT_TBL;
        }
    }
    xs
}

fn assign<T: Clone>(v: &mut Vec<T>, n: usize, d: T) {
    v.clear();
    v.resize(n, d);
}

impl DWSolver {
    pub fn init(&mut self, n_nodes: i32, n_links: i32, _groups: &XSectGroups, ctx: &SimulationContext) {
        self.n_nodes_ = n_nodes;
        self.n_links_ = n_links;
        self.groups_attached_ = true;

        // MIN_SURFAREA floor: a user value is in project area units.
        {
            let ucf_len = ucf::UCF[ucf::LENGTH as usize][ucf::get_unit_system(ctx.options.flow_units as i32) as usize];
            self.min_surf_area_ = if ctx.options.min_surf_area > 0.0 { ctx.options.min_surf_area / (ucf_len * ucf_len) } else { constants::MIN_SURFAREA };
        }

        let un = n_nodes as usize;
        let ul = n_links as usize;

        // Unsteady-friction cross-link stencil (issue #156).
        if self.unsteady_friction != 0 {
            assign(&mut self.uf_nb_up_, ul, -1);
            assign(&mut self.uf_nb_dn_, ul, -1);
            assign(&mut self.uf_sg_up_, ul, 0);
            assign(&mut self.uf_sg_dn_, ul, 0);
            let mut ccount = vec![0i32; un];
            let mut other = vec![0i32; un];
            let mut inc = vec![[-1i32, -1i32]; un];
            for j in 0..n_links {
                let ujj = j as usize;
                let is_cond = ctx.links.type_[ujj] == LinkType::Conduit;
                for nd in [ctx.links.node1[ujj], ctx.links.node2[ujj]] {
                    if nd < 0 || nd >= n_nodes {
                        continue;
                    }
                    let und = nd as usize;
                    if !is_cond {
                        other[und] += 1;
                        continue;
                    }
                    if ccount[und] < 2 {
                        inc[und][ccount[und] as usize] = j;
                    }
                    ccount[und] += 1;
                }
            }
            for j in 0..n_links {
                let ujj = j as usize;
                if ctx.links.type_[ujj] != LinkType::Conduit {
                    continue;
                }
                let n1 = ctx.links.node1[ujj];
                let n2 = ctx.links.node2[ujj];
                if n1 >= 0 && n1 < n_nodes {
                    let u1 = n1 as usize;
                    if ccount[u1] == 2 && other[u1] == 0 {
                        let k = if inc[u1][0] == j { inc[u1][1] } else { inc[u1][0] };
                        if k >= 0 {
                            self.uf_nb_up_[ujj] = k;
                            self.uf_sg_up_[ujj] = if ctx.links.node2[k as usize] == n1 { 1 } else { -1 };
                        }
                    }
                }
                if n2 >= 0 && n2 < n_nodes {
                    let u2 = n2 as usize;
                    if ccount[u2] == 2 && other[u2] == 0 {
                        let k = if inc[u2][0] == j { inc[u2][1] } else { inc[u2][0] };
                        if k >= 0 {
                            self.uf_nb_dn_[ujj] = k;
                            self.uf_sg_dn_[ujj] = if ctx.links.node1[k as usize] == n2 { 1 } else { -1 };
                        }
                    }
                }
            }
        }

        // Seepage possibility is static input.
        self.any_conduit_seep_ = false;
        self.losses_all_zero_ = true;
        {
            let cd = &ctx.link_subtypes.conduits;
            for r in 0..cd.count() {
                if cd.seep_rate[r as usize] > 0.0 {
                    self.any_conduit_seep_ = true;
                    break;
                }
            }
        }

        self.xnode_.resize(un);

        self.area1_.resize(ul, 0.0);
        self.area2_.resize(ul, 0.0);
        self.area_mid_.resize(ul, 0.0);
        self.hrad_mid_.resize(ul, 0.0);
        self.width_mid_.resize(ul, 0.0);
        self.depth1_.resize(ul, 0.0);
        self.depth2_.resize(ul, 0.0);
        self.depth_mid_.resize(ul, 0.0);

        self.velocity_.resize(ul, 0.0);
        self.froude_.resize(ul, 0.0);
        self.sigma_.resize(ul, 0.0);
        self.dqdh_.resize(ul, 0.0);
        self.new_flow_.resize(ul, 0.0);
        assign(&mut self.q1_, ul, 0.0);
        self.area_old_.resize(ul, 0.0);
        assign(&mut self.bypassed_, ul, 0);
        self.surf_area1_.resize(ul, 0.0);
        self.surf_area2_.resize(ul, 0.0);
        self.hrad1_.resize(ul, 0.0);
        self.width1_.resize(ul, 0.0);
        self.width2_.resize(ul, 0.0);
        self.h1_.resize(ul, 0.0);
        self.h2_.resize(ul, 0.0);
        self.fasnh_.resize(ul, 1.0);
        self.wcap_d1_.resize(ul, 0.0);
        self.wcap_d2_.resize(ul, 0.0);
        self.wcap_dm_.resize(ul, 0.0);

        // Conduit index list; DUMMY sections excluded (legacy isTrueConduit).
        self.conduit_idx_.clear();
        self.conduit_idx_.reserve(ul);
        for j in 0..n_links {
            let ujc = j as usize;
            if ctx.links.type_[ujc] == LinkType::Conduit && ctx.links.xsect_shape[ujc] != XsectShape::Dummy {
                self.conduit_idx_.push(j);
            }
        }
        self.n_conduits_ = self.conduit_idx_.len() as i32;

        // Per-link invariants
        self.is_open_.resize(ul, 0);
        self.is_force_main_.resize(ul, 0);
        self.has_losses_.resize(ul, 0);
        self.barrels_d_.resize(ul, 1.0);
        self.cached_length_.resize(ul, 0.0);
        self.inv_length_.resize(ul, 0.0);

        let links = &ctx.links;
        let cd = &ctx.link_subtypes.conduits;
        for ci in 0..self.n_conduits_ {
            let j = self.conduit_idx_[ci as usize];
            let uj = j as usize;
            let ucr = ctx.link_subtypes.conduit_row(j) as usize;
            let shape = links.xsect_shape[uj];

            self.is_open_[uj] = is_open_shape(shape) as u8;
            self.is_force_main_[uj] = (shape == XsectShape::ForceMain) as u8;
            self.has_losses_[uj] = (cd.loss_inlet[ucr] != 0.0 || cd.loss_outlet[ucr] != 0.0 || cd.loss_avg[ucr] != 0.0) as u8;
            self.barrels_d_[uj] = cmax(cd.barrels[ucr], 1) as f64;
            // legacy initLinks: Conduit.q1 = newFlow / barrels.
            self.q1_[uj] = links.flow[uj] / self.barrels_d_[uj];
            let mut len = cd.mod_length[ucr];
            if len <= 0.0 {
                len = cd.length[ucr];
            }
            self.cached_length_[uj] = len;
            self.inv_length_[uj] = if len > 0.0 { 1.0 / len } else { 0.0 };

            // Seed the conduit's initial mid-area from its initial flow depth.
            let xs0 = build_xsp(ctx, uj);
            self.area_mid_[uj] = xsect::get_a_of_y(&xs0, links.depth[uj]);
        }

        // Momentum category arrays
        self.category_.resize(ul, MomentumCategory::SkipDry);

        self.refresh_conduit_tile(ctx);
        self.build_conduit_node_csr(ctx);
        self.build_virtual_junction_pairs(ctx);

        // Anderson acceleration state arrays
        self.aa_y_prev_.resize(un, 0.0);
        self.aa_g_prev_.resize(un, 0.0);
        self.aa_r_prev_.resize(un, 0.0);
        self.aa_skip_.resize(un, 0);

        // TPA initialization (issue #156 Phase 5)
        if self.surcharge_method == SurchargeMethod::Tpa {
            let ucf_len2 = ucf::UCF[ucf::LENGTH as usize][ucf::get_unit_system(ctx.options.flow_units as i32) as usize];
            let a_fts = cmax(self.tpa_celerity / ucf_len2, 1.0);
            let ucn = self.n_conduits_ as usize;
            assign(&mut self.tpa_w_, ucn, 0.0);
            assign(&mut self.tpa_latch_, ucn, 0);
            assign(&mut self.tpa_latch_changed_, ucn, 0);
            for ci in 0..self.n_conduits_ {
                let uj = self.conduit_idx_[ci as usize] as usize;
                let af = ctx.links.xsect_a_full[uj];
                self.tpa_w_[ci as usize] = if af > 0.0 { GRAVITY * af / (a_fts * a_fts) } else { 0.0 };
            }
        }

        // Dynamic Preissmann Slot initialization
        if self.surcharge_method == SurchargeMethod::DynamicSlot {
            self.dps_config_.c_pT = ctx.options.dps_target_celerity * 3.28084;
            self.dps_config_.alpha = cmax(ctx.options.dps_alpha, 2.0);
            self.dps_config_.r = cmax(ctx.options.dps_decay_time, 0.001);
            self.dps_config_.c_pT_sq = self.dps_config_.c_pT * self.dps_config_.c_pT;

            self.dps_.resize(self.n_conduits_ as usize);
            self.sim_time_ = 0.0;

            for ci in 0..self.n_conduits_ {
                let j = self.conduit_idx_[ci as usize];
                let uj = j as usize;
                let uci = ci as usize;

                if self.is_open_[uj] != 0 {
                    self.dps_.P_hat_0[uci] = 1.0;
                    self.dps_.P_hat[uci] = 1.0;
                    self.dps_.P[uci] = 1.0;
                    continue;
                }

                let af = links.xsect_a_full[uj];
                let tw = links.xsect_w_max[uj];
                let l_d = if tw > 0.0 { af / tw } else { 0.0 };
                let c_g = if l_d > 0.0 { (GRAVITY * l_d).sqrt() } else { 1.0 };

                self.dps_.P_hat_0[uci] = self.dps_config_.c_pT / (self.dps_config_.alpha * c_g);
                self.dps_.P_hat_0[uci] = cmax(self.dps_.P_hat_0[uci], 1.0);

                self.dps_.P_hat[uci] = self.dps_.P_hat_0[uci];
                self.dps_.P[uci] = self.dps_.P_hat[uci];
                self.dps_.As[uci] = 0.0;
                self.dps_.hs[uci] = 0.0;
                self.dps_.hs_prev_iter[uci] = 0.0;

                self.dps_.T_s_target[uci] = if self.dps_config_.c_pT_sq > 0.0 { (GRAVITY * af) / self.dps_config_.c_pT_sq } else { 0.0 };

                self.dps_.surcharged[uci] = 0;
                self.dps_.t_s[uci] = 0.0;
            }
        }
    }

    /// Gather the timestep-invariant per-conduit fields into the conduit-dense tile.
    pub fn refresh_conduit_tile(&mut self, ctx: &SimulationContext) {
        let nc = self.n_conduits_ as usize;
        self.tile_uj_.resize(nc, 0);
        self.tile_crow_.resize(nc, 0);
        self.tile_n1_.resize(nc, 0);
        self.tile_n2_.resize(nc, 0);
        self.tile_inv1_elev_.resize(nc, 0.0);
        self.tile_inv2_elev_.resize(nc, 0.0);
        self.tile_z1_off_.resize(nc, 0.0);
        self.tile_z2_off_.resize(nc, 0.0);
        self.tile_y_full_.resize(nc, 0.0);
        self.tile_a_full_.resize(nc, 0.0);
        self.tile_r_full_.resize(nc, 0.0);
        self.tile_w_max_.resize(nc, 0.0);
        self.tile_length_.resize(nc, 0.0);
        self.tile_inv_length_.resize(nc, 0.0);
        self.tile_links_length_.resize(nc, 0.0);
        self.tile_beta_.resize(nc, 0.0);
        self.tile_q_max_.resize(nc, 0.0);
        self.tile_rough_factor_.resize(nc, 0.0);
        self.tile_barrels_d_.resize(nc, 0.0);
        self.tile_is_open_.resize(nc, 0);
        self.tile_is_force_main_.resize(nc, 0);
        self.tile_is_closed_.resize(nc, 0);
        self.tile_has_losses_.resize(nc, 0);
        self.tile_xsect_batch_shape_.resize(nc, 0);
        // value-initialised XsectShape is 0, CIRCULAR
        self.tile_shape_.resize(nc, XsectShape::Circular);
        self.tile_has_offset_.resize(nc, 0);
        self.tile_culvert_code_.resize(nc, 0);
        self.tile_slope_.resize(nc, 0.0);
        self.tile_q_limit_.resize(nc, 0.0);
        self.tile_loss_inlet_.resize(nc, 0.0);
        self.tile_loss_outlet_.resize(nc, 0.0);
        self.tile_loss_avg_.resize(nc, 0.0);
        self.tile_roughness_.resize(nc, 0.0);
        self.tile_fm_sbot_.resize(nc, 0.0);
        self.tile_fm_rbot_.resize(nc, 0.0);
        self.tile_has_flap_gate_.resize(nc, 0);
        self.tile_direction_.resize(nc, 0);
        assign(&mut self.tile_uj_to_ci_, self.n_links_ as usize, -1);

        let links = &ctx.links;
        let nodes = &ctx.nodes;
        let cd = &ctx.link_subtypes.conduits;
        for ci in 0..self.n_conduits_ {
            let uci = ci as usize;
            let j = self.conduit_idx_[uci];
            let uj = j as usize;
            let ucr = ctx.link_subtypes.conduit_row(j) as usize;
            let n1 = links.node1[uj];
            let n2 = links.node2[uj];
            let un1 = n1 as usize;
            let un2 = n2 as usize;

            self.tile_uj_[uci] = j;
            self.tile_crow_[uci] = ucr as i32;
            self.tile_n1_[uci] = n1;
            self.tile_n2_[uci] = n2;
            self.tile_inv1_elev_[uci] = nodes.invert_elev[un1];
            self.tile_inv2_elev_[uci] = nodes.invert_elev[un2];
            self.tile_z1_off_[uci] = links.offset1[uj];
            self.tile_z2_off_[uci] = links.offset2[uj];
            self.tile_y_full_[uci] = links.xsect_y_full[uj];
            self.tile_a_full_[uci] = links.xsect_a_full[uj];
            self.tile_r_full_[uci] = links.xsect_r_full[uj];
            self.tile_w_max_[uci] = links.xsect_w_max[uj];
            self.tile_length_[uci] = self.cached_length_[uj];
            self.tile_inv_length_[uci] = self.inv_length_[uj];
            self.tile_links_length_[uci] = if cd.true_length[ucr] > 0.0 { cd.true_length[ucr] } else { cd.length[ucr] };
            self.tile_beta_[uci] = cd.beta[ucr];
            self.tile_q_max_[uci] = cd.q_max[ucr];
            self.tile_rough_factor_[uci] = cd.rough_factor[ucr];
            self.tile_barrels_d_[uci] = self.barrels_d_[uj];
            self.tile_is_open_[uci] = self.is_open_[uj];
            self.tile_is_force_main_[uci] = self.is_force_main_[uj];
            self.tile_is_closed_[uci] = if links.is_closed[uj] != 0 { 1 } else { 0 };
            self.tile_has_losses_[uci] = self.has_losses_[uj];
            self.tile_xsect_batch_shape_[uci] = links.xsect_batch_shape[uj];
            self.tile_shape_[uci] = links.xsect_shape[uj];
            self.tile_has_offset_[uci] = if links.offset1[uj] > 0.0 || links.offset2[uj] > 0.0 { 1 } else { 0 };
            self.tile_culvert_code_[uci] = cd.culvert_code[ucr];
            self.tile_slope_[uci] = cd.slope[ucr];
            self.tile_q_limit_[uci] = links.q_limit[uj];
            self.tile_loss_inlet_[uci] = cd.loss_inlet[ucr];
            self.tile_loss_outlet_[uci] = cd.loss_outlet[ucr];
            self.tile_loss_avg_[uci] = cd.loss_avg[ucr];
            self.tile_roughness_[uci] = cd.roughness[ucr];
            self.tile_fm_sbot_[uci] = links.xsect_s_bot[uj];
            self.tile_fm_rbot_[uci] = links.xsect_r_bot[uj];
            self.tile_has_flap_gate_[uci] = if links.has_flap_gate[uj] != 0 { 1 } else { 0 };
            self.tile_direction_[uci] = links.direction[uj] as i8;
            self.tile_uj_to_ci_[uj] = ci;
        }
    }

    /// Node-to-incident-conduit CSR, entries in ascending conduit link order.
    pub fn build_conduit_node_csr(&mut self, ctx: &SimulationContext) {
        let links = &ctx.links;
        let nodes = &ctx.nodes;
        let un = self.n_nodes_ as usize;

        assign(&mut self.csr_row_, un + 1, 0);
        for ci in 0..self.n_conduits_ {
            let uj = self.conduit_idx_[ci as usize] as usize;
            let n1 = links.node1[uj];
            let n2 = links.node2[uj];
            if n1 >= 0 && n1 < self.n_nodes_ {
                self.csr_row_[n1 as usize + 1] += 1;
            }
            if n2 >= 0 && n2 < self.n_nodes_ {
                self.csr_row_[n2 as usize + 1] += 1;
            }
        }
        for i in 0..un {
            self.csr_row_[i + 1] += self.csr_row_[i];
        }

        let n_entries = self.csr_row_[un] as usize;
        assign(&mut self.csr_link_, n_entries, 0);
        assign(&mut self.csr_is_n2_, n_entries, 0);
        assign(&mut self.csr_other_outfall_, n_entries, 0);

        let mut cursor: Vec<i32> = self.csr_row_[..un].to_vec();
        for ci in 0..self.n_conduits_ {
            let j = self.conduit_idx_[ci as usize];
            let uj = j as usize;
            let n1 = links.node1[uj];
            let n2 = links.node2[uj];
            let of1 = n1 >= 0 && n1 < self.n_nodes_ && nodes.type_[n1 as usize] == NodeType::Outfall;
            let of2 = n2 >= 0 && n2 < self.n_nodes_ && nodes.type_[n2 as usize] == NodeType::Outfall;
            if n1 >= 0 && n1 < self.n_nodes_ {
                let k = cursor[n1 as usize] as usize;
                cursor[n1 as usize] += 1;
                self.csr_link_[k] = j;
                self.csr_is_n2_[k] = 0;
                self.csr_other_outfall_[k] = if of2 { 1 } else { 0 };
            }
            if n2 >= 0 && n2 < self.n_nodes_ {
                let k = cursor[n2 as usize] as usize;
                cursor[n2 as usize] += 1;
                self.csr_link_[k] = j;
                self.csr_is_n2_[k] = 1;
                self.csr_other_outfall_[k] = if of1 { 1 } else { 0 };
            }
        }
    }

    /// The virtual-junction pair table (empty for a model without virtual junctions).
    pub fn build_virtual_junction_pairs(&mut self, ctx: &SimulationContext) {
        let links = &ctx.links;
        let nodes = &ctx.nodes;
        let cd = &ctx.link_subtypes.conduits;

        self.vjunc_.clear();
        assign(&mut self.vj_pair_n1_, self.n_links_ as usize, -1);
        assign(&mut self.vj_pair_n2_, self.n_links_ as usize, -1);
        assign(&mut self.vj_wet_floor_, self.n_nodes_ as usize, 0.0);

        let mut row_of_node = vec![-1i32; self.n_nodes_ as usize];
        for i in 0..self.n_nodes_ {
            if (i as usize) < nodes.is_virtual.len() && nodes.is_virtual[i as usize] != 0 {
                row_of_node[i as usize] = self.vjunc_.len() as i32;
                let p = VJuncPair { node: i, ..VJuncPair::default() };
                self.vjunc_.push(p);
            }
        }
        if self.vjunc_.is_empty() {
            return;
        }

        for j in 0..self.n_links_ {
            let uj = j as usize;
            if links.type_[uj] != LinkType::Conduit {
                continue;
            }
            for at_node1 in [true, false] {
                let n = if at_node1 { links.node1[uj] } else { links.node2[uj] };
                if n < 0 || n >= self.n_nodes_ {
                    continue;
                }
                let r = row_of_node[n as usize];
                if r < 0 {
                    continue;
                }
                let p = &mut self.vjunc_[r as usize];
                if p.link_a < 0 {
                    p.link_a = j;
                } else if p.link_b < 0 {
                    p.link_b = j;
                }
                if at_node1 {
                    self.vj_pair_n1_[uj] = r;
                } else {
                    self.vj_pair_n2_[uj] = r;
                }
            }
        }

        for pi in 0..self.vjunc_.len() {
            let mut p = self.vjunc_[pi];
            if p.link_a < 0 || p.link_b < 0 {
                continue; // validation rejects these
            }
            let ua = p.link_a as usize;
            let ub = p.link_b as usize;

            if links.node2[ua] == p.node && links.node1[ub] == p.node {
                p.up_link = p.link_a;
                p.dn_link = p.link_b;
                p.through = 1;
            } else if links.node2[ub] == p.node && links.node1[ua] == p.node {
                p.up_link = p.link_b;
                p.dn_link = p.link_a;
                p.through = 1;
            } else {
                p.through = 0; // sag (both node2) or peak (both node1)
            }

            let ca = ctx.link_subtypes.conduit_row(p.link_a);
            let cb = ctx.link_subtypes.conduit_row(p.link_b);
            let la = if ca >= 0 { cd.length[ca as usize] } else { 0.0 };
            let lb = if cb >= 0 { cd.length[cb as usize] } else { 0.0 };
            p.lambda = 0.5 * (la + lb);

            {
                let xs = build_xsp(ctx, ua);
                let y_seed = K_VJ_WET_SEED_FRAC * links.xsect_y_full[ua];
                let w_seed = if y_seed > 0.0 { xsect::get_w_of_y(&xs, y_seed) } else { 0.0 };
                let barrels = if ca >= 0 { cd.barrels[ca as usize] } else { 1 };
                self.vj_wet_floor_[p.node as usize] = w_seed * p.lambda * barrels as f64;
            }
            self.vjunc_[pi] = p;
        }
    }

    /// Configure the dynamic-wave worker count (SWMM_DW_THREADS overrides).
    pub fn set_num_threads(&mut self, n: i32, warnings: Option<&mut Vec<Bytes>>) {
        if let Some(e) = std::env::var_os("SWMM_DW_THREADS") {
            let forced = crate::core::cstdlib::atoi(e.as_encoded_bytes());
            if forced >= 1 {
                self.num_threads_ = forced;
                if let Some(w) = warnings {
                    w.push(
                        format!(
                            "Dynamic-wave thread count forced to {} by the SWMM_DW_THREADS environment variable ([OPTIONS] THREADS ignored for dynamic wave).",
                            forced
                        )
                        .into_bytes(),
                    );
                }
                return;
            }
        }
        self.num_threads_ = thread_info::dw_threads(n, self.n_conduits_, warnings);
    }

    pub fn last_converged(&self) -> bool {
        self.last_converged_
    }

    pub fn is_initialized(&self) -> bool {
        !self.xnode_.sumdqdh.is_empty()
    }
}
