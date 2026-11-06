// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/fv/NetworkMeshData.hpp,
// FvClosureKernels.hpp (FvClosure) and ExplicitFvSolver.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026
// Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The explicit finite-volume network router (FLOW_ROUTING FV): the mesh, state and
// section-closure data the Router owns. The closure kernels and the exact section
// geometry are in `kernels`, the mesh builder in `mesh_builder`, the CPU backend's
// initialization half in `explicit_solver`.
//
// NOT YET TRANSLATED (phase 1e): the implicit and pressurized solvers (ImplicitNetworkSolver,
// ImplicitReachSolver, PressurizedHeadSolver, NativeGmres, NativeBlockIlu,
// NativeConduitCondensation, NativeLinearWorkspace, NonlinearPressureStorage), the
// checked advance (CheckedNetworkAdvance.hpp), INetworkSolver's step interface, and
// NetworkSolverFactory's GPU-plugin search: a 1D reference build has no plugin library,
// so makeNetworkSolver always returns the CPU ExplicitFvSolver (see Router::init_fv).
//
// Differences in form, never in result:
//   * FvGeometry's `XSectParams xs` borrows transect tables upstream; here it is
//     XSectParams<'static> holding only the scalars, with the table index in xs.transect.
//   * FvClosure's inline C arrays are fixed-size Rust arrays (boxed in FvGeometry so the
//     geometry vector stays cheap to move).

pub mod explicit_flux;
pub mod explicit_nodes;
pub mod explicit_solver;
pub mod explicit_step;
pub mod kernels;
pub mod mesh_builder;
pub mod native_gmres;

pub use explicit_solver::ExplicitFvSolver;

use crate::hydraulics::hyd_closure_kernels::CulvertCurve;
use crate::hydraulics::xsect_batch::XSectParams;

pub const K_NODE_VOL_SAMPLES: i32 = 129;
pub const K_CLOSURE_PANELS: usize = 128;
pub const K_CLOSURE_TABULATED: u8 = 0;
pub const K_CLOSURE_POLYNOMIAL: u8 = 1;

/// One section's closure table (fixed-size POD upstream).
#[allow(non_snake_case)]
#[derive(Debug, Clone)]
pub struct FvClosure {
    pub y_full: f64,
    pub dh: f64,
    pub inv_dh: f64,
    pub a_crown: f64,
    pub t_slot: f64,
    pub i1_crown: f64,
    pub r_full: f64,
    pub inv_a_crown_n: f64,
    pub c2: f64,
    pub c1: f64,
    pub p0: f64,
    pub p1: f64,
    pub kind: u8,
    pub A: [f64; K_CLOSURE_PANELS + 1],
    pub M: [f64; K_CLOSURE_PANELS + 1],
    pub I1: [f64; K_CLOSURE_PANELS + 1],
    pub R: [f64; K_CLOSURE_PANELS + 1],
    pub j_of_a: [i16; K_CLOSURE_PANELS + 1],
}

impl Default for FvClosure {
    fn default() -> Self {
        FvClosure {
            y_full: 0.0,
            dh: 0.0,
            inv_dh: 0.0,
            a_crown: 0.0,
            t_slot: 0.0,
            i1_crown: 0.0,
            r_full: 0.0,
            inv_a_crown_n: 0.0,
            c2: 0.0,
            c1: 0.0,
            p0: 0.0,
            p1: 0.0,
            kind: K_CLOSURE_TABULATED,
            A: [0.0; K_CLOSURE_PANELS + 1],
            M: [0.0; K_CLOSURE_PANELS + 1],
            I1: [0.0; K_CLOSURE_PANELS + 1],
            R: [0.0; K_CLOSURE_PANELS + 1],
            j_of_a: [0; K_CLOSURE_PANELS + 1],
        }
    }
}

/// One cross-section block shared by the cells of a conduit.
#[derive(Debug, Clone)]
pub struct FvGeometry {
    pub xs: XSectParams<'static>,
    pub y_full: f64,
    pub a_full: f64,
    pub w_max: f64,
    pub r_full: f64,
    pub barrel_scale: f64,
    pub y_crown: f64,
    pub t_slot: f64,
    pub a_crown: f64,
    pub i1_crown: f64,
    pub is_open: u8,
    pub barrels: i32,
    pub closure_tbl: Box<FvClosure>,
}

impl Default for FvGeometry {
    fn default() -> Self {
        FvGeometry {
            xs: XSectParams::default(),
            y_full: 0.0,
            a_full: 0.0,
            w_max: 0.0,
            r_full: 0.0,
            barrel_scale: 1.0,
            y_crown: 0.0,
            t_slot: 0.0,
            a_crown: 0.0,
            i1_crown: 0.0,
            is_open: 0,
            barrels: 1,
            closure_tbl: Box::default(),
        }
    }
}

/// The FV network mesh (NetworkMeshData).
#[derive(Debug, Clone, Default)]
pub struct NetworkMeshData {
    pub cell_geom: Vec<i32>,
    pub cell_conduit: Vec<i32>,
    pub cell_dx: Vec<f64>,
    pub cell_zb: Vec<f64>,
    pub cell_dzdx: Vec<f64>,
    pub cell_face0: Vec<i32>,
    pub cell_face1: Vec<i32>,
    pub cell_side0: Vec<i8>,
    pub cell_side1: Vec<i8>,
    pub face_cl: Vec<i32>,
    pub face_cr: Vec<i32>,
    pub face_node: Vec<i32>,
    pub face_geom: Vec<i32>,
    pub face_gate: Vec<u8>,
    pub face_culvert: Vec<i32>,
    pub face_zb: Vec<f64>,
    pub face_dx: Vec<f64>,
    pub face_dir_l: Vec<i8>,
    pub face_dir_r: Vec<i8>,
    pub face_virtual: Vec<u8>,
    pub face_vj_node: Vec<i32>,
    pub node_vj_face: Vec<i32>,
    pub conduit_cell_begin: Vec<i32>,
    pub conduit_cell_count: Vec<i32>,
    pub conduit_link: Vec<i32>,
    pub conduit_section: Vec<i32>,
    pub conduit_roughness: Vec<f64>,
    pub conduit_rough_factor: Vec<f64>,
    pub conduit_loss_inlet: Vec<f64>,
    pub conduit_loss_outlet: Vec<f64>,
    pub conduit_slope: Vec<f64>,
    pub conduit_culvert_code: Vec<i32>,
    pub conduit_culvert_curve: Vec<CulvertCurve>,
    pub conduit_culvert_mitered: Vec<u8>,
    pub chain_ptr: Vec<i32>,
    pub chain_cells: Vec<i32>,
    pub chain_dir: Vec<i8>,
    pub cell_chain: Vec<i32>,
    pub cell_chain_pos: Vec<i32>,
    pub struct_link: Vec<i32>,
    pub struct_n1: Vec<i32>,
    pub struct_n2: Vec<i32>,
    pub struct_is_dummy: Vec<u8>,
    pub node_dummy_drain: Vec<u8>,
    pub node_face_ptr: Vec<i32>,
    pub node_face_idx: Vec<i32>,
    pub node_face_sign: Vec<f64>,
    pub node_face_zb: Vec<f64>,
    pub node_invert: Vec<f64>,
    pub node_full_depth: Vec<f64>,
    pub node_ponded_area: Vec<f64>,
    pub node_sur_depth: Vec<f64>,
    pub node_can_pond: Vec<u8>,
    pub node_kind: Vec<u8>,
    pub node_area: Vec<f64>,
    pub node_vol_off: Vec<i32>,
    pub node_vol_dmax: Vec<f64>,
    pub node_vol_atop: Vec<f64>,
    pub node_vol_tbl: Vec<f64>,
    pub geom: Vec<FvGeometry>,
}

impl NetworkMeshData {
    pub fn n_cells(&self) -> i32 {
        self.cell_dx.len() as i32
    }
    pub fn n_faces(&self) -> i32 {
        self.face_cl.len() as i32
    }
    pub fn n_nodes(&self) -> i32 {
        self.node_invert.len() as i32
    }
    pub fn n_conduits(&self) -> i32 {
        self.conduit_cell_begin.len() as i32
    }
    pub fn n_chains(&self) -> i32 {
        if self.chain_ptr.is_empty() {
            0
        } else {
            self.chain_ptr.len() as i32 - 1
        }
    }

    /// Per-face section index derived from the adjacent cell (no-op when already sized).
    pub fn derive_face_geom(&mut self) {
        let nf = self.n_faces() as usize;
        if self.face_geom.len() == nf {
            return;
        }
        self.face_geom.clear();
        self.face_geom.resize(nf, 0);
        for f in 0..nf {
            let c = if self.face_cl[f] >= 0 { self.face_cl[f] } else { self.face_cr[f] };
            self.face_geom[f] = if c >= 0 && (c as usize) < self.cell_geom.len() { self.cell_geom[c as usize] } else { 0 };
        }
    }

    pub fn clear(&mut self) {
        *self = NetworkMeshData::default();
    }
}

pub const K_NODE_JUNCTION: u8 = 0;
pub const K_NODE_VIRTUAL: u8 = 1;
pub const K_NODE_STORAGE: u8 = 2;
pub const K_NODE_OUTFALL: u8 = 3;

/// The FV network state (NetworkStateData).
#[derive(Debug, Clone, Default)]
pub struct NetworkStateData {
    pub cell_a: Vec<f64>,
    pub cell_q: Vec<f64>,
    pub cell_h: Vec<f64>,
    pub node_volume: Vec<f64>,
    pub node_head: Vec<f64>,
    pub node_surf_area: Vec<f64>,
    pub node_overflow: Vec<f64>,
    pub node_overflow_vol: Vec<f64>,
    pub cell_phi: Vec<f64>,
    pub n_species: i32,
    pub cell_tpa: Vec<u8>,
}

impl NetworkStateData {
    pub fn resize(&mut self, n_cells: i32, n_nodes: i32, n_species_in: i32) {
        let nc = n_cells as usize;
        let nn = n_nodes as usize;
        fn a(v: &mut Vec<f64>, n: usize) {
            v.clear();
            v.resize(n, 0.0);
        }
        a(&mut self.cell_a, nc);
        a(&mut self.cell_q, nc);
        a(&mut self.cell_h, nc);
        a(&mut self.node_volume, nn);
        a(&mut self.node_head, nn);
        a(&mut self.node_surf_area, nn);
        a(&mut self.node_overflow, nn);
        a(&mut self.node_overflow_vol, nn);
        self.n_species = n_species_in;
        a(&mut self.cell_phi, nc * n_species_in as usize);
    }
}

