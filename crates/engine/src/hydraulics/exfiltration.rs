// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydraulics/Exfiltration.hpp and
// Exfiltration.cpp, with lidnode::bottomConductivity from src/engine/hydrology/LidNode.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Storage-node exfiltration: the per-node SoA and ExfilSolver::init.
//
// NOT YET TRANSLATED (run-time step function, with the routing step):
// ExfilSolver::computeAll.
//
// Differences in form, never in result:
//   * infil::grnampt_init is hydrology::infiltration::grnampt_init. lidnode::
//     bottomConductivity is translated here as a private helper, because LidNode's
//     run-time half is not in the crate yet; when it lands, call it instead (the
//     arithmetic is copied statement for statement).
//   * Members keep upstream's names and are `pub` (soa_ is private upstream) so the
//     gate's state dumper can read them.

use crate::core::context::SimulationContext;
use crate::core::unit_conversion as ucf;
use crate::data::lid_node::{GreenAmptState, LidNodeBottom, LidNodeLayerKind};
use crate::data::node::{NodeType, StorageShape};
use crate::data::storage_geometry::storage_shape_is_geometric;
use crate::data::table::table_lookup_ex;
use crate::hydrology::infiltration::grnampt_init;
use crate::hydrology::lid_node::layers;

#[derive(Debug, Clone, Default)]
pub struct ExfilSoA {
    pub count: i32,
    /// Storage node index
    pub node_idx: Vec<i32>,
    /// Bottom area (ft2)
    pub btm_area: Vec<f64>,
    /// Depth where bank seepage starts
    pub bank_min_depth: Vec<f64>,
    /// Depth at max bank area
    pub bank_max_depth: Vec<f64>,
    /// Maximum bank seepage area (ft2)
    pub bank_max_area: Vec<f64>,
    /// Bottom Green-Ampt state
    pub btm_ga: Vec<GreenAmptState>,
    /// Bank Green-Ampt state
    pub bank_ga: Vec<GreenAmptState>,
}

impl ExfilSoA {
    pub fn resize(&mut self, n: i32) {
        self.count = n;
        let un = n as usize;
        fn a<T: Clone>(v: &mut Vec<T>, n: usize, d: T) {
            v.clear();
            v.resize(n, d);
        }
        a(&mut self.node_idx, un, -1);
        a(&mut self.btm_area, un, 0.0);
        a(&mut self.bank_min_depth, un, 0.0);
        a(&mut self.bank_max_depth, un, 0.0);
        a(&mut self.bank_max_area, un, 0.0);
        self.btm_ga.resize(un, GreenAmptState::default());
        self.bank_ga.resize(un, GreenAmptState::default());
    }
}

/// lidnode::bottomConductivity (LidNode.cpp).
fn bottom_conductivity(c: &SimulationContext, r: i32) -> f64 {
    if r < 0 {
        return 0.0;
    }
    let st = &c.node_subtypes.storages;
    let ur = r as usize;
    if st.lid[ur].control >= 0 {
        if st.lid[ur].bottom != LidNodeBottom::Current {
            return 0.0;
        }
        let stack = layers(c, st.lid[ur].control);
        if let Some(last) = stack.last() {
            if last.kind == LidNodeLayerKind::Bottom {
                return last.params[0];
            }
        }
    }
    st.exfil_ksat[ur]
}

#[derive(Debug, Clone, Default)]
pub struct ExfilSolver {
    pub soa_: ExfilSoA,
}

impl ExfilSolver {
    pub fn init(&mut self, ctx: &mut SimulationContext) {
        let n_nodes = ctx.nodes.count();

        // First pass: count storage nodes that have exfiltration parameters.
        let mut n_exfil = 0i32;
        for i in 0..n_nodes {
            let ui = i as usize;
            let sr = ctx.node_subtypes.storage_row(i);
            let ksat = if sr >= 0 { bottom_conductivity(ctx, sr) } else { 0.0 };
            if ctx.nodes.type_[ui] == NodeType::Storage && ksat > 0.0 {
                n_exfil += 1;
            }
        }

        if n_exfil == 0 {
            self.soa_.count = 0;
            return;
        }

        // Second pass: populate ExfilSoA
        self.soa_.resize(n_exfil);

        const BIG: f64 = 1.0E10;
        let mut k = 0usize;
        for i in 0..n_nodes {
            let ui = i as usize;
            let sr = ctx.node_subtypes.storage_row(i);
            let st = &ctx.node_subtypes.storages;
            let exf_suction = if sr >= 0 { st.exfil_suction[sr as usize] } else { 0.0 };
            let exf_ksat = if sr >= 0 { bottom_conductivity(ctx, sr) } else { 0.0 };
            let exf_imd = if sr >= 0 { st.exfil_imd[sr as usize] } else { 0.0 };
            if ctx.nodes.type_[ui] != NodeType::Storage || exf_ksat <= 0.0 {
                continue;
            }

            let uk = k;
            self.soa_.node_idx[uk] = i;

            grnampt_init(&mut self.soa_.btm_ga[uk], exf_suction, exf_ksat, exf_imd, &ctx.options);
            grnampt_init(&mut self.soa_.bank_ga[uk], exf_suction, exf_ksat, exf_imd, &ctx.options);

            let curve_idx = if sr >= 0 { st.curve[sr as usize] } else { -1 };

            if curve_idx >= 0 {
                // TABULAR: storage shape given by a storage curve
                let curve = &ctx.tables.tables[curve_idx as usize];

                self.soa_.btm_area[uk] = table_lookup_ex(curve, 0.0);

                self.soa_.bank_min_depth[uk] = 0.0;
                self.soa_.bank_max_depth[uk] = 0.0;
                self.soa_.bank_max_area[uk] = 0.0;

                if !curve.x.is_empty() {
                    let mut alast = curve.y[0];
                    for ci in 1..curve.x.len() {
                        let d = curve.x[ci];
                        let a = curve.y[ci];

                        if a < alast {
                            break;
                        } else if a > alast {
                            self.soa_.bank_max_area[uk] = a;
                            self.soa_.bank_max_depth[uk] = d;
                        } else if self.soa_.bank_max_area[uk] == 0.0 {
                            self.soa_.bank_min_depth[uk] = d;
                        } else {
                            break;
                        }
                        alast = a;
                    }
                }

                // Legacy converts the TABULAR geometry from user to internal units here.
                {
                    let ucf_len = ucf::ucf(ucf::LENGTH, &ctx.options);
                    self.soa_.btm_area[uk] /= ucf_len * ucf_len;
                    self.soa_.bank_max_area[uk] /= ucf_len * ucf_len;
                    self.soa_.bank_min_depth[uk] /= ucf_len;
                    self.soa_.bank_max_depth[uk] /= ucf_len;
                }
            } else {
                // FUNCTIONAL: area = A * depth^B + C (raw user-unit coefficients).
                let a_coeff = if sr >= 0 { st.a[sr as usize] } else { 0.0 };
                let b_coeff = if sr >= 0 { st.b[sr as usize] } else { 0.0 };
                let c_coeff = if sr >= 0 { st.c[sr as usize] } else { 0.0 };
                let sshape = if sr >= 0 { st.shape[sr as usize] } else { StorageShape::Functional };

                let mut btm = c_coeff;
                if storage_shape_is_geometric(sshape) {
                    // Geometric shapes: the bottom area is simply c.
                } else if b_coeff == 0.0 {
                    btm += a_coeff;
                }
                self.soa_.btm_area[uk] = btm;

                self.soa_.bank_min_depth[uk] = 0.0;
                self.soa_.bank_max_depth[uk] = BIG;
                self.soa_.bank_max_area[uk] = BIG;
            }

            k += 1;
        }
    }

    pub fn state(&self) -> &ExfilSoA {
        &self.soa_
    }
}
