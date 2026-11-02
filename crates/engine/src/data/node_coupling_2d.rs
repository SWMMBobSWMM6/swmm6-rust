// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/NodeCoupling2DData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// (Upstream's file has no SPDX/license comment block; its doc comment carries the
// copyright and "@license Apache-2.0", reproduced above.)
//
// Sparse per-node [2D_NODE_COUPLING] rows keyed on the node index. Parameters are stored
// as authored; NaN means not authored. The unauthored NaN is built from the bits of
// std::numeric_limits<double>::quiet_NaN() (0x7FF8000000000000) rather than f64::NAN,
// whose bit pattern Rust does not promise. `type` is `type_`; unsetNodeCoupling2DParams
// is unset_node_coupling_2d_params. erase_row with a bad row is undefined in C++ and a
// panic here.

/// Coupling structure type of a node. NONE = not coupleable.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i8)]
pub enum NodeCoupling2DType {
    None = 0,
    Manhole = 1,
    Grate = 2,
    Curb = 3,
    /// Open-top storage / pond under the mesh.
    Open = 4,
    /// Tailwater + discharge mechanism, not an exchange law.
    Outfall = 5,
}

/// Optional KEY=VALUE parameters of a [2D_NODE_COUPLING] row, in write order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum NodeCoupling2DKey {
    /// Manhole diameter (length).
    D = 0,
    /// Weir crest length (length).
    Perimeter = 1,
    /// Orifice / open area (length^2).
    Area = 2,
    /// 1D chamber plan area above the rim (length^2).
    ChamberArea = 3,
    /// Open share of a covered lid, (0, 1].
    LidOpenFraction = 4,
    /// Crest policy (NodeCoupling2DCrest).
    Crest = 5,
    /// Coefficient preset (NodeCoupling2DCoeffs).
    Coeffs = 6,
    /// Free-weir coefficient.
    Cw = 7,
    /// Submerged-weir coefficient.
    Csw = 8,
    /// Orifice coefficient, surface to sewer.
    Coi = 9,
    /// Orifice coefficient, sewer to surface.
    Coo = 10,
    /// Per-node exchange law (NodeCoupling2DLaw).
    Law = 11,
    /// Outfall discharge spread radius (length).
    Spread = 12,
    /// Outfall tailwater relaxation time (s).
    TailwaterRelax = 13,
    Count = 14,
}

pub const K_NODE_COUPLING_2D_KEYS: i32 = NodeCoupling2DKey::Count as i32;

/// CREST policy values (stored in the CREST parameter slot).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum NodeCoupling2DCrest {
    NodeRim = 0,
    MeshBed = 1,
    Max = 2,
    Min = 3,
}

/// COEFFS preset values (stored in the COEFFS parameter slot).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum NodeCoupling2DCoeffs {
    Default = 0,
    Rubinato2017 = 1,
    Rubinato2017Lower = 2,
    Rubinato2017Upper = 3,
}

/// LAW values (stored in the LAW parameter slot).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum NodeCoupling2DLaw {
    Orifice = 0,
    WeirOrifice = 1,
}

pub type NodeCoupling2DParams = [f64; K_NODE_COUPLING_2D_KEYS as usize];

/// A parameter array with every slot unauthored (quiet NaN).
pub fn unset_node_coupling_2d_params() -> NodeCoupling2DParams {
    [f64::from_bits(0x7FF8_0000_0000_0000); K_NODE_COUPLING_2D_KEYS as usize]
}

/// Sparse [2D_NODE_COUPLING] rows (one per authored node).
#[derive(Debug, Clone, Default)]
pub struct NodeCoupling2DRows {
    /// Node index.
    pub node: Vec<i32>,
    pub type_: Vec<NodeCoupling2DType>,
    /// NaN = not authored.
    pub params: Vec<NodeCoupling2DParams>,
}

impl NodeCoupling2DRows {
    pub fn count(&self) -> i32 {
        self.node.len() as i32
    }

    /// First row for `node_idx`, or -1; linear scan.
    pub fn row_of(&self, node_idx: i32) -> i32 {
        for r in 0..self.node.len() {
            if self.node[r] == node_idx {
                return r as i32;
            }
        }
        -1
    }

    /// The existing row for `node_idx` (its type unchanged), or a new unauthored row of type `t`.
    pub fn ensure_row(&mut self, node_idx: i32, t: NodeCoupling2DType) -> i32 {
        let r = self.row_of(node_idx);
        if r >= 0 {
            return r;
        }
        self.node.push(node_idx);
        self.type_.push(t);
        self.params.push(unset_node_coupling_2d_params());
        self.count() - 1
    }

    pub fn erase_row(&mut self, r: i32) {
        let d = r as isize as usize;
        self.node.remove(d);
        self.type_.remove(d);
        self.params.remove(d);
    }

    /// Drop the deleted node's first row, then shift every remaining node key above `idx`.
    pub fn erase_node(&mut self, idx: i32) {
        let r = self.row_of(idx);
        if r >= 0 {
            self.erase_row(r);
        }
        for n in &mut self.node {
            if *n > idx {
                *n -= 1;
            }
        }
    }

    pub fn clear(&mut self) {
        self.node.clear();
        self.type_.clear();
        self.params.clear();
    }
}
