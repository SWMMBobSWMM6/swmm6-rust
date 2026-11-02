// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/NodeSubtypes.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson (first pass by tools/soa2rs.py, finished by hand). Not an
// OpenSWMM release.
//
// The dense node-subtype side-tables (storage, outfall, divider), one Vec per attribute
// in upstream's field order, joined to NodeData by `node_idx`, plus the reverse map
// `subtype_row`. Every method sizes, inserts and erases exactly the columns upstream's
// does, with the same seeded values. std::lower_bound is `partition_point(|&k| k < i)`,
// the same first position whose key is not less than `i`. A row index that is out of
// range is undefined in C++ and a panic here. A negative count converted to size_t
// (reserve, rebuild_index) is treated as zero, as in data/node.rs. erase_node's generic
// `drop` lambda is written as one key search per table, keeping its short-circuit
// order (storages, then outfalls, then dividers).

#![allow(clippy::all)]
use crate::core::options::Bytes;
use crate::data::lid_node::{LidNodeConfig, LidNodeState};
use crate::data::node::{DividerType, NodeData, NodeType, OutfallType, StorageShape};

// ============================================================================
// StorageData — side-table for STORAGE nodes
// ============================================================================

/// Dense SoA for storage-unit properties (one row per STORAGE node). Row r corresponds
/// to base node `node_idx[r]`.
#[derive(Debug, Clone, Default)]
pub struct StorageData {
    pub lid: Vec<LidNodeConfig>,
    pub lid_state: Vec<LidNodeState>,
    /// Base NodeData index this row belongs to (the join key).
    pub node_idx: Vec<i32>,

    /// Storage curve index into TableData (-1 = not tabulated).
    pub curve: Vec<i32>,
    /// Curve name for deferred resolution.
    pub curve_name: Vec<Bytes>,
    /// Surface-area relation. Selects how a/b/c below are interpreted.
    pub shape: Vec<StorageShape>,
    /// Area coefficient A (legacy a1).
    pub a: Vec<f64>,
    /// Area coefficient B (legacy a2).
    pub b: Vec<f64>,
    /// Area coefficient C (legacy a0). Baseline area at zero depth.
    pub c: Vec<f64>,
    /// Raw shape parameter 1 — major axis / base length L (geometric shapes only).
    pub p1: Vec<f64>,
    /// Raw shape parameter 2 — minor axis / base width W.
    pub p2: Vec<f64>,
    /// Raw shape parameter 3 — side slope Z (run/rise), or height for PARABOLOID.
    pub p3: Vec<f64>,
    /// Seepage rate (project units/day).
    pub seep_rate: Vec<f64>,
    /// Fraction of potential evaporation realized (0-1).
    pub evap_frac: Vec<f64>,
    /// Evaporation loss this timestep (ft3).
    pub evap_loss: Vec<f64>,
    /// Exfiltration loss this timestep (ft3).
    pub exfil_loss: Vec<f64>,
    /// Raw (uncapped) exfiltration rate this timestep (cfs).
    pub exfil_rate: Vec<f64>,
    /// Green-Ampt suction head for exfiltration.
    pub exfil_suction: Vec<f64>,
    /// Green-Ampt saturated conductivity for exfiltration.
    pub exfil_ksat: Vec<f64>,
    /// Green-Ampt initial moisture deficit for exfiltration (0-1).
    pub exfil_imd: Vec<f64>,
}

impl StorageData {
    pub fn count(&self) -> i32 {
        self.node_idx.len() as i32
    }

    pub fn clear(&mut self) {
        self.node_idx.clear();
        self.curve.clear();
        self.curve_name.clear();
        self.shape.clear();
        self.a.clear();
        self.b.clear();
        self.c.clear();
        self.p1.clear();
        self.p2.clear();
        self.p3.clear();
        self.seep_rate.clear();
        self.evap_frac.clear();
        self.evap_loss.clear();
        self.exfil_loss.clear();
        self.exfil_rate.clear();
        self.exfil_suction.clear();
        self.exfil_ksat.clear();
        self.exfil_imd.clear();
        self.lid.clear();
        self.lid_state.clear();
    }

    pub fn reserve(&mut self, n: i32) {
        let un = usize::try_from(n).unwrap_or(0);
        self.node_idx.reserve(un.saturating_sub(self.node_idx.len()));
        self.curve.reserve(un.saturating_sub(self.curve.len()));
        self.curve_name.reserve(un.saturating_sub(self.curve_name.len()));
        self.shape.reserve(un.saturating_sub(self.shape.len()));
        self.a.reserve(un.saturating_sub(self.a.len()));
        self.b.reserve(un.saturating_sub(self.b.len()));
        self.c.reserve(un.saturating_sub(self.c.len()));
        self.p1.reserve(un.saturating_sub(self.p1.len()));
        self.p2.reserve(un.saturating_sub(self.p2.len()));
        self.p3.reserve(un.saturating_sub(self.p3.len()));
        self.seep_rate.reserve(un.saturating_sub(self.seep_rate.len()));
        self.evap_frac.reserve(un.saturating_sub(self.evap_frac.len()));
        self.evap_loss.reserve(un.saturating_sub(self.evap_loss.len()));
        self.exfil_loss.reserve(un.saturating_sub(self.exfil_loss.len()));
        self.exfil_rate.reserve(un.saturating_sub(self.exfil_rate.len()));
        self.exfil_suction.reserve(un.saturating_sub(self.exfil_suction.len()));
        self.exfil_ksat.reserve(un.saturating_sub(self.exfil_ksat.len()));
        self.exfil_imd.reserve(un.saturating_sub(self.exfil_imd.len()));
        self.lid.reserve(un.saturating_sub(self.lid.len()));
        self.lid_state.reserve(un.saturating_sub(self.lid_state.len()));
    }

    /// Insert a default row at the sorted key position (std::lower_bound) and return
    /// it. Seeds FUNCTIONAL shape, curve = -1, zero coefficients/losses and default LID
    /// config/state.
    pub fn add_default(&mut self, i: i32) -> i32 {
        let p = self.node_idx.partition_point(|&k| k < i);
        self.node_idx.insert(p, i);
        self.curve.insert(p, -1);
        self.curve_name.insert(p, Bytes::new());
        // FUNCTIONAL is the historical default for a curve-less storage row: it is what
        // `curve < 0` has always meant, so existing callers see no behaviour change.
        self.shape.insert(p, StorageShape::Functional);
        self.a.insert(p, 0.0);
        self.b.insert(p, 0.0);
        self.c.insert(p, 0.0);
        self.p1.insert(p, 0.0);
        self.p2.insert(p, 0.0);
        self.p3.insert(p, 0.0);
        self.seep_rate.insert(p, 0.0);
        self.evap_frac.insert(p, 0.0);
        self.evap_loss.insert(p, 0.0);
        self.exfil_loss.insert(p, 0.0);
        self.exfil_rate.insert(p, 0.0);
        self.exfil_suction.insert(p, 0.0);
        self.exfil_ksat.insert(p, 0.0);
        self.exfil_imd.insert(p, 0.0);
        self.lid.insert(p, LidNodeConfig::default());
        self.lid_state.insert(p, LidNodeState::default());
        p as i32
    }

    /// Erase dense row `r` (not a base-node index) from every column.
    pub fn erase_at(&mut self, r: i32) {
        let p = r as usize;
        self.node_idx.remove(p);
        self.curve.remove(p);
        self.curve_name.remove(p);
        self.shape.remove(p);
        self.a.remove(p);
        self.b.remove(p);
        self.c.remove(p);
        self.p1.remove(p);
        self.p2.remove(p);
        self.p3.remove(p);
        self.seep_rate.remove(p);
        self.evap_frac.remove(p);
        self.evap_loss.remove(p);
        self.exfil_loss.remove(p);
        self.exfil_rate.remove(p);
        self.exfil_suction.remove(p);
        self.exfil_ksat.remove(p);
        self.exfil_imd.remove(p);
        self.lid.remove(p);
        self.lid_state.remove(p);
    }
}

// ============================================================================
// OutfallData — side-table for OUTFALL nodes
// ============================================================================

/// Dense SoA for outfall boundary-condition properties (one row per OUTFALL).
#[derive(Debug, Clone, Default)]
pub struct OutfallData {
    /// Base NodeData index this row belongs to (the join key).
    pub node_idx: Vec<i32>,

    /// Boundary condition type (FREE/NORMAL/FIXED/TIDAL/TIMESERIES).
    pub bc_type: Vec<OutfallType>,
    /// Fixed stage, or tidal/timeseries table index (per bc_type). For
    /// TIDAL/TIMESERIES an unresolved reference is -1 (never 0).
    pub param: Vec<f64>,
    /// Stage-data table name for TIDAL/TIMESERIES (empty otherwise).
    pub param_name: Vec<Bytes>,
    /// Flap gate present (0/1).
    pub has_flap_gate: Vec<u8>,
    /// Subcatchment index to route discharge to (-1 = none).
    pub route_to: Vec<i32>,
    /// Cached connected-conduit index (-1 = none).
    pub link_idx: Vec<i32>,
    /// Conduit offset at the outfall end.
    pub link_offset: Vec<f64>,
    /// Cached 2D surface head at the coupling point (sentinel -1e30).
    pub head_2d: Vec<f64>,
    /// Cached wet/dry ramp factor [0,1] for the 2D tailwater override; default 0.
    pub ramp_2d: Vec<f64>,
}

impl OutfallData {
    pub fn count(&self) -> i32 {
        self.node_idx.len() as i32
    }

    pub fn clear(&mut self) {
        self.node_idx.clear();
        self.bc_type.clear();
        self.param.clear();
        self.param_name.clear();
        self.has_flap_gate.clear();
        self.route_to.clear();
        self.link_idx.clear();
        self.link_offset.clear();
        self.head_2d.clear();
        self.ramp_2d.clear();
    }

    pub fn reserve(&mut self, n: i32) {
        let un = usize::try_from(n).unwrap_or(0);
        self.node_idx.reserve(un.saturating_sub(self.node_idx.len()));
        self.bc_type.reserve(un.saturating_sub(self.bc_type.len()));
        self.param.reserve(un.saturating_sub(self.param.len()));
        self.param_name.reserve(un.saturating_sub(self.param_name.len()));
        self.has_flap_gate.reserve(un.saturating_sub(self.has_flap_gate.len()));
        self.route_to.reserve(un.saturating_sub(self.route_to.len()));
        self.link_idx.reserve(un.saturating_sub(self.link_idx.len()));
        self.link_offset.reserve(un.saturating_sub(self.link_offset.len()));
        self.head_2d.reserve(un.saturating_sub(self.head_2d.len()));
        self.ramp_2d.reserve(un.saturating_sub(self.ramp_2d.len()));
    }

    /// Insert a default row at the sorted key position and return it. Seeds FREE
    /// boundary, no route/link references, no flap gate, head_2d = -1e30, zero ramp.
    pub fn add_default(&mut self, i: i32) -> i32 {
        let p = self.node_idx.partition_point(|&k| k < i);
        self.node_idx.insert(p, i);
        self.bc_type.insert(p, OutfallType::Free);
        self.param.insert(p, 0.0);
        self.param_name.insert(p, Bytes::new());
        self.has_flap_gate.insert(p, 0u8);
        self.route_to.insert(p, -1);
        self.link_idx.insert(p, -1);
        self.link_offset.insert(p, 0.0);
        self.head_2d.insert(p, -1.0e30);
        self.ramp_2d.insert(p, 0.0);
        p as i32
    }

    /// Erase dense row `r` (not a base-node index) from every column.
    pub fn erase_at(&mut self, r: i32) {
        let p = r as usize;
        self.node_idx.remove(p);
        self.bc_type.remove(p);
        self.param.remove(p);
        self.param_name.remove(p);
        self.has_flap_gate.remove(p);
        self.route_to.remove(p);
        self.link_idx.remove(p);
        self.link_offset.remove(p);
        self.head_2d.remove(p);
        self.ramp_2d.remove(p);
    }
}

// ============================================================================
// DividerData — side-table for DIVIDER nodes
// ============================================================================

/// Dense SoA for flow-divider properties (one row per DIVIDER node).
#[derive(Debug, Clone, Default)]
pub struct DividerData {
    /// Base NodeData index this row belongs to (the join key).
    pub node_idx: Vec<i32>,

    /// Diversion method (CUTOFF/OVERFLOW/TABULAR/WEIR).
    pub method: Vec<DividerType>,
    /// Cutoff flow for CUTOFF dividers.
    pub cutoff: Vec<f64>,
    /// Weir discharge coefficient for WEIR dividers.
    pub cd: Vec<f64>,
    /// Weir max depth for WEIR dividers.
    pub max_depth: Vec<f64>,
    /// Diversion curve index for TABULAR dividers (-1 = none).
    pub curve: Vec<i32>,
    /// Diversion link index (-1 = not set).
    pub link: Vec<i32>,
    /// Diversion link name (deferred resolution).
    pub link_name: Vec<Bytes>,
    /// Diversion curve name (deferred resolution, TABULAR only).
    pub curve_name: Vec<Bytes>,
}

impl DividerData {
    pub fn count(&self) -> i32 {
        self.node_idx.len() as i32
    }

    pub fn clear(&mut self) {
        self.node_idx.clear();
        self.method.clear();
        self.cutoff.clear();
        self.cd.clear();
        self.max_depth.clear();
        self.curve.clear();
        self.link.clear();
        self.link_name.clear();
        self.curve_name.clear();
    }

    pub fn reserve(&mut self, n: i32) {
        let un = usize::try_from(n).unwrap_or(0);
        self.node_idx.reserve(un.saturating_sub(self.node_idx.len()));
        self.method.reserve(un.saturating_sub(self.method.len()));
        self.cutoff.reserve(un.saturating_sub(self.cutoff.len()));
        self.cd.reserve(un.saturating_sub(self.cd.len()));
        self.max_depth.reserve(un.saturating_sub(self.max_depth.len()));
        self.curve.reserve(un.saturating_sub(self.curve.len()));
        self.link.reserve(un.saturating_sub(self.link.len()));
        self.link_name.reserve(un.saturating_sub(self.link_name.len()));
        self.curve_name.reserve(un.saturating_sub(self.curve_name.len()));
    }

    /// Insert a default row at the sorted key position and return it. Seeds CUTOFF
    /// method, unset curve/link references and zero parameters.
    pub fn add_default(&mut self, i: i32) -> i32 {
        let p = self.node_idx.partition_point(|&k| k < i);
        self.node_idx.insert(p, i);
        self.method.insert(p, DividerType::Cutoff);
        self.cutoff.insert(p, 0.0);
        self.cd.insert(p, 0.0);
        self.max_depth.insert(p, 0.0);
        self.curve.insert(p, -1);
        self.link.insert(p, -1);
        self.link_name.insert(p, Bytes::new());
        self.curve_name.insert(p, Bytes::new());
        p as i32
    }

    /// Erase dense row `r` (not a base-node index) from every column.
    pub fn erase_at(&mut self, r: i32) {
        let p = r as usize;
        self.node_idx.remove(p);
        self.method.remove(p);
        self.cutoff.remove(p);
        self.cd.remove(p);
        self.max_depth.remove(p);
        self.curve.remove(p);
        self.link.remove(p);
        self.link_name.remove(p);
        self.curve_name.remove(p);
    }
}

// ============================================================================
// NodeSubtypes — container + reverse map
// ============================================================================

/// The three node subtype side-tables plus the reverse index map. `subtype_row[i]` is
/// the row of node `i` within whichever side-table matches `nodes.type_[i]`, or -1 for
/// JUNCTION nodes.
#[derive(Debug, Clone, Default)]
pub struct NodeSubtypes {
    pub storages: StorageData,
    pub outfalls: OutfallData,
    pub dividers: DividerData,

    /// base node index → row in its subtype table (-1 for junctions).
    pub subtype_row: Vec<i32>,
}

/// First r in [0, keys.len()) with keys[r] == i (the search in erase_node's `drop`).
fn first_row_with_key(keys: &[i32], i: i32) -> Option<i32> {
    for r in 0..keys.len() as i32 {
        if keys[r as usize] == i {
            return Some(r);
        }
    }
    None
}

impl NodeSubtypes {
    pub fn clear(&mut self) {
        self.storages.clear();
        self.outfalls.clear();
        self.dividers.clear();
        self.subtype_row.clear();
    }

    /// Replace the reverse map with -1 then write each side-table row index, visiting
    /// storages, outfalls, dividers in that order (later tables overwrite duplicates).
    pub fn rebuild_index(&mut self, n_nodes: i32) {
        self.subtype_row.clear();
        self.subtype_row.resize(usize::try_from(n_nodes).unwrap_or(0), -1);
        for r in 0..self.storages.count() {
            self.subtype_row[self.storages.node_idx[r as usize] as usize] = r;
        }
        for r in 0..self.outfalls.count() {
            self.subtype_row[self.outfalls.node_idx[r as usize] as usize] = r;
        }
        for r in 0..self.dividers.count() {
            self.subtype_row[self.dividers.node_idx[r as usize] as usize] = r;
        }
    }

    /// Change one base node's type and maintain its dense side-table row. Returns the
    /// target row, or -1 for JUNCTION; a node already of type `t` with a row returns
    /// that row unchanged.
    pub fn set_node_type(&mut self, nodes: &mut NodeData, i: i32, t: NodeType) -> i32 {
        let ui = i as usize;
        if i >= self.subtype_row.len() as i32 {
            self.subtype_row.resize(ui + 1, -1);
        }

        // Contexts sized column by column may omit storage_area; unsized means no chamber.
        if t != NodeType::Junction && ui < nodes.storage_area.len() {
            nodes.storage_area[ui] = 0.0;
        }
        let old = nodes.type_[ui];
        let existing = self.subtype_row[ui];
        if old == t && existing >= 0 {
            return existing; // already this subtype with a row — idempotent.
        }

        let mut shifted = false;
        if existing >= 0 {
            // remove the old subtype row (real re-type/convert).
            match old {
                NodeType::Storage => self.storages.erase_at(existing),
                NodeType::Outfall => self.outfalls.erase_at(existing),
                NodeType::Divider => self.dividers.erase_at(existing),
                _ => {}
            }
            shifted = true;
        }

        nodes.type_[ui] = t;

        let mut row = -1;
        match t {
            NodeType::Storage => {
                row = self.storages.add_default(i);
                if row != self.storages.count() - 1 {
                    shifted = true;
                }
            }
            NodeType::Outfall => {
                row = self.outfalls.add_default(i);
                if row != self.outfalls.count() - 1 {
                    shifted = true;
                }
            }
            NodeType::Divider => {
                row = self.dividers.add_default(i);
                if row != self.dividers.count() - 1 {
                    shifted = true;
                }
            }
            NodeType::Junction => {} // no subtype row
        }

        if shifted {
            self.rebuild_index(nodes.count());
            return if t == NodeType::Junction { -1 } else { self.subtype_row[ui] };
        }
        self.subtype_row[ui] = row; // O(1) ascending-parse / end-insert path
        row
    }

    /// Drop a deleted base node's subtype row (the first match across storages,
    /// outfalls, dividers), decrement surviving keys above `i`, rebuild the map. Call
    /// after NodeData::erase_at(i).
    pub fn erase_node(&mut self, i: i32, n_after: i32) {
        // Erase the row whose node_idx == i (in whichever table holds it).
        if let Some(r) = first_row_with_key(&self.storages.node_idx, i) {
            self.storages.erase_at(r);
        } else if let Some(r) = first_row_with_key(&self.outfalls.node_idx, i) {
            self.outfalls.erase_at(r);
        } else if let Some(r) = first_row_with_key(&self.dividers.node_idx, i) {
            self.dividers.erase_at(r);
        }

        let shift = |keys: &mut Vec<i32>| {
            for k in keys.iter_mut() {
                if *k > i {
                    *k -= 1;
                }
            }
        };
        shift(&mut self.storages.node_idx);
        shift(&mut self.outfalls.node_idx);
        shift(&mut self.dividers.node_idx);

        self.rebuild_index(n_after);
    }

    /// Dense storage row of base node `i`, or -1.
    pub fn storage_row(&self, i: i32) -> i32 {
        if i >= 0 && i < self.subtype_row.len() as i32 {
            let r = self.subtype_row[i as usize];
            if r >= 0 && r < self.storages.count() && self.storages.node_idx[r as usize] == i {
                return r;
            }
        }
        -1
    }

    /// Dense outfall row of base node `i`, or -1.
    pub fn outfall_row(&self, i: i32) -> i32 {
        if i >= 0 && i < self.subtype_row.len() as i32 {
            let r = self.subtype_row[i as usize];
            if r >= 0 && r < self.outfalls.count() && self.outfalls.node_idx[r as usize] == i {
                return r;
            }
        }
        -1
    }

    /// Dense divider row of base node `i`, or -1.
    pub fn divider_row(&self, i: i32) -> i32 {
        if i >= 0 && i < self.subtype_row.len() as i32 {
            let r = self.subtype_row[i as usize];
            if r >= 0 && r < self.dividers.count() && self.dividers.node_idx[r as usize] == i {
                return r;
            }
        }
        -1
    }
}
