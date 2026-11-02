// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/LinkSubtypes.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson (first pass by tools/soa2rs.py, finished by hand). Not an
// OpenSWMM release.
//
// The five dense link-subtype side tables (conduit, pump, orifice, weir, outlet), each
// joined to the base LinkData row by `link_idx`, and the reverse map `subtype_row`.
// Each method visits exactly the columns upstream's does: ConduitData::erase_at() leaves
// the five gw_* columns unerased, as upstream does and documents. std::lower_bound is
// lower_bound() below, libstdc++'s algorithm step for step, so the insert position is
// the same even on keys that are not sorted. Upstream's generic lambda in erase_link()
// is the private SideTable trait. A C++ erase/insert/index past the end is undefined;
// here it panics.

#![allow(clippy::all)]
use crate::data::link::{LinkData, LinkType};

/// std::lower_bound over ints as libstdc++ writes it (__lower_bound): the first position
/// whose key is not less than `val`, found by the same halving sequence.
fn lower_bound(keys: &[i32], val: i32) -> usize {
    let mut first: usize = 0;
    let mut len: usize = keys.len();
    while len > 0 {
        let half = len >> 1;
        let middle = first + half;
        if keys[middle] < val {
            first = middle + 1;
            len = len - half - 1;
        } else {
            len = half;
        }
    }
    first
}

#[derive(Debug, Clone, Default)]
pub struct ConduitData {
    pub link_idx: Vec<i32>,
    pub roughness: Vec<f64>,
    pub length: Vec<f64>,
    pub true_length: Vec<f64>,
    pub slope: Vec<f64>,
    pub mod_length: Vec<f64>,
    pub barrels: Vec<i32>,
    pub beta: Vec<f64>,
    pub rough_factor: Vec<f64>,
    pub q_full: Vec<f64>,
    pub q_max: Vec<f64>,
    pub loss_inlet: Vec<f64>,
    pub loss_outlet: Vec<f64>,
    pub loss_avg: Vec<f64>,
    pub seep_rate: Vec<f64>,
    pub culvert_code: Vec<i32>,
    pub evap_loss_rate: Vec<f64>,
    pub seep_loss_rate: Vec<f64>,
    pub normal_flow_limited: Vec<u8>,
    pub inlet_control: Vec<u8>,
    pub full_state: Vec<i8>,
    pub gw_coupled: Vec<u8>,
    pub gw_head_rel: Vec<f64>,
    pub gw_kc: Vec<f64>,
    pub gw_dc: Vec<f64>,
    pub gw_gain_max: Vec<f64>,
}

impl ConduitData {
    pub fn count(&self) -> i32 {
        self.link_idx.len() as i32
    }

    pub fn clear(&mut self) {
        self.link_idx.clear();
        self.roughness.clear();
        self.length.clear();
        self.true_length.clear();
        self.slope.clear();
        self.mod_length.clear();
        self.barrels.clear();
        self.beta.clear();
        self.rough_factor.clear();
        self.q_full.clear();
        self.q_max.clear();
        self.loss_inlet.clear();
        self.loss_outlet.clear();
        self.loss_avg.clear();
        self.seep_rate.clear();
        self.culvert_code.clear();
        self.evap_loss_rate.clear();
        self.seep_loss_rate.clear();
        self.normal_flow_limited.clear();
        self.inlet_control.clear();
        self.full_state.clear();
        self.gw_coupled.clear();
        self.gw_head_rel.clear();
        self.gw_kc.clear();
        self.gw_dc.clear();
        self.gw_gain_max.clear();
    }

    pub fn reserve(&mut self, n: i32) {
        let un = usize::try_from(n).unwrap_or(0);
        self.link_idx.reserve(un.saturating_sub(self.link_idx.len()));
        self.roughness.reserve(un.saturating_sub(self.roughness.len()));
        self.length.reserve(un.saturating_sub(self.length.len()));
        self.true_length.reserve(un.saturating_sub(self.true_length.len()));
        self.slope.reserve(un.saturating_sub(self.slope.len()));
        self.mod_length.reserve(un.saturating_sub(self.mod_length.len()));
        self.barrels.reserve(un.saturating_sub(self.barrels.len()));
        self.beta.reserve(un.saturating_sub(self.beta.len()));
        self.rough_factor.reserve(un.saturating_sub(self.rough_factor.len()));
        self.q_full.reserve(un.saturating_sub(self.q_full.len()));
        self.q_max.reserve(un.saturating_sub(self.q_max.len()));
        self.loss_inlet.reserve(un.saturating_sub(self.loss_inlet.len()));
        self.loss_outlet.reserve(un.saturating_sub(self.loss_outlet.len()));
        self.loss_avg.reserve(un.saturating_sub(self.loss_avg.len()));
        self.seep_rate.reserve(un.saturating_sub(self.seep_rate.len()));
        self.culvert_code.reserve(un.saturating_sub(self.culvert_code.len()));
        self.evap_loss_rate.reserve(un.saturating_sub(self.evap_loss_rate.len()));
        self.seep_loss_rate.reserve(un.saturating_sub(self.seep_loss_rate.len()));
        self.normal_flow_limited.reserve(un.saturating_sub(self.normal_flow_limited.len()));
        self.inlet_control.reserve(un.saturating_sub(self.inlet_control.len()));
        self.full_state.reserve(un.saturating_sub(self.full_state.len()));
        self.gw_coupled.reserve(un.saturating_sub(self.gw_coupled.len()));
        self.gw_head_rel.reserve(un.saturating_sub(self.gw_head_rel.len()));
        self.gw_kc.reserve(un.saturating_sub(self.gw_kc.len()));
        self.gw_dc.reserve(un.saturating_sub(self.gw_dc.len()));
        self.gw_gain_max.reserve(un.saturating_sub(self.gw_gain_max.len()));
    }

    pub fn add_default(&mut self, i: i32) -> i32 {
        let p = lower_bound(&self.link_idx, i);
        self.link_idx.insert(p, i);
        self.roughness.insert(p, 0.01);
        self.length.insert(p, 0.0);
        self.true_length.insert(p, 0.0);
        self.slope.insert(p, 0.0);
        self.mod_length.insert(p, 0.0);
        self.barrels.insert(p, 1);
        self.beta.insert(p, 0.0);
        self.rough_factor.insert(p, 0.0);
        self.q_full.insert(p, 0.0);
        self.q_max.insert(p, 0.0);
        self.loss_inlet.insert(p, 0.0);
        self.loss_outlet.insert(p, 0.0);
        self.loss_avg.insert(p, 0.0);
        self.seep_rate.insert(p, 0.0);
        self.culvert_code.insert(p, 0);
        self.evap_loss_rate.insert(p, 0.0);
        self.seep_loss_rate.insert(p, 0.0);
        self.normal_flow_limited.insert(p, 0);
        self.inlet_control.insert(p, 0);
        self.full_state.insert(p, 0);
        self.gw_coupled.insert(p, 0);
        self.gw_head_rel.insert(p, 0.0);
        self.gw_kc.insert(p, 0.0);
        self.gw_dc.insert(p, 0.0);
        self.gw_gain_max.insert(p, 0.0);
        p as i32
    }

    pub fn erase_at(&mut self, r: i32) {
        let p = r as usize;
        self.link_idx.remove(p);
        self.roughness.remove(p);
        self.length.remove(p);
        self.true_length.remove(p);
        self.slope.remove(p);
        self.mod_length.remove(p);
        self.barrels.remove(p);
        self.beta.remove(p);
        self.rough_factor.remove(p);
        self.q_full.remove(p);
        self.q_max.remove(p);
        self.loss_inlet.remove(p);
        self.loss_outlet.remove(p);
        self.loss_avg.remove(p);
        self.seep_rate.remove(p);
        self.culvert_code.remove(p);
        self.evap_loss_rate.remove(p);
        self.seep_loss_rate.remove(p);
        self.normal_flow_limited.remove(p);
        self.inlet_control.remove(p);
        self.full_state.remove(p);
    }
}

#[derive(Debug, Clone, Default)]
pub struct PumpData {
    pub link_idx: Vec<i32>,
    pub curve: Vec<i32>,
    pub init_state: Vec<u8>,
    pub startup: Vec<f64>,
    pub shutoff: Vec<f64>,
    pub curve_type: Vec<i32>,
}

impl PumpData {
    pub fn count(&self) -> i32 {
        self.link_idx.len() as i32
    }

    pub fn clear(&mut self) {
        self.link_idx.clear();
        self.curve.clear();
        self.init_state.clear();
        self.startup.clear();
        self.shutoff.clear();
        self.curve_type.clear();
    }

    pub fn reserve(&mut self, n: i32) {
        let un = usize::try_from(n).unwrap_or(0);
        self.link_idx.reserve(un.saturating_sub(self.link_idx.len()));
        self.curve.reserve(un.saturating_sub(self.curve.len()));
        self.init_state.reserve(un.saturating_sub(self.init_state.len()));
        self.startup.reserve(un.saturating_sub(self.startup.len()));
        self.shutoff.reserve(un.saturating_sub(self.shutoff.len()));
        self.curve_type.reserve(un.saturating_sub(self.curve_type.len()));
    }

    pub fn add_default(&mut self, i: i32) -> i32 {
        let p = lower_bound(&self.link_idx, i);
        self.link_idx.insert(p, i);
        self.curve.insert(p, -1);
        self.init_state.insert(p, 0);
        self.startup.insert(p, 0.0);
        self.shutoff.insert(p, 0.0);
        self.curve_type.insert(p, -1);
        p as i32
    }

    pub fn erase_at(&mut self, r: i32) {
        let p = r as usize;
        self.link_idx.remove(p);
        self.curve.remove(p);
        self.init_state.remove(p);
        self.startup.remove(p);
        self.shutoff.remove(p);
        self.curve_type.remove(p);
    }
}

#[derive(Debug, Clone, Default)]
pub struct OrificeData {
    pub link_idx: Vec<i32>,
    pub orifice_type: Vec<f64>,
    pub cd: Vec<f64>,
    pub orate: Vec<f64>,
}

impl OrificeData {
    pub fn count(&self) -> i32 {
        self.link_idx.len() as i32
    }

    pub fn clear(&mut self) {
        self.link_idx.clear();
        self.orifice_type.clear();
        self.cd.clear();
        self.orate.clear();
    }

    pub fn reserve(&mut self, n: i32) {
        let un = usize::try_from(n).unwrap_or(0);
        self.link_idx.reserve(un.saturating_sub(self.link_idx.len()));
        self.orifice_type.reserve(un.saturating_sub(self.orifice_type.len()));
        self.cd.reserve(un.saturating_sub(self.cd.len()));
        self.orate.reserve(un.saturating_sub(self.orate.len()));
    }

    pub fn add_default(&mut self, i: i32) -> i32 {
        let p = lower_bound(&self.link_idx, i);
        self.link_idx.insert(p, i);
        self.orifice_type.insert(p, 0.0);
        self.cd.insert(p, 0.0);
        self.orate.insert(p, 0.0);
        p as i32
    }

    pub fn erase_at(&mut self, r: i32) {
        let p = r as usize;
        self.link_idx.remove(p);
        self.orifice_type.remove(p);
        self.cd.remove(p);
        self.orate.remove(p);
    }
}

#[derive(Debug, Clone, Default)]
pub struct WeirData {
    pub link_idx: Vec<i32>,
    pub weir_type: Vec<f64>,
    pub cd: Vec<f64>,
    pub end_contractions: Vec<f64>,
    pub crest_height: Vec<f64>,
    pub cd2: Vec<f64>,
    pub can_surcharge: Vec<u8>,
    pub road_width: Vec<f64>,
    pub road_surface: Vec<i8>,
    pub cd_curve: Vec<i32>,
}

impl WeirData {
    pub fn count(&self) -> i32 {
        self.link_idx.len() as i32
    }

    pub fn clear(&mut self) {
        self.link_idx.clear();
        self.weir_type.clear();
        self.cd.clear();
        self.end_contractions.clear();
        self.crest_height.clear();
        self.cd2.clear();
        self.can_surcharge.clear();
        self.road_width.clear();
        self.road_surface.clear();
        self.cd_curve.clear();
    }

    pub fn reserve(&mut self, n: i32) {
        let un = usize::try_from(n).unwrap_or(0);
        self.link_idx.reserve(un.saturating_sub(self.link_idx.len()));
        self.weir_type.reserve(un.saturating_sub(self.weir_type.len()));
        self.cd.reserve(un.saturating_sub(self.cd.len()));
        self.end_contractions.reserve(un.saturating_sub(self.end_contractions.len()));
        self.crest_height.reserve(un.saturating_sub(self.crest_height.len()));
        self.cd2.reserve(un.saturating_sub(self.cd2.len()));
        self.can_surcharge.reserve(un.saturating_sub(self.can_surcharge.len()));
        self.road_width.reserve(un.saturating_sub(self.road_width.len()));
        self.road_surface.reserve(un.saturating_sub(self.road_surface.len()));
        self.cd_curve.reserve(un.saturating_sub(self.cd_curve.len()));
    }

    pub fn add_default(&mut self, i: i32) -> i32 {
        let p = lower_bound(&self.link_idx, i);
        self.link_idx.insert(p, i);
        self.weir_type.insert(p, 0.0);
        self.cd.insert(p, 0.0);
        self.end_contractions.insert(p, 0.0);
        self.crest_height.insert(p, 0.0);
        self.cd2.insert(p, 0.0);
        self.can_surcharge.insert(p, 1);
        self.road_width.insert(p, 0.0);
        self.road_surface.insert(p, 0);
        self.cd_curve.insert(p, -1);
        p as i32
    }

    pub fn erase_at(&mut self, r: i32) {
        let p = r as usize;
        self.link_idx.remove(p);
        self.weir_type.remove(p);
        self.cd.remove(p);
        self.end_contractions.remove(p);
        self.crest_height.remove(p);
        self.cd2.remove(p);
        self.can_surcharge.remove(p);
        self.road_width.remove(p);
        self.road_surface.remove(p);
        self.cd_curve.remove(p);
    }
}

#[derive(Debug, Clone, Default)]
pub struct OutletData {
    pub link_idx: Vec<i32>,
    pub outlet_type: Vec<f64>,
    pub crest_height: Vec<f64>,
    pub coeff: Vec<f64>,
    pub expon: Vec<f64>,
    pub curve: Vec<i32>,
}

impl OutletData {
    pub fn count(&self) -> i32 {
        self.link_idx.len() as i32
    }

    pub fn clear(&mut self) {
        self.link_idx.clear();
        self.outlet_type.clear();
        self.crest_height.clear();
        self.coeff.clear();
        self.expon.clear();
        self.curve.clear();
    }

    pub fn reserve(&mut self, n: i32) {
        let un = usize::try_from(n).unwrap_or(0);
        self.link_idx.reserve(un.saturating_sub(self.link_idx.len()));
        self.outlet_type.reserve(un.saturating_sub(self.outlet_type.len()));
        self.crest_height.reserve(un.saturating_sub(self.crest_height.len()));
        self.coeff.reserve(un.saturating_sub(self.coeff.len()));
        self.expon.reserve(un.saturating_sub(self.expon.len()));
        self.curve.reserve(un.saturating_sub(self.curve.len()));
    }

    pub fn add_default(&mut self, i: i32) -> i32 {
        let p = lower_bound(&self.link_idx, i);
        self.link_idx.insert(p, i);
        self.outlet_type.insert(p, 0.0);
        self.crest_height.insert(p, 0.0);
        self.coeff.insert(p, 0.0);
        self.expon.insert(p, 0.0);
        self.curve.insert(p, -1);
        p as i32
    }

    pub fn erase_at(&mut self, r: i32) {
        let p = r as usize;
        self.link_idx.remove(p);
        self.outlet_type.remove(p);
        self.crest_height.remove(p);
        self.coeff.remove(p);
        self.expon.remove(p);
        self.curve.remove(p);
    }
}

/// What erase_link()'s generic `drop` lambda needs of a side table.
trait SideTable {
    fn rows(&self) -> i32;
    fn key(&self, r: usize) -> i32;
    fn erase_row(&mut self, r: i32);
}

macro_rules! side_table {
    ($t:ty) => {
        impl SideTable for $t {
            fn rows(&self) -> i32 {
                self.count()
            }
            fn key(&self, r: usize) -> i32 {
                self.link_idx[r]
            }
            fn erase_row(&mut self, r: i32) {
                self.erase_at(r)
            }
        }
    };
}
side_table!(ConduitData);
side_table!(PumpData);
side_table!(OrificeData);
side_table!(WeirData);
side_table!(OutletData);

#[derive(Debug, Clone, Default)]
pub struct LinkSubtypes {
    pub conduits: ConduitData,
    pub pumps: PumpData,
    pub orifices: OrificeData,
    pub weirs: WeirData,
    pub outlets: OutletData,
    pub subtype_row: Vec<i32>,
}

impl LinkSubtypes {
    pub fn build_default_rows(&mut self, l: &LinkData) {
        self.clear();
        let n = l.count();
        for i in 0..n {
            match l.type_[i as usize] {
                LinkType::Conduit => {
                    self.conduits.add_default(i);
                }
                LinkType::Pump => {
                    self.pumps.add_default(i);
                }
                LinkType::Orifice => {
                    self.orifices.add_default(i);
                }
                LinkType::Weir => {
                    self.weirs.add_default(i);
                }
                LinkType::Outlet => {
                    self.outlets.add_default(i);
                }
            }
        }
        self.rebuild_index(n);
    }

    pub fn clear(&mut self) {
        self.conduits.clear();
        self.pumps.clear();
        self.orifices.clear();
        self.weirs.clear();
        self.outlets.clear();
        self.subtype_row.clear();
    }

    pub fn rebuild_index(&mut self, n_links: i32) {
        let un = usize::try_from(n_links).unwrap_or(0);
        self.subtype_row.clear();
        self.subtype_row.resize(un, -1);
        fn idx(subtype_row: &mut [i32], keys: &[i32]) {
            let mut r: i32 = 0;
            while r < keys.len() as i32 {
                subtype_row[keys[r as usize] as usize] = r;
                r += 1;
            }
        }
        idx(&mut self.subtype_row, &self.conduits.link_idx);
        idx(&mut self.subtype_row, &self.pumps.link_idx);
        idx(&mut self.subtype_row, &self.orifices.link_idx);
        idx(&mut self.subtype_row, &self.weirs.link_idx);
        idx(&mut self.subtype_row, &self.outlets.link_idx);
    }

    pub fn set_link_type(&mut self, links: &mut LinkData, i: i32, t: LinkType) -> i32 {
        let ui = i as usize;
        if i >= self.subtype_row.len() as i32 {
            self.subtype_row.resize(ui + 1, -1);
        }

        let old = links.type_[ui];
        let existing = self.subtype_row[ui];
        if old == t && existing >= 0 {
            return existing; // already this subtype with a row — idempotent
        }

        let mut shifted = false;
        if existing >= 0 {
            // remove the old subtype row (real re-type/convert)
            match old {
                LinkType::Conduit => self.conduits.erase_at(existing),
                LinkType::Pump => self.pumps.erase_at(existing),
                LinkType::Orifice => self.orifices.erase_at(existing),
                LinkType::Weir => self.weirs.erase_at(existing),
                LinkType::Outlet => self.outlets.erase_at(existing),
            }
            shifted = true;
        }

        links.type_[ui] = t;

        // Upstream starts `row` at -1 for a type no case matches; a LinkType cannot be one.
        let row: i32 = match t {
            LinkType::Conduit => {
                let row = self.conduits.add_default(i);
                if row != self.conduits.count() - 1 {
                    shifted = true;
                }
                row
            }
            LinkType::Pump => {
                let row = self.pumps.add_default(i);
                if row != self.pumps.count() - 1 {
                    shifted = true;
                }
                row
            }
            LinkType::Orifice => {
                let row = self.orifices.add_default(i);
                if row != self.orifices.count() - 1 {
                    shifted = true;
                }
                row
            }
            LinkType::Weir => {
                let row = self.weirs.add_default(i);
                if row != self.weirs.count() - 1 {
                    shifted = true;
                }
                row
            }
            LinkType::Outlet => {
                let row = self.outlets.add_default(i);
                if row != self.outlets.count() - 1 {
                    shifted = true;
                }
                row
            }
        };

        if shifted {
            self.rebuild_index(links.count());
            return self.subtype_row[ui];
        }
        self.subtype_row[ui] = row; // O(1) ascending-parse / end-insert path
        row
    }

    pub fn erase_link(&mut self, i: i32, n_after: i32) {
        fn drop<T: SideTable>(tbl: &mut T, i: i32) -> bool {
            let mut r: i32 = 0;
            while r < tbl.rows() {
                if tbl.key(r as usize) == i {
                    tbl.erase_row(r);
                    return true;
                }
                r += 1;
            }
            false
        }
        let _ = drop(&mut self.conduits, i)
            || drop(&mut self.pumps, i)
            || drop(&mut self.orifices, i)
            || drop(&mut self.weirs, i)
            || drop(&mut self.outlets, i);

        fn shift(keys: &mut [i32], i: i32) {
            for k in keys.iter_mut() {
                if *k > i {
                    *k -= 1;
                }
            }
        }
        shift(&mut self.conduits.link_idx, i);
        shift(&mut self.pumps.link_idx, i);
        shift(&mut self.orifices.link_idx, i);
        shift(&mut self.weirs.link_idx, i);
        shift(&mut self.outlets.link_idx, i);

        self.rebuild_index(n_after);
    }

    // O(1) row helpers (validate link_idx[r] == i, like NodeSubtypes::*_row).
    pub fn conduit_row(&self, i: i32) -> i32 {
        self.row_in(i, &self.conduits.link_idx)
    }
    pub fn pump_row(&self, i: i32) -> i32 {
        self.row_in(i, &self.pumps.link_idx)
    }
    pub fn orifice_row(&self, i: i32) -> i32 {
        self.row_in(i, &self.orifices.link_idx)
    }
    pub fn weir_row(&self, i: i32) -> i32 {
        self.row_in(i, &self.weirs.link_idx)
    }
    pub fn outlet_row(&self, i: i32) -> i32 {
        self.row_in(i, &self.outlets.link_idx)
    }

    fn row_in(&self, i: i32, keys: &[i32]) -> i32 {
        if i >= 0 && i < self.subtype_row.len() as i32 {
            let r = self.subtype_row[i as usize];
            if r >= 0 && r < keys.len() as i32 && keys[r as usize] == i {
                return r;
            }
        }
        -1
    }
}
