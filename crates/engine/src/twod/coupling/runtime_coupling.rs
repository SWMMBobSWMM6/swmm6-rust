// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/coupling/RuntimeCoupling.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Runtime (C API) sources and forcings on the 2D surface: the whole header (wave 3 carried the
// data and the queries initialize reads; wave 4b adds the rest). Only the swmm_2d_* API fills
// the rows, and that API is not translated, so in every run the translation makes they hold
// no rows; the run path reaches active, get, rate, nextTime and expire. RuntimeReceipt is in
// data/boundary_data.rs.
//
// Differences in form, never in result:
//   * The unordered_map index members are HashMaps: find and insert use them, and clear(id)
//     iterates one only to set flags, so their order is not observable.
//   * `std::move` of a RuntimeReceipt leaves the source's doubles and empties its two vectors
//     (libstdc++'s moved-from vector), which `move_receipt` reproduces (replace, U550).
//   * rebuild's `cells.at(cell)` throws std::out_of_range for a cell outside the new size;
//     here it returns Err with libstdc++'s what() (the rows before it are already placed, as
//     upstream's are). The API that calls it is not translated.
//   * find returns the row's index (upstream a pointer into rows).
#![allow(clippy::all)]

use std::collections::HashMap;

use crate::core::cmath::cmin;
use crate::core::options::Bytes;
pub use crate::twod::data::boundary_data::RuntimeReceipt;

/// runtimeWaterHeatCapacity: J/(m3 K).
pub const RUNTIME_WATER_HEAT_CAPACITY: f64 = 4186000.0;

#[derive(Debug, Clone)]
pub struct RuntimeSource {
    pub id: Bytes,
    pub owner: Bytes,
    pub cell: i32,
    pub flow: f64,
    pub heat: f64,
    pub boundary_edge: i32,
    pub boundary_head: f64,
    pub until: f64,
    pub enabled: bool,
    pub concentrations: Vec<f64>,
    pub rates: Vec<f64>,
    pub receipt: RuntimeReceipt,
}

impl Default for RuntimeSource {
    fn default() -> Self {
        RuntimeSource {
            id: Bytes::new(),
            owner: Bytes::new(),
            cell: -1,
            flow: 0.0,
            heat: 0.0,
            boundary_edge: -1,
            boundary_head: f64::NAN,
            until: f64::INFINITY,
            enabled: true,
            concentrations: Vec::new(),
            rates: Vec::new(),
            receipt: RuntimeReceipt::default(),
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct RuntimeSources {
    pub rows: Vec<RuntimeSource>,
    pub cells: Vec<Vec<usize>>,
    pub index: HashMap<Bytes, HashMap<i32, usize>>,
}

/// `dst = std::move(src)` for a RuntimeReceipt: the doubles are copied, the vectors move
/// (the source's become empty).
fn move_receipt(src: &mut RuntimeReceipt) -> RuntimeReceipt {
    RuntimeReceipt {
        requested: src.requested,
        applied: src.applied,
        water_in: src.water_in,
        water_out: src.water_out,
        heat_requested: src.heat_requested,
        heat_applied: src.heat_applied,
        last_flow: src.last_flow,
        last_dt: src.last_dt,
        requested_mass: std::mem::take(&mut src.requested_mass),
        applied_mass: std::mem::take(&mut src.applied_mass),
    }
}

/// libstdc++'s std::vector::at range-check what().
fn range_check(n: usize, size: usize) -> Bytes {
    format!("vector::_M_range_check: __n (which is {}) >= this->size() (which is {})", n, size).into_bytes()
}

impl RuntimeSources {
    /// find(id, cell): the row's index.
    pub fn find(&self, id: &[u8], cell: i32) -> Option<usize> {
        let a = self.index.get(id)?;
        a.get(&cell).copied()
    }

    /// rebuild(count): the per-cell lists of enabled rows.
    pub fn rebuild(&mut self, count: i32) -> Result<(), Bytes> {
        self.cells.clear();
        self.cells.resize(count as usize, Vec::new());
        for k in 0..self.rows.len() {
            if self.rows[k].enabled {
                let c = self.rows[k].cell as usize; // size_type conversion of an int
                if c >= self.cells.len() {
                    return Err(range_check(c, self.cells.len()));
                }
                self.cells[c].push(k);
            }
        }
        Ok(())
    }

    /// replace(value): a new (id, cell) row is appended with zeroed mass receipts. For an
    /// existing row upstream moves the row's receipt into `value` and returns, so `value` is
    /// dropped and the row keeps its old values, its receipt's two mass vectors emptied (U550).
    pub fn replace(&mut self, mut value: RuntimeSource) {
        if let Some(k) = self.find(&value.id, value.cell) {
            value.receipt = move_receipt(&mut self.rows[k].receipt);
            return;
        }
        let n = self.rows.len();
        self.index.entry(value.id.clone()).or_default().insert(value.cell, n);
        value.receipt.requested_mass.clear();
        value.receipt.requested_mass.resize(value.concentrations.len(), 0.0);
        value.receipt.applied_mass.clear();
        value.receipt.applied_mass.resize(value.concentrations.len(), 0.0);
        self.rows.push(value);
    }

    /// prune(): drop disabled rows from every cell list (order kept).
    pub fn prune(&mut self) {
        let rows = &self.rows;
        for cell in self.cells.iter_mut() {
            cell.retain(|&k| rows[k].enabled);
        }
    }

    /// clear(id, cell): disable one (id, cell) row, or every row of `id` (cell < 0).
    pub fn clear(&mut self, id: &[u8], cell: i32) -> bool {
        let Some(found) = self.index.get(id) else {
            return false;
        };
        if cell >= 0 {
            let Some(k) = self.find(id, cell) else {
                return false;
            };
            self.rows[k].enabled = false;
            if (cell as usize) < self.cells.len() {
                let rows = &self.rows;
                self.cells[cell as usize].retain(|&k| rows[k].enabled);
            }
            return true;
        } else {
            let ks: Vec<usize> = found.values().copied().collect();
            for k in ks {
                self.rows[k].enabled = false;
            }
        }
        self.prune();
        true
    }

    /// syncReceipts(old): each row takes the receipt of its (id, cell) row in `old`.
    pub fn sync_receipts(&mut self, old: &RuntimeSources) {
        for row in self.rows.iter_mut() {
            if let Some(k) = old.find(&row.id, row.cell) {
                row.receipt = old.rows[k].receipt.clone();
            }
        }
    }

    pub fn active(&self) -> bool {
        self.rows.iter().any(|r| r.enabled)
    }

    /// nextTime(now): the earliest enabled `until` after now (infinity if none).
    pub fn next_time(&self, now: f64) -> f64 {
        let mut next = f64::INFINITY;
        for r in &self.rows {
            if r.enabled && r.until > now {
                next = cmin(next, r.until);
            }
        }
        next
    }

    /// expire(now): disable every enabled row whose `until` has passed; prune if any.
    pub fn expire(&mut self, now: f64) -> bool {
        let mut changed = false;
        for r in self.rows.iter_mut() {
            if r.enabled && r.until <= now {
                r.enabled = false;
                changed = true;
            }
        }
        if changed {
            self.prune();
        }
        changed
    }
}

#[derive(Debug, Clone)]
pub struct RuntimeForcing {
    pub channel: i32,
    pub cell: i32,
    pub mode: i32,
    pub rate: f64,
    pub until: f64,
    pub enabled: bool,
    pub owner: Bytes,
    pub concentrations: Vec<f64>,
    pub receipt: RuntimeReceipt,
}

impl Default for RuntimeForcing {
    fn default() -> Self {
        RuntimeForcing {
            channel: 0,
            cell: -1,
            mode: 1,
            rate: 0.0,
            until: f64::INFINITY,
            enabled: true,
            owner: Bytes::new(),
            concentrations: Vec::new(),
            receipt: RuntimeReceipt::default(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct RuntimeForcings {
    pub rows: Vec<RuntimeForcing>,
    pub cells: Vec<[i32; 4]>,
    pub uniform: [i32; 4],
    pub index: HashMap<i32, HashMap<i32, usize>>,
}

impl Default for RuntimeForcings {
    fn default() -> Self {
        RuntimeForcings { rows: Vec::new(), cells: Vec::new(), uniform: [-1, -1, -1, -1], index: HashMap::new() }
    }
}

impl RuntimeForcings {
    /// find(channel, cell): the row's index.
    pub fn find(&self, channel: i32, cell: i32) -> Option<usize> {
        let a = self.index.get(&channel)?;
        a.get(&cell).copied()
    }

    /// get(channel, cell): the cell's own enabled row, else the uniform one.
    pub fn get(&self, channel: i32, cell: i32) -> Option<&RuntimeForcing> {
        let ch = (channel - 1) as usize;
        let mut k = if cell >= 0 && (cell as usize) < self.cells.len() { self.cells[cell as usize][ch] } else { -1 };
        if k < 0 || !self.rows[k as usize].enabled {
            k = self.uniform[ch];
        }
        if k >= 0 && self.rows[k as usize].enabled {
            Some(&self.rows[k as usize])
        } else {
            None
        }
    }

    /// rate(channel, cell, native): the native rate, replaced (mode 1) or increased (mode 2)
    /// by a runtime forcing.
    pub fn rate(&self, channel: i32, cell: i32, native: f64) -> f64 {
        match self.get(channel, cell) {
            None => native,
            Some(row) => {
                if row.mode == 2 {
                    native + row.rate
                } else {
                    row.rate
                }
            }
        }
    }

    /// replace(row): a new (channel, cell) row is appended. For an existing one upstream moves
    /// the old row's receipt into `row` and drops `row`, so the old row keeps its values, its
    /// receipt's mass vectors emptied (U550).
    pub fn replace(&mut self, mut row: RuntimeForcing) {
        if let Some(k) = self.find(row.channel, row.cell) {
            row.receipt = move_receipt(&mut self.rows[k].receipt);
        } else {
            let n = self.rows.len();
            self.index.entry(row.channel).or_default().insert(row.cell, n);
            self.rows.push(row);
        }
    }

    /// rebuild(count): the per-cell and uniform row of each channel.
    pub fn rebuild(&mut self, count: i32) -> Result<(), Bytes> {
        self.cells.clear();
        self.cells.resize(count as usize, [-1, -1, -1, -1]);
        self.uniform = [-1, -1, -1, -1];
        for k in 0..self.rows.len() {
            if self.rows[k].enabled {
                let row = &self.rows[k];
                let ch = (row.channel - 1) as usize;
                if row.cell < 0 {
                    self.uniform[ch] = k as i32;
                } else {
                    let c = row.cell as usize;
                    if c >= self.cells.len() {
                        return Err(range_check(c, self.cells.len()));
                    }
                    self.cells[c][ch] = k as i32;
                }
            }
        }
        Ok(())
    }

    pub fn active(&self) -> bool {
        self.rows.iter().any(|r| r.enabled)
    }

    /// nextTime(now): the earliest enabled `until` after now (infinity if none).
    pub fn next_time(&self, now: f64) -> f64 {
        let mut time = f64::INFINITY;
        for row in &self.rows {
            if row.enabled && row.until > now {
                time = cmin(time, row.until);
            }
        }
        time
    }

    /// expire(now): disable every enabled row whose `until` has passed (no rebuild).
    pub fn expire(&mut self, now: f64) {
        for row in self.rows.iter_mut() {
            if row.enabled && row.until <= now {
                row.enabled = false;
            }
        }
    }
}
