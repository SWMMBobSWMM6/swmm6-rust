// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/2d/coupling/RuntimeCoupling.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 12 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Runtime (C API) sources and forcings on the 2D surface: the data and the queries the run path
// reads. Only the swmm_2d_* API fills them, and that API is not translated, so in every run the
// translation makes they hold no rows. RuntimeReceipt is in data/boundary_data.rs. The
// unordered_map index members are HashMaps: only find and insert use them (no iteration), so
// their order is not observable.

#![allow(clippy::all)]

use std::collections::HashMap;

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

impl RuntimeSources {
    pub fn active(&self) -> bool {
        self.rows.iter().any(|r| r.enabled)
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

    pub fn active(&self) -> bool {
        self.rows.iter().any(|r| r.enabled)
    }
}
