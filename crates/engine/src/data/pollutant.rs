// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/PollutantData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson (first pass by tools/soa2rs.py, finished by hand). Not an
// OpenSWMM release.
//
// The pollutant definitions, one Vec per attribute (structure of arrays).
// resize_pollutants() replaces every column with upstream's defaults. `snow_only` is a
// std::vector<bool> upstream and a Vec<bool> here; the per-object concentration arrays
// live in the node, link and subcatchment stores, not here.

#![allow(clippy::all)]
use crate::core::options::Bytes;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i8)]
pub enum MassUnits {
    MgPerL = 0,
    UgPerL = 1,
    CountsPerL = 2,
}

#[derive(Debug, Clone, Default)]
pub struct PollutantData {
    pub units: Vec<MassUnits>,
    pub mwt: Vec<f64>,
    pub k_decay: Vec<f64>,
    pub c_rain: Vec<f64>,
    pub c_gw: Vec<f64>,
    pub c_rdii: Vec<f64>,
    pub c_dwf: Vec<f64>,
    pub init_conc: Vec<f64>,
    pub co_pollut: Vec<i32>,
    pub co_frac: Vec<f64>,
    pub snow_only: Vec<bool>,
    pub comments: Vec<Bytes>,
}

// Upstream has no member initializers here, so Default is derived.

impl PollutantData {
    pub fn n_pollutants(&self) -> i32 {
        self.units.len() as i32
    }

    pub fn resize_pollutants(&mut self, n: i32) {
        let un = usize::try_from(n).unwrap_or(0);
        self.units.clear();
        self.units.resize(un, MassUnits::MgPerL);
        self.mwt.clear();
        self.mwt.resize(un, 1.0);
        self.k_decay.clear();
        self.k_decay.resize(un, 0.0);
        self.c_rain.clear();
        self.c_rain.resize(un, 0.0);
        self.c_gw.clear();
        self.c_gw.resize(un, 0.0);
        self.c_rdii.clear();
        self.c_rdii.resize(un, 0.0);
        self.c_dwf.clear();
        self.c_dwf.resize(un, 0.0);
        self.init_conc.clear();
        self.init_conc.resize(un, 0.0);
        self.co_pollut.clear();
        self.co_pollut.resize(un, -1);
        self.co_frac.clear();
        self.co_frac.resize(un, 0.0);
        self.snow_only.clear();
        self.snow_only.resize(un, false);
        self.comments.clear();
        self.comments.resize(un, Vec::new());
    }

    pub fn shrink_to_fit(&mut self) {
        self.units.shrink_to_fit();
        self.mwt.shrink_to_fit();
        self.k_decay.shrink_to_fit();
        self.c_rain.shrink_to_fit();
        self.c_gw.shrink_to_fit();
        self.c_rdii.shrink_to_fit();
        self.c_dwf.shrink_to_fit();
        self.init_conc.shrink_to_fit();
        self.co_pollut.shrink_to_fit();
        self.co_frac.shrink_to_fit();
        self.snow_only.shrink_to_fit();
        self.comments.shrink_to_fit();
    }
}
