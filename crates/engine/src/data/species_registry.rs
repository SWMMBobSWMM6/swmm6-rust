// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/SpeciesRegistry.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Single source of truth for all transported constituents (Unified Transport §4.1, T0a).
// Names and units are bytes, compared byte-exact (case-sensitive) as std::string == does.
// add() takes the name and units by value, as upstream moves them in.

use crate::core::options::Bytes;

/// Constituent kinds (master plan §4.1).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum SpeciesKind {
    /// legacy [POLLUTANTS] entry
    Pollutant = 0,
    /// __WATER_AGE__ (phase A1)
    ReservedAge = 1,
    /// __TEMPERATURE__ (phase H1)
    ReservedTemperature = 2,
    /// reactions component BULK species (R1)
    MsxBulk = 3,
    /// reactions component WALL species (R1)
    MsxWall = 4,
}

#[derive(Debug, Clone, Default)]
pub struct SpeciesRegistry {
    name_: Vec<Bytes>,
    kind_: Vec<SpeciesKind>,
    units_: Vec<Bytes>,
    n_pollutants_: i32,
}

// n_pollutants_ = 0 is Rust's default, so Default is derived.

impl SpeciesRegistry {
    /// Erase all species and reset the pollutant count (capacities retained).
    pub fn clear(&mut self) {
        self.name_.clear();
        self.kind_.clear();
        self.units_.clear();
        self.n_pollutants_ = 0;
    }

    /// Append a uniquely named species; -1 on an exact-name collision, with no mutation.
    pub fn add(&mut self, name: Bytes, kind: SpeciesKind, units: Bytes) -> i32 {
        if self.find(&name) >= 0 {
            return -1;
        }
        self.name_.push(name);
        self.kind_.push(kind);
        self.units_.push(units);
        if kind == SpeciesKind::Pollutant {
            self.n_pollutants_ += 1;
        }
        self.name_.len() as i32 - 1
    }

    /// First byte-exact species name, or -1; linear scan over all kinds.
    pub fn find(&self, name: &[u8]) -> i32 {
        for i in 0..self.name_.len() {
            if self.name_[i] == name {
                return i as i32;
            }
        }
        -1
    }

    pub fn count(&self) -> i32 {
        self.name_.len() as i32
    }

    pub fn pollutant_count(&self) -> i32 {
        self.n_pollutants_
    }

    /// POLLUTANT, MSX_BULK, RESERVED_AGE and RESERVED_TEMPERATURE entries.
    pub fn transported_count(&self) -> i32 {
        let mut n = 0;
        for &k in &self.kind_ {
            if k == SpeciesKind::Pollutant
                || k == SpeciesKind::MsxBulk
                || k == SpeciesKind::ReservedAge
                || k == SpeciesKind::ReservedTemperature
            {
                n += 1;
            }
        }
        n
    }

    pub fn name(&self, i: i32) -> &Bytes {
        &self.name_[i as usize]
    }

    pub fn kind(&self, i: i32) -> SpeciesKind {
        self.kind_[i as usize]
    }

    pub fn units(&self, i: i32) -> &Bytes {
        &self.units_[i as usize]
    }

    /// Drop entries from index n onward and recount pollutant kinds. Negative n is zero;
    /// n >= count() leaves storage unchanged.
    pub fn truncate_to(&mut self, n: i32) {
        let mut n = n;
        if n < 0 {
            n = 0;
        }
        if n >= self.count() {
            return;
        }
        let un = n as usize;
        self.name_.truncate(un);
        self.kind_.truncate(un);
        self.units_.truncate(un);
        self.n_pollutants_ = 0;
        for &k in &self.kind_ {
            if k == SpeciesKind::Pollutant {
                self.n_pollutants_ += 1;
            }
        }
    }
}
