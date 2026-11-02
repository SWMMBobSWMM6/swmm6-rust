// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/data/MsxSurfaceData.hpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// BW-MSX surface buildup / washoff / sweeping rows, parameters and state for the
// reactions component's species. Rust has no nested type declarations, so upstream's
// MsxSurfaceRows::Buildup / ::Washoff / ::Loading are MsxSurfaceRowsBuildup /
// MsxSurfaceRowsWashoff / MsxSurfaceRowsLoading. std::vector<bool> is Vec<bool>.
// MsxSurfaceRows::empty() is is_empty(). Every member initializer is Rust's default.

use crate::core::options::Bytes;

/// MsxSurfaceRows::Buildup
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MsxSurfaceRowsBuildup {
    pub landuse: Bytes,
    pub species: Bytes,
    /// 0=NONE,1=POW,2=EXP,3=SAT,4=EXT
    pub func_type: i32,
    pub c1: f64,
    pub c2: f64,
    pub c3: f64,
    /// EXT: time series name (resolved later)
    pub ts_name: Bytes,
    /// 0=PER_AREA, 1=PER_CURB
    pub normalizer: i32,
}

/// MsxSurfaceRows::Washoff
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MsxSurfaceRowsWashoff {
    pub landuse: Bytes,
    pub species: Bytes,
    /// 0=NONE,1=EXP,2=RC,3=EMC
    pub func_type: i32,
    pub coeff: f64,
    pub expon: f64,
    pub sweep_effic: f64,
    pub bmp_effic: f64,
}

/// MsxSurfaceRows::Loading
#[derive(Debug, Clone, Default, PartialEq)]
pub struct MsxSurfaceRowsLoading {
    pub subcatch: Bytes,
    pub species: Bytes,
    /// initial buildup, mass per unit area
    pub value: f64,
}

/// Rows parked by the section handlers for names that are not pollutants.
#[derive(Debug, Clone, Default)]
pub struct MsxSurfaceRows {
    pub buildup: Vec<MsxSurfaceRowsBuildup>,
    pub washoff: Vec<MsxSurfaceRowsWashoff>,
    pub loadings: Vec<MsxSurfaceRowsLoading>,
}

impl MsxSurfaceRows {
    pub fn is_empty(&self) -> bool {
        self.buildup.is_empty() && self.washoff.is_empty() && self.loadings.is_empty()
    }
}

/// Resolved MSX surface-quality parameters and state.
#[derive(Debug, Clone, Default)]
pub struct MsxSurfaceData {
    /// parse-time, by name
    pub rows: MsxSurfaceRows,

    pub n_landuses: i32,
    /// MSX species count at resolve time
    pub n_species: i32,
    pub n_subcatch: i32,

    // ---- parameters [lu * n_species + m]
    /// 0=NONE,1=POW,2=EXP,3=SAT,4=EXT
    pub bu_type: Vec<i32>,
    pub bu_c1: Vec<f64>,
    pub bu_c2: Vec<f64>,
    /// EXT: c3 = time-series index
    pub bu_c3: Vec<f64>,
    /// 0=PER_AREA, 1=PER_CURB
    pub bu_normalizer: Vec<i32>,
    /// derived (99.9 % of max)
    pub bu_max_days: Vec<f64>,
    /// 0=NONE,1=EXP,2=RC,3=EMC
    pub wo_type: Vec<i32>,
    pub wo_coeff: Vec<f64>,
    pub wo_expon: Vec<f64>,
    pub wo_sweep_effic: Vec<f64>,
    pub wo_bmp_effic: Vec<f64>,

    /// initial loading [sc * n_species + m] (mass per unit area)
    pub init_loading: Vec<f64>,

    /// per-species unit factor, concentration-mass → user mass (D-A27)
    pub mcf: Vec<f64>,

    // ---- state
    /// bu_idx(sc, lu, m), per normalizer unit
    pub buildup: Vec<f64>,
    /// [sc * n_species + m], conc units
    pub washoff_conc: Vec<f64>,
    /// Previous runoff step's washoff_conc.
    pub washoff_conc_old: Vec<f64>,
    /// [sc * n_species + m] ponded-water species mass (concentration mass)
    pub ponded_qual: Vec<f64>,
    /// [sc * n_species + m] run-on load rate (concentration mass / s)
    pub runon_qual_rate: Vec<f64>,
    /// Per species: rainfall concentration (0 = rain falls clean)
    pub c_rain: Vec<f64>,
    /// Per species: snow-only buildup flag
    pub snow_only: Vec<bool>,
    /// LID underdrain removal fraction [control * n_species + m]
    pub lid_rmvl: Vec<f64>,
    /// [node * n_species + m] underdrain species load (concentration mass / s)
    pub lid_drain_load: Vec<f64>,
    pub lid_drain_load_old: Vec<f64>,
    /// [sc * n_species + m] the species twin of SubcatchData::lid_drain_runon_qual
    pub lid_drain_runon_load: Vec<f64>,

    // ---- ledgers per species [m], user mass
    pub led_init_buildup: Vec<f64>,
    /// net accrual over the run
    pub led_buildup: Vec<f64>,
    /// delivered to the network
    pub led_runoff_load: Vec<f64>,
    pub led_sweeping: Vec<f64>,
    pub led_bmp_removal: Vec<f64>,
    /// rainfall concentration × rain volume
    pub led_wet_deposition: Vec<f64>,
    /// infiltrated from ponded water and LID units
    pub led_infil_loss: Vec<f64>,
    /// signed wet-pool reaction loss, user mass
    pub led_reacted: Vec<f64>,
    /// numerical addition after a negative RATE
    pub led_reaction_clamped: Vec<f64>,
    /// ponded mass written off when a surface dried
    pub led_ponded_lost: Vec<f64>,
    /// [sc * n_species + m], per-subcatchment total (report)
    pub led_subcatch_load: Vec<f64>,

    pub resolved: bool,
}

impl MsxSurfaceData {
    /// resolved && n_species > 0
    pub fn active(&self) -> bool {
        self.resolved && self.n_species > 0
    }

    /// lu*n_species+m, unchecked.
    pub fn pidx(&self, lu: i32, m: i32) -> usize {
        (lu as usize) * (self.n_species as usize) + (m as usize)
    }

    /// (sc*n_landuses+lu)*n_species+m, unchecked.
    pub fn bu_idx(&self, sc: i32, lu: i32, m: i32) -> usize {
        ((sc as usize) * (self.n_landuses as usize) + (lu as usize)) * (self.n_species as usize)
            + (m as usize)
    }

    /// sc*n_species+m, unchecked.
    pub fn sidx(&self, sc: i32, m: i32) -> usize {
        (sc as usize) * (self.n_species as usize) + (m as usize)
    }

    /// Replace resolved buildup/washoff parameter matrices. Preserves c_rain only when
    /// its length already equals nm.
    pub fn resize_params(&mut self, nlu: i32, nm: i32) {
        self.n_landuses = nlu;
        self.n_species = nm;
        let n = (nlu as usize) * (nm as usize);
        self.bu_type.clear();
        self.bu_type.resize(n, 0);
        self.bu_c1.clear();
        self.bu_c1.resize(n, 0.0);
        self.bu_c2.clear();
        self.bu_c2.resize(n, 0.0);
        self.bu_c3.clear();
        self.bu_c3.resize(n, 0.0);
        self.bu_normalizer.clear();
        self.bu_normalizer.resize(n, 0);
        self.bu_max_days.clear();
        self.bu_max_days.resize(n, 0.0);
        self.wo_type.clear();
        self.wo_type.resize(n, 0);
        self.wo_coeff.clear();
        self.wo_coeff.resize(n, 0.0);
        self.wo_expon.clear();
        self.wo_expon.resize(n, 0.0);
        self.wo_sweep_effic.clear();
        self.wo_sweep_effic.resize(n, 0.0);
        self.wo_bmp_effic.clear();
        self.wo_bmp_effic.resize(n, 0.0);
        self.mcf.clear();
        self.mcf.resize(nm as usize, 1.0);
        if self.c_rain.len() != nm as usize {
            self.c_rain.clear();
            self.c_rain.resize(nm as usize, 0.0);
        }
        self.snow_only.clear();
        self.snow_only.resize(nm as usize, false);
    }

    /// Replace subcatchment state and ledgers with zero arrays. Preserves init_loading
    /// only when its shape matches; lid_rmvl and lid_drain_load(_old) are left as they are.
    pub fn resize_state(&mut self, nsc: i32) {
        self.n_subcatch = nsc;
        let nm = self.n_species as usize;
        self.buildup.clear();
        self.buildup.resize((nsc as usize) * (self.n_landuses as usize) * nm, 0.0);
        self.washoff_conc.clear();
        self.washoff_conc.resize((nsc as usize) * nm, 0.0);
        self.washoff_conc_old.clear();
        self.washoff_conc_old.resize((nsc as usize) * nm, 0.0);
        self.ponded_qual.clear();
        self.ponded_qual.resize((nsc as usize) * nm, 0.0);
        self.runon_qual_rate.clear();
        self.runon_qual_rate.resize((nsc as usize) * nm, 0.0);
        self.lid_drain_runon_load.clear();
        self.lid_drain_runon_load.resize((nsc as usize) * nm, 0.0);
        if self.init_loading.len() != (nsc as usize) * nm {
            self.init_loading.clear();
            self.init_loading.resize((nsc as usize) * nm, 0.0);
        }
        self.led_init_buildup.clear();
        self.led_init_buildup.resize(nm, 0.0);
        self.led_buildup.clear();
        self.led_buildup.resize(nm, 0.0);
        self.led_runoff_load.clear();
        self.led_runoff_load.resize(nm, 0.0);
        self.led_sweeping.clear();
        self.led_sweeping.resize(nm, 0.0);
        self.led_bmp_removal.clear();
        self.led_bmp_removal.resize(nm, 0.0);
        self.led_wet_deposition.clear();
        self.led_wet_deposition.resize(nm, 0.0);
        self.led_infil_loss.clear();
        self.led_infil_loss.resize(nm, 0.0);
        self.led_ponded_lost.clear();
        self.led_ponded_lost.resize(nm, 0.0);
        self.led_reacted.clear();
        self.led_reacted.resize(nm, 0.0);
        self.led_reaction_clamped.clear();
        self.led_reaction_clamped.resize(nm, 0.0);
        self.led_subcatch_load.clear();
        self.led_subcatch_load.resize((nsc as usize) * nm, 0.0);
    }
}
