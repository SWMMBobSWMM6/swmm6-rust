// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/surface/DegreeDaySnow.hpp and
// DegreeDaySnow.cpp, reached through hydrology/Snow.hpp (which only includes it); the
// module is named after Snow.hpp and namespace snow, as SWMMEngine.hpp names them
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The degree-day snow solver's per-subcatchment store (SnowSoA, resize) and the solver
// shell with init, which surface::initializeNativeSnow calls from initHydrology.
//
// Differences of form:
//   * SnowModel (surface/SnowModel.hpp) is an abstract base with one implementation in a
//     1D build's subcatchment path; the virtual methods are inherent methods here.
//   * `double adc_imperv[10]` are fixed arrays; member initializers carry into a
//     hand-written Default. `std::vector<bool> external` is Vec<bool>.
//   * The solver member is public so the gate dumper can read it; upstream's is private.
//   * setMeltCoeffs was translated on 5 November 2026 for the runoff step.
//   * NOT YET TRANSLATED (run-time step functions, phase 1c): step, pack, unpack, both
//     execute overloads, both plowSnow overloads, meltElement,
//     rainMeltRate, and the file-static helpers (mixAge, getArealSnowCover,
//     getArealDepletion, ...).

#![allow(non_snake_case)]

/// Plowable, Impervious, Pervious
pub const N_SUBAREAS: i32 = 3;
pub const SNOW_PLOWABLE: i32 = 0;
pub const SNOW_IMPERV: i32 = 1;
pub const SNOW_PERV: i32 = 2;

/// Per-subcatchment snowpack state. Per-subarea arrays are flat [subcatch * 3 + subarea].
#[derive(Debug, Clone)]
pub struct SnowSoA {
    pub n_subcatch: i32,
    pub external: Vec<bool>,
    /// UEB melt plus throughfall, ft/s
    pub external_net: Vec<f64>,

    /// Snow water equivalent (ft)
    pub wsnow: Vec<f64>,
    /// Free water in pack (ft)
    pub fw: Vec<f64>,
    /// Cold content (ft water equiv)
    pub coldc: Vec<f64>,
    /// Antecedent temperature index (deg F)
    pub ati: Vec<f64>,
    /// Areal depletion index
    pub awe: Vec<f64>,
    /// Melt rate output (ft/sec)
    pub imelt: Vec<f64>,

    /// Base melt temperature (deg F)
    pub tbase: Vec<f64>,
    /// Degree-day melt factor (ft/deg-F/sec)
    pub dhm: Vec<f64>,
    /// Min melt coeff (winter solstice)
    pub dhmin: Vec<f64>,
    /// Max melt coeff (summer solstice)
    pub dhmax: Vec<f64>,
    /// Free water capacity fraction
    pub fwfrac: Vec<f64>,

    /// Fraction of total area for each subarea
    pub fArea: Vec<f64>,

    /// S2b — water age of the pack's water (seconds), complete-mix over wsnow + fw.
    pub age: Vec<f64>,
    /// pending DD-to-UEB plow heat, kJ/m2 of receiving subarea
    pub ueb_transfer_energy: Vec<f64>,
    /// driver ownership mask, true by default
    pub enabled: Vec<u8>,
    /// completed source footprints use actual donor/receiver areas
    pub conservative_plowing: bool,
    /// S2b — age of the water this step published as imelt (seconds).
    pub out_age: Vec<f64>,

    /// Snow depth for 100% cover (ft)
    pub si: Vec<f64>,
    /// Snow coverage area at start of new-snow ADC
    pub sba: Vec<f64>,
    /// Snow water equiv at end of new-snow ADC
    pub sbws: Vec<f64>,
    /// Computed areal snow coverage (0-1), per subarea
    pub asc: Vec<f64>,

    // Global ADC curves (shared by all subcatchments, 10 points each)
    pub adc_imperv: [f64; 10],
    pub adc_perv: [f64; 10],

    /// Plowable fraction of impervious area
    pub snn: Vec<f64>,

    /// Depth at which plowing begins (ft)
    pub weplow: Vec<f64>,
    /// Plowing fractions [subcatch * 5 + i]
    pub sfrac: Vec<f64>,
    /// Target subcatchment for plowed snow
    pub to_subcatch: Vec<i32>,

    /// ATI weighting factor
    pub tipm: f64,
    /// Negative melt ratio
    pub rnm: f64,
    /// Current snowmelt season factor (-1 to +1)
    pub season: f64,
    /// Cumulative snow plowed out of system (ft3)
    pub removed: f64,
    /// S2b — age of water arriving from the sky (seconds).
    pub precip_age: f64,
    /// S2b — whether to maintain age / out_age at all.
    pub track_age: bool,
}

impl Default for SnowSoA {
    fn default() -> Self {
        SnowSoA {
            n_subcatch: 0,
            external: Vec::new(),
            external_net: Vec::new(),
            wsnow: Vec::new(),
            fw: Vec::new(),
            coldc: Vec::new(),
            ati: Vec::new(),
            awe: Vec::new(),
            imelt: Vec::new(),
            tbase: Vec::new(),
            dhm: Vec::new(),
            dhmin: Vec::new(),
            dhmax: Vec::new(),
            fwfrac: Vec::new(),
            fArea: Vec::new(),
            age: Vec::new(),
            ueb_transfer_energy: Vec::new(),
            enabled: Vec::new(),
            conservative_plowing: false,
            out_age: Vec::new(),
            si: Vec::new(),
            sba: Vec::new(),
            sbws: Vec::new(),
            asc: Vec::new(),
            adc_imperv: [1.0; 10],
            adc_perv: [1.0; 10],
            snn: Vec::new(),
            weplow: Vec::new(),
            sfrac: Vec::new(),
            to_subcatch: Vec::new(),
            tipm: 0.5,
            rnm: 0.6,
            season: 0.0,
            removed: 0.0,
            precip_age: 0.0,
            track_age: false,
        }
    }
}

fn assign<T: Clone>(v: &mut Vec<T>, n: usize, d: T) {
    v.clear();
    v.resize(n, d);
}

impl SnowSoA {
    pub fn resize(&mut self, n: i32) {
        self.n_subcatch = n;
        assign(&mut self.enabled, n as usize, 1);
        assign(&mut self.external, n as usize, false);
        let total = (n * N_SUBAREAS) as usize;
        assign(&mut self.external_net, total, 0.0);
        assign(&mut self.ueb_transfer_energy, total, 0.0);
        let un = n as usize;

        assign(&mut self.wsnow, total, 0.0);
        assign(&mut self.fw, total, 0.0);
        assign(&mut self.coldc, total, 0.0);
        assign(&mut self.ati, total, 32.0);
        // 1.0, matching legacy snow_initSnowpack (snow.c:199).
        assign(&mut self.awe, total, 1.0);
        assign(&mut self.imelt, total, 0.0);
        assign(&mut self.tbase, total, 32.0);
        assign(&mut self.dhm, total, 0.0);
        assign(&mut self.dhmin, total, 0.0);
        assign(&mut self.dhmax, total, 0.0);
        assign(&mut self.fwfrac, total, 0.1);
        assign(&mut self.fArea, total, 0.0);
        assign(&mut self.age, total, 0.0);
        assign(&mut self.out_age, total, 0.0);
        assign(&mut self.si, total, 0.0);
        assign(&mut self.sba, total, 0.0);
        assign(&mut self.sbws, total, 0.0);
        assign(&mut self.asc, total, 1.0);
        assign(&mut self.snn, un, 0.0);
        assign(&mut self.weplow, un, 0.0);
        assign(&mut self.sfrac, (n * 5) as usize, 0.0);
        assign(&mut self.to_subcatch, un, -1);
    }
}

/// The degree-day snow solver (snow::SnowSolver is an alias of it).
#[derive(Debug, Clone, Default)]
pub struct DegreeDaySnow {
    pub soa_: SnowSoA,
}

pub type SnowSolver = DegreeDaySnow;

impl DegreeDaySnow {
    pub fn init(&mut self, n_subcatch: i32) {
        self.soa_.resize(n_subcatch);
    }

    /// SnowModel::initializeElements: init((n + N_SUBAREAS - 1) / N_SUBAREAS).
    pub fn initialize_elements(&mut self, n: i32) {
        self.init((n + N_SUBAREAS - 1) / N_SUBAREAS);
    }

    pub fn state_size(&self) -> i32 {
        12
    }

    pub fn state(&self) -> &SnowSoA {
        &self.soa_
    }

    pub fn state_mut(&mut self) -> &mut SnowSoA {
        &mut self.soa_
    }
}

// DegreeDaySnow::setMeltCoeffs (translated 5 November 2026 for the runoff step).
impl DegreeDaySnow {
    /// Seasonal melt coefficient interpolation (legacy snow.c snow_setMeltCoeffs):
    /// season = sin(0.0172615 * (day - 81)), period 364.000 days.
    pub fn set_melt_coeffs(&mut self, day_of_year: i32) {
        const K_SEASON_RAD: f64 = 0.0172615;
        let season = (K_SEASON_RAD * (day_of_year as f64 - 81.0)).sin(); // LIBM
        self.soa_.season = season;

        let total = self.soa_.n_subcatch * N_SUBAREAS;
        for i in 0..total as usize {
            self.soa_.dhm[i] = 0.5 * (self.soa_.dhmax[i] * (1.0 + season) + self.soa_.dhmin[i] * (1.0 - season));
        }
    }
}
