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
//   * Translated 5 November 2026 for the runoff step: setMeltCoeffs, the per-subcatchment
//     execute and plowSnow overloads (the ones stepRunoff calls), meltElement,
//     rainMeltRate and the file-static helpers (mixAge, getArealSnowCover,
//     getArealDepletion, updateColdContent, reduceColdContent, meltSnowpack,
//     meltThinPack, whose TrackAge template parameter is a bool argument). meltThinPack's
//     `#pragma omp simd` loop has no cross-lane operation; it runs as a plain loop.
//   * NOT YET TRANSLATED: step, pack, unpack and the scalar (broadcast) execute and
//     plowSnow overloads (not on the stepRunoff path).

#![allow(non_snake_case)]

use crate::core::cmath::{cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::unit_conversion as ucf;

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

// ============================================================================
// Run-time snow (DegreeDaySnow.cpp), translated 5 November 2026 for the runoff step
// ============================================================================

#[inline]
fn mix_age(v_have: f64, a_have: f64, v_in: f64, a_in: f64) -> f64 {
    let v = v_have + v_in;
    if !(v > 0.0) {
        return 0.0;
    }
    (v_have * a_have + v_in * a_in) / v
}

/// legacy getArealSnowCover (snow.c:691)
fn get_areal_snow_cover(adc: &[f64; 10], awesi: f64) -> f64 {
    if awesi <= 0.0 {
        return 0.0;
    }
    if awesi >= 0.9999 {
        return 1.0;
    }
    let m = (awesi * 10.0 + 0.00001) as i32;
    let asc1 = adc[m as usize];
    let asc2 = if m >= 9 { 1.0 } else { adc[(m + 1) as usize] };
    asc1 + (asc2 - asc1) / 0.1 * (awesi - 0.1 * (m as f32 as f64))
}

/// legacy getArealDepletion (snow.c:610)
fn get_areal_depletion(soa: &mut SnowSoA, ui: usize, subarea: i32, snowfall: f64, dt: f64) -> f64 {
    if subarea == SNOW_PLOWABLE {
        return 1.0;
    }
    let si = soa.si[ui];
    if si == 0.0 || soa.wsnow[ui] >= si {
        soa.awe[ui] = 1.0;
        return 1.0;
    }
    if soa.wsnow[ui] == 0.0 {
        soa.awe[ui] = 1.0;
        return 0.0;
    }
    let adc = if subarea == SNOW_PERV { soa.adc_perv } else { soa.adc_imperv };

    if snowfall > 0.0 {
        let mut awe = (soa.wsnow[ui] - snowfall * dt) / si;
        awe = cmax(awe, 0.0);
        let sba = get_areal_snow_cover(&adc, awe);
        let mut sbws = awe + (0.75 * snowfall * dt) / si;
        sbws = cmin(sbws, 1.0);
        soa.awe[ui] = awe;
        soa.sba[ui] = sba;
        soa.sbws[ui] = sbws;
        return 1.0;
    }

    let awe = soa.awe[ui];
    let sba = soa.sba[ui];
    let sbws = soa.sbws[ui];
    let awesi = soa.wsnow[ui] / si;
    if awesi < soa.awe[ui] {
        soa.awe[ui] = 1.0;
        get_areal_snow_cover(&adc, awesi)
    } else if awesi >= soa.sbws[ui] {
        1.0
    } else {
        sba + (1.0 - sba) / (sbws - awe) * (awesi - awe)
    }
}

/// legacy updateColdContent (snow.c:801)
fn update_cold_content(soa: &mut SnowSoA, ui: usize, temp: f64, asc: f64, snowfall: f64, dt: f64) {
    let mut ati = soa.ati[ui];
    let mut cc = soa.coldc[ui];
    if snowfall * 43200.0 > 0.02 {
        ati = temp;
    } else {
        let tipm = 1.0 - (1.0 - soa.tipm).powf(dt / (6.0 * 3600.0)); // LIBM
        ati += tipm * (temp - ati);
    }
    ati = cmin(ati, soa.tbase[ui]);
    cc += soa.rnm * soa.dhm[ui] * (ati - temp) * dt * asc;
    cc = cmax(cc, 0.0);
    let cc_max = soa.wsnow[ui] * 0.007 / 12.0 * (soa.tbase[ui] - ati);
    cc = cmin(cc, cc_max);
    soa.coldc[ui] = cc;
    soa.ati[ui] = ati;
}

/// legacy reduceColdContent (snow.c:854)
fn reduce_cold_content(soa: &mut SnowSoA, ui: usize, mut smelt: f64, cc_factor: f64) -> f64 {
    let mut cc = soa.coldc[ui];
    if smelt * cc_factor > cc {
        smelt -= cc / cc_factor;
        cc = 0.0;
    } else {
        cc -= smelt * cc_factor;
        smelt = 0.0;
    }
    soa.coldc[ui] = cc;
    smelt
}

/// legacy meltSnowpack (snow.c:729)
fn melt_snowpack(soa: &mut SnowSoA, ui: usize, temp: f64, rmelt: f64, asc: f64, snowfall: f64, dt: f64) -> f64 {
    let mut smelt;
    if rmelt > 0.0 {
        smelt = rmelt;
    } else if temp >= soa.tbase[ui] {
        smelt = soa.dhm[ui] * (temp - soa.tbase[ui]);
    } else {
        update_cold_content(soa, ui, temp, asc, snowfall, dt);
        return 0.0;
    }
    smelt *= asc;
    let cc_factor = dt * soa.rnm * asc;
    smelt = reduce_cold_content(soa, ui, smelt, cc_factor);
    soa.ati[ui] = soa.tbase[ui];
    smelt
}

/// meltThinPack<TrackAge>: every surface of an all-thin pack melts at once.
fn melt_thin_pack(soa: &mut SnowSoA, first: usize, dt: f64, track_age: bool) {
    for i in 0..N_SUBAREAS as usize {
        let k = first + i;
        let mut imelt = soa.imelt[k];
        imelt += (soa.wsnow[k] + soa.fw[k]) / dt;
        if track_age {
            soa.out_age[k] = soa.age[k];
            soa.age[k] = 0.0;
        }
        soa.wsnow[k] = 0.0;
        soa.fw[k] = 0.0;
        soa.coldc[k] = 0.0;
        soa.asc[k] = 0.0;
        soa.imelt[k] = 0.0 + imelt;
        if track_age && !(soa.imelt[k] > 0.0) {
            soa.out_age[k] = 0.0;
        }
    }
}

impl DegreeDaySnow {
    /// rainMeltRate (legacy getRainmelt, snow.c:775); can be negative below 32 F.
    pub fn rain_melt_rate(temp: f64, wind: f64, gamma: f64, ea: f64, rainfall: f64) -> f64 {
        let rainfall = rainfall * 43200.0;
        if rainfall > 0.02 {
            let uadj = 0.006 * wind;
            let t1 = temp - 32.0;
            let t2 = 7.5 * gamma * uadj;
            let t3 = 8.5 * uadj * (ea - 0.18);
            let smelt = t1 * (0.001167 + t2 + 0.007 * rainfall) + t3;
            return smelt / 43200.0;
        }
        0.0
    }

    /// execute (per-subcatchment rain and snow arrays): legacy snow_getSnowMelt.
    #[allow(clippy::too_many_arguments)]
    pub fn execute(&mut self, _ctx: &SimulationContext, dt: f64, temp: f64, wind: f64, rainfall: &[f64], snowfall: &[f64], gamma: f64, ea: f64) {
        let n = self.soa_.n_subcatch;
        if n == 0 {
            return;
        }
        for j in 0..n {
            let uj = j as usize;
            if self.soa_.enabled[uj] == 0 || (uj < self.soa_.external.len() && self.soa_.external[uj]) {
                continue;
            }
            let rain = rainfall[uj];
            let snow = snowfall[uj];
            let rmelt = Self::rain_melt_rate(temp, wind, gamma, ea, rain);
            let first = (j * N_SUBAREAS) as usize;
            if self.soa_.wsnow[first] <= 0.001 / 12.0 && self.soa_.wsnow[first + 1] <= 0.001 / 12.0 && self.soa_.wsnow[first + 2] <= 0.001 / 12.0 {
                let ta = self.soa_.track_age;
                melt_thin_pack(&mut self.soa_, first, dt, ta);
                continue;
            }
            for i in SNOW_PLOWABLE..=SNOW_PERV {
                let ui = (j * N_SUBAREAS + i) as usize;
                self.melt_element(ui, i, dt, temp, rain, snow, rmelt);
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn melt_element(&mut self, ui: usize, i: i32, dt: f64, temp: f64, rain: f64, snow: f64, rmelt: f64) {
        let soa = &mut self.soa_;
        let asc;
        let mut smelt;
        let mut imelt = soa.imelt[ui];

        if soa.wsnow[ui] <= 0.001 / 12.0 {
            asc = 0.0;
            smelt = 0.0;
            imelt += (soa.wsnow[ui] + soa.fw[ui]) / dt;
            if soa.track_age {
                soa.out_age[ui] = soa.age[ui];
                soa.age[ui] = 0.0;
            }
            soa.wsnow[ui] = 0.0;
            soa.fw[ui] = 0.0;
            soa.coldc[ui] = 0.0;
        } else {
            asc = get_areal_depletion(soa, ui, i, snow, dt);
            smelt = melt_snowpack(soa, ui, temp, rmelt, asc, snow, dt);

            let mut vmelt = smelt * dt;
            vmelt = cmin(vmelt, soa.wsnow[ui]);
            soa.wsnow[ui] -= vmelt;
            let rain_on_snow = rain * dt * asc;
            if soa.track_age && rain_on_snow > 0.0 {
                soa.age[ui] = mix_age(soa.wsnow[ui] + soa.fw[ui], soa.age[ui], rain_on_snow, soa.precip_age);
            }
            soa.fw[ui] += vmelt + rain_on_snow;
            vmelt = soa.fw[ui] - soa.fwfrac[ui] * soa.wsnow[ui];
            vmelt = cmax(vmelt, 0.0);
            soa.fw[ui] -= vmelt;
            smelt = vmelt / dt;
            if soa.track_age {
                soa.out_age[ui] = soa.age[ui];
            }
        }

        soa.asc[ui] = asc;
        soa.imelt[ui] = smelt + imelt;
        if soa.track_age && !(soa.imelt[ui] > 0.0) {
            soa.out_age[ui] = 0.0;
        }
    }

    /// plowSnow (per-subcatchment snowfall): accumulation, then plowing (legacy
    /// snow_plowSnow).
    pub fn plow_snow(&mut self, ctx: &SimulationContext, dt: f64, snowfall: &[f64]) {
        let n = self.soa_.n_subcatch;
        if n == 0 {
            return;
        }
        let soa = &mut self.soa_;

        for j in 0..n {
            let uj = j as usize;
            if soa.enabled[uj] == 0 || (uj < soa.external.len() && soa.external[uj]) {
                continue;
            }
            for k in SNOW_PLOWABLE..=SNOW_PERV {
                let idx = (j * N_SUBAREAS + k) as usize;
                if soa.fArea[idx] > 0.0 {
                    let add = snowfall[uj] * dt;
                    if soa.track_age {
                        let have = soa.wsnow[idx] + soa.fw[idx];
                        if have > 0.0 {
                            soa.age[idx] += dt;
                        }
                        if add > 0.0 {
                            soa.age[idx] = mix_age(have, soa.age[idx], add, soa.precip_age);
                        }
                    }
                    soa.wsnow[idx] += add;
                    soa.imelt[idx] = 0.0;
                }
            }
        }

        let ucf_land = ucf::ucf(ucf::LANDAREA, &ctx.options);
        for j in 0..n {
            let uj = j as usize;
            if soa.enabled[uj] == 0 || (uj < soa.external.len() && soa.external[uj]) {
                continue;
            }
            let plow_idx = (j * N_SUBAREAS + SNOW_PLOWABLE) as usize;
            if soa.fArea[plow_idx] <= 0.0 {
                continue;
            }
            if soa.wsnow[plow_idx] < soa.weplow[uj] {
                continue;
            }

            let exc = soa.wsnow[plow_idx];
            let sf = (j * 5) as usize;
            let mut sfrac_total = 0.0;
            let donor_age = if soa.track_age { soa.age[plow_idx] } else { 0.0 };

            if soa.sfrac[sf] > 0.0 {
                let lid_ft2 = if uj < ctx.subcatches.total_lid_area_ft2.len() { ctx.subcatches.total_lid_area_ft2[uj] } else { 0.0 };
                let mut area_ft2 = ctx.subcatches.area[uj] / ucf_land - lid_ft2;
                if area_ft2 < 0.0 {
                    area_ft2 = 0.0;
                }
                soa.removed += soa.sfrac[sf] * exc * soa.fArea[plow_idx] * area_ft2;
            }
            sfrac_total += soa.sfrac[sf];

            let imperv_idx = (j * N_SUBAREAS + SNOW_IMPERV) as usize;
            if soa.fArea[imperv_idx] > 0.0 {
                let f = soa.fArea[plow_idx] / soa.fArea[imperv_idx];
                let moved = soa.sfrac[sf + 1] * exc * f;
                if soa.track_age && moved > 0.0 {
                    soa.age[imperv_idx] = mix_age(soa.wsnow[imperv_idx] + soa.fw[imperv_idx], soa.age[imperv_idx], moved, donor_age);
                }
                soa.wsnow[imperv_idx] += moved;
                sfrac_total += soa.sfrac[sf + 1];
            }

            let perv_idx = (j * N_SUBAREAS + SNOW_PERV) as usize;
            if soa.fArea[perv_idx] > 0.0 {
                let f = soa.fArea[plow_idx] / soa.fArea[perv_idx];
                let moved = soa.sfrac[sf + 2] * exc * f;
                if soa.track_age && moved > 0.0 {
                    soa.age[perv_idx] = mix_age(soa.wsnow[perv_idx] + soa.fw[perv_idx], soa.age[perv_idx], moved, donor_age);
                }
                soa.wsnow[perv_idx] += moved;
                sfrac_total += soa.sfrac[sf + 2];
            }

            if dt > 0.0 {
                soa.imelt[plow_idx] = soa.sfrac[sf + 3] * exc / dt;
                if soa.track_age {
                    soa.out_age[plow_idx] = donor_age;
                }
            }
            sfrac_total += soa.sfrac[sf + 3];

            if soa.sfrac[sf + 4] > 0.0 && soa.to_subcatch[uj] >= 0 {
                let m = soa.to_subcatch[uj];
                if m < n {
                    let target_perv = (m * N_SUBAREAS + SNOW_PERV) as usize;
                    if soa.fArea[target_perv] > 0.0 {
                        let mut f = soa.fArea[plow_idx] / soa.fArea[target_perv];
                        if soa.conservative_plowing || ((m as usize) < soa.external.len() && soa.external[m as usize]) {
                            let actual_area = |sc: i32| -> f64 {
                                let lid = if (sc as usize) < ctx.subcatches.total_lid_area_ft2.len() { ctx.subcatches.total_lid_area_ft2[sc as usize] } else { 0.0 };
                                cmax(0., ctx.subcatches.area[sc as usize] / ucf::ucf(ucf::LANDAREA, &ctx.options) - lid)
                            };
                            let target_area = actual_area(m);
                            f *= if target_area > 0.0 { actual_area(j) / target_area } else { 0.0 };
                        }
                        let moved = soa.sfrac[sf + 4] * exc * f;
                        if soa.track_age && moved > 0.0 {
                            soa.age[target_perv] = mix_age(soa.wsnow[target_perv] + soa.fw[target_perv], soa.age[target_perv], moved, donor_age);
                        }
                        if soa.external[m as usize] && moved > 0.0 && exc > 0.0 {
                            let cold = soa.coldc[plow_idx] * soa.sfrac[sf + 4];
                            soa.ueb_transfer_energy[target_perv] -= cold * f * 333500.0 * 0.3048;
                            soa.coldc[plow_idx] -= cold;
                        }
                        soa.wsnow[target_perv] += moved;
                        if f > 0.0 {
                            sfrac_total += soa.sfrac[sf + 4];
                        }
                    }
                }
            }

            sfrac_total = cmin(sfrac_total, 1.0);
            soa.wsnow[plow_idx] = exc * (1.0 - sfrac_total);
            if soa.track_age && !(soa.wsnow[plow_idx] + soa.fw[plow_idx] > 0.0) {
                soa.age[plow_idx] = 0.0;
            }
        }
    }
}
