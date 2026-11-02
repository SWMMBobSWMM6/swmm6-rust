// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/surface/NativeSnow.cpp
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright (c) 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 2 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The adapter between the parsed [SNOWPACKS] data and the degree-day snow solver:
// initializeNativeSnow (called by SWMMEngine::initHydrology), hasUebSnow and the
// depth/volume helpers. The error is a std::string upstream (empty = success); here
// Bytes, empty on success.
//
// NOT YET TRANSLATED (run-time report helper, phase 1c): nativeSnowResults (it needs
// ueb::bulkTemperature).

use crate::core::cmath::cmin;
use crate::core::context::SimulationContext;
use crate::core::options::Bytes;
use crate::core::unit_conversion as ucf;
use crate::hydrology::snow::{SnowSoA, SnowSolver, N_SUBAREAS, SNOW_IMPERV, SNOW_PERV, SNOW_PLOWABLE};
use crate::hydrology::ueb_subcatchments::initialize_subcatchments;
use crate::hydrology::ueb_types::SnowKind;

pub fn has_ueb_snow(c: &SimulationContext) -> bool {
    if c.options.ignore_snow_melt {
        return false;
    }
    for &pack in &c.subcatches.snowpack {
        if pack >= 0 && c.snowpacks.config.model_for(pack) == SnowKind::Ueb {
            return true;
        }
    }
    false
}

pub fn native_snow_depth(s: &SnowSoA, sc: i32) -> f64 {
    let mut water = 0.0;
    for k in 0..3 {
        let i = (sc * 3 + k) as usize;
        water += s.wsnow[i] * s.fArea[i];
    }
    water
}

pub fn native_snow_volume(s: &SnowSoA, sc: i32, area: f64) -> f64 {
    let mut water = 0.0;
    for k in 0..3 {
        let i = (sc * 3 + k) as usize;
        water += (s.wsnow[i] + s.fw[i]) * s.fArea[i];
    }
    water * area
}

pub fn initialize_native_snow(ctx_: &mut SimulationContext, snow_: &mut SnowSolver) -> Bytes {
    snow_.init(ctx_.n_subcatches());
    {
        let soa = snow_.state_mut();
        let unit_sys_snow = ucf::get_unit_system(ctx_.options.flow_units as i32) as usize;
        let u = |q: i32| ucf::UCF[q as usize][unit_sys_snow];

        // Transfer global snowmelt parameters
        soa.tipm = ctx_.options.snow_ati_wt;
        soa.rnm = ctx_.options.snow_nrg_ratio;

        // Transfer ADC curves from options
        for k in 0..10 {
            soa.adc_imperv[k] = ctx_.options.adc_imperv[k];
            soa.adc_perv[k] = ctx_.options.adc_perv[k];
        }

        // Transfer per-subcatchment snowpack parameters
        for i in 0..ctx_.n_subcatches() {
            let ui = i as usize;
            let sp_idx = ctx_.subcatches.snowpack[ui];
            if sp_idx < 0 {
                continue;
            }
            let usp = sp_idx as usize;
            let sp = &ctx_.snowpacks;

            // Surface types: 0=PLOWABLE, 1=IMPERVIOUS, 2=PERVIOUS
            let surfaces: [Option<&[f64; 7]>; 3] = [
                if usp < sp.plowable.len() { Some(&sp.plowable[usp]) } else { None },
                if usp < sp.impervious.len() { Some(&sp.impervious[usp]) } else { None },
                if usp < sp.pervious.len() { Some(&sp.pervious[usp]) } else { None },
            ];

            for k in 0..N_SUBAREAS {
                let idx = (i * N_SUBAREAS + k) as usize;
                let p = match surfaces[k as usize] {
                    None => continue,
                    Some(p) => p,
                };
                // p[0]=Cmin, p[1]=Cmax, p[2]=Tbase, p[3]=FWF, p[4]=SD0, p[5]=FW0, p[6]=SNN0
                soa.dhmin[idx] = p[0] * u(ucf::TEMPERATURE) / u(ucf::RAINFALL);
                soa.dhmax[idx] = p[1] * u(ucf::TEMPERATURE) / u(ucf::RAINFALL);
                soa.tbase[idx] = p[2];
                if unit_sys_snow == 1 {
                    // SI: convert degC to degF
                    soa.tbase[idx] = (9.0 / 5.0) * soa.tbase[idx] + 32.0;
                }
                soa.fwfrac[idx] = p[3];
                soa.wsnow[idx] = p[4] / u(ucf::RAINDEPTH);
                soa.fw[idx] = cmin(p[5], p[3] * p[4]) / u(ucf::RAINDEPTH);
                soa.ati[idx] = soa.tbase[idx];
                // F6: the 7th field is SNN0 on the PLOWABLE row and SD100 on the other two
                // (snow.c:348-352).
                if k == SNOW_PLOWABLE {
                    soa.snn[ui] = p[6];
                } else {
                    soa.si[idx] = p[6] / u(ucf::RAINDEPTH);
                }
            }

            // Fractional area of each snow surface (legacy snow_initSnowpack,
            // snow.c:178-182).
            {
                let fimp = ctx_.subcatches.frac_imperv[ui];
                let snn = soa.snn[ui];
                let base = (i * N_SUBAREAS) as usize;
                soa.fArea[base + SNOW_PLOWABLE as usize] = snn * fimp;
                soa.fArea[base + SNOW_IMPERV as usize] = (1.0 - snn) * fimp;
                soa.fArea[base + SNOW_PERV as usize] = 1.0 - fimp;
                // Zero initial state where the surface has no area (snow.c:184-197).
                for k in 0..N_SUBAREAS {
                    let idx = base + k as usize;
                    if soa.fArea[idx] <= 0.0 {
                        soa.wsnow[idx] = 0.0;
                        soa.fw[idx] = 0.0;
                        soa.si[idx] = 0.0;
                    }
                }
            }

            // Transfer plowing/removal parameters
            if usp < sp.removal.len() {
                let r = &sp.removal[usp];
                let sf = (i * 5) as usize;
                soa.weplow[ui] = r[0] / u(ucf::RAINDEPTH);
                soa.sfrac[sf] = r[1]; // fraction removed
                soa.sfrac[sf + 1] = r[2]; // fraction to impervious
                soa.sfrac[sf + 2] = r[3]; // fraction to pervious
                soa.sfrac[sf + 3] = r[4]; // fraction immediate melt
                soa.sfrac[sf + 4] = r[5]; // fraction to other subcatchment

                // Resolve target subcatchment for plowed snow
                if usp < sp.removal_subcatch.len() && !sp.removal_subcatch[usp].is_empty() {
                    soa.to_subcatch[ui] = ctx_.subcatch_names.find(&sp.removal_subcatch[usp]);
                }
            }
        }
    }
    let error = initialize_subcatchments(ctx_, snow_.state_mut());
    if !error.is_empty() {
        return error;
    }
    for sc in 0..ctx_.n_subcatches() {
        ctx_.subcatches.snow_depth[sc as usize] = native_snow_depth(snow_.state(), sc);
    }
    Bytes::new()
}
