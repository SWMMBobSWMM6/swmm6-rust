// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/hydrology/LID.cpp (HydroCouple/openswmm.engine,
// branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the
// Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The LID process kernel (LID.cpp's file-static LidProcView, LegacyLidKernel, makeView,
// runUnitLegacy: legacy lidproc.c) and LIDSolver::execute, the call stepRunoff makes.
//
// Differences of form, never in result:
//   * LegacyLidKernel holds its unit through the pointers P (the parameter view), G (the
//     group) and U (the unit). Here the view is a field and the group is passed to each
//     method that reads or writes it.
//   * The pointer out-parameters of getOutflow (lidEvap, lidInfil, lidDrain) and of
//     getSurfaceOverflowRate (surfaceDepth) are &mut.
//   * saveLidReport writes through std::ofstream; here through std::fs::File, unbuffered
//     (the file's final content is the same). printf's %8.3f / %8.4f are fmt_f (Rust's
//     exact decimal formatting, rounding ties to even as glibc does; NaN printed as glibc
//     prints it). The report file's content is not part of any gate.

use crate::core::cmath::{cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::datetime;
use crate::hydrology::lid::LidRptFile;
use std::io::Write;
use crate::hydrology::infiltration::{self as infil, InfilModel};
use crate::hydrology::lid::{has_native_bottom, LIDGroupSoA, LIDSolver, LIDType};

pub const L_SURF: usize = 0;
pub const L_SOIL: usize = 1;
pub const L_STOR: usize = 2;
pub const L_PAVE: usize = 3;
pub const L_MAX: usize = 4;
pub const L_BIG: f64 = 1.0e10;
pub const L_ZERO: f64 = 1.0e-10;
pub const L_STOPTOL: f64 = 0.00328;
pub const L_MINFLOW: f64 = 2.3e-8;

/// Per-unit parameter view (ft, ft/s; legacy TLidProc after validateLidProc).
#[derive(Debug, Clone, Copy)]
pub struct LidProcView {
    pub type_: LIDType,
    pub surf_thick: f64,
    pub surf_void: f64,
    pub surf_alpha: f64,
    pub surf_side: f64,
    pub can_overflow: bool,
    pub pave_thick: f64,
    pub pave_void: f64,
    pub pave_imperv: f64,
    pub pave_ksat: f64,
    pub pave_clog: f64,
    pub pave_regen_days: f64,
    pub pave_regen_deg: f64,
    pub soil_thick: f64,
    pub soil_poros: f64,
    pub soil_fc: f64,
    pub soil_wp: f64,
    pub soil_ksat: f64,
    pub soil_kslope: f64,
    pub stor_thick: f64,
    pub stor_void: f64,
    pub stor_ksat: f64,
    pub stor_clog: f64,
    pub drain_coeff: f64,
    pub drain_expon: f64,
    pub drain_offset: f64,
    pub drain_delay: f64,
    pub drain_hopen: f64,
    pub drain_hclose: f64,
    pub dm_thick: f64,
    pub dm_alpha: f64,
    pub area: f64,
    pub full_width: f64,
    pub ucf_rainfall: f64,
    pub ucf_raindepth: f64,
}

/// legacy lidproc.c's file-scope variables and process functions for one unit.
#[derive(Debug, Clone, Copy)]
struct LegacyLidKernel {
    tstep: f64,
    evap_rate: f64,
    max_native_infil: f64,
    spatial_fraction: f64,
    spatial_limit: f64,
    spatial_exfil_ratio: f64,
    completed: bool,
    surface_inflow: f64,
    surface_infil: f64,
    surface_evap: f64,
    surface_outflow: f64,
    pave_evap: f64,
    pave_perc: f64,
    soil_evap: f64,
    soil_perc: f64,
    storage_inflow: f64,
    storage_exfil: f64,
    storage_evap: f64,
    storage_drain: f64,
    surface_volume: f64,
    pave_volume: f64,
    soil_volume: f64,
    storage_volume: f64,
    p: LidProcView,
    u: usize,
    old_runoff_days: f64,
}

impl LegacyLidKernel {
    fn new(p: LidProcView, u: usize) -> Self {
        LegacyLidKernel {
            tstep: 0.0,
            evap_rate: 0.0,
            max_native_infil: L_BIG,
            spatial_fraction: -1.0,
            spatial_limit: L_BIG,
            spatial_exfil_ratio: 0.0,
            completed: false,
            surface_inflow: 0.0,
            surface_infil: 0.0,
            surface_evap: 0.0,
            surface_outflow: 0.0,
            pave_evap: 0.0,
            pave_perc: 0.0,
            soil_evap: 0.0,
            soil_perc: 0.0,
            storage_inflow: 0.0,
            storage_exfil: 0.0,
            storage_evap: 0.0,
            storage_drain: 0.0,
            surface_volume: 0.0,
            pave_volume: 0.0,
            soil_volume: 0.0,
            storage_volume: 0.0,
            p,
            u,
            old_runoff_days: 0.0,
        }
    }

    // legacy getSurfaceOutflowRate
    fn get_surface_outflow_rate(&self, depth: f64) -> f64 {
        let delta = depth - self.p.surf_thick;
        if delta < 0.0 {
            return 0.0;
        }
        let mut outflow = self.p.surf_alpha * delta.powf(5.0 / 3.0) * self.p.full_width / self.p.area; // LIBM
        outflow = cmin(outflow, delta / self.tstep);
        if self.completed {
            outflow = cmin(outflow, cmax(0.0, self.surface_inflow + self.surface_volume / self.tstep - self.surface_evap - self.surface_infil));
        }
        outflow
    }

    // legacy getSurfaceOverflowRate (clamps the depth it is handed)
    fn get_surface_overflow_rate(&self, surface_depth: &mut f64) -> f64 {
        let delta = *surface_depth - self.p.surf_thick;
        if delta <= 0.0 {
            return 0.0;
        }
        *surface_depth = self.p.surf_thick;
        delta * self.p.surf_void / self.tstep
    }

    // legacy getPavementPermRate
    fn get_pavement_perm_rate(&mut self, g: &mut LIDGroupSoA) -> f64 {
        let u = self.u;
        let mut perm_reduction = 0.0;
        let clog_factor = self.p.pave_clog;
        let regen_days = self.p.pave_regen_days;
        if clog_factor > 0.0 {
            if regen_days > 0.0 && self.old_runoff_days >= g.next_regen_day[u] {
                g.vol_treated[u] *= 1.0 - self.p.pave_regen_deg;
                g.next_regen_day[u] += regen_days;
            }
            perm_reduction = g.vol_treated[u] / clog_factor;
            perm_reduction = cmin(perm_reduction, 1.0);
        }
        self.p.pave_ksat * (1.0 - perm_reduction)
    }

    // legacy getSoilPercRate
    fn get_soil_perc_rate(&self, theta: f64) -> f64 {
        if theta <= self.p.soil_fc {
            return 0.0;
        }
        let delta = self.p.soil_poros - theta;
        self.p.soil_ksat * (-delta * self.p.soil_kslope).exp() // LIBM
    }

    // legacy getStorageExfilRate
    fn get_storage_exfil_rate(&mut self, g: &LIDGroupSoA) -> f64 {
        if self.spatial_fraction < 0.0 && self.p.stor_ksat == 0.0 {
            return 0.0;
        }
        if self.spatial_fraction < 0.0 && self.max_native_infil == 0.0 {
            return 0.0;
        }
        let mut clog_factor = self.p.stor_clog;
        if clog_factor > 0.0 {
            clog_factor = g.wb_inflow[self.u] / clog_factor;
            clog_factor = cmin(clog_factor, 1.0);
        }
        let infil = self.p.stor_ksat * (1.0 - clog_factor);
        let native = cmin(infil, self.max_native_infil);
        if self.spatial_fraction < 0.0 {
            return native;
        }
        let inside = if has_native_bottom(self.p.type_, self.p.stor_thick, self.p.stor_ksat) {
            if self.p.stor_thick <= 0.0 {
                self.spatial_limit
            } else {
                cmin(infil, self.spatial_limit)
            }
        } else {
            0.0
        };
        let total = self.spatial_fraction * inside + (1.0 - self.spatial_fraction) * native;
        self.spatial_exfil_ratio = if total > 0.0 { self.spatial_fraction * inside / total } else { 0.0 };
        total
    }

    // legacy getStorageDrainRate
    fn get_storage_drain_rate(&self, g: &LIDGroupSoA, storage_depth: f64, soil_theta: f64, pave_depth: f64, surface_depth: f64) -> f64 {
        let p = &self.p;
        let mut head = storage_depth;
        let mut outflow = 0.0;
        if storage_depth >= p.stor_thick {
            if p.soil_thick > 0.0 {
                if soil_theta > p.soil_fc {
                    head += (soil_theta - p.soil_fc) / (p.soil_poros - p.soil_fc) * p.soil_thick;
                    if soil_theta >= p.soil_poros {
                        if p.pave_thick > 0.0 {
                            head += pave_depth;
                            if pave_depth >= p.pave_thick {
                                head += surface_depth;
                            }
                        } else {
                            head += surface_depth;
                        }
                    }
                }
            } else if p.pave_thick > 0.0 {
                head += pave_depth;
                if pave_depth >= p.pave_thick {
                    head += surface_depth;
                }
            }
        }
        if g.old_drain_flow[self.u] == 0.0 && head <= p.drain_hopen {
            return 0.0;
        }
        if g.old_drain_flow[self.u] > 0.0 && head <= p.drain_hclose {
            return 0.0;
        }
        head -= p.drain_offset;
        if head > L_ZERO {
            head *= p.ucf_raindepth;
            outflow = p.drain_coeff * head.powf(p.drain_expon); // LIBM
            outflow /= p.ucf_rainfall;
        }
        outflow
    }

    // legacy getDrainMatOutflow
    fn get_drain_mat_outflow(&self, depth: f64) -> f64 {
        let mut result = self.soil_perc;
        if self.p.dm_alpha > 0.0 {
            result = self.p.dm_alpha * depth.powf(5.0 / 3.0) * self.p.full_width / self.p.area * self.p.stor_void; // LIBM
        }
        result
    }

    // legacy getEvapRates
    fn get_evap_rates(&mut self, surface_vol: f64, pave_vol: f64, soil_vol: f64, storage_vol: f64, perv_frac: f64) {
        let mut avail_evap = self.evap_rate;
        self.surface_evap = cmin(avail_evap, surface_vol / self.tstep);
        self.surface_evap = cmax(0.0, self.surface_evap);
        avail_evap = cmax(0.0, avail_evap - self.surface_evap);
        avail_evap *= perv_frac;
        if self.surface_infil > 0.0 {
            self.pave_evap = 0.0;
            self.soil_evap = 0.0;
            self.storage_evap = 0.0;
        } else {
            self.pave_evap = cmin(avail_evap, pave_vol / self.tstep);
            avail_evap = cmax(0.0, avail_evap - self.pave_evap);
            self.soil_evap = cmin(avail_evap, soil_vol / self.tstep);
            avail_evap = cmax(0.0, avail_evap - self.soil_evap);
            self.storage_evap = cmin(avail_evap, storage_vol / self.tstep);
        }
    }

    // legacy roofFluxRates
    fn roof_flux_rates(&mut self, x: &[f64; L_MAX], f: &mut [f64; L_MAX]) {
        let mut surface_depth = x[L_SURF];
        self.surface_volume = if self.completed { surface_depth * self.p.surf_void } else { surface_depth };
        self.get_evap_rates(self.surface_volume, 0.0, 0.0, 0.0, 1.0);
        self.surface_infil = 0.0;
        if self.p.surf_alpha > 0.0 {
            self.surface_outflow = self.get_surface_outflow_rate(surface_depth);
        } else {
            self.get_surface_overflow_rate(&mut surface_depth);
        }
        self.storage_drain = cmin(self.p.drain_coeff / self.p.ucf_rainfall, self.surface_outflow);
        self.surface_outflow -= self.storage_drain;
        f[L_SURF] = self.surface_inflow - self.surface_evap - self.storage_drain - self.surface_outflow;
        if self.completed {
            f[L_SURF] /= self.p.surf_void;
        }
    }

    // legacy greenRoofFluxRates
    fn green_roof_flux_rates(&mut self, x: &[f64; L_MAX], f: &mut [f64; L_MAX]) {
        let p = self.p;
        let soil_thickness = p.soil_thick;
        let storage_thickness = p.stor_thick;
        let soil_porosity = p.soil_poros;
        let storage_void_frac = p.stor_void;
        let soil_field_cap = p.soil_fc;
        let soil_wilt_point = p.soil_wp;
        let surface_depth = x[L_SURF];
        let soil_theta = x[L_SOIL];
        let storage_depth = x[L_STOR];
        self.surface_volume = surface_depth * p.surf_void;
        self.soil_volume = soil_theta * soil_thickness;
        self.storage_volume = storage_depth * storage_void_frac;
        let mut avail_volume = self.soil_volume - soil_wilt_point * soil_thickness;
        self.get_evap_rates(self.surface_volume, 0.0, avail_volume, self.storage_volume, 1.0);
        if self.completed {
            self.surface_infil = cmin(self.surface_infil, cmax(0.0, self.surface_inflow + self.surface_volume / self.tstep - self.surface_evap));
        }
        if soil_theta >= soil_porosity {
            self.storage_evap = 0.0;
        }
        self.soil_perc = self.get_soil_perc_rate(soil_theta);
        avail_volume = (soil_theta - soil_field_cap) * soil_thickness;
        let mut max_rate = cmax(avail_volume, 0.0) / self.tstep - self.soil_evap;
        self.soil_perc = cmin(self.soil_perc, max_rate);
        self.soil_perc = cmax(self.soil_perc, 0.0);
        self.storage_exfil = 0.0;
        self.storage_drain = self.get_drain_mat_outflow(storage_depth);
        if soil_theta >= soil_porosity && storage_depth >= storage_thickness {
            max_rate = cmin(self.soil_perc, self.storage_drain);
            self.soil_perc = max_rate;
            self.storage_drain = max_rate;
            self.surface_infil = cmin(self.surface_infil, max_rate);
        } else {
            max_rate = storage_depth * storage_void_frac / self.tstep - self.storage_evap;
            if storage_depth >= storage_thickness {
                max_rate += self.soil_perc;
            }
            max_rate = cmax(max_rate, 0.0);
            self.storage_drain = cmin(self.storage_drain, max_rate);
            max_rate = (storage_thickness - storage_depth) * storage_void_frac / self.tstep + self.storage_drain + self.storage_evap;
            self.soil_perc = cmin(self.soil_perc, max_rate);
            max_rate = (soil_porosity - soil_theta) * soil_thickness / self.tstep + self.soil_perc + self.soil_evap;
            self.surface_infil = cmin(self.surface_infil, max_rate);
        }
        self.surface_outflow = self.get_surface_outflow_rate(surface_depth);
        f[L_SURF] = (self.surface_inflow - self.surface_evap - self.surface_infil - self.surface_outflow) / p.surf_void;
        f[L_SOIL] = (self.surface_infil - self.soil_evap - self.soil_perc) / p.soil_thick;
        f[L_STOR] = (self.soil_perc - self.storage_evap - self.storage_drain) / p.stor_void;
    }

    // legacy biocellFluxRates
    fn biocell_flux_rates(&mut self, g: &LIDGroupSoA, x: &[f64; L_MAX], f: &mut [f64; L_MAX]) {
        let p = self.p;
        let soil_thickness = p.soil_thick;
        let soil_porosity = p.soil_poros;
        let soil_field_cap = p.soil_fc;
        let soil_wilt_point = p.soil_wp;
        let storage_thickness = p.stor_thick;
        let storage_void_frac = p.stor_void;
        let surface_depth = x[L_SURF];
        let soil_theta = x[L_SOIL];
        let storage_depth = x[L_STOR];
        self.surface_volume = surface_depth * p.surf_void;
        self.soil_volume = soil_theta * soil_thickness;
        self.storage_volume = storage_depth * storage_void_frac;
        let mut avail_volume = self.soil_volume - soil_wilt_point * soil_thickness;
        self.get_evap_rates(self.surface_volume, 0.0, avail_volume, self.storage_volume, 1.0);
        if self.completed {
            self.surface_infil = cmin(self.surface_infil, cmax(0.0, self.surface_inflow + self.surface_volume / self.tstep - self.surface_evap));
        }
        if soil_theta >= soil_porosity {
            self.storage_evap = 0.0;
        }
        self.soil_perc = self.get_soil_perc_rate(soil_theta);
        avail_volume = (soil_theta - soil_field_cap) * soil_thickness;
        let mut max_rate = cmax(avail_volume, 0.0) / self.tstep - self.soil_evap;
        self.soil_perc = cmin(self.soil_perc, max_rate);
        self.soil_perc = cmax(self.soil_perc, 0.0);
        self.storage_exfil = self.get_storage_exfil_rate(g);
        self.storage_drain = 0.0;
        if p.drain_coeff > 0.0 {
            self.storage_drain = self.get_storage_drain_rate(g, storage_depth, soil_theta, 0.0, surface_depth);
        }
        if storage_thickness == 0.0 {
            self.storage_evap = 0.0;
            max_rate = cmin(self.soil_perc, self.storage_exfil);
            self.soil_perc = max_rate;
            self.storage_exfil = max_rate;
            max_rate = (soil_porosity - soil_theta) * soil_thickness / self.tstep + self.soil_perc + self.soil_evap;
            self.surface_infil = cmin(self.surface_infil, max_rate);
        } else if soil_theta >= soil_porosity && storage_depth >= storage_thickness {
            max_rate = self.storage_exfil + self.storage_drain;
            if self.soil_perc < max_rate {
                max_rate = self.soil_perc;
                if max_rate > self.storage_exfil {
                    self.storage_drain = max_rate - self.storage_exfil;
                } else {
                    self.storage_exfil = max_rate;
                    self.storage_drain = 0.0;
                }
            } else {
                self.soil_perc = max_rate;
            }
            self.surface_infil = cmin(self.surface_infil, max_rate);
        } else {
            max_rate = self.soil_perc - self.storage_evap + storage_depth * storage_void_frac / self.tstep;
            self.storage_exfil = cmin(self.storage_exfil, max_rate);
            self.storage_exfil = cmax(self.storage_exfil, 0.0);
            if self.storage_drain > 0.0 {
                max_rate = -self.storage_exfil - self.storage_evap;
                if storage_depth >= storage_thickness {
                    max_rate += self.soil_perc;
                }
                if p.drain_offset <= storage_depth {
                    max_rate += (storage_depth - p.drain_offset) * storage_void_frac / self.tstep;
                }
                max_rate = cmax(max_rate, 0.0);
                self.storage_drain = cmin(self.storage_drain, max_rate);
            }
            max_rate = self.storage_exfil + self.storage_drain + self.storage_evap + (storage_thickness - storage_depth) * storage_void_frac / self.tstep;
            self.soil_perc = cmin(self.soil_perc, max_rate);
            max_rate = (soil_porosity - soil_theta) * soil_thickness / self.tstep + self.soil_perc + self.soil_evap;
            self.surface_infil = cmin(self.surface_infil, max_rate);
        }
        self.surface_outflow = self.get_surface_outflow_rate(surface_depth);
        f[L_SURF] = (self.surface_inflow - self.surface_evap - self.surface_infil - self.surface_outflow) / p.surf_void;
        f[L_SOIL] = (self.surface_infil - self.soil_evap - self.soil_perc) / p.soil_thick;
        if storage_thickness == 0.0 {
            f[L_STOR] = 0.0;
        } else {
            f[L_STOR] = (self.soil_perc - self.storage_evap - self.storage_exfil - self.storage_drain) / p.stor_void;
        }
    }

    // legacy trenchFluxRates
    fn trench_flux_rates(&mut self, g: &LIDGroupSoA, x: &[f64; L_MAX], f: &mut [f64; L_MAX]) {
        let p = self.p;
        let storage_thickness = p.stor_thick;
        let storage_void_frac = p.stor_void;
        let surface_depth = x[L_SURF];
        let storage_depth = x[L_STOR];
        self.surface_volume = surface_depth * p.surf_void;
        self.soil_volume = 0.0;
        self.storage_volume = storage_depth * storage_void_frac;
        self.get_evap_rates(self.surface_volume, 0.0, 0.0, self.storage_volume, 1.0);
        if surface_depth > 0.0 {
            self.storage_evap = 0.0;
        }
        self.storage_inflow = self.surface_inflow + self.surface_volume / self.tstep;
        if self.completed {
            self.storage_inflow = cmax(0.0, self.storage_inflow - self.surface_evap);
        }
        self.storage_exfil = self.get_storage_exfil_rate(g);
        self.storage_drain = 0.0;
        if p.drain_coeff > 0.0 {
            self.storage_drain = self.get_storage_drain_rate(g, storage_depth, 0.0, 0.0, surface_depth);
        }
        let mut max_rate = self.storage_inflow - self.storage_evap + storage_depth * storage_void_frac / self.tstep;
        self.storage_exfil = cmin(self.storage_exfil, max_rate);
        self.storage_exfil = cmax(self.storage_exfil, 0.0);
        if self.storage_drain > 0.0 {
            max_rate = -self.storage_exfil - self.storage_evap;
            if storage_depth >= storage_thickness {
                max_rate += self.storage_inflow;
            }
            if p.drain_offset <= storage_depth {
                max_rate += (storage_depth - p.drain_offset) * storage_void_frac / self.tstep;
            }
            max_rate = cmax(max_rate, 0.0);
            self.storage_drain = cmin(self.storage_drain, max_rate);
        }
        max_rate = (storage_thickness - storage_depth) * storage_void_frac / self.tstep + self.storage_exfil + self.storage_evap + self.storage_drain;
        self.storage_inflow = cmin(self.storage_inflow, max_rate);
        self.surface_infil = self.storage_inflow;
        self.surface_outflow = self.get_surface_outflow_rate(surface_depth);
        f[L_SURF] = (self.surface_inflow - self.surface_evap - self.storage_inflow - self.surface_outflow) / p.surf_void;
        f[L_STOR] = (self.storage_inflow - self.storage_evap - self.storage_exfil - self.storage_drain) / p.stor_void;
        f[L_SOIL] = 0.0;
    }

    // legacy pavementFluxRates
    fn pavement_flux_rates(&mut self, g: &mut LIDGroupSoA, x: &[f64; L_MAX], f: &mut [f64; L_MAX]) {
        let p = self.p;
        let perv_frac = 1.0 - p.pave_imperv;
        let storage_inflow;
        let mut avail_volume;
        let mut max_rate;
        let pave_void_frac = p.pave_void * perv_frac;
        let pave_thickness = p.pave_thick;
        let soil_thickness = p.soil_thick;
        let soil_porosity = p.soil_poros;
        let soil_field_cap = p.soil_fc;
        let soil_wilt_point = p.soil_wp;
        let storage_thickness = p.stor_thick;
        let storage_void_frac = p.stor_void;
        let surface_depth = x[L_SURF];
        let pave_depth = x[L_PAVE];
        let soil_theta = x[L_SOIL];
        let storage_depth = x[L_STOR];
        self.surface_volume = surface_depth * p.surf_void;
        self.pave_volume = pave_depth * pave_void_frac;
        self.soil_volume = soil_theta * soil_thickness;
        self.storage_volume = storage_depth * storage_void_frac;
        avail_volume = self.soil_volume - soil_wilt_point * soil_thickness;
        self.get_evap_rates(self.surface_volume, self.pave_volume, avail_volume, self.storage_volume, perv_frac);
        if pave_depth >= pave_thickness || (soil_thickness > 0.0 && soil_theta >= soil_porosity) {
            self.storage_evap = 0.0;
        }
        self.surface_infil = self.surface_inflow + (self.surface_volume / self.tstep);
        if self.completed {
            self.surface_infil = cmax(0.0, self.surface_infil - self.surface_evap);
        }
        self.pave_perc = self.get_pavement_perm_rate(g) * perv_frac;
        self.surface_infil = cmin(self.surface_infil, self.pave_perc);
        max_rate = self.pave_volume / self.tstep + self.surface_infil - self.pave_evap;
        max_rate = cmax(max_rate, 0.0);
        self.pave_perc = cmin(self.pave_perc, max_rate);
        if soil_thickness > 0.0 {
            self.soil_perc = self.get_soil_perc_rate(soil_theta);
            avail_volume = (soil_theta - soil_field_cap) * soil_thickness;
            max_rate = cmax(avail_volume, 0.0) / self.tstep - self.soil_evap;
            self.soil_perc = cmin(self.soil_perc, max_rate);
            self.soil_perc = cmax(self.soil_perc, 0.0);
        } else {
            self.soil_perc = self.pave_perc;
        }
        self.storage_exfil = self.get_storage_exfil_rate(g);
        self.storage_drain = 0.0;
        if p.drain_coeff > 0.0 {
            self.storage_drain = self.get_storage_drain_rate(g, storage_depth, soil_theta, pave_depth, surface_depth);
        }
        if soil_thickness == 0.0 && storage_depth >= storage_thickness && pave_depth >= pave_thickness {
            max_rate = self.storage_evap + self.storage_drain + self.storage_exfil;
            if self.pave_perc > max_rate {
                self.pave_perc = max_rate;
            } else {
                self.storage_exfil = cmin(self.storage_exfil, self.pave_perc);
                self.storage_drain = self.pave_perc - self.storage_exfil;
            }
            self.soil_perc = self.pave_perc;
            self.surface_infil = cmin(self.surface_infil, self.pave_perc);
        } else if soil_thickness > 0.0 && storage_depth >= storage_thickness && soil_theta >= soil_porosity && pave_depth >= pave_thickness {
            max_rate = self.storage_exfil + self.storage_drain;
            if self.soil_perc < max_rate {
                max_rate = self.soil_perc;
            } else {
                max_rate = cmin(max_rate, self.pave_perc);
            }
            if max_rate > self.storage_exfil {
                self.storage_drain = max_rate - self.storage_exfil;
            } else {
                self.storage_exfil = max_rate;
                self.storage_drain = 0.0;
            }
            self.soil_perc = max_rate;
            self.pave_perc = max_rate;
            self.surface_infil = cmin(self.surface_infil, self.pave_perc);
        } else if soil_thickness > 0.0 && storage_depth >= storage_thickness && soil_theta >= soil_porosity {
            max_rate = self.storage_drain + self.storage_exfil;
            if self.soil_perc > max_rate {
                self.soil_perc = max_rate;
            } else {
                self.storage_exfil = cmin(self.storage_exfil, self.soil_perc);
                self.storage_drain = self.soil_perc - self.storage_exfil;
            }
            self.pave_perc = cmin(self.pave_perc, self.soil_perc);
            avail_volume = (pave_thickness - pave_depth) * pave_void_frac;
            max_rate = avail_volume / self.tstep + self.pave_perc + self.pave_evap;
            self.surface_infil = cmin(self.surface_infil, max_rate);
        } else if soil_thickness > 0.0 && pave_depth >= pave_thickness && soil_theta >= soil_porosity {
            self.pave_perc = cmin(self.pave_perc, self.soil_perc);
            self.soil_perc = self.pave_perc;
            self.surface_infil = cmin(self.surface_infil, self.pave_perc);
            max_rate = cmax(self.storage_volume / self.tstep + self.soil_perc - self.storage_evap, 0.0);
            self.storage_exfil = cmin(self.storage_exfil, max_rate);
        } else {
            max_rate = self.soil_perc - self.storage_evap + self.storage_volume / self.tstep;
            max_rate = cmax(0.0, max_rate);
            self.storage_exfil = cmin(self.storage_exfil, max_rate);
            if self.storage_drain > 0.0 {
                max_rate = -self.storage_exfil - self.storage_evap;
                if storage_depth >= storage_thickness {
                    max_rate += self.soil_perc;
                }
                if p.drain_offset <= storage_depth {
                    max_rate += (storage_depth - p.drain_offset) * storage_void_frac / self.tstep;
                }
                max_rate = cmax(max_rate, 0.0);
                self.storage_drain = cmin(self.storage_drain, max_rate);
            }
            avail_volume = (storage_thickness - storage_depth) * storage_void_frac;
            max_rate = avail_volume / self.tstep + self.storage_evap + self.storage_drain + self.storage_exfil;
            max_rate = cmax(max_rate, 0.0);
            if soil_thickness > 0.0 {
                self.soil_perc = cmin(self.soil_perc, max_rate);
                max_rate = (soil_porosity - soil_theta) * soil_thickness / self.tstep + self.soil_perc;
            }
            self.pave_perc = cmin(self.pave_perc, max_rate);
            avail_volume = (pave_thickness - pave_depth) * pave_void_frac;
            max_rate = avail_volume / self.tstep + self.pave_perc + self.pave_evap;
            self.surface_infil = cmin(self.surface_infil, max_rate);
        }
        self.surface_outflow = self.get_surface_outflow_rate(surface_depth);
        f[L_SURF] = self.surface_inflow - self.surface_evap - self.surface_infil - self.surface_outflow;
        if self.completed {
            f[L_SURF] /= p.surf_void;
        }
        f[L_PAVE] = (self.surface_infil - self.pave_evap - self.pave_perc) / pave_void_frac;
        if p.soil_thick > 0.0 {
            f[L_SOIL] = (self.pave_perc - self.soil_evap - self.soil_perc) / soil_thickness;
            storage_inflow = self.soil_perc;
        } else {
            f[L_SOIL] = 0.0;
            storage_inflow = self.pave_perc;
            self.soil_perc = 0.0;
        }
        f[L_STOR] = (storage_inflow - self.storage_evap - self.storage_exfil - self.storage_drain) / storage_void_frac;
    }

    // legacy swaleFluxRates
    fn swale_flux_rates(&mut self, x: &[f64; L_MAX], f: &mut [f64; L_MAX]) {
        let p = self.p;
        let mut depth = x[L_SURF];
        depth = cmin(depth, p.surf_thick);
        let d_store = 0.0;
        let mut slope = p.surf_side;
        let mut top_width = p.full_width;
        top_width = cmax(top_width, 0.5);
        let mut bot_width = top_width - 2.0 * slope * p.surf_thick;
        if bot_width < 0.5 {
            bot_width = 0.5;
            slope = 0.5 * (top_width - 0.5) / p.surf_thick;
        }
        let lid_area = p.area;
        let length = lid_area / top_width;
        let surf_width = bot_width + 2.0 * slope * depth;
        let surf_area = length * surf_width;
        let mut flow_area = (depth * (bot_width + slope * depth)) * p.surf_void;
        let volume = length * flow_area;
        let surf_inflow = self.surface_inflow * lid_area;
        self.surface_evap = self.evap_rate * surf_area;
        self.surface_evap = cmin(self.surface_evap, volume / self.tstep);
        self.storage_exfil = self.surface_infil * surf_area;
        let x_depth = depth - d_store;
        if x_depth <= L_ZERO {
            self.surface_outflow = 0.0;
        } else {
            flow_area -= (d_store * (bot_width + slope * d_store)) * p.surf_void;
            if flow_area < L_ZERO {
                self.surface_outflow = 0.0;
            } else {
                bot_width += 2.0 * d_store * slope;
                let mut hyd_radius = bot_width + 2.0 * x_depth * (1.0 + slope * slope).sqrt();
                hyd_radius = flow_area / hyd_radius;
                self.surface_outflow = p.surf_alpha * flow_area * hyd_radius.powf(2. / 3.); // LIBM
            }
        }
        let mut dvdt = surf_inflow - self.surface_evap - self.storage_exfil - self.surface_outflow;
        if depth == p.surf_thick && dvdt > 0.0 {
            self.surface_outflow += dvdt;
            dvdt = 0.0;
        }
        self.surface_evap /= lid_area;
        self.storage_exfil /= lid_area;
        self.surface_outflow /= lid_area;
        f[L_SURF] = dvdt / surf_area;
        f[L_SOIL] = 0.0;
        f[L_STOR] = 0.0;
        self.surface_volume = volume / lid_area;
        self.soil_volume = 0.0;
        self.storage_volume = 0.0;
    }

    // completed swales (2D path only; translated for completeness of getOutflow)
    fn completed_swale(&mut self, x: &mut [f64; L_MAX], f: &mut [f64; L_MAX]) {
        let p = self.p;
        let top = cmax(p.full_width, 0.5);
        let height = p.surf_thick;
        let mut side = p.surf_side;
        let mut bottom = top - 2.0 * side * height;
        if bottom < 0.5 {
            bottom = 0.5;
            side = 0.5 * (top - 0.5) / height;
        }
        let volume = |depth: f64| depth * (bottom + side * depth) * p.surf_void / top;
        let old_depth = x[L_SURF];
        let initial = volume(old_depth);
        let xc = *x;
        self.swale_flux_rates(&xc, f);
        let available = initial + self.surface_inflow * self.tstep;
        let demand = (self.surface_evap + self.storage_exfil + self.surface_outflow) * self.tstep;
        let scale = if demand > available && demand > 0.0 { available / demand } else { 1.0 };
        self.surface_evap *= scale;
        self.storage_exfil *= scale;
        self.surface_outflow *= scale;
        let mut fin = cmax(0.0, available - demand * scale);
        let excess = cmax(0.0, fin - volume(height));
        self.surface_outflow += excess / self.tstep;
        fin -= excess;
        let section = fin * top / p.surf_void;
        x[L_SURF] = 2.0 * section / (bottom + (bottom * bottom + 4.0 * side * section).sqrt());
        self.surface_volume = fin;
        f[L_SURF] = (x[L_SURF] - old_depth) / self.tstep;
    }

    // legacy barrelFluxRates
    fn barrel_flux_rates(&mut self, g: &LIDGroupSoA, x: &[f64; L_MAX], f: &mut [f64; L_MAX]) {
        let p = self.p;
        let storage_depth = x[L_STOR];
        self.surface_volume = 0.0;
        self.soil_volume = 0.0;
        self.storage_volume = storage_depth;
        self.surface_infil = 0.0;
        self.surface_outflow = 0.0;
        self.storage_drain = 0.0;
        if p.drain_delay == 0.0 || g.dry_time[self.u] >= p.drain_delay {
            let head = storage_depth - p.drain_offset;
            if head > 0.0 {
                self.storage_drain = self.get_storage_drain_rate(g, storage_depth, 0.0, 0.0, 0.0);
                let max_value = head / self.tstep;
                self.storage_drain = cmin(self.storage_drain, max_value);
            }
        }
        self.storage_inflow = self.surface_inflow;
        let max_value = (p.stor_thick - storage_depth) / self.tstep + self.storage_drain;
        self.storage_inflow = cmin(self.storage_inflow, max_value);
        self.surface_infil = self.storage_inflow;
        f[L_SURF] = self.surface_inflow - self.storage_inflow;
        if self.completed {
            f[L_SURF] /= p.surf_void;
        }
        f[L_STOR] = self.storage_inflow - self.storage_drain;
        f[L_SOIL] = 0.0;
    }

    fn flux_rates(&mut self, g: &mut LIDGroupSoA, x: &[f64; L_MAX], f: &mut [f64; L_MAX]) {
        match self.p.type_ {
            LIDType::BioCell | LIDType::RainGarden => self.biocell_flux_rates(g, x, f),
            LIDType::GreenRoof => self.green_roof_flux_rates(x, f),
            LIDType::InfilTrench => self.trench_flux_rates(g, x, f),
            LIDType::PermPavement => self.pavement_flux_rates(g, x, f),
            LIDType::RainBarrel => self.barrel_flux_rates(g, x, f),
            LIDType::RoofDiscon => self.roof_flux_rates(x, f),
            LIDType::VegSwale => self.swale_flux_rates(x, f),
        }
    }

    // legacy modpuls_solve
    #[allow(clippy::too_many_arguments)]
    fn modpuls_solve(&mut self, g: &mut LIDGroupSoA, x: &mut [f64; L_MAX], x_min: &[f64; L_MAX], x_max: &[f64; L_MAX], x_tol: &[f64; L_MAX], q_old: &[f64; L_MAX], q: &mut [f64; L_MAX], dt: f64, omega: f64) -> i32 {
        let mut steps = 1;
        let max_steps = 20;
        let x_old = *x;
        let mut x_prev = *x;
        while steps < max_steps {
            let mut can_stop = true;
            let xc = *x;
            self.flux_rates(g, &xc, q);
            for i in 0..L_MAX {
                x[i] = x_old[i] + (omega * q_old[i] + (1.0 - omega) * q[i]) * dt;
                x[i] = cmin(x[i], x_max[i]);
                x[i] = cmax(x[i], x_min[i]);
                if omega > 0.0 && (x[i] - x_prev[i]).abs() > x_tol[i] {
                    can_stop = false;
                }
                x_prev[i] = x[i];
            }
            if can_stop {
                return steps;
            }
            steps += 1;
        }
        0
    }

    // legacy lidproc_getOutflow
    #[allow(clippy::too_many_arguments)]
    fn get_outflow(&mut self, g: &mut LIDGroupSoA, inflow: f64, evap: f64, infil: f64, max_infil: f64, t_step: f64, infil_factor: f64, recovery_factor: f64, lid_evap: &mut f64, lid_infil: &mut f64, lid_drain: &mut f64) -> f64 {
        let u = self.u;
        let mut x = [0.0f64; L_MAX];
        let mut x_min = [0.0f64; L_MAX];
        let mut x_max = [0.0f64; L_MAX];
        let mut f_old = [0.0f64; L_MAX];
        let mut f = [0.0f64; L_MAX];
        let x_tol = [L_STOPTOL, L_STOPTOL, L_STOPTOL, L_STOPTOL];
        let mut omega = 0.0;

        self.tstep = t_step;
        self.evap_rate = evap;
        self.max_native_infil = max_infil;
        x[L_SURF] = g.surf_depth[u];
        x[L_SOIL] = g.soil_moist[u];
        x[L_STOR] = g.stor_depth[u];
        x[L_PAVE] = g.pave_depth[u];

        self.surface_volume = 0.0;
        self.pave_volume = 0.0;
        self.soil_volume = 0.0;
        self.storage_volume = 0.0;
        self.surface_inflow = inflow;
        self.surface_infil = 0.0;
        self.surface_evap = 0.0;
        self.surface_outflow = 0.0;
        self.pave_evap = 0.0;
        self.pave_perc = 0.0;
        self.soil_evap = 0.0;
        self.soil_perc = 0.0;
        self.storage_inflow = 0.0;
        self.storage_exfil = 0.0;
        self.storage_evap = 0.0;
        self.storage_drain = 0.0;
        f_old[L_SURF] = g.f_old_surf[u];
        f_old[L_SOIL] = g.f_old_soil[u];
        f_old[L_STOR] = g.f_old_stor[u];
        f_old[L_PAVE] = g.f_old_pave[u];
        for i in 0..L_MAX {
            f[i] = 0.0;
            x_min[i] = 0.0;
            x_max[i] = L_BIG;
        }

        if self.p.type_ == LIDType::PermPavement {
            self.surface_infil = 0.0;
        } else if g.soil_infil[u].Ks > 0.0 && !(self.spatial_fraction == 1.0 && self.p.type_ == LIDType::VegSwale) {
            let sd = g.surf_depth[u];
            self.surface_infil = infil::grnampt_get_infil(&mut g.soil_infil[u], self.surface_inflow, sd, self.tstep, InfilModel::ModGreenAmpt as i32, infil_factor, recovery_factor);
        } else {
            self.surface_infil = infil;
        }
        if self.spatial_fraction >= 0.0 && self.p.type_ == LIDType::VegSwale {
            let native = cmin(self.surface_infil, self.max_native_infil);
            self.surface_infil = self.spatial_fraction * self.spatial_limit + (1.0 - self.spatial_fraction) * native;
            self.spatial_exfil_ratio = if self.surface_infil > 0.0 { self.spatial_fraction * self.spatial_limit / self.surface_infil } else { 0.0 };
        }

        if self.completed && self.spatial_fraction < 0.0 && self.p.type_ == LIDType::VegSwale {
            self.surface_infil = cmin(self.surface_infil, self.max_native_infil);
        }

        if self.p.soil_thick > 0.0 {
            x_min[L_SOIL] = self.p.soil_wp;
            x_max[L_SOIL] = self.p.soil_poros;
        }
        if self.p.pave_thick > 0.0 {
            x_max[L_PAVE] = self.p.pave_thick;
        }
        if self.p.stor_thick > 0.0 {
            x_max[L_STOR] = self.p.stor_thick;
        }
        if self.p.type_ == LIDType::GreenRoof {
            x_max[L_STOR] = self.p.dm_thick;
        }
        if self.p.type_ == LIDType::VegSwale {
            omega = 0.5;
        }

        if self.completed && self.p.type_ == LIDType::VegSwale {
            self.completed_swale(&mut x, &mut f);
        } else {
            self.modpuls_solve(g, &mut x, &x_min, &x_max, &x_tol, &f_old, &mut f, t_step, omega);
        }

        if self.p.can_overflow || self.p.full_width == 0.0 {
            let mut xs = x[L_SURF];
            self.surface_outflow += self.get_surface_overflow_rate(&mut xs);
            x[L_SURF] = xs;
        }

        g.surf_depth[u] = x[L_SURF];
        g.pave_depth[u] = x[L_PAVE];
        g.soil_moist[u] = x[L_SOIL];
        g.stor_depth[u] = x[L_STOR];
        g.f_old_surf[u] = f[L_SURF];
        g.f_old_soil[u] = f[L_SOIL];
        g.f_old_stor[u] = f[L_STOR];
        g.f_old_pave[u] = f[L_PAVE];

        *lid_evap = self.surface_evap + self.pave_evap + self.soil_evap + self.storage_evap;
        *lid_infil = self.storage_exfil;
        *lid_drain = self.storage_drain;
        self.surface_outflow
    }
}

/// makeView: one unit's parameters (legacy TLidProc + TLidUnit area/width).
pub fn make_view(g: &LIDGroupSoA, u: usize) -> LidProcView {
    LidProcView {
        type_: g.type_,
        surf_thick: g.surf_store[u],
        surf_void: g.surf_void_frac[u],
        surf_alpha: g.surf_alpha[u],
        surf_side: g.surf_side_slope[u],
        can_overflow: g.can_overflow[u] != 0,
        pave_thick: g.pave_thick[u],
        pave_void: g.pave_void[u],
        pave_imperv: g.pave_imperv_frac[u],
        pave_ksat: g.pave_ksat[u],
        pave_clog: g.pave_clog_factor[u],
        pave_regen_days: g.pave_regen_days[u],
        pave_regen_deg: g.pave_regen_deg[u],
        soil_thick: g.soil_thick[u],
        soil_poros: g.soil_poros[u],
        soil_fc: g.soil_fc[u],
        soil_wp: g.soil_wp[u],
        soil_ksat: g.soil_ksat[u],
        soil_kslope: g.soil_kslope[u],
        stor_thick: g.stor_thick[u],
        stor_void: g.stor_void[u],
        stor_ksat: g.stor_ksat[u],
        stor_clog: g.stor_clog[u],
        drain_coeff: g.drain_coeff[u],
        drain_expon: g.drain_expon[u],
        drain_offset: g.drain_offset[u],
        drain_delay: g.drain_delay[u],
        drain_hopen: g.drain_hopen[u],
        drain_hclose: g.drain_hclose[u],
        dm_thick: g.drainmat_thick[u],
        dm_alpha: g.drainmat_alpha[u],
        area: if g.unit_area[u] > 0.0 { g.unit_area[u] } else { g.area[u] },
        full_width: if g.unit_area[u] > 0.0 { g.unit_width[u] } else { g.full_width[u] },
        ucf_rainfall: g.ucf_rainfall,
        ucf_raindepth: g.ucf_raindepth,
    }
}

/// runUnitLegacy: one unit's step (legacy lidproc_getOutflow + lidproc_saveResults' water
/// balance), on the ordinary 1D path (no spatial fraction, not completed).
#[allow(clippy::too_many_arguments)]
pub fn run_unit_legacy(g: &mut LIDGroupSoA, u: usize, inflow: f64, evap: f64, native_infil: f64, max_native_infil: f64, dt: f64, old_runoff_sec: f64, infil_factor: f64, recovery_factor: f64) {
    let view = make_view(g, u);
    let mut k = LegacyLidKernel::new(view, u);
    k.old_runoff_days = old_runoff_sec / 86400.0;
    k.spatial_fraction = -1.0;
    k.spatial_limit = L_BIG;
    k.completed = false;
    let completed = false;

    let mut lid_evap = 0.0;
    let mut lid_infil = 0.0;
    let mut lid_drain = 0.0;
    let lid_runoff = k.get_outflow(g, inflow, evap, native_infil, max_native_infil, dt, infil_factor, recovery_factor, &mut lid_evap, &mut lid_infil, &mut lid_drain);

    let total_evap = k.surface_evap + k.pave_evap + k.soil_evap + k.storage_evap;
    let total_volume = k.surface_volume + k.pave_volume + k.soil_volume + k.storage_volume;
    g.vol_treated[u] += k.surface_inflow * dt;
    g.wb_inflow[u] += k.surface_inflow * dt;
    g.wb_evap[u] += total_evap * dt;
    g.wb_infil[u] += k.storage_exfil * dt;
    g.wb_surf_flow[u] += k.surface_outflow * dt;
    g.wb_drain_flow[u] += k.storage_drain * dt;
    g.wb_final_vol[u] = total_volume;
    let is_dry = k.surface_inflow < L_MINFLOW && k.surface_outflow < L_MINFLOW && k.storage_drain < L_MINFLOW && k.storage_exfil < L_MINFLOW && total_evap < L_MINFLOW;
    g.is_wet[u] = if is_dry { 0 } else { 1 };

    g.surface_runoff[u] = lid_runoff;
    g.drain_flow[u] = lid_drain;
    g.evap_loss[u] = lid_evap * dt;
    g.infil_loss[u] = lid_infil * dt;
    g.spatial_infil_loss[u] = g.infil_loss[u] * k.spatial_exfil_ratio;
    g.wb_spatial_infil[u] += g.spatial_infil_loss[u];

    g.rpt_surf_infil[u] = k.surface_infil;
    g.rpt_pave_perc[u] = k.pave_perc;
    g.rpt_soil_perc[u] = k.soil_perc;

    g.in_surf[u] = k.surface_inflow;
    if completed {
        g.evap_surf[u] = k.surface_evap;
        g.evap_pave[u] = k.pave_evap;
        g.evap_soil[u] = k.soil_evap;
        g.evap_stor[u] = k.storage_evap;
    }
    g.in_pave[u] = if g.type_ == LIDType::PermPavement { k.surface_infil } else { 0.0 };
    g.in_soil[u] = if g.type_ == LIDType::PermPavement {
        k.pave_perc
    } else if g.soil_thick[u] > 0.0 {
        k.surface_infil
    } else {
        0.0
    };
    g.in_stor[u] = if g.type_ == LIDType::RainBarrel || g.type_ == LIDType::InfilTrench {
        k.storage_inflow
    } else if g.soil_thick[u] > 0.0 && g.stor_thick[u] > 0.0 {
        k.soil_perc
    } else if g.type_ == LIDType::PermPavement {
        k.pave_perc
    } else {
        0.0
    };
}

impl LIDSolver {
    /// LIDSolver::execute: every unit with area, in group order (legacy lid_getRunoff ->
    /// evalLidUnit), then its rain-barrel dry clock.
    pub fn execute(&mut self, ctx: &SimulationContext, dt: f64, _rainfall: f64, evap_rate: f64, source_owners: Option<&Vec<i32>>) -> Result<(), &'static str> {
        const MIN_RUNOFF: f64 = 2.31481e-8;
        if let Some(o) = source_owners {
            if o.len() != ctx.n_subcatches() as usize {
                return Err("Ordinary LID ownership requires a matching source map.");
            }
        }
        for gi in 0..self.groups_.len() {
            if self.groups_[gi].count == 0 {
                continue;
            }
            for uu in 0..self.groups_[gi].count as usize {
                let g = &mut self.groups_[gi];
                if g.area[uu] <= 0.0 {
                    continue;
                }
                let sc = g.subcatch_idx[uu];
                if let Some(o) = source_owners {
                    if sc >= 0 && o[sc as usize] >= 0 {
                        continue;
                    }
                }
                let usc = sc as usize;
                let e_dry = if ctx.options.evap_dry_only && g.subcatch_rain[uu] > 0.0 { 0.0 } else { evap_rate };
                let evap = if sc >= 0 { ctx.forcing.effective_evap_rate(usc, e_dry) } else { e_dry };
                g.evap_rate_unit[uu] = evap;
                let native = if sc >= 0 && usc < self.native_infil_.len() { self.native_infil_[usc] } else { 0.0 };
                let max_native = if sc >= 0 && usc < self.max_native_infil_.len() { self.max_native_infil_[usc] } else { L_BIG };
                let infil_factor = if sc >= 0 && usc < self.infil_factor_.len() { self.infil_factor_[usc] } else { 1.0 };
                let inflow = g.inflow[uu];
                run_unit_legacy(g, uu, inflow, evap, native, max_native, dt, self.old_runoff_sec_, infil_factor, self.recovery_factor_);
                if g.subcatch_rain[uu] > MIN_RUNOFF {
                    g.dry_time[uu] = 0.0;
                } else {
                    g.dry_time[uu] += dt;
                }
                if !self.groups_[gi].rpt_path[uu].is_empty() {
                    self.save_lid_report(ctx, gi as i32, uu, dt);
                }
            }
        }
        Ok(())
    }

    /// saveLidReport: one row of the unit's [LID_USAGE] report file (legacy createLidRptFile,
    /// initLidRptFile and the writing half of lidproc_saveResults).
    pub fn save_lid_report(&mut self, ctx: &SimulationContext, gi: i32, u: usize, dt: f64) {
        let key = (gi, u as i32);
        let g = &self.groups_[gi as usize];
        if !self.rpt_files_.contains_key(&key) {
            let mut st = Box::new(LidRptFile::default());
            use std::os::unix::ffi::OsStrExt;
            let path = std::path::Path::new(std::ffi::OsStr::from_bytes(&g.rpt_path[u]));
            st.f = std::fs::File::create(path).ok();
            st.failed = st.f.is_none();
            if !st.failed {
                const HEAD1: [&str; 14] = ["\n                    \t", "  Elapsed\t", "    Total\t", "    Total\t", "  Surface\t", " Pavement\t", "     Soil\t", "  Storage\t", "  Surface\t", "    Drain\t", "  Surface\t", " Pavement\t", "     Soil\t", "  Storage"];
                const HEAD2: [&str; 14] = ["\n                    \t", "     Time\t", "   Inflow\t", "     Evap\t", "    Infil\t", "     Perc\t", "     Perc\t", "    Exfil\t", "   Runoff\t", "  OutFlow\t", "    Level\t", "    Level\t", " Moisture\t", "    Level"];
                const UNITS1: [&str; 14] = ["\nDate        Time    \t", "    Hours\t", "    in/hr\t", "    in/hr\t", "    in/hr\t", "    in/hr\t", "    in/hr\t", "    in/hr\t", "    in/hr\t", "    in/hr\t", "   inches\t", "   inches\t", "  Content\t", "   inches"];
                const UNITS2: [&str; 14] = ["\nDate        Time    \t", "    Hours\t", "    mm/hr\t", "    mm/hr\t", "    mm/hr\t", "    mm/hr\t", "    mm/hr\t", "    mm/hr\t", "    mm/hr\t", "    mm/hr\t", "       mm\t", "       mm\t", "  Content\t", "       mm"];
                let is_us = g.ucf_raindepth == 12.0;
                let title: &[u8] = if ctx.title_notes.is_empty() { b"" } else { &ctx.title_notes[0] };
                let ci = g.control_idx[u];
                let si = g.subcatch_idx[u];
                let mut o: Vec<u8> = Vec::new();
                o.extend_from_slice(b"SWMM5 LID Report File\n");
                o.extend_from_slice(b"\nProject:  ");
                o.extend_from_slice(title);
                o.extend_from_slice(b"\nLID Unit: ");
                if ci >= 0 {
                    o.extend_from_slice(ctx.lid_names.name_of(ci));
                }
                o.extend_from_slice(b" in Subcatchment ");
                if si >= 0 {
                    o.extend_from_slice(ctx.subcatch_names.name_of(si));
                }
                o.extend_from_slice(b"\n");
                for c in HEAD1.iter().chain(HEAD2.iter()) {
                    o.extend_from_slice(c.as_bytes());
                }
                for c in if is_us { UNITS1.iter() } else { UNITS2.iter() } {
                    o.extend_from_slice(c.as_bytes());
                }
                o.extend_from_slice(b"\n----------- --------");
                for _ in 1..14 {
                    o.extend_from_slice(b"\t ---------");
                }
                if let Some(f) = st.f.as_mut() {
                    let _ = f.write_all(&o);
                }
            }
            self.rpt_files_.insert(key, st);
        }
        let g = &self.groups_[gi as usize];
        let r = self.rpt_files_.get_mut(&key).unwrap();
        if r.failed {
            return;
        }
        let total_evap = if dt > 0.0 { g.evap_loss[u] / dt } else { 0.0 };
        let stor_exfil = if dt > 0.0 { g.infil_loss[u] / dt } else { 0.0 };
        let (ucf_q, ucf_d) = (g.ucf_rainfall, g.ucf_raindepth);
        let rv: [f64; 12] = [
            g.in_surf[u] * ucf_q,
            total_evap * ucf_q,
            g.rpt_surf_infil[u] * ucf_q,
            g.rpt_pave_perc[u] * ucf_q,
            g.rpt_soil_perc[u] * ucf_q,
            stor_exfil * ucf_q,
            g.surface_runoff[u] * ucf_q,
            g.drain_flow[u] * ucf_q,
            g.surf_depth[u] * ucf_d,
            g.pave_depth[u] * ucf_d,
            if g.soil_thick[u] > 0.0 { g.soil_moist[u] } else { 0.0 },
            g.stor_depth[u] * ucf_d,
        ];
        let is_dry = g.in_surf[u] < L_MINFLOW && g.surface_runoff[u] < L_MINFLOW && g.drain_flow[u] < L_MINFLOW && stor_exfil < L_MINFLOW && total_evap < L_MINFLOW;

        if !is_dry && r.was_dry > 1 {
            if let Some(f) = r.f.as_mut() {
                let _ = f.write_all(&r.held);
            }
        }

        let new_sec = self.old_runoff_sec_ + dt;
        let when = ctx.options.start_date + new_sec / 86400.0;
        let (yy, mo, dd) = datetime::decode_date(when);
        let (hh, mi, ss) = datetime::decode_time(when);
        let ts = format!("{:02}/{:02}/{:04} {:02}:{:02}:{:02}", mo, dd, yy, hh, mi, ss);
        let mut row = format!("\n{:>20}\t {}\t {}\t {}\t {}\t {}\t {}\t {}\t{}\t {}\t {}\t {}\t {}\t {}", ts, fmt_f(new_sec / 3600.0, 3), fmt_f(rv[0], 3), fmt_f(rv[1], 4), fmt_f(rv[2], 3), fmt_f(rv[3], 3), fmt_f(rv[4], 3), fmt_f(rv[5], 3), fmt_f(rv[6], 3), fmt_f(rv[7], 3), fmt_f(rv[8], 3), fmt_f(rv[9], 3), fmt_f(rv[10], 3), fmt_f(rv[11], 3));
        row.truncate(511); // char row[512]
        r.held = row.into_bytes();

        if is_dry {
            if r.was_dry == 0 {
                if let Some(f) = r.f.as_mut() {
                    let _ = f.write_all(&r.held);
                }
            }
            r.was_dry += 1;
        } else {
            if let Some(f) = r.f.as_mut() {
                let _ = f.write_all(&r.held);
            }
            r.was_dry = 0;
        }
    }
}

/// printf("%8.<prec>f", v) as glibc prints it.
fn fmt_f(v: f64, prec: usize) -> String {
    if v.is_nan() {
        let s = if v.is_sign_negative() { "-nan" } else { "nan" };
        return format!("{:>8}", s);
    }
    format!("{:>8.*}", prec, v)
}
