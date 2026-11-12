// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/plugins/DefaultReportPlugin.cpp
// (DefaultReportPlugin::write_results) and src/engine/quality/MsxLidReactions.cpp (heldMass)
// (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026 Caleb
// Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// The legacy-format .rpt writer, part 2: write_results, the continuity blocks (RDII, runoff
// quantity and quality, groundwater, flow routing, quality routing, reaction species), the
// accuracy and time-step statistics, and the summary tables (subcatchment runoff and
// washoff, groundwater, LID performance, node depth, inflow, surcharge, flooding, storage,
// outfall loading, link flow, flow classification, conduit surcharge, pumping, streets,
// inlets, link pollutant loads). See core/report_plugin.rs for the conventions.
//
// The 2D blocks (12 November 2026, phase 2a wave 5): the 2D surface routing continuity, the
// 2D solver statistics and the 1D <-> 2D exchange reconciliation are compiled in both builds,
// as upstream compiles them, and print only when the 2D router has set
// mass_balance_2d.active (never in a 1D build). Upstream's #ifdef OPENSWMM_HAS_2D blocks (the
// 2D aquifer and aquifer quality continuity, the 2D aquifer clock, the 2D surface washoff
// summary and the 2D outfall tailwater summary) compile only with the `twod` feature; they
// read the router, which upstream reaches through ctx.twod_io.surface_router and the engine
// passes here. The aquifer blocks, the clock and the washoff summary need a phase-2b owner
// (the two-zone kernel, surface quality) that SurfaceRouter2D::initialize refuses, so a run
// cannot reach them; each is guarded with the crate's "STUB:" panic.

use crate::cfprintf;
use crate::core::cmath::{cclamp, cmax, cmin};
use crate::core::context::{SimulationContext, MAX_STATS, K_MAX_CONTROL_LOG};
use crate::core::cprintf::CFile;
use crate::core::datetime;
use crate::core::context::FileMode;
use crate::core::options::RoutingModel;
use crate::core::plugins::DefaultReportPlugin;
use crate::core::report_plugin::{elapsed_to_parts, flow_fmt, land_area_to_ft2, lt_str, nt_str, wr, FLOW_UNIT_WORDS};
use crate::core::unit_conversion as ucf;
use crate::data::gage::RainSource;
use crate::data::link::{LinkData, LinkType, XsectShape};
use crate::data::node::NodeType;
use crate::data::pollutant::MassUnits;

/// msxlid::heldMass: a node's LID-held mass of one reaction species.
fn msx_lid_held_mass(c: &SimulationContext, node: i32, species: i32) -> f64 {
    let row = c.node_subtypes.storage_row(node);
    let ns = c.reactions.n_species();
    if row < 0 || row >= c.reactions.lid.node_mass.len() as i32 || species < 0 || species >= ns {
        return 0.0;
    }
    let mut total = 0.0;
    let mass = &c.reactions.lid.node_mass[row as usize];
    let mut i = species as usize;
    while i < mass.len() {
        total += mass[i];
        i += ns as usize;
    }
    total
}

/// `p < v.size() ? v[p] : 0.0`
fn at(v: &[f64], p: i32) -> f64 {
    if p >= 0 && (p as usize) < v.len() {
        v[p as usize]
    } else {
        0.0
    }
}
fn atu(v: &[f64], s: usize) -> f64 {
    if s < v.len() {
        v[s]
    } else {
        0.0
    }
}

/// The relative continuity error of legacy's quality and loading balances.
fn qual_err(in_: f64, out_: f64) -> f64 {
    if (in_ - out_).abs() < 0.001 {
        0.0
    } else if in_ > 0.0 {
        100.0 * (1.0 - out_ / in_)
    } else if out_ > 0.0 {
        100.0 * (in_ / out_ - 1.0)
    } else {
        0.0
    }
}

impl DefaultReportPlugin {
    /// write_results(f, ctx): the simulation result sections.
    pub(crate) fn write_results(&mut self, f: &mut CFile, ctx: &SimulationContext, #[cfg(feature = "twod")] router: Option<&crate::twod::surface_router_2d::SurfaceRouter2D>) {
        let opt = &ctx.options;
        let mut fu = opt.flow_units as i32;
        if !(0..=5).contains(&fu) {
            fu = 0;
        }
        let ff = flow_fmt(fu);
        let fuw = FLOW_UNIT_WORDS[fu as usize];

        let du = ucf::DisplayUnits::from(opt);
        let si_report = du.unit_system == 1;
        let qcf = du.flow;
        let depth_vcf = du.raindepth;
        let len_ucf = du.length;
        let svol_ucf = du.volume;
        let len_word = du.length_word;
        let vcf = if si_report { 28.317 / 1.0e6 } else { 7.48 / 1.0e6 };
        let vol_word = du.mvol_word;
        let cont_us = if si_report { 1 } else { 0 };
        let cont_vcf1 = ucf::UCF[ucf::LENGTH as usize][cont_us] * ucf::UCF[ucf::LANDAREA as usize][cont_us];
        let cont_vcf2 = if si_report { 2.4466 / 86400.0 } else { 0.64632 / 86400.0 };
        let cont_mass = ucf::UCF[ucf::MASS as usize][cont_us];
        let load_word: &[u8] = if si_report { b"kg" } else { b"lbs" };
        let is_count = |p: i32| -> bool { p >= 0 && (p as usize) < ctx.pollutants.units.len() && ctx.pollutants.units[p as usize] == MassUnits::CountsPerL };
        let depth_word = du.depth_word;
        let _ = len_word;

        let has_rdii = !ctx.rdii_assigns.node_idx.is_empty();
        let mut has_gw = false;
        for j in 0..ctx.n_subcatches() {
            if ctx.subcatches.gw_aquifer[j as usize] >= 0 {
                has_gw = true;
                break;
            }
        }

        // ---- RDII Continuity (legacy report_writeRdiiError) ----
        if has_rdii && opt.rpt_continuity && !opt.ignore_rainfall && !opt.ignore_rdii && ctx.files.rdii_mode != FileMode::Use {
            let mb = &ctx.mass_balance;
            let sewer_rain = mb.rdii_rainfall;
            let rdii_prod = mb.rdii_produced;
            let ratio = if sewer_rain > 0.0 { rdii_prod / sewer_rain } else { 0.0 };
            let ucf1 = cont_vcf1;
            let ucf2 = cont_vcf2;
            cfprintf!(f, b"\n  **********************           Volume        Volume");
            cfprintf!(f, if si_report { b"\n  Rainfall Dependent I/I        hectare-m      10^6 ltr" as &[u8] } else { b"\n  Rainfall Dependent I/I        acre-feet      10^6 gal" });
            cfprintf!(f, b"\n  **********************        ---------     ---------");
            cfprintf!(f, b"\n  Sewershed Rainfall ......%14.3f%14.3f", sewer_rain * ucf1, sewer_rain * ucf2);
            cfprintf!(f, b"\n  RDII Produced ...........%14.3f%14.3f", rdii_prod * ucf1, rdii_prod * ucf2);
            cfprintf!(f, b"\n  RDII Ratio ..............%14.3f", ratio);
            wr(f, b"");
            wr(f, b"");
        }

        // ---- Runoff Quantity Continuity (legacy report_writeRunoffError) ----
        {
            let mb = &ctx.mass_balance;
            let mut total_area_ft2 = 0.0;
            for i in 0..ctx.n_subcatches() {
                total_area_ft2 += ctx.subcatches.area[i as usize] * land_area_to_ft2(fu);
            }
            let runoff_err_pct = mb.runoff_error() * 100.0;
            if ctx.n_subcatches() > 0 && total_area_ft2 > 0.0 && (opt.rpt_continuity || runoff_err_pct > 10.0) {
                cfprintf!(f, b"\n  **************************        Volume         Depth");
                cfprintf!(f, if si_report { b"\n  Runoff Quantity Continuity     hectare-m            mm" as &[u8] } else { b"\n  Runoff Quantity Continuity     acre-feet        inches" });
                cfprintf!(f, b"\n  **************************     ---------       -------");
                let row = |f: &mut CFile, label: &[u8], vol_ft3: f64| {
                    let af = vol_ft3 * cont_vcf1;
                    let depth_in = vol_ft3 / total_area_ft2 * depth_vcf;
                    cfprintf!(f, b"\n  %s%14.3f%14.3f", label, af, depth_in);
                };
                let has_snow = ctx.snowpack_names.size() > 0;
                if mb.runoff_init_store > 0.0 {
                    row(f, b"Initial LID Storage ......", mb.runoff_init_store);
                }
                if has_snow {
                    row(f, b"Initial Snow Cover .......", mb.runoff_init_snow);
                }
                row(f, b"Total Precipitation ......", mb.runoff_rainfall);
                if mb.runoff_runon > 0.0 {
                    row(f, b"Outfall Runon ............", mb.runoff_runon);
                }
                row(f, b"Evaporation Loss .........", mb.runoff_evap);
                row(f, b"Infiltration Loss ........", mb.runoff_infil);
                row(f, b"Surface Runoff ...........", mb.runoff_runoff);
                if mb.runoff_lid_drain > 0.0 {
                    row(f, b"LID Drainage .............", mb.runoff_lid_drain);
                }
                if has_snow {
                    row(f, b"Snow Removed .............", mb.runoff_snowremov);
                    row(f, b"Final Snow Cover .........", mb.runoff_final_snow);
                }
                row(f, b"Final Storage ............", mb.runoff_final_store);
                cfprintf!(f, b"\n  Continuity Error (%%) .....%14.3f", runoff_err_pct);
                wr(f, b"");
                wr(f, b"");
            }
        }

        // ---- Runoff Quality Continuity (legacy report_writeLoadingError) ----
        let ms_rq = &ctx.reactions.surface;
        let nm_rq = if ms_rq.active() && !opt.ignore_quality { ms_rq.n_species } else { 0 };
        if ctx.n_subcatches() > 0 && (ctx.n_pollutants() > 0 || nm_rq > 0) && !opt.ignore_quality {
            let np = ctx.n_pollutants();
            let ncol = np + nm_rq;
            let mb = &ctx.mass_balance;
            let logn = |x: f64| if x > 0.0 { x.log10() } else { x }; // LIBM
            let mut vals: Vec<[f64; 9]> = vec![[0.0; 9]; ncol as usize];
            let mut col_name: Vec<Vec<u8>> = vec![Vec::new(); ncol as usize];
            let mut col_counts = vec![false; ncol as usize];
            let mut col_mass = vec![false; ncol as usize];
            let mut max_err = 0.0f64;
            for p in 0..np {
                let in_ = at(&mb.qual_init_buildup, p) + at(&mb.qual_surface_buildup, p) + at(&mb.qual_wet_deposition, p);
                let out_ = at(&mb.qual_sweeping, p) + at(&mb.qual_infil_loss, p) + at(&mb.qual_bmp_removal, p) + at(&mb.qual_runoff_load, p) + at(&mb.qual_final_buildup, p);
                let err = qual_err(in_, out_);
                max_err = cmax(max_err, err);
                let v = &mut vals[p as usize];
                *v = [
                    at(&mb.qual_init_buildup, p),
                    at(&mb.qual_surface_buildup, p),
                    at(&mb.qual_wet_deposition, p),
                    at(&mb.qual_sweeping, p),
                    at(&mb.qual_infil_loss, p),
                    at(&mb.qual_bmp_removal, p),
                    at(&mb.qual_runoff_load, p),
                    at(&mb.qual_final_buildup, p),
                    err,
                ];
                let counts = (p as usize) < ctx.pollutants.units.len() && ctx.pollutants.units[p as usize] == MassUnits::CountsPerL;
                if counts {
                    for k in 0..8 {
                        v[k] = logn(v[k]);
                    }
                }
                col_name[p as usize] = ctx.pollutant_names.name_of(p).to_vec();
                col_counts[p as usize] = counts;
            }
            for m in 0..nm_rq {
                let um = m as usize;
                let final_bu = crate::transport::msx_surface::remaining_buildup(ctx, m);
                let in_ = at(&ms_rq.led_init_buildup, m) + at(&ms_rq.led_buildup, m) + at(&ms_rq.led_wet_deposition, m) + at(&ms_rq.led_reaction_clamped, m);
                let out_ = at(&ms_rq.led_sweeping, m) + at(&ms_rq.led_infil_loss, m) + at(&ms_rq.led_bmp_removal, m) + at(&ms_rq.led_runoff_load, m) + final_bu + at(&ms_rq.led_reacted, m);
                let err = qual_err(in_, out_);
                max_err = cmax(max_err, err);
                vals[np as usize + um] = [
                    at(&ms_rq.led_init_buildup, m),
                    at(&ms_rq.led_buildup, m),
                    at(&ms_rq.led_wet_deposition, m),
                    at(&ms_rq.led_sweeping, m),
                    at(&ms_rq.led_infil_loss, m),
                    at(&ms_rq.led_bmp_removal, m),
                    at(&ms_rq.led_runoff_load, m),
                    final_bu,
                    err,
                ];
                col_name[np as usize + um] = ctx.reactions.species_name[um].clone();
                col_mass[np as usize + um] = ms_rq.mcf[um] == 1.0;
            }
            if opt.rpt_continuity || max_err > 10.0 {
                const LABELS: [&[u8]; 9] = [
                    b"Initial Buildup ..........",
                    b"Surface Buildup ..........",
                    b"Wet Deposition ...........",
                    b"Sweeping Removal .........",
                    b"Infiltration Loss ........",
                    b"BMP Removal ..............",
                    b"Surface Runoff ...........",
                    b"Remaining Buildup ........",
                    b"Continuity Error (%) .....",
                ];
                let mut p1 = 0;
                while p1 < ncol {
                    let p2 = cmin(p1 + 5, ncol);
                    cfprintf!(f, b"\n  **************************");
                    for p in p1..p2 {
                        cfprintf!(f, b"%14s", col_name[p as usize]);
                    }
                    cfprintf!(f, b"\n  Runoff Quality Continuity ");
                    for p in p1..p2 {
                        let up = p as usize;
                        cfprintf!(f, b"%14s", if col_counts[up] { b"LogN" as &[u8] } else if col_mass[up] { b"mass" } else { load_word });
                    }
                    cfprintf!(f, b"\n  **************************");
                    for _p in p1..p2 {
                        cfprintf!(f, b"    ----------");
                    }
                    for k in 0..9 {
                        if k == 8 && (ctx.reactions.has_subcatch_reactions() || !ctx.reactions.buildup_laws.is_empty() || !ctx.reactions.lid_controls.is_empty()) {
                            for row in 0..2 {
                                cfprintf!(
                                    f,
                                    b"\n  %s",
                                    if row == 0 {
                                        if ctx.reactions.buildup_laws.is_empty() && ctx.reactions.lid_controls.is_empty() {
                                            b"Wet Pool Mass Reacted ...." as &[u8]
                                        } else {
                                            b"Surface Mass Reacted ....."
                                        }
                                    } else {
                                        b"Clamped to Zero (added) .."
                                    }
                                );
                                for p in p1..p2 {
                                    cfprintf!(f, b"%14.3f", if p < np { 0.0 } else { at(if row == 0 { &ms_rq.led_reacted } else { &ms_rq.led_reaction_clamped }, p - np) });
                                }
                            }
                        }
                        cfprintf!(f, b"\n  %s", LABELS[k]);
                        for p in p1..p2 {
                            cfprintf!(f, b"%14.3f", vals[p as usize][k]);
                        }
                    }
                    wr(f, b"");
                    wr(f, b"");
                    p1 += 5;
                }
            }
        }

        if opt.rpt_continuity && !ctx.reactions.booster_mass.is_empty() {
            cfprintf!(f, b"\n  Reaction Booster Injection (declared species mass units)\n");
            for species in 0..ctx.reactions.n_species() {
                let us = species as usize;
                cfprintf!(f, b"  %-24s %.10g %s\n", ctx.reactions.species_name[us], ctx.reactions.booster_mass[us] * 28.317, ctx.reactions.species_units[us]);
            }
            cfprintf!(f, b"\n");
        }

        // ---- Groundwater Continuity (legacy report_writeGwaterError) ----
        if has_gw && opt.rpt_continuity && !opt.ignore_groundwater {
            let mb = &ctx.mass_balance;
            let mut gw_area_ft2 = 0.0;
            for i in 0..ctx.n_subcatches() {
                let ui = i as usize;
                if ctx.subcatches.gw_aquifer[ui] >= 0 {
                    gw_area_ft2 += ctx.subcatches.area[ui] * land_area_to_ft2(fu);
                }
            }
            if gw_area_ft2 <= 0.0 {
                gw_area_ft2 = 1.0;
            }
            let ucf_vol = cont_vcf1;
            let ucf_dep = depth_vcf;
            cfprintf!(f, b"\n  **************************        Volume         Depth");
            cfprintf!(f, if si_report { b"\n  Groundwater Continuity         hectare-m            mm" as &[u8] } else { b"\n  Groundwater Continuity         acre-feet        inches" });
            cfprintf!(f, b"\n  **************************     ---------       -------");
            let gw_row = |f: &mut CFile, label: &[u8], vol_ft3: f64| {
                cfprintf!(f, b"\n  %s%14.3f%14.3f", label, vol_ft3 * ucf_vol, vol_ft3 / gw_area_ft2 * ucf_dep);
            };
            gw_row(f, b"Initial Storage ..........", mb.gw_init_storage);
            gw_row(f, b"Infiltration .............", mb.gw_infil);
            if mb.gw_infil_2d_recharge != 0.0 {
                gw_row(f, b"  of which from 2D .......", mb.gw_infil_2d_recharge);
            }
            gw_row(f, b"Upper Zone ET ............", mb.gw_upper_evap);
            gw_row(f, b"Lower Zone ET ............", mb.gw_lower_evap);
            gw_row(f, b"Deep Percolation .........", mb.gw_lower_perc);
            gw_row(f, b"Groundwater Flow .........", mb.gw_lateral_flow);
            gw_row(f, b"Final Storage ............", mb.gw_final_storage);
            let total_in = mb.gw_infil + mb.gw_init_storage;
            let total_out = mb.gw_upper_evap + mb.gw_lower_evap + mb.gw_lower_perc + mb.gw_lateral_flow + mb.gw_final_storage;
            let pct_error = if (total_in - total_out).abs() < 1.0 {
                0.0
            } else if total_in > 0.0 {
                100.0 * (1.0 - total_out / total_in)
            } else if total_out > 0.0 {
                100.0 * (total_in / total_out - 1.0)
            } else {
                0.0
            };
            cfprintf!(f, b"\n  Continuity Error (%%) .....%14.3f", pct_error);
            wr(f, b"");
            wr(f, b"");
        }

        // ---- Flow Routing Continuity (legacy report_writeFlowError) ----
        if ctx.n_nodes() > 0 && !opt.ignore_routing && (opt.rpt_continuity || ctx.mass_balance.routing_error() * 100.0 > 10.0) {
            let mb = &ctx.mass_balance;
            cfprintf!(f, b"\n  **************************        Volume        Volume");
            cfprintf!(f, if si_report { b"\n  Flow Routing Continuity        hectare-m      10^6 ltr" as &[u8] } else { b"\n  Flow Routing Continuity        acre-feet      10^6 gal" });
            cfprintf!(f, b"\n  **************************     ---------     ---------");
            let ucf1 = cont_vcf1;
            let ucf2 = cont_vcf2;
            let row = |f: &mut CFile, label: &[u8], vol_ft3: f64| {
                cfprintf!(f, b"\n  %s%14.3f%14.3f", label, vol_ft3 * ucf1, vol_ft3 * ucf2);
            };
            let rr = &mb.routing_report;
            row(f, b"Dry Weather Inflow .......", rr[0]);
            row(f, b"Wet Weather Inflow .......", rr[1]);
            row(f, b"Groundwater Inflow .......", rr[2]);
            if rr[3] != 0.0 {
                row(f, b"Conduit GW Inflow ........", rr[3]);
            }
            row(f, b"RDII Inflow ..............", rr[4]);
            row(f, b"External Inflow ..........", rr[5]);
            row(f, b"External Outflow .........", rr[8]);
            row(f, b"Flooding Loss ............", rr[6]);
            if rr[7] != 0.0 {
                row(f, b"2D Coupling Outflow ......", rr[7]);
            }
            row(f, b"Evaporation Loss .........", rr[9]);
            row(f, b"Exfiltration Loss ........", rr[10]);
            row(f, b"Initial Stored Volume ....", mb.routing_init_storage);
            row(f, b"Final Stored Volume ......", mb.routing_final_storage);
            {
                let mut slot_ft3 = 0.0;
                for &sv in &ctx.links.slot_volume {
                    slot_ft3 += sv;
                }
                if slot_ft3 > 0.0 {
                    row(f, b"Final Slot Storage .......", slot_ft3);
                }
            }
            {
                let mut tin = mb.routing_init_storage + rr[1] + rr[4] + rr[3];
                let mut tout = mb.routing_final_storage + rr[6] + rr[7] + rr[9] + rr[10];
                if rr[0] >= 0.0 {
                    tin += rr[0];
                } else {
                    tout -= rr[0];
                }
                if rr[2] >= 0.0 {
                    tin += rr[2];
                } else {
                    tout -= rr[2];
                }
                if rr[5] >= 0.0 {
                    tin += rr[5];
                } else {
                    tout -= rr[5];
                }
                if rr[8] >= 0.0 {
                    tout += rr[8];
                } else {
                    tin -= rr[8];
                }
                let pct = if (tin - tout).abs() < 1.0 {
                    0.0
                } else if tin.abs() > 0.0 {
                    100.0 * (1.0 - tout / tin)
                } else if tout.abs() > 0.0 {
                    100.0 * (tin / tout - 1.0)
                } else {
                    0.0
                };
                cfprintf!(f, b"\n  Continuity Error (%%) .....%14.3f", pct);
            }
        }

        // ---- 2D Surface Routing Continuity (SI volumes; a separate balance from the 1D one) ----
        if opt.rpt_continuity && ctx.mass_balance_2d.active {
            let mb2 = &ctx.mass_balance_2d;
            wr(f, b"");
            wr(f, b"");
            cfprintf!(f, b"\n  **************************        Volume        Volume");
            cfprintf!(f, b"\n  2D Surface Routing Continuity  cubic meters      10^6 ltr");
            cfprintf!(f, b"\n  **************************     ---------     ---------");
            let row2 = |f: &mut CFile, label: &[u8], vol_m3: f64| {
                cfprintf!(f, b"\n  %s%14.3f%14.3f", label, vol_m3, vol_m3 * 1.0e-3);
            };
            row2(f, b"Initial Stored Volume ....", mb2.init_storage);
            row2(f, b"Rainfall Inflow ..........", mb2.rainfall_in);
            row2(f, b"1D -> 2D Spill Inflow ....", mb2.coupling_1d_to_2d_in);
            row2(f, b"Outfall Inflow ...........", mb2.outfall_in);
            row2(f, b"Boundary Inflow ..........", mb2.boundary_in);
            if mb2.external_in != 0.0 {
                row2(f, b"Runtime Source Inflow ....", mb2.external_in);
            }
            if mb2.external_out != 0.0 {
                row2(f, b"Runtime Source Outflow ...", mb2.external_out);
            }
            row2(f, b"2D -> 1D Drain Outflow ...", mb2.coupling_2d_to_1d_out);
            row2(f, b"Outfall Withdrawal .......", mb2.outfall_out);
            row2(f, b"Boundary Outflow .........", mb2.boundary_out);
            row2(f, b"Evaporation Loss .........", mb2.evap_out);
            row2(f, b"Infiltration Loss ........", mb2.infil_out);
            // the SUBCATCH_AQUIFER share: delivered to the legacy aquifer, and still in flight
            if mb2.infil_to_aquifer > 0.0 || mb2.infil_aquifer_pending > 0.0 {
                row2(f, b"  to Aquifer (delivered) .", mb2.infil_to_aquifer);
                row2(f, b"  to Aquifer (in flight) .", mb2.infil_aquifer_pending);
            }
            row2(f, b"Final Stored Volume ......", mb2.final_storage);

            cfprintf!(f, b"\n  Continuity Error (%%) .....%14.3f", mb2.error() * 100.0);

            // OPENSWMM_HAS_2D: the two-zone [2D_AQUIFER] ledger and its species continuity
            // (ctx.twod_io.aquifer_state / aquifer_transport: the subsurface kernel, phase 2b)
            #[cfg(feature = "twod")]
            if let Some(r) = router {
                if r.subsurface_.active() {
                    panic!("STUB: the 2D Aquifer Continuity report blocks (DefaultReportPlugin::write_results; phase 2b)");
                }
            }

            // 2D Solver Statistics: cumulative marcher throughput, when populated
            if mb2.solver_nsteps >= 0 {
                wr(f, b"");
                wr(f, b"");
                cfprintf!(f, b"\n  *************************");
                cfprintf!(f, b"\n  2D Solver Statistics");
                cfprintf!(f, b"\n  *************************");
                let srow = |f: &mut CFile, label: &[u8], v: i64| {
                    cfprintf!(f, b"\n  %s%14ld", label, v);
                };
                srow(f, b"Internal Steps ...........", mb2.solver_nsteps);
                srow(f, b"Face-Kernel Evals ........", mb2.solver_nrhs);
                cfprintf!(f, b"\n  Avg Internal Step (s) ....%14.4f", mb2.solver_avg_h);
                cfprintf!(f, b"\n  Last Internal Step (s) ...%14.4f", mb2.solver_last_h);

                // marcher telemetry: active-fraction spread and LTS tier occupancy
                if mb2.solver_active_mean >= 0.0 {
                    cfprintf!(
                        f,
                        b"\n  Active Cells min/mean/max %8.1f /%5.1f /%5.1f  (%%)",
                        100.0 * mb2.solver_active_min,
                        100.0 * mb2.solver_active_mean,
                        100.0 * mb2.solver_active_max
                    );
                }
                if mb2.solver_n_tiers > 0 {
                    let mut total: i64 = 0;
                    for k in 0..mb2.solver_n_tiers as usize {
                        total += mb2.solver_tier_cells[k];
                    }
                    if total > 0 {
                        for k in 0..mb2.solver_n_tiers {
                            cfprintf!(f, b"\n  LTS Tier %d Occupancy (%%) .%14.1f", k, 100.0 * mb2.solver_tier_cells[k as usize] as f64 / total as f64);
                        }
                    }
                }
            }

            // OPENSWMM_HAS_2D: the 2D aquifer's independent clock (phase 2b)
            #[cfg(feature = "twod")]
            if let Some(r) = router {
                if r.subsurface_.independent_clock() {
                    panic!("STUB: the 2D Aquifer Clock report block (DefaultReportPlugin::write_results; phase 2b)");
                }
            }

            // 1D <-> 2D Exchange Reconciliation: the routing continuity with the exchange
            // treated as an internal transfer (U450: the error divides by the adjusted inflow).
            if mb2.coupling_1d_to_2d_in > 0.0 || mb2.coupling_2d_to_1d_out > 0.0 {
                let mb1 = &ctx.mass_balance;
                let land_vcf = du.landvol;
                let mvol_vcf = du.mvol;
                let m3_to_1d = 1.0 / 0.028316846592;
                let spill_1d = mb2.coupling_1d_to_2d_in * m3_to_1d;
                let drain_1d = mb2.coupling_2d_to_1d_out * m3_to_1d;

                wr(f, b"");
                wr(f, b"");
                cfprintf!(f, b"\n  **************************        Volume        Volume");
                cfprintf!(f, if si_report { b"\n  1D <-> 2D Exchange Reconcil.   hectare-m      10^6 ltr" as &[u8] } else { b"\n  1D <-> 2D Exchange Reconcil.   acre-feet      10^6 gal" });
                cfprintf!(f, b"\n  **************************     ---------     ---------");
                let rowx = |f: &mut CFile, label: &[u8], vol_1d: f64| {
                    cfprintf!(f, b"\n  %s%14.3f%14.3f", label, vol_1d * land_vcf, vol_1d * mvol_vcf);
                };
                rowx(f, b"1D -> 2D Spill ...........", spill_1d);
                rowx(f, b"2D -> 1D Drain ...........", drain_1d);
                rowx(f, b"Net 1D -> 2D .............", spill_1d - drain_1d);

                let in_adj = mb1.routing_dry_weather + mb1.routing_wet_weather + mb1.routing_gw_inflow + mb1.routing_link_gw_inflow + mb1.routing_rdii + (mb1.routing_external - drain_1d) + mb1.routing_init_storage;
                let out_adj = (mb1.routing_flooding + mb1.routing_coupling_out - spill_1d) + mb1.routing_outflow + mb1.routing_evap_loss + mb1.routing_seep_loss + mb1.routing_final_storage;
                let err_adj = if in_adj > 0.0 { (in_adj - out_adj) / in_adj } else { 0.0 };
                cfprintf!(f, b"\n  Flow Continuity w/ Exchange Internal (%%) %8.3f", err_adj * 100.0);
            }
        }

        wr(f, b"");
        wr(f, b"");

        // ---- Quality Routing Continuity (legacy report_writeQualError) ----
        if ctx.n_nodes() > 0 && !opt.ignore_routing && ctx.n_pollutants() > 0 && !opt.ignore_quality {
            let np = ctx.n_pollutants();
            let mb = &ctx.mass_balance;
            const LPER_FT3: f64 = 28.317;
            let mut vals: Vec<[f64; 12]> = vec![[0.0; 12]; np as usize];
            let mut max_err = 0.0f64;
            for p in 0..np {
                let raw: [f64; 11] = [
                    at(&mb.qual_routing_dw_in, p),
                    at(&mb.qual_routing_wet, p),
                    at(&mb.qual_routing_gw_in, p),
                    at(&mb.qual_routing_ii_in, p),
                    at(&mb.qual_routing_ex_in, p),
                    at(&mb.qual_routing_outflow, p),
                    at(&mb.qual_routing_flood, p),
                    at(&mb.qual_routing_seep, p),
                    at(&mb.qual_routing_reacted, p),
                    at(&mb.qual_routing_init, p),
                    at(&mb.qual_routing_final, p),
                ];
                let in_ = raw[0] + raw[1] + raw[2] + raw[3] + raw[4] + raw[9];
                let out_ = raw[6] + raw[5] + raw[8] + raw[7] + raw[10];
                let err = qual_err(in_, out_);
                if err.abs() > max_err.abs() {
                    max_err = err;
                }
                let mu = if (p as usize) < ctx.pollutants.units.len() { ctx.pollutants.units[p as usize] } else { MassUnits::MgPerL };
                let v = &mut vals[p as usize];
                for k in 0..11 {
                    let x = raw[k];
                    let y = if mu == MassUnits::CountsPerL {
                        let z = LPER_FT3 * x;
                        if z > 0.0 {
                            z.log10() // LIBM
                        } else {
                            z
                        }
                    } else {
                        let mut cf = LPER_FT3 * cont_mass;
                        if mu == MassUnits::UgPerL {
                            cf /= 1000.0;
                        }
                        x * cf
                    };
                    v[k] = y;
                }
                v[11] = err;
            }
            if opt.rpt_continuity || max_err > 10.0 {
                const LABELS: [&[u8]; 12] = [
                    b"Dry Weather Inflow .......",
                    b"Wet Weather Inflow .......",
                    b"Groundwater Inflow .......",
                    b"RDII Inflow ..............",
                    b"External Inflow ..........",
                    b"External Outflow .........",
                    b"Flooding Loss ............",
                    b"Exfiltration Loss ........",
                    b"Mass Reacted .............",
                    b"Initial Stored Mass ......",
                    b"Final Stored Mass ........",
                    b"Continuity Error (%) .....",
                ];
                let mut p1 = 0;
                while p1 < np {
                    let p2 = cmin(p1 + 5, np);
                    cfprintf!(f, b"\n  **************************");
                    for p in p1..p2 {
                        cfprintf!(f, b"%14s", ctx.pollutant_names.name_of(p));
                    }
                    cfprintf!(f, b"\n  Quality Routing Continuity");
                    for p in p1..p2 {
                        let counts = (p as usize) < ctx.pollutants.units.len() && ctx.pollutants.units[p as usize] == MassUnits::CountsPerL;
                        cfprintf!(f, b"%14s", if counts { b"LogN" as &[u8] } else { load_word });
                    }
                    cfprintf!(f, b"\n  **************************");
                    for _p in p1..p2 {
                        cfprintf!(f, b"    ----------");
                    }
                    for k in 0..12 {
                        cfprintf!(f, b"\n  %s", LABELS[k]);
                        for p in p1..p2 {
                            cfprintf!(f, b"%14.3f", vals[p as usize][k]);
                        }
                    }
                    wr(f, b"");
                    wr(f, b"");
                    p1 += 5;
                }
            }
        }

        // ---- Reaction Species Routing Continuity (review B14) ----
        if ctx.n_nodes() > 0 && !opt.ignore_routing && !opt.ignore_quality && ctx.reactions.msx_ledger.armed && ctx.reactions.n_species() > 0 {
            let rx = &ctx.reactions;
            let led = &rx.msx_ledger;
            let nsp = rx.n_species();
            let uns = nsp as usize;
            let mut bulk: Vec<usize> = Vec::new();
            for s in 0..uns {
                if rx.species_is_wall[s] == 0 {
                    bulk.push(s);
                }
            }
            const LPER_FT3: f64 = 28.317;
            const LABELS: [&[u8]; 10] = [
                b"Inflow (all sources) .....",
                b"External Outflow .........",
                b"Flooding Loss ............",
                b"Exfiltration Loss ........",
                b"Mass Reacted .............",
                b"Clamped to Zero (added) ..",
                b"Initial Stored Mass ......",
                b"Final Stored Mass ........",
                b"Continuity Error (%) .....",
                b"",
            ];
            let mut vals: Vec<[f64; 9]> = vec![[0.0; 9]; uns];
            let mut unit_word: Vec<Vec<u8>> = vec![Vec::new(); uns];
            for s in 0..uns {
                let mut fin = atu(&led.final_dry, s);
                for j in 0..ctx.n_nodes() {
                    fin += atu(&rx.msx_node_conc, j as usize * uns + s) * ctx.nodes.volume[j as usize] + msx_lid_held_mass(ctx, j, s as i32);
                }
                for j in 0..ctx.n_links() {
                    fin += atu(&rx.msx_link_conc, j as usize * uns + s) * ctx.links.volume[j as usize];
                }
                let raw: [f64; 8] = [atu(&led.inflow, s), atu(&led.outflow, s), atu(&led.flooding, s), atu(&led.seepage, s), atu(&led.reacted, s), atu(&led.clamped, s), atu(&led.init, s), fin];
                let in_ = raw[0] + raw[5] + raw[6];
                let out_ = raw[1] + raw[2] + raw[3] + raw[4] + raw[7];
                let err = qual_err(in_, out_);
                let mut u = if s < rx.species_units.len() { rx.species_units[s].clone() } else { Vec::new() };
                for ch in u.iter_mut() {
                    *ch = ch.to_ascii_uppercase();
                }
                let mut cf = LPER_FT3;
                if u.starts_with(b"MG") {
                    cf *= cont_mass;
                    unit_word[s] = load_word.to_vec();
                } else if u.starts_with(b"UG") {
                    cf *= cont_mass / 1000.0;
                    unit_word[s] = load_word.to_vec();
                } else {
                    unit_word[s] = b"conc*L".to_vec();
                }
                for k in 0..8 {
                    vals[s][k] = raw[k] * cf;
                }
                vals[s][8] = err;
            }
            let mut s1 = 0usize;
            while s1 < bulk.len() {
                let s2 = cmin(s1 + 5, bulk.len());
                cfprintf!(f, b"\n  **************************");
                for s in s1..s2 {
                    cfprintf!(f, b"%14s", rx.species_name[bulk[s]]);
                }
                cfprintf!(f, b"\n  Reaction Species Routing  ");
                for s in s1..s2 {
                    cfprintf!(f, b"%14s", unit_word[bulk[s]]);
                }
                cfprintf!(f, b"\n  **************************");
                for _s in s1..s2 {
                    cfprintf!(f, b"    ----------");
                }
                for k in 0..9 {
                    cfprintf!(f, b"\n  %s", LABELS[k]);
                    for s in s1..s2 {
                        cfprintf!(f, b"%14.3f", vals[bulk[s]][k]);
                    }
                }
                wr(f, b"");
                wr(f, b"");
                s1 += 5;
            }
        }

        // ---- Highest Continuity Errors (legacy report_writeMaxStats) ----
        let dw_like = opt.routing_model == RoutingModel::Dynwave || opt.routing_model == RoutingModel::Fv;
        if opt.rpt_flowstats && dw_like && ctx.n_links() > 0 {
            let mut errors: Vec<(i32, f64)> = Vec::new();
            for j in 0..ctx.n_nodes() {
                let uj = j as usize;
                if ctx.nodes.degree[uj] <= 0 {
                    continue;
                }
                let in_vol = ctx.nodes.stat_total_inflow_vol[uj];
                if in_vol <= 0.1 {
                    continue;
                }
                let out_vol = ctx.nodes.stat_total_outflow_vol[uj] + ctx.nodes.volume[uj];
                let x = 100.0 * (1.0 - out_vol / in_vol);
                if x.abs() > 1.0 {
                    errors.push((j, x));
                }
            }
            // std::stable_sort by |error| descending
            errors.sort_by(|a, b| {
                if a.1.abs() > b.1.abs() {
                    std::cmp::Ordering::Less
                } else if b.1.abs() > a.1.abs() {
                    std::cmp::Ordering::Greater
                } else {
                    std::cmp::Ordering::Equal
                }
            });
            if !errors.is_empty() {
                wr(f, b"*************************");
                wr(f, b"Highest Continuity Errors");
                wr(f, b"*************************");
                let count = cmin(5, errors.len() as i32);
                for i in 0..count as usize {
                    cfprintf!(f, b"\n  Node %s (%.2f%%)", ctx.node_names.name_of(errors[i].0), errors[i].1);
                }
            }
        }

        wr(f, b"");
        wr(f, b"");

        // ---- Flow statistics sections (FLOWSTATS) ----
        if opt.rpt_flowstats {
            // Time-Step Critical Elements
            if opt.routing_model == RoutingModel::Dynwave && ctx.n_links() > 0 && opt.variable_step != 0.0 {
                wr(f, b"***************************");
                wr(f, b"Time-Step Critical Elements");
                wr(f, b"***************************");
                let mut k = 0;
                for i in 0..MAX_STATS {
                    let ms = &ctx.max_courant_crit[i];
                    if ms.index < 0 {
                        continue;
                    }
                    k += 1;
                    if ms.obj_type == 0 {
                        cfprintf!(f, b"\n  Node %s", ctx.node_names.name_of(ms.index));
                    } else {
                        cfprintf!(f, b"\n  Link %s", ctx.link_names.name_of(ms.index));
                    }
                    cfprintf!(f, b" (%.2f%%)", ms.value);
                }
                if k == 0 {
                    cfprintf!(f, b"\n  None");
                }
            }

            wr(f, b"");
            wr(f, b"");

            // Highest Flow Instability Indexes
            if ctx.n_links() > 0 {
                wr(f, b"********************************");
                wr(f, b"Highest Flow Instability Indexes");
                wr(f, b"********************************");
                if ctx.max_flow_turns[0].index <= 0 {
                    cfprintf!(f, b"\n  All links are stable.");
                } else {
                    for i in 0..MAX_STATS {
                        let ms = &ctx.max_flow_turns[i];
                        if ms.index < 0 {
                            continue;
                        }
                        cfprintf!(f, b"\n  Link %s (%.0f)", ctx.link_names.name_of(ms.index), ms.value);
                    }
                }
            }

            wr(f, b"");
            wr(f, b"");

            // Most Frequent Nonconverging Nodes
            if ctx.n_nodes() > 0 && ctx.n_links() > 0 && opt.routing_model == RoutingModel::Dynwave {
                wr(f, b"*********************************");
                wr(f, b"Most Frequent Nonconverging Nodes");
                wr(f, b"*********************************");
                if ctx.max_non_converged[0].index <= 0 || ctx.max_non_converged[0].value < 0.00005 {
                    wr(f, b"Convergence obtained at all time steps.");
                } else {
                    for i in 0..MAX_STATS {
                        let ms = &ctx.max_non_converged[i];
                        if ms.index < 0 || ms.value <= 0.0 {
                            continue;
                        }
                        cfprintf!(f, b"\n  Node %s (%.2f%%)", ctx.node_names.name_of(ms.index), 100.0 * ms.value);
                    }
                }
            }

            wr(f, b"");
            wr(f, b"");

            // Routing Time Step Summary
            if ctx.n_links() > 0 && ctx.routing_stats.n_steps > 0 {
                wr(f, b"*************************");
                wr(f, b"Routing Time Step Summary");
                wr(f, b"*************************");
                let rs = &ctx.routing_stats;
                let t_total = rs.steady_time + rs.sum_step;
                let steady_pct = cmin(if t_total > 0.0 { 100.0 * rs.steady_time / t_total } else { 0.0 }, 100.0);
                if rs.n_steps > 0 {
                    cfprintf!(f, b"\n  Minimum Time Step           :  %7.2f sec", if rs.min_step < 1.0e30 { rs.min_step } else { 0.0 });
                    cfprintf!(f, b"\n  Average Time Step           :  %7.2f sec", rs.avg_step());
                    cfprintf!(f, b"\n  Maximum Time Step           :  %7.2f sec", rs.max_step);
                    cfprintf!(f, b"\n  %% of Time in Steady State   :  %7.2f", steady_pct);
                    let is_fv = opt.routing_model == RoutingModel::Fv;
                    if is_fv {
                        cfprintf!(f, b"\n  Average Substeps per Step   :  %7.2f", rs.computed_avg_iterations());
                        cfprintf!(f, b"\n  %% of Steps Not Converging   :      n/a");
                    } else {
                        cfprintf!(f, b"\n  Average Iterations per Step :  %7.2f", rs.computed_avg_iterations());
                        cfprintf!(f, b"\n  %% of Steps Not Converging   :  %7.2f", rs.pct_non_converged());
                    }
                    let n_bins = crate::core::context::N_TIME_BINS;
                    let mut binned: i64 = 0;
                    for i in 0..n_bins {
                        binned += rs.step_counts[i] as i64;
                    }
                    if opt.routing_model == RoutingModel::Dynwave && opt.variable_step > 0.0 && binned > 0 {
                        cfprintf!(f, b"\n  Time Step Frequencies       :");
                        for i in 0..n_bins {
                            let pct = 100.0 * rs.step_counts[i] as f64 / binned as f64;
                            cfprintf!(f, b"\n     %6.3f - %6.3f sec      :  %7.2f %%", rs.step_intervals[i], rs.step_intervals[i + 1], pct);
                        }
                    }
                }
            }

            // FV Solver Statistics (skipped on the -1 sentinel: not an FV run)
            if ctx.routing_stats.fv_nsteps >= 0 {
                let rs = &ctx.routing_stats;
                wr(f, b"");
                wr(f, b"");
                wr(f, b"*********************");
                wr(f, b"FV Solver Statistics");
                wr(f, b"*********************");
                let srow = |f: &mut CFile, label: &[u8], v: i64| {
                    cfprintf!(f, b"\n  %s%14ld", label, v);
                };
                srow(f, b"Explicit Substeps ........", rs.fv_nsteps as i64);
                srow(f, b"Face Flux Evaluations ....", rs.fv_nflux as i64);
                cfprintf!(f, b"\n  Avg Substep (s) ..........%14.6f", rs.fv_avg_h);
                cfprintf!(f, b"\n  Min Substep (s) ..........%14.6f", rs.fv_min_h);
                cfprintf!(f, b"\n  Last Substep (s) .........%14.6f", rs.fv_last_h);
                srow(f, b"LTS Macro Cycles Fired ...", rs.fv_macro_cycles as i64);
                srow(f, b"LTS Macro Cycles Rejected ", rs.fv_macro_rejected as i64);
                if rs.fv_active_mean >= 0.0 {
                    cfprintf!(f, b"\n  Active Faces min/mean/max %8.1f /%5.1f /%5.1f  (%%)", 100.0 * rs.fv_active_min, 100.0 * rs.fv_active_mean, 100.0 * rs.fv_active_max);
                }
                {
                    let tot = rs.fv_dt_argmin_pressurized as i64 + rs.fv_dt_argmin_band as i64 + rs.fv_dt_argmin_free as i64 + rs.fv_dt_argmin_node as i64;
                    if tot > 0 {
                        let pct = |v: i64| 100.0 * v as f64 / tot as f64;
                        cfprintf!(f, b"\n  dt Argmin Pressurized (%%) %14.1f", pct(rs.fv_dt_argmin_pressurized as i64));
                        cfprintf!(f, b"\n  dt Argmin Taper Band (%%) .%14.1f", pct(rs.fv_dt_argmin_band as i64));
                        cfprintf!(f, b"\n  dt Argmin Free (%%) .......%14.1f", pct(rs.fv_dt_argmin_free as i64));
                        cfprintf!(f, b"\n  dt Argmin Node Bound (%%) .%14.1f", pct(rs.fv_dt_argmin_node as i64));
                    }
                }
                if rs.fv_n_tiers > 0 {
                    let mut total: i64 = 0;
                    let mut k = 0;
                    while k < rs.fv_n_tiers && k < 8 {
                        total += rs.fv_tier_cells[k as usize] as i64;
                        k += 1;
                    }
                    if total > 0 {
                        let mut k = 0;
                        while k < rs.fv_n_tiers && k < 8 {
                            cfprintf!(f, b"\n  LTS Tier %d Occupancy (%%) .%14.1f", k, 100.0 * rs.fv_tier_cells[k as usize] as f64 / total as f64);
                            k += 1;
                        }
                    }
                }
            }

            // Slot Storage Summary (FV slot program R0)
            {
                let mut sum_slot_dt = 0.0;
                let mut sum_vol_dt = 0.0;
                for j in 0..ctx.links.stat_slot_vol_dt.len() {
                    sum_slot_dt += ctx.links.stat_slot_vol_dt[j];
                    sum_vol_dt += ctx.links.stat_vol_dt[j];
                }
                if sum_slot_dt > 0.0 {
                    wr(f, b"");
                    wr(f, b"");
                    wr(f, b"*********************");
                    wr(f, b"Slot Storage Summary");
                    wr(f, b"*********************");
                    cfprintf!(f, b"\n  Run Slot Share (%%) .......%14.2f", 100.0 * sum_slot_dt / cmax(sum_vol_dt, 1e-30));
                    cfprintf!(f, b"\n  Peak Slot Share (%%) ......%14.2f", 100.0 * ctx.routing_stats.slot_peak_share);
                    cfprintf!(f, b"\n  Hours Share Above 1%% .....%14.2f", ctx.routing_stats.slot_time_above_s / 3600.0);
                    let mut header = false;
                    for j in 0..ctx.n_links() {
                        let uj = j as usize;
                        if ctx.links.stat_peak_slot_share[uj] < 0.01 {
                            continue;
                        }
                        if !header {
                            cfprintf!(f, b"\n\n  Link                    Peak Share    Run Share   Hrs >1%%\n  ----------------------------------------------------------");
                            header = true;
                        }
                        let run_share = if ctx.links.stat_vol_dt[uj] > 0.0 { ctx.links.stat_slot_vol_dt[uj] / ctx.links.stat_vol_dt[uj] } else { 0.0 };
                        cfprintf!(f, b"\n  %-20s %9.2f%%   %9.2f%% %9.2f", ctx.link_names.name_of(j), 100.0 * ctx.links.stat_peak_slot_share[uj], 100.0 * run_share, ctx.links.stat_time_slot_above[uj] / 3600.0);
                    }
                }
            }

            wr(f, b"");
            wr(f, b"");

            // Virtual Junction Summary
            if !ctx.vj_diag.node_idx.is_empty() {
                wr(f, b"***********************");
                wr(f, b"Virtual Junction Summary");
                wr(f, b"***********************");
                cfprintf!(f, "\n                        Upstream         Downstream         Momentum Residual (cfs\u{b7}ft/s)".as_bytes());
                cfprintf!(f, b"\n  Name                 Conduit          Conduit               Maximum          Mean");
                cfprintf!(f, b"\n  ------------------------------------------------------------------------------------");
                let mut any_fed = false;
                for r in 0..ctx.vj_diag.node_idx.len() {
                    let ni = ctx.vj_diag.node_idx[r];
                    let ju = ctx.vj_diag.up_link[r];
                    let jd = ctx.vj_diag.dn_link[r];
                    let n = ctx.vj_diag.resid_n[r];
                    let mean = if n > 0 { ctx.vj_diag.resid_sum[r] / n as f64 } else { 0.0 };
                    let uni = ni as usize;
                    let max_lat = if ni >= 0 && uni < ctx.nodes.stat_max_lat_inflow.len() { ctx.nodes.stat_max_lat_inflow[uni] } else { 0.0 };
                    let ju_name: &[u8] = if ju >= 0 { ctx.link_names.name_of(ju) } else { b"*" };
                    let jd_name: &[u8] = if jd >= 0 { ctx.link_names.name_of(jd) } else { b"*" };
                    if max_lat.abs() > 0.0 {
                        cfprintf!(f, b"\n  %-20s %-16s %-16s %13.6f %13.6f   lateral %.3f %s", ctx.node_names.name_of(ni), ju_name, jd_name, ctx.vj_diag.resid_max[r], mean, max_lat * qcf, fuw);
                        any_fed = true;
                    } else {
                        cfprintf!(f, b"\n  %-20s %-16s %-16s %13.6f %13.6f", ctx.node_names.name_of(ni), ju_name, jd_name, ctx.vj_diag.resid_max[r], mean);
                    }
                }
                if any_fed {
                    cfprintf!(f, b"\n  (lateral: maximum lateral inflow at the node; its residual includes the added zero-momentum mass)");
                }
                wr(f, b"");
                wr(f, b"");
            }
        } // end rpt_flowstats

        // ---- Conduit Time Step Summary ([REPORT] LINK_STEPS) ----
        if opt.rpt_link_steps && ctx.n_links() > 0 {
            let mode = ctx.routing_stats.lstep_mode;
            wr(f, b"");
            wr(f, b"");
            wr(f, b"*************************");
            wr(f, b"Conduit Time Step Summary");
            wr(f, b"*************************");
            if mode == 0 {
                wr(f, b"Not available for this routing method.");
            } else {
                const NB: usize = LinkData::N_LSTEP_BINS as usize;
                let rs = opt.routing_step;
                if mode == 1 {
                    wr(f, b"Dynamic wave: the CFL time step each conduit allows (Courant factor");
                    wr(f, b"applied, capped at the routing step); all conduits advance at the");
                    wr(f, b"global routing step. Dry conduits set no limit and are not sampled.");
                } else {
                    wr(f, b"Finite volume: the local time step each conduit took (its finest");
                    wr(f, b"cell under local time stepping, else the global substep).");
                }
                let mut lbl: Vec<Vec<u8>> = Vec::with_capacity(NB);
                for b in 0..NB {
                    let hi = rs / (1i32 << b) as f64;
                    let mut s = if b < NB - 1 { crate::cformat!(b"%.3g-%.3g", 0.5 * hi, hi) } else { crate::cformat!(b"<=%.3g", hi) };
                    s.truncate(23); // snprintf into char[24]
                    lbl.push(s);
                }
                let mut rule = b"\n  ".to_vec();
                rule.extend(std::iter::repeat_n(b'-', 20 + 27 + 1 + 13 * NB + 11));
                let mut band = b" % of Time with Step in Range (sec) ".to_vec();
                let bw = 13 * NB - 1;
                let pad = if bw > band.len() { bw - band.len() } else { 0 };
                let mut nb = vec![b'-'; pad / 2];
                nb.extend_from_slice(&band);
                nb.extend(std::iter::repeat_n(b'-', pad - pad / 2));
                band = nb;
                cfprintf!(f, b"\n");
                cfprintf!(f, b"%s", rule);
                cfprintf!(f, b"\n  %-20s %8s %8s %8s  %s  %9s", b"", b"Minimum", b"Average", b"Maximum", band, b"Percent");
                cfprintf!(f, b"\n  %-20s %8s %8s %8s ", b"Conduit", b"Step", b"Step", b"Step");
                for b in 0..NB {
                    cfprintf!(f, b" %12s", lbl[b]);
                }
                cfprintf!(f, b"  %9s", b"Converged");
                cfprintf!(f, b"%s", rule);

                let l = &ctx.links;
                for j in 0..ctx.n_links() {
                    let uj = j as usize;
                    if l.type_[uj] != LinkType::Conduit || l.rpt_flag[uj] == 0 {
                        continue;
                    }
                    cfprintf!(f, b"\n  %-20s", ctx.link_names.name_of(j));
                    let tt = l.stat_lstep_time[uj];
                    if tt > 0.0 {
                        cfprintf!(f, b" %8.4f %8.4f %8.4f ", l.stat_lstep_min[uj], l.stat_lstep_dt_time[uj] / tt, l.stat_lstep_max[uj]);
                        for b in 0..NB {
                            cfprintf!(f, b" %12.2f", 100.0 * l.stat_lstep_bin_time[uj * NB + b] / tt);
                        }
                    } else {
                        cfprintf!(f, b" %8s %8s %8s ", b"-", b"-", b"-");
                        for _b in 0..NB {
                            cfprintf!(f, b" %12s", b"-");
                        }
                    }
                    if mode == 1 && l.stat_conv_total[uj] > 0 {
                        cfprintf!(f, b"  %9.2f", 100.0 * l.stat_conv_steps[uj] as f64 / l.stat_conv_total[uj] as f64);
                    } else {
                        cfprintf!(f, b"  %9s", b"-");
                    }
                }
                cfprintf!(f, b"%s", rule);
            }
            wr(f, b"");
        }

        // ---- Rainfall File Summary (legacy report_writeRainStats) ----
        {
            let mut any_rain_file = false;
            for g in 0..ctx.n_gages() {
                if ctx.gages.source[g as usize] == RainSource::FileRain {
                    any_rain_file = true;
                    break;
                }
            }
            if any_rain_file {
                wr(f, b"");
                wr(f, b"*********************");
                wr(f, b"Rainfall File Summary");
                wr(f, b"*********************");
                cfprintf!(f, b"\n  Station    First        Last         Recording   Periods    Periods    Periods");
                cfprintf!(f, b"\n  ID         Date         Date         Frequency  w/Precip    Missing    Malfunc.");
                cfprintf!(f, b"\n  -------------------------------------------------------------------------------");
                for g in 0..ctx.n_gages() {
                    let ug = g as usize;
                    if ctx.gages.source[ug] != RainSource::FileRain {
                        continue;
                    }
                    let mut d1 = b"***********".to_vec();
                    let mut d2 = b"***********".to_vec();
                    if ctx.gages.file_first_date[ug] > 0.0 {
                        let (y, m, d) = datetime::decode_date(ctx.gages.file_first_date[ug]);
                        d1 = crate::cformat!(b"%02d/%02d/%04d", m, d, y);
                        d1.truncate(15);
                    }
                    if ctx.gages.file_last_date[ug] > 0.0 {
                        let (y, m, d) = datetime::decode_date(ctx.gages.file_last_date[ug]);
                        d2 = crate::cformat!(b"%02d/%02d/%04d", m, d, y);
                        d2.truncate(15);
                    }
                    let sta: &[u8] = if ctx.gages.station_id[ug].is_empty() { ctx.gage_names.name_of(g) } else { &ctx.gages.station_id[ug] };
                    cfprintf!(f, b"\n  %-10s %10s   %-10s   %5d min    %6ld     %6ld     %6ld", sta, d1, d2, ctx.gages.interval_sec[ug] / 60, ctx.gages.file_periods_precip[ug], 0i64, 0i64);
                    wr(f, b"");
                }
                wr(f, b"");
            }
        }

        // ---- Control Actions Taken (legacy report_writeRuleAction) ----
        if opt.rpt_controls {
            wr(f, b"*********************");
            wr(f, b"Control Actions Taken");
            wr(f, b"*********************");
            if !ctx.control_log.is_empty() {
                for entry in &ctx.control_log {
                    let (yr, mo, dy) = datetime::decode_date(entry.date);
                    let (hr, mn, sc) = datetime::decode_time(entry.date);
                    let mut datebuf = crate::cformat!(b"%02d/%02d/%04d", mo, dy, yr);
                    datebuf.truncate(15);
                    let mut timebuf = crate::cformat!(b"%02d:%02d:%02d", hr, mn, sc);
                    timebuf.truncate(15);
                    let rule_name: &[u8] = if entry.rule_idx >= 0 && (entry.rule_idx as usize) < ctx.control_rule_names.len() { &ctx.control_rule_names[entry.rule_idx as usize] } else { b"Rule?" };
                    cfprintf!(f, b"\n  %11s: %8s Link %s setting changed to %6.2f by Control %s", datebuf, timebuf, ctx.link_names.name_of(entry.link_idx), entry.new_setting, rule_name);
                }
                if ctx.control_log_dropped > 0 {
                    cfprintf!(f, b"\n  (%zu further control actions not listed: the log is capped at %zu entries)", ctx.control_log_dropped, K_MAX_CONTROL_LOG);
                }
            }
            wr(f, b"");
            wr(f, b"");
        }

        // ---- Subcatchment Runoff Summary (legacy writeSubcatchRunoff) ----
        let runoff_rpt = ctx.n_subcatches() > 0 && (!opt.ignore_rainfall || (ctx.snowpack_names.size() > 0 && !opt.ignore_snow_melt) || (ctx.aquifer_names.size() > 0 && !opt.ignore_groundwater));
        if runoff_rpt {
            wr(f, b"***************************");
            wr(f, b"Subcatchment Runoff Summary");
            wr(f, b"***************************");
            cfprintf!(f, b"\n");
            cfprintf!(
                f,
                b"\n  ------------------------------------------------------------------------------------------------------------------------------\n                            Total      Total      Total      Total     Imperv       Perv      Total       Total     Peak  Runoff\n                           Precip      Runon       Evap      Infil     Runoff     Runoff     Runoff      Runoff   Runoff   Coeff\n  Subcatchment                 %2s         %2s         %2s         %2s         %2s         %2s         %2s    %8s      %3s\n  ------------------------------------------------------------------------------------------------------------------------------",
                depth_word,
                depth_word,
                depth_word,
                depth_word,
                depth_word,
                depth_word,
                depth_word,
                vol_word,
                fuw
            );
            for j in 0..ctx.n_subcatches() {
                let uj = j as usize;
                let a = ctx.subcatches.area[uj];
                if a <= 0.0 {
                    continue;
                }
                let area_ft2 = a / ucf::ucf(ucf::LANDAREA, opt);
                let sc = &ctx.subcatches;
                let precip_in = if area_ft2 > 0.0 { sc.stat_precip_vol[uj] * depth_vcf / area_ft2 } else { 0.0 };
                let evap_in = if area_ft2 > 0.0 { sc.stat_evap_vol[uj] * depth_vcf / area_ft2 } else { 0.0 };
                let infil_in = if area_ft2 > 0.0 { sc.stat_infil_vol[uj] * depth_vcf / area_ft2 } else { 0.0 };
                let imperv_in = if area_ft2 > 0.0 { sc.stat_imperv_vol[uj] * depth_vcf / area_ft2 } else { 0.0 };
                let perv_in = if area_ft2 > 0.0 { sc.stat_perv_vol[uj] * depth_vcf / area_ft2 } else { 0.0 };
                let runoff_in = if area_ft2 > 0.0 { sc.stat_runoff_vol[uj] * depth_vcf / area_ft2 } else { 0.0 };
                let runoff_vol = sc.stat_runoff_vol[uj] * vcf;
                let peak = sc.stat_max_runoff[uj] * qcf;
                let r = sc.stat_precip_vol[uj] + sc.stat_runon_vol[uj];
                let coeff = if r > 0.0 { sc.stat_runoff_vol[uj] / r } else { 0.0 };
                cfprintf!(f, b"\n  %-20s", ctx.subcatch_names.name_of(j));
                cfprintf!(f, b" %10.2f", precip_in);
                cfprintf!(f, b" %10.2f", sc.stat_runon_vol[uj] * depth_vcf / area_ft2);
                cfprintf!(f, b" %10.2f", evap_in);
                cfprintf!(f, b" %10.2f", infil_in);
                cfprintf!(f, b" %10.2f", imperv_in);
                cfprintf!(f, b" %10.2f", perv_in);
                cfprintf!(f, b" %10.2f", runoff_in);
                cfprintf!(f, b"%12.2f", runoff_vol);
                cfprintf!(f, b" %8.2f", peak);
                cfprintf!(f, b"%8.3f", coeff);
            }
            wr(f, b"");
        }

        wr(f, b"");
        wr(f, b"");

        // ---- Subcatchment Washoff Summary (legacy writeSubcatchLoads) ----
        let ms_wo = &ctx.reactions.surface;
        let nm_wo = if ms_wo.active() && ctx.n_subcatches() > 0 { ms_wo.n_species } else { 0 };
        if runoff_rpt && (ctx.n_pollutants() > 0 || nm_wo > 0) && !opt.ignore_quality {
            let ns = ctx.n_subcatches();
            let np_p = ctx.n_pollutants();
            let np = np_p + nm_wo;
            let col_name = |p: i32| -> &[u8] {
                if p < np_p {
                    ctx.pollutant_names.name_of(p)
                } else {
                    &ctx.reactions.species_name[(p - np_p) as usize]
                }
            };
            wr(f, b"****************************");
            wr(f, b"Subcatchment Washoff Summary");
            wr(f, b"****************************");
            cfprintf!(f, b"\n");
            cfprintf!(f, b" \n  ----------------------------------");
            for _p in 1..np {
                cfprintf!(f, b"--------------");
            }
            cfprintf!(f, b" \n                                ");
            for p in 0..np {
                cfprintf!(f, b"%14s", col_name(p));
            }
            cfprintf!(f, b" \n  Subcatchment                  ");
            for p in 0..np {
                if p >= np_p {
                    let mcf = ms_wo.mcf[(p - np_p) as usize];
                    cfprintf!(f, b"%14s", if mcf == 1.0 { b"mass" as &[u8] } else if si_report { b"kg" as &[u8] } else { b"lbs" as &[u8] });
                    continue;
                }
                let up = p as usize;
                let mu = if up < ctx.pollutants.units.len() { ctx.pollutants.units[up] } else { MassUnits::MgPerL };
                cfprintf!(f, b"%14s", if mu == MassUnits::CountsPerL { b"LogN" as &[u8] } else { load_word });
            }
            cfprintf!(f, b" \n  ----------------------------------");
            for _p in 1..np {
                cfprintf!(f, b"--------------");
            }
            let mut sys_loads = vec![0.0f64; np as usize];
            for j in 0..ns {
                let uj = j as usize;
                cfprintf!(f, b"\n  %-30s", ctx.subcatch_names.name_of(j));
                for p in 0..np {
                    let up = p as usize;
                    let load_lbs = if p < np_p {
                        let idx = uj * np_p as usize + up;
                        if idx < ctx.subcatches.total_load.len() {
                            ctx.subcatches.total_load[idx]
                        } else {
                            0.0
                        }
                    } else {
                        let k = ms_wo.sidx(j, p - np_p);
                        if k < ms_wo.led_subcatch_load.len() {
                            ms_wo.led_subcatch_load[k]
                        } else {
                            0.0
                        }
                    };
                    sys_loads[up] += load_lbs;
                    cfprintf!(f, b"%14.3f", if is_count(p) && load_lbs > 0.0 { load_lbs.log10() } else { load_lbs }); // LIBM
                }
            }
            cfprintf!(f, b" \n  ----------------------------------");
            for _p in 1..np {
                cfprintf!(f, b"--------------");
            }
            cfprintf!(f, b"\n  System                        ");
            for p in 0..np {
                let x = sys_loads[p as usize];
                cfprintf!(f, b"%14.3f", if is_count(p) && x > 0.0 { x.log10() } else { x }); // LIBM
            }
            wr(f, b"");
            wr(f, b"");
        }

        // ---- OPENSWMM_HAS_2D: 2D Surface Washoff Summary (surface quality, phase 2b) ----
        #[cfg(feature = "twod")]
        if let Some(sq) = &ctx.twod_io.surface_quality {
            if sq.borrow().active() && !opt.ignore_quality {
                panic!("STUB: the 2D Surface Washoff Summary (DefaultReportPlugin::write_results; phase 2b)");
            }
        }

        // ---- OPENSWMM_HAS_2D: 2D Outfall Tailwater Summary ----
        // Each outfall coupled to the 2D surface: how long the surface set its stage (engaged)
        // and drowned its crown (submerged), the highest stage, and the reverse volume.
        #[cfg(feature = "twod")]
        if let Some(r) = router {
            let ex = r.outfall_exchange();
            let mut header = false;
            for j in 0..ctx.n_nodes() {
                let (engaged, submerged, max_stage) = match crate::twod::surface_router_2d_end::outfall_tailwater(ex, j) {
                    Some(t) => t,
                    None => continue,
                };
                if !header {
                    wr(f, b"****************************");
                    wr(f, b"2D Outfall Tailwater Summary");
                    wr(f, b"****************************");
                    cfprintf!(f, b"\n");
                    cfprintf!(f, b"\n  ------------------------------------------------------------------------");
                    cfprintf!(f, b"\n                          Hours        Hours          Max       Reverse");
                    cfprintf!(f, b"\n                       Tailwater    Submerged    Tailwater        Volume");
                    cfprintf!(f, b"\n  Outfall Node            Engaged      (crown)    %9s    %10s", du.length_word, vol_word);
                    cfprintf!(f, b"\n  ------------------------------------------------------------------------");
                    header = true;
                }
                let uj = j as usize;
                let reverse_ft3 = if uj < ctx.outfall_exports.reverse_m3.len() { ctx.outfall_exports.reverse_m3[uj] / 0.028316846592 } else { 0.0 };
                cfprintf!(f, b"\n  %-20s %10.2f   %10.2f", ctx.node_names.name_of(j), engaged / 3600.0, submerged / 3600.0);
                if max_stage > -1.0e29 {
                    cfprintf!(f, b"   %10.3f", max_stage * len_ucf);
                } else {
                    cfprintf!(f, b"   %10s", b"-");
                }
                cfprintf!(f, b"    %10.3f", reverse_ft3 * vcf);
            }
            if header {
                wr(f, b"");
                wr(f, b"");
            }
        }

        // ---- Groundwater Summary (legacy writeGroundwater) ----
        if runoff_rpt && has_gw && !opt.ignore_groundwater {
            wr(f, b"*******************");
            wr(f, b"Groundwater Summary");
            wr(f, b"*******************");
            cfprintf!(f, b"\n");
            cfprintf!(f, b"\n  -----------------------------------------------------------------------------------------------------\n                                            Total    Total  Maximum  Average  Average    Final    Final\n                          Total    Total    Lower  Lateral  Lateral    Upper    Water    Upper    Water\n                          Infil     Evap  Seepage  Outflow  Outflow   Moist.    Table   Moist.    Table");
            cfprintf!(
                f,
                if si_report {
                    b"\n  Subcatchment               mm       mm       mm       mm      %3s                 m                 m" as &[u8]
                } else {
                    b"\n  Subcatchment               in       in       in       in      %3s                ft                ft"
                },
                fuw
            );
            cfprintf!(f, b"\n  -----------------------------------------------------------------------------------------------------");
            for j in 0..ctx.n_subcatches() {
                let uj = j as usize;
                let sc = &ctx.subcatches;
                if sc.gw_aquifer[uj] < 0 {
                    continue;
                }
                let area_ft2 = sc.area[uj] * land_area_to_ft2(fu);
                let infil_in = if area_ft2 > 0.0 { sc.stat_gw_infil_vol[uj] / area_ft2 * depth_vcf } else { 0.0 };
                let evap_in = if area_ft2 > 0.0 { (sc.stat_gw_upper_evap_vol[uj] + sc.stat_gw_lower_evap_vol[uj]) / area_ft2 * depth_vcf } else { 0.0 };
                let seep_in = if area_ft2 > 0.0 { sc.stat_gw_deep_perc_vol[uj] / area_ft2 * depth_vcf } else { 0.0 };
                let lat_in = if area_ft2 > 0.0 { sc.stat_gw_flow_vol[uj] / area_ft2 * depth_vcf } else { 0.0 };
                let max_flow = sc.stat_gw_max_flow[uj] * qcf;
                let seconds = sc.stat_gw_time[uj];
                let avg_theta = if seconds > 0.0 { sc.stat_gw_sum_theta[uj] / seconds } else { 0.0 };
                let avg_depth = if seconds > 0.0 { sc.stat_gw_sum_depth[uj] * len_ucf / seconds } else { 0.0 };
                let fin_theta = sc.stat_gw_final_theta[uj];
                let fin_depth = sc.stat_gw_final_depth[uj] * len_ucf;
                cfprintf!(f, b"\n  %-20s %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f", ctx.subcatch_names.name_of(j), infil_in, evap_in, seep_in, lat_in, max_flow, avg_theta, avg_depth, fin_theta, fin_depth);
            }
            wr(f, b"");
        }

        wr(f, b"");
        wr(f, b"");

        // ---- LID Performance Summary (legacy writeLidPerformance) ----
        if runoff_rpt && ctx.lid_usage.count() > 0 {
            let n_usage = ctx.lid_usage.count();
            let mut has_lids = false;
            for j in 0..n_usage {
                let uj = j as usize;
                if uj < ctx.lid_usage.wb_inflow.len() && (ctx.lid_usage.wb_inflow[uj] > 0.0 || ctx.lid_usage.wb_drain_flow[uj] > 0.0 || ctx.lid_usage.wb_surf_flow[uj] > 0.0) {
                    has_lids = true;
                    break;
                }
            }
            wr(f, b"********************");
            wr(f, b"LID Performance Summary");
            wr(f, b"********************");
            cfprintf!(f, b"\n");
            if !has_lids {
                wr(f, b"");
                wr(f, b"No LID performance data.");
            } else {
                cfprintf!(f, b"\n  -----------------------------------------------------------------\n                          Total    Evap   Infil  Surface   Drain\n                         Inflow    Loss    Loss    Runoff    Flow");
                cfprintf!(f, if si_report { b"\n  Control Group    Subcatch      mm      mm      mm       mm      mm" as &[u8] } else { b"\n  Control Group    Subcatch      in      in      in       in      in" });
                cfprintf!(f, b"\n  -----------------------------------------------------------------");
                for j in 0..n_usage {
                    let uj = j as usize;
                    let li = ctx.lid_usage.lid_index[uj];
                    let si = ctx.lid_usage.subcatch_index[uj];
                    let safe = |v: &[f64]| -> f64 {
                        if uj < v.len() {
                            v[uj] * depth_vcf
                        } else {
                            0.0
                        }
                    };
                    cfprintf!(f, b"\n  %-16s %-12s", ctx.lid_names.name_of(li), ctx.subcatch_names.name_of(si));
                    cfprintf!(f, b" %7.2f %7.2f %7.2f %8.2f %7.2f", safe(&ctx.lid_usage.wb_inflow), safe(&ctx.lid_usage.wb_evap), safe(&ctx.lid_usage.wb_infil), safe(&ctx.lid_usage.wb_surf_flow), safe(&ctx.lid_usage.wb_drain_flow));
                }
            }
            wr(f, b"");
            wr(f, b"");
        }

        // ---- Node Depth Summary (legacy writeNodeDepths) ----
        let route_rpt = ctx.n_links() > 0 && !opt.ignore_routing;
        if route_rpt {
            wr(f, b"******************");
            wr(f, b"Node Depth Summary");
            wr(f, b"******************");
            cfprintf!(f, b"\n");
            cfprintf!(f, b"\n  ---------------------------------------------------------------------------------\n                                 Average  Maximum  Maximum  Time of Max    Reported\n                                   Depth    Depth      HGL   Occurrence   Max Depth");
            cfprintf!(f, if si_report { b"\n  Node                 Type       Meters   Meters   Meters  days hr:min      Meters" as &[u8] } else { b"\n  Node                 Type         Feet     Feet     Feet  days hr:min        Feet" });
            cfprintf!(f, b"\n  ---------------------------------------------------------------------------------");
            let report_steps = cmax(ctx.routing_stats.report_steps as i64, 1i64);
            for j in 0..ctx.n_nodes() {
                let uj = j as usize;
                let nt = ctx.nodes.type_[uj] as i32;
                let max_d_int = ctx.nodes.stat_max_depth[uj];
                let avg_d = ctx.nodes.stat_sum_depth[uj] / report_steps as f64 * len_ucf;
                let max_d = max_d_int * len_ucf;
                let max_hgl = (ctx.nodes.invert_elev[uj] + max_d_int) * len_ucf;
                let rpt_max = ctx.nodes.stat_max_rpt_depth[uj];
                let (days, hrs, mins) = elapsed_to_parts(ctx.nodes.stat_max_depth_date[uj], ctx.options.report_start);
                cfprintf!(f, b"\n  %-20s", ctx.node_names.name_of(j));
                cfprintf!(f, b" %-9s", nt_str(nt));
                cfprintf!(f, b" %7.2f  %7.2f  %7.2f  %4d  %02d:%02d  %10.2f", avg_d, max_d, max_hgl, days, hrs, mins, rpt_max);
            }
        }

        wr(f, b"");
        wr(f, b"");

        // ---- Node Inflow Summary (legacy writeNodeFlows) ----
        if route_rpt {
            wr(f, b"*******************");
            wr(f, b"Node Inflow Summary");
            wr(f, b"*******************");
            cfprintf!(f, b"\n");
            cfprintf!(
                f,
                b"\n  -------------------------------------------------------------------------------------------------\n                                  Maximum  Maximum                  Lateral       Total        Flow\n                                  Lateral    Total  Time of Max      Inflow      Inflow     Balance\n                                   Inflow   Inflow   Occurrence      Volume      Volume       Error\n  Node                 Type           %3s      %3s  days hr:min    %8s    %8s     Percent\n  -------------------------------------------------------------------------------------------------",
                fuw,
                fuw,
                vol_word,
                vol_word
            );
            for j in 0..ctx.n_nodes() {
                let uj = j as usize;
                let nt = ctx.nodes.type_[uj] as i32;
                let max_lat = ctx.nodes.stat_max_lat_inflow[uj] * qcf;
                let max_tot = ctx.nodes.stat_max_total_inflow[uj] * qcf;
                let vol_lat = ctx.nodes.stat_lat_inflow_vol[uj] * vcf;
                let vol_tot = ctx.nodes.stat_total_inflow_vol[uj] * vcf;
                let (days, hrs, mins) = elapsed_to_parts(ctx.nodes.stat_max_inflow_date[uj], ctx.options.report_start);
                cfprintf!(f, b"\n  %-20s", ctx.node_names.name_of(j));
                cfprintf!(f, b" %-9s", nt_str(nt));
                cfprintf!(f, ff, max_lat);
                cfprintf!(f, ff, max_tot);
                cfprintf!(f, b"  %4d  %02d:%02d", days, hrs, mins);
                cfprintf!(f, b"%12.3g", vol_lat);
                cfprintf!(f, b"%12.3g", vol_tot);
                let in_v = ctx.nodes.stat_total_inflow_vol[uj];
                let out_v = ctx.nodes.stat_total_outflow_vol[uj] + ctx.nodes.volume[uj];
                if out_v.abs() < 1.0 {
                    cfprintf!(f, b"%12.3f %s", (in_v - out_v) * vcf * 1.0e6, if si_report { b"ltr" as &[u8] } else { b"gal" });
                } else {
                    cfprintf!(f, b"%12.3f", (in_v - out_v) / out_v * 100.0);
                }
            }
        }

        wr(f, b"");
        wr(f, b"");

        // ---- Node Surcharge Summary (legacy writeNodeSurcharge), DYNWAVE only ----
        if route_rpt && opt.routing_model as i32 == 2 {
            wr(f, b"**********************");
            wr(f, b"Node Surcharge Summary");
            wr(f, b"**********************");
            let mut n_written = 0;
            for j in 0..ctx.n_nodes() {
                let uj = j as usize;
                if ctx.nodes.type_[uj] == NodeType::Outfall {
                    continue;
                }
                if ctx.nodes.stat_time_surcharged[uj] == 0.0 {
                    continue;
                }
                let t = cmax(0.01, ctx.nodes.stat_time_surcharged[uj] / 3600.0);
                if n_written == 0 {
                    wr(f, b"");
                    wr(f, b"Surcharging occurs when water rises above the top of the highest conduit.");
                    cfprintf!(f, b"\n  ---------------------------------------------------------------------\n                                               Max. Height   Min. Depth\n                                   Hours       Above Crown    Below Rim");
                    cfprintf!(f, if si_report { b"\n  Node                 Type      Surcharged         Meters       Meters" as &[u8] } else { b"\n  Node                 Type      Surcharged           Feet         Feet" });
                    cfprintf!(f, b"\n  ---------------------------------------------------------------------");
                    n_written = 1;
                }
                let maxd = ctx.nodes.stat_max_depth[uj];
                let d1 = cmax(0.0, maxd + ctx.nodes.invert_elev[uj] - ctx.nodes.crown_elev[uj]);
                let d2 = cmax(0.0, ctx.nodes.full_depth[uj] - maxd);
                cfprintf!(f, b"\n  %-20s", ctx.node_names.name_of(j));
                cfprintf!(f, b" %-9s", nt_str(ctx.nodes.type_[uj] as i32));
                cfprintf!(f, b"  %9.2f      %9.3f    %9.3f", t, d1 * len_ucf, d2 * len_ucf);
            }
            if n_written == 0 {
                wr(f, b"");
                wr(f, b"No nodes were surcharged.");
            }
        }

        wr(f, b"");
        wr(f, b"");

        // ---- Node Flooding Summary (legacy writeNodeFlooding) ----
        if route_rpt {
            let mut any_flooding = false;
            for j in 0..ctx.n_nodes() {
                if ctx.nodes.stat_time_flooded[j as usize] > 0.0 {
                    any_flooding = true;
                    break;
                }
            }
            wr(f, b"*********************");
            wr(f, b"Node Flooding Summary");
            wr(f, b"*********************");
            if !any_flooding {
                wr(f, b"");
                wr(f, b"No nodes were flooded.");
            } else {
                cfprintf!(f, b"\n");
                cfprintf!(f, b"\n  Flooding refers to all water that overflows a node, whether it ponds or not.");
                cfprintf!(
                    f,
                    b"\n  --------------------------------------------------------------------------\n                                                             Total   Maximum\n                                 Maximum   Time of Max       Flood    Ponded\n                        Hours       Rate    Occurrence      Volume     Depth\n  Node                 Flooded       %3s   days hr:min    %8s    %6s\n  --------------------------------------------------------------------------",
                    fuw,
                    vol_word,
                    len_word
                );
                for j in 0..ctx.n_nodes() {
                    let uj = j as usize;
                    if ctx.nodes.stat_time_flooded[uj] <= 0.0 {
                        continue;
                    }
                    let hours = cmax(0.01, ctx.nodes.stat_time_flooded[uj] / 3600.0);
                    let max_rate = ctx.nodes.stat_max_overflow[uj] * qcf;
                    let vol_mgal = ctx.nodes.stat_vol_flooded[uj] * vcf;
                    let (days, hrs, mins) = elapsed_to_parts(ctx.nodes.stat_max_overflow_date[uj], ctx.options.report_start);
                    let mut ponded = 0.0;
                    if opt.routing_model as i32 == 2 {
                        ponded = ctx.nodes.stat_max_depth[uj] - ctx.nodes.full_depth[uj];
                        ponded = cmax(ponded, 0.0);
                    }
                    cfprintf!(f, b"\n  %-20s", ctx.node_names.name_of(j));
                    cfprintf!(f, b" %7.2f ", hours);
                    cfprintf!(f, ff, max_rate);
                    cfprintf!(f, b"   %4d  %02d:%02d", days, hrs, mins);
                    cfprintf!(f, b"%12.3f", vol_mgal);
                    cfprintf!(f, b" %9.3f", ponded * len_ucf);
                }
            }
        }

        wr(f, b"");
        wr(f, b"");

        // ---- LID Storage Moisture Profile (final per-cell profile) ----
        for r in 0..ctx.node_subtypes.storages.count() {
            let state = &ctx.node_subtypes.storages.lid_state[r as usize];
            if state.cells.is_empty() {
                continue;
            }
            let n = ctx.node_subtypes.storages.node_idx[r as usize];
            cfprintf!(f, b"\n\n  LID Storage Moisture Profile: %s\n  Layer      Bottom         Top     Moisture\n", ctx.node_names.name_of(n));
            for c in &state.cells {
                let wet = cclamp((ctx.nodes.depth[n as usize] - c.bottom) / (c.top - c.bottom), 0.0, 1.0);
                cfprintf!(f, b"  %5d %11.4f %11.4f %12.6f\n", c.layer, c.bottom * len_ucf, c.top * len_ucf, c.theta + (c.porosity - c.theta) * wet);
            }
        }

        // ---- Storage Volume Summary (legacy writeStorageVolumes) ----
        if route_rpt {
            let mut any_storage = false;
            for j in 0..ctx.n_nodes() {
                if ctx.nodes.type_[j as usize] == NodeType::Storage {
                    any_storage = true;
                    break;
                }
            }
            if any_storage {
                wr(f, b"**********************");
                wr(f, b"Storage Volume Summary");
                wr(f, b"**********************");
                cfprintf!(f, b"\n");
                cfprintf!(f, b"\n  ------------------------------------------------------------------------------------------------\n                         Average    Avg   Evap  Exfil     Maximum    Max    Time of Max    Maximum\n                          Volume   Pcnt   Pcnt   Pcnt      Volume   Pcnt     Occurrence    Outflow");
                cfprintf!(
                    f,
                    if si_report {
                        b"\n  Storage Unit           1000 m\xB3   Full   Loss   Loss     1000 m\xB3   Full    days hr:min        %3s" as &[u8]
                    } else {
                        b"\n  Storage Unit          1000 ft\xB3   Full   Loss   Loss    1000 ft\xB3   Full    days hr:min        %3s"
                    },
                    fuw
                );
                cfprintf!(f, b"\n  ------------------------------------------------------------------------------------------------");
                let report_steps = cmax(ctx.routing_stats.report_steps as i64, 1i64);
                for j in 0..ctx.n_nodes() {
                    let uj = j as usize;
                    if ctx.nodes.type_[uj] != NodeType::Storage {
                        continue;
                    }
                    let full_vol = ctx.nodes.full_volume[uj];
                    let avg_vol_ft3 = ctx.nodes.stat_sum_volume[uj] / report_steps as f64;
                    let max_vol_ft3 = ctx.nodes.stat_storage_max_volume[uj];
                    let pct_avg = if full_vol > 0.0 { avg_vol_ft3 / full_vol * 100.0 } else { 0.0 };
                    let pct_max = if full_vol > 0.0 { max_vol_ft3 / full_vol * 100.0 } else { 0.0 };
                    let avg_vol = avg_vol_ft3 * svol_ucf / 1000.0;
                    let max_vol = max_vol_ft3 * svol_ucf / 1000.0;
                    let (days, hrs, mins) = elapsed_to_parts(ctx.nodes.stat_storage_max_volume_date[uj], ctx.options.report_start);
                    let max_outflow = ctx.nodes.stat_storage_max_outflow[uj] * qcf;
                    cfprintf!(f, b"\n  %-20s", ctx.node_names.name_of(j));
                    cfprintf!(f, b"%10.3f  %5.1f  %5.1f  %5.1f  %10.3f  %5.1f", avg_vol, pct_avg, 0.0, 0.0, max_vol, pct_max);
                    cfprintf!(f, b"    %4d  %02d:%02d  ", days, hrs, mins);
                    cfprintf!(f, ff, max_outflow);
                }
                wr(f, b"");
                wr(f, b"");
            }
        }

        // ---- Outfall Loading Summary (legacy writeOutfallLoads) ----
        let mut any_outfall = false;
        let mut j = 0;
        while j < ctx.n_nodes() && !any_outfall {
            any_outfall = ctx.nodes.type_[j as usize] == NodeType::Outfall;
            j += 1;
        }
        if route_rpt && any_outfall {
            let np = ctx.n_pollutants();
            wr(f, b"***********************");
            wr(f, b"Outfall Loading Summary");
            wr(f, b"***********************");
            cfprintf!(f, b"\n");
            cfprintf!(f, b" \n  -----------------------------------------------------------");
            for _p in 0..np {
                cfprintf!(f, b"--------------");
            }
            cfprintf!(f, b" \n                         Flow       Avg       Max       Total");
            for _p in 0..np {
                cfprintf!(f, b"         Total");
            }
            cfprintf!(f, b" \n                         Freq      Flow      Flow      Volume");
            for p in 0..np {
                cfprintf!(f, b"%14s", ctx.pollutant_names.name_of(p));
            }
            cfprintf!(f, b" \n  Outfall Node           Pcnt       %3s       %3s    %8s", fuw, fuw, vol_word);
            for p in 0..np {
                cfprintf!(f, b"%14s", if ctx.pollutants.units[p as usize] == MassUnits::CountsPerL { b"LogN" as &[u8] } else if si_report { b"kg" as &[u8] } else { b"lbs" as &[u8] });
            }
            cfprintf!(f, b" \n  -----------------------------------------------------------");
            for _p in 0..np {
                cfprintf!(f, b"--------------");
            }

            let mut sys_vol = 0.0;
            let mut sys_flow_sum = 0.0;
            let mut sys_max_flow = 0.0;
            let mut sys_freq_sum = 0.0;
            let mut outfall_count = 0;
            let mut pol_totals = vec![0.0f64; np as usize];
            let total_steps = cmax(ctx.routing_stats.report_steps as i64, 1i64);

            for j in 0..ctx.n_nodes() {
                let uj = j as usize;
                if ctx.nodes.type_[uj] != NodeType::Outfall {
                    continue;
                }
                outfall_count += 1;
                let periods = ctx.nodes.stat_outfall_periods[uj] as i64;
                let avg_flow = if periods > 0 { ctx.nodes.stat_outfall_avg_flow[uj] * qcf / periods as f64 } else { 0.0 };
                let max_flow = ctx.nodes.stat_outfall_max_flow[uj] * qcf;
                let freq = 100.0 * periods as f64 / total_steps as f64;
                let vol = ctx.nodes.stat_total_inflow_vol[uj] * vcf;
                sys_freq_sum += freq;
                sys_flow_sum += avg_flow;
                sys_max_flow += max_flow;
                sys_vol += vol;
                cfprintf!(f, b"\n  %-20s", ctx.node_names.name_of(j));
                cfprintf!(f, b"%7.2f", freq);
                cfprintf!(f, b" ");
                cfprintf!(f, ff, avg_flow);
                cfprintf!(f, b" ");
                cfprintf!(f, ff, max_flow);
                cfprintf!(f, b"%12.3f", vol);
                let base = uj * np as usize;
                for p in 0..np {
                    let idx = base + p as usize;
                    let mu = ctx.pollutants.units[p as usize];
                    let mcf = if mu == MassUnits::MgPerL {
                        cont_mass
                    } else if mu == MassUnits::UgPerL {
                        cont_mass / 1000.0
                    } else {
                        1.0
                    };
                    let mut load = (if idx < ctx.nodes.stat_total_load.len() { ctx.nodes.stat_total_load[idx] } else { 0.0 }) * 28.317 * mcf;
                    pol_totals[p as usize] += load;
                    if mu == MassUnits::CountsPerL && load > 0.0 {
                        load = load.log10(); // LIBM
                    }
                    cfprintf!(f, b"%14.3f", load);
                }
            }

            cfprintf!(f, b" \n  -----------------------------------------------------------");
            for _p in 0..np {
                cfprintf!(f, b"--------------");
            }
            cfprintf!(f, b"\n  System              ");
            let sys_freq = if outfall_count > 0 { sys_freq_sum / outfall_count as f64 } else { 0.0 };
            cfprintf!(f, b"%7.2f ", sys_freq);
            cfprintf!(f, ff, sys_flow_sum);
            cfprintf!(f, b" ");
            let _ = sys_max_flow;
            cfprintf!(f, ff, ctx.routing_stats.max_outfall_flow * qcf);
            cfprintf!(f, b"%12.3f", sys_vol);
            for p in 0..np {
                let mut x = pol_totals[p as usize];
                if ctx.pollutants.units[p as usize] == MassUnits::CountsPerL && x > 0.0 {
                    x = x.log10(); // LIBM
                }
                cfprintf!(f, b"%14.3f", x);
            }
        }

        wr(f, b"");
        wr(f, b"");

        // ---- Link Flow Summary (legacy writeLinkFlows) ----
        if route_rpt {
            wr(f, b"********************");
            wr(f, b"Link Flow Summary");
            wr(f, b"********************");
            cfprintf!(f, b"\n");
            cfprintf!(
                f,
                b"\n  -----------------------------------------------------------------------------\n                                 Maximum  Time of Max   Maximum    Max/    Max/\n                                  |Flow|   Occurrence   |Veloc|    Full    Full\n  Link                 Type          %3s  days hr:min    %6s    Flow   Depth\n  -----------------------------------------------------------------------------",
                fuw,
                if si_report { b"m/sec" as &[u8] } else { b"ft/sec" }
            );
            for j in 0..ctx.n_links() {
                let uj = j as usize;
                let lt = ctx.links.type_[uj] as i32;
                let mf = ctx.links.stat_max_flow[uj] * qcf;
                let mv = ctx.links.stat_max_veloc[uj] * len_ucf;
                let fill = ctx.links.stat_max_filling[uj];
                let mut flow_ratio = 0.0;
                if lt == LinkType::Conduit as i32 {
                    let cr = ctx.link_subtypes.conduit_row(j);
                    let cd = &ctx.link_subtypes.conduits;
                    let mf_cfs = ctx.links.stat_max_flow[uj];
                    let qf = if cr >= 0 { cd.q_full[cr as usize] } else { 0.0 };
                    let barrels = cmax(if cr >= 0 { cd.barrels[cr as usize] } else { 1 }, 1);
                    flow_ratio = if qf > 0.0 { mf_cfs / qf / barrels as f64 } else { 0.0 };
                }
                let (days, hrs, mins) = elapsed_to_parts(ctx.links.stat_max_flow_date[uj], ctx.options.report_start);
                let shape = ctx.links.xsect_shape[uj];
                let pump = lt == LinkType::Pump as i32;
                let dummy = !pump && (shape == XsectShape::Dummy || lt == LinkType::Outlet as i32);
                cfprintf!(f, b"\n  %-20s", ctx.link_names.name_of(j));
                if dummy {
                    cfprintf!(f, b" DUMMY   ");
                } else if shape == XsectShape::Irregular {
                    cfprintf!(f, b" CHANNEL ");
                } else {
                    cfprintf!(f, b" %-7s ", lt_str(lt));
                }
                cfprintf!(f, ff, mf);
                cfprintf!(f, b"  %4d  %02d:%02d", days, hrs, mins);

                if lt == LinkType::Pump as i32 {
                    let mut q_full = 0.0;
                    let pr = ctx.link_subtypes.pump_row(j);
                    if pr >= 0 {
                        let ci = ctx.link_subtypes.pumps.curve[pr as usize];
                        if ci >= 0 && ci < ctx.tables.tables.len() as i32 {
                            let ty = &ctx.tables.tables[ci as usize].y;
                            if !ty.is_empty() {
                                // *std::max_element: the first largest by operator<
                                let mut largest = ty[0];
                                for &v in &ty[1..] {
                                    if largest < v {
                                        largest = v;
                                    }
                                }
                                q_full = largest / qcf;
                            }
                        }
                    }
                    if q_full > 0.0 {
                        cfprintf!(f, b"          ");
                        cfprintf!(f, b"  %6.2f", ctx.links.stat_max_flow[uj] / q_full);
                        continue;
                    }
                }
                if dummy || pump {
                    continue;
                }
                if lt == LinkType::Conduit as i32 {
                    if mv > 50.0 {
                        cfprintf!(f, b"    >50.00");
                    } else {
                        cfprintf!(f, b"   %7.2f", mv);
                    }
                    cfprintf!(f, b"  %6.2f", flow_ratio);
                } else {
                    cfprintf!(f, b"                  ");
                }
                let mut no_full_depth = ctx.links.xsect_y_full[uj] <= 0.0;
                if lt == LinkType::Orifice as i32 {
                    let orow = ctx.link_subtypes.orifice_row(j);
                    if orow >= 0 && ctx.link_subtypes.orifices.orifice_type[orow as usize] as f64 == 0.0 {
                        no_full_depth = true;
                    }
                }
                if no_full_depth {
                    cfprintf!(f, b"        ");
                } else {
                    cfprintf!(f, b"  %6.2f", fill);
                }
            }
        }

        wr(f, b"");
        wr(f, b"");

        // ---- Flow Classification Summary (legacy writeFlowClass) ----
        if route_rpt && (opt.routing_model == RoutingModel::Dynwave || opt.routing_model == RoutingModel::Fv) {
            wr(f, b"***************************");
            wr(f, b"Flow Classification Summary");
            wr(f, b"***************************");
            cfprintf!(f, b"\n\n  -------------------------------------------------------------------------------------\n                      Adjusted    ---------- Fraction of Time in Flow Class ---------- \n                       /Actual         Up    Down  Sub   Sup   Up    Down  Norm  Inlet \n  Conduit               Length    Dry  Dry   Dry   Crit  Crit  Crit  Crit  Ltd   Ctrl  \n  -------------------------------------------------------------------------------------");
            let total_secs = ctx.routing_stats.report_time;
            for j in 0..ctx.n_links() {
                let uj = j as usize;
                if ctx.links.type_[uj] != LinkType::Conduit {
                    continue;
                }
                if ctx.links.xsect_shape[uj] == XsectShape::Dummy {
                    continue;
                }
                let mut len_ratio = 1.0;
                let cr = ctx.link_subtypes.conduit_row(j);
                let cd = &ctx.link_subtypes.conduits;
                let l = if cr >= 0 { cd.length[cr as usize] } else { 0.0 };
                if l > 0.0 {
                    len_ratio = (if cr >= 0 { cd.mod_length[cr as usize] } else { 0.0 }) / l;
                }
                cfprintf!(f, b"\n  %-20s", ctx.link_names.name_of(j));
                cfprintf!(f, b"  %6.2f ", len_ratio);
                for c in 0..LinkData::N_FLOW_CLASSES as usize {
                    let idx = uj * LinkData::N_FLOW_CLASSES as usize + c;
                    let t = if idx < ctx.links.stat_flow_class.len() { ctx.links.stat_flow_class[idx] } else { 0.0 };
                    cfprintf!(f, b"  %4.2f", t / total_secs);
                }
                cfprintf!(f, b"  %4.2f", ctx.links.stat_norm_ltd[uj] / total_secs);
                cfprintf!(f, b"  %4.2f", ctx.links.stat_inlet_ctrl[uj] / total_secs);
            }
        }
        wr(f, b"");

        // ---- Conduit Surcharge Summary (legacy writeLinkSurcharge) ----
        if route_rpt {
            wr(f, b"*************************");
            wr(f, b"Conduit Surcharge Summary");
            wr(f, b"*************************");
            let mut n_written = 0;
            for j in 0..ctx.n_links() {
                let uj = j as usize;
                if ctx.links.type_[uj] != LinkType::Conduit || ctx.links.xsect_shape[uj] == XsectShape::Dummy {
                    continue;
                }
                let mut t: [f64; 5] = [
                    ctx.links.stat_time_surcharged[uj] / 3600.0,
                    ctx.links.stat_time_full_upstream[uj] / 3600.0,
                    ctx.links.stat_time_full_dnstream[uj] / 3600.0,
                    ctx.links.stat_time_full_both[uj] / 3600.0,
                    0.0,
                ];
                if t[0] + t[1] + t[2] + t[3] == 0.0 {
                    continue;
                }
                t[4] = ctx.links.stat_time_capacity_limited[uj] / 3600.0;
                for x in t.iter_mut() {
                    *x = cmax(0.01, *x);
                }
                if n_written == 0 {
                    wr(f, b"");
                    cfprintf!(f, b"\n  ----------------------------------------------------------------------------\n                                                           Hours        Hours \n                         --------- Hours Full --------   Above Full   Capacity\n  Conduit                Both Ends  Upstream  Dnstream   Normal Flow   Limited\n  ----------------------------------------------------------------------------");
                    n_written = 1;
                }
                cfprintf!(f, b"\n  %-20s", ctx.link_names.name_of(j));
                cfprintf!(f, b"    %8.2f  %8.2f  %8.2f  %8.2f     %8.2f", t[0], t[1], t[2], t[3], t[4]);
            }
            if n_written == 0 {
                wr(f, b"");
                wr(f, b"No conduits were surcharged.");
            }
        }

        wr(f, b"");
        wr(f, b"");

        // ---- Pumping Summary (legacy writePumpFlows) ----
        if route_rpt {
            let mut any_pump = false;
            for j in 0..ctx.n_links() {
                if ctx.links.type_[j as usize] == LinkType::Pump {
                    any_pump = true;
                    break;
                }
            }
            if any_pump {
                wr(f, b"***************");
                wr(f, b"Pumping Summary");
                wr(f, b"***************");
                cfprintf!(
                    f,
                    b"\n\n  ---------------------------------------------------------------------------------------------------------\n                                                  Min       Avg       Max     Total     Power    %% Time Off\n                        Percent   Number of      Flow      Flow      Flow    Volume     Usage    Pump Curve\n  Pump                 Utilized   Start-Ups       %3s       %3s       %3s  %8s     Kw-hr    Low   High\n  ---------------------------------------------------------------------------------------------------------",
                    fuw,
                    fuw,
                    fuw,
                    vol_word
                );
                let total_seconds = ctx.routing_stats.report_time;
                for j in 0..ctx.n_links() {
                    let uj = j as usize;
                    if ctx.links.type_[uj] != LinkType::Pump {
                        continue;
                    }
                    let max_flow = ctx.links.stat_max_flow[uj] * qcf;
                    let vol = ctx.links.stat_pump_volume[uj] * vcf;
                    let on_time = ctx.links.stat_pump_on_time[uj];
                    let pct_utilized = if total_seconds > 0.0 { on_time / total_seconds * 100.0 } else { 0.0 };
                    let start_ups = ctx.links.stat_pump_cycles[uj] as i32;
                    let avg_flow = if ctx.links.stat_pump_periods[uj] > 0 { ctx.links.stat_pump_sum_flow[uj] / ctx.links.stat_pump_periods[uj] as f64 * qcf } else { 0.0 };
                    let energy_kwh = ctx.links.stat_pump_energy[uj];
                    cfprintf!(f, b"\n  %-20s %8.2f  %10d %9.2f %9.2f %9.2f %9.3f %9.2f", ctx.link_names.name_of(j), pct_utilized, start_ups, 0.0, avg_flow, max_flow, vol, energy_kwh);
                    cfprintf!(
                        f,
                        b" %6.1f %6.1f",
                        if on_time > 0.0 { ctx.links.stat_pump_off_curve_low[uj] / on_time * 100.0 } else { 0.0 },
                        if on_time > 0.0 { ctx.links.stat_pump_off_curve_high[uj] / on_time * 100.0 } else { 0.0 }
                    );
                }
                wr(f, b"");
                wr(f, b"");
            }
        }

        // ---- Street tables ----
        let street_of_link = |j: i32| -> i32 {
            let uj = j as usize;
            if ctx.links.xsect_shape[uj] != XsectShape::StreetXsect {
                return -1;
            }
            let tt = ctx.links.xsect_curve[uj];
            if tt < 0 || tt >= ctx.transect_tables.len() as i32 {
                return -1;
            }
            let sname = &ctx.transect_tables[tt as usize].name;
            for s in 0..ctx.streets.count() {
                let a = &ctx.streets.names[s as usize];
                if a.len() != sname.len() {
                    continue;
                }
                if a.iter().zip(sname.iter()).all(|(x, y)| x.to_ascii_lowercase() == y.to_ascii_lowercase()) {
                    return s;
                }
            }
            -1
        };

        // Street Flow Summary (legacy writeStreetStats)
        if route_rpt && ctx.streets.count() > 0 {
            let mut header = false;
            let inv_len = 1.0 / len_ucf;
            for j in 0..ctx.n_links() {
                let si = street_of_link(j);
                if si < 0 {
                    continue;
                }
                let uj = j as usize;
                let su = si as usize;
                if !header {
                    wr(f, b"*******************");
                    wr(f, b"Street Flow Summary");
                    wr(f, b"*******************");
                    cfprintf!(
                        f,
                        b"\n\n  -------------------------------------------------------\n                          Peak   Maximum   Maximum\n                          Flow    Spread     Depth\n  Street Conduit           %-3s     %-5s     %-5s\n  -------------------------------------------------------",
                        fuw,
                        if si_report { b"m" as &[u8] } else { b"ft" },
                        if si_report { b"m" as &[u8] } else { b"ft" }
                    );
                    header = true;
                }
                let max_flow = ctx.links.stat_max_flow[uj];
                let max_depth = ctx.links.stat_max_filling[uj] * ctx.links.xsect_y_full[uj];
                let xs = crate::hydraulics::link::build_xsect_params(&ctx.links, uj, Some(&ctx.transect_tables));
                let sides = cmax(ctx.streets.sides[su], 1);
                let mut max_spread = crate::hydraulics::xsection::get_w_of_y(&xs, max_depth) / sides as f64;
                max_spread = cmin(max_spread, ctx.streets.t_crown[su] * inv_len);
                cfprintf!(f, b"\n  %-16s %9.3f %9.3f %9.3f", ctx.link_names.name_of(j), max_flow * qcf, max_spread * len_ucf, max_depth * len_ucf);
            }
            if header {
                wr(f, b"");
                wr(f, b"");
            }
        }

        // Street Inlet Flow Summary
        if route_rpt && ctx.inlet_usages.count() > 0 {
            let ni = ctx.inlet_usages.count();
            let has_stats = ctx.inlet_usages.stat_capture_vol.len() as i32 >= ni;
            let has_diag = ctx.inlet_diag.count() >= ni;
            wr(f, b"**************************");
            wr(f, b"Street Inlet Flow Summary");
            wr(f, b"**************************");
            if !has_stats {
                wr(f, b"");
                wr(f, b"No inlet statistics available.");
            } else {
                const FT3_TO_KGAL: f64 = 7.48052 / 1000.0;
                cfprintf!(
                    f,
                    b"\n\n  ------------------------------------------------------------------------------------------------------------------------------------------\n                                                              Peak      Peak       Avg.    Bypass      Back      Peak      Peak      Vol.      Vol.\n                                                              Flow   Capture    Capture      Flow      Flow   Capture    Bypass  Captured  Bypassed\n  Inlet Location        Inlet Design       Placement  Count     %-3s      Pcnt       Pcnt      Pcnt      Pcnt   / Inlet      %-3s  1000 Gal  1000 Gal\n  ------------------------------------------------------------------------------------------------------------------------------------------",
                    fuw,
                    fuw
                );
                for i in 0..ni {
                    let ui = i as usize;
                    let li = ctx.inlet_usages.link_index[ui];
                    let nh = ctx.inlet_usages.node_host[ui];
                    let di = ctx.inlet_usages.design_index[ui];
                    let mut host_name = b"?".to_vec();
                    if nh >= 0 {
                        host_name = ctx.node_names.name_of(nh).to_vec();
                        host_name.extend_from_slice(b" (node)");
                    } else if li >= 0 {
                        host_name = ctx.link_names.name_of(li).to_vec();
                    }
                    let inlet_name: &[u8] = if di >= 0 && di < ctx.inlets.count() { &ctx.inlets.names[di as usize] } else { b"?" };
                    let cap_vol = ctx.inlet_usages.stat_capture_vol[ui];
                    let byp_vol = ctx.inlet_usages.stat_bypass_vol[ui];
                    let peak = (if has_diag { ctx.inlet_diag.peak_flow[ui] } else { ctx.inlet_usages.stat_peak_flow[ui] }) * qcf;
                    let cap_kgal = cap_vol * FT3_TO_KGAL;
                    let byp_kgal = byp_vol * FT3_TO_KGAL;
                    let mut placement: &[u8] = b"ON-GRADE";
                    let mut n_inlets = 1;
                    let (mut pfc, mut afc, mut bpf, mut bff) = (0.0, 0.0, 0.0, 0.0);
                    let (mut peak_cap_per_inlet, mut peak_bypass) = (0.0, 0.0);
                    if has_diag {
                        let dg = &ctx.inlet_diag;
                        placement = if dg.is_sag[ui] != 0 { b"ON-SAG  " } else { b"ON-GRADE" };
                        n_inlets = cmax(dg.num_inlets[ui], 1);
                        let fp = dg.flow_periods[ui] as f64 / 100.0;
                        if fp > 0.0 {
                            let cp = dg.capture_periods[ui] as f64 / 100.0;
                            pfc = dg.peak_flow_capture[ui];
                            if cp > 0.0 {
                                afc = dg.avg_flow_capture[ui] / cp;
                                bpf = dg.bypass_freq[ui] / cp;
                            }
                            bff = dg.backflow_periods[ui] as f64 / fp;
                            let hl = if li >= 0 { li } else { dg.up_link[ui] };
                            if hl >= 0 && hl < ctx.n_links() {
                                let uh = hl as usize;
                                let max_flow = ctx.links.stat_max_flow[uh];
                                let si = street_of_link(hl);
                                let sides = if si >= 0 { cmax(ctx.streets.sides[si as usize], 1) } else { 1 };
                                peak_cap_per_inlet = (max_flow / sides as f64) * qcf * 0.01 * pfc / n_inlets as f64;
                                peak_bypass = max_flow * qcf * 0.01 * (100.0 - pfc);
                            }
                        }
                    }
                    cfprintf!(
                        f,
                        b"\n  %-20s  %-16s   %-9s  %5d  %8.3f  %8.2f  %9.2f  %8.2f  %8.2f  %8.3f  %8.3f  %8.3f  %8.3f",
                        host_name,
                        inlet_name,
                        placement,
                        n_inlets,
                        peak,
                        pfc,
                        afc,
                        bpf,
                        bff,
                        peak_cap_per_inlet,
                        peak_bypass,
                        cap_kgal,
                        byp_kgal
                    );
                }
            }
            wr(f, b"");
            wr(f, b"");
        }

        // ---- Link Pollutant Load Summary (legacy writeLinkLoads) ----
        if route_rpt && ctx.n_pollutants() > 0 && !opt.ignore_quality {
            let nl = ctx.n_links();
            let np = ctx.n_pollutants();
            const LT_PER_FT3: f64 = 28.317;
            let mut mass_cf = vec![0.0f64; np as usize];
            let mut unit_str: Vec<&[u8]> = vec![b""; np as usize];
            for p in 0..np {
                let up = p as usize;
                let mu = if up < ctx.pollutants.units.len() { ctx.pollutants.units[up] } else { MassUnits::MgPerL };
                unit_str[up] = if mu == MassUnits::CountsPerL { b"LogN" } else { load_word };
                mass_cf[up] = LT_PER_FT3
                    * (if mu == MassUnits::MgPerL {
                        cont_mass
                    } else if mu == MassUnits::UgPerL {
                        cont_mass / 1000.0
                    } else {
                        1.0
                    });
            }
            wr(f, b"***************************");
            wr(f, b"Link Pollutant Load Summary");
            wr(f, b"***************************");
            cfprintf!(f, b"\n");
            cfprintf!(f, b" \n  ----------------------------------");
            for _p in 1..np {
                cfprintf!(f, b"--------------");
            }
            cfprintf!(f, b" \n                                ");
            for p in 0..np {
                cfprintf!(f, b"%14s", ctx.pollutant_names.name_of(p));
            }
            cfprintf!(f, b" \n  Link                          ");
            for p in 0..np {
                cfprintf!(f, b"%14s", unit_str[p as usize]);
            }
            cfprintf!(f, b" \n  ----------------------------------");
            for _p in 1..np {
                cfprintf!(f, b"--------------");
            }
            for j in 0..nl {
                let uj = j as usize;
                cfprintf!(f, b"\n  %-30s", ctx.link_names.name_of(j));
                for p in 0..np {
                    let up = p as usize;
                    let idx = uj * np as usize + up;
                    let mut x = (if idx < ctx.links.stat_total_load.len() { ctx.links.stat_total_load[idx] } else { 0.0 }) * mass_cf[up];
                    if is_count(p) && x > 0.0 {
                        x = x.log10(); // LIBM
                    }
                    cfprintf!(f, if x < 10000.0 { b"%14.3f" as &[u8] } else { b"%14.3e" }, x);
                }
            }
            wr(f, b"");
            wr(f, b"");
        }
    }
}
