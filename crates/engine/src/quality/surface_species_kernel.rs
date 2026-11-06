// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/quality/SurfaceSpeciesKernel.hpp and
// SurfaceSpeciesKernel.cpp, and the block builders pollutantSurfaceMcf /
// pollutantSurfaceBlock in src/engine/core/SWMMEngine.cpp (HydroCouple/openswmm.engine,
// branch swmm6_rel, commit 88329d6), Copyright (c) 2026 Caleb Buahin, licensed under the
// Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 5 November 2026, by Claude (Anthropic) at the
// request of Robert Dickinson. Not an OpenSWMM release.
//
// Buildup, sweeping, washoff and the ponded complete-mix balance of one species family
// (surfacespecies::step), and the post-LID mass bookings (surfacespecies::bookPostLid).
//
// Differences of form, never in result:
//   * Upstream's Block is a set of pointers into the family's arrays. Only the pollutant
//     family is translated (stepSurfaceQuality / bookWashoffLoads build it with
//     pollutantSurfaceBlock): its arrays are named directly here (ctx.subcatches.conc,
//     ponded_qual, runon_qual_rate, total_load; ctx.pollutants.c_rain, snow_only,
//     co_pollut, co_frac; the qual_* ledgers; the sweep triple), the buildup store and the
//     land-use parameters are passed in. The pollutant block sets no reacted ponded mass,
//     ponded-reactor mask or LID-reactor mask, so those branches take their null paths.
//     The reactions species family (msxsurf::step / bookPostLid) is NOT YET TRANSLATED.

use crate::core::cmath::{cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::datetime;
use crate::core::unit_conversion as ucf;
use crate::data::pollutant::MassUnits;
use crate::data::table::table_lookup_cursor;
use crate::quality::landuse::{BuildupType, LanduseSolver, WashoffType};

/// pollutantSurfaceMcf: per-pollutant concentration-mass -> user-mass multipliers.
pub fn pollutant_surface_mcf(ctx: &SimulationContext) -> Vec<f64> {
    let np = ctx.n_pollutants();
    let mut mcf = vec![0.0f64; if np > 0 { np as usize } else { 0 }];
    let mass_ucf = ucf::ucf(ucf::MASS, &ctx.options);
    for p in 0..np.max(0) as usize {
        let mut m = mass_ucf;
        if p < ctx.pollutants.units.len() {
            match ctx.pollutants.units[p] {
                MassUnits::UgPerL => m = mass_ucf / 1000.0,
                MassUnits::CountsPerL => m = 1.0,
                _ => {}
            }
        }
        mcf[p] = m;
    }
    mcf
}

#[inline]
fn bu_idx(nlu: i32, n: i32, sc: i32, lu: i32, s: i32) -> usize {
    sc as usize * (nlu * n) as usize + lu as usize * n as usize + s as usize
}

/// surfacespecies::step for the pollutant family.
pub fn step_pollutants(ctx: &mut SimulationContext, ls: &LanduseSolver, buildup: &mut [f64], mcf: &[f64], dt_runoff: f64, old_runoff_time: f64) {
    let n = ctx.n_pollutants();
    if n <= 0 {
        return;
    }
    let np = n;
    let nlu = ctx.n_landuses();
    const MIN_RUNOFF_FTS: f64 = 2.31481e-8;
    let dt_days = dt_runoff / ucf::SEC_PER_DAY;
    const K_L_PER_FT3: f64 = 28.317;

    let ucf_landarea_sq = ucf::ucf(ucf::LANDAREA, &ctx.options);
    for i in 0..ctx.n_subcatches() {
        let ui = i as usize;
        let q = ctx.subcatches.runoff[ui];
        let area_ac = ctx.subcatches.area[ui];
        let area_ft2_sc = area_ac / ucf_landarea_sq;
        let runoff_fts = if ui < ctx.subcatches.subarea_runoff_rate.len() { ctx.subcatches.subarea_runoff_rate[ui] } else { 0.0 };

        let mut q_out1 = q;
        {
            let lid_ft2_sc = if ui < ctx.subcatches.total_lid_area_ft2.len() { ctx.subcatches.total_lid_area_ft2[ui] } else { 0.0 };
            if lid_ft2_sc > 0.0 {
                q_out1 += ctx.subcatches.rainfall[ui] * lid_ft2_sc;
                if area_ft2_sc == lid_ft2_sc {
                    q_out1 += ctx.subcatches.runon_rate[ui] * area_ft2_sc;
                }
            }
        }
        if ui < ctx.subcatches.washoff_vout1.len() {
            ctx.subcatches.washoff_vout1[ui] = q_out1 * dt_runoff;
        }

        for p in 0..np {
            let up = p as usize;
            let sq_idx = ui * np as usize + up;
            let mcf_p = mcf[up];
            let mut total_washoff_load = 0.0;
            let mut avg_bmp_effic = 0.0;

            for lu in 0..nlu {
                let cov_idx = ui * nlu as usize + lu as usize;
                let frac = if cov_idx < ctx.subcatches.coverage.len() { ctx.subcatches.coverage[cov_idx] / 100.0 } else { 0.0 };
                avg_bmp_effic += frac * ls.washoff_params[(lu * np + p) as usize].bmp_effic / 100.0;
                if frac <= 0.0 {
                    continue;
                }

                let k = (lu * np + p) as usize;
                let bp = ls.buildup_params[k];
                let wp = ls.washoff_params[k];
                let bu = bu_idx(nlu, n, i, lu, p);

                let norm = if bp.normalizer == 0 { frac * area_ac } else { frac * ctx.subcatches.curb_length[ui] };

                // street sweeping, after buildup, before washoff
                let apply_sweep = |ctx: &mut SimulationContext, buildup: &mut [f64]| {
                    let sk = ui * nlu as usize + lu as usize;
                    if sk >= ctx.subcatches.sweep_now.len() || ctx.subcatches.sweep_now[sk] == 0 {
                        return;
                    }
                    let removal_frac = ctx.landuses.sweep_removal[lu as usize];
                    let effic = wp.sweep_effic / 100.0;
                    let old_bu = buildup[bu];
                    let mut new_bu = old_bu * (1.0 - removal_frac * effic);
                    new_bu = cmin(old_bu, new_bu);
                    new_bu = cmax(0.0, new_bu);
                    buildup[bu] = new_bu;
                    ctx.mass_balance.qual_sweeping[up] += (old_bu - new_bu) * norm;
                };

                let subcatch_is_dry = runoff_fts < MIN_RUNOFF_FTS;
                let snow_only_blocked = up < ctx.pollutants.snow_only.len() && ctx.pollutants.snow_only[up] && (ui >= ctx.subcatches.snow_depth.len() || ctx.subcatches.snow_depth[ui] < 0.001 / 12.0);
                if bp.type_ != BuildupType::None && subcatch_is_dry && !snow_only_blocked {
                    let mass = buildup[bu];

                    if bp.type_ == BuildupType::External {
                        let max_bu = bp.coeff[0];
                        let sf = bp.coeff[1];
                        let ts_idx = bp.coeff[2] as i32;
                        let mut rate = 0.0;
                        if ts_idx >= 0 && ts_idx < ctx.tables.tables.len() as i32 {
                            let abs_t = datetime::add_seconds(ctx.options.start_date, old_runoff_time);
                            rate = sf * table_lookup_cursor(&mut ctx.tables.tables[ts_idx as usize], abs_t);
                        }
                        let mut new_mass = cmin(mass + rate * dt_days, max_bu);
                        new_mass = cmax(new_mass, 0.0);
                        let buildup_change = new_mass - mass;
                        if buildup_change != 0.0 {
                            buildup[bu] = new_mass;
                            ctx.mass_balance.qual_surface_buildup[up] += buildup_change * norm;
                        }
                        apply_sweep(ctx, buildup);
                        continue;
                    }

                    let mut days = 0.0;
                    let (c0, c1, c2) = (bp.coeff[0], bp.coeff[1], bp.coeff[2]);
                    if mass > 0.0 {
                        match bp.type_ {
                            BuildupType::Power => {
                                days = if c1 * c2 > 0.0 { (mass / c1).powf(1.0 / c2) } else { 0.0 }; // LIBM
                            }
                            BuildupType::Expon => {
                                days = if c0 * c1 > 0.0 && mass < c0 { -(1.0 - mass / c0).ln() / c1 } else { 0.0 }; // LIBM
                            }
                            BuildupType::Satur => {
                                days = if c0 > mass { mass * c2 / (c0 - mass) } else { bp.max_days };
                            }
                            _ => {}
                        }
                    }
                    days += dt_days;
                    let mut new_mass = 0.0;
                    if days > 0.0 {
                        match bp.type_ {
                            BuildupType::Power => {
                                new_mass = cmin(c1 * days.powf(c2), c0); // LIBM
                            }
                            BuildupType::Expon => {
                                new_mass = c0 * (1.0 - (-c1 * days).exp()); // LIBM
                            }
                            BuildupType::Satur => {
                                new_mass = if c2 + days > 0.0 { c0 * days / (c2 + days) } else { 0.0 };
                            }
                            _ => {}
                        }
                    }
                    new_mass = cmax(new_mass, mass);
                    let buildup_change = new_mass - mass;
                    buildup[bu] = new_mass;
                    if buildup_change > 0.0 {
                        ctx.mass_balance.qual_surface_buildup[up] += buildup_change * norm;
                    }
                }
                apply_sweep(ctx, buildup);

                if wp.type_ != WashoffType::None && runoff_fts >= MIN_RUNOFF_FTS && q > 0.0 {
                    let bu_v = buildup[bu];
                    let area_ft2 = area_ac / ucf::ucf(ucf::LANDAREA, &ctx.options);
                    let q_expon = runoff_fts * ucf::ucf(ucf::RAINFALL, &ctx.options);
                    let landuse_ft2 = frac * area_ft2;

                    let mut load = 0.0;
                    match wp.type_ {
                        WashoffType::Emc => {
                            load = wp.coeff * K_L_PER_FT3 * q * frac;
                        }
                        WashoffType::Expon => {
                            if bu_v > 0.0 && mcf_p > 0.0 && landuse_ft2 > 0.0 {
                                let mut wq = (wp.coeff / 3600.0) * q_expon.powf(wp.expon) * (bu_v * norm) / mcf_p; // LIBM
                                wq /= runoff_fts * landuse_ft2;
                                load = ((wq * q) * landuse_ft2) / area_ft2;
                            }
                        }
                        WashoffType::Rating => {
                            if area_ft2 > 0.0 {
                                let wq = wp.coeff * ucf::ucf(ucf::FLOW, &ctx.options).powf(wp.expon) * (runoff_fts * landuse_ft2).powf(wp.expon - 1.0); // LIBM
                                load = ((wq * q) * landuse_ft2) / area_ft2;
                            }
                        }
                        _ => {}
                    }

                    let avail = if mcf_p > 0.0 { (bu_v * norm) / mcf_p } else { 0.0 };
                    let max_load = if dt_runoff > 0.0 { avail / dt_runoff } else { 0.0 };
                    let mut capped = false;
                    if load > max_load && bp.type_ != BuildupType::None {
                        load = max_load;
                        capped = true;
                    }

                    let washed_mass = load * mcf_p * dt_runoff;
                    let has_buildup_fn = bp.type_ != BuildupType::None;
                    if !has_buildup_fn && washed_mass > bu_v * norm {
                        if up < ctx.mass_balance.qual_surface_buildup.len() {
                            ctx.mass_balance.qual_surface_buildup[up] += washed_mass;
                        }
                        buildup[bu] = 0.0;
                    } else if capped {
                        buildup[bu] = 0.0;
                    } else if norm > 0.0 {
                        let washed = washed_mass / norm;
                        buildup[bu] = cmax(buildup[bu] - washed, 0.0);
                    }

                    let bmp_removed = load * (wp.bmp_effic / 100.0);
                    if bmp_removed > 0.0 {
                        if up < ctx.mass_balance.qual_bmp_removal.len() {
                            ctx.mass_balance.qual_bmp_removal[up] += bmp_removed * dt_runoff * mcf_p;
                        }
                        load -= bmp_removed;
                    }

                    total_washoff_load += load;
                }
            } // end land use loop

            // Ponded quality (legacy findPondedLoads)
            let lid_ft2_pq = if ui < ctx.subcatches.total_lid_area_ft2.len() { ctx.subcatches.total_lid_area_ft2[ui] } else { 0.0 };
            let area_ft2_full_pq = area_ac / ucf::ucf(ucf::LANDAREA, &ctx.options);
            let non_lid_ft2_pq = area_ft2_full_pq - lid_ft2_pq;
            if sq_idx < ctx.subcatches.ponded_qual.len() && non_lid_ft2_pq > 0.0 {
                const L_PER_FT3: f64 = 28.317;
                let area_ft2_pq = non_lid_ft2_pq;
                let rrate = if ui < ctx.subcatches.rainfall.len() { ctx.subcatches.rainfall[ui] } else { 0.0 };

                let v_rain_pq = rrate * area_ft2_pq * dt_runoff;
                let v_runon_pq = ctx.subcatches.runon_inflow[ui] * dt_runoff;
                let v_infil_pq = ctx.subcatches.infil_loss[ui] * area_ft2_pq * dt_runoff;
                let v_outflow_pq = q * dt_runoff;

                let v_ponded_pq = ctx.subcatches.old_ponded_depth[ui] * area_ft2_pq;
                let v_inflow_pq = v_rain_pq + v_runon_pq + v_ponded_pq;

                let c_rain_pq = if up < ctx.pollutants.c_rain.len() { ctx.pollutants.c_rain[up] } else { 0.0 };
                let w_rain_pq = c_rain_pq * L_PER_FT3 * v_rain_pq;

                if w_rain_pq > 0.0 && up < ctx.mass_balance.qual_wet_deposition.len() {
                    ctx.mass_balance.qual_wet_deposition[up] += w_rain_pq * mcf_p;
                }

                let w_mass_pq = ctx.subcatches.ponded_qual[sq_idx];

                if v_inflow_pq <= 0.0 {
                    if ctx.subcatches.ponded_qual[sq_idx] != 0.0 && up < ctx.mass_balance.qual_ponded_lost.len() {
                        ctx.mass_balance.qual_ponded_lost[up] += ctx.subcatches.ponded_qual[sq_idx] * mcf_p;
                    }
                    ctx.subcatches.ponded_qual[sq_idx] = 0.0;
                } else {
                    let w_runon_pq = if sq_idx < ctx.subcatches.runon_qual_rate.len() { ctx.subcatches.runon_qual_rate[sq_idx] * dt_runoff } else { 0.0 };
                    let mut w_total_pq = w_mass_pq + w_rain_pq + w_runon_pq;
                    let c_ponded_pq = w_total_pq / v_inflow_pq;

                    let w_infil_pq = cmin(c_ponded_pq * v_infil_pq, w_total_pq);
                    w_total_pq -= w_infil_pq;
                    if w_infil_pq > 0.0 && up < ctx.mass_balance.qual_infil_loss.len() {
                        ctx.mass_balance.qual_infil_loss[up] += w_infil_pq * mcf_p;
                    }

                    let mut w_outflow_pq = cmin(c_ponded_pq * v_outflow_pq, w_total_pq);
                    w_total_pq -= w_outflow_pq;
                    let _ = w_total_pq;

                    let bmp_pq = avg_bmp_effic * w_outflow_pq;
                    if bmp_pq > 0.0 {
                        if up < ctx.mass_balance.qual_bmp_removal.len() {
                            ctx.mass_balance.qual_bmp_removal[up] += bmp_pq * mcf_p;
                        }
                        w_outflow_pq -= bmp_pq;
                    }

                    if w_outflow_pq > 0.0 && dt_runoff > 0.0 {
                        total_washoff_load += w_outflow_pq / dt_runoff;
                    }

                    ctx.subcatches.ponded_qual[sq_idx] = cmax(c_ponded_pq * ctx.subcatches.ponded_depth[ui] * area_ft2_pq, 0.0);
                }
            }

            // LID area wet deposition and runon (legacy findLidLoads)
            if ui < ctx.subcatches.total_lid_area_ft2.len() {
                let lid_ft2 = ctx.subcatches.total_lid_area_ft2[ui];
                if lid_ft2 > 0.0 && dt_runoff > 0.0 {
                    const L_PER_FT3_37: f64 = 28.317;
                    let rrate37 = if ui < ctx.subcatches.rainfall.len() { ctx.subcatches.rainfall[ui] } else { 0.0 };
                    let c_rain37 = if up < ctx.pollutants.c_rain.len() { ctx.pollutants.c_rain[up] } else { 0.0 };
                    let v_lid_rain = rrate37 * lid_ft2 * dt_runoff;
                    let w_lid_rain = c_rain37 * L_PER_FT3_37 * v_lid_rain;
                    if w_lid_rain > 0.0 {
                        total_washoff_load += w_lid_rain / dt_runoff;
                        if up < ctx.mass_balance.qual_wet_deposition.len() {
                            ctx.mass_balance.qual_wet_deposition[up] += w_lid_rain * mcf_p;
                        }
                    }
                    let full_ft2 = ctx.subcatches.area[ui] / ucf::ucf(ucf::LANDAREA, &ctx.options);
                    if (lid_ft2 - full_ft2).abs() < 1.0 {
                        let w_lid_runon = if sq_idx < ctx.subcatches.runon_qual_rate.len() { ctx.subcatches.runon_qual_rate[sq_idx] * dt_runoff } else { 0.0 };
                        if w_lid_runon > 0.0 {
                            total_washoff_load += w_lid_runon / dt_runoff;
                        }
                    }
                }
            }

            let mut conc = 0.0;
            if q_out1 > 0.0 && total_washoff_load > 0.0 {
                conc = total_washoff_load / q_out1 / K_L_PER_FT3;
            }
            ctx.subcatches.conc[sq_idx] = conc;
        } // end pollutant loop

        // Co-pollutant washoff (legacy findWashoffLoads, surfqual.c:505-521)
        if q_out1 > 0.0 {
            for p in 0..np {
                let up = p as usize;
                if up >= ctx.pollutants.co_pollut.len() {
                    continue;
                }
                let k = ctx.pollutants.co_pollut[up];
                if k < 0 || k >= np {
                    continue;
                }
                let frac = if up < ctx.pollutants.co_frac.len() { ctx.pollutants.co_frac[up] } else { 0.0 };
                if frac <= 0.0 {
                    continue;
                }
                let ip = ui * np as usize + up;
                let ik = ui * np as usize + k as usize;
                if ip >= ctx.subcatches.conc.len() || ik >= ctx.subcatches.conc.len() {
                    continue;
                }
                let add = frac * ctx.subcatches.conc[ik];
                if add <= 0.0 {
                    continue;
                }
                ctx.subcatches.conc[ip] += add;
                let mcf_cp = mcf[up];
                let w = add * q_out1 * 28.317;
                if up < ctx.mass_balance.qual_surface_buildup.len() {
                    ctx.mass_balance.qual_surface_buildup[up] += w * dt_runoff * mcf_cp;
                }
            }
        }
    }
}

/// surfacespecies::bookPostLid for the pollutant family.
pub fn book_post_lid_pollutants(ctx: &mut SimulationContext, mcf: &[f64], dt_runoff: f64, lid_infil_vol: &[f64]) {
    let np = ctx.n_pollutants();
    let ns = ctx.n_subcatches();
    if np <= 0 || ns <= 0 || ctx.options.ignore_quality {
        return;
    }
    const K_L_PER_FT3: f64 = 28.317;
    for i in 0..ns {
        let ui = i as usize;
        if ctx.subcatches.area[ui] <= 0.0 {
            continue;
        }
        let v_out1 = if ui < ctx.subcatches.washoff_vout1.len() { ctx.subcatches.washoff_vout1[ui] } else { 0.0 };
        let v_surf_out = ctx.subcatches.runoff[ui] * dt_runoff;
        let v_drain = if ui < ctx.subcatches.lid_drain_flow.len() { ctx.subcatches.lid_drain_flow[ui] * dt_runoff } else { 0.0 };
        let v_out2 = v_surf_out + v_drain;

        let lid_ft2 = if ui < ctx.subcatches.total_lid_area_ft2.len() { ctx.subcatches.total_lid_area_ft2[ui] } else { 0.0 };
        let has_lid = lid_ft2 > 0.0;
        let v_lost = if has_lid { cmax(v_out1 - v_out2, 0.0) } else { 0.0 };
        let v_infil = cmin(lid_infil_vol[ui], v_lost);

        let load_reaches_system = ctx.subcatches.outlet_node[ui] >= 0 || ctx.subcatches.outlet_subcatch[ui] == i;

        for p in 0..np as usize {
            let sq_idx = ui * np as usize + p;
            if sq_idx >= ctx.subcatches.conc.len() {
                continue;
            }
            let c_out = ctx.subcatches.conc[sq_idx] * K_L_PER_FT3;
            if c_out <= 0.0 {
                continue;
            }
            let mcf_p = mcf[p];
            if has_lid {
                let m_infil = c_out * v_infil * mcf_p;
                if m_infil > 0.0 && p < ctx.mass_balance.qual_infil_loss.len() {
                    ctx.mass_balance.qual_infil_loss[p] += m_infil;
                }
                let m_bmp = c_out * (v_lost - v_infil) * mcf_p;
                if m_bmp > 0.0 && p < ctx.mass_balance.qual_bmp_removal.len() {
                    ctx.mass_balance.qual_bmp_removal[p] += m_bmp;
                }
            }
            ctx.subcatches.total_load[sq_idx] += c_out * v_out2 * mcf_p;
            if load_reaches_system && p < ctx.mass_balance.qual_runoff_load.len() {
                ctx.mass_balance.qual_runoff_load[p] += c_out * v_surf_out * mcf_p;
            }
        }
    }
}
