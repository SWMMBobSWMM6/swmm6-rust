// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/components/HeatFluxModules/BedExchange.cpp and
// BedExchange.hpp (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6), Copyright 2026
// Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 6 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// [SEDIMENT_EXCHANGE] (phase H6b): the channel bed as a second thermal body coupled to the water
// and to the deep ground (relaxPair, the exact 2x2 linear solve), and the bed's solute pool
// (exchangePair). fillChannelThermalGeometry serves the C API's channel regions and is not
// translated.

use crate::core::context::SimulationContext;
use crate::core::engine_step_routing::StepError;
use crate::data::bed_zone::SedimentConfig;
use crate::data::heat_override::HeatElement;
use crate::data::table::table_tseries_lookup_cursor;
use crate::hydraulics::link as link_mod;
use crate::hydraulics::xsection as xsect;
use crate::transport::heat_fluxes::{K_CU_FT_TO_CU_M, K_SQ_FT_TO_SQ_M};
use crate::transport::heat_overrides::sediment_for;

/// heat::BedCoupling
#[derive(Debug, Clone, Copy, Default)]
pub struct BedCoupling {
    pub g_wb: f64,
    pub g_bg: f64,
    pub c_w: f64,
    pub c_b: f64,
    pub t_gr: f64,
}

impl BedCoupling {
    pub fn viable(&self) -> bool {
        self.c_w > 0.0 && self.c_b > 0.0
    }
}

/// heat::PairStep
#[derive(Debug, Clone, Copy, Default)]
pub struct PairStep {
    pub dt_w: f64,
    pub dt_b: f64,
}

/// heat::SolutePairStep
#[derive(Debug, Clone, Copy, Default)]
pub struct SolutePairStep {
    pub dc_w: f64,
    pub dc_b: f64,
}

fn phi1(z: f64) -> f64 {
    if z.abs() < 1e-8 {
        return 1.0 + 0.5 * z;
    }
    z.exp_m1() / z // LIBM
}

fn phi2(z: f64) -> f64 {
    if z.abs() < 1e-5 {
        return 0.5 + z / 6.0;
    }
    (z.exp_m1() - z) / (z * z) // LIBM
}

fn build_xsp(ctx: &SimulationContext, uk: usize) -> crate::hydraulics::xsect_batch::XSectParams<'static> {
    let links = &ctx.links;
    let mut xs = crate::hydraulics::xsect_batch::XSectParams::default();
    xs.type_ = link_mod::translate_shape(links.xsect_shape[uk]);
    xs.y_full = links.xsect_y_full[uk];
    xs.a_full = links.xsect_a_full[uk];
    xs.w_max = links.xsect_w_max[uk];
    xs.r_full = links.xsect_r_full[uk];
    xs.s_full = links.xsect_s_full[uk];
    xs.s_max = links.xsect_s_max[uk];
    xs.y_bot = links.xsect_y_bot[uk];
    xs.a_bot = links.xsect_a_bot[uk];
    xs.s_bot = links.xsect_s_bot[uk];
    xs.r_bot = links.xsect_r_bot[uk];
    xs
}

/// relaxPair: the water and the bed relaxed together over dt (exact for the linearized pair).
#[allow(clippy::too_many_arguments)]
pub fn relax_pair(g: &BedCoupling, t_w: f64, t_b: f64, j0: f64, j1: f64, h: f64, area_m2: f64, dt: f64) -> PairStep {
    let mut out = PairStep::default();
    if !g.viable() || !(dt > 0.0) {
        return out;
    }
    let mut djdt = 0.0;
    if h > 0.0 {
        let s = (j1 - j0) / h;
        if s > 0.0 && s.is_finite() {
            djdt = s;
        }
    }
    if g.g_bg == 0.0 && djdt == 0.0 && j0 == 0.0 {
        let total = g.c_w + g.c_b;
        let rate = g.g_wb / g.c_w + g.g_wb / g.c_b;
        let difference = (t_w - t_b) * (-rate * dt).exp_m1(); // LIBM
        let dw = (g.c_b / total) * difference;
        let db = -(g.c_w / g.c_b) * dw;
        if !dw.is_finite() || !db.is_finite() {
            return out;
        }
        out.dt_w = dw;
        out.dt_b = db;
        return out;
    }
    let a11 = -(area_m2 * djdt + g.g_wb) / g.c_w;
    let a12 = g.g_wb / g.c_w;
    let a21 = g.g_wb / g.c_b;
    let a22 = -(g.g_wb + g.g_bg) / g.c_b;
    let c1 = (-area_m2 * j0 + g.g_wb * (t_b - t_w)) / g.c_w;
    let c2 = (g.g_wb * (t_w - t_b) + g.g_bg * (g.t_gr - t_b)) / g.c_b;
    let tr = a11 + a22;
    let diff = a11 - a22;
    let disc = diff * diff + 4.0 * a12 * a21;
    let s = (if disc > 0.0 { disc } else { 0.0 }).sqrt();
    let alpha;
    let beta;
    if s * dt > 1e-8 {
        let l1 = 0.5 * (tr + s);
        let l2 = 0.5 * (tr - s);
        let p1 = phi1(l1 * dt);
        let p2 = phi1(l2 * dt);
        alpha = (p1 - p2) / s;
        beta = (p2 * l1 - p1 * l2) / s;
    } else {
        let l = 0.5 * tr;
        alpha = dt * phi2(l * dt);
        beta = 1.0;
    }
    let mc1 = a11 * c1 + a12 * c2;
    let mc2 = a21 * c1 + a22 * c2;
    let dw = dt * (alpha * mc1 + beta * c1);
    let db = dt * (alpha * mc2 + beta * c2);
    if !dw.is_finite() || !db.is_finite() {
        return out;
    }
    out.dt_w = dw;
    out.dt_b = db;
    out
}

/// exchangePair: the exact two-pool solute exchange.
pub fn exchange_pair(c_w: f64, c_b: f64, vol_w: f64, vol_b: f64, q_exch: f64, dt: f64) -> SolutePairStep {
    let mut out = SolutePairStep::default();
    if !(vol_w > 0.0) || !(vol_b > 0.0) || !(q_exch > 0.0) || !(dt > 0.0) {
        return out;
    }
    let vt = vol_w + vol_b;
    let mu = q_exch * vt / (vol_w * vol_b);
    let d0 = c_w - c_b;
    let dd = d0 * (-mu * dt).exp_m1(); // LIBM
    let share = dd / vt;
    if !share.is_finite() {
        return out;
    }
    out.dc_w = vol_b * share;
    out.dc_b = -vol_w * share;
    out
}

/// bedExchangeEnabled
pub fn bed_exchange_enabled(ctx: &SimulationContext) -> bool {
    ctx.options.heat_transport && ctx.heat_config.sediment_exchange
}

/// groundTemperature: the GLOBAL deep-ground temperature (series or constant).
pub fn ground_temperature(ctx: &mut SimulationContext) -> Result<f64, StepError> {
    let ts = ctx.heat_config.sediment.ground_ts_index;
    if ts >= 0 && ts < ctx.tables.count() as i32 {
        let d = ctx.current_date;
        return Ok(table_tseries_lookup_cursor(&mut ctx.tables[ts], d));
    }
    Ok(ctx.heat_config.sediment.ground_temp)
}

/// linkBedAreaM2: wetted perimeter x length x barrels, m2.
pub fn link_bed_area_m2(ctx: &SimulationContext, link: i32) -> f64 {
    let uj = link as usize;
    let cr = ctx.link_subtypes.conduit_row(link);
    if cr < 0 {
        return 0.0;
    }
    let ucr = cr as usize;
    let cd = &ctx.link_subtypes.conduits;
    let mut length = cd.length[ucr];
    if !(length > 0.0) {
        length = cd.mod_length[ucr];
    }
    if !(length > 0.0) {
        return 0.0;
    }
    let depth = ctx.links.depth[uj];
    if !(depth > 0.0) {
        return 0.0;
    }
    let xs = build_xsp(ctx, uj);
    let area = xsect::get_a_of_y(&xs, depth);
    let rhyd = xsect::get_r_of_y(&xs, depth);
    if !(area > 0.0) || !(rhyd > 0.0) {
        return 0.0;
    }
    let perim_ft = area / rhyd;
    let m2 = perim_ft * length * cd.barrels[ucr] as f64 * K_SQ_FT_TO_SQ_M;
    if m2 > 0.0 {
        m2
    } else {
        0.0
    }
}

/// linkFreeSurfaceFt2: top width x length x barrels for an open conduit, ft2.
pub fn link_free_surface_ft2(ctx: &SimulationContext, link: i32) -> f64 {
    let uj = link as usize;
    let cr = ctx.link_subtypes.conduit_row(link);
    if cr < 0 {
        return 0.0;
    }
    if !xsect::is_open(ctx.links.xsect_batch_shape[uj]) {
        return 0.0;
    }
    let ucr = cr as usize;
    let cd = &ctx.link_subtypes.conduits;
    let mut length = cd.length[ucr];
    if !(length > 0.0) {
        length = cd.mod_length[ucr];
    }
    let depth = ctx.links.depth[uj];
    if !(length > 0.0) || !(depth > 0.0) {
        return 0.0;
    }
    let xs = build_xsp(ctx, uj);
    let top_width = xsect::get_w_of_y(&xs, depth);
    if !(top_width > 0.0) {
        return 0.0;
    }
    top_width * length * cd.barrels[ucr] as f64
}

/// bedExchangeQ: the solute exchange flow, m3/s.
pub fn bed_exchange_q(cfg: &SedimentConfig, bed_m2: f64) -> f64 {
    if !(bed_m2 > 0.0) || !(cfg.bed_thickness > 0.0) {
        return 0.0;
    }
    cfg.solute_diffusivity * bed_m2 / cfg.bed_thickness + cfg.hyporheic_velocity * bed_m2
}

/// groundTempFor: a series (GLOBAL) wins over the element's constant.
pub fn ground_temp_for(ctx: &SimulationContext, sc: &SedimentConfig, t_gr_global: f64) -> f64 {
    if ctx.heat_config.sediment.ground_ts_index >= 0 {
        t_gr_global
    } else {
        sc.ground_temp
    }
}

/// bedCouplingFromContact
pub fn bed_coupling_from_contact(ctx: &SimulationContext, sc: &SedimentConfig, bed_m2: f64, vol_ft3: f64, t_gr: f64) -> BedCoupling {
    let mut g = BedCoupling::default();
    if !(bed_m2 > 0.0) || !(vol_ft3 > 0.0) {
        return g;
    }
    let rho_cp_s = sc.sed_density * sc.sed_specific_heat;
    let y_bed = sc.bed_thickness;
    let y_gr = sc.ground_depth;
    if !(y_bed > 0.0) {
        return g;
    }
    let k_eff = sc.thermal_diffusivity * rho_cp_s;
    let g_cond = k_eff * bed_m2 / y_bed;
    let g_adv = ctx.options.water_density * ctx.options.water_specific_heat * sc.hyporheic_velocity * bed_m2;
    g.g_wb = g_cond + g_adv;
    g.g_bg = if y_gr > 0.0 { k_eff * bed_m2 / y_gr } else { 0.0 };
    g.c_w = ctx.options.water_density * ctx.options.water_specific_heat * vol_ft3 * K_CU_FT_TO_CU_M;
    g.c_b = rho_cp_s * bed_m2 * y_bed;
    g.t_gr = t_gr;
    g
}

/// bedCouplingForLink
pub fn bed_coupling_for_link(ctx: &SimulationContext, link: i32, vol_ft3: f64, t_gr: f64) -> BedCoupling {
    if !(vol_ft3 > 0.0) {
        return BedCoupling::default();
    }
    let sc = sediment_for(ctx, &HeatElement::link(link));
    bed_coupling_from_contact(ctx, sc, link_bed_area_m2(ctx, link), vol_ft3, ground_temp_for(ctx, sc, t_gr))
}

/// seedBedTemperature: size the per-link bed temperatures and seed them once.
pub fn seed_bed_temperature(ctx: &mut SimulationContext) -> Result<(), StepError> {
    let nl = ctx.n_links();
    let sc = &ctx.heat_config.sediment;
    let t0 = if sc.has_initial_temp { sc.initial_temp } else { sc.ground_temp };
    let bed = &mut ctx.bed_state;
    if !bed.sized(nl) {
        bed.link_temp.clear();
        bed.link_temp.resize((if nl > 0 { nl } else { 0 }) as usize, t0);
        bed.seeded = false;
    }
    if !bed.seeded {
        for t in bed.link_temp.iter_mut() {
            *t = t0;
        }
        bed.seeded = true;
    }
    Ok(())
}

/// applyBedSoluteExchange: the per-link bed solute pools (rows offset..offset+n_arr of n_total),
/// for an array laid out [link * n_arr + s]. `which` names the array: the pollutant
/// concentrations (links.conc) or the reactions species (reactions.msx_link_conc).
pub fn apply_bed_solute_exchange(ctx: &mut SimulationContext, msx: bool, n_arr: i32, offset: i32, n_total: i32, dt: f64) -> Result<(), StepError> {
    if !bed_exchange_enabled(ctx) || !(dt > 0.0) {
        return Ok(());
    }
    if n_arr <= 0 || offset < 0 || n_total < offset + n_arr {
        return Ok(());
    }
    seed_bed_temperature(ctx)?;
    let nl_all = ctx.n_links();
    let need = (if nl_all > 0 { nl_all } else { 0 }) as usize * n_total as usize;
    if ctx.bed_state.n_species != n_total || ctx.bed_state.link_conc.len() != need {
        ctx.bed_state.link_conc.clear();
        ctx.bed_state.link_conc.resize(need, 0.0);
        ctx.bed_state.n_species = n_total;
    }
    let nl = ctx.n_links();
    let na = n_arr as usize;
    let nt = n_total as usize;
    let off = offset as usize;
    for j in 0..nl {
        let uj = j as usize;
        if uj < ctx.heat_state.mesh_contact_link.len() && ctx.heat_state.mesh_contact_link[uj] != 0 {
            continue;
        }
        if uj >= ctx.bed_state.link_temp.len() {
            break;
        }
        let vol_ft3 = ctx.links.volume[uj];
        if !(vol_ft3 > 0.0) {
            continue;
        }
        let bed_m2 = link_bed_area_m2(ctx, j);
        if !(bed_m2 > 0.0) {
            continue;
        }
        let sc = sediment_for(ctx, &HeatElement::link(j)).clone();
        let vol_w = vol_ft3 * K_CU_FT_TO_CU_M;
        let vol_b = bed_m2 * sc.bed_thickness;
        let q_exch = bed_exchange_q(&sc, bed_m2);
        if !(q_exch > 0.0) || !(vol_b > 0.0) {
            continue;
        }
        for s in 0..na {
            let wi = uj * na + s;
            let bi = uj * nt + off + s;
            let cw = if msx { ctx.reactions.msx_link_conc[wi] } else { ctx.links.conc[wi] };
            let cb = ctx.bed_state.link_conc[bi];
            let st = exchange_pair(cw, cb, vol_w, vol_b, q_exch, dt);
            if msx {
                ctx.reactions.msx_link_conc[wi] += st.dc_w;
            } else {
                ctx.links.conc[wi] += st.dc_w;
            }
            ctx.bed_state.link_conc[bi] += st.dc_b;
        }
    }
    Ok(())
}
