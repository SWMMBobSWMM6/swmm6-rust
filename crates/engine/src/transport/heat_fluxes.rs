// SPDX-License-Identifier: Apache-2.0
//
// Derived from OpenSWMM Engine src/engine/transport/components/HeatFluxModules/SurfaceExchange.cpp
// and .hpp, RadiativeExchange.cpp and .hpp, HeatFluxes.cpp and .hpp, SolarRadiation.hpp, and the
// solar-forcing half of src/engine/hydrology/Radiation.cpp (solarPosition, airMass,
// pressureFromElevation, birdClearSkyGHI, cloudShortwaveFactor, cloudLongwaveFactor,
// updateSolarForcing) (HydroCouple/openswmm.engine, branch swmm6_rel, commit 88329d6),
// Copyright 2026 Caleb Buahin, licensed under the Apache License, Version 2.0.
// Modified: translated from C++ to Rust, 6 November 2026, by Claude (Anthropic) at the request
// of Robert Dickinson. Not an OpenSWMM release.
//
// The heat flux modules (namespace openswmm::transport::heat): the surface energy balance
// (evaporation and conduction, Bowen ratio), the radiative balance (shortwave, back, atmospheric
// and land-cover longwave), their sum (netFluxOut), the semi-implicit relaxation (relaxT), the
// per-step solar forcing, and applyHeatFluxes for the LEGACY temperature mirror's nodes and
// links.
//
// Differences in form, never in result:
//   * updateSolarForcing can throw (resolveRadiation's error); here it returns
//     Err(StepError::Exception) for that throw. NOT YET TRANSLATED: the branch that resolves
//     the snow configuration's authored shortwave / longwave / met file (resolveRadiation)
//     returns StepError::NotTranslated.
//   * applyHeatFluxes' bed pair (SEDIMENT_EXCHANGE: bedCouplingForLink, relaxPair,
//     seedBedTemperature, groundTemperature) is in heat_bed.rs.

use crate::core::cmath::{cmax, cmin};
use crate::core::context::SimulationContext;
use crate::core::datetime;
use crate::core::engine_step_routing::StepError;
use crate::core::unit_conversion::get_unit_system;
use crate::data::heat::{RadiativeConfig, ShortwaveMode, SolarConfig};
use crate::data::heat_override::HeatElement;
use crate::data::table::table_tseries_lookup_cursor;
use crate::hydraulics::link as link_mod;
use crate::hydraulics::node::{get_surf_area, store_volume};
use crate::hydraulics::xsection as xsect;
use crate::transport::heat_bed;
use crate::transport::heat_overrides::radiative_for;

// ---------------------------------------------------------------------------
// SurfaceExchange
// ---------------------------------------------------------------------------

pub const K_BOWEN_COEFF: f64 = 0.061;
pub const K_PROBE_C: f64 = 1.0e-3;
pub const K_SQ_FT_TO_SQ_M: f64 = 0.09290304;
pub const K_CU_FT_TO_CU_M: f64 = 0.028316846592;
pub const K_MPH_TO_MS: f64 = 0.44704;

fn net_flux_out_surface(t_w: f64, t_air: f64, rh: f64, wind_ms: f64, a: f64, b: f64, rho: f64, p_ratio: f64) -> f64 {
    let je = latent_flux(t_w, t_air, rh, wind_ms, a, b, rho);
    je + sensible_flux(je, bowen_ratio(t_w, t_air, rh, p_ratio))
}

/// airTempCelsius: the climate state's air temperature (deg F) in deg C.
pub fn air_temp_celsius(ctx: &SimulationContext) -> f64 {
    (ctx.climate_state.temperature - 32.0) * 5.0 / 9.0
}

/// relaxT: the semi-implicit temperature change of one well-mixed body over dt.
#[allow(clippy::too_many_arguments)]
pub fn relax_t(j0: f64, j1: f64, h: f64, area_m2: f64, vol_m3: f64, dt: f64, rho: f64, cp: f64) -> f64 {
    let heat_capacity = rho * cp * vol_m3;
    if !(heat_capacity > 0.0) || !(area_m2 > 0.0) || !(dt > 0.0) {
        return 0.0;
    }
    let explicit_dt = -j0 * area_m2 * dt / heat_capacity;
    if !(h > 0.0) {
        return explicit_dt;
    }
    let djdt = (j1 - j0) / h;
    if !(djdt > 0.0) {
        return explicit_dt;
    }
    let k = area_m2 * djdt / heat_capacity;
    if !(k > 0.0) {
        return explicit_dt;
    }
    let d_t = (j0 / djdt) * (-k * dt).exp_m1(); // LIBM
    if d_t.is_finite() {
        d_t
    } else {
        explicit_dt
    }
}

/// equilibriumT
pub fn equilibrium_t(t0_c: f64, j0: f64, j1: f64, h: f64) -> f64 {
    if !(h > 0.0) {
        return t0_c;
    }
    let djdt = (j1 - j0) / h;
    if !(djdt > 0.0) {
        return t0_c;
    }
    t0_c - j0 / djdt
}

pub fn saturation_vapour_pressure(t_c: f64) -> f64 {
    0.61275 * (17.27 * t_c / (237.3 + t_c)).exp() // LIBM
}

pub fn latent_heat_of_vaporization(t_c: f64) -> f64 {
    1000.0 * (2499.0 - 2.36 * t_c)
}

pub fn wind_function(wind_ms: f64, a: f64, b: f64) -> f64 {
    a + b * cmax(wind_ms, 0.0)
}

pub fn evaporation_rate(t_water_c: f64, t_air_c: f64, humidity_pct: f64, wind_ms: f64, a: f64, b: f64) -> f64 {
    let es_w = saturation_vapour_pressure(t_water_c);
    let e_a = (humidity_pct / 100.0) * saturation_vapour_pressure(t_air_c);
    wind_function(wind_ms, a, b) * (es_w - e_a)
}

#[allow(clippy::too_many_arguments)]
pub fn latent_flux(t_water_c: f64, t_air_c: f64, humidity_pct: f64, wind_ms: f64, a: f64, b: f64, water_density: f64) -> f64 {
    let e = evaporation_rate(t_water_c, t_air_c, humidity_pct, wind_ms, a, b);
    water_density * latent_heat_of_vaporization(t_water_c) * e
}

pub fn bowen_ratio(t_water_c: f64, t_air_c: f64, humidity_pct: f64, pressure_ratio: f64) -> f64 {
    let es_w = saturation_vapour_pressure(t_water_c);
    let e_a = (humidity_pct / 100.0) * saturation_vapour_pressure(t_air_c);
    let deficit = es_w - e_a;
    if deficit.abs() < 1.0e-12 {
        return 0.0;
    }
    K_BOWEN_COEFF * pressure_ratio * (t_water_c - t_air_c) / deficit
}

pub fn sensible_flux(latent_flux: f64, bowen_ratio: f64) -> f64 {
    bowen_ratio * latent_flux
}

/// surfaceFluxOut: the surface energy balance, positive out of the water (W/m2).
pub fn surface_flux_out(ctx: &SimulationContext, elem: &HeatElement, t_w: f64) -> f64 {
    if !ctx.options.heat_transport || !ctx.heat_config.surface_exchange {
        return 0.0;
    }
    let t_air = ctx.forcing.element_air_temp_f(elem, ctx.climate_state.temperature);
    let rh = ctx.forcing.element_humidity(elem, ctx.climate_state.humidity);
    let wind = ctx.forcing.element_wind_mph(elem, ctx.climate_state.wind_speed);
    net_flux_out_surface(t_w, (t_air - 32.0) * 5.0 / 9.0, rh, wind * K_MPH_TO_MS, ctx.options.wind_func_coeff_a, ctx.options.wind_func_coeff_b, ctx.options.water_density, ctx.options.pressure_ratio)
}

// ---------------------------------------------------------------------------
// RadiativeExchange
// ---------------------------------------------------------------------------

pub const K_STEFAN_BOLTZMANN: f64 = 5.67e-8;
pub const K_KELVIN_OFFSET: f64 = 273.15;

fn kelvin4(t_c: f64) -> f64 {
    let k = t_c + K_KELVIN_OFFSET;
    k * k * k * k
}

pub fn net_shortwave(incoming_wm2: f64, albedo: f64, shade_factor: f64) -> f64 {
    (1.0 - albedo) * incoming_wm2 * cmax(0.0, 1.0 - shade_factor)
}

pub fn back_longwave(t_water_c: f64, emiss_water: f64) -> f64 {
    emiss_water * K_STEFAN_BOLTZMANN * kelvin4(t_water_c)
}

pub fn atmospheric_emissivity(e_a_kpa: f64, atm_emiss_coeff: f64, cloud_factor: f64) -> f64 {
    let clear = atm_emiss_coeff + 0.0027 * (cmax(0.0, e_a_kpa) * 1000.0).sqrt();
    if cloud_factor == 1.0 {
        return clear;
    }
    clear * cloud_factor
}

pub fn atmospheric_longwave(t_air_c: f64, humidity_pct: f64, atm_emiss_coeff: f64, lw_reflection: f64, sky_view: f64, cloud_factor: f64) -> f64 {
    let e_a = (humidity_pct / 100.0) * saturation_vapour_pressure(t_air_c);
    K_STEFAN_BOLTZMANN * kelvin4(t_air_c) * atmospheric_emissivity(e_a, atm_emiss_coeff, cloud_factor) * (1.0 - lw_reflection) * sky_view
}

pub fn land_cover_longwave(t_air_c: f64, emiss_landcover: f64, sky_view: f64) -> f64 {
    emiss_landcover * K_STEFAN_BOLTZMANN * (1.0 - sky_view) * kelvin4(t_air_c)
}

/// netRadiativeFluxOut (jin_wm2 < 0: the configured constant; incoming_longwave < 0: computed).
pub fn net_radiative_flux_out(t_water_c: f64, t_air_c: f64, humidity_pct: f64, cfg: &RadiativeConfig, jin_wm2: f64, cloud_factor: f64, incoming_longwave: f64) -> f64 {
    let jin = if jin_wm2 < 0.0 { cfg.shortwave_wm2 } else { jin_wm2 };
    let jsn = net_shortwave(jin, cfg.albedo, cfg.shade_factor);
    let jbr = back_longwave(t_water_c, cfg.emiss_water);
    let jan = if incoming_longwave >= 0.0 { incoming_longwave * (1.0 - cfg.lw_reflection) * cfg.sky_view } else { atmospheric_longwave(t_air_c, humidity_pct, cfg.atm_emiss_coeff, cfg.lw_reflection, cfg.sky_view, cloud_factor) };
    let t_lc_c = if cfg.landcover_temp.is_nan() { t_air_c } else { cfg.landcover_temp };
    let jlc = land_cover_longwave(t_lc_c, cfg.emiss_landcover, cfg.sky_view);
    jbr - jsn - jan - jlc
}

/// radiativeFluxOut(ctx, elem, t_w): enabled by [HEAT_FLUXES] RADIATIVE_EXCHANGE.
pub fn radiative_flux_out(ctx: &SimulationContext, elem: &HeatElement, t_w: f64) -> f64 {
    radiative_flux_out_enabled(ctx, elem, t_w, ctx.heat_config.radiative_exchange)
}

/// radiativeFluxOut(ctx, elem, t_w, enabled)
pub fn radiative_flux_out_enabled(ctx: &SimulationContext, elem: &HeatElement, t_w: f64, enabled: bool) -> f64 {
    if !ctx.options.heat_transport || !enabled {
        return 0.0;
    }
    let cc = &ctx.heat_config.cloud;
    let rc = radiative_for(ctx, elem);
    let t_air_c = (ctx.forcing.element_air_temp_f(elem, ctx.climate_state.temperature) - 32.0) * 5.0 / 9.0;
    let rh = ctx.forcing.element_humidity(elem, ctx.climate_state.humidity);
    let jin = ctx.forcing.element_shortwave(elem, ctx.heat_state.shortwave_now);
    net_radiative_flux_out(t_w, t_air_c, rh, rc, jin, cloud_longwave_factor(ctx.heat_state.cloud_now, cc.lw_cloud_k), if ctx.climate_state.radiation.shared_heat_forcing { ctx.climate_state.radiation.longwave } else { -1.0 })
}

// ---------------------------------------------------------------------------
// HeatFluxes
// ---------------------------------------------------------------------------

/// netFluxOut: every flux family, positive out of the water.
pub fn net_flux_out(ctx: &SimulationContext, elem: &HeatElement, t_w: f64) -> f64 {
    surface_flux_out(ctx, elem, t_w) + radiative_flux_out(ctx, elem, t_w)
}

/// applyHeatFluxes: the LEGACY mirror's nodes (storage surfaces) and links (open conduits,
/// and the bed pair under SEDIMENT_EXCHANGE).
pub fn apply_heat_fluxes(ctx: &mut SimulationContext, dt: f64) -> Result<(), StepError> {
    if !ctx.options.heat_transport {
        return Ok(());
    }
    let bed_on = heat_bed::bed_exchange_enabled(ctx);
    if !ctx.heat_config.surface_exchange && !ctx.heat_config.radiative_exchange && !bed_on {
        return Ok(());
    }
    if !(dt > 0.0) {
        return Ok(());
    }
    let mut t_gr = 0.0;
    if bed_on {
        heat_bed::seed_bed_temperature(ctx)?;
        t_gr = heat_bed::ground_temperature(ctx)?;
    }
    update_solar_forcing(ctx, dt)?;
    let rho = ctx.options.water_density;
    let cp = ctx.options.water_specific_heat;
    let unit_sys = get_unit_system(ctx.options.flow_units as i32);
    let step = |ctx: &SimulationContext, e: &HeatElement, t_w: f64, area_ft2: f64, vol_ft3: f64| -> f64 { relax_t(net_flux_out(ctx, e, t_w), net_flux_out(ctx, e, t_w + K_PROBE_C), K_PROBE_C, area_ft2 * K_SQ_FT_TO_SQ_M, vol_ft3 * K_CU_FT_TO_CU_M, dt, rho, cp) };

    // Nodes: storage surfaces only (getSurfArea is 0 elsewhere)
    let nn = ctx.n_nodes();
    for i in 0..nn {
        let ui = i as usize;
        if ui >= ctx.heat_state.node_temp.len() {
            break;
        }
        let vol_ft3 = store_volume(&ctx.nodes, i);
        if !(vol_ft3 > 0.0) {
            continue;
        }
        let area_ft2 = get_surf_area(&ctx.nodes, i, ctx.nodes.depth[ui], Some(&ctx.tables), unit_sys, Some(&ctx.node_subtypes));
        if !(area_ft2 > 0.0) {
            continue;
        }
        let d = step(ctx, &HeatElement::node(i), ctx.heat_state.node_temp[ui], area_ft2, vol_ft3);
        ctx.heat_state.node_temp[ui] += d;
    }

    // Links: open conduits (top width x length x barrels); the bed pair for every conduit
    let nl = ctx.n_links();
    for j in 0..nl {
        let uj = j as usize;
        if uj < ctx.heat_state.mesh_contact_link.len() && ctx.heat_state.mesh_contact_link[uj] != 0 {
            continue;
        }
        if uj >= ctx.heat_state.link_temp.len() {
            break;
        }
        let vol_ft3 = ctx.links.volume[uj];
        if !(vol_ft3 > 0.0) {
            continue;
        }
        let cr = ctx.link_subtypes.conduit_row(j);
        if cr < 0 {
            continue;
        }
        let mut area_ft2 = 0.0;
        if xsect::is_open(ctx.links.xsect_batch_shape[uj]) {
            let ucr = cr as usize;
            let cd = &ctx.link_subtypes.conduits;
            let mut length = cd.length[ucr];
            if !(length > 0.0) {
                length = cd.mod_length[ucr];
            }
            let depth = ctx.links.depth[uj];
            if length > 0.0 && depth > 0.0 {
                let xs = build_xsp(ctx, uj);
                let top_width = xsect::get_w_of_y(&xs, depth);
                if top_width > 0.0 {
                    area_ft2 = top_width * length * cd.barrels[ucr] as f64;
                }
            }
        }
        if bed_on && uj < ctx.bed_state.link_temp.len() {
            let g = heat_bed::bed_coupling_for_link(ctx, j, vol_ft3, t_gr);
            if g.viable() {
                let t_w = ctx.heat_state.link_temp[uj];
                let le = HeatElement::link(j);
                let s = heat_bed::relax_pair(&g, t_w, ctx.bed_state.link_temp[uj], net_flux_out(ctx, &le, t_w), net_flux_out(ctx, &le, t_w + K_PROBE_C), K_PROBE_C, area_ft2 * K_SQ_FT_TO_SQ_M, dt);
                ctx.heat_state.link_temp[uj] += s.dt_w;
                ctx.bed_state.link_temp[uj] += s.dt_b;
                continue;
            }
        }
        if !(area_ft2 > 0.0) {
            continue;
        }
        let d = step(ctx, &HeatElement::link(j), ctx.heat_state.link_temp[uj], area_ft2, vol_ft3);
        ctx.heat_state.link_temp[uj] += d;
    }
    Ok(())
}

/// HeatFluxes.cpp buildXsp: the link store's section parameters (no transect tables).
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

// ---------------------------------------------------------------------------
// Radiation.cpp: solar position, Bird clear sky, cloud, the per-step forcing
// ---------------------------------------------------------------------------

pub const K_SOLAR_CONSTANT: f64 = 1367.0;
pub const K_STD_PRESSURE_MB: f64 = 1013.25;
const K_PI: f64 = 3.14159265358979323846;
const K_DEG_RAD: f64 = K_PI / 180.0;
const K_RAD_DEG: f64 = 180.0 / K_PI;
const K_FT_TO_M: f64 = 0.3048;

/// climate::SolarPosition
#[derive(Debug, Clone, Copy)]
pub struct SolarPosition {
    pub declination_deg: f64,
    pub eq_time_min: f64,
    pub hour_angle_deg: f64,
    pub zenith_deg: f64,
    pub cos_zenith: f64,
    pub elevation_deg: f64,
    pub eccentricity: f64,
}

impl Default for SolarPosition {
    fn default() -> Self {
        SolarPosition { declination_deg: 0.0, eq_time_min: 0.0, hour_angle_deg: 0.0, zenith_deg: 90.0, cos_zenith: 0.0, elevation_deg: 0.0, eccentricity: 1.0 }
    }
}

/// solarPosition: Spencer (1971), NOAA's form.
pub fn solar_position(day_of_year: i32, hour_local: f64, lat_deg: f64, lon_deg: f64, tz_hours: f64) -> SolarPosition {
    let mut p = SolarPosition::default();
    let g = 2.0 * K_PI / 365.0 * (day_of_year as f64 - 1.0 + (hour_local - 12.0) / 24.0);
    // LIBM (cos, sin, acos)
    let c1 = g.cos();
    let s1 = g.sin();
    let c2 = (2.0 * g).cos();
    let s2 = (2.0 * g).sin();
    let c3 = (3.0 * g).cos();
    let s3 = (3.0 * g).sin();
    p.eq_time_min = 229.18 * (0.000075 + 0.001868 * c1 - 0.032077 * s1 - 0.014615 * c2 - 0.040849 * s2);
    let decl = 0.006918 - 0.399912 * c1 + 0.070257 * s1 - 0.006758 * c2 + 0.000907 * s2 - 0.002697 * c3 + 0.001480 * s3;
    p.declination_deg = decl * K_RAD_DEG;
    p.eccentricity = 1.00011 + 0.034221 * c1 + 0.001280 * s1 + 0.000719 * c2 + 0.000077 * s2;
    let time_offset_min = p.eq_time_min + 4.0 * lon_deg - 60.0 * tz_hours;
    let tst_min = hour_local * 60.0 + time_offset_min;
    let mut ha = tst_min / 4.0 - 180.0;
    while ha < -180.0 {
        ha += 360.0;
    }
    while ha > 180.0 {
        ha -= 360.0;
    }
    p.hour_angle_deg = ha;
    let lat = lat_deg * K_DEG_RAD;
    let cz = lat.sin() * decl.sin() + lat.cos() * decl.cos() * (ha * K_DEG_RAD).cos();
    let cz_clamped = cmax(-1.0, cmin(1.0, cz));
    p.zenith_deg = cz_clamped.acos() * K_RAD_DEG;
    p.elevation_deg = 90.0 - p.zenith_deg;
    p.cos_zenith = cmax(0.0, cz_clamped);
    p
}

/// airMass: Kasten & Young (1989).
pub fn air_mass(zenith_deg: f64) -> f64 {
    if zenith_deg >= 90.0 {
        return 0.0;
    }
    let z = zenith_deg;
    let denom = (z * K_DEG_RAD).cos() + 0.50572 * (96.07995 - z).powf(-1.6364); // LIBM
    if !(denom > 0.0) {
        return 0.0;
    }
    1.0 / denom
}

pub fn pressure_from_elevation(elevation_m: f64) -> f64 {
    K_STD_PRESSURE_MB * (1.0 - 2.25577e-5 * elevation_m).powf(5.25588) // LIBM
}

/// birdClearSkyGHI: Bird & Hulstrom (1981).
pub fn bird_clear_sky_ghi(pos: &SolarPosition, pressure_mb: f64, cfg: &SolarConfig) -> f64 {
    // LIBM throughout (exp, pow)
    let cz = pos.cos_zenith;
    if !(cz > 0.0) {
        return 0.0;
    }
    let am = air_mass(pos.zenith_deg);
    if !(am > 0.0) {
        return 0.0;
    }
    let amp = am * (pressure_mb / K_STD_PRESSURE_MB);
    let t_r = (-0.0903 * amp.powf(0.84) * (1.0 + amp - amp.powf(1.01))).exp();
    let xo = cfg.ozone_cm * am;
    let t_o = 1.0 - 0.1611 * xo * (1.0 + 139.48 * xo).powf(-0.3034) - 0.002715 * xo / (1.0 + 0.044 * xo + 0.0003 * xo * xo);
    let t_um = (-0.0127 * amp.powf(0.26)).exp();
    let xw = cfg.precip_water_cm * am;
    let t_w = 1.0 - 2.4959 * xw / ((1.0 + 79.034 * xw).powf(0.6828) + 6.385 * xw);
    let tau_a = 0.2758 * cfg.aod380 + 0.35 * cfg.aod500;
    let t_a = (-tau_a.powf(0.873) * (1.0 + tau_a - tau_a.powf(0.7088)) * am.powf(0.9108)).exp();
    const K_BA: f64 = 0.85;
    let t_aa = 1.0 - 0.1 * (1.0 - am + am.powf(1.06)) * (1.0 - t_a);
    let t_as = if t_aa > 0.0 { t_a / t_aa } else { 0.0 };
    let i_ext = K_SOLAR_CONSTANT * pos.eccentricity;
    let i_dn = 0.9662 * i_ext * t_r * t_o * t_um * t_w * t_a;
    let i_d = i_dn * cz;
    let denom = 1.0 - am + am.powf(1.02);
    let i_as = if denom > 0.0 { i_ext * cz * 0.79 * t_o * t_um * t_w * t_aa * (0.5 * (1.0 - t_r) + K_BA * (1.0 - t_as)) / denom } else { 0.0 };
    let rs = 0.0685 + (1.0 - K_BA) * (1.0 - t_as);
    let refl_denom = 1.0 - cfg.ground_albedo * rs;
    let ghi = if refl_denom > 0.0 { (i_d + i_as) / refl_denom } else { i_d + i_as };
    cmax(0.0, ghi)
}

pub fn cloud_shortwave_factor(cloud_fraction: f64, k: f64, n: f64) -> f64 {
    if !(cloud_fraction > 0.0) {
        return 1.0;
    }
    let c = cmin(1.0, cloud_fraction);
    cmax(0.0, 1.0 - k * c.powf(n)) // LIBM
}

pub fn cloud_longwave_factor(cloud_fraction: f64, k_lw: f64) -> f64 {
    if !(cloud_fraction > 0.0) {
        return 1.0;
    }
    let c = cmin(1.0, cloud_fraction);
    1.0 + k_lw * c * c
}

/// updateSolarForcing(ctx, seconds): this step's incoming shortwave and cloud fraction.
pub fn update_solar_forcing(ctx: &mut SimulationContext, seconds: f64) -> Result<(), StepError> {
    update_solar_forcing_mesh(ctx, seconds, false)
}

/// updateSolarForcing(ctx, seconds, mesh_enabled)
pub fn update_solar_forcing_mesh(ctx: &mut SimulationContext, _seconds: f64, mesh_enabled: bool) -> Result<(), StepError> {
    if !ctx.options.heat_transport || (!ctx.heat_config.radiative_exchange && !mesh_enabled) {
        ctx.heat_state.shortwave_now = -1.0;
        ctx.heat_state.cloud_now = 0.0;
        return Ok(());
    }
    {
        let snow = &ctx.snowpacks.config;
        if snow.shortwave.authored || snow.longwave.authored || !snow.met_file.original.is_empty() {
            return Err(StepError::NotTranslated("heat::updateSolarForcing (resolveRadiation from the snow configuration)"));
        }
    }
    ctx.climate_state.radiation.shared_heat_forcing = false;

    // Cloud fraction first
    let mut cloud = 0.0;
    {
        let cc = &ctx.heat_config.cloud;
        if cc.configured {
            if cc.use_timeseries {
                let ts = cc.ts_index;
                if ts >= 0 && ts < ctx.tables.count() as i32 {
                    let d = ctx.current_date;
                    cloud = table_tseries_lookup_cursor(&mut ctx.tables[ts], d);
                }
            } else {
                cloud = cc.fraction;
            }
        }
    }
    ctx.heat_state.cloud_now = cmax(0.0, cmin(1.0, cloud));

    // Incoming shortwave
    let mut jin = 0.0;
    let mode = ctx.heat_config.radiative.sw_mode;
    match mode {
        ShortwaveMode::Constant => {
            jin = ctx.heat_config.radiative.shortwave_wm2;
        }
        ShortwaveMode::Timeseries => {
            let ts = ctx.heat_config.radiative.sw_ts_index;
            if ts >= 0 && ts < ctx.tables.count() as i32 {
                let d = ctx.current_date;
                jin = table_tseries_lookup_cursor(&mut ctx.tables[ts], d);
            }
        }
        ShortwaveMode::Computed => {
            let sc = &ctx.heat_config.solar;
            let site = ctx.snowpacks.config.site.as_ref();
            if site.is_none() && (!sc.has_latitude || !sc.has_longitude) {
                jin = 0.0;
            } else {
                let doy = datetime::day_of_year(ctx.current_date) as f64;
                let day_frac = ctx.current_date - ctx.current_date.floor();
                let hour_local = day_frac * 24.0;
                let pos = solar_position(doy as i32, hour_local, site.map(|s| s.latitude).unwrap_or(sc.latitude_deg), site.map(|s| s.longitude).unwrap_or(sc.longitude_deg), site.map(|s| s.timezone).unwrap_or(sc.timezone_hours));
                let elev_m = match site {
                    Some(s) => s.elevation,
                    None => {
                        if sc.has_elevation {
                            sc.elevation_m
                        } else {
                            ctx.climate_state.elev * K_FT_TO_M
                        }
                    }
                };
                jin = bird_clear_sky_ghi(&pos, pressure_from_elevation(elev_m), sc);
            }
        }
    }
    let cc = &ctx.heat_config.cloud;
    jin *= cloud_shortwave_factor(ctx.heat_state.cloud_now, cc.sw_atten_k, cc.sw_atten_n);
    ctx.heat_state.shortwave_now = cmax(0.0, jin);
    Ok(())
}
