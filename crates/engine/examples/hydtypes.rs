// Rust side of gate/diff/hydtypes.cpp: the hydrology header types' initializers,
// constants, enum discriminants and inline functions, printed line for line as the C++
// harness prints them.
#![allow(clippy::excessive_precision, clippy::field_reassign_with_default)]
use std::mem::size_of;
use swmm6_rust::data::hydrology::*;
use swmm6_rust::hydrology::climate_types::*;
use swmm6_rust::hydrology::hourly_meteorology::*;
use swmm6_rust::hydrology::radiation::*;
use swmm6_rust::hydrology::ueb_types::*;

fn hx(v: f64) {
    print!(" {:016x}", v.to_bits());
}
fn ln(s: &str) {
    print!("{}", s);
}
fn nl() {
    println!();
}

fn d_parameters(p: &Parameters) {
    print!("Parameters {}", p.temperature as i32);
    hx(p.rain_temperature); hx(p.snow_temperature); hx(p.emissivity); hx(p.roughness);
    hx(p.measurement_height); hx(p.snow_density); hx(p.soil_density);
    hx(p.soil_heat_capacity); hx(p.active_soil_depth); hx(p.liquid_retention);
    hx(p.drainage_conductivity); hx(p.equilibrium_conductance); hx(p.snow_conductivity);
    hx(p.soil_conductivity); hx(p.low_frequency); hx(p.amplitude); hx(p.stability);
    hx(p.ground_albedo); hx(p.visible_albedo); hx(p.infrared_albedo);
    hx(p.new_snow_depth); hx(p.age_dust); hx(p.visible_age); hx(p.infrared_age);
    hx(p.max_step); nl();
}
fn d_site(s: &Site) {
    ln("Site"); hx(s.forest); hx(s.drift); hx(s.pressure); hx(s.ground_heat);
    hx(s.extinction_depth); nl();
}
fn d_state(s: &State) {
    ln("State"); hx(s.energy); hx(s.water); hx(s.age); hx(s.refreezing_depth);
    hx(s.total_refreezing_depth); print!(" {}", s.history.len()); nl();
}
fn d_fluxes(f: &Fluxes) {
    ln("Fluxes"); hx(f.rain); hx(f.snow); hx(f.liquid); hx(f.sublimation); hx(f.net_energy);
    hx(f.outgoing_melt_energy); hx(f.surface_temperature); hx(f.bulk_temperature);
    hx(f.albedo); hx(f.water_residual); hx(f.energy_residual);
    print!(" {} {} {}", f.steps, f.retries, f.diagnostics); nl();
}
fn d_lumped_state(s: &LumpedState) {
    d_state(&s.point); ln("LumpedState"); hx(s.water); hx(s.maximum); hx(s.cover); nl();
}
fn d_lumped_fluxes(f: &LumpedFluxes) {
    d_fluxes(&f.point); ln("LumpedFluxes"); hx(f.rain); hx(f.snow); hx(f.liquid);
    hx(f.throughfall); hx(f.sublimation); hx(f.active_fraction); hx(f.water_residual);
    hx(f.net_energy); hx(f.geometry_energy); hx(f.energy_residual); nl();
}
fn d_rad_source(r: &RadiationSource) {
    println!("RadiationSource {} {} {}", r.mode as i32, r.series.len(), r.authored as i32);
}
fn d_mesh_row(r: &MeshSnowRow) {
    println!("MeshSnowRow {} {} {} {} {} {} {}", r.pack.len(), r.forest.is_some() as i32,
             r.drift.is_some() as i32, r.ground_heat.is_some() as i32,
             r.extinction_depth.is_some() as i32, r.slope.is_some() as i32,
             r.aspect.is_some() as i32);
}
fn d_rad_state(r: &RadiationState) {
    ln("RadiationState"); hx(r.shortwave); hx(r.longwave); hx(r.cosine); hx(r.horizontal);
    hx(r.transmittance); hx(r.cloud); hx(r.temperature); hx(r.humidity); hx(r.wind);
    hx(r.pressure); hx(r.declination); hx(r.hour_angle); hx(r.latitude); hx(r.date);
    hx(r.seconds);
    print!(" {} {} {}", r.shortwave_source as i32, r.longwave_source as i32,
           r.shared_heat_forcing as i32);
    nl();
}
fn d_ma(m: &MovingAvg7) {
    ln("MovingAvg7");
    for v in m.ta { hx(v); }
    for v in m.tr { hx(v); }
    print!(" {} {}", m.front, m.count); nl();
}

fn main() {
    // ---- constants
    ln("const"); hx(K_SOLAR_CONSTANT); hx(K_STD_PRESSURE_MB); hx(MM_PER_INCH); nl();

    // ---- enum discriminants
    println!("enum TemperatureMethod {} {}", TemperatureMethod::ModifiedForceRestore as i32,
             TemperatureMethod::Equilibrium as i32);
    println!("enum SnowKind {} {}", SnowKind::DegreeDay as i32, SnowKind::Ueb as i32);
    println!("enum RadiationSourceMode {} {} {}", RadiationSourceMode::Computed as i32,
             RadiationSourceMode::Timeseries as i32, RadiationSourceMode::File as i32);
    println!("enum MeshSnowConfigMode {} {} {}", MeshSnowConfigMode::Auto as i32,
             MeshSnowConfigMode::Yes as i32, MeshSnowConfigMode::No as i32);
    println!("enum RadiationOrigin {} {} {} {}", RadiationOrigin::Computed as i32,
             RadiationOrigin::Timeseries as i32, RadiationOrigin::MetFile as i32,
             RadiationOrigin::LegacyHeat as i32);
    println!("enum EvapMethod {} {} {} {} {}", EvapMethod::Constant as i32,
             EvapMethod::Monthly as i32, EvapMethod::Timeseries as i32,
             EvapMethod::Temperature as i32, EvapMethod::Pan as i32);
    println!("sizeof enums {} {} {} {} {} {}", size_of::<TemperatureMethod>(),
             size_of::<SnowKind>(), size_of::<RadiationSourceMode>(),
             size_of::<MeshSnowConfigMode>(), size_of::<RadiationOrigin>(),
             size_of::<EvapMethod>());

    // ---- defaults, UEB point
    d_parameters(&Parameters::default());
    d_site(&Site::default());
    {
        let f = Forcing::default();
        ln("Forcing"); hx(f.seconds); hx(f.precipitation); hx(f.temperature); hx(f.wind);
        hx(f.humidity); hx(f.shortwave); hx(f.longwave); hx(f.cosine); hx(f.snow_catch);
        hx(f.snow_fraction); nl();
    }
    {
        let t = ThermalSample::default();
        ln("ThermalSample"); hx(t.seconds); hx(t.surface); hx(t.bulk); nl();
    }
    d_state(&State::default());
    d_fluxes(&Fluxes::default());
    d_lumped_state(&LumpedState::default());
    d_lumped_fluxes(&LumpedFluxes::default());
    {
        let h = HeldState::default();
        d_lumped_state(&h.accepted); d_lumped_state(&h.begin); d_lumped_state(&h.end);
        d_lumped_fluxes(&h.flux);
        ln("HeldState"); hx(h.duration); hx(h.elapsed); print!(" {}", h.pending() as i32); nl();
    }
    {
        let n = NativeState::default();
        ln("NativeState");
        for v in n.forcing { hx(v); }
        print!(" {} {} {} {} {} {} {}", n.columns.len(), n.last.len(), n.rain.len(),
               n.snow.len(), n.sample_end.len(), n.sample_seconds.len(), n.active as i32);
        hx(n.sublimation); hx(n.precipitation); hx(n.liquid); hx(n.throughfall);
        hx(n.net_energy); hx(n.geometry_energy); hx(n.plowed_energy); hx(n.initial_energy);
        hx(n.initial_water); hx(n.removed); hx(n.transferred_water); hx(n.last_seconds);
        hx(n.completed_seconds); nl();
    }

    // ---- defaults, configuration
    {
        let s = SnowpackConfig::default();
        print!("SnowpackConfig {}", s.model as i32); hx(s.initial_energy); nl();
        d_site(&s.site);
    }
    d_rad_source(&RadiationSource::default());
    {
        let c = ClimateSite::default();
        ln("ClimateSite"); hx(c.latitude); hx(c.longitude); hx(c.timezone); hx(c.elevation); nl();
    }
    {
        let c = SnowConfig::default();
        print!("SnowConfig {} {} {}", c.model as i32, c.model_authored as i32,
               c.division_temperature as i32);
        hx(c.step); hx(c.bc_a); hx(c.bc_c);
        for m in &c.monthly_range { print!(" {}", m.is_some() as i32); }
        print!(" {} {} {}", c.site.is_some() as i32, c.met_file.is_empty() as i32, c.packs.len());
        nl();
        d_parameters(&c.parameters);
        d_rad_source(&c.shortwave);
        d_rad_source(&c.longwave);
    }
    d_mesh_row(&MeshSnowRow::default());
    {
        let m = MeshSnowConfig::default();
        print!("MeshSnowConfig {}", m.mode as i32); hx(m.step);
        println!(" {} {} {}", m.defaults.len(), m.cells.len(), m.assigned() as i32);
    }

    // ---- defaults, climate
    d_rad_state(&RadiationState::default());
    {
        let s = SolarPosition::default();
        ln("SolarPosition"); hx(s.declination_deg); hx(s.eq_time_min); hx(s.hour_angle_deg);
        hx(s.zenith_deg); hx(s.cos_zenith); hx(s.elevation_deg); hx(s.eccentricity); nl();
    }
    d_ma(&MovingAvg7::default());
    {
        let c = ClimateState::default();
        d_rad_state(&c.radiation);
        ln("ClimateState"); hx(c.temperature); hx(c.temp_range); hx(c.evap_rate);
        hx(c.wind_speed); hx(c.humidity); hx(c.temperature_src); hx(c.file_tmin);
        hx(c.file_tmax); hx(c.wind_speed_src); hx(c.gamma); hx(c.ea); hx(c.latitude);
        hx(c.elev); hx(c.snow_divt); nl();
        ln("monthly_evap"); for v in c.monthly_evap { hx(v); } nl();
        ln("evaprate_ucf"); hx(c.evaprate_ucf); nl();
        ln("adjust_evap"); for v in c.adjust_evap { hx(v); } nl();
        ln("next_evap"); hx(c.next_evap_date); hx(c.next_evap_rate);
        print!(" {}", c.evap_ts_pos); nl();
        ln("adjust_temp"); for v in c.adjust_temp { hx(v); } nl();
        ln("adjust_rain"); for v in c.adjust_rain { hx(v); } nl();
        ln("adjust_hydcon"); for v in c.adjust_hydcon { hx(v); } nl();
        print!("evap_method {}", c.evap_method as i32); hx(c.infil_factor); nl();
        d_ma(&c.temp_ma);
        ln("recovery"); hx(c.recovery_factor);
        print!(" {} {} {} {}", c.temp_ts_index, c.evap_ts_index, c.humidity_ts_index,
               c.recovery_pat_index);
        nl();
        ln("subdaily"); hx(c.tmin_daily); hx(c.tmax_daily); hx(c.prev_tmax); hx(c.hrsr);
        hx(c.hrss); hx(c.hrday); hx(c.dhrdy); hx(c.dydif); hx(c.dtlong);
        print!(" {} {}", c.has_minmax as i32, c.last_temp_doy); nl();
    }
    {
        let r = HourlyMetRecord::default();
        ln("HourlyMetRecord"); hx(r.date); for v in r.values { hx(v); } nl();
        let w = SnowDailyWeather::default();
        ln("SnowDailyWeather"); hx(w.date); hx(w.tmin); hx(w.tmax); hx(w.previous_max);
        hx(w.wind); nl();
        let h = HourlyMeteorology::default();
        println!("HourlyMeteorology {} {}", h.is_empty() as i32, h.records().len());
    }

    // ---- HeldState::pending over a sweep (elapsed, duration)
    {
        let inf = f64::INFINITY;
        let nan = f64::NAN;
        let vs = [0.0, -0.0, 1.0, -1.0, 0.5, 1e-300, 4.9406564584124654e-324,
                  1.7976931348623157e308, inf, -inf, nan, -nan, 3600.0,
                  3600.0000000000005, 3599.9999999999995, 1e-5];
        ln("pending");
        for &e in &vs {
            for &d in &vs {
                let mut h = HeldState::default();
                h.elapsed = e;
                h.duration = d;
                print!("{}", h.pending() as i32);
            }
        }
        nl();
    }

    // ---- SnowConfig::model_for
    {
        let mut c = SnowConfig::default();
        let u = SnowpackConfig::default(); // default model Ueb
        let mut d = SnowpackConfig::default(); d.model = SnowKind::DegreeDay;
        c.packs.insert(0, u); c.packs.insert(2, d); c.packs.insert(5, u); c.packs.insert(-3, d);
        ln("modelFor DD");
        for p in -4..=7 { print!(" {}", c.model_for(p) as i32); }
        nl();
        c.model = SnowKind::Ueb;
        ln("modelFor UEB");
        for p in -4..=7 { print!(" {}", c.model_for(p) as i32); }
        nl();
        ln("packs order");
        for (k, v) in &c.packs { print!(" {}:{}", k, v.model as i32); }
        nl();
    }

    // ---- MeshSnowConfig::assigned and BTreeMap key order
    {
        let mut m = MeshSnowConfig::default();
        print!("assigned {}", m.assigned() as i32);
        m.defaults.insert(b"*".to_vec(), MeshSnowRow::default());
        print!(" {}", m.assigned() as i32);
        m.defaults.clear();
        print!(" {}", m.assigned() as i32);
        m.cells.insert(3, MeshSnowRow::default());
        print!(" {}", m.assigned() as i32);
        m.mode = MeshSnowConfigMode::No;
        println!(" {}", m.assigned() as i32);
        let keys: [&[u8]; 14] = [b"b", b"a", b"A", b"*", b"", b"\x80", b"\xff", b"ab", b"a\0b",
                                 b"a", b"Z", b"\x7f", b"~", b"aa"];
        for k in keys { m.defaults.insert(k.to_vec(), MeshSnowRow::default()); }
        ln("defaults order");
        for k in m.defaults.keys() {
            print!(" [");
            for ch in k { print!("{:02x}", ch); }
            print!("]");
        }
        nl();
        let cks = [7, -1, 0, 2147483647, -2147483647 - 1, 3, 12];
        for k in cks { m.cells.insert(k, MeshSnowRow::default()); }
        ln("cells order");
        for k in m.cells.keys() { print!(" {}", k); }
        nl();
    }

    // ---- data/hydrology.rs
    {
        let mut s = SnowpackStore::default();
        print!("SnowpackStore {}", s.count());
        s.names.push(b"A".to_vec()); s.names.push(b"B".to_vec()); s.names.push(b"C".to_vec());
        println!(" {} {} {} {} {} {} {} {}", s.count(), s.report_values.len(),
                 s.daily_weather.len(), s.plowable.len(), s.impervious.len(),
                 s.pervious.len(), s.removal.len(), s.removal_subcatch.len());
        print!("SnowpackStore.config {}", s.config.model as i32); hx(s.config.step); nl();
        println!("SnowpackStore.mesh {} {} {}", s.mesh_config.mode as i32,
                 s.mesh_config.assigned() as i32, s.meteorology.is_empty() as i32);
        let mut a = AquiferStore::default();
        print!("AquiferStore {}", a.count());
        for _ in 0..4 { a.names.push(b"x".to_vec()); }
        println!(" {}", a.count());
        let mut l = LidControlStore::default();
        print!("LidControlStore {}", l.count());
        for _ in 0..2 { l.names.push(b"y".to_vec()); }
        println!(" {}", l.count());
        let mut u = LidUsageStore::default();
        print!("LidUsageStore {}", u.count());
        for i in 0..5 { u.subcatch_index.push(i); }
        println!(" {}", u.count());
        let ns = [3, 0, 5, 2, 1, 6];
        for n in ns {
            if !u.wb_inflow.is_empty() {
                u.wb_inflow[0] = 7.25; u.wb_evap[0] = -1.5; *u.wb_final_vol.last_mut().unwrap() = 3.0;
            }
            u.resize_wb(n);
            print!("resize_wb {} {} {} {} {} {} {} {}", n, u.wb_inflow.len(),
                   u.wb_evap.len(), u.wb_infil.len(), u.wb_surf_flow.len(),
                   u.wb_drain_flow.len(), u.wb_init_vol.len(), u.wb_final_vol.len());
            for &v in &u.wb_inflow { hx(v); }
            for &v in &u.wb_evap { hx(v); }
            for &v in &u.wb_final_vol { hx(v); }
            nl();
        }
        println!("LidUsageStore other {} {} {}", u.area.len(), u.rpt_file.len(),
                 u.from_perv.len());
    }
}
