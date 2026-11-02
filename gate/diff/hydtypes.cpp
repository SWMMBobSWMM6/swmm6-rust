// Differential for the hydrology headers (UebSnow/UebLumped/UebHeld/UebConfig/
// MeshSnowConfig/UebSubcatchments, Radiation, Climate, HourlyMeteorology) and
// data/HydrologyData.hpp: every member initializer and constant as hex bits, enum
// discriminants, and each inline function (HeldState::pending, SnowConfig::modelFor,
// MeshSnowConfig::assigned, std::map key order, the four count()s, resize_wb) over a
// sweep of inputs. Rust side: crates/engine/examples/hydtypes.rs.
// Build: g++ -std=c++20 -O2 -ffp-contract=off -I$U -I$U/data -I$U/core hydtypes.cpp
//        (U = upstream src/engine)
#include "data/HydrologyData.hpp"
#include "hydrology/Climate.hpp"
#include "hydrology/Radiation.hpp"
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
using namespace openswmm;
namespace ueb = openswmm::surface::ueb;
namespace srf = openswmm::surface;
namespace cl = openswmm::climate;

static void hx(double v) {
    uint64_t u;
    std::memcpy(&u, &v, 8);
    std::printf(" %016llx", (unsigned long long)u);
}
static void ln(const char* s) { std::printf("%s", s); }
static void nl() { std::printf("\n"); }

static void dParameters(const ueb::Parameters& p) {
    std::printf("Parameters %d", (int)p.temperature);
    hx(p.rain_temperature); hx(p.snow_temperature); hx(p.emissivity); hx(p.roughness);
    hx(p.measurement_height); hx(p.snow_density); hx(p.soil_density);
    hx(p.soil_heat_capacity); hx(p.active_soil_depth); hx(p.liquid_retention);
    hx(p.drainage_conductivity); hx(p.equilibrium_conductance); hx(p.snow_conductivity);
    hx(p.soil_conductivity); hx(p.low_frequency); hx(p.amplitude); hx(p.stability);
    hx(p.ground_albedo); hx(p.visible_albedo); hx(p.infrared_albedo);
    hx(p.new_snow_depth); hx(p.age_dust); hx(p.visible_age); hx(p.infrared_age);
    hx(p.max_step); nl();
}
static void dSite(const ueb::Site& s) {
    ln("Site"); hx(s.forest); hx(s.drift); hx(s.pressure); hx(s.ground_heat);
    hx(s.extinction_depth); nl();
}
static void dState(const ueb::State& s) {
    ln("State"); hx(s.energy); hx(s.water); hx(s.age); hx(s.refreezing_depth);
    hx(s.total_refreezing_depth); std::printf(" %zu", s.history.size()); nl();
}
static void dFluxes(const ueb::Fluxes& f) {
    ln("Fluxes"); hx(f.rain); hx(f.snow); hx(f.liquid); hx(f.sublimation); hx(f.net_energy);
    hx(f.outgoing_melt_energy); hx(f.surface_temperature); hx(f.bulk_temperature);
    hx(f.albedo); hx(f.water_residual); hx(f.energy_residual);
    std::printf(" %d %d %d", f.steps, f.retries, f.diagnostics); nl();
}
static void dLumpedState(const ueb::LumpedState& s) {
    dState(s.point); ln("LumpedState"); hx(s.water); hx(s.maximum); hx(s.cover); nl();
}
static void dLumpedFluxes(const ueb::LumpedFluxes& f) {
    dFluxes(f.point); ln("LumpedFluxes"); hx(f.rain); hx(f.snow); hx(f.liquid);
    hx(f.throughfall); hx(f.sublimation); hx(f.active_fraction); hx(f.water_residual);
    hx(f.net_energy); hx(f.geometry_energy); hx(f.energy_residual); nl();
}
static void dRadSource(const srf::RadiationSource& r) {
    std::printf("RadiationSource %d %zu %d\n", (int)r.mode, r.series.size(), (int)r.authored);
}
static void dMeshRow(const srf::MeshSnowRow& r) {
    std::printf("MeshSnowRow %zu %d %d %d %d %d %d\n", r.pack.size(), (int)r.forest.has_value(),
                (int)r.drift.has_value(), (int)r.ground_heat.has_value(),
                (int)r.extinction_depth.has_value(), (int)r.slope.has_value(),
                (int)r.aspect.has_value());
}
static void dRadState(const cl::RadiationState& r) {
    ln("RadiationState"); hx(r.shortwave); hx(r.longwave); hx(r.cosine); hx(r.horizontal);
    hx(r.transmittance); hx(r.cloud); hx(r.temperature); hx(r.humidity); hx(r.wind);
    hx(r.pressure); hx(r.declination); hx(r.hour_angle); hx(r.latitude); hx(r.date);
    hx(r.seconds);
    std::printf(" %d %d %d", (int)r.shortwave_source, (int)r.longwave_source,
                (int)r.shared_heat_forcing);
    nl();
}
static void dMA(const cl::MovingAvg7& m) {
    ln("MovingAvg7");
    for (double v : m.ta) hx(v);
    for (double v : m.tr) hx(v);
    std::printf(" %d %d", m.front, m.count); nl();
}

int main() {
    // ---- constants
    ln("const"); hx(cl::kSolarConstant); hx(cl::kStdPressureMb); hx(cl::MM_PER_INCH); nl();

    // ---- enum discriminants
    std::printf("enum TemperatureMethod %d %d\n", (int)ueb::TemperatureMethod::ModifiedForceRestore,
                (int)ueb::TemperatureMethod::Equilibrium);
    std::printf("enum SnowKind %d %d\n", (int)srf::SnowKind::DegreeDay, (int)srf::SnowKind::Ueb);
    std::printf("enum RadiationSourceMode %d %d %d\n", (int)srf::RadiationSource::Mode::Computed,
                (int)srf::RadiationSource::Mode::Timeseries, (int)srf::RadiationSource::Mode::File);
    std::printf("enum MeshSnowConfigMode %d %d %d\n", (int)srf::MeshSnowConfig::Mode::Auto,
                (int)srf::MeshSnowConfig::Mode::Yes, (int)srf::MeshSnowConfig::Mode::No);
    std::printf("enum RadiationOrigin %d %d %d %d\n", (int)cl::RadiationOrigin::Computed,
                (int)cl::RadiationOrigin::Timeseries, (int)cl::RadiationOrigin::MetFile,
                (int)cl::RadiationOrigin::LegacyHeat);
    std::printf("enum EvapMethod %d %d %d %d %d\n", (int)cl::EvapMethod::CONSTANT,
                (int)cl::EvapMethod::MONTHLY, (int)cl::EvapMethod::TIMESERIES,
                (int)cl::EvapMethod::TEMPERATURE, (int)cl::EvapMethod::PAN);
    std::printf("sizeof enums %zu %zu %zu %zu %zu %zu\n", sizeof(ueb::TemperatureMethod),
                sizeof(srf::SnowKind), sizeof(srf::RadiationSource::Mode),
                sizeof(srf::MeshSnowConfig::Mode), sizeof(cl::RadiationOrigin),
                sizeof(cl::EvapMethod));

    // ---- defaults, UEB point
    dParameters(ueb::Parameters{});
    dSite(ueb::Site{});
    {
        ueb::Forcing f{};
        ln("Forcing"); hx(f.seconds); hx(f.precipitation); hx(f.temperature); hx(f.wind);
        hx(f.humidity); hx(f.shortwave); hx(f.longwave); hx(f.cosine); hx(f.snow_catch);
        hx(f.snow_fraction); nl();
    }
    {
        ueb::ThermalSample t{};
        ln("ThermalSample"); hx(t.seconds); hx(t.surface); hx(t.bulk); nl();
    }
    dState(ueb::State{});
    dFluxes(ueb::Fluxes{});
    dLumpedState(ueb::LumpedState{});
    dLumpedFluxes(ueb::LumpedFluxes{});
    {
        ueb::HeldState h{};
        dLumpedState(h.accepted); dLumpedState(h.begin); dLumpedState(h.end);
        dLumpedFluxes(h.flux);
        ln("HeldState"); hx(h.duration); hx(h.elapsed); std::printf(" %d", (int)h.pending()); nl();
    }
    {
        ueb::NativeState n{};
        ln("NativeState");
        for (double v : n.forcing) hx(v);
        std::printf(" %zu %zu %zu %zu %zu %zu %d", n.columns.size(), n.last.size(), n.rain.size(),
                    n.snow.size(), n.sample_end.size(), n.sample_seconds.size(), (int)n.active);
        hx(n.sublimation); hx(n.precipitation); hx(n.liquid); hx(n.throughfall);
        hx(n.net_energy); hx(n.geometry_energy); hx(n.plowed_energy); hx(n.initial_energy);
        hx(n.initial_water); hx(n.removed); hx(n.transferred_water); hx(n.last_seconds);
        hx(n.completed_seconds); nl();
    }

    // ---- defaults, configuration
    {
        srf::SnowpackConfig s{};
        std::printf("SnowpackConfig %d", (int)s.model); hx(s.initial_energy); nl();
        dSite(s.site);
    }
    dRadSource(srf::RadiationSource{});
    {
        srf::ClimateSite c{};
        ln("ClimateSite"); hx(c.latitude); hx(c.longitude); hx(c.timezone); hx(c.elevation); nl();
    }
    {
        srf::SnowConfig c{};
        std::printf("SnowConfig %d %d %d", (int)c.model, (int)c.model_authored,
                    (int)c.division_temperature);
        hx(c.step); hx(c.bc_a); hx(c.bc_c);
        for (const auto& m : c.monthly_range) std::printf(" %d", (int)m.has_value());
        std::printf(" %d %d %zu", (int)c.site.has_value(), (int)c.met_file.empty(), c.packs.size());
        nl();
        dParameters(c.parameters);
        dRadSource(c.shortwave);
        dRadSource(c.longwave);
    }
    dMeshRow(srf::MeshSnowRow{});
    {
        srf::MeshSnowConfig m{};
        std::printf("MeshSnowConfig %d", (int)m.mode); hx(m.step);
        std::printf(" %zu %zu %d\n", m.defaults.size(), m.cells.size(), (int)m.assigned());
    }

    // ---- defaults, climate
    dRadState(cl::RadiationState{});
    {
        cl::SolarPosition s{};
        ln("SolarPosition"); hx(s.declination_deg); hx(s.eq_time_min); hx(s.hour_angle_deg);
        hx(s.zenith_deg); hx(s.cos_zenith); hx(s.elevation_deg); hx(s.eccentricity); nl();
    }
    dMA(cl::MovingAvg7{});
    {
        cl::ClimateState c{};
        dRadState(c.radiation);
        ln("ClimateState"); hx(c.temperature); hx(c.temp_range); hx(c.evap_rate);
        hx(c.wind_speed); hx(c.humidity); hx(c.temperature_src); hx(c.file_tmin);
        hx(c.file_tmax); hx(c.wind_speed_src); hx(c.gamma); hx(c.ea); hx(c.latitude);
        hx(c.elev); hx(c.snow_divt); nl();
        ln("monthly_evap"); for (double v : c.monthly_evap) hx(v); nl();
        ln("evaprate_ucf"); hx(c.evaprate_ucf); nl();
        ln("adjust_evap"); for (double v : c.adjust_evap) hx(v); nl();
        ln("next_evap"); hx(c.next_evap_date); hx(c.next_evap_rate);
        std::printf(" %d", c.evap_ts_pos); nl();
        ln("adjust_temp"); for (double v : c.adjust_temp) hx(v); nl();
        ln("adjust_rain"); for (double v : c.adjust_rain) hx(v); nl();
        ln("adjust_hydcon"); for (double v : c.adjust_hydcon) hx(v); nl();
        std::printf("evap_method %d", (int)c.evap_method); hx(c.infil_factor); nl();
        dMA(c.temp_ma);
        ln("recovery"); hx(c.recovery_factor);
        std::printf(" %d %d %d %d", c.temp_ts_index, c.evap_ts_index, c.humidity_ts_index,
                    c.recovery_pat_index);
        nl();
        ln("subdaily"); hx(c.tmin_daily); hx(c.tmax_daily); hx(c.prev_tmax); hx(c.hrsr);
        hx(c.hrss); hx(c.hrday); hx(c.dhrdy); hx(c.dydif); hx(c.dtlong);
        std::printf(" %d %d", (int)c.has_minmax, c.last_temp_doy); nl();
    }
    {
        cl::HourlyMetRecord r{};
        ln("HourlyMetRecord"); hx(r.date); for (double v : r.values) hx(v); nl();
        cl::SnowDailyWeather w{};
        ln("SnowDailyWeather"); hx(w.date); hx(w.tmin); hx(w.tmax); hx(w.previous_max);
        hx(w.wind); nl();
        cl::HourlyMeteorology h;
        std::printf("HourlyMeteorology %d %zu\n", (int)h.empty(), h.records().size());
    }

    // ---- HeldState::pending over a sweep (elapsed, duration)
    {
        const double inf = std::numeric_limits<double>::infinity();
        const double nan = std::numeric_limits<double>::quiet_NaN();
        const double vs[] = {0.0, -0.0, 1.0, -1.0, 0.5, 1e-300, 4.9406564584124654e-324,
                             1.7976931348623157e308, inf, -inf, nan, -nan, 3600.0,
                             3600.0000000000005, 3599.9999999999995, 1e-5};
        ln("pending");
        for (double e : vs)
            for (double d : vs) {
                ueb::HeldState h;
                h.elapsed = e;
                h.duration = d;
                std::printf("%d", (int)h.pending());
            }
        nl();
    }

    // ---- SnowConfig::modelFor
    {
        srf::SnowConfig c;
        srf::SnowpackConfig u;                 // default model Ueb
        srf::SnowpackConfig d; d.model = srf::SnowKind::DegreeDay;
        c.packs[0] = u; c.packs[2] = d; c.packs[5] = u; c.packs[-3] = d;
        ln("modelFor DD");
        for (int p = -4; p <= 7; ++p) std::printf(" %d", (int)c.modelFor(p));
        nl();
        c.model = srf::SnowKind::Ueb;
        ln("modelFor UEB");
        for (int p = -4; p <= 7; ++p) std::printf(" %d", (int)c.modelFor(p));
        nl();
        ln("packs order");
        for (const auto& kv : c.packs) std::printf(" %d:%d", kv.first, (int)kv.second.model);
        nl();
    }

    // ---- MeshSnowConfig::assigned and std::map key order
    {
        srf::MeshSnowConfig m;
        std::printf("assigned %d", (int)m.assigned());
        m.defaults["*"] = srf::MeshSnowRow{};
        std::printf(" %d", (int)m.assigned());
        m.defaults.clear();
        std::printf(" %d", (int)m.assigned());
        m.cells[3] = srf::MeshSnowRow{};
        std::printf(" %d", (int)m.assigned());
        m.mode = srf::MeshSnowConfig::Mode::No;
        std::printf(" %d\n", (int)m.assigned());
        const std::string keys[] = {"b", "a", "A", "*", "", std::string("\x80", 1),
                                    std::string("\xff", 1), "ab", std::string("a\0b", 3),
                                    "a", "Z", std::string("\x7f", 1), "~", "aa"};
        for (const auto& k : keys) m.defaults[k] = srf::MeshSnowRow{};
        ln("defaults order");
        for (const auto& kv : m.defaults) {
            std::printf(" [");
            for (unsigned char ch : kv.first) std::printf("%02x", ch);
            std::printf("]");
        }
        nl();
        const int cks[] = {7, -1, 0, 2147483647, -2147483647 - 1, 3, 12};
        for (int k : cks) m.cells[k] = srf::MeshSnowRow{};
        ln("cells order");
        for (const auto& kv : m.cells) std::printf(" %d", kv.first);
        nl();
    }

    // ---- data/HydrologyData.hpp
    {
        SnowpackStore s;
        std::printf("SnowpackStore %d", s.count());
        s.names.push_back("A"); s.names.push_back("B"); s.names.push_back("C");
        std::printf(" %d %zu %zu %zu %zu %zu %zu %zu\n", s.count(), s.report_values.size(),
                    s.daily_weather.size(), s.plowable.size(), s.impervious.size(),
                    s.pervious.size(), s.removal.size(), s.removal_subcatch.size());
        std::printf("SnowpackStore.config %d", (int)s.config.model); hx(s.config.step); nl();
        std::printf("SnowpackStore.mesh %d %d %d\n", (int)s.mesh_config.mode,
                    (int)s.mesh_config.assigned(), (int)s.meteorology.empty());
        AquiferStore a;
        std::printf("AquiferStore %d", a.count());
        for (int i = 0; i < 4; ++i) a.names.push_back("x");
        std::printf(" %d\n", a.count());
        LidControlStore l;
        std::printf("LidControlStore %d", l.count());
        for (int i = 0; i < 2; ++i) l.names.push_back("y");
        std::printf(" %d\n", l.count());
        LidUsageStore u;
        std::printf("LidUsageStore %d", u.count());
        for (int i = 0; i < 5; ++i) u.subcatch_index.push_back(i);
        std::printf(" %d\n", u.count());
        const int ns[] = {3, 0, 5, 2, 1, 6};
        for (int n : ns) {
            if (!u.wb_inflow.empty()) {
                u.wb_inflow[0] = 7.25; u.wb_evap[0] = -1.5; u.wb_final_vol.back() = 3.0;
            }
            u.resize_wb(n);
            std::printf("resize_wb %d %zu %zu %zu %zu %zu %zu %zu", n, u.wb_inflow.size(),
                        u.wb_evap.size(), u.wb_infil.size(), u.wb_surf_flow.size(),
                        u.wb_drain_flow.size(), u.wb_init_vol.size(), u.wb_final_vol.size());
            for (double v : u.wb_inflow) hx(v);
            for (double v : u.wb_evap) hx(v);
            for (double v : u.wb_final_vol) hx(v);
            nl();
        }
        std::printf("LidUsageStore other %zu %zu %zu\n", u.area.size(), u.rpt_file.size(),
                    u.from_perv.size());
    }
    return 0;
}
