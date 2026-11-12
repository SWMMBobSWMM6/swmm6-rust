// The 2D output-side differential, upstream side (phase 2a wave 5). Runs a deck through the 2D
// reference engine's C API as upstream's CLI does (open, initialize, start(save = 1), the step
// loop, end, report, close, destroy), writing the .rpt and .out to the paths given, and prints:
//   * at every report period, the 2D surface block SWMMEngine::fillSurfaceSnapshot fills into
//     the snapshot (this harness defines fillSurfaceSnapshot; the library's postOutputSnapshot
//     calls it through its PLT, link with -Wl,-E, and it calls the library's own through
//     dlsym(RTLD_NEXT) and prints what it filled): "SNAP k t <sim_time bits>", the counts, and
//     one line per array with its length and FNV-1a 64 hash of the values' bits (every value,
//     as hex bits, with SNAP2D_FULL=1);
//   * after end(), "END <rc>" and every member of ctx.mass_balance_2d (SurfaceRouter2D::finalize
//     publishes the solver statistics there), as hex bits;
//   * "DONE <exit code>" after report and close (the CLI's code).
// Nothing in a CLI run of this build reads the snapshot's surface block (its reader is the HDF5
// 2D writer, not in the reference build), so this is the only check of fillSurfaceSnapshot.
// The Rust twin is crates/engine/examples/snap2d.rs. Build: gate/build.sh snap2d.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <climits>
#include <clocale>
#include <cmath>
#include <compare>
#include <concepts>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cwchar>
#include <cwctype>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <ios>
#include <iosfwd>
#include <istream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <new>
#include <numbers>
#include <optional>
#include <ostream>
#include <queue>
#include <ratio>
#include <semaphore>
#include <sstream>
#include <stdexcept>
#include <stop_token>
#include <streambuf>
#include <string>
#include <string_view>
#include <system_error>
#include <thread>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>
#include <omp.h>
#include <dlfcn.h>
#define private public
#define protected public
#include "openswmm/engine/openswmm_engine.h"
#include "core/SWMMEngine.hpp"
#undef private
#undef protected

static int g_k = 0;
static bool g_full = false;

static std::uint64_t fnv(const void* p, std::size_t n, std::uint64_t h = 1469598103934665603ULL) {
    const auto* b = static_cast<const unsigned char*>(p);
    for (std::size_t i = 0; i < n; ++i) {
        h ^= b[i];
        h *= 1099511628211ULL;
    }
    return h;
}
static unsigned long long bits(double x) {
    return std::bit_cast<unsigned long long>(x);
}
static void vd(const char* name, const std::vector<double>& v) {
    std::printf("  %s n %zu h %016llx\n", name, v.size(),
                static_cast<unsigned long long>(fnv(v.data(), v.size() * sizeof(double))));
    if (g_full)
        for (std::size_t i = 0; i < v.size(); ++i)
            std::printf("    %s[%zu] %016llx\n", name, i, bits(v[i]));
}
static void vu8(const char* name, const std::vector<uint8_t>& v) {
    std::printf("  %s n %zu h %016llx\n", name, v.size(),
                static_cast<unsigned long long>(fnv(v.data(), v.size())));
}

void openswmm::SWMMEngine::fillSurfaceSnapshot(SimulationSnapshot& snap) const noexcept {
    using Real = void (*)(const openswmm::SWMMEngine*, SimulationSnapshot&);
    static Real real = reinterpret_cast<Real>(
        dlsym(RTLD_NEXT, "_ZNK8openswmm10SWMMEngine19fillSurfaceSnapshotERNS_18SimulationSnapshotE"));
    if (!real) {
        std::fprintf(stderr, "snap2d: the library's fillSurfaceSnapshot was not found\n");
        std::_Exit(3);
    }
    real(this, snap);
    ++g_k;
    std::printf("SNAP %d t %016llx\n", g_k, bits(snap.sim_time));
    std::printf("  counts %d %d %d %d species %d buildup %d gw %d\n", snap.surface_tri_count,
                snap.surface_vert_count, snap.surface_quad_count, snap.surface_edge_stride,
                snap.surface_species_count, snap.surface_buildup_count, snap.gw2d_active ? 1 : 0);
    vu8("cell_nv", snap.surface_cell_nv);
    vd("depth", snap.surface_depth);
    vd("head", snap.surface_head);
    vd("grad_hx", snap.surface_grad_hx);
    vd("grad_hy", snap.surface_grad_hy);
    vd("grad_hx_lim", snap.surface_grad_hx_lim);
    vd("grad_hy_lim", snap.surface_grad_hy_lim);
    for (std::size_t k = 0; k < snap.surface_snow.size(); ++k)
        std::printf("  snow[%zu] n %zu\n", k, snap.surface_snow[k].size());
    vd("rainfall", snap.surface_rainfall);
    vd("coupling_flux", snap.surface_coupling_flux);
    vd("net_source", snap.surface_net_source);
    vd("infil_rate", snap.surface_infil_rate);
    vd("infil_cum", snap.surface_infil_cum);
    vd("rain_cum", snap.surface_rain_cum);
    vd("species_conc", snap.surface_species_conc);
    vd("buildup", snap.surface_buildup);
    vd("edge_flux", snap.surface_edge_flux);
    vd("vert_head", snap.surface_vert_head);
    vd("vert_depth", snap.surface_vert_depth);
    vd("face_vx", snap.surface_face_vx);
    vd("face_vy", snap.surface_face_vy);
    vd("continuity_err", snap.surface_continuity_err);
    vd("stat_max_depth", snap.surface_stat_max_depth);
    vd("stat_max_velocity", snap.surface_stat_max_velocity);
    vd("stat_max_cont_err", snap.surface_stat_max_cont_err);
}

static void mb2(const openswmm::SimulationContext& ctx) {
    const auto& m = ctx.mass_balance_2d;
    const std::pair<const char*, double> d[] = {
        {"init_storage", m.init_storage}, {"final_storage", m.final_storage},
        {"external_in", m.external_in}, {"external_out", m.external_out},
        {"rainfall_in", m.rainfall_in}, {"coupling_1d_to_2d_in", m.coupling_1d_to_2d_in},
        {"coupling_2d_to_1d_out", m.coupling_2d_to_1d_out}, {"outfall_in", m.outfall_in},
        {"outfall_out", m.outfall_out}, {"boundary_in", m.boundary_in},
        {"boundary_out", m.boundary_out}, {"evap_out", m.evap_out}, {"infil_out", m.infil_out},
        {"infil_to_aquifer", m.infil_to_aquifer}, {"infil_aquifer_pending", m.infil_aquifer_pending},
        {"aquifer_in", m.aquifer_in}, {"solver_avg_h", m.solver_avg_h},
        {"solver_last_h", m.solver_last_h}, {"solver_active_min", m.solver_active_min},
        {"solver_active_mean", m.solver_active_mean}, {"solver_active_max", m.solver_active_max}};
    for (const auto& [k, v] : d)
        std::printf("  mb2.%s %016llx\n", k, bits(v));
    std::printf("  mb2.active %d\n", m.active ? 1 : 0);
    std::printf("  mb2.solver_nsteps %ld\n", static_cast<long>(m.solver_nsteps));
    std::printf("  mb2.solver_nrhs %ld\n", static_cast<long>(m.solver_nrhs));
    std::printf("  mb2.solver_n_tiers %d\n", m.solver_n_tiers);
    for (int k = 0; k < 8; ++k)
        std::printf("  mb2.solver_tier_cells[%d] %ld\n", k, static_cast<long>(m.solver_tier_cells[k]));
}

int main(int argc, char** argv) {
    if (argc < 4)
        return 2;
    g_full = std::getenv("SNAP2D_FULL") != nullptr;
    const char* inp = argv[1];
    const char* rpt = argv[2];
    const char* out = argv[3];
    SWMM_Engine engine = swmm_engine_create();
    int err = swmm_engine_open(engine, inp, rpt, out, nullptr);
    if (err != SWMM_OK) {
        swmm_engine_destroy(engine);
        std::printf("DONE %d\n", err);
        return 0;
    }
    err = swmm_engine_initialize(engine);
    if (err != SWMM_OK) {
        swmm_engine_close(engine);
        swmm_engine_destroy(engine);
        std::printf("DONE %d\n", err);
        return 0;
    }
    err = swmm_engine_start(engine, out[0] != '\0' ? 1 : 0);
    if (err != SWMM_OK) {
        swmm_engine_close(engine);
        swmm_engine_destroy(engine);
        std::printf("DONE %d\n", err);
        return 0;
    }
    double elapsed = 0.0;
    while (true) {
        err = swmm_engine_step(engine, &elapsed);
        if (err != SWMM_OK || elapsed <= 0.0)
            break;
    }
    int result = err;
    auto fin = [&](int status) {
        if (status != SWMM_OK && result == SWMM_OK)
            result = status;
    };
    const int rc_end = swmm_engine_end(engine);
    fin(rc_end);
    std::printf("END %d\n", rc_end);
    mb2(static_cast<openswmm::SWMMEngine*>(engine)->context());
    fin(swmm_engine_report(engine));
    fin(swmm_engine_close(engine));
    swmm_engine_destroy(engine);
    std::printf("DONE %d\n", result);
    return 0;
}
