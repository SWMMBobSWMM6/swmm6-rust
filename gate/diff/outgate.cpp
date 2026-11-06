// The output-side gate, upstream side (gate (a): the writers without the runoff step).
// Runs a deck through the reference engine's C API exactly as upstream's CLI does
// (src/cli/main.cpp: open, initialize, start(save = 1), the step loop, end, report, close,
// destroy), writing the .rpt and .out to the paths given.
//
//   outgate deck.inp rpt out start
//       stop after start(): on success flush every stream and _Exit (the .rpt holds the
//       preamble, the .out the header, as after start()); on a failure take the CLI's own
//       failure path (close, destroy). Prints "rc <code>" and exits with it.
//   outgate deck.inp rpt out run K0 K1 [STRIDE]   (STRIDE -n: about n periods, spread evenly)
//       the whole run. This harness defines SWMMEngine::postOutputSnapshot (the library's
//       step() calls it through its PLT; link with -Wl,-E): for each call that writes a
//       report period (output due and the report date not before REPORT_START), counted
//       from 1, with K0 <= p <= K1 and (p - K0) % STRIDE == 0, it prints "PERIOD p" and the
//       whole state before calling the library's own postOutputSnapshot. After the step
//       loop it prints "ENDSTATE <step rc>" and the state end() starts from, then runs end,
//       report and close and prints "DONE <exit code>".
// The state is ctxdump::dump, rtdump, qualdump and outside(), the engine members the output
// side reads that the other dumpers do not print. The Rust twin is
// crates/engine/examples/outgate.rs.
// Build: gate/build.sh outgate
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
#include "hydraulics/fv/ExplicitFvSolver.hpp"
#include "plugins/DefaultOutputPlugin.hpp"
#include "hydraulics/TimestepController.hpp"
#include "core/DateTime.hpp"
#define CTXDUMP_QUALITY
#include "rtdump.hpp"

static int g_period = 0, g_k0 = 1, g_k1 = 0, g_stride = 1;
static bool g_run = false;

// The engine members the output side reads beyond the context, rtdump and qualdump. The
// Rust twin is load_outside in examples/outgate.rs.
static void outside(openswmm::SWMMEngine& e) {
    using namespace ctxdump;
    pf("engine.old_runoff_ms_", e.old_runoff_ms_);
    pf("engine.new_runoff_ms_", e.new_runoff_ms_);
    vf("engine.old_gw_rate_", e.old_gw_rate_);
    pi("engine.save_results_", e.save_results_);
    pi("engine.do_routing_", e.do_routing_);
    {
        std::vector<long long> t(e.hotstart_save_done_.begin(), e.hotstart_save_done_.end());
        vi("engine.hotstart_save_done_", t);
    }
    {
        auto f = [](const std::vector<float>& v) { return std::vector<double>(v.begin(), v.end()); };
        const auto& a = e.avg_;
        vf("engine.avg_.node_depth", f(a.node_depth));
        vf("engine.avg_.node_head", f(a.node_head));
        vf("engine.avg_.node_volume", f(a.node_volume));
        vf("engine.avg_.node_lat_inflow", f(a.node_lat_inflow));
        vf("engine.avg_.node_total_inflow", f(a.node_total_inflow));
        vf("engine.avg_.node_overflow", f(a.node_overflow));
        vf("engine.avg_.link_flow", f(a.link_flow));
        vf("engine.avg_.link_depth", f(a.link_depth));
        vf("engine.avg_.link_velocity", f(a.link_velocity));
        vf("engine.avg_.link_volume", f(a.link_volume));
        vf("engine.avg_.link_capacity", f(a.link_capacity));
        pi("engine.avg_.n_steps", a.n_steps);
    }
    vf("runoff_.soa_.area", e.runoff_.soa_.area);
    {
        const auto& g = e.groundwater_.state();
        vf("groundwater_.gw_flow", g.gw_flow);
        vf("groundwater_.bottom_elev", g.bottom_elev);
        vf("groundwater_.lower_depth", g.lower_depth);
        vf("groundwater_.theta", g.theta);
        vf("groundwater_.upper_evap", g.upper_evap);
        vf("groundwater_.lower_evap", g.lower_evap);
        vf("groundwater_.total_depth", g.total_depth);
        vf("groundwater_.porosity", g.porosity);
    }
    {
        const auto& sn = e.snow_.state();
        pi("snow_.n_subcatch", sn.n_subcatch);
        std::vector<long long> ex(sn.external.begin(), sn.external.end());
        vi("snow_.external", ex);
        vf("snow_.external_net", sn.external_net);
        vf("snow_.imelt", sn.imelt);
        vf("snow_.fArea", sn.fArea);
        std::vector<long long> en(sn.enabled.begin(), sn.enabled.end());
        vi("snow_.enabled", en);
    }
    pi("lid_.groups.n", e.lid_.numGroups());
    for (int t = 0; t < e.lid_.numGroups(); ++t) {
        const auto& g = e.lid_.group(t);
        pi("lid_.group.count", g.count);
        std::vector<long long> s(g.subcatch_idx.begin(), g.subcatch_idx.end());
        vi("lid_.group.subcatch_idx", s);
        std::vector<long long> c(g.control_idx.begin(), g.control_idx.end());
        vi("lid_.group.control_idx", c);
        vf("lid_.group.area", g.area);
        vf("lid_.group.evap_loss", g.evap_loss);
        vf("lid_.group.infil_loss", g.infil_loss);
        vf("lid_.group.wb_inflow", g.wb_inflow);
        vf("lid_.group.wb_evap", g.wb_evap);
        vf("lid_.group.wb_infil", g.wb_infil);
        vf("lid_.group.wb_surf_flow", g.wb_surf_flow);
        vf("lid_.group.wb_drain_flow", g.wb_drain_flow);
        vf("lid_.group.wb_init_vol", g.wb_init_vol);
        vf("lid_.group.wb_final_vol", g.wb_final_vol);
    }
    int np = -1;
    for (auto* p : e.plugins_.output_plugins_)
        if (auto* d = dynamic_cast<openswmm::DefaultOutputPlugin*>(p))
            np = d->n_periods_;
    pi("engine.out_plugin.n_periods_", np);
}

static void dumpall(openswmm::SWMMEngine& e) {
    ctxdump::dump(e.context());
    ctxdump::rtdump(e);
    ctxdump::qualdump(e);
    outside(e);
}

void openswmm::SWMMEngine::postOutputSnapshot(double dt_step) noexcept {
    using Real = void (*)(openswmm::SWMMEngine*, double);
    static Real real = reinterpret_cast<Real>(dlsym(RTLD_NEXT, "_ZN8openswmm10SWMMEngine18postOutputSnapshotEd"));
    if (!real) {
        std::fprintf(stderr, "outgate: the library's postOutputSnapshot was not found\n");
        std::_Exit(3);
    }
    if (g_run && hydraulics::TimestepController::output_due(ctx_)) {
        const double report_date =
            datetime::addSeconds(ctx_.options.start_date, (ctx_.next_report_ms + 1.0) / 1000.0);
        if (!(report_date < ctx_.options.report_start)) {
            ++g_period;
            if (g_period >= g_k0 && (g_k1 <= 0 || g_period <= g_k1) &&
                (g_period - g_k0) % g_stride == 0) {
                std::printf("PERIOD %d\n", g_period);
                dumpall(*this);
            }
        }
    }
    real(this, dt_step);
}

int main(int argc, char** argv) {
    if (argc < 5)
        return 2;
    const char* inp = argv[1];
    const char* rpt = argv[2];
    const char* out = argv[3];
    g_run = std::strcmp(argv[4], "run") == 0;
    if (g_run) {
        g_k0 = argc > 5 ? std::atoi(argv[5]) : 1;
        g_k1 = argc > 6 ? std::atoi(argv[6]) : 0;
        g_stride = argc > 7 ? std::atoi(argv[7]) : 1;
    }
    SWMM_Engine engine = swmm_engine_create();
    int err = swmm_engine_open(engine, inp, rpt, out, nullptr);
    if (err != SWMM_OK) {
        swmm_engine_destroy(engine);
        std::printf("rc %d\n", err);
        return err;
    }
    err = swmm_engine_initialize(engine);
    if (err != SWMM_OK) {
        swmm_engine_close(engine);
        swmm_engine_destroy(engine);
        std::printf("rc %d\n", err);
        return err;
    }
    err = swmm_engine_start(engine, out[0] != '\0' ? 1 : 0);
    if (err != SWMM_OK) {
        swmm_engine_close(engine);
        swmm_engine_destroy(engine);
        std::printf("rc %d\n", err);
        return err;
    }
    if (g_run && g_stride <= 0) {
        // STRIDE -n: about n periods, evenly spread over the run (the first one included)
        const int target = g_stride < 0 ? -g_stride : 100;
        const auto& o = static_cast<openswmm::SWMMEngine*>(engine)->context().options;
        const double n_est = (o.end_date - o.start_date) * 86400.0 / std::max(o.report_step, 1.0);
        g_stride = std::max(1, static_cast<int>(std::ceil(n_est / target)));
    }
    if (!g_run) {
        std::printf("rc 0\n");
        std::fflush(nullptr);
        std::_Exit(0);
    }
    double elapsed = 0.0;
    while (true) {
        err = swmm_engine_step(engine, &elapsed);
        if (err != SWMM_OK || elapsed <= 0.0)
            break;
    }
    std::printf("ENDSTATE %d\n", err);
    dumpall(*static_cast<openswmm::SWMMEngine*>(engine));
    std::fflush(stdout);
    int result = err;
    auto fin = [&](int status) {
        if (status != SWMM_OK && result == SWMM_OK)
            result = status;
    };
    fin(swmm_engine_end(engine));
    fin(swmm_engine_report(engine));
    fin(swmm_engine_close(engine));
    swmm_engine_destroy(engine);
    std::printf("DONE %d\n", result);
    std::fflush(stdout);
    return result;
}
