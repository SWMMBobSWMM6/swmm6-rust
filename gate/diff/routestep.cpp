// The routing-step cross test, upstream side. Opens a deck in the reference engine through
// its C API, initializes and starts it, then steps it with swmm_engine_step. This harness
// defines SWMMEngine::stepRouting itself: the library calls stepRouting through its PLT, so
// step() reaches this definition, which dumps the whole state, calls the library's own
// stepRouting (found with dlsym(RTLD_NEXT)) with the dt step() chose, and dumps again. So
// every dumped pair is exactly what upstream's step() hands its routing step and what it
// gets back.
//
// Output, for each stepRouting call c with K0 <= c <= K1:
//   BEFORE c <dt as hex bits>
//   <ctxdump::dump> <rtdump> <qualdump> <stepdump>
//   AFTER c
//   <the same dumps>            (or "throw <what>" when stepRouting throws)
// and "END c" when the run ends before K1. The run stops after call K1.
// The Rust twin is crates/engine/examples/routestep.rs. Compile with -Wl,-E (gate/build.sh)
// so the library's PLT call binds to this definition.
// Usage: routestep deck.inp K0 K1
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
#define CTXDUMP_QUALITY
#include "rtdump.hpp"

static int g_call = 0, g_k0 = 1, g_k1 = 1;

// Engine members stepRouting reads or writes that the other dumpers do not print. The Rust
// twin is stepdump in examples/routestep.rs.
static void stepdump(const openswmm::SWMMEngine& e) {
    ctxdump::pi("engine.last_step_steady_", e.last_step_steady_);
    ctxdump::pf("engine.new_rule_time_ms_", e.new_rule_time_ms_);
    ctxdump::pi("engine.checkpoint_enabled_", e.checkpoint_enabled_);
    ctxdump::pi("engine.ard_init_attempted_", e.ard_init_attempted_);
    ctxdump::pi("engine.rdii_iface_file_.open", e.rdii_iface_file_.isOpen());
    ctxdump::pi("engine.rdii_iface_file_.writing", e.rdii_iface_file_.isWriting());
    // read by assembleLateralInflows (filled by stepRunoff / the groundwater step)
    ctxdump::vf("engine.wet_q_interp_", e.wet_q_interp_);
    ctxdump::pi("engine.lid_drain_q_units_.n", (long long)e.lid_drain_q_units_.size());
    for (const auto& u : e.lid_drain_q_units_) {
        ctxdump::pi("engine.lid_drain_q_units_[].first", u.first);
        ctxdump::pf("engine.lid_drain_q_units_[].second", u.second);
    }
    ctxdump::vf("engine.gw_q_interp_", e.gw_q_interp_);
    {
        std::vector<long long> t(e.gw_q_node_.begin(), e.gw_q_node_.end());
        ctxdump::vi("engine.gw_q_node_", t);
    }
}

static void dumpall(openswmm::SWMMEngine& e) {
    ctxdump::dump(e.context());
    ctxdump::rtdump(e);
    ctxdump::qualdump(e);
    stepdump(e);
}

void openswmm::SWMMEngine::stepRouting(double dt_routing) {
    using Real = void (*)(openswmm::SWMMEngine*, double);
    static Real real = reinterpret_cast<Real>(dlsym(RTLD_NEXT, "_ZN8openswmm10SWMMEngine11stepRoutingEd"));
    if (!real) {
        std::fprintf(stderr, "routestep: the library's stepRouting was not found\n");
        std::_Exit(3);
    }
    ++g_call;
    const bool d = g_call >= g_k0 && g_call <= g_k1;
    if (d) {
        uint64_t u;
        std::memcpy(&u, &dt_routing, 8);
        std::printf("BEFORE %d %016llx\n", g_call, (unsigned long long)u);
        dumpall(*this);
    }
    try {
        real(this, dt_routing);
    } catch (const std::exception& ex) {
        if (d) {
            std::printf("AFTER %d\nthrow %s\n", g_call, ex.what());
        }
        if (g_call >= g_k1) {
            std::fflush(stdout);
            std::_Exit(0);
        }
        throw;
    }
    if (d) {
        std::printf("AFTER %d\n", g_call);
        dumpall(*this);
    }
    if (g_call >= g_k1) {
        std::fflush(stdout);
        std::_Exit(0);
    }
}

int main(int argc, char** argv) {
    if (argc < 4)
        return 2;
    g_k0 = std::atoi(argv[2]);
    g_k1 = std::atoi(argv[3]);
    SWMM_Engine e = swmm_engine_create();
    int rc = swmm_engine_open(e, argv[1], "/dev/null", "/dev/null", nullptr);
    if (rc == 0)
        rc = swmm_engine_initialize(e);
    if (rc == 0)
        rc = swmm_engine_start(e, 1);
    if (rc != 0) {
        std::printf("NOSTART %d\n", rc);
        return 0;
    }
    double elapsed = 0.0;
    while ((rc = swmm_engine_step(e, &elapsed)) == 0 && elapsed > 0.0) {
    }
    std::printf("END %d rc %d\n", g_call, rc);
    std::fflush(stdout);
    std::_Exit(0);
}
