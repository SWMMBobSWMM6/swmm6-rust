// The per-step state gate, upstream side: open a deck in the reference engine through its C
// API (the handle IS the SWMMEngine), then swmm_engine_initialize, swmm_engine_start(1) and
// swmm_engine_step as the CLI's run loop does, dumping the engine's whole state: the
// SimulationContext, the members open() sets (ctxdump, engdump), the router and hydraulic
// solvers (rtdump), the hydrology solvers (hydrodump), the quality members (qualdump) and the
// run lifecycle members (lifedump). The Rust twin is crates/engine/examples/ctxstep.rs.
// Compile with the engine's own layout defines (gate/build.sh builds it):
// -DNDEBUG -DSWMM_OS_LINUX -DSWMM_USE_OPENMP -fopenmp.
//
// Usage: ctxstep deck.inp init|start|step [N [hash|full K]]
//   init   open, initialize; print the return codes, then the full dump
//   start  open, initialize, start; the same
//   step   open, initialize, start (dumped in full, or hashed with "hash"), then up to N
//          steps (default 1). Each step prints "step k rc <rc> elapsed <bits> state <s>" and
//          then the full dump (no third argument), its FNV-1a 64 hash ("hash"), or the hash
//          except at step K ("full K"; K = 0 is the dump after start). The loop stops when a step fails or the run ends
//          (elapsed 0), as the CLI's does.
// A dump is only printed after a call returns; where upstream terminates (U1, U57) stdout
// is lost with the process, as before.
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
#include <codecvt>
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
#include <filesystem>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <istream>
#include <limits>
#include <locale>
#include <map>
#include <memory>
#include <mutex>
#include <new>
#include <numbers>
#include <numeric>
#include <optional>
#include <ostream>
#include <queue>
#include <ratio>
#include <semaphore>
#include <set>
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
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>
#include <omp.h>
#define private public
#define protected public
#include "openswmm/engine/openswmm_engine.h"
#include "core/SWMMEngine.hpp"
#include "hydraulics/fv/ExplicitFvSolver.hpp"
#define CTXDUMP_ENGINE
#define CTXDUMP_QUALITY
#include "ctxdump.hpp"
#include "rtdump.hpp"
#include "hydrodump.hpp"
#include "lifedump.hpp"
#include <fcntl.h>
#include <unistd.h>

static void full_dump(openswmm::SWMMEngine* eng) {
    ctxdump::dump(eng->context());
    ctxdump::engdump(*eng);
    ctxdump::rtdump(*eng);
    ctxdump::hydrodump(*eng);
    ctxdump::qualdump(*eng);
    ctxdump::lifedump(*eng);
}

// the dump's FNV-1a 64 hash: the dumpers print to stdout, so fd 1 goes to a scratch file
// while they run
static unsigned long long hashed_dump(openswmm::SWMMEngine* eng) {
    std::fflush(stdout);
    char path[] = "/tmp/ctxstepXXXXXX";
    const int fd = mkstemp(path);
    const int saved = dup(1);
    dup2(fd, 1);
    full_dump(eng);
    std::fflush(stdout);
    dup2(saved, 1);
    close(saved);
    unsigned long long h = 1469598103934665603ULL;
    lseek(fd, 0, SEEK_SET);
    unsigned char buf[65536];
    ssize_t n;
    while ((n = read(fd, buf, sizeof buf)) > 0)
        for (ssize_t i = 0; i < n; ++i) {
            h ^= buf[i];
            h *= 1099511628211ULL;
        }
    close(fd);
    unlink(path);
    return h;
}

static void dump_or_hash(openswmm::SWMMEngine* eng, bool hash, const char* what) {
    if (hash)
        std::printf("hash %s %016llx\n", what, hashed_dump(eng));
    else
        full_dump(eng);
}

int main(int argc, char** argv) {
    const char* mode = argc > 2 ? argv[2] : "init";
    const int nsteps = argc > 3 ? std::atoi(argv[3]) : 1;
    const char* how = argc > 4 ? argv[4] : "";
    const int full_at = (std::strcmp(how, "full") == 0 && argc > 5) ? std::atoi(argv[5]) : -1;
    const bool hashing = std::strcmp(how, "hash") == 0 || full_at >= 0;
    SWMM_Engine h = swmm_engine_create();
    int rc = swmm_engine_open(h, argv[1], "/dev/null", "/dev/null", nullptr);
    std::printf("open_rc %d\n", rc);
    auto* eng = static_cast<openswmm::SWMMEngine*>(h);
    if (rc == 0) {
        rc = swmm_engine_initialize(h);
        std::printf("init_rc %d\n", rc);
    }
    if (std::strcmp(mode, "init") == 0) {
        full_dump(eng);
        return 0;
    }
    if (rc == 0) {
        rc = swmm_engine_start(h, 1);
        std::printf("start_rc %d\n", rc);
    }
    if (std::strcmp(mode, "start") == 0) {
        full_dump(eng);
        return 0;
    }
    dump_or_hash(eng, hashing && full_at != 0, "start");
    for (int k = 1; rc == 0 && k <= nsteps; ++k) {
        double elapsed = -1.0;
        rc = swmm_engine_step(h, &elapsed);
        unsigned long long u;
        std::memcpy(&u, &elapsed, 8);
        std::printf("step %d rc %d elapsed %016llx state %d\n", k, rc, u, (int)eng->context().state);
        char what[32];
        std::snprintf(what, sizeof what, "%d", k);
        dump_or_hash(eng, hashing && k != full_at, what);
        if (elapsed == 0.0)
            break;
    }
    return 0;
}
