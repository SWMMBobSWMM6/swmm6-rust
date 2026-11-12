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
//   helper K NAME SIDE   (the step helpers' differential) open, initialize, start, run K
//          real steps, then write the context (SIDE.ctx, ctxdump::dump) and the engine
//          members the helpers read (SIDE.eng) and call one of step()'s helpers directly on
//          that state, printing its result and the context dump after it. The Rust twin
//          loads SIDE.ctx and SIDE.eng into its own engine (after its own initialize and
//          start, which the init and start gates hold identical) and calls its translation.
//          NAME "all" runs every helper on the same state (each in a forked child) and
//          prints "helper NAME" and the hash of the context after it, instead of the dump.
//          NAME: stats (updateStatistics(dt)), massbal (updateRoutingMassBalance(dt)),
//          finalstore (computeFinalStorage), finalqual (computeFinalQualityMassBalance),
//          nodetotals (accumulateNodeRoutingTotals(dt/2)), diverged (checkRoutingDiverged),
//          events (isBetweenEvents(current_date), isInSteadyState(0)), lateral
//          (assembleLateralInflows(dt); decks with street inlets are skipped, as computeAll
//          belongs to the routing step), forcings (resetStepMassBalance, applyForcings(dt)).
//          dt is the last step's (elapsed_ms - old_elapsed_ms) / 1000, or ROUTING_STEP at K 0.
//   runoff K M SIDE [full J]   (the runoff step's differential) open, initialize, start, run
//          K real steps, then (K > 0) write the context (SIDE.ctx), the hydrology solvers
//          (SIDE.hyd, hydrodump) and the engine members the runoff step reads (SIDE.eng);
//          then M times: save_lat_qual_state, resetStepMassBalance, applyForcings(dt) and
//          stepRunoff(dt) as step() calls them, printing "runoff m" and the hash of the runoff
//          dump after the first and the last call and after every call that took a runoff step (runoff_dump: the context, hydrodump and the runoff engine members;
//          the full dump at m = J), then advance the routing clock by hand
//          (old_elapsed_ms = elapsed_ms; elapsed_ms += 1000 dt, capped at the duration),
//          stopping once elapsed_ms reaches the duration. dt is the last step's, or
//          ROUTING_STEP at K 0. An exception stepRunoff throws prints "runoff m threw <what>"
//          and ends the loop. The Rust twin loads SIDE.* (K > 0) into its own started
//          engine and runs its translation the same way.
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
#ifdef OPENSWMM_HAS_2D
// the 2D build (gate/build.sh ctxstep2d): the 2D router's whole state too
#include "twoddump.hpp"
#endif
#include <fcntl.h>
#include <unistd.h>
#include <sys/wait.h>

static void full_dump(openswmm::SWMMEngine* eng) {
    // The writer thread updates the output plugin's period count (lifedump prints it) after
    // step() returns; drain it first, so the dump does not race it (the translation writes
    // in the step).
    eng->io_thread_.flush();
    ctxdump::dump(eng->context());
    ctxdump::engdump(*eng);
    ctxdump::rtdump(*eng);
    ctxdump::hydrodump(*eng);
    ctxdump::qualdump(*eng);
    ctxdump::lifedump(*eng);
#ifdef OPENSWMM_HAS_2D
    ctxdump::twoddump(*eng);
#endif
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

static unsigned long long hashed_ctx(openswmm::SWMMEngine* eng) {
    std::fflush(stdout);
    char path[] = "/tmp/ctxstepXXXXXX";
    const int fd = mkstemp(path);
    const int saved = dup(1);
    dup2(fd, 1);
    ctxdump::dump(eng->context());
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

static void wvec(FILE* f, const char* tag, const std::vector<double>& v) {
    std::fprintf(f, "%s %zu", tag, v.size());
    for (double x : v) {
        uint64_t u;
        std::memcpy(&u, &x, 8);
        std::fprintf(f, " %016llx", (unsigned long long)u);
    }
    std::fprintf(f, "\n");
}

static void wnum(FILE* f, const char* tag, double x) {
    uint64_t u;
    std::memcpy(&u, &x, 8);
    std::fprintf(f, "%s %016llx\n", tag, (unsigned long long)u);
}

// one helper on the engine's current state; 1 when the deck is skipped for it
static int run_helper(openswmm::SWMMEngine* eng, const std::string& n, double dt) {
    auto& ctx = eng->ctx_;
    if (n == "stats") {
        eng->updateStatistics(dt);
    } else if (n == "massbal") {
        eng->updateRoutingMassBalance(dt);
    } else if (n == "finalstore") {
        eng->computeFinalStorage();
    } else if (n == "finalqual") {
        eng->computeFinalQualityMassBalance();
    } else if (n == "nodetotals") {
        eng->accumulateNodeRoutingTotals(dt / 2.0);
    } else if (n == "diverged") {
        if (eng->router_.fvSolver()) {
            std::printf("diverged skipped (FV)\n");
            return 1;
        }
        std::printf("diverged rc %d\n", eng->checkRoutingDiverged());
    } else if (n == "events") {
        std::printf("between %d steady %d\n", eng->isBetweenEvents(ctx.current_date) ? 1 : 0,
                    eng->isInSteadyState(0) ? 1 : 0);
    } else if (n == "lateral") {
        if (eng->inlet_.soa_.count > 0) {
            std::printf("lateral skipped (inlets)\n");
            return 1;
        }
        eng->assembleLateralInflows(dt);
    } else if (n == "forcings") {
        eng->resetStepMassBalance();
        eng->applyForcings(dt);
    }
    return 0;
}

static const char* const kHelpers[] = {"stats", "massbal", "finalstore", "finalqual", "nodetotals",
                                       "diverged", "events", "lateral", "forcings"};

static int helper_mode(SWMM_Engine h, openswmm::SWMMEngine* eng, int k_steps, const char* name, const char* side) {
    int rc = 0;
    int k = 0;
    for (; rc == 0 && k < k_steps; ++k) {
        double elapsed = -1.0;
        rc = swmm_engine_step(h, &elapsed);
        if (elapsed == 0.0)
            break;
    }
    auto& ctx = eng->ctx_;
    std::printf("helper_at %d rc %d\n", k, rc);
    // the side files: the context, then the engine members the helpers read
    {
        std::string p = std::string(side) + ".ctx";
        std::fflush(stdout);
        FILE* f = std::fopen(p.c_str(), "w");
        const int saved = dup(1);
        dup2(fileno(f), 1);
        ctxdump::dump(ctx);
        std::fflush(stdout);
        dup2(saved, 1);
        close(saved);
        std::fclose(f);
    }
    {
        std::string p = std::string(side) + ".eng";
        FILE* f = std::fopen(p.c_str(), "w");
        std::fprintf(f, "helper_at %d rc %d\n", k, rc);
        std::fprintf(f, "next_event_ %d\n", eng->next_event_);
        std::fprintf(f, "fv_accepted %d\n", eng->quality_.fvTransportAccepted() ? 1 : 0);
        const auto& r = eng->runoff_.soa_;
        wvec(f, "area", r.area);
        wvec(f, "imperv_pct", r.imperv_pct);
        wvec(f, "frac_imperv0", r.frac_imperv0);
        wvec(f, "frac_imperv1", r.frac_imperv1);
        wvec(f, "depth_imperv0", r.depth_imperv0);
        wvec(f, "depth_imperv1", r.depth_imperv1);
        wvec(f, "depth_perv", r.depth_perv);
        const auto& sn = eng->snow_.soa_;
        std::vector<double> en(sn.enabled.begin(), sn.enabled.end());
        wvec(f, "snow_enabled", en);
        wvec(f, "snow_wsnow", sn.wsnow);
        wvec(f, "snow_fw", sn.fw);
        wvec(f, "snow_fArea", sn.fArea);
        wnum(f, "snow_removed", sn.removed);
        wvec(f, "buildup", eng->surface_quality_.buildup);
        wvec(f, "lid_drain_q_interp_", eng->lid_drain_q_interp_);
        wvec(f, "wet_q_interp_", eng->wet_q_interp_);
        wvec(f, "gw_q_interp_", eng->gw_q_interp_);
        std::vector<double> gqn(eng->gw_q_node_.begin(), eng->gw_q_node_.end());
        wvec(f, "gw_q_node_", gqn);
        std::vector<double> un, uq;
        for (auto& p : eng->lid_drain_q_units_) {
            un.push_back(p.first);
            uq.push_back(p.second);
        }
        wvec(f, "lid_units_node", un);
        wvec(f, "lid_units_q", uq);
        // the LID units' moisture state (lid_.storedVolume reads it; the groups' other
        // members are initHydrology's, the same on both sides)
        std::fprintf(f, "lid_groups %d\n", eng->lid_.numGroups());
        for (int t = 0; t < eng->lid_.numGroups(); ++t) {
            const auto& g = eng->lid_.group(t);
            std::fprintf(f, "lidg_count %d\n", g.count);
            wvec(f, "lidg_surf_depth", g.surf_depth);
            wvec(f, "lidg_soil_moist", g.soil_moist);
            wvec(f, "lidg_stor_depth", g.stor_depth);
            wvec(f, "lidg_pave_depth", g.pave_depth);
        }
        std::fprintf(f, "end\n");
        std::fclose(f);
    }
    const double dt = k > 0 ? (ctx.elapsed_ms - ctx.old_elapsed_ms) / 1000.0 : ctx.options.routing_step;
    if (std::strcmp(name, "all") != 0) {
        if (run_helper(eng, name, dt) == 0)
            ctxdump::dump(ctx);
        return 0;
    }
    // "all": every helper on the same state, each in a forked child; the dumps hashed
    for (const char* n : kHelpers) {
        std::printf("helper %s\n", n);
        std::fflush(stdout);
        const pid_t pid = fork();
        if (pid == 0) {
            if (run_helper(eng, n, dt) == 0)
                std::printf("hash %s %016llx\n", n, hashed_ctx(eng));
            std::fflush(stdout);
            _exit(0);
        }
        int st = 0;
        waitpid(pid, &st, 0);
    }
    return 0;
}

// the runoff step's dump: the context, the hydrology solvers and the engine members the
// runoff step writes (the router, quality-routing and lifecycle state are not loaded on the
// Rust side and the runoff step does not touch them, so they are left out)
static void runoff_dump(openswmm::SWMMEngine* eng) {
    ctxdump::dump(eng->ctx_);
    ctxdump::hydrodump(*eng);
    ctxdump::pf("old_runoff_time_", eng->old_runoff_time_);
    ctxdump::pf("new_runoff_time_", eng->new_runoff_time_);
    ctxdump::pf("old_runoff_ms_", eng->old_runoff_ms_);
    ctxdump::pf("prev_runoff_step_sec_", eng->prev_runoff_step_sec_);
    ctxdump::pf("new_runoff_ms_", eng->new_runoff_ms_);
    ctxdump::pi("has_runoff_", eng->has_runoff_ ? 1 : 0);
    ctxdump::pi("has_snow_", eng->has_snow_ ? 1 : 0);
    ctxdump::vf("wet_q_interp_", eng->wet_q_interp_);
    ctxdump::vf("lid_drain_q_interp_", eng->lid_drain_q_interp_);
    for (auto& q : eng->lid_drain_q_units_) {
        ctxdump::pi("lid_drain_q_units_.first", q.first);
        ctxdump::pf("lid_drain_q_units_.second", q.second);
    }
    ctxdump::vf("lid_vlidin_vol_", eng->lid_vlidin_vol_);
    ctxdump::vf("lid_qsurf_cfs_", eng->lid_qsurf_cfs_);
    ctxdump::vf("surface_quality_.buildup", eng->surface_quality_.buildup);
}

static unsigned long long hashed_runoff(openswmm::SWMMEngine* eng) {
    std::fflush(stdout);
    char path[] = "/tmp/ctxstepXXXXXX";
    const int fd = mkstemp(path);
    const int saved = dup(1);
    dup2(fd, 1);
    runoff_dump(eng);
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

// the side file of the runoff step: the engine members stepRunoff reads that neither the
// context dump nor hydrodump covers (the lifecycle's runoff clocks and work vectors, the
// surface quality buildup)
static void write_runoff_side(openswmm::SWMMEngine* eng, const char* side, int k, int rc) {
    {
        std::string p = std::string(side) + ".ctx";
        std::fflush(stdout);
        FILE* f = std::fopen(p.c_str(), "w");
        const int saved = dup(1);
        dup2(fileno(f), 1);
        ctxdump::dump(eng->ctx_);
        std::fflush(stdout);
        dup2(saved, 1);
        close(saved);
        std::fclose(f);
    }
    {
        std::string p = std::string(side) + ".hyd";
        std::fflush(stdout);
        FILE* f = std::fopen(p.c_str(), "w");
        const int saved = dup(1);
        dup2(fileno(f), 1);
        ctxdump::hydrodump(*eng);
        std::fflush(stdout);
        dup2(saved, 1);
        close(saved);
        std::fclose(f);
    }
    std::string p = std::string(side) + ".eng";
    FILE* f = std::fopen(p.c_str(), "w");
    std::fprintf(f, "runoff_at %d rc %d\n", k, rc);
    wnum(f, "old_runoff_time_", eng->old_runoff_time_);
    wnum(f, "new_runoff_time_", eng->new_runoff_time_);
    wnum(f, "old_runoff_ms_", eng->old_runoff_ms_);
    wnum(f, "prev_runoff_step_sec_", eng->prev_runoff_step_sec_);
    wnum(f, "new_runoff_ms_", eng->new_runoff_ms_);
    std::fprintf(f, "has_runoff_ %d\n", eng->has_runoff_ ? 1 : 0);
    std::fprintf(f, "has_snow_ %d\n", eng->has_snow_ ? 1 : 0);
    wvec(f, "wet_q_interp_", eng->wet_q_interp_);
    wvec(f, "lid_drain_q_interp_", eng->lid_drain_q_interp_);
    std::vector<double> un, uq;
    for (auto& q : eng->lid_drain_q_units_) {
        un.push_back(q.first);
        uq.push_back(q.second);
    }
    wvec(f, "lid_units_node", un);
    wvec(f, "lid_units_q", uq);
    wvec(f, "lid_vlidin_vol_", eng->lid_vlidin_vol_);
    wvec(f, "lid_qsurf_cfs_", eng->lid_qsurf_cfs_);
    wvec(f, "buildup", eng->surface_quality_.buildup);
    std::fprintf(f, "end\n");
    std::fclose(f);
}

static int runoff_mode(SWMM_Engine h, openswmm::SWMMEngine* eng, int k_steps, int m_calls, const char* side,
                       int full_at) {
    int rc = 0;
    int k = 0;
    for (; rc == 0 && k < k_steps; ++k) {
        double elapsed = -1.0;
        rc = swmm_engine_step(h, &elapsed);
        if (elapsed == 0.0)
            break;
    }
    auto& ctx = eng->ctx_;
    std::printf("runoff_at %d rc %d\n", k, rc);
    if (k_steps > 0)
        write_runoff_side(eng, side, k, rc);
    if (rc != 0)
        return 0;
    const double dt = k > 0 ? (ctx.elapsed_ms - ctx.old_elapsed_ms) / 1000.0 : ctx.options.routing_step;
    const double total = ctx.options.totalDurationMs();
    double last_new_ms = eng->new_runoff_ms_;
    for (int m = 1; m <= m_calls && ctx.elapsed_ms < total; ++m) {
        ctx.save_lat_qual_state();
        eng->resetStepMassBalance();
        eng->applyForcings(dt);
        try {
            eng->stepRunoff(dt);
        } catch (const std::exception& e) {
            std::printf("runoff %d threw %s\n", m, e.what());
            break;
        }
        // dumped after the first and last call and after every call that took a runoff step
        if (m == 1 || m == m_calls || m == full_at || eng->new_runoff_ms_ != last_new_ms) {
            std::printf("runoff %d\n", m);
            if (m == full_at)
                runoff_dump(eng);
            else
                std::printf("hash %d %016llx\n", m, hashed_runoff(eng));
        }
        last_new_ms = eng->new_runoff_ms_;
        ctx.old_elapsed_ms = ctx.elapsed_ms;
        ctx.elapsed_ms += 1000.0 * dt;
        if (ctx.elapsed_ms > total)
            ctx.elapsed_ms = total;
    }
    return 0;
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
    if (std::strcmp(mode, "helper") == 0)
        return rc == 0 ? helper_mode(h, eng, nsteps, argv[4], argv[5]) : 0;
    if (std::strcmp(mode, "runoff") == 0) {
        const int full_j = (argc > 7 && std::strcmp(argv[6], "full") == 0) ? std::atoi(argv[7]) : -1;
        return rc == 0 ? runoff_mode(h, eng, nsteps, std::atoi(argv[4]), argv[5], full_j) : 0;
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
