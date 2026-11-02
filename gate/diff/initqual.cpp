// Differential for SWMMEngine::initQuality, initGeometry and initMassBalance (core/SWMMEngine.cpp):
// open a deck in the reference engine through its C API (the handle IS the SWMMEngine), call
// one init leg directly, then dump the engine's SimulationContext, the members open() sets and
// the quality members initQuality sets. The Rust twin is crates/engine/examples/initqual.rs.
// Compile with the engine's own layout defines (gate/build.sh builds it):
// -DNDEBUG -DSWMM_OS_LINUX -DSWMM_USE_OPENMP -fopenmp.
//
// Usage: initqual deck.inp quality|geometry|massref|massbal [side-file]
//   quality   open, initQuality()
//   geometry  open, initGeometry()
//   massref   (run first, in its own process; prints nothing) open, swmm_engine_initialize, then
//             write the hydrology solver state initMassBalance reads to the side file
//   massbal   open, install the side file's hydrology state, initQuality(), initGeometry(),
//             initMassBalance() (init_modules' order: initHydrology fills those solvers first)
// An init leg runs only when open returns 0.
//
// Why massbal is not initMassBalance() alone: initMassBalance reads four hydrology solvers that
// initHydrology fills (runoff_.soa().area[j] for every subcatchment, lid_.storedVolume(),
// snow_.state(), groundwater_.state()). Right after open they are empty, and upstream reads
// past the end of runoff_.soa().area (SIGSEGV on any deck with a subcatchment). Its pollutant
// term also needs initQuality's mass_balance.resize_quality. So a reference run (massref, a
// separate process so that an abort inside upstream's own initialize, e.g. U1, cannot take the
// comparison with it) opens the same deck, runs swmm_engine_initialize and writes the members
// initMassBalance reads, and the LID groups' members that initQuality's
// msxlid::initializeUnits and initMassBalance's lid_.storedVolume() read; massbal installs them
// in this engine (rebuilding the LID groups) and the translation installs them in its
// stand-ins. Every other input is the open + initQuality + initGeometry state on both sides.
// With no complete side file, or arrays too short for the reads, both print "massbal skipped".
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
#define private public
#define protected public
#include "openswmm/engine/openswmm_engine.h"
#include "core/SWMMEngine.hpp"
#define CTXDUMP_ENGINE
#define CTXDUMP_QUALITY
#include "ctxdump.hpp"

static void wvec(FILE* f, const char* tag, const std::vector<double>& v) {
    std::fprintf(f, "%s %zu", tag, v.size());
    for (double x : v) {
        uint64_t u;
        std::memcpy(&u, &x, 8);
        std::fprintf(f, " %016llx", (unsigned long long)u);
    }
    std::fprintf(f, "\n");
}

static void wint(FILE* f, const char* tag, int v) { std::fprintf(f, "%s %d\n", tag, v); }

static bool rvec(FILE* f, const char* tag, std::vector<double>& v) {
    char t[64];
    std::size_t n = 0;
    if (std::fscanf(f, "%63s %zu", t, &n) != 2 || std::strcmp(t, tag) != 0)
        return false;
    v.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        unsigned long long u;
        if (std::fscanf(f, "%llx", &u) != 1)
            return false;
        std::memcpy(&v[i], &u, 8);
    }
    return true;
}

static bool rdint(FILE* f, const char* tag, int& v) {
    char t[64];
    return std::fscanf(f, "%63s %d", t, &v) == 2 && std::strcmp(t, tag) == 0;
}

// massref: the reference engine's hydrology solver state after its own initialize().
static void massref(const char* deck, const char* side) {
    SWMM_Engine ha = swmm_engine_create();
    int rca = swmm_engine_open(ha, deck, "/dev/null", "/dev/null", nullptr);
    int rci = rca == 0 ? swmm_engine_initialize(ha) : -1;
    auto* a = static_cast<openswmm::SWMMEngine*>(ha);
    const auto& ctx = a->ctx_;
    const auto nsc = static_cast<std::size_t>(ctx.n_subcatches());
    const auto& ra = a->runoff_.soa_.area;
    const auto& gw = a->groundwater_.soa_;
    const auto& sn = a->snow_.soa_;
    bool any_snow = false;
    for (std::size_t i = 0; i < nsc && i < ctx.subcatches.snowpack.size(); ++i)
        any_snow = any_snow || ctx.subcatches.snowpack[i] >= 0;
    const bool ok = rca == 0 && ra.size() >= nsc && gw.total_depth.size() >= nsc &&
                    gw.lower_depth.size() >= nsc && gw.theta.size() >= nsc && gw.porosity.size() >= nsc &&
                    (!any_snow || (sn.enabled.size() >= nsc && sn.fw.size() >= sn.wsnow.size() &&
                                   sn.fArea.size() >= sn.wsnow.size()));
    FILE* f = std::fopen(side, "w");
    std::fprintf(f, "ok %d rc %d\n", ok ? 1 : 0, rci);
    if (ok) {
        const double lv = a->lid_.storedVolume();
        uint64_t u;
        std::memcpy(&u, &lv, 8);
        wvec(f, "runoff_area", ra);
        std::fprintf(f, "lid_stored_volume %016llx\n", (unsigned long long)u);
        std::vector<double> en(sn.enabled.begin(), sn.enabled.end());
        wvec(f, "snow_enabled", en);
        wvec(f, "snow_wsnow", sn.wsnow);
        wvec(f, "snow_fw", sn.fw);
        wvec(f, "snow_fArea", sn.fArea);
        wvec(f, "gw_total_depth", gw.total_depth);
        wvec(f, "gw_lower_depth", gw.lower_depth);
        wvec(f, "gw_theta", gw.theta);
        wvec(f, "gw_porosity", gw.porosity);
        // the LID groups: the members msxlid::initializeUnits and storedVolume() read
        wint(f, "lid_groups", a->lid_.numGroups());
        for (int t = 0; t < a->lid_.numGroups(); ++t) {
            const auto& g = a->lid_.group(t);
            wint(f, "lidg_type", static_cast<int>(g.type));
            wint(f, "lidg_count", g.count);
            std::vector<double> sci(g.subcatch_idx.begin(), g.subcatch_idx.end());
            std::vector<double> cti(g.control_idx.begin(), g.control_idx.end());
            std::vector<double> cov(g.stor_covered.begin(), g.stor_covered.end());
            wvec(f, "lidg_subcatch_idx", sci);
            wvec(f, "lidg_control_idx", cti);
            wvec(f, "lidg_stor_covered", cov);
            wvec(f, "lidg_area", g.area);
            wvec(f, "lidg_surf_store", g.surf_store);
            wvec(f, "lidg_soil_thick", g.soil_thick);
            wvec(f, "lidg_stor_void", g.stor_void);
            wvec(f, "lidg_pave_void", g.pave_void);
            wvec(f, "lidg_pave_imperv_frac", g.pave_imperv_frac);
            wvec(f, "lidg_surf_void_frac", g.surf_void_frac);
            wvec(f, "lidg_surf_side_slope", g.surf_side_slope);
            wvec(f, "lidg_unit_width", g.unit_width);
            wvec(f, "lidg_surf_depth", g.surf_depth);
            wvec(f, "lidg_soil_moist", g.soil_moist);
            wvec(f, "lidg_stor_depth", g.stor_depth);
            wvec(f, "lidg_pave_depth", g.pave_depth);
        }
        std::fprintf(f, "end\n");
    }
    std::fclose(f);
}

// massbal: install the side file's hydrology state into this engine; false when the file is
// missing, incomplete or not ok (rc is -2 for a missing or incomplete file).
static bool install(openswmm::SWMMEngine* eng, const char* side, int& rc) {
    rc = -2;
    FILE* f = std::fopen(side, "r");
    if (!f)
        return false;
    int ok = 0, r = 0;
    char t1[8], t2[8];
    bool good = std::fscanf(f, "%7s %d %7s %d", t1, &ok, t2, &r) == 4;
    if (!good) {
        std::fclose(f);
        return false;
    }
    if (!ok) {
        rc = r;
        std::fclose(f);
        return false;
    }
    std::vector<double> area, en, wsnow, fw, farea, td, ld, th, po, lidv;
    int ng = 0;
    char tag[64];
    unsigned long long u = 0;
    good = rvec(f, "runoff_area", area) && std::fscanf(f, "%63s %llx", tag, &u) == 2 &&
           rvec(f, "snow_enabled", en) && rvec(f, "snow_wsnow", wsnow) && rvec(f, "snow_fw", fw) &&
           rvec(f, "snow_fArea", farea) && rvec(f, "gw_total_depth", td) && rvec(f, "gw_lower_depth", ld) &&
           rvec(f, "gw_theta", th) && rvec(f, "gw_porosity", po) && rdint(f, "lid_groups", ng);
    std::vector<openswmm::lid::LIDGroupSoA> groups(good ? static_cast<std::size_t>(ng) : 0);
    for (int t = 0; good && t < ng; ++t) {
        auto& g = groups[static_cast<std::size_t>(t)];
        int type = 0;
        std::vector<double> sci, cti, cov;
        good = rdint(f, "lidg_type", type) && rdint(f, "lidg_count", g.count) &&
               rvec(f, "lidg_subcatch_idx", sci) && rvec(f, "lidg_control_idx", cti) &&
               rvec(f, "lidg_stor_covered", cov) &&
               rvec(f, "lidg_area", g.area) &&
               rvec(f, "lidg_surf_store", g.surf_store) &&
               rvec(f, "lidg_soil_thick", g.soil_thick) &&
               rvec(f, "lidg_stor_void", g.stor_void) &&
               rvec(f, "lidg_pave_void", g.pave_void) &&
               rvec(f, "lidg_pave_imperv_frac", g.pave_imperv_frac) &&
               rvec(f, "lidg_surf_void_frac", g.surf_void_frac) &&
               rvec(f, "lidg_surf_side_slope", g.surf_side_slope) &&
               rvec(f, "lidg_unit_width", g.unit_width) &&
               rvec(f, "lidg_surf_depth", g.surf_depth) &&
               rvec(f, "lidg_soil_moist", g.soil_moist) &&
               rvec(f, "lidg_stor_depth", g.stor_depth) &&
               rvec(f, "lidg_pave_depth", g.pave_depth);
        g.type = static_cast<openswmm::lid::LIDType>(type);
        g.subcatch_idx.assign(sci.begin(), sci.end());
        g.control_idx.assign(cti.begin(), cti.end());
        g.stor_covered.assign(cov.begin(), cov.end());
    }
    good = good && std::fscanf(f, "%63s", tag) == 1 && std::strcmp(tag, "end") == 0;
    std::fclose(f);
    if (!good)
        return false;
    rc = r;
    eng->runoff_.soa_.area = area;
    eng->lid_.groups_ = groups; // initializeUnits and storedVolume read groups_ only
    eng->snow_.soa_.enabled.assign(en.size(), 0);
    for (std::size_t i = 0; i < en.size(); ++i)
        eng->snow_.soa_.enabled[i] = static_cast<unsigned char>(en[i]);
    eng->snow_.soa_.wsnow = wsnow;
    eng->snow_.soa_.fw = fw;
    eng->snow_.soa_.fArea = farea;
    eng->groundwater_.soa_.total_depth = td;
    eng->groundwater_.soa_.lower_depth = ld;
    eng->groundwater_.soa_.theta = th;
    eng->groundwater_.soa_.porosity = po;
    return true;
}

int main(int argc, char** argv) {
    const char* mode = argc > 2 ? argv[2] : "quality";
    if (std::strcmp(mode, "massref") == 0) {
        massref(argv[1], argv[3]);
        return 0;
    }
    SWMM_Engine h = swmm_engine_create();
    int rc = swmm_engine_open(h, argv[1], "/dev/null", "/dev/null", nullptr);
    std::printf("open_rc %d\n", rc);
    auto* eng = static_cast<openswmm::SWMMEngine*>(h);
    if (rc == 0) {
        if (std::strcmp(mode, "quality") == 0) {
            eng->initQuality();
        } else if (std::strcmp(mode, "geometry") == 0) {
            eng->initGeometry();
        } else if (std::strcmp(mode, "massbal") == 0) {
            int rci = -2;
            const bool ok = install(eng, argv[3], rci);
            std::printf("massbal side rc %d\n", rci);
            if (ok) {
                eng->initQuality();
                eng->initGeometry();
                eng->initMassBalance();
            } else {
                std::printf("massbal skipped\n");
            }
        }
    }
    ctxdump::dump(eng->context());
    ctxdump::engdump(*eng);
    ctxdump::qualdump(*eng);
    return 0;
}
