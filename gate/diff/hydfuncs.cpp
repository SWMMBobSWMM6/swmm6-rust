// Differential harness for hydraulics/Transect.cpp, Street.cpp, ForceMain.cpp,
// HydClosureKernels.hpp, Link.cpp, Culvert.cpp and math/FindRoot.cpp. Upstream's compiled
// engine library supplies every out-of-line function (it exports all C++ symbols); the inline
// ones (hydkernels::*, link::filledCircularOffsetBump) compile into this harness with the
// engine's flags.
//
//   hydfuncs gen FILE   write the inputs to FILE: every [TRANSECTS] row, [STREETS] row and
//                       SHAPE curve of every corpus and upstream test deck (read with
//                       upstream's own section handlers), synthetic transects / streets /
//                       shape curves (valid, degenerate, invalid, non-finite), setParams cases
//                       for every closed and open shape, and the sweep lists (depths, flows,
//                       slopes, Reynolds numbers, culvert heads, root-finder cases).
//   hydfuncs FILE       run every command of FILE in order and print the results as hex bits
//                       (ints in decimal). The Rust example crates/engine/examples/hydfuncs.rs
//                       runs the same commands and must print the same bytes.
//
// Build (from the repository root):
//   U=/root/work/swmm6src/upstream/src/engine; B=/root/work/swmm6src/build
//   g++ -std=c++20 -O3 -ffp-contract=off -mavx2 -mfma -DNDEBUG -DSWMM_OS_LINUX -I$U \
//       -I$B/include gate/diff/hydfuncs.cpp $B/src/engine/libopenswmm.engine.so \
//       -Wl,-rpath,$B/src/engine -o gate/diff/hydfuncs
#include "hydraulics/Transect.hpp"
#include "hydraulics/Street.hpp"
#include "hydraulics/ForceMain.hpp"
#include "hydraulics/HydClosureKernels.hpp"
#include "hydraulics/Link.hpp"
#include "hydraulics/Culvert.hpp"
#include "hydraulics/XSectBatch.hpp"
#include "math/FindRoot.hpp"
#include "data/InfraData.hpp"
#include "data/LinkData.hpp"
#include "data/TableData.hpp"
#include "core/SimulationContext.hpp"
#include "core/UnitConversion.hpp"
#include "input/InputReader.hpp"
#include "input/SectionRegistry.hpp"
#include "input/handlers/InfraHandler.hpp"
#include "input/handlers/LinksHandler.hpp"
#include "input/handlers/TablesHandler.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace openswmm;
using transect::TransectData;
static const int NT = transect::N_TRANSECT_TBL;

// ---------------------------------------------------------------- output
static std::string OUT;
static void hx(double v) {
    uint64_t u;
    std::memcpy(&u, &v, 8);
    char b[24];
    std::snprintf(b, sizeof b, "%016llx ", (unsigned long long)u);
    OUT += b;
}
static void hi(long long v) {
    OUT += std::to_string(v);
    OUT += ' ';
}
static void hs(const char* s) {
    OUT += s;
    OUT += ' ';
}
static void nl() {
    OUT += '\n';
    if (OUT.size() > (1u << 22)) {
        std::fwrite(OUT.data(), 1, OUT.size(), stdout);
        OUT.clear();
    }
}
static void harr(const double* a, size_t n) {
    for (size_t i = 0; i < n; ++i)
        hx(a[i]);
}
static void hname(const std::string& s) {
    if (s.empty()) {
        OUT += "~ ";
        return;
    }
    static const char* H = "0123456789abcdef";
    for (unsigned char c : s) {
        OUT += H[c >> 4];
        OUT += H[c & 15];
    }
    OUT += ' ';
}

// ---------------------------------------------------------------- gen helpers
static std::mt19937_64 rnd(20261102);
static double u01() { return (double)(rnd() >> 11) * 0x1.0p-53; }
static int ri(int lo, int hi) { return lo + (int)(rnd() % (uint64_t)(hi - lo + 1)); }
static void wd(FILE* f, double v) {
    uint64_t u;
    std::memcpy(&u, &v, 8);
    std::fprintf(f, " %016llx", (unsigned long long)u);
}
static void wname(FILE* f, const std::string& s) {
    if (s.empty()) {
        std::fprintf(f, " ~");
        return;
    }
    std::fprintf(f, " ");
    for (unsigned char c : s)
        std::fprintf(f, "%02x", c);
}
static void wlist(FILE* f, const char* name, const std::vector<double>& v) {
    std::fprintf(f, "L %s %zu", name, v.size());
    for (double x : v)
        wd(f, x);
    std::fprintf(f, "\n");
}

struct TRow {
    std::string name;
    double nl, nr, nc, xlb, xrb, xf, yf, lf;
    std::vector<double> st, el;
};
static void wtr(FILE* f, const TRow& r) {
    std::fprintf(f, "TR");
    wname(f, r.name);
    for (double v : {r.nl, r.nr, r.nc, r.xlb, r.xrb, r.xf, r.yf, r.lf})
        wd(f, v);
    std::fprintf(f, " %zu", r.st.size());
    for (double v : r.st)
        wd(f, v);
    std::fprintf(f, " %zu", r.el.size());
    for (double v : r.el)
        wd(f, v);
    std::fprintf(f, "\n");
}
static void wte(FILE* f, const std::vector<double>& ucf) {
    std::fprintf(f, "TE %zu", ucf.size());
    for (double v : ucf)
        wd(f, v);
    std::fprintf(f, "\n");
}
static void wst(FILE* f, const street::StreetParams& s) {
    std::fprintf(f, "ST");
    for (double v : {s.width, s.curb_height, s.slope, s.roughness, s.gutter_depression,
                     s.gutter_width})
        wd(f, v);
    std::fprintf(f, " %d", s.sides);
    for (double v : {s.back_width, s.back_slope, s.back_roughness})
        wd(f, v);
    std::fprintf(f, "\n");
}
static void wcu(FILE* f, const std::string& name, double yfull, double wmax, int n,
                const std::vector<double>& x, const std::vector<double>& y) {
    std::fprintf(f, "CU");
    wname(f, name);
    wd(f, yfull);
    wd(f, wmax);
    std::fprintf(f, " %d %zu", n, x.size());
    for (double v : x)
        wd(f, v);
    std::fprintf(f, " %zu", y.size());
    for (double v : y)
        wd(f, v);
    std::fprintf(f, "\n");
}

static int gen(const char* path) {
    FILE* f = std::fopen(path, "w");
    if (!f)
        return 1;
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    const double INF = std::numeric_limits<double>::infinity();
    const double UCF_US = ucf::Ucf[ucf::LENGTH][0], UCF_SI = ucf::Ucf[ucf::LENGTH][1];

    // ---- sweep lists
    wlist(f, "LD", {0.0, -0.0, 1e-6, 0.005, 0.01, std::nextafter(0.01, 1.0), 0.02, 0.05, 0.1,
                    0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.999, 1.0, 1.0001, 1.2, 2.0,
                    3.0, -0.5, NaN, INF});
    wlist(f, "LQV", {0.0, 1.0, -2.5, 30.0, NaN});
    wlist(f, "LV", {0.0, 2.3, -4.0, NaN});
    wlist(f, "LN", {0.013, 0.03, 0.0, -0.01});
    wlist(f, "LS", {0.001, -0.02, 0.0, 0.5});
    wlist(f, "LQF", {0.0, 1e-6, 0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0, 1.05, 1.2, 2.0, -1.0,
                     NaN, INF});
    {
        std::vector<double> cd;
        for (int k = 0; k <= 60; ++k)
            cd.push_back(k * 0.05);
        for (double v : {-0.1, 1e-9, 0.001, 0.94, 0.95, std::nextafter(0.95, 0.0),
                         std::nextafter(0.95, 1.0), 0.96, 0.999, 1.6, 1.7, 1.8, 2.5, 10.0, 100.0,
                         NaN, INF})
            cd.push_back(v);
        wlist(f, "CD", cd);
    }
    wlist(f, "CS", {0.0, 0.01, -0.02});
    // force main
    wlist(f, "FV", {0.0, -0.0, 1e-9, 1e-4, 0.01, 0.5, 1.0, 2.5, 5.0, 10.0, 30.0, -3.0, NaN, INF,
                    -INF});
    wlist(f, "FR", {0.0, -1.0, 1e-6, 0.01, 0.125, 0.25, 0.5, 1.0, 2.5, NaN, INF});
    wlist(f, "FP", {0.0, -10.0, 1e-6, 0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 100.0, 120.0,
                    130.0, 150.0, NaN, INF});
    wlist(f, "FE", {0.0, 1e-6, 0.0005, 0.005, 0.05, -0.001, NaN, INF});
    wlist(f, "FH", {0.0, -1.0, 0.01, 0.25, 1.0, NaN, INF});
    wlist(f, "FRE", {-5.0, 0.0, 1.0, 9.99, 10.0, 500.0, 1999.9, 2000.0, 2000.0000001, 2500.0,
                     3000.0, 3999.999, 4000.0, 4000.0001, 1e5, 1e8, 9.99e9, 1e10, 1e12, NaN,
                     INF});
    wlist(f, "FB", {0.0, 1e-5, 0.0005, 0.01, 100.0, 120.0, 140.0, -5.0, NaN, INF});
    wlist(f, "FY", {0.0, -1.0, 0.5, 1.0, 3.0, 10.0, NaN, INF});
    wlist(f, "FS", {0.0, -0.01, 1e-5, 0.001, 0.02, 1.0, NaN});
    wlist(f, "FL", {1.0, 0.0, -2.0, 0.5, 1.75, 10.0, NaN, INF});
    // culvertInflow kernel
    {
        std::vector<double> kh;
        for (int k = 0; k <= 40; ++k)
            kh.push_back(k * 0.1);
        for (double v : {-0.3, 1e-9, 0.95 * 3.0, 2.84999, 2.85001, 100.0, NaN, INF})
            kh.push_back(v);
        wlist(f, "KH", kh);
    }
    wlist(f, "KS", {0.0, 0.01, -0.05, 0.3});

    // ---- every [TRANSECTS], [STREETS] row and SHAPE curve in the corpus
    std::vector<std::string> decks;
    for (const char* root : {"/root/work/corpus", "/root/work/swmm6src/upstream/tests"}) {
        for (auto& e : std::filesystem::recursive_directory_iterator(root)) {
            if (e.is_regular_file() && e.path().extension() == ".inp")
                decks.push_back(e.path().string());
        }
    }
    std::sort(decks.begin(), decks.end());
    long nTr = 0, nSt = 0, nCu = 0;
    for (const auto& d : decks) {
        input::SectionRegistry reg;
        reg.register_builtin("TRANSECTS", input::handle_transects);
        reg.register_builtin("STREETS", input::handle_streets);
        reg.register_builtin("CURVES", input::handle_curves);
        input::InputReader rd(reg);
        SimulationContext ctx;
        try {
            rd.read(d, ctx);
        } catch (...) {
            continue;
        }
        const auto& ts = ctx.transects;
        if (ts.count() > 0) {
            std::fprintf(f, "TS\n");
            for (int t = 0; t < ts.count(); ++t) {
                auto u = (size_t)t;
                TRow r{ts.names[u],         ts.n_left[u],       ts.n_right[u],
                       ts.n_channel[u],     ts.x_left_bank[u],  ts.x_right_bank[u],
                       ts.x_factor[u],      ts.y_factor[u],     ts.length_factor[u],
                       ts.stations[u],      ts.elevations[u]};
                wtr(f, r);
                ++nTr;
            }
            wte(f, {UCF_US, UCF_SI});
        }
        const auto& ss = ctx.streets;
        for (int s = 0; s < ss.count(); ++s) {
            auto u = (size_t)s;
            for (double ucf_len : {UCF_US, UCF_SI}) {
                street::StreetParams sp;
                sp.width = ss.t_crown[u] / ucf_len;
                sp.curb_height = ss.h_curb[u] / ucf_len;
                sp.slope = ss.sx[u] / 100.0;
                sp.roughness = ss.n_road[u];
                sp.gutter_depression = ss.gutter_depres[u] / ucf_len;
                sp.gutter_width = ss.gutter_width[u] / ucf_len;
                sp.sides = ss.sides[u];
                sp.back_width = ss.back_width[u] / ucf_len;
                sp.back_slope = ss.back_slope[u] / 100.0;
                sp.back_roughness = ss.back_n[u];
                wst(f, sp);
                ++nSt;
            }
        }
        for (const auto& tb : ctx.tables.tables) {
            if (tb.type != TableType::CURVE_SHAPE)
                continue;
            for (double yf : {1.0, 3.5, 0.75}) {
                double wn = 0.0;
                for (double v : tb.y)
                    if (v > wn)
                        wn = v;
                wcu(f, tb.id, yf, wn * yf, (int)tb.x.size(), tb.x, tb.y);
                ++nCu;
            }
        }
    }
    std::fprintf(stderr, "corpus: %zu decks, %ld transect rows, %ld street cases, %ld shape cases\n",
                 decks.size(), nTr, nSt, nCu);

    // ---- synthetic transects
    auto randT = [&](int kind) {
        TRow r;
        r.name = "SYN" + std::to_string(kind);
        int n = ri(2, 30);
        double x = u01() * 50.0 - 10.0;
        double base = (rnd() % 3 == 0) ? 800.0 : 0.0;
        for (int i = 0; i < n; ++i) {
            if (i > 0 && rnd() % 6 != 0)
                x += 0.5 + u01() * 20.0; // else a vertical wall (repeated station)
            r.st.push_back(x);
            double e = base + u01() * 12.0;
            if (rnd() % 7 == 0 && i > 0)
                e = r.el.back(); // flat run
            r.el.push_back(e);
        }
        r.nc = 0.01 + u01() * 0.05;
        auto pickN = [&]() {
            switch (rnd() % 5) {
            case 0: return 0.0;
            case 1: return -0.02;
            case 2: return r.nc;
            default: return 0.02 + u01() * 0.1;
            }
        };
        r.nl = pickN();
        r.nr = pickN();
        auto pickBank = [&](bool left) {
            switch (rnd() % 4) {
            case 0: return r.st[(size_t)ri(0, n - 1)];
            case 1: return left ? r.st.front() : r.st.back();
            case 2: return u01() * 60.0 - 10.0;
            default: return left ? r.st[(size_t)std::min(1, n - 1)] : r.st[(size_t)std::max(0, n - 2)];
            }
        };
        r.xlb = pickBank(true);
        r.xrb = pickBank(false);
        const double xfs[] = {0.0, 1.0, 2.5, 0.5};
        r.xf = xfs[rnd() % 4];
        const double yfs[] = {0.0, 799.0, -3.0, 2.0};
        r.yf = yfs[rnd() % 4];
        const double lfs[] = {1.0, 0.0, -1.0, 1.3, 0.8};
        r.lf = lfs[rnd() % 5];
        return r;
    };
    for (int g = 0; g < 40; ++g) {
        std::fprintf(f, "TS\n");
        int rows = ri(1, 6);
        for (int k = 0; k < rows; ++k)
            wtr(f, randT(g * 10 + k));
        if (g % 4 == 0)
            wte(f, {UCF_US, UCF_SI, 0.0, NaN, 3.0});
        else
            wte(f, {UCF_US, UCF_SI});
    }
    {   // degenerate and invalid rows
        std::fprintf(f, "TS\n");
        auto T = [&](std::string nm, double nl, double nr, double nc, double xlb, double xrb,
                     double xf, double yf, double lf, std::vector<double> st,
                     std::vector<double> el) {
            wtr(f, TRow{nm, nl, nr, nc, xlb, xrb, xf, yf, lf, std::move(st), std::move(el)});
        };
        T("empty", 0.02, 0.02, 0.03, 0, 0, 1, 0, 1, {}, {});
        T("one", 0.02, 0.02, 0.03, 0, 0, 1, 0, 1, {1.0}, {2.0});
        T("two", 0.02, 0.02, 0.03, 0, 10, 1, 0, 1, {0.0, 10.0}, {5.0, 0.0});
        T("mismatch", 0.02, 0.02, 0.03, 0, 10, 1, 0, 1, {0.0, 5.0, 10.0}, {5.0, 0.0});
        T("flat", 0.02, 0.02, 0.03, 0, 10, 1, 0, 1, {0.0, 5.0, 10.0}, {1.0, 1.0, 1.0});
        T("n0", 0.02, 0.02, 0.0, 0, 10, 1, 0, 1, {0.0, 5.0, 10.0}, {5.0, 0.0, 5.0});
        T("nneg", 0.02, 0.02, -0.03, 0, 10, 1, 0, 1, {0.0, 5.0, 10.0}, {5.0, 0.0, 5.0});
        T("nnan", 0.02, 0.02, NaN, 0, 10, 1, 0, 1, {0.0, 5.0, 10.0}, {5.0, 0.0, 5.0});
        T("vnotch", 0.0, 0.0, 0.03, 0, 10, 1, 0, 1, {0.0, 5.0, 10.0}, {5.0, 0.0, 5.0});
        T("decr", 0.05, 0.06, 0.03, 10, 0, 1, 0, 1, {10.0, 5.0, 0.0}, {5.0, 0.0, 5.0});
        T("elnan", 0.05, 0.06, 0.03, 0, 10, 1, 0, 1, {0.0, 5.0, 10.0}, {5.0, NaN, 5.0});
        T("elnan0", 0.05, 0.06, 0.03, 0, 10, 1, 0, 1, {0.0, 5.0, 10.0}, {NaN, 0.0, 5.0});
        T("stnan", 0.05, 0.06, 0.03, 0, 10, 1, 0, 1, {0.0, NaN, 10.0}, {5.0, 0.0, 5.0});
        T("elinf", 0.05, 0.06, 0.03, 0, 10, 1, 0, 1, {0.0, 5.0, 10.0}, {INF, 0.0, 5.0});
        T("lfnan", 0.05, 0.06, 0.03, 0, 10, 1, 0, NaN, {0.0, 5.0, 10.0}, {5.0, 0.0, 5.0});
        T("banksnan", 0.05, 0.06, 0.03, NaN, NaN, 1, 0, 1, {0.0, 5.0, 10.0}, {5.0, 0.0, 5.0});
        T("walls", 0.05, 0.08, 0.03, 2, 8, 1, 0, 1, {0, 0, 2, 2, 5, 8, 8, 10, 10},
          {9, 6, 6, 3, 0, 3, 6, 6, 9});
        T("spike", 0.05, 0.08, 0.03, 2, 8, 1, 0, 1, {0, 2, 4, 5, 6, 8, 10}, {5, 1, 0, 7, 0, 1, 5});
        T("tiny", 0.05, 0.08, 0.03, 0, 1e-9, 1, 0, 1, {0, 1e-9, 2e-9}, {1e-12, 0, 1e-12});
        T("huge", 0.05, 0.08, 0.03, 0, 1e9, 1, 1e6, 1, {0, 1e9, 2e9}, {1e7, 0, 1e7});
        T(std::string("nul\0name", 8), 0.05, 0.08, 0.03, 0, 10, 1, 0, 1, {0.0, 5.0, 10.0}, {5.0, 0.0, 5.0});
        T(std::string(600, 'L'), 0.05, 0.08, 0.03, 0, 10, 1, 0, 1, {0.0, 5.0, 10.0},
          {5.0, 0.0, 5.0});
        wte(f, {UCF_US, UCF_SI, 0.0, NaN});
    }
    {   // uniquely named rows whose SWMM_TRACE_TRANSECT dump holds inf, -nan and subnormals
        std::fprintf(f, "TS\n");
        const double NNaN = -NaN;
        wtr(f, TRow{"dumpinf", 0.05, 0.06, 0.03, 0, 10, 1, 0, 1, {0.0, 5.0, 10.0}, {INF, 0.0, 5.0}});
        wtr(f, TRow{"dumpnegnan", 0.05, 0.06, 0.03, 0, 10, 1, 0, 1, {0.0, 5.0, 10.0}, {NNaN, 0.0, 5.0}});
        wtr(f, TRow{"dumpsub", 0.05, 0.06, 0.03, 0, 2e-160, 1, 0, 1, {0.0, 1e-160, 2e-160},
                    {1e-160, 0.0, 1e-160}});
        wte(f, {1.0});
    }

    // ---- synthetic streets
    auto S = [&](double w, double ch, double sl, double n, double gd, double gw, int sides,
                 double bw, double bs, double bn) {
        street::StreetParams sp;
        sp.width = w;
        sp.curb_height = ch;
        sp.slope = sl;
        sp.roughness = n;
        sp.gutter_depression = gd;
        sp.gutter_width = gw;
        sp.sides = sides;
        sp.back_width = bw;
        sp.back_slope = bs;
        sp.back_roughness = bn;
        wst(f, sp);
    };
    for (int sides : {1, 2, 0, 3}) {
        S(20, 0.5, 0.02, 0.016, 0.167, 2, sides, 0, 0, 0);
        S(20, 0.5, 0.02, 0.016, 0.167, 2, sides, 10, 0.04, 0.02);
        S(20, 0.5, 0.02, 0.016, 0, 0, sides, 10, 0.04, 0);
        S(20, 0.5, 0.0, 0.016, 0, 0, sides, 0, 0, 0);
        S(0, 0, 0, 0.016, 0, 0, sides, 0, 0, 0);
        S(20, 0.5, 0.02, 0.016, 0.167, 2, sides, 10, 0.5, 0.02); // backing above crown
        S(20, 0.5, 0.2, 0.016, 0.167, 2, sides, 1, 0.01, 0.02);  // crown above backing
        S(16, 0.5, 0.02, 0.0, 0.167, 2, sides, 5, 0.02, 0.02);   // zero road n
        S(16, 0.5, 0.02, 0.016, 0.167, 20, sides, 5, 0.02, 0.02); // gutter wider than road
        S(NaN, 0.5, 0.02, 0.016, 0.167, 2, sides, 5, 0.02, 0.02);
    }
    for (int k = 0; k < 60; ++k)
        S(5 + u01() * 30, u01(), u01() * 0.08, 0.01 + u01() * 0.03, (rnd() % 2) * u01() * 0.3,
          (rnd() % 2) * u01() * 3, 1 + (int)(rnd() % 2), (rnd() % 2) * u01() * 20, u01() * 0.1,
          (rnd() % 3) * 0.02);

    // ---- synthetic shape curves
    for (double yf : {1.0, 3.0, 0.0, -1.0, NaN, 1e-3, 1e3}) {
        std::vector<double> x = {0.0, 0.2, 0.5, 0.8, 1.0}, y = {0.4, 0.9, 1.0, 0.8, 0.0};
        wcu(f, "sc1", yf, yf, 5, x, y);
    }
    auto C = [&](const char* nm, std::vector<double> x, std::vector<double> y, int n = -1) {
        for (double yf : {1.0, 2.5})
            wcu(f, nm, yf, 0.0, n < 0 ? (int)x.size() : n, x, y);
    };
    C("nozero", {0.1, 0.5, 1.0}, {0.5, 1.0, 0.5});
    C("second_above1", {0.0, 1.5}, {0.5, 1.0});
    C("second_neg", {0.0, -0.1}, {1.0, 1.0});
    C("above1", {0.0, 0.5, 1.2, 1.5}, {0.2, 1.0, 0.7, 0.1});
    C("decrx", {0.0, 0.5, 0.4, 1.0}, {0.2, 1.0, 0.7, 0.1});
    C("negw", {0.0, 0.5, 0.7, 1.0}, {0.2, 1.0, -0.7, 0.1});
    C("negw2", {0.0, 0.5}, {0.2, -1.0});
    C("first1", {1.0, 1.0}, {1.0, 1.0});
    C("firstneg", {-0.1, 1.0}, {1.0, 1.0});
    C("firstnegw", {0.0, 1.0}, {-1.0, 1.0});
    C("one0", {0.0}, {1.0});
    C("one5", {0.5}, {1.0});
    C("none", {0.5, 1.0}, {1.0, 1.0}, 0);
    C("short", {0.0, 0.3, 0.6, 1.0}, {0.5, 1.0, 1.0, 0.5}, 2);
    C("dupx", {0.0, 0.3, 0.3, 0.7, 0.7, 1.0}, {0.5, 0.6, 1.0, 1.0, 0.4, 0.2});
    C("allzero", {0.0, 0.5, 1.0}, {0.0, 0.0, 0.0});
    C("rect", {0.0, 1.0}, {1.0, 1.0});
    C("tri", {0.0, 1.0}, {0.0, 1.0});
    C("short1", {0.0, 0.5}, {1.0, 1.0});
    C("nanx", {0.0, NaN, 1.0}, {1.0, 1.0, 1.0});
    C("nany", {0.0, 0.5, 1.0}, {1.0, NaN, 1.0});
    C("fine", {0.0, 0.019999, 0.02, 0.020001, 0.98, 0.99999, 1.0}, {0.1, 0.5, 0.6, 0.7, 0.9, 1.0, 0.3});
    for (int k = 0; k < 40; ++k) {
        int n = ri(2, 25);
        std::vector<double> x, y;
        double xx = (rnd() % 2) ? 0.0 : u01() * 0.3;
        for (int i = 0; i < n; ++i) {
            x.push_back(xx);
            y.push_back(u01() * 2.0);
            xx += (rnd() % 8 == 0) ? 0.0 : u01() * 2.0 / n;
        }
        wcu(f, "rc", 0.5 + u01() * 5, 1.0, n, x, y);
    }

    // ---- setParams cases for the link and culvert sweeps: type p0 p1 p2 p3 ucf full
    auto SH = [&](int t, double a, double b, double c, double d, double ucf_, int full) {
        std::fprintf(f, "SH %d", t);
        for (double v : {a, b, c, d, ucf_})
            wd(f, v);
        std::fprintf(f, " %d\n", full);
    };
    SH(0, 0, 0, 0, 0, 1.0, 0);
    SH(1, 3, 0, 0, 0, 1.0, 1);
    SH(1, 0.25, 0, 0, 0, 1.0, 0);
    SH(1, 0.9, 0, 0, 0, UCF_SI, 1);
    SH(2, 3, 0.5, 0, 0, 1.0, 0);
    SH(3, 4, 6, 0, 0, 1.0, 1);
    SH(3, 1, 10, 0, 0, 1.0, 0);
    SH(4, 4, 6, 0, 0, 1.0, 1);
    SH(4, 4, 6, 2, 0, 1.0, 0);
    SH(5, 5, 10, 2, 3, 1.0, 1);
    SH(6, 3, 6, 0, 0, 1.0, 0);
    SH(7, 3, 8, 0, 0, 1.0, 0);
    SH(8, 3, 8, 2, 0, 1.0, 0);
    SH(8, 3, 8, 0.5, 0, 1.0, 0);
    SH(9, 5, 6, 1, 0, 1.0, 0);
    SH(10, 5, 6, 4, 0, 1.0, 0);
    SH(11, 5, 6, 4, 0, 1.0, 0);
    SH(12, 2, 3, 0, 0, 1.0, 1);
    SH(12, 0, 1, 5, 0, 1.0, 0);
    SH(13, 2, 3, 0, 0, 1.0, 1);
    SH(14, 3, 4, 0, 0, 1.0, 1);
    SH(14, 0, 1, 50, 0, 1.0, 0);
    for (int t = 15; t <= 21; ++t)
        SH(t, 3, 0, 0, 0, 1.0, t == 15 ? 1 : 0);
    SH(24, 3, 130, 0, 0, 1.0, 1);

    // ---- single commands
    std::fprintf(f, "XS\nBX\nFB\nFM\nKC\nGC\nHP\nXP\n");

    // ---- root finder cases: fid x1 x2 rts xacc
    auto FR = [&](int fid, double x1, double x2, double rts, double xacc) {
        std::fprintf(f, "FR %d", fid);
        for (double v : {x1, x2, rts, xacc})
            wd(f, v);
        std::fprintf(f, "\n");
    };
    for (int fid = 0; fid <= 10; ++fid) {
        for (double xacc : {1e-3, 1e-8, 0.0, 0.5}) {
            FR(fid, 0.0, 3.0, 1.5, xacc);
            FR(fid, 3.0, 0.0, 1.5, xacc);
            FR(fid, -2.0, 4.0, 0.0, xacc);
            FR(fid, 1.0, 1.0, 1.0, xacc);
            FR(fid, 0.5, 2.5, 2.4, xacc);
            FR(fid, -10.0, 10.0, -9.0, xacc);
            FR(fid, 0.0, 1e-3, 5e-4, xacc);
            FR(fid, NaN, 2.0, 1.0, xacc);
            FR(fid, 0.0, INF, 1.0, xacc);
        }
        for (int k = 0; k < 10; ++k) {
            double a = u01() * 8 - 4, b = u01() * 8 - 4;
            FR(fid, a, b, a + (b - a) * u01(), 1e-6);
        }
    }
    std::fclose(f);
    return 0;
}

// ---------------------------------------------------------------- input
static std::map<std::string, std::vector<double>> LISTS;
static const std::vector<double>& LS_(const char* n) { return LISTS[n]; }

static double rdh(std::istringstream& is) {
    std::string s;
    is >> s;
    uint64_t u = std::stoull(s, nullptr, 16);
    double v;
    std::memcpy(&v, &u, 8);
    return v;
}
static std::string rdname(std::istringstream& is) {
    std::string s;
    is >> s;
    if (s == "~")
        return std::string();
    std::string r;
    for (size_t i = 0; i + 1 < s.size(); i += 2)
        r += (char)std::stoi(s.substr(i, 2), nullptr, 16);
    return r;
}
static std::vector<double> rdvec(std::istringstream& is) {
    size_t n;
    is >> n;
    std::vector<double> v;
    for (size_t i = 0; i < n; ++i)
        v.push_back(rdh(is));
    return v;
}

// ---------------------------------------------------------------- dumps
static void dumpLut(const xsect::LocateLut& L) {
    hx(L.t0);
    hx(L.scale);
    hi(L.j_last);
    for (int b = 0; b <= xsect::LocateLut::kBuckets; ++b)
        hi(L.lo[b]);
}
static void dumpTd(const TransectData& td) {
    hname(td.name);
    hx(td.n_left);
    hx(td.n_right);
    hx(td.n_channel);
    hx(td.length_factor);
    hi((long long)td.stations.size());
    harr(td.stations.data(), td.stations.size());
    hi((long long)td.elevations.size());
    harr(td.elevations.data(), td.elevations.size());
    hx(td.x_left_bank);
    hx(td.x_right_bank);
    hx(td.y_full);
    hx(td.w_max);
    hx(td.a_full);
    hx(td.r_full);
    hx(td.s_max);
    hx(td.a_max);
    harr(td.area_tbl, NT);
    harr(td.hrad_tbl, NT);
    harr(td.width_tbl, NT);
    dumpLut(td.area_lut);
}
// Fields, then whether each table pointer is null (0), points into `td` (1), or elsewhere (2).
static void dumpXs(const XSectParams& xs, const TransectData* td) {
    hi(xs.type);
    hi(xs.culvert_code);
    hi(xs.transect);
    for (double v : {xs.y_full, xs.w_max, xs.yw_max, xs.a_full, xs.r_full, xs.s_full, xs.s_max,
                     xs.y_bot, xs.a_bot, xs.s_bot, xs.r_bot})
        hx(v);
    hi(xs.transect_tbl_size);
    auto P = [&](const void* p, const void* mine) { hi(p == nullptr ? 0 : (p == mine ? 1 : 2)); };
    P(xs.area_tbl, td ? (const void*)td->area_tbl : nullptr);
    P(xs.hrad_tbl, td ? (const void*)td->hrad_tbl : nullptr);
    P(xs.width_tbl, td ? (const void*)td->width_tbl : nullptr);
    P(xs.area_lut, td ? (const void*)&td->area_lut : nullptr);
}

// ---------------------------------------------------------------- sweeps
static double base_of(double v) { return (std::isfinite(v) && v > 0.0) ? v : 1.0; }

static void cullines(const XSectParams& xs, const std::vector<int>& codes,
                     const std::vector<double>& slopes) {
    const double Yb = base_of(xs.y_full);
    const double Qb = base_of(xs.a_full * std::sqrt(32.2 * Yb));
    for (int code : codes) {
        for (double slope : slopes) {
            for (double q0 : {1e30, 0.5 * Qb, 0.0}) {
                hs("CUL");
                hi(code);
                hx(slope);
                hx(q0);
                for (double fr : LS_("CD")) {
                    double y = fr * Yb;
                    double dqdh = 7.0;
                    bool controls = true;
                    double q = culvert::getInflow(q0, y, xs, slope, code, dqdh, controls);
                    hx(q);
                    hx(dqdh);
                    hi(controls ? 1 : 0);
                }
                nl();
            }
        }
    }
}

static void linkSweep(const XSectParams& xs, int level) {
    const double Yb = base_of(xs.y_full);
    std::vector<double> D;
    for (double fr : LS_("LD"))
        D.push_back(fr * Yb);
    for (double flow : LS_("LQV")) {
        for (int b : {1, 0, 3, -2}) {
            hs("VEL");
            hx(flow);
            hi(b);
            for (double d : D)
                hx(link::getVelocity(xs, flow, d, b));
            nl();
        }
    }
    for (double v : LS_("LV")) {
        hs("FRO");
        hx(v);
        for (double d : D)
            hx(link::getFroude(xs, v, d));
        nl();
    }
    hs("CAP");
    for (double d : D)
        hx(link::getCapacity(xs, d));
    nl();
    const auto& LN = LS_("LN");
    const auto& LSl = LS_("LS");
    const size_t nn = level ? LN.size() : 2, ns = level ? LSl.size() : 2;
    for (size_t i = 0; i < nn; ++i) {
        for (size_t j = 0; j < ns; ++j) {
            double beta = -9, rf = -9, qf = -9;
            link::computeConveyance(LN[i], LSl[j], xs.s_full, beta, rf, qf);
            hs("CONV");
            hx(beta);
            hx(rf);
            hx(qf);
            nl();
            const double Qb = base_of(qf);
            for (double qmax : {-1.0, 0.0, 0.5 * Qb, Qb, std::numeric_limits<double>::quiet_NaN()}) {
                hs("DFF");
                hx(qmax);
                for (double fr : LS_("LQF"))
                    hx(link::getDepthFromFlow(xs, beta, fr * Qb, qmax));
                nl();
            }
        }
    }
    if (level)
        cullines(xs, [] {
            std::vector<int> c;
            for (int k = -1; k <= 58; ++k)
                c.push_back(k);
            return c;
        }(), LS_("CS"));
    else
        cullines(xs, {0, 1, 5, 12, 37, 50, 58}, {0.01});
}

static std::vector<TransectData> POOL;
static std::vector<XSectParams> SHXS;

static void tabulatedProbe(const TransectData& td, int idx) {
    for (int t : {22, 25, 23, 1}) {
        for (double yf : {td.y_full, 2.5}) {
            XSectParams xs{};
            xs.type = t;
            xs.a_bot = 0.123;
            link::applyTabulatedXSectParams(xs, td, yf, idx);
            hs("ATP");
            hi(t);
            hx(yf);
            dumpXs(xs, &td);
            nl();
        }
    }
}

// ---------------------------------------------------------------- root-finder test functions
static long long EVN;
static uint64_t EVH;
static void ev(double x) {
    uint64_t u;
    std::memcpy(&u, &x, 8);
    EVH = (EVH * 1099511628211ull) ^ u;
    ++EVN;
}
static double tf(int fid, double x, double* df) {
    ev(x);
    switch (fid) {
    case 0: *df = 3.0 * x * x - 2.0; return x * x * x - 2.0 * x - 5.0;
    case 1: *df = -std::sin(x) - 1.0; return std::cos(x) - x;
    case 2: *df = 2.0 * x; return x * x - 2.0;
    case 3: *df = std::exp(x); return std::exp(x) - 3.0;
    case 4: {
        double s = std::sqrt(1.0 + x * x);
        *df = 1.0 / (s * s * s);
        return x / s;
    }
    case 5: *df = 3.0 * (x - 1.0) * (x - 1.0); return (x - 1.0) * (x - 1.0) * (x - 1.0);
    case 6: *df = 0.0; return 0.0;
    case 7: *df = 0.0; return (x < 0.5) ? -1.0 : 1.0;
    case 8: *df = 0.5 / std::sqrt(x); return std::sqrt(x) - 1.0;
    case 9: *df = -1.0 / ((x - 0.3) * (x - 0.3)); return 1.0 / (x - 0.3);
    default: *df = 2.0 * x; return x * x + 1.0;
    }
}

// ---------------------------------------------------------------- run
static int run(const char* path) {
    std::ifstream in(path);
    if (!in)
        return 1;
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    std::string line;
    TransectStore ts;
    long cmd = 0;
    while (std::getline(in, line)) {
        std::istringstream is(line);
        std::string tag;
        is >> tag;
        ++cmd;
        if (tag == "L") {
            std::string name;
            is >> name;
            LISTS[name] = rdvec(is);
        } else if (tag == "TS") {
            ts = TransectStore{};
        } else if (tag == "TR") {
            ts.names.push_back(rdname(is));
            ts.comments.emplace_back();
            ts.n_left.push_back(rdh(is));
            ts.n_right.push_back(rdh(is));
            ts.n_channel.push_back(rdh(is));
            ts.x_left_bank.push_back(rdh(is));
            ts.x_right_bank.push_back(rdh(is));
            ts.x_left_encroachment.push_back(0.0);
            ts.x_right_encroachment.push_back(0.0);
            ts.x_factor.push_back(rdh(is));
            ts.y_factor.push_back(rdh(is));
            ts.length_factor.push_back(rdh(is));
            ts.stations.push_back(rdvec(is));
            ts.elevations.push_back(rdvec(is));
        } else if (tag == "TE") {
            std::vector<double> ucfs = rdvec(is);
            for (double u : ucfs) {
                for (int i = 0; i < ts.count(); ++i) {
                    TransectData td;
                    bool ret = transect::buildFromStore(ts, i, u, td);
                    hs("BFS");
                    hx(u);
                    hi(i);
                    hi(ret ? 1 : 0);
                    dumpTd(td);
                    nl();
                    if (ret) {
                        TransectData t2;
                        t2.name = td.name;
                        t2.n_left = td.n_left;
                        t2.n_right = td.n_right;
                        t2.n_channel = td.n_channel;
                        t2.length_factor = td.length_factor;
                        t2.stations = td.stations;
                        t2.elevations = td.elevations;
                        t2.x_left_bank = td.x_left_bank;
                        t2.x_right_bank = td.x_right_bank;
                        transect::buildTables(t2, false);
                        hs("BT0");
                        dumpTd(t2);
                        nl();
                    }
                    const int idx = (int)POOL.size();
                    POOL.push_back(td);
                    const TransectData& pt = POOL.back();
                    tabulatedProbe(pt, idx);
                    if (ret && pt.y_full > 0.0 && pt.a_full > 0.0 && std::isfinite(pt.y_full) &&
                        std::isfinite(pt.a_full)) {
                        XSectParams xs{};
                        xs.type = 22;
                        link::applyTabulatedXSectParams(xs, pt, pt.y_full, idx);
                        linkSweep(xs, 0);
                    }
                }
            }
        } else if (tag == "ST") {
            street::StreetParams sp;
            sp.width = rdh(is);
            sp.curb_height = rdh(is);
            sp.slope = rdh(is);
            sp.roughness = rdh(is);
            sp.gutter_depression = rdh(is);
            sp.gutter_width = rdh(is);
            is >> sp.sides;
            sp.back_width = rdh(is);
            sp.back_slope = rdh(is);
            sp.back_roughness = rdh(is);
            TransectData td;
            street::buildTransect(sp, td);
            hs("STR");
            dumpTd(td);
            nl();
            const int idx = (int)POOL.size();
            POOL.push_back(td);
            const TransectData& pt = POOL.back();
            tabulatedProbe(pt, idx);
            if (pt.y_full > 0.0 && pt.a_full > 0.0 && std::isfinite(pt.y_full) &&
                std::isfinite(pt.a_full)) {
                XSectParams xs{};
                xs.type = 25;
                link::applyTabulatedXSectParams(xs, pt, pt.y_full, idx);
                linkSweep(xs, 0);
            }
        } else if (tag == "CU") {
            std::string name = rdname(is);
            double yf = rdh(is), wm = rdh(is);
            int n;
            is >> n;
            std::vector<double> x = rdvec(is), y = rdvec(is);
            TransectData td;
            td.name = name;
            td.y_full = yf;
            td.a_full = 0.0;
            td.r_full = 0.0;
            td.w_max = wm;
            transect::buildCustomTables(td, yf, x.data(), y.data(), n);
            hs("CUS");
            dumpTd(td);
            nl();
            const int idx = (int)POOL.size();
            POOL.push_back(td);
            const TransectData& pt = POOL.back();
            tabulatedProbe(pt, idx);
            XSectParams xs{};
            xs.type = 23;
            link::applyTabulatedXSectParams(xs, pt, yf, idx);
            if (xs.y_full > 0.0 && xs.a_full > 0.0 && std::isfinite(xs.y_full) &&
                std::isfinite(xs.a_full) && pt.area_lut.scale != 0.0)
                linkSweep(xs, 0);
        } else if (tag == "SH") {
            int t, full;
            is >> t;
            double p[5] = {rdh(is), rdh(is), rdh(is), rdh(is), 0.0};
            double u = rdh(is);
            is >> full;
            XSectParams xs{};
            int rc = xsect::setParams(xs, t, p, u);
            hs("SHP");
            hi(t);
            hi(rc);
            dumpXs(xs, nullptr);
            nl();
            SHXS.push_back(xs);
            linkSweep(xs, full);
        } else if (tag == "XS") {
            hs("XSH");
            for (int s = 0; s <= 25; ++s)
                hi(link::translateShape(static_cast<XsectShape>(s)));
            nl();
        } else if (tag == "BX") {
            const int nP = (int)POOL.size();
            const int nS = (int)SHXS.size();
            const int n = nP + nS + 8;
            LinkData links;
            links.resize(n);
            const int shapes[] = {22, 23, 25, 1, 0, 24, 3};
            for (int i = 0; i < n; ++i) {
                auto u = (size_t)i;
                if (i < nS) {
                    const XSectParams& s = SHXS[u];
                    links.xsect_batch_shape[u] = s.type;
                    links.xsect_y_full[u] = s.y_full;
                    links.xsect_a_full[u] = s.a_full;
                    links.xsect_w_max[u] = s.w_max;
                    links.xsect_r_full[u] = s.r_full;
                    links.xsect_s_full[u] = s.s_full;
                    links.xsect_s_max[u] = s.s_max;
                    links.xsect_y_bot[u] = s.y_bot;
                    links.xsect_a_bot[u] = s.a_bot;
                    links.xsect_s_bot[u] = s.s_bot;
                    links.xsect_r_bot[u] = s.r_bot;
                    links.xsect_curve[u] = i - 2;
                } else {
                    links.xsect_batch_shape[u] = shapes[i % 7];
                    links.xsect_y_full[u] = i * 0.25;
                    links.xsect_a_full[u] = i * 0.5;
                    links.xsect_w_max[u] = i * 0.125;
                    links.xsect_r_full[u] = i * 0.0625;
                    links.xsect_s_full[u] = i * 1.5;
                    links.xsect_s_max[u] = i * 1.75;
                    links.xsect_y_bot[u] = i * 0.03125;
                    links.xsect_a_bot[u] = i * 2.0;
                    links.xsect_s_bot[u] = i * 3.0;
                    links.xsect_r_bot[u] = i * 4.0;
                    links.xsect_curve[u] = (i - nS) - 3;
                }
                for (int with : {0, 1}) {
                    XSectParams xs = link::buildXSectParams(links, u, with ? &POOL : nullptr);
                    const TransectData* td =
                        (xs.transect >= 0 && xs.transect < nP) ? &POOL[(size_t)xs.transect] : nullptr;
                    hs("BXP");
                    hi(i);
                    hi(with);
                    dumpXs(xs, td);
                    nl();
                }
            }
        } else if (tag == "FB") {
            LinkData links;
            links.resize(5 * 26);
            for (int i = 0; i < 5 * 26; ++i) {
                links.type[(size_t)i] = static_cast<LinkType>(i / 26);
                links.xsect_shape[(size_t)i] = static_cast<XsectShape>(i % 26);
                links.xsect_y_bot[(size_t)i] = i * 0.125 + 0.5;
            }
            hs("FBP");
            for (int i = 0; i < 5 * 26; ++i)
                for (bool r : {false, true})
                    hx(link::filledCircularOffsetBump(links, (size_t)i, r));
            nl();
        } else if (tag == "FM") {
            using forcemain::FrictionModel;
            const auto &FV = LS_("FV"), &FRr = LS_("FR"), &FP = LS_("FP");
            for (double v : FV)
                for (double r : FRr) {
                    hs("FHW");
                    for (double c : FP) {
                        hx(forcemain::getFricSlope_HW(v, r, c));
                        hx(hydkernels::fricSlopeHW(v, r, c));
                    }
                    nl();
                    hs("FDW");
                    for (double c : FP) {
                        hx(forcemain::getFricSlope_DW(v, r, c));
                        hx(hydkernels::fricSlopeDW(v, r, c));
                    }
                    nl();
                }
            for (double e : LS_("FE"))
                for (double h : LS_("FH")) {
                    hs("FFF");
                    for (double re : LS_("FRE"))
                        hx(forcemain::getFricFactor(e, h, re));
                    nl();
                }
            for (int m : {0, 1})
                for (double rb : LS_("FB")) {
                    for (double yf : LS_("FY")) {
                        hs("FEN");
                        for (double s : LS_("FS"))
                            hx(forcemain::getEquivN(static_cast<FrictionModel>(m), rb, yf, s, 0.013));
                        nl();
                    }
                    hs("FRF");
                    for (double lf : LS_("FL"))
                        hx(forcemain::getRoughFactor(static_cast<FrictionModel>(m), rb, lf));
                    nl();
                }
            const size_t nb = std::min({FV.size(), FRr.size(), FP.size()});
            for (int m : {0, 1})
                for (int cnt : {(int)nb, 0, -1, 3}) {
                    std::vector<double> out(nb, -7.0);
                    forcemain::batchFricSlope(FV.data(), FRr.data(), FP.data(), out.data(),
                                              static_cast<FrictionModel>(m), cnt);
                    hs("FBS");
                    harr(out.data(), out.size());
                    nl();
                }
        } else if (tag == "KC") {
            std::vector<hydkernels::CulvertCurve> curves;
            for (int c = 0; c <= 57; ++c) {
                culvert::CulvertCoeffs k = culvert::getCoeffs(c);
                curves.push_back({k.K, k.M, k.C, k.Y});
            }
            curves.push_back({0.01, 0.0, 0.03, 0.7});
            curves.push_back({0.01, 2.0, 0.0, 0.7});
            curves.push_back({0.01, 2.0, -0.03, 0.7});
            curves.push_back({-0.01, 2.0, 0.03, 0.7});
            curves.push_back({NaN, 2.0, 0.03, 0.7});
            curves.push_back({0.5, 0.667, 0.03, 0.3});
            for (const auto& cc : curves)
                for (bool mit : {false, true})
                    for (double s : LS_("KS")) {
                        for (double yf : {3.0, 0.0, NaN}) {
                            for (double q : {1e30, 5.0}) {
                                hs("KCI");
                                for (double h : LS_("KH")) {
                                    double dq = 7.0;
                                    hx(hydkernels::culvertInflow(q, h, yf, 7.0685835, s, cc, mit, dq));
                                    hx(dq);
                                }
                                nl();
                            }
                        }
                        hs("KCA");
                        for (double af : {0.0, -1.0, NaN, 1e-12}) {
                            double dq = 7.0;
                            hx(hydkernels::culvertInflow(1e30, 1.5, 3.0, af, s, cc, mit, dq));
                            hx(dq);
                        }
                        nl();
                    }
        } else if (tag == "GC") {
            hs("GCO");
            for (int c = -2; c <= 60; ++c) {
                culvert::CulvertCoeffs k = culvert::getCoeffs(c);
                for (double v : {k.form, k.K, k.M, k.C, k.Y})
                    hx(v);
            }
            nl();
        } else if (tag == "XP") {
            // Hand-made sections (y_full, a_full overridden after setParams): reaches the
            // a_full <= 0 arm of getCapacity and the BIG return of the submerged culvert flow.
            for (int t : {3, 1})
                for (double yf : {3.0, 0.0, -3.0, NaN})
                    for (double af : {7.0, 0.0, -1.0, NaN}) {
                        XSectParams xs{};
                        double p[5] = {3.0, 2.0, 0.0, 0.0, 0.0};
                        xsect::setParams(xs, t, p, 1.0);
                        xs.y_full = yf;
                        xs.a_full = af;
                        hs("XPP");
                        hi(t);
                        hx(yf);
                        hx(af);
                        for (double fr : LS_("LD")) {
                            const double d = fr * 3.0;
                            hx(link::getCapacity(xs, d));
                            hx(link::getVelocity(xs, 2.0, d, 1));
                            hx(link::getFroude(xs, 1.5, d));
                        }
                        nl();
                        for (int code : {1, 12, 30, 37}) {
                            hs("XPC");
                            hi(code);
                            for (double fr : LS_("CD")) {
                                double dq = 7.0;
                                bool ctl = true;
                                hx(culvert::getInflow(1e30, fr * 3.0, xs, 0.01, code, dq, ctl));
                                hx(dq);
                                hi(ctl ? 1 : 0);
                            }
                            nl();
                        }
                    }
        } else if (tag == "HP") {
            const auto& H = LS_("FV");
            for (double q : LS_("LQF")) {
                hs("HPW");
                for (double h1 : H)
                    for (double h2 : H)
                        hx(link::getHydPower(q, h1, h2));
                nl();
            }
        } else if (tag == "FR") {
            int fid;
            is >> fid;
            double x1 = rdh(is), x2 = rdh(is), rts0 = rdh(is), xacc = rdh(is);
            EVN = 0;
            EVH = 1469598103934665603ull;
            double rts = rts0;
            int n = findroot::newton(x1, x2, &rts, xacc, [fid](double x, double* f, double* df) {
                *f = tf(fid, x, df);
            });
            hs("FRN");
            hi(n);
            hx(rts);
            hi(EVN);
            hi((long long)EVH);
            EVN = 0;
            EVH = 1469598103934665603ull;
            double r = findroot::ridder(x1, x2, xacc, [fid](double x) {
                double d;
                return tf(fid, x, &d);
            });
            hs("FRR");
            hx(r);
            hi(EVN);
            hi((long long)EVH);
            nl();
        } else if (!tag.empty()) {
            std::fprintf(stderr, "unknown command %s at %ld\n", tag.c_str(), cmd);
            return 2;
        }
    }
    std::fwrite(OUT.data(), 1, OUT.size(), stdout);
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 3 && std::string(argv[1]) == "gen")
        return gen(argv[2]);
    if (argc == 2)
        return run(argv[1]);
    std::fprintf(stderr, "usage: hydfuncs gen FILE | hydfuncs FILE\n");
    return 1;
}
