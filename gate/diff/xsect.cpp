// Differential harness for the cross-section geometry core: hydraulics/xsect_tables.hpp,
// XSectKernels.hpp, XSection.cpp and XSectBatch.hpp/.cpp. Upstream's compiled engine
// library supplies every out-of-line function; the inline ones (XsectEval members called
// directly, xsect_tables) compile into this harness with the engine's flags.
//
//   xsect gen FILE   write the inputs to FILE: depth/area fractions (0, -0, subnormal,
//                    tiny, every table grid point and its two neighbours, a dense interior,
//                    near full, full, above full, negative, NaN, infinities, random),
//                    absolute values, flows, four synthetic transect tables, and ~300
//                    setParams cases (every shape code, typical/tiny/huge/zero/negative/
//                    NaN/inf dimensions, other units, the combinations setParams rejects,
//                    unknown codes).
//   xsect FILE       print, as hex bits (ints in decimal): every built-in table and bucket
//                    map; lookup/invLookup/locate on every table; the circular helpers;
//                    for each case setParams' return and fields, then getAofY, getWofY,
//                    getRofY, getYofA, getSofA, getRofA, getdSdA, getAofS, getYcrit,
//                    getAmax, isOpen over the sweeps; then the batch API (build,
//                    attachTransectTables, every compute*, bypass masks, team slices,
//                    every xsect_batch kernel) over a network made of the cases.
//
// A call on which upstream is undefined (a table read out of bounds after a non-finite
// int cast) or does not return (unbounded recursion for DUMMY/unknown codes, the
// power-function perimeter loop with a zero step, getYcrit's enumeration from INT_MIN)
// is not made; it prints "-", and the Rust side applies the same predicate.
//
// Build (from the repository root):
//   U=/root/work/swmm6src/upstream/src/engine; B=/root/work/swmm6src/build
//   g++ -std=c++20 -O3 -ffp-contract=off -mavx2 -mfma -DNDEBUG -DSWMM_OS_LINUX -I$U \
//       -I$B/include gate/diff/xsect.cpp $B/src/engine/libopenswmm.engine.so \
//       -Wl,-rpath,$B/src/engine -o gate/diff/xsect
#include "hydraulics/XSectBatch.hpp"
#include "hydraulics/XSectKernels.hpp"
#include "hydraulics/xsect_tables.hpp"
#include "hydraulics/Transect.hpp"
#include "core/SimulationContext.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

// Defined with external linkage in XSectBatch.cpp but not declared in its header.
namespace openswmm::xsect_batch {
void hydrad_rect_closed(const double*, const double*, const double*, double*, int);
void hydrad_rect_open(const double*, const double*, const double*, double*, int);
void width_rect_closed(const double*, const double*, const double*, double*, int);
} // namespace openswmm::xsect_batch

using namespace openswmm;
namespace xt = openswmm::xsect_tables;
namespace xb = openswmm::xsect_batch;

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
static void dash() { OUT += "- "; }
static void nl() { OUT += '\n'; }
static void harr(const double* a, int n) {
    for (int i = 0; i < n; ++i)
        hx(a[i]);
}

// ---------------------------------------------------------------- generator
static uint64_t g_s = 0x5eed2026ull;
static uint64_t rnd() {
    uint64_t z = (g_s += 0x9e3779b97f4a7c15ull);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
    return z ^ (z >> 31);
}
static double u01() { return (double)(rnd() >> 11) * 0x1.0p-53; }

static void wd(FILE* f, double v) {
    uint64_t u;
    std::memcpy(&u, &v, 8);
    std::fprintf(f, " %016llx", (unsigned long long)u);
}

static int gen(const char* path) {
    FILE* f = std::fopen(path, "w");
    if (!f)
        return 1;
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    const double INF = std::numeric_limits<double>::infinity();

    std::vector<double> F = {0.0,     -0.0,  5e-324, 1e-300, 1e-30,  1e-12,     1e-9,
                             1e-6,    1e-5,  1.0000001e-5,   9.999999e-6, 1e-4, 1e-3,
                             0.005,   0.01,  0.015,  0.02,   0.03,   0.04,      0.0399999,
                             0.0400001};
    for (int k = 1; k < 100; ++k)
        F.push_back(k / 100.0);
    for (int j = 0; j <= 50; ++j) {
        double g = j / 50.0;
        F.push_back(std::nextafter(g, -1.0));
        F.push_back(g);
        F.push_back(std::nextafter(g, 2.0));
    }
    for (double v : {1.0 - 1e-12, 1.01, 1.08, 1.5, 2.0, 3.0, 10.0, 100.0, -1e-6, -0.001, -0.02,
                     -0.5, -1.0, -100.0, NaN, INF, -INF, 1e300, 1e10, 1e6})
        F.push_back(v);
    for (int k = 0; k < 40; ++k)
        F.push_back(u01() * 1.2);
    std::fprintf(f, "F %zu", F.size());
    for (double v : F)
        wd(f, v);
    std::fprintf(f, "\n");

    std::vector<double> A = {0.0, 1e-8, 0.1, 0.5, 1.0, 2.5, 7.0, 30.0, 1e4, -2.0};
    std::fprintf(f, "A %zu", A.size());
    for (double v : A)
        wd(f, v);
    std::fprintf(f, "\n");

    std::vector<double> Q = {0.0,  -1.0, 1e-300, 1e-160, 1e-6,   1e-3, 0.01, 0.1, 0.5, 1.0, 2.0,
                             5.0,  10.0, 30.0,   100.0,  300.0,  1e3,  1e4,  1e6, 1e155, NaN, INF};
    std::fprintf(f, "Q %zu", Q.size());
    for (double v : Q)
        wd(f, v);
    std::fprintf(f, "\n");

    // ---- four synthetic transect tables (51 rows, normalized)
    const int N = transect::N_TRANSECT_TBL;
    auto putT = [&](int id, double yf, double wm, double af, double rf, double sf, double sm,
                    double ab, const std::vector<double>& ar, const std::vector<double>& hr,
                    const std::vector<double>& wt) {
        std::fprintf(f, "T %d", id);
        for (double v : {yf, wm, af, rf, sf, sm, ab})
            wd(f, v);
        for (double v : ar)
            wd(f, v);
        for (double v : hr)
            wd(f, v);
        for (double v : wt)
            wd(f, v);
        std::fprintf(f, "\n");
    };
    {   // T0: trapezoid, bottom 4 ft, side slope 2, 3 ft deep
        std::vector<double> ar(N), hr(N), wt(N);
        double Af = (4.0 + 2.0 * 3.0) * 3.0, Pf = 4.0 + 2.0 * 3.0 * std::sqrt(5.0), Rf = Af / Pf;
        for (int j = 0; j < N; ++j) {
            double y = 3.0 * j / 50.0, a = (4.0 + 2.0 * y) * y, p = 4.0 + 2.0 * y * std::sqrt(5.0);
            ar[j] = a / Af;
            hr[j] = (j == 0 ? 0.0 : a / p) / Rf;
            wt[j] = (4.0 + 4.0 * y) / 16.0;
        }
        double sf = Af * std::pow(Rf, 2.0 / 3.0);
        putT(0, 3.0, 16.0, Af, Rf, sf, sf, Af, ar, hr, wt);
    }
    {   // T1: the circular tables
        std::vector<double> ar(xt::A_Circ, xt::A_Circ + N), hr(xt::R_Circ, xt::R_Circ + N),
            wt(xt::W_Circ, xt::W_Circ + N);
        double af = 3.14159 * 1.0, rf = 0.5, sf = af * std::pow(rf, 2.0 / 3.0);
        putT(1, 2.0, 2.0, af, rf, sf, 1.08 * sf, 0.9756 * af, ar, hr, wt);
    }
    {   // T2: random monotone area with flat runs, random radius and width
        std::vector<double> ar(N), hr(N), wt(N);
        ar[0] = 0.0;
        for (int j = 1; j < N; ++j)
            ar[j] = ar[j - 1] + ((j % 7 == 3) ? 0.0 : u01());
        double top = ar[N - 1];
        for (int j = 0; j < N; ++j) {
            ar[j] /= top;
            hr[j] = 0.01 + u01();
            wt[j] = 0.2 + u01();
        }
        double sf = 4.0 * std::pow(0.6, 2.0 / 3.0);
        putT(2, 1.7, 5.0, 4.0, 0.6, sf, 1.1 * sf, 3.8, ar, hr, wt);
    }
    {   // T3: non-monotone area (no bucket map), declining radius top
        std::vector<double> ar(N), hr(N), wt(N);
        for (int j = 0; j < N; ++j) {
            ar[j] = j / 50.0 + 0.05 * std::sin((double)j);
            hr[j] = 1.2 - std::fabs(j / 50.0 - 0.8);
            wt[j] = 1.0 - j / 60.0;
        }
        double sf = 5.0 * std::pow(0.7, 2.0 / 3.0);
        putT(3, 2.5, 3.0, 5.0, 0.7, sf, sf, 5.0, ar, hr, wt);
    }

    // ---- setParams cases: type p0 p1 p2 p3 ucf tid
    auto P = [&](int t, double a, double b, double c, double d, double ucf, int tid) {
        std::fprintf(f, "P %d", t);
        for (double v : {a, b, c, d, ucf})
            wd(f, v);
        std::fprintf(f, " %d\n", tid);
    };
    struct B4 {
        double p[4];
    };
    auto bases = [&](int t) -> std::vector<B4> {
        switch (t) {
        case 0: return {{{0, 0, 0, 0}}, {{3, 4, 5, 6}}};
        case 1: return {{{3, 0, 0, 0}}, {{0.25, 0, 0, 0}}};
        case 2: return {{{3, 0.5, 0, 0}}, {{3, 0, 0, 0}}, {{3, 2.9, 0, 0}}, {{3, 3, 0, 0}}, {{3, 4, 0, 0}}};
        case 3: return {{{4, 6, 0, 0}}, {{1, 10, 0, 0}}};
        case 4: return {{{4, 6, 0, 0}}, {{4, 6, 1, 0}}, {{4, 6, 2, 0}}};
        case 5: return {{{5, 10, 2, 3}}, {{5, 10, 0, 0}}, {{5, 0, 1, 1}}, {{5, 10, 0.5, 0}}};
        case 6: return {{{3, 6, 0, 0}}, {{2, 0.5, 0, 0}}};
        case 7: return {{{3, 8, 0, 0}}, {{1, 0.2, 0, 0}}};
        case 8: return {{{3, 8, 2, 0}}, {{3, 8, 0.5, 0}}, {{3, 8, 1, 0}}, {{3, 8, 0, 0}}};
        case 9: return {{{5, 6, 1, 0}}, {{5, 6, 4.9, 0}}, {{5, 0, 1, 0}}, {{5, 6, 0, 0}}, {{1, 6, 3, 0}}};
        case 10: return {{{5, 6, 4, 0}}, {{5, 6, 0, 0}}, {{5, 6, 3, 0}}, {{5, 0, 1, 0}}, {{1, 6, 3, 0}}};
        case 11: return {{{5, 6, 4, 0}}, {{5, 6, 0, 0}}, {{5, 0, 1, 0}}, {{1, 6, 3, 0}}};
        case 12:
        case 13:
            return {{{2, 3, 0, 0}},  {{0, 1, 5, 0}},    {{7, 0, 0, 0}},   {{0, 1, 23, 0}},
                    {{0, 1, 1.9, 0}}, {{0, 1, 24, 0}},  {{0, 1, 1e12, 0}}, {{0, 1, INF, 0}},
                    {{2, -1, 0, 0}},  {{0, 0, 0, 0}}};
        case 14:
            return {{{3, 4, 0, 0}},   {{0, 1, 50, 0}},  {{102, 0, 0, 0}}, {{0, 1, 102, 0}},
                    {{0, 1, 103, 0}}, {{3, -1, 0, 0}},  {{0, 1, 1e12, 0}}};
        case 24: return {{{3, 130, 0, 0}}, {{1, 0.01, 0, 0}}};
        default:
            if (t >= 15 && t <= 21)
                return {{{3, 0, 0, 0}}, {{0.5, 0, 0, 0}}};
            return {{{3, 4, 5, 6}}};
        }
    };
    std::vector<int> types;
    for (int t = 0; t <= 25; ++t)
        types.push_back(t);
    for (int t : {-1, 26, 100})
        types.push_back(t);
    for (int t : types) {
        const bool tab = (t == 22 || t == 23 || t == 25);
        std::vector<B4> bs = bases(t);
        if (tab) {
            for (int tid : {0, 1, 2, 3, -1, 99})
                P(t, 0, 0, 0, 0, 1.0, tid);
            continue;
        }
        for (auto& b : bs)
            P(t, b.p[0], b.p[1], b.p[2], b.p[3], 1.0, -1);
        const double* p = bs[0].p;
        P(t, p[0], p[1], p[2], p[3], 12.0, -1);
        P(t, p[0], p[1], p[2], p[3], 0.3048, -1);
        P(t, p[0] * 1e-6, p[1] * 1e-6, p[2], p[3], 1.0, -1);
        P(t, p[0] * 1e5, p[1] * 1e5, p[2], p[3], 1.0, -1);
        P(t, 0.0, p[1], p[2], p[3], 1.0, -1);
        P(t, -p[0], p[1], p[2], p[3], 1.0, -1);
        P(t, NaN, p[1], p[2], p[3], 1.0, -1);
        P(t, INF, p[1], p[2], p[3], 1.0, -1);
        P(t, 0.0, 0.0, 0.0, 0.0, 1.0, -1);
    }
    std::fclose(f);
    return 0;
}

// ---------------------------------------------------------------- input
struct TRec {
    double yf, wm, af, rf, sf, sm, ab;
    std::vector<double> ar, hr, wt;
};
struct PRec {
    int type;
    double p[4];
    double ucf;
    int tid;
};
static std::vector<double> Fv, Av, Qv;
static std::vector<TRec> Tv;
static std::vector<PRec> Pv;

static double rdh(std::istringstream& is) {
    std::string s;
    is >> s;
    uint64_t u = std::stoull(s, nullptr, 16);
    double v;
    std::memcpy(&v, &u, 8);
    return v;
}
static int load(const char* path) {
    std::ifstream in(path);
    if (!in)
        return 1;
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream is(line);
        std::string tag;
        is >> tag;
        if (tag == "F" || tag == "A" || tag == "Q") {
            size_t n;
            is >> n;
            auto& v = tag == "F" ? Fv : tag == "A" ? Av : Qv;
            for (size_t i = 0; i < n; ++i)
                v.push_back(rdh(is));
        } else if (tag == "T") {
            int id;
            is >> id;
            TRec t;
            t.yf = rdh(is);
            t.wm = rdh(is);
            t.af = rdh(is);
            t.rf = rdh(is);
            t.sf = rdh(is);
            t.sm = rdh(is);
            t.ab = rdh(is);
            const int N = transect::N_TRANSECT_TBL;
            for (int j = 0; j < N; ++j)
                t.ar.push_back(rdh(is));
            for (int j = 0; j < N; ++j)
                t.hr.push_back(rdh(is));
            for (int j = 0; j < N; ++j)
                t.wt.push_back(rdh(is));
            Tv.push_back(t);
        } else if (tag == "P") {
            PRec p;
            is >> p.type;
            for (int k = 0; k < 4; ++k)
                p.p[k] = rdh(is);
            p.ucf = rdh(is);
            is >> p.tid;
            Pv.push_back(p);
        }
    }
    return 0;
}

// ---------------------------------------------------------------- predicates
static bool isTabS(int t) { return t >= 15 && t <= 21; }
static bool isRecType(int t) { return t == 0 || t < 0 || t > 25; } // DUMMY or unknown
static bool isTabulated(int t) { return t == 22 || t == 23 || t == 25; }

// powerfunc_getPofY(xs, y) terminates (and in a bounded number of steps).
static bool pofy_ok(const XSectParams& xs, double y) {
    if (!(y > 0.0))
        return true;
    double dy1 = 0.02 * xs.y_full;
    if (std::isnan(dy1))
        return true;
    if (!(dy1 > 0.0))
        return false;
    if (std::isinf(dy1))
        return true;
    return std::isfinite(y) && y / dy1 <= 1.0e4;
}
static double pf_yofa(const XSectParams& xs, double a) {
    return std::pow(a / xs.r_bot, 1.0 / (xs.s_bot + 1.0));
}
static bool pf_sane(const XSectParams& xs) {
    return xs.y_full > 0.0 && std::isfinite(xs.y_full) && xs.a_full > 0.0 &&
           std::isfinite(xs.a_full) && xs.r_bot > 0.0 && std::isfinite(xs.r_bot) &&
           xs.s_bot > 0.0 && std::isfinite(xs.s_bot);
}
static bool ok_RofY(const XSectParams& xs, double y) {
    if (xs.type == 8)
        return pofy_ok(xs, y);
    return true;
}
static bool ok_RofA(const XSectParams& xs, double a) {
    if (isRecType(xs.type))
        return a <= 0.0;
    if (xs.type == 8) // powerfunc_getRofA's guard is `a <= 0`: NaN reaches the loop
        return a <= 0.0 || pofy_ok(xs, pf_yofa(xs, a));
    return true;
}
static bool ok_SofA(const XSectParams& xs, double a) {
    if (isRecType(xs.type))
        return a <= 0.0;
    if (xs.type == 8)
        return a == 0.0 || ok_RofA(xs, a);
    return true;
}
static bool tab_ok(double alpha, int n) {
    double v = alpha / (1.0 / static_cast<double>(n - 1));
    return v > -1.0 && v < 2147483648.0;
}
static bool ok_dSdA(const XSectParams& xs, double a) {
    int t = xs.type;
    if (t == 1 || t == 24) {
        double alpha = a / xs.a_full;
        if (alpha <= 1.0e-30 || alpha < 0.04)
            return true;
        return tab_ok(alpha, 51);
    }
    if (isTabS(t))
        return tab_ok(a / xs.a_full, 51);
    if (t == 3 || t == 4 || t == 9 || t == 10 || t == 11 || t == 5 || t == 6)
        return true;
    double alpha = a / xs.a_full;
    double alpha1 = alpha - 0.001;
    double alpha2 = alpha + 0.001;
    if (alpha1 < 0.0)
        alpha1 = 0.0;
    double a1 = alpha1 * xs.a_full;
    double a2 = alpha2 * xs.a_full;
    return ok_SofA(xs, a2) && ok_SofA(xs, a1);
}
static bool ok_AofS(const XSectParams& xs, double s) {
    if (s <= 0.0)
        return true;
    int t = xs.type;
    if (t == 0)
        return true;
    if (t < 0 || t > 25)
        return false;
    if (t == 8)
        return pf_sane(xs);
    return true;
}
static bool ok_Ycrit(const XSectParams& xs, double q) {
    if (q <= 0.0)
        return true;
    double q2g = q * q / xsect::GRAVITY;
    if (q2g == 0.0)
        return true;
    int t = xs.type;
    if (t == 0 || t == 3 || t == 4 || t == 6 || t == 7 || t == 8)
        return true;
    double y0 = 1.01 * std::pow(q2g / xs.y_full, 0.25);
    if (y0 >= xs.y_full)
        y0 = 0.97 * xs.y_full;
    double r = xs.a_full / (xsect::PI / 4.0 * xs.y_full * xs.y_full);
    if (r >= 0.5 && r <= 2.0) {
        double dy = xs.y_full / 25;
        double v = y0 / dy;
        return v >= -100.0 && v <= 100.0;
    }
    return true;
}
static bool batch_ok(const XSectParams& xs) { return xs.type != 8 || pf_sane(xs); }
static double base_of(double v) { return (std::isfinite(v) && v > 0.0) ? v : 1.0; }

// ---------------------------------------------------------------- tables
struct TL {
    const char* name;
    const double* t;
    int n;
    int lut; // LutId or -1
};
static std::vector<TL> tableList() {
    const auto& T = xsect::hostTables();
    return {
        {"A_Arch", T.A_Arch, T.N_A_Arch, 6},
        {"A_Baskethandle", T.A_Baskethandle, T.N_A_Baskethandle, -1},
        {"A_Circ", T.A_Circ, T.N_A_Circ, -1},
        {"A_Egg", T.A_Egg, T.N_A_Egg, -1},
        {"A_HorizEllipse", T.A_HorizEllipse, T.N_A_HorizEllipse, 4},
        {"A_Horseshoe", T.A_Horseshoe, T.N_A_Horseshoe, -1},
        {"A_VertEllipse", T.A_VertEllipse, T.N_A_VertEllipse, 5},
        {"R_Arch", T.R_Arch, T.N_R_Arch, -1},
        {"R_Baskethandle", T.R_Baskethandle, T.N_R_Baskethandle, -1},
        {"R_Circ", T.R_Circ, T.N_R_Circ, -1},
        {"R_Egg", T.R_Egg, T.N_R_Egg, -1},
        {"R_HorizEllipse", T.R_HorizEllipse, T.N_R_HorizEllipse, -1},
        {"R_Horseshoe", T.R_Horseshoe, T.N_R_Horseshoe, -1},
        {"R_VertEllipse", T.R_VertEllipse, T.N_R_VertEllipse, -1},
        {"S_BasketHandle", T.S_BasketHandle, T.N_S_BasketHandle, 13},
        {"S_Catenary", T.S_Catenary, T.N_S_Catenary, 11},
        {"S_Circ", T.S_Circ, T.N_S_Circ, 7},
        {"S_Egg", T.S_Egg, T.N_S_Egg, 8},
        {"S_Gothic", T.S_Gothic, T.N_S_Gothic, 10},
        {"S_Horseshoe", T.S_Horseshoe, T.N_S_Horseshoe, 9},
        {"S_SemiCirc", T.S_SemiCirc, T.N_S_SemiCirc, 14},
        {"S_SemiEllip", T.S_SemiEllip, T.N_S_SemiEllip, 12},
        {"W_Arch", T.W_Arch, T.N_W_Arch, -1},
        {"W_BasketHandle", T.W_BasketHandle, T.N_W_BasketHandle, -1},
        {"W_Catenary", T.W_Catenary, T.N_W_Catenary, -1},
        {"W_Circ", T.W_Circ, T.N_W_Circ, -1},
        {"W_Egg", T.W_Egg, T.N_W_Egg, -1},
        {"W_Gothic", T.W_Gothic, T.N_W_Gothic, -1},
        {"W_HorizEllipse", T.W_HorizEllipse, T.N_W_HorizEllipse, -1},
        {"W_Horseshoe", T.W_Horseshoe, T.N_W_Horseshoe, -1},
        {"W_SemiCirc", T.W_SemiCirc, T.N_W_SemiCirc, -1},
        {"W_SemiEllip", T.W_SemiEllip, T.N_W_SemiEllip, -1},
        {"W_VertEllipse", T.W_VertEllipse, T.N_W_VertEllipse, -1},
        {"Y_BasketHandle", T.Y_BasketHandle, T.N_Y_BasketHandle, -1},
        {"Y_Catenary", T.Y_Catenary, T.N_Y_Catenary, 1},
        {"Y_Circ", T.Y_Circ, T.N_Y_Circ, -1},
        {"Y_Egg", T.Y_Egg, T.N_Y_Egg, -1},
        {"Y_Gothic", T.Y_Gothic, T.N_Y_Gothic, 0},
        {"Y_Horseshoe", T.Y_Horseshoe, T.N_Y_Horseshoe, -1},
        {"Y_SemiCirc", T.Y_SemiCirc, T.N_Y_SemiCirc, 3},
        {"Y_SemiEllip", T.Y_SemiEllip, T.N_Y_SemiEllip, 2},
    };
}
static void dumpLut(const xsect::LocateLut& L) {
    hx(L.t0);
    hx(L.scale);
    hi(L.j_last);
    for (int b = 0; b <= xsect::LocateLut::kBuckets; ++b)
        hi(L.lo[b]);
}

static void dumpXs(const XSectParams& xs) {
    hi(xs.type);
    hi(xs.culvert_code);
    hi(xs.transect);
    for (double v : {xs.y_full, xs.w_max, xs.yw_max, xs.a_full, xs.r_full, xs.s_full, xs.s_max,
                     xs.y_bot, xs.a_bot, xs.s_bot, xs.r_bot})
        hx(v);
}

static void dumpGroups(const char* tag, const XSectGroups& G) {
    hs("G");
    hs(tag);
    hi(G.numGroups());
    nl();
    for (int gi = 0; gi < G.numGroups(); ++gi) {
        const ShapeGroup& g = G.group(gi);
        hs("g");
        hi(gi);
        hi((int)g.shape);
        hi(g.count);
        hi((long long)g.link_idx.size());
        for (int v : g.link_idx)
            hi(v);
        for (const auto* col : {&g.y_full, &g.a_full, &g.r_full, &g.s_full, &g.w_max,
                                &g.inv_y_full, &g.y_bot, &g.a_bot, &g.s_bot, &g.r_bot}) {
            hi((long long)col->size());
            for (double v : *col)
                hx(v);
        }
        for (const auto* tv : {&g.area_tables, &g.hrad_tables, &g.width_tables}) {
            hi((long long)tv->size());
            for (const double* p : *tv) {
                if (!p)
                    hs("n");
                else {
                    hx(p[0]);
                    hx(p[17]);
                    hx(p[50]);
                }
            }
        }
        hi(g.transect_tbl_size);
        hi((long long)g.buf_d.size());
        hi((long long)g.buf_r.size());
        hi((long long)g.buf_r2.size());
        nl();
    }
    hs("find");
    for (int s = 0; s <= 25; ++s) {
        const ShapeGroup* p = G.findGroup(static_cast<XSectShape>(s));
        hi(p ? p->count : -1);
    }
    nl();
}

static const double SENT = -12345.678;

static void prline(const char* tag, int vi, const std::vector<double>& v) {
    hs(tag);
    hi(vi);
    for (double x : v)
        hx(x);
    nl();
}

// ---------------------------------------------------------------- run
static int run(const char* path) {
    if (load(path))
        return 1;
    const int N = transect::N_TRANSECT_TBL;
    const int nT = (int)Tv.size();

    // ---- tables and bucket maps
    auto TLv = tableList();
    for (const auto& t : TLv) {
        hs("TBL");
        hs(t.name);
        hi(t.n);
        harr(t.t, t.n);
        nl();
    }
    hs("TBL Amax");
    harr(xt::Amax, 26);
    nl();
    hs("TBL Ellipse");
    hi(xt::NumCodesEllipse);
    harr(xt::MinorAxis_Ellipse, 23);
    harr(xt::MajorAxis_Ellipse, 23);
    harr(xt::Afull_Ellipse, 23);
    harr(xt::Rfull_Ellipse, 23);
    nl();
    hs("TBL Arch");
    hi(xt::NumCodesArch);
    harr(xt::Yfull_Arch, 102);
    harr(xt::Wmax_Arch, 102);
    harr(xt::Afull_Arch, 102);
    harr(xt::Rfull_Arch, 102);
    nl();
    {
        const auto& T = xsect::hostTables();
        hs("TBL hostAmax");
        harr(T.Amax, 26);
        nl();
        for (int id = 0; id < (int)xsect::LutId::COUNT; ++id) {
            hs("LUT");
            hi(id);
            dumpLut(*T.lut(static_cast<xsect::LutId>(id)));
            nl();
        }
    }
    hs("CONST");
    for (double v : {xsect::TINY, xsect::PI, xsect::GRAVITY, xsect::RECT_ALFMAX,
                     xsect::RECT_TRIANG_ALFMAX, xsect::RECT_ROUND_ALFMAX})
        hx(v);
    nl();

    // ---- lookup helpers on every table
    std::vector<double> FA = Fv;
    FA.insert(FA.end(), Av.begin(), Av.end());
    for (const auto& t : TLv) {
        hs("L");
        hs(t.name);
        for (double x : FA)
            hx(xsect::lookup(x, t.t, t.n));
        nl();
        hs("I");
        hs(t.name);
        for (double y : FA)
            hx(xsect::invLookup(y, t.t, t.n, nullptr));
        nl();
        if (t.lut >= 0) {
            const xsect::LocateLut* L = xsect::hostTables().lut(static_cast<xsect::LutId>(t.lut));
            hs("IL");
            hs(t.name);
            for (double y : FA)
                hx(xsect::invLookup(y, t.t, t.n, L));
            nl();
        }
        hs("LOC");
        hs(t.name);
        for (double y : FA)
            hi(xsect::locate(y, t.t, t.n - 1));
        nl();
    }
    {
        const auto& E = xsect::hostEval();
        hs("YC");
        for (double x : FA)
            hx(xsect::getYcircular(x));
        nl();
        hs("SC");
        for (double x : FA)
            hx(xsect::getScircular(x));
        nl();
        hs("AC");
        for (double x : FA)
            hx(E.getAcircular(x));
        nl();
        hs("TA");
        for (double x : FA)
            hx(E.getThetaOfAlpha(x));
        nl();
        hs("TP");
        for (double x : FA)
            hx(E.getThetaOfPsi(x));
        nl();
    }

    // ---- transect LUTs
    std::vector<xsect::LocateLut> tluts(nT);
    for (int t = 0; t < nT; ++t) {
        xsect::build_invlookup_lut(tluts[t], Tv[t].ar.data(), N);
        hs("TLUT");
        hi(t);
        dumpLut(tluts[t]);
        nl();
    }

    // ---- per-element cases
    std::vector<XSectParams> net;
    std::vector<int> curve;
    for (size_t c = 0; c < Pv.size(); ++c) {
        const PRec& P = Pv[c];
        XSectParams xs;
        int ret = xsect::setParams(xs, P.type, P.p, P.ucf);
        if (isTabulated(P.type) && P.tid >= 0 && P.tid < nT) {
            const TRec& T = Tv[P.tid];
            xs.y_full = T.yf;
            xs.w_max = T.wm;
            xs.a_full = T.af;
            xs.r_full = T.rf;
            xs.s_full = T.sf;
            xs.s_max = T.sm;
            xs.a_bot = T.ab;
            xs.area_tbl = T.ar.data();
            xs.hrad_tbl = T.hr.data();
            xs.width_tbl = T.wt.data();
            xs.transect_tbl_size = N;
            xs.area_lut = &tluts[P.tid];
        }
        hs("C");
        hi((long long)c);
        hi(ret);
        dumpXs(xs);
        nl();

        const double Yb = base_of(xs.y_full), Ab = base_of(xs.a_full), Sb = base_of(xs.s_full);
        double qb = xs.a_full * std::sqrt(xsect::GRAVITY * xs.y_full);
        const double Qb = base_of(qb);
        std::vector<double> D, AA, SS, QQ;
        for (double f : Fv) {
            D.push_back(f * Yb);
            AA.push_back(f * Ab);
            SS.push_back(f * Sb);
            QQ.push_back(f * Qb);
        }
        D.insert(D.end(), Av.begin(), Av.end());
        AA.insert(AA.end(), Av.begin(), Av.end());
        SS.insert(SS.end(), Av.begin(), Av.end());
        QQ.insert(QQ.end(), Qv.begin(), Qv.end());

        hs("AofY");
        for (double y : D)
            hx(xsect::getAofY(xs, y));
        nl();
        hs("WofY");
        for (double y : D)
            hx(xsect::getWofY(xs, y));
        nl();
        hs("RofY");
        for (double y : D)
            ok_RofY(xs, y) ? hx(xsect::getRofY(xs, y)) : dash();
        nl();
        hs("YofA");
        for (double a : AA)
            hx(xsect::getYofA(xs, a));
        nl();
        hs("SofA");
        for (double a : AA)
            ok_SofA(xs, a) ? hx(xsect::getSofA(xs, a)) : dash();
        nl();
        hs("RofA");
        for (double a : AA)
            ok_RofA(xs, a) ? hx(xsect::getRofA(xs, a)) : dash();
        nl();
        hs("dSdA");
        for (double a : AA)
            ok_dSdA(xs, a) ? hx(xsect::getdSdA(xs, a)) : dash();
        nl();
        hs("AofS");
        for (double s : SS)
            ok_AofS(xs, s) ? hx(xsect::getAofS(xs, s)) : dash();
        nl();
        hs("Ycrit");
        for (double q : QQ)
            ok_Ycrit(xs, q) ? hx(xsect::getYcrit(xs, q)) : dash();
        nl();
        hs("Amax");
        hx(xsect::getAmax(xs));
        hi(xsect::isOpen(P.type) ? 1 : 0);
        nl();

        if (batch_ok(xs)) {
            net.push_back(xs);
            curve.push_back(isTabulated(P.type) ? P.tid : -1);
        }
    }
    for (int t : {-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                  21, 22, 23, 24, 25, 26, 100}) {
        hs("open");
        hi(t);
        hi(xsect::isOpen(t) ? 1 : 0);
        nl();
    }

    // ---- batch: network of the batch-safe cases
    const int n = (int)net.size();
    SimulationContext ctx;
    ctx.links.resize(n);
    for (int i = 0; i < n; ++i)
        ctx.links.xsect_curve[(size_t)i] = curve[(size_t)i];
    for (int t = 0; t < nT; ++t) {
        transect::TransectData td;
        td.y_full = Tv[t].yf;
        td.w_max = Tv[t].wm;
        td.a_full = Tv[t].af;
        td.r_full = Tv[t].rf;
        for (int j = 0; j < N; ++j) {
            td.area_tbl[j] = Tv[t].ar[j];
            td.hrad_tbl[j] = Tv[t].hr[j];
            td.width_tbl[j] = Tv[t].wt[j];
        }
        ctx.transect_tables.push_back(td);
    }
    std::vector<double> by(n), ba(n);
    for (int i = 0; i < n; ++i) {
        by[i] = base_of(net[i].y_full);
        ba[i] = base_of(net[i].a_full);
    }
    std::vector<double> Fb;
    for (size_t i = 0; i < Fv.size(); ++i) {
        double f = Fv[i];
        bool keep = std::isnan(f) || (std::isfinite(f) && std::fabs(f) <= 100.0);
        if (keep && (i % 3 == 0 || i < 21))
            Fb.push_back(f);
    }
    const int nv = (int)Fb.size();
    auto dvec = [&](int vi) {
        std::vector<double> d(n);
        for (int i = 0; i < n; ++i)
            d[i] = Fb[(size_t)(vi % nv)] * by[i];
        return d;
    };
    auto avec = [&](int vi) {
        std::vector<double> a(n);
        for (int i = 0; i < n; ++i)
            a[i] = Fb[(size_t)(vi % nv)] * ba[i];
        return a;
    };
    std::vector<std::vector<uint8_t>> masks(4, std::vector<uint8_t>(n));
    for (int i = 0; i < n; ++i) {
        masks[0][i] = (i % 3 == 0);
        masks[1][i] = 1;
        masks[2][i] = 0;
        masks[3][i] = ((i * 7) % 5 < 2);
    }

    hs("NET");
    hi(n);
    hi(nv);
    nl();
    XSectGroups G;
    G.build(net.data(), n);
    dumpGroups("build", G);
    {   // before attach: per-link groups have no table vectors and are left unwritten
        auto d = dvec(5);
        std::vector<double> o1(n, SENT), o2(n, SENT), o3(n, SENT);
        G.computeAreas(d.data(), o1.data(), n);
        G.computeHydRad(d.data(), o2.data(), n);
        G.computeWidths(d.data(), o3.data(), n);
        prline("preA", 5, o1);
        prline("preR", 5, o2);
        prline("preW", 5, o3);
    }
    G.attachTransectTables(ctx);
    dumpGroups("attach", G);
    G.attachTransectTables(ctx);
    dumpGroups("attach2", G);

    for (int vi = 0; vi < nv; ++vi) {
        auto d = dvec(vi), d2 = dvec(vi + 1), dm = dvec(vi + 2), a = avec(vi);
        std::vector<double> o1(n, SENT), o2(n, SENT);
        G.computeAreas(d.data(), o1.data(), n);
        prline("bA", vi, o1);
        std::fill(o1.begin(), o1.end(), SENT);
        G.computeHydRad(d.data(), o1.data(), n);
        prline("bR", vi, o1);
        std::fill(o1.begin(), o1.end(), SENT);
        G.computeWidths(d.data(), o1.data(), n);
        prline("bW", vi, o1);
        std::fill(o1.begin(), o1.end(), SENT);
        G.computeAreaAndHydRad(d.data(), o1.data(), o2.data(), n);
        prline("bAR_a", vi, o1);
        prline("bAR_r", vi, o2);
        std::fill(o1.begin(), o1.end(), SENT);
        G.computeSectionFactors(a.data(), o1.data(), n);
        prline("bSF", vi, o1);
        std::fill(o1.begin(), o1.end(), SENT);
        G.computeDepthsFromArea(a.data(), o1.data(), n);
        prline("bDA", vi, o1);

        auto triple = [&](const char* tag) {
            std::vector<double> a1(n, SENT), a2(n, SENT), am(n, SENT), h1(n, SENT), hm(n, SENT);
            G.computeAreaHydRadTriple(d.data(), d2.data(), dm.data(), a1.data(), a2.data(),
                                      am.data(), h1.data(), hm.data(), n);
            std::vector<double> w1(n, SENT), w2(n, SENT), wm(n, SENT);
            G.computeWidthsTriple(d.data(), d2.data(), dm.data(), w1.data(), w2.data(),
                                  wm.data(), n);
            std::string t(tag);
            prline((t + "_a1").c_str(), vi, a1);
            prline((t + "_a2").c_str(), vi, a2);
            prline((t + "_am").c_str(), vi, am);
            prline((t + "_h1").c_str(), vi, h1);
            prline((t + "_hm").c_str(), vi, hm);
            prline((t + "_w1").c_str(), vi, w1);
            prline((t + "_w2").c_str(), vi, w2);
            prline((t + "_wm").c_str(), vi, wm);
        };
        triple("T");
        if (vi % 5 == 0) {
            for (int m = 0; m < 4; ++m) {
                G.setBypassMask(masks[m].data());
                char tag[16];
                std::snprintf(tag, sizeof tag, "M%d", m);
                triple(tag);
                std::fill(o1.begin(), o1.end(), SENT);
                G.computeAreas(d.data(), o1.data(), n); // unmasked
                prline((std::string(tag) + "_A").c_str(), vi, o1);
            }
            G.setBypassMask(nullptr);
            triple("M_off");
        }
        if (vi % 7 == 0) {
            for (int nt : {2, 3, 5, 0, -3}) {
                for (int withMask : {0, 1}) {
                    if (withMask)
                        G.setBypassMask(masks[3].data());
                    std::vector<double> a1(n, SENT), a2(n, SENT), am(n, SENT), h1(n, SENT),
                        hm(n, SENT), w1(n, SENT), w2(n, SENT), wm(n, SENT);
                    int workers = nt < 1 ? 1 : nt;
                    for (int tid = 0; tid < workers; ++tid) {
                        G.computeAreaHydRadTripleTeam(d.data(), d2.data(), dm.data(), a1.data(),
                                                      a2.data(), am.data(), h1.data(), hm.data(),
                                                      n, tid, nt);
                        G.computeWidthsTripleTeam(d.data(), d2.data(), dm.data(), w1.data(),
                                                  w2.data(), wm.data(), n, tid, nt);
                    }
                    if (withMask)
                        G.setBypassMask(nullptr);
                    char tag[24];
                    std::snprintf(tag, sizeof tag, "TM%d_%d", nt, withMask);
                    std::string t(tag);
                    prline((t + "_a1").c_str(), vi, a1);
                    prline((t + "_a2").c_str(), vi, a2);
                    prline((t + "_am").c_str(), vi, am);
                    prline((t + "_h1").c_str(), vi, h1);
                    prline((t + "_hm").c_str(), vi, hm);
                    prline((t + "_w1").c_str(), vi, w1);
                    prline((t + "_w2").c_str(), vi, w2);
                    prline((t + "_wm").c_str(), vi, wm);
                }
            }
        }
    }
    {
        XSectGroups G2;
        G2.build(net.data(), n / 2);
        dumpGroups("half", G2);
        G2.setBypassMask(masks[0].data());
        G2.build(net.data(), 0);
        dumpGroups("zero", G2);
        G2.build(nullptr, -5);
        dumpGroups("neg", G2);
    }

    // ---- every batch kernel over the network's columns
    std::vector<double> cy(n), ca(n), cr(n), cs(n), cw(n), cyb(n), csb(n), crb(n);
    for (int i = 0; i < n; ++i) {
        cy[i] = net[i].y_full;
        ca[i] = net[i].a_full;
        cr[i] = net[i].r_full;
        cs[i] = net[i].s_full;
        cw[i] = net[i].w_max;
        cyb[i] = net[i].y_bot;
        csb[i] = net[i].s_bot;
        crb[i] = net[i].r_bot;
    }
    std::vector<const double*> ta(n), tr(n);
    for (int i = 0; i < n; ++i) {
        int k = i % 5;
        ta[i] = k < nT ? Tv[k].ar.data() : nullptr;
        int k2 = (i + 1) % 5;
        tr[i] = k2 < nT ? Tv[k2].hr.data() : nullptr;
    }
    const auto& HT = xsect::hostTables();
    struct NT {
        const char* name;
        const double* t;
        int n;
    };
    const NT Atab[] = {{"A_Circ", HT.A_Circ, HT.N_A_Circ},           {"A_Egg", HT.A_Egg, HT.N_A_Egg},
                       {"A_Horseshoe", HT.A_Horseshoe, HT.N_A_Horseshoe},
                       {"A_Baskethandle", HT.A_Baskethandle, HT.N_A_Baskethandle},
                       {"A_HorizEllipse", HT.A_HorizEllipse, HT.N_A_HorizEllipse},
                       {"A_VertEllipse", HT.A_VertEllipse, HT.N_A_VertEllipse},
                       {"A_Arch", HT.A_Arch, HT.N_A_Arch}};
    const NT Ytab[] = {{"Y_Gothic", HT.Y_Gothic, HT.N_Y_Gothic},
                       {"Y_Catenary", HT.Y_Catenary, HT.N_Y_Catenary},
                       {"Y_SemiEllip", HT.Y_SemiEllip, HT.N_Y_SemiEllip},
                       {"Y_SemiCirc", HT.Y_SemiCirc, HT.N_Y_SemiCirc}};
    const int Ylut[] = {0, 1, 2, 3};
    const NT Rtab[] = {{"R_Egg", HT.R_Egg, HT.N_R_Egg},
                       {"R_Horseshoe", HT.R_Horseshoe, HT.N_R_Horseshoe},
                       {"R_Baskethandle", HT.R_Baskethandle, HT.N_R_Baskethandle},
                       {"R_HorizEllipse", HT.R_HorizEllipse, HT.N_R_HorizEllipse},
                       {"R_VertEllipse", HT.R_VertEllipse, HT.N_R_VertEllipse},
                       {"R_Arch", HT.R_Arch, HT.N_R_Arch}};
    const NT Wtab[] = {{"W_Egg", HT.W_Egg, HT.N_W_Egg},
                       {"W_Horseshoe", HT.W_Horseshoe, HT.N_W_Horseshoe},
                       {"W_Gothic", HT.W_Gothic, HT.N_W_Gothic},
                       {"W_Catenary", HT.W_Catenary, HT.N_W_Catenary},
                       {"W_SemiEllip", HT.W_SemiEllip, HT.N_W_SemiEllip},
                       {"W_BasketHandle", HT.W_BasketHandle, HT.N_W_BasketHandle},
                       {"W_SemiCirc", HT.W_SemiCirc, HT.N_W_SemiCirc},
                       {"W_HorizEllipse", HT.W_HorizEllipse, HT.N_W_HorizEllipse},
                       {"W_VertEllipse", HT.W_VertEllipse, HT.N_W_VertEllipse},
                       {"W_Arch", HT.W_Arch, HT.N_W_Arch}};
    for (int vi = 0; vi < nv; vi += 3) {
        auto d = dvec(vi);
        std::vector<double> o(n, SENT), o2(n, SENT);
        auto K = [&](const char* tag) {
            prline(tag, vi, o);
            std::fill(o.begin(), o.end(), SENT);
        };
        xb::area_circular(d.data(), cy.data(), ca.data(), o.data(), n);
        K("k_area_circular");
        xb::area_rect(d.data(), cw.data(), o.data(), n);
        K("k_area_rect");
        xb::area_trapezoidal(d.data(), cyb.data(), csb.data(), o.data(), n);
        K("k_area_trapezoidal");
        xb::area_triangular(d.data(), csb.data(), o.data(), n);
        K("k_area_triangular");
        xb::area_parabolic(d.data(), crb.data(), o.data(), n);
        K("k_area_parabolic");
        xb::area_powerfunc(d.data(), csb.data(), crb.data(), o.data(), n);
        K("k_area_powerfunc");
        for (const auto& t : Atab) {
            xb::area_tabulated(d.data(), cy.data(), ca.data(), t.t, t.n, o.data(), n);
            K((std::string("k_area_tabulated_") + t.name).c_str());
        }
        for (int j = 0; j < 4; ++j) {
            const auto& t = Ytab[j];
            xb::area_inv_tabulated(d.data(), cy.data(), ca.data(), t.t, t.n, nullptr, o.data(), n);
            K((std::string("k_area_inv_") + t.name).c_str());
            xb::area_inv_tabulated(d.data(), cy.data(), ca.data(), t.t, t.n,
                                   HT.lut(static_cast<xsect::LutId>(Ylut[j])), o.data(), n);
            K((std::string("k_area_inv_lut_") + t.name).c_str());
        }
        xb::perlink_tabulated(d.data(), cy.data(), ca.data(), ta.data(), N, o.data(), n);
        K("k_perlink");
        xb::perlink_tabulated_pair(d.data(), cy.data(), ca.data(), cr.data(), ta.data(), tr.data(),
                                   N, o.data(), o2.data(), n);
        prline("k_pair_b", vi, o2);
        std::fill(o2.begin(), o2.end(), SENT);
        K("k_pair_a");
        xb::hydrad_circular(d.data(), cy.data(), cr.data(), o.data(), n);
        K("k_hydrad_circular");
        xb::area_hydrad_circular(d.data(), cy.data(), ca.data(), cr.data(), o.data(), o2.data(), n);
        prline("k_ahc_r", vi, o2);
        std::fill(o2.begin(), o2.end(), SENT);
        K("k_ahc_a");
        xb::hydrad_trapezoidal(d.data(), cyb.data(), csb.data(), crb.data(), o.data(), n);
        K("k_hydrad_trapezoidal");
        xb::hydrad_triangular(d.data(), csb.data(), crb.data(), o.data(), n);
        K("k_hydrad_triangular");
        xb::hydrad_rect(d.data(), cw.data(), o.data(), n);
        K("k_hydrad_rect");
        xb::hydrad_rect_closed(d.data(), cw.data(), ca.data(), o.data(), n);
        K("k_hydrad_rect_closed");
        xb::hydrad_rect_open(d.data(), cw.data(), csb.data(), o.data(), n);
        K("k_hydrad_rect_open");
        for (const auto& t : Rtab) {
            xb::hydrad_tabulated(d.data(), cy.data(), cr.data(), t.t, t.n, o.data(), n);
            K((std::string("k_hydrad_tabulated_") + t.name).c_str());
        }
        xb::width_circular(d.data(), cy.data(), cw.data(), o.data(), n);
        K("k_width_circular");
        xb::width_trapezoidal(d.data(), cyb.data(), csb.data(), o.data(), n);
        K("k_width_trapezoidal");
        xb::width_triangular(d.data(), csb.data(), o.data(), n);
        K("k_width_triangular");
        xb::width_rect(cw.data(), o.data(), n);
        K("k_width_rect");
        xb::width_rect_closed(d.data(), cy.data(), cw.data(), o.data(), n);
        K("k_width_rect_closed");
        for (const auto& t : Wtab) {
            xb::width_tabulated(d.data(), cy.data(), cw.data(), t.t, t.n, o.data(), n);
            K((std::string("k_width_tabulated_") + t.name).c_str());
        }
        // count 0 and negative: no writes
        xb::area_circular(d.data(), cy.data(), ca.data(), o.data(), 0);
        xb::perlink_tabulated(d.data(), cy.data(), ca.data(), ta.data(), N, o.data(), -4);
        K("k_count0");
    }

    std::fwrite(OUT.data(), 1, OUT.size(), stdout);
    return 0;
}

int main(int argc, char** argv) {
    if (argc == 3 && std::string(argv[1]) == "gen")
        return gen(argv[2]);
    if (argc == 2)
        return run(argv[1]);
    std::fprintf(stderr, "usage: xsect gen FILE | xsect FILE\n");
    return 2;
}
