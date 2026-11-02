// Differential harness for hydraulics/XSectLookup.hpp, with Transect.hpp's TransectData
// defaults. Upstream's own headers are compiled; nothing is copied.
//
//   xslookup gen FILE   build the tables and queries and write them to FILE: the 50 legacy
//                       tables of xsect_tables.hpp, then generated ones of 1..300 entries
//                       (increasing, with duplicates, flat, decreasing, random, with a NaN,
//                       with a declining top, dyadic). Queries cover 0, -0, negatives, every
//                       grid point and table value and their neighbours, midpoints, values
//                       above range, NaN and infinities, and random draws.
//   xslookup FILE       call every inline function on every table and query; print each
//                       result's bits in hex (ints in decimal), one line per query.
//
// Inputs on which upstream is undefined (static_cast<int> of a value that does not fit
// an int, or a table read past its end) are not called; the line prints "-" there, and
// the Rust side applies the same predicate. lookup_fast is not called: the translation
// covers the default (SWMM_XSECT_FAST_LOOKUP off) build only.
//
// Build: g++ -std=c++20 -O2 -ffp-contract=off -I<upstream>/src/engine/hydraulics
#include "XSectLookup.hpp"
#include "Transect.hpp"
#include "xsect_tables.hpp"
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

using namespace openswmm;
using xsect::LocateLut;

static void hx(double v) {
    uint64_t u;
    std::memcpy(&u, &v, 8);
    std::printf("%016llx ", (unsigned long long)u);
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

struct Case {
    std::vector<double> t, q;
};

static void add_queries(Case& c) {
    const std::vector<double>& t = c.t;
    const int n = (int)t.size();
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<double>& q = c.q;
    const double fixed[] = {0.0,   -0.0,   1.0,    std::nextafter(1.0, 0.0), std::nextafter(1.0, 2.0),
                            1.5,   2.0,    10.0,   1e6,   1e300, -1e-300, -1e-12, -0.001, -0.5,
                            -1.0,  -1e6,   NaN,    inf,   -inf,  0.25,   0.5,    0.75,   1e-300,
                            5e-324, 0.01, 0.02,   0.04,  0.06};
    for (double v : fixed) q.push_back(v);
    if (n >= 2) {
        const double delta = 1.0 / (double)(n - 1);
        q.push_back(-0.5 * delta);
        q.push_back(-0.999 * delta);
        for (int i = 0; i < n; ++i) {
            double g1 = (double)i / (double)(n - 1), g2 = i * delta;
            q.push_back(g1);
            q.push_back(g2);
            q.push_back(std::nextafter(g1, -inf));
            q.push_back(std::nextafter(g1, inf));
            q.push_back((i + 0.5) * delta);
        }
    }
    // y-type queries: the table's own values, their neighbours, midpoints, beyond range
    double lo = inf, hi = -inf;
    for (int j = 0; j < n; ++j) {
        double v = t[j];
        q.push_back(v);
        q.push_back(std::nextafter(v, -inf));
        q.push_back(std::nextafter(v, inf));
        if (j + 1 < n) q.push_back((t[j] + t[j + 1]) / 2.0);
        if (v == v) {
            if (v < lo) lo = v;
            if (v > hi) hi = v;
        }
    }
    if (lo <= hi) {
        q.push_back(lo - 1.0);
        q.push_back(hi + 1.0);
        double span = hi - lo;
        for (int k = 0; k < 30; ++k) q.push_back(lo - 0.1 * span + u01() * 1.2 * span);
    }
    for (int k = 0; k < 30; ++k) q.push_back(u01());
    for (int k = 0; k < 20; ++k) q.push_back(-0.1 + u01() * 1.3);
}

static void gen(const char* path) {
    std::vector<Case> cs;
    using namespace openswmm::xsect_tables;
#define T(name) do { Case c; c.t.assign(name, name + sizeof(name) / sizeof(double)); cs.push_back(c); } while (0)
    T(Amax); T(A_Circ); T(R_Circ); T(Y_Circ); T(S_Circ); T(W_Circ); T(A_Egg); T(R_Egg); T(Y_Egg);
    T(S_Egg); T(W_Egg); T(A_Horseshoe); T(R_Horseshoe); T(Y_Horseshoe); T(S_Horseshoe);
    T(W_Horseshoe); T(Y_Gothic); T(S_Gothic); T(W_Gothic); T(Y_Catenary); T(S_Catenary);
    T(W_Catenary); T(Y_SemiEllip); T(S_SemiEllip); T(W_SemiEllip); T(A_Baskethandle);
    T(R_Baskethandle); T(Y_BasketHandle); T(S_BasketHandle); T(W_BasketHandle); T(Y_SemiCirc);
    T(S_SemiCirc); T(W_SemiCirc); T(A_HorizEllipse); T(R_HorizEllipse); T(W_HorizEllipse);
    T(A_VertEllipse); T(R_VertEllipse); T(W_VertEllipse); T(A_Arch); T(R_Arch); T(W_Arch);
    T(MinorAxis_Ellipse); T(MajorAxis_Ellipse); T(Afull_Ellipse); T(Rfull_Ellipse);
    T(Yfull_Arch); T(Wmax_Arch); T(Afull_Arch); T(Rfull_Arch);
#undef T
    const int sizes[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 20, 21, 26, 51, 64, 100, 200, 254, 255, 256, 300};
    for (int n : sizes) {
        for (int kind = 0; kind < 10; ++kind) {
            Case c;
            std::vector<double>& t = c.t;
            t.resize(n);
            switch (kind) {
            case 0: // strictly increasing from 0
                t[0] = 0.0;
                for (int j = 1; j < n; ++j) t[j] = t[j - 1] + 0.001 + u01();
                break;
            case 1: // nondecreasing with duplicates, on a 1/1024 grid
                t[0] = 0.0;
                for (int j = 1; j < n; ++j) t[j] = t[j - 1] + (u01() < 0.4 ? 0.0 : (double)(1 + (rnd() % 64)) / 1024.0);
                break;
            case 2: // nondecreasing from a negative offset, duplicates
                t[0] = -5.0 * u01();
                for (int j = 1; j < n; ++j) t[j] = t[j - 1] + (u01() < 0.3 ? 0.0 : u01());
                break;
            case 3: // flat
                for (int j = 0; j < n; ++j) t[j] = 0.75;
                break;
            case 4: // random, non-monotone
                for (int j = 0; j < n; ++j) t[j] = -1.0 + 2.0 * u01();
                break;
            case 5: // decreasing
                t[0] = 1.0;
                for (int j = 1; j < n; ++j) t[j] = t[j - 1] - 0.5 * u01();
                break;
            case 6: // increasing with one NaN
                t[0] = 0.0;
                for (int j = 1; j < n; ++j) t[j] = t[j - 1] + u01();
                t[rnd() % n] = std::numeric_limits<double>::quiet_NaN();
                break;
            case 7: // increasing, declining top two rows (section-factor style)
                t[0] = 0.0;
                for (int j = 1; j < n; ++j) t[j] = t[j - 1] + 0.01 + u01();
                if (n >= 3) {
                    t[n - 2] = t[n - 3] + 0.25 * u01();
                    t[n - 1] = t[n - 3] - 0.01 - 0.25 * u01();
                }
                break;
            case 8: // dyadic: t[j] = j/(n-1) by the same divide the lookup uses
                for (int j = 0; j < n; ++j) t[j] = n > 1 ? (double)j / (double)(n - 1) : 0.0;
                break;
            case 9: // geometry-like: 0..1 increasing, convex, plus a long duplicate run at the top
                for (int j = 0; j < n; ++j) {
                    double x = n > 1 ? (double)j / (double)(n - 1) : 0.0;
                    t[j] = x * x;
                }
                for (int j = (3 * n) / 4; j < n; ++j) t[j] = t[(3 * n) / 4];
                break;
            }
            cs.push_back(c);
        }
    }
    for (Case& c : cs) add_queries(c);
    FILE* f = std::fopen(path, "wb");
    int nt = (int)cs.size();
    std::fwrite(&nt, 4, 1, f);
    for (Case& c : cs) {
        int n = (int)c.t.size(), nq = (int)c.q.size();
        std::fwrite(&n, 4, 1, f);
        std::fwrite(c.t.data(), 8, n, f);
        std::fwrite(&nq, 4, 1, f);
        std::fwrite(c.q.data(), 8, nq, f);
    }
    std::fclose(f);
}

// ---------------------------------------------------------------- the run
// True when upstream's lookup_exact(x, t, n) is defined: x/delta converts to an int,
// and every table read stays inside the n entries.
static bool ok_lookup(double x, int n) {
    const double delta = 1.0 / static_cast<double>(n - 1);
    const double q = x / delta;
    if (!(q > -1.0 && q < 2147483647.0)) return false;
    int i = static_cast<int>(q);
    if (i >= n - 1) return true;
    if (i < 2 && i + 2 > n - 1) return false;
    return true;
}
// True when lut_bucket's truncating cast is defined.
static bool ok_bucket(const LocateLut& L, double y) {
    const double v = (y - L.t0) * L.scale;
    return v > -2147483648.0 && v < 2147483647.0;
}
static void dump_lut(const LocateLut& L) {
    hx(L.t0);
    hx(L.scale);
    std::printf("%d ", L.j_last);
    for (int b = 0; b <= LocateLut::kBuckets; ++b) std::printf("%02x", L.lo[b]);
    std::printf(" ");
}

static void run(const char* path) {
    // Transect.hpp: TransectData's member initializers
    {
        transect::TransectData td;
        std::printf("TD %d %zu ", transect::N_TRANSECT_TBL, td.name.size());
        hx(td.n_left); hx(td.n_right); hx(td.n_channel); hx(td.length_factor);
        std::printf("%zu %zu ", td.stations.size(), td.elevations.size());
        hx(td.x_left_bank); hx(td.x_right_bank); hx(td.y_full); hx(td.w_max); hx(td.a_full);
        hx(td.r_full); hx(td.s_max); hx(td.a_max);
        for (int j = 0; j < transect::N_TRANSECT_TBL; ++j) { hx(td.area_tbl[j]); hx(td.hrad_tbl[j]); hx(td.width_tbl[j]); }
        dump_lut(td.area_lut);
        std::printf("%d\n", LocateLut::kBuckets);
    }
    const double NaN = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const double params[] = {0.0, -1.0, NaN, 0.37, 1.0, 2.5, inf};

    FILE* f = std::fopen(path, "rb");
    int nt;
    if (std::fread(&nt, 4, 1, f) != 1) return;
    long total = 0;
    for (int tix = 0; tix < nt; ++tix) {
        int n, nq;
        if (std::fread(&n, 4, 1, f) != 1) return;
        std::vector<double> t(n);
        if (std::fread(t.data(), 8, n, f) != (size_t)n) return;
        if (std::fread(&nq, 4, 1, f) != 1) return;
        std::vector<double> q(nq);
        if (std::fread(q.data(), 8, nq, f) != (size_t)nq) return;

        // jLast values: 0, 1, 2, n/2, n-1 (each once, within the table)
        std::vector<int> jl;
        for (int v : {0, 1, 2, n / 2, n - 1})
            if (v >= 0 && v <= n - 1) {
                bool dup = false;
                for (int w : jl) dup = dup || (w == v);
                if (!dup) jl.push_back(v);
            }
        std::vector<LocateLut> L(jl.size());
        std::printf("T %d %d ", tix, n);
        for (size_t k = 0; k < jl.size(); ++k) {
            xsect::build_locate_lut(L[k], t.data(), jl[k]);
            dump_lut(L[k]);
        }
        LocateLut Lnull;
        Lnull.scale = 9.0; // build must reset it
        xsect::build_locate_lut(Lnull, nullptr, n - 1);
        dump_lut(Lnull);
        LocateLut Linv, Linv_null;
        xsect::build_invlookup_lut(Linv, t.data(), n);
        Linv_null.t0 = 3.0;
        xsect::build_invlookup_lut(Linv_null, nullptr, n);
        dump_lut(Linv);
        dump_lut(Linv_null);
        int neff = n;
        if (n >= 4 && t[n - 3] > t[n - 1]) neff = n - 2;
        std::printf("%d\n", neff);
        const LocateLut& Lfull = L.back(); // built for jLast == n-1

        for (double x : q) {
            hx(x);
            if (ok_lookup(x, n)) hx(xsect::lookup_exact(x, t.data(), n)); else std::printf("- ");
            for (double p : params) {
                double xn = xsect::norm_x(x, p);
                hx(xn);
                if (ok_lookup(xn, n)) hx(xsect::norm_lookup(x, p, t.data(), n)); else std::printf("- ");
            }
            std::printf("| ");
            for (size_t k = 0; k < jl.size(); ++k) {
                const int j = jl[k];
                std::printf("%d %d %d %d ", xsect::locate_bisect(x, t.data(), j),
                            xsect::locate_maybe_lut(x, t.data(), j, &L[k]),
                            xsect::locate_maybe_lut(x, t.data(), j, nullptr),
                            xsect::locate_maybe_lut(x, t.data(), j, &Lfull));
                if (L[k].scale > 0.0) std::printf("%d ", xsect::locate_lut(x, t.data(), j, L[k])); else std::printf("- ");
                if (ok_bucket(L[k], x)) std::printf("%d ", xsect::lut_bucket(L[k], x)); else std::printf("- ");
            }
            if (n >= 4) std::printf("%d ", xsect::locate_maybe_lut(x, t.data(), neff - 1, &Linv));
            std::printf("\n");
            ++total;
        }
    }
    std::fclose(f);
    std::fprintf(stderr, "%d tables, %ld queries\n", nt, total);
}

int main(int argc, char** argv) {
    if (argc == 3 && std::strcmp(argv[1], "gen") == 0) {
        gen(argv[2]);
        return 0;
    }
    if (argc == 2) {
        run(argv[1]);
        return 0;
    }
    std::fprintf(stderr, "usage: xslookup gen FILE | xslookup FILE\n");
    return 2;
}
