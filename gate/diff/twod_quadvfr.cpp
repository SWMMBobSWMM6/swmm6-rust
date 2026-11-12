// Differential harness for src/engine/2d/mesh/QuadVfr.hpp (the two-plane quad VFR):
// quadVfrPrecompute over every 4-vertex elevation pattern from {0,1,2} (all three diagonal cases,
// ties) on square, reversed and skewed quads plus degenerate (collinear) and random quads; then
// quadWetFraction, quadMeanDepthFromEta, quadMeanDepthSlope, quadDEtaDMeanDepth and
// quadEtaFromMeanDepth (with the eta -> depth -> eta round trip) on each, at stages and depths
// around every vertex elevation and branch point. Prints every input and result as hex bits; the
// Rust twin is crates/engine/examples/twod_quadvfr.rs.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
//
// Build: as gate/diff/twod_vfr.cpp (gate/twodkern.sh), with twod_quadvfr in place of twod_vfr.
#include "2d/mesh/QuadVfr.hpp"
#include "twodgen.hpp"

#include <algorithm>
#include <utility>

using namespace openswmm::twoD;
using namespace twodgen;

static const std::vector<double> EPS = {0.0, -0.0, 1e-6, 1e-3, 0.01, 0.05, 0.2, 0.5, 0.99, 1.0, -0.1};

static void quad_funcs(Out& o, Gen& g, const double* zs_in, double A1, double A2, int nstage) {
    double zs[6];
    for (int k = 0; k < 6; ++k) zs[k] = opq(zs_in[k]);
    A1 = opq(A1);
    A2 = opq(A2);
    std::vector<double> st = {-INF, INF};
    for (int k = 0; k < 6; ++k) {
        st.push_back(zs[k]);
        st.push_back(std::nextafter(zs[k], INF));
        st.push_back(std::nextafter(zs[k], -INF));
        st.push_back(zs[k] + 0.25);
        st.push_back(zs[k] - 0.25);
    }
    const double zw = canon((A1 * (zs[0] + zs[1] + zs[2]) / 3.0 + A2 * (zs[3] + zs[4] + zs[5]) / 3.0) / (A1 + A2));
    st.push_back(zw);
    for (int r = 0; r < nstage; ++r) {
        const double eta = opq(r < static_cast<int>(st.size()) ? st[r] : g.uni(zs[0] - 1.0, zs[5] + 1.0));
        const double eps = opq(g.pick(EPS));
        auto pre = [&](const char* tag) {
            o.tag(tag);
            for (int k = 0; k < 6; ++k) o.hx(zs[k]);
            o.hx(A1);
            o.hx(A2);
        };
        pre("qwf"); o.hx(eta); o.sep(); o.hx(quadWetFraction(zs, A1, A2, eta)); o.nl();
        pre("qmd"); o.hx(eta); o.hx(eps);
        const double d = quadMeanDepthFromEta(zs, A1, A2, eta, eps);
        o.sep();
        o.hx(d);
        o.hx(quadEtaFromMeanDepth(zs, A1, A2, opq(d), eps));
        o.nl();
        pre("qsl"); o.hx(eta); o.hx(eps); o.sep(); o.hx(quadMeanDepthSlope(zs, A1, A2, eta, eps)); o.nl();
        pre("qde"); o.hx(eta); o.hx(eps); o.sep(); o.hx(quadDEtaDMeanDepth(zs, A1, A2, eta, eps)); o.nl();
        // the inverse at depths around its branch points
        const double ztop = zs[2] > zs[5] ? zs[2] : zs[5];
        const double tw = canon(ztop - zw);
        const double dsel[] = {0.0, -0.0, -1.0, DMIN, 1e-9, 1e-3, tw, std::nextafter(tw, -INF),
                               std::nextafter(tw, INF), g.uni(0.0, 3.0), INF};
        const double md = opq(dsel[g.below(11)]);
        pre("qef"); o.hx(md); o.hx(eps); o.sep(); o.hx(quadEtaFromMeanDepth(zs, A1, A2, md, eps)); o.nl();
    }
}

static void precompute(Out& o, Gen& g, const double* x, const double* y, const double* z, int nstage) {
    double xs[4], ys[4], zz[4], zs[6], A1 = -7.0, A2 = -7.0;
    for (int k = 0; k < 4; ++k) {
        xs[k] = opq(x[k]);
        ys[k] = opq(y[k]);
        zz[k] = opq(z[k]);
    }
    const int kase = quadVfrPrecompute(xs, ys, zz, zs, A1, A2);
    o.tag("pre");
    for (int k = 0; k < 4; ++k) { o.hx(xs[k]); o.hx(ys[k]); o.hx(zz[k]); }
    o.sep();
    for (int k = 0; k < 6; ++k) o.hx(zs[k]);
    o.hx(A1); o.hx(A2); o.hi(kase); o.nl();
    quad_funcs(o, g, zs, A1, A2, nstage);
}

int main() {
    Out o;
    Gen g(0x2d0f1a7e5eed0002ull);
    const double SQX[4] = {0, 1, 1, 0}, SQY[4] = {0, 0, 1, 1};
    const double RVX[4] = {0, 0, 1, 1}, RVY[4] = {0, 1, 1, 0};
    const double SKX[4] = {0, 2.5, 3.1, 0.4}, SKY[4] = {0, 0.2, 1.7, 1.3};
    const double LNX[4] = {0, 1, 2, 3}, LNY[4] = {0, 1, 2, 3.0000001}; // nearly collinear
    const double* XS[4] = {SQX, RVX, SKX, LNX};
    const double* YS[4] = {SQY, RVY, SKY, LNY};
    const std::vector<double> ZV = {0.0, 1.0, 2.0};
    for (int q = 0; q < 4; ++q)
        grid({3, 3, 3, 3}, [&](const std::vector<std::size_t>& i) {
            const double z[4] = {ZV[i[0]], ZV[i[1]], ZV[i[2]], ZV[i[3]]};
            precompute(o, g, XS[q], YS[q], z, 40);
        });
    // distinct elevations in every cyclic order (all 24 permutations of 4 levels)
    const double L4[4] = {10.0, 10.5, 11.75, 13.0};
    int perm[4] = {0, 1, 2, 3};
    do {
        const double z[4] = {L4[perm[0]], L4[perm[1]], L4[perm[2]], L4[perm[3]]};
        for (int q = 0; q < 3; ++q) precompute(o, g, XS[q], YS[q], z, 40);
    } while (std::next_permutation(perm, perm + 4));
    // random quads: jittered rectangles at a datum, random elevations (some ties, some flat)
    const std::vector<double> SPZ = {0.0, 1.0, 1e-9, 1e-10, 100.0};
    for (int k = 0; k < 6000; ++k) {
        const double w = g.uni(0.1, 50.0), h = g.uni(0.1, 50.0), x0 = g.uni(-1e3, 1e3), y0 = g.uni(-1e3, 1e3);
        double x[4] = {x0, x0 + w, x0 + w, x0}, y[4] = {y0, y0, y0 + h, y0 + h};
        for (int v = 0; v < 4; ++v) {
            x[v] += g.uni(-0.2, 0.2) * w;
            y[v] += g.uni(-0.2, 0.2) * h;
        }
        const double base = g.uni(-10.0, 500.0);
        double z[4];
        const int mode = static_cast<int>(g.below(4));
        for (int v = 0; v < 4; ++v) {
            switch (mode) {
            case 0: z[v] = base + g.uni(0.0, 3.0); break;
            case 1: z[v] = base + g.pick(SPZ); break;
            case 2: z[v] = g.draw(SPZ, -5.0, 5.0); break;
            default: z[v] = base + g.uni(0.0, 1e-8); break;
            }
        }
        if (g.below(2) == 0) { // reverse orientation
            std::swap(x[1], x[3]);
            std::swap(y[1], y[3]);
            std::swap(z[1], z[3]);
        }
        precompute(o, g, x, y, z, 8);
    }
    // raw closure data not from the precompute: zero sub-areas, unsorted triples, wild values
    for (int k = 0; k < 3000; ++k) {
        double zs[6];
        for (int v = 0; v < 6; ++v) zs[v] = g.draw(SPZ, -5.0, 5.0);
        if (g.below(4) != 0) {
            vfrSort3(zs[0], zs[1], zs[2]);
            vfrSort3(zs[3], zs[4], zs[5]);
        }
        const std::vector<double> SA = {0.0, -0.0, 1.0, 0.5, 1e-12, 1e6};
        const double A1 = g.draw(SA, 0.0, 100.0), A2 = g.draw(SA, 0.0, 100.0);
        quad_funcs(o, g, zs, A1, A2, 4);
    }
    o.flush();
    std::fprintf(stderr, "twod_quadvfr: %ld lines\n", o.lines);
    return 0;
}
