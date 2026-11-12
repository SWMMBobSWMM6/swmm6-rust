// Differential harness for src/engine/2d/mesh/VfrClosure.hpp (the planar-bed triangle VFR):
// vfrSort3, vfrWetFraction, vfrMeanDepthFromEtaExact, vfrStageAtWetFraction,
// vfrMeanDepthFromEta, vfrEtaFromMeanDepth, vfrDryEta, vfrDEtaDMeanDepth and the
// eta -> depth -> eta round trip, over structured bed triples (flat, near-flat at the 1e-9
// relief guard, coincident vertices, large datums) and seeded random ones, with stages and depths
// at and around every branch point (z1, z2, z3, zbar, the eps switch stage, h_at_z2, z3 - zbar).
// Prints every input and result as hex bits; the Rust twin is crates/engine/examples/twod_vfr.rs.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
//
// Build (the 2D reference build's flags and defines, from build2d's flags.make; gate/twodkern.sh):
//   U=/root/work/swmm6src/upstream; B=/root/work/swmm6src/build2d
//   g++ -ffp-contract=off -O3 -DNDEBUG -std=c++20 -flto=auto -fno-fat-lto-objects -fPIC \
//       -fvisibility=default -mavx2 -mfma -fopenmp -DOPENSWMM_HAS_2D=1 -DSWMM_OS_LINUX \
//       -DSWMM_USE_OPENMP -I$U/include -I$U/src/engine -I$B/include -Igate/diff \
//       gate/diff/twod_vfr.cpp -o gate/diff/twod_vfr
#include "2d/mesh/VfrClosure.hpp"
#include "twodgen.hpp"

using namespace openswmm::twoD;
using namespace twodgen;

struct Tri { double a, b, c; };

static std::vector<Tri> structured() {
    std::vector<Tri> t = {
        {0, 0, 0}, {1, 1, 1}, {-0.0, 0.0, 0.0}, {0, 0, 1}, {0, 1, 1}, {0, 0.5, 1}, {0, 0.25, 1},
        {0, 1e-10, 2e-10}, {0, 5e-10, 1e-9}, {0, 1e-9, 1}, {0, 9.99e-10, 1}, {0, 1 - 1e-9, 1},
        {0, 1 - 5e-10, 1}, {100.0, 100.3, 101.7}, {-5, -2, 3}, {1000.0, 1000.0 + 1e-6, 1000.0 + 3e-6},
        {12345.678, 12345.679, 12345.68}, {0, 2, 3}, {-1e-3, 0, 1e-3}, {0, 1e-6, 1e-3},
        {0.1, 0.2, 0.3}, {1, 2, 1000}, {-0.0, 1e-300, 2e-300}, {0, 1e-9, 2e-9},
        {5, 5, 5 + 1e-9}, {5, 5 + 1e-9, 5 + 1e-9}, {DMIN, 2 * DMIN, 3 * DMIN},
        // unsorted (outside the documented precondition; still deterministic)
        {3, 1, 2}, {1, 0, 0.5},
    };
    return t;
}

static std::vector<double> stages(Tri z) {
    z = {opq(z.a), opq(z.b), opq(z.c)}; // no folding of the derived values on either side
    const double zbar = canon((z.a + z.b + z.c) / 3.0);
    std::vector<double> e = {z.a, z.b, z.c, zbar, canon(0.5 * (z.a + z.b)), canon(0.5 * (z.b + z.c)), z.a - 1.0,
                             z.c + 1.0, z.a - 1e-12, z.c + 1e-12, -INF, INF};
    for (double v : {z.a, z.b, z.c}) {
        e.push_back(std::nextafter(v, INF));
        e.push_back(std::nextafter(v, -INF));
        e.push_back(std::nextafter(std::nextafter(v, INF), INF));
    }
    return e;
}

static const std::vector<double> EPS = {0.0, -0.0, 1e-6, 1e-4, 0.01, 0.05, 0.1, 0.3, 0.9, 1.0, -0.01, 0.999999};

static void line_tri(Out& o, const char* tag, const Tri& z) {
    o.tag(tag);
    o.hx(z.a);
    o.hx(z.b);
    o.hx(z.c);
}

// all of one triple's functions at one (eta, eps), and the depths derived from them
static void one(Out& o, Tri z, double eta, double eps) {
    const double z1 = opq(z.a), z2 = opq(z.b), z3 = opq(z.c);
    eta = opq(eta);
    eps = opq(eps);
    line_tri(o, "wf", z); o.hx(eta); o.sep(); o.hx(vfrWetFraction(z1, z2, z3, eta)); o.nl();
    line_tri(o, "mdx", z); o.hx(eta); o.sep(); o.hx(vfrMeanDepthFromEtaExact(z1, z2, z3, eta)); o.nl();
    line_tri(o, "md", z); o.hx(eta); o.hx(eps);
    const double d = vfrMeanDepthFromEta(z1, z2, z3, eta, eps);
    o.sep();
    o.hx(d);
    // round trip eta -> depth -> eta
    o.hx(vfrEtaFromMeanDepth(z1, z2, z3, opq(d), eps));
    o.nl();
    line_tri(o, "de", z); o.hx(eta); o.hx(eps); o.sep(); o.hx(vfrDEtaDMeanDepth(z1, z2, z3, eta, eps)); o.nl();
}

static void per_eps(Out& o, Tri z, double eps) {
    const double z1 = opq(z.a), z2 = opq(z.b), z3 = opq(z.c);
    eps = opq(eps);
    line_tri(o, "st", z); o.hx(eps); o.sep(); o.hx(vfrStageAtWetFraction(z1, z2, z3, eps)); o.nl();
    line_tri(o, "dry", z); o.hx(eps); o.sep(); o.hx(vfrDryEta(z1, z2, z3, eps)); o.nl();
    // depths at every branch point of the inverse
    z = {opq(z.a), opq(z.b), opq(z.c)};
    const double zbar = (z.a + z.b + z.c) / 3.0, relief = z.c - z.a, d21 = z.b - z.a;
    std::vector<double> ds = {0.0, -0.0, -1.0, DMIN, 1e-12, 1e-6, 0.5, z.c - zbar, d21 * d21 / (3.0 * relief),
                              1e3, INF, -INF};
    if (eps > 0.0 && relief >= kVfrFlatRelief) {
        const double es = vfrStageAtWetFraction(z.a, z.b, z.c, eps);
        ds.push_back(vfrMeanDepthFromEtaExact(z.a, z.b, z.c, es));
    }
    const std::size_t n0 = ds.size();
    for (double& d : ds) d = canon(d);
    for (std::size_t i = 4; i < n0; ++i) {
        ds.push_back(std::nextafter(ds[i], INF));
        ds.push_back(std::nextafter(ds[i], -INF));
    }
    for (double d : ds) {
        d = opq(d);
        line_tri(o, "ef", z); o.hx(d); o.hx(eps); o.sep(); o.hx(vfrEtaFromMeanDepth(z1, z2, z3, d, eps)); o.nl();
    }
}

int main() {
    Out o;
    Gen g(0x2d0f1a7e5eed0001ull);
    // vfrSort3 over a grid of triples with ties and signed zeros
    const std::vector<double> SZ = {-1.0, 0.0, -0.0, 1e-10, 1.0, 2.0, 1.0 + 1e-12, 5.0, 100.0};
    grid({SZ.size(), SZ.size(), SZ.size()}, [&](const std::vector<std::size_t>& i) {
        double a = opq(SZ[i[0]]), b = opq(SZ[i[1]]), c = opq(SZ[i[2]]);
        o.tag("sort"); o.hx(a); o.hx(b); o.hx(c); o.sep();
        vfrSort3(a, b, c);
        o.hx(a); o.hx(b); o.hx(c); o.nl();
    });
    // structured triples: every stage x every eps
    for (Tri z : structured()) {
        for (double eps : EPS) {
            per_eps(o, z, eps);
            for (double eta : stages(z)) one(o, z, eta, eps);
        }
    }
    // random triples (sorted, as the solver calls them), a few stages and eps each
    const std::vector<double> SPZ = {0.0, -0.0, 1.0, 1e-9, 1e-10, 100.0, -3.0};
    for (int k = 0; k < 6000; ++k) {
        Tri z;
        switch (g.below(4)) {
        case 0: z = {g.uni(-10, 100), g.uni(-10, 100), g.uni(-10, 100)}; break;
        case 1: { const double b = g.uni(-50, 500); z = {b, b + g.pos_ebits(1023 - 35, 1023 + 2), b + g.pos_ebits(1023 - 35, 1023 + 2)}; } break;
        case 2: z = {g.draw(SPZ, -10, 10), g.draw(SPZ, -10, 10), g.draw(SPZ, -10, 10)}; break;
        default: { const double b = g.uni(0, 2000); z = {b, b + g.uni(0, 1e-8), b + g.uni(0, 3)}; } break;
        }
        if (g.below(16) != 0) vfrSort3(z.a, z.b, z.c);
        const std::vector<double> st = stages(z);
        for (int r = 0; r < 4; ++r) {
            double eps = g.below(3) == 0 ? g.pick(EPS) : g.uni(0, 0.5);
            double eta;
            switch (g.below(3)) {
            case 0: eta = g.pick(st); break;
            case 1: eta = g.uni(z.a - 1.0, z.c + 1.0); break;
            default: eta = g.nudge(g.pick(st)); break;
            }
            if (r == 0) per_eps(o, z, eps);
            one(o, z, eta, eps);
        }
    }
    o.flush();
    std::fprintf(stderr, "twod_vfr: %ld lines\n", o.lines);
    return 0;
}
