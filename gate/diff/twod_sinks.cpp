// Differential harness for the inline sink ramps of src/engine/2d/solver/SurfaceFluxCalculator.hpp:
// evapSink and infilSink, over a grid of their branch points (nonpositive rate or depth, depth at
// and around dry_depth) and seeded random draws. Prints every input and result as hex bits; the
// Rust twin is crates/engine/examples/twod_sinks.rs.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
//
// Build: as gate/diff/twod_vfr.cpp (gate/twodkern.sh), with twod_sinks in place of twod_vfr.
#include "2d/solver/SurfaceFluxCalculator.hpp"
#include "twodgen.hpp"

using namespace openswmm::twoD;
using namespace twodgen;
using V = std::vector<double>;

static void args(Out& o, std::initializer_list<double> a) {
    for (double v : a) o.hx(v);
    o.sep();
}

int main() {
    Out o;
    Gen g(0x2d0f1a7e5eed0006ull);
    const int N = 50000;
    const V RT = {0.0, -0.0, 1e-7, 2.7e-8, -1e-7, 1.0, DMIN};
    const V DP = {0.0, -0.0, 1e-4, std::nextafter(1e-4, 0.0), std::nextafter(1e-4, INF), 5e-5, 1e-12, 1.0, -1.0, DMIN};
    const V DD = {1e-4, 0.0, -0.0, 1e-3, 1e-12, -1e-4};
    grid({RT.size(), DP.size(), DD.size()}, [&](const std::vector<std::size_t>& i) {
        const double a = opq(RT[i[0]]), b = opq(DP[i[1]]), c = opq(DD[i[2]]);
        o.tag("ev"); args(o, {a, b, c}); o.hx(evapSink(a, b, c)); o.nl();
        o.tag("in"); args(o, {a, b, c}); o.hx(infilSink(a, b, c)); o.nl();
    });
    for (int k = 0; k < N; ++k) {
        const double a = g.draw(RT, 0.0, 1e-5), b = g.draw(DP, 0.0, 2e-4), c = g.draw(DD, 0.0, 1e-3);
        o.tag("ev"); args(o, {a, b, c}); o.hx(evapSink(opq(a), opq(b), opq(c))); o.nl();
        o.tag("in"); args(o, {a, b, c}); o.hx(infilSink(opq(a), opq(b), opq(c))); o.nl();
    }
    o.flush();
    std::fprintf(stderr, "twod_sinks: %ld lines\n", o.lines);
    return 0;
}
