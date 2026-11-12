// Differential harness for src/engine/2d/coupling/PondModel.hpp (the above-rim pond of a
// 2D-coupled junction): area, volume, depthForVolume (glibc cbrt and a safeguarded Newton) with
// the volume -> depth round trip, and equalizeCap, over grids of their branch points (x at and
// around the band, band zero and negative for the legacy step, nonpositive volumes and areas,
// unknown windows) and seeded random draws. Prints every input and result as hex bits; the Rust
// twin is crates/engine/examples/twod_pond.rs.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
//
// Build: as gate/diff/twod_vfr.cpp (gate/twodkern.sh), with twod_pond in place of twod_vfr.
#include "2d/coupling/PondModel.hpp"
#include "twodgen.hpp"

using namespace openswmm::twoD::pond;
using namespace twodgen;
using V = std::vector<double>;

static void args(Out& o, std::initializer_list<double> a) {
    for (double v : a) o.hx(v);
    o.sep();
}

int main() {
    Out o;
    Gen g(0x2d0f1a7e5eed0009ull);
    const int N = 60000;
    const V XX = {0.0, -0.0, -1.0, 0.05, std::nextafter(0.05, 0.0), std::nextafter(0.05, INF), 0.025, 1e-12, 1.0, DMIN, 0.1};
    const V AQ = {10.0, 0.0, -0.0, -5.0, 1.0, 1e6, 1e-9};
    const V BD = {0.05, 0.0, -0.0, -0.05, 0.1, 1e-9, 1.0};
    grid({XX.size(), AQ.size(), BD.size()}, [&](const std::vector<std::size_t>& i) {
        const double x = opq(XX[i[0]]), a = opq(AQ[i[1]]), b = opq(BD[i[2]]);
        o.tag("pa"); args(o, {x, a, b}); o.hx(area(x, a, b)); o.hx(volume(x, a, b)); o.nl();
        o.tag("pd"); args(o, {x, a, b}); o.hx(depthForVolume(x, a, b)); o.nl(); // x as a volume
    });
    for (int k = 0; k < N; ++k) {
        const double a = g.draw(AQ, 0.0, 100.0), b = g.draw(BD, 0.0, 0.2);
        const double x = g.below(2) == 0 ? g.draw(XX, -0.01, 0.3) : canon(b * g.uni(0.0, 1.2));
        o.tag("pa"); args(o, {x, a, b});
        const double ar = area(opq(x), opq(a), opq(b)), vo = volume(opq(x), opq(a), opq(b));
        o.hx(ar); o.hx(vo);
        o.hx(depthForVolume(opq(vo), opq(a), opq(b))); // round trip
        o.nl();
        const double v = g.below(2) == 0 ? g.draw(XX, -0.01, 3.0) : canon(0.5 * a * b * g.uni(0.0, 1.1));
        o.tag("pd"); args(o, {v, a, b}); o.hx(depthForVolume(opq(v), opq(a), opq(b))); o.nl();
    }
    const V DH = {0.0, -0.0, 1.0, -1.0, 1e-12, -0.3};
    const V WN = {0.0, -0.0, 1.0, 30.0, -1.0, 1e-9};
    const V A1 = {0.0, -0.0, 1.0, 12.5, -2.0, 1e-9};
    grid({DH.size(), WN.size(), A1.size()}, [&](const std::vector<std::size_t>& i) {
        const double a = opq(DH[i[0]]), b = opq(WN[i[1]]), c = opq(A1[i[2]]);
        o.tag("eq"); args(o, {a, b, c}); o.hx(equalizeCap(a, b, c)); o.nl();
    });
    for (int k = 0; k < N; ++k) {
        const double a = g.draw(DH, -2.0, 2.0), b = g.draw(WN, 0.0, 60.0), c = g.draw(A1, 0.0, 50.0);
        o.tag("eq"); args(o, {a, b, c}); o.hx(equalizeCap(opq(a), opq(b), opq(c))); o.nl();
    }
    o.flush();
    std::fprintf(stderr, "twod_pond: %ld lines\n", o.lines);
    return 0;
}
