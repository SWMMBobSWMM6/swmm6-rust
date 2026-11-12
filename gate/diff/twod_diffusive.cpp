// Differential harness for src/engine/2d/solver/DiffusiveKernels.hpp (MOMENTUM_EQUATION
// DIFFUSIVE_WAVE): faceDischarge and cellDiffusiveDt (both glibc cbrt), over grids of their
// branch points (dry, zero slope, the slope floor s_eps and its neighbours, both flow directions)
// and seeded random draws. Prints every input and result as hex bits; the Rust twin is
// crates/engine/examples/twod_diffusive.rs.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
//
// Build: as gate/diff/twod_vfr.cpp (gate/twodkern.sh), with twod_diffusive in place of twod_vfr.
#include "2d/solver/DiffusiveKernels.hpp"
#include "twodgen.hpp"

using namespace openswmm::twoD::diffusive;
using namespace twodgen;
using V = std::vector<double>;

static void args(Out& o, std::initializer_list<double> a) {
    for (double v : a) o.hx(v);
    o.sep();
}

int main() {
    Out o;
    Gen g(0x2d0f1a7e5eed0005ull);
    const int N = 60000;
    const V HF = {0.0, -0.0, 1e-12, 1e-3, 0.3, 1.0, 4.0, -1.0, DMIN, 1e-200};
    const V SL = {0.0, -0.0, 1e-6, -1e-6, std::nextafter(1e-6, 0.0), std::nextafter(1e-6, INF), 1e-3, -1e-3, 0.2, -0.2};
    const V NF = {0.0, -0.0, 0.03, 0.1, -0.02, 1e-9};
    const V SE = {1e-6, 0.0, 1e-4, 1e-9};
    grid({HF.size(), SL.size(), NF.size(), SE.size()}, [&](const std::vector<std::size_t>& i) {
        const double a = opq(HF[i[0]]), b = opq(SL[i[1]]), c = opq(NF[i[2]]), d = opq(SE[i[3]]);
        o.tag("fq"); args(o, {a, b, c, d}); o.hx(faceDischarge(a, b, c, d)); o.nl();
    });
    for (int k = 0; k < N; ++k) {
        const double a = g.draw(HF, 0.0, 3.0), b = g.draw(SL, -0.05, 0.05), c = g.draw(NF, 0.01, 0.1), d = g.draw(SE, 0.0, 1e-4);
        o.tag("fq"); args(o, {a, b, c, d}); o.hx(faceDischarge(opq(a), opq(b), opq(c), opq(d))); o.nl();
    }
    const V AL = {0.25, 0.5, 1.0, 0.0};
    const V LC = {1.0, 10.0, 0.0, 1e-3};
    const V SA = {0.0, 1e-6, std::nextafter(1e-6, 0.0), std::nextafter(1e-6, INF), 1e-3, 0.2, -0.1};
    grid({AL.size(), LC.size(), HF.size(), NF.size(), SA.size(), SE.size()}, [&](const std::vector<std::size_t>& i) {
        const double a = opq(AL[i[0]]), b = opq(LC[i[1]]), c = opq(HF[i[2]]), d = opq(NF[i[3]]), e = opq(SA[i[4]]),
                     f = opq(SE[i[5]]);
        o.tag("dt"); args(o, {a, b, c, d, e, f}); o.hx(cellDiffusiveDt(a, b, c, d, e, f)); o.nl();
    });
    for (int k = 0; k < N; ++k) {
        const double a = g.draw(AL, 0.1, 1.0), b = g.draw(LC, 0.1, 100.0), c = g.draw(HF, 0.0, 3.0);
        const double d = g.draw(NF, 0.01, 0.1), e = g.draw(SA, 0.0, 0.1), f = g.draw(SE, 0.0, 1e-4);
        o.tag("dt"); args(o, {a, b, c, d, e, f});
        o.hx(cellDiffusiveDt(opq(a), opq(b), opq(c), opq(d), opq(e), opq(f)));
        o.nl();
    }
    o.flush();
    std::fprintf(stderr, "twod_diffusive: %ld lines\n", o.lines);
    return 0;
}
