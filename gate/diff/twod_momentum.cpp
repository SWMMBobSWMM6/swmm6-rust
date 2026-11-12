// Differential harness for src/engine/2d/solver/MomentumPhysics2D: the class's member defaults
// and MomentumPhysics2D::rotate (std::cos / std::sin), called in the 2D reference library itself
// (build2d's libopenswmm.engine.so exports it), so the comparison is against upstream's own
// compiled code. Each call's vectors start with the unit pairs (1, 0) and (0, 1), so cos and
// sin of every angle are visible directly; angles cover zero, signed zero, multiples and
// fractions of pi, tiny, huge (argument reduction) and non-finite values and seeded random ones.
// Prints every input and result as hex bits; the Rust twin is
// crates/engine/examples/twod_momentum.rs.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
//
// Build: as gate/diff/twod_vfr.cpp (gate/twodkern.sh), with twod_momentum in place of twod_vfr,
// plus /root/work/swmm6src/build2d/src/engine/libopenswmm.engine.so
// -Wl,-rpath,/root/work/swmm6src/build2d/src/engine.
#include "2d/solver/MomentumPhysics2D.hpp"
#include "twodgen.hpp"

using namespace openswmm::twoD;
using namespace twodgen;
using V = std::vector<double>;

int main() {
    Out o;
    Gen g(0x2d0f1a7e5eed0007ull);
    {
        MomentumPhysics2D mp;
        o.tag("def");
        o.sep();
        o.hx(mp.last_limit); o.hi(static_cast<long long>(mp.substeps)); o.hi(static_cast<long long>(mp.retries));
        o.hx(MomentumPhysics2D::kappa);
        o.hi(static_cast<long long>(mp.viscosity.size() + mp.strain.size() + mp.vorticity.size() + mp.wall_distance.size()));
        o.nl();
    }
    const double PI = 3.141592653589793;
    const V AN = {0.0, -0.0, PI, -PI, PI / 2, PI / 4, 2 * PI, 1e-8, -1e-8, 1e-300, DMIN, 1.0, -1.0, 7.292e-5,
                  1e6, 1e22, -1e22, 1.7e308, INF, -INF, QNAN, 0.5, 100.0, 3.0e-3};
    const V QV = {0.0, -0.0, 1.0, -1.0, 2.5, 1e-12, -3.0};
    auto call = [&](double angle, int n) {
        V qx = {1.0, 0.0}, qy = {0.0, 1.0};
        for (int k = 0; k < n; ++k) {
            qx.push_back(g.draw(QV, -5.0, 5.0));
            qy.push_back(g.draw(QV, -5.0, 5.0));
        }
        o.tag("rot");
        o.hx(angle);
        for (std::size_t k = 0; k < qx.size(); ++k) { o.hx(qx[k]); o.hx(qy[k]); }
        o.sep();
        for (std::size_t k = 0; k < qx.size(); ++k) { qx[k] = opq(qx[k]); qy[k] = opq(qy[k]); }
        MomentumPhysics2D::rotate(qx, qy, opq(angle));
        for (std::size_t k = 0; k < qx.size(); ++k) { o.hx(qx[k]); o.hx(qy[k]); }
        o.nl();
    };
    for (double a : AN) call(a, 6);
    for (int k = 0; k < 60000; ++k) {
        double a;
        switch (g.below(5)) {
        case 0: a = g.uni(-2 * PI, 2 * PI); break;
        case 1: a = g.ebits(1023 - 30, 1023 + 3); break;   // small to ~16
        case 2: a = g.ebits(1, 2046); break;                // anywhere (argument reduction)
        case 3: a = g.nudge(g.pick(AN)); break;
        default: a = g.uni(-1e-3, 1e-3); break;              // Coriolis-sized f*dt
        }
        call(a, static_cast<int>(g.below(9)));
    }
    // empty vectors and mismatched lengths are not called (upstream indexes qy by qx's size)
    o.flush();
    std::fprintf(stderr, "twod_momentum: %ld lines\n", o.lines);
    return 0;
}
