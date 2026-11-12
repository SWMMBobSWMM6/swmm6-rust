// Differential harness for src/engine/2d/coupling/ExchangeLaw.hpp (COUPLING_LAW WEIR_ORIFICE):
// phi, dphi, smoothstep, softmin and exchange. phi/dphi run over a grid around the
// regularization head eps; softmin over ratios from 1 down to underflow of r^16 and s^17;
// exchange over random structures (crest, weir length, orifice area, coefficients including the
// defaults) and head pairs in every regime: surface -> sewer free (h1 <= z), inside the
// submergence band, fully submerged, sewer -> surface onto a dry and a wet surface, equal heads,
// heads at the crest. Prints every input and result as hex bits; the Rust twin is
// crates/engine/examples/twod_exchange.rs.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
//
// Build: as gate/diff/twod_vfr.cpp (gate/twodkern.sh), with twod_exchange in place of twod_vfr.
#include "2d/coupling/ExchangeLaw.hpp"
#include "twodgen.hpp"

using namespace openswmm::twoD::exchange;
using namespace twodgen;
using V = std::vector<double>;

static void args(Out& o, std::initializer_list<double> a) {
    for (double v : a) o.hx(v);
    o.sep();
}

int main() {
    Out o;
    Gen g(0x2d0f1a7e5eed0008ull);
    const int N = 60000;
    // ---- phi, dphi
    const V XS = {0.0, -0.0, -1.0, 0.02, std::nextafter(0.02, 0.0), std::nextafter(0.02, INF), 0.01, 1e-12, 1.0, 4.0, DMIN, 1e-3};
    const V ES = {0.02, 0.0, 1e-3, 0.1, 1e-9, -0.02, 1.0};
    grid({XS.size(), ES.size()}, [&](const std::vector<std::size_t>& i) {
        const double x = opq(XS[i[0]]), e = opq(ES[i[1]]);
        o.tag("phi"); args(o, {x, e}); o.hx(phi(x, e)); o.hx(dphi(x, e)); o.nl();
    });
    for (int k = 0; k < N; ++k) {
        const double x = g.draw(XS, -0.05, 0.5), e = g.draw(ES, 1e-4, 0.1);
        o.tag("phi"); args(o, {x, e}); o.hx(phi(opq(x), opq(e))); o.hx(dphi(opq(x), opq(e))); o.nl();
    }
    // ---- smoothstep
    const V TS = {0.0, -0.0, 1.0, 0.5, -0.5, 1.5, 1e-12, std::nextafter(1.0, 0.0), std::nextafter(1.0, INF), DMIN};
    for (double t : TS) { o.tag("ss"); args(o, {t}); o.hx(smoothstep(opq(t))); o.nl(); }
    for (int k = 0; k < N; ++k) {
        const double t = g.draw(TS, -0.2, 1.2);
        o.tag("ss"); args(o, {t}); o.hx(smoothstep(opq(t))); o.nl();
    }
    // ---- softmin
    const V AB = {0.0, -0.0, -1.0, 1.0, 2.0, 0.5, 1e-3, 1e-20, 1e20, 1e-300, 1e300, DMIN, 0.9999999};
    grid({AB.size(), AB.size()}, [&](const std::vector<std::size_t>& i) {
        const double a = opq(AB[i[0]]), b = opq(AB[i[1]]);
        double wa = -7, wb = -7;
        const double q = softmin(a, b, wa, wb);
        o.tag("sm"); args(o, {a, b}); o.hx(q); o.hx(wa); o.hx(wb); o.nl();
    });
    for (int k = 0; k < N; ++k) {
        const double a = g.draw(AB, 0.0, 2.0);
        const double b = g.below(3) == 0 ? canon(a * g.uni(0.5, 2.0)) : g.draw(AB, 0.0, 2.0);
        double wa = -7, wb = -7;
        const double q = softmin(opq(a), opq(b), wa, wb);
        o.tag("sm"); args(o, {a, b}); o.hx(q); o.hx(wa); o.hx(wb); o.nl();
    }
    // ---- exchange(geo, c, h2, h1)
    const V ZC = {0.0, 10.0, -0.0, 100.25};
    const V LW = {0.0, 1.0, 3.14, 0.5};
    const V AO = {0.0, 0.1, 0.785, 1e-4};
    const V CO = {0.54, 0.36, 0.167, 0.16, 0.0, 1.0};
    const V HE = {0.02, 1e-3, 0.1, 1e-9};
    for (int k = 0; k < 3 * N; ++k) {
        Geometry geo;
        geo.z_crest = g.draw(ZC, -5.0, 50.0);
        geo.l_weir = g.draw(LW, 0.0, 5.0);
        geo.a_orif = g.draw(AO, 0.0, 2.0);
        Coefficients c;
        if (g.below(3) != 0) {
            c.c_w = g.draw(CO, 0.0, 1.0);
            c.c_sw = g.draw(CO, 0.0, 1.0);
            c.c_oi = g.draw(CO, 0.0, 1.0);
            c.c_oo = g.draw(CO, 0.0, 1.0);
            c.h_eps = g.draw(HE, 1e-4, 0.1);
        }
        const double z = geo.z_crest;
        // heads relative to the crest, in every regime
        const V RH = {0.0, -0.0, 0.01, 0.02, 0.05, 0.1, 0.5, 1.0, -0.5, 1e-12, -1e-12, 2.0};
        double h2, h1;
        switch (g.below(6)) {
        case 0: h2 = canon(z + g.pick(RH)); h1 = canon(z + g.pick(RH)); break;
        case 1: h2 = canon(z + g.uni(0.0, 2.0)); h1 = canon(z + g.uni(-1.0, 2.0)); break;
        case 2: { const double hs = g.uni(1e-3, 1.0); h2 = canon(z + hs); h1 = canon(z + hs * g.uni(0.0, 0.12)); } break; // band
        case 3: h2 = canon(z + g.uni(-1.0, 2.0)); h1 = h2; break;
        case 4: h2 = g.draw(RH, -1.0, 60.0); h1 = g.draw(RH, -1.0, 60.0); break;
        default: { const double hd = g.uni(1e-3, 1.0); h1 = canon(z + hd); h2 = canon(z + hd * g.uni(-0.5, 0.15)); } break;
        }
        h2 = opq(h2);
        h1 = opq(h1);
        Geometry gq{opq(geo.z_crest), opq(geo.l_weir), opq(geo.a_orif)};
        Coefficients cq{opq(c.c_w), opq(c.c_sw), opq(c.c_oi), opq(c.c_oo), opq(c.h_eps)};
        const Result r = exchange(gq, cq, h2, h1);
        o.tag("ex");
        args(o, {geo.z_crest, geo.l_weir, geo.a_orif, c.c_w, c.c_sw, c.c_oi, c.c_oo, c.h_eps, h2, h1});
        o.hx(r.q); o.hx(r.g); o.nl();
    }
    o.flush();
    std::fprintf(stderr, "twod_exchange: %ld lines\n", o.lines);
    return 0;
}
