// Differential harness for src/engine/2d/solver/SweKernels.hpp (MOMENTUM_EQUATION FULL_SWE):
// physicalFlux, waveSpeeds, hllcFlux, faceFlux, faceFluxReconBed, faceFluxRecon, bjLimiter,
// frictionUpdate (glibc cbrt), criticalDepth (glibc cbrt) and depthFromInvariantAndDischarge.
// Face states cover wet-wet, wet-dry, dry-wet and dry-dry sides, rest states (lake at rest over
// a step), supercritical flow both ways (every HLLC branch), the 1e-14 contact-denominator guard
// and the h_dry velocity cut-off; normals come from a fixed list and from random unit vectors.
// Prints every input and result as hex bits; the Rust twin is crates/engine/examples/twod_swe.rs.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
//
// Build: as gate/diff/twod_vfr.cpp (gate/twodkern.sh), with twod_swe in place of twod_vfr.
#include "2d/solver/SweKernels.hpp"
#include "twodgen.hpp"

using namespace openswmm::twoD::swe;
using namespace twodgen;
using V = std::vector<double>;

static void args(Out& o, std::initializer_list<double> a) {
    for (double v : a) o.hx(v);
    o.sep();
}

static const V HS = {0.0, -0.0, 1e-12, 1e-6, 0.5, 1.0, 2.0, -1.0, DMIN};
static const V US = {0.0, -0.0, 1.0, -1.0, 3.0, -3.0, 10.0, -10.0, 1e-12};
static const V NXS = {1.0, 0.0, -1.0, -0.0, 0.6, -0.8, 0.7071067811865476, -0.7071067811865476};
static const V NYS = {0.0, 1.0, -0.0, -1.0, 0.8, 0.6, 0.7071067811865476, 0.7071067811865476};

static FaceSide side(Gen& g) {
    FaceSide s;
    s.h = g.draw(HS, 0.0, 3.0);
    s.un = g.draw(US, -8.0, 8.0);
    s.ut = g.draw(US, -8.0, 8.0);
    return s;
}

static void normal(Gen& g, double& nx, double& ny) {
    if (g.below(2) == 0) {
        const std::size_t k = g.below(NXS.size());
        nx = NXS[k];
        ny = NYS[k];
    } else {
        nx = g.uni(-1.0, 1.0);
        ny = std::sqrt(1.0 - nx * nx);
        if (g.below(2) == 0) ny = -ny;
    }
}

static void face_flux_line(Out& o, const char* tag, std::initializer_list<double> in, const FaceFlux& f, double a,
                           double b, double c, double d, bool r) {
    o.tag(tag);
    args(o, in);
    o.hx(f.mass); o.hx(f.mx); o.hx(f.my); o.hx(f.sstar);
    o.hx(a); o.hx(b); o.hx(c); o.hx(d); o.hi(r ? 1 : 0);
    o.nl();
}

int main() {
    Out o;
    Gen g(0x2d0f1a7e5eed0004ull);
    const int N = 40000;

    // ---- physicalFlux, waveSpeeds, hllcFlux over a grid of sides, then random sides
    const V HG = {0.0, 1e-12, 0.5, 2.0, -1.0};
    const V UG = {0.0, -0.0, 1.0, -1.0, 10.0, -10.0};
    const V TG = {0.0, 2.0};
    grid({HG.size(), UG.size(), TG.size(), HG.size(), UG.size(), TG.size()}, [&](const std::vector<std::size_t>& i) {
        FaceSide L, R;
        L.h = opq(HG[i[0]]); L.un = opq(UG[i[1]]); L.ut = opq(TG[i[2]]);
        R.h = opq(HG[i[3]]); R.un = opq(UG[i[4]]); R.ut = opq(TG[i[5]]);
        double sl = -7, sr = -7, fh = -7, fn = -7, ft = -7, ss = -7;
        waveSpeeds(L, R, sl, sr);
        hllcFlux(L, R, fh, fn, ft, ss);
        o.tag("hllc"); args(o, {L.h, L.un, L.ut, R.h, R.un, R.ut});
        o.hx(sl); o.hx(sr); o.hx(fh); o.hx(fn); o.hx(ft); o.hx(ss); o.nl();
    });
    for (int k = 0; k < 2 * N; ++k) {
        FaceSide L = side(g), R = side(g);
        if (g.below(8) == 0) R.h = L.h; // equal depths: the contact denominator degenerates
        L.h = opq(L.h); L.un = opq(L.un); L.ut = opq(L.ut);
        R.h = opq(R.h); R.un = opq(R.un); R.ut = opq(R.ut);
        double pf = -7, pn = -7, pt = -7, sl = -7, sr = -7, fh = -7, fn = -7, ft = -7, ss = -7;
        physicalFlux(L, pf, pn, pt);
        waveSpeeds(L, R, sl, sr);
        hllcFlux(L, R, fh, fn, ft, ss);
        o.tag("hllc"); args(o, {L.h, L.un, L.ut, R.h, R.un, R.ut});
        o.hx(sl); o.hx(sr); o.hx(fh); o.hx(fn); o.hx(ft); o.hx(ss);
        o.hx(pf); o.hx(pn); o.hx(pt);
        o.nl();
    }

    // ---- faceFlux(etaL, hL, qxL, qyL, etaR, hR, qxR, qyR, nx, ny, h_dry, out, corr...)
    const V ETA = {0.0, 1.0, 1.5, 2.0, -0.0, 1e-12};
    const V QX = {0.0, -0.0, 0.5, -0.5, 3.0, -3.0};
    const V HD = {1e-6, 0.0, 1e-3, 1e-12};
    // lake at rest over steps and dry beds: a grid
    grid({ETA.size(), HS.size(), ETA.size(), HS.size(), 4}, [&](const std::vector<std::size_t>& i) {
        const double etaL = opq(ETA[i[0]]), hL = opq(HS[i[1]]), etaR = opq(ETA[i[2]]), hR = opq(HS[i[3]]);
        const double nx = opq(NXS[i[4]]), ny = opq(NYS[i[4]]), hd = opq(1e-6);
        FaceFlux f{-7, -7, -7, -7};
        double a = -7, b = -7, c = -7, d = -7;
        const bool r = faceFlux(etaL, hL, opq(0.0), opq(0.0), etaR, hR, opq(0.0), opq(0.0), nx, ny, hd, f, a, b, c, d);
        face_flux_line(o, "ff", {etaL, hL, 0.0, 0.0, etaR, hR, 0.0, 0.0, nx, ny, hd}, f, a, b, c, d, r);
    });
    for (int k = 0; k < 2 * N; ++k) {
        const double zL = g.draw(ETA, -2.0, 5.0), zR = g.below(4) == 0 ? zL : g.draw(ETA, -2.0, 5.0);
        const double hL = g.draw(HS, 0.0, 3.0), hR = g.draw(HS, 0.0, 3.0);
        const double etaL = g.below(16) == 0 ? g.draw(ETA, -2.0, 5.0) : canon(zL + hL);
        const double etaR = g.below(16) == 0 ? g.draw(ETA, -2.0, 5.0) : canon(zR + hR);
        const double qxL = g.draw(QX, -5.0, 5.0), qyL = g.draw(QX, -5.0, 5.0);
        const double qxR = g.draw(QX, -5.0, 5.0), qyR = g.draw(QX, -5.0, 5.0);
        double nx, ny;
        normal(g, nx, ny);
        const double hd = g.draw(HD, 0.0, 1e-2);
        FaceFlux f{-7, -7, -7, -7};
        double a = -7, b = -7, c = -7, d = -7;
        const bool r = faceFlux(opq(etaL), opq(hL), opq(qxL), opq(qyL), opq(etaR), opq(hR), opq(qxR), opq(qyR), opq(nx),
                                opq(ny), opq(hd), f, a, b, c, d);
        face_flux_line(o, "ff", {etaL, hL, qxL, qyL, etaR, hR, qxR, qyR, nx, ny, hd}, f, a, b, c, d, r);
    }

    // ---- faceFluxReconBed and faceFluxRecon
    for (int k = 0; k < 2 * N; ++k) {
        const double zL = g.draw(ETA, -2.0, 5.0), zR = g.below(4) == 0 ? zL : g.draw(ETA, -2.0, 5.0);
        const double hLc = g.draw(HS, 0.0, 3.0), hRc = g.draw(HS, 0.0, 3.0);
        const double zLf = g.below(3) == 0 ? zL : zL + g.uni(-0.2, 0.2);
        const double zRf = g.below(3) == 0 ? zR : zR + g.uni(-0.2, 0.2);
        const double etaLf = g.below(16) == 0 ? g.draw(ETA, -2.0, 5.0) : canon(zL + hLc + g.uni(-0.1, 0.1));
        const double etaRf = g.below(16) == 0 ? g.draw(ETA, -2.0, 5.0) : canon(zR + hRc + g.uni(-0.1, 0.1));
        const double uxL = g.draw(US, -5.0, 5.0), uyL = g.draw(US, -5.0, 5.0);
        const double uxR = g.draw(US, -5.0, 5.0), uyR = g.draw(US, -5.0, 5.0);
        double nx, ny;
        normal(g, nx, ny);
        const double hd = g.draw(HD, 0.0, 1e-2);
        FaceFlux f{-7, -7, -7, -7};
        double a = -7, b = -7, c = -7, d = -7;
        bool r = faceFluxReconBed(opq(etaLf), opq(uxL), opq(uyL), opq(zL), opq(hLc), opq(etaRf), opq(uxR), opq(uyR),
                                  opq(zR), opq(hRc), opq(nx), opq(ny), opq(hd), f, a, b, c, d, opq(zLf), opq(zRf));
        face_flux_line(o, "frb", {etaLf, uxL, uyL, zL, hLc, etaRf, uxR, uyR, zR, hRc, nx, ny, hd, zLf, zRf}, f, a, b, c,
                       d, r);
        FaceFlux f2{-7, -7, -7, -7};
        a = b = c = d = -7;
        r = faceFluxRecon(opq(etaLf), opq(uxL), opq(uyL), opq(zL), opq(hLc), opq(etaRf), opq(uxR), opq(uyR), opq(zR),
                          opq(hRc), opq(nx), opq(ny), opq(hd), f2, a, b, c, d);
        face_flux_line(o, "frc", {etaLf, uxL, uyL, zL, hLc, etaRf, uxR, uyR, zR, hRc, nx, ny, hd}, f2, a, b, c, d, r);
    }

    // ---- bjLimiter(w, wf, wmin, wmax)
    const V BW = {0.0, -0.0, 1.0, 1e-14, -1e-14, 2e-14, -2e-14, 0.5, -1.0, 2.0};
    grid({BW.size(), BW.size(), BW.size(), BW.size()}, [&](const std::vector<std::size_t>& i) {
        const double a = opq(BW[i[0]]), b = opq(BW[i[1]]), c = opq(BW[i[2]]), d = opq(BW[i[3]]);
        o.tag("bj"); args(o, {a, b, c, d}); o.hx(bjLimiter(a, b, c, d)); o.nl();
    });
    for (int k = 0; k < N; ++k) {
        const double w = g.draw(BW, -2.0, 2.0);
        const double wf = g.below(4) == 0 ? w + g.uni(-3e-14, 3e-14) : g.draw(BW, -3.0, 3.0);
        const double wmin = g.below(4) == 0 ? g.draw(BW, -3.0, 3.0) : w - g.uni(0.0, 1.0);
        const double wmax = g.below(4) == 0 ? g.draw(BW, -3.0, 3.0) : w + g.uni(0.0, 1.0);
        o.tag("bj"); args(o, {w, wf, wmin, wmax}); o.hx(bjLimiter(opq(w), opq(wf), opq(wmin), opq(wmax))); o.nl();
    }

    // ---- frictionUpdate(qx, qy, h, n, dt)
    const V NM = {0.0, 0.02, 0.05, 0.1, -0.03};
    const V DTS = {0.0, -0.0, 0.1, 1.0, 10.0, -1.0};
    for (int k = 0; k < 2 * N; ++k) {
        double qx = g.draw(QX, -5.0, 5.0), qy = g.draw(QX, -5.0, 5.0);
        const double h = g.draw(HS, 0.0, 3.0), n = g.draw(NM, 0.0, 0.1), dt = g.draw(DTS, 0.0, 30.0);
        o.tag("fr"); args(o, {qx, qy, h, n, dt});
        qx = opq(qx);
        qy = opq(qy);
        frictionUpdate(qx, qy, opq(h), opq(n), opq(dt));
        o.hx(qx); o.hx(qy); o.nl();
    }

    // ---- criticalDepth(q)
    for (int k = 0; k < N; ++k) {
        const double q = g.draw(QX, -10.0, 10.0);
        o.tag("cd"); args(o, {q}); o.hx(criticalDepth(opq(q))); o.nl();
    }

    // ---- depthFromInvariantAndDischarge(r, q, h_interior)
    const V RV = {0.0, 1.0, 5.0, 10.0, -1.0, 2.0 * std::sqrt(9.80665), 6.264, 20.0};
    const V QV = {0.0, -0.0, -1.0, 1.0, -0.5, 0.5, -5.0, 5.0, 1e-9};
    const V HV = {0.0, -0.0, 1.0, 0.5, 2.0, 1e-6, -1.0, 1e-12};
    grid({RV.size(), QV.size(), HV.size()}, [&](const std::vector<std::size_t>& i) {
        const double a = opq(RV[i[0]]), b = opq(QV[i[1]]), c = opq(HV[i[2]]);
        o.tag("dri"); args(o, {a, b, c}); o.hx(depthFromInvariantAndDischarge(a, b, c)); o.nl();
    });
    for (int k = 0; k < 2 * N; ++k) {
        // mostly physical: an interior state (h, u), r = u + 2 sqrt(g h), a nearby discharge
        const double h = g.draw(HV, 0.0, 3.0), u = g.draw(US, -3.0, 3.0);
        const double r = g.below(4) == 0 ? g.draw(RV, -5.0, 20.0) : canon(u + 2.0 * std::sqrt(9.80665 * h));
        const double q = g.below(4) == 0 ? g.draw(QV, -5.0, 5.0) : canon(u * h * g.uni(0.5, 1.5));
        o.tag("dri"); args(o, {r, q, h}); o.hx(depthFromInvariantAndDischarge(opq(r), opq(q), opq(h))); o.nl();
    }

    o.flush();
    std::fprintf(stderr, "twod_swe: %ld lines\n", o.lines);
    return 0;
}
