// Differential harness for the scalar kernels of src/engine/2d/solver/InertialKernels.hpp:
// momentumBlendFraction (glibc log, expm1), qMagnitude, etaDepthScalar, volumeFromEtaScalar,
// etaDepthQuadScalar, volumeFromEtaQuadScalar, faceFlowDepth, faceDepthFromEta,
// faceFlowDepthVfr, inertialFaceUpdate (glibc cbrt), inertialAdvection, froudeCap, cellCflDt and
// positivityScale. Each runs over a grid of its branch points and special values (zero, signed
// zero, tiny, the thresholds 1e-9, 1e-12, 1e-30 and their neighbours, both flow directions) where
// its arity allows, and over seeded random draws mixing specials, physical ranges, wide
// magnitudes and non-finite values. Prints every input and result as hex bits; the Rust twin is
// crates/engine/examples/twod_inertial.rs.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
//
// Build: as gate/diff/twod_vfr.cpp (gate/twodkern.sh), with twod_inertial in place of twod_vfr.
#include "2d/solver/InertialKernels.hpp"
#include "twodgen.hpp"

using namespace openswmm::twoD;
using namespace openswmm::twoD::inertial;
using namespace twodgen;
using V = std::vector<double>;

static void args(Out& o, std::initializer_list<double> a) {
    for (double v : a) o.hx(v);
    o.sep();
}

int main() {
    Out o;
    Gen g(0x2d0f1a7e5eed0003ull);
    const int N = 40000;

    // ---- momentumBlendFraction(theta, dt, reference_step)
    const V TH = {1.0, 0.5, 0.9, 0.99, 0.999999, 1e-3, 0.0, -0.0, -0.5, std::nextafter(1.0, 0.0), 2.0, DMIN};
    const V DT = {0.0, -0.0, 1e-12, 1e-9, 0.01, 1.0, 60.0, 3600.0, -1.0, DMIN};
    const V RS = {0.0, -0.0, 1.0, 0.5, 1e-6, 1e6, -1.0, DMIN};
    grid({TH.size(), DT.size(), RS.size()}, [&](const std::vector<std::size_t>& i) {
        const double a = opq(TH[i[0]]), b = opq(DT[i[1]]), c = opq(RS[i[2]]);
        o.tag("mbf"); args(o, {a, b, c}); o.hx(momentumBlendFraction(a, b, c)); o.nl();
    });
    for (int k = 0; k < N; ++k) {
        const double a = g.draw(TH, 0.0, 1.0), b = g.draw(DT, 0.0, 120.0), c = g.draw(RS, 0.0, 120.0);
        o.tag("mbf"); args(o, {a, b, c}); o.hx(momentumBlendFraction(opq(a), opq(b), opq(c))); o.nl();
    }

    // ---- qMagnitude(x, y)
    const V QS = {0.0, -0.0, 1.0, -1.0, 3.0, 4.0, 1e-160, 1e160, 1e-300, DMIN, 0.1, -2.5};
    grid({QS.size(), QS.size()}, [&](const std::vector<std::size_t>& i) {
        const double a = opq(QS[i[0]]), b = opq(QS[i[1]]);
        o.tag("qm"); args(o, {a, b}); o.hx(qMagnitude(a, b)); o.nl();
    });
    for (int k = 0; k < N; ++k) {
        const double a = g.draw(QS, -10.0, 10.0), b = g.draw(QS, -10.0, 10.0);
        o.tag("qm"); args(o, {a, b}); o.hx(qMagnitude(opq(a), opq(b))); o.nl();
    }

    // ---- etaDepthScalar / volumeFromEtaScalar (triangle closure, both FLAT and VFR)
    const V AR = {1.0, 0.0, -0.0, 1e-30, std::nextafter(1e-30, INF), std::nextafter(1e-30, 0.0), 25.0, 1e4};
    const V ZZ = {0.0, -0.0, 1.0, 1e-9, 100.0, -3.0, 0.5};
    const V EP = {0.0, -0.0, 1e-6, 0.01, 0.1, 0.5, 1.0};
    const V VV = {0.0, -0.0, -1.0, 1e-12, 1.0, 50.0, DMIN};
    for (int k = 0; k < N; ++k) {
        const double area = g.draw(AR, 0.0, 500.0), cz = g.draw(ZZ, -5.0, 50.0);
        const double base = g.uni(-5.0, 50.0);
        double z1, z2, z3;
        if (g.below(3) == 0) {
            z1 = g.draw(ZZ, -5.0, 5.0); z2 = g.draw(ZZ, -5.0, 5.0); z3 = g.draw(ZZ, -5.0, 5.0);
        } else {
            z1 = base + g.uni(0.0, 2.0); z2 = base + g.uni(0.0, 2.0); z3 = base + g.uni(0.0, 2.0);
        }
        const bool vfr = g.below(2) == 0;
        const double eps = g.draw(EP, 0.0, 0.3);
        const double vol = g.draw(VV, 0.0, 2.0 * area);
        const double eta = g.below(2) == 0 ? g.uni(z1 - 0.5, z3 + 2.0) : g.draw(ZZ, -5.0, 50.0);
        double e = -7.0, d = -7.0;
        etaDepthScalar(opq(area), opq(cz), opq(z1), opq(z2), opq(z3), vfr, opq(eps), opq(vol), e, d);
        o.tag("eds"); args(o, {area, cz, z1, z2, z3, vfr ? 1.0 : 0.0, eps, vol}); o.hx(e); o.hx(d); o.nl();
        o.tag("vfe"); args(o, {area, cz, z1, z2, z3, vfr ? 1.0 : 0.0, eps, eta});
        o.hx(volumeFromEtaScalar(opq(area), opq(cz), opq(z1), opq(z2), opq(z3), vfr, opq(eps), opq(eta)));
        o.nl();
    }

    // ---- etaDepthQuadScalar / volumeFromEtaQuadScalar
    for (int k = 0; k < N / 2; ++k) {
        double x[4] = {0, 1, 1, 0}, y[4] = {0, 0, 1, 1}, zq[4], zs[6], A1, A2;
        const double base = g.uni(-5.0, 50.0), sc = g.uni(0.1, 30.0);
        for (int v = 0; v < 4; ++v) { x[v] *= sc; y[v] *= sc; zq[v] = base + g.uni(0.0, 2.0); }
        if (g.below(5) == 0) zq[g.below(4)] = zq[g.below(4)];
        quadVfrPrecompute(x, y, zq, zs, A1, A2);
        const double area = g.below(8) == 0 ? g.pick(AR) : A1 + A2;
        const double cz = (zq[0] + zq[1] + zq[2] + zq[3]) / 4.0;
        const bool vfr = g.below(2) == 0;
        const double eps = g.draw(EP, 0.0, 0.3);
        const double vol = g.draw(VV, 0.0, 2.0 * area);
        const double eta = g.below(2) == 0 ? g.uni(base - 0.5, base + 3.0) : g.draw(ZZ, -5.0, 50.0);
        for (int v = 0; v < 6; ++v) zs[v] = opq(zs[v]);
        double e = -7.0, d = -7.0;
        etaDepthQuadScalar(opq(area), opq(cz), zs, opq(A1), opq(A2), vfr, opq(eps), opq(vol), e, d);
        o.tag("edq"); args(o, {area, cz, zs[0], zs[1], zs[2], zs[3], zs[4], zs[5], A1, A2, vfr ? 1.0 : 0.0, eps, vol});
        o.hx(e); o.hx(d); o.nl();
        o.tag("vfq"); args(o, {area, cz, zs[0], zs[1], zs[2], zs[3], zs[4], zs[5], A1, A2, vfr ? 1.0 : 0.0, eps, eta});
        o.hx(volumeFromEtaQuadScalar(opq(area), opq(cz), zs, opq(A1), opq(A2), vfr, opq(eps), opq(eta)));
        o.nl();
    }

    // ---- faceFlowDepth, faceDepthFromEta, faceFlowDepthVfr
    const V ET = {0.0, -0.0, 1.0, 1.0 + 1e-9, 1.0 + 5e-10, 2.0, 0.5, -1.0, 1e-12, 3.0};
    grid({ET.size(), ET.size(), ET.size()}, [&](const std::vector<std::size_t>& i) {
        const double a = opq(ET[i[0]]), b = opq(ET[i[1]]), c = opq(ET[i[2]]);
        o.tag("ffd"); args(o, {a, b, c}); o.hx(faceFlowDepth(a, b, c)); o.nl();
        o.tag("fde"); args(o, {a, b, c}); o.hx(faceDepthFromEta(a, b, c)); o.nl();
    });
    for (int k = 0; k < N; ++k) {
        const double zlo = g.draw(ET, -2.0, 10.0);
        const int m = static_cast<int>(g.below(4));
        const double zhi = m == 0 ? zlo + g.uni(0.0, 2e-9) : m == 1 ? g.draw(ET, -2.0, 10.0) : zlo + g.uni(0.0, 3.0);
        const double sel[] = {zlo, zhi, std::nextafter(zlo, INF), std::nextafter(zhi, INF), std::nextafter(zhi, -INF)};
        const double eta = g.below(3) == 0 ? sel[g.below(5)] : g.uni(zlo - 1.0, zhi + 1.0);
        const double eta2 = g.draw(ET, zlo - 1.0, zhi + 1.0);
        o.tag("fde"); args(o, {eta, zlo, zhi}); o.hx(faceDepthFromEta(opq(eta), opq(zlo), opq(zhi))); o.nl();
        o.tag("ffv"); args(o, {eta, eta2, zlo, zhi}); o.hx(faceFlowDepthVfr(opq(eta), opq(eta2), opq(zlo), opq(zhi))); o.nl();
        o.tag("ffd"); args(o, {eta, eta2, zlo}); o.hx(faceFlowDepth(opq(eta), opq(eta2), opq(zlo))); o.nl();
    }

    // ---- inertialFaceUpdate(q, qhat, hf, dt, slope, n2, q_mag, adv)
    const V HF = {1e-3, 1.0, 0.0, -0.0, 1e-12, 1e-6, 5.0, DMIN, 1e-300, -1.0};
    const V SL = {0.0, -0.0, 1e-3, -1e-3, 0.05, -0.05, 1e-12, 1.0};
    const V QQ = {0.0, -0.0, 0.5, -0.5, 2.0, 1e-9, -10.0};
    const V NN = {0.0, 0.0004, 0.0025, 0.01, 1.0};
    for (int k = 0; k < 2 * N; ++k) {
        const double q = g.draw(QQ, -5.0, 5.0), qhat = g.draw(QQ, -5.0, 5.0), hf = g.draw(HF, 0.0, 3.0);
        const double dt = g.draw(DT, 0.0, 10.0), slope = g.draw(SL, -0.1, 0.1), n2 = g.draw(NN, 0.0, 0.01);
        const double qmag = g.draw(QQ, 0.0, 5.0), adv = g.below(2) == 0 ? 0.0 : g.draw(QQ, -1.0, 1.0);
        o.tag("ifu"); args(o, {q, qhat, hf, dt, slope, n2, qmag, adv});
        o.hx(inertialFaceUpdate(opq(q), opq(qhat), opq(hf), opq(dt), opq(slope), opq(n2), opq(qmag), opq(adv)));
        o.nl();
    }
    // the defaulted adv = 0.0
    for (int k = 0; k < N / 4; ++k) {
        const double q = g.draw(QQ, -5.0, 5.0), qhat = g.draw(QQ, -5.0, 5.0), hf = g.draw(HF, 0.0, 3.0);
        const double dt = g.draw(DT, 0.0, 10.0), slope = g.draw(SL, -0.1, 0.1), n2 = g.draw(NN, 0.0, 0.01);
        const double qmag = g.draw(QQ, 0.0, 5.0);
        o.tag("ifu0"); args(o, {q, qhat, hf, dt, slope, n2, qmag});
        o.hx(inertialFaceUpdate(opq(q), opq(qhat), opq(hf), opq(dt), opq(slope), opq(n2), opq(qmag)));
        o.nl();
    }

    // ---- inertialAdvection(q_f, unL, hL, unR, hR, inv_dx_normal)
    const V UN = {0.0, -0.0, 1.0, -1.0, 0.3, -0.3, 1e-12, -1e-12};
    for (int k = 0; k < N; ++k) {
        const double qf = g.draw(QQ, -5.0, 5.0), unL = g.draw(UN, -3.0, 3.0), hL = g.draw(HF, 0.0, 3.0);
        const double unR = g.draw(UN, -3.0, 3.0), hR = g.draw(HF, 0.0, 3.0), idx = g.draw(AR, 0.01, 2.0);
        o.tag("adv"); args(o, {qf, unL, hL, unR, hR, idx});
        o.hx(inertialAdvection(opq(qf), opq(unL), opq(hL), opq(unR), opq(hR), opq(idx)));
        o.nl();
    }

    // ---- froudeCap(q, hf, fr_max) (fr_max < 0 included: there std::clamp's hi < lo)
    const V FR = {1.0, 0.0, -0.0, 0.9, 2.0, -1.0, 1e-9};
    grid({QQ.size(), HF.size(), FR.size()}, [&](const std::vector<std::size_t>& i) {
        const double a = opq(QQ[i[0]]), b = opq(HF[i[1]]), c = opq(FR[i[2]]);
        o.tag("fc"); args(o, {a, b, c}); o.hx(froudeCap(a, b, c)); o.nl();
    });
    for (int k = 0; k < N; ++k) {
        const double q = g.draw(QQ, -20.0, 20.0), hf = g.draw(HF, 0.0, 3.0), fr = g.draw(FR, 0.0, 3.0);
        o.tag("fc"); args(o, {q, hf, fr}); o.hx(froudeCap(opq(q), opq(hf), opq(fr))); o.nl();
    }

    // ---- cellCflDt(alpha, lchar, h, speed)
    const V AL = {0.5, 0.9, 1.0, 0.0, -0.5};
    const V LC = {1.0, 10.0, 0.0, 1e-3, 1e4};
    const V HH = {0.0, -0.0, 1e-30, 1.0e-25, 1.0198e-25, 1.02e-25, 1.0, 1e-6, -1.0};
    const V SP = {0.0, -0.0, 1e-12, std::nextafter(1e-12, INF), 1.0, 5.0, -1.0};
    grid({AL.size(), LC.size(), HH.size(), SP.size()}, [&](const std::vector<std::size_t>& i) {
        const double a = opq(AL[i[0]]), b = opq(LC[i[1]]), c = opq(HH[i[2]]), d = opq(SP[i[3]]);
        o.tag("cfl"); args(o, {a, b, c, d}); o.hx(cellCflDt(a, b, c, d)); o.nl();
    });
    for (int k = 0; k < N; ++k) {
        const double a = g.draw(AL, 0.1, 1.0), b = g.draw(LC, 0.1, 100.0), c = g.draw(HH, 0.0, 3.0), d = g.draw(SP, 0.0, 5.0);
        o.tag("cfl"); args(o, {a, b, c, d}); o.hx(cellCflDt(opq(a), opq(b), opq(c), opq(d))); o.nl();
    }

    // ---- positivityScale(V, outflow_m3s, dt, beta)
    const V PV = {0.0, -0.0, -1.0, 1.0, 10.0, 1e-12, DMIN};
    const V PO = {0.0, -0.0, -1.0, 1.0, 0.1, 100.0, 1e-12};
    const V PD = {0.0, -0.0, 1.0, 0.1, 10.0, -1.0};
    const V PB = {1.0, 0.0, 0.5, 0.95};
    grid({PV.size(), PO.size(), PD.size(), PB.size()}, [&](const std::vector<std::size_t>& i) {
        const double a = opq(PV[i[0]]), b = opq(PO[i[1]]), c = opq(PD[i[2]]), d = opq(PB[i[3]]);
        o.tag("pos"); args(o, {a, b, c, d}); o.hx(positivityScale(a, b, c, d)); o.nl();
    });
    for (int k = 0; k < N; ++k) {
        const double a = g.draw(PV, -1.0, 20.0), b = g.draw(PO, 0.0, 20.0), c = g.draw(PD, 0.0, 5.0), d = g.draw(PB, 0.0, 1.0);
        o.tag("pos"); args(o, {a, b, c, d}); o.hx(positivityScale(opq(a), opq(b), opq(c), opq(d))); o.nl();
    }

    o.flush();
    std::fprintf(stderr, "twod_inertial: %ld lines\n", o.lines);
    return 0;
}
