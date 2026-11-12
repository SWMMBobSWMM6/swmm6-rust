// The run-time 1D-2D exchange differential, upstream side (phase 2a wave 4b). Runs a coupled
// deck on the 2D reference (build2d) through its C API, and interposes on the exchange calls
// the library makes through its PLT (the harness is linked with -Wl,-E, so the library's calls
// bind to the definitions below, which find the library's own with dlsym(RTLD_NEXT)):
//   Q   computeNodeCouplingQ        (the marcher's live junction exchange)
//   G   computeNodeCouplingDQdh1d   (inside computeCouplingConductances)
//   CC  SurfaceRouter2D::computeCouplingConductances   (stepRouting's DW callback)
//   OB  updateOutfallBoundaries     \
//   TW  OutfallExchange::recordTailwater   } updateOutfallsPreRouting
//   PR  OutfallExchange::prepare     |
//   TP  SurfaceRouter2D::updateTrackingPonds /
//   AP  OutfallExchange::apply      (advancePostRouting)
//   EQ  SurfaceRouter2D::publishEqualizeAreas (coAdvanceStep, before the advance)
//   DW  dynwave::DWSolver::initNodeStates (each Picard iteration: the coupled junctions' pond
//       area, COUPLING_POND_MODE / COUPLING_RIM_BAND, and the node flow reset)
//   CA  SurfaceRouter2D::coAdvanceStep: its junction ledger (exch -> coupling_volume -> queue;
//       coupling_cum_applied_ is not compared: accumulateMassBalance, wave 4a's, also adds to it)
// Each interposer captures the inputs the call reads (the router's, state's, context's
// arrays, whole), calls upstream's compiled function, and prints one record: the inputs, then
// the outputs, every double as its 64-bit pattern in hex. Before each coAdvanceStep, on a copy
// of the live state, it also calls the three functions upstream never calls (SY
// scatterCouplingToYdot, AD accumulateOutfallDischargeStep, IN injectAccumulatedExchange), and
// for every Q and G record a few variants of the same call (provisional volumes, a head
// offset, the equalize flag and the law flipped), so the paths the run does not take are
// checked on its real states. At the end, TQ: OutfallExchange::tailwater for every node.
// The Rust twin, crates/engine/examples/twod_exch.rs, opens, initializes and starts the deck in
// the translation, reads these records on stdin, loads each record's inputs into its own
// engine, calls the translation and compares the outputs bit for bit.
//
// Record format: "R <kind> <seq> <step>", then lines "i <name> <n> <hex>..." (inputs) and
// "o <name> <n> <hex>..." (outputs; "t <what>" when the call throws), then ".".
// Usage: twod_exch deck.inp [max_steps]   (env TWOD_EXCH_QEVERY=n: one Q record in n)
// Build (gate/twodexch.sh): build2d's compile line, -DOPENSWMM_HAS_2D=1, -Wl,-E -ldl.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cassert>
#include <cctype>
#include <cerrno>
#include <climits>
#include <clocale>
#include <cmath>
#include <compare>
#include <concepts>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cwchar>
#include <cwctype>
#include <deque>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <ios>
#include <iosfwd>
#include <istream>
#include <limits>
#include <locale>
#include <map>
#include <memory>
#include <mutex>
#include <new>
#include <numbers>
#include <numeric>
#include <optional>
#include <ostream>
#include <queue>
#include <ratio>
#include <semaphore>
#include <set>
#include <sstream>
#include <stdexcept>
#include <stop_token>
#include <streambuf>
#include <string>
#include <string_view>
#include <system_error>
#include <thread>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>
#include <omp.h>
#include <dlfcn.h>
#define private public
#define protected public
#include "openswmm/engine/openswmm_engine.h"
#include "core/SWMMEngine.hpp"
#include "2d/SurfaceRouter2D.hpp"
#include "2d/coupling/NodeCoupling.hpp"
#include "2d/coupling/OutfallExchange.hpp"

using namespace openswmm;
using namespace openswmm::twoD;

static SWMMEngine* g_eng = nullptr;
static long g_seq = 0;
static int g_step = 0;
static long g_qcount = 0, g_qevery = 1;
static bool g_on = false; // records only between start and end

template <class F> static F realfn(const char* sym) {
    void* p = dlsym(RTLD_NEXT, sym);
    if (!p) {
        std::fprintf(stderr, "twod_exch: %s not found\n", sym);
        std::_Exit(3);
    }
    return reinterpret_cast<F>(p);
}

// One record, assembled in memory (inputs are captured before a call that changes them).
struct Rec {
    std::string s;
    void head(const char* kind) {
        char b[96];
        std::snprintf(b, sizeof b, "R %s %ld %d\n", kind, ++g_seq, g_step);
        s += b;
    }
    void arr(char tag, const char* name, const double* p, std::size_t n) {
        char b[48];
        std::snprintf(b, sizeof b, "%c %s %zu", tag, name, n);
        s += b;
        for (std::size_t i = 0; i < n; ++i) {
            uint64_t u;
            std::memcpy(&u, &p[i], 8);
            std::snprintf(b, sizeof b, " %016llx", (unsigned long long)u);
            s += b;
        }
        s += '\n';
    }
    void arr(char tag, const char* name, const std::vector<double>& v) { arr(tag, name, v.data(), v.size()); }
    void iarr(char tag, const char* name, const std::vector<int>& v) {
        std::vector<double> d(v.begin(), v.end());
        arr(tag, name, d);
    }
    void sc(char tag, const char* name, double v) { arr(tag, name, &v, 1); }
    void thrown(const char* what) {
        s += "t ";
        s += what;
        s += '\n';
    }
    void flush() {
        s += ".\n";
        std::fwrite(s.data(), 1, s.size(), stdout);
        s.clear();
    }
};

static SurfaceRouter2D& router() { return g_eng->surface_router_; }

// The node arrays the exchange reads.
static void nodes_in(Rec& r, const NodeData& n) {
    r.arr('i', "nhead", n.head);
    r.arr('i', "ndepth", n.depth);
    r.arr('i', "ninv", n.invert_elev);
    r.arr('i', "nfull", n.full_depth);
}

// The state arrays and attachments the junction exchange reads.
static void state_in(Rec& r, const SurfaceStateData& st) {
    r.arr('i', "depth", st.depth);
    r.arr('i', "head", st.head);
    r.arr('i', "vol", st.volume);
    r.sc('i', "win", st.coupling_window_s);
    r.sc('i', "hc", st.coupling_head_cells ? 1.0 : 0.0);
    if (st.coupling_head_cells)
        r.iarr('i', "hcells", *st.coupling_head_cells);
    r.sc('i', "pa", st.node_pond_area_m2 ? 1.0 : 0.0);
    if (st.node_pond_area_m2)
        r.arr('i', "pond", *st.node_pond_area_m2);
}

// The k-th junction point as buildCouplingPoints made it (the k-th non-outfall point of
// coupling_points_): initialize moves a vertex-coupled junction point to its lowest stencil cell
// (vertex_idx -1), so the exchange's vertex branches run only in the variants below.
static const CouplingPoint* original_point(int k) {
    int j = -1;
    for (const auto& p : router().coupling_points_)
        if (!p.is_outfall && ++j == k)
            return &p;
    return nullptr;
}

static int point_index(const CouplingPoint& cp) {
    const auto& v = router().node_coupling_points_;
    if (v.empty() || &cp < v.data() || &cp >= v.data() + v.size())
        return -1;
    return static_cast<int>(&cp - v.data());
}

using QFn = double (*)(const CouplingPoint&, const MeshData&, const SurfaceStateData&, const NodeData&,
                       const SolverOptions2D&, const double*, double);
using GFn = double (*)(const CouplingPoint&, const MeshData&, const SurfaceStateData&, const NodeData&,
                       const SolverOptions2D&);
static QFn realQ() {
    static QFn f = realfn<QFn>("_ZN8openswmm4twoD20computeNodeCouplingQERKNS0_13CouplingPointERKNS0_8MeshDataERKNS0_16SurfaceStateDataERKNS_8NodeDataERKNS0_15SolverOptions2DEPKdd");
    return f;
}
static GFn realG() {
    static GFn f = realfn<GFn>("_ZN8openswmm4twoD25computeNodeCouplingDQdh1dERKNS0_13CouplingPointERKNS0_8MeshDataERKNS0_16SurfaceStateDataERKNS_8NodeDataERKNS0_15SolverOptions2DE");
    return f;
}

namespace openswmm::twoD {

double computeNodeCouplingQ(const CouplingPoint& cp, const MeshData& mesh, const SurfaceStateData& state,
                            const NodeData& nodes, const SolverOptions2D& opts, const double* prov,
                            double h_off) noexcept {
    const double q = realQ()(cp, mesh, state, nodes, opts, prov, h_off);
    if (!g_on || (g_qcount++ % g_qevery) != 0)
        return q;
    Rec r;
    r.head("Q");
    const int k = point_index(cp);
    r.sc('i', "k", k);
    r.sc('i', "hoff", h_off);
    r.sc('i', "prov", prov ? 1.0 : 0.0);
    if (prov)
        r.arr('i', "provv", prov, state.volume.size());
    state_in(r, state);
    nodes_in(r, nodes);
    r.sc('o', "q", q);
    // Variants on the same state: provisional volumes (the state's, a quarter of them, none),
    // a head offset either way, the equalize flag and the law flipped.
    std::vector<double> quarter(state.volume.size()), zero(state.volume.size(), 0.0);
    for (std::size_t i = 0; i < quarter.size(); ++i)
        quarter[i] = state.volume[i] * 0.25;
    CouplingPoint eq = cp, law = cp;
    eq.equalize = !cp.equalize;
    law.law = cp.law == 1 ? 0 : 1;
    double v[7] = {realQ()(cp, mesh, state, nodes, opts, state.volume.data(), h_off),
                   realQ()(cp, mesh, state, nodes, opts, quarter.data(), h_off),
                   realQ()(cp, mesh, state, nodes, opts, zero.data(), h_off),
                   realQ()(cp, mesh, state, nodes, opts, prov, h_off + 0.05),
                   realQ()(cp, mesh, state, nodes, opts, prov, h_off - 0.05),
                   realQ()(eq, mesh, state, nodes, opts, prov, h_off),
                   realQ()(law, mesh, state, nodes, opts, prov, h_off)};
    r.arr('o', "qv", v, 7);
    // The vertex branches: the point at its original vertex (live and provisional heads,
    // either law); -12345 when the junction was coupled at a cell.
    double w[4] = {-12345.0, -12345.0, -12345.0, -12345.0};
    if (const auto* o = original_point(k); o && o->vertex_idx >= 0) {
        CouplingPoint vx = cp, vl;
        vx.vertex_idx = o->vertex_idx;
        vl = vx;
        vl.law = law.law;
        w[0] = realQ()(vx, mesh, state, nodes, opts, prov, h_off);
        w[1] = realQ()(vx, mesh, state, nodes, opts, state.volume.data(), h_off);
        w[2] = realQ()(vl, mesh, state, nodes, opts, prov, h_off);
        w[3] = realQ()(vx, mesh, state, nodes, opts, quarter.data(), h_off);
    }
    r.arr('o', "qw", w, 4);
    r.flush();
    return q;
}

double computeNodeCouplingDQdh1d(const CouplingPoint& cp, const MeshData& mesh,
                                 const SurfaceStateData& state, const NodeData& nodes,
                                 const SolverOptions2D& opts) noexcept {
    const double g = realG()(cp, mesh, state, nodes, opts);
    if (!g_on)
        return g;
    Rec r;
    r.head("G");
    r.sc('i', "k", point_index(cp));
    state_in(r, state);
    nodes_in(r, nodes);
    r.sc('o', "g", g);
    CouplingPoint law = cp;
    law.law = cp.law == 1 ? 0 : 1;
    r.sc('o', "gv", realG()(law, mesh, state, nodes, opts));
    double w[2] = {-12345.0, -12345.0};
    if (const auto* o = original_point(point_index(cp)); o && o->vertex_idx >= 0) {
        CouplingPoint vx = cp;
        vx.vertex_idx = o->vertex_idx;
        w[0] = realG()(vx, mesh, state, nodes, opts);
        vx.law = law.law;
        w[1] = realG()(vx, mesh, state, nodes, opts);
    }
    r.arr('o', "gw", w, 2);
    r.flush();
    return g;
}

void SurfaceRouter2D::computeCouplingConductances(const SimulationContext& ctx,
                                                  std::vector<std::pair<int, double>>& out) const {
    using F = void (*)(const SurfaceRouter2D*, const SimulationContext&, std::vector<std::pair<int, double>>&);
    static F real = realfn<F>("_ZNK8openswmm4twoD15SurfaceRouter2D27computeCouplingConductancesERKNS_17SimulationContextERSt6vectorISt4pairIidESaIS7_EE");
    const std::size_t n0 = out.size();
    real(this, ctx, out);
    if (!g_on)
        return;
    // A const call: the inputs are unchanged; the nested G records come first.
    Rec r;
    r.head("CC");
    r.sc('i', "n0", static_cast<double>(n0));
    state_in(r, state_);
    nodes_in(r, ctx.nodes);
    std::vector<double> pairs;
    for (std::size_t i = n0; i < out.size(); ++i) {
        pairs.push_back(out[i].first);
        pairs.push_back(out[i].second);
    }
    r.arr('o', "pairs", pairs);
    r.flush();
}

void updateOutfallBoundaries(const std::vector<CouplingPoint>& cps, const MeshData& mesh,
                             const SurfaceStateData& state, SimulationContext& ctx,
                             const SolverOptions2D& opts, double dt) {
    using F = void (*)(const std::vector<CouplingPoint>&, const MeshData&, const SurfaceStateData&,
                       SimulationContext&, const SolverOptions2D&, double);
    static F real = realfn<F>("_ZN8openswmm4twoD23updateOutfallBoundariesERKSt6vectorINS0_13CouplingPointESaIS2_EERKNS0_8MeshDataERKNS0_16SurfaceStateDataERNS_17SimulationContextERKNS0_15SolverOptions2DEd");
    Rec r;
    if (g_on) {
        r.head("OB");
        r.sc('i', "all", &cps == &router().coupling_points_ ? 1.0 : 0.0);
        r.sc('i', "dt", dt);
        state_in(r, state);
        r.arr('i', "h2d", ctx.node_subtypes.outfalls.head_2d);
        r.arr('i', "r2d", ctx.node_subtypes.outfalls.ramp_2d);
    }
    real(cps, mesh, state, ctx, opts, dt);
    if (g_on) {
        r.arr('o', "h2d", ctx.node_subtypes.outfalls.head_2d);
        r.arr('o', "r2d", ctx.node_subtypes.outfalls.ramp_2d);
        r.flush();
    }
}

static void ofx_points_in(Rec& r, char tag, const OutfallExchange& x) {
    std::vector<double> res, tot;
    for (const auto& p : x.points_) {
        res.insert(res.end(), p.reserved.begin(), p.reserved.end());
        tot.push_back(p.total);
    }
    r.arr(tag, "reserved", res);
    r.arr(tag, "total", tot);
}

void OutfallExchange::recordTailwater(const SimulationContext& ctx, double dt) {
    using F = void (*)(OutfallExchange*, const SimulationContext&, double);
    static F real = realfn<F>("_ZN8openswmm4twoD15OutfallExchange15recordTailwaterERKNS_17SimulationContextEd");
    Rec r;
    if (g_on) {
        r.head("TW");
        r.sc('i', "dt", dt);
        r.arr('i', "h2d", ctx.node_subtypes.outfalls.head_2d);
        r.arr('i', "r2d", ctx.node_subtypes.outfalls.ramp_2d);
        r.arr('i', "eng", tw_engaged_s_);
        r.arr('i', "sub", tw_submerged_s_);
        r.arr('i', "max", tw_max_);
    }
    real(this, ctx, dt);
    if (g_on) {
        r.arr('o', "eng", tw_engaged_s_);
        r.arr('o', "sub", tw_submerged_s_);
        r.arr('o', "max", tw_max_);
        r.flush();
    }
}

void OutfallExchange::prepare(const MeshData& mesh, const SurfaceStateData& state,
                              const SolverOptions2D& opts, SimulationContext& ctx, double dt) {
    using F = void (*)(OutfallExchange*, const MeshData&, const SurfaceStateData&, const SolverOptions2D&,
                       SimulationContext&, double);
    static F real = realfn<F>("_ZN8openswmm4twoD15OutfallExchange7prepareERKNS0_8MeshDataERKNS0_16SurfaceStateDataERKNS0_15SolverOptions2DERNS_17SimulationContextEd");
    Rec r;
    if (g_on) {
        r.head("PR");
        r.sc('i', "dt", dt);
        r.arr('i', "vol", state.volume);
        r.sc('i', "pdt", prepared_dt_);
        r.sc('i', "prep", prepared_ ? 1.0 : 0.0);
        r.arr('i', "nres", node_reserved_);
        r.arr('i', "limit", ctx.surface_outfall_link_limit);
        ofx_points_in(r, 'i', *this);
    }
    real(this, mesh, state, opts, ctx, dt);
    if (g_on) {
        r.sc('o', "pdt", prepared_dt_);
        r.sc('o', "prep", prepared_ ? 1.0 : 0.0);
        r.arr('o', "nres", node_reserved_);
        r.arr('o', "limit", ctx.surface_outfall_link_limit);
        ofx_points_in(r, 'o', *this);
        r.flush();
    }
}

bool OutfallExchange::apply(const MeshData& mesh, SurfaceStateData& state, const SolverOptions2D& opts,
                            SimulationContext& ctx, double dt, std::vector<double>& accumulated,
                            const std::vector<double>* node_conc, std::vector<double>* cell_accumulated) {
    using F = bool (*)(OutfallExchange*, const MeshData&, SurfaceStateData&, const SolverOptions2D&,
                       SimulationContext&, double, std::vector<double>&, const std::vector<double>*,
                       std::vector<double>*);
    static F real = realfn<F>("_ZN8openswmm4twoD15OutfallExchange5applyERKNS0_8MeshDataERNS0_16SurfaceStateDataERKNS0_15SolverOptions2DERNS_17SimulationContextEdRSt6vectorIdSaIdEEPKSE_PSE_");
    Rec r;
    auto outs = [&](char tag) {
        r.arr(tag, "vol", state.volume);
        r.arr(tag, "depth", state.depth);
        r.arr(tag, "head", state.head);
        r.arr(tag, "capp", state.coupling_applied);
        r.arr(tag, "pend", pending_);
        r.arr(tag, "acc", accumulated);
        r.sc(tag, "mbin", ctx.mass_balance_2d.outfall_in);
        r.sc(tag, "mbout", ctx.mass_balance_2d.outfall_out);
        r.sc(tag, "mbfs", ctx.mass_balance_2d.final_storage);
        r.sc(tag, "prep", prepared_ ? 1.0 : 0.0);
        r.arr(tag, "runoff", state.transport.cell_runoff_vol);
        if (cell_accumulated)
            r.arr(tag, "cacc", *cell_accumulated);
    };
    if (g_on) {
        r.head("AP");
        r.sc('i', "dt", dt);
        r.sc('i', "ns", state.transport.n_species);
        r.sc('i', "nc", node_conc ? 1.0 : 0.0);
        if (node_conc)
            r.arr('i', "conc", *node_conc);
        r.sc('i', "ca", cell_accumulated ? 1.0 : 0.0);
        r.arr('i', "inflow", ctx.nodes.inflow);
        r.arr('i', "outflow", ctx.nodes.outflow);
        r.arr('i', "nres", node_reserved_);
        ofx_points_in(r, 'i', *this);
        outs('i');
    }
    bool changed = false;
    try {
        changed = real(this, mesh, state, opts, ctx, dt, accumulated, node_conc, cell_accumulated);
    } catch (const std::exception& e) {
        if (g_on) {
            r.thrown(e.what());
            outs('o');
            r.sc('o', "recv", static_cast<double>(recv_.size()));
            r.flush();
        }
        throw;
    }
    if (g_on) {
        r.sc('o', "changed", changed ? 1.0 : 0.0);
        outs('o');
        r.sc('o', "recv", static_cast<double>(recv_.size()));
        r.flush();
    }
    return changed;
}

void SurfaceRouter2D::updateTrackingPonds(SimulationContext& ctx, double dt) {
    using F = void (*)(SurfaceRouter2D*, SimulationContext&, double);
    static F real = realfn<F>("_ZN8openswmm4twoD15SurfaceRouter2D19updateTrackingPondsERNS_17SimulationContextEd");
    Rec r;
    if (g_on) {
        r.head("TP");
        r.sc('i', "dt", dt);
        r.arr('i', "depth", state_.depth);
        r.arr('i', "parea", ctx.coupled_pond_area);
        r.arr('i', "pband", ctx.coupled_pond_band);
        r.arr('i', "ndepth", ctx.nodes.depth);
        r.arr('i', "nhead", ctx.nodes.head);
    }
    real(this, ctx, dt);
    if (g_on) {
        r.arr('o', "parea", ctx.coupled_pond_area);
        r.arr('o', "ndepth", ctx.nodes.depth);
        r.arr('o', "nhead", ctx.nodes.head);
        r.flush();
    }
}

void SurfaceRouter2D::publishEqualizeAreas(const SimulationContext& ctx, double window) {
    using F = void (*)(SurfaceRouter2D*, const SimulationContext&, double);
    static F real = realfn<F>("_ZN8openswmm4twoD15SurfaceRouter2D20publishEqualizeAreasERKNS_17SimulationContextEd");
    Rec r;
    if (g_on) {
        r.head("EQ");
        r.sc('i', "win", window);
        r.arr('i', "ndepth", ctx.nodes.depth);
        r.arr('i', "parea", ctx.coupled_pond_area);
        r.arr('i', "pband", ctx.coupled_pond_band);
        r.arr('i', "pond", node_pond_area_m2_);
    }
    real(this, ctx, window);
    if (g_on) {
        r.sc('o', "pa", state_.node_pond_area_m2 ? 1.0 : 0.0);
        r.sc('o', "pa_target", state_.node_pond_area_m2 == &node_pond_area_m2_ ? 1.0 : 0.0);
        r.sc('o', "win", state_.coupling_window_s);
        r.arr('o', "pond", node_pond_area_m2_);
        r.flush();
    }
}

// The three functions upstream never calls, on a copy of the live state before a batch.
static void window_helpers(SurfaceRouter2D& rt, SimulationContext& ctx, double dt) {
    const auto& mesh = rt.mesh_;
    const auto& opts = rt.options_;
    const auto& st = rt.state_;
    const auto& pts = rt.node_coupling_points_;
    const auto& all = rt.coupling_points_;
    const auto nc = st.volume.size();
    // SY: each junction point's live Q, and +/-1e-3 m3/s, scattered onto a zero ydot.
    for (std::size_t k = 0; k < pts.size(); ++k) {
        Rec r;
        r.head("SY");
        r.sc('i', "k", static_cast<double>(k));
        state_in(r, st);
        nodes_in(r, ctx.nodes);
        const double qs[3] = {realQ()(pts[k], mesh, st, ctx.nodes, opts, nullptr, 0.0), 1.0e-3, -1.0e-3};
        r.arr('o', "qs", qs, 3);
        for (double q : qs) {
            std::vector<double> ydot(nc, 0.0);
            scatterCouplingToYdot(mesh, st, pts[k], q, ydot.data());
            r.arr('o', "ydot", ydot);
        }
        // the vertex branch: the point at its original vertex
        if (const auto* o = original_point(static_cast<int>(k)); o && o->vertex_idx >= 0) {
            CouplingPoint vx = pts[k];
            vx.vertex_idx = o->vertex_idx;
            for (double q : qs) {
                std::vector<double> ydot(nc, 0.0);
                scatterCouplingToYdot(mesh, st, vx, q, ydot.data());
                r.arr('o', "ydotv", ydot);
            }
        }
        r.flush();
    }
    // AD: one routing step's outfall discharge accumulation, on the full budget and on a
    // quarter of it (so withdrawals clamp), with and without the sample row.
    if (!all.empty()) {
        Rec r;
        r.head("AD");
        r.sc('i', "dt", dt);
        r.arr('i', "vol", st.volume);
        r.arr('i', "head", st.head);
        r.arr('i', "inflow", ctx.nodes.inflow);
        r.arr('i', "outflow", ctx.nodes.outflow);
        for (int pass = 0; pass < 2; ++pass) {
            std::vector<double> accum(all.size(), 0.0), budget(st.volume), sample(all.size(), -7.0);
            if (pass == 1)
                for (auto& b : budget)
                    b *= 0.25;
            const int clamped = accumulateOutfallDischargeStep(all, mesh, st, ctx, opts, dt, accum, budget,
                                                               pass == 0 ? sample.data() : nullptr);
            r.sc('o', "clamped", clamped);
            r.arr('o', "accum", accum);
            r.arr('o', "budget", budget);
            r.arr('o', "sample", sample);
            // IN: the accumulated volumes injected over the window, both signs.
            for (double sign : {1.0, -1.0}) {
                SurfaceStateData cp = st;
                std::fill(cp.coupling_flux.begin(), cp.coupling_flux.end(), 0.0);
                std::vector<double> acc2(accum);
                for (std::size_t k = 0; k < acc2.size(); ++k)
                    if (acc2[k] == 0.0)
                        acc2[k] = 1.0e-3 * static_cast<double>(k + 1); // a junction slot too
                injectAccumulatedExchange(all, mesh, cp, acc2, dt, sign, nullptr);
                r.arr('o', "flux", cp.coupling_flux);
            }
        }
        r.flush();
    }
}

void SurfaceRouter2D::coAdvanceStep(SimulationContext& ctx, double dt, double t) {
    using F = void (*)(SurfaceRouter2D*, SimulationContext&, double, double);
    static F real = realfn<F>("_ZN8openswmm4twoD15SurfaceRouter2D13coAdvanceStepERNS_17SimulationContextEdd");
    Rec r;
    if (g_on && dt > 0.0) {
        window_helpers(*this, ctx, dt);
        r.head("CA");
        r.sc('i', "dt", dt);
        r.sc('i', "common", common_stage_ ? 1.0 : 0.0);
        r.sc('i', "ns", state_.transport.n_species);
        r.arr('i', "queue", ctx.nodes.coupling_queue);
        r.arr('i', "cvol", ctx.nodes.coupling_volume);
        r.sc('i', "cdr", ctx.coupling_delivery_remaining);
    }
    real(this, ctx, dt, t);
    if (g_on && dt > 0.0) {
        r.arr('i', "exch", solver_ ? solver_->last_coupling_exchange() : std::vector<double>{});
        r.arr('o', "queue", ctx.nodes.coupling_queue);
        r.arr('o', "cvol", ctx.nodes.coupling_volume);
        r.sc('o', "cdr", ctx.coupling_delivery_remaining);
        r.flush();
    }
}

} // namespace openswmm::twoD

void openswmm::dynwave::DWSolver::initNodeStates(SimulationContext& ctx) {
    using F = void (*)(DWSolver*, SimulationContext&);
    static F real = realfn<F>("_ZN8openswmm7dynwave8DWSolver14initNodeStatesERNS_17SimulationContextE");
    Rec r;
    if (g_on) {
        r.head("DW");
        r.sc('i', "pond", ctx.options.allow_ponding ? 1.0 : 0.0);
        r.arr('i', "cnode", std::vector<double>(ctx.coupled_node.begin(), ctx.coupled_node.end()));
        r.arr('i', "parea", ctx.coupled_pond_area);
        r.arr('i', "pband", ctx.coupled_pond_band);
        r.arr('i', "ndepth", ctx.nodes.depth);
        r.arr('i', "lat", ctx.nodes.lat_flow);
        r.arr('i', "loss", ctx.nodes.losses);
    }
    real(this, ctx);
    if (g_on) {
        r.arr('o', "surf", xnode_.new_surf_area);
        r.arr('o', "dqdh", xnode_.sumdqdh);
        r.arr('o', "conv", std::vector<double>(xnode_.converged.begin(), xnode_.converged.end()));
        r.arr('o', "inflow", ctx.nodes.inflow);
        r.arr('o', "outflow", ctx.nodes.outflow);
        r.flush();
    }
}

int main(int argc, char** argv) {
    if (argc < 2)
        return 2;
    const int max_steps = argc > 2 ? std::atoi(argv[2]) : -1;
    if (const char* e = std::getenv("TWOD_EXCH_QEVERY"))
        g_qevery = std::max(1L, std::atol(e));
    SWMM_Engine h = swmm_engine_create();
    g_eng = reinterpret_cast<SWMMEngine*>(h);
    int rc = swmm_engine_open(h, argv[1], "/dev/null", "/dev/null", nullptr);
    if (rc == 0)
        rc = swmm_engine_initialize(h);
    if (rc == 0)
        rc = swmm_engine_start(h, 0);
    std::printf("START %d active %d points %zu junctions %zu\n", rc, router().active_ ? 1 : 0,
                router().coupling_points_.size(), router().node_coupling_points_.size());
    if (rc != 0 || !router().active_) {
        std::fflush(stdout);
        std::_Exit(0);
    }
    g_on = true;
    double elapsed = 0.0;
    for (;;) {
        if (max_steps >= 0 && g_step >= max_steps)
            break;
        ++g_step;
        rc = swmm_engine_step(h, &elapsed);
        if (rc != 0 || elapsed <= 0.0)
            break;
    }
    g_on = false;
    // TQ: each node's tailwater record.
    {
        Rec r;
        r.head("TQ");
        std::vector<double> v;
        for (int n = 0; n < g_eng->ctx_.n_nodes(); ++n) {
            double a = -1, b = -2, c = -3;
            const bool ok = router().outfallExchange().tailwater(n, a, b, c);
            v.insert(v.end(), {ok ? 1.0 : 0.0, a, b, c});
        }
        r.arr('o', "tw", v);
        r.flush();
    }
    std::printf("END %d rc %d records %ld\n", g_step, rc, g_seq);
    std::fflush(stdout);
    std::_Exit(0);
}
