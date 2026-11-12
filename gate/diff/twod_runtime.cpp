// RuntimeCoupling.hpp differential (phase 2a wave 4b): the runtime surface sources and forcings,
// compiled from upstream's own header. A fixed pseudo-random sequence of operations (replace,
// clear, rebuild with out-of-range cells, prune, expire, nextTime, syncReceipts, active, get,
// rate) runs on both stores; after each operation the harness prints the operation's result
// and the whole store (every row with its receipt, the per-cell lists or channel slots, and
// find() for every key), doubles as hex bits. The Rust twin is
// crates/engine/examples/twod_runtime.rs (twod/coupling/runtime_coupling.rs); the outputs must
// be byte-identical.
// Build (gate/twodexch.sh): build2d's compile line; header-only, no library.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
#include <cinttypes>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include "2d/coupling/RuntimeCoupling.hpp"

using namespace openswmm::twoD;

static uint64_t g_state = 0x9E3779B97F4A7C15ull;
static uint64_t rnd() { // xorshift64*
    g_state ^= g_state >> 12;
    g_state ^= g_state << 25;
    g_state ^= g_state >> 27;
    return g_state * 0x2545F4914F6CDD1Dull;
}
static int ri(int n) { return static_cast<int>(rnd() % static_cast<uint64_t>(n)); }
static double rd() { return static_cast<double>(ri(4001) - 2000) / 16.0; }
static double runtil() { return ri(4) == 0 ? std::numeric_limits<double>::infinity() : static_cast<double>(ri(400)) / 4.0; }

static void hx(double v) {
    uint64_t u;
    std::memcpy(&u, &v, 8);
    std::printf(" %016" PRIx64, u);
}
static void vec(const char* n, const std::vector<double>& v) {
    std::printf(" %s%zu", n, v.size());
    for (double x : v) hx(x);
}
static void receipt(const RuntimeReceipt& r) {
    std::printf(" rc");
    hx(r.requested); hx(r.applied); hx(r.water_in); hx(r.water_out);
    hx(r.heat_requested); hx(r.heat_applied); hx(r.last_flow); hx(r.last_dt);
    vec("rm", r.requested_mass);
    vec("am", r.applied_mass);
}
static const char* kIds[] = {"a", "b", "c", "d"};

static void dump(const RuntimeSources& s) {
    for (const auto& r : s.rows) {
        std::printf("row %s %s %d", r.id.c_str(), r.owner.c_str(), r.cell);
        hx(r.flow); hx(r.heat);
        std::printf(" %d", r.boundary_edge);
        hx(r.boundary_head); hx(r.until);
        std::printf(" %d", r.enabled ? 1 : 0);
        vec("c", r.concentrations);
        vec("r", r.rates);
        receipt(r.receipt);
        std::printf("\n");
    }
    std::printf("cells %zu:", s.cells.size());
    for (const auto& c : s.cells) {
        std::printf(" [");
        for (auto k : c) std::printf(" %zu", k);
        std::printf(" ]");
    }
    std::printf("\nfind");
    for (const char* id : kIds)
        for (int cell = -1; cell < 9; ++cell) {
            const auto* r = s.find(id, cell);
            std::printf(" %td", r ? r - s.rows.data() : -1);
        }
    std::printf("\n");
}

static void dump(const RuntimeForcings& f) {
    for (const auto& r : f.rows) {
        std::printf("frow %d %d %d", r.channel, r.cell, r.mode);
        hx(r.rate); hx(r.until);
        std::printf(" %d %s", r.enabled ? 1 : 0, r.owner.c_str());
        vec("c", r.concentrations);
        receipt(r.receipt);
        std::printf("\n");
    }
    std::printf("fcells %zu:", f.cells.size());
    for (const auto& c : f.cells) std::printf(" %d,%d,%d,%d", c[0], c[1], c[2], c[3]);
    std::printf(" uniform %d,%d,%d,%d\nffind", f.uniform[0], f.uniform[1], f.uniform[2], f.uniform[3]);
    auto& m = const_cast<RuntimeForcings&>(f);
    for (int ch = 1; ch <= 4; ++ch)
        for (int cell = -1; cell < 8; ++cell) {
            const auto* r = m.find(ch, cell);
            std::printf(" %td", r ? r - f.rows.data() : -1);
        }
    std::printf("\n");
}

static RuntimeReceipt rreceipt() {
    RuntimeReceipt r;
    r.requested = rd(); r.applied = rd(); r.water_in = rd(); r.water_out = rd();
    r.heat_requested = rd(); r.heat_applied = rd(); r.last_flow = rd(); r.last_dt = rd();
    for (int k = ri(3); k > 0; --k) r.requested_mass.push_back(rd());
    for (int k = ri(3); k > 0; --k) r.applied_mass.push_back(rd());
    return r;
}

int main() {
    RuntimeSources s, snap;
    RuntimeForcings f;
    for (int step = 0; step < 6000; ++step) {
        const int op = ri(16);
        std::printf("op %d %d:", step, op);
        switch (op) {
        case 0: case 1: case 2: {
            RuntimeSource v;
            v.id = kIds[ri(4)];
            v.owner = ri(2) ? "own" : "";
            v.cell = ri(9);
            v.flow = rd(); v.heat = rd();
            v.boundary_edge = ri(3) - 1;
            if (ri(2)) v.boundary_head = rd();
            v.until = runtil();
            v.enabled = ri(5) != 0;
            for (int k = ri(3); k > 0; --k) { v.concentrations.push_back(rd()); v.rates.push_back(rd()); }
            if (ri(2)) v.receipt = rreceipt();
            std::printf(" replace %s %d\n", v.id.c_str(), v.cell);
            s.replace(std::move(v));
            break;
        }
        case 3: {
            const char* id = ri(5) == 4 ? "zz" : kIds[ri(4)];
            const int cell = ri(10) - 1;
            std::printf(" clear %s %d -> %d\n", id, cell, s.clear(id, cell) ? 1 : 0);
            break;
        }
        case 4: {
            const int count = 3 + ri(6);
            std::printf(" rebuild %d ->", count);
            try {
                s.rebuild(count);
                std::printf(" ok\n");
            } catch (const std::out_of_range& e) {
                std::printf(" out_of_range %s\n", e.what());
            }
            break;
        }
        case 5: std::printf(" prune\n"); s.prune(); break;
        case 6: {
            const double now = static_cast<double>(ri(400)) / 4.0;
            std::printf(" expire"); hx(now);
            std::printf(" -> %d\n", s.expire(now) ? 1 : 0);
            break;
        }
        case 7: {
            const double now = static_cast<double>(ri(400)) / 4.0;
            std::printf(" nextTime"); hx(now); std::printf(" ->"); hx(s.nextTime(now)); std::printf(" active %d\n", s.active() ? 1 : 0);
            break;
        }
        case 8: std::printf(" snapshot\n"); snap = s; break;
        case 9: std::printf(" syncReceipts\n"); s.syncReceipts(snap); break;
        case 10: case 11: {
            RuntimeForcing r;
            r.channel = 1 + ri(4);
            r.cell = ri(8) - 1;
            r.mode = ri(3);
            r.rate = rd();
            r.until = runtil();
            r.enabled = r.mode != 0;
            r.owner = ri(2) ? "own" : "";
            for (int k = ri(3); k > 0; --k) r.concentrations.push_back(rd());
            if (ri(2)) r.receipt = rreceipt();
            std::printf(" freplace %d %d\n", r.channel, r.cell);
            f.replace(std::move(r));
            break;
        }
        case 12: {
            const int count = 3 + ri(6);
            std::printf(" frebuild %d ->", count);
            try {
                f.rebuild(count);
                std::printf(" ok\n");
            } catch (const std::out_of_range& e) {
                std::printf(" out_of_range %s\n", e.what());
            }
            break;
        }
        case 13: {
            const double now = static_cast<double>(ri(400)) / 4.0;
            std::printf(" fexpire"); hx(now); f.expire(now); std::printf("\n");
            break;
        }
        case 14: {
            const double now = static_cast<double>(ri(400)) / 4.0;
            std::printf(" fnextTime"); hx(now); std::printf(" ->"); hx(f.nextTime(now)); std::printf(" active %d\n", f.active() ? 1 : 0);
            break;
        }
        default: {
            std::printf(" fget");
            for (int ch = 1; ch <= 4; ++ch)
                for (int cell = -1; cell < 9; ++cell) {
                    const auto* r = f.get(ch, cell);
                    std::printf(" %td", r ? r - f.rows.data() : -1);
                    hx(f.rate(ch, cell, 0.5 * cell));
                }
            std::printf("\n");
            break;
        }
        }
        dump(s);
        dump(f);
    }
    return 0;
}
