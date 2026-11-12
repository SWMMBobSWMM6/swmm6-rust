// Shared input generator and printer for the 2D kernel differentials (gate/diff/twod_*.cpp).
// Its Rust twin is crates/engine/examples/common/twodgen.rs; both draw the same doubles from
// one xorshift64* stream with integer and exact floating-point operations only, so both sides
// see identical input bits. Every output line carries its inputs and results as hex bits.
// Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
#pragma once
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

namespace twodgen {

// An opaque copy: the compiler cannot see the value, so it cannot constant-fold a kernel (or a
// libm call, which GCC folds with MPFR's correctly rounded result, not glibc's) on grid inputs.
inline double opq(double x) {
    asm volatile("" : "+x"(x));
    return x;
}

inline double bitsd(uint64_t b) { double d; std::memcpy(&d, &b, 8); return d; }
inline uint64_t dbits(double d) { uint64_t b; std::memcpy(&b, &d, 8); return b; }

const double INF = std::numeric_limits<double>::infinity();
const double QNAN = bitsd(0x7ff8000000000000ull);  // +NaN
const double NQNAN = bitsd(0xfff8000000000000ull); // -NaN (the x86 default NaN)
const double DMIN = bitsd(1);                      // smallest subnormal
const double NMIN = std::numeric_limits<double>::min();
const double DMAX = std::numeric_limits<double>::max();

// The generator's own arithmetic on two drawn values can meet two NaNs of opposite sign, and then
// g++ and rustc may return different ones (PLAN.md, hazards). Every derived input goes through
// canon(), which makes any NaN +NaN, so both sides always see the same input bits.
inline double canon(double x) { return std::isnan(x) ? QNAN : x; }

struct Gen {
    uint64_t s;
    explicit Gen(uint64_t seed) : s(seed) {}
    uint64_t nx() { s ^= s >> 12; s ^= s << 25; s ^= s >> 27; return s * 0x2545f4914f6cdd1dull; }
    uint64_t below(uint64_t n) { return nx() % n; }
    double u01() { return static_cast<double>(nx() >> 11) * 0x1.0p-53; }
    double uni(double lo, double hi) { return canon(lo + (hi - lo) * u01()); }
    // random sign, biased exponent in [lo, hi], random significand
    double ebits(int lo, int hi) {
        const uint64_t r = nx();
        const int e = lo + static_cast<int>(r % static_cast<uint64_t>(hi - lo + 1));
        return bitsd(((r >> 63) << 63) | (static_cast<uint64_t>(e) << 52) | (nx() & ((1ull << 52) - 1)));
    }
    double pos_ebits(int lo, int hi) { return std::fabs(ebits(lo, hi)); }
    double pick(const std::vector<double>& v) { return v[below(v.size())]; }
    // x moved by -4..4 ulps (nextafter)
    double nudge(double x) {
        const int k = static_cast<int>(below(9)) - 4;
        for (int i = 0; i < k; ++i) x = std::nextafter(x, INF);
        for (int i = 0; i > k; --i) x = std::nextafter(x, -INF);
        return x;
    }
    double nonfinite() {
        switch (below(4)) {
        case 0: return INF;
        case 1: return -INF;
        case 2: return QNAN;
        default: return NQNAN;
        }
    }
    // One argument value: a special (maybe nudged by a few ulps), a uniform draw over the
    // physical range [lo, hi], a random magnitude near unity or anywhere in the finite range,
    // and (1 in 64) a non-finite value.
    double draw(const std::vector<double>& sp, double lo, double hi) {
        if (below(64) == 0) return nonfinite();
        switch (below(10)) {
        case 0: case 1: case 2: return pick(sp);
        case 3: case 4: return nudge(pick(sp));
        case 5: case 6: case 7: return uni(lo, hi);
        case 8: return ebits(1023 - 45, 1023 + 12);
        default: return ebits(1, 2046);
        }
    }
};

// ---------------------------------------------------------------- output
struct Out {
    std::string buf;
    long lines = 0;
    void tag(const char* t) { buf += t; buf += ' '; }
    void hx(double v) {
        char b[20];
        std::snprintf(b, sizeof b, "%016llx ", static_cast<unsigned long long>(dbits(v)));
        buf += b;
    }
    void hi(long long v) { buf += std::to_string(v); buf += ' '; }
    void sep() { buf += "| "; } // inputs | results
    void nl() {
        buf += '\n';
        ++lines;
        if (buf.size() > (1u << 22)) flush();
    }
    void flush() { std::fwrite(buf.data(), 1, buf.size(), stdout); buf.clear(); }
};

// The cartesian product of index ranges, for grids: calls f(idx) for every combination.
template <class F>
void grid(const std::vector<std::size_t>& sizes, F f) {
    std::vector<std::size_t> idx(sizes.size(), 0);
    for (std::size_t s : sizes) if (s == 0) return;
    while (true) {
        f(idx);
        std::size_t k = 0;
        while (k < sizes.size()) {
            if (++idx[k] < sizes[k]) break;
            idx[k] = 0;
            ++k;
        }
        if (k == sizes.size()) return;
    }
}

} // namespace twodgen
