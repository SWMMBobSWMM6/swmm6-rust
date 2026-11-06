// Differential harness for the two C library hazards of the FV step (6 Nov 2026):
//   * NativeGmres.cpp's `dot`, which accumulates in x87 extended precision (long double), copied
//     here verbatim, against core::ext80::dot_ext80;
//   * glibc's cbrt (FvKernels.hpp frictionUpdate's std::cbrt) against core::cmath::glibc_cbrt.
// Both sides draw the same inputs from one xorshift64* stream: vectors of 1..64 doubles over
// wide exponent ranges, near-cancelling pairs, subnormals and zeros; and single values for
// cbrt. Prints every result as hex bits. Rust side: crates/engine/examples/ext80dot.rs.
// Build: g++ -std=c++20 -O3 -mavx2 -mfma -ffp-contract=off gate/diff/ext80dot.cpp -o gate/diff/ext80dot
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>
static uint64_t st = 0x9e3779b97f4a7c15ull;
static uint64_t nx() { st ^= st >> 12; st ^= st << 25; st ^= st >> 27; return st * 0x2545f4914f6cdd1dull; }
static double mk(int lo, int hi) {
  uint64_t r = nx();
  int e = lo + (int)(r % (uint64_t)(hi - lo + 1));
  uint64_t b = ((r >> 63) << 63) | ((uint64_t)e << 52) | (nx() & ((1ull << 52) - 1));
  double d; std::memcpy(&d, &b, 8); return d;
}
static double dot(const double* a, const double* b, std::size_t n) {
  long double sum = 0.0;
  for (std::size_t i = 0; i < n; ++i)
    sum += static_cast<long double>(a[i]) * b[i];
  return static_cast<double>(sum);
}
static void hx(double v) { uint64_t u; std::memcpy(&u, &v, 8); std::printf("%016llx\n", (unsigned long long)u); }
int main() {
  std::vector<double> a, b;
  for (int t = 0; t < 300000; ++t) {
    const int kind = (int)(nx() % 6);
    const std::size_t n = 1 + nx() % 64;
    a.assign(n, 0.0); b.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
      switch (kind) {
      case 0: a[i] = mk(1023 - 30, 1023 + 30); b[i] = mk(1023 - 30, 1023 + 30); break;
      case 1: a[i] = mk(1, 2046); b[i] = mk(1023 - 200, 1023 + 200); break;
      case 2: a[i] = mk(1023 - 5, 1023 + 5); b[i] = (i % 2) ? -a[i - 1] * (1.0 + std::ldexp((double)(nx() % 8), -50)) : mk(1023 - 5, 1023 + 5); break;
      case 3: { uint64_t s = nx() & ((1ull << 52) - 1); std::memcpy(&a[i], &s, 8); b[i] = mk(1023 + 500, 1023 + 1000); } break;
      case 4: a[i] = (nx() % 4 == 0) ? 0.0 : mk(1023 - 60, 1023 + 60); b[i] = (nx() % 4 == 0) ? -0.0 : mk(1023 - 60, 1023 + 60); break;
      default: a[i] = mk(1, 400); b[i] = mk(1, 400); break;
      }
    }
    hx(dot(a.data(), b.data(), n));
  }
  for (int t = 0; t < 300000; ++t) {
    double x = (t % 3 == 0) ? mk(1, 2046) : mk(1023 - 40, 1023 + 40);
    if (t % 97 == 0) { uint64_t s = nx() & ((1ull << 52) - 1); std::memcpy(&x, &s, 8); }
    hx(std::cbrt(x));
  }
  return 0;
}
