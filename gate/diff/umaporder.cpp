// Iteration order of libstdc++'s std::unordered_map under insertion, for the two hashes the
// translation emulates (core::std_unordered_map): std::hash<int> and openswmm::CiHash.
// 200 pseudo-random insertion sequences each, up to ~300 keys (several rehashes), printed as
// the insertion index of each element in begin()..end() order. Rust twin:
// crates/engine/examples/umaporder.rs. Build: g++ -std=c++20 -O1 -I<upstream>/src/engine
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>
#include "core/StringCase.hpp"
int main() {
    unsigned s = 12345;
    auto rnd = [&] { s = s * 1103515245u + 12345u; return (s >> 8) & 0xffffff; };
    for (int trial = 0; trial < 200; ++trial) {
        int n = trial % 7 == 0 ? (int)(rnd() % 300) : (int)(rnd() % 40);
        std::unordered_map<int, int> mi;
        std::vector<int> keys;
        for (int i = 0; i < n; ++i) {
            int k = (int)(rnd() % (trial % 3 ? 60 : 100000));
            if (mi.find(k) == mi.end()) { mi[k] = (int)keys.size(); keys.push_back(k); }
        }
        std::printf("I %zu:", keys.size());
        for (auto& kv : mi) std::printf(" %d", kv.second);
        std::printf("\n");
        std::unordered_map<std::string, int, openswmm::CiHash, openswmm::CiEqual> ms;
        std::vector<std::string> sk;
        for (int i = 0; i < n; ++i) {
            std::string k = "uh" + std::to_string(rnd() % 500);
            if (i % 5 == 0) k[0] = 'U';
            if (ms.find(k) == ms.end()) { ms[k] = (int)sk.size(); sk.push_back(k); }
        }
        std::printf("S %zu:", sk.size());
        for (auto& kv : ms) std::printf(" %d", kv.second);
        std::printf("\n");
    }
}
