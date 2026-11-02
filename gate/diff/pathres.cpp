// Differential for core/PathResolver.cpp (openswmm::io) and the libstdc++ path members it
// uses, against the reference library. Reads "a<TAB>b" lines; prints every result.
// Build like ctxread (without the handler file); Rust side: examples/pathres.rs.
#include "core/PathResolver.hpp"
#include <filesystem>
#include <iostream>
#include <string>
namespace fs = std::filesystem;
static void p(const char* l, const std::string& s) { std::cout << l << " " << s.size() << ":" << s << "\n"; }
int main() {
    std::string line;
    while (std::getline(std::cin, line)) {
        auto t = line.find('\t');
        std::string a = line.substr(0, t), b = t == std::string::npos ? "" : line.substr(t + 1);
        p("parentDir", openswmm::io::parentDir(a));
        p("normalise", openswmm::io::normaliseSeparators(a));
        p("isAbs", openswmm::io::isAbsolutePath(a) ? "1" : "0");
        p("resolve", openswmm::io::resolveRelative(a, b));
        auto r = openswmm::io::makeRelative(a, b);
        p("mkrel.path", r.path); p("mkrel.class", std::to_string((int)r.classification));
        p("mkrel.up", std::to_string(r.up_levels)); p("mkrel.warn", r.warning);
        auto r2 = openswmm::io::makeRelative(a, b, 1);
        p("mkrel1.path", r2.path); p("mkrel1.warn", r2.warning);
        fs::path pa(a), pb(b);
        p("norm", pa.lexically_normal().generic_string());
        p("normnative", pa.lexically_normal().native());
        p("gen", pa.generic_string());
        p("join", (pa / pb).native());
        p("joinnorm", (pa / pb).lexically_normal().native());
        p("rel", pa.lexically_relative(pb).native());
        p("prox", pa.lexically_proximate(pb).native());
        p("parent", pa.parent_path().native());
        fs::path rf = pa; rf.remove_filename(); p("rmfn", rf.native());
        p("fname", pa.filename().native());
        p("hasfn", pa.has_filename() ? "1" : "0");
        p("hasrel", pa.has_relative_path() ? "1" : "0");
    }
}
