// Parse a deck with upstream's own InputReader and section handlers, taken from the
// reference engine library itself (libopenswmm.engine.so exports them), then dump the
// context. argv[2] (optional) is a comma-separated list of section tags to enable, or ALL.
// argv[3] (optional) "resolve": after reading, run upstream's
// openswmm::input::resolve_cross_references(ctx) from the library, then dump.
#include "input/InputReader.hpp"
#include "input/SectionRegistry.hpp"
#include "input/PostParseResolver.hpp"
#include "ctxdump.hpp"
#include <string>
#include <set>
namespace openswmm::input { void register_reference_handlers(SectionRegistry&, const std::set<std::string>&, bool all); }
int main(int argc, char** argv) {
    using namespace openswmm;
    std::set<std::string> on; bool all = argc > 2 && std::string(argv[2]) == "ALL";
    if (argc > 2 && !all) { std::string s = argv[2]; size_t p = 0; while (p <= s.size()) { size_t q = s.find(',', p); if (q == std::string::npos) q = s.size(); if (q > p) on.insert(s.substr(p, q - p)); p = q + 1; } }
    const bool resolve = argc > 3 && std::string(argv[3]) == "resolve";
    input::SectionRegistry reg;
    input::register_reference_handlers(reg, on, all);
    input::InputReader rd(reg);
    SimulationContext ctx;
    rd.read(argv[1], ctx);
    if (resolve)
        input::resolve_cross_references(ctx);
    ctxdump::dump(ctx);
}
