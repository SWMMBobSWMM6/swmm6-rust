// The [INITIAL_QUALITY] MSX branches of the post-parse resolver, which need a configured
// reactions component (set at engine init, so ctxread never reaches them): parse a deck,
// mark reactions configured with species CL2 and NH3 and two .rxn per-element seeds (node
// 0 NH3 0.5, matching; link 0 CL2 9.0, conflicting), resolve with upstream's library, dump.
#include "input/InputReader.hpp"
#include "input/SectionRegistry.hpp"
#include "input/PostParseResolver.hpp"
#include "ctxdump.hpp"
#include <set>
#include <string>
namespace openswmm::input { void register_reference_handlers(SectionRegistry&, const std::set<std::string>&, bool all); }
int main(int argc, char** argv) {
    using namespace openswmm;
    input::SectionRegistry reg;
    input::register_reference_handlers(reg, {}, true);
    input::InputReader rd(reg);
    SimulationContext ctx;
    rd.read(argv[1], ctx);
    ctx.reactions.configured = true;
    ctx.reactions.species_name = {"CL2", "NH3"};
    ctx.reactions.init_elem_is_link = {0, 1};
    ctx.reactions.init_elem_idx = {0, 0};
    ctx.reactions.init_elem_species = {1, 0};
    ctx.reactions.init_elem_value = {0.5, 9.0};
    input::resolve_cross_references(ctx);
    ctxdump::dump(ctx);
}
