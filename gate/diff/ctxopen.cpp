// Open (and optionally initialize) a deck in the reference engine through its C API, then
// dump the engine's own SimulationContext and the SWMMEngine members open() sets. The C API
// handle is the SWMMEngine object, and context() is an inline accessor, so this harness must
// be compiled with the engine's own defines (layout must match): -DNDEBUG -DSWMM_OS_LINUX
// -DSWMM_USE_OPENMP -fopenmp. The Rust twin is crates/engine/examples/ctxopen.rs.
// Usage: ctxopen deck.inp [open|lenient|metadata|init]
//   open      swmm_engine_open(deck, "/dev/null", "/dev/null", NULL)   (the default)
//   lenient   the same after swmm_engine_set_lenient_open(e, 1)
//   metadata  swmm_engine_open_metadata(deck)
//   init      open, then swmm_engine_initialize
#include "openswmm/engine/openswmm_engine.h"
#include "core/SWMMEngine.hpp"
#define CTXDUMP_ENGINE
#include "ctxdump.hpp"
#include <cstring>
int main(int argc, char** argv) {
    const char* mode = argc > 2 ? argv[2] : "open";
    SWMM_Engine e = swmm_engine_create();
    if (std::strcmp(mode, "lenient") == 0)
        swmm_engine_set_lenient_open(e, 1);
    int rc = std::strcmp(mode, "metadata") == 0 ? swmm_engine_open_metadata(e, argv[1])
                                                : swmm_engine_open(e, argv[1], "/dev/null", "/dev/null", nullptr);
    std::printf("open_rc %d\n", rc);
    if (rc == 0 && std::strcmp(mode, "init") == 0)
        std::printf("init_rc %d\n", swmm_engine_initialize(e));
    auto* eng = static_cast<openswmm::SWMMEngine*>(e);
    ctxdump::dump(eng->context());
    ctxdump::engdump(*eng);
    return 0;
}
