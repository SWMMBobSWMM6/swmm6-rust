// Open (and optionally initialize) a deck in the reference engine through its C API, then
// dump the engine's own SimulationContext. The C API handle is the SWMMEngine object, and
// context() is an inline accessor, so this harness must be compiled with the engine's own
// defines (layout must match): -DNDEBUG -DSWMM_OS_LINUX -DSWMM_USE_OPENMP -fopenmp.
// Usage: ctxopen deck.inp [open|init]
#include "openswmm/engine/openswmm_engine.h"
#include "core/SWMMEngine.hpp"
#include "ctxdump.hpp"
#include <cstring>
int main(int argc, char** argv) {
    SWMM_Engine e = swmm_engine_create();
    int rc = swmm_engine_open(e, argv[1], "/dev/null", "/dev/null", nullptr);
    std::printf("open_rc %d\n", rc);
    if (rc == 0 && argc > 2 && std::strcmp(argv[2], "init") == 0)
        std::printf("init_rc %d\n", swmm_engine_initialize(e));
    ctxdump::dump(static_cast<openswmm::SWMMEngine*>(e)->context());
    return 0;
}
