// Open a deck in the reference engine through its C API, run the private
// SWMMEngine::initHydraulics(), then dump the engine's SimulationContext, the members
// open() sets, and everything initHydraulics builds (the Router and its solvers, the
// inflow / RDII / structure / exfiltration / inlet / interface-file solvers, the control
// rules). The C API handle IS the SWMMEngine; `#define private public` (after every
// standard header, so their access specifiers are untouched) lets the harness call the
// private method and read the private members. Compile with the engine's own layout
// defines: -DNDEBUG -DSWMM_OS_LINUX -DSWMM_USE_OPENMP -fopenmp.
// The Rust twin is crates/engine/examples/inithyd.rs.
// Usage: inithyd deck.inp
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
#include <fstream>
#include <functional>
#include <initializer_list>
#include <ios>
#include <iosfwd>
#include <istream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <new>
#include <numbers>
#include <optional>
#include <ostream>
#include <queue>
#include <ratio>
#include <semaphore>
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
#include <utility>
#include <variant>
#include <vector>
#include <omp.h>
#define private public
#define protected public
#include "openswmm/engine/openswmm_engine.h"
#include "core/SWMMEngine.hpp"
#include "hydraulics/fv/ExplicitFvSolver.hpp"
#define CTXDUMP_ENGINE
#include "ctxdump.hpp"
#include "rtdump.hpp"
int main(int argc, char** argv) {
    if (argc < 2)
        return 2;
    SWMM_Engine e = swmm_engine_create();
    int rc = swmm_engine_open(e, argv[1], "/dev/null", "/dev/null", nullptr);
    std::printf("open_rc %d\n", rc);
    auto* eng = static_cast<openswmm::SWMMEngine*>(e);
    if (rc == 0)
        eng->initHydraulics();
    ctxdump::dump(eng->context());
    ctxdump::engdump(*eng);
    ctxdump::rtdump(*eng);
    return 0;
}
