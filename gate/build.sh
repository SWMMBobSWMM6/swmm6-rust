#!/bin/bash
# Build the gate harnesses: the C++ ones that link the reference library (they include the
# generated gate/diff/ctxdump.hpp, rtdump.hpp, hydrodump.hpp, so rebuild them after
# tools/ctxdump.py runs; ctxstep also includes lifedump.hpp), then the translation's examples. The harnesses read the engine's
# own classes, so they use the engine's layout defines (-DNDEBUG -DSWMM_OS_LINUX
# -DSWMM_USE_OPENMP -fopenmp). Usage: gate/build.sh [harness ...]   (default: all)
# The differentials compiled from upstream source (tokdump, xsect, hydfuncs, ...) carry their
# own build lines in their headers and do not depend on the generated dumpers.
set -e
H=$(cd "$(dirname "$0")/.." && pwd); cd "$H"
U=/root/work/swmm6src/upstream; B=/root/work/swmm6src/build
CXX="g++ -std=c++20 -O3 -ffp-contract=off -mavx2 -mfma -DNDEBUG -DSWMM_OS_LINUX -DSWMM_USE_OPENMP -fopenmp"
INC="-I$U/src/engine -I$U/include -I$B/include -Igate/diff"
LIB="$B/src/engine/libopenswmm.engine.so -Wl,-rpath,$B/src/engine"
build() { # name, sources...
  local n=$1; shift
  echo "g++ $n"; $CXX $INC "$@" -o "gate/diff/$n" $LIB
}
ALL="ctxread ctxopen msxiq pathres inithyd inithydro initqual ctxstep routestep outgate umaporder steptrace"
# The 2D harnesses (not in ALL; `gate/build.sh ctxopen2d`): compiled with the 2D layout define
# against the 2D reference library (build2d, OPENSWMM_HAS_2D=1), and their Rust twins built with
# the twod feature into target2d, so the 1D examples in target/ are not rebuilt.
B2=/root/work/swmm6src/build2d
INC2D="-I$U/src/engine -I$U/include -I$B2/include -Igate/diff"
LIB2D="$B2/src/engine/libopenswmm.engine.so -Wl,-rpath,$B2/src/engine"
build2d() { # name, sources...
  local n=$1; shift
  echo "g++ $n (2D)"; $CXX -DOPENSWMM_HAS_2D=1 $INC2D "$@" -o "gate/diff/$n" $LIB2D
}
ONED=0; TWOD=0
for h in ${@:-$ALL}; do
  case $h in
    ctxopen2d) TWOD=1 ;;
    *) ONED=1 ;;
  esac
  case $h in
    ctxopen2d) build2d ctxopen2d gate/diff/ctxopen2d.cpp ;;
    ctxread) build ctxread gate/diff/ctxread.cpp gate/diff/ctxread_handlers.cpp ;;
    msxiq) build msxiq gate/diff/msxiq.cpp gate/diff/ctxread_handlers.cpp ;;
    ctxopen|pathres|inithyd|inithydro|initqual|ctxstep) build $h gate/diff/$h.cpp ;;
    # -Wl,-E exports the harness's own SWMMEngine::stepRouting, so the library's PLT call
    # from step() binds to it (see routestep.cpp)
    routestep) build routestep gate/diff/routestep.cpp -Wl,-E -ldl ;;
    # outgate defines SWMMEngine::postOutputSnapshot the same way (see outgate.cpp)
    outgate) build outgate gate/diff/outgate.cpp -Wl,-E -ldl ;;
    umaporder) echo "g++ umaporder"; g++ -std=c++20 -O1 -I$U/src/engine gate/diff/umaporder.cpp -o gate/diff/umaporder ;;
    steptrace) echo "cc steptrace"; cc -O2 gate/steptrace.c -I$U/include/openswmm/engine -I$B/include -o gate/steptrace $LIB ;;
    *) echo "unknown harness $h"; exit 2 ;;
  esac
done
if [ $ONED = 1 ]; then cargo build --release --features gate --examples; fi
if [ $TWOD = 1 ]; then CARGO_TARGET_DIR="$H/target2d" cargo build --release --features twod --example ctxopen2d; fi
