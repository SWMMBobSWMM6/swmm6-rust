#!/bin/bash
# Phase 2a wave 2: the G4 differentials of the 2D leaf numerics, in one command. Builds every
# gate/diff/twod_*.cpp harness from upstream's own 2D headers with the 2D reference build's exact
# compile line (build2d's flags.make: Release, g++ 13, -ffp-contract=off -mavx2 -mfma, OpenMP,
# OPENSWMM_HAS_2D=1; twod_momentum links build2d's library and calls its compiled
# MomentumPhysics2D::rotate), builds the Rust twins (examples/twod_*.rs, feature twod), runs both
# sides and compares their stdout byte for byte. Where `cmp` finds a difference, gate/twodcmp.py
# classifies it: the harness passes only if every differing line differs in nothing but the sign
# of a NaN result on a line with a non-finite input (the x86 / LLVM operand-order hazard of
# PLAN.md), and the count of such lines is reported.
#
# Usage: gate/twodkern.sh [harness ...]   (default: all; names without the twod_ prefix work)
#   OUT=dir keeps the outputs in dir (default: a temporary directory, removed at the end).
# The 1D harness builds (gate/build.sh) are unchanged.
# Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
set -e
H=$(cd "$(dirname "$0")/.." && pwd); cd "$H"
U=/root/work/swmm6src/upstream; B=/root/work/swmm6src/build2d
CXX="g++ -ffp-contract=off -O3 -DNDEBUG -std=c++20 -flto=auto -fno-fat-lto-objects -fPIC -fvisibility=default -mavx2 -mfma -fopenmp"
DEF="-DOPENSWMM_HAS_2D=1 -DSWMM_OS_LINUX -DSWMM_USE_OPENMP"
INC="-I$U/include -I$U/src/engine -I$B/include -Igate/diff"
ALL="vfr quadvfr inertial swe diffusive sinks momentum exchange pond"
NAMES=""
for h in ${@:-$ALL}; do NAMES="$NAMES twod_${h#twod_}"; done

for n in $NAMES; do
  [ -f gate/diff/$n.cpp ] || { echo "unknown harness $n"; exit 2; }
  EXTRA=""
  [ $n = twod_momentum ] && EXTRA="$B/src/engine/libopenswmm.engine.so -Wl,-rpath,$B/src/engine"
  echo "g++ $n"
  $CXX $DEF $INC gate/diff/$n.cpp -o gate/diff/$n $EXTRA
done
EX=""; for n in $NAMES; do EX="$EX --example $n"; done
CARGO_BUILD_JOBS=${CARGO_BUILD_JOBS:-1} cargo build --release --features twod $EX

if [ -n "$OUT" ]; then mkdir -p "$OUT"; T=$OUT; else T=$(mktemp -d); trap 'rm -rf "$T"' EXIT; fi
fail=0
printf '%-16s %10s  %s\n' harness lines result
for n in $NAMES; do
  gate/diff/$n > "$T/$n.cpp.txt" 2>/dev/null
  target/release/examples/$n > "$T/$n.rs.txt" 2>/dev/null
  lines=$(wc -l < "$T/$n.cpp.txt")
  if cmp -s "$T/$n.cpp.txt" "$T/$n.rs.txt"; then
    printf '%-16s %10s  %s\n' $n $lines identical
  elif r=$(python3 gate/twodcmp.py "$T/$n.cpp.txt" "$T/$n.rs.txt" 0); then
    printf '%-16s %10s  %s\n' $n $lines "$(echo "$r" | tail -1)"
  else
    printf '%-16s %10s  %s\n' $n $lines DIFFER
    python3 gate/twodcmp.py "$T/$n.cpp.txt" "$T/$n.rs.txt" 3 || true
    fail=1
  fi
done
exit $fail
