#!/bin/bash
# Phase 2a wave 4b: the run-time 1D-2D exchange differential. Builds gate/diff/twod_exch.cpp
# with the 2D reference's compile line (build2d's flags.make: Release, g++ 13,
# -ffp-contract=off -mavx2 -mfma, OpenMP, OPENSWMM_HAS_2D=1) against build2d's library, linked
# with -Wl,-E so the library's PLT calls to the exchange functions bind to the harness's
# interposers, and the Rust twin (examples/twod_exch.rs, features "twod gate", into target2d).
# Each deck runs to its end on the reference; every exchange call's inputs and upstream's
# outputs stream into the twin, which runs the translation on the same inputs and compares the
# outputs bit for bit (see the two files' headers). First, gate/diff/twod_runtime.cpp against
# examples/twod_runtime.rs: RuntimeCoupling.hpp on a fixed random sequence of operations.
#
# Usage: gate/twodexch.sh [deck.inp ...]   (paths are made absolute)
#   default decks: gate/decks/twodc, twodx, twodi, twod, upstream examples/2d_complete_example.inp;
#   ALL2D=1 adds every upstream deck with a [2D_OPTIONS] section.
#   NOBUILD=1 skips the builds; KEEP=file keeps the per-deck results; SHOW=n differing records
#   shown per deck (default 3).
# Results per deck: identical, differ, nostart (the reference stops before stepping, or the mesh
# is inactive), stub (the translation refuses a phase-2b path in initialize/start), noresult.
# Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
set -e
H=$(cd "$(dirname "$0")/.." && pwd); cd "$H"
U=/root/work/swmm6src/upstream; B=/root/work/swmm6src/build2d
CXX="g++ -ffp-contract=off -O3 -DNDEBUG -std=c++20 -flto=auto -fno-fat-lto-objects -fPIC -fvisibility=default -mavx2 -mfma -fopenmp"
DEF="-DOPENSWMM_HAS_2D=1 -DSWMM_OS_LINUX -DSWMM_USE_OPENMP"
INC="-I$U/include -I$U/src/engine -I$B/include -Igate/diff"
if [ -z "$NOBUILD" ]; then
  echo "g++ twod_exch twod_runtime"
  nice $CXX $DEF $INC gate/diff/twod_exch.cpp -o gate/diff/twod_exch $B/src/engine/libopenswmm.engine.so -Wl,-rpath,$B/src/engine -Wl,-E -ldl
  nice $CXX $DEF $INC gate/diff/twod_runtime.cpp -o gate/diff/twod_runtime
  CARGO_BUILD_JOBS=${CARGO_BUILD_JOBS:-1} CARGO_TARGET_DIR="$H/target2d" nice cargo build --release --features "twod gate" --example twod_exch --example twod_runtime
fi
C_BIN="$H/gate/diff/twod_exch"; R_BIN="$H/target2d/release/examples/twod_exch"
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT
# RuntimeCoupling.hpp: the operation-sequence differential (no deck)
"$H/gate/diff/twod_runtime" > "$W/rt_c.txt"; "$H/target2d/release/examples/twod_runtime" > "$W/rt_r.txt"
if cmp -s "$W/rt_c.txt" "$W/rt_r.txt"; then echo "twod_runtime: identical, $(wc -l < "$W/rt_c.txt") lines"
else echo "twod_runtime: DIFFER at $(cmp "$W/rt_c.txt" "$W/rt_r.txt" | head -1)"; fi
if [ $# -gt 0 ]; then for d in "$@"; do realpath "$d"; done > "$W/decks.txt"
else
  { find "$H/gate/decks/twodc" "$H/gate/decks/twodx" "$H/gate/decks/twodi" "$H/gate/decks/twod" -name '*.inp' 2>/dev/null
    echo "$U/examples/2d_complete_example.inp"
    [ -n "$ALL2D" ] && grep -rl --include='*.inp' '^\[2D_OPTIONS\]' "$U/tests" "$U/examples" "$U/docs" 2>/dev/null | grep -v '^'"$U"'/examples/2d_complete_example.inp$'
  } | sort -u > "$W/decks.txt"
fi
: > "$W/res.txt"; mkdir -p "$W/bak"
while IFS= read -r d; do
  # back up any files the deck SAVEs to, and restore them after the run
  python3 "$H/gate/savefiles.py" "$d" > "$W/saves.txt" 2>/dev/null || : > "$W/saves.txt"; k=0
  while IFS= read -r f; do k=$((k+1)); rm -f "$W/bak/$k"; [ -f "$f" ] && cp -p "$f" "$W/bak/$k"; done < "$W/saves.txt"
  ( cd "$(dirname "$d")" && timeout 3600 nice "$C_BIN" "$d" 2>/dev/null | TWOD_EXCH_SHOW=${SHOW:-3} timeout 3600 nice "$R_BIN" "$d" > "$W/r.txt" 2> "$W/err.txt" ) 2>/dev/null || true
  k=0; while IFS= read -r f; do k=$((k+1)); if [ -f "$W/bak/$k" ]; then cp -p "$W/bak/$k" "$f"; else rm -f "$f"; fi; done < "$W/saves.txt"
  res=$(head -1 "$W/r.txt")
  if grep -q panicked "$W/err.txt"; then res="panic $(grep -m1 -A1 panicked "$W/err.txt" | tail -1 | head -c 160)"
  elif [ -z "$res" ]; then res="noresult"; fi
  echo "$res $d" >> "$W/res.txt"
  grep '^cov ' "$W/r.txt" >> "$W/cov.txt" || true
  case "$res" in identical*|nostart*|stub*) ;; *) echo "== $d"; sed -n '2,200p' "$W/r.txt" | grep -v '^cov ' | head -60 ;; esac
done < "$W/decks.txt"
for k in identical differ nostart stub terminated noend noresult panic; do
  printf '%s %s\n' "$k" "$(grep -c "^$k " "$W/res.txt")"
done | tr '\n' ' '; echo "($(wc -l < "$W/decks.txt") decks)"
# total records per kind over the identical decks, and the coverage of the paths
grep '^identical ' "$W/res.txt" | sed 's/^identical [0-9]* (\([^)]*\)).*/\1/' | tr ' ' '\n' | paste - - | \
  awk '{n[$1]+=$2} END {printf "records:"; for (k in n) printf " %s %d", k, n[k]; printf "\n"}'
sed 's/^cov //' "$W/cov.txt" 2>/dev/null | tr ' ' '\n' | grep '=' | \
  awk -F= '{n[$1]+=$2; d[$1]++} END {printf "coverage (records, decks):"; for (k in n) printf " %s %d %d,", k, n[k], d[k]; printf "\n"}'
grep -v '^identical \|^nostart \|^stub ' "$W/res.txt" | head -20 || true
[ -n "$KEEP" ] && cp "$W/res.txt" "$KEEP"
exit 0
