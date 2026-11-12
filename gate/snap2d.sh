#!/bin/bash
# The 2D output-side differential (phase 2a wave 5): each deck through gate/diff/snap2d (the 2D
# reference library, its fillSurfaceSnapshot interposed) and examples/snap2d (the translation,
# twod + gate, target2d), in the deck's own folder with its SAVE files backed up and restored
# (gate/withsaves.sh), comparing stdout: the 2D surface block of every report period's snapshot
# (hashed per array), ctx.mass_balance_2d after end() (SurfaceRouter2D::finalize) and the exit
# code the CLI would return. A deck where the translation stops at code not translated yet
# ("STUB <what>", everything printed before it identical) counts apart.
#   gate/snap2d.sh LIST [OUTDIR]   (LIST: one deck path per line, absolute)
# Build first: gate/build.sh snap2d. SNAP2D_TIMEOUT=s sets the per-side limit (default 1800).
# Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
H=$(cd "$(dirname "$0")/.." && pwd)
LIST=$1; W=${2:-$(mktemp -d)}; mkdir -p "$W"; W=$(cd "$W" && pwd)   # absolute: the runs are in the decks' folders
C_BIN="$H/gate/diff/snap2d"; R_BIN="${SNAP2D_EXAMPLES:-$H/target2d/release/examples}/snap2d"
T=${SNAP2D_TIMEOUT:-1800}
n=0; same=0; differ=0; stub=0; snaps=0
while IFS= read -r d; do
  [ -z "$d" ] && continue
  n=$((n+1))
  OMP="env OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OMP_DYNAMIC=false"
  "$H/gate/withsaves.sh" "$d" timeout "$T" nice $OMP "$C_BIN" "$d" "$W/c.rpt" "$W/c.out" > "$W/c.txt" 2> /dev/null
  "$H/gate/withsaves.sh" "$d" timeout "$T" nice $OMP "$R_BIN" "$d" "$W/r.rpt" "$W/r.out" > "$W/r.txt" 2> /dev/null
  if cmp -s "$W/c.txt" "$W/r.txt"; then
    same=$((same+1)); snaps=$((snaps + $(grep -c '^SNAP ' "$W/c.txt")))
  elif tail -1 "$W/r.txt" | grep -q '^STUB ' && [ "$(head -n -1 "$W/r.txt")" = "$(head -n "$(($(wc -l < "$W/r.txt") - 1))" "$W/c.txt")" ]; then
    stub=$((stub+1)); echo "stub $d: $(tail -1 "$W/r.txt")"
  else
    differ=$((differ+1)); echo "differ $d: $(diff "$W/c.txt" "$W/r.txt" | head -3 | tr '\n' ' ')"
    [ $differ = 1 ] && cp "$W/c.txt" "$W/first_c.txt" && cp "$W/r.txt" "$W/first_r.txt"
  fi
done < "$LIST"
echo "snap2d: $n decks: $same identical ($snaps report-period snapshots), $differ differ, $stub stopped at a stub"
