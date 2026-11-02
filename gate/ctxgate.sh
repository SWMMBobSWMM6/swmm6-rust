#!/bin/bash
# Parse every corpus and upstream test deck with upstream (reference library) and with the
# translation, dump the whole context, compare. Usage: gate/ctxgate.sh [TAGS|ALL] [resolve]
# With "resolve", both sides run resolve_cross_references after reading, then dump.
# Usage: gate/ctxgate.sh open [lenient|metadata]
# With "open", both sides run the whole SWMMEngine::open (upstream through its C API,
# swmm_engine_open, in gate/diff/ctxopen; the translation in examples/ctxopen), then dump
# the context and the engine members open() sets. open_rc and the exit code are compared too.
# "lenient" sets swmm_engine_set_lenient_open first; "metadata" opens through
# swmm_engine_open_metadata.
H=$(cd "$(dirname "$0")/.." && pwd); W=$(mktemp -d)
if [ "${1:-}" = open ]; then
  T=open; R=${2:-}; C_BIN="$H/gate/diff/ctxopen"; R_BIN="$H/target/release/examples/ctxopen"
else
  T=${1:-ALL}; R=${2:-}; C_BIN="$H/gate/diff/ctxread"; R_BIN="$H/target/release/examples/ctxread"
fi
same=0; diff=0; : > "$W/differ.txt"
{ find /root/work/corpus -name '*.inp'; find /root/work/swmm6src/upstream/tests -name '*.inp'; find "$H/gate/decks" -name '*.inp'; } | sort > "$W/decks.txt"
while IFS= read -r d; do
  if [ "$T" = open ]; then A=("$d" ${R:-open}); else A=("$d" "$T" $R); fi
  ( cd "$(dirname "$d")" && timeout 120 "$C_BIN" "${A[@]}" > "$W/c.txt" 2>/dev/null; echo "exit $?" >> "$W/c.txt"; timeout 120 "$R_BIN" "${A[@]}" > "$W/r.txt" 2>/dev/null; echo "exit $?" >> "$W/r.txt" )
  if cmp -s "$W/c.txt" "$W/r.txt"; then same=$((same+1)); else diff=$((diff+1)); echo "$d" >> "$W/differ.txt"; fi
done < "$W/decks.txt"
echo "ctxgate $T $R: $same identical, $diff differ"; head -${CTXGATE_SHOW:-5} "$W/differ.txt"
[ -n "$CTXGATE_KEEP" ] && cp "$W/differ.txt" "$CTXGATE_KEEP"; rm -rf "$W"
