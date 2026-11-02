#!/bin/bash
# Parse every corpus and upstream test deck with upstream (reference library) and with the
# translation, dump the whole context, compare. Usage: gate/ctxgate.sh [TAGS|ALL] [resolve]
# With "resolve", both sides run resolve_cross_references after reading, then dump.
T=${1:-ALL}; R=${2:-}; H=$(cd "$(dirname "$0")/.." && pwd); W=$(mktemp -d)
same=0; diff=0; : > "$W/differ.txt"
{ find /root/work/corpus -name '*.inp'; find /root/work/swmm6src/upstream/tests -name '*.inp'; find "$H/gate/decks" -name '*.inp'; } | sort > "$W/decks.txt"
while IFS= read -r d; do
  ( cd "$(dirname "$d")" && timeout 120 "$H/gate/diff/ctxread" "$d" "$T" $R > "$W/c.txt" 2>/dev/null; echo "exit $?" >> "$W/c.txt"; timeout 120 "$H/target/release/examples/ctxread" "$d" "$T" $R > "$W/r.txt" 2>/dev/null; echo "exit $?" >> "$W/r.txt" )
  if cmp -s "$W/c.txt" "$W/r.txt"; then same=$((same+1)); else diff=$((diff+1)); echo "$d" >> "$W/differ.txt"; fi
done < "$W/decks.txt"
echo "ctxgate $T $R: $same identical, $diff differ"; head -${CTXGATE_SHOW:-5} "$W/differ.txt"
[ -n "$CTXGATE_KEEP" ] && cp "$W/differ.txt" "$CTXGATE_KEEP"; rm -rf "$W"
