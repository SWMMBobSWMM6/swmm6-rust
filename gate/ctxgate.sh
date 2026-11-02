#!/bin/bash
# Parse every corpus and upstream test deck with upstream (reference library) and with the
# translation, dump the whole context, compare. Usage: gate/ctxgate.sh [TAGS|ALL]
T=${1:-ALL}; H=$(cd "$(dirname "$0")/.." && pwd); W=$(mktemp -d)
same=0; diff=0; : > "$W/differ.txt"
{ find /root/work/corpus -name '*.inp'; find /root/work/swmm6src/upstream/tests -name '*.inp'; } | sort > "$W/decks.txt"
while IFS= read -r d; do
  ( cd "$(dirname "$d")" && timeout 120 "$H/gate/diff/ctxread" "$d" "$T" > "$W/c.txt" 2>/dev/null; timeout 120 "$H/target/release/examples/ctxread" "$d" "$T" > "$W/r.txt" 2>/dev/null )
  if cmp -s "$W/c.txt" "$W/r.txt"; then same=$((same+1)); else diff=$((diff+1)); echo "$d" >> "$W/differ.txt"; fi
done < "$W/decks.txt"
echo "ctxgate $T: $same identical, $diff differ"; head -5 "$W/differ.txt"; rm -rf "$W"
