#!/bin/bash
# The translation's command-line program (target/release/openswmm) for
# /home/claude/swmm6-ports/gate/gate.py, which runs `ENGINE model.inp report.rpt output.out`
# in the model's folder: the files the deck SAVEs to ([FILES], gate/savefiles.py) are backed
# up first and restored after, so a SAVE in the corpus never changes a file another deck USEs.
# The exit status is the program's (bash reports a child killed by signal N as 128 + N, as
# the reference's was recorded). OPENSWMM_BIN overrides the program.
H=$(cd "$(dirname "$0")/.." && pwd)
BIN=${OPENSWMM_BIN:-$H/target/release/openswmm}
d="$PWD/$1"
B=$(mktemp -d); python3 "$H/gate/savefiles.py" "$d" > "$B/saves.txt" 2>/dev/null; k=0
while IFS= read -r f; do k=$((k+1)); [ -f "$f" ] && cp -p "$f" "$B/$k"; done < "$B/saves.txt"
"$BIN" "$@"; rc=$?
k=0; while IFS= read -r f; do k=$((k+1)); if [ -f "$B/$k" ]; then cp -p "$B/$k" "$f"; else rm -f "$f"; fi; done < "$B/saves.txt"
rm -rf "$B"; exit $rc
