#!/bin/bash
# Run a command on a corpus deck in the deck's own folder with the files the deck SAVEs to
# ([FILES], gate/savefiles.py) backed up first and restored after, as ctxgate.sh does.
#   gate/withsaves.sh DECK CMD [ARGS...]    (DECK is an absolute path; CMD runs in its folder)
H=$(cd "$(dirname "$0")/.." && pwd); d=$1; shift
B=$(mktemp -d); python3 "$H/gate/savefiles.py" "$d" > "$B/saves.txt" 2>/dev/null; k=0
while IFS= read -r f; do k=$((k+1)); [ -f "$f" ] && cp -p "$f" "$B/$k"; done < "$B/saves.txt"
( cd "$(dirname "$d")" && "$@" ); rc=$?
k=0; while IFS= read -r f; do k=$((k+1)); if [ -f "$B/$k" ]; then cp -p "$B/$k" "$f"; else rm -f "$f"; fi; done < "$B/saves.txt"
rm -rf "$B"; exit $rc
