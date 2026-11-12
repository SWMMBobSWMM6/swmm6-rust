#!/usr/bin/env python3
# Compare a 2D kernel differential's two outputs (gate/diff/twod_*.cpp against
# examples/twod_*.rs) line by line, after `cmp` has found them different. Each line is a tag and
# hex-bit tokens, the inputs, then "|" and the results. Reports, per tag, the lines that differ:
#   input    : the inputs differ (the two generators disagree; a harness bug)
#   nan-sign : every differing result is a NaN on both sides, differing only in the sign bit, and
#              some input is non-finite (the x86/LLVM operand-order hazard of PLAN.md: two NaNs
#              of opposite sign meet in a commutative operation)
#   nan-sign-finite : the same, with all inputs finite (to be looked at)
#   nan      : every differing result is a NaN on both sides (payload differs)
#   value    : anything else (a real difference)
# with the first few examples of each. Exit status 0 when the only differences are nan-sign.
# Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
import sys
from collections import defaultdict


def nonfinite(tok):
    try:
        b = int(tok, 16)
    except ValueError:
        return False
    return len(tok) == 16 and (b >> 52) & 0x7FF == 0x7FF


def isnan(tok):
    try:
        b = int(tok, 16)
    except ValueError:
        return False
    return len(tok) == 16 and (b >> 52) & 0x7FF == 0x7FF and (b & ((1 << 52) - 1)) != 0


def main():
    a_path, b_path = sys.argv[1], sys.argv[2]
    shown = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    counts = defaultdict(lambda: defaultdict(int))
    examples = defaultdict(list)
    total = 0
    with open(a_path) as fa, open(b_path) as fb:
        for n, (la, lb) in enumerate(zip(fa, fb), 1):
            total += 1
            if la == lb:
                continue
            ta, tb = la.split(), lb.split()
            tag = ta[0] if ta else "?"
            kind = "value"
            ia = ta.index("|") if "|" in ta else len(ta)
            ib = tb.index("|") if "|" in tb else len(tb)
            if ta[:ia] != tb[:ib]:
                kind = "input"
            elif len(ta) == len(tb):
                diffs = [(x, y) for x, y in zip(ta, tb) if x != y]
                if all(isnan(x) and isnan(y) for x, y in diffs):
                    if all(int(x, 16) ^ int(y, 16) == 1 << 63 for x, y in diffs):
                        kind = "nan-sign" if any(nonfinite(x) for x in ta[1:ia]) else "nan-sign-finite"
                    else:
                        kind = "nan"
            counts[tag][kind] += 1
            if len([e for e in examples[(tag, kind)]]) < shown:
                examples[(tag, kind)].append((n, la.rstrip(), lb.rstrip()))
        extra_a = sum(1 for _ in fa)
        extra_b = sum(1 for _ in fb)
    print(f"{total} lines compared")
    if extra_a or extra_b:
        print(f"length differs: {extra_a} more lines in {a_path}, {extra_b} more in {b_path}")
    for tag in sorted(counts):
        for kind, c in sorted(counts[tag].items()):
            print(f"  {tag:10s} {kind:9s} {c} lines")
            for n, la, lb in examples[(tag, kind)]:
                print(f"    line {n}\n      C++  {la}\n      Rust {lb}")
    # pass: identical, or differing only in NaN signs on lines with a non-finite input
    bad = sum(c for tag in counts for kind, c in counts[tag].items() if kind != "nan-sign")
    nsign = sum(c for tag in counts for kind, c in counts[tag].items() if kind == "nan-sign")
    if bad or extra_a or extra_b:
        print("DIFFER")
        sys.exit(1)
    print(f"identical except {nsign} lines whose results differ only in a NaN's sign (non-finite inputs)")
    sys.exit(0)


if __name__ == "__main__":
    main()
