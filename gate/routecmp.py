#!/usr/bin/env python3
"""Compare the routing-step cross test's two streams (gate/ctxgate.sh route).

  routecmp.py upstream.txt translation.txt [loaded.txt]

upstream.txt is gate/diff/routestep's output (BEFORE c / AFTER c blocks), translation.txt is
examples/routestep's (the same BEFORE blocks copied, its own AFTER blocks), loaded.txt the
loader's round trip (examples/routestep ... load). Prints one line:
  identical <calls>                     every AFTER block equal byte for byte
  nottranslated <call> <what>           the translation reached a path it does not translate
  loaddiff <call> <label>               the loader's round trip differs (before stepping)
  differ <call> <label> <index>         the first AFTER block that differs, its first differing
                                        line's label (and the element index for vector lines)
  nostart <rest>                        upstream did not start the run
  nosteps                               upstream ran no routing step in the window
"""
import sys


def blocks(text):
    """{call: (before_lines, after_lines)} and the order of calls."""
    out, order, cur, side = {}, [], None, None
    for ln in text.split(b'\n'):
        if ln.startswith(b'BEFORE '):
            cur = int(ln.split()[1]); order.append(cur); out[cur] = ([], []); side = 0
            continue
        if ln.startswith(b'AFTER '):
            side = 1
            continue
        if ln.startswith(b'END ') or ln.startswith(b'NOSTART'):
            cur = None
            continue
        if cur is not None:
            out[cur][side].append(ln)
    return out, order


def first_diff(a, b):
    """(label, element index) of the first differing line; a vector line's label is the
    previous header line."""
    last = b'?'
    for i in range(max(len(a), len(b))):
        x = a[i] if i < len(a) else b'<end>'
        y = b[i] if i < len(b) else b'<end>'
        if not x.startswith(b' '):
            last = x.split(b' ')[0]
        if x != y:
            if x.startswith(b' ') and y.startswith(b' '):
                xs, ys = x.split(), y.split()
                k = next((k for k in range(max(len(xs), len(ys))) if k >= len(xs) or k >= len(ys) or xs[k] != ys[k]), -1)
                return last.decode(errors='replace'), k
            return (x.split(b' ')[0] if not x.startswith(b' ') else last).decode(errors='replace'), -1
    return None, -1


def main():
    up = open(sys.argv[1], 'rb').read()
    tr = open(sys.argv[2], 'rb').read()
    if up.startswith(b'NOSTART'):
        print('nostart', up.split(b'\n')[0][8:].decode(errors='replace'))
        return
    ub, order = blocks(up)
    if not order:
        print('nosteps')
        return
    if len(sys.argv) > 3:
        lb, _ = blocks(open(sys.argv[3], 'rb').read())
        for c in order:
            if c not in lb or lb[c][0] != ub[c][0]:
                lab, k = first_diff(ub[c][0], lb.get(c, ([], []))[0])
                print('loaddiff', c, lab, k)
                return
    tb, _ = blocks(tr)
    for c in order:
        if c not in tb:
            print('differ', c, '<missing>', -1)
            return
        ta = tb[c][1]
        if len(ta) >= 1 and ta[0].startswith(b'nottranslated '):
            print('nottranslated', c, ta[0][14:].decode(errors='replace'))
            return
        if ta != ub[c][1]:
            lab, k = first_diff(ub[c][1], ta)
            print('differ', c, lab, k)
            return
    print('identical', len(order))


if __name__ == '__main__':
    main()
