#!/usr/bin/env python3
"""Run decks through two OpenSWMM 6 command-line programs and compare what each leaves:
the exit status, the .out bytes and the .rpt with the wall-clock lines masked.

Usage: gate/twobuilds.py A_BIN B_BIN LIST OUTDIR [-j N] [--timeout S]
  LIST is a file with one deck path per line. Each run goes through gate/openswmm_saves.sh
  (the deck's SAVE files are backed up and restored, OpenMP as the reference ran), in the
  deck's own folder, writing its .rpt/.out under OUTDIR. Prints one line per deck
  (same | differ: what) and a summary; OUTDIR/results.tsv keeps the per-deck record
  (deck, status, exit A, exit B, .out sizes).

Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson.
"""
import argparse, concurrent.futures as cf, hashlib, os, re, subprocess, sys

H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLOCK = re.compile(rb'(Analysis begun on|Analysis ended on|Total elapsed time)')


def run(bin_, deck, rpt, out, timeout):
    env = dict(os.environ, OPENSWMM_BIN=bin_)
    try:
        p = subprocess.run([os.path.join(H, 'gate', 'openswmm_saves.sh'), os.path.basename(deck), rpt, out],
                           cwd=os.path.dirname(os.path.abspath(deck)), env=env, timeout=timeout,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return p.returncode
    except subprocess.TimeoutExpired:
        return 'timeout'


def read(p):
    try:
        return open(p, 'rb').read()
    except OSError:
        return None


def masked(b):
    if b is None:
        return None
    return b'\n'.join(l for l in b.split(b'\n') if not CLOCK.search(l))


def one(i, deck, a, b, outdir, timeout):
    tag = '%04d_%s' % (i, hashlib.sha1(deck.encode()).hexdigest()[:8])
    res = {}
    for side, bin_ in (('a', a), ('b', b)):
        rpt = os.path.join(outdir, tag + '.' + side + '.rpt')
        out = os.path.join(outdir, tag + '.' + side + '.out')
        for f in (rpt, out):
            if os.path.exists(f):
                os.remove(f)
        res[side] = (run(bin_, deck, rpt, out, timeout), read(rpt), read(out))
    (ra, rpa, oa), (rb, rpb, ob) = res['a'], res['b']
    diffs = []
    if ra != rb:
        diffs.append('exit %s/%s' % (ra, rb))
    if oa != ob:
        diffs.append('out')
    if masked(rpa) != masked(rpb):
        diffs.append('rpt')
    for side in 'ab':
        for ext in ('rpt', 'out'):
            p = os.path.join(outdir, tag + '.' + side + '.' + ext)
            if not diffs and os.path.exists(p):
                os.remove(p)
    return deck, ('same' if not diffs else 'differ: ' + ', '.join(diffs)), ra, rb, \
        len(oa) if oa is not None else -1, len(ob) if ob is not None else -1, tag


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('a'); ap.add_argument('b'); ap.add_argument('list'); ap.add_argument('outdir')
    ap.add_argument('-j', type=int, default=2); ap.add_argument('--timeout', type=int, default=1800)
    o = ap.parse_args()
    os.makedirs(o.outdir, exist_ok=True)
    decks = [l.strip() for l in open(o.list) if l.strip()]
    res_path = os.path.join(o.outdir, 'results.tsv')
    done = {}
    if os.path.exists(res_path):  # resume: keep the decks already compared
        for l in open(res_path):
            f = l.rstrip('\n').split('\t')
            if len(f) >= 2:
                done[f[0]] = f[1]
    n_same = sum(1 for d in decks if done.get(d) == 'same')
    todo = [(i, d) for i, d in enumerate(decks) if d not in done]
    with open(res_path, 'a') as log, cf.ThreadPoolExecutor(o.j) as ex:
        futs = [ex.submit(one, i, d, o.a, o.b, o.outdir, o.timeout) for i, d in todo]
        for f in cf.as_completed(futs):
            deck, st, ra, rb, la, lb, tag = f.result()
            n_same += st == 'same'
            log.write('\t'.join(map(str, (deck, st, ra, rb, la, lb, tag))) + '\n'); log.flush()
            if st != 'same':
                print(st, deck, tag, flush=True)
    n_done = len(decks)  # every deck now has a line in results.tsv
    print('twobuilds: %d same, %d differ, of %d' % (n_same, n_done - n_same, len(decks)), flush=True)


if __name__ == '__main__':
    main()
