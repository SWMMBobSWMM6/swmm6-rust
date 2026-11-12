#!/usr/bin/env python3
"""Run decks through two OpenSWMM 6 command-line programs and compare what each leaves:
the exit status, the .out bytes and the .rpt with the wall-clock lines masked (and, with
--stderr, the bytes each program writes to stderr).

Usage: gate/twobuilds.py A_BIN B_BIN LIST OUTDIR [-j N] [--timeout S] [--stderr] [--summary]
  LIST is a file with one deck path per line. Each run goes through gate/openswmm_saves.sh
  (the deck's SAVE files are backed up and restored, OpenMP as the reference ran), in the
  deck's own folder, writing its .rpt/.out under OUTDIR. Prints one line per deck that is not
  the same (differ: what | stopped: why) and a summary; OUTDIR/results.tsv keeps the per-deck
  record (deck, status, exit A, exit B, .out sizes, tag, 2D flag) and makes a rerun resume.
  A deck counts as "stopped" when a program exits with status 3 and says on stderr that it
  reached code not translated yet (the translation's command-line program,
  src/bin/openswmm.rs); its files are then not compared. --summary prints the counts per deck
  category (the corpus, upstream's tests, examples and docs, each folder of gate/decks) and,
  within each, how many decks carry a [2D_*] section and how many ran to completion.
  A is the reference, B the translation. On stderr the PID in bash's notice of a program
  killed by a signal (the U1/U57 aborts) is masked.

Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson;
--stderr, the stopped status and --summary added 12 November 2026 (phase 2a wave 5).
"""
import argparse, collections, concurrent.futures as cf, hashlib, os, re, subprocess, sys

H = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLOCK = re.compile(rb'(Analysis begun on|Analysis ended on|Total elapsed time)')
STUB = b'not translated yet: '
TWOD = re.compile(rb'^\s*\[2D_', re.M | re.I)
WRAPPER_PID = re.compile(rb'^(.*openswmm_saves\.sh: line \d+:\s+)\d+ ', re.M)


def run(bin_, deck, rpt, out, err, timeout):
    env = dict(os.environ, OPENSWMM_BIN=bin_)
    with open(err, 'wb') as ef:
        # its own process group, so a timeout ends the program and not only the wrapper
        p = subprocess.Popen(['nice', os.path.join(H, 'gate', 'openswmm_saves.sh'), os.path.basename(deck), rpt, out],
                             cwd=os.path.dirname(os.path.abspath(deck)), env=env,
                             stdout=subprocess.DEVNULL, stderr=ef, start_new_session=True)
        try:
            return p.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(p.pid, 9)
            p.wait()
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


def errmasked(b):
    """stderr without the PID bash prints when the program dies on a signal (the wrapper's
    "openswmm_saves.sh: line N: PID Aborted ..." notice)."""
    if b is None:
        return None
    return WRAPPER_PID.sub(rb'\1PID ', b)


def stub_of(rc, err):
    if rc == 3 and err and STUB in err:
        return err.split(STUB, 1)[1].split(b'\n', 1)[0].decode('latin-1')
    return None


def one(i, deck, a, b, outdir, timeout, with_stderr):
    tag = '%04d_%s' % (i, hashlib.sha1(deck.encode()).hexdigest()[:8])
    res = {}
    for side, bin_ in (('a', a), ('b', b)):
        rpt = os.path.join(outdir, tag + '.' + side + '.rpt')
        out = os.path.join(outdir, tag + '.' + side + '.out')
        err = os.path.join(outdir, tag + '.' + side + '.err')
        for f in (rpt, out, err):
            if os.path.exists(f):
                os.remove(f)
        rc = run(bin_, deck, rpt, out, err, timeout)
        res[side] = (rc, read(rpt), read(out), read(err))
    (ra, rpa, oa, ea), (rb, rpb, ob, eb) = res['a'], res['b']
    diffs = []
    stopped = stub_of(rb, eb) or stub_of(ra, ea)
    if stopped is None:
        if ra != rb:
            diffs.append('exit %s/%s' % (ra, rb))
        if oa != ob:
            diffs.append('out')
        if masked(rpa) != masked(rpb):
            diffs.append('rpt')
        if with_stderr and errmasked(ea) != errmasked(eb):
            diffs.append('stderr')
    st = 'stopped: ' + stopped if stopped is not None else ('same' if not diffs else 'differ: ' + ', '.join(diffs))
    if st == 'same' or stopped is not None:
        for side in 'ab':
            for ext in ('rpt', 'out', 'err'):
                p = os.path.join(outdir, tag + '.' + side + '.' + ext)
                if os.path.exists(p):
                    os.remove(p)
    try:
        twod = 1 if TWOD.search(open(deck, 'rb').read()) else 0
    except OSError:
        twod = 0
    return deck, st, ra, rb, len(oa) if oa is not None else -1, len(ob) if ob is not None else -1, tag, twod


def category(deck):
    p = deck.replace('\\', '/')
    m = re.search(r'/gate/decks/([^/]+)/', p)
    if m:
        return 'gate/decks/' + m.group(1)
    if '/gate/decks/' in p:
        return 'gate/decks'
    for k in ('tests', 'examples', 'docs'):
        if re.search(r'/(upstream|up)/%s/' % k, p):
            return 'upstream ' + k
    if '/corpus/' in p:
        return 'corpus'
    return 'other'


def summary(rows):
    cats = collections.OrderedDict()
    for deck, st, ra, rb, twod in sorted(rows, key=lambda r: category(r[0])):
        c = cats.setdefault(category(deck), collections.Counter())
        k = st.split(':')[0]
        c[k] += 1
        c['decks'] += 1
        c['2d'] += twod
        if k == 'same' and str(ra) == '0':
            c['complete'] += 1
            c['complete 2d'] += twod
        if k == 'same' and str(ra) != '0':
            c['same, exit %s' % ra] += 1
    tot = collections.Counter()
    for name, c in cats.items():
        tot.update(c)
        extra = ', '.join('%s %d' % (k, v) for k, v in sorted(c.items()) if k.startswith('same, exit'))
        print('  %-24s %5d decks (%d with [2D_*]): same %d (complete runs %d, %d of them 2D%s), differ %d, stopped %d'
              % (name, c['decks'], c['2d'], c['same'], c['complete'], c['complete 2d'],
                 '; ' + extra if extra else '', c['differ'], c['stopped']))
    print('  %-24s %5d decks (%d with [2D_*]): same %d (complete runs %d, %d of them 2D), differ %d, stopped %d'
          % ('all', tot['decks'], tot['2d'], tot['same'], tot['complete'], tot['complete 2d'], tot['differ'], tot['stopped']))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('a'); ap.add_argument('b'); ap.add_argument('list'); ap.add_argument('outdir')
    ap.add_argument('-j', type=int, default=2); ap.add_argument('--timeout', type=int, default=1800)
    ap.add_argument('--stderr', action='store_true'); ap.add_argument('--summary', action='store_true')
    o = ap.parse_args()
    o.outdir = os.path.abspath(o.outdir)
    o.a, o.b = os.path.abspath(o.a), os.path.abspath(o.b)  # the runs are in the decks' folders
    os.makedirs(o.outdir, exist_ok=True)
    decks = [l.strip() for l in open(o.list) if l.strip()]
    res_path = os.path.join(o.outdir, 'results.tsv')
    done = {}
    if os.path.exists(res_path):  # resume: keep the decks already compared
        for l in open(res_path):
            f = l.rstrip('\n').split('\t')
            if len(f) >= 2:
                done[f[0]] = f
    todo = [(i, d) for i, d in enumerate(decks) if d not in done]
    with open(res_path, 'a') as log, cf.ThreadPoolExecutor(o.j) as ex:
        futs = [ex.submit(one, i, d, o.a, o.b, o.outdir, o.timeout, o.stderr) for i, d in todo]
        for f in cf.as_completed(futs):
            r = f.result()
            deck, st, tag = r[0], r[1], r[6]
            fields = list(map(str, r))
            done[deck] = fields
            log.write('\t'.join(fields) + '\n'); log.flush()
            if st != 'same':
                print(st, deck, tag, flush=True)
    rows = []
    for d in decks:
        f = done[d]
        rows.append((d, f[1], f[2], f[3], int(f[7]) if len(f) > 7 else 0))
    n = collections.Counter(r[1].split(':')[0] for r in rows)
    print('twobuilds: %d same, %d differ, %d stopped, of %d' % (n['same'], n['differ'], n['stopped'], len(decks)), flush=True)
    if o.summary:
        summary(rows)


if __name__ == '__main__':
    main()
