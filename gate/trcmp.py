"""Compare two step traces (steptrace format) bit for bit and report the first departure:
which step, at what elapsed time, which quantity and which object, with both values."""
import struct, sys

FIELDS = ('subcatch runoff', 'node depth', 'node volume', 'link flow', 'link depth')

def frames(path):
    with open(path, 'rb') as f:
        assert f.read(4) == b'S6TR', path + ': not a step trace'
        ver, ns, nn, nl = struct.unpack('<4i', f.read(16))
        sizes = (ns, nn, nn, nl, nl)
        width = 8 * (1 + sum(sizes))
        while True:
            rec = f.read(width)
            if len(rec) < width: return (ns, nn, nl), sizes
            yield rec, sizes

def main(a, b):
    ga, gb = frames(a), frames(b)
    step = 0
    while True:
        try: ra, sizes = next(ga)
        except StopIteration as e: ra = None
        try: rb, _ = next(gb)
        except StopIteration: rb = None
        if ra is None or rb is None:
            if ra is None and rb is None:
                print(f'identical: {step} steps'); return 0
            print(f'step counts differ: one trace ends after {step} steps'); return 1
        step += 1
        if ra == rb: continue
        va = struct.unpack(f'<{len(ra)//8}d', ra); vb = struct.unpack(f'<{len(rb)//8}d', rb)
        if va[0] != vb[0]:
            print(f'step {step}: elapsed time differs {va[0]!r} vs {vb[0]!r}'); return 1
        k = 1
        for name, n in zip(FIELDS, sizes):
            for j in range(n):
                if ra[8*(k+j):8*(k+j+1)] != rb[8*(k+j):8*(k+j+1)]:
                    print(f'step {step} (t={va[0]:.9f} d): {name} [{j}] {va[k+j]!r} vs {vb[k+j]!r}'); return 1
            k += n

if __name__ == '__main__':
    sys.exit(main(sys.argv[1], sys.argv[2]))
