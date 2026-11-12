#!/usr/bin/env python3
"""Generate the run-time exchange gate decks in gate/decks/twodx/ (phase 2a wave 4b).

gate/decks/twodc varies one [2D_OPTIONS] key at a time on the coupled example; these decks add
the exchange paths none of them reaches (see gate/twodexch.sh's coverage line):

  outfall_backflow(_gated)   surface water at a 1D outfall's cell stands above the network, so
                             the outfall draws water off the mesh (OutfallExchange withdrawals,
                             the reserved donor budget, the tailwater engaged and drowned); the
                             gated copy shuts the reverse flow
  outfall_backflow_spread    the same with [2D_NODE_COUPLING] OUTFALL SPREAD and TAILWATER_RELAX
  outfall_vertex(_stencil, _vertex, _vfr)
                             the outfall coupled at a mesh vertex (the deepest stencil cell, the
                             sampling cells of COUPLING_HEAD STENCIL_WET_MEAN / VERTEX, VFR)
  multi_point(_weir, _stencil)
                             one junction coupled at two vertices and a triangle (the ledger sums
                             points per node), under ORIFICE, WEIR_ORIFICE with tracking ponds,
                             and STENCIL_WET_MEAN heads with the equalize cap

Run from anywhere: `python3 gate/twodx_gen.py`.
Written for swmm6-rust, 12 November 2026, by Claude (Anthropic) for Robert Dickinson. The decks
derive from OpenSWMM's examples/2d_complete_example.inp (Copyright 2026 Caleb Buahin, Apache
License, Version 2.0) through gate/twodc_gen.py; the header of each deck says so.
"""
import os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import twodc_gen as g  # noqa: E402

OUT = os.path.join(g.H, 'gate', 'decks', 'twodx')
HEADER = g.HEADER.replace('gate/twodc_gen.py', 'gate/twodx_gen.py')


def write(name, what, secs):
    p = os.path.join(OUT, name + '.inp')
    with open(p, 'w', encoding='utf-8', newline='\n') as f:
        f.write(HEADER.format(name=name, what=what))
        f.write(g.join(secs))


def base_backflow(gated):
    """The outfall deck with the outfall's cell flooded from its south edge, a junction below
    the outfall and a second outfall below that, so surface water runs back into OUT3."""
    s = g.base_outfall(gated)
    # JX below OUT3 (reverse fall), drained by CY to a lower free outfall OUT4
    s = g.put(s, 'JUNCTIONS', g.get(s, 'JUNCTIONS').replace(
        'JX      101.5      1.0       0.0        0.0       0.0',
        'JX      99.5       1.0       0.0        0.0       0.0'))
    s = g.put(s, 'OUTFALLS', g.get(s, 'OUTFALLS').rstrip('\n') + '\nOUT4    99.0       FREE   *          NO\n\n')
    s = g.put(s, 'CONDUITS', g.get(s, 'CONDUITS').rstrip('\n') +
              '\nCY      JX        OUT4     30     0.013      0         0          0         0\n\n')
    s = g.put(s, 'XSECTIONS', g.get(s, 'XSECTIONS').rstrip('\n') + '\nCY      CIRCULAR  0.30   0      0      0      1\n\n')
    s = g.put(s, 'COORDINATES', g.get(s, 'COORDINATES').rstrip('\n') + '\nJX       5.0       -3.0\nOUT3     2.0       2.0\nOUT4     5.0       -6.0\n\n')
    # a tidal stage on T0's south edge floods the outfall's cell for two hours
    ts = g.get(s, 'TIMESERIES') + ('\nBACK_TS     01/01/2024   00:00  100.00\nBACK_TS                  01:00  100.00\n'
                                   'BACK_TS                  01:30  100.90\nBACK_TS                  03:30  100.90\n'
                                   'BACK_TS                  04:00  100.00\nBACK_TS                  06:00  100.00\n')
    s = g.put(s, 'TIMESERIES', ts)
    bc = g.get(s, '2D_BOUNDARY_CONDITIONS').replace(
        '   0     2      NORMAL_FLOW       0.0050          *         *           ;; T0 south edge (v0,v1)',
        '   0     2      TS_STAGE          BACK_TS         *         back        ;; T0 south edge (v0,v1): flood')
    s = g.put(s, '2D_BOUNDARY_CONDITIONS', bc)
    return s


def vertex_outfall(s):
    """Move OUT3's coupling from triangle T0 to vertex v1 (south mid: T0, T2, T3)."""
    m = g.get(s, '2D_TRIANGLE_NODE_MAP')
    m = '\n'.join(l for l in m.split('\n') if 'OUT3' not in l) + '\n'
    s = g.put(s, '2D_TRIANGLE_NODE_MAP', m)
    v = g.get(s, '2D_VERTEX_NODE_MAP').rstrip('\n') + '\n   1       OUT3       0.65   0.5       ;; south-mid vertex -> outfall\n\n'
    return g.put(s, '2D_VERTEX_NODE_MAP', v)


def multi_point(s):
    """J1 coupled at v4 (as in the example), at v3 and at triangle T4 too."""
    v = g.get(s, '2D_VERTEX_NODE_MAP').rstrip('\n') + '\n   3       J1         0.60   1.0       ;; west-mid vertex -> J1 again\n\n'
    s = g.put(s, '2D_VERTEX_NODE_MAP', v)
    t = g.get(s, '2D_TRIANGLE_NODE_MAP').rstrip('\n') + '\n   4                 J1          0.55   0.8       ;; NW-lower cell -> J1 again\n\n'
    return g.put(s, '2D_TRIANGLE_NODE_MAP', t)


def opt(s, kv):
    o = g.get(s, '2D_OPTIONS')
    for k, v in kv:
        o = g.set_opt(o, k, v)
    return g.put(s, '2D_OPTIONS', o)


def main():
    os.makedirs(OUT, exist_ok=True)
    for f in os.listdir(OUT):
        if f.endswith('.inp'):
            os.remove(os.path.join(OUT, f))
    n = 0
    write('outfall_backflow', 'surface water runs back into a 1D outfall (withdrawals)', base_backflow(False)); n += 1
    write('outfall_backflow_gated', 'the same with a flap gate on the outfall', base_backflow(True)); n += 1
    s = base_backflow(False)
    s = g.put(s, '2D_NODE_COUPLING', 'OUT3  OUTFALL SPREAD=12 TAILWATER_RELAX=60\n\n', after='2D_TRIANGLE_NODE_MAP')
    write('outfall_backflow_spread', 'backflow with OUTFALL SPREAD=12 TAILWATER_RELAX=60', s); n += 1
    s = vertex_outfall(base_backflow(False))
    write('outfall_vertex', 'the backflow outfall coupled at a mesh vertex', s); n += 1
    write('outfall_vertex_stencil', 'vertex outfall, COUPLING_HEAD STENCIL_WET_MEAN',
          opt(vertex_outfall(base_backflow(False)), [('COUPLING_HEAD', 'STENCIL_WET_MEAN')])); n += 1
    write('outfall_vertex_vertex', 'vertex outfall, COUPLING_HEAD VERTEX',
          opt(vertex_outfall(base_backflow(False)), [('COUPLING_HEAD', 'VERTEX')])); n += 1
    write('outfall_vertex_vfr', 'vertex outfall, CELL_CLOSURE VFR',
          opt(vertex_outfall(base_backflow(False)), [('CELL_CLOSURE', 'VFR')])); n += 1
    write('multi_point', 'J1 coupled at two vertices and a triangle', multi_point(g.base_spill())); n += 1
    write('multi_point_weir', 'J1 at three points, COUPLING_LAW WEIR_ORIFICE, COUPLING_POND_MODE TRACKING',
          opt(multi_point(g.base_spill()), [('COUPLING_LAW', 'WEIR_ORIFICE'), ('COUPLING_POND_MODE', 'TRACKING')])); n += 1
    write('multi_point_stencil', 'J1 at three points, COUPLING_HEAD STENCIL_WET_MEAN, COUPLING_EQUALIZE YES',
          opt(multi_point(g.base_spill()), [('COUPLING_HEAD', 'STENCIL_WET_MEAN'), ('COUPLING_EQUALIZE', 'YES')])); n += 1
    print('twodx_gen: %d decks in %s' % (n, OUT))


if __name__ == '__main__':
    main()
