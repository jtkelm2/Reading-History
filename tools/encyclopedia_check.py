"""Check an encyclopedia export (tools/encyclopedia_export.py) against the runtime lookup.

At random generic vertices of every stratum (random parameter, random symmetry image), count
  table         distinct histories realizable with reverse folds at the vertex, from Table.valid over
                all rotations (src/lookup.py with the exact cells), as the solver sees them;
  encyclopedia  (stratum, symmetry, card) triples with the vertex on the stratum's domain and the
                card realizable there.
Both sides are compared as sets of histories (hidden directions and steps, in the vertex's own
frame), modulo the symmetries of the vertex: every history must appear exactly once.
Usage: python3 tools/encyclopedia_check.py EXPORT.json [samples per stratum] [seed]"""
import sys, json, math, random
from fractions import Fraction as F
from common import ROOT
sys.path.insert(0, str(ROOT / 'src')); sys.path.insert(0, str(ROOT / 'tools'))
import lookup, cells, rigid
from encyclopedia_export import G, act, inverse

D = json.load(open(sys.argv[1]))
n_samp = int(sys.argv[2]) if len(sys.argv) > 2 else 5
rnd = random.Random(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
T = lookup.Table(str(ROOT / 'tables'), 'reverse', str(ROOT / 'cells'))
St = D['strata']
for s in St:
    if s['kind'] == 'line':
        s['P0e'] = tuple(F(x) for x in s['P0_exact']); s['P1e'] = tuple(F(x) for x in s['P1_exact'])


TP = 2 * math.pi
def dirs_of(Q): return [float(x) * math.pi for x in (F(0), Q[0], Q[0] + Q[1], 1 + Q[1])]
def iso(g, A):
    """The isometry T_g on directions, in the frame with ray directions A."""
    r, m = g
    return (lambda x: (x - A[r]) % TP) if m == 0 else (lambda x: (A[r] - x) % TP)
def canon(Q, steps_dirs):
    """A history at Q (list of steps, each a list of directions in Q's frame), modulo Q's symmetries."""
    A = dirs_of(Q); keys = []
    for g in [g for g in G if act(g, Q) == Q]:
        f = iso(g, A)
        keys.append(tuple(tuple(sorted(round(f(x), 5) % round(TP, 5) for x in S)) for S in steps_dirs))
    return min(keys)


def table_set(Q):
    A = dirs_of(Q); hs = set()
    for r in range(4):
        a = [(A[(i + r) % 4] - A[r]) % TP for i in range(4)]
        for t, h in T.valid(a, None):
            d = [a[x[1]] if x[0] == 'v' else h[x[1]] for x in t['rays']]
            hs.add(canon(Q, [[(d[j] + A[r]) % TP for j in S] for S in t['steps']]))
    return hs


def enc_set(Q):
    A = dirs_of(Q); hs = set(); n = 0
    for s in St:
        seen = set()
        for g in G:
            R = act(g, Q)
            if R in seen: continue
            seen.add(R)
            if s['kind'] == 'point':
                if R != tuple(F(x) for x in s['point_exact']): continue
                tt = 0.5
            else:
                (x0, y0), (x1, y1) = s['P0e'], s['P1e']
                if (x1 - x0) * (R[1] - y0) != (y1 - y0) * (R[0] - x0): continue
                tq = (R[0] - x0) / (x1 - x0) if x1 != x0 else (R[1] - y0) / (y1 - y0)
                if not (0 < tq <= F(s["tmax"])) or tq == 1: continue
                tt = float(tq)
            sp = next((v for x, xe, v in s.get('special', []) if F(xe) == tq), None) if s['kind'] == 'line' else None
            for i, c in enumerate(s['cards']):
                ok = sp[i] if sp is not None else any(lo < tt < hi and v for lo, hi, v in c['iv'])
                if not ok: continue
                n += 1
                h = tuple(c['g']); q = act(inverse(h), R); a = dirs_of(q)
                lams = next(l for uu, l in lookup.branches(c, a) if list(uu) == c['u']) if c['scaffold'] else None
                d = cells.replay(rigid.NUM, c, a, lams)
                assert d is not None, 'a realizable card does not replay'
                to_rep = iso(h, a); rr, mm = g; AR = A[rr]
                to_q = (lambda y: (y + AR) % TP) if mm == 0 else (lambda y: (AR - y) % TP)
                hs.add(canon(Q, [[to_q(to_rep(d[j])) for j in S] for S in c['steps']]))
    return hs, n


bad = tot = 0
for s in St:
    if s['kind'] != 'line': continue
    for _ in range(n_samp):
        tq = F(rnd.randrange(1, 10 ** 6), 10 ** 6) * F(s['tmax'])
        P = (s['P0e'][0] + tq * (s['P1e'][0] - s['P0e'][0]), s['P0e'][1] + tq * (s['P1e'][1] - s['P0e'][1]))
        Q = act(rnd.choice(G), P)
        a, (b, nb) = table_set(Q), enc_set(Q); tot += 1
        if a != b or nb != len(b):
            bad += 1
            if bad <= 5: print('disagree at', [float(x) * 180 for x in Q], 'table', len(a), 'encyclopedia', len(b), 'cards', nb)
special = {tuple(F(x) for x in s['point_exact']) for s in St if s['kind'] == 'point'} | {(F(p[2]), F(p[3])) for p in D['pins']}
for s in St:                                                  # special points of the lines
    if s['kind'] != 'line': continue
    for x, xe, v in s['special']:
        tq = F(xe); P = (s['P0e'][0] + tq * (s['P1e'][0] - s['P0e'][0]), s['P0e'][1] + tq * (s['P1e'][1] - s['P0e'][1]))
        special.add(min((act(g, P) for g in G), key=lambda q: (q[0] > q[1] or q[1] > F(1, 2), q)))
sb = dup = 0
for Q in sorted(special):
    a, (b, nb) = table_set(Q), enc_set(Q)
    if a != b:
        sb += 1
        if sb <= 8: print('special vertex', [float(x) * 180 for x in Q], 'table', len(a), 'encyclopedia', len(b), 'only table', len(a - b), 'only encyclopedia', len(b - a))
    dup += nb - len(b)
print('special vertices checked', len(special), 'disagree', sb, '; cards equivalent to another card at their vertex', dup)
print('generic vertices checked', tot, 'disagree', bad)
