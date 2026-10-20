"""Check exact cell trees (src/cells.py) against the numeric computation.
1. Witnesses: at the stored vertex of every leaf, the numeric computation gives the leaf's value,
   and evaluating the tree there reaches a leaf with that value.
2. Corpus: at every rotation of every corpus vertex shape of the tree's degree, for every type,
   the tree's value equals the numeric computation.
Usage: python3 tools/verify_cells.py CELLS.json TABLE.json [--witnesses=PATH] [pattern.fold ...]"""
import sys, json, math
from common import fh, resolve_tol
sys.path.insert(0, str(fh.__file__).rsplit('/', 1)[0])
import cells, rigid

WP = next((x[12:] for x in sys.argv if x.startswith('--witnesses=')), None)
sys.argv = [x for x in sys.argv if not x.startswith('--witnesses=')]
C = json.load(open(sys.argv[1])); T = json.load(open(sys.argv[2]))
if WP: C['witnesses'] = json.load(open(WP))['witnesses']
types = T['types'][:len(C['trees'])]; n = T['meta']['n']
den = C['meta'].get('den', 4)
bad = []; nw = 0
for i, (t, tree, ws) in enumerate(zip(types, C['trees'], C['witnesses'])):
    for w in ws:
        b, val = w[0], w[1]; u = tuple(w[2]) if len(w) > 2 else None
        a = [0.0] + [math.pi * x for x in b]; nw += 1
        num = cells.decide(rigid.NUM, t, a, None, u); ev = cells.evaluate(tree, a, den=den, u=u)
        if num != val or ev != val: bad.append(('witness', i, b, u, val, num, ev))
print('witnesses: %d checked, %d disagree' % (nw, len(bad)))
seen, nc, badc = set(), 0, 0
for path in sys.argv[3:]:
    P = resolve_tol(path)
    for p in P.interior:
        if len(P.H[p]) != n: continue
        D = sorted(d for d, _ in P.H[p]); secs = tuple(round(P.deg(D[(j + 1) % n] - D[j]), 4) for j in range(n))
        if secs in seen: continue
        seen.add(secs)
        for r in range(n):
            a = [fh.nrm(D[(j + r) % n] - D[r]) for j in range(n)]
            for i, (t, tree) in enumerate(zip(types, C['trees'])):
                for u in (cells.branch_keys(t) if 'scaffold' in t else [None]):
                    num = cells.decide(rigid.NUM, t, a, None, u); ev = cells.evaluate(tree, a, den=den, u=u); nc += 1
                    if num != ev:
                        badc += 1
                        if badc <= 5: print('corpus disagreement', secs, r, i, u, num, ev)
print('corpus: %d (vertex rotation, type[, branch]) pairs checked over %d shapes, %d disagree' % (nc, len(seen), badc))
for x in bad[:5]: print(x)
assert not bad and not badc
print('all agree')
