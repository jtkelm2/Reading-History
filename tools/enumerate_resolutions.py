"""List every minimal-cost geometry of a pattern within the search space. A geometry is the
choice of hidden directions at the vertices together with the paths of the cancelled segments.
Usage: python3 tools/enumerate_resolutions.py pattern.fold [--max=N]"""
import sys
from common import ROOT, fh, resolve_tol, geometries

args = [a for a in sys.argv[1:] if not a.startswith('--')]
cap = int(next((a[6:] for a in sys.argv[1:] if a.startswith('--max=')), 1000))
P = resolve_tol(args[0]); res = fh.solve(P, log=lambda s: None)
assert res['status'] == 'ok', res['status']
G, complete = geometries(P, res, cap)
sels = {}
for g in G: sels.setdefault(g['selection'], []).append(g)
print('cost %d: %d minimal geometries%s, %d choices of hidden directions' % (
    res['cost'], len(G), '' if complete else ' (INCOMPLETE: cap or time budget reached)', len(sels)))
common = set.intersection(*[set(s) for s in sels]) if sels else set()
print('the choices differ at vertices', sorted({p for s in sels for p, _ in s} - {p for p, _ in common}))
for i, g in enumerate(sorted(G, key=lambda g: g['length'])):
    print('%3d  length %.4f  %s' % (i, g['length'], {p: R for p, R in g['selection'] if (p, R) not in common}))
