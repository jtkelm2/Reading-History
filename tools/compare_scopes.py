"""Folds and forced order of one pattern under each notion of 'fixed':
  directions  the hidden directions at the vertices are fixed, the segments' paths are free;
  found       the geometry the solver found (directions and paths);
  shortest    the minimal geometry whose cancelled segments are shortest (needs all geometries);
  minimal     every minimal resolution (visible creases only).
Usage: python3 tools/compare_scopes.py pattern.fold [--budget=SECONDS per analysis]"""
import sys, time
import networkx as nx
from common import fh, resolve_tol, geometries

args = [a for a in sys.argv[1:] if not a.startswith('--')]
budget = float(next((a[9:] for a in sys.argv[1:] if a.startswith('--budget=')), 1800))
P = resolve_tol(args[0])
res = fh.solve(P, log=lambda s: None)
assert res['status'] == 'ok', res['status']
print('cost', res['cost'], 'lower bound', res['lower_bound'], flush=True)
out = {}

def blocks(folds):
    """The partition of the visible creases induced by the folds, and crease -> block."""
    B = [frozenset(k for k in ks if isinstance(k, int)) for ks in folds.values()]
    B = [b for b in B if b]
    return B, {k: b for b in B for k in b}

def show(name, r, **kw):
    t0 = time.time()
    folds, lt, val, paths = fh.analyse(P, r, budget=budget, **kw)
    out[name] = (folds, lt)
    vis = [x for x, ks in folds.items() if any(isinstance(k, int) for k in ks)]
    G = nx.DiGraph(); G.add_nodes_from(folds); G.add_edges_from(lt)
    print('%-10s folds %d (with visible creases %d), relations %d, Hasse %d, longest chain %d, undecided queries %d, %.0fs' % (
        name, len(folds), len(vis), len(lt), len(fh.reduce(list(folds), lt)), nx.dag_longest_path_length(G) + 1,
        r['analysis_unknown'], time.time() - t0), flush=True)

show('directions', res, pin=False)
show('found', res)
G, complete = geometries(P, res, budget=budget)
print('minimal geometries: %d%s' % (len(G), '' if complete else ' (INCOMPLETE)'), flush=True)
if complete:
    g = min(G, key=lambda g: g['length'])
    show('shortest', dict(res, m=g['m'], choice=g['choice']))
show('minimal', res, scope='minimal')

# consistency: 'minimal' ranges over a superset of the resolutions of 'directions', so its
# partition of the visible creases refines that of 'directions', and its order is contained in it
(fd, ld), (fm, lm) = out['directions'], out['minimal']
Bd, bd = blocks(fd); Bm, bm = blocks(fm)
assert all(any(b <= c for c in Bd) for b in Bm), 'minimal folds do not refine the directions folds'
vis = lambda folds, x: [k for k in folds[x] if isinstance(k, int)]
Ld = {(bd[vis(fd, x)[0]], bd[vis(fd, y)[0]]) for x, y in ld if vis(fd, x) and vis(fd, y)}
for x, y in lm:
    assert (bd[vis(fm, x)[0]], bd[vis(fm, y)[0]]) in Ld, 'a minimal relation fails under directions'
print('consistency checks passed', flush=True)
