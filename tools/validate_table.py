"""Compare the lookup tables with the solver's own local search, vertex by vertex.

At every interior vertex of degree covered by the tables, and every number k of hidden rays the
tables cover, compare the set of local histories:
  table     src/lookup.py, anchored types (no step making only hidden rays), replayed;
  search    Pattern.options (k = 0: histories(), complete; k >= 1: lazy_histories), restricted
            likewise. Scaffold types are compared by tools/validate_scaffold.py.
Hidden creases are compared by direction (degrees, 4 decimals) and by first/cancelling crease.
Each vertex shape (sector angles up to rotation and reflection) is compared once over all
patterns. Prints one JSON line per pattern.
Usage: python3 tools/validate_table.py TABLEDIR pattern.fold [...]"""
import sys, json, time
from common import fh, resolve_tol
sys.path.insert(0, str(fh.__file__).rsplit('/', 1)[0])
import lookup


def canon(P, h):
    return tuple(frozenset(k if isinstance(k, int) else ('h', round(P.deg(k[0][1]), 4), k[1]) for k in st) for st in h)


def hidden_only(h):
    return any(all(not isinstance(k, int) and k[1] == 1 for k in st) for st in h)


T = lookup.Table(sys.argv[1])
kmax = {}
for m in T.meta: kmax[m['n']] = max(kmax.get(m['n'], -1), m['k'])
seen = set(); failed = False
for path in sys.argv[2:]:
    P = resolve_tol(path); opt = lookup.options_for(T, P)
    out = {'file': path.rsplit('/', 1)[-1], 'vertices': 0, 'compared': 0, 'equal': 0, 'only_search': 0, 'only_table': 0,
           'search_lazy': 0, 'examples': []}
    t0 = time.time()
    for p in P.interior:
        n = len(P.H[p])
        if n not in kmax: continue
        D = sorted(d for d, _ in P.H[p])
        secs = [round(P.deg(D[(i + 1) % n] - D[i]), 4) for i in range(n)]
        shape = min(tuple(s_[i:] + s_[:i]) for s_ in (secs, secs[::-1]) for i in range(n))
        if shape in seen: continue                      # each vertex shape once (up to symmetry)
        seen.add(shape); out['vertices'] += 1
        for k in range(kmax[n] + 1):
            old = P.__class__.options(P, p, k, 0)
            if old and old[0][1] == 'LAZY': out['search_lazy'] += 1; continue
            A = {canon(P, h) for _, hs, *_ in old for h in hs if not hidden_only(h)}
            B = {canon(P, h) for _, hs, *_ in opt(p, k, 0) for h in hs if not hidden_only(h)}   # anchored types only
            out['compared'] += 1
            if A == B: out['equal'] += 1
            else:
                out['only_search'] += len(A - B); out['only_table'] += len(B - A)
                if len(out['examples']) < 3:
                    out['examples'].append({'vertex': p, 'k': k, 'rays_deg': [P.deg(d) for d, _ in P.H[p]],
                                            'only_search': [str(x) for x in list(A - B)[:2]], 'only_table': [str(x) for x in list(B - A)[:2]]})
    out['seconds'] = round(time.time() - t0, 1)
    print(json.dumps(out), flush=True)
    failed = failed or out['equal'] != out['compared']
print('all equal' if not failed else 'DIFFERENCES FOUND')
sys.exit(1 if failed else 0)
