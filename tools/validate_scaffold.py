"""Compare the scaffold types of the lookup tables with the solver's own local search, vertex by vertex.

At every interior vertex whose degree the tables cover, and every number k of hidden rays they
cover, compare the local histories that have a scaffold step (a step making only hidden rays):
  table     src/lookup.py, scaffold types (*_s.json), replayed at the vertex;
  search    Pattern.options (lazy_histories: scaffold lines from bisector images), restricted to
            histories whose scaffold lines are all pinned (tools/pinning.py, free_dims = 0).
The search's lines are a special case, so search ⊆ table is expected. Every table history is also
checked to be pinned by the numeric measure, independently of the elimination that admitted it.
Each vertex shape is compared once over all patterns. Prints one JSON line per pattern.
Usage: python3 tools/validate_scaffold.py TABLEDIR pattern.fold [...]"""
import sys, json, time
from common import fh, resolve_tol
sys.path.insert(0, str(fh.__file__).rsplit('/', 1)[0])
import lookup, pinning


def canon(P, h):
    return tuple(frozenset(k if isinstance(k, int) else ('h', round(P.deg(k[0][1]), 4), k[1]) for k in st) for st in h)


T = lookup.Table(sys.argv[1])
kmax = {}
for m in T.meta:
    if m.get('scaffold'): kmax[m['n']] = max(kmax.get(m['n'], -1), m['k'])
seen = set()
for path in sys.argv[2:]:
    P = resolve_tol(path); opt = lookup.options_for(T, P)
    out = {'file': path.rsplit('/', 1)[-1], 'vertices': 0, 'compared': 0, 'search_subset': 0, 'search': 0, 'table': 0,
           'only_search': 0, 'table_unpinned': 0, 'search_lazy': 0, 'examples': []}
    t0 = time.time()
    for p in P.interior:
        n = len(P.H[p])
        if n not in kmax: continue
        D = sorted(d for d, _ in P.H[p])
        secs = [round(P.deg(D[(i + 1) % n] - D[i]), 4) for i in range(n)]
        shape = min(tuple(s_[i:] + s_[:i]) for s_ in (secs, secs[::-1]) for i in range(n))
        if shape in seen: continue
        seen.add(shape); out['vertices'] += 1
        for k in range(kmax[n] + 1):
            old = P.__class__.options(P, p, k, 0)
            if old and old[0][1] == 'LAZY': out['search_lazy'] += 1; continue
            A = {canon(P, h) for _, hs, *_ in old for h in hs if pinning.scaffold_steps(h) and pinning.free_dims(P, p, h) == 0}
            Bh = [h for _, hs, *_ in opt(p, k, 0) for h in hs if pinning.scaffold_steps(h)]
            B = {canon(P, h) for h in Bh}
            out['table_unpinned'] += sum(1 for h in Bh if pinning.free_dims(P, p, h) != 0)
            out['compared'] += 1; out['search'] += len(A); out['table'] += len(B)
            if A <= B: out['search_subset'] += 1
            else:
                out['only_search'] += len(A - B)
                if len(out['examples']) < 3:
                    out['examples'].append({'vertex': p, 'k': k, 'rays_deg': [P.deg(d) for d, _ in P.H[p]],
                                            'only_search': [str(x) for x in list(A - B)[:2]]})
    out['seconds'] = round(time.time() - t0, 1)
    print(json.dumps(out), flush=True)
