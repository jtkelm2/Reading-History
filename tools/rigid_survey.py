"""How many of the tables' local histories are rigidly realizable (src/rigid.py)?

For every vertex shape of the given patterns (each shape once, up to rotation and reflection) and
every number k of hidden rays the tables cover, count the histories the tables give at the vertex
and those that are rigidly realizable.  Prints one JSON line per shape, then a summary.
Usage: python3 tools/rigid_survey.py TABLEDIR [--reverse] pattern.fold [...]
--reverse also allows reverse folds."""
import sys, json, time
from common import fh, resolve_tol
sys.path.insert(0, str(fh.__file__).rsplit('/', 1)[0])
import lookup, rigid

REV = '--reverse' in sys.argv
sys.argv = [a for a in sys.argv if a != '--reverse']
T = lookup.Table(sys.argv[1])
kmax = {}
for m in T.meta: kmax[m['n']] = max(kmax.get(m['n'], -1), m['k'])
seen, rows = set(), []
for path in sys.argv[2:]:
    P = resolve_tol(path); opt = lookup.options_for(T, P)
    for p in P.interior:
        n = len(P.H[p])
        if n not in kmax: continue
        D = sorted(d for d, _ in P.H[p])
        secs = [round(P.deg(D[(i + 1) % n] - D[i]), 4) for i in range(n)]
        shape = min(tuple(s_[i:] + s_[:i]) for s_ in (secs, secs[::-1]) for i in range(n))
        if shape in seen: continue
        seen.add(shape)
        dirs = {e: d for d, e in P.H[p]}
        row = {'file': path.rsplit('/', 1)[-1], 'vertex': p, 'sectors': list(shape), 'by_k': {}}
        t0 = time.monotonic()
        for k in range(kmax[n] + 1):
            hs = [h for _, H, _ in opt(p, k, 0) for h in H]
            ok = 0
            for h in hs:
                steps = [[dirs[x] if isinstance(x, int) else x[0][1] for x in st] for st in h]
                ok += rigid.realizable(steps, reverse=REV)[0]
            if hs: row['by_k'][k] = [len(hs), ok]
        row['seconds'] = round(time.monotonic() - t0, 2)
        rows.append(row); print(json.dumps(row), flush=True)

tot = {}
for r in rows:
    for k, (a, b) in r['by_k'].items(): t = tot.setdefault(k, [0, 0]); t[0] += a; t[1] += b
first = lambda r, i: min((int(k) for k, v in r['by_k'].items() if v[i] > 0), default=None)
print(json.dumps({'shapes': len(rows), 'histories_by_k (total, rigid)': tot,
                  'shapes whose least k rises under rigidity': [r['sectors'] for r in rows if first(r, 0) != first(r, 1)],
                  'shapes with no rigid history within the tables': [r['sectors'] for r in rows if r['by_k'] and first(r, 1) is None]}))
