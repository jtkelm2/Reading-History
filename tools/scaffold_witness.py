"""Do pinned scaffold lines repair a pattern that the physical tables cannot resolve globally?

Options at each vertex: the tables filtered by the exact cells (realizable with reverse folds), at
every depth the tables cover; with --scaffolds also the histories of the solver's own local search
(Pattern.options) that have a scaffold step (a step making only hidden rays; the search takes its line
from bisector images), are realizable with reverse folds, and whose scaffold lines are all pinned
(tools/pinning.py, free_dims = 0), within the same depths.

  --need=V1,V2   check that only these vertices can choose consistently (all options in the model).
  --witness      look for a global resolution: start from the solution with the unfiltered tables,
                 pin each vertex whose chosen history is physical to its hidden directions, and drop
                 the pins in each unsat core until sat (a witness, verified against every assertion)
                 or a core without pins (unsat over the whole option set).
  --arith=3      z3's dense difference-logic solver (default; the constraints are difference
                 constraints on integer times); --arith=default for z3's own choice.
  --seconds=S    per check (default 240).   --out=PATH  write the witness as JSON.
Usage: python3 tools/scaffold_witness.py FILE [--scaffolds] [--need=...|--witness] [...]"""
import sys, math, time, json
from common import fh, resolve_tol, ROOT
sys.path.insert(0, str(ROOT / 'src'))
import lookup, rigid, pinning
z3 = fh.z3

path = sys.argv[1]
fl = dict((a[2:].split('=', 1) + ['1'])[:2] for a in sys.argv[2:])
if fl.get('arith', '3') != 'default': z3.set_param('smt.arith.solver', int(fl.get('arith', '3')))
t0 = time.monotonic()
P = resolve_tol(path)
N = P.N; U = 2 * math.pi / N if N else None
deg = lambda x: round(math.degrees(fh.nrm(x)), 4)
search = P.options
tables = str(ROOT / 'tables')
T = lookup.Table(tables, 'reverse', str(ROOT / 'cells'))
kmax = {}
for m in T.meta: kmax[m['n']] = max(kmax.get(m['n'], -1), m['k'])
tab = lookup.options_for(T, P)


def steps_of(p, h):
    dirs = {e: d for d, e in P.H[p]}
    return [[dirs[kk] if kk in dirs else kk[0][1] for kk in S] for S in h]


def describe(p, h):
    """Rays made / erased per step, in degrees."""
    out = []
    for S, ds in zip(h, steps_of(p, h)):
        mk = sorted(deg(d) for kk, d in zip(S, ds) if not (isinstance(kk, tuple) and kk[1] == 2))
        er = sorted(deg(d) for kk, d in zip(S, ds) if isinstance(kk, tuple) and kk[1] == 2)
        out.append('%s/%s' % (mk, er) if er else '%s' % mk)
    return ' ; '.join(out)


opts, src, added = {}, {}, 0
for p in P.interior:
    n = len(P.H[p]); ol = []
    for k in range(kmax.get(n, -1) + 1):
        ol += list(tab(p, k))
        if 'scaffolds' in fl:
            for R, hs, col in search(p, k, 0):
                if hs == 'LAZY': continue
                ok = [h for h in hs if pinning.scaffold_steps(h) and rigid.realizable_fast(steps_of(p, h), reverse=True)[0]
                      and pinning.free_dims(P, p, h) == 0]
                if ok: ol.append((R, ok, col)); src[p, len(ol) - 1] = 'scaffold'; added += len(ok)
    opts[p] = ol
if 'solve' in fl:                       # the solver's own search (core-guided, minimum cost) on these options
    byk = {}
    for p, ol in opts.items():
        for R, hs, col in ol: byk.setdefault((p, len(R)), []).append((R, hs, col))
    P.options = lambda p, k, c=0: [] if c else byk.get((p, k), [])
    res = fh.solve(P, log=lambda s: None, budget=float(fl.get('seconds', 240)))
    out = {k: res[k] for k in ('status', 'cost', 'lower_bound', 'certified_minimal', 'core') if k in res}
    if res['status'] == 'ok':
        sid = {id(opts[p][i][0]) for p, i in src}
        out['scaffold_vertices'] = sorted(p for (p, i), b in res['model'].sel.items()
                                          if z3.is_true(res['m'].eval(b)) and id(res['opts'][p][i][0]) in sid)
    print(json.dumps(dict(out, seconds=round(time.monotonic() - t0, 1)))); sys.exit()
M = fh.Model(P, opts)
print(json.dumps({'file': path.rsplit('/', 1)[-1], 'options': sum(len(o) for o in opts.values()),
                  'scaffold_histories': added, 'seconds': round(time.monotonic() - t0, 1)}), flush=True)
inv = {v: a for a, v in M.track.items()}
M.s.set('timeout', int(1000 * float(fl.get('seconds', 240))))

if 'need' in fl:
    need = [int(x) for x in fl['need'].split(',')]
    r = M.s.check([inv[p] for p in need])
    print(json.dumps({'need': need, 'result': str(r), 'seconds': round(time.monotonic() - t0, 1)}), flush=True)
    if r == z3.sat: sys.exit()
    if r == z3.unsat: print('core', sorted(M.track[c] for c in M.s.unsat_core()))
    sys.exit()

# --witness
P2 = resolve_tol(path); P2.options = lookup.options_for(lookup.Table(tables), P2)
res = fh.solve(P2, log=lambda s: None, budget=300)
assert res['status'] == 'ok', res['status']
key = lambda R: tuple(sorted(round(x, 6) for x in R))
pins = {}
for (p, i), b in res['model'].sel.items():
    if not z3.is_true(res['m'].eval(b)): continue
    R, hs, _ = res['opts'][p][i]
    h = next(h for h in hs if z3.is_true(res['m'].eval(res['model'].history(h), model_completion=True)))
    if not rigid.realizable_fast(steps_of(p, h), reverse=True)[0]: continue
    lits = [M.sel[p, j] for j, o in enumerate(opts[p]) if key(o[0]) == key(R)]
    if lits:
        a = z3.Bool('pin_%d' % p); M.s.add(z3.Implies(a, z3.Or(lits))); pins[p] = a
while True:
    r = M.s.check(list(M.track) + list(pins.values()))
    print(json.dumps({'pins': sorted(pins), 'result': str(r), 'seconds': round(time.monotonic() - t0, 1)}), flush=True)
    if r != z3.unsat: break
    core = [str(c) for c in M.s.unsat_core()]
    drop = [int(c[4:]) for c in core if c.startswith('pin_')]
    if not drop:
        print('unsat over the whole option set; core', sorted(int(c[4:]) for c in core)); break
    for p in drop: del pins[p]
if r == z3.sat:
    m = M.s.model()
    bad = [a for a in M.s.assertions() if not z3.is_true(m.eval(a, model_completion=True))]
    bad += [a for a in list(M.track) + list(pins.values()) if not z3.is_true(m.eval(a, model_completion=True))]
    assert not bad, 'model does not satisfy %d assertions' % len(bad)
    W = []
    for p in sorted(opts):
        for i, o in enumerate(opts[p]):
            if z3.is_true(m.eval(M.sel[p, i])):
                h = next(h for h in o[1] if z3.is_true(m.eval(M.history(h), model_completion=True)))
                st = steps_of(p, h)
                W.append({'vertex': p, 'degree': len(P.H[p]), 'rays': [deg(d) for d, _ in P.H[p]],
                          'source': src.get((p, i), 'table'), 'hidden': sorted(deg(x) for x in o[0]),
                          'history': describe(p, h), 'rigid': rigid.realizable_fast(st)[0],
                          'reverse': rigid.realizable_fast(st, reverse=True)[0], 'free_dims': int(pinning.free_dims(P, p, h))})
    assert all(w['reverse'] for w in W) and all(w['free_dims'] == 0 for w in W)
    sc = [w['vertex'] for w in W if w['source'] == 'scaffold']
    print(json.dumps({'witness': 'verified against every assertion', 'hidden_rays': sum(len(w['hidden']) for w in W),
                      'scaffold_vertices': sc}))
    for w in W:
        if w['source'] == 'scaffold': print('  v%d (degree %d, rays %s): %s' % (w['vertex'], w['degree'], w['rays'], w['history']))
    if 'out' in fl: json.dump({'file': path.rsplit('/', 1)[-1], 'vertices': W}, open(fl['out'], 'w'), indent=1)
