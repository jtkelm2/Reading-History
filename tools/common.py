"""Shared helpers for the tools: import the solver from src/ and resolve repo paths."""
import sys, math, pathlib
ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))
import fold_history as fh  # noqa: E402
z3 = fh.z3


def resolve_tol(path, tol=None):
    """The pattern, read with the first tolerance at which every interior vertex is flat."""
    for t in ((tol,) if tol else (1e-6, 1e-4, 5e-4, 2e-3)):
        P = fh.Pattern(str(path), t)
        if all(fh.maps(sorted(d for d, _ in P.H[p])) is not None for p in P.interior): break
    return P


def segments(P, M, m, choice):
    """The cancelled segments of model m, as a set of edges between rounded points."""
    E = set()
    for p, R in choice.items():
        for r in R:
            pts, _ = fh.walk(P, M, m, P.V[p], r, (p, r), p, None)
            for a, b in zip(pts, pts[1:]):
                E.add(frozenset(((round(a[0], 6), round(a[1], 6)), (round(b[0], 6), round(b[1], 6)))))
    return E


def geometries(P, res, cap=1000, budget=1800):
    """Minimal-cost geometries within the search space: dicts with the model m, the hidden
    directions (choice), a hashable selection, and the total length of the segments.  Stops after
    cap geometries or budget seconds; the second value says whether the list is complete."""
    import time
    M, opts = res['model'], res['opts']; s = M.s; out = []; t0 = time.time(); r = None
    s.push(); s.add(M.cost() == res['cost'])
    while len(out) < cap:
        s.set('timeout', max(1, int(1000 * (budget - (time.time() - t0)))))
        r = s.check(list(M.track))
        if r != z3.sat: break
        m = s.model()
        if fh.refine(P, M, opts, m, [s]): continue
        on = [(p, i) for (p, i), b in M.sel.items() if z3.is_true(m.eval(b))]
        choice = {p: opts[p][i][0] for p, i in on if opts[p][i][0]}
        E = segments(P, M, m, choice)
        out.append({'m': m, 'choice': choice, 'length': sum(math.dist(*tuple(e)) for e in E),
                    'selection': frozenset((p, tuple(P.deg(r) for r in R)) for p, R in choice.items())})
        s.add(z3.Not(z3.And([M.sel[k] for k in on] + [M.geometry(choice, m)])))
    s.pop()
    return out, r == z3.unsat


def resolve(path, tol=None, geometry='found'):
    """Pattern, solve result, and (folds, forced order, times, cancelled paths) for one file.
    geometry='shortest' analyses the minimal geometry whose cancelled segments are shortest."""
    P = resolve_tol(path, tol)
    res = fh.solve(P, log=lambda s: None)
    assert res['status'] == 'ok', res['status']
    if geometry == 'shortest':
        G, complete = geometries(P, res)
        assert complete, 'the minimal geometries were not all enumerated'
        g = min(G, key=lambda g: g['length'])
        res = dict(res, m=g['m'], choice=g['choice'])
    return P, res, fh.analyse(P, res)


def linear_extensions(nodes, lt):
    idx = {x: i for i, x in enumerate(nodes)}; n = len(nodes)
    pred = [0] * n
    for a, b in lt: pred[idx[b]] |= 1 << idx[a]
    cnt = [0] * (1 << n); cnt[0] = 1
    for m in range(1 << n):
        if cnt[m]:
            for v in range(n):
                if not m >> v & 1 and pred[v] & m == pred[v]: cnt[m | 1 << v] += cnt[m]
    return cnt[-1]
