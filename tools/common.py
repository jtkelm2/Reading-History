"""Shared helpers for the tools: import the solver from src/ and resolve repo paths."""
import sys, math, pathlib
ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))
import fold_history as fh  # noqa: E402
from fold_history import geometries, segments  # noqa: E402
z3 = fh.z3


def resolve_tol(path, tol=None):
    """The pattern, read with the first tolerance at which every interior vertex is flat."""
    for t in ((tol,) if tol else (1e-6, 1e-4, 5e-4, 2e-3)):
        P = fh.Pattern(str(path), t)
        if all(fh.maps(sorted(d for d, _ in P.H[p])) is not None for p in P.interior): break
    return P


def resolve(path, tol=None, geometry='found', table=None):
    """Pattern, solve result, and (folds, forced order, times, cancelled paths) for one file.
    geometry='shortest' analyses the minimal geometry whose cancelled segments are shortest.
    table=DIR takes the local histories from the lookup tables in DIR (src/lookup.py)."""
    P = resolve_tol(path, tol)
    if table:
        import lookup
        P.options = lookup.options_for(lookup.Table(table), P)
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
