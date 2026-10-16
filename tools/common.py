"""Shared helpers for the tools: import the solver from src/ and resolve repo paths."""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))
import fold_history as fh  # noqa: E402


def resolve(path, tol=None):
    """Pattern, solve result, and (folds, forced order, times, cancelled paths) for one file."""
    tols = (tol,) if tol else (1e-6, 1e-4, 5e-4, 2e-3)
    for t in tols:
        P = fh.Pattern(str(path), t)
        if all(fh.maps(sorted(d for d, _ in P.H[p])) is not None for p in P.interior): break
    res = fh.solve(P, log=lambda s: None)
    assert res['status'] == 'ok', res['status']
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
