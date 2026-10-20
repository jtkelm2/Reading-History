"""Local freedom of the scaffold lines of a local history (keys as in Pattern.options).

A scaffold step makes only hidden rays.  Give each scaffold step a line parameter; every other step's
line is fixed by an anchor (a visible ray, or a hidden ray being erased), and new hidden rays are
preimages of the step's line.  The residuals of (L) and of flatness are piecewise linear in the
parameters.  free_dims = (number of scaffold steps) - (rank of their Jacobian at the history):
0 means every scaffold line is pinned (locally isolated), > 0 means a continuum of valid lines.
Computed numerically (finite differences, eps = 1e-6); a reported free direction is confirmed by a
finite step along it that keeps every residual 0."""
import math
import numpy as np
from common import fh

PSI = 0.1234567     # a generic rotation: no ray sits on the cut at angle 0, so small perturbations
                    # never change the sector the folding maps start from


def _maps(S):
    """Folding maps of the sorted state S without the closure check, and the closure defect (mod 2pi)."""
    out = [(1, 0.0)]
    for c in S[1:]:
        e, o = out[-1]; out.append((-e, fh.nrm(2 * e * c + o)))
    e, o = out[-1]
    v = 2 * e * S[0] + o
    return out, (v + math.pi) % (2 * math.pi) - math.pi


def _replay(P, p, h, delta):
    """Residuals of the history h at vertex p with scaffold step i's line moved by delta[i]."""
    dirs = {e: d for d, e in P.H[p]}
    steps = [[(kk, fh.nrm((dirs[kk] if kk in dirs else kk[0][1]) + PSI)) for kk in T] for T in h]
    d, S, res = {}, [], []
    for i, T in enumerate(steps):
        Sd = sorted(x for _, x in S)
        phis = None
        if Sd:
            phis, defect = _maps(Sd); res.append(defect)
        img = (lambda x: fh.apply(phis[fh.sector_of(x, Sd)], x)) if Sd else (lambda x: x)
        cur = lambda kk, x: d[kk[0][1]] if isinstance(kk, tuple) and kk[0][1] in d else x
        anchor = next(((kk, x) for kk, x in T if not isinstance(kk, tuple) or kk[1] == 2), None)
        lam = img(cur(*anchor)) if anchor else img(T[0][1]) + delta.get(i, 0.0)
        for kk, x in T:
            if isinstance(kk, tuple) and kk[1] == 1:
                if Sd:
                    f = phis[fh.sector_of(x, Sd)]
                    cands = [fh.inverse_apply(f, lam), fh.inverse_apply(f, fh.nrm(lam + math.pi))]
                else:
                    cands = [fh.nrm(lam), fh.nrm(lam + math.pi)]
                d[kk[0][1]] = min(cands, key=lambda c: abs((c - x + math.pi) % (2 * math.pi) - math.pi))
        for kk, x in T:
            v = img(cur(kk, x)) - lam
            res.append((v + math.pi / 2) % math.pi - math.pi / 2)
        for kk, x in T:                     # rays tagged: a crease id 0 must not equal a direction 0.0
            if not isinstance(kk, tuple): S.append((('v', kk), x))
            elif kk[1] == 1: S.append((('h', kk[0][1]), d[kk[0][1]]))
            else: S = [(a, y) for a, y in S if a != ('h', kk[0][1])]
    return res


def scaffold_steps(h):
    return [i for i, T in enumerate(h) if all(isinstance(kk, tuple) and kk[1] == 1 for kk in T)]


def free_dims(P, p, h, eps=1e-6):
    sc = scaffold_steps(h)
    if not sc: return 0
    r0 = _replay(P, p, h, {})
    assert max(map(abs, r0)) < 1e-6, 'the history does not replay'
    cols = []
    for i in sc:
        r = _replay(P, p, h, {i: eps})
        cols.append([(a - b) / eps for a, b in zip(r, r0)])
    J = np.array(cols).T
    rank = np.linalg.matrix_rank(J, tol=1e-3) if J.size else 0
    if rank < len(sc):
        v = np.linalg.svd(J)[2][-1]
        r = _replay(P, p, h, {i: 1e-3 * x for i, x in zip(sc, v)})
        assert max(map(abs, r)) < 1e-8, 'free direction not confirmed'
    return len(sc) - rank
