"""Rigid realizability of a local history (flatvertex.tex, Definitions 2.1, 2.2, 3.1, 3.2).

A layered state is a flat state S with, at every point u of the image circle outside the image
points E_S, a strict total order (bottom to top) on the faces whose image contains u, locally
constant in u, satisfying (O1)-(O3). A rigid step folds a step D along its line l by a half-turn:
it is given by a labelling (which pieces move) and a choice of the up component of the circle
minus the two rays of l. With senses (reverse folds), every moving piece is lifted over its own
side; a moving run keeps the common sense or is split at one spine into halves with opposite
senses. A history is realizable if some choices make every step rigid.

Orders are stored on the arcs between consecutive image points. All geometric decisions go
through an arithmetic object (src/geom.py): Num() for floats in radians, or Sym(...) for exact
affine forms whose signs an oracle decides; the same code serves both.

Usage: python3 src/rigid.py   (checks the three histories of flatvertex.tex, Section 4)
"""
import itertools, math
from functools import cmp_to_key
import geom

NUM = geom.Num()


def arcs_of(B):
    """Consecutive arcs (a, b) of the circle cut at the sorted points B; one full arc if B is empty."""
    if not B: return [(None, None)]
    return [(B[i], B[(i + 1) % len(B)]) for i in range(len(B))]


def inside(ar, u, arc):
    a, b = arc
    return True if a is None else geom.between(ar, u, a, b)


class State:
    """S: sorted directions; phis: folding maps; arcs: [(a, b, order)] covering the circle."""
    def __init__(self, ar, S, arcs):
        self.ar, self.S = ar, S
        self.phis = geom.maps(ar, S) if S else [(1, ar.ZERO)]
        self.arcs = arcs
    def face_of(self, x):
        return geom.sector_of(self.ar, x, self.S) if self.S else 0
    def image(self, x):
        return geom.apply(self.ar, self.phis[self.face_of(x)], x)
    def E(self):
        return geom.uniq(self.ar, [self.image(s) for s in self.S])
    def order_at(self, u):
        for a, b, o in self.arcs:
            if inside(self.ar, u, (a, b)): return o
        raise ValueError('no arc contains %r' % (u,))


def initial(ar=NUM):
    return State(ar, [], [(None, None, (0,))])


def img_arc(ar, f, a, b):
    """Image of the open arc (a, b) under x -> e x + o, as a counterclockwise arc."""
    e, o = f
    if a is None: return (None, None)
    x, y = geom.apply(ar, f, a), geom.apply(ar, f, b)
    return (x, y) if e == 1 else (y, x)


def runs(R, D, mu):
    """The moving runs of a step: maximal cyclic sequences of consecutive moving pieces (indices
    into the pieces (R[k], R[k+1])). Consecutive pieces of a run meet at a point of S not in D."""
    n = len(R); out = []
    starts = [k for k in range(n) if mu[k] and not mu[(k - 1) % n]]
    for k0 in starts:
        run, k = [], k0
        while mu[k] and (not run or k != k0):
            run.append(k); k = (k + 1) % n
        out.append(run)
    return out


def pieces_and_labels(st, D, mu0):
    ar = st.ar
    D = [geom.nrm(ar, d) for d in D]
    R = geom.uniq(ar, list(st.S) + D)
    inD = lambda x: any(geom.eq(ar, x, d) for d in D)
    n = len(R); mu, cur = [], mu0
    for k in range(n):
        mu.append(cur)
        if inD(R[(k + 1) % n]): cur ^= 1
    return D, R, inD, mu, cur == mu0


def sense_options(st, D, mu0):
    """Sense assignments allowed with labelling mu0: each moving run either keeps the common sense
    0, or is split at one spine (a point between two consecutive pieces of the run) into two halves
    with opposite senses (either half first). Yields a list sense[k] over all pieces."""
    D, R, inD, mu, ok = pieces_and_labels(st, D, mu0)
    n = len(R)
    per_run = []
    for run in runs(R, D, mu):
        opts = [None] + [(j, o) for j in range(1, len(run)) for o in (0, 1)]
        per_run.append((run, opts))
    for choice in itertools.product(*[opts for _, opts in per_run]):
        sense = [0] * n
        for (run, _), c in zip(per_run, choice):
            if c is None: continue
            j, o = c
            for i, k in enumerate(run): sense[k] = o if i < j else 1 - o
        yield sense


def step(st, D, mu0, up, sense=None, why=None):
    """The layered state after the step D with labelling mu0, up side up and senses (default: all
    pieces sense 0, i.e. the rigid half-turn of flatvertex Definition 2.2), or None.
    A moving piece with sense 0 is lifted over the up side, one with sense 1 over the other side.
    If why is a list, the reason for a rejection is appended to it."""
    fail = lambda r: (why.append(r) if why is not None else None)
    ar, S = st.ar, st.S
    D, R, inD, mu, ok = pieces_and_labels(st, D, mu0)
    if not ok: return None                                       # |D| odd: not a step
    n = len(R)
    pieces = [(R[k], R[(k + 1) % n]) for k in range(n)]
    parent = [st.face_of(R[k]) for k in range(n)]
    lam = geom.modpi(ar, st.image(D[0]))
    L2 = [lam, geom.nrm(ar, lam + ar.PI)]
    uparc = (L2[0], L2[1]) if up == 0 else (L2[1], L2[0])
    downarc = (uparc[1], uparc[0])
    if sense is None: sense = [0] * n
    U = [uparc if not sense[k] else downarc for k in range(n)]           # the side each piece is lifted over
    sig = lambda x: geom.nrm(ar, 2 * lam - x)
    phi_piece = [img_arc(ar, st.phis[parent[k]], *pieces[k]) for k in range(n)]
    def psi_piece(k):
        a, b = phi_piece[k]
        return (a, b) if not mu[k] else (sig(b), sig(a))
    psi = [psi_piece(k) for k in range(n)]
    ES = st.E()
    # (F') moving against staying, (F'') moving pieces of opposite senses against each other
    for a, b in arcs_of(geom.uniq(ar, ES + L2)):
        u = geom.mid(ar, a, b)
        P = [k for k in range(n) if inside(ar, u, phi_piece[k])]
        mov = [k for k in P if mu[k]]; sta = [k for k in P if not mu[k]]
        if not mov: continue
        pos = {f: i for i, f in enumerate(st.order_at(u))}
        for p in mov:
            lifted = inside(ar, u, U[p])
            for q in sta:
                if lifted and not pos[parent[q]] < pos[parent[p]]: return fail("(F')")
                if not lifted and not pos[parent[p]] < pos[parent[q]]: return fail("(F')")
            for q in mov:
                if sense[q] != sense[p] and lifted and not pos[parent[q]] < pos[parent[p]]: return fail("(F'')")
    # new order, in the coordinates of Psi
    S2 = geom.uniq(ar, [x for x in S if not inD(x)] + [d for d in D if not any(geom.eq(ar, d, x) for x in S)])
    new = State(ar, S2, None)
    G = [new.face_of(R[k]) for k in range(n)]
    fine = []
    for a, b in arcs_of(geom.uniq(ar, ES + [sig(x) for x in ES] + L2)):
        u = geom.mid(ar, a, b)
        P = [k for k in range(n) if inside(ar, u, psi[k])]
        sta = [k for k in P if not mu[k]]
        below = [k for k in P if mu[k] and inside(ar, u, U[k])]         # arrived from underneath
        above = [k for k in P if mu[k] and not inside(ar, u, U[k])]     # arrived over the top
        if sta:
            pos = {f: i for i, f in enumerate(st.order_at(u))}
            sta.sort(key=lambda k: pos[parent[k]])                       # rule (i)
        if below or above:
            pos = {f: i for i, f in enumerate(st.order_at(sig(u)))}
            below.sort(key=lambda k: -pos[parent[k]]); above.sort(key=lambda k: -pos[parent[k]])  # rule (ii)
        seq = below + sta + above                                        # rule (iii)
        order = tuple(G[k] for k in seq)
        if len(set(order)) != len(order): raise AssertionError('Psi is not injective on a face')
        fine.append((a, b, order))
    # g with Phi_S2 = g o Psi (checked on every piece)
    def psi_map(k):
        e, o = st.phis[parent[k]]
        return (e, o) if not mu[k] else (-e, geom.nrm(ar, 2 * lam - o))
    gs = []
    for k in range(n):
        ep, op = psi_map(k); e2, o2 = new.phis[G[k]]
        gs.append((e2 * ep, geom.nrm(ar, o2 - e2 * ep * op)))
    eg, og = gs[0]
    if any(e != eg or not geom.eq(ar, o, og) for e, o in gs): raise AssertionError('Phi_S2 is not g o Psi')
    g = lambda y: geom.nrm(ar, eg * y + og)
    moved = []
    for a, b, o in fine:
        if a is None: moved.append((None, None, o)); continue
        moved.append((g(a), g(b), o) if eg == 1 else (g(b), g(a), o))
    # extend to a layer order: constant across breakpoints that are not image points of S2
    E2 = new.E()
    isE = lambda x: any(geom.eq(ar, x, v) for v in E2)
    if moved[0][0] is not None:
        moved.sort(key=cmp_to_key(lambda s, t: ar.sgn(s[0] - t[0])))
        m = len(moved)
        start = next((i for i in range(m) if isE(moved[i][0])), None)
        if start is None:                                                # no image points: one order
            if any(o != moved[0][2] for _, _, o in moved): return fail('no extension')
            arcs = [(None, None, moved[0][2])]
        else:
            arcs, i = [], start
            for _ in range(m):
                a, b, o = moved[i]
                if isE(a): arcs.append([a, b, o])
                else:
                    if o != arcs[-1][2]: return fail('no extension')
                    arcs[-1][1] = b
                i = (i + 1) % m
            arcs = [tuple(x) for x in arcs]
    else:
        arcs = moved
    new.arcs = arcs
    r = layered_ok(new, True)
    return new if r is True else fail(r)


def layered_ok(st, reason=False):
    """(O1)-(O3) for a state whose orders are stored on the arcs between its image points.
    With reason, returns True or the name of the first condition that fails."""
    no = lambda r: r if reason else False
    ar, S = st.ar, st.S
    if not S: return True
    m = len(S)
    by_v = []                                                            # [(v, rays)]
    for idx, s in enumerate(S):
        v = st.image(s)
        jb = (idx - 1) % m                                               # face before s
        eps = '-' if st.phis[jb][0] == 1 else '+'
        hit = next((g for g in by_v if geom.eq(ar, g[0], v)), None)
        if hit is None: by_v.append((v, [((jb, idx), eps)]))
        else: hit[1].append(((jb, idx), eps))
    for v, rays in by_v:
        minus = next(o for a, b, o in st.arcs if a is not None and geom.eq(ar, b, v))
        plus = next(o for a, b, o in st.arcs if a is not None and geom.eq(ar, a, v))
        theta = set(minus) & set(plus)
        if [f for f in minus if f in theta] != [f for f in plus if f in theta]: return no('(O1)')
        side = {'-': minus, '+': plus}
        for (t, e) in rays:
            pos = {f: i for i, f in enumerate(side[e])}
            if t[0] not in pos or t[1] not in pos: return no('taco missing')
            lo, hi = sorted((pos[t[0]], pos[t[1]]))
            if any(lo < pos[c] < hi for c in theta if c in pos): return no('(O3)')
        for i in range(len(rays)):
            for j in range(i + 1, len(rays)):
                (t1, e1), (t2, e2) = rays[i], rays[j]
                if e1 != e2: continue
                pos = {f: k for k, f in enumerate(side[e1])}
                a1, a2 = sorted((pos[t1[0]], pos[t1[1]])); b1, b2 = sorted((pos[t2[0]], pos[t2[1]]))
                if (a1 < b1 < a2) != (a1 < b2 < a2): return no('(O2)')
    return True


def realizable(steps, reverse=False, ar=NUM):
    """Some choice of labelling and up side (and, if reverse, of senses: reverse folds) at every
    step makes all steps rigid. steps: list of lists of directions (rays made or erased by each
    step), in the number type of ar. Returns (True, choices) or (False, None); a choice is
    (mu0, up, split) with split True when some run is split at a spine."""
    def dfs(st, i, ch):
        if i == len(steps): return list(ch)
        for mu0 in (0, 1):
            senses = list(sense_options(st, steps[i], mu0)) if reverse else [None]
            for up in (0, 1):
                for sense in senses:
                    nst = step(st, steps[i], mu0, up, sense)
                    if nst is not None:
                        r = dfs(nst, i + 1, ch + [(mu0, up, bool(sense and any(sense)))])
                        if r is not None: return r
        return None
    r = dfs(initial(ar), 0, [])
    return (r is not None, r)



# ---------------------------------------------------------------- geometry once, choices combinatorially
class Skeleton:
    """The geometry of a step D from a flat state S, for both labelings: everything the step
    needs except the layer orders. Only this part consults the arithmetic; given a skeleton, a step
    for a choice of labeling, up side and senses is a computation on orders (comb_step)."""


def skeleton(ar, S, D):
    st = State(ar, S, None)
    E = st.E(); arcsE = arcs_of(E)
    def arc_index(u):
        for i, arc in enumerate(arcsE):
            if inside(ar, u, arc): return i
        raise ValueError('no arc')
    D = [geom.nrm(ar, d) for d in D]
    R = geom.uniq(ar, list(S) + D)
    inD = [any(geom.eq(ar, R[k], d) for d in D) for k in range(len(R))]
    n = len(R); base, cur = [], 0
    for k in range(n):
        base.append(cur)
        if inD[(k + 1) % n]: cur ^= 1
    if cur != 0: return None
    sk = Skeleton(); sk.n = n
    pieces = [(R[k], R[(k + 1) % n]) for k in range(n)]
    parent = [st.face_of(R[k]) for k in range(n)]
    lam = geom.modpi(ar, st.image(D[0]))
    L2 = [lam, geom.nrm(ar, lam + ar.PI)]
    inup0 = lambda u: inside(ar, u, (L2[0], L2[1]))
    sig = lambda x: geom.nrm(ar, 2 * lam - x)
    phi_piece = [img_arc(ar, st.phis[parent[k]], *pieces[k]) for k in range(n)]
    sk.parent = parent
    sk.F = []
    for a, b in arcs_of(geom.uniq(ar, E + L2)):
        u = geom.mid(ar, a, b)
        P = [k for k in range(n) if inside(ar, u, phi_piece[k])]
        if P: sk.F.append((P, arc_index(u), inup0(u)))
    fine = arcs_of(geom.uniq(ar, E + [sig(x) for x in E] + L2))
    S2 = geom.uniq(ar, [x for x in S if not any(geom.eq(ar, x, d) for d in D)] + [d for d in D if not any(geom.eq(ar, d, x) for x in S)])
    st2 = State(ar, S2, None)
    sk.S2 = S2
    sk.G = [st2.face_of(R[k]) for k in range(n)]
    E2 = st2.E(); sk.nE2 = max(1, len(E2))
    isE = lambda x: any(geom.eq(ar, x, v) for v in E2)
    e2index = lambda x: next(i for i, v in enumerate(E2) if geom.eq(ar, x, v))
    sk.mu, sk.N, sk.compress, sk.senses = [], [], [], []
    for L in (0, 1):
        mu = [m ^ L for m in base]; sk.mu.append(mu)
        psi = [phi_piece[k] if not mu[k] else (sig(phi_piece[k][1]), sig(phi_piece[k][0])) for k in range(n)]
        rows = []
        for a, b in fine:
            u = geom.mid(ar, a, b)
            P = [k for k in range(n) if inside(ar, u, psi[k])]
            rows.append((P, arc_index(u), arc_index(sig(u)) if any(mu[k] for k in P) else None, inup0(u)))
        sk.N.append(rows)
        # g with Phi_S2 = g o Psi, and where each fine arc lands among the arcs of E2
        def psi_map(k):
            e, o = st.phis[parent[k]]
            return (e, o) if not mu[k] else (-e, geom.nrm(ar, 2 * lam - o))
        gs = []
        for k in range(n):
            ep, op = psi_map(k); e2, o2 = st2.phis[sk.G[k]]
            gs.append((e2 * ep, geom.nrm(ar, o2 - e2 * ep * op)))
        eg, og = gs[0]
        if any(e != eg or not geom.eq(ar, o, og) for e, o in gs): raise AssertionError('Phi_S2 is not g o Psi')
        g = lambda y: geom.nrm(ar, eg * y + og)
        if not E2 or fine[0][0] is None:
            sk.compress.append([(j, 0) for j in range(len(fine))])
        else:
            starts = [(g(a) if eg == 1 else g(b), j) for j, (a, b) in enumerate(fine)]
            starts.sort(key=cmp_to_key(lambda s1, s2: ar.sgn(s1[0] - s2[0])))
            first = next(i for i, (x, j) in enumerate(starts) if isE(x))
            comp, cur = [], None
            for t in range(len(starts)):
                x, j = starts[(first + t) % len(starts)]
                if isE(x): cur = e2index(x)
                comp.append((j, cur))
            sk.compress.append(comp)
        sk.senses.append(list(_senses(R, mu)))
    # (O1)-(O3) data for S2: rays grouped by image point, with their tacos and sides
    m = len(S2); groups = {}
    for idx, s_ in enumerate(S2):
        i = e2index(st2.image(s_)); jb = (idx - 1) % m
        groups.setdefault(i, []).append(((jb, idx), '-' if st2.phis[jb][0] == 1 else '+'))
    sk.V = [((i - 1) % len(E2), i, rays) for i, rays in sorted(groups.items())] if E2 else []
    return sk


def _senses(R, mu):
    n = len(R); per_run = []
    for run in runs(R, None, mu):
        per_run.append((run, [None] + [(j, o) for j in range(1, len(run)) for o in (0, 1)]))
    for choice in itertools.product(*[opts for _, opts in per_run]):
        sense = [0] * n
        for (run, _), c in zip(per_run, choice):
            if c is None: continue
            j, o = c
            for i, k in enumerate(run): sense[k] = o if i < j else 1 - o
        yield tuple(sense)


def comb_step(sk, orders, L, up, sense):
    """The orders after the step with labeling L, up side up and senses, or None (a string: the
    failing condition, when why is wanted)."""
    mu, parent, G = sk.mu[L], sk.parent, sk.G
    for P, ai, inup in sk.F:                                           # (F'), (F'')
        mov = [k for k in P if mu[k]]
        if not mov: continue
        pos = {f: i for i, f in enumerate(orders[ai])}
        sta = [k for k in P if not mu[k]]
        for p in mov:
            lifted = inup ^ up ^ sense[p]
            for q in sta:
                if lifted and not pos[parent[q]] < pos[parent[p]]: return "(F')"
                if not lifted and not pos[parent[p]] < pos[parent[q]]: return "(F')"
            for q in mov:
                if sense[q] != sense[p] and lifted and not pos[parent[q]] < pos[parent[p]]: return "(F'')"
    fine = []
    for P, ai, ais, inup in sk.N[L]:                                    # rules (i)-(iii)
        sta = [k for k in P if not mu[k]]
        if sta:
            pos = {f: i for i, f in enumerate(orders[ai])}
            sta.sort(key=lambda k: pos[parent[k]])
        below = [k for k in P if mu[k] and (inup ^ up ^ sense[k])]
        above = [k for k in P if mu[k] and not (inup ^ up ^ sense[k])]
        if below or above:
            pos = {f: i for i, f in enumerate(orders[ais])}
            below.sort(key=lambda k: -pos[parent[k]]); above.sort(key=lambda k: -pos[parent[k]])
        o = tuple(G[k] for k in below + sta + above)
        if len(set(o)) != len(o): raise AssertionError('Psi is not injective on a face')
        fine.append(o)
    out = [None] * sk.nE2                                               # extension to a layer order
    for j, i in sk.compress[L]:
        if out[i] is None: out[i] = fine[j]
        elif out[i] != fine[j]: return 'no extension'
    for minus_i, plus_i, rays in sk.V:                                  # (O1)-(O3)
        minus, plus = out[minus_i], out[plus_i]
        theta = set(minus) & set(plus)
        if [f for f in minus if f in theta] != [f for f in plus if f in theta]: return '(O1)'
        side = {'-': minus, '+': plus}
        for (t, e) in rays:
            pos = {f: i for i, f in enumerate(side[e])}
            if t[0] not in pos or t[1] not in pos: return 'taco missing'
            lo, hi = sorted((pos[t[0]], pos[t[1]]))
            if any(lo < pos[c] < hi for c in theta if c in pos): return '(O3)'
        for i in range(len(rays)):
            for j in range(i + 1, len(rays)):
                (t1, e1), (t2, e2) = rays[i], rays[j]
                if e1 != e2: continue
                pos = {f: k for k, f in enumerate(side[e1])}
                a1, a2 = sorted((pos[t1[0]], pos[t1[1]])); b1, b2 = sorted((pos[t2[0]], pos[t2[1]]))
                if (a1 < b1 < a2) != (a1 < b2 < a2): return '(O2)'
    return tuple(out)


def skeletons(steps, ar=NUM):
    """The skeletons of all steps of a history, or None if some step is not a step."""
    S, out = [], []
    for D in steps:
        sk = skeleton(ar, S, D)
        if sk is None: return None
        out.append(sk); S = sk.S2
    return out


def realizable_sk(sks, reverse=False):
    """realizable() from precomputed skeletons; the search over choices is combinatorial, and
    states (step, orders) already reached are not explored again."""
    seen = set()
    def dfs(i, orders, ch):
        if i == len(sks): return list(ch)
        if (i, orders) in seen: return None
        seen.add((i, orders))
        sk = sks[i]
        for L in (0, 1):
            for up in (0, 1):
                for sense in (sk.senses[L] if reverse else [(0,) * sk.n]):
                    r = comb_step(sk, orders, L, up, sense)
                    if isinstance(r, tuple):
                        res = dfs(i + 1, r, ch + [(L, up, any(sense))])
                        if res is not None: return res
        return None
    r = dfs(0, ((0,),), [])
    return (r is not None, r)


def realizable_fast(steps, reverse=False, ar=NUM):
    sks = skeletons(steps, ar)
    if sks is None: return (False, None)
    return realizable_sk(sks, reverse)

if __name__ == '__main__':
    u = math.pi / 4
    H = {'H1': [(0, 4), (1, 7), (3, 5)], 'H2': [(1, 5), (3, 7), (0, 2, 4, 6), (2, 6)], 'H3': [(1, 5), (3, 7), (0, 4)]}
    for name, h in H.items():
        for rev in (False, True):
            ok, ch = realizable([[k * u for k in st] for st in h], reverse=rev)
            print(name, 'with reverse folds' if rev else 'half-turns only   ', 'realizable' if ok else 'not realizable', ch or '')
