"""Rigid realizability of a local history (flatvertex.tex, Definitions 2.1 and 2.2).

A layered state is a flat state S with, at every point u of the image circle outside the image
points E_S, a strict total order (bottom to top) on the faces whose image contains u, locally
constant in u, satisfying (O1)-(O3). A rigid step folds a step D along its line l by a half-turn:
it is given by a labelling (which pieces move) and a choice of the up component of the circle
minus the two rays of l. (F) requires the moving layers to lie above the staying ones over the up
side and below them over the down side. The new order is given by rules (i)-(iii) and must
extend to a layered state. A history is rigidly realizable if some choices make every step rigid.

Orders are stored on the arcs between consecutive image points. Directions are floats, compared
with fold_history.eq / between (tolerance fold_history.EPS).

Usage: python3 src/rigid.py   (checks the three histories of flatvertex.tex, Section 3)
"""
import math, sys
import fold_history as fh

PI, TWO_PI = math.pi, 2 * math.pi


def uniq(xs):
    out = []
    for x in sorted(fh.nrm(x) for x in xs):
        if not out or not fh.eq(x, out[-1]): out.append(x)
    if len(out) > 1 and fh.eq(out[0], out[-1]): out.pop()
    return out


def arcs_of(B):
    """Consecutive arcs (a, b) of the circle cut at the sorted points B; one full arc if B is empty."""
    if not B: return [(None, None)]
    return [(B[i], B[(i + 1) % len(B)]) for i in range(len(B))]


def mid(a, b):
    if a is None: return 0.0
    return fh.nrm(a + ((b - a) % TWO_PI or TWO_PI) / 2)


def inside(u, arc):
    a, b = arc
    return True if a is None else fh.between(u, a, b)


class State:
    """S: sorted directions; phis: folding maps; arcs: [(a, b, order)] covering the circle."""
    def __init__(self, S, arcs):
        self.S = S
        self.phis = fh.maps(S) if S else [(1, 0.0)]
        self.arcs = arcs
    def face_of(self, x):
        return fh.sector_of(x, self.S) if self.S else 0
    def image(self, x):
        return fh.apply(self.phis[self.face_of(x)], x)
    def E(self):
        return uniq(self.image(s) for s in self.S)
    def order_at(self, u):
        for a, b, o in self.arcs:
            if inside(u, (a, b)): return o
        raise ValueError('no arc contains %r' % u)


def initial():
    return State([], [(None, None, (0,))])


def img_arc(f, a, b):
    """Image of the open arc (a, b) under x -> e x + o, as a counterclockwise arc."""
    e, o = f
    if a is None: return (None, None)
    x, y = fh.apply(f, a), fh.apply(f, b)
    return (x, y) if e == 1 else (y, x)


def step(st, D, mu0, up):
    """The layered state after the rigid step D with labelling mu0 and up side up, or None."""
    S = st.S
    D = [fh.nrm(d) for d in D]
    R = uniq(list(S) + D)
    inD = lambda x: any(fh.eq(x, d) for d in D)
    n = len(R)
    pieces = [(R[k], R[(k + 1) % n]) for k in range(n)]
    parent = [st.face_of(R[k]) for k in range(n)]
    mu, cur = [], mu0
    for k in range(n):
        mu.append(cur)
        if inD(R[(k + 1) % n]): cur ^= 1
    if cur != mu0: return None                                   # |D| odd: not a step
    lam = st.image(D[0]) % PI
    L2 = [lam, fh.nrm(lam + PI)]
    uparc = (L2[0], L2[1]) if up == 0 else (L2[1], L2[0])
    sig = lambda x: fh.nrm(2 * lam - x)
    phi_piece = [img_arc(st.phis[parent[k]], *pieces[k]) for k in range(n)]
    def psi_piece(k):
        a, b = phi_piece[k]
        return (a, b) if not mu[k] else (sig(b), sig(a))
    psi = [psi_piece(k) for k in range(n)]
    ES = st.E()
    # (F)
    for a, b in arcs_of(uniq(ES + L2)):
        u = mid(a, b)
        P = [k for k in range(n) if inside(u, phi_piece[k])]
        mov = [k for k in P if mu[k]]; sta = [k for k in P if not mu[k]]
        if not (mov and sta): continue
        pos = {f: i for i, f in enumerate(st.order_at(u))}
        on_up = inside(u, uparc)
        for p in mov:
            for q in sta:
                if on_up and not pos[parent[q]] < pos[parent[p]]: return None
                if not on_up and not pos[parent[p]] < pos[parent[q]]: return None
    # new order, in the coordinates of Psi
    S2 = uniq([x for x in S if not inD(x)] + [d for d in D if not any(fh.eq(d, x) for x in S)])
    new = State(S2, None)
    G = [new.face_of(R[k]) for k in range(n)]
    fine = []
    for a, b in arcs_of(uniq(ES + [sig(x) for x in ES] + L2)):
        u = mid(a, b)
        P = [k for k in range(n) if inside(u, psi[k])]
        sta = [k for k in P if not mu[k]]; mov = [k for k in P if mu[k]]
        if sta:
            pos = {f: i for i, f in enumerate(st.order_at(u))}
            sta.sort(key=lambda k: pos[parent[k]])                       # rule (i)
        if mov:
            pos = {f: i for i, f in enumerate(st.order_at(sig(u)))}
            mov.sort(key=lambda k: -pos[parent[k]])                      # rule (ii): reversed
        seq = (mov + sta) if inside(u, uparc) else (sta + mov)           # rule (iii)
        order = tuple(G[k] for k in seq)
        if len(set(order)) != len(order): raise AssertionError('Psi is not injective on a face')
        fine.append((a, b, order))
    # g with Phi_S2 = g o Psi (checked on every piece)
    def psi_map(k):
        e, o = st.phis[parent[k]]
        return (e, o) if not mu[k] else (-e, fh.nrm(2 * lam - o))
    gs = []
    for k in range(n):
        ep, op = psi_map(k); e2, o2 = new.phis[G[k]]
        gs.append((e2 * ep, fh.nrm(o2 - e2 * ep * op)))
    eg, og = gs[0]
    if any(e != eg or not fh.eq(o, og) for e, o in gs): raise AssertionError('Phi_S2 is not g o Psi')
    g = lambda y: fh.nrm(eg * y + og)
    moved = []
    for a, b, o in fine:
        if a is None: moved.append((None, None, o)); continue
        moved.append((g(a), g(b), o) if eg == 1 else (g(b), g(a), o))
    # extend to a layer order: constant across breakpoints that are not image points of S2
    E2 = new.E()
    isE = lambda x: any(fh.eq(x, v) for v in E2)
    if moved[0][0] is not None:
        moved.sort(key=lambda t: t[0])
        m = len(moved)
        start = next((i for i in range(m) if isE(moved[i][0])), None)
        if start is None:                                                # no image points: one order
            if any(o != moved[0][2] for _, _, o in moved): return None
            arcs = [(None, None, moved[0][2])]
        else:
            arcs, i = [], start
            for _ in range(m):
                a, b, o = moved[i]
                if isE(a): arcs.append([a, b, o])
                else:
                    if o != arcs[-1][2]: return None
                    arcs[-1][1] = b
                i = (i + 1) % m
            arcs = [tuple(x) for x in arcs]
    else:
        arcs = moved
    new.arcs = arcs
    return new if layered_ok(new) else None


def layered_ok(st):
    """(O1)-(O3) for a state whose orders are stored on the arcs between its image points."""
    S = st.S
    if not S: return True
    m = len(S)
    by_v = {}
    for idx, s in enumerate(S):
        v = st.image(s)
        key = next((w for w in by_v if fh.eq(w, v)), v)
        jb = (idx - 1) % m                                               # face before s
        eps = '-' if st.phis[jb][0] == 1 else '+'
        by_v.setdefault(key, []).append(((jb, idx), eps))
    for v, rays in by_v.items():
        minus = next(o for a, b, o in st.arcs if a is not None and fh.eq(b, v))
        plus = next(o for a, b, o in st.arcs if a is not None and fh.eq(a, v))
        theta = set(minus) & set(plus)
        if [f for f in minus if f in theta] != [f for f in plus if f in theta]: return False        # (O1)
        side = {'-': minus, '+': plus}
        for (t, e) in rays:
            pos = {f: i for i, f in enumerate(side[e])}
            if t[0] not in pos or t[1] not in pos: return False
            lo, hi = sorted((pos[t[0]], pos[t[1]]))
            if any(lo < pos[c] < hi for c in theta if c in pos): return False                    # (O3)
        for i in range(len(rays)):
            for j in range(i + 1, len(rays)):
                (t1, e1), (t2, e2) = rays[i], rays[j]
                if e1 != e2: continue
                pos = {f: k for k, f in enumerate(side[e1])}
                a1, a2 = sorted((pos[t1[0]], pos[t1[1]])); b1, b2 = sorted((pos[t2[0]], pos[t2[1]]))
                if (a1 < b1 < a2) != (a1 < b2 < a2): return False                                # (O2)
    return True


def realizable(steps, limit=None):
    """Some choice of labelling and up side at every step makes all steps rigid.
    steps: list of lists of directions (rays made or erased by each step).
    Returns (True, choices) or (False, None)."""
    def dfs(st, i, ch):
        if i == len(steps): return list(ch)
        for mu0 in (0, 1):
            for up in (0, 1):
                nst = step(st, steps[i], mu0, up)
                if nst is not None:
                    r = dfs(nst, i + 1, ch + [(mu0, up)])
                    if r is not None: return r
        return None
    r = dfs(initial(), 0, [])
    return (r is not None, r)


if __name__ == '__main__':
    u = PI / 4
    H = {'H1': [(0, 4), (1, 7), (3, 5)], 'H2': [(1, 5), (3, 7), (0, 2, 4, 6), (2, 6)], 'H3': [(1, 5), (3, 7), (0, 4)]}
    for name, h in H.items():
        ok, ch = realizable([[k * u for k in st] for st in h])
        print(name, 'rigidly realizable' if ok else 'not rigidly realizable', ch or '')
