"""Exact realizability cells of history types (lookup tables, refined for rigidity).

For a history type t of degree n (src/lookup.py), the vertex is a point b = (b_1, ..., b_{n-1}),
its directions in units of pi (b_0 = 0), with 0 < b_1 < ... < b_{n-1} < 2 and Kawasaki's
condition b_1 - b_2 + ... + b_{n-1} = 1. Replaying t and deciding its realizability (src/rigid.py)
consults only signs of affine forms in b. We run that computation symbolically: an exact oracle
(exact polygon arithmetic over the rationals) decides each sign from the constraints gathered so far; when they do not
decide it, the computation is split into the feasible signs and restarted in each part. The result
is a decision tree: inner nodes are affine forms with a branch per sign, leaves say, for every vertex of the cell, whether t is valid
(by (L), flatness and its cyclic order) and realizable: 2 rigidly, 1 only with reverse folds,
0 not at all. The
leaves partition the vertices into polyhedral cells (relatively open, of every dimension).

Condition (C) is not replayed: realizability implies it (flatvertex, Lemma 3.3(3)).

Usage:
  python3 src/cells.py build TABLE.json OUT.json [--seconds=S] [--first=N] [--start=N] [--witnesses=PATH]
  python3 src/cells.py merge OUT.json PART1.json PART2.json ...   (parts of one table, in order)
Witnesses (one vertex per leaf, for tools/verify_cells.py) go to PATH, or into OUT.
"""
import itertools, json, math, sys, time
from fractions import Fraction as Fr
import geom, rigid


def scaffold_den(types):
    """4 times the lcm, over the scaffold types, of the products of their pivots: the scaffold lines
    are affine in b with denominators dividing these products (src/lookup.py, eliminate)."""
    L = 1
    for t in types:
        if 'scaffold' in t:
            prod = 1
            for j, c, d in t['scaffold']['pivots']: prod *= abs(d[j])
            L = L * prod // math.gcd(L, prod)
    return 4 * L


def branch_keys(t):
    """The branches of a scaffold type: one integer u_j in [0, |pivot j|) per scaffold line."""
    sizes = [None] * len(t['scaffold']['steps'])
    for j, c, d in t['scaffold']['pivots']: sizes[j] = abs(d[j])
    return list(itertools.product(*[range(x) for x in sizes]))


def lines(ar, t, a, u):
    """The scaffold lines of branch u at the vertex a, in ar's number type, by back substitution
    through the pivot rows (lookup.branches): {step: direction mod pi of the step's first ray}.
    Later rows substitute the unreduced lines (affine in a), so that branch labels are continuous."""
    sc = t['scaffold']; steps, piv = sc['steps'], sc['pivots']; lam = [None] * len(steps)
    for j, c, d in reversed(piv):
        x = u[j] * ar.PI
        for q, cq in enumerate(c):
            if cq: x = x - cq * a[q + 1]
        for q in range(len(steps)):
            if q != j and d[q]: x = x - d[q] * lam[q]
        lam[j] = x.div(d[j]) if isinstance(x, geom.Lin) else x / d[j]
    return {st: geom.modpi(ar, lam[j]) for j, st in enumerate(steps)}


def _norm(P, L):
    if L < 0: P, L = [-x for x in P], -L
    g = math.gcd(L, *P)
    return (tuple(x // g for x in P), L // g)


class PolyOracle:
    """Signs of affine forms on a relatively open convex region of dimension <= 2, kept as the
    vertices of its closure (exact rational points (P, L): b_i = P_i / L), in cyclic order.
    f > 0 on the region iff f >= 0 at every vertex and f is not 0 at all of them; f = 0 iff it is 0
    at every vertex (the vertices span the region's affine hull); otherwise f takes both signs and
    the region splits into three nonempty parts. Splitting is exact clipping."""
    def __init__(self, m, vertices):
        self.m = m
        self.V = [list(vertices)]
        self.cache, self.memo, self.added = {}, {}, [[]]
        self.pts = [None]; self.checks = 0
    canon = staticmethod(lambda x: Oracle_canon(x))
    @staticmethod
    def val(x, v):
        P, L = v
        return sum(c * p for c, p in zip(x.c, P)) + x.k * L
    def sign(self, x):
        key, f = Oracle_canon(x)
        if key is None: return (x.k > 0) - (x.k < 0)
        if key in self.cache: return f * self.cache[key]
        xn = geom.Lin(key[0], key[1]); self.checks += 1
        vals = [self.val(xn, v) for v in self.V[-1]]
        pos, neg = any(v > 0 for v in vals), any(v < 0 for v in vals)
        if pos and neg: raise geom.Need(xn, [-1, 0, 1])
        sg = 1 if pos else (-1 if neg else 0)
        self.remember(key, sg); return f * sg
    def point(self):
        """The centroid of the vertices: a point of the region's relative interior."""
        if self.pts[-1] is None:
            V = self.V[-1]; S = [Fr(0)] * self.m
            for P, L in V:
                for i in range(self.m): S[i] += Fr(P[i], L)
            S = [x / len(V) for x in S]
            L = 1
            for x in S: L = L * x.denominator // math.gcd(L, x.denominator)
            self.pts[-1] = ([int(x * L) for x in S], L)
        return self.pts[-1]
    def remember(self, key, s):
        self.cache[key] = s; self.added[-1].append(('s', key))
    def recall(self, key): return self.memo.get(key)
    def store(self, key, value):
        if key not in self.memo: self.memo[key] = value; self.added[-1].append(('m', key))
    def push(self, x, s):
        self.V.append(clip(self.V[-1], x, s)); self.pts.append(None); self.added.append([])
        key, f = Oracle_canon(x)
        self.remember(key, f * s)
    def pop(self):
        self.V.pop(); self.pts.pop()
        for kind, key in self.added.pop():
            del (self.cache if kind == 's' else self.memo)[key]


def Oracle_canon(x):
    """(key, f): x scaled to coprime integers with a positive leading coefficient, and the sign f
    of the scale, so that sign(x) = f * sign(key)."""
    lead = next((c for c in x.c if c), None)
    if lead is None: return None, 1
    g = math.gcd(*x.c, x.k)
    f = 1 if lead > 0 else -1
    return (tuple(f * c // g for c in x.c), f * x.k // g), f


def _cross(A, B, fa, fb):
    """The point of the segment AB where an affine form with values fa, fb (numerators) vanishes."""
    (PA, LA), (PB, LB) = A, B
    return _norm([fa * pb - fb * pa for pa, pb in zip(PA, PB)], fa * LB - fb * LA)


def clip(V, x, s):
    """The closure of (region with vertices V) intersected with x < 0, x = 0 or x > 0 (s = -1, 0, 1)."""
    vals = [PolyOracle.val(x, v) for v in V]
    out = []
    n = len(V)
    for i in range(n):
        A, B = V[i], V[(i + 1) % n]; fa, fb = vals[i], vals[(i + 1) % n]
        if s == 0:
            if fa == 0: out.append(A)
        elif fa * s >= 0: out.append(A)
        if (fa < 0 < fb or fb < 0 < fa) and n > 1: out.append(_cross(A, B, fa, fb))
    res = []
    for v in out:
        if v not in res: res.append(v)
    return res


def pieces(t):
    """The pieces of a type: for each assignment of integers to a basis of its required congruences
    (with Kawasaki's condition b_1 - b_2 + ... = 1), the closure of the part of that affine
    subspace where 0 < b_1 < ... < b_{n-1} < 2, if it is nonempty with the subspace's dimension.
    Returns (basis rows, [(values, vertices)])."""
    import itertools
    n = t['n']; m = n - 1
    kaw = tuple((-1) ** i for i in range(m))                     # b_1 - b_2 + b_3 - ...
    rows, basis = [kaw] + [tuple(v) for v in t['req']], []
    def rank(R):
        M = [[Fr(x) for x in r] for r in R]; rk = 0
        for c in range(m):
            piv = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
            if piv is None: continue
            M[rk], M[piv] = M[piv], M[rk]
            for i in range(len(M)):
                if i != rk and M[i][c] != 0:
                    f = M[i][c] / M[rk][c]; M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
            rk += 1
        return rk
    for r in rows:
        if rank(basis + [r]) > len(basis): basis.append(r)
    others = [r for r in rows if r not in basis]
    ranges = [range(1, 2)] + [range(sum(2 * c for c in r if c < 0), sum(2 * c for c in r if c > 0) + 1) for r in basis[1:]]
    order = [geom.Lin(tuple(geom.DEN if j == 0 else 0 for j in range(m)))] + \
            [geom.Lin(tuple((geom.DEN if j == i + 1 else (-geom.DEN if j == i else 0)) for j in range(m))) for i in range(m - 1)] + \
            [geom.Lin(tuple(-geom.DEN if j == m - 1 else 0 for j in range(m)), 2 * geom.DEN)]   # all > 0
    out = []
    for vals in itertools.product(*ranges):
        # solve basis . b = vals: reduced row echelon form over the rationals
        M = [[Fr(x) for x in r] + [Fr(v)] for r, v in zip(basis, vals)]; piv = []; rk = 0
        for c in range(m):
            p_ = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
            if p_ is None: continue
            M[rk], M[p_] = M[p_], M[rk]; M[rk] = [x / M[rk][c] for x in M[rk]]
            for i in range(len(M)):
                if i != rk and M[i][c] != 0:
                    f = M[i][c]; M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
            piv.append(c); rk += 1
        free = [c for c in range(m) if c not in piv]
        def point(tv):
            b = [Fr(0)] * m
            for j, c in enumerate(free): b[c] = Fr(tv[j])
            for i, c in enumerate(piv): b[c] = M[i][m] - sum(M[i][f_] * b[f_] for f_ in free)
            L = 1
            for x in b: L = L * x.denominator // math.gcd(L, x.denominator)
            return _norm([int(x * L) for x in b], L), b
        d = len(free)
        corners = {0: [()], 1: [(0,), (2,)], 2: [(0, 0), (2, 0), (2, 2), (0, 2)]}[d]
        V = [point(c)[0] for c in corners]
        # the other required congruences must hold (integers) on the subspace
        b0 = point((0,) * d)[1]
        if any((sum(Fr(c) * x for c, x in zip(r, b0))).denominator != 1 for r in others): continue
        for f in order:
            V = clip(V, f, 1) if len(V) > 1 else V
            if not V: break
        if not V: continue
        if any(all(PolyOracle.val(f, v) == 0 for v in V) for f in order): continue   # a strict order fails
        if d == 0 and any(PolyOracle.val(f, V[0]) <= 0 for f in order): continue
        dim = 0 if len(V) == 1 else (1 if len(V) == 2 or _collinear(V) else 2)
        if dim != d: continue
        out.append((list(vals[1:]), V))
    return [list(r) for r in basis[1:]], out


def _collinear(V):
    (P0, L0), (P1, L1) = V[0], V[1]
    u = [Fr(a, L1) - Fr(b, L0) for a, b in zip(P1, P0)]
    for P, L in V[2:]:
        w = [Fr(a, L) - Fr(b, L0) for a, b in zip(P, P0)]
        for i in range(len(u)):
            for j in range(i + 1, len(u)):
                if u[i] * w[j] - u[j] * w[i] != 0: return False
    return True


def replay(ar, t, a, lams=None):
    """lookup.replay in any arithmetic, without (C): the hidden directions, or None."""
    rays, N = t['rays'], len(t['rays'])
    d = [a[r[1]] if r[0] == 'v' else None for r in rays]
    status = [0] * N
    for i, T in enumerate(t['steps']):
        S = [j for j in range(N) if status[j] == 1]
        Sd = geom.sort(ar, [d[j] for j in S])
        phis = geom.maps(ar, Sd)
        if phis is None: return None
        anchor = next((j for j in T if rays[j][0] == 'v' or status[j] == 1), None)
        if anchor is None:                                   # a scaffold step: its first ray from its line
            anchor = T[0]; x = lams[i]
            if S:
                p = max((q for q in S if q < anchor), default=max(S)); q = min((q for q in S if q > anchor), default=min(S))
                c = [y for y in (x, x + ar.PI) if geom.between(ar, y, d[p], d[q])]
                if len(c) != 1: return None
                x = c[0]
            d[anchor] = x
        lam = geom.apply(ar, phis[geom.sector_of(ar, d[anchor], Sd)], d[anchor]) if S else d[anchor]
        for j in T:
            if j != anchor and status[j] == 0 and rays[j][0] == 'h':
                if not S: d[j] = geom.nrm(ar, lam + ar.PI); continue
                p = max((i for i in S if i < j), default=max(S)); q = min((i for i in S if i > j), default=min(S))
                f = phis[geom.sector_of(ar, d[p], Sd)]
                c = [x for x in (geom.inverse_apply(ar, f, lam), geom.inverse_apply(ar, f, geom.nrm(ar, lam + ar.PI)))
                     if geom.between(ar, x, d[p], d[q])]
                if len(c) != 1: return None
                d[j] = c[0]
        imgs = [geom.modpi(ar, geom.apply(ar, phis[geom.sector_of(ar, d[j], Sd)], d[j])) if S else geom.modpi(ar, d[j]) for j in T]
        if any(not geom.eq(ar, geom.modpi(ar, x - imgs[0]), ar.ZERO) and not geom.eq(ar, geom.modpi(ar, x - imgs[0]), ar.PI) for x in imgs): return None
        for j in T: status[j] = 1 if status[j] == 0 else 2
        if geom.maps(ar, geom.sort(ar, [d[j] for j in range(N) if status[j] == 1])) is None: return None
    for j in range(N - 1):                                         # the cyclic order of the type
        if ar.sgn(d[j + 1] - d[j]) <= 0: return None
    if ar.sgn(ar.TWO_PI - d[-1]) <= 0: return None
    return d


def decide(ar, t, a, tid=None, u=None):
    """At the vertex a (directions in ar's number type): 2 if t is valid and rigidly realizable,
    1 if valid and realizable with reverse folds only, 0 otherwise; for a scaffold type, with the
    scaffold lines of branch u. All arithmetic happens in the replay and the skeletons; the search
    over choices is combinatorial."""
    O = getattr(ar, 'oracle', None)
    key = ('replay', tid, u)
    if O is not None and O.recall(key) is not None: d = O.recall(key)
    else:
        d = replay(ar, t, a, lines(ar, t, a, u) if u is not None else None)
        if O is not None: O.store(key, d if d is not None else False)
    if d is None or d is False: return 0
    sks = rigid.skeletons([[d[j] for j in T] for T in t['steps']], ar)
    if sks is None: return 0
    if rigid.realizable_sk(sks, False)[0]: return 2
    return 1 if rigid.realizable_sk(sks, True)[0] else 0


def build_tree(t, deadline):
    n = t['n']; m = n - 1
    basis, ps = pieces(t)
    stats = {'leaves': 0, 'nodes': 0, 'checks': 0, 'pieces': len(ps)}; witnesses = []
    tree = {'basis': basis, 'pieces': {}}
    for (vals, V), u in itertools.product(ps, branch_keys(t) if 'scaffold' in t else [None]):
        O = PolyOracle(m, V); ar = geom.Sym(m, O)
        a = [ar.ZERO] + [ar.var(i) for i in range(m)]
        def explore():
            if time.monotonic() > deadline: raise TimeoutError
            try:
                r = decide(ar, t, a, 0, u); stats['leaves'] += 1
                P, L = O.point()                                         # one vertex in the cell
                witnesses.append([[p / L for p in P], r] + ([list(u)] if u is not None else []))
                return r
            except geom.Need as e:
                stats['nodes'] += 1; ch = {}
                for s_ in e.signs:
                    O.push(e.form, s_); ch[s_] = explore(); O.pop()
                vals_ = list(ch.values())
                if all(isinstance(v, int) for v in vals_) and len(set(vals_)) == 1:
                    return vals_[0]                                  # all branches agree: no split needed
                return {'f': [str(c) for c in e.form.c], 'k': str(e.form.k), 'c': {str(s_): v for s_, v in ch.items()}}
        tree['pieces'][piece_key(vals, u)] = explore()
        stats['checks'] += O.checks
    return tree, stats, witnesses


def piece_key(vals, u=None):
    return ','.join(map(str, vals)) + ('' if u is None else '|' + ','.join(map(str, u)))


def evaluate(tree, a, tol=1e-7, den=4, u=None):
    """The value of tree at the vertex with directions a (radians, a[0] = 0), for a scaffold type
    with the scaffold lines of branch u: its piece is given by the values of the basis congruences;
    0 if it lies in no piece. den: the denominator of the forms (the cells file's meta)."""
    b = [x / math.pi for x in a[1:]]
    vals = []
    for r in tree['basis']:
        v = sum(c * x for c, x in zip(r, b))
        if abs(v - round(v)) > 1e-6: return 0
        vals.append(int(round(v)))
    node = tree['pieces'].get(piece_key(vals, u))
    if node is None: return 0
    while isinstance(node, dict):
        v = (sum(int(c) * x for c, x in zip(node['f'], b)) + int(node['k'])) / den
        sg = 0 if abs(v) < tol else (1 if v > 0 else -1)
        node = node['c'].get(str(sg))
        if node is None: return None                                # outside the explored region
    return node


def main(argv):
    if argv[0] == 'merge':                                  # merge OUT1 OUT2 ... into OUT (types in order)
        out, parts = argv[1], argv[2:]
        D = [json.load(open(x)) for x in parts]
        meta = dict(D[0]['meta'])                           # the first part was built from the whole table
        meta['trees'] = sum(len(x['trees']) for x in D)
        for k in ('seconds', 'leaves', 'tree_nodes', 'sign_queries'): meta[k] = round(sum(x['meta'][k] for x in D), 1)
        meta['max_type_seconds'] = max(x['meta']['max_type_seconds'] for x in D)
        meta['complete'] = meta['trees'] == meta['types'] and D[-1]['meta']['complete']
        json.dump({'meta': meta, 'trees': [t for x in D for t in x['trees']], 'stats': [r for x in D for r in x['stats']]}, open(out, 'w'))
        print(json.dumps(meta)); return
    if argv[0] == 'build':
        src, out = argv[1], argv[2]
        secs = float(next((x[10:] for x in argv if x.startswith('--seconds=')), 3600))
        d = json.load(open(src)); t0 = time.monotonic(); deadline = t0 + secs
        first = int(next((x[8:] for x in argv if x.startswith('--first=')), 0))
        start = int(next((x[8:] for x in argv if x.startswith('--start=')), 0))
        geom.set_den(scaffold_den(d['types']))         # over the whole table, so that parts agree
        if first: d['types'] = d['types'][:first]
        if start: d['types'] = d['types'][start:]
        trees, rows, wits, complete = [], [], [], True
        for i, t in enumerate(d['types']):
            t1 = time.monotonic()
            try: tree, st, w = build_tree(t, deadline)
            except TimeoutError: complete = False; break
            trees.append(tree); wits.append(w); st['seconds'] = round(time.monotonic() - t1, 3); rows.append(st)
        meta = {'n': d['meta']['n'], 'k': d['meta']['k'], 'types': len(d['types']), 'trees': len(trees),
                'complete': complete and len(trees) == len(d['types']), 'seconds': round(time.monotonic() - t0, 1),
                'leaves': sum(r['leaves'] for r in rows), 'tree_nodes': sum(r['nodes'] for r in rows),
                'sign_queries': sum(r['checks'] for r in rows), 'max_type_seconds': max((r['seconds'] for r in rows), default=0),
                'den': geom.DEN}
        wpath = next((x[12:] for x in argv if x.startswith('--witnesses=')), None)
        json.dump({'meta': meta, 'trees': trees, 'stats': rows} if wpath else {'meta': meta, 'trees': trees, 'stats': rows, 'witnesses': wits}, open(out, 'w'))
        if wpath: json.dump({'witnesses': wits}, open(wpath, 'w'))
        print(json.dumps(meta), flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
