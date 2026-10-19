"""Exact realizability cells of history types (lookup tables, refined for rigidity).

For a history type t of degree n (src/lookup.py), the vertex is a point b = (b_1, ..., b_{n-1}),
its directions in units of pi (b_0 = 0), with 0 < b_1 < ... < b_{n-1} < 2 and Kawasaki's
condition b_1 - b_2 + ... + b_{n-1} = 1. Replaying t and deciding its realizability (src/rigid.py)
consults only signs of affine forms in b. We run that computation symbolically: an exact oracle
(z3 over the rationals) decides each sign from the constraints gathered so far; when they do not
decide it, the computation is split into the feasible signs and restarted in each part. The result
is a decision tree: inner nodes are affine forms with a branch per sign, leaves say, for every vertex of the cell, whether t is valid
(by (L), flatness and its cyclic order) and realizable: 2 rigidly, 1 only with reverse folds,
0 not at all. The
leaves partition the vertices into polyhedral cells (relatively open, of every dimension).

Condition (C) is not replayed: realizability implies it (flatvertex, Lemma 3.3(3)).

Usage:
  python3 src/cells.py build TABLE.json OUT.json [--seconds=S]
"""
import json, math, sys, time
from fractions import Fraction as Fr
import z3
import geom, rigid


class Oracle:
    def __init__(self, m):
        self.m = m
        self.B = [z3.Real('b%d' % (i + 1)) for i in range(m)]
        self.s = z3.Solver()
        prev = z3.RealVal(0)
        for x in self.B: self.s.add(x > prev); prev = x
        self.s.add(prev < 2)
        self.s.add(z3.Sum([(1 if i % 2 == 0 else -1) * x for i, x in enumerate(self.B)]) == 1)   # Kawasaki
        self.cache = [{}]; self.model = [None]; self.checks = 0
    def z(self, x):
        D = geom.DEN
        return z3.Sum([z3.RealVal(0)] + [z3.Q(c, D) * self.B[i] for i, c in enumerate(x.c) if c]) + z3.Q(x.k, D)
    @staticmethod
    def canon(x):
        """(key, f): key is x scaled to coprime integers with a positive leading coefficient, f the
        sign of the scale, so that sign(x) = f * sign(key)."""
        lead = next((c for c in x.c if c), None)
        if lead is None: return None, 1
        g = math.gcd(*x.c, x.k)
        f = 1 if lead > 0 else -1
        return (tuple(f * c // g for c in x.c), f * x.k // g), f
    def rel(self, x, s):
        e = self.z(x)
        return e < 0 if s < 0 else (e == 0 if s == 0 else e > 0)
    def feasible(self, x, s):
        self.s.push(); self.s.add(self.rel(x, s)); self.checks += 1
        r = self.s.check(); self.s.pop()
        if r == z3.unknown: raise RuntimeError('z3 unknown')
        return r == z3.sat
    def sign(self, x):
        key, f = self.canon(x)
        if key is None: return (x.k > 0) - (x.k < 0)
        for d in reversed(self.cache):
            if key in d: return f * d[key]
        xn = geom.Lin(key[0], key[1])
        signs = [s for s in (-1, 0, 1) if self.feasible(xn, s)]
        if len(signs) == 1: self.cache[-1][key] = signs[0]; return f * signs[0]
        raise geom.Need(xn, signs)
    def point(self):
        """One point of the current region (a z3 model), as integers P over a common denominator L."""
        if self.model[-1] is None:
            self.checks += 1; assert self.s.check() == z3.sat
            M = self.s.model(); pt = []
            for b in self.B:
                v = M.eval(b, model_completion=True); pt.append(Fr(v.numerator_as_long(), v.denominator_as_long()))
            L = 1
            for v in pt: L = L * v.denominator // math.gcd(L, v.denominator)
            self.model[-1] = ([int(v * L) for v in pt], L)
        return self.model[-1]
    def push(self, x, s):
        self.s.push(); self.s.add(self.rel(x, s))
        key, _ = self.canon(x)
        self.cache.append({key: s}); self.model.append(None)
    def pop(self):
        self.s.pop(); self.cache.pop(); self.model.pop()


def replay(ar, t, a):
    """lookup.replay in any arithmetic, without (C): the hidden directions, or None."""
    rays, N = t['rays'], len(t['rays'])
    d = [a[r[1]] if r[0] == 'v' else None for r in rays]
    status = [0] * N
    for T in t['steps']:
        S = [j for j in range(N) if status[j] == 1]
        Sd = geom.sort(ar, [d[j] for j in S])
        phis = geom.maps(ar, Sd)
        if phis is None: return None
        anchor = next(j for j in T if rays[j][0] == 'v' or status[j] == 1)
        lam = geom.apply(ar, phis[geom.sector_of(ar, d[anchor], Sd)], d[anchor]) if S else d[anchor]
        for j in T:
            if status[j] == 0 and rays[j][0] == 'h':
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


def decide(ar, t, a):
    """At the vertex a (directions in ar's number type): 2 if t is valid and rigidly realizable,
    1 if valid and realizable with reverse folds only, 0 otherwise. All arithmetic happens in the
    replay and the skeletons; the search over choices is combinatorial."""
    d = replay(ar, t, a)
    if d is None: return 0
    sks = rigid.skeletons([[d[j] for j in T] for T in t['steps']], ar)
    if sks is None: return 0
    if rigid.realizable_sk(sks, False)[0]: return 2
    return 1 if rigid.realizable_sk(sks, True)[0] else 0


def build_tree(t, deadline):
    n = t['n']; m = n - 1
    O = Oracle(m); ar = geom.Sym(m, O)
    a = [ar.ZERO] + [ar.var(i) for i in range(m)]
    stats = {'leaves': 0, 'nodes': 0}
    def explore():
        if time.monotonic() > deadline: raise TimeoutError
        try:
            r = decide(ar, t, a); stats['leaves'] += 1
            return r
        except geom.Need as e:
            stats['nodes'] += 1; ch = {}
            for s in e.signs:
                O.push(e.form, s); ch[s] = explore(); O.pop()
            vals = list(ch.values())
            if all(isinstance(v, int) for v in vals) and len(set(vals)) == 1:
                return vals[0]                                      # all branches agree: no split needed
            return {'f': [str(c) for c in e.form.c], 'k': str(e.form.k), 'c': {str(s): v for s, v in ch.items()}}
    tree = explore()
    stats['checks'] = O.checks
    return tree, stats


def evaluate(tree, a, tol=1e-7):
    """The leaf of tree at the vertex with directions a (radians, a[0] = 0)."""
    b = [x / math.pi for x in a[1:]]
    while isinstance(tree, dict):
        v = (sum(int(c) * x for c, x in zip(tree['f'], b)) + int(tree['k'])) / geom.DEN
        s = 0 if abs(v) < tol else (1 if v > 0 else -1)
        tree = tree['c'].get(str(s))
        if tree is None: return None                               # outside the explored region
    return tree


def main(argv):
    if argv[0] == 'build':
        src, out = argv[1], argv[2]
        secs = float(next((x[10:] for x in argv if x.startswith('--seconds=')), 3600))
        d = json.load(open(src)); t0 = time.monotonic(); deadline = t0 + secs
        trees, rows, complete = [], [], True
        for i, t in enumerate(d['types']):
            t1 = time.monotonic()
            try: tree, st = build_tree(t, deadline)
            except TimeoutError: complete = False; break
            trees.append(tree); st['seconds'] = round(time.monotonic() - t1, 3); rows.append(st)
        meta = {'n': d['meta']['n'], 'k': d['meta']['k'], 'types': len(d['types']), 'trees': len(trees),
                'complete': complete and len(trees) == len(d['types']), 'seconds': round(time.monotonic() - t0, 1),
                'leaves': sum(r['leaves'] for r in rows), 'tree_nodes': sum(r['nodes'] for r in rows),
                'z3_checks': sum(r['checks'] for r in rows), 'max_type_seconds': max((r['seconds'] for r in rows), default=0)}
        json.dump({'meta': meta, 'trees': trees, 'stats': rows}, open(out, 'w'))
        print(json.dumps(meta), flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
