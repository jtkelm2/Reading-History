"""Lookup tables of local histories (first build).

Prepopulation.  For a degree n and a number k of hidden rays, every *history type* valid at some
flat vertex of degree n is enumerated: a cyclic order of the n visible rays and k hidden rays, and
a sequence of single-line steps.  The visible directions are symbolic (a_0 = 0, a_1 < ... < a_{n-1}),
so the depth-first search over steps keeps a prefix exactly when z3 finds some flat vertex and
some hidden directions satisfying it: (P), (L) by Theorem 11, flatness of every state, (C), and the
cyclic order.  Steps whose rays are all new hidden rays are excluded, so every hidden direction is
an integer combination of visible directions modulo pi (no free hidden directions).  By rotation
symmetry the first step contains visible ray 0; the runtime tries every rotation.

Each type records the congruences c.a = 0 (mod pi) it needs, with hidden directions substituted:
the equation part of its cell.  The table indexes types by these congruence vectors.

Runtime.  At a vertex, for each rotation of its rays, evaluate the congruence vectors, take the
types whose required congruences hold, and verify each by replaying its steps at the vertex (the
inequalities and (C)).  No step search happens at runtime.  Results are cached per vertex shape.

Usage:
  python3 src/lookup.py build N K OUTDIR [--seconds=S]     prepopulate one (degree, hidden rays)
  python3 src/lookup.py stats OUTDIR                        summarise the tables
"""
import itertools, json, math, os, sys, time
import z3
import fold_history as fh

PI = math.pi
TWO_PI = 2 * math.pi
ZPI = z3.RealVal(repr(PI))


class Lin:
    """sum(coef * var) + const, vars hashable keys."""
    __slots__ = ('coef', 'const')
    def __init__(self, coef=None, const=0.0): self.coef, self.const = dict(coef or {}), const
    def __add__(self, o):
        o = o if isinstance(o, Lin) else Lin(const=o); c = dict(self.coef)
        for i, a in o.coef.items(): c[i] = c.get(i, 0) + a
        return Lin({i: a for i, a in c.items() if a}, self.const + o.const)
    def __rmul__(self, k): return Lin({i: k * a for i, a in self.coef.items() if k * a}, k * self.const)
    def __sub__(self, o): return self + (-1) * (o if isinstance(o, Lin) else Lin(const=o))
    def z3(self, V): return z3.Sum([z3.RealVal(repr(self.const))] + [a * V[i] for i, a in self.coef.items()])


def vec(lin, n):
    """The congruence vector of a mod-pi linear form in the visible directions a_1..a_{n-1},
    sign-normalised; None if it is identically zero."""
    v = [lin.coef.get(('a', i), 0) for i in range(1, n)]
    if not any(v): return None
    if next(x for x in v if x) < 0: v = [-x for x in v]
    return tuple(v)


# ---------------------------------------------------------------- prepopulation
class Builder:
    def __init__(self, n, k, deadline, check_ms=20000, scaffold=False):
        """scaffold: enumerate instead the types with at least one scaffold step (a step making only
        hidden rays), keeping those whose scaffold lines are all pinned (see eliminate)."""
        self.n, self.k, self.deadline, self.check_ms, self.scaffold = n, k, deadline, check_ms, scaffold
        self.types, self.nodes, self.checks, self.unknown, self.complete = [], 0, 0, 0, True
        self.unpinned = 0

    def placements(self):
        for comp in itertools.product(range(self.k + 1), repeat=self.n):
            if sum(comp) == self.k: yield comp

    def run(self):
        for comp in self.placements():
            if time.monotonic() > self.deadline: self.complete = False; break
            self.search(comp)
        return self

    def search(self, comp):
        n, k = self.n, self.k
        rays = []
        for g in range(n):
            rays.append(('v', g)); rays += [('h', i) for i in range(sum(comp[:g]), sum(comp[:g + 1]))]
        self.rays, self.comp = rays, comp
        V = {('a', i): z3.Real('a%d' % i) for i in range(1, n)}
        V.update({('h', i): z3.Real('h%d' % i) for i in range(k)})
        self.V = V
        pos = [Lin() if r == ('v', 0) else Lin({r if r[0] == 'h' else ('a', r[1]): 1}) for r in rays]
        self.pos = pos
        s = z3.Solver(); s.set('timeout', self.check_ms); self.s = s; self.nm = 0
        N = len(rays)
        for j in range(N):                                   # strict cyclic order, a_0 = 0
            nxt = pos[j + 1] if j + 1 < N else Lin(const=TWO_PI)
            s.add((nxt - pos[j]).z3(V) > 0)
        alt = Lin()                                           # the vertex is flat (Kawasaki)
        for i in range(n): alt = alt + ((-1) ** i) * pos[rays.index(('v', i))]
        self.cong_z3(alt, s)
        # mod-pi expressions in the visible directions; hidden ones are filled in when made
        self.vis = [Lin({('a', r[1]): 1}) if r[0] == 'v' and r[1] else (Lin() if r[0] == 'v' else None) for r in rays]
        self.dfs({j: 0 for j in range(N)}, [], [], ())

    def cong_z3(self, e, s=None):
        self.nm += 1; m = z3.Int('m%d' % self.nm)
        c = e.z3(self.V) == m * ZPI
        if s is not None: s.add(c)
        return c

    def wedge(self, status, x, y):
        N = len(self.rays)
        return [(x + d) % N for d in range(1, (y - x) % N) if status[(x + d) % N] == 1]

    def theorem11(self, P, status, x, y):
        """The mod-pi value Theorem 11 gives for y from x across the present rays between them,
        P the position map (Lin per ray)."""
        W = self.wedge(status, x, y); sv = Lin()
        for w in W: sv = P[w] - sv
        return (2 * sv - P[x]) if len(W) % 2 else (P[x] + 2 * sv)

    def collision(self, S, T):
        """Condition (C) as z3 constraints (as in exact_local.Vertex.collision)."""
        pos, V = self.pos, self.V
        maps = [(1, Lin())]
        for j in S[1:]:
            e, o = maps[-1]; maps.append((-e, o + (2 * e) * pos[j]))
        def sector(t):
            for i in range(len(S)):
                a, b = S[i], S[(i + 1) % len(S)]
                if t == a or (a < t < b if a < b else (t > a or t < b)): return i
        e, o = maps[sector(T[0])]
        rho = e * pos[T[0]] + o
        R = sorted(set(S) | set(T)); cur, mv = 0, []
        for i in range(len(R)):
            mv.append(cur)
            if R[(i + 1) % len(R)] in T: cur ^= 1
        cs = []
        for rr in (rho, rho + PI):
            st = {0: [], 1: []}
            for i in range(len(R)):
                a = pos[R[i]]; b = pos[R[i + 1]] if i + 1 < len(R) else pos[R[0]] + TWO_PI
                e, o = maps[sector(R[i])]
                u = e * (rr - o) - a
                self.nm += 1; q = z3.Int('q%d' % self.nm); w = u.z3(V) + 2 * ZPI * q
                cs += [w >= 0, w < 2 * ZPI]
                st[mv[i]].append(z3.And(w > 0, w < (b - a).z3(V)))
            cs.append(z3.Or(z3.Not(z3.Or(st[0])), z3.Not(z3.Or(st[1]))))
        return cs

    def candidates_T(self, status):
        """Steps from the current state, filtered combinatorially: even size, at least one anchor
        (a visible ray not yet made, or a hidden ray erased) unless scaffold steps are enumerated,
        at most one new ray per sector of the state (a sector has angle <= pi, so it holds at most
        one preimage of a line), a line from the empty state, and visible ray 0 in the first step
        that contains a visible ray."""
        N, rays = len(self.rays), self.rays
        S = [j for j in range(N) if status[j] == 1]
        cand = [j for j in range(N) if status[j] == 0 or (status[j] == 1 and rays[j][0] == 'h')]
        first_visible = not any(status[j] for j in range(N) if rays[j][0] == 'v')
        def rot_ok(T):
            return not first_visible or ('v', 0) in [rays[t] for t in T] or all(rays[t][0] == 'h' for t in T)
        if not S:
            for x, y in itertools.combinations(cand, 2):
                if ('v', 0) in (rays[x], rays[y]) or (self.scaffold and rays[x][0] == rays[y][0] == 'h'): yield [x, y]
            return
        sec = {}
        for t in cand:
            if status[t] == 0:
                for i in range(len(S)):
                    a, b = S[i], S[(i + 1) % len(S)]
                    if (a < t < b) if a < b else (t > a or t < b): sec[t] = i; break
        for size in range(2, len(cand) + 1, 2):
            for T in itertools.combinations(cand, size):
                new = [t for t in T if status[t] == 0]
                if len({sec[t] for t in new}) < len(new): continue
                if not self.scaffold and not any((rays[t][0] == 'v') or status[t] == 1 for t in T): continue
                if not rot_ok(T): continue
                yield list(T)

    def dfs(self, status, hist, req, scaf):
        """req: the congruences (Lin, = 0 mod pi) so far; scaf: the indices of the scaffold steps."""
        if time.monotonic() > self.deadline: self.complete = False; return
        self.nodes += 1
        N, rays = len(self.rays), self.rays
        if all(status[j] == 1 if rays[j][0] == 'v' else status[j] == 2 for j in range(N)):
            t = {'n': self.n, 'k': self.k, 'rays': [list(r) for r in rays], 'steps': [list(t) for t in hist]}
            if not scaf:
                t['req'] = sorted(set(v for v in (vec(l, self.n) for l in req) if v is not None))
            else:
                e = eliminate(req, self.n, scaf)
                if e is None: self.unpinned += 1; return
                t['req'], t['scaffold'] = e
            self.types.append(t)
            return
        if self.scaffold and not scaf and sum(1 for j in range(N) if rays[j][0] == 'h' and status[j] == 0) < 2:
            return                                             # no scaffold step can follow
        for T in self.candidates_T(status):
            S = [j for j in range(N) if status[j] == 1]
            # symbolic mod-pi values (for the index): fill new hidden rays from an anchor; a scaffold
            # step has none, and the direction of its first ray becomes a new variable
            vis = list(self.vis)
            anchor = next((t for t in T if vis[t] is not None and (rays[t][0] == 'v' or status[t] == 1)), None)
            scaf2 = scaf
            if anchor is None:
                anchor = T[0]; vis[anchor] = Lin({('l', len(hist)): 1}); scaf2 = scaf + (len(hist),)
            newreq = []; ok = True
            for y in T:
                if y == anchor: continue
                val = self.theorem11(vis, status, anchor, y)
                if vis[y] is None: vis[y] = val
                else:
                    d = vis[y] - val
                    if d.coef: newreq.append(d)
            S2 = sorted([j for j in S if j not in T] + [j for j in T if status[j] == 0])
            alt = Lin()
            for i, j in enumerate(S2): alt = alt + ((-1) ** i) * vis[j]
            if alt.coef: newreq.append(alt)
            # exact constraints in z3
            cs = []
            x = T[0]
            for y in T[1:]:
                W = self.wedge(status, x, y); sv = Lin()
                for w in W: sv = self.pos[w] - sv
                e = (self.pos[x] + self.pos[y] - 2 * sv) if len(W) % 2 else (self.pos[y] - self.pos[x] - 2 * sv)
                cs.append(self.cong_z3(e))
            alt = Lin()
            for i, j in enumerate(S2): alt = alt + ((-1) ** i) * self.pos[j]
            cs.append(self.cong_z3(alt))
            if S: cs += self.collision(S, T)
            self.s.push(); self.s.add(cs); self.checks += 1
            r = self.s.check()
            if r == z3.unknown: self.unknown += 1; self.complete = False
            if r == z3.sat:
                st2 = dict(status)
                for j in T: st2[j] = 1 if status[j] == 0 else 2
                saved = self.vis; self.vis = vis
                self.dfs(st2, hist + [T], req + newreq, scaf2)
                self.vis = saved
            self.s.pop()
            if time.monotonic() > self.deadline: self.complete = False; return


def eliminate(req, n, scaf):
    """Eliminate the scaffold lines from the congruences req (Lin in ('a', i) and ('l', step),
    each = 0 mod pi), by unimodular integer row operations: Euclid's algorithm on each line's column
    leaves at most one row with a nonzero entry there (its pivot row), which is set aside. Rows left
    with no line are the congruences in the visible directions alone: by back substitution through
    the pivot rows (triangular in the lines), the vertices with some lines satisfying req are exactly
    those satisfying them. None if some line has no pivot (it is not pinned: a continuum of lines).
    Returns (req vectors, {'steps': scaf, 'pivots': [[j, c, d]]}): pivot row j reads
    c.a + d.lam = 0 (mod pi) with d[j] != 0 and d[i] = 0 for i < j."""
    m, s = n - 1, len(scaf)
    rem = [[l.coef.get(('a', i), 0) for i in range(1, n)] + [l.coef.get(('l', st), 0) for st in scaf] for l in req]
    assert all(isinstance(x, int) for r in rem for x in r)
    piv = []
    for j in range(s):
        col = m + j
        while True:
            nz = sorted((r for r in rem if r[col]), key=lambda r: abs(r[col]))
            if len(nz) <= 1: break
            p = nz[0]
            for r in nz[1:]:
                q = r[col] // p[col]
                for i in range(m + s): r[i] -= q * p[i]
        nz = [r for r in rem if r[col]]
        if not nz: return None
        rem.remove(nz[0]); piv.append([j, nz[0][:m], nz[0][m:]])
    out = set()
    for r in rem:
        v = r[:m]
        if any(v):
            if next(x for x in v if x) < 0: v = [-x for x in v]
            out.add(tuple(v))
    return sorted(out), {'steps': list(scaf), 'pivots': piv}


def build(n, k, outdir, seconds, scaffold=False):
    t0 = time.monotonic()
    B = Builder(n, k, t0 + seconds, scaffold=scaffold).run()
    os.makedirs(outdir, exist_ok=True)
    meta = {'n': n, 'k': k, 'types': len(B.types), 'nodes': B.nodes, 'checks': B.checks, 'unknown': B.unknown,
            'complete': B.complete, 'seconds': round(time.monotonic() - t0, 1), 'z3': z3.get_version_string()}
    if scaffold: meta.update(scaffold=True, unpinned=B.unpinned)
    with open(os.path.join(outdir, 'deg%d_k%d%s.json' % (n, k, '_s' if scaffold else '')), 'w') as f:
        json.dump({'meta': meta, 'types': B.types}, f)
    return meta


# ---------------------------------------------------------------- runtime
class Table:
    def __init__(self, outdir, rigid=False, cells=None):
        """rigid: False, or 'halfturn' / 'reverse': keep only histories realizable at the vertex with
        rigid half-turns, or with half-turns and reverse folds (src/rigid.py)."""
        self.rigid = rigid
        self.by_n, self.by_nk = {}, {}; self.meta = []
        self.by_file = {}
        for fn in sorted(os.listdir(outdir)):
            if not fn.endswith('.json'): continue
            d = json.load(open(os.path.join(outdir, fn))); self.meta.append(d['meta'])
            self.by_file[fn] = d['types']
            for t in d['types']:
                t['req'] = [tuple(v) for v in t['req']]
                self.by_n.setdefault(t['n'], []).append(t)
                self.by_nk.setdefault((t['n'], t['k']), []).append(t)
        self.vectors = {n: sorted({v for t in ts for v in t['req']}) for n, ts in self.by_n.items()}
        self.cache = {}; self.stats = {'lookups': 0, 'cache_hits': 0, 'replays': 0, 'not_rigid': 0, 'tree_lookups': 0, 'searched': 0}
        if cells:                                       # exact realizability cells (src/cells.py), per table file
            for fn in sorted(os.listdir(outdir)):
                cf = os.path.join(cells, fn)
                if not fn.endswith('.json') or not os.path.exists(cf): continue
                C = json.load(open(cf)); T = json.load(open(os.path.join(outdir, fn)))['types']
                if not C['meta'].get('complete') or len(C['trees']) != len(T): continue
                for t, tr in zip(self.by_file[fn], C['trees']): t['tree'] = tr; t['den'] = C['meta'].get('den', 4)

    def realizable(self, t, a, h, u=None):
        """Realizability of a valid type at the vertex (a scaffold type: with the lines of branch u):
        from its exact cells when available (2 rigidly, 1 only with reverse folds), else by search
        (counted in stats['searched'])."""
        if 'tree' in t:
            v = _tree_value(t['tree'], a, t['den'], u); self.stats['tree_lookups'] += 1
            if v is not None: return v == 2 or (v == 1 and self.rigid == 'reverse')
        self.stats['searched'] += 1
        return rigid_ok(t, a, h, self.rigid == 'reverse')

    def valid(self, a, k=None):
        """Valid (type, hidden directions) at the vertex with visible directions a (a[0] = 0), with
        k hidden rays (all k if None)."""
        n = len(a); key = (tuple(round(x, 9) for x in a), k)
        self.stats['lookups'] += 1
        if key in self.cache: self.stats['cache_hits'] += 1; return self.cache[key]
        sat = set()
        for v in self.vectors.get(n, []):
            x = sum(c * a[i + 1] for i, c in enumerate(v))
            if abs((x + PI / 2) % PI - PI / 2) < 1e-6: sat.add(v)
        out = []
        for t in (self.by_n.get(n, []) if k is None else self.by_nk.get((n, k), [])):
            if all(v in sat for v in t['req']):
                seen = set()
                for u, lams in (branches(t, a) if 'scaffold' in t else [(None, None)]):
                    self.stats['replays'] += 1
                    h = replay(t, a, lams)
                    if h is None: continue
                    hk = tuple(round(x, 9) for x in h)
                    if hk in seen: continue
                    seen.add(hk)
                    if self.rigid and not self.realizable(t, a, h, u): self.stats['not_rigid'] += 1; continue
                    out.append((t, h))
        self.cache[key] = out
        return out


def _tree_value(tree, a, den=4, u=None):
    import cells
    return cells.evaluate(tree, a, den=den, u=u)


def rigid_ok(t, a, h, reverse=False):
    import rigid
    d = [a[r[1]] if r[0] == 'v' else h[r[1]] for r in t['rays']]
    return rigid.realizable_fast([[d[j] for j in T] for T in t['steps']], reverse=reverse)[0]


def branches(t, a):
    """The scaffold lines of a scaffold type at the vertex a, by back substitution through the
    pivot rows of its elimination: each line once per solution mod pi. Yields (u, {step: line mod pi}),
    u the branch (u_j in [0, |pivot j|) per line, as in cells.lines). Later rows substitute the
    unreduced lines, which are affine in a: reducing them mod pi first would make the labels u jump
    where a line crosses 0 = pi, and disagree with the exact cells."""
    sc = t['scaffold']; steps, piv = sc['steps'], sc['pivots']
    lam, us = [None] * len(steps), [None] * len(steps)
    def rec(i):
        if i < 0:
            yield tuple(us), {st: lam[j] % PI for j, st in enumerate(steps)}; return
        j, c, d = piv[i]
        rhs = -(sum(ci * a[q + 1] for q, ci in enumerate(c)) + sum(d[q] * lam[q] for q in range(len(steps)) if q != j and d[q]))
        for u in range(abs(d[j])):
            lam[j] = (rhs + PI * u) / d[j]; us[j] = u
            yield from rec(i - 1)
        lam[j] = us[j] = None
    yield from rec(len(piv) - 1)


def replay(t, a, lams=None):
    """Check type t at the vertex a; return the hidden directions (by hidden index), or None.
    lams: for a scaffold type, {step index: direction mod pi of the step's first ray}."""
    rays, N = t['rays'], len(t['rays'])
    d = [a[r[1]] if r[0] == 'v' else None for r in rays]
    status = [0] * N
    for i, T in enumerate(t['steps']):
        S = [j for j in range(N) if status[j] == 1]
        Sd = sorted(d[j] for j in S)
        phis = fh.maps(Sd)
        if phis is None: return None
        anchor = next((j for j in T if rays[j][0] == 'v' or status[j] == 1), None)
        if anchor is None:                                   # a scaffold step: its first ray from its line
            anchor = T[0]; x = lams[i]
            if S:
                p = max((q for q in S if q < anchor), default=max(S)); q = min((q for q in S if q > anchor), default=min(S))
                c = [y for y in (x, x + PI) if fh.between(y, d[p], d[q])]
                if len(c) != 1: return None
                x = c[0]
            d[anchor] = x
        lam = fh.apply(phis[fh.sector_of(d[anchor], Sd)], d[anchor]) if S else d[anchor]
        for j in T:
            if j != anchor and status[j] == 0 and rays[j][0] == 'h':
                if not S: d[j] = fh.nrm(lam + PI); continue
                p = max((i for i in S if i < j), default=max(S)); q = min((i for i in S if i > j), default=min(S))
                f = phis[fh.sector_of(d[p], Sd)]
                c = [x for x in (fh.inverse_apply(f, lam), fh.inverse_apply(f, fh.nrm(lam + PI))) if fh.between(x, d[p], d[q])]
                if len(c) != 1: return None
                d[j] = c[0]
        imgs = [fh.apply(phis[fh.sector_of(d[j], Sd)], d[j]) % PI for j in T] if S else [d[j] % PI for j in T]
        if any(abs((x - imgs[0] + PI / 2) % PI - PI / 2) > 1e-6 for x in imgs): return None
        if S and not fh.collision_free(Sd, [d[j] for j in T], phis, imgs[0]): return None
        for j in T: status[j] = 1 if status[j] == 0 else 2
        if fh.maps(sorted(d[j] for j in range(N) if status[j] == 1)) is None: return None
    order = [d[j] for j in range(N)]                       # the cyclic order of the type
    if any(not (0 <= order[j] < TWO_PI) for j in range(N)): return None
    if any(not (order[j + 1] - order[j] > 1e-7) for j in range(N - 1)) or TWO_PI - order[-1] <= 1e-7: return None
    return [d[j] for j in range(N) if rays[j][0] == 'h']


def options_for(table, P):
    """A replacement for Pattern.options: options at p with k hidden rays from the table (no collapses)."""
    memo = {}
    def options(p, k, c=0):
        if c: return []
        if (p, k) not in memo:                         # per depth, on demand
            items = P.H[p]; n = len(items); D = [x for x, _ in items]; E = [e for _, e in items]
            byR = {}
            for r in range(n):
                a = [fh.nrm(D[(i + r) % n] - D[r]) for i in range(n)]
                for t, hdirs in table.valid(a, k):
                    H = [fh.nrm(h + D[r]) for h in hdirs]
                    R = tuple(sorted(H, key=lambda x: round(x, 9)))
                    rk = tuple(round(x, 9) for x in R)
                    o = byR.setdefault(rk, [R, set()])
                    canon = {round(x, 9): x for x in o[0]}      # one float per hidden direction per option
                    made = set(); steps = []
                    for T in t['steps']:
                        st = []
                        for j in T:
                            ray = t['rays'][j]
                            if ray[0] == 'v': st.append(E[(ray[1] + r) % n])
                            else: st.append(((p, canon[round(H[ray[1]], 9)]), 2 if j in made else 1))
                        made |= set(T); steps.append(frozenset(st))
                    o[1].add(tuple(steps))
            memo[p, k] = [(R, sorted(hs, key=str), 0) for R, hs in byR.values()]
        return memo[p, k]
    return options


def main(argv):
    if argv[0] == 'build':
        n, k, outdir = int(argv[1]), int(argv[2]), argv[3]
        secs = float(next((a[10:] for a in argv if a.startswith('--seconds=')), 3600))
        print(json.dumps(build(n, k, outdir, secs, '--scaffold' in argv)), flush=True)
    elif argv[0] == 'stats':
        T = Table(argv[1])
        for m in T.meta: print(json.dumps(m))
        print({n: (len(ts), len(T.vectors[n])) for n, ts in T.by_n.items()})


if __name__ == '__main__':
    main(sys.argv[1:])
