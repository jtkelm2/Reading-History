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
    def __init__(self, n, k, deadline, check_ms=20000):
        self.n, self.k, self.deadline, self.check_ms = n, k, deadline, check_ms
        self.types, self.nodes, self.checks, self.unknown, self.complete = [], 0, 0, 0, True

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
        self.dfs({j: 0 for j in range(N)}, [], [])

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
        (a visible ray not yet made, or a hidden ray erased), at most one new ray per sector of the
        state (a sector has angle <= pi, so it holds at most one preimage of a line), a line from
        the empty state, and visible ray 0 in the first step."""
        N, rays = len(self.rays), self.rays
        S = [j for j in range(N) if status[j] == 1]
        cand = [j for j in range(N) if status[j] == 0 or (status[j] == 1 and rays[j][0] == 'h')]
        if not S:
            for x, y in itertools.combinations(cand, 2):
                if ('v', 0) in (rays[x], rays[y]): yield [x, y]
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
                if not any((rays[t][0] == 'v') or status[t] == 1 for t in T): continue
                yield list(T)

    def dfs(self, status, hist, req):
        if time.monotonic() > self.deadline: self.complete = False; return
        self.nodes += 1
        N, rays = len(self.rays), self.rays
        if all(status[j] == 1 if rays[j][0] == 'v' else status[j] == 2 for j in range(N)):
            self.types.append({'n': self.n, 'k': self.k, 'rays': [list(r) for r in rays], 'steps': [list(t) for t in hist],
                               'req': sorted(set(req))})
            return
        for T in self.candidates_T(status):
            S = [j for j in range(N) if status[j] == 1]
            # symbolic mod-pi values (for the index): fill new hidden rays from an anchor
            vis = list(self.vis); anchor = next(t for t in T if vis[t] is not None and (rays[t][0] == 'v' or status[t] == 1))
            newreq = []; ok = True
            for y in T:
                if y == anchor: continue
                val = self.theorem11(vis, status, anchor, y)
                if vis[y] is None: vis[y] = val
                else:
                    v = vec(vis[y] - val, self.n)
                    if v is not None: newreq.append(v)
            S2 = sorted([j for j in S if j not in T] + [j for j in T if status[j] == 0])
            alt = Lin()
            for i, j in enumerate(S2): alt = alt + ((-1) ** i) * vis[j]
            v = vec(alt, self.n)
            if v is not None: newreq.append(v)
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
                self.dfs(st2, hist + [T], req + newreq)
                self.vis = saved
            self.s.pop()
            if time.monotonic() > self.deadline: self.complete = False; return


def build(n, k, outdir, seconds):
    t0 = time.monotonic()
    B = Builder(n, k, t0 + seconds).run()
    os.makedirs(outdir, exist_ok=True)
    meta = {'n': n, 'k': k, 'types': len(B.types), 'nodes': B.nodes, 'checks': B.checks, 'unknown': B.unknown,
            'complete': B.complete, 'seconds': round(time.monotonic() - t0, 1), 'z3': z3.get_version_string()}
    with open(os.path.join(outdir, 'deg%d_k%d.json' % (n, k)), 'w') as f:
        json.dump({'meta': meta, 'types': B.types}, f)
    return meta


# ---------------------------------------------------------------- runtime
class Table:
    def __init__(self, outdir, rigid=False):
        """rigid: keep only histories that are rigidly realizable at the vertex (src/rigid.py)."""
        self.rigid = rigid
        self.by_n = {}; self.meta = []
        for fn in sorted(os.listdir(outdir)):
            if not fn.endswith('.json'): continue
            d = json.load(open(os.path.join(outdir, fn))); self.meta.append(d['meta'])
            for t in d['types']:
                t['req'] = [tuple(v) for v in t['req']]
                self.by_n.setdefault(t['n'], []).append(t)
        self.vectors = {n: sorted({v for t in ts for v in t['req']}) for n, ts in self.by_n.items()}
        self.cache = {}; self.stats = {'lookups': 0, 'cache_hits': 0, 'replays': 0, 'not_rigid': 0}

    def valid(self, a):
        """Valid (type, hidden directions) at the vertex with visible directions a (a[0] = 0)."""
        n = len(a); key = tuple(round(x, 9) for x in a)
        self.stats['lookups'] += 1
        if key in self.cache: self.stats['cache_hits'] += 1; return self.cache[key]
        sat = set()
        for v in self.vectors.get(n, []):
            x = sum(c * a[i + 1] for i, c in enumerate(v))
            if abs((x + PI / 2) % PI - PI / 2) < 1e-6: sat.add(v)
        out = []
        for t in self.by_n.get(n, []):
            if all(v in sat for v in t['req']):
                self.stats['replays'] += 1
                h = replay(t, a)
                if h is not None and self.rigid and not rigid_ok(t, a, h): self.stats['not_rigid'] += 1; continue
                if h is not None: out.append((t, h))
        self.cache[key] = out
        return out


def rigid_ok(t, a, h):
    import rigid
    d = [a[r[1]] if r[0] == 'v' else h[r[1]] for r in t['rays']]
    return rigid.realizable([[d[j] for j in T] for T in t['steps']])[0]


def replay(t, a):
    """Check type t at the vertex a; return the hidden directions (by hidden index), or None."""
    rays, N = t['rays'], len(t['rays'])
    d = [a[r[1]] if r[0] == 'v' else None for r in rays]
    status = [0] * N
    for T in t['steps']:
        S = [j for j in range(N) if status[j] == 1]
        Sd = sorted(d[j] for j in S)
        phis = fh.maps(Sd)
        if phis is None: return None
        anchor = next(j for j in T if rays[j][0] == 'v' or status[j] == 1)
        lam = fh.apply(phis[fh.sector_of(d[anchor], Sd)], d[anchor]) if S else d[anchor]
        for j in T:
            if status[j] == 0 and rays[j][0] == 'h':
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
        if p not in memo:
            items = P.H[p]; n = len(items); D = [x for x, _ in items]; E = [e for _, e in items]
            byk = {}
            for r in range(n):
                a = [fh.nrm(D[(i + r) % n] - D[r]) for i in range(n)]
                for t, hdirs in table.valid(a):
                    H = [fh.nrm(h + D[r]) for h in hdirs]
                    R = tuple(sorted(H, key=lambda x: round(x, 9)))
                    rk = tuple(round(x, 9) for x in R)
                    o = byk.setdefault(t['k'], {}).setdefault(rk, [R, set()])
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
            memo[p] = {kk: [(R, sorted(hs, key=str), 0) for R, hs in d.values()] for kk, d in byk.items()}
        return memo[p].get(k, [])
    return options


def main(argv):
    if argv[0] == 'build':
        n, k, outdir = int(argv[1]), int(argv[2]), argv[3]
        secs = float(next((a[10:] for a in argv if a.startswith('--seconds=')), 3600))
        print(json.dumps(build(n, k, outdir, secs)), flush=True)
    elif argv[0] == 'stats':
        T = Table(argv[1])
        for m in T.meta: print(json.dumps(m))
        print({n: (len(ts), len(T.vectors[n])) for n, ts in T.by_n.items()})


if __name__ == '__main__':
    main(sys.argv[1:])
