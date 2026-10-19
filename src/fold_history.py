"""Structure identification by local fold histories, for crease patterns in FOLD format.

Usage: python3 src/fold_history.py [--out=PREFIX] [--budget=SECONDS] [--tol=RADIANS] [--triage]
                                   [--geometry=shortest|found] [--table=DIR] pattern.fold
--table=DIR takes local histories from the lookup tables in DIR (src/lookup.py) instead of
searching for them at each vertex; vertices not covered by the tables get no options.
--rigid[=reverse|halfturn] (with --table) keeps only local histories realizable with rigid half-turns
and reverse folds (default), or with half-turns only (src/rigid.py).
The analysis fixes the minimal geometry whose cancelled segments are shortest (default), or the one
the solver found.

Moves (shadow model: a flat folded state is a map that reflects across every crease):
  single-line fold   new creases at a vertex are preimages of one fold line under the current
                     folding maps, an even number of them, and no straddling sectors collide;
  cancellation       a fold lands on an existing crease and erases it (Lemma 4); the erased
                     segment is a hidden pair, traced across the sheet (it bends at earlier creases,
                     passes later ones, and ends at the boundary or at another vertex's hidden ray);
  collapse           several creases at once with any image lines; costs (lines - 1).
Each interior vertex needs a local history of such moves; one global order of creases must be
consistent with all of them.  Cost: collapses first, then hidden pairs.
Directions are real numbers, snapped to multiples of 2pi/N when some N <= 64 fits.
High-degree vertices are encoded symbolically (wedge-product and parity conditions) and checked
lazily for collisions.  Vertices get more hidden pairs or collapses only when they appear in an
unsatisfiable core.  Output: the folds (creases forced to be simultaneous in every solution of the
chosen resolution), the forced order between folds, and its Hasse diagram.
"""
import json, math, sys, itertools, shutil, subprocess, time
import z3

TWO_PI = 2 * math.pi
TOL = 1e-6          # input clustering tolerance (radians)
EPS = 1e-7          # equality of derived directions (after snapping)
MAXK = 4            # hidden pairs at one vertex
MAXBEND = 2
MAXN = 64           # try to snap to multiples of 2pi/N for N <= MAXN


def nrm(x):
    x %= TWO_PI
    return 0.0 if TWO_PI - x < EPS else x

def eq(a, b): return abs((a - b + math.pi) % TWO_PI - math.pi) < EPS

def linekey(x):
    y = x % math.pi
    if math.pi - y < EPS: y = 0.0
    return round(y / EPS / 10)

def between(x, a, b):
    if eq(a, b): return not eq(x, a)
    da, db = (x - a) % TWO_PI, (b - a) % TWO_PI
    return EPS < da < db - EPS


def maps(S):
    if not S: return [(1, 0.0)]
    out = [(1, 0.0)]
    for c in S[1:]:
        s, o = out[-1]; out.append((-s, nrm(2 * s * c + o)))
    s, o = out[-1]
    s2, o2 = -s, nrm(2 * s * S[0] + o)
    return out if s2 == 1 and eq(o2, 0.0) else None

def sector_of(x, S):
    if not S: return 0
    m = len(S)
    for j in range(m):
        if eq(x, S[j]) or between(x, S[j], S[(j + 1) % m]): return j
    raise ValueError(x)

def apply(f, x): return nrm(f[0] * x + f[1])
def inverse_apply(f, y): return nrm(f[0] * (y - f[1]))

def collision_free(S, switch, phis, line):
    R = sorted(set(S) | set(switch)); m = len(R)
    sw = lambda x: any(eq(x, y) for y in switch)
    mv, cur = [], 0
    for k in range(m):
        mv.append(cur)
        if sw(R[(k + 1) % m]): cur ^= 1
    if cur: return False
    for rho in (line, nrm(line + math.pi)):
        seen = set()
        for k in range(m):
            a, b = R[k], R[(k + 1) % m]
            if between(inverse_apply(phis[sector_of(a, S)], rho), a, b): seen.add(mv[k])
        if len(seen) > 1: return False
    return True

def histories(items, limit=20000):
    """Single-line local histories of the visible rays items = [(direction, key)], with no hidden
    rays: every weak order whose classes, in order, are steps satisfying (P), (L) and (C)."""
    by = dict((k, d) for d, k in items); res = []
    def rec(S, rem, hist):
        if len(res) >= limit: return
        if not rem: res.append(tuple(hist)); return
        phis = maps(S)
        if phis is None: return
        cand = {}
        for k in rem:
            img = apply(phis[sector_of(by[k], S)], by[k])
            cand.setdefault(linekey(img), (img % math.pi, []))[1].append(k)
        for key, (line, ks) in cand.items():
            for size in range(2, len(ks) + 1, 2):
                for T in itertools.combinations(ks, size):
                    dirs = [by[k] for k in T]
                    if not collision_free(S, dirs, phis, line): continue
                    hist.append(frozenset(T)); rec(sorted(S + dirs), rem - set(T), hist); hist.pop()
    rec([], frozenset(k for _, k in items), [])
    return res


class Pattern:
    def __init__(self, path, tol=TOL):
        d = json.load(open(path))
        self.V = V = [tuple(map(float, p[:2])) for p in d['vertices_coords']]
        A = [a.upper() for a in (d.get('edges_assignment') or ['U'] * len(d['edges_vertices']))]
        self.crease = [(u, v, a) for (u, v), a in zip(d['edges_vertices'], A) if a in 'MVU']
        self.bound = [(u, v) for (u, v), a in zip(d['edges_vertices'], A) if a == 'B']
        raw = {}
        for i, (u, v, _) in enumerate(self.crease):
            for p, q in ((u, v), (v, u)):
                raw.setdefault(p, []).append((math.atan2(V[q][1] - V[p][1], V[q][0] - V[p][0]) % TWO_PI, i))
        # cluster line directions (mod pi) with tolerance, then try to snap to a cyclic group
        lines = sorted(t % math.pi for hs in raw.values() for t, _ in hs)
        reps, cl = [], []
        for t in lines:
            if cl and t - cl[-1] < tol: cl.append(t)
            else:
                if cl: reps.append(sum(cl) / len(cl))
                cl = [t]
        reps.append(sum(cl) / len(cl))
        if len(reps) > 1 and math.pi - reps[-1] + reps[0] < tol: reps.pop()
        self.N = None
        base = reps[0]
        for N in range(2, MAXN + 1, 2):
            u = TWO_PI / N
            if all(abs(((r - base) / u) - round((r - base) / u)) * u < max(tol, 1e-9) for r in reps):
                self.N = N; break
        def snap(t):
            if self.N:
                u = TWO_PI / self.N
                return nrm(base + round((t - base) / u) * u)
            r = min(reps, key=lambda r: abs((t % math.pi) - r))
            return nrm(r + (math.pi if abs(t - r) > math.pi / 2 and abs(t - r - math.pi) < math.pi / 2 else 0)) \
                if abs((t % math.pi) - r) < tol else nrm(t)
        self.H = {p: sorted((snap(t), e) for t, e in hs) for p, hs in raw.items()}
        self.ndirs = len(reps)
        onb = {x for e in self.bound for x in e}
        self.interior = sorted(p for p in self.H if p not in onb)
        self.cdir = {}
        for p, hs in self.H.items():
            for t, e in hs:
                if self.crease[e][0] == p: self.cdir[e] = t

    def deg(self, r): return round(math.degrees(r) % 360, 3)

    def trace(self, o, r, skip_v=None, skip_c=None):
        V = self.V; dv = (math.cos(r), math.sin(r))
        def hit(a, b):
            e = (b[0] - a[0], b[1] - a[1]); den = dv[0] * e[1] - dv[1] * e[0]
            if abs(den) < 1e-12: return None
            w = (a[0] - o[0], a[1] - o[1])
            return (w[0] * e[1] - w[1] * e[0]) / den, (w[0] * dv[1] - w[1] * dv[0]) / den
        tb = math.inf
        for u, v in self.bound:
            h = hit(V[u], V[v])
            if h and h[0] > 1e-9 and -1e-9 <= h[1] <= 1 + 1e-9: tb = min(tb, h[0])
        if tb == math.inf: return None
        scale = max(1.0, max(abs(c) for q in V for c in q))
        ev = {}
        for w in range(len(V)):
            if w == skip_v: continue
            q = (V[w][0] - o[0], V[w][1] - o[1]); t = q[0] * dv[0] + q[1] * dv[1]
            if 1e-7 * scale < t < tb - 1e-7 * scale and abs(q[0] * dv[1] - q[1] * dv[0]) < 1e-7 * scale:
                ev[round(t / scale, 9)] = ('v', w, [(e, eq(abs((dd - r) % math.pi), math.pi / 2)) for dd, e in self.H.get(w, [])])
        for i, (u, v, _) in enumerate(self.crease):
            if i == skip_c or (skip_v is not None and skip_v in (u, v)): continue
            h = hit(V[u], V[v])
            if h is None:
                q = (V[u][0] - o[0], V[u][1] - o[1])
                if abs(q[0] * dv[1] - q[1] * dv[0]) < 1e-7 * scale:
                    ts = sorted([q[0] * dv[0] + q[1] * dv[1], (V[v][0] - o[0]) * dv[0] + (V[v][1] - o[1]) * dv[1]])
                    if ts[1] > 1e-7 * scale and ts[0] < tb - 1e-7 * scale: ev[round(max(ts[0], 0) / scale + 1e-6, 9)] = ('a', i)
                continue
            t, s = h
            if 1e-7 * scale < t < tb - 1e-7 * scale and 1e-7 < s < 1 - 1e-7:
                zd = self.cdir[i]
                ev[round(t / scale, 9)] = ('x', i, eq(abs((zd - r) % math.pi), math.pi / 2),
                                           (o[0] + t * dv[0], o[1] + t * dv[1]), zd)
        return [ev[k] for k in sorted(ev)], (o[0] + tb * dv[0], o[1] + tb * dv[1])

    def options(self, p, k, c=0):
        """Options at p with k hidden rays and c collapse steps: [(R, histories or 'LAZY', collapse cost)],
        R the tuple of hidden directions."""
        reals = [(d, e) for d, e in self.H[p]]
        if k == 0 and c == 0:
            lim = LAZY_T if self.N else LAZY_T_REAL    # the symbolic encoding is slow without a unit
            hs = histories(reals, limit=lim + 1)
            if len(hs) > lim: return [((), 'LAZY', 0)]
            return [((), hs, 0)] if hs else []
        res = lazy_histories(p, reals, k, c)
        out = []
        for (R, gc), hs in sorted(res.items(), key=lambda kv: kv[0][1]):
            if any(self.trace(self.V[p], r) is None for r in R): continue
            # rename hidden keys to ((p, r), i) with r the float direction used in R
            def ren(key):
                if isinstance(key, tuple) and len(key) == 2 and isinstance(key[0], tuple):
                    r = [x for x in R if abs(round(math.degrees(x), 6) - key[0][1]) < 1e-5][0]
                    return ((p, r), key[1])
                return key
            out.append((R, [tuple(frozenset(ren(x) for x in st) for st in h) for h in hs], gc))
        return out


class Model:
    def __init__(self, P, opts):
        self.P, self.opts = P, opts
        self.s = z3.Solver(); self.T = {}
        self.sel = {(p, i): z3.Bool('s_%d_%d' % (p, i)) for p, ol in opts.items() for i in range(len(ol))}
        self.track = {}
        for p, ol in opts.items():
            a = z3.Bool('trk_%d' % p); self.track[a] = p
            alts = []
            for i, (R, hs, *_) in enumerate(ol):
                c = [z3.Or([self.history(h) for h in hs])] if hs != 'LAZY' else [self.symbolic(p)]
                c += [self.path(P.V[p], r, 0, self.t(((p, r), 1)), self.t(((p, r), 2)), p, None) for r in R]
                self.s.add(z3.Implies(self.sel[p, i], z3.And(c)))
            self.s.add(z3.Implies(a, z3.PbEq([(self.sel[p, i], 1) for i in range(len(ol))], 1)))
            self.s.add(z3.PbLe([(self.sel[p, i], 1) for i in range(len(ol))], 1))

    def symbolic(self, p):
        """Single-line local histories at p (real creases only), without enumeration:
        (P) every set of creases made at one time at p has even size;
        (L) if x and y are made together, the product of reflections in the creases strictly
            inside the ccw wedge x->y that are earlier than x maps x onto the line of y.
        The product of reflections in w_1 < ... < w_k is x -> 2s - x (k odd) or x + 2s (k even),
        with s = w_k - w_(k-1) + ... (alternating sum), so (L) is linear given the earlier set.
        The collision test (C) is checked lazily on models."""
        items = self.P.H[p]; n = len(items); cs = []
        N = self.P.N
        if N and N & (N - 1) == 0:                     # Z/N as bit-vectors: arithmetic mod N is free
            B = N.bit_length() - 1
            R = lambda x: z3.BitVecVal(round(x / (TWO_PI / N)) % N, B)
            for _, x in items:
                k = z3.Int('par_%d_%d' % (p, x))
                cs.append(z3.Sum([z3.If(self.t(z) == self.t(x), 1, 0) for _, z in items]) == 2 * k)
            for i in range(n):
                for j in range(i + 1, n):
                    (dx, x), (dy, y) = items[i], items[j]
                    sv, par = z3.BitVecVal(0, B), z3.BoolVal(False)
                    for dw, w in items[i + 1:j]:
                        e = self.t(w) < self.t(x)
                        sv = z3.If(e, R(dw) - sv, sv); par = z3.If(e, z3.Not(par), par)
                    expr = z3.If(par, 2 * sv - R(dx) - R(dy), R(dx) + 2 * sv - R(dy))
                    line0 = expr & z3.BitVecVal(N // 2 - 1, B) == 0
                    cs.append(z3.Implies(self.t(x) == self.t(y), line0))
            return z3.And(cs)
        if N:                                          # exact integer arithmetic in units of 2pi/N
            PI = z3.IntVal(N // 2); eps = z3.IntVal(0)
            R = lambda x: z3.IntVal(round(x / (TWO_PI / N)) % N)
        else:
            PI = z3.RealVal(repr(round(math.pi, 12))); eps = z3.RealVal('0.000001')
            R = lambda x: z3.RealVal(repr(round(x, 12)))
        for _, x in items:
            k = z3.Int('par_%d_%d' % (p, x))
            cs.append(z3.Sum([z3.If(self.t(z) == self.t(x), 1, 0) for _, z in items]) == 2 * k)
        for i in range(n):
            for j in range(i + 1, n):
                (dx, x), (dy, y) = items[i], items[j]
                sv, par = (z3.IntVal(0) if self.P.N else z3.RealVal(0)), z3.BoolVal(False)
                for dw, w in items[i + 1:j]:
                    e = self.t(w) < self.t(x)
                    sv = z3.If(e, R(dw) - sv, sv); par = z3.If(e, z3.Not(par), par)
                m = z3.Int('m_%d_%d_%d' % (p, x, y))
                expr = z3.If(par, 2 * sv - R(dx) - R(dy), R(dx) + 2 * sv - R(dy)) - m * PI
                cs.append(z3.Implies(self.t(x) == self.t(y), z3.And(expr >= -eps, expr <= eps)))
        return z3.And(cs)

    def t(self, k):
        if k not in self.T: self.T[k] = z3.Int('t_%s' % (k,))
        return self.T[k]

    def history(self, h):
        cs, reps = [], []
        for step in h:
            ks = sorted(step, key=str); reps.append(ks[0])
            cs += [self.t(k) == self.t(ks[0]) for k in ks[1:]]
        return z3.And(cs + [self.t(a) < self.t(b) for a, b in zip(reps, reps[1:])])

    def joins(self, w, r2, h1, h2):
        return [z3.And(self.sel[w, i], h1 == self.t(((w, Rw[j]), 1)), h2 == self.t(((w, Rw[j]), 2)))
                for i, (Rw, *_) in enumerate(self.opts.get(w, [])) for j in range(len(Rw)) if eq(Rw[j], r2)]

    def path(self, o, r, depth, h1, h2, skip_v, skip_c, m=None):
        """The cancelled segment leaving o in direction r, with creases at times h1 < h2, has a
        valid path.  Given a model m, only the branch that m takes (its first true one): this pins
        the segment's path, bends and endpoint."""
        P = self.P
        tr = P.trace(o, r, skip_v, skip_c)
        if tr is None: return z3.BoolVal(False)
        terms, pre = [], []
        for e in tr[0]:
            if e[0] == 'x':
                if e[2]: continue
                z = self.t(e[1])
                if depth < MAXBEND:
                    terms.append(z3.And(pre + [z < h1, self.path(e[3], nrm(2 * e[4] - r - math.pi), depth + 1, h1, h2, None, e[1], m)]))
                pre = pre + [z > h2]
            elif e[0] == 'v':
                terms += [z3.And(pre + [j]) for j in self.joins(e[1], nrm(r + math.pi), h1, h2)]
                pre = pre + [self.t(z) > h2 for z, perp in e[2] if not perp]
            else:
                pre = pre + [self.t(e[1]) > h2]
        terms.append(z3.And(pre))
        if m is None: return z3.Or(terms)
        return next((c for c in terms if z3.is_true(m.eval(c, model_completion=True))), z3.BoolVal(False))

    def geometry(self, choice, m):
        """The cancelled segments of model m, pinned: every hidden ray takes the path it takes in m."""
        return z3.And([self.path(self.P.V[p], r, 0, self.t(((p, r), 1)), self.t(((p, r), 2)), p, None, m)
                       for p, R in choice.items() for r in R])

    def cost(self):
        return z3.Sum([z3.IntVal(0)] + [z3.If(self.sel[p, i], len(o[0]) + W_COLLAPSE * o[2], 0) for p, ol in self.opts.items() for i, o in enumerate(ol)])


def first_failure(items, tv):
    """items: [(direction, key)] real creases at a vertex; tv: key -> time.  Simulate the steps in
    time order.  Returns None, or a learned conflict:
      ('pair', x, y, wedge, earlier)  x, y made together but their images differ; this depends only
                                      on which creases strictly inside the ccw wedge x->y are earlier;
      ('parity', T)                   an odd number of creases made at one time;
      ('exact', A, T)                 anything else (collision): earlier set A, step T."""
    S, done = [], []
    for t in sorted({tv[k] for _, k in items}):
        T = [(d, k) for d, k in items if tv[k] == t]
        if len(T) % 2: return ('parity', [k for _, k in T])
        phis = maps(S)
        img = [(linekey(apply(phis[sector_of(d, S)], d)), d, k) for d, k in T]
        bad = [z for z in img if z[0] != img[0][0]]
        if bad:
            (_, dx, x), (_, dy, y) = img[0], bad[0]
            wedge = [k for d, k in items if between(d, dx, dy)]
            return ('pair', x, y, wedge, [w for w in wedge if w in done])
        line = apply(phis[sector_of(T[0][0], S)], T[0][0]) % math.pi
        S2 = sorted(S + [d for d, _ in T])
        if not collision_free(S, [d for d, _ in T], phis, line) or maps(S2) is None:
            return ('exact', list(done), [k for _, k in T])
        S = S2; done = done + [k for _, k in T]
    return None

def learn(M, b, items, f):
    t = M.t
    if f[0] == 'pair':
        _, x, y, wedge, earlier = f
        return z3.Implies(b, z3.Or([t(x) != t(y)] + [(t(w) < t(x)) != (w in earlier) for w in wedge]))
    if f[0] == 'parity':
        T = f[1]; rest = [e for _, e in items if e not in T]
        return z3.Implies(b, z3.Or([t(x) != t(T[0]) for x in T[1:]] + [t(c) == t(T[0]) for c in rest]))
    _, A, T = f; rest = [e for _, e in items if e not in A and e not in T]
    return z3.Implies(b, z3.Or([t(a) >= t(T[0]) for a in A] + [t(x) != t(T[0]) for x in T[1:]] + [t(c) <= t(T[0]) for c in rest]))

KINDS = {}
def refine(P, M, opts, m, sinks):
    """Check every selected lazy vertex in model m; add learned clauses to all solvers in sinks."""
    n = 0
    for (p, i), b in M.sel.items():
        if opts[p][i][1] == 'LAZY' and z3.is_true(m.eval(b)):
            items = [(d, e) for d, e in P.H[p]]
            tv = {e: m.eval(M.t(e), model_completion=True).as_long() for _, e in items}
            f = first_failure(items, tv)
            if f:
                KINDS[f[0]] = KINDS.get(f[0], 0) + 1
                c = learn(M, b, items, f)
                for s_ in sinks: s_.add(c)
                n += 1
    return n

W_COLLAPSE = 100
LAZY_T = 2000
LAZY_T_REAL = 20000
MAXC = 1

def levels(kmax=MAXK):
    return [(c, k) for c in range(MAXC + 1) for k in range(kmax + 1)]

def solve(P, log=print, budget=600):
    """Core-guided expansion over the levels (c, k) = (collapse steps, hidden rays), c <= MAXC and
    k <= MAXK, in lexicographic order.  Each vertex starts at its first level with an option; a
    vertex in an unsat core is given the options of its next level.  The cost is then minimised over
    the options given so far.  lower_bound (each vertex at its first level, a collapse step costing
    at least 1) is a lower bound over the whole search space, so cost == lower_bound certifies that
    the minimum found is the minimum over the search space."""
    t0 = time.monotonic()
    LV = levels(); level, kmin, opts, ladder = {}, {}, {}, {}
    for p in P.interior:
        for li, (c, k) in enumerate(LV):
            o = P.options(p, k, c)
            if o:
                kmin[p] = li; level[p] = 0; opts[p] = list(o); ladder[p] = LV[li:]
                break
        else:
            return {'status': 'local fail', 'first': p}
    lower = sum(k + W_COLLAPSE * c for c, k in (LV[kmin[p]] for p in P.interior))
    left = lambda: max(1, int(1000 * (budget - (time.monotonic() - t0))))   # ms; every z3 call is bounded
    rounds = 0; lazy_rounds = [0]; rebuild = True
    while True:
        if rebuild:
            rounds += 1; M = Model(P, opts); rebuild = False
        M.s.set('timeout', left())
        r = M.s.check(list(M.track))
        if r == z3.sat:
            learned = refine(P, M, opts, M.s.model(), [M.s])
            if lazy_rounds[0] % 50 == 0: log('    lazy %d learned %d %s (%.0fs)' % (lazy_rounds[0], learned, KINDS, time.monotonic() - t0))
            if not learned: break
            lazy_rounds[0] += 1
            if time.monotonic() - t0 > budget: return {'status': 'timeout', 'rounds': rounds}
            continue
        if r != z3.sat and r != z3.unsat:
            return {'status': 'timeout' if time.monotonic() - t0 >= budget - 1 else 'unknown', 'rounds': rounds}
        core = [M.track[c] for c in M.s.unsat_core()]
        grown = []
        for p in core:
            seq = ladder[p]
            if level[p] + 1 < len(seq):
                level[p] += 1; c, k = seq[level[p]]
                opts[p] = opts[p] + list(P.options(p, k, c)); grown.append(p)
        rebuild = True
        log('  round %d: core %s, expanded %s (%.0fs)' % (rounds, core, grown, time.monotonic() - t0))
        if not grown: return {'status': 'global unsat', 'core': core, 'rounds': rounds}
        if time.monotonic() - t0 > budget: return {'status': 'timeout', 'rounds': rounds}
    o = z3.Optimize(); o.add(M.s.assertions()); o.add(list(M.track))
    c = M.cost(); o.minimize(c)
    while True:
        o.set('timeout', left())
        if o.check() != z3.sat: return {'status': 'timeout', 'rounds': rounds}
        m = o.model()
        if not refine(P, M, opts, m, [o, M.s]): break
        lazy_rounds[0] += 1
    ch = {p: opts[p][i][0] for (p, i), b in M.sel.items() if opts[p][i][0] and z3.is_true(m.eval(b))}
    col = sorted(p for (p, i), b in M.sel.items() if opts[p][i][2] and z3.is_true(m.eval(b)))
    cv = m.eval(c).as_long()
    gsum = sum(opts[p][i][2] for (p, i), b in M.sel.items() if z3.is_true(m.eval(b)))
    return {'status': 'ok', 'model': M, 'opts': opts, 'choice': ch, 'hidden': cv % W_COLLAPSE, 'collapses': col, 'excess_lines': gsum, 'm': m,
            'cost': cv, 'lower_bound': lower, 'certified_minimal': cv == lower,
            'rounds': rounds, 'lazy_refinements': lazy_rounds[0], 'expanded': sorted(p for p in level if level[p] > 0),
            'defects': sorted(p for p in kmin if kmin[p] > 0)}

def lazy_histories(p, reals, k, ctarget=0, limit=20000):
    """Histories at vertex p using exactly k hidden pairs, which are introduced only where a fold
    line has a preimage ray carrying no real crease.  reals: [(direction, key)].
    Returns {R: [history]} with R a tuple of hidden directions."""
    realdir = {key: d for d, key in reals}
    bis = []
    D = [d for d, _ in reals]
    for a, b in itertools.combinations(D, 2): bis.append(nrm((a + b) / 2))
    out = {}; count = [0]
    def hkey(x): return (p, round(math.degrees(x), 6))
    def rec(S, kinds, rem, alive, used, hist, ncol=0):
        # S: sorted ray directions present; kinds: dir-> ('real', key) or ('hid', hkey)
        if count[0] >= limit: return
        if not rem and not alive:
            if used == k and ncol == ctarget:
                R = tuple(sorted(hd for hd in used_dirs))
                out.setdefault((R, gcost[0]), []).append(tuple(hist)); count[0] += 1
            return
        phis = maps(S)
        if phis is None: return
        if ncol < ctarget:                  # a collapse: several real creases at once, any lines
            rl = sorted(rem, key=str)
            for size in range(2, len(rl) + 1, 2):
                for T in itertools.combinations(rl, size):
                    S2 = sorted(S + [realdir[x] for x in T])
                    if maps(S2) is None: continue
                    kl = len({linekey(apply(phis[sector_of(realdir[x], S)], realdir[x])) for x in T})
                    gcost[0] += max(1, kl - 1)
                    hist.append(frozenset(T)); rec(S2, kinds, rem - set(T), alive, used, hist, ncol + 1); hist.pop()
                    gcost[0] -= max(1, kl - 1)
        lines = {}
        for key in rem:
            d = realdir[key]; lines.setdefault(linekey(apply(phis[sector_of(d, S)], d)), apply(phis[sector_of(d, S)], d) % math.pi)
        for x in alive:
            v = apply(phis[sector_of(x, S)], x); lines.setdefault(linekey(v), v % math.pi)
        if used < k and not S:
            for b in bis: lines.setdefault(linekey(b), b % math.pi)
        elif used < k:
            for b in bis:
                v = apply(phis[sector_of(b, S)], b) if not any(eq(b, y) for y in S) else None
                if v is not None: lines.setdefault(linekey(v), v % math.pi)
        remdir = [(realdir[key], key) for key in rem]
        for lk, psi in lines.items():
            cand = []                       # (direction, kind, key)
            m = len(S)
            secs = [(S[j], S[(j + 1) % m], phis[j]) for j in range(m)] if S else [(0.0, 0.0, (1, 0.0))]
            for a, b, f in secs:
                for rho in (psi, nrm(psi + math.pi)):
                    x = inverse_apply(f, rho)
                    if (between(x, a, b) if S else True):
                        rk = [key for d, key in remdir if eq(d, x)]
                        if rk: cand.append((x, 'real', rk[0]))
                        elif used < k and not any(eq(x, y) for y in used_dirs): cand.append((x, 'new', None))
            for x in alive:
                if linekey(apply(phis[sector_of(x, S)], x)) == lk: cand.append((x, 'cancel', None))
            # dedupe by direction
            seen, c2 = [], []
            for c in cand:
                if not any(eq(c[0], y) for y in seen): seen.append(c[0]); c2.append(c)
            cand = c2
            for size in range(2, len(cand) + 1, 2):
                for T in itertools.combinations(cand, size):
                    nnew = sum(1 for c in T if c[1] == 'new')
                    if used + nnew > k: continue
                    if not any(c[1] == 'real' for c in T) and not any(c[1] == 'cancel' for c in T) and nnew < 2: continue
                    dirs = [c[0] for c in T]
                    if not collision_free(S, dirs, phis, psi): continue
                    S2 = [y for y in S if not any(c[1] == 'cancel' and eq(c[0], y) for c in T)]
                    S2 = sorted(S2 + [c[0] for c in T if c[1] != 'cancel'])
                    kinds2 = dict(kinds); step = []
                    alive2 = [y for y in alive if not any(c[1] == 'cancel' and eq(c[0], y) for c in T)]
                    rem2 = set(rem)
                    for c in T:
                        if c[1] == 'real': step.append(c[2]); rem2.discard(c[2])
                        elif c[1] == 'new': step.append((hkey(c[0]), 1)); alive2.append(c[0]); used_dirs.append(c[0])
                        else: step.append((hkey(c[0]), 2))
                    hist.append(frozenset(step))
                    rec(S2, kinds2, frozenset(rem2), alive2, used + nnew, hist, ncol)
                    hist.pop()
                    for c in T:
                        if c[1] == 'new': used_dirs.remove(c[0])
    used_dirs = []; gcost = [0]
    rec([], {}, frozenset(key for _, key in reals), [], 0, [])
    return out



# ---------------------------------------------------------------- analysis of the chosen resolution
class Checker:
    """Satisfiability over a set of resolutions, with lazy vertices refined until the model is valid.
    scope='geometry': the resolutions with the geometry of res, that is the option chosen at every
    vertex and, if pin, the path of every cancelled segment.  scope='minimal': every resolution of
    minimal cost among the options in res['opts']."""
    def __init__(self, P, res, pin=True, budget=3600, scope='geometry'):
        self.P, self.M, self.opts = P, res['model'], res['opts']
        self.deadline = time.monotonic() + budget; self.unknown = 0
        M, m = self.M, res['m']
        self.s = M.s
        if scope == 'minimal':
            g = z3.FreshBool('minimal'); self.s.add(z3.Implies(g, M.cost() == res['cost']))
            self.base = list(M.track) + [g]
        else:
            self.base = list(M.track) + [b if z3.is_true(m.eval(b)) else z3.Not(b) for b in M.sel.values()]
        if pin and scope != 'minimal':
            g = z3.FreshBool('geometry'); self.s.add(z3.Implies(g, M.geometry(res['choice'], m))); self.base.append(g)
        self.n = 0

    def check(self, extra=None):
        lits = list(self.base)
        if extra is not None:
            self.n += 1; a = z3.FreshBool('q'); self.s.add(z3.Implies(a, extra)); lits.append(a)
        while True:
            self.s.set('timeout', max(1, int(1000 * (self.deadline - time.monotonic()))))
            r = self.s.check(lits)
            if r == z3.unsat: return None
            if r != z3.sat: self.unknown += 1; return 'unknown'
            m = self.s.model()
            if not refine(self.P, self.M, self.opts, m, [self.s]): return m


def walk(P, M, m, o, r, k1, skip_v, skip_c, depth=0):
    val = lambda k: m.eval(M.t(k), model_completion=True).as_long()
    h1 = val((k1, 1))
    ev, end = P.trace(o, r, skip_v, skip_c)
    for e in ev:
        if e[0] == 'x' and not e[2] and val(e[1]) < h1 and depth < MAXBEND:
            pts, how = walk(P, M, m, e[3], nrm(2 * e[4] - r - math.pi), k1, None, e[1], depth + 1)
            return [o] + pts, how
        if e[0] == 'v' and any(z3.is_true(m.eval(j)) for j in M.joins(e[1], nrm(r + math.pi), M.t((k1, 1)), M.t((k1, 2)))):
            return [o, P.V[e[1]]], 'joins vertex %d' % e[1]
    return [o, end], 'boundary'


def analyse(P, res, log=print, pin=True, budget=3600, scope='geometry'):
    """Folds and forced order over the resolutions of Checker(scope).  With scope='minimal' only
    visible creases are analysed, since hidden ones need not exist in every minimal resolution.
    A query that z3 cannot decide within the budget counts as not proven, so the folds may be finer
    and the order smaller, never wrong; res['analysis_unknown'] is the number of such queries."""
    C = Checker(P, res, pin, budget, scope); M = C.M
    m = C.check()
    if not isinstance(m, z3.ModelRef): raise RuntimeError('analysis: the resolutions are not confirmed to exist')
    keys = list(range(len(P.crease)))
    if scope != 'minimal': keys += [((p, r), i) for p, R in res['choice'].items() for r in R for i in (1, 2)]
    val = {k: m.eval(M.t(k), model_completion=True).as_long() for k in keys}
    par = {k: k for k in keys}
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for p in P.interior:
        its = [e for _, e in P.H[p]] + ([((p, r), i) for r in res['choice'].get(p, ()) for i in (1, 2)] if scope != 'minimal' else [])
        for a, b in itertools.combinations(its, 2):
            if val[a] == val[b] and find(a) != find(b) and C.check(M.t(a) != M.t(b)) is None:
                par[find(a)] = find(b)
    folds = {}
    for k in keys: folds.setdefault(find(k), []).append(k)
    reps = list(folds)
    cand = {(x, y) for x in reps for y in reps if val[x] < val[y]}; lt = set()
    while cand:
        x, y = cand.pop()
        m2 = C.check(M.t(y) <= M.t(x))
        if m2 is None: lt.add((x, y))
        elif not isinstance(m2, str):
            v2 = {k: m2.eval(M.t(k), model_completion=True).as_long() for k in reps}
            cand = {(a, b) for a, b in cand if v2[a] < v2[b]}
    mp = m if scope != 'minimal' else res['m']           # the drawn segments are those of res
    paths = {(p, r): walk(P, M, mp, P.V[p], r, (p, r), p, None) for p, R in res['choice'].items() for r in R}
    res['analysis_unknown'] = C.unknown
    return folds, lt, val, paths


def segments(P, M, m, choice):
    """The cancelled segments of model m, as a set of edges between rounded points."""
    E = set()
    for p, R in choice.items():
        for r in R:
            pts, _ = walk(P, M, m, P.V[p], r, (p, r), p, None)
            for a, b in zip(pts, pts[1:]):
                E.add(frozenset(((round(a[0], 6), round(a[1], 6)), (round(b[0], 6), round(b[1], 6)))))
    return E


def geometries(P, res, cap=1000, budget=1800):
    """Minimal-cost geometries within the search space: dicts with the model m, the hidden
    directions (choice), a hashable selection, and the total length of the segments.  Stops after
    cap geometries or budget seconds; the second value says whether the list is complete."""
    M, opts = res['model'], res['opts']; s = M.s; out = []; t0 = time.monotonic(); r = None
    s.push(); s.add(M.cost() == res['cost'])
    while len(out) < cap:
        s.set('timeout', max(1, int(1000 * (budget - (time.monotonic() - t0)))))
        r = s.check(list(M.track))
        if r != z3.sat: break
        m = s.model()
        if refine(P, M, opts, m, [s]): continue
        on = [(p, i) for (p, i), b in M.sel.items() if z3.is_true(m.eval(b))]
        choice = {p: opts[p][i][0] for p, i in on if opts[p][i][0]}
        E = segments(P, M, m, choice)
        out.append({'m': m, 'choice': choice, 'length': sum(math.dist(*tuple(e)) for e in E),
                    'selection': frozenset((p, tuple(P.deg(r) for r in R)) for p, R in choice.items())})
        s.add(z3.Not(z3.And([M.sel[k] for k in on] + [M.geometry(choice, m)])))
    s.pop()
    return out, r == z3.unsat


def reduce(nodes, lt):
    succ = {x: {b for a, b in lt if a == x} for x in nodes}
    return {(a, b) for a, b in lt if not any(b in succ[c] for c in succ[a] if c != b)}


def report(P, folds, lt, paths, collapsed, out):
    import networkx as nx
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    V = P.V; reps = list(folds)
    where = {k: x for x, ks in folds.items() for k in ks}
    G = nx.DiGraph(); G.add_nodes_from(reps); G.add_edges_from(lt)
    depth = {}
    for x in nx.topological_sort(G): depth[x] = max([depth[a] + 1 for a in G.predecessors(x)], default=0)
    D = max(depth.values()) if depth else 0
    hidden_only = {x for x in reps if not any(isinstance(k, int) for k in folds[x])}
    def anchor(x):
        best = None
        for k in folds[x]:
            pts = [V[q] for q in P.crease[k][:2]] if isinstance(k, int) else paths[k[0]][0][:2]
            L = math.dist(*pts)
            if best is None or L > best[0]: best = (L, ((pts[0][0] + pts[1][0]) / 2, (pts[0][1] + pts[1][1]) / 2))
        return best[1]
    # left-right mirror symmetry of folds, if any
    xs = [q[0] for q in V]; X = min(xs) + max(xs)
    vkey = {(round(q[0], 6), round(q[1], 6)): i for i, q in enumerate(V)}
    vm = {i: vkey.get((round(X - q[0], 6) + 0.0, round(q[1], 6))) for i, q in enumerate(V)}
    seg = {frozenset(c[:2]): i for i, c in enumerate(P.crease)}
    hk = [k for k in where if not isinstance(k, int)]
    def mk(k):
        if isinstance(k, int):
            u, v, _ = P.crease[k]; return seg.get(frozenset((vm[u], vm[v])))
        (p, r), i = k
        for (q, s_), j in hk:
            if q == vm[p] and j == i and eq(s_, nrm(math.pi - r)): return ((q, s_), j)
    mirror = {}
    for x, ks in folds.items():
        im = {where.get(mk(k)) for k in ks}
        if len(im) == 1 and None not in im: mirror[x] = im.pop()
    symmetric = len(mirror) == len(reps) and all((mirror[a], mirror[b]) in lt for a, b in lt)
    order = sorted(reps, key=lambda x: (depth[x], round(anchor(x)[0], 3), anchor(x)[1]))
    label, n = {}, 0
    for x in order:
        if x in label: continue
        n += 1; label[x] = str(n); y = mirror.get(x) if symmetric else None
        if y is not None and y != x:
            if anchor(x)[0] > anchor(y)[0]: label[x], label[y] = str(n) + "'", str(n)
            else: label[y] = str(n) + "'"
    hasse = reduce(reps, lt)
    cmap = plt.get_cmap('viridis')
    col = {x: matplotlib.colors.to_hex(cmap(0.92 * depth[x] / max(D, 1))) for x in reps}
    txt = {x: ('white' if depth[x] < 0.45 * D else 'black') for x in reps}
    fig, ax = plt.subplots(figsize=(13, 13.6))
    for u, v in P.bound: ax.plot([V[u][0], V[v][0]], [V[u][1], V[v][1]], 'k-', lw=1)
    for i, (u, v, a) in enumerate(P.crease):
        ax.plot([V[u][0], V[v][0]], [V[u][1], V[v][1]], '-', color=col[where[i]], lw=2.4, solid_capstyle='round')
    drawn = set()
    for pts, how in paths.values():
        for a, b in zip(pts, pts[1:]):
            k = frozenset(((round(a[0], 4), round(a[1], 4)), (round(b[0], 4), round(b[1], 4))))
            if k not in drawn:
                drawn.add(k); ax.plot([a[0], b[0]], [a[1], b[1]], ls=(0, (4, 3)), color='crimson', lw=1.4)
    for p in collapsed:
        ax.plot(*V[p], 'o', ms=14, mfc='none', mec='darkorange', mew=2.2)
    for x in reps:
        cx, cy = anchor(x)
        ax.text(cx, cy, label[x], fontsize=8.5, ha='center', va='center', weight='bold', color=txt[x],
                bbox=dict(boxstyle='round,pad=0.12', fc=col[x], ec='crimson' if x in hidden_only else 'k', lw=0.8))
    ax.set_aspect('equal'); ax.axis('off')
    sm = plt.cm.ScalarMappable(cmap=matplotlib.colors.ListedColormap([cmap(0.92 * i / max(D, 1)) for i in range(D + 1)]),
                               norm=matplotlib.colors.Normalize(-0.5, D + 0.5))
    cb = plt.colorbar(sm, ax=ax, orientation='horizontal', fraction=0.03, pad=0.02, ticks=range(D + 1))
    cb.set_label('depth in the forced order')
    ax.set_title('Folds (creases forced to be made together), labelled as in the Hasse diagram.\n'
                 'Dashed crimson: hidden (cancelled) creases.  Orange rings: vertices made by a collapse.')
    plt.savefig(out + '_folds.png', dpi=110, bbox_inches='tight'); plt.close(fig)
    pos = {}
    for lev in range(D + 1):
        layer = sorted((x for x in reps if depth[x] == lev), key=lambda x: (anchor(x)[0], anchor(x)[1]))
        span = (max(xs) - min(xs)) or 1.0
        xl = [(anchor(x)[0] - min(xs)) / span * 16 for x in layer]
        for _ in range(500):
            moved = False
            for i in range(len(xl) - 1):
                if xl[i + 1] - xl[i] < 0.8:
                    d = (0.8 - (xl[i + 1] - xl[i])) / 2 + 1e-3; xl[i] -= d; xl[i + 1] += d; moved = True
            if not moved: break
        for x, px in zip(layer, xl): pos[x] = (px, -lev * 1.25)
    lines = ['digraph H {', 'splines=true; overlap=false; outputorder=edgesfirst;',
             'node [shape=circle, fixedsize=true, width=0.4, fontsize=10, style=filled, penwidth=0.8];',
             'edge [arrowsize=0.45, penwidth=0.7, color="#555555"];']
    for x in reps:
        lines.append('"%s" [pos="%.1f,%.1f!", fillcolor="%s", fontcolor="%s"%s];' % (
            label[x], 50 * pos[x][0], 72 * pos[x][1], col[x], txt[x], ', color="crimson", penwidth=1.6' if x in hidden_only else ''))
    for a, b in sorted(hasse, key=lambda e: (label[e[0]], label[e[1]])):
        lines.append('"%s" -> "%s";' % (label[a], label[b]))
    lines.append('}')
    open(out + '_hasse.dot', 'w').write('\n'.join(lines))
    if shutil.which('neato'):
        subprocess.run(['neato', '-n2', '-Tpng', '-Gdpi=100', out + '_hasse.dot', '-o', out + '_hasse.png'], check=True)
    else:                                   # no Graphviz: the same layout, with straight arrows
        fig, ax = plt.subplots(figsize=(max(6, 0.75 * 16), max(3, 1.1 * (D + 1))))
        for a, b in hasse:
            ax.annotate('', xy=pos[b], xytext=pos[a], arrowprops=dict(arrowstyle='-|>', color='#555555', lw=0.7,
                        shrinkA=9, shrinkB=9, mutation_scale=8))
        for x in reps:
            ax.text(*pos[x], label[x], ha='center', va='center', fontsize=8, color=txt[x], zorder=3,
                    bbox=dict(boxstyle='circle,pad=0.3', fc=col[x], ec='crimson' if x in hidden_only else 'k',
                              lw=1.6 if x in hidden_only else 0.8))
        ax.set_xlim(min(p[0] for p in pos.values()) - 0.6, max(p[0] for p in pos.values()) + 0.6)
        ax.set_ylim(-D * 1.25 - 0.6, 0.6); ax.axis('off')
        plt.savefig(out + '_hasse.png', dpi=100, bbox_inches='tight'); plt.close(fig)
    with open(out + '_folds.txt', 'w') as f:
        f.write('label depth  creases (id+assignment; hN@deg#1/2 = first/cancelling hidden crease from vertex N)\n')
        for x in order:
            f.write('%-4s %3d   %s\n' % (label[x], depth[x], ' '.join(
                ('%d%s' % (k, P.crease[k][2]) if isinstance(k, int) else 'h%d@%g#%d' % (k[0][0], P.deg(k[0][1]), k[1])) for k in folds[x])))
        f.write('\nHasse edges (earlier -> later):\n')
        for a, b in sorted(hasse, key=lambda e: (depth[e[0]], label[e[0]], label[e[1]])):
            f.write('%s -> %s\n' % (label[a], label[b]))
    return {'folds': len(reps), 'hidden_only': len(hidden_only), 'relations': len(lt), 'hasse_edges': len(hasse),
            'depth': D, 'mirror_symmetric': symmetric, 'isolated': sorted(label[x] for x in reps if G.degree(x) == 0)}


def main(argv):
    out, budget, tols, analyse_it, geometry, table, rigid_only = None, 300.0, (1e-6, 1e-4, 5e-4, 2e-3), True, 'shortest', None, False
    args = []
    for a in argv:
        if a.startswith('--geometry='): geometry = a[11:]; continue
        if a.startswith('--table='): table = a[8:]; continue
        if a.startswith('--rigid'): rigid_only = a[8:] or 'reverse'; continue
        if a.startswith('--out='): out = a[6:]
        elif a.startswith('--budget='): budget = float(a[9:])
        elif a.startswith('--tol='): tols = (float(a[6:]),)
        elif a == '--triage': analyse_it = False
        else: args.append(a)
    path = args[0]; t0 = time.monotonic()
    summary = {'file': path.split('/')[-1]}
    for tol in tols:
        P = Pattern(path, tol)
        notflat = [p for p in P.interior if maps(sorted(t for t, _ in P.H[p])) is None]
        if not notflat: break
    summary.update(tol=tol, N=P.N, creases=len(P.crease), interior=len(P.interior))
    if notflat:
        summary.update(status='not flat-foldable', vertices=notflat); print(json.dumps(summary)); return
    if table:                         # local histories from the lookup tables only (src/lookup.py)
        import lookup
        tl = time.monotonic(); T = lookup.Table(table, rigid_only); summary['table_load_s'] = round(time.monotonic() - tl, 2)
        opt = lookup.options_for(T, P)
    else: opt = P.options
    spent = [0.0]                     # time spent producing local histories, in either mode
    def timed(p, k, c=0):
        t1 = time.monotonic(); r = opt(p, k, c); spent[0] += time.monotonic() - t1; return r
    P.options = timed
    res = solve(P, log=lambda s: None if analyse_it is False else print(s, flush=True), budget=budget)
    summary['status'] = res['status']; summary['local_s'] = round(spent[0], 2)
    if table:
        summary['table'] = T.stats
        if res['status'] == 'local fail': summary['first_degree'] = len(P.H[res['first']])
    for k in ('hidden', 'excess_lines', 'collapses', 'cost', 'lower_bound', 'certified_minimal', 'rounds', 'lazy_refinements', 'expanded', 'core', 'first'):
        if k in res: summary[k] = res[k]
    if res['status'] == 'ok' and analyse_it:
        if geometry == 'shortest' and res['choice']:
            G, complete = geometries(P, res, budget=budget)
            summary.update(minimal_geometries=len(G), geometries_complete=complete)
            if complete:
                g = min(G, key=lambda g: g['length'])
                res = dict(res, m=g['m'], choice=g['choice'])
            else: summary['geometry'] = 'found (the minimal geometries were not all enumerated)'
        folds, lt, val, paths = analyse(P, res, budget=budget)
        summary['analysis_unknown'] = res['analysis_unknown']
        for (p, r), (pts, how) in sorted(paths.items(), key=lambda kv: (kv[0][0], kv[0][1])):
            print('  hidden from %d at %.2f deg: %s, %s' % (p, P.deg(r), ' -> '.join('(%.3f, %.3f)' % q for q in pts), how))
        summary.update(report(P, folds, lt, paths, res['collapses'], out or path.split('/')[-1].rsplit('.', 1)[0]))
    summary['t'] = round(time.monotonic() - t0, 1)
    print(json.dumps(summary))


if __name__ == '__main__':
    main(sys.argv[1:])
