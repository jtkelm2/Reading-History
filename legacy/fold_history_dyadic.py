"""Structure identification by local fold histories, for crease patterns in FOLD format.

Usage: python3 fold_history.py [--out=PREFIX] [--pick=K] [--tol=RADIANS] pattern.fold

Model (single-line folds, shadow model of the companion paper):
  * Directions are integers mod N, measured from the first crease (rotation invariant).
  * At each interior vertex v, a local history is a sequence of steps.  A step creates a set T of
    half-edges at once (one fold) and is valid iff
      (K) the creases already present at v are flat-foldable there (product of reflections = id);
      (L) every half-edge of T is a preimage of one fold line under the current folding maps
          of the sectors at v (Lemma 3 with any number of reflections; replaces Assumption 8);
      (P) |T| is even (paper switches between moving and staying only at T; Lemma 5);
      (C) no two uncut sectors straddle the same ray of the fold line with one moving and one
          staying (the collision argument of Prop. 20).
    Every crease at v must be created by exactly one step.  Pairings are not chosen: creases
    made in one step are simultaneous, which is all the global order needs.
  * A vertex with no history is a defect.  It may carry one or two hidden pairs: two coincident
    creases, the second cancelling the first (Lemma 4).  A hidden crease runs straight until it
    reaches the boundary, joins the opposite hidden ray of another vertex, or bends across an
    earlier crease (at most MAXBEND bends, every continuation traced).  A crease it crosses
    obliquely without bending, or runs along, is later than both hidden creases.
  * Global: integer times t(e) per crease; each vertex takes one (hidden set, history); the
    selected histories impose equalities within steps and strict order between steps.
  * Output: every feasible hidden selection; for the chosen one, the folds (creases forced to
    be simultaneous in every solution), the forced order between folds, and its Hasse diagram.
"""
import json, math, sys, itertools, subprocess
import z3

TOL = 1e-6
MAXK = 12
MAX_HIDDEN = 2
MAXBEND = 2


# ---------------------------------------------------------------- geometry
def load(path):
    d = json.load(open(path))
    V = [tuple(map(float, p[:2])) for p in d['vertices_coords']]
    A = [a.upper() for a in (d.get('edges_assignment') or ['U'] * len(d['edges_vertices']))]
    crease = [(u, v, a) for (u, v), a in zip(d['edges_vertices'], A) if a in 'MVU']
    bound = [(u, v) for (u, v), a in zip(d['edges_vertices'], A) if a == 'B']
    return V, crease, bound


def find_unit(rel):
    for k in range(2, MAXK + 1):
        N = 2 ** k
        if all(abs(t * N / (2 * math.pi) - round(t * N / (2 * math.pi))) * 2 * math.pi / N < TOL for t in rel):
            return N
    return None


class Pattern:
    def __init__(self, path):
        self.V, self.crease, self.bound = load(path)
        V = self.V
        ang = {}
        for i, (u, v, _) in enumerate(self.crease):
            for p, q in ((u, v), (v, u)):
                ang.setdefault(p, []).append((math.atan2(V[q][1] - V[p][1], V[q][0] - V[p][0]), i))
        self.theta0 = ang[self.crease[0][0]][0][0]
        rel = [(t - self.theta0) % (2 * math.pi) for hs in ang.values() for t, _ in hs]
        self.N = N = find_unit(rel)
        if N is None: raise SystemExit('directions are not dyadic multiples of a common unit')
        self.H = {p: sorted((round(((t - self.theta0) % (2 * math.pi)) * N / (2 * math.pi)) % N, e) for t, e in hs)
                  for p, hs in ang.items()}
        onb = {x for e in self.bound for x in e}
        self.interior = sorted(p for p in self.H if p not in onb)
        self.cdir = [self.H[u][[e for _, e in self.H[u]].index(i)][0] for i, (u, v, _) in enumerate(self.crease)]

    def deg(self, r):
        return round((math.degrees(self.theta0) + r * 360 / self.N) % 360, 4)

    def unit(self, r):
        th = self.theta0 + 2 * math.pi * r / self.N
        return math.cos(th), math.sin(th)

    def trace(self, o, r, skip_v=None, skip_c=None):
        """Events along the ray from point o in direction r, nearest first, and the boundary end:
        ('x', crease, perpendicular, point, crease direction)  interior crossing,
        ('v', vertex, [(crease, perpendicular)])             passage through a vertex,
        ('a', crease)                                         running along a crease."""
        V, N = self.V, self.N; dv = self.unit(r)
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
        ev = {}
        for w in range(len(V)):
            if w == skip_v: continue
            q = (V[w][0] - o[0], V[w][1] - o[1]); t = q[0] * dv[0] + q[1] * dv[1]
            if 1e-7 < t < tb - 1e-7 and abs(q[0] * dv[1] - q[1] * dv[0]) < 1e-7:
                ev[round(t, 9)] = ('v', w, [(e, (dd - r) % (N // 2) == N // 4) for dd, e in self.H.get(w, [])])
        for i, (u, v, _) in enumerate(self.crease):
            if i == skip_c or (skip_v is not None and skip_v in (u, v)): continue
            h = hit(V[u], V[v])
            if h is None:
                q = (V[u][0] - o[0], V[u][1] - o[1])
                if abs(q[0] * dv[1] - q[1] * dv[0]) < 1e-7:
                    ts = sorted([q[0] * dv[0] + q[1] * dv[1], (V[v][0] - o[0]) * dv[0] + (V[v][1] - o[1]) * dv[1]])
                    if ts[1] > 1e-7 and ts[0] < tb - 1e-7: ev[round(max(ts[0], 0) + 1e-6, 9)] = ('a', i)
                continue
            t, s = h
            if 1e-7 < t < tb - 1e-7 and 1e-7 < s < 1 - 1e-7:
                ev[round(t, 9)] = ('x', i, (self.cdir[i] - r) % (N // 2) == N // 4,
                                   (o[0] + t * dv[0], o[1] + t * dv[1]), self.cdir[i])
        return [ev[k] for k in sorted(ev)], (o[0] + tb * dv[0], o[1] + tb * dv[1])


# ---------------------------------------------------------------- local histories
def maps(S, N):
    """Folding map of each sector (S[j], S[j+1]) as x -> s*x + o, or None if not flat-foldable."""
    if not S: return [(1, 0)]
    out = [(1, 0)]
    for c in S[1:]:
        s, o = out[-1]; out.append((-s, (2 * s * c + o) % N))
    s, o = out[-1]
    return out if (-s, (2 * s * S[0] + o) % N) == (1, 0) else None


def between(x, a, b, N):
    return x != a if a == b else 0 < (x - a) % N < (b - a) % N


def sector_of(x, S, N):
    if not S: return 0
    m = len(S)
    for j in range(m):
        if x == S[j] or between(x, S[j], S[(j + 1) % m], N): return j


def collision_free(S, switch, phis, N, line):
    R = sorted(set(S) | set(switch)); m = len(R)
    mv, cur = [], 0
    for k in range(m):
        mv.append(cur)
        if R[(k + 1) % m] in switch: cur ^= 1
    if cur: return False
    for rho in (line, line + N // 2):
        seen = set()
        for k in range(m):
            a, b = R[k], R[(k + 1) % m]
            s, o = phis[sector_of(a, S, N)]
            if between((s * (rho - o)) % N, a, b, N): seen.add(mv[k])
        if len(seen) > 1: return False
    return True


def histories(items, N):
    """items: (direction, key, kind), kind in real/hid1/hid2; hid2 key (h, 2) cancels hid1 key (h, 1).
    Returns the valid histories, each a tuple of steps (frozensets of keys)."""
    by = {k: (d, kind) for d, k, kind in items}; res = []
    def rec(S, rem, hist):
        if not rem: res.append(tuple(hist)); return
        phis = maps(S, N)
        if phis is None: return
        cand = {}
        for k in rem:
            d, kind = by[k]
            if kind == 'hid2' and (k[0], 1) in rem: continue
            s, o = phis[sector_of(d, S, N)]
            cand.setdefault((s * d + o) % (N // 2), []).append(k)
        for line, ks in cand.items():
            for size in range(2, len(ks) + 1, 2):
                for T in itertools.combinations(ks, size):
                    dirs = [by[k][0] for k in T]
                    if len(set(dirs)) < size or not collision_free(S, dirs, phis, N, line): continue
                    S2 = sorted((set(S) | {by[k][0] for k in T if by[k][1] != 'hid2'})
                                - {by[k][0] for k in T if by[k][1] == 'hid2'})
                    hist.append(frozenset(T)); rec(S2, rem - set(T), hist); hist.pop()
    rec([], frozenset(k for _, k, _ in items), [])
    return res


def local_options(P, p):
    base = [(d, e, 'real') for d, e in P.H[p]]
    hs = histories(base, P.N)
    if hs: return [((), hs)]
    dirs = {d for d, _ in P.H[p]}
    for k in range(1, MAX_HIDDEN + 1):
        out = []
        for R in itertools.combinations([r for r in range(P.N) if r not in dirs], k):
            items = base + [x for r in R for x in ((r, ((p, r), 1), 'hid1'), (r, ((p, r), 2), 'hid2'))]
            hs = histories(items, P.N)
            if hs: out.append((R, hs))
        if out: return out
    return []


# ---------------------------------------------------------------- global model
class Model:
    def __init__(self, P, opts):
        self.P, self.opts = P, opts
        self.s = z3.Solver(); self.T = {}
        self.sel = {(p, i): z3.Bool('s_%d_%d' % (p, i)) for p, ol in opts.items() for i in range(len(ol))}
        for p, ol in opts.items():
            for i, (R, hs) in enumerate(ol):
                b = self.sel[p, i]
                self.s.add(z3.Implies(b, z3.Or([self.history(h) for h in hs])))
                for r in R:
                    h1, h2 = self.t(((p, r), 1)), self.t(((p, r), 2))
                    self.s.add(z3.Implies(b, self.path(P.V[p], r, 0, h1, h2, p, None)))
            self.s.add(z3.PbEq([(self.sel[p, i], 1) for i in range(len(ol))], 1))

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
        return [z3.And(self.sel[w, i], h1 == self.t(((w, r2), 1)), h2 == self.t(((w, r2), 2)))
                for i, (Rw, _) in enumerate(self.opts.get(w, [])) if r2 in Rw]

    def path(self, o, r, depth, h1, h2, skip_v, skip_c):
        P, N = self.P, self.P.N
        tr = P.trace(o, r, skip_v, skip_c)
        if tr is None: return z3.BoolVal(False)
        terms, pre = [], []
        for e in tr[0]:
            if e[0] == 'x':
                if e[2]: continue                         # right angle: no constraint
                z = self.t(e[1])
                if depth < MAXBEND:
                    terms.append(z3.And(pre + [z < h1, self.path(e[3], (2 * e[4] - r - N // 2) % N,
                                                                  depth + 1, h1, h2, None, e[1])]))
                pre = pre + [z > h2]
            elif e[0] == 'v':
                terms += [z3.And(pre + [j]) for j in self.joins(e[1], (r + N // 2) % N, h1, h2)]
                pre = pre + [self.t(z) > h2 for z, perp in e[2] if not perp]
            else:
                pre = pre + [self.t(e[1]) > h2]
        return z3.Or(terms + [z3.And(pre)])

    def selections(self, limit=1000):
        s = self.s; out = []; s.push()
        while len(out) < limit and s.check() == z3.sat:
            m = s.model()
            ch = {p: self.opts[p][i][0] for (p, i), b in self.sel.items()
                  if self.opts[p][i][0] and z3.is_true(m.eval(b))}
            out.append(ch)
            s.add(z3.Or([z3.Not(self.sel[p, i]) for (p, i) in self.sel if self.opts[p][i][0] == ch.get(p)]))
        s.pop(); return out

    def fix(self, ch):
        for (p, i), b in self.sel.items(): self.s.add(b == (self.opts[p][i][0] == ch.get(p, ())))

    def value(self, m, k): return m.eval(self.t(k), model_completion=True).as_long()

    def walk(self, m, o, r, depth, k1, skip_v, skip_c):
        """The hidden polyline actually used in model m, and how it ends."""
        P, N = self.P, self.P.N
        h1, h2 = self.value(m, (k1, 1)), self.value(m, (k1, 2))
        ev, end = P.trace(o, r, skip_v, skip_c)
        for e in ev:
            if e[0] == 'x' and not e[2] and self.value(m, e[1]) < h1:
                pts, how = self.walk(m, e[3], (2 * e[4] - r - N // 2) % N, depth + 1, k1, None, e[1])
                return [o] + pts, how
            if e[0] == 'v':
                r2 = (r + N // 2) % N; w = e[1]
                if (w, r2) != k1 and ((w, r2), 1) in self.T and any(z3.is_true(m.eval(j)) for j in self.joins(w, r2, self.t((k1, 1)), self.t((k1, 2)))):
                    return [o, P.V[w]], 'joins vertex %d' % w
        return [o, end], 'boundary'


# ---------------------------------------------------------------- analysis of one resolution
def analyse(P, opts, ch):
    M = Model(P, opts); M.fix(ch); s = M.s
    assert s.check() == z3.sat
    m = s.model()
    keys = list(range(len(P.crease))) + [((p, r), i) for p, R in ch.items() for r in R for i in (1, 2)]
    val = {k: M.value(m, k) for k in keys}
    # folds: creases forced to be simultaneous
    par = {k: k for k in keys}
    def find(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for p in P.interior:
        its = [e for _, e in P.H[p]] + [((p, r), i) for r in ch.get(p, ()) for i in (1, 2)]
        for a, b in itertools.combinations(its, 2):
            if val[a] == val[b] and find(a) != find(b) and s.check(M.t(a) != M.t(b)) == z3.unsat:
                par[find(a)] = find(b)
    folds = {}
    for k in keys: folds.setdefault(find(k), []).append(k)
    reps = list(folds)
    # forced order between folds (each counter-model discards every relation it violates)
    cand = {(x, y) for x in reps for y in reps if val[x] < val[y]}; lt = set()
    while cand:
        x, y = cand.pop()
        if s.check(M.t(y) <= M.t(x)) == z3.unsat: lt.add((x, y))
        else:
            m2 = s.model(); v2 = {k: M.value(m2, k) for k in reps}
            cand = {(a, b) for a, b in cand if v2[a] < v2[b]}
    paths = {(p, r): M.walk(m, P.V[p], r, 0, (p, r), p, None) for p, R in ch.items() for r in R}
    return folds, lt, val, paths


def reduce(nodes, lt):
    succ = {x: {b for a, b in lt if a == x} for x in nodes}
    return {(a, b) for a, b in lt if not any(b in succ[c] for c in succ[a] if c != b)}


# ---------------------------------------------------------------- output
def main(path, out='out', pick=None):
    P = Pattern(path); N = P.N
    print('%s: %d creases, %d interior vertices, N = %d (unit %.4g deg)' % (path, len(P.crease), len(P.interior), N, 360 / N))
    opts = {p: local_options(P, p) for p in P.interior}
    defects = [p for p in P.interior if opts[p] and opts[p][0][0]]
    nohope = [p for p in P.interior if not opts[p]]
    nh = {p: sum(len(h) for _, h in opts[p]) for p in P.interior}
    print('local histories per vertex: min %d, max %d; defects %d; unresolvable %d'
          % (min(nh.values()), max(nh.values()), len(defects), len(nohope)))
    for p in defects:
        print('  defect %3d at (%.3f, %.3f): hidden options %s' % (p, *P.V[p], [tuple(P.deg(r) for r in R) for R, _ in opts[p]]))
    if nohope: print('no resolution: vertices', nohope); return
    M = Model(P, opts)
    sels = M.selections()
    print('feasible hidden selections:', len(sels))
    if not sels: return
    common = {p: R for p, R in sels[0].items() if all(c.get(p) == R for c in sels)}
    print('  common to all:', {p: tuple(P.deg(r) for r in R) for p, R in common.items()})
    for i, c in enumerate(sels):
        print('  %d: differs at' % i, {p: tuple(P.deg(r) for r in R) for p, R in c.items() if p not in common})
    # choose the selection with least total hidden length unless told otherwise
    res = []
    for i, c in enumerate(sels if pick is None else [sels[pick]]):
        folds, lt, val, paths = analyse(P, opts, c)
        L = sum(math.dist(a, b) for pts, _ in paths.values() for a, b in zip(pts, pts[1:]))
        res.append((L, i if pick is None else pick, c, folds, lt, val, paths))
        print('  selection %d: hidden length %.4f, %d folds, %d forced relations' % (res[-1][1], L, len(folds), len(lt)))
    L, i, c, folds, lt, val, paths = min(res, key=lambda x: x[0])
    print('analysed selection', i)
    for (p, r), (pts, how) in sorted(paths.items()):
        print('  hidden from %d at %.1f deg: %s, %s' % (p, P.deg(r), ' -> '.join('(%.3f, %.3f)' % q for q in pts), how))
    report(P, folds, lt, val, paths, out)


def report(P, folds, lt, val, paths, out):
    import networkx as nx
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    V, N = P.V, P.N; reps = list(folds)
    where = {k: x for x, ks in folds.items() for k in ks}
    G = nx.DiGraph(); G.add_nodes_from(reps); G.add_edges_from(lt)
    depth = {}
    for x in nx.topological_sort(G): depth[x] = max([depth[a] + 1 for a in G.predecessors(x)], default=0)
    D = max(depth.values())
    hidden_only = {x for x in reps if not any(isinstance(k, int) for k in folds[x])}
    # where to put a label: midpoint of the longest crease (or hidden segment) of the fold
    def anchor(x):
        best = None
        for k in folds[x]:
            pts = [V[q] for q in P.crease[k][:2]] if isinstance(k, int) else paths[k[0]][0][:2]
            L = math.dist(*pts)
            if best is None or L > best[0]: best = (L, ((pts[0][0] + pts[1][0]) / 2, (pts[0][1] + pts[1][1]) / 2))
        return best[1]
    # left-right mirror map on folds (vertices, creases and hidden rays), if the pattern has one
    xs = [q[0] for q in V]; X = min(xs) + max(xs)
    vkey = {(round(q[0], 6), round(q[1], 6)): i for i, q in enumerate(V)}
    vm = {i: vkey.get((round(X - q[0], 6) + 0.0, round(q[1], 6))) for i, q in enumerate(V)}
    seg = {frozenset(c[:2]): i for i, c in enumerate(P.crease)}
    def mk(k):
        if isinstance(k, int):
            u, v, _ = P.crease[k]; return seg.get(frozenset((vm[u], vm[v])))
        (p, r), i = k
        th = math.pi - math.radians(P.deg(r))
        return ((vm[p], round(((th - P.theta0) % (2 * math.pi)) * N / (2 * math.pi)) % N), i)
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
    # ---- crease pattern coloured by depth
    fig, ax = plt.subplots(figsize=(13, 13.6))
    for u, v in P.bound: ax.plot([V[u][0], V[v][0]], [V[u][1], V[v][1]], 'k-', lw=1)
    for i, (u, v, a) in enumerate(P.crease):
        ax.plot([V[u][0], V[v][0]], [V[u][1], V[v][1]], '-', color=col[where[i]], lw=2.4, solid_capstyle='round')
    drawn = set()
    for pts, how in paths.values():
        for a, b in zip(pts, pts[1:]):
            k = frozenset(((round(a[0], 4), round(a[1], 4)), (round(b[0], 4), round(b[1], 4))))
            if k in drawn: continue
            drawn.add(k); ax.plot([a[0], b[0]], [a[1], b[1]], ls=(0, (4, 3)), color='crimson', lw=1.4)
    for x in reps:
        cx, cy = anchor(x)
        ax.text(cx, cy, label[x], fontsize=8.5, ha='center', va='center', weight='bold', color=txt[x],
                bbox=dict(boxstyle='round,pad=0.12', fc=col[x], ec='crimson' if x in hidden_only else 'k', lw=0.8))
    ax.set_aspect('equal'); ax.axis('off')
    sm = plt.cm.ScalarMappable(cmap=matplotlib.colors.ListedColormap([cmap(0.92 * i / max(D, 1)) for i in range(D + 1)]),
                               norm=matplotlib.colors.Normalize(-0.5, D + 0.5))
    cb = plt.colorbar(sm, ax=ax, orientation='horizontal', fraction=0.03, pad=0.02, ticks=range(D + 1))
    cb.set_label('depth in the forced order (longest chain of forced predecessors)')
    ax.set_title('Folds: creases forced to be made together, labelled as in the Hasse diagram.\n'
                 'Dashed crimson: hidden (cancelled) creases.  A prime marks the mirror image.')
    plt.savefig(out + '_folds.png', dpi=110, bbox_inches='tight'); plt.close(fig)
    # ---- Hasse diagram: y = depth, x = position of the fold on the sheet (spread within a layer)
    pos = {}
    for lev in range(D + 1):
        layer = sorted((x for x in reps if depth[x] == lev), key=lambda x: (anchor(x)[0], anchor(x)[1]))
        xsl = [anchor(x)[0] * 16 for x in layer]
        for _ in range(200):                        # push apart, symmetric about the sheet's centre
            moved = False
            for i in range(len(xsl) - 1):
                if xsl[i + 1] - xsl[i] < 0.8:
                    d = (0.8 - (xsl[i + 1] - xsl[i])) / 2 + 1e-3; xsl[i] -= d; xsl[i + 1] += d; moved = True
            if not moved: break
        for x, px in zip(layer, xsl): pos[x] = (px, -lev * 1.25)
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
    subprocess.run(['neato', '-n2', '-Tpng', '-Gdpi=100', out + '_hasse.dot', '-o', out + '_hasse.png'], check=True)
    iso = sorted(label[x] for x in reps if G.degree(x) == 0)
    print('%d folds (%d entirely hidden), %d forced relations, %d Hasse edges, depth %d; mirror-symmetric order: %s'
          % (len(reps), len(hidden_only), len(lt), len(hasse), D, symmetric))
    print('incomparable with every other fold:', iso)
    with open(out + '_folds.txt', 'w') as f:
        f.write('label depth  creases (id+assignment; hN@deg#1/2 = first/cancelling hidden crease from vertex N)\n')
        for x in order:
            f.write('%-4s %3d   %s\n' % (label[x], depth[x], ' '.join(
                ('%d%s' % (k, P.crease[k][2]) if isinstance(k, int) else 'h%d@%g#%d' % (k[0][0], P.deg(k[0][1]), k[1])) for k in folds[x])))
        f.write('\nHasse edges (earlier -> later):\n')
        for a, b in sorted(hasse, key=lambda e: (depth[e[0]], label[e[0]], label[e[1]])):
            f.write('%s -> %s\n' % (label[a], label[b]))


if __name__ == '__main__':
    args = sys.argv[1:]; out = 'out'; pick = None
    while args and args[0].startswith('--'):
        a = args.pop(0)
        if a.startswith('--out='): out = a[6:]
        elif a.startswith('--pick='): pick = int(a[7:])
        elif a.startswith('--tol='): TOL = float(a[6:])
    main(args[0], out, pick)
