"""Structure identification by reflection, for crease patterns in FOLD format.

Usage: python3 reflection_analysis.py [--tol=RADIANS] pattern.fold [more.fold ...]

Edges assigned M, V or U are creases, B edges are the boundary, and F edges are
taken to be reference lines and ignored. The steps follow the companion paper:
  1. find the coarsest unit 2*pi/N (N a power of two) making all directions integral;
  2. at each interior vertex, list the perfect pairings of half-edges (straight,
     d' = d + N/2, or reflected across another crease there, d' = 2c - d) whose
     local order is acyclic;
  3. a vertex with none is a defect. It may take one or two hidden pairs of
     cancelling creases. Each hidden pair runs straight to the edge of the sheet;
     every crease it crosses at an angle other than 90 degrees, or runs along,
     must have been folded after both of its chains;
  4. search the per-vertex options for combinations whose order is acyclic,
     using as few hidden pairs as possible, and report every such combination.
Vertices with no collinear pair of creases are listed separately: the first crease
made through a point is straight there, so such a vertex needs a cancelled crease
or a move made all at once.
"""
import json, math, sys, itertools

TOL = 1e-6        # radians
MAXK = 12         # try N = 4, 8, ..., 2**MAXK
MAX_HIDDEN = 2    # hidden pairs allowed at one vertex
MAX_SOLUTIONS = 10000


def load(path):
    d = json.load(open(path))
    V = [tuple(map(float, p[:2])) for p in d['vertices_coords']]
    A = d.get('edges_assignment') or ['U'] * len(d['edges_vertices'])
    crease, bound, ref = [], [], []
    for (u, v), a in zip(d['edges_vertices'], A):
        a = a.upper()
        (crease if a in 'MVU' else bound if a == 'B' else ref).append((u, v, a))
    return V, crease, bound, ref


def find_unit(angles):
    for k in range(2, MAXK + 1):
        N = 2 ** k
        if all(abs(t * N / (2 * math.pi) - round(t * N / (2 * math.pi))) * 2 * math.pi / N < TOL for t in angles):
            return N
    return None


def acyclic(edges):
    succ, indeg = {}, {}
    for a, b in edges:
        succ.setdefault(a, set()).add(b); succ.setdefault(b, set())
    for a in succ:
        indeg.setdefault(a, 0)
        for b in succ[a]: indeg[b] = indeg.get(b, 0) + 1
    stack = [x for x in succ if indeg[x] == 0]; seen = 0
    while stack:
        x = stack.pop(); seen += 1
        for y in succ[x]:
            indeg[y] -= 1
            if indeg[y] == 0: stack.append(y)
    return seen == len(succ)


def hidden_ray(V, crease, bound, p, theta):
    """Trace a hidden segment from vertex p straight to the edge.
    Returns (end point, ids of creases it crosses obliquely or runs along) or None."""
    o = V[p]; d = (math.cos(theta), math.sin(theta)); tb = math.inf
    for u, v, _ in bound:
        a, b = V[u], V[v]; e = (b[0] - a[0], b[1] - a[1]); den = d[0] * e[1] - d[1] * e[0]
        if abs(den) < 1e-12: continue
        w = (a[0] - o[0], a[1] - o[1])
        t = (w[0] * e[1] - w[1] * e[0]) / den; s = (w[0] * d[1] - w[1] * d[0]) / den
        if t > 1e-9 and -1e-9 <= s <= 1 + 1e-9: tb = min(tb, t)
    if tb == math.inf: return None
    eps = 1e-6 * max(1.0, tb); crossed = set()
    for i, (u, v, _) in enumerate(crease):
        if p in (u, v): continue
        a, b = V[u], V[v]; e = (b[0] - a[0], b[1] - a[1]); L = math.hypot(*e)
        w = (a[0] - o[0], a[1] - o[1]); den = d[0] * e[1] - d[1] * e[0]
        if abs(den) < 1e-9 * L:                       # parallel: does it run along the segment?
            if abs(w[0] * d[1] - w[1] * d[0]) < 1e-6 * max(1.0, L):
                ta = w[0] * d[0] + w[1] * d[1]; tb2 = (b[0] - o[0]) * d[0] + (b[1] - o[1]) * d[1]
                if max(ta, tb2) > eps and min(ta, tb2) < tb - eps: crossed.add(i)
            continue
        t = (w[0] * e[1] - w[1] * e[0]) / den; s = (w[0] * d[1] - w[1] * d[0]) / den
        if eps < t < tb - eps and -1e-9 <= s <= 1 + 1e-9:
            if abs(d[0] * e[0] + d[1] * e[1]) > 1e-9 * L:   # not perpendicular
                crossed.add(i)
    return (o[0] + tb * d[0], o[1] + tb * d[1]), crossed


def analyze(path):
    V, crease, bound, ref = load(path)
    out = {'file': path.split('/')[-1], 'creases': len(crease), 'reference lines': len(ref)}
    H = {}
    for i, (u, v, _) in enumerate(crease):
        for p, q in ((u, v), (v, u)):
            H.setdefault(p, []).append((math.atan2(V[q][1] - V[p][1], V[q][0] - V[p][0]) % (2 * math.pi), i))
    N = find_unit([t for hs in H.values() for t, _ in hs])
    out['N'] = N
    if N is None:
        out['directions not dyadic (deg mod 180)'] = sorted({round(math.degrees(t) % 180, 3) for hs in H.values() for t, _ in hs
            if abs(t * 2 ** MAXK / (2 * math.pi) - round(t * 2 ** MAXK / (2 * math.pi))) > 1e-3})
        return out
    half = N // 2
    D = {p: [(round(t * N / (2 * math.pi)) % N, e) for t, e in sorted(hs)] for p, hs in H.items()}
    on_boundary = {x for u, v, _ in bound for x in (u, v)}
    interior = [p for p in D if p not in on_boundary]
    xy = lambda p: tuple(round(c, 3) for c in V[p])
    straight = lambda a, b: (a + half - b) % N == 0
    mirrors = lambda d, a, b, skip: [g for g in range(len(d)) if g not in skip and (2 * d[g] - a - b) % N == 0]
    admissible = lambda d, i, j: straight(d[i], d[j]) or (d[i] != d[j] and bool(mirrors(d, d[i], d[j], (i, j))))
    out['interior vertices'] = len(interior)
    out['no collinear pair'] = [xy(p) for p in interior
        if not any(straight(a, b) for (a, _), (b, _) in itertools.combinations(D[p], 2))]

    def perfect(d, n_real):
        """Perfect pairings; hidden half-edges (index >= n_real) must pair with real ones."""
        res = []
        def rec(rem, cur):
            if not rem: res.append(cur); return
            i = rem[0]
            for j in rem[1:]:
                if i >= n_real and j >= n_real: continue
                if admissible(d, i, j): rec([x for x in rem if x not in (i, j)], cur + [(i, j)])
        rec(list(range(len(d))), []); return res

    def relations(d, pairs, label):
        return [(label[i], {label[g] for g in mirrors(d, d[i], d[j], (i, j))} - {label[i]})
                for i, j in pairs if not straight(d[i], d[j])]

    def some_acyclic(rel, fixed=()):
        if not all(s for _, s in rel): return False
        return any(acyclic(list(fixed) + [(w, l) for (l, _), w in zip(rel, c)])
                   for c in itertools.product(*[sorted(s) for _, s in rel]))

    def local_ok(d, pairs):
        label = {}
        for n, (i, j) in enumerate(pairs): label[i] = label[j] = n
        return some_acyclic(relations(d, pairs, label))

    # ---- local stage: options per vertex.  An option = (hidden list, pairs), where a
    # hidden entry is (direction, id1, id2, end point, crossed crease ids).
    next_id = [len(crease)]
    rays = {}
    def ray(p, r):
        if (p, r) not in rays: rays[(p, r)] = hidden_ray(V, crease, bound, p, 2 * math.pi * r / N)
        return rays[(p, r)]

    def hidden_options(p, k):
        d = [x for x, _ in D[p]]; n = len(d); res = []
        free = [r for r in range(N) if r not in d and ray(p, r) is not None]
        for rs in itertools.combinations(free, k):
            d2 = d + [r for r in rs for _ in (0, 1)]
            for s in perfect(d2, n):
                mate = {i: j for pr in s for i, j in (pr, pr[::-1])}
                if any(mate[n + 2 * t] > mate[n + 2 * t + 1] for t in range(k)): continue   # copies interchangeable
                if not local_ok(d2, s): continue
                hid = []
                for r in rs:
                    end, crossed = ray(p, r)
                    hid.append((r, next_id[0], next_id[0] + 1, end, crossed)); next_id[0] += 2
                res.append((hid, s))
        return res

    fixed_choice, variable = {}, {}
    status = {}
    for p in interior:
        d = [x for x, _ in D[p]]
        o = [([], s) for s in perfect(d, len(d)) if local_ok(d, s)]
        if len(o) == 1: fixed_choice[p] = o[0]; status[p] = 'forced'
        elif o: variable[p] = {0: o}; status[p] = 'ambiguous'
        else: variable[p] = {}; status[p] = 'defect'
    defects = [p for p in interior if status[p] == 'defect']
    out['ambiguous vertices (xy, options)'] = [(xy(p), len(variable[p][0])) for p in interior if status[p] == 'ambiguous']

    # ---- global stage: depth-first search with an incremental cycle check,
    # by increasing total number of hidden pairs.
    def edges_of(p, hid):
        return [e for _, e in D[p]] + [x for h in hid for x in (h[1], h[2])]
    def dirs_of(p, hid):
        return [x for x, _ in D[p]] + [h[0] for h in hid for _ in (0, 1)]

    def order_of(assign):
        """Union-find over paired creases, then (definite relations, multi-witness relations)."""
        parent = {}
        def find(x):
            parent.setdefault(x, x)
            while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
            return x
        for p, (hid, pairs) in assign.items():
            E = edges_of(p, hid)
            for i, j in pairs: parent[find(E[i])] = find(E[j])
        fixed, multi = set(), []
        for p, (hid, pairs) in assign.items():
            E, d = edges_of(p, hid), dirs_of(p, hid)
            for l, s in relations(d, pairs, {i: find(E[i]) for i in range(len(d))}):
                if not s: return None, None, None
                if len(s) == 1: fixed.add((next(iter(s)), l))
                else: multi.append((l, s))
            for r, h1, h2, end, crossed in hid:          # crossed creases are later than both hidden chains
                for z in crossed:
                    fixed.add((find(h1), find(z))); fixed.add((find(h2), find(z)))
        return fixed, multi, find

    def full_ok(assign):
        fixed, multi, find = order_of(assign)
        if fixed is None: return None
        for c in itertools.product(*[sorted(s) for _, s in multi]):
            order = fixed | {(w, l) for (l, _), w in zip(multi, c)}
            if acyclic(list(order)): return order, find
        return None

    def options_with_budget(p, kmax):
        ks = [0] if status[p] == 'ambiguous' else range(1, kmax + 1)
        for k in ks:
            if k not in variable[p]:
                variable[p][k] = hidden_options(p, k)
        return [(k, o) for k in ks for o in variable[p].get(k, [])]

    order_vars = sorted(variable, key=lambda p: (status[p] != 'defect', p))
    solutions, used_budget = [], None
    if not defects and not variable:
        res = full_ok(fixed_choice)
        if res: solutions.append((dict(fixed_choice), res)); used_budget = 0
    for budget in range(len(defects), MAX_HIDDEN * len(defects) + 1):
        kmax = min(MAX_HIDDEN, budget - len(defects) + 1)
        opts = {p: options_with_budget(p, kmax) for p in variable}
        if any(not opts[p] for p in defects): break
        found = []
        def dfs(i, assign, spent):
            if len(found) >= MAX_SOLUTIONS: return
            if i == len(order_vars):
                if spent != budget: return
                res = full_ok(assign)
                if res: found.append((dict(assign), res))
                return
            p = order_vars[i]
            rest_min = sum(1 for q in order_vars[i + 1:] if status[q] == 'defect')
            for k, o in opts[p]:
                if spent + k + rest_min > budget: continue
                assign[p] = o
                fixed, multi, _ = order_of(assign)
                if fixed is not None and acyclic(list(fixed)): dfs(i + 1, assign, spent + k)
                del assign[p]
        dfs(0, dict(fixed_choice), 0)
        if found: solutions, used_budget = found, budget; break
    out['defects'] = [xy(p) for p in defects]
    out['hidden pairs needed'] = used_budget
    out['resolutions'] = len(solutions)
    if not solutions: return out

    def describe(assign):
        rows = []
        for p in interior:
            hid, pairs = assign[p]
            if status[p] == 'forced': continue
            d = dirs_of(p, hid)
            rows.append((xy(p), status[p], [(h[0], tuple(round(c, 1) for c in h[3]), len(h[4])) for h in hid],
                         sorted((d[i], d[j]) for i, j in pairs)))
        return rows
    if len(solutions) == 1:
        assign, (order, find) = solutions[0]
        used = {find(e) for p, (hid, _) in assign.items() for e in edges_of(p, hid)} | {find(e) for e in range(len(crease))}
        out['chains'] = len(used)
        out['order relations'] = len(order)
        idx = {c: k for k, c in enumerate(sorted(used))}
        if len(used) <= 24:
            pred = [0] * len(used)
            for a, b in order:
                if a in idx and b in idx: pred[idx[b]] |= 1 << idx[a]
            cnt = [0] * (1 << len(used)); cnt[0] = 1
            for m in range(1 << len(used)):
                if cnt[m]:
                    for v in range(len(used)):
                        if not m >> v & 1 and pred[v] & m == pred[v]: cnt[m | 1 << v] += cnt[m]
            out['linear extensions'] = cnt[-1]
        out['resolution (xy, status, hidden: direction/end/crossings, pairs)'] = describe(assign)
    else:
        varies = {}
        for assign, _ in solutions:
            for row in describe(assign): varies.setdefault(row[0], set()).add(repr(row[2:]))
        out['vertices that differ between resolutions (xy, distinct choices)'] = [(k, len(v)) for k, v in varies.items() if len(v) > 1]
        out['vertices fixed across resolutions'] = [k for k, v in varies.items() if len(v) == 1]
    return out


if __name__ == '__main__':
    args = sys.argv[1:]
    if args and args[0].startswith('--tol='): TOL = float(args.pop(0)[6:])
    for path in args:
        r = analyze(path)
        print('\n' + r.pop('file'))
        for k, v in r.items(): print('  %-44s %s' % (k, v))
