"""Export the degree-4 lookup tables as an encyclopedia of strata (JSON for docs/encyclopedia).

Coordinates. A flat degree-4 vertex is (x, y) = (theta_1, theta_2) in units of pi, in the open unit
square; theta_3 = 1 - x, theta_4 = 1 - y, and b = (a_1, a_2, a_3)/pi = (x, x + y, 1 + y).
Symmetry. g = (r, m) relabels a vertex: ray r becomes ray 0, counterclockwise (m = 0) or mirrored
(m = 1); on (x, y) this is the dihedral group of the square, with fundamental triangle
0 <= x <= y <= 1/2 (theta_1 <= theta_2 <= 90 deg).
Strata. Every (type, piece[, branch]) of a degree-4 cells tree lives on the affine subspace of its
piece: a line or a point of the square. A stratum is an orbit of these under the group, shown on a
representative. Its cards are the histories on it up to symmetry: duplicates (rotated or mirrored
types, and images under the stratum's stabilizer) are merged, after checking that their
realizability agrees.
Realizability along a line is computed exactly from the trees (src/cells.py): every split form is
affine in the line parameter t, so the leaves are intervals with rational endpoints.
Usage: python3 tools/encyclopedia_export.py OUT.json [patterns.fold ...]   (patterns: corpus pins)"""
import sys, json, math, glob, os, itertools
from fractions import Fraction as Fr
from common import fh, resolve_tol, ROOT
sys.path.insert(0, str(ROOT / 'src'))
import cells, rigid, lookup

F = Fr
SAMPLES = [F(k, 97) for k in (13, 29, 41, 50, 61, 77, 89)]          # generic parameters for keys


# ---------------------------------------------------------------- the group on (x, y)
def act(g, p):
    """Image of the vertex p = (x, y) under g = (r, m)."""
    x, y = p; s = [x, y, 1 - x, 1 - y]                                # s[i] = theta_{i+1}
    r, m = g
    if m == 0: return (s[r % 4], s[(r + 1) % 4])                      # theta_{r+1}, theta_{r+2}
    return (s[(r - 1) % 4], s[(r - 2) % 4])                           # theta_r, theta_{r-1}

G = [(r, m) for m in (0, 1) for r in range(4)]


def inverse(g):
    probe = [(F(1, 7), F(2, 9)), (F(3, 11), F(1, 5)), (F(1, 3), F(5, 13))]
    for h in G:
        if all(act(h, act(g, p)) == p for p in probe): return h


def line_through(p, q):
    """Normalised integer (a, b, c) with a x + b y = c."""
    a, b = q[1] - p[1], p[0] - q[0]; c = a * p[0] + b * p[1]
    L = 1
    for v in (a, b, c): L = L * v.denominator // math.gcd(L, v.denominator)
    a, b, c = int(a * L), int(b * L), int(c * L)
    g_ = math.gcd(math.gcd(abs(a), abs(b)), abs(c)) or 1
    a, b, c = a // g_, b // g_, c // g_
    if (a, b) < (0, 0) or (a == 0 and b < 0): a, b, c = -a, -b, -c
    if a < 0: a, b, c = -a, -b, -c
    return (a, b, c)


def clip(line, poly):
    """Segment of the line a x + b y = c inside the convex polygon (closed); None if empty or a point."""
    a, b, c = line
    if b != 0: p0, d = (F(0), F(c, b)), (F(1), F(-a, b))
    else: p0, d = (F(c, a), F(0)), (F(0), F(1))
    lo, hi = None, None
    n = len(poly)
    for i in range(n):                                               # half-planes of a ccw polygon
        P, Q = poly[i], poly[(i + 1) % n]
        ex, ey = Q[0] - P[0], Q[1] - P[1]
        f0 = ex * (p0[1] - P[1]) - ey * (p0[0] - P[0]); f1 = ex * d[1] - ey * d[0]   # >= 0 inside
        if f1 == 0:
            if f0 < 0: return None
            continue
        tt = -f0 / f1
        if f1 > 0: lo = tt if lo is None or tt > lo else lo
        else: hi = tt if hi is None or tt < hi else hi
    if lo is None or hi is None or hi <= lo: return None
    A = (p0[0] + lo * d[0], p0[1] + lo * d[1]); B = (p0[0] + hi * d[0], p0[1] + hi * d[1])
    return tuple(sorted([A, B]))


SQUARE = [(F(0), F(0)), (F(1), F(0)), (F(1), F(1)), (F(0), F(1))]
TRI = [(F(0), F(0)), (F(1, 2), F(1, 2)), (F(0), F(1, 2))]


def in_tri(p): return 0 <= p[0] <= p[1] <= F(1, 2)


# ---------------------------------------------------------------- strata of the pieces
def subspace(basis, vals):
    """The piece's subspace in (x, y): ('line', (a, b, c)) or ('point', (x, y))."""
    rows = []
    for r, v in zip(basis, vals):                                    # r . (x, x + y, 1 + y) = v
        rows.append([F(r[0] + r[1]), F(r[1] + r[2]), F(v - r[2])])
    M = [list(r) for r in rows]; rk = 0; piv = []
    for c in range(2):
        p = next((i for i in range(rk, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[rk], M[p] = M[p], M[rk]; M[rk] = [x / M[rk][c] for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] != 0: f = M[i][c]; M[i] = [a - f * b for a, b in zip(M[i], M[rk])]
        piv.append(c); rk += 1
    assert all(r[0] == 0 and r[1] == 0 and r[2] == 0 for r in M[rk:]), 'inconsistent piece'
    if rk == 2: return ('point', (M[0][2], M[1][2]))
    assert rk == 1
    a, b, c = M[0]
    p = (c / a, F(0)) if a != 0 else (F(0), c / b)
    q = (p[0] - b, p[1] + a)
    return ('line', line_through(p, q))


def image(kind, obj, g):
    if kind == 'point': return act(g, obj)
    a, b, c = obj
    p = (F(c, a), F(0)) if a else (F(0), F(c, b))
    q = (p[0] - b, p[1] + a)
    return line_through(act(g, p), act(g, q))


def representative(kind, obj):
    imgs = {g: image(kind, obj, g) for g in G}
    if kind == 'point':
        rep = min((v for v in imgs.values() if in_tri(v)), key=lambda v: (v[0], v[1]))
        return rep
    def score(L):
        s = clip(L, TRI); length = 0 if s is None else (s[1][0] - s[0][0]) ** 2 + (s[1][1] - s[0][1]) ** 2
        return (-length, L)
    return min(set(imgs.values()), key=score)


# ---------------------------------------------------------------- realizability along a line
def walk(node, alpha_beta, lo, hi, den, out, zeros=None):
    """Leaves of tree node on the open interval (lo, hi) of t; alpha_beta(f, k) -> (alpha, beta) with
    form value (alpha + beta t)/den. Appends (lo, hi, value) (value None: outside the explored region)."""
    if hi <= lo: return
    if not isinstance(node, dict): out.append((lo, hi, node)); return
    al, be = alpha_beta(node['f'], node['k'])
    ch = node['c']
    if be == 0:
        s = (al > 0) - (al < 0); walk(ch.get(str(s)), alpha_beta, lo, hi, den, out); return
    t0 = -al / be
    if zeros is not None and lo < t0 < hi: zeros.add(t0)
    neg, pos = ((lo, min(hi, t0)), (max(lo, t0), hi)) if be > 0 else ((max(lo, t0), hi), (lo, min(hi, t0)))
    walk(ch.get('-1'), alpha_beta, neg[0], neg[1], den, out, zeros)
    walk(ch.get('1'), alpha_beta, pos[0], pos[1], den, out, zeros)


def at_point(node, b, den):
    while isinstance(node, dict):
        v = sum(F(int(c)) * x for c, x in zip(node['f'], b)) + int(node['k'])
        node = node['c'].get(str((v > 0) - (v < 0)))
    return node


def merge(iv):
    out = []
    for lo, hi, v in sorted(iv):
        if out and out[-1][2] == v and out[-1][1] == lo: out[-1] = (out[-1][0], hi, v)
        else: out.append((lo, hi, v))
    return out


# ---------------------------------------------------------------- main
def main():
    out_path = sys.argv[1]
    types = []
    for f in sorted(glob.glob(str(ROOT / 'tables' / 'deg4_k*.json'))):
        name = os.path.basename(f)
        T = json.load(open(f))['types']; C = json.load(open(str(ROOT / 'cells' / name)))
        assert C['meta']['complete'] and len(C['trees']) == len(T), name
        for i, (t, tr) in enumerate(zip(T, C['trees'])):
            t['tree'] = tr; t['den'] = C['meta'].get('den', 4); t['id'] = '%s#%d' % (name[:-5], i)
            types.append(t)
    print('types', len(types), file=sys.stderr)

    strata = {}                                                   # (kind, rep) -> data
    for t in types:
        tr = t['tree']
        for key, node in tr['pieces'].items():
            vs, _, us = key.partition('|')
            vals = [int(x) for x in vs.split(',') if x != '']
            u = tuple(int(x) for x in us.split(',')) if us else None
            kind, obj = subspace(tr['basis'], vals)
            rep = representative(kind, obj)
            S = strata.setdefault((kind, rep), {'kind': kind, 'rep': rep, 'items': []})
            S['items'].append((t, node, u, [g for g in G if image(kind, obj, g) == rep]))
    print('strata', len(strata), file=sys.stderr)

    pins, pin_points = {}, set()
    for path in sys.argv[2:]:
        P_ = resolve_tol(path)
        for p in P_.interior:
            if len(P_.H[p]) != 4: continue
            D = sorted(d for d, _ in P_.H[p])
            s = [D[1] - D[0], D[2] - D[1]]
            q = (F(round(s[0] / math.pi * 2 ** 20), 2 ** 20), F(round(s[1] / math.pi * 2 ** 20), 2 ** 20))
            q = tuple(x.limit_denominator(720) if abs(float(x.limit_denominator(720)) - float(x)) < 1e-6 else x for x in q)
            c = min((v for v in (act(g, q) for g in G) if in_tri(v)), key=lambda v: (v[0], v[1]), default=None)
            if c is None: continue
            pin_points.add(c)
            pins.setdefault(c, set()).add(os.path.basename(path)[:-5])
    out_strata, checks = [], {'merged': 0, 'merge_disagree': 0, 'undecided_intervals': 0, 'C_samples': 0, 'C_fail': 0, 'formula_checks': 0, 'formula_fail': 0}
    for (kind, rep), S in sorted(strata.items(), key=lambda kv: (kv[0][0], str(kv[0][1]))):
        if kind == 'line':
            seg = clip(rep, SQUARE); P0, P1 = seg
            P = lambda tt: (P0[0] + tt * (P1[0] - P0[0]), P0[1] + tt * (P1[1] - P0[1]))
        # the stabilizer of the line: elements fixing it pointwise merge histories; an element reversing
        # it halves the domain to t in [0, 1/2], where a history and its image are both kept
        tmax = F(1)
        if kind == 'line':
            stab = [g for g in G if image('line', rep, g) == rep]
            if any(act(g, P0) == P1 for g in stab): tmax = F(1, 2)
        cards = {}
        classes = []
        for t, node, u, gs in S['items']:
            if kind == 'line':
                cl_ = {}
                for g in gs: cl_.setdefault(act(inverse(g), P0), []).append(g)
                classes += [(t, node, u, c_) for c_ in cl_.values()]
            else: classes.append((t, node, u, gs))
        for t, node, u, gs in classes:
            best = None
            for g in gs:
                gi = inverse(g)                                       # rep frame -> type frame
                # realizability intervals in the rep parameter
                if kind == 'line':
                    q0, q1 = act(gi, P(F(0))), act(gi, P(F(1)))       # type frame, affine in t
                    b0 = (q0[0], q0[0] + q0[1], 1 + q0[1]); b1 = (q1[0], q1[0] + q1[1], 1 + q1[1])
                    def ab(f, k, b0=b0, b1=b1):
                        al = sum(F(int(c)) * x for c, x in zip(f, b0)) + int(k)
                        be = sum(F(int(c)) * (y - x) for c, x, y in zip(f, b0, b1))
                        return al, be
                    zs = set(); iv = []; walk(node, ab, F(0), F(1), t['den'], iv, zs); iv = merge(iv)
                    pts = [act(gi, P(s)) for s in SAMPLES]
                else:
                    zs = set(); q = act(gi, rep); iv = [(F(0), F(1), at_point(node, (q[0], q[0] + q[1], 1 + q[1]), t['den']))]
                    pts = [q]
                # canonical key: structure in the rep frame + hidden directions at the samples
                r, m = g
                order = [(i if m == 0 else -i) for i in range(len(t['rays']))]
                vpos = [j for j, ry in enumerate(t['rays']) if ry == ['v', r]][0]
                N = len(t['rays'])
                pos = {j: ((j - vpos) % N if m == 0 else (vpos - j) % N) for j in range(N)}
                kinds = [None] * N
                for j, ry in enumerate(t['rays']):
                    kinds[pos[j]] = 'v%d' % ((ry[1] - r) % 4 if m == 0 else (r - ry[1]) % 4) if ry[0] == 'v' else 'h'
                steps = tuple(tuple(sorted(pos[j] for j in T)) for T in t['steps'])
                dirs = []
                for q in pts:
                    a = [0.0, q[0] * math.pi, (q[0] + q[1]) * math.pi, (1 + q[1]) * math.pi]
                    lams = None
                    if 'scaffold' in t:
                        lams = next((l for uu, l in lookup.branches(t, a) if uu == u), None)
                    d = cells.replay(rigid.NUM, t, a, lams) if (lams is not None or 'scaffold' not in t) else None
                    if d is None: dirs.append(None); continue
                    ar = a[r]
                    dirs.append(tuple(round(((x - ar) if m == 0 else (ar - x)) % (2 * math.pi), 6)
                                      for j, x in sorted(((pos[j], d[j]) for j in range(N) if t['rays'][j][0] == 'h'))))
                k = (tuple(kinds), steps, tuple(dirs) if any(dirs) else (t['id'], u))
                if best is None or str(k) < str(best[0]): best = (k, g, iv, zs, node)
            k, g, iv, zs, node = best
            iv = [(lo, min(hi, tmax), v) for lo, hi, v in iv if lo < tmax]
            if k in cards:
                checks['merged'] += 1
                if [(lo, hi, v) for lo, hi, v in cards[k]['iv']] != iv: checks['merge_disagree'] += 1
                continue
            checks['undecided_intervals'] += sum(1 for _, _, v in iv if v is None)
            cards[k] = {'t': t, 'u': u, 'g': g, 'iv': iv, 'kinds': k[0], 'steps': k[1], 'zeros': zs, 'node': node}

        # (C) on the cards' realizable intervals: lookup.replay (with (C)) at each interval's midpoint
        for c in cards.values():
            gi = inverse(c['g'])
            for lo, hi, v in c['iv']:
                if not v: continue
                q = act(gi, P((lo + hi) / 2) if kind == 'line' else rep)
                a = [0.0, float(q[0]) * math.pi, float(q[0] + q[1]) * math.pi, float(1 + q[1]) * math.pi]
                lams = next((l for uu, l in lookup.branches(c['t'], a) if uu == c['u']), None) if 'scaffold' in c['t'] else None
                checks['C_samples'] += 1
                if lookup.replay(c['t'], a, lams) is None: checks['C_fail'] += 1

        # breakpoints, dominance, summary function
        cl = list(cards.values())
        bps = sorted({F(0), tmax} | {x for c in cl for lo, hi, _ in c['iv'] for x in (lo, hi)})
        def val(c, tt):
            for lo, hi, v in c['iv']:
                if lo < tt < hi: return v or 0
            return 0
        mids = [(x + y) / 2 for x, y in zip(bps, bps[1:])]
        def exact(c, tq):
            q = act(inverse(c['g']), P(tq)) if kind == 'line' else act(inverse(c['g']), rep)
            v = at_point(c['node'], (q[0], q[0] + q[1], 1 + q[1]), c['t']['den'])
            if v is None: checks['undecided_points'] = checks.get('undecided_points', 0) + 1
            return v or 0
        # physical: realizable on an interval, or at one of its own split points (or the midpoint)
        for c in cl:
            own = [z for z in c['zeros'] if 0 < z <= tmax and z < 1] + ([tmax] if tmax < 1 else [])
            c['physical'] = any(v for _, _, v in c['iv']) or (kind == 'line' and any(exact(c, z) for z in own))
        phys = [c for c in cl if c['physical']]
        spts = set()
        if kind == 'line':
            for c in phys:
                spts |= {z for z in c['zeros'] if 0 < z <= tmax and z < 1} | {x for lo, hi, _ in c["iv"] for x in (lo, hi) if 0 < x <= tmax and x < 1}
            for pq in pin_points:
                for g in G:
                    R = act(g, pq)
                    if rep[0] * R[0] + rep[1] * R[1] == rep[2]:
                        tq = (R[0] - P0[0]) / (P1[0] - P0[0]) if P1[0] != P0[0] else (R[1] - P0[1]) / (P1[1] - P0[1])
                        if 0 < tq <= tmax and tq < 1: spts.add(tq)
            if tmax < 1: spts.add(tmax)
        spts = sorted(spts)
        prof = {id(c): [val(c, z) for z in mids] + [exact(c, z) for z in spts] for c in phys}
        for c in cl:
            c['dominated'] = c['physical'] and any(
                d is not c and d['t']['k'] < c['t']['k'] and all(x >= y for x, y in zip(prof[id(d)], prof[id(c)])) for d in phys)
        summary = []
        for (x, y), z in zip(zip(bps, bps[1:]), mids):
            ks = [c['t']['k'] for c in cl if val(c, z) >= 1]
            summary.append([float(x), float(y), min(ks) if ks else None, any(val(c, z) == 2 for c in cl)])
        sm = []
        for row in summary:
            if sm and sm[-1][2:] == row[2:] and sm[-1][1] == row[0]: sm[-1][1] = row[1]
            else: sm.append(row)

        # folded images in the fundamental triangle, with their rep parameters
        folds = []
        if kind == 'line':
            for g in G:
                L = image('line', rep, g); s = clip(L, TRI)
                if s is None: continue
                gi = inverse(g)
                ts = []
                for Q in s:
                    R = act(gi, Q)                                   # back on the rep line
                    tt = (R[0] - P0[0]) / (P1[0] - P0[0]) if P1[0] != P0[0] else (R[1] - P0[1]) / (P1[1] - P0[1])
                    ts.append(tt if tt <= tmax else 1 - tt)
                assert all(0 <= x <= tmax for x in ts), 'fold outside the domain'
                f_ = [[float(s[0][0]), float(s[0][1])], [float(s[1][0]), float(s[1][1])], [float(ts[0]), float(ts[1])]]
                if f_ not in folds: folds.append(f_)
        else:
            folds.append([[float(rep[0]), float(rep[1])]])

        def scaffold_formula(c):
            """Each scaffold line of the card in the rep frame: lambda = sum coef_i A_i + const (mod pi),
            exact, as (coefs over A_1..A_3, const) with Fractions in units of pi."""
            t, u, (r, m) = c['t'], c['u'], c['g']
            if 'scaffold' not in t: return None
            gi = inverse(c['g'])
            q = act(gi, P(F(1, 2)) if kind == 'line' else rep)
            a_num = [F(0), q[0], q[0] + q[1], 1 + q[1]]               # type frame, units of pi
            # A_{sigma(i)} (rep frame, in [0, 2)) for the type's visible rays at this point
            def repdir(x): return ((x - a_num[r]) if m == 0 else (a_num[r] - x)) % 2
            sig = [((i - r) % 4 if m == 0 else (r - i) % 4) for i in range(4)]
            eps = 1 if m == 0 else -1
            # a_i = eps * (A_sig(i) - A_sig(0)) + 2 n_i exactly; forms over (A_1, A_2, A_3, 1)
            forms = []
            for i in range(4):
                f = [F(0)] * 4
                if sig[i]: f[sig[i] - 1] += eps
                if sig[0]: f[sig[0] - 1] -= eps
                val_ = eps * (repdir(a_num[i]) - repdir(a_num[0]))
                f[3] += a_num[i] - val_
                forms.append(f)
            sc = t['scaffold']; steps_, piv = sc['steps'], sc['pivots']
            lam = [None] * len(steps_)
            for j, cc, d in reversed(piv):
                x = [F(0), F(0), F(0), F(u[j])]
                for qi, cq in enumerate(cc):
                    x = [xv - cq * fv for xv, fv in zip(x, forms[qi + 1])]
                for qq in range(len(steps_)):
                    if qq != j and d[qq]: x = [xv - d[qq] * lv for xv, lv in zip(x, lam[qq])]
                lam[j] = [xv / d[j] for xv in x]
            res = []
            for j, st in enumerate(steps_):
                l = [eps * (lv - fr) for lv, fr in zip(lam[j], forms[r])]
                res.append({'step': st, 'coef': [str(v) for v in l[:3]], 'const': str(l[3] % 1)})
            return res

        out_cards = []
        for c in sorted(cl, key=lambda c: (c['t']['k'], len(c['t']['steps']), c['t']['id'])):
            t = c['t']
            if not c['physical']: continue
            fm = scaffold_formula(c)
            if fm:                                                    # check the formula against the replay
                gi = inverse(c['g']); r, m = c['g']
                for lo, hi, v in c['iv']:
                    if not v: continue
                    tt = (lo + hi) / 2
                    q = act(gi, P(tt) if kind == 'line' else rep)
                    a = [0.0, float(q[0]) * math.pi, float(q[0] + q[1]) * math.pi, float(1 + q[1]) * math.pi]
                    lams = next(l for uu, l in lookup.branches(t, a) if uu == c['u'])
                    d = cells.replay(rigid.NUM, t, a, lams)
                    Q = P(tt) if kind == 'line' else rep
                    A = [0, Q[0], Q[0] + Q[1], 1 + Q[1]]
                    for f in fm:
                        anchor = t['steps'][f['step']][0]
                        actual = (((d[anchor] - a[r]) if m == 0 else (a[r] - d[anchor])) / math.pi) % 1
                        pred = float(sum(F(cf) * A[i + 1] for i, cf in enumerate(f['coef'])) + F(f['const'])) % 1
                        checks['formula_checks'] += 1
                        if min(abs(actual - pred), 1 - abs(actual - pred)) > 1e-7: checks['formula_fail'] += 1
            out_cards.append({'id': t['id'], 'k': t['k'], 'rays': t['rays'], 'steps': t['steps'],
                              'scaffold': t.get('scaffold'), 'u': list(c['u']) if c['u'] else None,
                              'g': list(c['g']), 'iv': [[float(lo), float(hi), v] for lo, hi, v in c['iv']],
                              'iv_exact': [[str(lo), str(hi)] for lo, hi, v in c['iv']],
                              'physical': c['physical'], 'dominated': c['dominated'],
                              'formula': fm})
        ordered = sorted(phys, key=lambda c: (c['t']['k'], len(c['t']['steps']), c['t']['id']))
        special = [[float(tq), str(tq), [exact(c, tq) for c in ordered]] for tq in spts] if kind == 'line' else []
        rec = {'kind': kind, 'cards': out_cards, 'summary': sm, 'folds': folds, 'special': special,
               'nonphysical': sum(1 for c in cl if not c['physical'])}
        if kind == 'line':
            rec.update(tmax=float(tmax), line=list(rep), P0=[float(P0[0]), float(P0[1])], P1=[float(P1[0]), float(P1[1])],
                       P0_exact=[str(P0[0]), str(P0[1])], P1_exact=[str(P1[0]), str(P1[1])])
        else:
            rec.update(point=[float(rep[0]), float(rep[1])], point_exact=[str(rep[0]), str(rep[1])])
        out_strata.append(rec)

    json.dump({'strata': out_strata, 'pins': [[float(k[0]) * 180, float(k[1]) * 180, str(k[0]), str(k[1]), sorted(v)] for k, v in sorted(pins.items())], 'checks': checks},
              open(out_path, 'w'), separators=(',', ':'))
    print(json.dumps(checks), file=sys.stderr)
    print('strata', len(out_strata), 'lines', sum(s['kind'] == 'line' for s in out_strata),
          'cards (physical)', sum(len(s['cards']) for s in out_strata), 'nonphysical', sum(s['nonphysical'] for s in out_strata),
          'shown (physical, not dominated)', sum(c['physical'] and not c['dominated'] for s in out_strata for c in s['cards']),
          'pins', len(pins), file=sys.stderr)


if __name__ == '__main__':
    main()
