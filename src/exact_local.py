"""Exact search for local histories with cancellations at one vertex, with no candidate fold lines.

Fix the cyclic order of the visible and hidden rays. By Theorem 2, condition (L) for a step is then
a linear congruence modulo pi in the hidden directions. Flatness of each state is one more
congruence, and condition (C) is a Boolean combination of linear inequalities once every angle
is normalised modulo 2pi by an integer that the inequalities determine. So, for a given number k
of hidden rays, a depth-first search over step sequences decides whether the vertex has a local
history. The search runs once for every placement of the k hidden rays in the gaps between
visible rays, with an incremental z3 problem over the real hidden directions that prunes every
infeasible prefix.

Search space, for a given k:
  * k hidden rays, each created by one step and erased by a later one;
  * hidden directions pairwise distinct and distinct from the visible ones, at least DELTA apart;
  * a step is any even set of rays: visible rays not yet made, hidden rays not yet made, and
    hidden rays present (erased).  Visible rays are never erased.
Congruences hold to within EPS_EQ, and (C) treats a preimage within DELTA of a piece's boundary
as not inside the piece.  Every history found is re-checked in floating point with the
functions of fold_history (maps, collision_free).

Usage: python3 src/exact_local.py pattern.fold [vertex ...] [--kmax=K] [--seconds=S per k]
"""
import itertools, math, sys, time
import z3
import fold_history as fh

PI = math.pi
DELTA = 1e-6
EPS_EQ = 1e-9
ZPI = z3.RealVal(repr(PI))


class Undecided(Exception):
    """The deadline passed, or z3 returned unknown, before the search finished."""


class Lin:
    """A linear term sum(coef * x_i) + const in the hidden directions x_i."""
    def __init__(self, coef=None, const=0.0): self.coef, self.const = dict(coef or {}), const
    def __add__(self, o):
        o = o if isinstance(o, Lin) else Lin(const=o); c = dict(self.coef)
        for i, a in o.coef.items(): c[i] = c.get(i, 0) + a
        return Lin({i: a for i, a in c.items() if a}, self.const + o.const)
    def __rmul__(self, k): return Lin({i: k * a for i, a in self.coef.items() if k * a}, k * self.const)
    def __sub__(self, o): return self + (-1) * (o if isinstance(o, Lin) else Lin(const=o))
    def z3(self, X): return z3.Sum([z3.RealVal(repr(self.const))] + [a * X[i] for i, a in self.coef.items()])
    def val(self, x): return self.const + sum(a * x[i] for i, a in self.coef.items())
    def is_const(self): return not self.coef


class Vertex:
    def __init__(self, A, k):
        """A: sorted visible directions in [0, 2pi); k: number of hidden rays."""
        self.A, self.k = A, k

    def placements(self):
        n = len(self.A)
        for comp in itertools.product(range(self.k + 1), repeat=n):
            if sum(comp) == self.k: yield comp

    def search(self, comp, limit=1, deadline=math.inf):
        """Histories for one placement: comp[g] hidden rays in the gap after visible ray g."""
        A, n = self.A, len(self.A)
        self.rays = []                                    # ('v', g) or ('h', i), counterclockwise
        for g in range(n):
            self.rays.append(('v', g))
            self.rays += [('h', i) for i in range(sum(comp[:g]), sum(comp[:g + 1]))]
        X = [z3.Real('x%d' % i) for i in range(self.k)]
        self.X = X
        s = z3.Solver(); self.s = s; self.nm = 0
        pos = []
        for r in self.rays: pos.append(Lin(const=A[r[1]]) if r[0] == 'v' else Lin({r[1]: 1}))
        self.pos = pos
        top = Lin(const=A[0] + 2 * PI)
        for j, r in enumerate(self.rays):                  # cyclic order, with margin DELTA
            if r[0] == 'h':
                nxt = pos[j + 1] if j + 1 < len(self.rays) else top
                s.add((pos[j] - pos[j - 1]).z3(X) >= DELTA, (nxt - pos[j]).z3(X) >= DELTA)
        N = len(self.rays)
        status = {j: 0 for j in range(N)}                 # 0 not made, 1 present, 2 erased (hidden)
        self.found, self.nodes, self.deadline = [], 0, deadline
        self.dfs(status, [], limit)
        return self.found

    # ---------------------------------------------------------------- constraints
    def cong(self, e, cs):
        """e == 0 (mod pi).  Constant terms are decided now; returns False if one fails."""
        if e.is_const():
            return abs(e.const - round(e.const / PI) * PI) < 1e-7
        self.nm += 1; m = z3.Int('m%d' % self.nm); d = e.z3(self.X) - m * ZPI
        cs += [d <= EPS_EQ, d >= -EPS_EQ]
        return True

    def step_constraints(self, status, T):
        """Constraints for step T (ray indices) from the state given by status, or None."""
        N, pos, X = len(self.rays), self.pos, self.X
        S = [j for j in range(N) if status[j] == 1]
        cs = []
        x = T[0]
        for y in T[1:]:                                    # (L), by Theorem 2
            W = [(x + d) % N for d in range(1, (y - x) % N) if status[(x + d) % N] == 1]
            sv = Lin()
            for w in W: sv = pos[w] - sv
            e = (pos[x] + pos[y] - 2 * sv) if len(W) % 2 else (pos[y] - pos[x] - 2 * sv)
            if not self.cong(e, cs): return None
        S2 = sorted([j for j in S if j not in T] + [j for j in T if status[j] == 0])
        alt = Lin()                                        # flatness of the new state
        for i, j in enumerate(S2): alt = alt + ((-1) ** i) * pos[j]
        if not self.cong(alt, cs): return None
        if S: cs += self.collision(S, T)
        return cs

    def collision(self, S, T):
        """Condition (C): no two pieces of different colour straddle the same ray of the fold line."""
        pos, X, N = self.pos, self.X, len(self.rays)
        top = 2 * PI
        maps = [(1, Lin())]
        for j in S[1:]:
            e, o = maps[-1]; maps.append((-e, o + (2 * e) * pos[j]))
        def sector(t):                                      # sector of S containing ray t, or starting at it
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
                a = pos[R[i]]; b = pos[R[i + 1]] if i + 1 < len(R) else pos[R[0]] + top
                e, o = maps[sector(R[i])]
                u = e * (rr - o) - a                         # preimage, relative to the piece start
                self.nm += 1; q = z3.Int('q%d' % self.nm); w = u.z3(X) + 2 * ZPI * q
                cs += [w >= 0, w < 2 * ZPI]
                st[mv[i]].append(z3.And(w > DELTA, w < (b - a).z3(X) - DELTA))
            cs.append(z3.Or(z3.Not(z3.Or(st[0])), z3.Not(z3.Or(st[1]))))
        return cs

    # ---------------------------------------------------------------- search
    def dfs(self, status, hist, limit):
        if len(self.found) >= limit: return
        self.nodes += 1
        N = len(self.rays)
        if all(status[j] == 1 if self.rays[j][0] == 'v' else status[j] == 2 for j in range(N)):
            self.found.append((self.s.model(), list(hist))); return
        cand = [j for j in range(N) if status[j] == 0 or (status[j] == 1 and self.rays[j][0] == 'h')]
        for size in range(2, len(cand) + 1, 2):
            for T in itertools.combinations(cand, size):
                cs = self.step_constraints(status, list(T))
                if cs is None: continue
                self.s.push(); self.s.add(cs)
                left = self.deadline - time.monotonic()
                if left <= 0: raise Undecided()
                if left < 1e6: self.s.set('timeout', max(1, int(1000 * left)))
                r = self.s.check()
                if r == z3.unknown: raise Undecided()
                if r == z3.sat:
                    st2 = dict(status)
                    for j in T: st2[j] = 1 if status[j] == 0 else 2
                    hist.append(T); self.dfs(st2, hist, limit); hist.pop()
                self.s.pop()
                if len(self.found) >= limit: return

    def directions(self, m):
        return [self.A[r[1]] if r[0] == 'v' else float(m.eval(self.X[r[1]]).as_fraction()) % (2 * PI)
                for r in self.rays]

    def verify(self, m, hist):
        """Re-check a history in floating point with fold_history's maps and collision_free."""
        D = self.directions(m); S = []
        for T in hist:
            phis = fh.maps(sorted(S))
            if phis is None: return 'state not flat before %s' % (T,)
            Sd = sorted(S)
            imgs = [fh.apply(phis[fh.sector_of(D[j], Sd)], D[j]) % PI for j in T]
            if any(abs((a - imgs[0] + PI / 2) % PI - PI / 2) > 1e-6 for a in imgs): return 'not one line at %s' % (T,)
            if not fh.collision_free(Sd, [D[j] for j in T], phis, imgs[0]): return 'collision at %s' % (T,)
            S = [x for x in S if not any(fh.eq(x, D[j]) for j in T)] + [D[j] for j in T if not any(fh.eq(D[j], x) for x in S)]
        return None if fh.maps(sorted(S)) is not None else 'final state not flat'

    def describe(self, m, hist):
        D = self.directions(m); out = []; made = set()
        for T in hist:
            out.append('{' + ', '.join(('%.4f' % math.degrees(D[j])) + ('' if self.rays[j][0] == 'v' else (' erased' if j in made else ' hidden'))
                                       for j in T) + '}')
            made |= set(T)
        return ' -> '.join(out)


def min_history(A, kmax=4, log=print, seconds=math.inf, kmin=0):
    """Smallest k in [kmin, kmax] with a local history, and one history: (k, (text, float check)).
    (None, None) if there is none up to kmax; ('undecided', k) if the search for k ran out of
    seconds (each k gets its own allowance)."""
    for k in range(kmin, kmax + 1):
        V = Vertex(A, k); nodes = 0; deadline = time.monotonic() + seconds
        try:
            for comp in V.placements():
                f = V.search(comp, deadline=deadline); nodes += V.nodes
                if f:
                    m, hist = f[0]; bad = V.verify(m, hist)
                    return k, (V.describe(m, hist), bad)
        except Undecided:
            log('  k=%d: undecided after %.0f s' % (k, seconds)); return 'undecided', k
        log('  k=%d: no history (%d search nodes)' % (k, nodes))
    return None, None


def main(argv):
    kmax, secs, args = 4, math.inf, []
    for a in argv:
        if a.startswith('--kmax='): kmax = int(a[7:])
        elif a.startswith('--seconds='): secs = float(a[10:])
        else: args.append(a)
    path, verts = args[0], [int(v) for v in args[1:]]
    for tol in (1e-6, 1e-4, 5e-4, 2e-3):
        P = fh.Pattern(path, tol)
        if all(fh.maps(sorted(d for d, _ in P.H[p])) is not None for p in P.interior): break
    for p in verts or P.interior:
        A = sorted(d for d, _ in P.H[p])
        print('vertex %d at (%.4f, %.4f), rays %s' % (p, *P.V[p], ', '.join('%.4f' % math.degrees(a) for a in A)), flush=True)
        k, h = min_history(A, kmax, log=lambda s: print(s, flush=True), seconds=secs)
        if k == 'undecided': print('  undecided at k=%d' % h, flush=True); continue
        print('  minimum %s hidden rays%s' % (k if k is not None else '> %d' % kmax,
              '' if h is None else ': ' + h[0] + ('' if h[1] is None else '   [float check FAILED: %s]' % h[1])), flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
