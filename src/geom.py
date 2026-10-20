"""Arithmetic on directions for the layered model, numeric or exact.

Every geometric decision of src/rigid.py goes through an arithmetic object `ar`:
  Num()      directions are floats in radians, compared with the tolerance fold_history.EPS;
  Sym(...)   directions are affine forms with rational coefficients in the variables b_1..b_m
             (directions in units of pi), and every sign is decided by an Oracle that holds a
             context of linear constraints; when the context does not decide a sign it raises
             Need(form, signs).
The functions below (nrm, eq, between, uniq, maps, ...) are written once, for both.
"""
import math
from fractions import Fraction as Fr
from functools import cmp_to_key

EPS = 1e-7


class Num:
    """Floats in radians."""
    PI, TWO_PI, ZERO = math.pi, 2 * math.pi, 0.0
    def sgn(self, x): return 0 if abs(x) < EPS else (1 if x > 0 else -1)
    def floor_div(self, x, period):
        q = math.floor(x / period)
        if period - (x - q * period) < EPS: q += 1                  # snap just below a period to it
        return q
    def half(self, x): return x / 2


class Need(Exception):
    """The oracle's context does not decide the sign of form; signs lists the feasible ones."""
    def __init__(self, form, signs): self.form, self.signs = form, signs


DEN = 4                          # Lin numerators are integers over DEN (midpoints halve once at most);
                                 # 4 times the denominators of scaffold lines for tables with them


def set_den(d):
    global DEN
    DEN = d


class Lin:
    """An affine form (sum(c_i b_i) + k) / DEN with integer numerators (directions in units of pi)."""
    __slots__ = ('c', 'k')
    def __init__(self, c, k=0): self.c, self.k = tuple(c), k
    def __add__(self, o):
        if not isinstance(o, Lin): return Lin(self.c, self.k + _num(o))
        return Lin(tuple(a + b for a, b in zip(self.c, o.c)), self.k + o.k)
    __radd__ = __add__
    def __neg__(self): return Lin(tuple(-a for a in self.c), -self.k)
    def __sub__(self, o): return self + (-o if isinstance(o, Lin) else -Fr(o))
    def __rsub__(self, o): return (-self) + o
    def __mul__(self, t):
        assert isinstance(t, int), t
        return Lin(tuple(a * t for a in self.c), self.k * t)
    __rmul__ = __mul__
    def half(self):
        assert all(a % 2 == 0 for a in self.c) and self.k % 2 == 0, 'denominator DEN too small'
        return Lin(tuple(a // 2 for a in self.c), self.k // 2)
    def div(self, q):
        assert all(a % q == 0 for a in self.c) and self.k % q == 0, 'denominator DEN too small'
        return Lin(tuple(a // q for a in self.c), self.k // q)
    def is_const(self): return not any(self.c)
    def num_at(self, P, L):
        """The numerator of the value at the point b_i = P_i / L, over DEN * L."""
        return sum(c * x for c, x in zip(self.c, P)) + self.k * L
    def __repr__(self):
        t = ['%d*b%d' % (c, i + 1) for i, c in enumerate(self.c) if c] + ([str(self.k)] if self.k or not any(self.c) else [])
        return '(' + ' + '.join(t) + ')/%d' % DEN


def _num(x):
    """A rational constant (in units of pi) as a numerator over DEN."""
    v = Fr(x) * DEN
    assert v.denominator == 1, x
    return int(v)


class Sym:
    """Exact affine forms; signs decided by an oracle."""
    def __init__(self, m, oracle):
        self.m, self.oracle = m, oracle
        self.ZERO = Lin((0,) * m); self.PI = Lin((0,) * m, DEN); self.TWO_PI = Lin((0,) * m, 2 * DEN)
    def var(self, i): return Lin(tuple(DEN if j == i else 0 for j in range(self.m)))
    def const(self, k): return Lin((0,) * self.m, _num(k))
    def sgn(self, x):
        if not isinstance(x, Lin): x = self.const(x)
        return self.oracle.sign(x)
    def floor_div(self, x, period):
        """The integer q with 0 <= x - q period < period, decided by sign queries."""
        p = period.k                                   # period is a constant Lin: PI or TWO_PI
        key = ('fd', x.c, x.k, p)
        q = self.oracle.recall(key)
        if q is not None: return q
        q = self._floor_div(x, p)
        self.oracle.store(key, q)
        return q
    def _floor_div(self, x, p):
        P, L = self.oracle.point()
        q = x.num_at(P, L) // (p * L)
        lo, hi = self.sgn(x - Lin((0,) * len(x.c), q * p)), self.sgn(x - Lin((0,) * len(x.c), (q + 1) * p))
        assert lo >= 0 and hi < 0, 'oracle inconsistent'
        return q
    def half(self, x): return x.half()


# ---------------------------------------------------------------- functions for both
def nrm(ar, x): return x - ar.floor_div(x, ar.TWO_PI) * ar.TWO_PI
def modpi(ar, x): return x - ar.floor_div(x, ar.PI) * ar.PI
def eq(ar, x, y): return ar.sgn(nrm(ar, x - y)) == 0
def between(ar, x, a, b):
    """x strictly inside the counterclockwise arc from a to b (the full circle minus a if a = b)."""
    if eq(ar, a, b): return not eq(ar, x, a)
    da, db = nrm(ar, x - a), nrm(ar, b - a)
    return ar.sgn(da) > 0 and ar.sgn(db - da) > 0
def sort(ar, xs): return sorted(xs, key=cmp_to_key(lambda x, y: ar.sgn(x - y)))
def uniq(ar, xs):
    out = []
    for x in sort(ar, [nrm(ar, x) for x in xs]):
        if not out or not eq(ar, x, out[-1]): out.append(x)
    if len(out) > 1 and eq(ar, out[0], out[-1]): out.pop()
    return out
def mid(ar, a, b):
    if a is None: return ar.ZERO
    d = nrm(ar, b - a)
    if ar.sgn(d) == 0: d = ar.TWO_PI
    return nrm(ar, a + ar.half(d))
def maps(ar, S):
    """Folding maps (e, o) of the sorted state S, x -> e x + o; None if S is not flat."""
    if not S: return [(1, ar.ZERO)]
    out = [(1, ar.ZERO)]
    for c in S[1:]:
        s, o = out[-1]; out.append((-s, nrm(ar, 2 * s * c + o)))
    s, o = out[-1]
    return out if -s == 1 and eq(ar, nrm(ar, 2 * s * S[0] + o), ar.ZERO) else None
def apply(ar, f, x): return nrm(ar, f[0] * x + f[1])
def inverse_apply(ar, f, y): return nrm(ar, f[0] * (y - f[1]))
def sector_of(ar, x, S):
    if not S: return 0
    m = len(S)
    for j in range(m):
        if eq(ar, x, S[j]) or between(ar, x, S[j], S[(j + 1) % m]): return j
    raise ValueError('no sector')
