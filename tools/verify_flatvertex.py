"""Recompute every computed claim of flatvertex.tex with src/rigid.py and assert it.
Usage: python3 tools/verify_flatvertex.py"""
import itertools, math, sys
from common import fh
sys.path.insert(0, str(fh.__file__).rsplit('/', 1)[0])
import rigid

Q = math.pi / 4                       # the unit of Section 4
E = math.pi / 8                       # the unit of the example vertex A'
D_ = math.radians
deg = lambda x: '%g' % round(math.degrees(x) % 360, 4)


def choices(st, D, reverse):
    for mu0 in (0, 1):
        for up in (0, 1):
            for sense in (list(rigid.sense_options(st, D, mu0)) if reverse else [None]):
                yield mu0, up, sense


def reachable(steps, reverse):
    """All layered states reachable by the steps, with the choices that reach them."""
    states = [(rigid.initial(), [])]
    for D in steps:
        states = [(ns, ch + [(mu0, up, sense)]) for st, ch in states for mu0, up, sense in choices(st, D, reverse)
                  for ns in [rigid.step(st, D, mu0, up, sense)] if ns is not None]
    return states


def key(st, unit):
    name = ['(%g,%g)' % (round(st.S[j] / unit, 6), round((st.S[(j + 1) % len(st.S)] / unit) % (2 * math.pi / unit) or 2 * math.pi / unit, 6))
            for j in range(len(st.S))]
    return tuple((round(a / unit, 6), round(b / unit, 6), tuple(name[f] for f in o)) for a, b, o in st.arcs if a is not None)


# ---- Section 4: H1, H2, H3
H = {'H1': [(0, 4), (1, 7), (3, 5)], 'H2': [(1, 5), (3, 7), (0, 2, 4, 6), (2, 6)], 'H3': [(1, 5), (3, 7), (0, 4)]}
res = {(n, r): rigid.realizable([[k * Q for k in st] for st in h], reverse=r)[0] for n, h in H.items() for r in (False, True)}
print('H1, H2, H3 rigidly realizable:', [res[n, False] for n in H], ' with reverse folds:', [res[n, True] for n in H])
assert [res[n, False] for n in H] == [True, False, False] and [res[n, True] for n in H] == [True, False, False]

paper = {'f2f1f6f5f4f7f8f3', 'f3f8f7f4f5f6f1f2', 'f4f7f8f3f2f1f6f5', 'f5f6f1f2f3f8f7f4',
         'f6f5f2f1f8f3f4f7', 'f7f4f3f8f1f2f5f6', 'f8f3f4f7f6f5f2f1', 'f1f2f5f6f7f4f3f8'}
S3 = reachable([[k * Q for k in st] for st in H['H2'][:3]], False)
orders = {''.join('f%d' % (f + 1) for f in next(o for a, b, o in st.arcs if len(o) == 8)) for st, _ in S3}
print('H2, steps 1-3: %d of 64 choices give layered states; orders as in the proof of the theorem: %s' % (len(S3), orders == paper))
assert len(S3) == 64 and orders == paper

# ---- Section 3: fold in half, then a second line (inside and outside reverse folds)
for ang, kind in ((45, 'inside'), (120, 'outside')):
    steps = [[0.0, math.pi], [D_(ang), D_(360 - ang)]]
    rig, rev = reachable(steps, False), reachable(steps, True)
    krig = {key(st, D_(1)) for st, _ in rig}
    new = {}
    for st, ch in rev:
        k = key(st, D_(1))
        if k not in krig: new.setdefault(k, []).append(ch)
    print('fold in half, then the %d deg line: %d layer states with rigid steps, %d with reverse folds' % (ang, len(krig), len(krig) + len(new)))
    for k, chs in new.items():
        # which run was split at step 2: the one with joint 0 (the tip) or the one with joint 180 (the base)
        R = rigid.geom.uniq(rigid.NUM, [0.0, math.pi, D_(ang), D_(360 - ang)])
        runs_split = set()
        for ch in chs:
            mu0, up, sense = ch[1]
            for k2 in range(len(R)):
                nxt = (k2 + 1) % len(R)
                if sense and sense[k2] != sense[nxt] and not any(fh.eq(R[nxt], d) for d in steps[1]):
                    runs_split.add('tip (joint 0)' if fh.eq(R[nxt], 0.0) else 'base (joint 180)')
        print('   new:', ' ; '.join('over (%g,%g): %s' % (a, b, ' < '.join(o)) for a, b, o in k if len(o) == 4), '| split run:', runs_split)
        assert runs_split == ({'base (joint 180)'} if kind == 'inside' else {'tip (joint 0)'})
    assert len(krig) == 4 and len(new) == 2
    # why the other kind of reverse fold is rejected
    st1 = rigid.step(rigid.initial(), steps[0], 0, 0)
    reasons = set()
    for mu0, up, sense in choices(st1, steps[1], True):
        R = rigid.geom.uniq(rigid.NUM, list(st1.S) + steps[1])
        if not sense or not any(sense): continue
        split_at = {R[(k2 + 1) % len(R)] for k2 in range(len(R)) if sense[k2] != sense[(k2 + 1) % len(R)]
                    and not any(fh.eq(R[(k2 + 1) % len(R)], d) for d in steps[1])}
        other = 0.0 if kind == 'inside' else math.pi            # the split giving the other kind
        if any(fh.eq(x, other) for x in split_at):
            why = []
            if rigid.step(st1, steps[1], mu0, up, sense, why) is None: reasons.update(why)
    print('   the %s reverse fold here is rejected by: %s' % ('outside' if kind == 'inside' else 'inside', sorted(reasons)))

# ---- Section 3: the vertex A' = {3, 8, 10, 13} pi/8 with the history ({8, 0}, {3, 13}, {10, 0})
K = [[8 * E, 0.0], [3 * E, 13 * E], [10 * E, 0.0]]
r0, r1 = rigid.realizable(K, False)[0], rigid.realizable(K, True)
print("A': rigidly realizable %s, with reverse folds %s" % (r0, r1[0]))
assert not r0 and r1[0]
stacks = set()
for st, ch in reachable(K[:2], False):
    names = ['(%g,%g)' % (round(st.S[j] / E), round(st.S[(j + 1) % 4] / E) or 16) for j in range(4)]
    for a, b, o in st.arcs:
        if a is not None and len(o) == 4: stacks.add(tuple(names[f] for f in o))
print("A': after steps 1-2, %d different stacks over the arc where all four faces overlap:" % len(stacks))
for s_ in sorted(stacks):
    mv = [i for i, f in enumerate(s_) if f in ('(8,13)', '(13,16)')]
    print('   ', ' < '.join(s_), '  faces (8,13), (13,16) at positions', mv)
    assert mv in ([1, 2], [0, 3])
# an explicit realization: half-turns at steps 1 and 2, then a split at the joint 13
found = None
for st, ch in reachable(K[:2], False):
    names = ['(%g,%g)' % (round(st.S[j] / E), round(st.S[(j + 1) % 4] / E) or 16) for j in range(4)]
    for mu0, up, sense in choices(st, K[2], True):
        if sense and any(sense):
            ns = rigid.step(st, K[2], mu0, up, sense)
            if ns is not None and found is None:
                o4 = next(o for a, b, o in st.arcs if a is not None and len(o) == 4)
                found = (' < '.join(names[f] for f in o4), ns)
print("A': a realization: steps 1-2 rigid, giving", found[0], '; step 3 splits the run (10,16) at 13')
names = ['(%g,%g)' % (round(found[1].S[j] / E), round(found[1].S[(j + 1) % 4] / E) or 16) for j in range(4)]
for a, b, o in found[1].arcs:
    if a is not None and len(o) > 1: print('    final stack over (%g,%g):' % (round(a / E, 4), round(b / E, 4)), ' < '.join(names[f] for f in o))
print('all assertions hold')
