"""Name the folds of the running pattern (examples/example1_running.fold) as in the paper."""
import math


def label_running(P, folds):
    V = P.V; names = {}
    near = lambda p, q: abs(p[0] - q[0]) < 0.02 and abs(p[1] - q[1]) < 0.02
    for x, ks in folds.items():
        segs = [(V[P.crease[k][0]], V[P.crease[k][1]]) for k in ks if isinstance(k, int)]
        pts = [p for s in segs for p in s]
        horiz = all(abs(a[1] - b[1]) < 1e-6 for a, b in segs)
        vert = all(abs(a[0] - b[0]) < 1e-6 for a, b in segs)
        touches = lambda q: any(near(p, q) for p in pts)
        if horiz:
            names[x] = {6: 'A^{+}', 4: 'A^{0}', 2: 'A^{-}'}[round(segs[0][0][1])]
        elif vert and all(abs(a[0] - 1.313) < 0.01 for a, _ in segs): names[x] = 'E'
        elif any(abs(p[0] - 1.313) < 0.01 for p in pts): names[x] = "E'"
        elif touches((4, 0)): names[x] = 'C'
        elif touches((5.107, 0)): names[x] = "C'"
        elif touches((6, 0)): names[x] = 'F'
        elif touches((8, 2)) and touches((8, 6)): names[x] = 'I'
        else:
            corner = (8, 8) if touches((8, 8)) else (8, 0)
            a, b = next(s for s in segs if near(s[0], corner) or near(s[1], corner))
            ang = math.degrees(math.atan2(abs(a[1] - b[1]), abs(a[0] - b[0])))
            i = min(((1, 22.5), (2, 45.0), (3, 67.5)), key=lambda t: abs(t[1] - ang))[0]
            names[x] = 'P_%d^{%s}' % (i, '+' if corner == (8, 8) else '-')
    assert len(set(names.values())) == len(folds), names
    return names
