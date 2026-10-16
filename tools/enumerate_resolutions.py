"""List every minimal-cost choice of hidden rays (within the search space) for a pattern.
Usage: python3 tools/enumerate_resolutions.py examples/komatsu_gentle_dragon.fold"""
import sys
from common import fh
z3 = fh.z3
P = fh.Pattern(sys.argv[1]); res = fh.solve(P, log=lambda s: None)
M, opts = res['model'], res['opts']; cost = M.cost(); best = res['m'].eval(cost).as_long()
s = M.s; s.add(cost == best); found = []
while s.check(list(M.track)) == z3.sat:
    m = s.model()
    if fh.refine(P, M, opts, m, [s]): continue
    lits = [M.sel[p, i] for (p, i), b in M.sel.items() if opts[p][i][0] and z3.is_true(m.eval(b))]
    found.append({p: [P.deg(r) for r in opts[p][i][0]] for (p, i), b in M.sel.items() if opts[p][i][0] and z3.is_true(m.eval(b))})
    s.add(z3.Or([z3.Not(l) for l in lits]))
common = {p for p in found[0] if all(f.get(p) == found[0][p] for f in found)}
print('%d minimal selections (cost %d); they differ at vertices %s' % (len(found), best, sorted({p for f in found for p in f} - common)))
for f in found: print({p: f[p] for p in sorted(f) if p not in common})
