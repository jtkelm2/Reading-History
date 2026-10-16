"""Compare, at every vertex that needs hidden rays or a collapse, the first level found by the
solver's local search (hidden-only steps restricted to bisector images) with the exact minimum
number of hidden rays (src/exact_local.py, any directions).  One JSON line per vertex.
Usage: python3 tools/compare_local.py pattern.fold [--seconds=S per k] [--kmax=K]"""
import sys, json, math
from common import fh, resolve_tol
sys.path.insert(0, str(fh.__file__).rsplit('/', 1)[0])
import exact_local as el

args = [a for a in sys.argv[1:] if not a.startswith('--')]
opt = dict(a[2:].split('=') for a in sys.argv[1:] if a.startswith('--'))
secs, kmax = float(opt.get('seconds', 60)), int(opt.get('kmax', 4))
P = resolve_tol(args[0]); name = args[0].rsplit('/', 1)[-1]
for p in P.interior:
    first = next(((c, k) for c, k in fh.levels() if P.options(p, k, c)), None)
    if first == (0, 0): continue
    A = sorted(d for d, _ in P.H[p])
    top = first[1] - 1 if first and first[0] == 0 else kmax      # search below the solver's level
    k, h = el.min_history(A, top, log=lambda s: None, seconds=secs)
    out = {'file': name, 'vertex': p, 'degree': len(A), 'rays_deg': [round(math.degrees(a), 4) for a in A],
           'solver_first_level': first, 'searched_up_to': top}
    if k == 'undecided': out.update(exact='undecided', at_k=h)
    elif k is None: out.update(exact='none up to %d' % top)
    else: out.update(exact=k, history=h[0], float_check=h[1] or 'ok')
    print(json.dumps(out), flush=True)
