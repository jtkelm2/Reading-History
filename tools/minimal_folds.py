"""Are the earliest folds straight lines?  For one pattern: solve, analyse (scope 'minimal' by
default: what every minimal resolution forces, visible creases only), and for every fold count
the distinct lines among its creases.  Prints one JSON line.
Usage: python3 tools/minimal_folds.py pattern.fold [--scope=minimal|geometry] [--budget=SECONDS]"""
import sys, json, math
import networkx as nx
from common import fh, resolve_tol

args = [a for a in sys.argv[1:] if not a.startswith('--')]
opt = dict(a[2:].split('=') for a in sys.argv[1:] if a.startswith('--'))
scope, budget = opt.get('scope', 'minimal'), float(opt.get('budget', 600))
P = resolve_tol(args[0]); out = {'file': args[0].rsplit('/', 1)[-1], 'scope': scope}
res = fh.solve(P, log=lambda s: None, budget=budget)
out['status'] = res['status']
if res['status'] == 'ok':
    folds, lt, val, paths = fh.analyse(P, res, budget=budget, scope=scope)
    V = P.V

    def line_of(k):
        """The line through crease k, as (direction mod pi, signed distance from the origin)."""
        (x0, y0), (x1, y1) = V[P.crease[k][0]], V[P.crease[k][1]]
        th = math.atan2(y1 - y0, x1 - x0) % math.pi
        return th, -math.sin(th) * x0 + math.cos(th) * y0

    def nlines(ks):
        L = []
        for k in ks:
            th, d = line_of(k)
            if not any(min(abs(th - a), math.pi - abs(th - a)) < 1e-6 and abs(d - (b if abs(th - a) < 1 else -b)) < 1e-6 * max(1, abs(d)) for a, b in L):
                L.append((th, d))
        return len(L)

    vis = {x: [k for k in ks if isinstance(k, int)] for x, ks in folds.items()}
    vis = {x: ks for x, ks in vis.items() if ks}
    lt_v = {(a, b) for a, b in lt if a in vis and b in vis}
    G = nx.DiGraph(); G.add_nodes_from(vis); G.add_edges_from(lt_v)
    minimal = [x for x in vis if G.in_degree(x) == 0 and G.out_degree(x) > 0]   # earliest, not isolated
    isolated = [x for x in vis if G.degree(x) == 0]
    out.update(folds=len(vis), relations=len(lt_v), unknown=res['analysis_unknown'],
               lines_minimal=sorted(nlines(vis[x]) for x in minimal),
               lines_other=sorted(nlines(vis[x]) for x in vis if x not in minimal and x not in isolated),
               lines_isolated=sorted(nlines(vis[x]) for x in isolated))
print(json.dumps(out), flush=True)
