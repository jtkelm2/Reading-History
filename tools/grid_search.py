"""Breadth-first search over layered states of a single vertex, with every fold line on the grid of
multiples of 2pi/(2N) (so all creases lie on the 2N grid directions). Steps are those of
flatvertex (any even subset of the preimages of a line, creating or erasing any rays, visible or
not), each tried with both labelings, both up sides and (with --reverse) all sense assignments.
Usage: python3 tools/grid_search.py N target_ray_indices... [--reverse] [--depth=D] [--seconds=S]
(run from the repository root; e.g. N = 8 for the 22.5-degree grid, rays in units of 22.5 degrees)"""
import sys, math, time, itertools, collections, json
sys.path.insert(0, 'src')
import rigid, geom
NUM = rigid.NUM

args = [a for a in sys.argv[1:] if not a.startswith('--')]
N = int(args[0]); M = 2 * N                    # M grid directions, N grid lines
A = frozenset(int(x) for x in args[1:])
rev = '--reverse' in sys.argv
depth = int(next((a[8:] for a in sys.argv if a.startswith('--depth=')), 12))
secs = float(next((a[10:] for a in sys.argv if a.startswith('--seconds=')), 1200))
U = 2 * math.pi / M
dirs = [i * U for i in range(M)]

def preimages(S):
    """For each grid line j: the grid rays whose image under the folding map of S lies on it."""
    st = rigid.State(NUM, [dirs[i] for i in sorted(S)], None)
    if S and st.phis is None: return None
    out = collections.defaultdict(list)
    for x in range(M):
        v = st.image(dirs[x]) % math.pi
        j = round(v / U) % N
        if abs(v - j * U) < 1e-7 or abs(v - j * U - math.pi) < 1e-7: out[j].append(x)
    return out

sk_cache = {}
def skel(S, D):
    key = (S, D)
    if key not in sk_cache:
        sk_cache[key] = rigid.skeleton(NUM, [dirs[i] for i in sorted(S)], [dirs[i] for i in sorted(D)])
    return sk_cache[key]

t0 = time.monotonic()
start = (frozenset(), ((0,),))
parent = {start: None}
frontier = [start]; found = None; level = 0
while frontier and level < depth and found is None and time.monotonic() - t0 < secs:
    nxt = []
    for state in frontier:
        S, orders = state
        pre = preimages(S)
        if pre is None: continue
        for j, rays in pre.items():
            for size in range(2, len(rays) + 1, 2):
                for D in itertools.combinations(rays, size):
                    D = frozenset(D)
                    sk = skel(S, D)
                    if sk is None: continue
                    S2 = frozenset(i for i in range(M) if any(abs(((dirs[i] - x + math.pi) % (2 * math.pi)) - math.pi) < 1e-7 for x in sk.S2))
                    for L in (0, 1):
                        for up in (0, 1):
                            for sense in (sk.senses[L] if rev else [(0,) * sk.n]):
                                r = rigid.comb_step(sk, orders, L, up, sense)
                                if not isinstance(r, tuple): continue
                                ns = (S2, r)
                                if ns in parent: continue
                                parent[ns] = (state, sorted(D), L, up, any(sense))
                                if S2 == A: found = ns; break
                                nxt.append(ns)
                            if found: break
                        if found: break
                    if found: break
                if found: break
            if found: break
        if found or time.monotonic() - t0 > secs: break
    level += 1
    print(json.dumps({'level': level, 'frontier': len(nxt), 'states': len(parent), 'seconds': round(time.monotonic() - t0, 1)}), flush=True)
    frontier = nxt
if found:
    path, s = [], found
    while parent[s] is not None:
        prev, D, L, up, split = parent[s]; path.append((D, split, sorted(s[0]))); s = prev
    path.reverse()
    print('FOUND in %d steps:' % len(path))
    for D, split, S2 in path:
        made = [d for d in D if d in S2]; erased = [d for d in D if d not in S2]
        print('  step D=%s  made %s  erased %s%s  -> state %s' % (D, made, erased, '  (reverse fold)' if split else '', S2))
else:
    print('NOT FOUND: frontier exhausted' if not frontier else 'NOT FOUND within the depth/time limit', 'states explored', len(parent))
