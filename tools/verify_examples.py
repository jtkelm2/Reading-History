"""Recompute the numbers the paper states for its two worked examples.
Usage: python3 tools/verify_examples.py"""
from common import ROOT, fh, resolve, linear_extensions, geometries
from running_labels import label_running

P, res, (folds, lt, val, paths) = resolve(ROOT / 'examples/example1_running.fold')
reps = list(folds); names = label_running(P, folds); H = fh.reduce(reps, lt)
print('running: hidden rays', res['hidden'], 'folds', len(reps), 'relations', len(lt), 'Hasse', len(H),
      'linear extensions', linear_extensions(reps, lt))
print('  Hasse edges:', sorted((names[a], names[b]) for a, b in H))
print('  incomparable:', [names[x] for x in reps if not any(x in e for e in lt)])
assert (len(reps), len(lt), len(H)) == (15, 30, 19)
G, complete = geometries(P, res)
print('  minimal geometries', len(G), 'complete', complete, 'certified minimal', res['certified_minimal'])
assert complete and len(G) == 1 and res['certified_minimal']

P, res, (folds, lt, val, paths) = resolve(ROOT / 'examples/example2_body.fold')
reps = list(folds)
print('body: hidden', {p: [P.deg(r) for r in R] for p, R in res['choice'].items()}, 'folds', len(reps),
      'relations', len(lt), 'Hasse', len(fh.reduce(reps, lt)), 'linear extensions', linear_extensions(reps, lt))
assert (len(reps), len(lt), linear_extensions(reps, lt)) == (6, 9, 24)
G, complete = geometries(P, res)
print('  minimal geometries', len(G), 'complete', complete, 'certified minimal', res['certified_minimal'])
assert complete and len(G) == 1 and res['certified_minimal']
print('ok')
