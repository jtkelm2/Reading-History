"""Regenerate the LaTeX fragments that are pasted into paper/reading_history.tex:
the appendix table (running pattern), the corpus table rows, and the crossover TikZ.
Usage: python3 tools/paper_tables.py   (writes into build/)"""
import json, math
from common import ROOT, fh, resolve
from running_labels import label_running

out = ROOT / 'build'; out.mkdir(exist_ok=True)

# appendix: local histories at every interior vertex of the running pattern
P, res, (folds, lt, val, paths) = resolve(ROOT / 'examples/example1_running.fold')
names = label_running(P, folds); where = {k: x for x, ks in folds.items() for k in ks}
u = math.pi / 32; rows = []
for p in P.interior:
    R = res['choice'].get(p, ())
    hs = next(o[1] for o in res['opts'][p] if o[0] == R)
    D = sorted([round(d / u) % 64 for d, _ in P.H[p]] + [round(r / u) % 64 for r in R for _ in (0, 1)])
    seqs = [['/'.join(sorted({names[where[k]] for k in st})) for st in h] for h in hs]
    rows.append((P.V[p], D, seqs))
lines = ['$(%.3f, %.3f)$ & %s & %s \\\\' % (x, y, ', '.join(map(str, D)),
         '; '.join(' $\\to$ '.join('$' + s + '$' for s in q) for q in seqs))
         for (x, y), D, seqs in sorted(rows, key=lambda r: (-r[0][1], r[0][0]))]
(out / 'appendix_rows.tex').write_text('\n'.join(lines) + '\n')

# corpus rows from results/corpus_results.jsonl
rows = [json.loads(l) for l in open(ROOT / 'results/corpus_results.jsonl')]
AUTHORS = {'traditional': 'Traditional', 'traditionaloripa': 'Traditional', 'ku': 'Ku', 'komatsu': 'Komatsu',
           'komatsuoripa': 'Komatsu', 'tanakaoripa': 'Tanaka', 'tanaka': 'Tanaka', 'tachioripa': 'Tachi',
           'furutaoripa': 'Furuta', 'kei': 'Kei', 'boxhard': 'boxhard', 'bernhayes': 'Bern--Hayes', 'miuraoripa': 'Miura'}
def short(f):
    n, a, *t = f.rsplit('.', 1)[0].split('_')
    t = ' '.join(t).replace('(Pokemon)', '').replace('45˚', '45$^\\circ$').strip()
    return '%s & %s, %s' % (n, AUTHORS.get(a, a), t)
groups = [('S', 'no cancellations'), ('H', 'cancellations at vertices with no single-line history'),
          ('X', 'cancellations also at vertices that have one'), ('C', 'collapses'), ('U', 'unresolved within the time budget')]
def cls(d):
    if d['status'] != 'ok': return 'U'
    if d.get('collapses'): return 'C'
    if d.get('rounds', 1) > 1: return 'X'
    return 'H' if d.get('hidden', 0) else 'S'
unit = lambda d: ('%g$^\\circ$' % (360 / d['N'])) if d.get('N') else '---'
lines = []
for c, title in groups:
    g = [d for d in rows if cls(d) == c]
    lines.append('\\midrule\\multicolumn{6}{l}{\\emph{%s} (%d)}\\\\' % (title, len(g)))
    for d in g:
        star = '' if c == 'U' or d.get('certified_minimal') else '$^*$'
        lines.append('%s & %s & %s & %s & %s \\\\' % (short(d['file']), d.get('interior') if d.get('interior') is not None else '',
                     unit(d), (str(d.get('hidden', '')) + star) if c != 'U' else '', len(d.get('collapses') or []) if c != 'U' else ''))
(out / 'corpus_rows.tex').write_text('\n'.join(lines) + '\n')

# crossover TikZ
Q = fh.Pattern(str(ROOT / 'examples/bernhayes_45_crossover.fold'))
xs = [q[0] for q in Q.V]; ys = [q[1] for q in Q.V]; s = 7.0 / max(max(xs) - min(xs), max(ys) - min(ys))
c = lambda p: ((Q.V[p][0] - min(xs)) * s, (Q.V[p][1] - min(ys)) * s)
L = ['\\draw[sheet] (%.3f,%.3f) -- (%.3f,%.3f);' % (*c(u_), *c(v_)) for u_, v_ in Q.bound]
L += ['\\draw[crease] (%.3f,%.3f) -- (%.3f,%.3f);' % (*c(u_), *c(v_)) for u_, v_, _ in Q.crease]
L += ['\\fill (%.3f,%.3f) circle (1.6pt);' % c(p) for p in Q.interior]
(out / 'crossover.tex').write_text('\n'.join(L) + '\n')
print('wrote', sorted(p.name for p in out.iterdir()))
