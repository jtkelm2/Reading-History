"""Write results/corpus_results.md from results/corpus_results.jsonl.
Usage: python3 tools/corpus_table.py"""
import json
from common import ROOT

rows = [json.loads(l) for l in open(ROOT / 'results/corpus_results.jsonl')]
L = ['Produced by `tools/run_corpus.sh` (170 s per pattern, one core, z3 4.16) and `tools/corpus_table.py`.',
     '"Certified" means the cost equals the lower bound, so the resolution is minimal over the whole',
     'search space; otherwise it is minimal over the options given during core-guided expansion.', '',
     '| Pattern | Interior vertices | Unit | Status | Hidden rays | Collapsed vertices | Cost | Lower bound | Certified | Rounds | Time (s) |',
     '|---|---|---|---|---|---|---|---|---|---|---|']
for d in rows:
    ok = d.get('status') == 'ok'
    L.append('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (
        d['file'].rsplit('.', 1)[0], d.get('interior', ''), ('%g°' % (360 / d['N'])) if d.get('N') else ('real' if 'N' in d else ''),
        d['status'], d.get('hidden', '') if ok else '', len(d.get('collapses') or []) if ok else '',
        d.get('cost', '') if ok else '', d.get('lower_bound', '') if ok else '',
        ('yes' if d.get('certified_minimal') else 'no') if ok else '', d.get('rounds', ''), d.get('t', '')))
(ROOT / 'results/corpus_results.md').write_text('\n'.join(L) + '\n')
print('%d patterns, %d resolved, %d certified' % (len(rows), sum(d['status'] == 'ok' for d in rows),
      sum(bool(d.get('certified_minimal')) for d in rows)))
