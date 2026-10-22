"""Summarise corpus runs (JSON lines from tools/run_corpus.sh) as the Markdown table of
docs/lookup_tables.md: per run, how many patterns are resolved, fail locally or globally, or time
out, and the median and largest time of the resolved patterns.
Usage: python3 tools/lookup_summary.py RUN.jsonl [...]"""
import sys, json, os
from collections import Counter

FLAGS = {'search': 'local search (no tables)', 'tables_anchored': '`--table` (anchored types only)',
         'tables_dl': '`--table --dl`', 'rigid_anchored': '`--rigid=halfturn --cells` (anchored types only)',
         'physical_anchored': '`--rigid=reverse --cells` (anchored types only)',
         'physical': '`--rigid=reverse --cells`', 'physical_dl': '`--rigid=reverse --cells --dl`'}
runs = []
for path in sys.argv[1:]:
    name = os.path.basename(path)[:-6]
    runs.append((name, {json.loads(l)['file']: json.loads(l) for l in open(path) if l.strip()}))
print('| run | options | resolved | local fail | global unsat | timeout | time of resolved: median / max |')
print('|---|---|---|---|---|---|---|')
for name, R in runs:
    c = Counter(d['status'] for d in R.values())
    ts = sorted(d.get('t', 0) for d in R.values() if d['status'] == 'ok')
    med = (ts[(len(ts) - 1) // 2] + ts[len(ts) // 2]) / 2 if ts else 0
    print('| `%s` | %s | %d | %d | %d | %d | %.1f / %.1f s |' % (name, FLAGS.get(name, ''), c['ok'], c['local fail'], c['global unsat'], c['timeout'], med, ts[-1] if ts else 0))
