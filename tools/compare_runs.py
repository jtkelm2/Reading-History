"""Compare two corpus runs (JSON lines from tools/run_corpus.sh), e.g. the local search against
the lookup tables.  Prints a per-pattern table (Markdown) and a summary.
Usage: python3 tools/compare_runs.py BASE.jsonl NEW.jsonl [BASE_NAME NEW_NAME]"""
import json, sys, statistics

A = {json.loads(l)['file']: json.loads(l) for l in open(sys.argv[1]) if l.strip()}
B = {json.loads(l)['file']: json.loads(l) for l in open(sys.argv[2]) if l.strip()}
na, nb = (sys.argv[3], sys.argv[4]) if len(sys.argv) > 4 else ('base', 'new')


def why(d):
    if d.get('status') == 'ok': return 'ok'
    if d.get('status') == 'local fail':
        deg = d.get('first_degree')
        return 'no local history (degree %s%s)' % (deg, ', not tabulated' if deg and deg > 6 else '') if deg else 'no local history'
    return d.get('status', '?')


def cell(d):
    if d.get('status') != 'ok': return why(d)
    return 'cost %s%s, %.1f s' % (d['cost'], '' if d.get('certified_minimal') else '*', d.get('t', 0))


print('| Pattern | %s | %s |' % (na, nb)); print('|---|---|---|')
for f in A:
    print('| %s | %s | %s |' % (f.rsplit('.', 1)[0], cell(A[f]), cell(B.get(f, {'status': 'missing'}))))
ok = lambda D: {f for f, d in D.items() if d.get('status') == 'ok'}
both = ok(A) & ok(B)
print()
print('resolved: %s %d, %s %d, both %d, only %s %d, only %s %d' % (na, len(ok(A)), nb, len(ok(B)), len(both),
      na, len(ok(A) - ok(B)), nb, len(ok(B) - ok(A))))
same = sum(A[f]['cost'] == B[f]['cost'] for f in both)
print('on patterns resolved by both: same cost %d of %d; %s cheaper %d; %s cheaper %d' % (
    same, len(both), na, sum(A[f]['cost'] < B[f]['cost'] for f in both), nb, sum(B[f]['cost'] < A[f]['cost'] for f in both)))
for nm, D in ((na, A), (nb, B)):
    ts = [D[f].get('t', 0) for f in both]
    if ts: print('%s time on those: total %.1f s, median %.2f s, max %.1f s' % (nm, sum(ts), statistics.median(ts), max(ts)))
reasons = {}
for f, d in B.items():
    if d.get('status') != 'ok': reasons.setdefault(why(d), []).append(f.rsplit('.', 1)[0])
print('%s failures by reason:' % nb)
for r, fs in sorted(reasons.items()): print('  %s: %d  %s' % (r, len(fs), ', '.join(fs)))
