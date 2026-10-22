"""Build the degree-4 encyclopedia page: python3 docs/encyclopedia/build.py EXPORT.json OUT.html
EXPORT.json comes from tools/encyclopedia_export.py; the page inlines it with core.js."""
import sys, json, os
here = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(sys.argv[1]))
for s in D['strata']:
    for k in ('P0_exact', 'P1_exact', 'point_exact'): s.pop(k, None)
    for c in s['cards']:
        c.pop('iv_exact', None); c.pop('id', None)
D.pop('checks', None)
page = open(os.path.join(here, 'page.html')).read()
core = open(os.path.join(here, 'core.js')).read()
data = json.dumps(D, separators=(',', ':')).replace('</', '<\\/')
out = page.replace('__CORE__', core).replace('__DATA__', data)
open(sys.argv[2], 'w').write(out)
print('wrote', sys.argv[2], len(out) // 1024, 'KB')
