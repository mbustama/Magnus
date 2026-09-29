"""Compare two paper_figure_cache.json files section by section (from issue #92).
usage: python diff_cache.py COMMITTED.json RECOMPUTED.json"""
import json, sys, math

SKIP = {'fingerprint', 'measured', 'machine', 'seconds'}


def leaves(x, path=''):
    if isinstance(x, dict):
        for k, v in x.items():
            if k not in SKIP:
                yield from leaves(v, path + '/' + str(k))
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from leaves(v, '%s[%d]' % (path, i))
    else:
        yield path, x


a, b = (json.load(open(f)) for f in sys.argv[1:3])
for sec in sorted(set(a) | set(b)):
    if sec not in a or sec not in b:
        print('%-40s only in %s' % (sec, 'committed' if sec in a else 'recomputed'))
        continue
    la, lb = dict(leaves(a[sec])), dict(leaves(b[sec]))
    if la.keys() != lb.keys():
        print('%-40s STRUCTURE differs (%d vs %d leaves)' % (sec, len(la), len(lb)))
        continue
    worst = []
    for p, va in la.items():
        vb = lb[p]
        if isinstance(va, (int, float)) and isinstance(vb, (int, float)) and not isinstance(va, bool):
            if va != vb and not (math.isnan(va) and math.isnan(vb)):
                worst.append((abs(vb - va), abs(vb - va)/max(abs(va), 1e-300), p, va, vb))
        elif va != vb:
            worst.append((math.inf, math.inf, p, va, vb))
    if not worst:
        print('%-40s identical' % sec)
        continue
    worst.sort(key=lambda t: -t[1])
    print('%-40s %d of %d leaves differ; largest relative changes:' % (sec, len(worst), len(la)))
    for d, r, p, va, vb in worst[:3]:
        print('      %-60s %r -> %r  (rel %.1e)' % (p[:60], va, vb, r))
