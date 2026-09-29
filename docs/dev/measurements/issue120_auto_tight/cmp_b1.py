"""B1 before vs after #122: probabilities, engines and warnings per case, configuration and tolerance."""
import json, sys, numpy as np
S = '/tmp/claude-1000/-home-mbustamante-Research-magnus/93bacad5-e69e-4b42-9e80-36335d33ac78/scratchpad'
def load(f):
    out = {}
    for line in open(f):
        r = json.loads(line)
        key = (r['kind'], r['name'], json.dumps(r['E_gev']))
        for row in r['rows']:
            for cfg, v in row.items():
                if cfg != 'tol' and isinstance(v, dict):
                    out[key + (row['tol'], cfg)] = v
    return out
a, b = load(S + '/tight_b1.jsonl'), load(S + '/tight_b1_post122.jsonl')
print('entries: before %d, after %d, common %d' % (len(a), len(b), len(set(a) & set(b))))
same_P = diffP = 0
for k in sorted(set(a) & set(b)):
    x, y = a[k], b[k]
    if x.get('P') is None or y.get('P') is None:
        continue
    if np.array_equal(np.asarray(x['P'], float), np.asarray(y['P'], float)):
        same_P += 1
    else:
        diffP += 1
        print('  P differs: %s %s tol %g %s  max |dP| %.2e' % (k[1], k[2][:30], k[3], k[4], np.abs(np.asarray(x['P'], float) - np.asarray(y['P'], float)).max()))
    if sorted(x.get('warnings') or []) != sorted(y.get('warnings') or []) or x.get('engine') != y.get('engine') or x.get('certified') != y.get('certified'):
        print('  changed: %s %s tol %g %s: %s %s -> %s %s' % (k[1], k[2][:30], k[3], k[4], x.get('engine'), x.get('warnings'), y.get('engine'), y.get('warnings')))
print('probabilities identical: %d, differing: %d' % (same_P, diffP))
