"""Which calls the proposed rule decides differently from main, in the suite and the notebooks."""
import json, sys, collections, os
S = '/tmp/claude-1000/-home-mbustamante-Research-magnus/93bacad5-e69e-4b42-9e80-36335d33ac78/scratchpad'
for f in ('route120_tests.jsonl', 'route120_nb.jsonl'):
    p = os.path.join(S, f)
    if not os.path.exists(p):
        print(f, 'missing'); continue
    rows = [json.loads(l) for l in open(p) if l.strip()]
    ch = [r for r in rows if r['old'] != r['new']]
    loose = [r for r in rows if r['tol'] >= 1e-6]
    print('\n%s: %d decisions, %d at >= 1e-6 (changed %d), %d below 1e-6, CHANGED %d'
          % (f, len(rows), len(loose), sum(r['old'] != r['new'] for r in loose), len(rows) - len(loose), len(ch)))
    by = collections.defaultdict(list)
    for r in ch:
        by[r['where'].split(' (')[0]].append(r)
    for w, rs in sorted(by.items()):
        print('  %s: %d changed' % (w, len(rs)))
        seen = set()
        for r in rs:
            k = (r['tol'], r['p'], round(r['phase'], 1), r['n_energies'], r['batched'])
            if k in seen: continue
            seen.add(k)
            print('      tol %.0e p=%d %s phase %.1f (limit %.0f) floor %d  %d energies%s  %s -> %s'
                  % (r['tol'], r['p'], r['method'], r['phase'], r['max_phase'], r['n_floor'], r['n_energies'],
                     ' batched' if r['batched'] else '', 'ladder' if r['old'] else 'hybrid', 'ladder' if r['new'] else 'hybrid'))
