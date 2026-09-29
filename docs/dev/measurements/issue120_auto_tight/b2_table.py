"""#120 B2: the proposed rule applied to B1's records, under each reading of 'tol'."""
import sys, json, numpy as np
sys.path.insert(0, '/home/mbustamante/Research/magnus/resources/handover_auto_tight')
import tight_measure as tm
R = tm.refs()
def err_lim(r, x, tol):
    keys = [(r['name'], r['E_gev'])] if r['kind'] == 'point' else [(r['name'][:-5], E) for E in r['E_gev']]
    out = [(abs(P - R[k][1]), max(tol*1e-2 + tol*abs(R[k][1]), abs(R[k][0] - R[k][1]))) for P, k in zip(x['P'], keys) if k in R]
    return max(out, key=lambda e: e[0]/e[1]) if out else None
def limit(tol, p, cap):
    return min(cap, 1e4*(tol/1e-6)**(1.0/p))
for path in sys.argv[1:]:
    recs = [json.loads(l) for l in open(path)]
    print('\n##### %s' % path.split('/')[-1])
    for reading in ('min(rtol, atol)', 'rtol'):
        for cap in (2000.0, np.inf):
            print('\n--- tol = %s, cap %s' % (reading, cap))
            print('%-4s %-6s %6s %7s %6s %7s %9s %8s %6s %5s' % ('p', 'rtol', 'limit', 'routed', 'silent', 'H_silent', 'newwarn', 'median', 'max', 'n>1.2'))
            for p in (4, 8):
                for rt in tm.TOLS:
                    tol = rt*1e-2 if reading.startswith('min') else rt
                    lim = limit(tol, p, cap)
                    routed = silent = hsilent = newwarn = 0; ratios = []; notes = []
                    for r in recs:
                        if not r['eligible']:
                            continue
                        row = {x['tol']: x for x in r['rows']}[rt]
                        L, H = row.get('L%d/1' % p), row['H%d' % p]
                        ph = r['route_1e-3']['phase']
                        if ph > lim or L is None or 'error' in L:
                            continue
                        routed += 1
                        el, eh = err_lim(r, L, rt), err_lim(r, H, rt)
                        lw, hw = bool(tm.HONEST & set(L['warnings'])), bool(tm.HONEST & set(H['warnings']))
                        if el and not lw and el[0] > el[1]:
                            silent += 1; notes.append('silent %.2fx %s %s' % (el[0]/el[1], r['name'], str(r['E_gev'])[:8]))
                        if eh and not hw and eh[0] > eh[1]:
                            hsilent += 1
                        if lw and not hw:
                            newwarn += 1; notes.append('newwarn %s %s phase %.0f' % (r['name'], str(r['E_gev'])[:8], ph))
                        ratios.append(L['t']/H['t'])
                    rs = np.array(ratios)
                    print('%-4d %-6g %6.0f %7d %6d %7d %9d %8.2f %6.2f %5d  %s' % (p, rt, lim, routed, silent, hsilent, newwarn, np.median(rs), rs.max(), (rs > 1.2).sum(), '; '.join(notes[:4])))
