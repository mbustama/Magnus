import os, sys, json, subprocess, time, statistics as st
S = sys.argv[1] + '/b8'
def load():
    return os.getloadavg()[0]
def fresh(v):
    t = time.perf_counter()
    out = subprocess.run([sys.executable, S + '/fresh.py', S + '/listing_%s.py' % v], capture_output=True, text=True, check=True).stdout
    wall = time.perf_counter() - t
    return dict(json.loads(out.strip().splitlines()[-1]), wall=wall)
res = {'load_start': load()}
for v in ('printed', 'default_tol'):
    res[v + '_first'] = fresh(v)                      # untimed-quality first run: cache fill, reported apart
    runs = [fresh(v) for _ in range(5)]
    res[v + '_fresh'] = {k: st.median(r[k] for r in runs) for k in ('wall', 'imports', 'calls')}
    res[v + '_fresh_all'] = runs
    out = subprocess.run([sys.executable, S + '/warm.py', S + '/listing_%s.py' % v, '5'], capture_output=True, text=True, check=True).stdout
    res[v + '_warm'] = json.loads(out.strip().splitlines()[-1])
    res['load_after_' + v] = load()
out = subprocess.run([sys.executable, S + '/warm.py', S + '/listing_hybrid.py', '3'], capture_output=True, text=True, check=True).stdout
res['hybrid_warm'] = json.loads(out.strip().splitlines()[-1])
res['hybrid_fresh'] = fresh('hybrid')
res['load_end'] = load()
json.dump(res, open(S + '/b8_results.json', 'w'), indent=1)
for v in ('printed', 'default_tol'):
    f, w = res[v + '_fresh'], res[v + '_warm']
    print('%-12s first run: wall %.2f s (calls %.2f) | fresh, median of 5: wall %.2f s = start-up+imports %.2f + calls %.2f | warm, median of 5: %.3f s | engines %s'
          % (v, res[v + '_first']['wall'], res[v + '_first']['calls'], f['wall'], f['wall'] - f['calls'], f['calls'], st.median(w['warm']), sorted(set(e for _, e in w['engines']))))
hw = res['hybrid_warm']
print('hybrid (old default route) at 1e-12, order 8: warm median of 3 %.2f s, fresh wall %.2f s, engines %s' % (st.median(hw['warm']), res['hybrid_fresh']['wall'], sorted(set(e for _, e in hw['engines']))))
print('load: start %.2f, after printed %.2f, after default %.2f, end %.2f' % (res['load_start'], res['load_after_printed'], res['load_after_default_tol'], res['load_end']))
