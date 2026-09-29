import sys, time, json, numpy as np
import magnus.oscprob as oscprob
src = compile(open(sys.argv[1]).read(), 'listing', 'exec'); n = int(sys.argv[2])
exec(src, {})
ts = []
for _ in range(n):
    t = time.perf_counter(); exec(src, {}); ts.append(time.perf_counter() - t)
names = ('osc_prob_2nu_matter_exp_density', 'osc_prob_3nu_matter_exp_density', 'osc_prob_4nu_matter_exp_density', 'osc_prob_5nu_matter_exp_density')
engines = []
for name in names:
    f = getattr(oscprob, name)
    def wrap(*a, _f=f, _n=name, **k):
        info = {}; out = _f(*a, strategy_info=info, **k); engines.append((_n.split('_')[2], info.get('engine'))); return out
    setattr(oscprob, name, wrap)
exec(src, {})
print(json.dumps(dict(warm=ts, engines=engines)))
