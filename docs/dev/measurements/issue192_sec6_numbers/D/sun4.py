from sun3 import *
import json
res={}
t0=time.time()
# ring spacing check at b=0.3, 100 GeV
bs=np.linspace(0.3-5.1e-4, 0.3+5.1e-4, 41); ps=[P(100.0,b) for b in bs]
print('100GeV within one ring spacing of 0.3: %.4f .. %.4f' % (min(ps), max(ps)), flush=True)
for b0 in (0.15, 0.3, 0.45):
    bs=np.linspace(b0, b0+2e-3, 41); ps=[P(1000.0,b) for b in bs]
    print('1TeV core b0=%.2f: P %.4f..%.4f half-range %.4f' % (b0, min(ps), max(ps), (max(ps)-min(ps))/2), flush=True)
for E in (300.0, 1e3, 1e4, 5e4):
    bs=np.arange(0.5, 0.995, 0.002); ps=np.array([P(E,b) for b in bs])
    i=ps.argmin(); bf=np.linspace(bs[i]-0.003, bs[i]+0.003, 25); pf=np.array([P(E,b) for b in bf]); j=pf.argmin()
    print('E=%g deepest ring P=%.4f at b=%.4f (coarse %.4f at %.3f)' % (E, pf[j], bf[j], ps[i], bs[i]), flush=True)
    if E==300.0: print('  averaged_varying at that b: %.4f; its own min over coarse grid: ' % Pvary(E,bf[j]), flush=True)
    res[E]=[bs.tolist(), ps.tolist()]
json.dump(res, open('sun4.json','w'))
print('t', time.time()-t0)
