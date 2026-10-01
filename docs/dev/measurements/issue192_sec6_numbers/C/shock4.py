import numpy as np, warnings, time
import magnus.oscprob as oscprob, magnus.globaldefs as gd, magnus.avgprob as avgprob
from shockdef import *
osc = gd.load_nufit_params('NuFIT 6.1')
Ls = np.linspace(10200.0, 80000.0, 4000)*gd.UNIT_KM
E = 15.0*gd.UNIT_MEV
def call(wkm, declared, **kw):
    ne = sn_shock_ne(wkm/7e4)
    extra = dict(t_breakpoints=edges_for(wkm)) if declared else {}
    with warnings.catch_warnings(record=True) as c:
        warnings.simplefilter('always')
        t=time.time()
        P = oscprob.osc_prob_matter_std_potential(3, ne, E, Ls[-1], osc, L0=10000.0*gd.UNIT_KM, average=True, nu_i=gd.NUE, nu_f=gd.NUE, density_is_of_number_of_electrons=True, **extra, **kw)
    msgs = [(x.category.__name__, str(x.message)) for x in c]
    return float(P), msgs, time.time()-t
for wkm in (70.0, 0.07):
    for dec in (False, True):
        P, m, t = call(wkm, dec)
        print('w=%g declared=%s P=%.4f (%.0fs)' % (wkm, dec, P, t), sorted({a for a,_ in m}))
        for a,b in m:
            if 'standard error' in b: print('   ', b[b.find('largest'):b.find('largest')+70])
P, m, t = call(0.07, True, average_spread=0.05, average_n_samples=161)
print('spread .05 n161 P=%.4f (%.0fs)' % (P, t))
for a,b in m:
    if 'standard error' in b: print('   ', b[b.find('largest'):b.find('largest')+70])
