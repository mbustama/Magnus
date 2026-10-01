import numpy as np, warnings, time, sys, json
import magnus.oscprob as oscprob
import magnus.globaldefs as gd
from shockdef import *
osc = gd.load_nufit_params('NuFIT 6.1')
ne_shock = sn_shock_ne(1e-6)
w = 0.07*gd.UNIT_KM
edges = [r + s*w/2 for r in (12348.0*gd.UNIT_KM, 30323.0*gd.UNIT_KM) for s in (-1.0, 1.0)]
Ls = np.linspace(10200.0, 80000.0, 4000)*gd.UNIT_KM
run = dict(L0=10000.0*gd.UNIT_KM, t_breakpoints=edges, rtol=1.0e-8, atol=1.0e-10, nu_i=gd.NUE, nu_f=gd.NUE, density_is_of_number_of_electrons=True)
E = 15.0*gd.UNIT_MEV
info = {}
t=time.time()
with warnings.catch_warnings(record=True) as c:
    warnings.simplefilter('always')
    P3 = oscprob.osc_prob_matter_std_potential(3, ne_shock, E, Ls, osc, strategy_info=info, **run)
print('engine', info.get('engine'), 'time %.0f' % (time.time()-t))
for x in c: print('  W', x.category.__name__, str(x.message)[:200])
np.save('P3_sharp.npy', np.asarray(P3))
