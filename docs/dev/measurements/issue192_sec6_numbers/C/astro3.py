import numpy as np, warnings
import magnus.oscprob as oscprob, magnus.earth as earth, magnus.globaldefs as gd
osc = gd.load_nufit_params('NuFIT 6.1')
E = np.logspace(3, 7, 60)*gd.UNIT_GEV
eps = dict(eps_ee=0.10, eps_em=0.05, eps_mt=0.03)
L_in = earth.distance_traveled_inside_earth(-1.0)
out={}
for ns in (32, None):
    kw = {} if ns is None else dict(n_slabs=ns)
    with warnings.catch_warnings(record=True) as C:
        warnings.simplefilter('always')
        out[ns] = np.asarray(oscprob.osc_prob_3nu_earth_nsi(E, costhz=-1.0, L=L_in*gd.UNIT_KM, rtol=1e-6, atol=1e-8, **kw, **osc, **eps))
    print(ns, sorted({w.category.__name__ for w in C}))
    for w in C[:2]: print('   ', str(w.message)[:250])
print('max change', np.abs(out[32]-out[None]).max())
L8 = 1e8*gd.UNIT_KM
with warnings.catch_warnings(record=True) as C2:
    warnings.simplefilter('always')
    P8 = oscprob.osc_prob_3nu_vacuum(E, L8, average=True, **osc)
print(str(C2[0].message))
