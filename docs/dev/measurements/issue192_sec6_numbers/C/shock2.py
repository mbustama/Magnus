import numpy as np, warnings, time
warnings.simplefilter('ignore')
import magnus.oscprob as oscprob, magnus.globaldefs as gd
from shockdef import *
osc = gd.load_nufit_params('NuFIT 6.1')
Ls = np.linspace(10200.0, 80000.0, 4000)*gd.UNIT_KM
E = 15.0*gd.UNIT_MEV
eps = dict(eps_ee=0.10, eps_em=0.05, eps_et=0.0, eps_mm=0.0, eps_mt=0.03, eps_tt=0.0)
osc4 = dict(osc, s14=0.10**0.5, s24=0.10**0.5, s34=0.0, D41=1.0, d14=0.0, d24=0.0)
osc5 = dict(osc4, s15=0.06**0.5, s25=0.06**0.5, s35=0.0, D51=1.7, d15=0.0, d35=0.0)
res={}
beyond = Ls/KM > R_FORWARD_KM
for wkm in (0.07, 70.0):
    ne = sn_shock_ne(wkm/7e4)
    run = dict(L0=L0, t_breakpoints=edges_for(wkm), rtol=1e-8, atol=1e-10, nu_i=gd.NUE, nu_f=gd.NUE, density_is_of_number_of_electrons=True)
    t=time.time()
    r = dict(P3=oscprob.osc_prob_matter_std_potential(3, ne, E, Ls, osc, **run),
             NSI=oscprob.osc_prob_matter_nsi(3, ne, E, Ls, osc, eps, **run),
             P4=oscprob.osc_prob_matter_std_potential(4, ne, E, Ls, osc4, **run),
             P5=oscprob.osc_prob_matter_std_potential(5, ne, E, Ls, osc5, **run),
             P2=oscprob.osc_prob_matter_std_potential(2, ne, E, Ls, dict(sth=osc['s12'], Dm2=osc['D21']), **run))
    r = {k: np.asarray(v, dtype=float) for k,v in r.items()}
    res[wkm]=r
    print('w=%g km (%.0f s):' % (wkm, time.time()-t), ' '.join('%s mean beyond FS %.4f' % (k, v[beyond].mean()) for k,v in r.items()), 'P2 min %.4f' % r['P2'].min())
    np.savez('scan_%g.npz' % wkm, **r)
before = Ls/KM < R_CONTACT_KM - 35
for k in res[0.07]:
    a,b=res[0.07][k],res[70.0][k]
    print(k, 'bitwise equal before contact (r<%g): %s; first differing L %.1f km' % (R_CONTACT_KM-35, np.array_equal(a[before],b[before]), Ls[np.argmax(a!=b)]/KM))
# without breakpoints, 3nu sharp
ne = sn_shock_ne(1e-6)
with warnings.catch_warnings(record=True) as c:
    warnings.simplefilter('always')
    Pnb = np.asarray(oscprob.osc_prob_matter_std_potential(3, ne, E, Ls, osc, L0=L0, rtol=1e-8, atol=1e-10, nu_i=gd.NUE, nu_f=gd.NUE, density_is_of_number_of_electrons=True), float)
print('no-bp warnings', sorted({x.category.__name__ for x in c}))
d = np.abs(Pnb - res[0.07]['P3'])
print('no-bp max err beyond FS %.4f, before FS %.2e' % (d[beyond].max(), d[~beyond].max()))
