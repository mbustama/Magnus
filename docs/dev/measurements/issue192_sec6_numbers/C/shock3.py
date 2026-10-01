import numpy as np, warnings, time
import magnus.oscprob as oscprob, magnus.globaldefs as gd
from shockdef import *
osc = gd.load_nufit_params('NuFIT 6.1')
Ls = np.linspace(10200.0, 80000.0, 4000)*gd.UNIT_KM
E = 15.0*gd.UNIT_MEV
eps = dict(eps_ee=0.10, eps_em=0.05, eps_et=0.0, eps_mm=0.0, eps_mt=0.03, eps_tt=0.0)
osc4 = dict(osc, s14=0.10**0.5, s24=0.10**0.5, s34=0.0, D41=1.0, d14=0.0, d24=0.0)
osc5 = dict(osc4, s15=0.06**0.5, s25=0.06**0.5, s35=0.0, D51=1.7, d15=0.0, d35=0.0)
warnings.simplefilter('ignore')
for wkm in (0.07, 70.0):
    ne = sn_shock_ne(wkm/7e4); old = np.load('scan_%g.npz' % wkm)
    run = dict(L0=L0, t_breakpoints=edges_for(wkm), rtol=1e-8, atol=1e-10, nu_i=gd.NUE, nu_f=gd.NUE, density_is_of_number_of_electrons=True, max_n_slabs=200000)
    t=time.time()
    r = dict(P3=oscprob.osc_prob_matter_std_potential(3, ne, E, Ls, osc, **run),
             NSI=oscprob.osc_prob_matter_nsi(3, ne, E, Ls, osc, eps, **run),
             P4=oscprob.osc_prob_matter_std_potential(4, ne, E, Ls, osc4, **run),
             P5=oscprob.osc_prob_matter_std_potential(5, ne, E, Ls, osc5, **run))
    print('w=%g (%.0fs)' % (wkm, time.time()-t), {k: '%.1e' % np.abs(np.asarray(v,float)-old[k]).max() for k,v in r.items()})
