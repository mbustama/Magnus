import numpy as np
from scipy.optimize import fsolve
import magnus.oscprob as oscprob, magnus.globaldefs as gd
E, L, rho = 2.0*gd.UNIT_GEV, 1300.0*gd.UNIT_KM, 3.0
kw = dict(nu_i=gd.NUMU, nu_f=gd.NUE, density_matter_is_in_g_per_cm3=True)
eps = dict(eps_ee=0.10, eps_em=0.05, eps_mt=0.03)
def A(d): return np.array([oscprob.osc_prob_3nu_matter_constant_density(E,L,rho,nubar=nb,dCP=d,**kw) for nb in (False,True)])
def B(d): return np.array([oscprob.osc_prob_3nu_matter_nsi_constant_density(E,L,rho,nubar=nb,dCP=d,**eps,**kw) for nb in (False,True)])
for g in [(-0.8614,-0.9642),(-0.5753,-0.3035),(-0.0186,0.0257),(0.2392,0.5992)]:
    s = fsolve(lambda x: A(x[0]*np.pi)-B(x[1]*np.pi), g, xtol=1e-12)
    print(s, A(s[0]*np.pi))
