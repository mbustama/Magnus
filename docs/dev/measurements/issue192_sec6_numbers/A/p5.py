import numpy as np, warnings
warnings.simplefilter('ignore')
import magnus.globaldefs as gd, magnus.oscprob as o
E, L = 1.0*gd.UNIT_GEV, 1000.0*gd.UNIT_KM
Es = np.logspace(-1, 1, 200)*gd.UNIT_GEV
Ls = np.linspace(1, 1300, 200)*gd.UNIT_KM
ex = {2: dict(sth=0.5557, Dm2=7.537e-5), 3: {}, 4: dict(s14=0.3, s24=0.3, D41=1.0), 5: dict(s14=0.3, s24=0.3, s15=0.2, s25=0.2, D41=1.0, D51=1.7)}
for n in (2,3,4,5):
    for env in ('vacuum','matter_constant_density'):
        f = getattr(o, f'osc_prob_{n}nu_{env}')
        kw = dict(ex[n])
        if env!='vacuum': kw.update(rho=3.0, density_matter_is_in_g_per_cm3=True)
        Pb = f(Es, L, **kw); Pp = np.array([f(e, L, **kw) for e in Es])
        Pb2 = f(E, Ls, **kw); Pp2 = np.array([f(E, l, **kw) for l in Ls])
        print(n, env, 'E-batch maxdiff', np.abs(Pb-Pp).max(), 'L-batch maxdiff', np.abs(Pb2-Pp2).max())
