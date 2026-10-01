import numpy as np
import magnus.oscprob as oscprob
import magnus.globaldefs as gd
import magnus.hamiltonians as H
osc = gd.load_nufit_params('NuFIT 6.1')
n, width = 24, 250.0*gd.UNIT_KM
rho = np.where(np.arange(n) % 2 == 0, 2.0, 8.0)
edges = np.arange(n + 1)*width
def castle_wall(l):
    k = np.searchsorted(edges, l, side='right') - 1
    return rho[np.clip(k, 0, n - 1)]
E = np.logspace(-0.7, 1.7, 400)*gd.UNIT_GEV
P, Pbar = [oscprob.osc_prob_matter_std_potential(
    3, castle_wall, E, n*width, osc,
    t_breakpoints=edges[1:-1], nubar=nubar,
    nu_i=gd.NUMU, nu_f=gd.NUE,
    density_matter_is_in_g_per_cm3=True,
    rtol=1.0e-8, atol=1.0e-10)
    for nubar in (False, True)]
EG = E/gd.UNIT_GEV
profs = {'castle': rho, 'serrated': np.tile(np.linspace(2,8,6),4),
         'random': np.random.default_rng(20260801).permutation(rho), 'uniform': np.full(n,5.0)}
curves = {}
for name, r in profs.items():
    def f(l, r=r):
        k = np.searchsorted(edges, l, side='right') - 1
        return r[np.clip(k, 0, n - 1)]
    for nb in (False, True):
        curves[name, nb] = np.asarray(oscprob.osc_prob_matter_std_potential(3, f, E, n*width, osc,
            t_breakpoints=edges[1:-1], nubar=nb, nu_i=gd.NUMU, nu_f=gd.NUE,
            density_matter_is_in_g_per_cm3=True, rtol=1e-8, atol=1e-10))
print('listing castle == loop castle', np.max(np.abs(np.asarray(P)-curves['castle',False])), np.max(np.abs(np.asarray(Pbar)-curves['castle',True])))
for nb in (False, True):
    ref = curves['uniform', nb]
    for name in profs:
        p = curves[name, nb]; k = np.argmax(p)
        d = np.abs(p-ref)
        print(nb, name, 'peak %.4f at %.3f GeV; max|P-Pu| %.4f at %.3f' % (p[k], EG[k], d.max(), EG[np.argmax(d)]))
# local peaks of castle below 1.5 GeV
for nb in (False, True):
    p = curves['castle', nb]; u = curves['uniform', nb]
    m = EG < 1.5
    k = np.argmax(np.where(m, p, -1))
    print('castle low-E peak', nb, '%.4f at %.4f GeV, uniform there %.4f ratio %.3f' % (p[k], EG[k], u[k], p[k]/u[k]))
    # also local maxima list
    idx = [i for i in range(1,len(p)-1) if p[i]>p[i-1] and p[i]>p[i+1] and EG[i]<1.5]
    print('  local maxima', [(round(EG[i],3), round(p[i],4)) for i in idx])
pk_nu = curves['castle',False][EG<1.5].max(); pk_nb = curves['castle',True][EG<1.5].max()
print('ratio nubar/nu low peak', pk_nb/pk_nu)
# oscillation lengths in matter at 0.46 GeV
import inspect
