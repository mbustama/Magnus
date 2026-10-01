import numpy as np
import magnus.oscprob as oscprob
import magnus.hamiltonians as ham
import magnus.earth as earth
import magnus.globaldefs as gd

osc = gd.load_nufit_params('NuFIT 6.1')
E = np.logspace(3, 7, 60)*gd.UNIT_GEV   # 1 TeV-10 PeV
f_src = np.array([1/3, 2/3, 0.0])       # pion decay
# A source 100 Mpc away: far enough for every pair
# of eigenvalues to have decohered at every energy,
# which the wrappers decide from the baseline
L_src = 100.0*3.0857e19*gd.UNIT_KM

def at_earth(P):
    # The source fractions, zero for a sterile state,
    # contracted with the averaged matrix; the active
    # fractions renormalized to one
    f = np.zeros(len(P)); f[:3] = f_src
    f = (f @ P)[:3]
    return f/f.sum()

# Standard: the averaged matrix at each energy
P_std = oscprob.osc_prob_3nu_vacuum(E, L_src,
 average=True, **osc)                   # (60, 3, 3)

# LIV: the n = 1 operator, aligned with the PMNS
# angles, sized to equal the vacuum term at 100 TeV
b3 = osc['D31']/(2*(100.0*gd.UNIT_TEV)**2)
P_liv = oscprob.osc_prob_3nu_vacuum_liv(E, L_src,
 average=True, sxi12=osc['s12'], sxi23=osc['s23'],
 sxi13=osc['s13'], dxiCP=0.0, b1=0.0, b2=0.0,
 b3=b3, Lambda=1.0, n_liv=1, **osc)

# NSI through the Earth: the flux arrives decohered,
# then crosses the Earth coherently, so the two legs
# compose as probability matrices.  The ladder starts
# at 32 slabs: the matter term alone winds more than
# pi across a coarser one
eps = dict(eps_ee=0.10, eps_em=0.05, eps_mt=0.03)
L_in = earth.distance_traveled_inside_earth(-1.0)
P_in = oscprob.osc_prob_3nu_earth_nsi(E, costhz=-1.0,
 L=L_in*gd.UNIT_KM, n_slabs=32, rtol=1e-6, atol=1e-8,
 **osc, **eps)
P_nsi = np.asarray(P_std) @ np.asarray(P_in)

# 3+1: one sterile state, sin^2 of both active-
# sterile angles 0.1, Dm41^2 = 1 eV^2
P_4 = oscprob.osc_prob_4nu_vacuum(E, L_src,
 average=True, s14=np.sqrt(0.1), s24=np.sqrt(0.1),
 D41=1.0, **osc)

# Pseudo-Dirac: the second mass state is paired
# with a partner split by 1e-13 eV^2
U = ham.pmns_mixing_matrix(osc['s12'], osc['s23'],
 osc['s13'], osc['dCP'])
m2 = np.array([0.0, osc['D21'], osc['D31']])
P_pd = oscprob.osc_prob_pseudo_dirac_vacuum(E,
 L_src, {1: 1.0e-13}, average=True, **osc)

cases = dict(std=P_std, liv=P_liv, nsi=P_nsi,
             four=P_4, pd=P_pd)
f_earth = {k: np.array([at_earth(P) for P in v])
           for k, v in cases.items()}
for k, v in f_earth.items():
    print(k, 'min', v.min(0).round(4), 'max', v.max(0).round(4), 'range of each', (v.max(0)-v.min(0)).round(5))
print('b3', b3)
std0 = f_earth['std'][0]
print('std', std0, 'nsi max shift', np.abs(f_earth['nsi']-f_earth['std']).max())
# sterile loss
P4a = np.asarray(P_4); Ppd = np.asarray(P_pd)
f = np.zeros(4); f[:3] = f_src; print('3+1 sterile share', (f@P4a[0])[3:].sum(), 1-(f@P4a[0])[:3].sum())
f = np.zeros(Ppd.shape[-1]); f[:3] = f_src; print('pd sterile share', 1-(f@Ppd[0])[:3].sum(), Ppd.shape)
print('pd fractions', f_earth['pd'][0], 'flat', np.ptp(f_earth['pd'],0))
liv = f_earth['liv']; Pliv=np.asarray(P_liv)
Em = E/gd.UNIT_TEV
print('LIV f_e: first %.4f max %.4f at %.0f TeV min %.4f at %.0f TeV last %.4f' % (liv[0,0], liv[:,0].max(), Em[liv[:,0].argmax()], liv[:,0].min(), Em[liv[:,0].argmin()], liv[-1,0]))
print('LIV f_mu: first %.4f min %.4f at %.0f TeV last %.4f' % (liv[0,1], liv[:,1].min(), Em[liv[:,1].argmin()], liv[-1,1]))
for a in range(3): print('LIV Paa', a, 'max %.4f min %.4f' % (Pliv[:,a,a].max(), Pliv[:,a,a].min()))
for a in range(3): print('std Paa', a, '%.4f' % np.asarray(P_std)[0,a,a])
# 3+1 renorm with f/f.sum done
np.save('liv.npy', liv)
nsi=f_earth['nsi']
for i in (0,10,20,30,59): print('nsi E=%.1f TeV' % Em[i], nsi[i].round(4), 'shift', (nsi[i]-std0).round(4))
