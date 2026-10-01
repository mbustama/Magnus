import numpy as np, warnings
warnings.simplefilter('ignore')
from scipy.integrate import solve_ivp
import magnus.oscprob as oscprob, magnus.globaldefs as gd, magnus.solarmodels as sm
import magnus.matter as matter, magnus.hamiltonians as hamiltonians
MODEL='BS05-AGS-OP'
NUFIT_NO = gd.load_nufit_params('NuFIT 6.1', 'NO')
ne_bs05 = sm.electron_density_profile(MODEL)
ENERGY = 5.0*gd.UNIT_MEV
L0, L1 = 0.0, gd.L_SCALE_SUN
print('L1/RS', L1/(gd.SUN_RADIUS*gd.UNIT_KM))
params2 = {'sth': NUFIT_NO['s12'], 'Dm2': NUFIT_NO['D21']}
hvac2 = hamiltonians.hamiltonian_2nu_vacuum_energy_independent(params2['sth'], params2['Dm2'])
e00 = np.diag([1.0, 0.0])
VCC = matter.vcc_func_from_rho_func(ne_bs05, 0.0, 1.0, 0.5, nubar=False, density_matter_is_in_g_per_cm3=False, density_is_of_number_of_electrons=True)
def H_of_l(l): return (1.0/ENERGY)*hvac2 + np.asarray(VCC(l))[..., None, None]*e00
def exact_U_many(H_func, l0, Ls, dim):
    def rhs(l, y): return (-1j*np.asarray(H_func(l)) @ y.reshape(dim, dim)).ravel()
    sol = solve_ivp(rhs, (float(l0), float(Ls[-1])), np.eye(dim, dtype=complex).ravel(), rtol=1e-12, atol=1e-14, method='DOP853', t_eval=Ls)
    return np.array([sol.y[:, i].reshape(dim, dim) for i in range(len(Ls))])
def to_P(U): return np.swapaxes(np.asarray(U).real**2 + np.asarray(U).imag**2, -1, -2)
L_OSC = 4.0*np.pi*ENERGY/params2['Dm2']
avg = float(np.asarray(oscprob.osc_prob_2nu_sun(ENERGY, L1, L0, **params2, density_profile=MODEL, average=True))[0][0])
avgd = float(np.asarray(oscprob.osc_prob_2nu_sun(ENERGY, L1, L0, **params2, density_profile=MODEL, average=True, average_initial_state='decohered'))[0][0])
print('average=True flavor %.6f decohered %.6f' % (avg, avgd))
for n in (6,12,24,48):
    Ls_w = np.linspace(L1 - n*L_OSC, L1, 20*n + 1)
    P_w = np.array([to_P(U) for U in exact_U_many(H_of_l, L0, Ls_w, 2)])
    m = P_w[:,0,0].mean(); print(n, '%.6f bias vs decoh %.4f vs flavor %.4f' % (m, m-avgd, m-avg))
TH = np.arcsin(params2['sth'])
def c2m(l, e):
    x = 2.0*e*float(np.asarray(VCC(l)))/params2['Dm2']; return np.cos(np.arctan2(np.sin(2*TH), np.cos(2*TH)-x))
worst=0
for Em in (1.,2.,5.,8.,10.,15.,20.):
    e=Em*gd.UNIT_MEV
    got=np.asarray(oscprob.osc_prob_2nu_sun(e, L1, L0, **params2, density_profile=MODEL, average=True, average_initial_state='decohered'))[0][0]
    worst=max(worst, abs(got-(0.5+0.5*c2m(L0,e)*c2m(L1,e))))
print('worst decohered vs textbook %.1e' % worst)
