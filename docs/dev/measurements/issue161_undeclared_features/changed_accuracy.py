import sys, pickle, warnings, numpy as np
warnings.simplefilter('ignore')
sys.path.insert(0, 'src')
import magnus.globaldefs as gd, magnus.hamiltonians as hm, magnus.matter as matter
from scipy.integrate import solve_ivp
KM = gd.UNIT_KM; GEV = gd.UNIT_GEV
p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']; std = {k: p[k] for k in ('s12','s23','s13','dCP','D21','D31')}
hv3 = np.asarray(hm.hamiltonian_3nu_vacuum_energy_independent(**std)); e00 = np.diag([1., 0, 0])
Vf = lambda rho: np.sqrt(2)*gd.GF*rho*gd.UNIT_G_PER_CM3/((gd.MASS_PROTON + gd.MASS_NEUTRON)/2)*0.5
Es = np.geomspace(0.5, 20, 12)*GEV
def dop(Hl, L, splits):
    psi = np.eye(3, dtype=complex); pts = [0.0] + sorted(s for s in splits if 0 < s < L) + [L]
    for a, b in zip(pts[:-1], pts[1:]):
        s = solve_ivp(lambda l, y: (-1j*Hl(l) @ y.reshape(3, 3)).ravel(), (a, b), psi.ravel(), method='DOP853', rtol=1e-12, atol=1e-14)
        psi = s.y[:, -1].reshape(3, 3)
    return (np.abs(psi)**2).T
def std_ref(rho, E, L, splits):
    vcc = matter.vcc_func_from_rho_func(rho, density_matter_is_in_g_per_cm3=True)
    return dop(lambda l: hv3/E + float(vcc(l))*e00, L, splits)
step = lambda x: np.where(np.asarray(x, float) < 2000*KM, 3.0, 8.0)
sp5 = lambda x: 3 + 50*np.exp(-((np.asarray(x, float) - 1000*KM)/(5*KM))**2)
sp50 = lambda x: 3 + 50*np.exp(-((np.asarray(x, float) - 1000*KM)/(50*KM))**2)
castle = lambda E: (lambda l: hv3/E + Vf(3.3 if (l < 3000*KM or l > 9000*KM) else 11.5)*e00)
refs = {
 'std step magnus': std_ref(step, 3*GEV, 5000*KM, [2000*KM]),
 'std spike5 single': std_ref(sp5, 3*GEV, 4000*KM, [960*KM, 1040*KM]),
 'std spike5 escan': np.array([std_ref(sp5, E, 4000*KM, [960*KM, 1040*KM]) for E in Es[:6]]),
 'std spike5 lscan': np.array([std_ref(sp5, 3*GEV, L, [960*KM, 1040*KM]) for L in np.linspace(0.2, 1, 6)*4000*KM]),
 'std spike50 single': std_ref(sp50, 3*GEV, 4000*KM, [600*KM, 1400*KM]),
 'std spike50 escan': np.array([std_ref(sp50, E, 4000*KM, [600*KM, 1400*KM]) for E in Es[:6]]),
 'std spike50 lscan': np.array([std_ref(sp50, 3*GEV, L, [600*KM, 1400*KM]) for L in np.linspace(0.2, 1, 6)*4000*KM]),
 'raw castle single': dop(castle(3*GEV), 6200*KM, [3000*KM, 9000*KM]),
 'raw castle escan': np.array([dop(castle(E), 6200*KM, [3000*KM, 9000*KM]) for E in Es[:5]]),
 'raw castle tight': dop(castle(3*GEV), 6200*KM, [3000*KM, 9000*KM]),
 'raw castle nocumul': np.array([dop(castle(3*GEV), L, [3000*KM, 9000*KM]) for L in np.linspace(1000, 11000, 4)*KM]),
 'raw castle strict': dop(castle(3*GEV), 6200*KM, [3000*KM, 9000*KM]),
}
a = pickle.load(open(sys.argv[1], 'rb')); b = pickle.load(open(sys.argv[2], 'rb'))
tol = {'raw castle tight': 1e-7}
for k, ref in refs.items():
    em = np.abs(a[k][0] - ref).max(); eb = np.abs(b[k][0] - ref).max()
    print('%-22s error main %.1e -> branch %.1e  (tolerance %.0e)%s' % (k, em, eb, tol.get(k, 1e-3), '' if eb <= em else '   <-- WORSE'))
