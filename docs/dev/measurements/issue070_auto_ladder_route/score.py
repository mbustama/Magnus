import time, warnings, numpy as np
from scipy.integrate import solve_ivp
warnings.simplefilter('ignore')
import magnus.globaldefs as gd, magnus.oscprob as op, magnus.hamiltonians as hm, magnus.matter as matter
osc = gd.load_nufit_params('NuFIT 6.1')
def run(L_km, rho0, h_km, E_gev):
    L = L_km*gd.UNIT_KM; E = E_gev*gd.UNIT_GEV
    kw = dict(L=L, L0=0.0, rho_central=rho0, l_scale=h_km*gd.UNIT_KM, density_matter_is_in_g_per_cm3=True)
    f = lambda **k: np.asarray(op.osc_prob_3nu_matter_exp_density(E, **kw, **osc, **k))
    # the same Hamiltonian, integrated by DOP853: build it through the package's own H
    info = {}
    hv = np.asarray(hm.hamiltonian_3nu_vacuum(E, **osc), dtype=complex)
    ne0 = rho0*gd.UNIT_G_PER_CM3/(0.5*(gd.MASS_PROTON + gd.MASS_NEUTRON))*0.5
    per = matter.VCC_func(0.0, lambda l: 1.0)
    def H(l): h = hv.copy(); h[0, 0] += per*ne0*np.exp(-l/(h_km*gd.UNIT_KM)); return h
    def rhs(x, y):
        U = y.reshape(3, 3, 2); dU = -1j*H(x) @ (U[..., 0] + 1j*U[..., 1]); return np.stack([dU.real, dU.imag], -1).ravel()
    y0 = np.stack([np.eye(3), np.zeros((3, 3))], -1).ravel()
    sol = solve_ivp(rhs, (0.0, L), y0, method='DOP853', rtol=1e-12, atol=1e-12)
    U = sol.y[:, -1].reshape(3, 3, 2); U = U[..., 0] + 1j*U[..., 1]
    Pref = np.abs(U.T)**2          # P[a, b] = |U_ba|^2, rows = initial flavor
    out = {}
    for name, k in (('auto 1e-3', {}), ('magnus 1e-3', dict(strategy='magnus')), ('magnus 1e-4', dict(strategy='magnus', rtol=1e-4, atol=1e-4))):
        P = f(**k); out[name] = np.max(np.abs(P - Pref))
    chk = np.max(np.abs(f(rtol=1e-9, atol=1e-9, strategy='magnus') - Pref))
    return out, chk
cases = [(25, 3e3, 10, 0.002), (25, 3e3, 10, 0.02), (25, 3e3, 10, 0.2), (250, 3e3, 100, 0.2), (2500, 3e3, 1000, 0.2),
         (2500, 3e3, 1000, 0.02), (25000, 3e3, 10000, 0.02)]
for c in cases:
    out, chk = run(*c)
    print('L=%6g km E=%6g GeV | ' % (c[0], c[3]) + ' | '.join('%s err %.1e' % kv for kv in out.items()) + ' | (ladder 1e-9 vs DOP853: %.1e)' % chk, flush=True)
