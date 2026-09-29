import time, warnings, numpy as np
from scipy.integrate import solve_ivp
import magnus.oscprob as op, magnus.globaldefs as gd, magnus.adiabatic as ad
osc = gd.load_nufit_params('NuFIT 6.1')
s14 = s24 = np.sqrt(0.10); s15 = s25 = np.sqrt(0.06)
common = dict(L=25.0*gd.UNIT_KM, L0=0.0, rho_central=3.e3, l_scale=10.0*gd.UNIT_KM,
              density_matter_is_in_g_per_cm3=True, nu_i=gd.NUE, nu_f=gd.NUE)
curves = [('2nu', op.osc_prob_2nu_matter_exp_density, 0.0005, 0.05, dict(sth=osc['s12'], Dm2=osc['D21'])),
          ('3nu', op.osc_prob_3nu_matter_exp_density, 0.002, 0.2, dict(osc)),
          ('4nu', op.osc_prob_4nu_matter_exp_density, 2.0, 20.0, dict(osc, s14=s14, s24=s24, D41=1.0)),
          ('5nu', op.osc_prob_5nu_matter_exp_density, 2.0, 20.0, dict(osc, s14=s14, s15=s15, s24=s24, s25=s25, D41=1.0, D51=1.7))]
for order in (8, None):
    for name, fn, lo, hi, ex in curves:
        E = np.logspace(np.log10(lo), np.log10(hi), 4)*gd.UNIT_GEV
        kw = dict(rtol=1e-12, atol=1e-14) if order is None else dict(rtol=1e-12, atol=1e-14, magnus_exp_order=order)
        info = {}
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            t = time.perf_counter(); fn(E, **common, **kw, **ex, strategy_info=info); dt = time.perf_counter() - t
        print('order %s %s: engine %s, declined %s, warnings %s, %.3f s' % (order, name, info.get('engine'), dict(info.get('declined', [])).get('hybrid'), sorted({x.category.__name__ for x in w}), dt))
# the worst point, against DOP853 on the captured Hamiltonian
cap = {}
real = ad.hybrid_propagator
def hp(H, *a, **k):
    cap['H'] = H; return real(H, *a, **k)
ad.hybrid_propagator = hp
op.osc_prob_3nu_matter_exp_density(np.array([0.002])*gd.UNIT_GEV, **common, **osc, rtol=1e-3, atol=1e-3, strategy='hybrid')
ad.hybrid_propagator = real
H = cap['H']
t = time.perf_counter()
sol = solve_ivp(lambda l, y: (-1j*np.asarray(H(l)) @ y.reshape(3, 3)).ravel(), (0.0, 25.0*gd.UNIT_KM), np.eye(3, dtype=complex).ravel(), rtol=1e-13, atol=1e-15, method='DOP853')
Pref = abs(sol.y[0, -1])**2; print('DOP853 %.2f s, P_ref %.15f' % (time.perf_counter() - t, Pref))
for order in (8, 4):
    info = {}
    P = float(np.asarray(op.osc_prob_3nu_matter_exp_density(np.array([0.002, 0.2])*gd.UNIT_GEV, **common, **osc, rtol=1e-12, atol=1e-14, magnus_exp_order=order, strategy_info=info)).ravel()[0])
    print('order %d at 2 MeV: engine %s, |P - ref| = %.2e (tol %.2e)' % (order, info.get('engine'), abs(P - Pref), 1e-14 + 1e-12*Pref))
