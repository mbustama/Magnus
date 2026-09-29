"""Listing 1 at the default tolerance under the #70 routing: engine, time, error vs the tight run."""
import sys, time, json
import numpy as np
import magnus.globaldefs as gd
import magnus.oscprob as oscprob

osc = gd.load_nufit_params('NuFIT 6.1')
s14 = s24 = np.sqrt(0.10); s15 = s25 = np.sqrt(0.06)
def energies(lo, hi):
    return np.logspace(np.log10(lo), np.log10(hi), 140)*gd.UNIT_GEV
base = dict(L=25.0*gd.UNIT_KM, L0=0.0, rho_central=3.e3, l_scale=10.0*gd.UNIT_KM,
            density_matter_is_in_g_per_cm3=True, nu_i=gd.NUE, nu_f=gd.NUE)
curves = {
    '2nu': (oscprob.osc_prob_2nu_matter_exp_density, energies(0.0005, 0.05), dict(sth=osc['s12'], Dm2=osc['D21'])),
    '3nu': (oscprob.osc_prob_3nu_matter_exp_density, energies(0.002, 0.2), dict(**osc)),
    '4nu': (oscprob.osc_prob_4nu_matter_exp_density, energies(2.0, 20.0), dict(**osc, s14=s14, s24=s24, D41=1.0)),
    '5nu': (oscprob.osc_prob_5nu_matter_exp_density, energies(2.0, 20.0),
            dict(**osc, s14=s14, s15=s15, s24=s24, s25=s25, D41=1.0, D51=1.7)),
}
mode = sys.argv[1]   # 'ref' or 'run'
tol = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-3
if mode == 'ref':
    out = {}
    for k, (f, E, p) in curves.items():
        out[k] = np.asarray(f(E, **base, **p, rtol=1e-12, atol=1e-14, strategy='hybrid'))
    np.savez('l1_ref.npz', **out); print('ref saved')
else:
    ref = np.load('l1_ref.npz')
    for k, (f, E, p) in curves.items():   # warm the kernels on a 2-point call
        f(E[:2], **base, **p, rtol=tol, atol=tol)
    tot = 0.0
    for k, (f, E, p) in curves.items():
        info = {}
        t = time.perf_counter()
        P = np.asarray(f(E, **base, **p, rtol=tol, atol=tol, strategy_info=info))
        dt = time.perf_counter() - t; tot += dt
        err = float(np.max(np.abs(P - ref[k])))
        print('%s t=%.3f s err=%.2e route=%s' % (k, dt, err, info.get('route', info.get('engine'))))
    print('total=%.3f' % tot)
