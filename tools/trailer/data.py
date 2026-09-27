r"""Computes, with Magnus, every number the trailer's new scenes show, and checks them.

Writes ``build/opening.npz``, ``build/adiabatic.npz`` and ``build/diagrams.json``, and prints
the checks the storyboard cites.  A few minutes on one core, most of it the adiabatic scene's
240 separate calls (one per point, so that each is answered by the engine ``strategy='auto'``
picks for it)::

    nice -n 19 python tools/trailer/data.py [opening] [code] [adiabatic] [diagrams]
"""
import collections
import inspect
import json
import sys
import time
import warnings

import numpy as np

from common import (BUILD, REPO, OPENING_E_GEV, OPENING_L_KM, opening_rho, ADIAB_E, ADIAB_TH,
                    ADIAB_DM2, ADIAB_S, ADIAB_W, DELTA, C2, S2, v_over_delta, adiab_rho_func,
                    reference_pee)
import magnus.oscprob as op
import magnus.globaldefs as gd
import magnus.matter as matter
import magnus.adiabatic as ad
from magnus.hamiltonians import hamiltonians3nu as h3, hamiltonians4nu as h4, hamiltonians5nu as h5

sys.path.insert(0, str(REPO / 'notebooks'))
P3 = gd.OSC_PARAMS_NU_FIT_6_1_SK_NO
ANGLES = ('s12', 's23', 's13', 'dCP', 'D21', 'D31')


def opening():
    """P(nu_mu -> nu_e) along 0-10,000 km at 3 GeV, in vacuum and in the varying-density track."""
    E, L = OPENING_E_GEV * gd.UNIT_GEV, OPENING_L_KM * gd.UNIT_KM
    Pv = op.osc_prob_3nu_vacuum(E, L, **{k: P3[k] for k in ANGLES}, nu_i=gd.NUMU, nu_f=gd.NUE)
    Pm = op.osc_prob_matter_std_potential(3, opening_rho, E, L, P3, nu_i=gd.NUMU, nu_f=gd.NUE,
                                          density_matter_is_in_g_per_cm3=True, rtol=1e-4, atol=1e-4)
    np.savez(BUILD / 'opening.npz', L=OPENING_L_KM, rho=opening_rho(L), Pv=np.ravel(Pv),
             Pm=np.ravel(Pm), E=OPENING_E_GEV)
    print('opening: max P vacuum %.3f, matter %.3f' % (np.max(Pv), np.max(Pm)))


def code():
    """The code moment's user Hamiltonian runs without a warning and matches the shipped wrapper."""
    H_vac = np.asarray(h3.hamiltonian_3nu_vacuum_energy_independent(*(P3[k] for k in ANGLES)))
    V = matter.vcc_func_from_rho_func(opening_rho, density_matter_is_in_g_per_cm3=True)
    P_e = matter.matter_potential_projector(3)

    def H(E, x):                                   # as on screen
        return H_vac / E + np.multiply.outer(V(x), P_e)

    energy = np.full(5, OPENING_E_GEV * gd.UNIT_GEV)
    L = np.linspace(2000, 10000, 5) * gd.UNIT_KM
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        P = op.osc_prob_energy_baseline(H, energy, L, nu_i=1, nu_f=0)
    Q = op.osc_prob_matter_std_potential(3, opening_rho, energy, L, P3, nu_i=1, nu_f=0,
                                         density_matter_is_in_g_per_cm3=True)
    names = sorted({w.category.__name__ for w in caught})
    print('code moment: max |P - wrapper| = %.1e, warnings: %s' % (np.max(np.abs(P - Q)), names or 'none'))


def adiabatic():
    """P_ee along the scene's path, one call per point under strategy='auto', checked against the
    independent reference; the instantaneous levels; the resonances; Magnus's patch window."""
    Lend = 3000.0 * ADIAB_S
    rho = adiab_rho_func()
    Ls = np.linspace(Lend / 240, Lend, 240)
    P, engines = [], []
    t0 = time.time()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        for L_ in Ls:
            info = {}
            P.append(float(np.ravel(op.osc_prob_matter_std_potential(
                2, rho, ADIAB_E, L_ / DELTA, dict(sth=np.sin(ADIAB_TH), Dm2=ADIAB_DM2), nu_i=0, nu_f=0,
                density_matter_is_in_g_per_cm3=True, strategy_info=info))[0]))
            engines.append(info.get('engine'))
    P = np.array(P)
    ref = np.array([reference_pee(L_) for L_ in Ls[23::24]])
    H0 = 0.5 * np.array([[-C2, S2], [S2, C2]])

    def H(x):
        return H0 + np.multiply.outer(v_over_delta(x), np.diag([1.0, 0.0]))

    wins = ad.find_nonadiabatic_windows(H, 0.0, Lend, n_probe=6400)[0]
    x = np.linspace(0, Lend, 20001)
    np.savez(BUILD / 'adiabatic.npz', x=x, V=v_over_delta(x), lam=np.linalg.eigvalsh(H(x)), Ls=Ls, P=P,
             engines=np.array(engines), res=x[np.where(np.diff(np.sign(v_over_delta(x) - C2)))[0]],
             wins=np.array(wins), c2=C2, s2=S2, km=1 / DELTA / gd.UNIT_KM, S=ADIAB_S, W=ADIAB_W)
    print('adiabatic: %.0f s, engines %s, warnings %s, max |P - reference| at %d points = %.1e, windows %s'
          % (time.time() - t0, dict(collections.Counter(engines)),
             sorted({w.category.__name__ for w in caught}) or 'none', len(ref),
             np.max(np.abs(P[23::24] - ref)), [tuple(np.round(w, 1)) for w in wins]))


def diagrams():
    """The slab ladder at one point, a 2000-energy Earth spectrum, and |U|^2 for 2-5 flavors."""
    import gen_profile_benchmarks as gpb
    import magnus.earth as earth
    out = {}
    kw = dict(nu_i=gd.NUMU, nu_f=gd.NUE, density_matter_is_in_g_per_cm3=True, strategy='magnus')
    E, L = OPENING_E_GEV * gd.UNIT_GEV, 1e4 * gd.UNIT_KM
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        out['ladder'] = [(n, float(np.ravel(op.osc_prob_matter_std_potential(
            3, opening_rho, E, L, P3, rtol=None, atol=None, n_slabs=n, **kw))[0]))
            for n in (2, 3, 4, 6, 9, 14, 21)]
        out['ladder_converged'] = float(np.ravel(op.osc_prob_matter_std_potential(
            3, opening_rho, E, L, P3, rtol=1e-9, atol=1e-9, **kw))[0])
    En = np.geomspace(1, 30, 2000) * gd.UNIT_GEV
    osc = gd.load_nufit_params('NuFIT 6.1')
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        Pe = np.ravel(op.osc_prob_3nu_earth(En, costhz=-0.8, L=earth.distance_traveled_inside_earth(-0.8) * gd.UNIT_KM,
                                            nu_i=gd.NUMU, nu_f=gd.NUE, **osc))
    out['fast'] = dict(E=(En / gd.UNIT_GEV).tolist(), P=Pe.tolist(),
                       warnings=sorted({w.category.__name__ for w in caught}))
    U = {}
    th = np.arcsin(np.sqrt(0.307))
    U[2] = np.array([[np.cos(th)**2, np.sin(th)**2], [np.sin(th)**2, np.cos(th)**2]])
    p3 = gpb.osc_params(3)
    U[3] = np.abs(np.asarray(h3.mixing_matrix_3x3(p3['s12'], p3['s23'], p3['s13'], p3['dCP'])))**2
    for d, fn in ((4, h4.mixing_matrix_4x4), (5, h5.mixing_matrix_5x5)):
        pp = gpb.osc_params(d)
        names = [('dCP' if a == 'd13' else a) for a in inspect.signature(fn).parameters]
        U[d] = np.abs(np.asarray(fn(*[pp[a] for a in names if a in pp])))**2
    out['U'] = {d: U[d].tolist() for d in U}
    (BUILD / 'diagrams.json').write_text(json.dumps(out))
    print('diagrams: ladder %s, converged %.5f; Earth spectrum warnings %s; |U|^2 rows sum to 1: %s'
          % ([(n, round(p, 5)) for n, p in out['ladder']], out['ladder_converged'], out['fast']['warnings'] or 'none',
             all(np.allclose(np.sum(U[d], axis=1), 1) for d in U)))


STEPS = {'opening': opening, 'code': code, 'adiabatic': adiabatic, 'diagrams': diagrams}

if __name__ == '__main__':
    BUILD.mkdir(exist_ok=True)
    for name in sys.argv[1:] or list(STEPS):
        STEPS[name]()
