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


def paper():
    """The remade paper scenes.  What the paper already computed is read from
    notebooks/paper_figure_cache.json (read only: this never writes it); the rest is computed
    here with the paper's own settings (notebook 28's cells, cited per block)."""
    import magnus.earth as earth
    cache = json.loads((REPO / 'notebooks' / 'paper_figure_cache.json').read_text())
    val = lambda k: cache[k]['value']                                  # noqa: E731
    out = {}
    OSC = gd.load_nufit_params('NuFIT 6.1')
    # cell 42: the oscillograms, 170 directions x 200 energies each
    out['osc_cz'] = np.linspace(-1.0, -0.05, 170)
    out['osc_e_gev'] = np.logspace(np.log10(2.0), np.log10(60.0), 200)
    out['osc_e_tev'] = np.logspace(np.log10(1.0), np.log10(30.0), 200)
    for k, tag in (('osc_3nu', 'Earth_3_nu'), ('osc_3p1', 'Earth_3_1'), ('osc_3p2', 'Earth_3_2')):
        out[k] = np.asarray(val('oscillogram_' + tag))
    out['osc_chord_km'] = np.array([earth.distance_traveled_inside_earth(float(c)) for c in out['osc_cz']])
    # cell 62: the long-range force in the Sun
    lr = val('solar_long_range')
    out['lr_e_mev'] = np.logspace(np.log10(0.1), np.log10(20.0), 70)
    out['lr_std'], out['lr_1'], out['lr_01'] = (np.asarray(lr[k]) for k in ('std', '1', '0.1'))
    # cell 66: the Sun in the nu_e channel, P against impact parameter at nine energies
    energies = (0.01, 10.0, 30.0, 100.0, 300.0, 1.0e3, 3.0e3, 1.0e4, 5.0e4)
    out['sun_e_gev'] = np.array(energies)
    discs = []
    for e in energies:
        P = np.asarray(val('solar_disk_%s' % ('%g' % e).replace('.', 'p').replace('+', '')), float)
        grid = np.linspace(0.0, 1.0, len(P))
        P[np.isnan(P)] = P[np.isfinite(P)][-1]
        discs.append(np.interp(np.linspace(0, 1, 2001), grid, P))
    out['sun_P'] = np.array(discs)
    # cell 69: the jet, P(nu_e -> nu_e) at Earth
    out['jet_e_tev'] = np.geomspace(0.1, 1.0e4, 161)
    for k in ('smooth', 'turbulent', 'stepped'):
        out['jet_' + k] = np.asarray(val('jet_' + k)['P_earth'])[:, gd.NUE, gd.NUE]
    # cells 86-88: geoneutrinos, four production points
    out['geo_e_mev'] = 1.0 / np.linspace(1.0 / 1.8, 1.0 / 3.3, 22500)
    for k, depth, cz, from_detector in (('local', 10.0, 0.078, True), ('far_crust', 20.0, -0.272, False),
                                        ('mantle', 1000.0, -0.552, False), ('core', 2800.0, -0.872, False)):
        out['geo_' + k] = np.asarray(val('geo_energy_' + k))
        depths = (dict(source_depth=1.4, detector_depth=depth) if from_detector
                  else dict(source_depth=depth, detector_depth=1.4))
        out['geo_L_' + k] = earth.distance_traveled_inside_earth(cz, **depths)          # km
    out['geo_vac_avg'] = float(op.osc_prob_3nu_vacuum(2.5 * gd.UNIT_MEV, 1.0e8 * gd.UNIT_KM, average=True,
                                                      nubar=True, nu_i=gd.NUE, nu_f=gd.NUE, **OSC))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # cell 44: from Fermilab to four sites
        E = np.logspace(np.log10(0.3), np.log10(10.0), 400)
        out['fnal_e_gev'] = E
        for site in ('snolab', 'homestake', 'cern', 'south_pole'):
            a, b = earth.loc_coords_dms['fermilab'], earth.loc_coords_dms[site]
            out['fnal_L_' + site] = earth.chord_length_inside_earth(a['lat'], a['lon'], b['lat'], b['lon'])
            out['fnal_' + site] = np.asarray(op.osc_prob_3nu_earth(E * gd.UNIT_GEV, loc_ini='fermilab', loc_fin=site,
                                                                   nu_i=gd.NUMU, nu_f=gd.NUE, **OSC))
        # CP violation: the bi-probability ellipses at Fermilab -> Homestake, 2.5 GeV, NO and IO
        dcp = np.linspace(0, 2 * np.pi, 121)
        out['cp_dcp'] = dcp
        for order in ('NO', 'IO'):
            p = gd.load_nufit_params('NuFIT 6.1', ordering=order)
            nu, nb = [], []
            for d in dcp:
                q = dict(p, dCP=float(d))
                kw = dict(loc_ini='fermilab', loc_fin='homestake', nu_i=gd.NUMU, nu_f=gd.NUE, **q)
                nu.append(float(np.ravel(op.osc_prob_3nu_earth(2.5 * gd.UNIT_GEV, **kw))[0]))
                nb.append(float(np.ravel(op.osc_prob_3nu_earth(2.5 * gd.UNIT_GEV, nubar=True, **kw))[0]))
            out['cp_%s_nu' % order], out['cp_%s_nubar' % order] = np.array(nu), np.array(nb)
        # cell 27: standard, NSI and LIV through the Earth at cos(theta_z) = -0.9
        cz = -0.9
        L = earth.distance_traveled_inside_earth(cz) * gd.UNIT_KM
        E = np.logspace(0, np.log10(40.0), 260) * gd.UNIT_GEV
        EPS = dict(eps_ee=0.10, eps_em=0.05 + 0.0j, eps_et=0.0j, eps_mm=0.0, eps_mt=0.03 + 0.0j, eps_tt=0.0)
        LIV = dict(b1=0.0, b2=0.0, b3=np.pi / L, Lambda=1.0, sxi12=0.0, sxi23=1.0 / np.sqrt(2.0), sxi13=0.0,
                   dxiCP=0.0, n_liv=0)
        kw = dict(costhz=cz, L=L, nu_i=gd.NUMU, nu_f=gd.NUMU, rtol=1e-6, atol=1e-6, **OSC)
        out['bsm_e_gev'] = E / gd.UNIT_GEV
        out['bsm_std'] = np.asarray(op.osc_prob_3nu_earth(E, **kw))
        out['bsm_nsi'] = np.asarray(op.osc_prob_3nu_earth_nsi(E, **kw, **EPS))
        out['bsm_liv'] = np.asarray(op.osc_prob_3nu_earth_liv(E, **kw, **LIV))
        # cell 84: flavor ratios at Earth from a pion-decay source, 100 TeV, 100 Mpc
        E_tri, L_src = 100.0 * gd.UNIT_TEV, 100.0 * 3.0857e19 * gd.CONV_KM_TO_INV_EV
        src = np.array([1 / 3, 2 / 3, 0.0])
        liv_a = dict(sxi12=OSC['s12'], sxi23=OSC['s23'], sxi13=OSC['s13'], dxiCP=0.0, b1=0.0, b2=0.0,
                     Lambda=1.0, n_liv=1)
        b3u = float(OSC['D31']) / (2.0 * E_tri**2)

        def frac(P, n=3):
            f = np.einsum('a,ab->b', np.pad(src, (0, P.shape[0] - 3)), P)[:3]
            return f / f.sum()

        ratios = np.concatenate([[0.0], np.logspace(-3.0, 3.0, 121)])
        out['tri_liv'] = np.array([frac(np.asarray(op.osc_prob_3nu_vacuum_liv(
            E_tri, L_src, average=True, b3=r * b3u, **liv_a, **OSC))) for r in ratios])
        s2 = np.linspace(0.0, 0.3, 61)
        out['tri_sterile'] = np.array([frac(np.asarray(op.osc_prob_4nu_vacuum(
            E_tri, L_src, average=True, s14=np.sqrt(x), s24=np.sqrt(x), s34=0.0, D41=1.0, **OSC))) for x in s2])
        out['tri_vac'] = out['tri_liv'][0]
        # cell 48: a buried body swept by a beam over 1500 km of crust
        from magnus import matter as mt
        L0, Es, alpha = 1500.0, np.linspace(25.0, 150.0, 400) * gd.UNIT_MEV, np.linspace(-15.0, 15.0, 160)
        ne = lambda rho: mt.num_density_e_func(0.0, lambda _: rho, electron_fraction=0.5,     # noqa: E731
                                               ratio_number_neutrons_to_protons=1.0,
                                               density_matter_is_in_g_per_cm3=True)
        ne_c, ne_b = ne(3.3), ne(10.0)
        ckw = dict(osc_params=OSC, L0=0.0, nu_i=gd.NUE, nu_f=gd.NUE, nubar=True, rtol=1e-6, atol=1e-6,
                   density_is_of_number_of_electrons=True)
        P0 = np.asarray(op.osc_prob_matter_std_potential(3, ne_c, Es, L0 * gd.UNIT_KM, **ckw))
        dP = np.zeros((len(Es), len(alpha)))
        R, D0 = 125.0, 750.0
        for j, a in enumerate(alpha):
            miss = abs(D0 * np.sin(np.radians(a)))
            if miss >= R:
                continue
            h, mid = np.sqrt(R**2 - miss**2), D0 * np.cos(np.radians(a))
            lo, hi = mid - h, mid + h
            prof = lambda l, lo=lo, hi=hi: ne_c + (ne_b - ne_c) * ((np.asarray(l) / gd.UNIT_KM >= lo) &  # noqa: E731
                                                                   (np.asarray(l) / gd.UNIT_KM <= hi))
            dP[:, j] = np.asarray(op.osc_prob_matter_std_potential(3, prof, Es, L0 * gd.UNIT_KM,
                                  t_breakpoints=np.array([lo, hi]) * gd.UNIT_KM, **ckw)) - P0
        out['cav_e_mev'], out['cav_alpha'], out['cav_dP'] = Es / gd.UNIT_MEV, alpha, dP
    np.savez(BUILD / 'paper.npz', **out)
    # the continents of the paper's globes: the LAND polygons of notebook 28 (Figure 3e), [lon, lat]
    import ast
    nb = json.loads((REPO / 'notebooks' / '28_magnus_paper_figures.ipynb').read_text())
    src = next(''.join(c['source']) for c in nb['cells'] if c['cell_type'] == 'code' and 'LAND = [' in ''.join(c['source']))
    line = next(ln for ln in src.splitlines() if ln.startswith('LAND = ['))
    (BUILD / 'land.json').write_text(json.dumps(ast.literal_eval(line[len('LAND = '):])))
    print('paper: oscillograms %s; Fermilab chords %s km; cavity |dP| max %.3f; flavor triangle %d + %d points'
          % (out['osc_3nu'].shape, [int(out['fnal_L_' + s]) for s in ('snolab', 'homestake', 'cern', 'south_pole')],
             np.abs(out['cav_dP']).max(), len(out['tri_liv']), len(out['tri_sterile'])))


STEPS = {'opening': opening, 'code': code, 'adiabatic': adiabatic, 'diagrams': diagrams, 'paper': paper}

if __name__ == '__main__':
    BUILD.mkdir(exist_ok=True)
    for name in sys.argv[1:] or list(STEPS):
        STEPS[name]()
