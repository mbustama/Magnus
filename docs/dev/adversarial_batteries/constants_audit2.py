# -*- coding: utf-8 -*-
"""Provenance for the seven remaining unaudited constants (robustness programme, item 3).

``constants_audit.py`` covered the two the brief ranked first.  These are the rest:

    adiabatic:  min_threshold (1e-6)  n_probe0 (200)  n_points0 (201)
                patch_atol (1e-7)     n_slabs0 (400)
    oscprob:    growth_factor_n_slabs (1.5)   min_n_tpts_per_slab (2)

**The population spans single points, baseline scans AND energy scans**, deliberately.  The
``threshold0`` audit swept a fixed baseline, concluded a tolerance-derived rule was safe, and was
refuted by an energy scan that got 20x worse -- see ``adiabatic.THRESHOLD0_PROVENANCE``.  A
population that does not contain the workload you are about to change is not evidence about it.

Each reference solution is computed once per workload and reused across every value of every
constant, so the ``solve_ivp`` cost is paid once rather than once per sweep point.

Run:  python constants_audit2.py [which ...]     (which in {probe, points, minthr, patch,
                                                  slabs0, growth, mintpts})
"""

import sys
import time
import warnings

import numpy as np

import harness as H
import magnus.adiabatic as ad
import magnus.oscprob as oscprob
from battery2 import bump_profile, ne_res_for

TOL = 1.0e-3
L0, L1 = 0.0, 1.0*H.L_SCALE
SPAN = L1 - L0


def workloads():
    """(label, d, energies, baselines) spanning point / baseline scan / energy scan."""
    p2, p3 = H.params_for(2), H.params_for(3)
    ner2 = ne_res_for(2, p2, 10.0e6)
    profiles = [
        ('solar', H.solar_ne()),
        ('multi-resonance', H.modulated_ne(amp=0.9, n_cycles=6.0, span=SPAN)),
        ('bump w=1e-2', bump_profile(ner2, 0.45*SPAN, 1e-2*SPAN)),
    ]
    out = []
    for pname, ne in profiles:
        for d, params in ((2, p2), (3, p3)):
            out.append(('%s d=%d point' % (pname, d), ne, d, params,
                        np.array([10.0e6]), np.array([L1])))
            out.append(('%s d=%d L-scan N=8' % (pname, d), ne, d, params,
                        np.full(8, 10.0e6), np.linspace(0.2*L1, L1, 8)))
            out.append(('%s d=%d E-scan N=8' % (pname, d), ne, d, params,
                        np.linspace(10.0e6, 100.0e6, 8), np.full(8, L1)))
    return out


def reference(ne, d, params, energies, baselines):
    """solve_ivp per (energy, baseline) point.  Computed once; reused for every sweep value."""
    out = []
    for e, L in zip(energies, baselines):
        H_of_l = H.H_factory(d, params, H.vcc_of(ne), float(e))
        out.append(H.P_of(H.exact_U(H_of_l, L0, float(L), d)))
    return np.array(out)


def score(ne, d, params, energies, baselines, Pref, **call_kwargs):
    e0 = energies[0]
    same_e = bool(np.all(energies == e0))
    E_arg = float(e0) if same_e else energies
    L_arg = float(baselines[0]) if len(baselines) == 1 else baselines
    t0 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P = np.asarray(oscprob.osc_prob_matter_std_potential(
            d, ne, E_arg, L_arg, params, L0=L0,
            density_is_of_number_of_electrons=True, **call_kwargs))
    dt = time.time() - t0
    return H.maxabs(P.reshape(Pref.shape) - Pref), dt


def sweep(name, values, apply_value, wls, refs):
    """apply_value(v) -> (call_kwargs, teardown) ; teardown() restores any monkeypatch."""
    print('\n## %s   default marked *' % name)
    header = ''.join('%13s' % ('%g' % v) for v in values)
    print('%-28s %s' % ('workload', header))
    worst = {v: 0.0 for v in values}
    for (label, ne, d, params, es, Ls), Pref in zip(wls, refs):
        row, times = [], []
        for v in values:
            kwargs, teardown = apply_value(v)
            try:
                err, dt = score(ne, d, params, es, Ls, Pref, **kwargs)
            finally:
                teardown()
            row.append(err)
            times.append(dt)
            worst[v] = max(worst[v], err)
        print('%-28s %s' % (label, ''.join('%13.2e' % e for e in row)))
    print('%-28s %s' % ('WORST over workloads', ''.join('%13.2e' % worst[v] for v in values)))
    return worst


def main():
    which = sys.argv[1:] or ['probe', 'points', 'minthr', 'patch', 'slabs0', 'growth',
                             'mintpts']
    wls = workloads()
    print('# Seven remaining constants.  %d workloads (point / L-scan / E-scan x 3 profiles '
          'x d=2,3), tolerance %.0e' % (len(wls), TOL))
    print('# references: solve_ivp/DOP853 rtol=1e-12, computed once each', flush=True)
    refs = []
    for i, (label, ne, d, params, es, Ls) in enumerate(wls):
        refs.append(reference(ne, d, params, es, Ls))
        print('   reference %d/%d done (%s)' % (i + 1, len(wls), label), flush=True)

    def noop():
        return None

    def hybrid_sweep(key):
        def apply(v):
            orig = ad.hybrid_propagator

            def patched(*a, **k):
                k[key] = v
                return orig(*a, **k)
            ad.hybrid_propagator = patched
            oscprob.adiabatic.hybrid_propagator = patched

            def teardown():
                ad.hybrid_propagator = orig
                oscprob.adiabatic.hybrid_propagator = orig
            return {}, teardown
        return apply

    def patch_sweep(key):
        """n_slabs0 / patch_atol live in _local_evolution_operator, unreachable from above."""
        def apply(v):
            orig = ad._local_evolution_operator

            def patched(*a, **k):
                k[key] = v
                return orig(*a, **k)
            ad._local_evolution_operator = patched

            def teardown():
                ad._local_evolution_operator = orig
            return {}, teardown
        return apply

    def plain(key):
        return lambda v: ({key: v}, noop)

    jobs = {
        'probe':   ('adiabatic n_probe0 (default 200*)', [50, 100, 200, 400, 800],
                    hybrid_sweep('n_probe0')),
        'points':  ('adiabatic n_points0 (default 201*)', [51, 101, 201, 401, 801],
                    hybrid_sweep('n_points0')),
        'minthr':  ('adiabatic min_threshold (default 1e-6*)', [1e-4, 1e-5, 1e-6, 1e-7, 1e-8],
                    hybrid_sweep('min_threshold')),
        'patch':   ('adiabatic patch_atol (default 1e-7*)', [1e-5, 1e-6, 1e-7, 1e-8, 1e-9],
                    patch_sweep('patch_atol')),
        'slabs0':  ('adiabatic n_slabs0 (default 400*)', [100, 200, 400, 800, 1600],
                    patch_sweep('n_slabs0')),
        'growth':  ('oscprob growth_factor_n_slabs (default 1.5*)', [1.2, 1.5, 2.0, 3.0],
                    plain('growth_factor_n_slabs')),
        'mintpts': ('oscprob min_n_tpts_per_slab (default 2*)', [2, 4, 8],
                    plain('min_n_tpts_per_slab')),
    }
    for key in which:
        name, values, apply_value = jobs[key]
        sweep(name, values, apply_value, wls, refs)


if __name__ == '__main__':
    main()
