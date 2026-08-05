# -*- coding: utf-8 -*-
"""When the hybrid strategy declines, how good is what answers instead?

Found incidentally while auditing ``patch_atol``: on a multi-resonance **energy scan**, forcing
the hybrid path to decline left the energy-batched separable engine answering **2.08e-02** wrong
against a requested 1e-3.  That is a property of the fallback, not of the constant that exposed
it, and it had never been measured.

``strategy='auto'`` promises the best available answer, and the dispatch order is fixed:
hybrid -> interaction picture -> separable -> cumulative -> general ladder.  Each engine decides
only whether it *applies*, never whether it is the most accurate of those that do.  This asks
what that costs.

Method: for each workload, get the oracle, then force each engine that applies to answer alone
and score it.  Engines are forced with the same ``_ENGINES_DISABLED`` mechanism
``cross_check_strategies`` uses, so "which engine answered" is read rather than assumed.

Reports, per workload: every engine's error, which one ``'auto'`` actually picks, and the
penalty for that choice against the best engine that applied.

Run:  python fallback_quality.py [n_random]     # the synthetic population
      python fallback_quality.py --physical     # the physically-motivated one
"""

import sys
import warnings

import numpy as np

import harness as H
import magnus.oscprob as oscprob
from battery2 import bump_profile, ne_res_for
from battery3 import noisy_ne

TOL = 1.0e-3
L0, L1 = 0.0, 1.0*H.L_SCALE
SPAN = L1 - L0

# label -> (strategy, cumulative, engines to switch off so this one is reached)
FORCE = {
    'hybrid':     ('hybrid', False, ('ip_exp', 'separable')),
    'ip_exp':     ('magnus', False, ('hybrid', 'separable')),
    'separable':  ('magnus', False, ('hybrid', 'ip_exp')),
    'magnus':     ('magnus', False, ('hybrid', 'ip_exp', 'separable')),
    'cumulative': ('magnus', True, ('hybrid', 'ip_exp', 'separable')),
}


def run(label, ne, d, params, energies, baselines, forced=None, l0=L0):
    """One engine's answer, or None if it declined.  `forced=None` means an ordinary call."""
    same_e = bool(np.all(energies == energies[0]))
    E_arg = float(energies[0]) if same_e else energies
    L_arg = float(baselines[0]) if len(baselines) == 1 else baselines
    kw, disabled = {}, ()
    if forced is not None:
        strategy, cumulative, disabled = FORCE[forced]
        kw = dict(strategy=strategy, cumulative=cumulative)
    info = {}
    try:
        with oscprob._engine_probe(disabled), warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            P = np.asarray(oscprob.osc_prob_matter_std_potential(
                d, ne, E_arg, L_arg, params, L0=l0,
                density_is_of_number_of_electrons=True, strategy_info=info, **kw))
    except Exception as exc:                       # noqa: BLE001 -- a decline, not a failure
        return None, type(exc).__name__, []
    if (forced is not None) and info.get('engine') != forced:
        return None, 'declined (%s answered)' % info.get('engine'), []
    return (P.reshape(len(energies), d, d), info.get('engine'),
            sorted({w.category.__name__ for w in caught}))


def workloads(n_random):
    p2, p3 = H.params_for(2), H.params_for(3)
    ner2 = ne_res_for(2, p2, 10.0e6)
    profiles = [('solar', H.solar_ne()),
                ('multi-resonance', H.modulated_ne(amp=0.9, n_cycles=6.0, span=SPAN)),
                ('noisy', noisy_ne()),
                ('bump w=1e-2', bump_profile(ner2, 0.45*SPAN, 1e-2*SPAN))]
    rng = np.random.default_rng(90210)
    for i in range(n_random):
        profiles.append(('random Fourier #%d' % i,
                         H.fourier_ne(rng, n_modes=int(rng.integers(2, 9)), span=SPAN,
                                      base_ratio=float(10.0**rng.uniform(-2.5, -1.0)),
                                      amp=float(rng.uniform(0.2, 0.9)))))
    out = []
    for pname, ne in profiles:
        for d, params in ((2, p2), (3, p3)):
            out += shapes(pname, ne, d, params, L0, L1, 10.0e6, 100.0e6)
    return out


def shapes(pname, ne, d, params, l0, l1, e_lo, e_hi):
    """The three request shapes that route differently.  The E-scan is the one that exposed
    the fallback-quality defect, so a population without one is not evidence about it."""
    return [('%s d=%d point' % (pname, d), ne, d, params,
             np.array([e_lo]), np.array([l1]), l0),
            ('%s d=%d E-scan N=8' % (pname, d), ne, d, params,
             np.linspace(e_lo, e_hi, 8), np.full(8, l1), l0),
            ('%s d=%d L-scan N=8' % (pname, d), ne, d, params,
             np.full(8, e_lo), np.linspace(l0 + 0.2*(l1 - l0), l1, 8), l0)]


def physical_workloads():
    """The same three shapes, over the physically-motivated population.

    Each family carries its own trajectory and its own energy band -- a supernova ray is 1e4 to
    8e4 km and 15-45 MeV, an Earth chord is a few thousand km -- so the span cannot stay a
    module global here the way it can for the solar-only population.

    **Scope limit, stated rather than hidden.**  The E-scan needs eight *separate* oracle
    integrations (eight different energies at one baseline), and on a supernova ray at 15 MeV
    one of those costs ~170 s at d = 3 -- about 2.7 h for the SN families alone.  So the SN
    families here are run at the top of their energy band, where the oracle is cheapest, and
    without the E-scan.  Their accuracy at every energy in the band is covered instead by
    ``warn_fp.py --physical``, whose oracle reads one integration at eight stops.
    """
    import physical_profiles as pp
    out = []
    for f in pp.families():
        sn = f['kind'] in ('sn_shock', 'sn_turbulence')
        for d, params in ((2, H.params_for(2)), (3, H.params_for(3))):
            es = f['energies']
            rows = shapes(f['label'], f['ne'], d, params, f['l0'], f['l1'],
                          es[-1] if sn else es[0], es[-1])
            out += [r for r in rows if not (sn and 'E-scan' in r[0])]
    return out


def main():
    physical = '--physical' in sys.argv[1:]
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    n_random = int(args[0]) if args else 3
    wls = physical_workloads() if physical else workloads(n_random)
    print('# Fallback quality: every applicable engine scored on the same request')
    print('# population: %s' % ('PHYSICAL' if physical else 'synthetic'))
    print('# %d workloads, tolerance %.0e\n' % (len(wls), TOL))
    print('%-30s %-11s %s' % ('workload', "auto picks", ''.join('%12s' % e for e in FORCE)))
    rows = []
    for (label, ne, d, params, es, Ls, l0) in wls:
        try:
            if len(es) > 1 and bool(np.all(es == es[0])):
                # A baseline scan is ONE integration read at eight stops, not eight
                # integrations.  On the physical population that is the difference between
                # affordable and not: an SN ray at 15 MeV costs ~170 s per solve at d = 3.
                Hf = H.H_factory(d, params, H.vcc_of(ne), float(es[0]))
                Pref = np.array([H.P_of(U) for U in H.exact_U_many(Hf, l0, Ls, d)])
            else:
                Pref = np.array([H.P_of(H.exact_U(H.H_factory(d, params, H.vcc_of(ne),
                                                              float(e)), l0, float(L), d))
                                 for e, L in zip(es, Ls)])
        except Exception:                          # noqa: BLE001
            continue
        P_auto, eng_auto, warns_auto = run(label, ne, d, params, es, Ls, l0=l0)
        if P_auto is None:
            continue
        errs = {}
        for name in FORCE:
            P, who, _ = run(label, ne, d, params, es, Ls, forced=name, l0=l0)
            errs[name] = H.maxabs(P - Pref) if P is not None else None
        err_auto = H.maxabs(P_auto - Pref)
        available = {k: v for k, v in errs.items() if v is not None}
        best = min(available.values()) if available else err_auto
        rows.append(dict(label=label, engine=eng_auto, err_auto=err_auto, errs=errs,
                         best=best, penalty=err_auto/max(best, 1e-300),
                         silent=bool(err_auto > TOL and not warns_auto)))
        print('%-30s %-11s %s' % (label, eng_auto,
                                  ''.join('%12s' % ('%.2e' % errs[e] if errs[e] is not None
                                                    else '-') for e in FORCE)), flush=True)

    print('\n=== SUMMARY ===')
    print('workloads scored                      : %d' % len(rows))
    bad = [r for r in rows if r['penalty'] > 10.0]
    print("'auto' more than 10x worse than the best engine that applied: %d" % len(bad))
    for r in sorted(bad, key=lambda r: -r['penalty']):
        best_engine = min((k for k, v in r['errs'].items() if v is not None),
                          key=lambda k: r['errs'][k])
        print('   %-30s auto=%-11s %.2e   best=%-11s %.2e   %.0fx%s'
              % (r['label'], r['engine'], r['err_auto'], best_engine, r['best'],
                 r['penalty'], '   SILENT' if r['silent'] else ''))
    outside = [r for r in rows if r['err_auto'] > TOL]
    print("\n'auto' outside tolerance              : %d" % len(outside))
    print("  of those, silent                    : %d"
          % len([r for r in outside if r['silent']]))
    for name in FORCE:
        applied = [r for r in rows if r['errs'][name] is not None]
        if applied:
            worst = max(r['errs'][name] for r in applied)
            print('  %-11s applied to %3d workloads, worst error %.2e'
                  % (name, len(applied), worst))
    # One output file per run: two jobs writing the same path produced a stale read last
    # session that looked exactly like a real test failure.
    out = 'fallback_rows_physical.npy' if physical else 'fallback_rows.npy'
    np.save(out, np.array(rows, dtype=object), allow_pickle=True)
    print('\nrows saved to %s' % out)


if __name__ == '__main__':
    main()
