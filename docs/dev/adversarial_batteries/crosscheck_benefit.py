# -*- coding: utf-8 -*-
"""Does the weak-band cross-check earn its cost? (robustness programme, item 6)

``CROSS_CHECK_MARGIN`` decides when ``strategy='auto'`` verifies a window-free hybrid result
against the general Magnus ladder.  The cost is measured (about 9 % of calls, 2x each).  This
measures the **benefit**, which is the harder half and the one it would be easy to assume:

  * how often do the two engines actually disagree beyond the requested tolerance?
  * when they disagree, **which one is right**?  Scored against ``solve_ivp``, because a
    disagreement alone does not say who is wrong, and a check that discarded the better answer
    would be worse than no check.

Run with the margin forced to 0.0 so that *every* window-free result is verified, which is the
population the constant selects a subset of.

Run:  python crosscheck_benefit.py [n_cases]
"""

import sys
import warnings

import numpy as np

import harness as H
import magnus.adiabatic as ad
import magnus.oscprob as oscprob

TOL = 1.0e-3
BOUND = 2*TOL


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    rng = np.random.default_rng(20260805)
    seen = agreed = disagreed = 0
    hybrid_better = general_better = 0
    rows = []

    for _ in range(n):
        d = int(rng.choice([2, 3]))
        energy = float(10.0**rng.uniform(7.3, 8.3))
        span = float(rng.uniform(0.3, 1.0))*H.L_SCALE
        ne = H.fourier_ne(rng, n_modes=int(rng.integers(2, 9)), span=span,
                          base_ratio=float(10.0**rng.uniform(-2.5, -1.0)),
                          amp=float(rng.uniform(0.2, 0.95)))
        L = float(rng.uniform(0.2, 1.0))*span
        params = H.params_for(d)
        try:
            H_of_l = H.H_factory(d, params, H.vcc_of(ne), energy)
            info = {}
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                U, _, certified = ad.hybrid_propagator(H_of_l, 0.0, L, rtol=TOL, atol=TOL,
                                                       info=info)
        except Exception:                          # noqa: BLE001
            continue
        if not certified or info.get('n_windows', 0) != 0:
            continue
        margin = ad.GAMMA_TO_ERROR*info.get('gamma_max', 0.0)/BOUND
        seen += 1

        ok, spread = oscprob._verify_against_general_path(
            H_of_l, 0.0, L, TOL, TOL, 4, 'gl', U)
        if ok:
            agreed += 1
            rows.append((margin, spread, None, None))
            continue

        disagreed += 1
        # Who was right?  Only worth an oracle when they disagree.
        try:
            Pref = H.P_of(H.exact_U(H_of_l, 0.0, L, d))
        except Exception:                          # noqa: BLE001
            continue
        P_hyb = np.swapaxes(U.real**2 + U.imag**2, -1, -2)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            P_gen = np.asarray(oscprob.osc_prob(H_of_l, 0.0, L, rtol=TOL, atol=TOL))
        e_h, e_g = H.maxabs(P_hyb - Pref), H.maxabs(P_gen - Pref)
        if e_g < e_h:
            general_better += 1
        else:
            hybrid_better += 1
        rows.append((margin, spread, e_h, e_g))

    print('window-free certified results examined : %d' % seen)
    print('  agreed with the general ladder       : %d (%.0f%%)'
          % (agreed, 100.0*agreed/max(seen, 1)))
    print('  DISAGREED beyond the tolerance       : %d (%.0f%%)'
          % (disagreed, 100.0*disagreed/max(seen, 1)))
    if disagreed:
        print('\n  of the disagreements, scored against solve_ivp:')
        print('    the general ladder was closer      : %d   <-- the check helps'
              % general_better)
        print('    the adiabatic answer was closer    : %d   <-- the check would HURT'
              % hybrid_better)
        print('\n  %-10s %-12s %-12s %-12s' % ('margin', 'spread', 'err hybrid', 'err general'))
        for m, s, eh, eg in sorted([r for r in rows if r[2] is not None],
                                   key=lambda r: -r[1]):
            print('  %-10.4f %-12.3e %-12.3e %-12.3e%s'
                  % (m, s, eh, eg, '   general better' if eg < eh else '   hybrid better'))

    covered = [r for r in rows if r[0] > oscprob.CROSS_CHECK_MARGIN]
    caught = [r for r in covered if r[2] is not None]
    print('\n  at the shipped CROSS_CHECK_MARGIN = %.2f:' % oscprob.CROSS_CHECK_MARGIN)
    print('    window-free results verified       : %d of %d (%.0f%%)'
          % (len(covered), seen, 100.0*len(covered)/max(seen, 1)))
    print('    disagreements caught               : %d of %d' % (len(caught), disagreed))


if __name__ == '__main__':
    main()
