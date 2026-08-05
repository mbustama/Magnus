# -*- coding: utf-8 -*-
"""Oracle-free invariants, swept over a profile matrix (robustness programme, item 4).

These encode *"the answer must not depend on which door you came in"*.  They need no
``solve_ivp``, so they are cheap enough to run in CI -- and the dispatch seams they cover are
exactly where every defect in ``FINDINGS_ADVERSARIAL_VALIDATION.md`` lived.

This script **measures** the invariants; ``tests/test_invariants.py`` asserts on them.  The two
are separate on purpose: the brief's own warning is that some of these disagreements are
*correct* -- ``auto`` and ``magnus`` are different methods and ``auto`` is usually the better
one -- so the tolerance per invariant has to come from a measured distribution rather than from
a guess.  Run this, read the table, then set the test's bounds from it.

Invariants, and what each would catch:

  I1  strategy 'auto' / 'magnus' / 'hybrid' agree            -- a dispatch seam moving answers
  I2  cumulative True / False / 'auto' agree                 -- the second dispatch seam
  I3  a scan agrees with the same points computed singly     -- the batching seams
  I4  shuffled baselines, unshuffled after                   -- ordering assumptions (exact)
  I5  n_jobs=2 agrees with n_jobs=1                          -- parallel warm-start state (exact)
  I6  U(0->L2) == U(L1->L2) U(0->L1)                         -- time ordering (operator level)
  I7  probability rows and columns sum to 1                  -- unitarity, on the output itself

Run:  python invariants.py
"""

import sys
import time

import numpy as np

import harness as H
import magnus.adiabatic as ad
import magnus.oscprob as oscprob
from battery2 import bump_profile, ne_res_for
from battery3 import noisy_ne

L0, L1 = 0.0, 1.0*H.L_SCALE
SPAN = L1 - L0
QUIET = dict(density_is_of_number_of_electrons=True)


def profiles():
    p2 = H.params_for(2)
    ner2 = ne_res_for(2, p2, 10.0e6)
    return [
        ('solar exponential', H.solar_ne()),
        ('multi-resonance', H.modulated_ne(amp=0.9, n_cycles=6.0, span=SPAN)),
        ('noisy', noisy_ne()),
        ('sinusoid span/7', H.sine_ne(SPAN/7.0, base_ratio=3.0e-2)),
        ('resolvable bump w=1e-2', bump_profile(ner2, 0.45*SPAN, 1e-2*SPAN)),
    ]


def call(ne, d, params, energy, Ls, **kw):
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P = oscprob.osc_prob_matter_std_potential(
            d, ne, energy, Ls if np.size(Ls) > 1 else float(np.ravel(Ls)[0]), params,
            L0=L0, **QUIET, **kw)
    return np.asarray(P).reshape(np.size(Ls), d, d)


def main():
    rows = []
    t0 = time.time()
    for label, ne in profiles():
        for d in (2, 3):
            params = H.params_for(d)
            for energy in (10.0e6, 50.0e6):
                for N in (1, 8, 30):
                    Ls = np.linspace(0.05*L1, L1, N) if N > 1 else np.array([L1])
                    tag = '%s d=%d E=%.0fMeV N=%d' % (label, d, energy/1e6, N)
                    r = dict(tag=tag, label=label, d=d, energy=energy, N=N)

                    P_auto = call(ne, d, params, energy, Ls)
                    P_mag = call(ne, d, params, energy, Ls, strategy='magnus')
                    P_hyb = call(ne, d, params, energy, Ls, strategy='hybrid')
                    r['I1 auto-magnus'] = H.maxabs(P_auto - P_mag)
                    r['I1 auto-hybrid'] = H.maxabs(P_auto - P_hyb)

                    if N > 1:
                        P_cT = call(ne, d, params, energy, Ls, cumulative=True)
                        P_cF = call(ne, d, params, energy, Ls, cumulative=False)
                        r['I2 cumT-cumF'] = H.maxabs(P_cT - P_cF)
                        r['I2 auto-cumT'] = H.maxabs(P_auto - P_cT)

                        singly = np.array([call(ne, d, params, energy, np.array([L]))[0]
                                           for L in Ls])
                        r['I3 scan-singly'] = H.maxabs(P_auto - singly)

                        perm = np.random.default_rng(4).permutation(N)
                        P_sh = call(ne, d, params, energy, Ls[perm])
                        back = np.empty_like(P_sh)
                        back[perm] = P_sh
                        r['I4 shuffled'] = H.maxabs(P_auto - back)

                        P_par = call(ne, d, params, energy, Ls, n_jobs=2)
                        r['I5 n_jobs'] = H.maxabs(P_auto - P_par)

                    # I6: composition, at the operator level (probabilities do not compose)
                    H_of_l = H.H_factory(d, params, H.vcc_of(ne), energy)
                    lm = 0.5*(L0 + L1)
                    U_full, _, _ = ad.hybrid_propagator(H_of_l, L0, L1)
                    U_a, _, _ = ad.hybrid_propagator(H_of_l, L0, lm)
                    U_b, _, _ = ad.hybrid_propagator(H_of_l, lm, L1)
                    r['I6 compose hybrid'] = H.maxabs(U_full - U_b @ U_a)
                    # A CONVERGED grid on both sides, not one slab.  The first version of this
                    # row used a single slab over the whole path and reported a "disagreement"
                    # of 1.5 -- which was the harness comparing two unconverged integrations,
                    # not the package failing to compose.  The instrument was wrong before the
                    # code was, which is the pattern this whole exercise keeps reproducing.
                    def V(a, b, n=4000):
                        edges = np.linspace(a, b, n + 1)
                        chain = oscprob.compute_evolution_operator_multiple_slabs(
                            H_of_l, np.column_stack([edges[:-1], edges[1:]]), 2, 6)
                        out = chain[0]
                        for U in chain[1:]:
                            out = U @ out
                        return out
                    V_full, V_a, V_b = V(L0, L1), V(L0, lm), V(lm, L1)
                    r['I6 compose magnus'] = H.maxabs(V_full - V_b @ V_a)

                    # I7: unitarity, read off the probabilities the user receives
                    r['I7 unitarity'] = max(H.maxabs(P_auto.sum(axis=1) - 1.0),
                                            H.maxabs(P_auto.sum(axis=2) - 1.0))
                    rows.append(r)
                    print('  %-46s done (%.0f s)' % (tag, time.time() - t0), flush=True)

    keys = [k for k in rows[0] if k[0] == 'I'] + ['I2 cumT-cumF', 'I2 auto-cumT',
                                                  'I3 scan-singly', 'I4 shuffled', 'I5 n_jobs']
    keys = sorted({k for r in rows for k in r if k.startswith('I')})
    print('\n=== INVARIANT SUMMARY over %d configurations ===' % len(rows))
    print('%-22s %11s %11s %11s   %s' % ('invariant', 'median', 'p90', 'max', 'worst case'))
    for k in keys:
        vals = np.array([r[k] for r in rows if k in r])
        worst = max((r for r in rows if k in r), key=lambda r: r[k])
        print('%-22s %11.3e %11.3e %11.3e   %s'
              % (k, np.median(vals), np.percentile(vals, 90), vals.max(), worst['tag']))
    np.save('invariant_rows.npy', np.array(rows, dtype=object), allow_pickle=True)
    print('\nrows saved to invariant_rows.npy')


if __name__ == '__main__':
    sys.exit(main())
