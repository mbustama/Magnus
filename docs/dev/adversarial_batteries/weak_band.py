# -*- coding: utf-8 -*-
"""Where is the hybrid path's self-certification actually weak? (robustness programme, item 6)

Every remaining silent miss lives below the N = 25 seam, on the hybrid path.  The proposal is to
have ``strategy='auto'`` verify itself there against an independent engine -- but doing that on
every call would double the cost of the package's most common request, so the trigger has to be
targeted.  This measures what to target.

For each case it records what the hybrid propagator knew about itself (did a window open, how
big was gamma against the certification bound) alongside the error it actually made, so the
trigger condition can be read off the data instead of guessed.

The hypothesis to test, from ``FINDINGS_ADVERSARIAL_VALIDATION.md`` section 3.2: agreement
between refinement levels carries no information when **no window opened**, because both levels
are then pure adiabatic transport converging to the same limit.  The gamma rule guards that, but
``GAMMA_TO_ERROR`` is itself good only to ~2x, so a result certified at 0.95 of the bound is a
much closer call than one at 0.01.

Run:  python weak_band.py [n_cases]
"""

import sys
import warnings

import numpy as np

import harness as H
import magnus.adiabatic as ad
import magnus.oscprob as oscprob

TOL = 1.0e-3
BOUND = TOL + TOL          # atol + rtol, the certification budget


def cases(n, seed=20260805):
    """Random smooth profiles in the weak band: N < 25, where the hybrid path answers."""
    rng = np.random.default_rng(seed)
    for _ in range(n):
        d = int(rng.choice([2, 3]))
        energy = float(10.0**rng.uniform(7.3, 8.3))
        span = float(rng.uniform(0.3, 1.0))*H.L_SCALE
        ne = H.fourier_ne(rng, n_modes=int(rng.integers(2, 9)), span=span,
                          base_ratio=float(10.0**rng.uniform(-2.5, -1.0)),
                          amp=float(rng.uniform(0.2, 0.95)))
        N = int(rng.choice([1, 3, 8, 16, 24]))
        Ls = np.linspace(0.05*span, span, N) if N > 1 else np.array([span])
        yield dict(d=d, energy=energy, span=span, ne=ne, N=N, Ls=Ls,
                   params=H.params_for(d))


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    rows = []
    for i, c in enumerate(cases(n)):
        H_of_l = H.H_factory(c['d'], c['params'], H.vcc_of(c['ne']), c['energy'])
        try:
            Pref = np.array([H.P_of(U)
                             for U in H.exact_U_many(H_of_l, 0.0, c['Ls'], c['d'])])
        except Exception:                      # noqa: BLE001
            continue

        info = {}
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always')
                P = np.asarray(oscprob.osc_prob_matter_std_potential(
                    c['d'], c['ne'], c['energy'],
                    c['Ls'] if c['N'] > 1 else float(c['Ls'][0]), c['params'], L0=0.0,
                    density_is_of_number_of_electrons=True, strategy_info=info)
                ).reshape(c['N'], c['d'], c['d'])
        except ValueError:
            # The random Fourier sum went negative somewhere; the package rejecting it is
            # correct validation of the generator's output, not a defect.  Same rate the
            # original fuzzer saw (7 of 250).
            continue
        err = H.maxabs(P - Pref)

        # What the hybrid propagator knew about itself, per point, independently of dispatch.
        gmax, nwin = 0.0, 0
        for L in c['Ls']:
            hinfo = {}
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                ad.hybrid_propagator(H_of_l, 0.0, float(L), rtol=TOL, atol=TOL, info=hinfo)
            gmax = max(gmax, hinfo.get('gamma_max', 0.0))
            nwin = max(nwin, hinfo.get('n_windows', 0))

        rows.append(dict(i=i, d=c['d'], N=c['N'], engine=info['engine'], err=err,
                         gamma_max=gmax, n_windows=nwin,
                         margin=ad.GAMMA_TO_ERROR*gmax/BOUND,
                         warned=bool(caught)))
        if (i + 1) % 20 == 0:
            print('   ... %d/%d' % (i + 1, n), flush=True)

    scored = [r for r in rows if np.isfinite(r['gamma_max'])]
    hyb = [r for r in scored if r['engine'] == 'hybrid']
    silent = [r for r in scored if r['err'] > TOL and not r['warned']]
    print('\nscored %d;  hybrid answered %d;  silent misses %d'
          % (len(scored), len(hyb), len(silent)))

    print('\n--- the hypothesis: is a silent miss a no-window case? ---')
    for tag, pop in (('all hybrid-answered', hyb),
                     ('  of which no window opened', [r for r in hyb if r['n_windows'] == 0]),
                     ('  of which a window opened', [r for r in hyb if r['n_windows'] > 0])):
        if not pop:
            print('%-30s n=0' % tag)
            continue
        e = np.array([r['err'] for r in pop])
        s = [r for r in pop if r['err'] > TOL and not r['warned']]
        print('%-30s n=%-3d median err %.2e  max %.2e  silent %d'
              % (tag, len(pop), np.median(e), e.max(), len(s)))

    print('\n--- trigger cost/benefit: fraction of calls a margin cut would verify ---')
    print('%-14s %10s %10s %12s' % ('margin >', 'triggered', 'of hybrid', 'silent caught'))
    nwin0 = [r for r in hyb if r['n_windows'] == 0]
    sil0 = [r for r in nwin0 if r['err'] > TOL and not r['warned']]
    for f in (0.0, 0.01, 0.05, 0.1, 0.3, 0.5):
        t = [r for r in nwin0 if r['margin'] > f]
        c = [r for r in sil0 if r['margin'] > f]
        print('%-14.2f %10d %9.0f%% %11d/%d'
              % (f, len(t), 100.0*len(t)/max(len(hyb), 1), len(c), len(sil0)))

    if silent:
        print('\nsilent misses in detail:')
        print('%4s %3s %4s %-11s %10s %10s %8s' %
              ('i', 'd', 'N', 'engine', 'err', 'gamma_max', 'margin'))
        for r in sorted(silent, key=lambda r: -r['err']):
            print('%4d %3d %4d %-11s %10.3e %10.3e %8.3f'
                  % (r['i'], r['d'], r['N'], r['engine'], r['err'],
                     r['gamma_max'], r['margin']))
    np.save('weak_band_rows.npy', np.array(rows, dtype=object), allow_pickle=True)


if __name__ == '__main__':
    main()
