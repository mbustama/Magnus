# -*- coding: utf-8 -*-
"""What would an aliasing check COST?  Measured before it is built, not after.

The proposal: warn when a baseline scan samples an oscillation more coarsely than the
oscillation itself, so the returned array is aliased and meaningless *as a curve* -- which is
wrong regardless of what the caller intended, and so is warnable without guessing intent.  The
criterion is Nyquist: the local oscillation length is :math:`2\\pi/\\Delta\\lambda` from the
smallest gap in the Hamiltonian's spectrum; compare it against the scan spacing.

Cost is the question that decides whether it can be on by default, so it is measured first.
Three variants, cheapest to dearest:

  * ``probe8`` / ``probe16`` / ``probe64`` -- evaluate H at k points along the trajectory,
    diagonalize, take the smallest gap.  Standalone cost, nothing reused.
  * ``per_point`` -- diagonalize at every scan point, the naive implementation.

Reported as a fraction of the scan it would be attached to.  Anything under ~1 % is free in
practice; anything above ~5 % has to be opt-in or must reuse eigenvalues the engine already
computed (``adiabatic`` diagonalizes at lines 661/864/918 and ``oscprob`` at 2955/3726, so
reuse is plausible but is a different piece of work).

**Machine discipline**: this is a timing measurement on a machine that may be shared, so it
alternates the variants, carries a control the change cannot touch, and reads **minima** --
interference only ever adds time, so the fastest observation is closest to the true cost.

Run:  python alias_cost.py [rounds]
"""

import sys
import time
import warnings

import numpy as np

import harness as H
import magnus.oscprob as op
import physical_profiles as pp

D = 3
N_SCAN = (2, 8, 32)


def alias_probe(H_of_l, l0, l1, k):
    """The whole proposed check: k samples, smallest spectral gap, shortest oscillation."""
    ls = np.linspace(float(l0), float(l1), k)
    shortest = np.inf
    for l in ls:
        lam = np.linalg.eigvalsh(np.asarray(H_of_l(l)))
        gap = float(np.min(np.diff(np.sort(lam))))
        if gap > 0.0:
            shortest = min(shortest, 2.0*np.pi/gap)
    return shortest


def timeit(fn, repeats=1):
    t0 = time.perf_counter()
    for _ in range(repeats):
        fn()
    return (time.perf_counter() - t0)/repeats


def main():
    rounds = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    params = H.params_for(D)
    fams, seen = [], set()
    for f in pp.families():
        if f['kind'] not in seen:
            seen.add(f['kind'])
            fams.append(f)

    def control():
        return op.osc_prob_matter_std_potential(
            D, 0.05*H.NE0, 30.0e6, np.linspace(0.2, 1.0, 8)*H.L_SCALE, params, L0=0.0,
            density_is_of_number_of_electrons=True)

    acc, ctrl = {}, []
    print('# Cost of the proposed aliasing check, as a fraction of the scan it guards')
    print('# %d profiles x N in %s, %d alternating rounds, minima reported\n'
          % (len(fams), N_SCAN, rounds))
    for rnd in range(rounds):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            t0 = timeit(control)
            for f in fams:
                energy = f['energies'][0]
                H_of_l = H.H_factory(D, params, H.vcc_of(f['ne']), energy)
                for N in N_SCAN:
                    Ls = np.linspace(f['l0'] + 0.2*(f['l1'] - f['l0']), f['l1'], N)

                    def scan(f=f, energy=energy, Ls=Ls):
                        return op.osc_prob_matter_std_potential(
                            D, f['ne'], energy, Ls, params, L0=f['l0'],
                            density_is_of_number_of_electrons=True)

                    key = (f['kind'], N)
                    for name, fn in (
                            ('scan', scan),
                            ('probe8', lambda h=H_of_l, f=f: alias_probe(h, f['l0'], f['l1'], 8)),
                            ('probe16', lambda h=H_of_l, f=f: alias_probe(h, f['l0'], f['l1'], 16)),
                            ('probe64', lambda h=H_of_l, f=f: alias_probe(h, f['l0'], f['l1'], 64)),
                            ('per_point', lambda h=H_of_l, f=f, N=N:
                                alias_probe(h, f['l0'], f['l1'], max(N, 2)))):
                        acc.setdefault(key, {}).setdefault(name, []).append(timeit(fn))
            t1 = timeit(control)
        ctrl.append(t1/t0)
        print('  round %d/%d, control drift %.2fx' % (rnd + 1, rounds, t1/t0), flush=True)

    print('\ncontrol drift: %s' % ', '.join('%.2fx' % c for c in ctrl))
    print('\n%-16s %4s %12s %10s %10s %10s %10s'
          % ('profile', 'N', 'scan ms', 'probe8', 'probe16', 'probe64', 'per_point'))
    worst = {}
    for (kind, N) in sorted(acc):
        row = acc[(kind, N)]
        s = float(np.min(row['scan']))
        cells = []
        for name in ('probe8', 'probe16', 'probe64', 'per_point'):
            frac = float(np.min(row[name]))/s
            worst[name] = max(worst.get(name, 0.0), frac)
            cells.append('%9.3f%%' % (100.0*frac))
        print('%-16s %4d %12.1f %s' % (kind, N, 1e3*s, ' '.join(cells)))

    print('\n=== VERDICT ===')
    for name in ('probe8', 'probe16', 'probe64', 'per_point'):
        v = ('free in practice' if worst[name] < 0.01 else
             'acceptable by default' if worst[name] < 0.05 else
             'too expensive to be on by default')
        print('%-10s worst overhead %7.3f%%   -> %s' % (name, 100.0*worst[name], v))
    return 0


if __name__ == '__main__':
    sys.exit(main())
