# -*- coding: utf-8 -*-
"""Would an aliasing warning fire usefully, or on everything?  Measured before it is built.

`FINDINGS_ROBUSTNESS_PROGRAMME.md` §5c: a warning whose false-positive rate cannot be
established is not shipped.  This establishes it for the proposed aliasing check *before* any
of it is written, because there is a specific way for the idea to fail that has nothing to do
with the implementation.

The criterion is Nyquist and it is objectively right: to represent an oscillation of wavelength
:math:`\\lambda` a scan needs spacing below :math:`\\lambda/2`, and above that the returned array
cannot represent the oscillation -- any plot or interpolation of it is an artefact, regardless of
what the caller intended.

**The risk is not that it fires wrongly, but that it fires on everything.**  A solar trajectory
is 398 oscillation lengths and a supernova ray about 4700, so a scan would need thousands of
points to be Nyquist-adequate.  If realistic scans are essentially all undersampled, a warning
that says so is noise however correct it is, and the right vehicle for the information is a
``strategy_info`` statistic that costs nothing and nags nobody.

Reports, over the physical families x realistic scan sizes:

  * the fraction of scans the criterion would fire on -- **the number that decides it**;
  * how many points Nyquist would actually require, per family;
  * how often the cheap 8-point probe disagrees with a 4096-point reference about the verdict,
    which is the only genuine false-positive channel once the criterion is fixed.

Run:  python alias_fp.py
"""

import sys

import numpy as np

import harness as H
import physical_profiles as pp

D = 3
N_SCAN = (2, 4, 8, 16, 32, 64, 128, 256, 1024)


def fastest_oscillation(H_of_l, l0, l1, k):
    """Shortest oscillation length on the path: 2*pi / the LARGEST eigenvalue spread.

    Aliasing is set by the fastest oscillation, which comes from the largest spread -- not from
    the smallest adjacent gap, which is the *longest* oscillation and is what the first version
    of this measured.  The largest spread tracks the matter potential and varies smoothly, which
    is why 8 samples suffice; the smallest gap has a sharp minimum at a resonance and would not.
    """
    ls = np.linspace(float(l0), float(l1), k)
    out = np.inf
    for l in ls:
        lam = np.linalg.eigvalsh(np.asarray(H_of_l(l)))
        spread = float(lam.max() - lam.min())
        if spread > 0.0:
            out = min(out, 2.0*np.pi/spread)
    return out


def main():
    params = H.params_for(D)
    print('# Would an aliasing warning fire usefully, or on everything?')
    print('# Nyquist: a scan is aliased when its spacing exceeds lambda_min / 2\n')
    print('%-15s %8s %12s %14s   %s'
          % ('family', 'E/MeV', 'lambda_min', 'N for Nyquist', 'fires at N = ' + str(N_SCAN)))
    fired = total = 0
    disagree = 0
    seen, rows = set(), []
    for f in pp.families():
        if f['kind'] in seen:
            continue
        seen.add(f['kind'])
        energy = f['energies'][0]
        H_of_l = H.H_factory(D, params, H.vcc_of(f['ne']), energy)
        lam8 = fastest_oscillation(H_of_l, f['l0'], f['l1'], 8)
        lam_ref = fastest_oscillation(H_of_l, f['l0'], f['l1'], 4096)
        span = f['l1'] - f['l0']
        n_needed = int(np.ceil(2.0*span/lam_ref)) + 1
        marks = []
        for N in N_SCAN:
            spacing = span/max(N - 1, 1)
            v8 = spacing > 0.5*lam8
            vref = spacing > 0.5*lam_ref
            marks.append('Y' if v8 else '.')
            fired += int(v8)
            total += 1
            disagree += int(v8 != vref)
        rows.append((f['kind'], n_needed))
        print('%-15s %8.0f %12.4e %14d   %s'
              % (f['kind'], energy/1e6, lam_ref, n_needed, ' '.join(marks)))

    print('\n=== THE DECIDING NUMBER ===')
    print('scans the criterion would fire on : %d of %d = %.0f%%' % (fired, total,
                                                                     100.0*fired/total))
    print('8-point probe disagrees with the 4096-point reference on the verdict: %d of %d'
          % (disagree, total))
    print('\npoints a scan would need to be Nyquist-adequate:')
    for kind, n in rows:
        print('   %-15s %10d' % (kind, n))
    print('\nverdict: %s' % (
        'fires on essentially everything -- ship a strategy_info statistic, NOT a warning'
        if fired/total > 0.9 else
        'fires selectively -- a warning is defensible, subject to the disagreement count'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
