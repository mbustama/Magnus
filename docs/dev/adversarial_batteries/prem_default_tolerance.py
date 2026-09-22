# -*- coding: utf-8 -*-
"""What the default tolerance actually delivers on an Earth crossing.

index.rst and methodology.rst both used to say "about 5e-4 on Earth crossings,
verified against 1e-7-tolerance references", and neither named a measurement.
Nothing in the repository produced that figure.  This is the sweep that replaced
it: the same entry point at the default rtol=atol=1e-3 against itself at 1e-7,
max |dP| over the whole 3x3 probability matrix, across eight chords from grazing
to core-crossing and six energies from 0.5 to 20 GeV.

Read it for what it is.  The reference is Magnus at a tighter tolerance, so this
measures a self-convergence gap and not an error against an independent truth --
which is what the claim it replaces asserted, and is why the number belongs
beside diagnostics.rst's oracle-scored populations rather than instead of them.

Result on 2026-09-22: median 9.24e-07, p90 1.19e-04, worst 2.22e-03 at
costhz = -1.00 and 0.5 GeV -- outside the 1e-3 that was requested.  A spread of
three orders, which is why the pages now give the shape rather than one number.
"""

import warnings

import numpy as np

import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.oscprob as op

COS = [-1.0, -0.95, -0.9, -0.8, -0.6, -0.4, -0.2, -0.05]
ENERGIES = np.array([0.5, 1.0, 2.0, 5.0, 10.0, 20.0])


def main():
    warnings.simplefilter('ignore')
    rows = []
    for costhz in COS:
        # distance_traveled_inside_earth returns km; every osc_prob baseline is
        # in eV^-1, and passing the raw value converges on a chord a few meters
        # long without complaining.
        L = earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM
        E = ENERGIES*gd.UNIT_GEV
        ref = np.asarray(op.osc_prob_3nu_earth(E, costhz=costhz, L=L,
                                               rtol=1e-7, atol=1e-7))
        got = np.asarray(op.osc_prob_3nu_earth(E, costhz=costhz, L=L))
        d = np.abs(got - ref)
        rows.extend((costhz, e, d[k].max()) for k, e in enumerate(ENERGIES))
        print('costhz %+5.2f  chord %7.0f km   max|dP| = %.3e'
              % (costhz, L/gd.UNIT_KM, d.max()), flush=True)

    v = np.array([r[2] for r in rows])
    print('\n%d configurations (%d chords x %d energies, full 3x3 matrix)'
          % (len(v), len(COS), len(ENERGIES)))
    print('  median %.2e   p90 %.2e   max %.2e'
          % (np.median(v), np.percentile(v, 90), v.max()))
    print('  worst: costhz %+.2f at %.1f GeV' % rows[int(np.argmax(v))][:2])
    print('  above 5e-4: %.1f%%   above the requested 1e-3: %.1f%%'
          % (100.0*(v > 5e-4).mean(), 100.0*(v > 1e-3).mean()))


if __name__ == '__main__':
    main()
