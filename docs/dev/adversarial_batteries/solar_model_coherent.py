# -*- coding: utf-8 -*-
r"""Coherent probabilities through the tabulated solar models: right, loud, or both?

The tables are interpolated linearly in log n_e, so the profile has a kink at every row, and
the Bahcall tables (four significant figures) step through the flat core.  The averaged route
does not care; a coherent probability over most of the Sun can.  This measures three things,
for P_ee from the centre, on every case where the default call warned plus a sample of those
where it did not:

  default    what osc_prob_3nu_sun returns with nothing but density_profile set, which engine
             answered, and what it warned;
  reference  strategy='magnus' with a breakpoint at every table row (so no slab straddles a
             kink) at a fixed 800 000 and 1 600 000 slabs; the two must agree before the
             reference is quoted;
  recipe     the call docs/source/solar_models.rst recommends: rows as t_breakpoints,
             strategy='magnus', n_slabs=200_000, max_n_slabs=10_000_000.

Run from the repository root:  python docs/dev/adversarial_batteries/solar_model_coherent.py
"""

import warnings

import numpy as np

import magnus.globaldefs as gd
import magnus.oscprob as op
from magnus import solarmodels as sm

R = gd.SUN_RADIUS*gd.UNIT_KM

# (model, energy [MeV], L [R_sun]).  The first five are every coherent call that warned on the
# 12 models x {1, 5, 10, 20} MeV x 0.9 R_sun grid, plus the one found at 0.95 R_sun; the rest
# certified on the hybrid engine and are here to check that certifying meant being right.
CASES = [('BP04', 10, 0.9), ('BP04', 20, 0.9), ('BS05-OP', 10, 0.9), ('BS05-OP', 20, 0.9),
         ('BS05-AGS-OP', 10, 0.95),
         ('BS05-AGS-OP', 10, 0.9), ('B16-GS98', 20, 0.9), ('BP2000', 20, 0.9),
         ('B23-GS98', 20, 0.9), ('BP04', 5, 0.9), ('B16-GS98', 1, 0.9)]


def _call(record, **kw):
    info = {}
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        P = float(np.asarray(op.osc_prob_3nu_sun(strategy_info=info, **kw)))
    names = sorted({x.category.__name__.replace('Warning', '') for x in w})
    return P, (info.get('engine'), names) if record else names


def main():
    print('%-12s %4s %5s | %-14s %-8s %-36s | %-11s | %-9s %-9s %s' % (
        'model', 'MeV', 'L/R', 'reference', 'ref agr', 'default: engine, warnings', 'default err',
        'recipe err', '', 'recipe warnings'))
    for m, e, f in CASES:
        kw = dict(energy=e*gd.UNIT_MEV, L=f*R, L0=0.0, nu_i=0, nu_f=0, density_profile=m)
        rows = sm.load_solar_model(m)['r_over_r_sun']*R
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ref = [float(np.asarray(op.osc_prob_3nu_sun(strategy='magnus', n_slabs=n,
                                                        max_n_slabs=n, t_breakpoints=rows, **kw)))
                   for n in (800_000, 1_600_000)]
        P, (engine, warned) = _call(True, **kw)
        P_r, (_, warned_r) = _call(True, strategy='magnus', t_breakpoints=rows, n_slabs=200_000,
                                   max_n_slabs=10_000_000, **kw)
        print('%-12s %4d %5.2f | %.10f %.1e %-8s %-27s | %.1e     | %.1e   %s' % (
            m, e, f, ref[1], abs(ref[1] - ref[0]), engine, ','.join(warned) or '-',
            abs(P - ref[1]), abs(P_r - ref[1]), ','.join(warned_r) or '-'), flush=True)


if __name__ == '__main__':
    main()
