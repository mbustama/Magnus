# -*- coding: utf-8 -*-
"""Attribution for the silent misses the physical population produced.

P1 of the physical-profile brief says a silent miss "must be attributed: which engine answered,
which mechanism, and whether ``t_breakpoints`` cures it".  ``warn_fp.py --physical`` found two
in 195 configurations:

  * ``BS05(AGS,OP) cubic``  d = 3, E = 100 MeV, N = 1, err = 1.444e-03
  * ``SN shock w=1e-03``    d = 3, E =  15 MeV, N = 1, err = 1.095e-03

Both are marginal -- 1.4x and 1.1x the requested tolerance, against the 0.54 and 2.9e-02 the
adversarial constructions produced -- but marginal is not the same as absent, and a silent miss
is the one failure mode this programme treats as unacceptable.

This asks, for each: which engine answered and what did the others say; what does
``strategy_info`` report; is it the same at neighbouring energies (a knife-edge or a plateau);
and does declaring ``t_breakpoints`` cure it.

Run:  python attribute_physical.py
"""

import sys
import warnings

import numpy as np

import harness as H
import magnus.oscprob as oscprob
import physical_profiles as pp
from fallback_quality import FORCE

TOL = 1e-3


def _one(ne, d, params, energy, l0, l1, kw=None, forced=None):
    """One call, returning (P, engine, warnings) or (None, reason, [])."""
    kw = dict(kw or {})
    disabled = ()
    if forced is not None:
        strategy, cumulative, disabled = FORCE[forced]
        kw.update(strategy=strategy, cumulative=cumulative)
    info = {}
    try:
        with oscprob._engine_probe(disabled), warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            P = np.asarray(oscprob.osc_prob_matter_std_potential(
                d, ne, energy, l1, params, L0=l0,
                density_is_of_number_of_electrons=True, strategy_info=info, **kw))
    except Exception as exc:                       # noqa: BLE001 -- a decline
        return None, type(exc).__name__, [], {}
    if forced is not None and info.get('engine') != forced:
        return None, 'declined (%s answered)' % info.get('engine'), [], info
    return P, info.get('engine'), sorted({w.category.__name__ for w in caught}), info


def attribute(label, family_label, d, energy, breakpoints=None):
    fam = next(f for f in pp.families() if f['label'] == family_label)
    ne, l0, l1 = fam['ne'], fam['l0'], fam['l1']
    params = H.params_for(d)
    Href = H.H_factory(d, params, H.vcc_of(ne), energy)
    Pref = H.P_of(H.exact_U(Href, l0, l1, d))

    print('\n' + '=' * 78)
    print('%s   d=%d  E=%.0f MeV  single point' % (label, d, energy/1e6))
    print('=' * 78)

    P, engine, warns, info = _one(ne, d, params, energy, l0, l1)
    err = H.maxabs(P - Pref)
    print("  'auto'      engine=%-11s err=%.3e   warns=%s"
          % (engine, err, ','.join(warns) or 'NONE  <-- silent'))
    print('  strategy_info: %s'
          % {k: v for k, v in sorted(info.items()) if k not in ('engine',)})

    print('\n  every engine forced:')
    for name in FORCE:
        Pf, who, wf, _ = _one(ne, d, params, energy, l0, l1, forced=name)
        if Pf is None:
            print('    %-11s %s' % (name, who))
        else:
            print('    %-11s err=%.3e   warns=%s'
                  % (name, H.maxabs(Pf - Pref), ','.join(wf) or '-'))

    # Knife-edge or plateau?  A single point that happens to land just outside is a different
    # problem from a band that is outside.
    print('\n  neighbouring energies (is it a knife-edge?):')
    for scale in (0.7, 0.85, 1.0, 1.15, 1.3):
        e = energy*scale
        Pe, eng, we, _ = _one(ne, d, params, e, l0, l1)
        ref = H.P_of(H.exact_U(H.H_factory(d, params, H.vcc_of(ne), e), l0, l1, d))
        ee = H.maxabs(Pe - ref)
        print('    E=%6.1f MeV  err=%.3e  %-11s %s%s'
              % (e/1e6, ee, eng, ','.join(we) or 'silent',
                 '   <-- outside' if ee > TOL else ''))

    # Does declaring the structure cure it?
    if breakpoints is not None:
        Pb, engb, wb, _ = _one(ne, d, params, energy, l0, l1,
                               kw={'t_breakpoints': breakpoints})
        print('\n  with t_breakpoints declared:')
        print('    engine=%-11s err=%.3e  (was %.3e, %.0fx better)  warns=%s'
              % (engb, H.maxabs(Pb - Pref), err,
                 err/max(H.maxabs(Pb - Pref), 1e-300), ','.join(wb) or '-'))

    # For the tabulated families: does a tighter request fix it?  A silent miss that a tighter
    # tolerance cures is a calibration problem; one it does not cure is structural.
    print('\n  tighter request (rtol=atol=1e-6):')
    Pt, engt, wt, _ = _one(ne, d, params, energy, l0, l1, kw={'rtol': 1e-6, 'atol': 1e-6})
    print('    engine=%-11s err=%.3e  warns=%s'
          % (engt, H.maxabs(Pt - Pref), ','.join(wt) or '-'))


def main():
    attribute('SILENT MISS 1: a real solar model, cubic-spline interpolated',
              'BS05(AGS,OP) cubic', 3, 100.0e6)
    attribute('SILENT MISS 2: supernova shock, 70 km fronts',
              'SN shock w=1e-03', 3, 15.0e6,
              breakpoints=pp.sn_shock_breakpoints(1.0e-3))
    # The linear interpolation of the same model is clean; worth showing side by side, because
    # it localises the defect to the interpolant rather than to the model.
    attribute('CONTROL: the same solar model, linearly interpolated',
              'BS05(AGS,OP) linear', 3, 100.0e6)
    return 0


if __name__ == '__main__':
    sys.exit(main())
