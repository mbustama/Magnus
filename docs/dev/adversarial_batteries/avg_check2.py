# -*- coding: utf-8 -*-
"""The averaging check of ``avg_check.py``, applied to the supernova families.

§13.17 showed the 5 MeV solar "silent miss" is mostly phase, and phase is unobservable: the
trajectory is 398 oscillation lengths and everything downstream averages it away.  The same
argument applies to a supernova -- the neutrino's phase at the stellar surface is not knowable
to anything like the precision that would make it observable, and detector energy resolution
finishes the job.  So the supernova numbers reported in §13.4 need the same test before they can
be called defects:

  * ``SN shock w=1e-3``  d = 3, 15 MeV -- the SILENT miss, 1.095e-03 instantaneous;
  * ``SN shock w=1e-6``  d = 3, 15 MeV -- the LOUD case, error up to 0.203.  If that averages
    away too, the headline supernova number is a phase artefact as well;
  * ``SN turbulence``    d = 3 -- broadband roughness, up to 1.39e-02.

The window is six oscillation lengths of the *shortest* relevant oscillation, which at d = 3 is
set by D31 rather than D21.

Run:  python avg_check2.py
"""
import sys
import time
import warnings

import numpy as np

import harness as H
import magnus.globaldefs as gd
import magnus.oscprob as op
import physical_profiles as pp

TOL = 1e-3


def check(label, fam_label, d, energy, n_pts):
    fam = next(f for f in pp.families() if f['label'] == fam_label)
    l0, l1 = fam['l0'], fam['l1']
    params = H.params_for(d)
    p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
    dm2 = p['D31'] if d >= 3 else p['D21']
    l_osc = 4.0*np.pi*energy/dm2
    Ls = np.linspace(max(l1 - 6.0*l_osc, l0 + 0.5*(l1 - l0)), l1, n_pts)
    t0 = time.time()
    Hf = H.H_factory(d, params, H.vcc_of(fam['ne']), energy)
    Pref = np.array([H.P_of(U) for U in H.exact_U_many(Hf, l0, Ls, d)])
    got, warned = [], False
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        for L in Ls:
            got.append(np.asarray(op.osc_prob_matter_std_potential(
                d, fam['ne'], energy, float(L), params, L0=l0,
                density_is_of_number_of_electrons=True)))
        warned = len(caught) > 0
    got = np.array(got)
    inst = float(np.max(np.abs(got[-1] - Pref[-1])))
    inst_worst = float(np.max(np.abs(got - Pref)))
    avg = float(np.max(np.abs(got.mean(axis=0) - Pref.mean(axis=0))))
    print('\n%s  (d=%d, %.0f MeV, %d points over %.1f osc lengths, %.0f s)'
          % (label, d, energy/1e6, n_pts, 6.0, time.time() - t0))
    print('   trajectory / oscillation length      : %.0f' % ((l1 - l0)/l_osc))
    print('   instantaneous at nominal baseline    : %.3e %s'
          % (inst, '(outside)' if inst > TOL else ''))
    print('   worst instantaneous over the window  : %.3e' % inst_worst)
    print('   error in the AVERAGE                 : %.3e %s   <-- the observable'
          % (avg, '(OUTSIDE)' if avg > TOL else '(inside)'))
    print('   averaging reduces the error by       : %.0fx' % (inst_worst/max(avg, 1e-30)))
    print('   any warning raised                   : %s' % warned)
    return dict(label=label, inst=inst, inst_worst=inst_worst, avg=avg, warned=warned)


def main():
    print('# Does the supernova error survive phase averaging?  tolerance %.0e' % TOL)
    rows = [
        check('SN shock w=1e-3  (the SILENT miss)', 'SN shock w=1e-03', 3, 15.0e6, 81),
        check('SN shock w=1e-6  (the LOUD 0.203)', 'SN shock w=1e-06', 3, 15.0e6, 81),
        check('SN turbulence C*=0.1', 'SN turbulence C*=0.1', 3, 45.0e6, 25),
    ]
    print('\n=== VERDICT ===')
    for r in rows:
        verdict = ('still a defect on the averaged observable' if r['avg'] > TOL
                   else 'phase artefact -- inside tolerance once averaged')
        print('%-38s avg %.3e  %s%s' % (r['label'], r['avg'], verdict,
                                        '' if r['warned'] else '   [SILENT]'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
