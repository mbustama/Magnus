# -*- coding: utf-8 -*-
"""Which features does the averaged route escalate?  (issue #60; SUDDEN_TRANSFER_THRESHOLD)

averaged_probabilities_adiabatic searches for non-adiabatic windows once, on 200 probes.  A
front narrower than their spacing is never examined, and the fully adiabatic answer used to
come back silently -- 0.04 on a supernova shock ray where the averaged probability is 0.37 to
0.59.  Such a front can be seen cheaply: one half of a probe interval carries nearly all of the
interval's change (adiabatic._concentrated_intervals, the first stage of the resolution test).
But a solar-model table interpolated in log-density shows the same shape at every grid point of
its core, where nothing happens, so seeing is not enough.

THE QUANTITY.  For each concentrated interval [l_a, l_b], the most probability an instantaneous
change from H(l_a) to H(l_b) could move between levels,

    max_{i != j} |<v_i(l_a)|v_j(l_b)>|^2 ,

an upper bound on what a monotone passage between the two positions transfers.  The averaged
route escalates when the largest such bound on a path exceeds SUDDEN_TRANSFER_THRESHOLD.

METHOD.  The largest bound, per path, over the populations below.  Whether the averaged answer
was right today comes from the issue-#60 work, where each case was scored against a decohered
reference built without the averaging engine: the resolved instantaneous probability, fronts
declared, averaged over energy and final position (for a solar chord, over impact parameter
and energy, with decohered injection and detection at the chord ends).

  A.  BS05(AGS,OP), cubic and linear, d = 2 and 3, 1-30 MeV.
  B.  Issue #60's shock ray at 15 MeV: fronts 0.07 to 2000 km wide, three placements each.
  C.  Supernova turbulence (the physical-profile families), d = 3, 5-30 MeV.
  D.  Earth crust with its layer edges undeclared, d = 3, 5-30 MeV.

Run from the repository root:  python docs/dev/adversarial_batteries/sudden_transfer_sweep.py
"""
import sys
import warnings

import numpy as np

import harness as H
import physical_profiles as pp
import magnus.globaldefs as gd
from magnus import avgprob, adiabatic

warnings.simplefilter('ignore')
KM = gd.UNIT_KM


def largest_bound(H_func, l0, l1, n_probe=200):
    ls, flagged, _ = adiabatic._concentrated_intervals(H_func, float(l0), float(l1), n_probe)
    return max([avgprob._sudden_transfer(H_func, ls[i], ls[i + 1]) for i in flagged], default=0.0)


def shock_H(w_km, dr_km, energy):
    """Issue #60's shock ray, fronts w_km wide, shifted by dr_km."""
    m_n = 0.5*(gd.MASS_PROTON + gd.MASS_NEUTRON)
    rc, rf = 12348.0 + dr_km, 30323.0 + dr_km

    def smooth(u):
        u = np.clip(np.asarray(u, dtype=float), 0.0, 1.0)
        return u*u*(3.0 - 2.0*u)

    def rare(r, rs):
        u = np.clip(1.0 - np.asarray(r, dtype=float)/rs, 0.0, 1.0)
        return np.exp((0.28 - 0.69*np.log(rs))*np.arcsin(u)**1.1)

    def ne(l):
        r = np.asarray(l, dtype=float)/KM
        f = 1.0 + smooth((rf + 0.5*w_km - r)/w_km)*(10.0*rare(r, rf) - 1.0)
        f = f*(1.0 + smooth((rc + 0.5*w_km - r)/w_km)*1.5)
        out = 1.0e14*r**(-2.4)*f*gd.UNIT_G_PER_CM3/m_n*0.5
        return out[()] if np.ndim(out) == 0 else out
    return H.H_factory(3, H.params_for(3), H.vcc_of(ne), energy), 1.0e4*KM, 8.0e4*KM


def main():
    fams = {f['label']: f for f in pp.families()}
    print('A. BS05 solar model')
    for label in [k for k in fams if k.startswith('BS05')]:
        f = fams[label]
        b = [largest_bound(H.H_factory(d, H.params_for(d), H.vcc_of(f['ne']), E), f['l0'], f['l1'])
             for d in (2, 3) for E in (1e6, 5e6, 10e6, 30e6)]
        print('   %-24s largest bound %.1e' % (label, max(b)))
    print('B. issue #60 shock ray, 15 MeV')
    for w in (0.07, 0.7, 7.0, 20.0, 70.0, 200.0, 700.0, 2000.0):
        b = [largest_bound(*shock_H(w, dr, 15.0e6)) for dr in (0.0, 117.0, 233.0)]
        print('   fronts %7.2f km          largest bounds %s' % (w, ', '.join('%.1e' % x for x in b)))
    for tag, prefix in (('C. supernova turbulence', 'SN turbulence'), ('D. Earth crust, edges undeclared', 'Earth crust')):
        print(tag)
        for label in [k for k in fams if k.startswith(prefix)]:
            f = fams[label]
            b = [largest_bound(H.H_factory(3, H.params_for(3), H.vcc_of(f['ne']), E), f['l0'], f['l1'])
                 for E in (5e6, 10e6, 30e6)]
            print('   %-24s largest bounds %s' % (label, ', '.join('%.1e' % x for x in b)))
    print('\nthreshold: %.0e' % avgprob.SUDDEN_TRANSFER_THRESHOLD)
    return 0


if __name__ == '__main__':
    sys.exit(main())
