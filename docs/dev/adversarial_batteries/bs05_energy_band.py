# -*- coding: utf-8 -*-
"""Does the BS05 silent miss survive at energies a solar neutrino actually has?

``warn_fp.py --physical`` found a silent miss on the cubic-spline interpolation of the real
BS05(AGS,OP) solar model at **100 MeV**, and ``attribute_physical.py`` found a worse one, still
silent, at 130 MeV (6.39e-03).  Neither is a solar-neutrino energy: the pp and 8B spectra stop
around 15-20 MeV.

The population uses 30 and 100 MeV on solar profiles for one reason -- ``solve_ivp`` is what made
the first version of ``warn_fp.py`` unfinishable, and it is cheapest at high energy -- so the
band that matters most is the one the population could least afford.  **A finding at an energy
nobody computes is not a finding about reachability**, so this sweeps down into the real band and
reports where the miss starts and stops.

The oracle is verified at each energy rather than assumed: at low energy the vacuum term
``h_vac/E`` grows and the integration gets harder, which is exactly where a quietly
non-converged oracle would invent an error.

Run:  python bs05_energy_band.py
"""

import sys
import warnings

import numpy as np

import harness as H
import magnus.oscprob as oscprob
import physical_profiles as pp

TOL = 1e-3
ENERGIES = (5.0e6, 10.0e6, 15.0e6, 20.0e6, 30.0e6, 50.0e6, 100.0e6)
SOLAR_MAX = 20.0e6
"""float: Top of the solar-neutrino band, for the verdict below.  The 8B spectrum ends near
15 MeV and hep near 18.8; 20 MeV is a generous ceiling."""


def one(ne, d, params, energy, l0, l1):
    info = {}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        P = np.asarray(oscprob.osc_prob_matter_std_potential(
            d, ne, energy, l1, params, L0=l0,
            density_is_of_number_of_electrons=True, strategy_info=info))
    return P, info, sorted({w.category.__name__ for w in caught})


def main():
    print('# Does the BS05 silent miss reach the solar-neutrino energy band?')
    print('# tolerance %.0e; trajectory = 1 solar scale height\n' % TOL)
    print('%-22s %3s %8s %12s %9s %10s  %s'
          % ('profile', 'd', 'E/MeV', 'err', 'certified', 'oracle mv', 'warns'))
    silent = []
    for kind in ('cubic', 'linear'):
        fam = next(f for f in pp.families() if f['label'] == 'BS05(AGS,OP) ' + kind)
        ne, l0, l1 = fam['ne'], fam['l0'], fam['l1']
        for d in (2, 3):
            params = H.params_for(d)
            for energy in ENERGIES:
                H_of_l = H.H_factory(d, params, H.vcc_of(ne), energy)
                try:
                    Pref = H.P_of(H.exact_U(H_of_l, l0, l1, d))
                    P, info, warns = one(ne, d, params, energy, l0, l1)
                except Exception as exc:               # noqa: BLE001
                    print('%-22s %3d %8.1f  %s' % (kind, d, energy/1e6, type(exc).__name__))
                    continue
                err = H.maxabs(P - Pref)
                moved, _ = H.oracle_converged(H_of_l, l0, l1, d, max(err, TOL))
                tag = ''
                if err > TOL:
                    tag = '   <-- SILENT' if not warns else '   <-- warned'
                    if not warns:
                        silent.append((err, kind, d, energy))
                print('%-22s %3d %8.1f %12.3e %9s %10.1e  %s%s'
                      % (kind, d, energy/1e6, err, info.get('certified'), moved,
                         ','.join(warns) or '-', tag), flush=True)

    # The question is "does ANY silent miss fall in the solar band", not "is the LARGEST one in
    # it".  The first version of this reported the worst miss and then tested that one row for
    # band membership, which answered a different question and got the verdict backwards: the
    # largest miss sits at 100 MeV, but there are two more at 5 MeV.
    print('\n=== VERDICT ===')
    if not silent:
        print('no silent miss anywhere in %s MeV' % ([e/1e6 for e in ENERGIES],))
        return 0
    print('silent misses (no warning, outside %.0e):' % TOL)
    for err, kind, d, energy in sorted(silent, key=lambda r: -r[0]):
        print('    %-8s d=%d E=%6.1f MeV  err=%.3e%s'
              % (kind, d, energy/1e6, err, '   <-- solar band' if energy <= SOLAR_MAX else ''))
    in_band = [r for r in silent if r[3] <= SOLAR_MAX]
    print('\n%d of %d silent misses are at or below %.0f MeV' %
          (len(in_band), len(silent), SOLAR_MAX/1e6))
    print('reaches a real solar-neutrino energy: %s'
          % ('YES' if in_band else 'NO -- only above the solar band'))
    if in_band:
        print('  worst in band: %s d=%d E=%.0f MeV err=%.3e'
              % (in_band[0][1], in_band[0][2], in_band[0][3]/1e6, in_band[0][0]))
    return 0


if __name__ == '__main__':
    sys.exit(main())
