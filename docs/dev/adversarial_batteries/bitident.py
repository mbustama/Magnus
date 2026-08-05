# -*- coding: utf-8 -*-
"""Battery 1.1 -- the bit-identity set.

These nine configurations must be BIT-IDENTICAL to main; anything else is a regression to
be justified.  Run under main's src and under the branch's src, dump raw probabilities,
diff element-wise.

Caveat carried from the handover: commit 9c7945a (the gamma sweep) landed AFTER the last
check and touches the hybrid path, so the single-point and sub-threshold rows may
legitimately have moved.  This script establishes WHICH moved; the write-up justifies each.
"""

import sys
import warnings

import numpy as np

warnings.simplefilter('ignore')

import magnus.globaldefs as gd            # noqa: E402
import magnus.matter as matter            # noqa: E402
import magnus.oscprob as oscprob          # noqa: E402

L_SCALE = gd.L_SCALE_SUN
R_SUN = gd.SUN_RADIUS*gd.UNIT_KM
NE0 = gd.NUM_DENSITY_E_SUN_CENTRAL
E = 50.0e6
LM = 0.5*R_SUN

_p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
P3 = {k: _p[k] for k in ('s12', 's23', 's13', 'dCP', 'D21', 'D31')}
P2 = {'sth': _p['s12'], 'Dm2': _p['D21']}

solar = matter.exp_density_profile(NE0, L_SCALE)
Ls60 = np.linspace(0.05*LM, LM, 60)
Ls8 = np.linspace(0.05*LM, LM, 8)


def std(*a, **k):
    return np.asarray(oscprob.osc_prob_matter_std_potential(
        *a, L0=0.0, density_is_of_number_of_electrons=True, **k))


CASES = {}
CASES["1 strategy='magnus' scan (N=60)"] = lambda: std(2, solar, E, Ls60, P2,
                                                       strategy='magnus')
CASES['2 sub-threshold scan (N=8)'] = lambda: std(2, solar, E, Ls8, P2)
CASES['3 single point'] = lambda: std(2, solar, E, LM, P2)
CASES['4 vacuum scan (N=60)'] = lambda: std(3, 0.0, E, Ls60, P3)
CASES['5 constant-density scan (N=60)'] = lambda: std(3, 0.05*NE0, E, Ls60, P3)
CASES['6 energy scan at fixed baseline'] = lambda: std(
    2, solar, np.linspace(10e6, 100e6, 40), LM, P2)
_E_earth = np.linspace(1.0, 20.0, 12)*gd.UNIT_GEV
CASES['7 osc_prob_earth PREM'] = lambda: np.asarray(oscprob.osc_prob_3nu_earth(
    _E_earth, costhz=-0.8, L=np.full(12, 2.0*6371.0*0.8)*gd.UNIT_KM,
    validate_input=False))
CASES['8 average=True'] = lambda: std(2, solar, E, Ls60, P2, average=True)
CASES['9 explicit cumulative=False'] = lambda: np.asarray(
    oscprob.osc_prob_energy_baseline(
        (lambda enu, l: (1.0/enu)*np.asarray(
            oscprob.hamiltonians.hamiltonian_2nu_vacuum_energy_independent(
                P2['sth'], P2['Dm2']), dtype=complex)
         + np.asarray(matter.vcc_func_from_rho_func(
             solar, 0.0, 1.0, 0.5, nubar=False, density_matter_is_in_g_per_cm3=False,
             density_is_of_number_of_electrons=True)(l))[..., None, None]
         * np.diag([1.0, 0.0]).astype(complex)),
        E, Ls60, 0.0, cumulative=False))
# Extra rows the handover's list does not name but the branch could plausibly move:
CASES['10 3nu solar scan N=60 (cumulative)'] = lambda: std(3, solar, E, Ls60, P3)
CASES['11 nubar scan N=60'] = lambda: std(2, solar, E, Ls60, P2, nubar=True)
# Row 12 MOVES BY DESIGN, and exists because rows 1-11 cannot see the change that matters.
# Row 7 is the only Earth case and it is a *scan*, which the separable engine answers without
# ever reaching magnus_expansion_multislab -- so the palindromic mirror leaves it bit-identical
# and the battery would report "0 of 11 moved" while the optimisation was live and moving
# single points.  A single Earth point takes the general Magnus ladder, where the mirror does
# fire.  Measured movement across 15 (costhz, energy) configurations: worst 8.6e-15 relative,
# typical 2e-15.  Set magnus.magnus.USE_PALINDROME = False to recover the pre-mirror numbers.
_c_earth = -0.8
CASES['12 osc_prob_earth PREM single point (MOVES)'] = lambda: np.asarray(
    oscprob.osc_prob_3nu_earth(2.0*gd.UNIT_GEV, costhz=_c_earth,
                               L=2.0*6371.0*abs(_c_earth)*gd.UNIT_KM,
                               validate_input=False))

out = {}
for name, fn in CASES.items():
    try:
        out[name] = np.asarray(fn(), dtype=float)
    except Exception as ex:                # noqa: BLE001
        out[name] = np.array([float('nan')])
        print('  %s RAISED %s: %s' % (name, type(ex).__name__, str(ex)[:120]),
              file=sys.stderr)

np.savez(sys.argv[1], **{k: v for k, v in out.items()})
print('wrote %s (%d cases)' % (sys.argv[1], len(out)))
