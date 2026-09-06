# -*- coding: utf-8 -*-
r"""sterile_projector_check.py

Does the sterile block of the matter projector follow the Earth's layered
composition, and what does it cost when it does not?

WHY THIS EXISTS AS A SCRIPT.  The number that motivated issue #47 -- an error of
about 0.4 in probability near the sterile matter resonance on a core-crossing
chord -- lived only in a prose assessment, and could not be reproduced from it
without rebuilding the comparison by hand.  Rebuilding it by hand went wrong
twice, once on a mixing angle and once on the baseline's units, both of which
return perfectly unitary and perfectly wrong answers.  This is that comparison,
runnable::

    python3 notebooks/sterile_projector_check.py

WHAT IT COMPARES.  Three arms, differing *only* in the matter projector.  All
three use the same low-level solver, the same layered-`Y_e` density, and the same
PREM layer edges as slab boundaries:

  truth    P(l) built from r(l) = (1 - Y_e(l)) / Y_e(l), the same composition the
           density uses.  On PREM this is both the "position-dependent P" and the
           "per-layer P" of issue #47: the slab edges already sit on every Y_e
           boundary, so the two coincide.
  scalar   P built once from r = 1.0, which is what the Earth wrappers used
           before the fix.
  wrapper  `osc_prob_4nu_earth` at its default, to show which of the two the
           shipped entry point agrees with.

WHAT THE RESIDUAL MEANS.  The wrapper agrees with the truth arm to about 2e-05,
and that figure is flat from rtol 1e-6 through 1e-10, so it is not discretisation.
It is a small construction difference: this script's truth arm passes a fixed
ratio to `vcc_func_from_rho_func` for the mass-density conversion, where the
wrapper threads the resolved r(l) through it as well.  Four orders below the
defect being measured, so it does not affect the comparison -- but it is the hand
-built arm that is approximate here, not the shipped one.

BASELINE UNITS.  `gd.EARTH_RADIUS` is in kilometres and `osc_prob` wants inverse
eV; the chord length is converted here.  Passing kilometres yields the identity
matrix -- unitary, converged and meaningless.
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

import pathlib
import sys
import warnings

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]/'src'))

import magnus.globaldefs as gd                            # noqa: E402
import magnus.hamiltonians as hams                        # noqa: E402
import magnus.matter as matter                            # noqa: E402
import magnus.oscprob as oscprob                          # noqa: E402

OSC = dict(s12=np.sqrt(0.307), s23=np.sqrt(0.55), s13=np.sqrt(0.022),
           dCP=1.36*np.pi, D21=7.42e-5, D31=2.51e-3)
STERILE = dict(s14=0.15, s24=0.10, s34=0.0, D41=1.0)
COSTHZ, ENERGY = -1.0, 2.4e12                             # vertical chord, 2.4 TeV
L = 2.0*gd.EARTH_RADIUS*gd.UNIT_KM                        # km -> eV^-1


def _pieces():
    """The density and the resolved neutron-to-proton ratio the wrappers use."""
    out = oscprob._earth_composition(
        COSTHZ, None, None, None, None, None, None, 'sterile_projector_check',
        num_flavors=4)
    rho_func, r_resolved = out if isinstance(out, tuple) else (out, 1.0)
    vcc = matter.vcc_func_from_rho_func(rho_func, 0.0, 1.0, 0.5, False, False, True)
    return vcc, r_resolved


def _probability(vcc, projector):
    """P through the general ladder, projector the only thing that varies."""
    # The signature interleaves each angle with its CP phase -- s12, s23, s13,
    # dCP, s14, d14, s24, d24, s34, then the splittings.  Grouping the angles
    # together puts a phase in a sine slot; the package catches it, but only
    # because a phase happens to exceed one.
    hvac = np.asarray(hams.hamiltonian_4nu_vacuum_energy_independent(
        OSC['s12'], OSC['s23'], OSC['s13'], OSC['dCP'],
        STERILE['s14'], 0.0, STERILE['s24'], 0.0, STERILE['s34'],
        OSC['D21'], OSC['D31'], STERILE['D41']))

    def H_of_l(l):
        v = np.asarray(vcc(l))
        P = projector(l) if callable(projector) else projector
        return hvac/ENERGY + v[..., None, None]*P if v.ndim else hvac/ENERGY + v*P

    info = {}
    with oscprob._engine_probe(info=info):
        P = np.asarray(oscprob.osc_prob(H_of_l, 0.0, L, rtol=1e-6, atol=1e-8,
                                        validate_input=False))
    return P, info.get('engine')


def main():
    warnings.simplefilter('ignore')
    vcc, r_resolved = _pieces()
    print('resolved r is %s' % ('a callable r(l)' if callable(r_resolved)
                                else 'the scalar %r' % (r_resolved,)))

    P_truth, e1 = _probability(vcc, matter.matter_potential_projector(4, r_resolved))
    P_scalar, e2 = _probability(vcc, matter.matter_potential_projector(4, 1.0))
    P_wrap = np.asarray(oscprob.osc_prob_4nu_earth(
        ENERGY, costhz=COSTHZ, L=L, **OSC, d14=0.0, d24=0.0, **STERILE,
        rtol=1e-6, atol=1e-8, validate_input=False))

    def worst(a, b):
        return float(np.max(np.abs(np.asarray(a)[:3, :3] - np.asarray(b)[:3, :3])))

    print('engines: truth=%s scalar=%s' % (e1, e2))
    print('  scalar (r = 1.0)   vs truth : %.4e   <- the defect' % worst(P_scalar, P_truth))
    print('  shipped wrapper    vs truth : %.4e' % worst(P_wrap, P_truth))
    print('  shipped wrapper    vs scalar: %.4e' % worst(P_wrap, P_scalar))
    print('\nThe wrapper agrees with whichever arm its default projector matches.')


if __name__ == '__main__':
    main()
