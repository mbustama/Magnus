# -*- coding: utf-8 -*-
"""How small can a gap be before eigh cannot tell it from zero?  (issue #59)

_point_adiabaticity scores a resonance candidate by gamma = |<v_j|dH/dl|v_k>| / gap^2, and
meant to score an exact crossing as infinite -- it returned inf when ``gap == 0``.  But an
exact crossing sampled in floating point does not give a gap of bitwise zero.  Bisection puts
the candidate within an ulp of the crossing and eigh returns the two eigenvalues to within a
few eps*max|lambda| of each other, so the gap comes out at round-off -- 2e-27 eV on the call of
issue #59 -- and when the two levels are decoupled the coupling is exactly zero, so gamma is
0/(2e-27)^2 = 0, no window opens, and adiabatic transport carries a state through the crossing
onto the wrong level.  Whether a given call escaped depended on the last bit.

THE QUANTITY.  The gap at a candidate in units of eigh's resolution,

    u  =  gap / (eps * max|lambda|) ,

which does not depend on the units H is written in.  DEGENERACY_ULPS must sit above u at every
exact crossing (so that each one is treated as the degeneracy it is) and below u at every
candidate that is not a crossing (so that nothing else changes).

METHOD.  Every candidate find_resonance_candidates returns, classified by whether the sorted
pair swaps identity across it: an exact crossing hands the eigenvector of index j just before
the candidate to index k just after.  The Hamiltonian is the one the hybrid engine itself is
handed, captured at the call.

  A.  Exact crossings: the three 4nu Sun wrappers and the 5nu one at their defaults (sterile
      decoupled), nu and nubar, 1-300 MeV, four values of n_n/n_p; the 3nu Sun wrapper with
      s13 = 0, where nu_3 decouples and crosses the nu_e level above ~70 MeV; and a two-level
      exact crossing written in four sets of units, 1e-20 to 1e+10.
  B.  Everything else those calls return, plus ordinary configurations -- 2nu, 3nu and coupled
      4nu (s14 = 0.1, D41 = 1 eV^2) Sun calls over the same energies.

Run from the repository root:  python docs/dev/adversarial_batteries/degenerate_gap_sweep.py
"""
import sys
import warnings

import numpy as np

import magnus.globaldefs as gd
import magnus.oscprob as op
from magnus import adiabatic

warnings.simplefilter('ignore')
EPS = np.finfo(float).eps
L_HALF = 0.5*gd.SUN_RADIUS*gd.UNIT_KM
P0 = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']


def captured(call):
    """The (H_func, l0, l1) the hybrid engine is handed by ``call``."""
    cap = {}
    orig = adiabatic.hybrid_propagator

    def spy(H, l0, l1, *a, **k):
        cap.setdefault('h', (H, l0, l1))
        return orig(H, l0, l1, *a, **k)
    adiabatic.hybrid_propagator = spy
    try:
        call()
    finally:
        adiabatic.hybrid_propagator = orig
    return cap.get('h')


def classify(H, l0, l1, n_probe=200):
    """(u, is_exact_crossing) for every candidate of H on [l0, l1]."""
    out = []
    delta = 1.0e-5*(l1 - l0)
    for c in adiabatic.find_resonance_candidates(H, l0, l1, n_probe=n_probe):
        lam = np.linalg.eigvalsh(np.asarray(H(c['l']), dtype=complex))
        u = c['gap']/(EPS*np.max(np.abs(lam))) if np.max(np.abs(lam)) > 0 else np.inf
        a, b = max(c['l'] - delta, l0), min(c['l'] + delta, l1)
        _, Va = np.linalg.eigh(np.asarray(H(a), dtype=complex))
        _, Vb = np.linalg.eigh(np.asarray(H(b), dtype=complex))
        swapped = abs(np.vdot(Va[:, c['j']], Vb[:, c['k']]))**2 > 0.5
        out.append((u, bool(swapped)))
    return out


def solar_calls():
    Es = np.geomspace(1.0, 300.0, 25)*gd.UNIT_MEV
    for f in (op.osc_prob_4nu_sun, op.osc_prob_4nu_sun_nsi, op.osc_prob_4nu_sun_liv,
              op.osc_prob_5nu_sun):
        for r in (1.0, 0.5, 0.3, 0.15):
            for nubar in (False, True):
                for E in Es:
                    yield ('%s r=%.2f' % (f.__name__, r),
                           lambda f=f, E=E, r=r, nubar=nubar: f(
                               energy=E, L=L_HALF, L0=0.0, nubar=nubar, strategy='hybrid',
                               ratio_number_neutrons_to_protons=r))
    for E in np.geomspace(60.0, 1000.0, 12)*gd.UNIT_MEV:
        yield ('3nu_sun s13=0', lambda E=E: op.osc_prob_3nu_sun(
            energy=E, L=L_HALF, L0=0.0, s13=0.0, strategy='hybrid'))
    for E in Es:
        for nubar in (False, True):
            kw = dict(energy=E, L=L_HALF, L0=0.0, nubar=nubar, strategy='hybrid')
            yield ('2nu_sun', lambda kw=kw: op.osc_prob_2nu_sun(**kw, sth=P0['s12'],
                                                                 Dm2=P0['D21']))
            yield ('3nu_sun', lambda kw=kw: op.osc_prob_3nu_sun(**kw))
            yield ('4nu_sun coupled', lambda kw=kw: op.osc_prob_4nu_sun(**kw, s14=0.1,
                                                                         D41=1.0))


def two_level(scale):
    """Two decoupled levels crossing at l = 0.37, written in units where |H| ~ scale."""
    def H(l):
        return scale*np.diag([1.0 + l, 2.0 - 2.0*l]).astype(complex)
    return H, 0.0, 1.0


def main():
    rows = {}
    for name, call in solar_calls():
        h = captured(call)
        if h is None:
            continue
        rows.setdefault(name, []).extend(classify(*h))
    for scale in (1.0e-20, 1.0e-10, 1.0, 1.0e10):
        rows.setdefault('two-level, |H| ~ %.0e' % scale, []).extend(classify(*two_level(scale)))

    A = [u for v in rows.values() for u, x in v if x]
    B = [u for v in rows.values() for u, x in v if not x]
    print('%-30s %8s %22s %8s %22s' % ('source', 'A: n', 'A: max u', 'B: n', 'B: min u'))
    for name, v in rows.items():
        a = [u for u, x in v if x]
        b = [u for u, x in v if not x]
        print('%-30s %8d %22s %8d %22s' % (
            name, len(a), '%.3g' % max(a) if a else '-', len(b), '%.3g' % min(b) if b else '-'))
    print('\nA  exact crossings    n = %6d   max u = %.3g' % (len(A), max(A)))
    print('B  everything else    n = %6d   min u = %.3g' % (len(B), min(B)))
    print('   B below 1e3: %d   below 1e6: %d   below 1e9: %d'
          % tuple(sum(u < t for u in B) for t in (1e3, 1e6, 1e9)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
