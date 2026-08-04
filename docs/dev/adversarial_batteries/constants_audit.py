# -*- coding: utf-8 -*-
"""Provenance for the unaudited calibration constants (robustness programme, item 2).

The brief's finding: ``GAMMA_TO_ERROR`` was set from five configurations that all sat in one
corner of the parameter space, and of the four constants measured properly in the previous
session, **three were wrong**.  The rest have the same provenance -- "it has always been that"
-- and this measures the two the brief ranks first.

Method, and it is the part that matters: for each constant, identify **the regime in which the
constant is consulted**, sample that regime, and report a distribution rather than one number.
The ``GAMMA_TO_ERROR`` mistake was reading a maximum over a population that included rows the
constant never decides.

Sub-tests:

  1  ``fd_step_frac = 1e-6`` against the step-size trade-off.  A central difference has
     truncation error ~ h^2 and cancellation error ~ eps/h, so there is an optimum and nothing
     has ever checked that 1e-6 is near it.  Scored against the **analytic** dH/dl, available
     here because H(l) = h_vac/E + C*ne(l)*proj and every profile below has a closed-form
     ne'(l) -- no reference differentiation, so the reference carries no error of its own.
  2  the same constant end to end: does moving it change the answer, or only the derivative?
  3  ``threshold0 = 0.1``, as a cost-against-accuracy trade, swept against the requested
     tolerance to test the brief's hypothesis that the right value is a rule, not a constant.

Run:  python constants_audit.py [1|2|3 ...]
"""

import sys
import time
import warnings

import numpy as np

import harness as H
import magnus.adiabatic as ad
from battery2 import _C_VCC, ne_res_for

L0, L1 = 0.0, 1.0*H.L_SCALE
SPAN = L1 - L0


# ------------------------------------------------------------------ profiles with analytic ne'
def exp_profile():
    ls = H.L_SCALE

    def ne(l):
        return H.scalarize(H.NE0*np.exp(-np.asarray(l, float)/ls))

    def dne(l):
        return -H.NE0*np.exp(-np.asarray(l, float)/ls)/ls
    return 'solar exponential', ne, dne


def sine_profile(period, base):
    def ne(l):
        return H.scalarize(base*(1.0 + 0.9*np.sin(2.0*np.pi*np.asarray(l, float)/period)))

    def dne(l):
        return base*0.9*(2.0*np.pi/period)*np.cos(2.0*np.pi*np.asarray(l, float)/period)
    return 'sinusoid, period span/7', ne, dne


def gauss_profile(ner, lc, w):
    a, b = 0.30*ner, (3.0 - 0.30)*ner

    def ne(l):
        x = np.asarray(l, float)
        return H.scalarize(a + b*np.exp(-0.5*((x - lc)/w)**2))

    def dne(l):
        x = np.asarray(l, float)
        return -b*(x - lc)/w**2*np.exp(-0.5*((x - lc)/w)**2)
    return 'gaussian bump w=1e-2 span', ne, dne


def analytic_profiles():
    p2 = H.params_for(2)
    ner2 = ne_res_for(2, p2, 10.0e6)
    return [exp_profile(),
            sine_profile(SPAN/7.0, 3.0e-2*H.NE0),
            gauss_profile(ner2, 0.45*SPAN, 1e-2*SPAN)]


# ------------------------------------------------------------------ 1: the step-size trade-off
def sub1():
    print('## 2.1  fd_step_frac against the central-difference step-size trade-off')
    print('    reference: the ANALYTIC dH/dl (no reference differentiation)')
    fracs = [10.0**k for k in range(-12, -1)]
    print('%-28s %3s %10s %s' % ('profile', 'd', 'E [MeV]',
                                 ''.join('%10s' % ('1e%d' % k) for k in range(-12, -1))))
    best_overall = {}
    for label, ne, dne in analytic_profiles():
        for d in (2, 3):
            params = H.params_for(d)
            for energy in (10.0e6, 50.0e6):
                H_of_l = H.H_factory(d, params, H.vcc_of(ne), energy)
                proj = np.zeros((d, d), dtype=complex)
                proj[0, 0] = 1.0
                # Sample the path where the diagnostics are actually evaluated: the probe grid.
                ls = np.linspace(L0, L1, 200)[1:-1][::7]
                # Normalised by the LARGEST |dH/dl| anywhere on the path, not per position.
                # Per-position normalisation divides ~0 by ~0 wherever the profile is flat --
                # which on a narrow Gaussian is almost everywhere -- and reported a relative
                # error of exactly 1.0 at every step size, i.e. the instrument saturated before
                # the code did.  A single global scale is what "how accurate is this
                # derivative" means for a diagnostic that compares gamma across the path.
                exact = np.array([_C_VCC*dne(float(l))*proj for l in ls])
                ref = max(float(np.max(np.abs(exact))), 1e-300)
                errs = []
                for frac in fracs:
                    h = SPAN*frac
                    got = np.array([ad._dH_dl(H_of_l, float(l), h, (L0, L1)) for l in ls])
                    errs.append(float(np.max(np.abs(got - exact))/ref))
                k_best = int(np.argmin(errs))
                best_overall.setdefault(label, []).append(fracs[k_best])
                print('%-28s %3d %10.0f %s   best 1e%d'
                      % (label, d, energy/1e6,
                         ''.join('%10.1e' % e for e in errs), int(np.log10(fracs[k_best]))))
    print('\n    optimum by profile: %s'
          % ', '.join('%s -> %s' % (k, sorted({'1e%d' % int(np.log10(f)) for f in v}))
                      for k, v in best_overall.items()))
    print()


# ------------------------------------------------------------------ 2: does it reach the answer?
def sub2():
    print('## 2.2  fd_step_frac end to end: hybrid_propagator error against solve_ivp')
    fracs = [1e-10, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]
    print('%-28s %3s %9s %s' % ('profile', 'd', 'E [MeV]',
                                ''.join('%11s' % ('%.0e' % f) for f in fracs)))
    for label, ne, _ in analytic_profiles():
        for d in (2, 3):
            params = H.params_for(d)
            for energy in (10.0e6, 50.0e6):
                H_of_l = H.H_factory(d, params, H.vcc_of(ne), energy)
                Pref = H.P_of(H.exact_U(H_of_l, L0, L1, d))
                row = []
                for frac in fracs:
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore')
                        U, win, cert = ad.hybrid_propagator(H_of_l, L0, L1, rtol=1e-3, atol=1e-3,
                                                            fd_step_frac=frac)
                    row.append('%9.2e%s' % (H.maxabs(H.P_of(U) - Pref), '' if cert else '*'))
                print('%-28s %3d %9.0f %s' % (label, d, energy/1e6, ' '.join(row)))
    print('    (* = uncertified)\n')


# ------------------------------------------------------------------ 3: threshold0
def sub3():
    print('## 2.3  threshold0: cost against accuracy, at three requested tolerances')
    print('    the regime it governs: it is the FIRST threshold tried, so it decides how many')
    print('    refinement iterations are needed before a window opens -- not, since the gamma')
    print('    rule was added, whether the answer may be certified.')
    thresholds = [1.0, 0.3, 0.1, 0.03, 0.01, 1e-3]
    for tol in (1e-2, 1e-3, 1e-5):
        print('\n  requested rtol = atol = %.0e' % tol)
        print('  %-28s %3s %9s %s' % ('profile', 'd', 'E [MeV]',
                                      ''.join('%20s' % ('t0=%.3g' % t) for t in thresholds)))
        for label, ne, _ in analytic_profiles():
            for d in (2, 3):
                params = H.params_for(d)
                for energy in (10.0e6,):
                    H_of_l = H.H_factory(d, params, H.vcc_of(ne), energy)
                    Pref = H.P_of(H.exact_U(H_of_l, L0, L1, d))
                    row = []
                    for t0 in thresholds:
                        t = time.time()
                        with warnings.catch_warnings():
                            warnings.simplefilter('ignore')
                            U, win, cert = ad.hybrid_propagator(
                                H_of_l, L0, L1, rtol=tol, atol=tol, threshold0=t0)
                        dt = time.time() - t
                        row.append('%8.2e%s %4dw%5.2fs'
                                   % (H.maxabs(H.P_of(U) - Pref), '' if cert else '*',
                                      len(win), dt))
                    print('  %-28s %3d %9.0f %s'
                          % (label, d, energy/1e6, ' '.join(row)), flush=True)
    print('\n    (* = uncertified;  Nw = windows opened;  seconds are RELATIVE only -- the box')
    print('     is loaded, so read the ratios across a row, never the absolute times)\n')


if __name__ == '__main__':
    for w in (sys.argv[1:] or ['1', '2', '3']):
        {'1': sub1, '2': sub2, '3': sub3}[w]()
