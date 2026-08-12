#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
r"""Independent re-derivation of the d = 3 ``SEV_TOL`` gate.

``expmkernels.SEV_TOL`` is the conditioning above which a 3x3 is handed back to
``eigh`` instead of being solved in closed form.  It gates on
:math:`m = \mathrm{tr}(X^2)/6` of the traceless part of :math:`K`, and its
docstring states the calibration behind the value 1e4:

    no cell at spectral scale <= 1e2 (m <= 1.1e3) is worse than ``eigh`` by more
    than 2e-13 absolute, while cells at scale 1e3 (m >= 1.1e5) reach 131x and at
    1e5 reach 7440x.

``docs/dev/HANDOVER_OVERHEAD.md`` records that this calibration is one of two
claims **neither** max-effort review checked -- both were taken on the word of
the commit that introduced them.  This script is that check.

What it does
------------

Sweeps the two axes that matter *together*, which is the whole point: the damage
needs a clustered spectrum **and** a large norm, and neither axis alone reaches
it.  Sweeping norms at generic separation, or separations at norm 1, is exactly
how the original unqualified claim came to be believed.

For each (spectral scale, separation) cell it builds a Hermitian K with that
spectrum in a random unitary basis, computes ``exp(-iK)`` three ways -- the
closed-form kernel, ``eigh``, and ``scipy.linalg.expm`` as an outside referee --
and reports the closed form's error against the referee as a ratio to ``eigh``'s.

Then it reports the smallest m at which the closed form is meaningfully worse,
and compares that to the gate.

Usage
-----

::

    python3 docs/dev/calibrate_sev_tol.py                 # the standard grid
    python3 docs/dev/calibrate_sev_tol.py --bases 16      # more random bases

It prints a table and a verdict.  It asserts nothing: this is a measurement, and
what to do about the number is a judgement about how much margin the gate should
carry.

UNRESOLVED, 2026-08-12
----------------------

This grid disagrees with the calibration the gate is set from, and the
disagreement is **not** settled:

* here, cells the 1e4 gate ADMITS reach 5.1e-13 absolute against ``eigh``'s
  6e-15 -- past the 2e-13 budget the SEV_TOL docstring claims for admitted
  cells -- with the smallest offender at m = 4.4e3.
* ``HANDOVER_OVERHEAD.md`` records the original calibration finding the first
  unsafe cell at 1.1e5, an 11x margin, and flags it as checked by neither
  max-effort review.
* ``tests/test_expm_backend.py::test_sev_tol_sits_inside_its_calibrated_window``
  encodes the original: it asserts m(scale 1e2) < SEV_TOL < m(scale 1e3), so
  lowering the gate to 1e3 **fails that test** -- 1e3 would decline cells the
  original measured as safe.

So the two do not merely differ in margin, they differ about whether a
scale-1e2 cell is safe at all.  Lowering the gate was considered and reverted
for exactly that reason: it is not a conservative tweak if it contradicts a
standing calibration, it is picking a side.

What would settle it, in rough order of cost: compare the spectrum families --
this script uses ``[-s, -s(1-d), s]`` while the test uses ``_GRID_SHAPES
['double-low']``, and they need not probe the same corner -- then identify which
cell here produced the 5.1e-13 and whether its m really sits below 1e4, then
decide whether the 2e-13 budget in the docstring is the right claim to hold the
gate to.

Worth knowing before spending that time: instrumenting the kernel across a PREM
chord, a constant-density call, a 60-energy Earth scan and a solar profile, the
severity actually reached is m <~ 10.  A Magnus slab has ||Omega|| <~ pi by
construction, so the ladder cannot produce a badly conditioned exponential and
neither candidate gate ever fires in ordinary use.  This is a question about
whether a documented guarantee is true, not about numbers users are getting.
"""

import argparse
import os
import sys

import numpy as np
import scipy.linalg as sla

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, os.pardir, 'src'))

import magnus.expmkernels as expmkernels             # noqa: E402


def traceless_m(K):
    """The quantity the gate is keyed on: m = tr(X^2)/6, X the traceless part."""
    X = K - np.eye(K.shape[-1])*np.trace(K)/K.shape[-1]
    return float(np.trace(X @ X).real/6.0)


def hermitian_with_spectrum(eigenvalues, rng):
    """A Hermitian matrix with the given spectrum, in a random unitary basis.

    The basis matters.  The closed form is not basis-invariant in finite
    precision, so a calibration done in one basis says nothing about another --
    which is why the handover re-checked in eight.
    """
    A = rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3))
    Q, _ = np.linalg.qr(A)
    return Q @ np.diag(eigenvalues).astype(complex) @ Q.conj().T


def measure(scale, separation, rng):
    """One cell: closed-form error against expm, as a ratio to eigh's."""
    # Two eigenvalues separated by `separation` relative to the scale, and a
    # third far from both -- clustering is a property of a *pair*.
    lam = np.array([-scale, -scale*(1.0 - separation), scale])
    K = hermitian_with_spectrum(lam, rng)

    ref = sla.expm(-1j*K)                       # outside referee

    w, V = np.linalg.eigh(K)
    via_eigh = (V*np.exp(-1j*w)) @ V.conj().T

    closed, _, _ = expmkernels.expm_herm_stack(K[None, ...].copy())
    closed = closed[0]

    err_closed = float(np.max(np.abs(closed - ref)))
    err_eigh = float(np.max(np.abs(via_eigh - ref)))
    ratio = err_closed/err_eigh if err_eigh > 0 else np.inf
    return traceless_m(K), err_closed, err_eigh, ratio


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--bases', type=int, default=8,
                    help='random unitary bases per cell (default 8)')
    ap.add_argument('--seed', type=int, default=20260811)
    ap.add_argument('--unsafe-ratio', type=float, default=10.0,
                    help='ratio to eigh above which a cell is flagged (secondary)')
    ap.add_argument('--abs-budget', type=float, default=2.0e-13,
                    help='absolute error the gate was calibrated against (primary)')
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    scales = [1.0e0, 1.0e1, 1.0e2, 1.0e3, 1.0e4, 1.0e5]
    separations = [1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0e-1, 5.0e-1]

    print('SEV_TOL = %.3g   (gate on m = tr(X^2)/6)' % expmkernels.SEV_TOL)
    print('numba kernels: %s' % expmkernels.HAVE_NUMBA)
    print('%d random bases per cell, worst taken\n' % args.bases)
    print('%-10s %-10s %-12s %-12s %-12s %s'
          % ('scale', 'separation', 'm', 'err closed', 'err eigh', 'ratio'))

    unsafe_m = []          # by ratio to eigh -- the secondary reading
    over_budget_m = []     # by absolute error -- what the gate was calibrated on
    safe_max_ratio = 0.0
    safe_max_abs = 0.0
    for scale in scales:
        for sep in separations:
            worst = None
            for _ in range(args.bases):
                m, ec, ee, ratio = measure(scale, sep, rng)
                if worst is None or ratio > worst[3]:
                    worst = (m, ec, ee, ratio)
            m, ec, ee, ratio = worst
            flag = ''
            if ratio > args.unsafe_ratio:
                unsafe_m.append(m)
                flag = '  <-- ratio'
            if ec > args.abs_budget:
                over_budget_m.append(m)
                flag = ('  <-- OVER BUDGET, BELOW THE GATE' if m < expmkernels.SEV_TOL
                        else '  <-- over budget')
            if m < expmkernels.SEV_TOL:
                safe_max_ratio = max(safe_max_ratio, ratio)
                safe_max_abs = max(safe_max_abs, ec)
            print('%-10.0e %-10.0e %-12.4g %-12.3e %-12.3e %8.1fx%s'
                  % (scale, sep, m, ec, ee, ratio, flag))

    print()
    print('--- the criterion the gate was calibrated against: absolute error ---')
    print('The docstring of SEV_TOL states a budget, not a ratio: "no cell at spectral')
    print('scale <= 1e2 is worse than eigh by more than 2e-13 absolute".  A ratio reading')
    print('is the wrong test and will condemn a gate that is doing its job -- a cell whose')
    print('errors are 9e-14 against 6e-15 is 15x worse and still far inside the budget.')
    print()
    print('worst absolute error among cells the gate ADMITS (m < %.3g): %.3e'
          % (expmkernels.SEV_TOL, safe_max_abs))
    print('budget                                                     : %.3e'
          % args.abs_budget)
    if safe_max_abs <= args.abs_budget:
        print('VERDICT: the gate HOLDS.  Everything it admits is inside the budget,')
        print('         by a factor of %.1f.' % (args.abs_budget/max(safe_max_abs, 1e-300)))
    else:
        print('VERDICT: the gate does NOT hold.  A cell below SEV_TOL exceeds the')
        print('         absolute budget it was calibrated to respect.')
    if over_budget_m:
        lo = min(over_budget_m)
        print('smallest m over budget : %.4g   (gate %.3g, margin %.2fx)'
              % (lo, expmkernels.SEV_TOL, lo/expmkernels.SEV_TOL))
        print('the handover records 1.100e5, an 11x margin, from the original')
        print('calibration -- compare with the number above.')

    print()
    print('--- secondary: ratio to eigh ---')
    print('Useful for seeing where the closed form starts to degrade at all, but it')
    print('amplifies as both errors approach machine epsilon, so it is not a verdict.')
    if unsafe_m:
        lo = min(unsafe_m)
        print('smallest m with ratio > %.0fx : %.4g  (%s the gate)'
              % (args.unsafe_ratio, lo, 'below' if lo < expmkernels.SEV_TOL else 'above'))
        print('worst ratio among admitted cells: %.1fx' % safe_max_ratio)
    else:
        print('no cell reached %.0fx of eigh on this grid; the grid may not span the'
              % args.unsafe_ratio)
        print('corner the gate exists for.  Widen the scales before concluding anything.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
