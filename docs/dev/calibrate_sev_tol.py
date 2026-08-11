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
                    help='ratio to eigh above which a cell counts as unsafe')
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    scales = [1.0e0, 1.0e1, 1.0e2, 1.0e3, 1.0e4, 1.0e5]
    separations = [1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2, 1.0e-1, 5.0e-1]

    print('SEV_TOL = %.3g   (gate on m = tr(X^2)/6)' % expmkernels.SEV_TOL)
    print('numba kernels: %s' % expmkernels.HAVE_NUMBA)
    print('%d random bases per cell, worst taken\n' % args.bases)
    print('%-10s %-10s %-12s %-12s %-12s %s'
          % ('scale', 'separation', 'm', 'err closed', 'err eigh', 'ratio'))

    unsafe_m = []
    safe_max_ratio = 0.0
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
                flag = '  <-- unsafe'
                if m < expmkernels.SEV_TOL:
                    flag = '  <-- UNSAFE AND BELOW THE GATE'
            elif m < expmkernels.SEV_TOL:
                safe_max_ratio = max(safe_max_ratio, ratio)
            print('%-10.0e %-10.0e %-12.4g %-12.3e %-12.3e %8.1fx%s'
                  % (scale, sep, m, ec, ee, ratio, flag))

    print()
    if unsafe_m:
        lo = min(unsafe_m)
        print('smallest m with ratio > %.0fx : %.4g' % (args.unsafe_ratio, lo))
        print('gate                        : %.4g' % expmkernels.SEV_TOL)
        print('margin                      : %.2fx' % (lo/expmkernels.SEV_TOL))
        print()
        if lo <= expmkernels.SEV_TOL:
            print('VERDICT: the gate does NOT keep the damage out of reach -- a cell')
            print('         below SEV_TOL is worse than eigh by more than %.0fx.'
                  % args.unsafe_ratio)
        else:
            print('VERDICT: gate holds.  Everything it admits is within %.1fx of eigh;'
                  % safe_max_ratio)
            print('         the first unsafe cell sits %.2fx above it.'
                  % (lo/expmkernels.SEV_TOL))
            print('         The handover records 1.100e5, an 11x margin, from the')
            print('         original calibration -- compare with the number above.')
    else:
        print('VERDICT: no cell reached %.0fx of eigh anywhere on this grid, so the')
        print('         grid does not span the corner the gate exists for.  Widen the')
        print('         scales before concluding the gate is unnecessary.'
              % ())
    return 0


if __name__ == '__main__':
    sys.exit(main())
