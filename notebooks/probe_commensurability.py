# -*- coding: utf-8 -*-
r"""probe_commensurability.py

Is a timing taken today comparable with the ones frozen in
``external_profile_benchmarks.json`` on 2026-08-10?

Figure 11 puts every series on one time axis and its caption claims they were timed in
one process on one machine.  Appending new Magnus curves measured in a *later* session
only keeps that claim honest if the machine still times the same way.  This probe asks
that question directly, and cheaply, before any sweep is run:

  1. The interleaved control workload that produced ``control_ratio`` in the stored file.
     Same nine-fold interleave, same 180x180 matmul.
  2. One stored Magnus point, re-timed: order 4 at the tightest tolerance available for
     three flavours, scored against the *stored* reference so the comparison is of speed
     alone and not of accuracy.

Nothing is written.  Run it on an idle machine::

    python notebooks/probe_commensurability.py

A re-timed point within a few percent of the stored one means the sessions are
commensurable and the new order-6 and order-8 series can be appended without disturbing
anything already in the file.  A large drift means the honest route is a full re-run,
which is a decision for the author rather than for this script.
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

import json
import pathlib
import platform
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import gen_profile_benchmarks as gpb                      # noqa: E402

HERE = pathlib.Path(__file__).parent
STORE = HERE/'external_profile_benchmarks.json'


def main():
    stored = json.loads(STORE.read_text())
    print('stored machine : %s' % stored['machine'])
    print('this machine   : %s' % platform.platform())
    print()

    # --- 1. the control, interleaved exactly as the generator does it ---------------
    best = {'a': np.inf, 'b': np.inf}
    for _ in range(9):
        for k in ('a', 'b'):
            best[k] = min(best[k], gpb.timed(gpb.control, repeat=1))
    ratio = best['a']/best['b']
    print('control ratio   stored %.3f   now %.3f   (drift %+.1f%%)'
          % (stored['control_ratio'], ratio,
             100.0*(ratio - stored['control_ratio'])/stored['control_ratio']))
    print('control block   %.3f ms  (absolute, for reference)' % (1.0e3*best['a']))
    print()

    # --- 2. one stored Magnus point, re-timed --------------------------------------
    case = next(c for c in stored['cases'] if c['flavours'] == 3)
    mg = next(s for s in case['series'] if s['name'] == 'Magnus')
    pt = min(mg['points'], key=lambda p: p['rtol'])       # the tightest, so the slowest
    prof = gpb.exponential_profile()
    ref = np.array(case['reference'])

    # magnus_points() recomputes the DOP853 referee, which the stored file already
    # carries.  Time the same call it would, and score against the stored reference.
    per_ne = gpb.matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)

    def ne_of(x):
        return prof['vcc'](x)/per_ne

    E = prof['energies']
    osc = gpb.osc_params(3)

    def call(r=pt['rtol']):
        return np.asarray(gpb.oscprob.osc_prob_matter_std_potential(
            3, ne_of, E, prof['baseline'], osc, L0=0.0,
            nu_i=gpb.gd.NUMU, nu_f=gpb.gd.NUMU,
            density_is_of_number_of_electrons=True, rtol=r, atol=r*1.0e-2,
            strategy='magnus'))

    P = call()
    us = 1.0e6*gpb.timed(call)/len(E)
    err = float(np.max(np.abs(P - ref)))
    d_t = 100.0*(us - pt['us_per_probability'])/pt['us_per_probability']
    print('Magnus order 4, 3nu, rtol = %.0e' % pt['rtol'])
    print('  us/probability  stored %10.2f   now %10.2f   (drift %+.1f%%)'
          % (pt['us_per_probability'], us, d_t))
    print('  max |dP|        stored %10.3e   now %10.3e   (must match: same code, '
          'same reference)' % (pt['max_abs_error'], err))
    print()
    verdict = 'COMMENSURABLE' if abs(d_t) < 10.0 else 'DRIFTED -- ask the author'
    print('verdict: %s' % verdict)


if __name__ == '__main__':
    main()
