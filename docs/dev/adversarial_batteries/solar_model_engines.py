# -*- coding: utf-8 -*-
r"""Which engine answers a Sun wrapper on each tabulated solar model, and what it warns.

P_ee from the centre to 0.9 R_sun, through osc_prob_3nu_sun, for the exponential fit and all
twelve models at 1, 5, 10 and 20 MeV, averaged and coherent.  Measured on 2026-09-23, when the
solar models were added:

  averaged   every call on the phase-average engine, none escalated, no warning;
  coherent   44 of 48 certified on the hybrid engine (three with MagnusConvergenceWarning);
             BP04 and BS05-OP at 10 and 20 MeV fell to the slab ladder and raised
             ToleranceNotAchievedWarning.  solar_model_coherent.py checks those answers
             against a converged reference.

Run from the repository root:  python docs/dev/adversarial_batteries/solar_model_engines.py
"""
import time
import warnings

import numpy as np

import magnus.globaldefs as gd
import magnus.oscprob as op
from magnus import solarmodels as sm

R = gd.SUN_RADIUS*gd.UNIT_KM
models = ['exp'] + list(sm.SOLAR_MODELS)
for avg in (True, False):
    for m in models:
        for e in (1.0, 5.0, 10.0, 20.0):
            info = {}
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter('always')
                t = time.time()
                P = float(np.asarray(op.osc_prob_3nu_sun(energy=e*gd.UNIT_MEV, L=0.9*R, L0=0.0, nu_i=0, nu_f=0,
                                                         density_profile=m, average=avg, strategy_info=info)))
                dt = time.time() - t
            ws = sorted({x.category.__name__.replace('Warning', '') for x in w})
            extra = {k: info[k] for k in ('escalated', 'resolved', 'certified') if k in info}
            print('%-4s %-14s %5.1f MeV  P=%.8f  %-8s %-18s cert=%-5s %5.1fs  %s %s' % (
                'avg' if avg else 'coh', m, e, P, info.get('engine'), info.get('family'), info.get('certified'),
                dt, ws or '', extra or ''), flush=True)
