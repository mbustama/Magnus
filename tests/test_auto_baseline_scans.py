# -*- coding: utf-8 -*-
"""Baseline scans under ``strategy='auto'`` at tight tolerances (issue #125).

Below ``AUTO_LADDER_MIN_TOLERANCE`` a single-energy scan of two to seven baselines used to go to
the hybrid strategy, since the threshold at which it yields to the cumulative scan,
``HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS = 8``, was set at loose tolerances.  There the cumulative
scan took 0.06 to 0.18 of the time at the median, with no silent miss the hybrid did not also
make (``docs/dev/measurements/issue125_baseline_scans/``), so below 1e-6 the threshold is
``HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS_TIGHT = 2``.  At 1e-6 and looser nothing changes.
"""

import warnings

import numpy as np
from scipy.integrate import solve_ivp

import magnus.globaldefs as gd
import magnus.hamiltonians as hm
import magnus.matter as mt
import magnus.oscprob as op

KM = gd.UNIT_KM
OSC = {k: gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT'][k]
       for k in ('s12', 's23', 's13', 'dCP', 'D21', 'D31')}
E = 0.02*gd.UNIT_GEV
LS = np.array([800.0, 1600.0, 2500.0])*KM


def rho(l):
    return 3.0 + 2.0*np.exp(-np.asarray(l, dtype=float)/(1000.0*KM))


def engine(L, tol, **kw):
    info = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        P = op.osc_prob_matter_std_potential(3, rho, E, L, OSC, L0=0.0, rtol=tol, atol=tol,
                                             density_matter_is_in_g_per_cm3=True,
                                             strategy_info=info, **kw)
    return np.asarray(P), info.get('engine')


def dop853(L):
    hv = np.asarray(hm.hamiltonian_3nu_vacuum_energy_independent(**OSC))/E
    per = mt.VCC_func(0.0, lambda l: 1.0)*gd.UNIT_G_PER_CM3/gd.ATOMIC_MASS_UNIT*0.5

    def rhs(x, y):
        H = hv + np.diag([per*rho(x), 0.0, 0.0])
        return (-1j*H @ y.reshape(3, 3)).ravel()
    s = solve_ivp(rhs, (0.0, float(L[-1])), np.eye(3, dtype=complex).ravel(), method='DOP853',
                  rtol=1e-12, atol=1e-13, t_eval=L)
    return np.array([(np.abs(s.y[:, i].reshape(3, 3))**2).T for i in range(len(L))])


def test_a_short_baseline_scan_at_a_tight_tolerance_takes_the_cumulative_scan():
    P, eng = engine(LS, 1e-9)
    assert eng == 'cumulative'                       # 'hybrid' before issue #125
    assert np.max(np.abs(P - dop853(LS))) < 1e-8


def test_the_loose_tolerance_and_the_other_routes_are_unchanged():
    assert engine(LS, 1e-3)[1] == 'cumulative'       # already so at 1e-6 and looser
    assert engine(LS[-1], 1e-9)[1] != 'cumulative'   # a single point is not a scan
    assert engine(LS, 1e-9, strategy='hybrid')[1] == 'hybrid'
    assert engine(LS, 1e-9, cumulative=False)[1] != 'cumulative'


def test_a_coherent_solar_scan_at_the_default_tolerance_keeps_the_hybrid():
    """The threshold of 8 still applies at 1e-6 and looser, where it was measured."""
    osc = gd.load_nufit_params('NuFIT 6.1')
    R = gd.SUN_RADIUS*KM
    for tol, expected in ((1e-3, 'hybrid'), (1e-9, 'cumulative')):
        info = {}
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            op.osc_prob_3nu_sun(0.01*gd.UNIT_GEV, np.array([0.05, 0.15, 0.3])*R, 0.0, **osc,
                                nu_i=gd.NUE, nu_f=gd.NUE, rtol=tol, atol=tol, strategy_info=info)
        assert info.get('engine') == expected
