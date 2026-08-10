# -*- coding: utf-8 -*-
r"""gen_solar_nusquids.py

The nuSQuIDS half of notebook 25's solar section, pre-generated.

    python notebooks/gen_solar_nusquids.py > notebooks/external_solar_nusquids.json

**Averaged, not instantaneous.**  Over the ray from the core to the surface the accumulated
phase is ~13 000 radians, so a single nuSQuIDS evaluation samples a probability that turns over
thousands of times; neighbouring energies land anywhere between 0.15 and 0.9.  Comparing that
curve against Mag(nu)s's ``average=True`` would be comparing two different quantities.

**The shape of the calculation: one dense curve, then a window.**  P_ee is evaluated on a
single log-uniform energy grid spanning the whole range, in one ``EvolveState``; the averaged
observable is then a window mean over that curve, in post-processing.  nuSQuIDS evolves the
entire energy vector as one ODE system, so this costs about one evolution rather than one per
target, and the window width becomes a free parameter that can be changed without re-running
anything.  Both curves are stored: the raw one is what shows *why* the instantaneous
probability is not the observable, and the windowed one is what compares against Mag(nu)s.

WHY THE EARLIER VERSION RETURNED <P_ee> OF 2.55, 58.9 AND 3.19: **the solver tolerance was far
too loose for the accumulated phase**, and nothing in the output said so.  It swept 1e-2 to
1e-4, and measured on this grid the answer is not a probability until 1e-6:

    tol     seconds   P_ee range        |1 - sum over flavours|
    1e-4      165     0.0057 .. 2.8336         4.4e-16
    1e-5      260     0.0055 .. 0.9931         3.3e-16      (another flavour out of range)
    1e-6      587     0.0055 .. 0.9928         3.3e-16      physical
    1e-7     1085     0.0055 .. 0.9928         3.3e-16      physical

**The unitarity check cannot catch this, and it is worth understanding why.**  The obvious
diagnostic -- sum the flavour probabilities and assert they are 1 -- passes at 3e-16 on output
containing P_ee = 2.83, and on looser settings still it passes on values running from -19 to
+45.  nuSQuIDS evolves SU(3) coefficients in which the identity component *is* the trace, so
the flavour sum is conserved by construction however badly the traceless part is integrated.
The check that bites is each probability lying in [0, 1], which is what ``main`` asserts below.

Two other things were wrong and are fixed here, but neither was the cause, and both were
suspected before being tested: the twelve-cluster energy grid (a log-uniform grid is still the
right way to drive it, and cheaper, but the clustering did not break anything), and the track
below.  They are recorded because a plausible diagnosis that measurement rejects is worth as
much as the one it accepts.

THE TRACK WAS ALSO WRONG, and it was wrong in a way that still converged.  ``SunASnu.Track``
takes an **impact parameter**, not a production radius: measured, ``GetFinalX()`` is
``2 sqrt(R^2 - b^2)``, the full chord across the Sun at that impact parameter.  Passing
``0.05 R_sun`` as though it were a production radius therefore propagated a neutrino across
2 R_sun of Sun, entering and leaving through the surface, instead of outward from the
production point through 0.95 R_sun.  The two-argument form is ``Track(x_initial, b)``, so the
radial path a solar neutrino actually takes is ``b = 0`` starting at ``x = R + r``.  R is read
off the body itself rather than assumed, because nuSQuIDS's solar radius is its own.

Cost is dominated by the requested solver tolerance -- measured 10.9 s per energy at 1e-6,
1.39 s at 1e-4, 0.54 s at 1e-3 -- so the tolerance is swept here and becomes nuSQuIDS's dial
on the speed/accuracy panel.

Conventions are matched as in ``gen_profile_benchmarks.py``: nuSQuIDS reads the *same*
BS2005-AGS,OP file Mag(nu)s uses, via a temporary copy carrying all twelve columns (its reader
needs uniform columns, and it indexes them by position, so a trimmed file makes it spline the
wrong ones; its own default model path does not exist in this install).
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import json
import os
import platform
import sys
import tempfile
import time

import numpy as np

import magnus.globaldefs as gd

TARGETS = np.logspace(np.log10(0.1), np.log10(20.0), 40)*gd.UNIT_MEV
SPREAD = 0.05
# The grid runs a little past the targets at both ends, so that the +/-SPREAD window of the
# first and last target is covered by real nodes rather than by the edge of the grid.
GRID_PAD = 1.0 + 2.0*SPREAD
N_GRID = 200
# Measured on this grid: 165 s at 1e-4, 260 s at 1e-5, 587 s at 1e-6, 1085 s at 1e-7.  The
# first two return values OUTSIDE [0, 1] -- P_ee = 2.83 at 1e-4 -- so the two loose settings
# are kept deliberately, as the evidence for the cliff described below, and not because
# anything is read off them.
TOLERANCES = (1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7)
PRODUCTION_R = 0.05          # in units of the solar radius
TABLE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     '..', 'docs', 'dev', 'adversarial_batteries', 'bs05_agsop.dat')


def energy_grid():
    """One log-uniform grid over the whole range -- see the note on the clustered grid."""
    return np.logspace(np.log10(TARGETS[0]/GRID_PAD),
                       np.log10(TARGETS[-1]*GRID_PAD), N_GRID)


def solar_table_for_nusquids():
    """A temporary copy of the same file, in the form nuSQuIDS's reader accepts."""
    rows = []
    with open(TABLE) as fh:
        for line in fh:
            f = line.split()
            if len(f) == 12:
                try:
                    rows.append([float(x) for x in f])
                except ValueError:
                    continue
    handle = tempfile.NamedTemporaryFile('w', suffix='.dat', delete=False)
    np.savetxt(handle.name, np.array(rows), fmt='%.10e')
    handle.close()
    return handle.name


def probability_curve(tolerance, osc, table_path):
    """P_ee at every node of one log-uniform grid, from one evolution.

    Returns the raw curve and the wall time.  The averaging is deliberately NOT done
    here: a window mean over this curve is post-processing, so the window can be changed
    without paying for another evolution.
    """
    import nuSQuIDS as nsq
    units = nsq.Const()

    grid = energy_grid()                       # already sorted: logspace ascending
    solver = nsq.nuSQUIDS(grid*units.eV, 3, nsq.NeutrinoType.neutrino, False)
    body = nsq.SunASnu(table_path)
    solver.Set_Body(body)
    # b = 0 and start at R + r: the radial path outward from the production point.  The
    # single-argument form takes an IMPACT PARAMETER, which is a chord across the whole
    # Sun -- see the note at the top.  R is read off the body, not assumed.
    r_sun = nsq.SunASnu.Track(0.0, 0.0).GetFinalX()/2.0
    solver.Set_Track(nsq.SunASnu.Track(r_sun*(1.0 + PRODUCTION_R), 0.0))
    solver.Set_MixingAngle(0, 1, np.arcsin(osc['s12']))
    solver.Set_MixingAngle(0, 2, np.arcsin(osc['s13']))
    solver.Set_MixingAngle(1, 2, np.arcsin(osc['s23']))
    solver.Set_CPPhase(0, 2, osc['dCP'])
    solver.Set_SquareMassDifference(1, osc['D21'])
    solver.Set_SquareMassDifference(2, osc['D31'])
    solver.Set_rel_error(tolerance)
    solver.Set_abs_error(tolerance)
    state = np.zeros((len(grid), 3))
    state[:, 0] = 1.0
    solver.Set_initial_state(state, nsq.Basis.flavor)

    t0 = time.perf_counter()
    solver.EvolveState()
    elapsed = time.perf_counter() - t0

    p_all = np.array([[solver.EvalFlavorAtNode(f, i) for f in range(3)]
                      for i in range(len(grid))])
    # Both checks are kept, and the pair is the point: `unitarity` is the one the obvious
    # diagnostic would use and it passes on nonsense, while `physical` is the one that
    # actually detects it.  Storing both lets the notebook show that.
    unitarity = float(np.max(np.abs(p_all.sum(axis=1) - 1.0)))
    physical = bool(p_all.min() >= -1.0e-6 and p_all.max() <= 1.0 + 1.0e-6)
    return grid, p_all[:, 0], unitarity, physical, float(p_all.min()), \
        float(p_all.max()), elapsed


DENSE_TARGETS = np.array([1.0, 5.0, 15.0])*gd.UNIT_MEV
DENSE_SAMPLES = 300


def dense_check(osc, table_path, tolerance=1.0e-6):
    r"""Does the window mean actually converge to the adiabatic average, and at what cost?

    The sweep above spreads its nodes over the whole range, which leaves only three or
    four inside each +/-5% window.  A mean of four samples of a probability that turns
    over tens of thousands of times within the window is a Monte-Carlo estimate with a
    standard error near 0.14 -- so those <P_ee> values scatter across most of [0, 1] and
    say nothing.  That is not nuSQuIDS being wrong; it is what estimating an average from
    four samples costs.

    This puts DENSE_SAMPLES nodes inside the window at a few energies instead, so the
    estimate is good enough to compare against Mag(nu)s.  The standard error is reported
    with it, because an average quoted without one cannot be checked.
    """
    import nuSQuIDS as nsq
    units = nsq.Const()

    offsets = np.linspace(-SPREAD, SPREAD, DENSE_SAMPLES)
    grid = np.sort(np.concatenate([t*(1.0 + offsets) for t in DENSE_TARGETS]))
    solver = nsq.nuSQUIDS(grid*units.eV, 3, nsq.NeutrinoType.neutrino, False)
    solver.Set_Body(nsq.SunASnu(table_path))
    r_sun = nsq.SunASnu.Track(0.0, 0.0).GetFinalX()/2.0
    solver.Set_Track(nsq.SunASnu.Track(r_sun*(1.0 + PRODUCTION_R), 0.0))
    solver.Set_MixingAngle(0, 1, np.arcsin(osc['s12']))
    solver.Set_MixingAngle(0, 2, np.arcsin(osc['s13']))
    solver.Set_MixingAngle(1, 2, np.arcsin(osc['s23']))
    solver.Set_CPPhase(0, 2, osc['dCP'])
    solver.Set_SquareMassDifference(1, osc['D21'])
    solver.Set_SquareMassDifference(2, osc['D31'])
    solver.Set_rel_error(tolerance)
    solver.Set_abs_error(tolerance)
    state = np.zeros((len(grid), 3))
    state[:, 0] = 1.0
    solver.Set_initial_state(state, nsq.Basis.flavor)

    t0 = time.perf_counter()
    solver.EvolveState()
    elapsed = time.perf_counter() - t0

    p_all = np.array([[solver.EvalFlavorAtNode(f, i) for f in range(3)]
                      for i in range(len(grid))])
    rows = []
    for target in DENSE_TARGETS:
        sel = (grid >= target*(1.0 - SPREAD)) & (grid <= target*(1.0 + SPREAD))
        q = p_all[sel, 0]
        rows.append(dict(energy_ev=float(target), n_samples=int(sel.sum()),
                         mean=float(q.mean()), std=float(q.std(ddof=1)),
                         stderr=float(q.std(ddof=1)/np.sqrt(sel.sum()))))
        print('    %6.2f MeV: <P_ee> = %.4f +/- %.4f  (%d samples)'
              % (target/gd.UNIT_MEV, rows[-1]['mean'], rows[-1]['stderr'],
                 rows[-1]['n_samples']), file=sys.stderr)
    return dict(tolerance=tolerance, samples_per_target=DENSE_SAMPLES,
                seconds_total=elapsed,
                physical=bool(p_all.min() >= -1e-6 and p_all.max() <= 1.0 + 1e-6),
                points=rows)


def window_average(grid, p_ee):
    """<P_ee> at each target: the mean of the curve over +/-SPREAD around it.

    A window in *energy* is what a finite energy resolution means.  It is not the same
    device as Mag(nu)s's adiabatic average, which is the L/E -> infinity limit, and the
    two agree only because the window here is many oscillations wide at every target --
    the accumulated phase is thousands of radians, so +/-5% spans hundreds of turns.
    """
    out = np.empty(len(TARGETS))
    counts = np.empty(len(TARGETS), dtype=int)
    for k, target in enumerate(TARGETS):
        sel = (grid >= target*(1.0 - SPREAD)) & (grid <= target*(1.0 + SPREAD))
        counts[k] = int(sel.sum())
        out[k] = float(np.mean(p_ee[sel]))
    return out, counts


def main():
    osc = gd.load_nufit_params('NuFIT 6.1', 'NO')
    table_path = solar_table_for_nusquids()
    out = {'note': 'nuSQuIDS SunASnu on the same BS2005-AGS,OP file as Magnus; P_ee on '
                   'one log-uniform grid of %d nodes, then a +/-%.0f%% window mean'
                   % (N_GRID, 100*SPREAD),
           'machine': platform.platform(),
           'osc_params': 'NuFIT 6.1 NO',
           'production_radius_over_rsun': PRODUCTION_R,
           'n_grid': N_GRID,
           'spread': SPREAD,
           'energy_ev': [float(e) for e in TARGETS],
           'series': []}
    for tol in TOLERANCES:
        grid, p_ee, unitarity, physical, p_lo, p_hi, elapsed = probability_curve(
            tol, osc, table_path)
        avg, counts = window_average(grid, p_ee)
        out['series'].append({'tolerance': tol,
                              'seconds_total': elapsed,
                              'seconds_per_target': elapsed/len(TARGETS),
                              'unitarity': unitarity,
                              'physical': physical,
                              'p_min': p_lo, 'p_max': p_hi,
                              'nodes_per_window': [int(c) for c in counts],
                              'grid_energy_ev': [float(e) for e in grid],
                              'P_ee_instantaneous': [float(x) for x in p_ee],
                              'P_ee': [float(x) for x in avg]})
        print('  tol %.0e: %.1f s, %d nodes, |1-sum| = %.2e, P in %.4f..%.4f  %s'
              % (tol, elapsed, len(grid), unitarity, p_lo, p_hi,
                 'physical' if physical else 'NOT A PROBABILITY'), file=sys.stderr)

    # The loose settings are recorded on purpose -- they are the evidence for the cliff.
    # What must not happen is the file being written when even the TIGHTEST setting is
    # unphysical, because then there is no curve in it worth plotting and the last
    # version of this script shipped exactly that.
    if not out['series'][-1]['physical']:
        raise SystemExit('tightest tolerance %.0e is still unphysical (P in %.4f..%.4f) '
                         '-- not writing' % (TOLERANCES[-1],
                                             out['series'][-1]['p_min'],
                                             out['series'][-1]['p_max']))
    json.dump(out, sys.stdout, indent=1)


if __name__ == '__main__':
    main()
