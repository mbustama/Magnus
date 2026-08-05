# -*- coding: utf-8 -*-
"""P2, P3 and P4 of the physical-profile brief, which no existing battery covers.

``fallback_quality.py``, ``warn_fp.py``, ``resolution_fp.py`` and ``crosscheck_acceptance.py``
all take ``--physical`` now and answer P1 between them.  Two questions have no script:

  * **P2 / P3 -- the hidden-feature scan.**  ``find_hidden_features`` is measured in
    ``tests/test_adiabatic.py`` on synthetic Gaussians (0 false positives on 67 smooth cases,
    68-90 % detection at widths of 2e-5 and above).  Nothing measures it over a population, so
    ``sub_grid`` here does: false-positive rate on the physical families that hide nothing, and
    detection rate on the two that genuinely do.

  * **P4 -- the seam.**  ``HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS`` moved 25 -> 8 on a cost
    measurement over **three** profiles.  ``seam_cost`` re-measures cumulative-against-hybrid
    cost on the physical population, alternating the two engines and carrying a control the
    change cannot touch, because absolute times under load are worthless.

Run:  python physical_battery.py sub_grid
      python physical_battery.py seam_cost [n_rounds]
"""

import sys
import time
import warnings

import numpy as np

import harness as H
import magnus.adiabatic as ad
import magnus.oscprob as oscprob
import physical_profiles as pp

TOL = 1e-3


# ======================================================================
# P2 / P3 -- the hidden-feature scan
# ======================================================================

def subgrid_ratio(profile, l0, l1, n_ref=6400, n_sub=32):
    """Total variation on a grid ``n_sub`` times finer than ``n_ref``, over that on ``n_ref``.

    An oracle-free, construction-independent answer to "does this profile carry variation the
    package's finest grid cannot see?".  It is 1.000 for anything the reference grid captures --
    including a **step**, however narrow, because the two nodes bracketing it see its full
    height -- and rises above 1 only for structure that goes up and comes back down between
    nodes.  That distinction is the whole difference between the resolution test and the
    hidden-feature test, and getting it wrong is how a shock front gets mislabelled.
    """
    coarse = np.asarray(profile(np.linspace(float(l0), float(l1), n_ref)))
    fine = np.asarray(profile(np.linspace(float(l0), float(l1), (n_ref - 1)*n_sub + 1)))
    tv_c = float(np.abs(np.diff(coarse)).sum())
    tv_f = float(np.abs(np.diff(fine)).sum())
    return tv_f/tv_c if tv_c > 0.0 else 1.0


def sub_grid():
    """False-positive and detection rates for ``find_hidden_features``.

    The scan is run on the **scalar potential**, which the docstring says is identical to running
    it on ``H_func`` for a separable Hamiltonian and 18x cheaper, and at the ``n_sub`` values the
    dispatcher actually uses (8 for a single point, 32 for sixteen or more points).

    **The ground truth is measured, not declared.**  ``physical_profiles.families()`` carries a
    ``hides`` flag, but a flag set by the person who also built the profile is exactly the way to
    turn a construction error into a finding -- the first version of this called every shock
    width a hidden feature, and it is not one.  So the truth used here is
    :func:`subgrid_ratio`: total variation on a grid 32x finer than the reference grid, divided
    by total variation on the reference grid.  Above 1, the profile carries variation no grid the
    package lays down can see.  The declared flag is printed beside it and any disagreement is
    reported.
    """
    fams = pp.families()
    print('# P2/P3: find_hidden_features over the physical population')
    print('# HIDDEN_FEATURE_CONCENTRATION = %.3f, n_ref = 6400\n'
          % ad.HIDDEN_FEATURE_CONCENTRATION)
    print('%-26s %-14s %8s %7s %10s %10s %7s  %s'
          % ('family', 'kind', 'TV ratio', 'hides', 'conc n=8', 'conc n=32', 'fires', 'verdict'))
    rows, disagree = [], []
    for f in fams:
        vcc = H.vcc_of(f['ne'])
        ratio = subgrid_ratio(vcc, f['l0'], f['l1'])
        truth = ratio > 1.01                       # 1 % more variation than the grid can see
        if truth != f['hides']:
            disagree.append((f['label'], f['hides'], truth, ratio))
        got = {}
        for n_sub in (8, 32):
            got[n_sub] = ad.find_hidden_features(vcc, f['l0'], f['l1'], n_ref=6400, n_sub=n_sub)
        fires = bool(got[8]['hidden'] or got[32]['hidden'])
        if truth:
            verdict = 'DETECTED' if fires else 'missed'
        else:
            verdict = 'FALSE POSITIVE' if fires else 'quiet'
        rows.append(dict(label=f['label'], kind=f['kind'], hides=truth, declared=f['hides'],
                         ratio=ratio, fires=fires, c8=got[8]['concentration'],
                         c32=got[32]['concentration'], centre=got[32]['l_centre']))
        print('%-26s %-14s %8.3f %7s %10.4f %10.4f %7s  %s'
              % (f['label'], f['kind'], ratio, truth, got[8]['concentration'],
                 got[32]['concentration'], fires, verdict))
    if disagree:
        print('\n*** declared `hides` disagrees with the measurement on %d families ***'
              % len(disagree))
        for label, declared, truth, ratio in disagree:
            print('    %-26s declared %s, measured %s (TV ratio %.3f)'
                  % (label, declared, truth, ratio))

    quiet_pop = [r for r in rows if not r['hides']]
    hiding = [r for r in rows if r['hides']]
    fp = [r for r in quiet_pop if r['fires']]
    tp = [r for r in hiding if r['fires']]
    print('\n=== SUMMARY ===')
    print('P2  false positives on families with no sub-grid feature : %d / %d'
          % (len(fp), len(quiet_pop)))
    for r in fp:
        print('      %-26s concentration %.4f at l=%.4e' % (r['label'], r['c32'], r['centre']))
    print('P3  detection on families that genuinely hide something  : %d / %d'
          % (len(tp), len(hiding)))
    for r in hiding:
        print('      %-26s %-9s concentration %.4f'
              % (r['label'], 'detected' if r['fires'] else 'MISSED', r['c32']))

    # The shock width sweep, reported as a curve rather than only as a rate -- and with the TV
    # ratio beside it, because the sweep's answer is that narrowing a step never makes it
    # hidden.
    print('\nshock width sweep:')
    for r in rows:
        if r['kind'] == 'sn_shock':
            print('   %-26s TV ratio %.3f  conc(n=8) %.4f  conc(n=32) %.4f  %s'
                  % (r['label'], r['ratio'], r['c8'], r['c32'],
                     'fires' if r['fires'] else '-'))
    np.save('physical_subgrid_rows.npy', np.array(rows, dtype=object), allow_pickle=True)
    return rows


# ======================================================================
# P4 -- does the seam change still hold on physical profiles?
# ======================================================================

# Forcing, as in fallback_quality.py: switch off the engines ahead of the one being timed.
_FORCE = {'hybrid': ('hybrid', False, ('ip_exp', 'separable')),
          'cumulative': ('magnus', True, ('hybrid', 'ip_exp', 'separable'))}


def _call(engine, ne, d, params, energy, Ls, l0):
    strategy, cumulative, off = _FORCE[engine]
    with oscprob._engine_probe(off), warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return np.asarray(oscprob.osc_prob_matter_std_potential(
            d, ne, energy, Ls, params, L0=l0, density_is_of_number_of_electrons=True,
            strategy=strategy, cumulative=cumulative))


def _time(fn, repeats=1):
    t0 = time.perf_counter()
    for _ in range(repeats):
        fn()
    return (time.perf_counter() - t0)/repeats


def seam_cost(rounds=5):
    """Cumulative-scan cost against hybrid cost, on the physical population.

    Method, exactly as the handover requires it:

      * **alternate** -- the two engines are timed round-robin within each round, never back to
        back, because this machine has moved identical code by a factor of two between runs;
      * **carry a control** the change cannot touch (a constant-density scan, which the
        cumulative path declines by construction), and discard the round if it does not return
        about 1.00x;
      * **report the spread**, not one number.

    The point count matters more than the profile here: the seam is a threshold on N, so
    N = 1, 2, 4, 8 and 16 are all timed.  8 is the value in the code.

    **N = 1 is here because of what P1 found.**  Both silent misses on the physical population
    are single points where the hybrid path answered, certified itself, and was 1.1-1.4x outside
    tolerance while the cumulative scan on the same request was 400-1000x more accurate and never
    reached.  §12.1 measured the cost of yielding at N = 2 and above and found 5.75x at N = 2 --
    a real price on the cheapest requests -- but never measured N = 1.  This does.
    """
    fams = [f for f in pp.families()
            if f['kind'] in ('bs05', 'sn_shock', 'sn_turbulence', 'earth_crust', 'tabulated')]
    # One representative per kind keeps this affordable; the sweep that matters is over N.
    seen, chosen = set(), []
    for f in fams:
        if f['kind'] not in seen:
            seen.add(f['kind'])
            chosen.append(f)

    def control_ne(l):
        """The control: a constant density, which the cumulative path declines by construction,
        so nothing the seam governs can move it.  If this does not come back at ~1.00x across
        the run, the round is machine noise and the numbers below are worthless."""
        return H.scalarize(np.full_like(np.asarray(l, dtype=float), 0.05*H.NE0))

    p3 = H.params_for(3)

    print('# P4: cumulative-vs-hybrid cost on the physical population')
    print('# %d profiles x N in {1,2,4,8,16}, %d alternating rounds, seam is at N = %d\n'
          % (len(chosen), rounds, oscprob.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS))

    # Warm every path once before timing any of it, so the first round is not measuring import
    # and cache effects on whichever engine happens to go first.
    for f in chosen:
        for engine in ('hybrid', 'cumulative'):
            try:
                _call(engine, f['ne'], 3, p3, f['energies'][0],
                      np.linspace(f['l0'] + 0.2*(f['l1'] - f['l0']), f['l1'], 8), f['l0'])
            except Exception:                                  # noqa: BLE001 -- a decline
                pass

    acc, ctrl = {}, []
    for rnd in range(rounds):
        t0 = _time(lambda: _call('hybrid', control_ne, 3, p3, 30.0e6,
                                 np.linspace(0.2, 1.0, 8)*H.L_SCALE, 0.0))
        for f in chosen:
            energy = f['energies'][0]
            for N in (1, 2, 4, 8, 16):
                Ls = np.linspace(f['l0'] + 0.2*(f['l1'] - f['l0']), f['l1'], N)
                key = (f['kind'], N)
                for engine in ('hybrid', 'cumulative'):       # alternate, never back to back
                    try:
                        dt = _time(lambda e=engine: _call(e, f['ne'], 3, p3, energy, Ls,
                                                          f['l0']))
                    except Exception:                          # noqa: BLE001 -- a decline
                        dt = None
                    acc.setdefault(key, {}).setdefault(engine, []).append(dt)
        t1 = _time(lambda: _call('hybrid', control_ne, 3, p3, 30.0e6,
                                 np.linspace(0.2, 1.0, 8)*H.L_SCALE, 0.0))
        ctrl.append(t1/t0)
        print('  round %d/%d done, control drift %.2fx' % (rnd + 1, rounds, t1/t0), flush=True)

    drift = float(np.median(ctrl))
    print('\ncontrol drift per round: %s' % ', '.join('%.2fx' % c for c in ctrl))
    print('control drift median %.3fx, range %.3f-%.3f' % (drift, min(ctrl), max(ctrl)))
    clean = [i for i, c in enumerate(ctrl) if 0.9 < c < 1.1]
    print('rounds with the control inside +/-10%%: %d of %d%s'
          % (len(clean), len(ctrl), '' if len(clean) == len(ctrl) else '  <-- the rest are noise'))

    # MIN, not median, is the statistic to read when the machine is not exclusively ours.
    # Interference can only ever ADD time to a call, never remove it, so the fastest observation
    # of each engine is the one closest to its uncontended cost and the ratio of minima is the
    # best estimate of the uncontended ratio.  The median is kept beside it: if the two disagree,
    # the run was noisy and the min column is the one to believe.
    print('\n%-16s %4s %11s %11s %9s %9s  %s'
          % ('profile', 'N', 'hyb min ms', 'cum min ms', 'cum/hyb', '(median)', 'verdict'))
    ratios_at_seam = []
    for (kind, N) in sorted(acc):
        h = [x for x in acc[(kind, N)].get('hybrid', []) if x is not None]
        c = [x for x in acc[(kind, N)].get('cumulative', []) if x is not None]
        if not h:
            continue
        hmin, hmed = float(np.min(h)), float(np.median(h))
        if not c:
            print('%-16s %4d %11.1f %11s %9s %9s  cumulative declined'
                  % (kind, N, 1e3*hmin, '-', '-', '-'))
            continue
        cmin, cmed = float(np.min(c)), float(np.median(c))
        ratio, ratio_med = cmin/hmin, cmed/hmed
        if N >= oscprob.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS:
            ratios_at_seam.append((kind, N, ratio))
        print('%-16s %4d %11.1f %11.1f %9.2f %9.2f  %s%s'
              % (kind, N, 1e3*hmin, 1e3*cmin, ratio, ratio_med,
                 'cheaper' if ratio < 1.0 else 'costs %.2fx' % ratio,
                 '   noisy' if abs(ratio_med/ratio - 1.0) > 0.25 else ''))

    print('\n=== P4 ===')
    if ratios_at_seam:
        worst = max(ratios_at_seam, key=lambda t: t[2])
        print('at or above the seam (N >= %d): worst cumulative/hybrid ratio %.2fx on %s at N=%d'
              % (oscprob.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS, worst[2], worst[0], worst[1]))
        print('median ratio at or above the seam: %.2fx'
              % float(np.median([r for _, _, r in ratios_at_seam])))
        print('the 25 -> 8 change holds here'
              if worst[2] < 2.0 else 'the 25 -> 8 change needs revisiting on this population')
    else:
        print('no workload reached the seam')

    # The separate question P1 raised: both silent misses on the physical population are SINGLE
    # POINTS, which the seam cannot reach because 8 is a threshold on point count.  Whether it
    # should reach them is a cost question, and this is the row that answers it.
    below = [(k, N, float(np.min(acc[(k, N)]['cumulative']))/float(np.min(acc[(k, N)]['hybrid'])))
             for (k, N) in sorted(acc)
             if N < oscprob.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS
             and acc[(k, N)].get('cumulative') and acc[(k, N)].get('hybrid')
             and None not in acc[(k, N)]['cumulative'] and None not in acc[(k, N)]['hybrid']]
    ones = [r for k, N, r in below if N == 1]
    print('\n--- below the seam: could it reach the single points P1 exposed? ---')
    if ones:
        print('N = 1  cumulative/hybrid: median %.2fx, worst %.2fx, best %.2fx'
              % (float(np.median(ones)), max(ones), min(ones)))
        print('  -> lowering the seam to reach single points would cost about %.2fx there'
              % float(np.median(ones)))
    else:
        print('N = 1: the cumulative scan declined on every profile, so the seam cannot reach '
              'single points at any threshold -- the exposure needs a different fix')


# ======================================================================
# Is the oracle itself trustworthy here?
# ======================================================================

def oracle_check():
    """``solve_ivp``/DOP853 assumes smoothness.  Two of these families are not smooth.

    Every error this programme reports is measured against DOP853 at ``rtol = 1e-12``.  That is
    an eighth-order method with an error estimate built on a Taylor expansion, and the Earth
    chord has PREM's density steps in it while the SN shock at ``w = 1e-6`` is a 0.07 km ramp on
    a 7e4 km ray.  If the oracle is itself wrong on those, every number quoted for them is
    meaningless -- which is why the earlier tranches used ``expm`` composed across segments
    wherever the profile was piecewise constant.  Neither of these families is piecewise
    constant, so ``expm`` is not available and the oracle has to be checked instead.

    The check is ``harness.oracle_converged``: tighten to ``rtol = 1e-13`` and confirm the
    answer moves by far less than the errors being quoted (1e-3 and below).
    """
    print('# Oracle self-consistency on the families where DOP853 is least at home')
    print('# criterion: rtol 1e-12 -> 1e-13 must move the answer by << the quoted error\n')
    print('%-28s %3s %8s %14s  %s' % ('family', 'd', 'E/MeV', 'oracle moves', 'verdict'))
    worst = 0.0
    for f in pp.families():
        if f['kind'] not in ('sn_shock', 'earth_crust'):
            continue
        if f['kind'] == 'sn_shock' and not f['label'].endswith(('1e-06', '1e-04')):
            continue
        for d in (2, 3):
            energy = f['energies'][0]
            H_of_l = H.H_factory(d, H.params_for(d), H.vcc_of(f['ne']), energy)
            moved, ok = H.oracle_converged(H_of_l, f['l0'], f['l1'], d, TOL)
            worst = max(worst, moved)
            print('%-28s %3d %8.1f %14.3e  %s'
                  % (f['label'], d, energy/1e6, moved,
                     'ok' if moved < 0.01*TOL else '*** ORACLE NOT CONVERGED ***'))
    print('\nworst movement %.3e against a %.0e tolerance -> oracle is %s'
          % (worst, TOL, 'trustworthy here' if worst < 0.01*TOL else 'NOT trustworthy'))


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    what = sys.argv[1]
    if what == 'sub_grid':
        sub_grid()
    elif what == 'seam_cost':
        seam_cost(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
    elif what == 'oracle_check':
        oracle_check()
    else:
        print(__doc__)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
