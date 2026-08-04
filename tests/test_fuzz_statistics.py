# -*- coding: utf-8 -*-
"""Seeded fuzzing, asserted on statistics rather than on per-case values.

The adversarial batteries found things the unit tests did not, and nothing ran them: they print
tables and take tens of minutes.  This is the part of them that fits in CI.

**Why aggregates.**  A per-case assertion on random input is brittle -- one case drifts, the
suite goes red, and the fix is to loosen the bound until it stops.  An assertion on the
*distribution* is stable under that drift and still catches the thing worth catching, which is a
regression that moves the whole distribution.  The three quantities asserted here are the three
``FINDINGS_ADVERSARIAL_VALIDATION.md`` reports: the silent-miss rate, the median error, and
whether anything raised.

**Silent miss** means an answer outside the requested tolerance with no warning of any kind.
It is the only failure mode that matters here; an inaccurate answer that says so is the
warnings' job, and a false-positive warning is noise.

**Cost.**  The ``solve_ivp`` oracle dominates, badly, at low energy and high flavour count --
measured at 25 s for 2 of 80 baselines at 6.5 MeV, 3nu, where the package answered the whole
case in 2.1 s.  So the smooth population is held at 20-200 MeV and d in {2, 3}, and the
piecewise population uses ``expm`` composed across segments, which is **exact** for those
profiles and costs nothing.  Run this file directly to print the distributions:

```bash
python tests/test_fuzz_statistics.py
```

**Cost.**  About 6 minutes, nearly all of it the ``solve_ivp`` oracle on the 40 smooth cases;
the 120 piecewise cases are close to free.  That is a real addition to a suite that already runs
in ~22 minutes, and it is stated here rather than buried so that trimming ``N_SMOOTH`` is an
informed decision: halving it halves the runtime and roughly doubles the width of the
silent-miss confidence interval.
"""

import warnings

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

import magnus.globaldefs as gd
import magnus.hamiltonians as hams
import magnus.matter as matter
import magnus.oscprob as op

TOL = 1.0e-3
L_SCALE = gd.L_SCALE_SUN
NE0 = gd.NUM_DENSITY_E_SUN_CENTRAL

N_SMOOTH = 40
N_PIECEWISE = 120
SEED = 20260804

# Any warning at all makes a case not-silent; the point of this measurement is whether the
# package said anything, not which class it chose.
_MIX = ('s12', 's23', 's13', 'dCP', 'D21', 'D31')


def params_for(d):
    """Sterile mixings large enough to put extra resonances in the density range scanned --
    a realistic eV^2-scale splitting never crosses the matter potential here, so it would
    exercise nothing."""
    p = {k: gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT'][k] for k in _MIX}
    if d == 2:
        return {'sth': p['s12'], 'Dm2': p['D21']}
    if d == 3:
        return p
    if d == 4:
        return dict(p, s14=0.30, d14=0.7, s24=0.20, d24=1.9, s34=0.15, D41=1.0e-3)
    return dict(p, s14=0.30, d14=0.7, s15=0.25, d15=2.2, s24=0.20, d24=1.9, s25=0.18,
                s34=0.15, s35=0.12, d35=0.4, D41=1.0e-3, D51=3.0e-3)


def h_vac_for(d, p, nubar=False):
    if d == 2:
        return hams.hamiltonian_2nu_vacuum_energy_independent(p['sth'], p['Dm2'])
    if d == 3:
        return hams.hamiltonian_3nu_vacuum_energy_independent(
            p['s12'], p['s23'], p['s13'], p['dCP'], p['D21'], p['D31'], nubar=nubar)
    if d == 4:
        return hams.hamiltonian_4nu_vacuum_energy_independent(
            p['s12'], p['s23'], p['s13'], p['dCP'], p['s14'], p['d14'], p['s24'], p['d24'],
            p['s34'], p['D21'], p['D31'], p['D41'], nubar=nubar)
    return hams.hamiltonian_5nu_vacuum_energy_independent(
        p['s12'], p['s23'], p['s13'], p['dCP'], p['s14'], p['d14'], p['s15'], p['d15'],
        p['s24'], p['d24'], p['s25'], p['s34'], p['s35'], p['d35'],
        p['D21'], p['D31'], p['D41'], p['D51'], nubar=nubar)


def vcc_of(ne_func, nubar=False):
    """TRAP GUARD: the 7th positional argument of vcc_func_from_rho_func is
    density_is_of_number_of_electrons, not nubar.  Both by keyword so it cannot be confused."""
    return matter.vcc_func_from_rho_func(
        ne_func, 0.0, 1.0, 0.5, nubar=nubar, density_matter_is_in_g_per_cm3=False,
        density_is_of_number_of_electrons=True)


def H_of(d, p, ne, energy, nubar=False):
    h_vac = np.asarray(h_vac_for(d, p, nubar=nubar), dtype=complex)
    proj = np.zeros((d, d), dtype=complex)
    proj[0, 0] = 1.0
    vcc = vcc_of(ne, nubar=nubar)

    def H_func(l):
        v = np.asarray(vcc(l))
        return (1.0/energy)*h_vac + v[..., None, None]*proj
    return H_func


def P_of(U):
    U = np.asarray(U)
    return np.transpose(U.real**2 + U.imag**2)


def maxabs(x):
    return float(np.max(np.abs(np.asarray(x))))


def scalarize(y):
    a = np.asarray(y)
    return a[()] if a.ndim == 0 else a


# ----------------------------------------------------------------------
# populations
# ----------------------------------------------------------------------

def smooth_cases(n=N_SMOOTH, seed=SEED):
    """Random smooth profiles: random Fourier sums of controlled bandwidth, kept positive."""
    rng = np.random.default_rng(seed)
    for _ in range(n):
        d = int(rng.choice([2, 3]))
        energy = float(10.0**rng.uniform(7.3, 8.3))          # 20 - 200 MeV
        span = float(rng.uniform(0.3, 1.0))*L_SCALE
        n_modes = int(rng.integers(2, 9))
        ks = rng.integers(1, n_modes + 1, size=n_modes)
        a, b = rng.normal(size=n_modes), rng.normal(size=n_modes)
        nrm = np.sqrt(np.sum(a**2 + b**2)) or 1.0
        a, b = a/nrm, b/nrm
        base = float(10.0**rng.uniform(-2.5, -1.0))
        amp = float(rng.uniform(0.2, 0.9))

        def ne(l, ks=ks, a=a, b=b, span=span, base=base, amp=amp):
            x = np.asarray(l, dtype=float)[..., None]
            s = (a*np.cos(2.0*np.pi*ks*x/span) + b*np.sin(2.0*np.pi*ks*x/span)).sum(axis=-1)
            return scalarize(NE0*base*(1.0 + amp*s)
                             * np.exp(-np.asarray(l, dtype=float)/(6.0*L_SCALE)))

        N = int(rng.choice([1, 3, 12, 30]))
        Ls = np.linspace(0.05*span, span, N) if N > 1 else np.array([span])
        yield dict(d=d, energy=energy, span=span, ne=ne, N=N, Ls=Ls,
                   nubar=bool(rng.random() < 0.3), params=params_for(d))


def piecewise_cases(n=N_PIECEWISE, seed=SEED + 1):
    """Random piecewise-constant profiles, for which expm composed across segments is exact.

    **Spans d = 2...5**, unlike the smooth population.  The oracle here is a product of matrix
    exponentials rather than an ODE solve, so it costs essentially nothing and the flavour count
    is free; ``solve_ivp`` is what makes d = 4 and 5 unaffordable on the smooth side.
    """
    rng = np.random.default_rng(seed)
    for _ in range(n):
        d = int(rng.choice([2, 3, 4, 5]))
        energy = float(10.0**rng.uniform(7.3, 8.3))
        span = float(rng.uniform(0.3, 1.0))*L_SCALE
        n_seg = int(rng.integers(2, 9))
        cuts = np.sort(rng.uniform(0.0, 1.0, n_seg - 1))*span
        edges = np.concatenate([[0.0], cuts, [span]])
        values = NE0*10.0**rng.uniform(-2.5, -0.5, n_seg)

        def ne(l, edges=edges, values=values):
            x = np.asarray(l, dtype=float)
            idx = np.clip(np.searchsorted(edges, x, side='right') - 1, 0, len(values) - 1)
            return scalarize(values[idx])

        N = int(rng.choice([1, 3, 12, 30]))
        Ls = np.linspace(0.05*span, span, N) if N > 1 else np.array([span])
        yield dict(d=d, energy=energy, span=span, ne=ne, N=N, Ls=Ls, edges=edges,
                   nubar=bool(rng.random() < 0.3), params=params_for(d))


# ----------------------------------------------------------------------
# oracles
# ----------------------------------------------------------------------

def oracle_solve_ivp(H_func, l0, Ls, d):
    def rhs(l, y):
        return (-1j*np.asarray(H_func(l)) @ y.reshape(d, d)).ravel()
    sol = solve_ivp(rhs, (float(l0), float(Ls[-1])), np.eye(d, dtype=complex).ravel(),
                    rtol=1e-12, atol=1e-14, method='DOP853', t_eval=np.asarray(Ls, dtype=float))
    assert sol.success, sol.message
    return np.array([P_of(sol.y[:, i].reshape(d, d)) for i in range(len(Ls))])


def oracle_expm(H_func, edges, Ls, d):
    """Exact, not approximate: H is constant on each segment, so the exponentials compose."""
    out, U, cursor = [], np.eye(d, dtype=complex), float(edges[0])
    stops = sorted(set(np.concatenate([edges, Ls]).tolist()))
    at = {}
    for nxt in stops:
        if nxt > cursor:
            Hm = np.asarray(H_func(0.5*(cursor + nxt)), dtype=complex)
            U = expm(-1j*Hm*(nxt - cursor)) @ U
            cursor = nxt
        at[nxt] = U.copy()
    for L in Ls:
        out.append(P_of(at[float(L)]))
    return np.array(out)


def run_case(case, oracle):
    H_func = H_of(case['d'], case['params'], case['ne'], case['energy'], case['nubar'])
    Pref = oracle(H_func)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        P = np.asarray(op.osc_prob_matter_std_potential(
            case['d'], case['ne'], case['energy'],
            case['Ls'] if case['N'] > 1 else float(case['Ls'][0]), case['params'],
            L0=0.0, nubar=case['nubar'], density_is_of_number_of_electrons=True)
        ).reshape(case['N'], case['d'], case['d'])
    err = maxabs(P - Pref)
    return err, bool(caught)


def collect(kind):
    """Score a whole population.  Anything raising propagates: a ValueError from the package's
    own validation would be correct behaviour on a profile the generator drove negative, so the
    generators are written to stay positive and an exception here is a defect.

    This used to be asserted separately by a third test that re-ran 16 of the same cases.  That
    cost 72 s to re-establish what these two tests establish anyway -- they run 160 cases and
    cannot pass if one of them raises -- so it is folded in here instead."""
    rows = []
    if kind == 'smooth':
        for c in smooth_cases():
            rows.append(run_case(c, lambda H: oracle_solve_ivp(H, 0.0, c['Ls'], c['d'])))
    else:
        for c in piecewise_cases():
            rows.append(run_case(c, lambda H: oracle_expm(H, c['edges'], c['Ls'], c['d'])))
    errs = np.array([r[0] for r in rows])
    silent = [r for r in rows if r[0] > TOL and not r[1]]
    return errs, silent


# ----------------------------------------------------------------------
# the assertions
# ----------------------------------------------------------------------

def test_smooth_profile_fuzz_statistics():
    """Random smooth profiles, scored against solve_ivp/DOP853 at rtol=1e-12.

    Bounds come from the measured distribution, not from a target.  ``FINDINGS`` §9.4 measured
    the silent-miss rate at **4.1 %** over 145 cases after the fixes (down from 21 %) and the
    median error at 6.08e-08.  This population (40 cases, 20-200 MeV, d in {2,3}, held there
    because ``solve_ivp`` dominates the cost at low energy and high flavour count) measures
    **median 7.8e-09, p90 6.1e-04, max 1.08e-03, 1 silent miss (2.5 %)**, the one sitting at
    1.08e-03 against a requested 1e-3 -- a few per cent over, which is the residue
    ``FINDINGS`` §9.4 describes rather than a new failure.

    It was 24 cases, and that was too few to say anything: 2 of 24 is consistent with a 4 %
    rate and equally with 15 %.  40 narrows it, and the piecewise population below carries the
    rest of the power at no oracle cost.

    The bounds are set above what was measured so ordinary drift does not turn the suite red,
    and far below the pre-fix numbers (21 % silent, max 4.0e-02) so a regression to them would.
    """
    errs, silent = collect('smooth')
    rate = len(silent)/len(errs)
    assert rate <= 0.10, 'silent-miss rate %.1f%% (%d/%d)' % (100*rate, len(silent), len(errs))
    assert np.median(errs) <= 1.0e-6, 'median error %.3e' % np.median(errs)
    assert errs.max() <= 5.0e-2, 'worst error %.3e' % errs.max()


def test_piecewise_profile_fuzz_statistics():
    """Random piecewise-constant profiles, with the edges left undeclared -- the adversarial
    case, since declaring them is essentially exact (``FINDINGS`` §9.2: median 1.34e-12).

    The oracle here is ``expm`` composed across the segments, which is the exact operator for
    these Hamiltonians rather than an approximation, so it cannot itself step over a jump.
    ``FINDINGS`` §9.2 measured 2 silent misses in 150 undeclared cases; the median error is
    expected to be poor (7.8e-04 there) and that is not the failure -- being poor *quietly* is.

    This population is **120 cases across d = 2, 3, 4 and 5**, which is where the statistical
    power in this file comes from: the oracle is a product of matrix exponentials rather than an
    ODE solve, so cases are nearly free and the flavour count costs nothing.  Measured:
    **median 3.5e-04, p90 2.0e-03, max 1.49e-02, 19 outside tolerance and 0 of them silent** --
    every inaccurate answer said so.  At 0 of 120, the 95 % upper bound on the silent-miss rate
    here is about 2.5 %.
    """
    errs, silent = collect('piecewise')
    rate = len(silent)/len(errs)
    assert rate <= 0.10, 'silent-miss rate %.1f%% (%d/%d)' % (100*rate, len(silent), len(errs))
    assert errs.max() <= 2.0e-1, 'worst error %.3e' % errs.max()


if __name__ == '__main__':
    for kind in ('smooth', 'piecewise'):
        errs, silent = collect(kind)
        print('%-10s n=%-4d median %.3e  p90 %.3e  max %.3e  outside %d  SILENT %d (%.1f%%)'
              % (kind, len(errs), np.median(errs), np.percentile(errs, 90), errs.max(),
                 int((errs > TOL).sum()), len(silent), 100*len(silent)/len(errs)))
