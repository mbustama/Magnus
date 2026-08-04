# -*- coding: utf-8 -*-
"""Battery 3 -- the routing seams.

Every threshold is a discontinuity, and discontinuities are where behaviour hides.

PASS CRITERIA, stated before running:
  (D1) accuracy must not jump DISCONTINUOUSLY across N = 25.  A large accuracy
       discontinuity at the seam is a defect even if both sides are inside tolerance.
       Quantified as err(N=24)/err(N=26) and its inverse; flagged above 10x.
  (D2) every answer meets 1e-3 or warns (silent miss = failure).
  (D3) escape hatches behave: strategy='magnus'/'hybrid' and cumulative=False/True do what
       they promise; cumulative=True RAISES when inapplicable, never falls back silently.
  (D4) a position-independent Hamiltonian must NEVER reach the cumulative scan, at every N
       and every d.
  (D5) t_slab_edges -> cumulative declines; t_breakpoints -> cumulative accepts AND the
       breakpoints are genuinely present in the grid it builds.
"""

import sys
import warnings

import numpy as np

import harness as H
import harness6 as H6
import magnus.oscprob as oscprob

TOL = 1e-3
E_DEF = 50.0e6
L_MAX = 0.5*H.R_SUN
NS = [1, 2, 3, 8, 24, 25, 26, 100]

WARN_OK = ('ToleranceNotAchievedWarning', 'MagnusConvergenceWarning',
           'HybridCertificationWarning')


# ---------------------------------------------------------------- profiles
def castle_wall_ne(n_walls=20, l_ini=None, l_fin=None, lo=None, hi=None):
    l_ini = l_ini if l_ini is not None else 0.05*L_MAX
    l_fin = l_fin if l_fin is not None else L_MAX
    lo = lo if lo is not None else 0.02*H.NE0
    hi = hi if hi is not None else 0.20*H.NE0

    def ne(l):
        x = np.asarray(l, dtype=float)
        u = (x - l_ini)/(l_fin - l_ini)
        idx = np.floor(np.clip(u, 0.0, 1.0 - 1e-15)*n_walls)
        return H.scalarize(np.where((x < l_ini) | (x > l_fin), lo,
                                   np.where(idx % 2 == 0, lo, hi)))
    edges = np.linspace(l_ini, l_fin, n_walls + 1)
    return ne, edges


def noisy_ne(seed=3, n_modes=25, amp=0.6):
    rng = np.random.default_rng(seed)
    ks = np.arange(1, n_modes + 1)
    ph = rng.uniform(0, 2*np.pi, n_modes)
    w = rng.normal(size=n_modes)/ks
    w = w/np.abs(w).sum()

    def ne(l):
        x = np.asarray(l, dtype=float)[..., None]
        s = (w*np.sin(2*np.pi*ks*x/L_MAX + ph)).sum(axis=-1)
        return H.scalarize(0.1*H.NE0*np.exp(-np.asarray(l, float)/(3*H.L_SCALE))
                           * (1.0 + amp*s))
    return ne


PROFILES = {}


def build_profiles():
    cw, cw_edges = castle_wall_ne()
    PROFILES['solar exponential'] = (H.solar_ne(), {})
    PROFILES['multi-resonance'] = (H.modulated_ne(amp=0.9, n_cycles=6.0, span=L_MAX), {})
    PROFILES['castle wall + breakpoints'] = (cw, {'t_breakpoints': cw_edges})
    PROFILES['noisy'] = (noisy_ne(), {})


# ---------------------------------------------------------------- spy
class Spy:
    """Records which internal path answered."""

    def __init__(self):
        self.cum = 0
        self.hyb = 0

    def __enter__(self):
        self._c, self._h = oscprob._osc_prob_cumulative_scan, oscprob._osc_prob_hybrid_dispatch

        def cum(*a, **k):
            self.cum += 1
            return self._c(*a, **k)

        def hyb(*a, **k):
            r = self._h(*a, **k)
            if r is not NotImplemented:
                self.hyb += 1
            return r

        oscprob._osc_prob_cumulative_scan = cum
        oscprob._osc_prob_hybrid_dispatch = hyb
        return self

    def __exit__(self, *e):
        oscprob._osc_prob_cumulative_scan = self._c
        oscprob._osc_prob_hybrid_dispatch = self._h
        return False


def run(ne, kw, N, d=2, strategy='auto', cumulative=None, energy=E_DEF, **extra):
    p = H.params_for(d)
    Ls = np.linspace(0.05*L_MAX, L_MAX, N) if N > 1 else np.array([L_MAX])
    call_kw = dict(kw)
    call_kw.update(extra)
    if cumulative is not None:
        call_kw['cumulative'] = cumulative
    with Spy() as sp, H.Caught() as c:
        P = np.asarray(oscprob.osc_prob_matter_std_potential(
            d, ne, energy, Ls if N > 1 else float(Ls[0]), p, L0=0.0,
            density_is_of_number_of_electrons=True, strategy=strategy, **call_kw))
    P = P.reshape((N, d, d))
    H_of_l = H6.H_family('std', d, energy, H.vcc_of(ne), p)
    Pref = np.array([H.P_of(U) for U in H.exact_U_many(H_of_l, 0.0, Ls, d)])
    err = H.maxabs(P - Pref)
    return err, c.names, sp.cum, sp.hyb


# ---------------------------------------------------------------- sub-tests
def sub1():
    print('## 3.1  N across the seam (HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS = 25), '
          'per profile')
    print('%-28s %s' % ('profile', ''.join('%12s' % ('N=%d' % n) for n in NS)))
    results = {}
    for name, (ne, kw) in PROFILES.items():
        errs, paths = [], []
        for N in NS:
            err, wn, cum, hyb = run(ne, kw, N)
            errs.append(err)
            paths.append('C' if cum else ('H' if hyb else 'g'))
            silent = err > TOL and not any(x in WARN_OK for x in wn)
            if silent:
                paths[-1] += '!'
        results[name] = (errs, paths)
        print('%-28s %s' % (name, ''.join('%12.2e' % e for e in errs)))
        print('%-28s %s' % ('  path (C/H/g, !=silent miss)',
                            ''.join('%12s' % p for p in paths)))
    print('\n  seam check (criterion D1): err(N=24) vs err(N=26)')
    for name, (errs, _) in results.items():
        e24, e26 = errs[NS.index(24)], errs[NS.index(26)]
        ratio = max(e24, e26)/max(min(e24, e26), 1e-300)
        print('    %-28s %.2e -> %.2e   %.1fx %s'
              % (name, e24, e26, ratio, 'DISCONTINUITY' if ratio > 10 else ''))
    print()


def sub2():
    print('## 3.2  the 2 <= N < 25 band on multi-resonance profiles specifically')
    print('    (the handover: "commit 8 should have improved this band too; confirm it did")')
    for cyc in (4.0, 6.0, 10.0):
        ne = H.modulated_ne(amp=0.9, n_cycles=cyc, span=L_MAX)
        row = []
        for N in (2, 5, 8, 16, 24, 25, 40):
            err, wn, cum, hyb = run(ne, {}, N)
            silent = err > TOL and not any(x in WARN_OK for x in wn)
            row.append('%9.2e%s%s' % (err, 'C' if cum else ('H' if hyb else 'g'),
                                      '!' if silent else ' '))
        print('  n_cycles=%-5.1f  %s' % (cyc, ' '.join(row)))
    print('  (columns: N = 2, 5, 8, 16, 24, 25, 40;  C=cumulative H=hybrid g=general, '
          '!=silent miss)\n')


def sub3():
    print('## 3.3  escape hatches at every N')
    ne = H.solar_ne()
    for N in (1, 2, 24, 25, 100):
        out = []
        for strat in ('auto', 'magnus', 'hybrid'):
            err, wn, cum, hyb = run(ne, {}, N, strategy=strat)
            out.append('%s:%8.2e[%s]' % (strat, err, 'C' if cum else ('H' if hyb else 'g')))
        for cval in (False, True):
            try:
                err, wn, cum, hyb = run(ne, {}, N, cumulative=cval)
                out.append('cum=%s:%8.2e[%s]' % (cval, err, 'C' if cum else 'g'))
            except Exception as ex:      # noqa: BLE001
                out.append('cum=%s:RAISED(%s)' % (cval, type(ex).__name__))
        print('  N=%-4d %s' % (N, '  '.join(out)))
    print('\n  cumulative=True on requests it cannot serve (must RAISE, never fall back):')
    p = H.params_for(2)
    Ls = np.linspace(0.05*L_MAX, L_MAX, 30)
    bad = [
        ('differing energies', dict(energy=np.linspace(10e6, 50e6, 30), L=Ls)),
        ('a baseline behind L0', dict(energy=E_DEF, L=np.r_[-1.0, Ls[1:]])),
        ('t_slab_edges given', dict(energy=E_DEF, L=Ls,
                                    t_slab_edges=np.column_stack(
                                        [np.linspace(0, L_MAX, 11)[:-1],
                                         np.linspace(0, L_MAX, 11)[1:]]))),
    ]
    for name, kw in bad:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                oscprob.osc_prob_energy_baseline(
                    lambda e, l: H6.H_family('std', 2, e, H.vcc_of(ne), p)(l),
                    kw.pop('energy'), kw.pop('L'), 0.0, cumulative=True, **kw)
            print('    %-24s NO RAISE   <-- criterion D3 FAILURE' % name)
        except Exception as ex:          # noqa: BLE001
            print('    %-24s raises %s  ok' % (name, type(ex).__name__))
    print()


def sub4():
    print('## 3.4  position-INDEPENDENT H must never reach the cumulative scan '
          '(every N, every d)')
    bad = 0
    for d in (2, 3, 4, 5):
        for N in (1, 2, 24, 25, 100, 400):
            for kind, rho in (('vacuum', 0.0), ('constant density', 0.05*H.NE0)):
                p = H.params_for(d)
                Ls = np.linspace(0.05*L_MAX, L_MAX, N) if N > 1 else np.array([L_MAX])
                with Spy() as sp, warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    oscprob.osc_prob_matter_std_potential(
                        d, rho, E_DEF, Ls if N > 1 else float(Ls[0]), p, L0=0.0,
                        density_is_of_number_of_electrons=True)
                if sp.cum:
                    bad += 1
                    print('    d=%d N=%-4d %-17s REACHED CUMULATIVE  <-- D4 FAILURE'
                          % (d, N, kind))
    print('    %d violations across 48 (d, N, kind) combinations%s'
          % (bad, '' if bad else '  -- ok'))
    print()


def sub5():
    print('## 3.5  t_slab_edges (cumulative must decline) and t_breakpoints '
          '(must be HONOURED, not dropped)')
    cw, edges = castle_wall_ne()
    Ls = np.linspace(0.05*L_MAX, L_MAX, 40)
    p = H.params_for(2)

    # t_slab_edges under cumulative='auto' -> must NOT take the cumulative scan
    pairs = np.column_stack([np.linspace(0, L_MAX, 201)[:-1], np.linspace(0, L_MAX, 201)[1:]])
    with Spy() as sp, warnings.catch_warnings():
        warnings.simplefilter('ignore')
        oscprob.osc_prob_matter_std_potential(
            2, cw, E_DEF, Ls, p, L0=0.0, density_is_of_number_of_electrons=True,
            t_slab_edges=pairs)
    print('    t_slab_edges + cumulative="auto": cumulative used = %d  %s'
          % (sp.cum, 'ok' if sp.cum == 0 else '<-- D5 FAILURE'))

    # t_breakpoints: confirm they are genuinely present in the grid the scan builds
    n_acc = 500
    grid, out_idx = oscprob._cumulative_scan_grid(Ls, 0.0, n_acc, edges)
    interior = edges[(edges > 0.0) & (edges < Ls[-1])]
    present = np.array([np.min(np.abs(grid - b)) == 0.0 for b in interior])
    print('    _cumulative_scan_grid: %d/%d interior breakpoints present exactly  %s'
          % (present.sum(), len(interior), 'ok' if present.all() else '<-- D5 FAILURE'))
    on_edge = np.array([grid[i] == L for i, L in zip(out_idx, Ls)])
    print('    every requested baseline lands exactly on a slab edge: %d/%d  %s'
          % (on_edge.sum(), len(Ls), 'ok' if on_edge.all() else '<-- D5 FAILURE'))

    # and end-to-end: does passing them actually improve the castle wall?
    for name, kw in (('without breakpoints', {}), ('with breakpoints',
                                                   {'t_breakpoints': edges})):
        err, wn, cum, hyb = run(cw, kw, 40)
        print('    castle wall N=40 %-20s err=%.3e  path=%s  warns=%s'
              % (name, err, 'C' if cum else ('H' if hyb else 'g'), ','.join(wn) or '-'))
    print()


if __name__ == '__main__':
    build_profiles()
    print('# Battery 3 -- routing seams.  E=%.0f MeV, L_max=%.2f R_sun, tol=%.0e\n'
          % (E_DEF/1e6, L_MAX/H.R_SUN, TOL))
    for w in (sys.argv[1:] or ['1', '2', '3', '4', '5']):
        {'1': sub1, '2': sub2, '3': sub3, '4': sub4, '5': sub5}[w]()
