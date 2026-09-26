# -*- coding: utf-8 -*-
r"""Tests of the energy-batched engine on grids with breakpoints (issue #71).

The breakpoints are inserted into every refinement level's grid, so at small slab counts a
nominal x1.5 step barely refines the real grid: on the core-crossing PREM chord, 4 -> 6 slabs is
20 -> 22 points.  Two such grids agree without having converged, and the engine accepted that:
at two flavors and rtol = atol = 1e-3, two of twelve energies came back 2.5e-3 off, with no
warning.  The per-point ladder already refused such agreements (``MIN_EFFECTIVE_REFINEMENT``,
see ``test_tolerance.py``).  The batched engine now applies the same bound on grids with
breakpoints, and grows the real grid rather than the nominal slab count, so that it does not
spend levels on steps the bound would refuse.

The reference is ``notebooks/prem_chord_reference.json``: extended precision, Richardson-
extrapolated, self-converged to 1e-12.
"""

import json
import pathlib
import sys
import warnings

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.magnus as mg
import magnus.oscprob as op

NOTEBOOKS = pathlib.Path(__file__).resolve().parents[1]/'notebooks'
sys.path.insert(0, str(NOTEBOOKS))


@pytest.fixture(scope='module')
def chord():
    import gen_profile_benchmarks as gpb
    import prem_chord_common as pcc
    ch = pcc.chord()
    refs = json.loads((NOTEBOOKS/'prem_chord_reference.json').read_text())
    per_ne = gpb.matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
    return dict(
        ne=lambda x: ch['vcc'](x)/per_ne, L=ch['baseline'], bp=ch['edges'][1:-1],
        E=np.asarray(refs['energy_ev'], dtype=float), osc=gpb.osc_params,
        ref={int(c['flavours']): np.asarray(c['reference']) for c in refs['cases']})


def scan(chord, d, tol, **kw):
    r"""P(nu_mu -> nu_mu) over the chord's energies, the slab counts of the levels computed,
    and the tolerance warnings raised."""
    seen = []
    real = mg.evolution_operators_from_samples

    def spy(At, widths, *args, **kwargs):
        seen.append(len(widths))
        return real(At, widths, *args, **kwargs)

    mg.evolution_operators_from_samples = spy
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            P = np.asarray(op.osc_prob_matter_std_potential(
                d, chord['ne'], chord['E'], chord['L'], chord['osc'](d), L0=0.0,
                nu_i=gd.NUMU, nu_f=gd.NUMU, density_is_of_number_of_electrons=True,
                t_breakpoints=chord['bp'], rtol=tol, atol=tol, **kw))
    finally:
        mg.evolution_operators_from_samples = real
    levels = [n for i, n in enumerate(seen) if i == 0 or n != seen[i - 1]]
    names = {w.category.__name__ for w in caught
             if issubclass(w.category, op.ToleranceNotAchievedWarning)}
    return P, levels, names


@pytest.mark.parametrize('strategy', ['auto', 'magnus'])
def test_the_prem_chord_scan_meets_its_tolerance(chord, strategy):
    r"""The case the issue is about, on the default route and the explicit one."""
    tol = 1e-3
    P, _, names = scan(chord, 2, tol, strategy=strategy)
    ref = chord['ref'][2]
    assert not names
    assert np.all(np.abs(P - ref) <= tol + tol*np.abs(ref))


@pytest.mark.parametrize('growth', [1.5, 1.1])
def test_every_level_refines_the_real_grid(chord, growth):
    r"""Points of consecutive levels (slabs + 1) differ by at least the bound, also when the
    user's growth factor is below it: the ladder then grows by the bound instead of refusing
    every agreement until it runs out of levels."""
    _, levels, names = scan(chord, 2, 1e-3, strategy='magnus', growth_factor_n_slabs=growth)
    points = np.asarray(levels) + 1
    assert len(points) >= 2 and not names
    assert np.all(points[1:]/points[:-1] >= op.MIN_EFFECTIVE_REFINEMENT)


def test_a_level_cut_short_by_the_cap_warns(chord):
    r"""With ``max_n_slabs=20`` the last step is 30 -> 37 points, x1.23: below the bound, so
    its agreement cannot certify, and the caps end the ladder with a warning."""
    _, levels, names = scan(chord, 2, 1e-3, strategy='magnus', max_n_slabs=20)
    points = np.asarray(levels) + 1
    assert points[-1]/points[-2] < op.MIN_EFFECTIVE_REFINEMENT
    assert 'ToleranceNotAchievedWarning' in names


def test_without_breakpoints_a_capped_step_still_counts():
    r"""The bound is not applied without breakpoints, where nominal and real steps coincide:
    a smooth profile keeps today's results and warnings bit for bit, and a last step cut
    short by the cap (108 -> 115 slabs here) certifies as before."""
    rho = lambda l: 3e3*np.exp(-np.asarray(l, dtype=float)/(100.0*gd.UNIT_KM))
    osc = gd.load_nufit_params('NuFIT 6.1')
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        op.osc_prob_matter_std_potential(
            3, rho, np.linspace(0.1, 0.4, 40)*gd.UNIT_GEV, 250.0*gd.UNIT_KM, osc, L0=0.0,
            density_matter_is_in_g_per_cm3=True, rtol=1e-6, atol=1e-6, strategy='magnus',
            max_n_slabs=115)
    assert not any(issubclass(w.category, op.ToleranceNotAchievedWarning) for w in caught)


@pytest.mark.parametrize('d', [2, 4])
def test_the_earth_wrapper_with_a_layered_electron_fraction_meets_its_tolerance(d):
    r"""The Earth wrappers declare the layer crossings themselves and use a layered Y_e; at
    four flavors the matter term then follows the composition (a position-resolved
    projector).  Two flavors missed on 6 of 12 energies here, by up to 2.3e-3, silently.
    The reference is the per-point ladder at 1e-11, with the batched engine disabled."""
    import gen_profile_benchmarks as gpb
    from magnus import earth
    fn = {2: op.osc_prob_2nu_earth, 4: op.osc_prob_4nu_earth}[d]
    energies = np.linspace(2.0, 12.0, 12)*gd.UNIT_GEV
    kw = dict(costhz=-0.9, L=earth.distance_traveled_inside_earth(-0.9)*gd.UNIT_KM,
              **gpb.osc_params(d), **({} if d == 2 else dict(nu_i=gd.NUMU, nu_f=gd.NUMU)))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with op._engine_probe(disabled=('separable',)):
            ref = np.asarray(fn(energies, rtol=1e-11, atol=1e-11, strategy='magnus', **kw))
    tol = 1e-3
    info = {}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        P = np.asarray(fn(energies, rtol=tol, atol=tol, strategy_info=info, **kw))
    assert info['engine'] == 'separable'
    assert not any(issubclass(w.category, op.ToleranceNotAchievedWarning) for w in caught)
    assert np.all(np.abs(P - ref) <= tol + tol*np.abs(ref))


def test_the_grid_point_count_matches_the_grid():
    r"""``_n_grid_points`` counts the grid without building it, including breakpoints that
    fall exactly on a node of the uniform grid, which add no point."""
    rng = np.random.default_rng(71)
    L0, L = 0.0, 1.0
    for _ in range(500):
        n = int(rng.integers(1, 60))
        bp = rng.uniform(L0, L, int(rng.integers(0, 20)))
        on_nodes = np.linspace(L0, L, n + 1)[rng.integers(0, n + 1, int(rng.integers(0, 4)))]
        raw = np.concatenate([bp, on_nodes, bp[:2], [L0, L]])
        bp_in = np.unique(raw)
        bp_in = bp_in[(bp_in > L0) & (bp_in < L)]
        grid = np.unique(np.concatenate([np.linspace(L0, L, n + 1), bp_in]))
        assert op._n_grid_points(L0, L, n, bp_in) == len(grid)
