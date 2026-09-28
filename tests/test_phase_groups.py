# -*- coding: utf-8 -*-
"""Issue #111: the energy-batched engine groups the energies of a scan by their own phase.

With eV-scale sterile splittings the phase across an Earth chord is ~1e4 rad and varies as
1/E, so one slab grid sized for the fastest energy did several times the slab-work of one
ladder per energy, and batching lost to one call per energy.  The engine now splits such a
scan into groups of energies of similar phase, each on its own grid, where its cost model
(BATCHED_PHASE_GROUPING) says that pays; a scan it does not split is computed exactly as
before.
"""
import warnings

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.magnus as mg
import magnus.oscprob as op

L_CHORD = 6371.0*gd.UNIT_KM
STERILE = dict(s14=0.1, s24=0.1, D41=1.0)


@pytest.fixture
def slab_work(monkeypatch):
    """Total slab-energies the batched engine exponentiates (sum over levels of active
    energies times slabs), recorded from the samples it hands to the Magnus kernel."""
    work = []
    real = mg.evolution_operators_from_samples

    def spy(At, widths, *args, **kwargs):
        if At.ndim == 5:                      # (nE, n_slabs, m, d, d): the batched engine
            work.append(At.shape[0]*At.shape[1])
        return real(At, widths, *args, **kwargs)
    monkeypatch.setattr(mg, 'evolution_operators_from_samples', spy)
    return work


def _no_grouping(monkeypatch):
    """Disable the split without touching anything else (a no-op on versions without it)."""
    cost = dict(getattr(op, 'BATCHED_PHASE_GROUPING', {}), margin=np.inf)
    monkeypatch.setattr(op, 'BATCHED_PHASE_GROUPING', cost, raising=False)


def test_sterile_scan_does_less_slab_work_than_one_shared_grid(monkeypatch, slab_work):
    energies = np.linspace(1.0, 40.0, 6)*gd.UNIT_GEV
    call = dict(costhz=-0.5, L=L_CHORD, rtol=1e-3, atol=1e-3, **STERILE)
    info = {}
    P = np.asarray(op.osc_prob_4nu_earth(energies, strategy_info=info, **call))
    assert info['engine'] == 'separable'
    grouped = sum(slab_work)
    slab_work.clear()
    _no_grouping(monkeypatch)
    P_one_grid = np.asarray(op.osc_prob_4nu_earth(energies, **call))
    one_grid = sum(slab_work)
    # The phase spans a factor 40 across the scan; on one grid every energy pays for the
    # fastest.  Measured: 0.19 of the work at 40 energies over 1-40 GeV.
    assert grouped < 0.5*one_grid, (grouped, one_grid)
    # Both are within the requested tolerance of each other and of the per-point ladder.
    with op._engine_probe(disabled=('separable',)):
        P_pt = np.array([np.asarray(op.osc_prob_4nu_earth(e, **call)) for e in energies])
    assert np.allclose(P, P_pt, rtol=1e-3, atol=1e-3)
    assert np.allclose(P_one_grid, P_pt, rtol=1e-3, atol=1e-3)


def test_a_scan_the_model_does_not_split_is_bit_identical(monkeypatch):
    # Three flavors, no steriles: 40 rad of phase, seeds of a few slabs, nothing to split.
    energies = np.linspace(0.5, 5.0, 12)*gd.UNIT_GEV
    P = np.asarray(op.osc_prob_3nu_earth(energies, costhz=-0.5, L=L_CHORD))
    _no_grouping(monkeypatch)
    assert np.array_equal(P, np.asarray(op.osc_prob_3nu_earth(energies, costhz=-0.5, L=L_CHORD)))


def test_phase_groups_split_by_seed_and_only_when_it_pays():
    groups = op._phase_groups
    # Seeds falling as 1/E over a factor 10, four flavors: several groups, each a run of
    # consecutive seeds, covering every energy once, in increasing order of seed.
    seeds = np.ceil(4000.0/np.linspace(0.5, 5.0, 40)*0.5)
    g = groups(seeds, 4, 1.5, 1e-3)
    assert g is not None and len(g) > 1
    idx = np.concatenate(g)
    assert sorted(idx.tolist()) == list(range(40))
    tops = [seeds[x].max() for x in g]
    assert tops == sorted(tops)
    assert all(seeds[x].max() <= seeds[y].min() for x, y in zip(g, g[1:]))
    # Equal seeds: one grid.
    assert groups(np.full(40, 4000.0), 4, 1.5, 1e-3) is None
    # Below the tolerance floor (a group's coarse first levels can agree by chance): one grid.
    assert groups(seeds, 4, 1.5, 1e-6) is None
    # Where the tolerance needs many levels anyway, the split has nothing to gain.
    assert groups(seeds, 4, 1.5, 1e-10) is None
    # Seeds too small for a split to pay its fixed cost, at any flavor count.
    assert groups(np.ceil(40.0/np.linspace(0.5, 5.0, 40)*0.5), 2, 1.5, 1e-3) is None
    assert groups(np.ceil(4.0/np.linspace(0.5, 5.0, 40)*0.5), 5, 1.5, 1e-3) is None
    # Two flavors are cheap per slab, so the same seeds split less than at four.
    g2 = groups(seeds, 2, 1.5, 1e-3)
    assert g2 is None or len(g2) < len(g)


def test_warnings_stay_once_per_call_when_a_scan_is_split():
    energies = np.linspace(1.0, 40.0, 6)*gd.UNIT_GEV
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        op.osc_prob_4nu_earth(energies, costhz=-0.5, L=L_CHORD, rtol=1e-3, atol=1e-3,
                              max_n_slabs=300, **STERILE)
    tol = [w for w in caught if issubclass(w.category, op.ToleranceNotAchievedWarning)]
    assert len(tol) == 1
