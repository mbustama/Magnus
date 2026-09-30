# -*- coding: utf-8 -*-
r"""Issue #71: the energy-batched ``'gl'`` ladder refuses an agreement reached while an energy's
own slabs are still wide (``BATCHED_GL_MAX_SLAB_NORM``), and the per-row slab norms it uses
come from ``magnus._row_slab_norms`` without changing anything else."""
import warnings

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.magnus as mg
import magnus.oscprob as op

P3 = gd.load_nufit_params('NuFIT 6.1')
P2 = dict(sth=P3['s12'], Dm2=P3['D31'])


def _exp_rho(rho0, scale_km):
    return lambda l: rho0*np.exp(-np.asarray(l, dtype=float)/(scale_km*gd.UNIT_KM))


def _scan(d, rho, L_km, E_gev, **kw):
    params = P3 if d == 3 else P2
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return np.asarray(op.osc_prob_matter_std_potential(
            d, rho, np.asarray(E_gev)*gd.UNIT_GEV, L_km*gd.UNIT_KM, params, L0=0.0,
            density_matter_is_in_g_per_cm3=True, strategy='magnus', integration_method='gl',
            **kw))


@pytest.mark.parametrize('d, L_km, scale_km, n_wrong_before', [(2, 250, 100, 15), (3, 250, 100, 3)])
def test_gl_scan_on_an_exponential_profile_is_within_tolerance(d, L_km, scale_km, n_wrong_before):
    """W1/W4 of the #71 pool: before the gate, 15 of 40 (2nu) and 3 of 40 (3nu) energies were
    returned outside rtol=atol=1e-3 (up to 19x) by two coarse levels that agreed by chance."""
    E = np.linspace(0.1, 0.4, 40)
    tol = 1e-3
    P = _scan(d, _exp_rho(3e3, scale_km), L_km, E, rtol=tol, atol=tol)
    R = _scan(d, _exp_rho(3e3, scale_km), L_km, E, rtol=1e-9, atol=1e-9)
    over = (np.abs(P - R) > tol + tol*np.abs(R)).reshape(len(E), -1).any(axis=1)
    assert not over.any(), '%d of %d energies outside the tolerance (was %d)' % (
        over.sum(), len(E), n_wrong_before)


def test_gate_result_does_not_depend_on_the_energy_chunking():
    """The per-row norms of several chunks are concatenated in order."""
    E = np.linspace(0.1, 0.4, 40)
    kw = dict(rtol=1e-3, atol=1e-3)
    P = _scan(3, _exp_rho(3e3, 100), 250, E, **kw)
    old = op.BATCH_WORKING_ENTRIES
    op.BATCH_WORKING_ENTRIES = 1500          # a few energies per chunk
    try:
        Pc = _scan(3, _exp_rho(3e3, 100), 250, E, **kw)
    finally:
        op.BATCH_WORKING_ENTRIES = old
    assert np.array_equal(P, Pc)


def test_row_slab_norms_are_the_per_row_maximum_and_leave_the_warning_unchanged():
    rng = np.random.default_rng(71)
    A = rng.standard_normal((6, 5, 3, 3)) + 1j*rng.standard_normal((6, 5, 3, 3))
    Om = 0.4*(A - np.conj(np.swapaxes(A, -1, -2)))
    with mg._deferred_slab_norm() as seen0:
        U0 = mg._expm_stack(Om, warn_wide=True)
    with mg._deferred_slab_norm() as seen1:
        with mg._row_slab_norms() as rows:
            U1 = mg._expm_stack(Om, warn_wide=True)
    assert mg._ROW_NORM_SINK.get() is None
    assert np.array_equal(U0, U1)
    assert seen0 == seen1                     # MagnusConvergenceWarning sees the same value
    assert len(rows) == 1 and rows[0].shape == (6,)
    lam = np.linalg.eigvalsh(1j*Om)
    expect = np.max(np.abs(lam), axis=(-1, -2))
    assert np.allclose(rows[0], expect, rtol=1e-12, atol=0.0)
    # The gate keeps the full norm; the warning measures the traceless part (issue #155 §1).
    traceless = np.max(np.abs(lam - lam.mean(axis=-1, keepdims=True)))
    assert np.isclose(seen1[0], traceless, rtol=1e-12, atol=0.0)
    assert seen1[0] <= float(rows[0].max()) + 1e-12
    # nothing is collected for a constant A, without warn_wide, and the sink is restored
    # after an exception
    with mg._row_slab_norms() as rows:
        mg._expm_stack(Om, warn_wide=True, A_is_const=True)
        mg._expm_stack(Om, warn_wide=False)
    assert rows == []
    with pytest.raises(RuntimeError):
        with mg._row_slab_norms():
            raise RuntimeError
    assert mg._ROW_NORM_SINK.get() is None
