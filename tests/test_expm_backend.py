# -*- coding: utf-8 -*-
r"""Tests for the matrix-exponential backends and the switch that selects them.

``EXPM_BACKEND`` chooses between the compiled Cayley-Hamilton kernel of
:mod:`magnus.expmkernels` and ``numpy.linalg.eigh``.  Four things are worth
stating about what is tested here, because each is a way this could have shipped
broken:

* **Assert that the kernel fired.**  Comparing the two backends proves nothing
  if the fast path silently declined -- both sides then run ``eigh`` and agree to
  0.0.  Every equivalence test below goes through :func:`_kernel_calls`, which
  counts entries into the kernel and asserts the count is non-zero.  Exactly
  this trap was hit before, by the palindromic-profile tests; see
  ``test_palindrome.py``.

* **Degeneracy is the risk, and a plausible wrong answer is the failure mode.**
  An earlier prototype of this kernel put the conjugate of
  :math:`\det X`'s cross term on the wrong factor, ran 6x faster, and returned
  matrices that were still nearly unitary and wrong by O(1).
  :func:`test_determinant_cross_term_conjugate` pins that term against
  ``np.linalg.det``, and both misplacements are wrong by O(1) rather than
  subtly.

* **Eigenvalue accuracy and exponential accuracy are decoupled.**  The
  closed-form eigenvalues degrade to ~1e-9 at a degeneracy, because
  :math:`\arccos` has infinite derivative where a repeated root sits, and
  :math:`\exp(-iK)` stays accurate to 1e-16 there anyway, because interpolation
  error is *second* order in the displacement of a coalescing node.  Both halves
  are asserted -- the sloppy eigenvalues included -- so that a future change
  which tightens one at the other's expense is visible.

* **The two backends must read the same triangle.**  ``eigh`` reads the lower
  one (``UPLO`` defaults to ``'L'``), and :func:`magnus.magnus._expm_stack`
  admits input that is anti-Hermitian only to 1e-12 relative.  A kernel reading
  the upper triangle therefore exponentiates a *different matrix* on such input,
  and the backends diverge by ~2e-12 -- which is what the first version of the
  kernel did.  :func:`test_backends_agree_on_near_hermitian_input` is what keeps
  them fed identically.
"""

import warnings

import numpy as np
import pytest

import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.oscprob as op
from magnus import expmkernels as ek
from magnus import magnus as mg


requires_numba = pytest.mark.skipif(not ek.HAVE_NUMBA,
                                    reason="numba is not installed")


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

@pytest.fixture
def restore_backend():
    r"""Puts ``EXPM_BACKEND`` back, whatever a test does to it."""
    saved = mg.EXPM_BACKEND
    yield
    mg.EXPM_BACKEND = saved


class _kernel_calls:
    r"""Context manager counting entries into the compiled kernel.

    ``with _kernel_calls() as n: ...; assert n.count`` is how every equivalence
    test here establishes that the fast path actually ran.  Without it, a
    declined kernel makes both sides of the comparison identical and the test
    passes by running nothing.
    """

    def __init__(self):
        self.count = 0
        self.matrices = 0

    def __enter__(self):
        self._saved = ek.expm_herm_stack

        def counted(K):
            self.count += 1
            self.matrices += int(np.prod(K.shape[:-2])) if K.ndim > 2 else 1
            return self._saved(K)

        ek.expm_herm_stack = counted
        return self

    def __exit__(self, *exc):
        ek.expm_herm_stack = self._saved
        return False


def _herm(shape, rng, scale=1.0):
    r"""A random Hermitian stack, exactly Hermitian (bitwise symmetric)."""
    M = rng.normal(size=shape) + 1j*rng.normal(size=shape)
    H = (M + np.conj(np.swapaxes(M, -1, -2)))/2.0
    return scale*H


def _scipy_expm(K):
    r"""Reference :math:`\exp(-iK)`, matrix by matrix, from scipy."""
    import scipy.linalg as sla
    d = K.shape[-1]
    flat = np.asarray(K).reshape(-1, d, d)
    return np.array([sla.expm(-1j*w) for w in flat]).reshape(K.shape)


def _unitarity(U):
    d = U.shape[-1]
    Uh = np.conj(np.swapaxes(U, -1, -2))
    return np.max(np.abs(Uh @ U - np.eye(d)))


# --------------------------------------------------------------------------
# the switch itself
# --------------------------------------------------------------------------

def test_default_backend_is_auto():
    r"""'auto' is the default, and it is the value that cannot fail."""
    assert mg.EXPM_BACKEND == 'auto'
    assert mg.valid_expm_backends == ['auto', 'numba', 'eigh']


@pytest.mark.parametrize('bad', ['numpy', 'Numba', 'EIGH', '', 'scipy', 0])
def test_invalid_backend_name_raises(bad):
    with pytest.raises(ValueError, match='expm_backend'):
        mg._expm_stack(-1j*np.eye(3), expm_backend=bad)


def test_explicit_numba_without_numba_raises(monkeypatch):
    r"""Asked for by name, a missing numba is an error, not a silent downgrade.

    The caller who writes ``'numba'`` wanted to know the compiled path was
    running; answering with ``eigh`` and saying nothing would defeat the point of
    naming it.  ``'auto'`` is the value that promises to work anywhere.
    """
    monkeypatch.setattr(ek, 'HAVE_NUMBA', False)
    with pytest.raises(ValueError, match='numba is not installed'):
        mg._expm_stack(-1j*np.eye(3), expm_backend='numba')


def test_auto_falls_back_without_numba(monkeypatch):
    r"""'auto' still answers, correctly, on an install with no numba."""
    monkeypatch.setattr(ek, 'HAVE_NUMBA', False)
    U = mg._expm_stack(-1j*0.3*np.eye(3), expm_backend='auto')
    assert np.max(np.abs(U - np.exp(-0.3j)*np.eye(3))) < 1e-15


@requires_numba
def test_eigh_backend_does_not_call_the_kernel():
    r"""The opt-out has to actually opt out, or comparisons against it are void."""
    with _kernel_calls() as n:
        mg._expm_stack(-1j*_herm((16, 3, 3), np.random.default_rng(0)),
                       expm_backend='eigh')
    assert n.count == 0


@requires_numba
@pytest.mark.parametrize('d', [4, 5])
def test_dim_four_and_five_delegate_to_eigh(d):
    r"""No closed form for a 4x4/5x5 Hermitian eigenproblem, so 4nu/5nu use eigh.

    They stay correct; they are simply not accelerated.  Asserted rather than
    assumed, because the alternative -- a kernel that quietly handled dimension 4
    with the 3x3 formulas -- would be wrong by O(1).
    """
    assert not ek.supports_dim(d)
    K = _herm((8, d, d), np.random.default_rng(d))
    with _kernel_calls() as n:
        U = mg._expm_stack(-1j*K, expm_backend='numba')
    assert n.count == 0
    assert np.max(np.abs(U - _scipy_expm(K))) < 1e-13


@requires_numba
def test_module_switch_reaches_oscprob(restore_backend):
    r"""``EXPM_BACKEND`` is the control that reaches the physics wrappers."""
    mg.EXPM_BACKEND = 'eigh'
    with _kernel_calls() as n:
        op.osc_prob_3nu_vacuum(1.0*gd.UNIT_GEV, 1300.0*gd.UNIT_KM,
                               validate_input=False)
    assert n.count == 0
    mg.EXPM_BACKEND = 'auto'
    with _kernel_calls() as n:
        op.osc_prob_3nu_vacuum(1.0*gd.UNIT_GEV, 1300.0*gd.UNIT_KM,
                               validate_input=False)
    assert n.count > 0


@requires_numba
def test_per_call_parameter_overrides_the_module_switch(restore_backend):
    mg.EXPM_BACKEND = 'eigh'
    K = _herm((4, 3, 3), np.random.default_rng(1))
    with _kernel_calls() as n:
        mg._expm_stack(-1j*K, expm_backend='auto')
    assert n.count == 1


@requires_numba
def test_expm_backend_threads_through_the_public_entry_points():
    r"""The three magnus.py entry points pass the parameter down."""
    H = _herm((3, 3), np.random.default_rng(2), scale=1e-13)
    edges = np.stack([np.linspace(0, 1e13, 6)[:-1],
                      np.linspace(0, 1e13, 6)[1:]], axis=1)
    calls = [
        lambda b: mg.magnus_expansion(lambda t: -1j*H, 0.0, 1e13, order=4,
                                      expm_backend=b),
        lambda b: mg.magnus_expansion_multislab(lambda t: -1j*H, edges,
                                                order=4, expm_backend=b),
        lambda b: mg.evolution_operators_from_samples(
            np.broadcast_to(-1j*H, (5, 2, 3, 3)),
            np.full(5, 2e12), order=4, expm_backend=b),
    ]
    for call in calls:
        with _kernel_calls() as n:
            got = call('numba')
        assert n.count > 0
        with _kernel_calls() as n:
            ref = call('eigh')
        assert n.count == 0
        assert np.max(np.abs(np.asarray(got) - np.asarray(ref))) < 1e-13


# --------------------------------------------------------------------------
# the algebra: degeneracies, which are the whole risk
# --------------------------------------------------------------------------

@requires_numba
@pytest.mark.parametrize('name,K', [
    ('zero matrix, d=3', np.zeros((3, 3))),
    ('zero matrix, d=2', np.zeros((2, 2))),
    ('H = I, three-fold degenerate', 1.7*np.eye(3)),
    ('H = -I', -1.7*np.eye(3)),
    ('diag(1, 1, 2), two-fold', np.diag([1.0, 1.0, 2.0])),
    ('diag(1, 2, 2), two-fold', np.diag([1.0, 2.0, 2.0])),
    ('diag(0, 0, 1)', np.diag([0.0, 0.0, 1.0])),
    ('H = I, d=2', 0.9*np.eye(2)),
    ('diag(1, 1), d=2', np.diag([1.0, 1.0])),
])
def test_exact_degeneracies(name, K):
    r"""The cases an eigenvector route died on, and this one must not.

    A coincident eigenvalue makes an eigenvector null space multi-dimensional and
    every row cross-product vanish; the prototype that preceded this kernel
    raised ``ZeroDivisionError`` on the very first of these.  Cayley-Hamilton has
    no eigenvectors to extract, and needs no confluent (Hermite) form either: a
    Hermitian matrix is never defective, so matching :math:`\exp` on the
    *distinct* eigenvalues is already exact.
    """
    K = np.asarray(K, dtype=complex)
    U, lam, _sev = ek.expm_herm_stack(K[None])
    assert np.all(np.isfinite(U)), name
    assert np.max(np.abs(U[0] - _scipy_expm(K))) < 1e-14, name
    assert _unitarity(U[0]) < 1e-14, name
    # 1e-7, not 1e-14, and deliberately: an exactly two-fold degenerate spectrum
    # is where the closed-form eigenvalues are at their worst, and diag(1, 1, 2)
    # comes out 4.1e-09 off LAPACK.  The exponential above is unaffected at
    # 2.5e-16, which is the whole point of
    # test_eigenvalues_are_sloppy_at_a_degeneracy_and_the_exponential_is_not.
    # The eigenvalues leave this module only for the slab-width convergence
    # warning, which compares them against a threshold of order 1, so 1e-9 there
    # is immaterial.  Tightening this assertion would be asserting something the
    # scheme does not provide and does not need.
    assert np.max(np.abs(np.sort(lam[0]) - np.linalg.eigvalsh(K))) < 1e-7, name


@requires_numba
@pytest.mark.parametrize('split', [1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12,
                                   1e-14, 1e-16, 0.0])
def test_near_degeneracy_sweep_is_the_msw_resonance_case(split):
    r"""Two eigenvalues genuinely approaching each other, in a generic basis.

    This is the shape of an MSW resonance, where a level crossing is avoided by a
    gap that shrinks but does not close, and it is the reason resonance
    robustness was asked for.  The error must be flat across the sweep: a scheme
    with a tolerance-based near-degenerate branch has a notch at its crossover,
    and this one has no branch to place.

    Both the exactly-zero splitting and the 1e-16 one are included; the first
    takes the guarded path, the second does not, and they must agree.
    """
    rng = np.random.default_rng(31)
    V = np.linalg.qr(_herm((3, 3), rng) + 1j*np.eye(3))[0]
    K = V @ np.diag([0.3, 0.3 + split, 1.1]).astype(complex) @ V.conj().T
    K = (K + K.conj().T)/2.0
    U, _, _sev = ek.expm_herm_stack(K[None])
    assert np.max(np.abs(U[0] - _scipy_expm(K))) < 1e-14
    assert _unitarity(U[0]) < 1e-14
    assert np.max(np.abs(U[0] - mg._expm_stack(-1j*K, expm_backend='eigh'))) < 1e-13


@requires_numba
@pytest.mark.parametrize('split', [1e-2, 1e-5, 1e-8, 1e-11, 1e-14, 0.0])
def test_all_three_eigenvalues_near_degenerate(split):
    r"""The case that drives the second divided difference's denominator to zero."""
    rng = np.random.default_rng(32)
    V = np.linalg.qr(_herm((3, 3), rng) + 1j*np.eye(3))[0]
    K = V @ np.diag([0.7, 0.7 + split, 0.7 + 2*split]).astype(complex) @ V.conj().T
    K = (K + K.conj().T)/2.0
    U, _, _sev = ek.expm_herm_stack(K[None])
    assert np.max(np.abs(U[0] - _scipy_expm(K))) < 1e-14
    assert _unitarity(U[0]) < 1e-14


@requires_numba
def test_eigenvalues_are_sloppy_at_a_degeneracy_and_the_exponential_is_not():
    r"""Pins the decoupling the kernel's accuracy actually rests on.

    The closed-form eigenvalues lose most of their digits as two roots coalesce
    -- :math:`\arccos` has infinite derivative at the ends of its range, which is
    exactly where a repeated root sits -- and that does *not* propagate to
    :math:`\exp(-iK)`, because the interpolation error at a node displaced by
    :math:`\delta` from a double eigenvalue is :math:`O(\delta^2)`, not
    :math:`O(\delta)`.

    Both halves are asserted, the bad one included: if a later change made the
    eigenvalues sharp, this test should be re-derived rather than merely relaxed,
    because the second assertion is the one carrying the accuracy claim.
    """
    rng = np.random.default_rng(33)
    V = np.linalg.qr(_herm((3, 3), rng) + 1j*np.eye(3))[0]
    K = V @ np.diag([0.3, 0.3, 1.1]).astype(complex) @ V.conj().T
    K = (K + K.conj().T)/2.0
    U, lam, _sev = ek.expm_herm_stack(K[None])
    eig_err = np.max(np.abs(np.sort(lam[0]) - np.linalg.eigvalsh(K)))
    exp_err = np.max(np.abs(U[0] - _scipy_expm(K)))
    assert eig_err > 1e-12, (
        "eigenvalues came out sharper than expected (%.2e); the decoupling "
        "argument in expmkernels' docstring should be revisited" % eig_err)
    assert exp_err < 1e-15, exp_err


@requires_numba
def test_determinant_cross_term_conjugate():
    r"""The exact term an earlier prototype got wrong, pinned against numpy.

    ``2*Re(p01*p12*conj(p02))``.  Moving the conjugate to either of the other two
    factors leaves a determinant that is wrong by O(1) -- and, downstream, an
    exponential that is still nearly unitary, still plausible, and wrong by O(1).
    That is why this is asserted directly on the determinant rather than only
    through the exponential.
    """
    rng = np.random.default_rng(34)
    worst_ok = 0.0
    worst_bad = 0.0
    for _ in range(400):
        X = _herm((3, 3), rng)
        X = X - np.trace(X).real/3.0*np.eye(3)          # traceless part
        a, b, c = X[0, 0].real, X[1, 1].real, X[2, 2].real
        p01, p02, p12 = X[0, 1], X[0, 2], X[1, 2]
        n01, n02, n12 = abs(p01)**2, abs(p02)**2, abs(p12)**2
        base = a*b*c - a*n12 - b*n02 - c*n01
        good = base + 2.0*(p01*p12*np.conj(p02)).real
        ref = np.linalg.det(X).real
        worst_ok = max(worst_ok, abs(good - ref))
        for bad in (base + 2.0*(p01*np.conj(p12)*p02).real,
                    base + 2.0*(np.conj(p01)*p12*p02).real):
            worst_bad = max(worst_bad, abs(bad - ref))
    assert worst_ok < 1e-12, worst_ok
    assert worst_bad > 1e-2, (
        "the misplaced-conjugate variants were not distinguishable here, so "
        "this test would not have caught the historical bug")


@requires_numba
def test_eigenvalues_are_returned_ascending():
    r"""The kernel stores its roots in a known order instead of sorting them.

    Asserted, not enforced, so a slip in the trigonometric formula surfaces as a
    failure rather than being tidied away by a sort.  Includes exact ties and
    trace offsets, where the ordering argument is tightest.
    """
    rng = np.random.default_rng(35)
    stacks = [_herm((2000, 3, 3), rng),
              _herm((2000, 3, 3), rng, scale=1e-8),
              _herm((200, 3, 3), rng) + 1e6*np.eye(3),
              np.broadcast_to(np.diag([1.0, 1.0, 2.0]).astype(complex),
                              (4, 3, 3)),
              np.zeros((4, 3, 3), dtype=complex),
              _herm((2000, 2, 2), rng)]
    for K in stacks:
        K = np.ascontiguousarray(K, dtype=complex)
        _, lam, _sev = ek.expm_herm_stack(K)
        assert np.all(np.diff(lam, axis=-1) >= 0.0)
        assert np.max(np.abs(lam - np.linalg.eigvalsh(K))) < 1e-7


# --------------------------------------------------------------------------
# backend equivalence
# --------------------------------------------------------------------------

@requires_numba
@pytest.mark.parametrize('d', [2, 3, 4, 5])
@pytest.mark.parametrize('n', [1, 2, 7, 108, 4096])
def test_backends_agree_on_random_hermitian_stacks(d, n):
    rng = np.random.default_rng(1000*d + n)
    Om = -1j*_herm((n, d, d), rng)
    with _kernel_calls() as spy:
        Uk = mg._expm_stack(Om, expm_backend='numba')
    Ue = mg._expm_stack(Om, expm_backend='eigh')
    if ek.supports_dim(d):
        assert spy.count == 1, "the kernel declined; this comparison is void"
    assert np.max(np.abs(Uk - Ue)) < 1e-14
    assert _unitarity(Uk) < 1e-13
    assert Uk.shape == Om.shape


@requires_numba
@pytest.mark.parametrize('d', [2, 3])
def test_backends_agree_on_near_hermitian_input(d):
    r"""Both routes must read the *same* triangle of a not-quite-Hermitian input.

    :func:`magnus.magnus._expm_stack` accepts input that is anti-Hermitian only
    to 1e-12 relative, and ``eigh`` reads the lower triangle (``UPLO='L'``).  A
    kernel reading the upper one exponentiates a different matrix on such input
    and diverges by ~2e-12 -- large enough to matter, small enough to be mistaken
    for ordinary rounding, which is what makes it worth a test of its own.
    """
    rng = np.random.default_rng(2000 + d)
    worst = 0.0
    for _ in range(200):
        M = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d))
        Om = (M - M.conj().T)/2.0
        # A *Hermitian* perturbation, so the whole of it is anti-Hermiticity
        # violation and its size is known rather than merely bounded: for
        # anti-Hermitian Om, ``Om + Om^H`` is exactly ``2*eps*D``.  Sized to land
        # just inside the 1e-12 gate -- overshoot it and _expm_stack routes to the
        # scipy fallback instead, both backends run identical code, and the test
        # passes while proving nothing.  That is what a first version of this did.
        D = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d))
        D = (D + D.conj().T)/2.0
        D = D/np.max(np.abs(D))
        Om = Om + (0.25e-12*np.max(np.abs(Om)))*D
        with _kernel_calls() as spy:
            Uk = mg._expm_stack(Om, expm_backend='numba')
        assert spy.count == 1, (
            "input fell through the anti-Hermiticity gate to the scipy "
            "fallback, so this comparison is void")
        Ue = mg._expm_stack(Om, expm_backend='eigh')
        worst = max(worst, np.max(np.abs(Uk - Ue)))
    assert worst < 1e-14, (
        "backends diverged by %.2e on input at the Hermiticity tolerance; "
        "check that the kernel reads the lower triangle, as eigh does" % worst)


@requires_numba
@pytest.mark.parametrize('scale', [1e-150, 1e-8, 1.0, 10.0, 100.0, 1e3, 1e4])
def test_kernel_is_no_worse_than_eigh_at_any_norm(scale):
    r"""Both routes lose digits linearly in :math:`\lVert K \rVert`; neither cliffs.

    Magnus slabs normally satisfy :math:`\lVert\Omega\rVert \lesssim 1`, but a
    constant Hamiltonian over a long baseline is exponentiated in one slab, and
    then the norm is the whole accumulated phase and reaches thousands.  The
    kernel must not be *worse* than the backend it replaces there.
    """
    rng = np.random.default_rng(int(-np.log10(scale)) + 50)
    K = _herm((32, 3, 3), rng)
    K *= scale/np.max(np.abs(K))
    ref = _scipy_expm(K)
    err_k = np.max(np.abs(ek.expm_herm_stack(np.ascontiguousarray(K))[0] - ref))
    err_e = np.max(np.abs(mg._expm_stack(-1j*K, expm_backend='eigh') - ref))
    tol = max(1e-15, 1e-14*scale)
    assert err_k < tol, (scale, err_k)
    assert err_k <= 10.0*max(err_e, 1e-17), (
        "kernel error %.2e against eigh's %.2e at ||K||=%.0e" % (err_k, err_e, scale))


@requires_numba
def test_compiled_kernel_matches_the_same_source_uncompiled():
    r"""Guards against numba miscompiling the kernel rather than the maths being wrong.

    The two are expected to agree *bitwise*: ``fastmath`` is off precisely so
    that the compiler may not reassociate the cancellations the stability
    argument depends on, and a difference here would mean it did.
    """
    rng = np.random.default_rng(36)
    for d, core, kern in ((2, ek._ch2_core, ek._ch2_kernel),
                          (3, ek._ch3_core, ek._ch3_kernel)):
        K = np.ascontiguousarray(_herm((128, d, d), rng), dtype=complex)
        o_py, l_py = np.empty_like(K), np.empty((128, d))
        o_jit, l_jit = np.empty_like(K), np.empty((128, d))
        core(K, o_py, l_py)
        kern(K, o_jit, l_jit)
        assert np.array_equal(o_py, o_jit), d
        assert np.array_equal(l_py, l_jit), d


@requires_numba
@pytest.mark.parametrize('name,K', [
    ('zero matrix, d=3', np.zeros((3, 3))),
    ('H = I, three-fold degenerate', 1.7*np.eye(3)),
    ('diag(1, 1, 2), two-fold', np.diag([1.0, 1.0, 2.0])),
    ('clamp at u = -1', np.diag([-2.0, 1.0, 1.0])),
    ('zero matrix, d=2', np.zeros((2, 2))),
    ('H = I, d=2', 0.9*np.eye(2)),
])
def test_uncompiled_core_handles_degeneracies_too(name, K):
    r"""The *un-jitted* source, on the degenerate inputs, not only on random ones.

    Written because of what the coverage report showed rather than what it said.
    The guards inside :func:`_ch3_core` -- the ``m <= 0`` branch for X = 0 and the
    two ``arccos`` clamps -- came back uncovered even though
    :func:`test_exact_degeneracies` exercises exactly those cases, because that
    test goes through the compiled kernel and coverage cannot see lines executed
    inside numba.  So the branches were tested and the *Python* form of them was
    not, and the Python form is what an install without numba runs.  The only
    other test to call the cores directly uses random Hermitian matrices, which
    are never degenerate.
    """
    K = np.ascontiguousarray(np.asarray(K, dtype=complex)[None])
    d = K.shape[-1]
    core = ek._ch2_core if d == 2 else ek._ch3_core
    kern = ek._ch2_kernel if d == 2 else ek._ch3_kernel
    o_py, l_py = np.empty_like(K), np.empty((1, d))
    o_jit, l_jit = np.empty_like(K), np.empty((1, d))
    core(K, o_py, l_py)
    kern(K, o_jit, l_jit)
    assert np.all(np.isfinite(o_py)), name
    assert np.array_equal(o_py, o_jit), name
    assert np.array_equal(l_py, l_jit), name
    assert np.max(np.abs(o_py[0] - _scipy_expm(K[0]))) < 1e-14, name


@requires_numba
def test_sinc_compiled_matches_python_including_zero():
    for x in [0.0, 1e-300, 1e-8, 1e-4 - 1e-18, 1e-4, 1e-4 + 1e-18, 0.5,
              1.0, np.pi, -np.pi, 10.0, -1e-5]:
        assert ek._sinc(x) == ek._sinc_py(x), x
    assert ek._sinc_py(0.0) == 1.0


@requires_numba
def test_real_valued_hermitian_input_is_accepted():
    r"""A real Hermitian K is a valid input; the output buffer must still be complex."""
    K = np.array([[[1.0, 0.5], [0.5, 2.0]]])
    U, lam, _sev = ek.expm_herm_stack(K)
    assert np.max(np.abs(U[0] - _scipy_expm(K[0]))) < 1e-14
    assert np.max(np.abs(lam[0] - np.linalg.eigvalsh(K[0]))) < 1e-14


@requires_numba
@pytest.mark.parametrize('shape', [(3, 3), (1, 3, 3), (4, 5, 3, 3),
                                   (2, 3, 4, 2, 2)])
def test_leading_axes_are_preserved(shape):
    r"""oscprob batches energy in front of the slab axis, so shape must survive."""
    rng = np.random.default_rng(37)
    Om = -1j*_herm(shape, rng)
    Uk = mg._expm_stack(Om, expm_backend='numba')
    Ue = mg._expm_stack(Om, expm_backend='eigh')
    assert Uk.shape == shape
    assert np.max(np.abs(Uk - Ue)) < 1e-14


@requires_numba
def test_non_contiguous_input_is_handled():
    r"""A sliced or transposed stack must not be read through the wrong strides."""
    rng = np.random.default_rng(38)
    big = -1j*_herm((16, 3, 3), rng)
    view = big[::3]
    assert not view.flags['C_CONTIGUOUS']
    Uk = mg._expm_stack(view, expm_backend='numba')
    Ue = mg._expm_stack(np.ascontiguousarray(view), expm_backend='eigh')
    assert np.max(np.abs(Uk - Ue)) < 1e-14


@requires_numba
def test_slab_norm_warning_still_fires_on_the_kernel_path():
    r"""The convergence warning reads the eigenvalues, which both routes return."""
    K = 40.0*np.eye(3)
    with pytest.warns(mg.MagnusConvergenceWarning):
        mg._expm_stack(-1j*K, warn_wide=True, expm_backend='numba')


@requires_numba
def test_non_anti_hermitian_input_still_reaches_scipy():
    r"""The general fallback is untouched by the backend switch."""
    M = np.array([[1.0, 2.0], [0.0, 0.5]], dtype=complex)   # not anti-Hermitian
    U = mg._expm_stack(M, expm_backend='numba')
    import scipy.linalg as sla
    assert np.max(np.abs(U - sla.expm(M))) < 1e-13


# --------------------------------------------------------------------------
# physics level
# --------------------------------------------------------------------------

def _physics_cases():
    E = 1.0*gd.UNIT_GEV
    L = 1300.0*gd.UNIT_KM
    return {
        '3nu vacuum': lambda: op.osc_prob_3nu_vacuum(E, L, validate_input=False),
        '3nu matter constant': lambda: op.osc_prob_3nu_matter_constant_density(
            E, L, 3.0*gd.UNIT_G_PER_CM3, validate_input=False),
        '3nu earth chord': lambda: op.osc_prob_3nu_earth(
            E, costhz=-1.0, L=earth.distance_traveled_inside_earth(-1.0),
            validate_input=False),
        '3nu earth energy scan': lambda: op.osc_prob_3nu_earth(
            np.logspace(-1, 1, 25)*gd.UNIT_GEV, costhz=-0.85,
            L=earth.distance_traveled_inside_earth(-0.85), validate_input=False),
        # strategy='magnus' deliberately: at the default 'auto' the Sun is
        # answered by the adiabatic engine, which never reaches _expm_stack, so
        # the case would compare two identical runs of code this test is not
        # about.  See test_adiabatic_solar_path_is_unaffected for that half.
        '3nu sun': lambda: op.osc_prob_3nu_sun(
            10.0*gd.UNIT_MEV, 0.3*gd.SUN_RADIUS*gd.UNIT_KM, 0.0,
            strategy='magnus', validate_input=False),
        '3nu NSI resonance': lambda: op.osc_prob_3nu_matter_nsi_constant_density(
            E, L, 5.0*gd.UNIT_G_PER_CM3, eps_ee=0.1, eps_em=0.05j, eps_et=0.0,
            eps_mm=0.0, eps_mt=0.0, eps_tt=0.0, validate_input=False),
        '3nu earth NSI': lambda: op.osc_prob_3nu_earth_nsi(
            E, costhz=-0.7, L=earth.distance_traveled_inside_earth(-0.7),
            eps_ee=0.3, eps_em=0.1j, eps_et=0.0, eps_mm=0.0, eps_mt=0.0,
            eps_tt=0.0, validate_input=False),
        '2nu matter constant': lambda: op.osc_prob_2nu_matter_constant_density(
            E, L, 3.0*gd.UNIT_G_PER_CM3, 0.55, 2.5e-3, validate_input=False),
        '4nu vacuum': lambda: op.osc_prob_4nu_vacuum(
            E, L, s14=0.1, d14=0.0, s24=0.1, d24=0.0, s34=0.1, D41=1.0,
            validate_input=False),
    }


# Per-case tolerance, because the cases differ by four orders of magnitude in how
# many exponentials they chain, and a single flat number would either be vacuous
# for the short ones or wrong for the long one.
#
# The solar call at strategy='magnus' chains 33,575 slab exponentials.  Each
# agrees between backends to ~1e-16, and an ordered product of N of them drifts by
# up to N*eps = 7.4e-12; the observed gap is 3.0e-12, i.e. within that bound and
# not a defect in either backend.  Every other case here chains fewer than a
# hundred and holds 1e-12 comfortably.
_PHYSICS_TOL = {'3nu sun': 1e-10}


@requires_numba
@pytest.mark.parametrize('name', list(_physics_cases()))
def test_probabilities_agree_between_backends(name, restore_backend):
    r"""The claim a user cares about: switching backend does not move the physics."""
    fn = _physics_cases()[name]
    tol = _PHYSICS_TOL.get(name, 1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        mg.EXPM_BACKEND = 'eigh'
        Pe = np.asarray(fn(), dtype=float)
        mg.EXPM_BACKEND = 'auto'
        with _kernel_calls() as spy:
            Pk = np.asarray(fn(), dtype=float)
    if name != '4nu vacuum':
        assert spy.count > 0, "kernel declined; this comparison is void"
    assert np.max(np.abs(Pk - Pe)) < tol, name
    assert np.max(np.abs(np.sum(Pk, axis=-1) - 1.0)) < tol, name


@requires_numba
def test_adiabatic_solar_path_is_unaffected(restore_backend):
    r"""The engine that does not exponentiate must not notice the switch at all.

    At the default strategy a solar profile is answered adiabatically, without
    reaching :func:`magnus.magnus._expm_stack`.  The two backends must then agree
    *exactly*, not merely closely: any difference would mean the switch had
    reached code that has no business depending on it.
    """
    def run():
        return np.asarray(op.osc_prob_3nu_sun(
            10.0*gd.UNIT_MEV, 0.3*gd.SUN_RADIUS*gd.UNIT_KM, 0.0,
            validate_input=False), dtype=float)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        mg.EXPM_BACKEND = 'eigh'
        Pe = run()
        mg.EXPM_BACKEND = 'auto'
        with _kernel_calls() as spy:
            Pk = run()
    assert spy.count == 0, "the solar default path reached the kernel after all"
    assert np.array_equal(Pk, Pe)


# ----------------------------------------------------------------------
# The separation-by-scale grid.
#
# This is the test whose absence let a false claim ship.  The kernel was
# verified against random spectra at many norms, and separately at many
# eigenvalue separations at norm ~1 -- and was up to 7440x worse than eigh
# where those two conditions hold *together*, which no single-axis sweep
# reaches.  The CHANGELOG and implementation_details.rst both asserted
# "the same order or slightly better at every norm from 1 to 1e5" on that
# evidence.
#
# arccos has infinite derivative at u = +/-1, so a clustered spectrum turns
# rounding in u into an eigenvalue error ~sqrt(eps)*||K||, which only matters
# once ||K|| is large.  magnus.expmkernels.SEV_TOL gates on the scale and
# hands the rest to eigh; these tests pin that the gate is in the right place
# and that the hole stays shut.
#
# No high-precision oracle is needed: eigh IS the reference the kernel has to
# match, so comparing the two directly is both sufficient and free.  Before the
# gate the two differed by 2.7e-07 on the worst cell; after it, by 1e-13.
# ----------------------------------------------------------------------

_GRID_SHAPES = {
    # A double root at the bottom of the spectrum drives u -> +1, at the top
    # u -> -1, and the two behave differently: at u = +1 the coincident pair
    # comes out of cos(-2pi/3) and cos(+2pi/3), bit-identical by symmetry, so
    # exact degeneracy is harmless there and not at the other end.
    'double-low': lambda S, d: [0.0, d, S],
    'double-high': lambda S, d: [0.0, S - d, S],
    'triple': lambda S, d: [0.0, d, 2.0*d],
    'unclustered': lambda S, d: [-0.5*S, 0.17*S, S],
}
_GRID_SCALES = [1.0, 1e1, 1e2, 1e3, 1e4, 1e5]
# Spanning exact degeneracy, the sqrt(eps) danger zone, and full separation.
_GRID_SEPS = [0.0, 1e-16, 1e-12, 1e-9, 1.49e-8, 1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]


def _grid_matrix(spec, rng):
    r"""A Hermitian 3x3 with exactly the requested spectrum, in a generic basis."""
    Q = np.linalg.qr(rng.normal(size=(3, 3)) + 1j*rng.normal(size=(3, 3)))[0]
    K = Q @ np.diag(np.asarray(spec, dtype=complex)) @ Q.conj().T
    return np.ascontiguousarray((K + K.conj().T)/2.0)


@requires_numba
@pytest.mark.parametrize('shape', list(_GRID_SHAPES))
def test_kernel_tracks_eigh_across_separation_and_scale(shape):
    r"""The two-condition grid: clustered spectrum crossed with large norm.

    Neither axis alone finds anything -- that is the point.  A sweep over norms
    at generic separation gives ratios of 0.4-2.5, and a sweep over separations
    at norm 1 gives 0.4-1.0.  Only the product was bad.
    """
    mk = _GRID_SHAPES[shape]
    worst = 0.0
    where = None
    for S in _GRID_SCALES:
        for f in _GRID_SEPS:
            K = _grid_matrix(mk(S, f*S), np.random.default_rng(20260809))
            Om = -1j*K
            d = np.max(np.abs(mg._expm_stack(Om, expm_backend='numba')
                              - mg._expm_stack(Om, expm_backend='eigh')))
            if d > worst:
                worst, where = d, (S, f)
    assert worst < 1e-12, (
        "%s: backends differ by %.3e at scale %.0e, separation %.0e*scale. "
        "Before magnus.expmkernels.SEV_TOL existed this reached 2.7e-07; if it "
        "has regressed, check that the gate still sends large-norm 3x3 matrices "
        "to eigh." % (shape, worst, where[0], where[1]))


@requires_numba
def test_expm_stack_is_accurate_in_absolute_terms_across_the_grid():
    r"""Absolute error, not the ratio to eigh -- the ratio is the wrong yardstick.

    A ratio of 19.6x at 8.7e-14 absolute is harmless and a ratio of 7440x at
    6.7e-08 is not, and a ratio test cannot tell them apart.  This asserts the
    quantity that matters.

    Note this goes through :func:`magnus.magnus._expm_stack`, not the raw kernel.
    That is the contract, not an evasion: ``expm_herm_stack`` is the bare closed
    form and *returns* its conditioning severity so the caller can decide, which
    is why the gate lives in the caller.  Written against the raw kernel this
    same assertion fails at 8.7e-08 -- see the companion test below, which pins
    exactly that.
    """
    # The tolerance has to grow with the norm, because every method's error does:
    # exp(-iK) carries a phase of size ||K||, so an eps-relative error in that
    # phase is eps*||K|| absolute.  At scale 1e5, eigh measures 3.0e-11 and scipy
    # is the same order, so a flat 1e-12 would be asserting better-than-possible.
    # It keeps its teeth regardless: the ungated kernel reached 2.7e-07 here,
    # three orders above the loosest tolerance below.
    over = []
    for shape, mk in _GRID_SHAPES.items():
        for S in _GRID_SCALES:
            for f in (0.0, 1.49e-8, 1e-4, 1e-2):
                K = _grid_matrix(mk(S, f*S), np.random.default_rng(7))
                U = mg._expm_stack(-1j*K, expm_backend='numba')
                err = np.max(np.abs(U - _scipy_expm(K)))
                tol = max(1.0e-13, 2.0e-15*S)
                if err > tol:
                    over.append((shape, S, f, err, tol))
    assert not over, over


@requires_numba
def test_the_raw_kernel_flags_every_cell_where_it_is_inaccurate():
    r"""The severity has to be trustworthy, since the gate is built on it.

    For every grid cell, either the bare kernel is accurate or it reports a
    severity above :data:`magnus.expmkernels.SEV_TOL`.  A cell that is wrong
    *and* unflagged would be a silently wrong answer reaching any caller that
    used the documented contract correctly.
    """
    unflagged_and_wrong = []
    for shape, mk in _GRID_SHAPES.items():
        for S in _GRID_SCALES:
            for f in (0.0, 1.49e-8, 1e-4, 1e-2):
                K = _grid_matrix(mk(S, f*S), np.random.default_rng(7))
                U, _, sev = ek.expm_herm_stack(K[None])
                err = np.max(np.abs(U[0] - _scipy_expm(K)))
                if err > 1e-12 and sev <= ek.SEV_TOL:
                    unflagged_and_wrong.append((shape, S, f, err, sev))
    assert not unflagged_and_wrong, unflagged_and_wrong


@requires_numba
def test_the_gate_does_not_fire_on_physical_spectra():
    r"""If the gate fired on real work the speed-up would quietly evaporate.

    A Magnus slab has ``||Omega|| <~ pi`` by construction, and an ordinary 3nu
    constant-density call measures ``||K|| ~ 4``, so nothing physical should be
    declined.  Measured over a PREM chord scan and a solar slab chain, the
    declined fraction is 0.00%; this pins the principle on the cheap cases.
    """
    rng = np.random.default_rng(4)
    # hierarchical 3nu: (0, dm21^2, dm31^2)/2E, at accumulated phases that occur
    for phase in (1.0, 10.0, 100.0):
        spec = np.array([0.0, 7.5e-5/2.51e-3, 1.0])*phase
        K = _grid_matrix(spec, rng)
        _, _, sev = ek.expm_herm_stack(K[None])
        assert sev <= ek.SEV_TOL, (phase, sev)
    # and a Magnus slab, whose norm the convergence warning bounds near pi
    K = _herm((64, 3, 3), rng, scale=1.0)
    _, _, sev = ek.expm_herm_stack(np.ascontiguousarray(K))
    assert sev <= ek.SEV_TOL, sev


@requires_numba
def test_sev_tol_sits_inside_its_calibrated_window():
    r"""Pins the constant itself, since it was chosen by measurement.

    Cells at spectral scale 1e2 are safe and must stay on the kernel; cells at
    1e3 are not and must be declined.  ``m = tr(X^2)/6`` for those two scales
    brackets the constant, so a well-meaning edit that moves it an order of
    magnitude either way fails here rather than silently reopening the hole or
    silently costing the speed-up.
    """
    rng = np.random.default_rng(5)
    safe = _grid_matrix(_GRID_SHAPES['double-low'](1e2, 1e-8*1e2), rng)
    unsafe = _grid_matrix(_GRID_SHAPES['double-low'](1e3, 1e-8*1e3), rng)
    _, _, sev_safe = ek.expm_herm_stack(safe[None])
    _, _, sev_unsafe = ek.expm_herm_stack(unsafe[None])
    assert sev_safe < ek.SEV_TOL < sev_unsafe, (sev_safe, ek.SEV_TOL, sev_unsafe)


@requires_numba
@pytest.mark.parametrize('d', [1, 4, 5, 6])
def test_public_kernel_entry_point_refuses_unsupported_dimensions(d):
    r"""``expm_herm_stack`` must consult ``supports_dim``, not fall through to 3x3.

    It did not, and every dimension other than 2 went to the 3x3 kernel: at d=4
    that returned no exception, an error of 2.4 against ``scipy.linalg.expm``, a
    unitarity violation of 11.3, and uninitialised memory in the fourth
    eigenvalue; at d=1 it indexed ``K[i,2,1]`` with numba's bounds checking off
    and segfaulted the interpreter.  Only ``_expm_stack``'s own guard kept the
    library safe, and this function is public.
    """
    K = np.ascontiguousarray(_herm((2, d, d), np.random.default_rng(d)),
                             dtype=complex)
    with pytest.raises(ValueError, match='supports_dim'):
        ek.expm_herm_stack(K)
