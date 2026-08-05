# -*- coding: utf-8 -*-
r"""Tests for the palindromic-profile optimisation.

The mechanism evaluates :math:`A` on the first half of a symmetric slab chain and derives the
rest by reversal.  Two things are worth stating about what is tested here, because both were
mistakes made while building it:

* **Assert that the gate fired.**  A test that compares the mirrored route against the plain one
  proves nothing if the gate quietly declined -- both sides then run the same code and agree to
  0.0.  That happened: an early version of these checks reported exact agreement at odd slab
  counts, which looked like success and was the gate refusing.  Every accuracy test below asserts
  :func:`magnus.magnus._mirror_applies` first.

* **Exercise odd slab counts.**  The middle slab of an odd chain has no mirror partner.  The
  prototype this was built from allocated the output with ``np.empty`` and never wrote that slab,
  returning uninitialised memory -- wrong by 7.1e-01 at 31 slabs, while every even count was
  exact to 1e-15.
"""

import numpy as np
import pytest

from magnus import magnus as mg


L_CHORD = 5.8115714094e13          # a PREM chord's length, in eV^-1
V0 = 4.0e-13


def _edges(n, length=L_CHORD):
    r"""A uniform chain of ``n`` slabs, built the way a caller would."""
    a = (length/n)*np.arange(n + 1, dtype=float)
    return np.stack([a[:-1], a[1:]], axis=1)


def _h_vac():
    return np.array([[0.0, 1.0e-13, 2.0e-14],
                     [1.0e-13, 3.0e-13, 5.0e-14],
                     [2.0e-14, 5.0e-14, 7.0e-13]], dtype=complex)


def _A_symmetric(l):
    r""":math:`A(l) = A(L - l)` to machine precision."""
    l = np.asarray(l, dtype=float)
    v = V0*(1.0 + 0.5*np.cos(2.0*np.pi*l/L_CHORD))
    e00 = np.zeros((3, 3))
    e00[0][0] = 1.0
    return -1j*(_h_vac() + v[..., None, None]*e00)


def _A_monotonic(l):
    r"""Solar-like: strictly decreasing, so emphatically not a palindrome."""
    l = np.asarray(l, dtype=float)
    v = V0*np.exp(-3.0*l/L_CHORD)
    e00 = np.zeros((3, 3))
    e00[0][0] = 1.0
    return -1j*(_h_vac() + v[..., None, None]*e00)


def _span(edges):
    return (edges[0, 0], edges[-1, 1])


# --------------------------------------------------------------- the predicate

@pytest.mark.parametrize('values, expected', [
    ([1.0, 2.0, 1.0], True),
    ([1.0, 2.0, 3.0], False),
    ([1.0, 2.0, 2.0, 1.0], True),
    ([1.0, 2.0, 3.0, 1.0], False),
    ([5.0], True),                       # too short to be anything else
    ([], True),
])
def test_palindromic_decides_on_exact_equality(values, expected):
    assert mg.palindromic(np.array(values, dtype=float)) is expected


def test_palindromic_of_nothing_is_true():
    assert mg.palindromic() is True


def test_palindromic_requires_every_array_to_agree():
    r"""Symmetric widths across an asymmetric profile is not a palindrome."""
    assert mg.palindromic(np.array([1.0, 2.0, 1.0]),
                          np.array([3.0, 4.0, 5.0])) is False


def test_palindromic_is_exact_not_approximate():
    r"""One ulp of asymmetry is asymmetry.

    A tolerance here would return a different answer for a nearly-symmetric profile depending on
    how nearly symmetric it is, which is a silent accuracy change keyed on an invisible property.
    """
    w = np.array([1.0, 2.0, 1.0])
    w[2] = np.nextafter(w[2], 2.0)
    assert mg.palindromic(w) is False


# -------------------------------------------------------------------- the gate

def test_gate_needs_a_declaration():
    r"""Palindromic widths alone must NOT arm the mirror.

    This is the trap the gate exists for: a monotonic profile on a uniform grid has palindromic
    widths, and mirroring it is wrong by 3.3e-01.  Symmetry of the grid says nothing about
    symmetry of the profile.
    """
    e = _edges(64)
    assert mg._mirror_applies(e, e[:, 1] - e[:, 0], None) is False


def test_gate_declines_a_sub_range_of_the_declared_interval():
    r"""Engines call this layer on blocks of a chain; a sub-range of a symmetric profile is not
    itself symmetric, so the declared span must be matched, not merely overlapped."""
    e = _edges(64)
    w = e[:, 1] - e[:, 0]
    full = _span(e)
    assert mg._mirror_applies(e, w, full) is True
    assert mg._mirror_applies(e[:32], w[:32], full) is False
    assert mg._mirror_applies(e[16:48], w[16:48], full) is False


def test_gate_declines_an_uneven_grid():
    r"""The width bound admits floating-point noise, not an actually uneven grid."""
    e = _edges(64)
    w = e[:, 1] - e[:, 0]
    e2 = e.copy()
    e2[0, 1] = e2[0, 0] + 1.5*w[0]        # one slab half again as wide
    e2[1, 0] = e2[0, 1]
    assert mg._mirror_applies(e2, e2[:, 1] - e2[:, 0], _span(e2)) is False


def test_gate_fires_across_the_practical_range_of_slab_counts():
    r"""Regression on the bound's scaling.

    ``edges[:, 1] - edges[:, 0]`` returns a relative asymmetry of about ``0.30*n*eps``, so a
    fixed ulp bound passes at small counts and fails at large ones -- an earlier version admitted
    7 slabs and refused 31, which would have shipped the optimisation as a near-total no-op.
    """
    for n in (2, 3, 7, 31, 32, 63, 64, 127, 128, 129, 256, 511, 512, 1024):
        e = _edges(n)
        assert mg._mirror_applies(e, e[:, 1] - e[:, 0], _span(e)) is True, n


def test_use_palindrome_switch_disarms_the_gate(monkeypatch):
    e = _edges(64)
    w = e[:, 1] - e[:, 0]
    monkeypatch.setattr(mg, 'USE_PALINDROME', False)
    assert mg._mirror_applies(e, w, _span(e)) is False


# ------------------------------------------------------------------- the route

@pytest.mark.parametrize('n_slabs', [2, 3, 7, 31, 32, 63, 64, 127, 128, 129])
@pytest.mark.parametrize('order, method', [
    (1, 'gl'), (2, 'gl'), (4, 'gl'), (6, 'gl'),
    (2, 'trapezoid'), (4, 'trapezoid'), (6, 'trapezoid'),
    (4, 'simpson'), (6, 'simpson'),
])
def test_mirrored_route_agrees_with_the_plain_one(n_slabs, order, method):
    r"""Agreement is to rounding, not bitwise, and that is inherent.

    The mirrored slab's nodes are reached as ``(L - b) + h*s`` on one route and ``a + h*s`` on the
    other: different floating-point expressions for the same real number.
    """
    e = _edges(n_slabs)
    assert mg._mirror_applies(e, e[:, 1] - e[:, 0], _span(e)) is True, "gate declined"
    kw = dict(order=order, integration_method=method, n_tpts_per_slab=41,
              validate_input=False)
    mirrored = mg.magnus_expansion_multislab(_A_symmetric, e, symmetric_over=_span(e), **kw)
    plain = mg.magnus_expansion_multislab(_A_symmetric, e, **kw)
    assert np.max(np.abs(mirrored - plain)) < 1e-13


@pytest.mark.parametrize('n_slabs', [64, 65])
def test_the_mirror_halves_the_hamiltonian_evaluations(n_slabs):
    r"""The saving is halved calls to the caller's Hamiltonian, and nothing else."""
    seen = []

    def counted(l):
        seen.append(np.size(l))
        return _A_symmetric(l)

    e = _edges(n_slabs)
    kw = dict(order=4, integration_method='gl', validate_input=False)

    seen.clear()
    mg.magnus_expansion_multislab(counted, e, symmetric_over=_span(e), **kw)
    mirrored = max(seen)

    seen.clear()
    mg.magnus_expansion_multislab(counted, e, **kw)
    plain = max(seen)

    # Odd counts evaluate the unpaired middle slab too, so the ratio is just under two.
    assert plain/mirrored >= 1.95


def test_the_middle_slab_of_an_odd_chain_is_computed_not_left_uninitialised():
    r"""Regression for the prototype's odd-count bug.

    ``Om[:m]`` and ``Om[n-m:]`` with ``m = n//2`` skip index ``m`` when ``n`` is odd, so the
    middle slab came back as whatever was in memory -- 7.1e-01 wrong at 31 slabs, while every
    even count was exact.  Checked per slab, so the middle one cannot hide behind the others.
    """
    for n in (3, 7, 31, 63, 129):
        e = _edges(n)
        assert mg._mirror_applies(e, e[:, 1] - e[:, 0], _span(e)) is True, n
        kw = dict(order=4, integration_method='gl', validate_input=False)
        mirrored = mg.magnus_expansion_multislab(_A_symmetric, e,
                                                 symmetric_over=_span(e), **kw)
        plain = mg.magnus_expansion_multislab(_A_symmetric, e, **kw)
        per_slab = np.max(np.abs(mirrored - plain), axis=(1, 2))
        assert per_slab[n//2] < 1e-13, (n, per_slab[n//2])
        assert np.max(per_slab) < 1e-13


def test_a_false_declaration_is_wrong_and_that_is_why_it_is_not_user_facing():
    r"""Pins the hazard that keeps ``symmetric_over`` out of the public API.

    The declaration is not checked and cannot cheaply be: testing it needs the evaluations it
    exists to avoid.  Declared of a monotonic profile, the mirror is wrong by ~1e-01 -- so this
    parameter is set by entry points that know the geometry, never by a caller's assertion.
    """
    e = _edges(128)
    kw = dict(order=4, integration_method='gl', validate_input=False)
    mirrored = mg.magnus_expansion_multislab(_A_monotonic, e, symmetric_over=_span(e), **kw)
    plain = mg.magnus_expansion_multislab(_A_monotonic, e, **kw)
    assert np.max(np.abs(mirrored - plain)) > 1e-3


def test_no_declaration_leaves_the_shipped_path_bit_identical():
    r"""The feature must be inert unless asked for: same call, same bits."""
    e = _edges(128)
    for order, method in ((2, 'gl'), (4, 'gl'), (4, 'trapezoid')):
        kw = dict(order=order, integration_method=method, n_tpts_per_slab=41,
                  validate_input=False)
        a = mg.magnus_expansion_multislab(_A_monotonic, e, **kw)
        b = mg.magnus_expansion_multislab(_A_monotonic, e, symmetric_over=None, **kw)
        assert np.array_equal(a, b)


# ------------------------------------------------------- the Earth entry points

def _chord(costhz):
    from magnus import earth, globaldefs as gd
    return earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM


def test_earth_chord_symmetry_declares_only_the_whole_chord():
    r"""A chord is symmetric over its full length and over no shorter prefix."""
    from magnus import oscprob
    c = -0.9
    L = _chord(c)
    assert oscprob._earth_chord_symmetry(c, L) == (0.0, L)
    assert oscprob._earth_chord_symmetry(c, [L, L]) == (0.0, L)
    assert oscprob._earth_chord_symmetry(c, 0.5*L) is None
    assert oscprob._earth_chord_symmetry(c, [L, 0.9*L]) is None
    assert oscprob._earth_chord_symmetry(0.5, L) is None      # up-going: no chord
    assert oscprob._earth_chord_symmetry(None, L) is None


def test_prem_crossings_are_exactly_symmetric_about_the_chord_midpoint():
    r"""What licenses the declaration, checked rather than assumed.

    The two roots of the crossing quadratic come out as ``d/2 +/- s``, so the crossings mirror
    exactly -- not nearly.  If this ever stops holding, the Earth slab grid stops being
    palindromic and the mirror silently stops firing, which is a performance regression that no
    accuracy test would catch.
    """
    from magnus import earth
    for costhz in (-0.05, -0.1, -0.3, -0.5, -0.7, -0.9, -0.99, -1.0):
        d = earth.distance_traveled_inside_earth(costhz)
        tb = earth.prem_layer_edges_along_chord(costhz)
        if len(tb) == 0:
            continue
        assert np.max(np.abs(tb + tb[::-1] - d)) == 0.0, costhz


@pytest.mark.parametrize('costhz', [-0.2, -0.5, -0.9])
def test_an_earth_single_point_takes_the_mirrored_path(costhz, monkeypatch):
    r"""The wiring works end to end, asserted by watching the gate rather than the answer."""
    from magnus import magnus as mgcore, oscprob, globaldefs as gd
    seen = []
    original = mgcore._mirror_applies
    monkeypatch.setattr(mgcore, '_mirror_applies',
                        lambda e, w, s: seen.append(bool(original(e, w, s))) or seen[-1])
    oscprob.osc_prob_3nu_earth(2.0*gd.UNIT_GEV, costhz=costhz, L=_chord(costhz),
                               nu_i=1, nu_f=0, validate_input=False)
    assert seen and all(seen), seen


def test_a_shorter_baseline_does_not_take_the_mirrored_path(monkeypatch):
    r"""The chord is not symmetric over half of itself, and the span check is what notices."""
    from magnus import magnus as mgcore, oscprob, globaldefs as gd
    seen = []
    original = mgcore._mirror_applies
    monkeypatch.setattr(mgcore, '_mirror_applies',
                        lambda e, w, s: seen.append(bool(original(e, w, s))) or seen[-1])
    oscprob.osc_prob_3nu_earth(2.0*gd.UNIT_GEV, costhz=-0.9, L=0.5*_chord(-0.9),
                               nu_i=1, nu_f=0, validate_input=False)
    assert seen and not any(seen), seen


def test_a_solar_profile_never_arms_the_mirror(monkeypatch):
    r"""Monotonic profiles must not reach this path; a false positive here is a correctness bug,
    not a missed optimisation."""
    from magnus import magnus as mgcore, oscprob, matter, globaldefs as gd
    seen = []
    original = mgcore._mirror_applies
    monkeypatch.setattr(mgcore, '_mirror_applies',
                        lambda e, w, s: seen.append(bool(original(e, w, s))) or seen[-1])
    profile = matter.exp_density_profile(gd.NUM_DENSITY_E_SUN_CENTRAL, gd.L_SCALE_SUN)
    params = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
    oscprob.osc_prob_matter_std_potential(
        2, profile, 50.0e6, 0.5*gd.SUN_RADIUS*gd.UNIT_KM,
        {'sth': params['s12'], 'Dm2': params['D21']},
        L0=0.0, density_is_of_number_of_electrons=True)
    assert not any(seen), seen


@pytest.mark.parametrize('costhz', [-0.2, -0.5, -0.8, -0.9, -1.0])
def test_the_earth_answer_moves_only_at_rounding(costhz, monkeypatch):
    r"""The mirror changes the Earth path's bits, and this pins by how much.

    Worst measured across 15 (costhz, energy) configurations is 8.6e-15 relative; the bound here
    is 1e-12, loose enough not to be flaky and tight enough that a real regression trips it.
    """
    from magnus import magnus as mgcore, oscprob, globaldefs as gd
    L = _chord(costhz)
    kw = dict(costhz=costhz, L=L, nu_i=1, nu_f=0, validate_input=False)
    with_mirror = float(oscprob.osc_prob_3nu_earth(2.0*gd.UNIT_GEV, **kw))
    monkeypatch.setattr(mgcore, 'USE_PALINDROME', False)
    without = float(oscprob.osc_prob_3nu_earth(2.0*gd.UNIT_GEV, **kw))
    assert abs(with_mirror - without) <= 1e-12*max(abs(without), 1e-3)
