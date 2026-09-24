# -*- coding: utf-8 -*-
r"""The paper figure cache's key survives a change of machine.

`notebooks/paper_figure_cache.json` stores what the paper's figures are drawn from, keyed
on a hash of the configuration.  That hash used to include raw float64 bytes of computed
profiles -- the output of `np.exp`, and of the power behind `np.logspace`.  NumPy dispatches
those through SIMD paths chosen by CPU feature, so the last unit in the last place is not
reproducible between machines, and neither was the key.

It bit exactly as you would expect.  The same commit, byte for byte, hit the cache on one
GitHub runner and missed on the next; with MAGNUS_PAPER_CACHE_ONLY set that is a failed
build rather than a slow one, so `main` went red on a merge that changed nothing.

Floats now enter the hash at FINGERPRINT_DIGITS significant figures.  A float64 carries
fifteen to seventeen, so twelve is clear of last-place noise and far below anything a
configuration change could hide in.  These tests hold both ends of that: a one-ULP
difference must NOT move the key, and a change a thousand times larger must.
"""

import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT/'notebooks'))


@pytest.fixture(scope='module')
def hashes():
    """The two hash functions, lifted from the generator rather than reimplemented.

    They are defined inside a notebook cell -- a string literal in make_notebooks.py --
    so they cannot be imported. Executing the definitions out of that source is what keeps
    this test measuring the real hash: a copy here would pass while the notebook's own
    version stayed broken, which is the whole failure mode being guarded.
    """
    import hashlib

    source = (ROOT/'notebooks'/'make_notebooks.py').read_text()
    namespace = {'np': np, 'hashlib': hashlib}
    for anchor in ('FINGERPRINT_DIGITS = 12', 'def _hashable', 'def _python_scalars',
                   'class _Rounded', 'def _rounded', 'def fingerprint', 'def _shallow_fingerprint',
                   'def legacy_fingerprint', 'def earlier_fingerprints'):
        assert anchor in source, (
            'make_notebooks.py no longer contains %r; the paper cache hash has been '
            'restructured and this test needs updating' % anchor)
        start = source.index(anchor)
        exec(compile(source[start:source.index('\n\n\n', start)], '<fingerprint>', 'exec'),
             namespace)
    return namespace['fingerprint'], namespace['legacy_fingerprint'], namespace['_python_scalars']


# Values at the magnitude the real keys carry: a matter potential is ~1e-13 eV, and a
# hash that is only tested on numbers near one would not exercise the scaling at all.
#
# No exact zero.  `np.nextafter(0.0, inf)` is 5e-324, which is not a last-place difference
# but an infinite relative one, and no significant-digit rounding absorbs it -- nor should
# it.  The quantity this guards is a computed profile, whose samples are nonzero.
SAMPLES = np.array([1.23456789012345e-13, 4.2e-14, 9.87e-13, -5.5e-13])


def test_one_ulp_does_not_move_the_key(hashes):
    """The failure that turned a latent problem into a red main."""
    fingerprint, _, _ = hashes
    nudged = np.nextafter(SAMPLES, np.inf)
    assert fingerprint(SAMPLES) == fingerprint(nudged), (
        'a one-ULP difference changes the cache key, so the key is machine-dependent '
        'again and the same commit can hit on one runner and miss on another')


def test_the_old_hash_really_did_move(hashes):
    """So the test above is measuring the fix and not a property floats had anyway."""
    _, legacy, _ = hashes
    assert legacy(SAMPLES) != legacy(np.nextafter(SAMPLES, np.inf))


def test_a_real_change_still_moves_the_key(hashes):
    """Tolerance bought at the twelfth digit, and no further."""
    fingerprint, _, _ = hashes
    for rel in (1.0e-9, 1.0e-6, 1.0e-3):
        assert fingerprint(SAMPLES) != fingerprint(SAMPLES*(1.0 + rel)), (
            'a relative change of %g leaves the cache key alone, so a moved '
            'configuration would be served stale values' % rel)


def test_plain_floats_are_quantized_too(hashes):
    """Tolerances and baselines reach the key as Python floats, not only as arrays."""
    fingerprint, _, _ = hashes
    x = 3000.0*1.0e9
    assert fingerprint(x) == fingerprint(np.nextafter(x, np.inf))
    assert fingerprint(x) != fingerprint(x*(1.0 + 1.0e-9))


def test_every_stored_fingerprint_is_a_full_hash():
    """A truncated label would compare unequal forever and recompute on every build.

    Written after an aborted migration that scraped keys from an error message, which
    reports only the first twelve characters.
    """
    import json
    blob = json.loads((ROOT/'notebooks'/'paper_figure_cache.json').read_text())
    short = {name: entry['fingerprint'] for name, entry in blob.items()
             if isinstance(entry, dict) and 'fingerprint' in entry
             and len(entry['fingerprint']) != 64}
    assert not short, 'cache entries carry truncated fingerprints: %s' % short


def test_a_numpy_scalar_hashes_as_the_python_scalar_it_holds(hashes):
    """NumPy 2 writes repr(np.float64(0.3)) as 'np.float64(0.3)', where NumPy 1 wrote '0.3'.

    A key holding one, such as a sorted dict of mixing parameters, therefore hit the cache under
    NumPy 1 and missed on CI under NumPy 2.  Every NumPy scalar in a key enters as the Python
    scalar it holds, whose repr is what NumPy 1 printed, so no stored key moves.
    """
    fingerprint, _, python_scalars = hashes
    plain = [('D41', 1.0), ('n', 3), ('s14', 0.31622776601683794), ('sterile', True)]
    held = [('D41', np.float64(1.0)), ('n', np.int64(3)),
            ('s14', np.sqrt(np.float64(0.10))), ('sterile', np.bool_(True))]
    assert fingerprint(held) == fingerprint(plain)
    assert fingerprint({'a': np.float64(0.5)}) == fingerprint({'a': 0.5})
    converted = python_scalars(held)
    assert all(type(v) in (float, int, bool) for _, v in converted)
    assert repr(converted) == repr(plain)


def test_a_nested_float_is_quantized_too(hashes):
    """A float inside a dict or a list of pairs reached the hash through repr, digit for digit.

    The long-range couplings of notebook 28's solar_long_range key are computed, and their
    last place differed between two GitHub runners: the section passed on the branch and
    missed on main.  One ulp inside a dict must not move the key; a change a thousand times
    larger still must.
    """
    fingerprint, _, _ = hashes
    g = 8.641658082371646e-54
    base = fingerprint('lri', {'1': g, '0.1': 3.193173940586491e-53}, [('s12', 0.55)])
    ulp = fingerprint('lri', {'1': float(np.nextafter(g, np.inf)), '0.1': 3.193173940586491e-53},
                      [('s12', 0.55)])
    assert ulp == base
    assert fingerprint('lri', {'1': g*(1 + 1e-9), '0.1': 3.193173940586491e-53},
                       [('s12', 0.55)]) != base


def test_an_entry_stored_under_the_shallow_hash_is_recognised():
    """Quantizing nested floats changed the key of every section that nests one.  Those entries
    are relabelled, not recomputed: the key the same configuration had before is among the
    earlier ones."""
    import hashlib

    source = (ROOT/'notebooks'/'make_notebooks.py').read_text()
    ns = {'np': np, 'hashlib': hashlib}
    for anchor in ('FINGERPRINT_DIGITS = 12', 'def _hashable', 'def _python_scalars',
                   'class _Rounded', 'def _rounded', 'def fingerprint', 'def _shallow_fingerprint',
                   'def legacy_fingerprint', 'def earlier_fingerprints'):
        start = source.index(anchor)
        exec(compile(source[start:source.index('\n\n\n', start)], '<fingerprint>', 'exec'), ns)
    parts = ('lri', [0.1, 1.0], {'1': 8.641658082371646e-54}, sorted({'s12': 0.55}.items()))
    assert ns['fingerprint'](*parts) != ns['_shallow_fingerprint'](*parts)
    assert ns['_shallow_fingerprint'](*parts) in ns['earlier_fingerprints'](*parts)
    assert ns['fingerprint'](*parts) not in ns['earlier_fingerprints'](*parts)
