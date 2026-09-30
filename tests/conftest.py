# -*- coding: utf-8 -*-
"""Pytest configuration: make the magnus package importable without installation, and skip
the tests that need a source checkout when there is none.

Some tests check the repository rather than the library -- that the documentation, the
notebooks, the paper's assets and the CI workflows agree with the code.  The sdist ships
``src/`` and ``tests/`` but not ``docs/``, ``notebooks/``, ``resources/`` or ``.github/``, so
from an unpacked sdist those tests could only fail (issue #164 §1).  They are marked
``checkout_only`` and skipped, with the reason, wherever ``docs/`` and ``notebooks/`` are
missing.  In a checkout nothing is skipped.
"""

import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# From an unpacked sdist with src/ removed, the tests run against the installed package.
if os.path.isdir(os.path.join(ROOT, 'src')):
    sys.path.insert(0, os.path.join(ROOT, 'src'))

CHECKOUT = all(os.path.isdir(os.path.join(ROOT, d)) for d in ('docs', 'notebooks'))
"""bool: Whether the tests run from a source checkout, which has docs/ and notebooks/."""

CHECKOUT_REASON = ("needs a source checkout: reads files the sdist does not ship "
                   "(docs/, notebooks/, resources/, .github/)")


def skip_module_without_checkout():
    """For a test module that reads checkout-only files when it is imported."""
    if not CHECKOUT:
        pytest.skip(CHECKOUT_REASON, allow_module_level=True)


def pytest_configure(config):
    config.addinivalue_line(
        'markers', 'checkout_only: needs a source checkout; skipped from an unpacked sdist')


def pytest_collection_modifyitems(config, items):
    if CHECKOUT:
        return
    skip = pytest.mark.skip(reason=CHECKOUT_REASON)
    for item in items:
        if 'checkout_only' in item.keywords:
            item.add_marker(skip)
