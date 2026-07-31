# -*- coding: utf-8 -*-
"""Tests of the version resolver (magnus.version).

The version is written down in exactly one place, the ``version`` field of
``pyproject.toml``, and reaches users through two quite different routes:
``importlib.metadata`` for an installed distribution, and a direct read of
``pyproject.toml`` for a source checkout on ``sys.path``.  Only the first
of those runs in a normal test session, so a coverage run found the entire
fallback unexecuted -- including the branch that decides whether an
installed user sees the real number or ``0.0.0+unknown``.

That failure mode has already happened once here: the lookup used the
*import* package name ``magnus`` instead of the distribution name
``magnuspy``, which raises PackageNotFoundError and falls through to a
``pyproject.toml`` that is absent from an installed wheel.  A source
checkout looked perfectly fine while every installed user would have
reported ``0.0.0+unknown``.  These tests exercise both routes so the two
cannot silently disagree again.
"""

import re
from pathlib import Path

import magnus.version as ver

PYPROJECT = Path(__file__).resolve().parents[1] / 'pyproject.toml'


def declared_version():
    """The version as literally written in pyproject.toml."""
    text = PYPROJECT.read_text(encoding='utf-8')
    match = re.search(r'^\s*version\s*=\s*["\']([^"\']+)["\']', text, re.MULTILINE)
    assert match is not None, "pyproject.toml has no version field"
    return match.group(1)


def test_version_from_pyproject_matches_the_declared_version():
    """The source-checkout route. This is the branch that never runs when
    the package is installed, which is exactly why it needs a test."""
    assert ver._version_from_pyproject() == declared_version()


def test_resolved_version_matches_the_declared_version():
    """Whichever route was taken in this session, the answer must be the
    number in pyproject.toml -- never the placeholder."""
    assert ver.__version__ == declared_version()
    assert ver.__version__ != '0.0.0+unknown'


def test_metadata_lookup_uses_the_distribution_name_not_the_import_name():
    """Regression test for the naming trap. `magnuspy` is the distribution,
    `magnus` the import package; querying the latter raises
    PackageNotFoundError. The test is meaningful only when the package is
    actually installed, so it skips rather than lies when it is not."""
    import importlib.metadata as md
    import pytest

    try:
        md.version('magnuspy')
    except md.PackageNotFoundError:
        pytest.skip("magnuspy is not installed in this environment")

    assert ver._version_from_metadata() == declared_version()


def test_metadata_lookup_returns_empty_when_the_distribution_is_missing(monkeypatch):
    """The not-installed case: PackageNotFoundError must become an empty
    string so the caller falls through to pyproject.toml, rather than
    propagating out of an import."""
    import importlib.metadata as md

    def raise_not_found(_name):
        raise md.PackageNotFoundError(_name)

    monkeypatch.setattr(md, 'version', raise_not_found)

    assert ver._version_from_metadata() == ''


def test_pyproject_read_failure_returns_empty(monkeypatch):
    """A checkout without pyproject.toml (a vendored copy, say) must yield
    '' so the caller reaches the '0.0.0+unknown' placeholder, instead of
    raising OSError at import time and taking the whole package down."""
    def raise_oserror(*args, **kwargs):
        raise OSError("no such file")

    monkeypatch.setattr(Path, 'read_text', raise_oserror)

    assert ver._version_from_pyproject() == ''


def test_pyproject_without_a_version_field_returns_empty(monkeypatch):
    """The file exists but carries no version line: still '', not a crash
    on a None match."""
    monkeypatch.setattr(Path, 'read_text', lambda *a, **k: "[project]\nname = 'magnuspy'\n")

    assert ver._version_from_pyproject() == ''
