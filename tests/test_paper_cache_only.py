# -*- coding: utf-8 -*-
r"""With MAGNUS_PAPER_CACHE_ONLY set, notebook 28 reads its cache and never recomputes.

Continuous integration rebuilds every notebook on each push, and the paper's notebook stores
what its figures are drawn from in `notebooks/paper_figure_cache.json`.  Both READMEs promise
that the variable forbids recomputing anything: a section whose configuration moved stops the
build and names itself, so that it is re-measured on a quiet machine and the cache committed.

Only the scan and timing sections kept that promise.  `cached`, the helper behind about thirty
others, printed "configuration moved, recomputing", computed on the runner, and wrote a cache
nobody would commit (issue #63).  These tests hold the contract on the helper itself, lifted
out of the generator the way `test_paper_cache_key_is_portable.py` lifts the hash: a copy here
would pass while the notebook's own version stayed broken.
"""

import hashlib
import json
import os
import pathlib
import platform
import time

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
ANCHORS = ('FINGERPRINT_DIGITS = 12', 'def _hashable', 'def _python_scalars', 'def fingerprint',
           'def legacy_fingerprint', 'def cache_miss', 'def write_cache', 'def cached(section')


@pytest.fixture
def notebook(tmp_path, monkeypatch):
    """Notebook 28's cache helpers, writing to a scratch cache file."""
    monkeypatch.delenv('MAGNUS_PAPER_REDO', raising=False)
    source = (ROOT/'notebooks'/'make_notebooks.py').read_text()
    ns = dict(np=np, hashlib=hashlib, json=json, os=os, time=time, platform=platform,
              pathlib=pathlib)
    for anchor in ANCHORS:
        assert anchor in source, (
            'make_notebooks.py no longer contains %r; the paper cache has been restructured '
            'and this test needs updating' % anchor)
        start = source.index(anchor)
        exec(compile(source[start:source.index('\n\n\n', start)], '<notebook 28>', 'exec'), ns)
    ns['MP_CACHE'] = tmp_path/'paper_figure_cache.json'
    ns['CACHE_ONLY'] = False
    return ns


@pytest.fixture
def compute():
    calls = []

    def run():
        calls.append(1)
        return {'P': [0.5, 0.25]}

    run.calls = calls
    return run


def test_without_the_variable_a_miss_is_recomputed_and_stored(notebook, compute):
    assert notebook['cached']('solar_3_1', ('config', 1.0), compute) == {'P': [0.5, 0.25]}
    assert len(compute.calls) == 1
    assert 'solar_3_1' in json.loads(notebook['MP_CACHE'].read_text())


def test_a_hit_is_read_back_under_the_variable(notebook, compute):
    notebook['cached']('solar_3_1', ('config', 1.0), compute)
    notebook['CACHE_ONLY'] = True
    assert notebook['cached']('solar_3_1', ('config', 1.0), compute) == {'P': [0.5, 0.25]}
    assert len(compute.calls) == 1


def test_a_missing_section_stops_the_build(notebook, compute):
    notebook['CACHE_ONLY'] = True
    with pytest.raises(RuntimeError, match="section 'solar_3_1'.*MAGNUS_PAPER_CACHE_ONLY"):
        notebook['cached']('solar_3_1', ('config', 1.0), compute)
    assert not compute.calls
    assert not notebook['MP_CACHE'].exists()


def test_a_moved_configuration_stops_the_build_and_keeps_the_stored_entry(notebook, compute):
    notebook['cached']('solar_3_1', ('config', 1.0), compute)
    before = notebook['MP_CACHE'].read_text()
    notebook['CACHE_ONLY'] = True
    with pytest.raises(RuntimeError, match="section 'solar_3_1'"):
        notebook['cached']('solar_3_1', ('config', 2.0), compute)
    assert len(compute.calls) == 1
    assert notebook['MP_CACHE'].read_text() == before


def test_the_variable_wins_over_a_request_to_redo(notebook, compute, monkeypatch):
    """Asking to recompute and forbidding it at once stops, as the timing section does."""
    notebook['cached']('solar_3_1', ('config', 1.0), compute)
    monkeypatch.setenv('MAGNUS_PAPER_REDO', '1')
    notebook['CACHE_ONLY'] = True
    with pytest.raises(RuntimeError, match="section 'solar_3_1'"):
        notebook['cached']('solar_3_1', ('config', 1.0), compute)
    assert len(compute.calls) == 1
