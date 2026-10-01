# -*- coding: utf-8 -*-
"""Warnings that say "Shown once per session" are shown once per session (issue #205).

Python's default filter remembers a message per call site and forgets it whenever the filters
change, which every ``warnings.catch_warnings()`` block does -- including the package's own.
One run of the paper's cavity example printed the same ``MagnusConvergenceWarning`` four
times.  The first test runs in a fresh interpreter, since this suite has long since shown the
warning by the time it gets there.
"""

import subprocess
import sys
import textwrap
import warnings

import numpy as np
import pytest

import magnus.magnus as mm

_REPRO = textwrap.dedent('''
    import warnings, numpy as np
    import magnus.magnus as mm
    At = np.broadcast_to(-1j*np.array([[0, 50.0], [50.0, 0]]), (1, 2, 2, 2)).copy()
    def call():
        mm.evolution_operators_from_samples(At, np.array([1.0]), 4, 'gl')
    call(); call()
    with warnings.catch_warnings():       # any change to the filters, as the package makes
        warnings.simplefilter('ignore', DeprecationWarning)
    call()
''')

_AT = np.broadcast_to(-1j*np.array([[0, 50.0], [50.0, 0]]), (1, 2, 2, 2)).copy()


def _call():
    mm.evolution_operators_from_samples(_AT, np.array([1.0]), 4, 'gl')


def test_shown_once_across_a_change_of_filters():
    err = subprocess.run([sys.executable, '-c', _REPRO], capture_output=True, text=True,
                         check=True).stderr
    assert err.count('MagnusConvergenceWarning') == 1, err      # 2 before issue #205


def test_a_callers_filter_still_wins():
    """'always' still shows every one, 'error' still raises, 'ignore' still hides."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        _call()
        _call()
    assert sum(issubclass(w.category, mm.MagnusConvergenceWarning) for w in caught) == 2
    with warnings.catch_warnings():
        warnings.simplefilter('error', mm.MagnusConvergenceWarning)
        with pytest.raises(mm.MagnusConvergenceWarning):
            _call()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('ignore')
        _call()
    assert not caught
