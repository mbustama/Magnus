# -*- coding: utf-8 -*-
"""Concurrent calls from several threads of one process (issue #153).

Four pieces of per-call state lived in module globals that a context manager swapped in and
out: the two slab-norm sinks of ``magnus.magnus`` and the engine trace and disabled set of
``magnus.oscprob``.  Two threads in a batched Earth energy scan overwrote each other's sink --
``ValueError: need at least one array to concatenate``, or a different refinement path and a
different answer -- and a thread asking for ``strategy='magnus'`` disabled the other engines for
every thread while it ran.  They are context variables now, so each thread sees its own.

``n_jobs`` uses processes, not threads, and was never affected.
"""

import threading
import warnings

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.oscprob as op

K = gd.UNIT_KM
E3 = np.array([1.0, 2.0, 3.0])*1e9


def _run_concurrently(fn, n_threads, n_rounds):
    """Every result of n_threads x n_rounds concurrent calls: an array, or the exception."""
    results = []
    for _ in range(n_rounds):
        out = [None]*n_threads

        def run(i):
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    out[i] = np.asarray(fn())
            except Exception as err:          # noqa: BLE001 -- collected and asserted on
                out[i] = err

        threads = [threading.Thread(target=run, args=(i,)) for i in range(n_threads)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        results.extend(out)
    return results


@pytest.mark.parametrize('fn', [
    lambda: op.osc_prob_3nu_earth(E3, costhz=-0.7, L=8000.0*K),
    lambda: op.osc_prob_5nu_earth(E3, costhz=-0.7, L=8000.0*K),
    lambda: op.osc_prob_3nu_earth_nsi(E3, costhz=-0.7, L=8000.0*K, eps_ee=0.1),
], ids=['3nu_earth', '5nu_earth', '3nu_earth_nsi'])
def test_concurrent_earth_scans_match_the_serial_answer_bit_for_bit(fn):
    """On main, 14 to 20 of 100 such calls raised and 7 more returned different numbers."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ref = np.asarray(fn())
    results = _run_concurrently(fn, n_threads=8, n_rounds=10)
    errors = [r for r in results if isinstance(r, Exception)]
    assert not errors, repr(errors[0])
    assert all(np.array_equal(r, ref) for r in results)


def test_a_strategy_in_one_thread_does_not_reach_another():
    """strategy='magnus' disables the other engines for its own call only, and strategy_info
    reports the trace of its own call.  On main, 30 of 40 such calls differed from serial."""
    R = gd.SUN_RADIUS*K

    def sun(strategy):
        info = {}
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            P = np.asarray(op.osc_prob_3nu_sun(1.0e7, R, 0.0, strategy=strategy,
                                               strategy_info=info))
        return P, info['engine'], len(info['trace'])

    ref = {s: sun(s) for s in ('auto', 'magnus')}
    assert ref['auto'][1] != ref['magnus'][1]      # otherwise the test proves nothing
    for _ in range(5):
        out = {}

        def run(key, strategy):
            out[key] = sun(strategy)

        threads = [threading.Thread(target=run, args=((i, s), s))
                   for i in range(2) for s in ('auto', 'magnus')]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
        for (_, s), (P, engine, n_trace) in out.items():
            assert engine == ref[s][1] and n_trace == ref[s][2]
            assert np.array_equal(P, ref[s][0])
