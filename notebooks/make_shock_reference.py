# -*- coding: utf-8 -*-
r"""Freezes the ``solve_ivp`` ground truth used by ``14_magnus_supernova_shock.ipynb``.

The notebook compares Mag(nu)s against a tight-tolerance ODE solution of the same
Schroedinger equation.  That solution is a **constant of the physics configuration** -- it
does not depend on Mag(nu)s at all -- so recomputing it on every run spends a quarter of an
hour re-deriving a number that cannot have changed.  Measured, notebook 14 took 933 s, a
third of the whole suite and its largest single cost, and the oracle was nearly all of it.

So it is computed once, here, and stored in ``shock_reference.json`` as hexadecimal floats,
which round-trip exactly: a reader gets the bits it was computed from rather than a decimal
rendering of them.

**The definitions are not copied.**  An earlier version of this script transcribed the
shock profile by hand and got it wrong -- inventing a power-law rarefaction in place of the
Fogli et al. form, and the wrong density normalisation -- which `solve_ivp` reported as an
overflow rather than as a wrong answer, and which would otherwise have frozen a reference
for a profile the notebook does not use.  So the notebook's own cells are executed here,
straight out of ``make_notebooks.py``, up to the point where the reference is needed.  There
is one definition of the physics and this reads it.

**What is frozen is the oracle, never the thing under test.**  Every Mag(nu)s number in
that notebook is still computed live; only its reference is stored.  The one risk that
introduces is a stale oracle outliving a change to the profile, so the file also carries a
fingerprint -- the electron density sampled along the ray -- and the notebook refuses a
reference whose fingerprint does not match what it just built.

Run ``python notebooks/make_shock_reference.py`` after any change to the shock profile, the
energy, the ray, or the sampled baselines.
"""

import json
import pathlib
import sys

import numpy as np
from scipy.integrate import solve_ivp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'src'))
sys.path.insert(0, str(HERE))

OUT = HERE/'shock_reference.json'
NOTEBOOK = '14_magnus_supernova_shock.ipynb'
WIDTHS = (1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
N_FINGERPRINT = 40


def notebook_namespace():
    r"""Executes the notebook's own definition cells and returns their namespace.

    Stops at the cell defining ``measure``, which is the first one that wants the
    reference this script exists to produce.
    """
    import make_notebooks

    cells = make_notebooks.books[NOTEBOOK].cells
    ns = {'__name__': '__notebook__'}
    for cell in cells:
        if cell.cell_type != 'code':
            continue
        src = cell.source
        if 'def measure(' in src:
            break
        exec(compile(src, '<%s>' % NOTEBOOK, 'exec'), ns)
    for name in ('sn_shock_ne', 'make_H', 'L0', 'Ls'):
        if name not in ns:
            raise SystemExit(
                '%s did not define %r; the notebook has been restructured and this '
                'script needs updating' % (NOTEBOOK, name))
    return ns


def hexed(a):
    a = np.asarray(a, dtype=float)
    return [float(x).hex() for x in a.ravel()], list(a.shape)


def main():
    import matplotlib
    matplotlib.use('Agg')                     # the cells import pyplot

    ns = notebook_namespace()
    sn_shock_ne, make_H = ns['sn_shock_ne'], ns['make_H']
    L0, Ls = ns['L0'], np.asarray(ns['Ls'], dtype=float)
    fingerprint_l = np.linspace(float(L0), float(Ls[-1]), N_FINGERPRINT)

    store = {
        'note': ('solve_ivp DOP853 rtol=1e-12 atol=1e-14; produced by '
                 'notebooks/make_shock_reference.py, do not edit by hand'),
        'Ls': hexed(Ls)[0],
        'fingerprint_l': hexed(fingerprint_l)[0],
        'cases': {},
    }

    for width in WIDTHS:
        ne = sn_shock_ne(width)
        H = make_H(ne)

        def rhs(l, y, H=H):
            return (-1j*np.asarray(H(l)) @ y.reshape(3, 3)).ravel()

        sol = solve_ivp(rhs, (float(L0), float(Ls[-1])),
                        np.eye(3, dtype=complex).ravel(),
                        rtol=1e-12, atol=1e-14, method='DOP853', t_eval=Ls)
        if not sol.success:
            raise SystemExit('solve_ivp failed at w=%.0e: %s' % (width, sol.message))
        U = np.array([sol.y[:, i].reshape(3, 3) for i in range(len(Ls))])
        P = np.swapaxes(U.real**2 + U.imag**2, -1, -2)
        flat, shape = hexed(P)
        store['cases']['%.0e' % width] = {
            'P': flat, 'shape': shape,
            'fingerprint_ne': hexed(np.asarray(ne(fingerprint_l), dtype=float))[0],
        }
        print('  frozen w=%.0e  %s' % (width, tuple(shape)), flush=True)

    OUT.write_text(json.dumps(store, indent=1))
    print('wrote %s (%.0f kB)' % (OUT.name, OUT.stat().st_size/1024))


if __name__ == '__main__':
    main()
