r"""Broadcasting helpers shared by the Hamiltonian builders (issue #155 §2).

Every builder takes its one position- or energy-dependent argument -- ``energy``, ``VCC`` or
``l`` -- as a number or as an array, the way the matter builders always have: a number returns
one ``(d, d)`` matrix, and an array of shape ``s`` returns a stack of shape ``s + (d, d)``.  That
is what lets a user's ``H_func`` built from them take the engine's vectorized path (see
:class:`magnus.magnus.ScalarHamiltonianWarning`).

Routine listings
----------------

    * stacked - An argument as a float array with two trailing axes
    * over_positions - A position-independent matrix repeated once per position
"""

import numpy as np

#: The scalar types a builder meets at every quadrature node; a lookup of the type here skips
#: ``np.ndim`` on the hot path, so scalar calls cost what they did before arrays were accepted.
SCALARS = frozenset((float, int, np.float64, np.float32, np.int64, np.int32))


def stacked(x):
    r"""``x`` as a float array with two trailing axes, ready to multiply a ``(d, d)`` matrix."""
    return np.asarray(x, dtype=float)[..., None, None]


def over_positions(l, H):
    r"""``H`` repeated once per position in ``l``, for a builder whose matrix ignores ``l``.

    A scalar ``l`` returns ``H`` itself.  An array returns a writable stack of shape
    ``broadcast(shape(l), H.shape[:-2]) + (d, d)``, so an energy stack and a position stack of
    the same length pair up entry by entry.
    """
    if type(l) in SCALARS or np.ndim(l) == 0:
        return H
    shape = np.broadcast_shapes(np.shape(l), np.shape(H)[:-2]) + np.shape(H)[-2:]
    return np.broadcast_to(H, shape).copy()
