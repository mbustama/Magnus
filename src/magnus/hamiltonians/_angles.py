# -*- coding: utf-8 -*-
r"""_angles.py

One guard, shared by every vacuum-Hamiltonian builder: a mixing-angle sine has to lie
in :math:`[-1, 1]`.

Every builder computes its cosines as :math:`\sqrt{1 - s^2}`.  Handed a sine outside
the unit interval, NumPy returns ``nan`` with only a ``RuntimeWarning``, and the
builder goes on to return a Hamiltonian full of ``nan`` -- which ``osc_prob`` then
propagates into ``nan`` probabilities.  Nothing raises, and a caller who is not
watching stderr sees a result-shaped object that is not a result.

**The mistake this catches is a slot error, not a physics error.**  Nobody types
``s25 = 3.79`` on purpose.  They pass eighteen positional arguments to
``hamiltonian_5nu_vacuum_energy_independent``, whose signature interleaves each angle
with its CP phase --- ``s14, d14, s15, d15, s24, d24, s25, s34, s35, d35`` --- group
the angles together as any reader would expect, and a *phase* lands in a sine slot.
That is exactly how this was found.

Written once here rather than inlined in each of the four builders.  The same structure
copied into several files, all agreeing with each other, is how this package's
four-flavour matter term stayed wrong through a max-effort review; see
:func:`magnus.matter.matter_potential_projector`.

.. versionadded:: 1.0.0

Routine listings
----------------

    * validate_sines - Raise if any mixing-angle sine is out of range
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import numpy as np


def validate_sines(source_func_name: str, **sines) -> None:
    r"""Raises :class:`ValueError` if any named sine is outside :math:`[-1, 1]`.

    Parameters
    ----------
    source_func_name : str
        Name of the calling builder, for the message.
    \**sines
        The mixing-angle sines, by parameter name, so the message can say which one is
        at fault rather than that one of eighteen arguments is.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        Naming the offending parameter and its value.
    """
    for name, value in sines.items():
        arr = np.asarray(value, dtype=float)
        if np.all(np.isfinite(arr)) and np.all(np.abs(arr) <= 1.0):
            continue
        # Imported here, not at module scope: magnus.globaldefs imports magnus.magnus,
        # and this module is imported by the builders, so a top-level import would put
        # a needless cycle in the package's import graph for a string.
        import magnus.globaldefs as gd
        bad = float(arr) if arr.ndim == 0 else float(
            arr.ravel()[np.argmax(~np.isfinite(arr.ravel())
                                  | (np.abs(arr.ravel()) > 1.0))])
        raise ValueError(
            gd.ERROR_MSG_NO_COLOR + " hamiltonians." + source_func_name + ": " + name
            + " is the sine of a mixing angle and must lie in [-1, 1]; got "
            + repr(bad) + ". If several angles were passed positionally, check the"
            " order: these signatures interleave each angle with its CP phase, so"
            " grouping the angles together silently puts a phase in a sine slot.")


__all__ = [
    'validate_sines',
]
