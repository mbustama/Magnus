# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
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

It is also the one place the ``angles`` convention is interpreted.  Callers may state
their mixing angles as ``'sin'`` (the default, and what :func:`magnus.globaldefs.
load_nufit_params` returns), ``'sin2'`` (the sine *squared*, which is what global fits
report), ``'rad'`` or ``'deg'``; everything downstream of :func:`resolve` sees sines and
radians and needs to know nothing about it.

Routine listings
----------------

    * validate_sines - Raise if any mixing-angle sine is out of range
    * validate_convention - Raise unless ``angles`` is one of the four accepted values
    * resolve - Convert mixing angles to sines and CP phases to radians
    * from_sines - The inverse: state stored sines in the caller's convention
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import numpy as np


def validate_sines(source_func_name: str, **sines) -> None:
    r"""Raises :class:`ValueError` if any named sine is outside :math:`[-1, 1]`.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    source_func_name : str
        Module-qualified name of the calling function, for the message -- ``resolve`` is
        reached from ``oscprobstd`` as well as from the builders, so the module is the
        caller's to state rather than this module's to assume.
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
            gd.ERROR_MSG_NO_COLOR + " " + source_func_name + ": " + name
            + " is the sine of a mixing angle and must lie in [-1, 1]; got "
            + repr(bad) + ". If several angles were passed positionally, check the"
            " order: these signatures interleave each angle with its CP phase, so"
            " grouping the angles together silently puts a phase in a sine slot.")


# theta_13 is the smallest mixing angle anyone measures, at about 8.5 degrees.  A
# parameter set whose every angle is under a degree is therefore not a small-angle
# study, it is sines that were handed to `angles='deg'` -- sin(theta) for the known
# angles runs 0.15 to 0.85, which read as degrees is fifty times too small.
#
# A warning rather than an error, for the same reason as matter.py's density pair:
# the threshold sits far below anything currently modelled, but that is a statement
# about the mixing people currently study, not a law.
IMPLAUSIBLE_MIXING_ANGLE_DEG = 1.0
r"""float: Below this, in degrees, a whole parameter set is probably sines. [deg]"""


def validate_convention(source_func_name: str, angles: str) -> str:
    r"""Raises :class:`ValueError` unless ``angles`` is one of
    :data:`magnus.globaldefs.ANGLE_CONVENTIONS`.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    source_func_name : str
        Name of the calling function, for the message.
    angles : str
        The convention to check.

    Returns
    -------
    str
        ``angles`` unchanged, so callers can write ``angles = validate_convention(...)``.

    Raises
    ------
    ValueError
        If ``angles`` is not one of the four accepted values.
    """
    import magnus.globaldefs as gd
    if angles in gd.ANGLE_CONVENTIONS:
        return angles
    raise ValueError(
        gd.ERROR_MSG_NO_COLOR + " " + source_func_name + ": angles must be one of "
        + ", ".join(repr(c) for c in gd.ANGLE_CONVENTIONS) + "; got " + repr(angles)
        + ".  'sin' (the default) is the sine of the angle, 'sin2' is the sine"
        " SQUARED -- which is what global fits report -- 'rad' is the angle itself"
        " in radians, and 'deg' the angle in degrees.")


def _warn_if_angles_are_probably_sines(source_func_name, values):
    r"""Warns when every angle declared to be in degrees is too small to be one.

    Called only for ``angles='deg'``.  See :data:`IMPLAUSIBLE_MIXING_ANGLE_DEG`.

    .. versionadded:: 1.0.0
    """
    import warnings

    import magnus.globaldefs as gd

    finite = [float(np.max(np.abs(np.asarray(v, dtype=float))))
              for v in values if np.all(np.isfinite(np.asarray(v, dtype=float)))]
    if not finite:
        return
    largest = max(finite)
    if largest == 0.0 or largest >= IMPLAUSIBLE_MIXING_ANGLE_DEG:
        return

    warnings.warn(
        gd.WARNING_MSG_NO_COLOR + " " + source_func_name + ": angles='deg' was given, but"
        " the largest mixing angle passed is " + format(largest, '.4g') + " degrees, and"
        " theta_13 -- the smallest angle anyone measures -- is about 8.5.  These are most"
        " likely sines: sin(theta) for the known angles runs 0.15 to 0.85, which read as"
        " degrees is about fifty times too small.  The call will return a converged,"
        " unitary, entirely wrong probability rather than an error.  Either drop"
        " angles='deg' (its default, 'sin', is what load_nufit_params returns), or pass"
        " the angles themselves.",
        # 4, not matter.py's 3: this chain is one frame deeper -- warn, this function,
        # resolve, the builder -- so 4 is what attributes it to the builder's caller.
        gd.MixingAngleConventionWarning, stacklevel=4)


def resolve(source_func_name: str, angles: str, sines: dict, phases: dict = None):
    r"""Converts mixing angles to sines and CP phases to radians.

    The one place any of the four conventions is interpreted.  Written here rather than
    inlined in each builder for the reason given at the top of this module: the same
    structure copied into several files, all agreeing with each other, is how this
    package's four-flavour matter term stayed wrong through a max-effort review.

    **Phases follow the convention only for** ``'deg'``.  A CP phase has no sine
    representation, so under ``'sin'`` and ``'sin2'`` it stays in radians; under
    ``'rad'`` radians are already what it is in.  Only ``'deg'`` converts it.

    **``'sin2'`` cannot express a negative sine**, since the square discards the sign and
    the root returns the non-negative branch.  ``'sin'`` accepts the full :math:`[-1, 1]`,
    and maximal mixing at :math:`s = -1` is a case the suite pins, so a caller who needs a
    negative sine has to state it in ``'sin'``, ``'rad'`` or ``'deg'``.  A negative value
    under ``'sin2'`` is therefore rejected rather than squared away: it means the caller
    believes they are in ``'sin'``.

    The ``'sin'`` path returns the values it was given, unconverted -- it is already the
    internal convention -- so the default costs nothing and is bit-identical to not having
    this function at all.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    source_func_name : str
        Name of the calling function, for any message.
    angles : str
        One of :data:`magnus.globaldefs.ANGLE_CONVENTIONS`.
    sines : dict
        Mixing-angle parameters, by name, in the ``angles`` convention.
    phases : dict, optional
        CP-phase parameters, by name, in the ``angles`` convention.

    Returns
    -------
    tuple of dict
        ``(sines, phases)``, the first as sines of the angles and the second in radians,
        with the same keys they were given.

    Raises
    ------
    ValueError
        If ``angles`` is not accepted, or a value is out of range for its convention.
    """
    validate_convention(source_func_name, angles)
    phases = {} if phases is None else phases

    if angles == 'sin':
        validate_sines(source_func_name, **sines)
        return dict(sines), dict(phases)

    if angles == 'sin2':
        # sin^2 cannot be negative, and a negative here is not a small slip: it means the
        # caller believes they are in 'sin', where a negative sine is perfectly legal.
        _validate_range(source_func_name, sines, 0.0, 1.0, 'sin2',
                        "the SQUARE of the sine of a mixing angle, so it must lie in"
                        " [0, 1]; a negative value means the convention is 'sin' rather"
                        " than 'sin2'")
        return ({name: np.sqrt(np.asarray(value, dtype=float))
                 for name, value in sines.items()}, dict(phases))

    if angles == 'rad':
        # 2*pi, not pi/2: a mixing angle outside the first quadrant is unconventional
        # rather than wrong, but a value above 2*pi is degrees in a radians slot.
        _validate_range(source_func_name, sines, -2.0*np.pi, 2.0*np.pi, 'rad',
                        "an angle in radians and must lie in [-2*pi, 2*pi]; a larger"
                        " value is degrees passed as radians")
        converted = {name: np.sin(np.asarray(value, dtype=float))
                     for name, value in sines.items()}
        return converted, dict(phases)

    _validate_range(source_func_name, sines, -360.0, 360.0, 'deg',
                    "an angle in degrees and must lie in [-360, 360]")
    _warn_if_angles_are_probably_sines(source_func_name, list(sines.values()))
    return ({name: np.sin(np.radians(np.asarray(value, dtype=float)))
             for name, value in sines.items()},
            {name: np.radians(np.asarray(value, dtype=float))
             for name, value in phases.items()})


def from_sines(angles: str, sines: dict, phases: dict = None):
    r"""The inverse of :func:`resolve`: states stored sines in the caller's convention.

    Needed wherever the package supplies values the caller did not: the predefined parameter
    sets and the NuFit tables are stored as sines, and handing a sine back to someone working
    in degrees produces a parameter set in **two** conventions at once.  Omitting one angle
    from a ``angles='deg'`` call would otherwise fill it from the defaults as a sine, and the
    builder would then read 0.1499 as 0.1499 degrees -- a converged, unitary, wrong answer,
    and only for the parameters the caller happened not to pass.

    One function rather than a conversion at each site, for the reason at the top of this
    module: two copies of the same arithmetic agreeing with each other is how a wrong
    four-flavour matter term survived a max-effort review.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    angles : str
        One of :data:`magnus.globaldefs.ANGLE_CONVENTIONS`.
    sines : dict
        Mixing-angle sines, by name, as stored.
    phases : dict, optional
        CP phases in radians, by name, as stored.

    Returns
    -------
    tuple of dict
        ``(angles, phases)`` restated in the requested convention.
    """
    phases = {} if phases is None else phases
    if angles == 'sin':
        return dict(sines), dict(phases)
    if angles == 'sin2':
        return ({k: np.asarray(v, dtype=float)**2 for k, v in sines.items()}, dict(phases))
    theta = {k: np.arcsin(np.asarray(v, dtype=float)) for k, v in sines.items()}
    if angles == 'rad':
        return theta, dict(phases)
    return ({k: np.degrees(v) for k, v in theta.items()},
            {k: np.degrees(np.asarray(v, dtype=float)) for k, v in phases.items()})


def _validate_range(source_func_name, values, low, high, convention, what):
    r"""Raises :class:`ValueError` if any named value falls outside ``[low, high]``."""
    for name, value in values.items():
        arr = np.asarray(value, dtype=float)
        if np.all(np.isfinite(arr)) and np.all(arr >= low) and np.all(arr <= high):
            continue
        import magnus.globaldefs as gd
        flat = arr.ravel()
        bad = float(arr) if arr.ndim == 0 else float(
            flat[np.argmax(~np.isfinite(flat) | (flat < low) | (flat > high))])
        raise ValueError(
            gd.ERROR_MSG_NO_COLOR + " " + source_func_name + ": with angles="
            + repr(convention) + ", " + name + " is " + what + "; got " + repr(bad) + ".")


__all__ = [
    'IMPLAUSIBLE_MIXING_ANGLE_DEG',
    'from_sines',
    'resolve',
    'validate_convention',
    'validate_sines',
]
