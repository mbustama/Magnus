# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
r"""solarmodels.py

Standard solar models, tabulated by their authors, as profiles along a radial path.

The ``osc_prob_*_sun*`` wrappers describe the Sun by an exponential fit to its electron density
unless told otherwise.  Twelve published standard solar models ship with the package, and naming
one through ``density_profile`` puts that model's profile in its place: the electron density, from
the tabulated mass density and hydrogen mass fraction, and, for the wrappers with sterile states,
the neutron-to-proton ratio at every radius.

======================== ==================================================================
Name                     Model
======================== ==================================================================
``'BP2000'``             Bahcall, Pinsonneault & Basu (2001)
``'BP04'``               Bahcall & Pinsonneault (2004)
``'BS05-OP'``            Bahcall, Serenelli & Basu (2005), GS98 composition
``'BS05-AGS-OP'``        Bahcall, Serenelli & Basu (2005), AGS05 composition
``'B16-GS98'``           Vinyoles et al. (2017), GS98 composition
``'B16-AGSS09met'``      Vinyoles et al. (2017), AGSS09met composition
``'B23-GS98'``           Herrera & Serenelli (2023), GS98 composition
``'B23-AGSS09'``         Herrera & Serenelli (2023), AGSS09 composition
``'B23-C11'``            Herrera & Serenelli (2023), C11 composition
``'B23-AAG21'``          Herrera & Serenelli (2023), AAG21 composition
``'B23-MB22m'``          Herrera & Serenelli (2023), MB22 meteoritic composition
``'B23-MB22p'``          Herrera & Serenelli (2023), MB22 photospheric composition
======================== ==================================================================

Names are matched without regard to case.  :func:`solar_model_info` gives each model's full
reference, the source it was taken from, the date, its terms of use, and the radial range it
covers; the tables themselves, in ``magnus/data/solar_models/``, carry the same in their headers,
together with the authors' original header.  They hold three columns of the original -- radius,
mass density and hydrogen mass fraction -- copied as written.

**Outside the table.**  Only the B23 tables start at the centre; the others begin between
0.0005 and 0.0065 :math:`R_\odot`.  Some stop well short of the surface: BP2000 and BP04 end at
0.95 :math:`R_\odot`, the BS05 models at 0.98.  Below the first row the density is held at its
first value, since the core is flat.  Past the last row it
continues along the logarithmic slope of the last tabulated interval, so that it keeps falling
rather than stopping at a constant or dropping to zero.  That is a continuation, not a model of
the solar atmosphere: from 0.95 :math:`R_\odot` it reaches about :math:`1.6 \times 10^{-3}`
g cm\ :sup:`-3` at the surface, where the B16 and B23 models, which are tabulated there, give
:math:`1.7 \times 10^{-7}`.  The wrappers' ``stop_at_table_edge`` refuses a probability past the
last row instead.  The composition is held at its first and last tabulated values.

**Between rows.**  The density is interpolated linearly in its logarithm, so the profile has a
kink at every row, and the Bahcall tables, printed to four significant figures, step through
the flat core.  Phase-averaged probabilities do not notice.  A coherent probability over most
of the Sun can: the hybrid engine may then not certify it and hand it to the slab ladder,
which warns when it runs out of slabs.  The rows are where to put ``t_breakpoints`` (see
:doc:`/solar_models`).

Routine listings
----------------

    * available_solar_models - The names of the models that ship with the package
    * canonical_name - A model's name as the package spells it
    * solar_model_info - Reference, source, terms and radial range of one model
    * load_solar_model - The tabulated radius, mass density and hydrogen fraction of one model
    * table_edge - Distance from the centre of the last tabulated radius
    * electron_density_profile - Electron number density along a radial path
    * neutron_to_proton_ratio_profile - Neutron-to-proton ratio along a radial path

.. versionadded:: 1.2.0
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import functools
from importlib import resources
from typing import Callable, Dict, Tuple

import numpy as np

import magnus.globaldefs as gd
from magnus import _validate as _v


EXPONENTIAL = 'exp'
r"""str: The ``density_profile`` value that selects the exponential fit to the Sun's electron
density, which is the default of every Sun wrapper.

.. versionadded:: 1.2.0
"""

SOLAR_MODELS = ('BP2000', 'BP04', 'BS05-OP', 'BS05-AGS-OP', 'B16-GS98', 'B16-AGSS09met',
                'B23-GS98', 'B23-AGSS09', 'B23-C11', 'B23-AAG21', 'B23-MB22m', 'B23-MB22p')
r"""tuple of str: The standard solar models that ship with the package, oldest first.

.. versionadded:: 1.2.0
"""

_BY_KEY = {name.lower(): name for name in SOLAR_MODELS}

# Electrons and neutrons per atomic mass unit of the solar plasma (issue #168).  The tables give
# the hydrogen mass fraction X; the rest, 1 - X, is counted as helium-4, and the metals with it.
# With the atomic masses of 1H and 4He [u], hydrogen carries 1/1.00783 electrons per atomic mass
# unit and helium 2/4.00260, each helium atom as many neutrons as electrons.  Against the full
# composition of BS05(OP) -- 1H, 4He, 3He, 12C, 14N, 16O at their atomic masses -- this is within
# 0.002% at the centre and 0.07% at worst, near 0.3 R_sun where 3He (0.66 electrons per u) peaks;
# the metals, at about 0.496, move it by 0.01%.  The textbook rho N_A (1 + X)/2 is 0.4 to 0.8%
# higher, and the mean free-nucleon mass used until 1.2.0 put n_e 0.1 to 0.4% lower.
_MASS_H1_U = 1.00782503
_MASS_HE4_U = 4.00260325


def _per_atomic_mass_unit(X):
    """Electrons and neutrons per atomic mass unit for a hydrogen mass fraction X."""
    helium = (1.0 - X)*2.0/_MASS_HE4_U
    return X/_MASS_H1_U + helium, helium


def available_solar_models() -> Tuple[str, ...]:
    r"""The names of the standard solar models that ship with the package.

    .. versionadded:: 1.2.0

    Returns
    -------
    tuple of str
        The names accepted by ``density_profile``, oldest model first.

    Examples
    --------
    .. jupyter-execute::

        from magnus import solarmodels

        print(solarmodels.available_solar_models())
    """
    return SOLAR_MODELS


def canonical_name(name: str) -> str:
    r"""A solar model's name as the package spells it; the match ignores case.

    .. versionadded:: 1.2.0

    Parameters
    ----------
    name : str
        A model name, in any case, such as ``'b16-gs98'``.

    Returns
    -------
    str
        The name as listed in :data:`SOLAR_MODELS`.

    Raises
    ------
    ValueError
        If the name is not one of the shipped models; the message lists them.
    """
    # A string: str(None) and str(1) used to be looked up, and refused as unknown models
    # (issue #160 §3).
    if not isinstance(name, str):
        raise _v.InputTypeError(gd.ERROR_MSG_NO_COLOR + " solarmodels: a solar model is named by a "
            "string, such as " + repr(SOLAR_MODELS[0]) + "; got " + type(name).__name__ + " " +
            repr(name) + ".")
    key = name.strip().lower()
    if key not in _BY_KEY:
        raise ValueError(gd.ERROR_MSG_NO_COLOR + " solarmodels: unknown solar model " + repr(name) +
            ".  The models that ship with the package are " + ", ".join(SOLAR_MODELS) +
            "; the default, " + repr(EXPONENTIAL) + ", is the exponential fit.")
    return _BY_KEY[key]


def _file_name(name: str) -> str:
    return canonical_name(name).lower().replace('-', '_') + '.dat'


@functools.lru_cache(maxsize=None)
def _parse(name: str) -> Tuple[Dict[str, str], np.ndarray]:
    """The metadata fields and the numeric table of one model, read once."""
    text = (resources.files('magnus') / 'data' / 'solar_models' / _file_name(name)).read_text(
        encoding='utf-8')
    meta, rows = {}, []
    for line in text.splitlines():
        if line.startswith('#'):
            body = line[1:].strip()
            for key in ('Reference', 'Source', 'Original', 'Retrieved', 'SHA-256', 'Terms',
                        'density_column'):
                if body.startswith(key + ':') and key not in meta:
                    meta[key] = body[len(key) + 1:].strip()
        elif line.strip():
            rows.append([float(t) for t in line.split()])
    table = np.array(rows)
    if table.ndim != 2 or table.shape[1] != 3 or not np.all(np.diff(table[:, 0]) > 0.0):
        raise RuntimeError(gd.ERROR_MSG_NO_COLOR + " solarmodels: the table for " + name +
            " is damaged: it must hold three columns with the radius strictly increasing.")
    return meta, table


def load_solar_model(name: str) -> Dict[str, np.ndarray]:
    r"""The tabulated radius, mass density and hydrogen mass fraction of one standard solar model.

    .. versionadded:: 1.2.0

    Parameters
    ----------
    name : str
        One of :data:`SOLAR_MODELS`, in any case.

    Returns
    -------
    dict
        ``'r_over_r_sun'``, the radius in units of the solar radius; ``'rho_g_per_cm3'``, the
        mass density in g cm\ :sup:`-3`; and ``'x_hydrogen'``, the hydrogen mass fraction.  Each
        is an array with one entry per tabulated row, the radius strictly increasing.  Where the
        authors tabulate the logarithm of the density, it is exponentiated here.
    """
    name = canonical_name(name)
    meta, table = _parse(name)
    rho = 10.0**table[:, 1] if meta.get('density_column') == 'log10_rho_g_per_cm3' else table[:, 1]
    return {'r_over_r_sun': table[:, 0].copy(), 'rho_g_per_cm3': np.array(rho),
            'x_hydrogen': table[:, 2].copy()}


def solar_model_info(name: str) -> Dict[str, object]:
    r"""Reference, source, terms of use and radial range of one standard solar model.

    .. versionadded:: 1.2.0

    Parameters
    ----------
    name : str
        One of :data:`SOLAR_MODELS`, in any case.

    Returns
    -------
    dict
        ``'name'``; ``'reference'``, the paper or data release to cite; ``'source'``, where the
        table was taken from; ``'original'``, the file it was taken from; ``'retrieved'``;
        ``'sha256'`` of that original file; ``'terms'``, the authors' terms of use;
        ``'rows'``; and ``'r_min'`` and ``'r_max'``, the first and last tabulated radius in units
        of the solar radius.

    Examples
    --------
    .. jupyter-execute::

        from magnus import solarmodels

        info = solarmodels.solar_model_info('B16-GS98')
        print(info['reference'])
        print('tabulated from %.4f to %.4f R_sun' % (info['r_min'], info['r_max']))
    """
    name = canonical_name(name)
    meta, table = _parse(name)
    return {'name': name, 'reference': meta.get('Reference', ''), 'source': meta.get('Source', ''),
            'original': meta.get('Original', ''), 'retrieved': meta.get('Retrieved', ''),
            'sha256': meta.get('SHA-256', '').split()[0] if meta.get('SHA-256') else '',
            'terms': meta.get('Terms', ''), 'rows': int(table.shape[0]),
            'r_min': float(table[0, 0]), 'r_max': float(table[-1, 0])}


def table_edge(name: str) -> float:
    r"""Distance from the centre of the Sun to a model's last tabulated radius.

    .. versionadded:: 1.2.0

    Parameters
    ----------
    name : str
        One of :data:`SOLAR_MODELS`, in any case.

    Returns
    -------
    float
        The distance [:math:`\text{eV}^{-1}`], in the units of the wrappers' ``L``.
    """
    return solar_model_info(name)['r_max']*gd.SUN_RADIUS*gd.UNIT_KM


def _positions(name: str) -> np.ndarray:
    return load_solar_model(name)['r_over_r_sun']*gd.SUN_RADIUS*gd.UNIT_KM


def electron_density_profile(name: str) -> Callable:
    r"""Electron number density of a standard solar model along a radial path.

    The density is :math:`n_e = \rho N_A\,[X/m_{\rm H} + 2(1 - X)/m_{\rm He}]`, with :math:`X`
    the hydrogen mass fraction and :math:`m_{\rm H} = 1.00783`, :math:`m_{\rm He} = 4.00260`
    the atomic masses of :math:`^1`\ H and :math:`^4`\ He in atomic mass units; everything
    heavier than hydrogen is counted as helium-4.  Against a model's full composition this is
    within 0.07%; the textbook :math:`\rho N_A (1 + X)/2`, which Bahcall's electron-density
    tables follow, is 0.4 to 0.8% higher.  It is interpolated linearly in :math:`\ln n_e`.  Below the first
    tabulated radius it holds its first value; past the last, it continues along the logarithmic
    slope of the last interval (see the module notes).

    .. versionadded:: 1.2.0

    Parameters
    ----------
    name : str
        One of :data:`SOLAR_MODELS`, in any case.

    Returns
    -------
    Callable
        ``n_e(l)``, with ``l`` the distance from the centre [:math:`\text{eV}^{-1}`], scalar or
        array, returning the number density in the package's natural units -- the form the
        wrappers take with ``density_is_of_number_of_electrons=True``.

    Examples
    --------
    .. jupyter-execute::

        import magnus.globaldefs as gd
        from magnus import solarmodels

        ne = solarmodels.electron_density_profile('B16-GS98')
        per_cm3 = gd.N_AV*gd.UNIT_PER_CM3
        for r in (0.0, 0.5, 0.9):
            print('r = %.1f R_sun: n_e = %6.2f N_A per cm^3'
                  % (r, ne(r*gd.SUN_RADIUS*gd.UNIT_KM)/per_cm3))
    """
    model = load_solar_model(name)
    x = _positions(name)
    electrons, _ = _per_atomic_mass_unit(model['x_hydrogen'])
    log_ne = np.log(model['rho_g_per_cm3']*gd.UNIT_G_PER_CM3/gd.ATOMIC_MASS_UNIT*electrons)
    slope = (log_ne[-1] - log_ne[-2])/(x[-1] - x[-2])
    x_first, x_last, first, last = x[0], x[-1], log_ne[0], log_ne[-1]

    def ne(l):
        s = np.asarray(l, dtype=float)
        inside = np.interp(s, x, log_ne)
        out = np.exp(np.where(s > x_last, last + slope*(s - x_last),
                              np.where(s < x_first, first, inside)))
        return out[()] if np.ndim(out) == 0 else out

    return ne


def neutron_to_proton_ratio_profile(name: str) -> Callable:
    r"""Neutron-to-proton ratio of a standard solar model along a radial path.

    With the hydrogen mass fraction :math:`X` and everything heavier counted as helium-4, as in
    :func:`electron_density_profile`, :math:`n_n/n_p = [2(1 - X)/m_{\rm He}]/[X/m_{\rm H} +
    2(1 - X)/m_{\rm He}]`, the ratio the sterile states' matter term needs.  :math:`X` is
    interpolated linearly, and held at its first and last tabulated values outside the table.

    .. versionadded:: 1.2.0

    Parameters
    ----------
    name : str
        One of :data:`SOLAR_MODELS`, in any case.

    Returns
    -------
    Callable
        ``r(l)``, with ``l`` the distance from the centre [:math:`\text{eV}^{-1}`], scalar or
        array.
    """
    model = load_solar_model(name)
    x, X = _positions(name), model['x_hydrogen']

    def ratio(l):
        electrons, neutrons = _per_atomic_mass_unit(np.interp(np.asarray(l, dtype=float), x, X))
        out = neutrons/electrons
        return out[()] if np.ndim(out) == 0 else out

    return ratio


__all__ = [
    'EXPONENTIAL',
    'SOLAR_MODELS',
    'available_solar_models',
    'canonical_name',
    'solar_model_info',
    'load_solar_model',
    'table_edge',
    'electron_density_profile',
    'neutron_to_proton_ratio_profile',
]
