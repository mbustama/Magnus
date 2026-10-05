# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
r"""cli.py

Command-line calculator for Mag$\nu$s: computes a single neutrino
oscillation probability (or probability matrix) from the command line,
without writing any Python. Wraps the same ``osc_prob_{2,3,4,5}nu_*``
functions used by the Python API (see :py:mod:`magnus.oscprob`
and :doc:`/cli`), dispatching to the right one based on ``--flavors``,
``--environment``, ``--scenario`` and ``--density-profile``.

Installed as the ``magnus`` console script (see ``[project.scripts]``
in pyproject.toml) and also runnable as ``python -m magnus``.  Its one
subcommand, ``prob``, is the default: ``magnus --flavors 3 ...`` runs
``magnus prob --flavors 3 ...`` (issue #138).

Routine listings
----------------

    * main - Entry point: parses argv, dispatches, prints the result
    * build_parser - Builds the argparse.ArgumentParser
    * SUBCOMMANDS - The subcommands build_parser registers; the first is the default
    * FLAVOR_NAME_TO_INDEX - Maps flavor names (e, mu, tau, s, s1, s2)
           to their globaldefs index
    * ENERGY_UNITS, LENGTH_UNITS - Unit-name to :math:`\text{eV}` / :math:`\text{eV}^{-1}`
           conversion factors
    * FLAVOR_LABELS - Row and column names used by the printed table
    * ALWAYS_FORWARD - Numerics keywords forwarded to every wrapper
    * DENSITY_PROFILES - Values ``--density-profile`` accepts: 'constant', 'exp', and
           every tabulated solar model in :mod:`magnus.solarmodels`
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import argparse
import math
import inspect
import json
import sys
import warnings
from typing import Optional

import numpy as np

import magnus.oscprob as oscprob
import magnus.globaldefs as gd
import magnus.solarmodels as solarmodels
from magnus.version import __version__


ENERGY_UNITS = {
    'eV': 1.0, 'keV': gd.UNIT_KEV, 'MeV': gd.UNIT_MEV,
    'GeV': gd.UNIT_GEV, 'TeV': gd.UNIT_TEV, 'PeV': gd.UNIT_PEV,
}

LENGTH_UNITS = {
    'eV-1': 1.0, 'km': gd.UNIT_KM, 'cm': gd.UNIT_CM,
}

FLAVOR_NAME_TO_INDEX = {
    'e': gd.NUE, 'mu': gd.NUMU, 'tau': gd.NUTAU,
    's': gd.NUS, 's1': gd.NUS1, 's2': gd.NUS2,
}

FLAVOR_LABELS = {
    2: ['0', '1'],
    3: ['nu_e', 'nu_mu', 'nu_tau'],
    4: ['nu_e', 'nu_mu', 'nu_tau', 'nu_s'],
    5: ['nu_e', 'nu_mu', 'nu_tau', 'nu_s1', 'nu_s2'],
}

# Refinement/numerics kwargs that every osc_prob_* wrapper accepts via **kwargs
# even where they are not explicit named parameters (see the "layer contract" in
# docs/source/architecture.rst).  One exception: 'strategy' is only taken by the
# position-dependent wrappers, so main drops it for vacuum and constant density.
# No logging keyword belongs here -- verbose reaches the wrappers because every
# one of them declares it.
ALWAYS_FORWARD = {'magnus_exp_order', 'n_jobs', 'integration_method', 'rtol', 'atol',
                  'strategy'}


# 'constant' and 'exp' describe a profile by its shape; the rest name a tabulated standard
# solar model (see magnus.solarmodels), and so apply to --environment sun only.
DENSITY_PROFILES = ('constant', 'exp') + solarmodels.SOLAR_MODELS


def _density_profile(value: str) -> str:
    r"""argparse type= callback: the canonical spelling of a density profile, in any case.

    A value it does not recognize is returned unchanged, so that argparse's own ``choices``
    check rejects it and lists what is accepted.
    """
    key = value.strip().lower()
    for name in DENSITY_PROFILES:
        if name.lower() == key:
            return name
    return value


def _flavor_index(value: str) -> int:
    r"""argparse type= callback: accepts an int string or a flavor name."""
    try:
        return int(value)
    except ValueError:
        pass
    key = value.strip().lower()
    if key not in FLAVOR_NAME_TO_INDEX:
        raise argparse.ArgumentTypeError(
            f"invalid flavor {value!r}; expected an integer index or one of "
            f"{sorted(FLAVOR_NAME_TO_INDEX)}")
    return FLAVOR_NAME_TO_INDEX[key]


def _finite_float(value: str) -> float:
    r"""argparse type= callback: a finite number.

    ``float()`` accepts 'nan' and 'inf', and a NaN baseline or phase used to print a table of
    NaN with exit code 0 (issue #160 §13).
    """
    try:
        x = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"invalid number {value!r}") from None
    if not math.isfinite(x):
        raise argparse.ArgumentTypeError(f"must be a finite number, not {value!r}")
    return x


def _precision(value: str) -> int:
    r"""argparse type= callback: a digit count from 0 to 17, what a float64 can show.

    -1 used to print the header and then a raw traceback; 100 printed digits beyond float64.
    """
    try:
        n = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"invalid int value {value!r}") from None
    if not 0 <= n <= 17:
        raise argparse.ArgumentTypeError(f"must be from 0 to 17 (the digits a double carries), "
                                         f"not {n}")
    return n


SUBCOMMANDS = ('prob',)
r"""tuple of str: Module-level constant

The subcommands :func:`build_parser` registers.  The first, ``prob``, is the default:
a command line that names none runs it (issue #138).

.. versionadded:: 1.2.0
"""


def _with_default_subcommand(argv):
    r"""Prepend the default subcommand unless ``argv`` names one or asks for top-level help.

    ``magnus --flavors 3 ...`` becomes ``magnus prob --flavors 3 ...``.  An empty command
    line, ``-h``/``--help`` and ``-V``/``--version`` are left to the top-level parser.

    .. versionadded:: 1.2.0

    Parameters
    ----------
    argv : list of str
        The command-line arguments, without the program name.

    Returns
    -------
    list of str
        ``argv``, with ``'prob'`` in front where no subcommand was given.
    """
    argv = list(argv)
    if argv and argv[0] not in SUBCOMMANDS + ('-h', '--help', '-V', '--version'):
        argv.insert(0, SUBCOMMANDS[0])
    return argv


class _HelpWithDefaultSubcommand(argparse.Action):
    r"""``magnus -h``: the top-level help, followed by the help of the default subcommand.

    ``prob`` runs when no subcommand is named, so its options are the ones a user of
    ``magnus --help`` is looking for; printing only the subcommand list would send them to
    ``magnus prob --help`` for every option (issue #138).

    .. versionadded:: 1.2.0
    """

    def __call__(self, parser, namespace, values, option_string=None):
        parser.print_help()
        print()
        parser.default_subparser.print_help()
        parser.exit()


def build_parser() -> argparse.ArgumentParser:
    r"""Builds the ``magnus`` command-line argument parser.

    .. versionadded:: 1.0.0

    Returns
    -------
    argparse.ArgumentParser
        The top-level parser, with the ``prob`` subcommand attached (also kept as its
        ``default_subparser`` attribute).  Its epilog says that ``prob`` is the default
        subcommand, which :func:`main` supplies when none is given, and its ``-h`` prints
        the help of ``prob`` after its own.

    .. versionchanged:: 1.2.0
       The epilog names ``prob`` as the default subcommand, and ``magnus -h`` also prints
       the options of ``prob`` (issue #138).

    .. versionchanged:: 1.2.0
       The --magnus-exp-order help states that 'gl' takes the even orders from 2 to 8 (issue
       #160 §5).
    """
    parser = argparse.ArgumentParser(
        prog='magnus',
        description="Magνs: neutrino oscillation probabilities via the Magnus expansion.",
        epilog="prob is the default subcommand: 'magnus ...' runs 'magnus prob ...', and its "
               "options follow.",
        add_help=False)
    parser.add_argument('-h', '--help', action=_HelpWithDefaultSubcommand, nargs=0,
        help="show this help message, followed by that of 'prob', and exit")
    parser.add_argument('-V', '--version', action='version', version=f'magnus {__version__}')
    sub = parser.add_subparsers(dest='command', required=True)

    p = sub.add_parser('prob', help='Compute a single oscillation probability (matrix or channel).')
    parser.default_subparser = p

    g_env = p.add_argument_group('Environment')
    g_env.add_argument('--flavors', type=int, choices=[2, 3, 4, 5], default=3,
        help='Number of neutrino flavors (default: 3).')
    g_env.add_argument('--environment', choices=['vacuum', 'matter', 'earth', 'sun'], default='vacuum',
        help='Propagation environment (default: vacuum).')
    g_env.add_argument('--scenario', choices=['std', 'nsi', 'liv'], default='std',
        help="Physics scenario on top of the environment: 'std' (Standard Model), "
             "'nsi' (non-standard interactions) or 'liv' (Lorentz-invariance violation). "
             "'nsi' is not available with --environment vacuum. Default: std.")
    g_env.add_argument('--density-profile', type=_density_profile, choices=DENSITY_PROFILES,
        default=None, metavar='PROFILE',
        help="Matter density profile.  With --environment matter: 'constant' (the default; "
             "requires --rho) or 'exp' (requires --rho-central and --l-scale).  With "
             "--environment sun: 'exp' (the default) or a tabulated standard solar model, "
             "named in any case: " + ', '.join(solarmodels.SOLAR_MODELS) + '.')
    g_env.add_argument('--nubar', action='store_true',
        help='Compute the probability for antineutrinos instead of neutrinos. No effect '
             'with --flavors 2 --environment vacuum --scenario std, where there is no CP '
             'phase and no matter, so the two probabilities are equal.')

    g_kin = p.add_argument_group('Energy and baseline')
    g_kin.add_argument('--energy', type=_finite_float, required=True, help='Neutrino energy.')
    g_kin.add_argument('--energy-unit', choices=list(ENERGY_UNITS), default='GeV',
        help='Unit of --energy (default: GeV).')
    g_kin.add_argument('--baseline', type=_finite_float, default=None,
        help='Baseline (final position). Required for vacuum, matter and sun. For earth, '
             'give it with --costhz, or omit it when --detector-depth or --source-depth is '
             'given (the baseline is then computed) or when --loc-ini and --loc-fin are '
             'given. Must be omitted with --detector-depth.')
    g_kin.add_argument('--l0', type=_finite_float, default=0.0,
        help='Initial position (used by --environment sun and --density-profile exp). '
             'Default: 0.0.')
    g_kin.add_argument('--baseline-unit', choices=list(LENGTH_UNITS), default='km',
        help='Unit of --baseline, --l0, --l-scale, --source-depth and --detector-depth '
             '(default: km).')

    g_mat = p.add_argument_group('Matter (--environment matter)')
    g_mat.add_argument('--rho', type=_finite_float, default=None,
        help='Matter density (constant profile).')
    g_mat.add_argument('--rho-central', type=_finite_float, default=None,
        help='Matter density at the center of the profile, l=0 (exponential profile).')
    g_mat.add_argument('--l-scale', type=_finite_float, default=None,
        help='Length scale of the exponential density decrease (exponential profile).')
    g_mat.add_argument('--density-unit', choices=['g/cm3', 'natural'], default='g/cm3',
        help='Unit of --rho/--rho-central: g/cm3 (converted internally) or natural units '
             '(eV^4). Default: g/cm3.')
    g_mat.add_argument('--ratio-n-to-p', type=_finite_float, default=1.0,
        help='Ratio of the number of neutrons to protons in matter. It sets the '
             'neutral-current potential of the sterile states, so it has no effect at 2 or '
             '3 flavors. Default: 1.0.')
    g_mat.add_argument('--electron-fraction', type=_finite_float, default=0.5,
        help='Y_e, the number of electrons per atomic mass unit of the matter, so that '
             'n_e = rho N_A Y_e. Default: 0.5.')

    g_earth = p.add_argument_group('Earth (--environment earth)')
    g_earth.add_argument('--costhz', type=_finite_float, default=None,
        help='Cosine of the neutrino zenith angle.')
    g_earth.add_argument('--loc-ini', default=None,
        help='Initial location name (e.g., fermilab); see magnus.earth.loc_coords_dms. '
             'Must be given together with --loc-fin, as an alternative to --costhz.')
    g_earth.add_argument('--loc-fin', default=None,
        help='Final location name; see --loc-ini.')
    # Both depths are read in --baseline-unit, so they are stated in the same unit as the
    # baseline they replace.  Naming --detector-depth makes --baseline unnecessary rather
    # than optional: the library raises if both arrive, so the branch below stops it here
    # with a message that names the flags rather than the parameters.
    g_earth.add_argument('--detector-depth', type=_finite_float, default=0.0,
        help='Depth of the detector below the surface, in --baseline-unit. The zenith '
             'angle is measured at the detector, so a buried one also sees downward-going '
             'neutrinos (--costhz > 0) through its overburden. Computes the baseline, so '
             '--baseline must be omitted. Default: 0 (a detector on the surface).')
    g_earth.add_argument('--source-depth', type=_finite_float, default=0.0,
        help="Depth of the neutrino's entry point below the surface, in --baseline-unit. "
             'Default: 0 (entry at the surface).')

    g_sun = p.add_argument_group('Sun (--environment sun)')
    g_sun.add_argument('--stop-at-table-edge', action='store_true',
        help='With a tabulated solar model, do not extrapolate past its last row: a baseline '
             'that ends beyond it returns nan, with a warning.  Without it, the density '
             "continues the last row's logarithmic slope.  Not available with 'exp', which "
             'has no table.')

    g_osc = p.add_argument_group('Standard oscillation parameters (2-flavor)')
    g_osc.add_argument('--angles', default='sin', choices=list(gd.ANGLE_CONVENTIONS),
        help="Convention for every mixing angle below (--sth, --s12, ..., --sxi, ...): 'sin' "
             "(default) takes sines, 'sin2' squared sines (the form global fits report), 'rad' "
             "radians and 'deg' degrees. Under 'deg' the CP phases (--dcp, --d14, ...) are "
             "read as degrees too; otherwise they stay in radians.")
    g_osc.add_argument('--sth', type=_finite_float, default=None,
        help='Mixing angle theta, in the convention set by --angles (required for --flavors 2).')
    g_osc.add_argument('--dm2', type=_finite_float, default=None, dest='Dm2',
        help='Mass-squared difference Delta m^2 [eV^2] (required for --flavors 2).')

    g_osc3 = p.add_argument_group('Standard oscillation parameters (3+ flavors)')
    g_osc3.add_argument('--s12', type=_finite_float, default=None, help='Mixing angle theta_12, per --angles. Default: NuFIT 6.1.')
    g_osc3.add_argument('--s23', type=_finite_float, default=None, help='Mixing angle theta_23, per --angles. Default: NuFIT 6.1.')
    g_osc3.add_argument('--s13', type=_finite_float, default=None, help='Mixing angle theta_13, per --angles. Default: NuFIT 6.1.')
    g_osc3.add_argument('--dcp', type=_finite_float, default=None, dest='dCP',
        help='delta_CP [radian, or degree with --angles deg]. Default: NuFIT 6.1.')
    g_osc3.add_argument('--dm21', type=_finite_float, default=None, dest='D21',
        help='Mass-squared difference Delta m^2_21 [eV^2]. Default: NuFIT 6.1.')
    g_osc3.add_argument('--dm31', type=_finite_float, default=None, dest='D31',
        help='Mass-squared difference Delta m^2_31 [eV^2]. Default: NuFIT 6.1.')
    g_osc3.add_argument('--osc-params-set', default='OSC_PARAMS_DEFAULT',
        dest='default_osc_params_set_name',
        choices=sorted(gd.OSC_PARAMS_PREDEFINED),
        # `metavar` hides the enumeration without weakening it: argparse still rejects a
        # name that is not in `choices`, and still prints the full list in the error it
        # raises.  What it stops is printing all of them in the usage line and again in
        # the option's own entry.  With five sets that was informative; with one per NuFIT
        # release, ordering and SK variant it is some 2500 characters of help, twice over.
        metavar='NAME',
        help='Predefined set used to fill in any of s12/s23/s13/dCP/D21/D31 left unspecified: '
             'one per NuFIT release, mass ordering (..._NO or ..._IO) and, from release 4.0 on, '
             'inclusion (..._SK_) or exclusion (..._NOSK_) of Super-Kamiokande atmospheric '
             'data.  OSC_PARAMS_DEFAULT is NuFIT 6.1 SK NO.  The full list is '
             'globaldefs.OSC_PARAMS_PREDEFINED; an unknown name prints it.')

    g_osc4 = p.add_argument_group('Additional sterile mixing (4+ flavors)')
    g_osc4.add_argument('--s14', type=_finite_float, default=0.0, help='Mixing angle theta_14, per --angles. Default: 0.0.')
    g_osc4.add_argument('--d14', type=_finite_float, default=0.0, help='delta_14 [radian, or degree with --angles deg]. Default: 0.0.')
    g_osc4.add_argument('--s24', type=_finite_float, default=0.0, help='Mixing angle theta_24, per --angles. Default: 0.0.')
    g_osc4.add_argument('--d24', type=_finite_float, default=0.0, help='delta_24 [radian, or degree with --angles deg]. Default: 0.0.')
    g_osc4.add_argument('--s34', type=_finite_float, default=0.0, help='Mixing angle theta_34, per --angles. Default: 0.0.')
    g_osc4.add_argument('--dm41', type=_finite_float, default=0.0, dest='D41',
        help='Mass-squared difference Delta m^2_41 [eV^2]. Default: 0.0.')

    g_osc5 = p.add_argument_group('Additional sterile mixing (5 flavors)')
    g_osc5.add_argument('--s15', type=_finite_float, default=0.0, help='Mixing angle theta_15, per --angles. Default: 0.0.')
    g_osc5.add_argument('--d15', type=_finite_float, default=0.0, help='delta_15 [radian, or degree with --angles deg]. Default: 0.0.')
    g_osc5.add_argument('--s25', type=_finite_float, default=0.0, help='Mixing angle theta_25, per --angles. Default: 0.0.')
    g_osc5.add_argument('--s35', type=_finite_float, default=0.0, help='Mixing angle theta_35, per --angles. Default: 0.0.')
    g_osc5.add_argument('--d35', type=_finite_float, default=0.0, help='delta_35 [radian, or degree with --angles deg]. Default: 0.0.')
    g_osc5.add_argument('--dm51', type=_finite_float, default=0.0, dest='D51',
        help='Mass-squared difference Delta m^2_51 [eV^2]. Default: 0.0.')

    g_nsi = p.add_argument_group('NSI parameters (--scenario nsi)')
    g_nsi.add_argument('--eps-aa', type=_finite_float, default=0.0, help='2-flavor diagonal NSI coupling.')
    g_nsi.add_argument('--eps-ab', type=_finite_float, default=0.0, help='2-flavor off-diagonal NSI coupling.')
    g_nsi.add_argument('--eps-ee', type=_finite_float, default=0.0, help='Diagonal NSI coupling of nu_e.')
    g_nsi.add_argument('--eps-em', type=_finite_float, default=0.0, help='Off-diagonal (e-mu) NSI coupling.')
    g_nsi.add_argument('--eps-et', type=_finite_float, default=0.0, help='Off-diagonal (e-tau) NSI coupling.')
    g_nsi.add_argument('--eps-mm', type=_finite_float, default=0.0, help='Diagonal NSI coupling of nu_mu.')
    g_nsi.add_argument('--eps-mt', type=_finite_float, default=0.0, help='Off-diagonal (mu-tau) NSI coupling.')
    g_nsi.add_argument('--eps-tt', type=_finite_float, default=0.0, help='Diagonal NSI coupling of nu_tau.')
    g_nsi.add_argument('--eps-es', type=_finite_float, default=0.0, help='(4nu) Off-diagonal (e-s) NSI coupling.')
    g_nsi.add_argument('--eps-ms', type=_finite_float, default=0.0, help='(4nu) Off-diagonal (mu-s) NSI coupling.')
    g_nsi.add_argument('--eps-ts', type=_finite_float, default=0.0, help='(4nu) Off-diagonal (tau-s) NSI coupling.')
    g_nsi.add_argument('--eps-ss', type=_finite_float, default=0.0, help='(4nu) Diagonal NSI coupling of nu_s.')
    g_nsi.add_argument('--eps-es1', type=_finite_float, default=0.0, help='(5nu) Off-diagonal (e-s1) NSI coupling.')
    g_nsi.add_argument('--eps-es2', type=_finite_float, default=0.0, help='(5nu) Off-diagonal (e-s2) NSI coupling.')
    g_nsi.add_argument('--eps-ms1', type=_finite_float, default=0.0, help='(5nu) Off-diagonal (mu-s1) NSI coupling.')
    g_nsi.add_argument('--eps-ms2', type=_finite_float, default=0.0, help='(5nu) Off-diagonal (mu-s2) NSI coupling.')
    g_nsi.add_argument('--eps-ts1', type=_finite_float, default=0.0, help='(5nu) Off-diagonal (tau-s1) NSI coupling.')
    g_nsi.add_argument('--eps-ts2', type=_finite_float, default=0.0, help='(5nu) Off-diagonal (tau-s2) NSI coupling.')
    g_nsi.add_argument('--eps-s1s1', type=_finite_float, default=0.0, help='(5nu) Diagonal NSI coupling of nu_s1.')
    g_nsi.add_argument('--eps-s1s2', type=_finite_float, default=0.0, help='(5nu) Off-diagonal (s1-s2) NSI coupling.')
    g_nsi.add_argument('--eps-s2s2', type=_finite_float, default=0.0, help='(5nu) Diagonal NSI coupling of nu_s2.')

    g_liv = p.add_argument_group('LIV parameters (--scenario liv)')
    g_liv.add_argument('--sxi', type=_finite_float, default=0.0, help='2-flavor LIV mixing angle xi, per --angles.')
    g_liv.add_argument('--sxi12', type=_finite_float, default=0.0, help='LIV mixing angle xi_12, per --angles.')
    g_liv.add_argument('--sxi23', type=_finite_float, default=0.0, help='LIV mixing angle xi_23, per --angles.')
    g_liv.add_argument('--sxi13', type=_finite_float, default=0.0, help='LIV mixing angle xi_13, per --angles.')
    g_liv.add_argument('--dxicp', type=_finite_float, default=None, dest='dxiCP',
        help='(3/4/5nu) LIV CP-violation phase of the 1-3 rotation [radian, or degree with --angles deg].')
    g_liv.add_argument('--dxi13', type=_finite_float, default=None,
        help='(4/5nu) Deprecated alias of --dxicp; issues a FutureWarning.')
    g_liv.add_argument('--sxi14', type=_finite_float, default=0.0, help='(4/5nu) LIV mixing angle xi_14, per --angles.')
    g_liv.add_argument('--dxi14', type=_finite_float, default=0.0, help='(4/5nu) LIV CP-violation phase of the 1-4 rotation [radian, or degree with --angles deg].')
    g_liv.add_argument('--sxi24', type=_finite_float, default=0.0, help='(4/5nu) LIV mixing angle xi_24, per --angles.')
    g_liv.add_argument('--dxi24', type=_finite_float, default=0.0, help='(4/5nu) LIV CP-violation phase of the 2-4 rotation [radian, or degree with --angles deg].')
    g_liv.add_argument('--sxi34', type=_finite_float, default=0.0, help='(4/5nu) LIV mixing angle xi_34, per --angles.')
    g_liv.add_argument('--sxi15', type=_finite_float, default=0.0, help='(5nu) LIV mixing angle xi_15, per --angles.')
    g_liv.add_argument('--dxi15', type=_finite_float, default=0.0, help='(5nu) LIV CP-violation phase of the 1-5 rotation [radian, or degree with --angles deg].')
    g_liv.add_argument('--sxi25', type=_finite_float, default=0.0, help='(5nu) LIV mixing angle xi_25, per --angles.')
    g_liv.add_argument('--sxi35', type=_finite_float, default=0.0, help='(5nu) LIV mixing angle xi_35, per --angles.')
    g_liv.add_argument('--dxi35', type=_finite_float, default=0.0, help='(5nu) LIV CP-violation phase of the 3-5 rotation [radian, or degree with --angles deg].')
    g_liv.add_argument('--b1', type=_finite_float, default=0.0, help='LIV eigenvalue b1 [eV].')
    g_liv.add_argument('--b2', type=_finite_float, default=0.0, help='LIV eigenvalue b2 [eV].')
    g_liv.add_argument('--b3', type=_finite_float, default=0.0, help='LIV eigenvalue b3 [eV].')
    g_liv.add_argument('--b4', type=_finite_float, default=0.0, help='LIV eigenvalue b4 [eV].')
    g_liv.add_argument('--b5', type=_finite_float, default=0.0, help='LIV eigenvalue b5 [eV].')
    g_liv.add_argument('--liv-lambda', type=_finite_float, default=1.0, dest='Lambda',
        help='LIV energy scale Lambda [eV]. Default: 1.0.')
    g_liv.add_argument('--n-liv', type=int, default=0,
        help='Power of the energy dependence of the LIV operator. Default: 0.')

    g_chan = p.add_argument_group('Channel selection')
    g_chan.add_argument('--nu-i', type=_flavor_index, default=None,
        help='Initial flavor (index or name: e, mu, tau, s, s1, s2). If given with --nu-f, '
             'prints a single probability instead of the full matrix.')
    g_chan.add_argument('--nu-f', type=_flavor_index, default=None,
        help='Final flavor; see --nu-i.')

    g_num = p.add_argument_group('Advanced numerics')
    g_num.add_argument('--magnus-exp-order', type=int, default=4, dest='magnus_exp_order',
        help='Highest order of the Magnus expansion (1-10; an even order from 2 to 8 with '
             'the default --integration-method gl). Default: 4.')
    g_num.add_argument('--integration-method', choices=['gl', 'trapezoid', 'simpson'], default='gl',
        help="Quadrature method. 'gl' (Gauss-Legendre collocation) needs only 1-4 Hamiltonian "
             "evaluations per slab and matches its quadrature order to the expansion order, so "
             "it is both the fastest and the most accurate for a smooth Hamiltonian. "
             "'trapezoid'/'simpson' instead sample a uniform grid of points in each slab, "
             "starting from 100 points per slab (the library default, which the CLI cannot "
             "change). At a declared breakpoint, they sample each side of a jump separately. "
             "Default: gl.")
    g_num.add_argument('--rtol', type=_finite_float, default=1.e-3,
        help='Relative tolerance on the agreement between successive refinement levels: a stopping rule, not a guaranteed accuracy. Default: 1e-3.')
    g_num.add_argument('--atol', type=_finite_float, default=1.e-3,
        help='Absolute tolerance on the same agreement; see --rtol. Default: 1e-3.')
    g_num.add_argument('--n-jobs', type=int, default=1, dest='n_jobs',
        help='Number of parallel joblib workers. Default: 1.')
    g_num.add_argument('--strategy', choices=['auto', 'hybrid', 'magnus'], default='auto',
        help="How to propagate a position-dependent Hamiltonian. 'magnus' uses only the "
             "Magnus-expansion engines. 'hybrid' also tries adiabatic transport with a "
             "Magnus patch at each non-adiabatic window and warns if it cannot certify the "
             "result. 'auto' sends a smooth profile whose estimated accumulated phase is at "
             "most 1e4 rad (at a tolerance of 1e-6 or looser) to the Magnus ladder. "
             "Otherwise, it tries hybrid and falls back to the Magnus engines without a "
             "warning, except for an undeclared density jump. Vacuum and constant-density "
             "environments accept only 'auto'. Default: auto.")
    g_num.add_argument('--verbose', type=int, default=0, choices=[0, 1, 2],
        help='Verbosity level. Default: 0.')

    g_out = p.add_argument_group('Output')
    g_out.add_argument('--json', action='store_true', help='Print the result as JSON instead of a table.')
    g_out.add_argument('--precision', type=_precision, default=4,
        help='Decimal digits shown in the table and in the single-channel value; ignored '
             'with --json. Default: 4.')

    return parser


def _std_osc_kwargs(flavors: int, args: argparse.Namespace) -> dict:
    r"""Collects the standard oscillation parameters for one flavor count.

    Returns only the keywords the ``osc_prob_{flavors}nu_*`` wrappers declare at that
    count: the two-flavor pair at two flavors, otherwise the three-flavor set together
    with ``default_osc_params_set_name``, extended with a fourth state's mixing at four
    or more and a fifth state's at five.  Values go through unconverted, in whichever
    convention ``--angles`` names.

    Parameters
    ----------
    flavors : int
        Number of neutrino flavors (2, 3, 4 or 5).
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    dict
        Keyword arguments for the wrapper, under the library's parameter names.
    """
    if flavors == 2:
        return {'sth': args.sth, 'Dm2': args.Dm2, 'angles': args.angles}
    kw = {'s12': args.s12, 's23': args.s23, 's13': args.s13, 'dCP': args.dCP,
          'D21': args.D21, 'D31': args.D31, 'angles': args.angles,
          'default_osc_params_set_name': args.default_osc_params_set_name}
    if flavors >= 4:
        kw.update({'s14': args.s14, 'd14': args.d14, 's24': args.s24, 'd24': args.d24,
                   's34': args.s34, 'D41': args.D41})
    if flavors == 5:
        kw.update({'s15': args.s15, 'd15': args.d15, 's25': args.s25, 's35': args.s35,
                   'd35': args.d35, 'D51': args.D51})
    return kw


def _nsi_kwargs(flavors: int, args: argparse.Namespace) -> dict:
    r"""Collects the non-standard-interaction couplings for one flavor count.

    Returns the ``eps_*`` keywords the NSI wrappers declare at that count: the
    two-flavor pair at two flavors, otherwise the active-sector block, plus the
    couplings to one sterile state at four flavors and to two at five.  Every coupling
    defaults to 0.0 in the parser, so an unset one is forwarded as zero rather than
    omitted.

    Parameters
    ----------
    flavors : int
        Number of neutrino flavors (2, 3, 4 or 5).
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    dict
        Keyword arguments for the wrapper, under the library's parameter names.
    """
    if flavors == 2:
        return {'eps_aa': args.eps_aa, 'eps_ab': args.eps_ab}
    kw = {'eps_ee': args.eps_ee, 'eps_em': args.eps_em, 'eps_et': args.eps_et,
          'eps_mm': args.eps_mm, 'eps_mt': args.eps_mt, 'eps_tt': args.eps_tt}
    if flavors == 4:
        kw.update({'eps_es': args.eps_es, 'eps_ms': args.eps_ms, 'eps_ts': args.eps_ts,
                   'eps_ss': args.eps_ss})
    if flavors == 5:
        kw.update({'eps_es1': args.eps_es1, 'eps_es2': args.eps_es2, 'eps_ms1': args.eps_ms1,
                   'eps_ms2': args.eps_ms2, 'eps_ts1': args.eps_ts1, 'eps_ts2': args.eps_ts2,
                   'eps_s1s1': args.eps_s1s1, 'eps_s1s2': args.eps_s1s2, 'eps_s2s2': args.eps_s2s2})
    return kw


def _liv_kwargs(flavors: int, args: argparse.Namespace) -> dict:
    r"""Collects the Lorentz-violating operator's parameters for one flavor count.

    Returns the mixing angles, the CP phases and the eigenvalues the LIV wrappers
    declare at that count.  Every count from three up takes ``dxiCP``; four and five add
    the sterile sector's angles, phases and eigenvalues.  ``--dxi13``, the former name of
    ``--dxicp`` at four and five flavors, is still accepted, with a :class:`FutureWarning`.  ``Lambda``
    and ``n_liv`` go in at every count.  An eigenvalue above the flavor count is
    accepted by the parser and dropped here.

    Parameters
    ----------
    flavors : int
        Number of neutrino flavors (2, 3, 4 or 5).
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    dict
        Keyword arguments for the wrapper, under the library's parameter names.
    """
    if flavors == 2:
        return {'sxi': args.sxi, 'b1': args.b1, 'b2': args.b2, 'Lambda': args.Lambda,
                'n_liv': args.n_liv}
    kw = {'sxi12': args.sxi12, 'sxi23': args.sxi23, 'sxi13': args.sxi13,
          'b1': args.b1, 'b2': args.b2, 'b3': args.b3, 'Lambda': args.Lambda, 'n_liv': args.n_liv}
    kw['dxiCP'] = 0.0 if args.dxiCP is None else args.dxiCP
    if flavors >= 4 and args.dxi13 is not None:
        warnings.warn("magnus: '--dxi13' is deprecated; use '--dxicp', which means the same.",
                      FutureWarning, stacklevel=2)
        kw['dxiCP'] = args.dxi13
    if flavors >= 4:
        kw.update({'sxi14': args.sxi14, 'dxi14': args.dxi14,
                   'sxi24': args.sxi24, 'dxi24': args.dxi24, 'sxi34': args.sxi34, 'b4': args.b4})
    if flavors == 5:
        kw.update({'sxi15': args.sxi15, 'dxi15': args.dxi15, 'sxi25': args.sxi25,
                   'sxi35': args.sxi35, 'dxi35': args.dxi35, 'b5': args.b5})
    return kw


def _env_kwargs(environment: str, density_profile: str, args: argparse.Namespace,
                 baseline_ev: Optional[float], l0_ev: float) -> dict:
    r"""Builds the environment's keywords, and refuses combinations it cannot serve.

    Returns what the chosen ``osc_prob_*`` wrapper needs: the baseline alone for vacuum;
    density, composition and unit flag for matter, with either one density or a central
    density and a scale height by profile; the trajectory for earth, with both depths
    converted to natural units; and the two positions and the density profile for sun.
    Every invalid combination of flags exits here rather than downstream, so the message
    names the flags the user typed rather than the library's parameters.

    .. versionchanged:: 1.2.0
       The sun branch forwards the density profile and ``--stop-at-table-edge``, and a
       solar model or ``--stop-at-table-edge`` outside ``--environment sun`` is refused.

    Parameters
    ----------
    environment : str
        'vacuum', 'matter', 'earth' or 'sun'.
    density_profile : str
        One of :data:`DENSITY_PROFILES`, already resolved from its default: 'constant' or
        'exp' for matter, 'exp' or a solar model for sun.  Ignored for vacuum and earth,
        unless it names a solar model.
    args : argparse.Namespace
        Parsed command-line arguments.
    baseline_ev : float or None
        Baseline [eV^-1], or None if ``--baseline`` was not given.
    l0_ev : float
        Initial position [eV^-1].

    Returns
    -------
    dict
        Keyword arguments for the wrapper, under the library's parameter names.

    Raises
    ------
    SystemExit
        If the flags given do not describe a trajectory this environment can
        propagate.  The message names the flags at fault, and the exit code is 1.
    """
    if environment != 'sun':
        if density_profile in solarmodels.SOLAR_MODELS:
            raise SystemExit(f"magnus prob: --density-profile {density_profile} is a solar "
                              f"model, so it applies to --environment sun only.")
        if args.stop_at_table_edge:
            raise SystemExit("magnus prob: --stop-at-table-edge applies to --environment sun "
                              "with a tabulated solar model only.")
    if environment == 'vacuum':
        return {'L': baseline_ev}
    if environment == 'matter':
        kw = {
            'L': baseline_ev,
            'ratio_number_neutrons_to_protons': args.ratio_n_to_p,
            'electron_fraction': args.electron_fraction,
            'density_matter_is_in_g_per_cm3': (args.density_unit == 'g/cm3'),
        }
        if density_profile == 'constant':
            if args.rho is None:
                raise SystemExit("magnus prob: --rho is required for --environment matter "
                                  "--density-profile constant.")
            kw['rho'] = args.rho
        else:
            if args.rho_central is None or args.l_scale is None:
                raise SystemExit("magnus prob: --rho-central and --l-scale are required for "
                                  "--environment matter --density-profile exp.")
            kw['L0'] = l0_ev
            kw['rho_central'] = args.rho_central
            kw['l_scale'] = args.l_scale * LENGTH_UNITS[args.baseline_unit]
        return kw
    if environment == 'earth':
        using_locations = bool(args.loc_ini and args.loc_fin)
        scale = LENGTH_UNITS[args.baseline_unit]
        source_depth = args.source_depth*scale
        detector_depth = args.detector_depth*scale
        buried = bool(source_depth or detector_depth)
        if buried and using_locations:
            raise SystemExit("magnus prob: --loc-ini/--loc-fin fix a surface-to-surface "
                              "chord, so neither --source-depth nor --detector-depth applies "
                              "to them. Use --costhz with the depths instead.")
        if detector_depth and baseline_ev is not None:
            raise SystemExit("magnus prob: --detector-depth says where the trajectory ends "
                              "and so does --baseline. Give one or the other.")
        if not using_locations:
            if args.costhz is None:
                raise SystemExit("magnus prob: --environment earth requires either --costhz "
                                  "(with --baseline, --detector-depth or --source-depth) or "
                                  "both --loc-ini and --loc-fin.")
            if baseline_ev is None and not buried:
                raise SystemExit("magnus prob: --costhz fixes the direction of the chord but "
                                  "not its length: give --baseline, or --detector-depth or "
                                  "--source-depth, from which the length is computed. "
                                  "--loc-ini/--loc-fin compute it without --costhz.")
        return {'costhz': args.costhz, 'loc_ini': args.loc_ini, 'loc_fin': args.loc_fin,
                'L': baseline_ev, 'source_depth': source_depth,
                'detector_depth': detector_depth}
    if environment == 'sun':
        if baseline_ev is None:  # pragma: no cover - pre-empted, see below
            # Unreachable from the command line: main() rejects a missing --baseline for
            # every environment except earth before it calls this function, so a solar run
            # without one has already exited.  Kept because this function is the one that
            # knows what the sun branch needs, and a future caller reaching it by another
            # route should still get a clear message rather than a KeyError downstream.
            raise SystemExit("magnus prob: --baseline is required for --environment sun.")
        if density_profile == 'constant':
            raise SystemExit("magnus prob: --density-profile constant is not available with "
                              "--environment sun; use 'exp' or a solar model.")
        if args.stop_at_table_edge and density_profile == 'exp':
            raise SystemExit("magnus prob: --stop-at-table-edge needs a tabulated solar model "
                              "(--density-profile); the exponential profile has no last row.")
        return {'L': baseline_ev, 'L0': l0_ev, 'density_profile': density_profile,
                'stop_at_table_edge': args.stop_at_table_edge}
    raise AssertionError(environment)  # pragma: no cover


def _wrapper_name(flavors: int, environment: str, scenario: str, density_profile: str) -> str:
    r"""Names the :mod:`magnus.oscprob` wrapper that serves a requested combination.

    Composes ``osc_prob_{flavors}nu_*`` from the environment, the physics scenario and,
    for matter, the density profile.  Non-standard interactions in vacuum is the one
    combination with no wrapper, since the couplings scale a matter potential that
    vacuum does not have, and it exits here rather than failing later as a missing
    attribute.

    Parameters
    ----------
    flavors : int
        Number of neutrino flavors (2, 3, 4 or 5).
    environment : str
        'vacuum', 'matter', 'earth' or 'sun'.
    scenario : str
        'std', 'nsi' or 'liv'.
    density_profile : str
        'constant' or 'exp'; read only when ``environment`` is 'matter'.

    Returns
    -------
    str
        Name of a function in :mod:`magnus.oscprob`.

    Raises
    ------
    SystemExit
        If ``scenario`` is 'nsi' and ``environment`` is 'vacuum'.
    """
    if environment == 'vacuum':
        if scenario == 'nsi':
            raise SystemExit("magnus prob: --scenario nsi is not available with --environment "
                              "vacuum (NSI couplings scale the matter potential, which vacuum "
                              "has none of); use --environment matter/earth/sun instead.")
        suffix = '_liv' if scenario == 'liv' else ''
        return f'osc_prob_{flavors}nu_vacuum{suffix}'
    if environment == 'matter':
        density_suffix = 'constant_density' if density_profile == 'constant' else 'exp_density'
        scenario_infix = {'std': '', 'nsi': '_nsi', 'liv': '_liv'}[scenario]
        return f'osc_prob_{flavors}nu_matter{scenario_infix}_{density_suffix}'
    if environment in ('earth', 'sun'):
        scenario_suffix = {'std': '', 'nsi': '_nsi', 'liv': '_liv'}[scenario]
        return f'osc_prob_{flavors}nu_{environment}{scenario_suffix}'
    raise AssertionError(environment)  # pragma: no cover


def _call(fn, candidate_kwargs: dict):
    r"""Calls fn with only the keys it actually accepts explicitly (plus the
    universally-forwarded refinement/numerics kwargs)."""
    sig = inspect.signature(fn)
    explicit_names = {n for n, par in sig.parameters.items()
                       if par.kind != inspect.Parameter.VAR_KEYWORD}
    kwargs = {k: v for k, v in candidate_kwargs.items()
              if (k in ALWAYS_FORWARD) or (k in explicit_names)}
    return fn(**kwargs)


def _format_table(P: np.ndarray, flavors: int, precision: int) -> str:
    r"""Renders a probability matrix as a labeled fixed-width text table.

    Rows are the initial flavor and columns the final one, both named by
    :data:`FLAVOR_LABELS`; at two flavors the system is abstract, so the labels are 0
    and 1.  Each column is as wide as the longest label plus two, or ``precision + 4``,
    whichever is larger, so the entries stay aligned at any precision.

    Parameters
    ----------
    P : np.ndarray
        Probability matrix, shape ``(flavors, flavors)``.
    flavors : int
        Number of neutrino flavors (2, 3, 4 or 5).
    precision : int
        Decimal digits shown for each entry.

    Returns
    -------
    str
        The table, newline-joined and without a trailing newline.
    """
    labels = FLAVOR_LABELS[flavors]
    width = max(len(l) for l in labels) + 2
    width = max(width, precision + 4)
    header = ' ' * (max(len(l) for l in labels) + 2) + ''.join(
        f'{l:>{width}}' for l in labels)
    lines = [header]
    for i, row_label in enumerate(labels):
        row = ''.join(f'{P[i, j]:>{width}.{precision}f}' for j in range(len(labels)))
        lines.append(f'{row_label:<{max(len(l) for l in labels) + 2}}{row}')
    return '\n'.join(lines)


class _ReadRecorder:
    r"""Wraps an argparse namespace and records which attributes are read."""

    def __init__(self, namespace):
        object.__setattr__(self, '_ns', namespace)
        object.__setattr__(self, 'read', set())

    def __getattr__(self, name):
        self.read.add(name)
        return getattr(self._ns, name)


# Flags that every run uses, or that main() reads for itself rather than through a builder.
_ALWAYS_USED = frozenset(('command', 'flavors', 'environment', 'scenario', 'density_profile',
    'energy', 'energy_unit', 'baseline', 'baseline_unit', 'nubar', 'nu_i', 'nu_f', 'verbose',
    'magnus_exp_order', 'n_jobs', 'integration_method', 'rtol', 'atol', 'json', 'precision',
    'help', 'version'))


def _unread_flags(prob: argparse.ArgumentParser, args: argparse.Namespace,
                  used: _ReadRecorder) -> list:
    r"""The flags set on the command line that no part of this run read, by their option names.

    A flag counts as set when its value differs from its default.

    .. versionadded:: 1.2.0
    """
    out = []
    for action in prob._actions:
        dest = action.dest
        if dest in _ALWAYS_USED or dest in used.read or not action.option_strings:
            continue
        if getattr(args, dest, action.default) != action.default:
            out.append(action.option_strings[-1])
    return out


def main(argv=None) -> int:
    r"""Entry point for the ``magnus`` console script / ``python -m magnus``.

    A command line that names no subcommand runs ``prob``, the only one:
    ``magnus --flavors 3 ...`` is ``magnus prob --flavors 3 ...``.

    .. versionadded:: 1.0.0

    .. versionchanged:: 1.2.0
       ``prob`` is optional (issue #138).

    Parameters
    ----------
    argv : list of str, optional
        Arguments to parse instead of ``sys.argv[1:]`` (mainly for testing).

    Returns
    -------
    int
        Always 0.  Every failure leaves through :class:`SystemExit` instead, so
        this function has no non-zero return.

    Raises
    ------
    SystemExit
        Exit code 2 for an argparse error or a value the library rejects, and 1 for a
        combination of flags that ``_env_kwargs`` or ``_wrapper_name`` refuses.
    """
    parser = build_parser()
    args = parser.parse_args(_with_default_subcommand(sys.argv[1:] if argv is None else argv))
    prob = parser.default_subparser

    # The deprecated name and the new one together: the old name used to win silently, while
    # the Python API refuses the pair (issue #160 §13).
    if args.dxi13 is not None and args.dxiCP is not None:
        prob.error("give --dxicp or its deprecated alias --dxi13, not both.")

    # The builders below read, from `args`, exactly the flags this run uses.  Anything the
    # user set and none of them read would be ignored in silence -- --rho in vacuum, an NSI
    # coupling without --scenario nsi, a sterile angle at three flavors -- so it is refused
    # instead, by name (issue #160 §13).
    used = _ReadRecorder(args)

    flavors = args.flavors
    environment = args.environment
    scenario = args.scenario

    if flavors == 2 and (args.sth is None or args.Dm2 is None):
        parser.error("--sth and --dm2 are both required for --flavors 2.")

    energy_ev = args.energy * ENERGY_UNITS[args.energy_unit]
    baseline_ev = None if args.baseline is None else args.baseline * LENGTH_UNITS[args.baseline_unit]
    l0_ev = args.l0 * LENGTH_UNITS[args.baseline_unit]

    if environment != 'earth' and baseline_ev is None:
        parser.error(f"--baseline is required for --environment {environment}.")

    # The default depends on the environment: a uniform slab for matter, the exponential
    # profile for the Sun (what the sun wrappers have always used).
    density_profile = args.density_profile
    if density_profile is None:
        density_profile = 'exp' if environment == 'sun' else 'constant'

    fn_name = _wrapper_name(flavors, environment, scenario, density_profile)
    fn = getattr(oscprob, fn_name)

    candidate = {'energy': energy_ev}
    candidate.update(_env_kwargs(environment, density_profile, used, baseline_ev, l0_ev))
    candidate.update(_std_osc_kwargs(flavors, used))
    if scenario == 'nsi':
        candidate.update(_nsi_kwargs(flavors, used))
    elif scenario == 'liv':
        candidate.update(_liv_kwargs(flavors, used))
    candidate.update({
        'nubar': args.nubar, 'nu_i': args.nu_i, 'nu_f': args.nu_f,
        'validate_input': True, 'verbose': args.verbose,
        'magnus_exp_order': args.magnus_exp_order, 'n_jobs': args.n_jobs,
        'integration_method': args.integration_method, 'rtol': args.rtol, 'atol': args.atol,
    })
    # `strategy` selects how a *position-dependent* Hamiltonian is propagated, so it is only
    # forwarded where the Hamiltonian actually depends on position.  Vacuum and constant-density
    # environments have no such dependence and their wrappers forward unknown keywords all the
    # way down to the Magnus core, which would reject it.
    if environment in ('earth', 'sun') or (environment == 'matter'
                                           and density_profile == 'exp'):
        candidate['strategy'] = used.strategy
    if 'L0' in candidate:
        used.l0

    ignored = _unread_flags(prob, args, used)
    if ignored:
        prob.error(", ".join(ignored) + (" does" if len(ignored) == 1 else " do") +
                   " not apply to --flavors " + str(flavors) + " --environment " + environment +
                   " --scenario " + scenario + (" --density-profile " + density_profile
                   if environment in ('matter', 'sun') else "") + ", and would be ignored.")

    try:
        P = _call(fn, candidate)
    except ValueError as error:
        # The library validates its own inputs and raises; surface that as a clean CLI error
        # (exit code 2, like any other argument problem) rather than a raw traceback.  The
        # library names its own arguments; the flag the user typed is named instead.
        message = str(error)
        for action in prob._actions:
            if action.option_strings and action.dest not in ('help', 'version'):
                message = message.replace(': ' + action.dest + ' must',
                                          ': ' + action.option_strings[-1] + ' must')
        parser.error(message)

    if args.json:
        payload = {
            'function': fn_name, 'flavors': flavors, 'environment': environment,
            'scenario': scenario, 'nubar': args.nubar,
            'energy_eV': energy_ev, 'baseline_eV-1': baseline_ev,
            'probability': np.asarray(P).tolist(),
        }
        print(json.dumps(payload, indent=2))
        return 0

    print(f"Magνs {__version__} -- {fn_name}")
    label = f"E = {args.energy:g} {args.energy_unit}"
    if baseline_ev is not None:
        label += f", L = {args.baseline:g} {args.baseline_unit}"
    if density_profile in solarmodels.SOLAR_MODELS:
        label += f", {density_profile} solar model"
    if args.nubar:
        label += ", antineutrinos"
    print(label)
    print()
    if args.nu_i is not None and args.nu_f is not None:
        print(f"P = {float(P):.{args.precision}f}")
    else:
        print(_format_table(np.asarray(P), flavors, args.precision))
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
